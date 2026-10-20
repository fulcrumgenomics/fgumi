//! `WriteConsolidatedRun` — the sink of the consolidation tail: writes each
//! job's compressed blocks to the job's merged run and reports the run to
//! `SpillWrite`.
//!
//! A `Detached` sink on the `sort-io` driver. It receives
//! `SpillBlockCompressConsolidation`'s blocks in ordinal order (the
//! compressor's ordered output), so one job's blocks arrive contiguously and
//! one file is open at a time. On a job's first block it creates the output the
//! job started with ([`ConsolidationRegistry::output_of`]) — `create_new`, so
//! it never truncates a live run — and writes the codec magic; it appends every
//! block; on the last it writes the trailer, flushes, and reports the run's
//! size through [`ConsolidationRegistry::complete`], which unparks
//! `SpillWrite`.

use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

use fgumi_pipeline_core::external_wait;
use fgumi_pipeline_core::step::{DetachedGroup, Step, StepCtx, StepKind, StepOutcome, StepProfile};
use fgumi_sort::{SpillCodec, spill_magic, spill_trailer};

use super::consolidation::ConsolidationRegistry;
use super::protocol::SpillBlockEvent;
use super::spill_block_compress::SpillBufPool;

/// Blocks written per `try_run`, whichever of this and [`BYTES_PER_CALL`] is
/// reached first (as `SpillWrite` batches its writes).
const BLOCKS_PER_CALL: usize = 64;

/// Compressed bytes written per `try_run`.
const BYTES_PER_CALL: u64 = 8 * 1024 * 1024;

/// The merged run being written: its job, path and writer.
struct OpenRun {
    job_id: u32,
    path: PathBuf,
    writer: BufWriter<File>,
}

/// `Detached` sink writing consolidated runs (see the module docs).
pub struct WriteConsolidatedRun {
    current: Option<OpenRun>,
    codec: SpillCodec,
    registry: Arc<ConsolidationRegistry>,
    /// Where written blocks go back for reuse (the tail's shared pool).
    buf_pool: Option<Arc<SpillBufPool>>,
}

impl WriteConsolidatedRun {
    /// A writer of `codec` runs reporting completions to `registry`.
    #[must_use]
    pub fn new(codec: SpillCodec, registry: Arc<ConsolidationRegistry>) -> Self {
        Self { current: None, codec, registry, buf_pool: None }
    }

    /// Return each written block's buffer to `pool` (the consolidation tail's
    /// shared pool, which its compressor draws output buffers from).
    #[must_use]
    pub fn with_buf_pool(mut self, pool: Arc<SpillBufPool>) -> Self {
        self.buf_pool = Some(pool);
        self
    }

    /// Write one event. Returns the compressed bytes written.
    fn write(&mut self, event: SpillBlockEvent) -> io::Result<u64> {
        let SpillBlockEvent::Block { file_id: job_id, is_last_in_file, bytes, .. } = event else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "WriteConsolidatedRun received a non-block event",
            ));
        };
        if self.current.is_none() {
            let path = self.registry.output_of(job_id).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "WriteConsolidatedRun: no output recorded for consolidation job {job_id}"
                    ),
                )
            })?;
            // Run-file I/O blocks on the filesystem — minutes on a busy parallel
            // filesystem — which is an external wait, not a wedge (#1068).
            let writer = external_wait(|| {
                let file =
                    OpenOptions::new().write(true).create_new(true).open(&path).map_err(|e| {
                        io::Error::new(
                            e.kind(),
                            format!("failed to create consolidated run {}: {e}", path.display()),
                        )
                    })?;
                let mut writer = BufWriter::with_capacity(256 * 1024, file);
                writer.write_all(spill_magic(self.codec))?;
                Ok::<_, io::Error>(writer)
            })?;
            self.current = Some(OpenRun { job_id, path, writer });
        }
        let open = self.current.as_mut().expect("opened above");
        if open.job_id != job_id {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "WriteConsolidatedRun: a block of job {job_id} arrived while job {} was open",
                    open.job_id
                ),
            ));
        }
        external_wait(|| open.writer.write_all(&bytes))?;
        let written = bytes.len() as u64;
        if let Some(pool) = &self.buf_pool {
            pool.give(bytes);
        }
        if is_last_in_file {
            let OpenRun { job_id, path, mut writer } = self.current.take().expect("open run");
            // Dropping the writer closes the file, which can block on a parallel
            // filesystem (NFS, Lustre) as long as the writes, so it is bracketed too.
            external_wait(|| {
                writer.write_all(spill_trailer(self.codec))?;
                writer.flush()?;
                drop(writer);
                Ok::<_, io::Error>(())
            })?;
            let len = external_wait(|| std::fs::metadata(&path))?.len();
            self.registry.complete(job_id, len);
        }
        Ok(written)
    }
}

impl Step for WriteConsolidatedRun {
    type Input = SpillBlockEvent;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "WriteConsolidatedRun",
            kind: StepKind::Detached,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn detached_group(&self) -> DetachedGroup {
        DetachedGroup::Shared(crate::sort::SORT_IO_GROUP)
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let (mut blocks, mut bytes) = (0usize, 0u64);
        while blocks < BLOCKS_PER_CALL && bytes < BYTES_PER_CALL {
            let Some(event) = ctx.input.pop() else { break };
            bytes += self.write(event)?;
            blocks += 1;
        }
        if blocks > 0 {
            return Ok(StepOutcome::Progress);
        }
        if ctx.input.is_drained() {
            if let Some(open) = &self.current {
                return Err(io::Error::other(format!(
                    "WriteConsolidatedRun: input drained with consolidated run {} unfinished",
                    open.path.display()
                )));
            }
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}
