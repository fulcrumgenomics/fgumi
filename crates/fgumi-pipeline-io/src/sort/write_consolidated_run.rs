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
//! block; on the last it writes the trailer and, for an indexed run, the
//! run-index footer, flushes, and reports where the run's frames end through
//! [`ConsolidationRegistry::complete`], which unparks `SpillWrite`.

use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

use fgumi_pipeline_core::external_wait;
use fgumi_pipeline_core::step::{DetachedGroup, Step, StepCtx, StepKind, StepOutcome, StepProfile};
use fgumi_sort::{RunIndex, SpillCodec, spill_magic, spill_trailer};

use super::consolidation::ConsolidationRegistry;
use super::protocol::SpillBlockEvent;
use super::spill_block_compress::SpillBufPool;
use super::spill_write::{fold_block_frames, write_index_footer};

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
    /// Bytes written so far, the codec magic included: the file position the
    /// next block lands at.
    bytes_written: u64,
    /// The run's index so far, at its job's stride: the first block's frames,
    /// extended by each later block's; `None` for good once a block carried
    /// no frames (as `SpillWrite` does for a phase-1 run).
    index: Option<RunIndex>,
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
        let SpillBlockEvent::Block {
            file_id: job_id,
            is_last_in_file,
            bytes,
            frames,
            index_stride,
            ..
        } = event
        else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "WriteConsolidatedRun received a non-block event",
            ));
        };
        let mut first_block = false;
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
            let bytes_written = spill_magic(self.codec).len() as u64;
            self.current = Some(OpenRun { job_id, path, writer, bytes_written, index: None });
            first_block = true;
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
        // The merged run is indexed at the stride its job was stamped with
        // (the live runs' stride when it started), folded block by block, so
        // its index is never held denser than that while the job runs.
        fold_block_frames(&mut open.index, first_block, frames, open.bytes_written, index_stride);
        external_wait(|| open.writer.write_all(&bytes))?;
        open.bytes_written += bytes.len() as u64;
        let written = bytes.len() as u64;
        if let Some(pool) = &self.buf_pool {
            pool.give(bytes);
        }
        if is_last_in_file {
            let OpenRun { job_id, path, mut writer, bytes_written, index } =
                self.current.take().expect("open run");
            // The merged run's index, at its job's stride, goes into its
            // footer, after the frames; `SpillWrite` reads it back at
            // `AllAnnounced`, brought to the final stride.
            let frame_bytes = bytes_written + spill_trailer(self.codec).len() as u64;
            // Dropping the writer closes the file, which can block on a parallel
            // filesystem (NFS, Lustre) as long as the writes, so it is bracketed too.
            let index_footer = external_wait(|| {
                writer.write_all(spill_trailer(self.codec))?;
                let index_footer = write_index_footer(&mut writer, index.as_ref())?;
                writer.flush()?;
                drop(writer);
                Ok::<_, io::Error>(index_footer)
            })
            .map_err(|e| {
                io::Error::new(
                    e.kind(),
                    format!("failed to write consolidated run {}: {e}", path.display()),
                )
            })?;
            self.registry.complete(job_id, frame_bytes, index_footer);
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

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    /// A test block: its length and its frames' `(comp_len, raw_len)` spans,
    /// or none.
    type TestBlock = (usize, Option<Vec<(u32, u32)>>);

    /// A block's index: one frame per `(comp_len, raw_len)` (comp lengths are
    /// the bytes the frame spans in the block), offsets from the block's start.
    fn frames(spans: &[(u32, u32)]) -> RunIndex {
        let mut index = RunIndex::new(8);
        for (i, &(comp, raw)) in spans.iter().enumerate() {
            index.push_frame(comp, raw, 1, &(i as u64).to_le_bytes());
        }
        index
    }

    /// Feed job 0's blocks (each `bytes` long, with `frames` or none, stamped
    /// `index_stride`) through `WriteConsolidatedRun` and return the index
    /// read back from the merged run's footer, after checking the completion
    /// names the run's frame bytes and whether it has a footer.
    fn write_job(blocks: Vec<(usize, Option<RunIndex>)>, index_stride: u32) -> Option<RunIndex> {
        let dir = tempfile::tempdir().unwrap();
        let registry = ConsolidationRegistry::new();
        let path = dir.path().join("merged.zsp");
        registry.begin(0, path.clone());
        let body_end = spill_magic(SpillCodec::Zstd).len() as u64
            + blocks.iter().map(|(len, _)| *len as u64).sum::<u64>();
        let mut writer = WriteConsolidatedRun::new(SpillCodec::Zstd, Arc::clone(&registry));
        let n = blocks.len();
        for (i, (len, frames)) in blocks.into_iter().enumerate() {
            writer
                .write(SpillBlockEvent::Block {
                    ordinal: i as u64,
                    file_id: 0,
                    key_kind: fgumi_sort::SpillKeyKind::Coordinate,
                    is_last_in_file: i + 1 == n,
                    records_ingested_so_far: 0,
                    bytes: vec![0u8; len],
                    frames,
                    index_stride,
                })
                .unwrap();
        }
        let (job_id, frame_bytes, index_footer) =
            registry.take_completed().expect("the job completed");
        assert_eq!((job_id, frame_bytes), (0, body_end));
        let slot =
            fgumi_sort::open_spill_slot(&path, 0, fgumi_sort::FramesEnd::At(frame_bytes)).unwrap();
        assert_eq!(slot.len(), body_end, "the slot's frames exclude the footer");
        assert_eq!(slot.has_index_footer(), index_footer);
        slot.read_run_index().unwrap()
    }

    /// The index rule both run writers share (`fold_block_frames`): every
    /// block indexed gives one index over the whole run, its offsets rebased
    /// past the codec magic and the earlier blocks, written as the run's
    /// footer; a block without frames drops the run's index for good, and the
    /// run has no footer. Columns: blocks (length, frame spans or none),
    /// expected entry `(offset, end)`s read back from the footer (`None`: no
    /// footer).
    #[rstest]
    #[case::frames_frames(
        vec![(30, Some(vec![(10, 100), (20, 200)])), (15, Some(vec![(15, 150)]))],
        Some(vec![(0, 10), (10, 30), (30, 45)])
    )]
    #[case::frames_frames_none(
        vec![(10, Some(vec![(10, 100)])), (5, Some(vec![(5, 50)])), (7, None)],
        None
    )]
    #[case::none_then_frames(vec![(10, None), (5, Some(vec![(5, 50)]))], None)]
    #[case::single_block(vec![(12, Some(vec![(4, 40), (8, 80)]))], Some(vec![(0, 4), (4, 12)]))]
    fn consolidated_run_index_follows_its_blocks(
        #[case] blocks: Vec<TestBlock>,
        #[case] expected: Option<Vec<(u64, u64)>>,
    ) {
        let magic = spill_magic(SpillCodec::Zstd).len() as u64;
        let blocks = blocks.into_iter().map(|(len, f)| (len, f.as_deref().map(frames))).collect();
        let index = write_job(blocks, 1);
        let got = index.map(|ix| (0..ix.len()).map(|i| (ix.offset(i), ix.end(i))).collect());
        let expected: Option<Vec<(u64, u64)>> =
            expected.map(|e| e.into_iter().map(|(o, end)| (o + magic, end + magic)).collect());
        assert_eq!(got, expected);
    }

    /// The merged run's footer is at the stride its job's blocks were
    /// stamped with, not dense: its entries group the frames, a short entry
    /// carried across blocks. Columns: the stamped stride, expected entry
    /// `(offset, end, raw)`s read back from the footer.
    #[rstest]
    #[case::dense(1, vec![(0, 10, 100), (10, 30, 200), (30, 45, 150), (45, 50, 50)])]
    #[case::stride_two(2, vec![(0, 30, 300), (30, 50, 200)])]
    #[case::stride_four(4, vec![(0, 50, 500)])]
    fn consolidated_run_footer_is_at_the_stamped_stride(
        #[case] stride: u32,
        #[case] expected: Vec<(u64, u64, u32)>,
    ) {
        let magic = spill_magic(SpillCodec::Zstd).len() as u64;
        let blocks = vec![
            (30, Some(frames(&[(10, 100), (20, 200)]))),
            (20, Some(frames(&[(15, 150), (5, 50)]))),
        ];
        let index = write_job(blocks, stride).expect("an indexed run");
        assert_eq!(index.stride(), stride);
        assert_eq!(index.covered_frames(), 4);
        let got: Vec<(u64, u64, u32)> =
            (0..index.len()).map(|i| (index.offset(i), index.end(i), index.raw_len(i))).collect();
        let expected: Vec<(u64, u64, u32)> =
            expected.into_iter().map(|(o, end, raw)| (o + magic, end + magic, raw)).collect();
        assert_eq!(got, expected);
    }
}
