//! `ReadBgzfBlocks` source step + `read_bam(path)` convenience helper.
//!
//! Reads raw BGZF blocks from a file (no decompression) and emits them as
//! `BgzfBlock` items with monotonically increasing `batch_serial`. The
//! header bytes are NOT skipped here — they pass through as part of the
//! first block(s); `FindBamBoundaries` strips them downstream. Same model
//! as the legacy `pipeline/bam.rs` block-reading thread.
//!
//! `Serial` + `Affinity::Reader`. Worker 0 is the only worker that ever
//! attempts the source's mutex; other workers `Skip` this step in
//! dispatch — no `try_lock` thrash. Mirrors legacy's "sticky read on
//! T0" pattern (`pipeline/bam.rs:3541` + `base.rs:4360-4392`).
//! Worker 0 still round-robins through every other step in the chain;
//! the priority-restart logic in the driver keeps it returning to the
//! source whenever the source has more data and downstream isn't
//! backpressured.

use std::collections::VecDeque;
use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use fgumi_bam_io::PipelineReaderOpts;
use fgumi_bgzf::reader::read_raw_blocks;
use noodles::sam::Header;
use parking_lot::Mutex;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::BgzfBlock;

/// Legacy default blocks-per-batch (kept for compatibility with the
/// original `read_bam(..., DEFAULT_BLOCKS_PER_BATCH, ...)` callers).
/// **New code should use [`crate::pipeline::steps::tuning::BamPipelineTuning::auto_tuned`]**
/// which scales the value with thread count (16/32/48/64) to match the
/// legacy pipeline's `auto_tuned` defaults.
pub const DEFAULT_BLOCKS_PER_BATCH: usize = 16;

/// `Exclusive + sticky` source step that reads raw BGZF blocks from a
/// file. Each `try_run` call reads up to `blocks_per_batch` blocks and
/// emits them one at a time as `BgzfBlock` items.
pub struct ReadBgzfBlocks {
    /// Held inside `Mutex<Option<...>>` because the `Step` trait requires
    /// `Clone` (which the runtime never invokes for `Exclusive` steps but
    /// the bound is uniform). Clone panics; the worker that owns the
    /// step holds the `&mut self` and never contends.
    reader: Arc<Mutex<Option<Box<dyn io::Read + Send>>>>,
    blocks_per_batch: usize,
    next_serial: u64,
    /// Pending block batch when emission is mid-flight (held across
    /// retries until each block is pushed). FIFO so that consecutive
    /// `pop_front` calls preserve the read-order serials assigned in
    /// `try_run`.
    pending: VecDeque<BgzfBlock>,
    held: HeldSlot<Unpushed<BgzfBlock>>,
    output_byte_limit: u64,
    finished: Arc<AtomicBool>,
}

impl ReadBgzfBlocks {
    #[must_use]
    pub fn new(
        reader: Box<dyn io::Read + Send>,
        blocks_per_batch: usize,
        output_byte_limit: u64,
    ) -> Self {
        Self {
            reader: Arc::new(Mutex::new(Some(reader))),
            blocks_per_batch: blocks_per_batch.max(1),
            next_serial: 0,
            pending: VecDeque::new(),
            held: HeldSlot::new(),
            output_byte_limit,
            finished: Arc::new(AtomicBool::new(false)),
        }
    }
}

impl Step for ReadBgzfBlocks {
    type Input = ();
    type Outputs = OrderedBytesSingle<BgzfBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ReadBgzfBlocks",
            kind: StepKind::Serial,
            // Worker 0 (the Affinity::Reader target) drives the source
            // sticky — fills the downstream queue in tight bursts before
            // yielding to round-robin. Mirrors legacy `bam.rs:3541` +
            // `base.rs:4360-4392` sticky read.
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn affinity(&self) -> Affinity {
        // Worker 0 is the dedicated reader; other workers Skip this step.
        // Matches legacy `bam.rs:3541` (`is_reader: thread_id == 0`).
        Affinity::Reader
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // 1. Drain the held slot first.
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // 2. Drain pending blocks (one per call iteration).
        if let Some(block) = self.pending.pop_front() {
            match ctx.outputs.push(block) {
                Ok(()) => return Ok(StepOutcome::Progress),
                Err(unpushed) => {
                    self.held.put(unpushed);
                    return Ok(StepOutcome::Progress);
                }
            }
        }

        if self.finished.load(Ordering::Acquire) {
            return Ok(StepOutcome::Finished);
        }

        // 3. Read up to `blocks_per_batch` raw BGZF blocks.
        let raw_blocks = {
            let mut guard = self.reader.lock();
            let reader =
                guard.as_mut().expect("ReadBgzfBlocks: reader missing — was clone() called?");
            read_raw_blocks(reader.as_mut(), self.blocks_per_batch)?
        };

        if raw_blocks.is_empty() {
            self.finished.store(true, Ordering::Release);
            return Ok(StepOutcome::Finished);
        }

        // Assign serials in read order; `pop_front` then yields them in
        // the same order without any reversal trick.
        for raw in raw_blocks {
            let serial = self.next_serial;
            self.next_serial += 1;
            self.pending.push_back(BgzfBlock {
                batch_serial: serial,
                uncompressed_size: u32::try_from(raw.uncompressed_size()).unwrap_or(u32::MAX),
                bytes: raw.data,
            });
        }

        // Emit one block from the freshly-filled queue this iteration so
        // we don't return `NoProgress` while having work in hand.
        if let Some(block) = self.pending.pop_front() {
            match ctx.outputs.push(block) {
                Ok(()) => Ok(StepOutcome::Progress),
                Err(unpushed) => {
                    self.held.put(unpushed);
                    Ok(StepOutcome::Progress)
                }
            }
        } else {
            Ok(StepOutcome::NoProgress)
        }
    }
}

/// Build a [`ReadBgzfBlocks`] step from an already-prepared reader + header.
///
/// The reader MUST be positioned at byte 0 of the BAM stream (i.e. include
/// the BAM header bytes). `FindBamBoundaries::new()` strips them downstream.
/// Used by both [`read_bam`] (file path) and [`read_bam_stdin`] (stdin via
/// `create_bam_reader_for_pipeline_with_opts`'s tee-replay reader).
#[must_use]
pub fn read_bam_from_reader(
    reader: Box<dyn io::Read + Send>,
    header: Header,
    blocks_per_batch: usize,
    output_byte_limit: u64,
) -> (ReadBgzfBlocks, Header) {
    (ReadBgzfBlocks::new(reader, blocks_per_batch, output_byte_limit), header)
}

/// Convenience helper: open a BAM file, parse its header, return the
/// `(step, header)` pair. Header parsing reads the first BGZF block(s)
/// once, then the file is reopened from byte 0 for the source's raw-block
/// reads. The header bytes pass through the chain a second time (the
/// `FindBamBoundaries` step strips them); the double-read of header
/// bytes is small overhead since headers are typically a few KB.
///
/// **Limitation (Phase 3):** this helper requires a seekable file —
/// stdin/pipe inputs use [`read_bam_stdin`] instead, which buffers the
/// header bytes via a `TeeReader`/`ChainedReader` pair so the source's
/// raw-block reads start at byte 0 of the replayed stream.
/// A Phase 4 follow-up will share the post-header reader between header
/// parsing and raw-block reads (per the design doc's read-once contract).
///
/// # Errors
///
/// Returns I/O errors from file open or BAM-header parse.
pub fn read_bam<P: AsRef<Path>>(
    path: P,
    opts: PipelineReaderOpts,
    blocks_per_batch: usize,
    output_byte_limit: u64,
) -> io::Result<(ReadBgzfBlocks, Header)> {
    let path = path.as_ref();
    // Parse header via the high-level reader (consumes the first BGZF
    // block(s) but the resulting reader is dropped right after).
    let (_, header) = fgumi_bam_io::create_raw_bam_reader_with_opts(path, 1, opts)
        .map_err(|e| io::Error::other(format!("create_raw_bam_reader_with_opts: {e}")))?;

    // Reopen from byte 0 for the source's raw-block reads. Wrap in a
    // 2 MiB `BufReader` to amortize disk-I/O syscalls — every BGZF
    // block read involves several short `read(2)` calls (header + extra
    // fields + body); a raw `File` would syscall on each. Matches
    // the legacy spill-merger's `BufReader::with_capacity(2 * 1024 *
    // 1024, file)` (`crates/fgumi-sort/src/external.rs:552`).
    let file = File::open(path)?;
    let reader: Box<dyn io::Read + Send> =
        Box::new(io::BufReader::with_capacity(2 * 1024 * 1024, file));
    Ok(read_bam_from_reader(reader, header, blocks_per_batch, output_byte_limit))
}

/// Stdin counterpart to [`read_bam`]. Uses
/// [`fgumi_bam_io::create_bam_reader_for_pipeline_with_opts`] to parse the
/// BAM header off stdin via a `TeeReader`, then replays the buffered header
/// bytes ahead of the remaining stdin stream so the source emits BGZF
/// blocks starting at byte 0 (same shape as the file path).
///
/// # Errors
///
/// Returns I/O errors from stdin read or BAM-header parse.
pub fn read_bam_stdin(
    opts: PipelineReaderOpts,
    blocks_per_batch: usize,
    output_byte_limit: u64,
) -> io::Result<(ReadBgzfBlocks, Header)> {
    let (reader, header) =
        fgumi_bam_io::create_bam_reader_for_pipeline_with_opts(Path::new("-"), opts).map_err(
            |e| io::Error::other(format!("create_bam_reader_for_pipeline_with_opts: {e}")),
        )?;
    Ok(read_bam_from_reader(reader, header, blocks_per_batch, output_byte_limit))
}

/// Path-aware dispatcher: routes to [`read_bam_stdin`] when `path` is a
/// stdin sentinel (`-` or `/dev/stdin`) and to [`read_bam`] otherwise.
///
/// Used by every command's typed-step pipeline as the single entry point
/// for opening a BAM source. Replaces the legacy "reject stdin on the new
/// path" branches with a uniform dispatch.
///
/// # Errors
///
/// Returns I/O errors from file open, stdin read, or BAM-header parse.
pub fn read_bam_auto<P: AsRef<Path>>(
    path: P,
    opts: PipelineReaderOpts,
    blocks_per_batch: usize,
    output_byte_limit: u64,
) -> io::Result<(ReadBgzfBlocks, Header)> {
    if fgumi_bam_io::is_stdin_path(path.as_ref()) {
        read_bam_stdin(opts, blocks_per_batch, output_byte_limit)
    } else {
        read_bam(path, opts, blocks_per_batch, output_byte_limit)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn profile_advertises_serial_reader_byordinal() {
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let header = noodles::sam::Header::default();
        let writer = fgumi_bam_io::create_raw_bam_writer(&path, &header, 1, 1).unwrap();
        writer.finish().unwrap();

        let (step, _hdr) =
            read_bam(&path, PipelineReaderOpts::default(), DEFAULT_BLOCKS_PER_BATCH, 1024 * 1024)
                .unwrap();
        let profile = step.profile();
        assert_eq!(profile.name, "ReadBgzfBlocks");
        assert_eq!(profile.kind, StepKind::Serial);
        assert!(profile.sticky);
        assert_eq!(step.affinity(), Affinity::Reader);
        assert_eq!(profile.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
        assert!(matches!(profile.output_queues[0], QueueSpec::ByteBounded { .. }));
    }

    /// Round-trip a known BAM through [`read_bam_from_reader`] — the same
    /// reader-from-`Box<dyn Read + Send>` shape that [`read_bam_stdin`]
    /// produces. Drives the step to completion against the in-memory
    /// reader and asserts that the emitted BGZF blocks reconstruct the
    /// original BAM byte-for-byte (sans the BGZF EOF marker, which
    /// `read_raw_blocks` filters out).
    #[test]
    fn read_bam_from_reader_round_trips_bytes() {
        const BGZF_EOF_LEN: usize = 28;

        // Write a non-empty BAM file. `create_raw_bam_writer` emits the
        // BAM header + BGZF EOF marker.
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let header = noodles::sam::Header::default();
        let writer = fgumi_bam_io::create_raw_bam_writer(&path, &header, 1, 1).unwrap();
        writer.finish().unwrap();

        let on_disk = std::fs::read(&path).unwrap();
        assert!(!on_disk.is_empty(), "BAM file should contain header + EOF block");

        // Build the source from a `Cursor<Vec<u8>>` reader (stand-in for
        // the chained-stdin reader shape) and drain the held queue.
        let cursor = std::io::Cursor::new(on_disk.clone());
        let reader: Box<dyn io::Read + Send> = Box::new(cursor);
        let (mut step, _hdr) =
            read_bam_from_reader(reader, header, DEFAULT_BLOCKS_PER_BATCH, 1024 * 1024);

        // Drain the source via the same `read_raw_blocks` call its
        // `try_run` uses. The concatenated block bytes equal everything
        // in the file except the trailing 28-byte EOF block.
        let mut collected = Vec::with_capacity(on_disk.len());
        let mut last_serial: Option<u64> = None;
        loop {
            let raw_blocks = {
                let mut guard = step.reader.lock();
                let reader = guard.as_mut().expect("reader present");
                fgumi_bgzf::reader::read_raw_blocks(reader.as_mut(), DEFAULT_BLOCKS_PER_BATCH)
                    .unwrap()
            };
            if raw_blocks.is_empty() {
                break;
            }
            for raw in raw_blocks {
                let serial = step.next_serial;
                step.next_serial += 1;
                if let Some(prev) = last_serial {
                    assert_eq!(serial, prev + 1, "serials must be monotonic");
                }
                last_serial = Some(serial);
                collected.extend_from_slice(&raw.data);
            }
        }
        let expected = &on_disk[..on_disk.len() - BGZF_EOF_LEN];
        assert_eq!(collected, expected, "concatenated blocks must equal source bytes minus EOF");
        assert!(last_serial.is_some(), "should have read at least one block");
    }
}
