//! `ReadBgzfBlocks` source step + `read_bam(path)` convenience helper.
//!
//! Reads raw BGZF blocks from a file (no decompression) and emits them as
//! `BgzfBlock` items with monotonically increasing `batch_serial`. The
//! header bytes are NOT skipped here — they pass through as part of the
//! first block(s); `FindBamBoundaries` strips them downstream.

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

use crate::types::BgzfBlock;
use fgumi_pipeline_core::{
    Unpushed,
    held::HeldSlot,
    outputs::OrderedBytesSingle,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// Legacy default blocks-per-batch.
pub const DEFAULT_BLOCKS_PER_BATCH: usize = 16;

/// `Exclusive + sticky` source step that reads raw BGZF blocks from a file.
pub struct ReadBgzfBlocks {
    reader: Arc<Mutex<Option<Box<dyn io::Read + Send>>>>,
    blocks_per_batch: usize,
    next_serial: u64,
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
            sticky: true,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn affinity(&self) -> Affinity {
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

        for raw in raw_blocks {
            let serial = self.next_serial;
            self.next_serial += 1;
            self.pending.push_back(BgzfBlock {
                batch_serial: serial,
                uncompressed_size: u32::try_from(raw.uncompressed_size()).map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "ReadBgzfBlocks: BGZF uncompressed_size out of range: {}",
                            raw.uncompressed_size()
                        ),
                    )
                })?,
                bytes: raw.data,
            });
        }

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
/// `(step, header)` pair.
///
/// The file is deliberately opened twice: first via
/// `create_raw_bam_reader_with_opts` to parse and return the `Header`, then
/// re-opened with `File::open` and a fresh `BufReader` at offset 0 so the raw
/// BGZF stream — including the header blocks — is emitted in full. Those header
/// blocks are stripped downstream by `FindBamBoundaries`. Do not "optimize" this
/// by reusing the first reader's position: skipping the header blocks corrupts
/// the raw-block stream.
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
    let (_, header) = fgumi_bam_io::create_raw_bam_reader_with_opts(path, 1, opts)
        .map_err(|e| io::Error::other(format!("create_raw_bam_reader_with_opts: {e}")))?;

    let file = File::open(path)?;
    let reader: Box<dyn io::Read + Send> =
        Box::new(io::BufReader::with_capacity(2 * 1024 * 1024, file));
    Ok(read_bam_from_reader(reader, header, blocks_per_batch, output_byte_limit))
}

/// Stdin counterpart to [`read_bam`].
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

    #[test]
    fn read_bam_from_reader_round_trips_bytes() {
        const BGZF_EOF_LEN: usize = 28;

        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let header = noodles::sam::Header::default();
        let writer = fgumi_bam_io::create_raw_bam_writer(&path, &header, 1, 1).unwrap();
        writer.finish().unwrap();

        let on_disk = std::fs::read(&path).unwrap();
        assert!(!on_disk.is_empty(), "BAM file should contain header + EOF block");

        let cursor = std::io::Cursor::new(on_disk.clone());
        let reader: Box<dyn io::Read + Send> = Box::new(cursor);
        let (mut step, _hdr) =
            read_bam_from_reader(reader, header, DEFAULT_BLOCKS_PER_BATCH, 1024 * 1024);

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
