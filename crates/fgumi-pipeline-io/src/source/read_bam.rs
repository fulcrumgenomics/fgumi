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

use fgumi_bam_io::PipelineReaderOpts;
use fgumi_bgzf::reader::read_raw_blocks;
use noodles::sam::Header;

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

/// `Serial + sticky` source step that reads raw BGZF blocks from a file.
///
/// The reader and the finished flag are plain owned fields, not `Arc`/atomics:
/// this is a `Serial` step, so the runtime drives a single shared instance and
/// never calls `new_worker_copy` on it (only `Parallel` steps are cloned per
/// worker). There is no second owner to share them with.
pub struct ReadBgzfBlocks {
    reader: Option<Box<dyn io::Read + Send>>,
    blocks_per_batch: usize,
    next_serial: u64,
    pending: VecDeque<BgzfBlock>,
    held: HeldSlot<Unpushed<BgzfBlock>>,
    output_byte_limit: u64,
    finished: bool,
}

impl ReadBgzfBlocks {
    #[must_use]
    pub fn new(
        reader: Box<dyn io::Read + Send>,
        blocks_per_batch: usize,
        output_byte_limit: u64,
    ) -> Self {
        Self {
            reader: Some(reader),
            blocks_per_batch: blocks_per_batch.max(1),
            next_serial: 0,
            pending: VecDeque::new(),
            held: HeldSlot::new(),
            output_byte_limit,
            finished: false,
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

        if self.finished {
            return Ok(StepOutcome::Finished);
        }

        // 3. Read up to `blocks_per_batch` raw BGZF blocks. The reader is taken
        // only at end of stream, so `None` here means `try_run` was called again
        // after it already returned `Finished`.
        let raw_blocks = {
            let reader = self
                .reader
                .as_mut()
                .expect("ReadBgzfBlocks: try_run called after the source reported Finished");
            read_raw_blocks(reader.as_mut(), self.blocks_per_batch)?
        };

        if raw_blocks.is_empty() {
            self.finished = true;
            // Release the reader (and its 2 MiB BufReader) as soon as the stream
            // is drained rather than holding it for the rest of the run.
            self.reader = None;
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
                index: None,
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
    let reader = build_source_reader(file, path, &opts)?;
    Ok(read_bam_from_reader(reader, header, blocks_per_batch, output_byte_limit))
}

/// The plain sequential / async-prefetch reader (`--read-streams 1` and every
/// non-sort command). `verify_crc` does not apply here — `ReadBgzfBlocks`
/// forwards compressed bytes without decoding them.
fn sequential_or_async_reader(
    file: File,
    path: &Path,
    async_reader: bool,
) -> Box<dyn io::Read + Send> {
    if async_reader {
        log::info!("async read enabled: spawning fgumi-prefetch thread for {}", path.display());
        Box::new(fgumi_bam_io::prefetch_reader::PrefetchReader::from_file(file))
    } else {
        Box::new(io::BufReader::with_capacity(2 * 1024 * 1024, file))
    }
}

/// Pick the reader for a regular-file source: the concurrent
/// [`fgumi_bam_io::scatter_reader::ScatterReader`] when `--read-streams` asks
/// for more than one stream and the file is seekable, otherwise the sequential
/// / async reader. The scatter-vs-fallback decision (and its warning) is shared
/// with the sort chain's reader via [`fgumi_bam_io::scatter_reader::decide_reader`]
/// so the two copies cannot drift; only the fallback reader differs (a 2 MiB
/// `BufReader` here vs a bare `File` there).
fn build_source_reader(
    file: File,
    path: &Path,
    opts: &PipelineReaderOpts,
) -> io::Result<Box<dyn io::Read + Send>> {
    match fgumi_bam_io::scatter_reader::decide_reader(file, opts.read_streams, path, "reader")? {
        fgumi_bam_io::scatter_reader::ScatterDecision::Scatter(scatter) => Ok(scatter),
        fgumi_bam_io::scatter_reader::ScatterDecision::Fallback(file) => {
            Ok(sequential_or_async_reader(file, path, opts.async_reader))
        }
    }
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
    // stdin is not seekable, so concurrent positional reads cannot apply. Don't
    // fail the run over a perf knob — warn (only on an explicit `Fixed(n>1)`;
    // the `Auto` default falls back silently) and read sequentially.
    fgumi_bam_io::scatter_reader::warn_read_streams_unavailable(
        opts.read_streams,
        "stdin",
        "is not a seekable regular file",
    );
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
    use rstest::rstest;

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

    // ---------------------------------------------------------------------
    // Driving the step through a real pipeline
    // ---------------------------------------------------------------------

    /// Sink that records every `BgzfBlock` the source emits, in arrival order.
    struct BlockSink {
        received: std::sync::Arc<parking_lot::Mutex<Vec<BgzfBlock>>>,
    }

    impl Step for BlockSink {
        type Input = BgzfBlock;
        type Outputs = ();

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "BlockSink",
                kind: StepKind::Serial,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(block) => {
                    self.received.lock().push(block);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// Write `record_count` records to a temp BAM and return its `path` plus the
    /// on-disk bytes.
    fn temp_bam(record_count: usize) -> (tempfile::TempPath, Vec<u8>) {
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let header = noodles::sam::Header::default();
        let mut writer = fgumi_bam_io::create_raw_bam_writer(&path, &header, 1, 1).unwrap();
        for i in 0..record_count {
            let name = format!("q{i}");
            let bytes = fgumi_raw_bam::testutil::make_bam_bytes(
                0,
                i32::try_from(i).unwrap(),
                0,
                name.as_bytes(),
                &[],
                10,
                -1,
                -1,
                &[],
            );
            writer.write_raw_record(&bytes).unwrap();
        }
        writer.finish().unwrap();
        let on_disk = std::fs::read(&path).unwrap();
        (path, on_disk)
    }

    /// Run a `ReadBgzfBlocks -> BlockSink` pipeline and return the emitted blocks.
    fn drive(path: &Path, blocks_per_batch: usize, threads: usize) -> Vec<BgzfBlock> {
        drive_with_opts(path, blocks_per_batch, threads, PipelineReaderOpts::default())
    }

    /// `drive` with explicit reader opts (used to exercise `--read-streams`).
    fn drive_with_opts(
        path: &Path,
        blocks_per_batch: usize,
        threads: usize,
        opts: PipelineReaderOpts,
    ) -> Vec<BgzfBlock> {
        let received = std::sync::Arc::new(parking_lot::Mutex::new(Vec::new()));
        let (source, _hdr) = read_bam(path, opts, blocks_per_batch, 1024 * 1024).unwrap();
        let sink = BlockSink { received: std::sync::Arc::clone(&received) };

        let builder = fgumi_pipeline_core::builder::Pipeline::builder();
        builder.chain(source).chain(sink).into_sink_marker();
        let pipeline = builder.build().unwrap();
        pipeline
            .run(fgumi_pipeline_core::builder::PipelineConfig { threads, ..Default::default() })
            .unwrap();

        // Returned in ARRIVAL order, deliberately unsorted: the step declares
        // `BranchOrdering::ByItemOrdinal`, so the sink must already see blocks in
        // serial order. Sorting here would normalize out-of-order delivery and let
        // an ordering regression pass.
        std::mem::take(&mut *received.lock())
    }

    #[rstest]
    #[case::single_block_batches(1)]
    #[case::default_batching(DEFAULT_BLOCKS_PER_BATCH)]
    #[case::larger_than_the_file(1024)]
    fn try_run_emits_every_block_with_dense_serials(#[case] blocks_per_batch: usize) {
        const BGZF_EOF_LEN: usize = 28;
        let (path, on_disk) = temp_bam(64);
        let blocks = drive(&path, blocks_per_batch, 1);

        assert!(!blocks.is_empty(), "a non-empty BAM must yield at least one block");

        // Serials are dense and start at zero.
        for (i, block) in blocks.iter().enumerate() {
            assert_eq!(block.batch_serial, i as u64, "serial {i} must be dense");
            // `bytes` is the raw *compressed* block (this step does not inflate),
            // so it is not the same length as `uncompressed_size`; only that the
            // declared inflated size is populated is checkable here.
            assert!(block.uncompressed_size > 0, "block {i} must declare an inflated size");
            assert!(!block.bytes.is_empty(), "block {i} must carry its compressed bytes");
        }

        // Concatenating the payloads reproduces the file minus its BGZF EOF block,
        // which is what `FindBamBoundaries` downstream expects to receive.
        let concatenated: Vec<u8> = blocks.iter().flat_map(|b| b.bytes.clone()).collect();
        assert_eq!(concatenated, on_disk[..on_disk.len() - BGZF_EOF_LEN]);
    }

    #[rstest]
    #[case::sequential(fgumi_bam_io::ReadStreams::Fixed(1))]
    #[case::four_streams(fgumi_bam_io::ReadStreams::Fixed(4))]
    #[case::auto(fgumi_bam_io::ReadStreams::Auto)]
    fn read_streams_emit_identical_blocks(#[case] read_streams: fgumi_bam_io::ReadStreams) {
        // Routing a seekable file through the scatter reader (Fixed(4)/Auto)
        // must yield byte-identical blocks to the plain sequential reader
        // (Fixed(1)). 2000 records keep the file comfortably larger than the
        // 64-record fixture while staying fast to build.
        let (path, _) = temp_bam(2000);
        let baseline = drive(&path, 4, 1);
        let opts = PipelineReaderOpts { read_streams, ..PipelineReaderOpts::default() };
        let actual = drive_with_opts(&path, 4, 1, opts);
        assert_eq!(actual.len(), baseline.len(), "block count must not depend on read-streams");
        for (a, b) in actual.iter().zip(baseline.iter()) {
            assert_eq!(a.batch_serial, b.batch_serial);
            assert_eq!(a.bytes, b.bytes, "block bytes must be identical across read-streams");
        }
    }

    #[test]
    fn try_run_emits_the_same_blocks_regardless_of_thread_count() {
        let (path, _) = temp_bam(64);
        let one = drive(&path, 4, 1);
        let many = drive(&path, 4, 4);
        assert_eq!(one.len(), many.len(), "block count must not depend on threads");
        for (a, b) in one.iter().zip(many.iter()) {
            assert_eq!(a.batch_serial, b.batch_serial);
            assert_eq!(a.bytes, b.bytes);
        }
    }

    #[test]
    fn try_run_on_a_header_only_bam_still_emits_the_header_block() {
        const BGZF_EOF_LEN: usize = 28;
        let (path, on_disk) = temp_bam(0);
        let blocks = drive(&path, DEFAULT_BLOCKS_PER_BATCH, 1);
        let concatenated: Vec<u8> = blocks.iter().flat_map(|b| b.bytes.clone()).collect();
        assert_eq!(concatenated, on_disk[..on_disk.len() - BGZF_EOF_LEN]);
    }

    // ---------------------------------------------------------------------
    // Constructor + dispatch
    // ---------------------------------------------------------------------

    #[rstest]
    #[case::zero_clamps_to_one(0, 1)]
    #[case::one_stays_one(1, 1)]
    #[case::larger_is_preserved(32, 32)]
    fn new_clamps_blocks_per_batch_to_at_least_one(
        #[case] requested: usize,
        #[case] expected: usize,
    ) {
        let reader: Box<dyn io::Read + Send> = Box::new(io::Cursor::new(Vec::new()));
        let step = ReadBgzfBlocks::new(reader, requested, 1024);
        assert_eq!(step.blocks_per_batch, expected);
    }

    #[test]
    fn read_bam_auto_routes_a_regular_path_to_the_file_reader() {
        let (path, _) = temp_bam(4);
        // `read_bam_auto` on a non-stdin path must behave exactly like `read_bam`:
        // same header, and a step that reads the same file.
        let (_, auto_hdr) =
            read_bam_auto(&path, PipelineReaderOpts::default(), 4, 1024 * 1024).unwrap();
        let (_, direct_hdr) =
            read_bam(&path, PipelineReaderOpts::default(), 4, 1024 * 1024).unwrap();
        assert_eq!(auto_hdr, direct_hdr);
        assert!(!fgumi_bam_io::is_stdin_path(AsRef::<Path>::as_ref(&path)));
    }

    #[rstest]
    #[case::dash("-")]
    #[case::dev_stdin("/dev/stdin")]
    fn stdin_sentinels_are_recognised(#[case] sentinel: &str) {
        // Guards the branch condition in `read_bam_auto` without consuming the
        // process's real stdin, which a test must not do.
        assert!(fgumi_bam_io::is_stdin_path(Path::new(sentinel)));
    }

    #[test]
    fn read_bam_errors_on_a_missing_file() {
        // `ReadBgzfBlocks` is not `Debug`, so inspect the variant directly rather
        // than via `expect_err`.
        let result = read_bam(
            Path::new("/nonexistent/definitely/not/here.bam"),
            PipelineReaderOpts::default(),
            4,
            1024,
        );
        match result {
            Ok(_) => panic!("a missing file must not open"),
            Err(e) => assert!(!e.to_string().is_empty()),
        }
    }
}
