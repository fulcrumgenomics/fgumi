//! The chain-native input trio `PlanInputReads → PreadInputSlices →
//! FrameBgzfBlocks`, driven through real pipelines.

use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Duration;

use fgumi_bam_io::PipelineReaderOpts;
use fgumi_bam_io::pread::test_sources::TimedSource;
use fgumi_bam_io::pread::{FILL_BYTES, PositionalSource, ReadStreamsPolicy, SliceBufferPool};
use fgumi_pipeline_core::builder::{Pipeline, PipelineConfig};
use fgumi_pipeline_core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use fgumi_pipeline_core::{PipelineError, PoolPlacement};
use rstest::rstest;

use super::frame_bgzf_blocks::{AFTER_FRAME as frame_hook, FrameBgzfBlocks};
use super::plan_input_reads::PlanInputReads;
use crate::pread::{AFTER_READ as pread_hook, InputLedger, PreadSlices};
use crate::types::BgzfBlock;

/// Sink that records every block in arrival order.
struct BlockSink {
    got: Arc<parking_lot::Mutex<Vec<BgzfBlock>>>,
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
            Some(b) => {
                self.got.lock().push(b);
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }
}

/// A BAM of `records` stored-block (incompressible-size) records, so the file
/// spans several 4 MiB fills: ~270 bytes per record.
pub(crate) fn bam_fixture(records: usize) -> (tempfile::TempDir, PathBuf) {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("in.bam");
    let header = noodles::sam::Header::default();
    let mut writer = fgumi_bam_io::create_raw_bam_writer(&path, &header, 1, 0).unwrap();
    for i in 0..records {
        let name = format!("q{:08x}", i.wrapping_mul(2_654_435_761));
        let bytes = fgumi_raw_bam::testutil::make_bam_bytes(
            0,
            i32::try_from(i % 1_000_000).unwrap(),
            0,
            name.as_bytes(),
            &[],
            150,
            -1,
            -1,
            &[],
        );
        writer.write_raw_record(&bytes).unwrap();
    }
    writer.finish().unwrap();
    (dir, path)
}

/// The sequential oracle: `ReadBgzfBlocks` over a plain reader.
pub(crate) fn read_with_read_bgzf_blocks(path: &Path) -> Vec<BgzfBlock> {
    let got = Arc::new(parking_lot::Mutex::new(Vec::new()));
    let (source, _hdr) =
        super::read_bam::read_bam(path, PipelineReaderOpts::default(), 16, 1 << 20).expect("open");
    let builder = Pipeline::builder();
    builder.chain(source).chain(BlockSink { got: Arc::clone(&got) }).into_sink_marker();
    builder.build().unwrap().run(PipelineConfig { threads: 2, ..Default::default() }).unwrap();
    std::mem::take(&mut *got.lock())
}

/// Run the native trio over `source` and return the blocks in arrival order.
pub(crate) fn read_with_native_trio_from(
    source: Arc<dyn PositionalSource>,
    policy: Arc<ReadStreamsPolicy>,
    threads: usize,
) -> Result<Vec<BgzfBlock>, PipelineError> {
    let len = source.byte_len().unwrap();
    let got = Arc::new(parking_lot::Mutex::new(Vec::new()));
    let ledger = Arc::new(InputLedger::default());
    let eligible = fgumi_pipeline_core::runtime::parallel_hosts(
        PoolPlacement::ExcludeReader,
        Some(0),
        threads,
        None,
    )
    .clone_count();
    let plan =
        PlanInputReads::from_source(source, len, policy, Arc::clone(&ledger), eligible, 1 << 20);
    let pread =
        PreadSlices::input(SliceBufferPool::new(2 * threads + 8), Arc::clone(&ledger), 8 << 20);
    let frame = FrameBgzfBlocks::new(ledger, 1 << 20);
    let builder = Pipeline::builder();
    builder
        .chain(plan)
        .chain(pread)
        .chain(frame)
        .chain(BlockSink { got: Arc::clone(&got) })
        .into_sink_marker();
    builder.build().unwrap().run(PipelineConfig { threads, ..Default::default() })?;
    Ok(std::mem::take(&mut *got.lock()))
}

fn open_source(path: &Path) -> Arc<dyn PositionalSource> {
    Arc::new(std::fs::File::open(path).unwrap())
}

/// The native trio over the file at `path`.
pub(crate) fn read_with_native_trio(
    path: &Path,
    policy: Arc<ReadStreamsPolicy>,
    threads: usize,
) -> Vec<BgzfBlock> {
    read_with_native_trio_from(open_source(path), policy, threads).expect("native read")
}

fn assert_same_blocks(got: &[BgzfBlock], oracle: &[BgzfBlock]) {
    assert_eq!(got.len(), oracle.len(), "block count");
    for (i, (a, b)) in got.iter().zip(oracle).enumerate() {
        assert_eq!(a.batch_serial, i as u64, "dense serials in arrival order");
        assert_eq!(
            (a.batch_serial, a.uncompressed_size, &a.bytes[..]),
            (b.batch_serial, b.uncompressed_size, &b.bytes[..]),
            "block {i}"
        );
    }
}

/// The native trio emits exactly the blocks `ReadBgzfBlocks` emits (same
/// bytes, same `uncompressed_size`, `batch_serial` dense from 0, in order),
/// for every stream count and thread count.
#[rstest]
fn native_reads_emit_identical_blocks(
    #[values(1usize, 2, 4)] threads: usize,
    #[values(1usize, 4, 8)] streams: usize,
) {
    let (_dir, path) = bam_fixture(40_000);
    let oracle = read_with_read_bgzf_blocks(&path);
    assert!(std::fs::metadata(&path).unwrap().len() > 2 * FILL_BYTES as u64, "several fills");
    let got = read_with_native_trio(&path, ReadStreamsPolicy::fixed(streams), threads);
    assert_same_blocks(&got, &oracle);
}

/// A BAM truncated inside its last block fails the run.
#[test]
fn truncated_bam_fails_the_native_path() {
    let (_dir, path) = bam_fixture(2_000);
    let len = std::fs::metadata(&path).unwrap().len();
    std::fs::OpenOptions::new().write(true).open(&path).unwrap().set_len(len - 37).unwrap();
    let err = read_with_native_trio_from(open_source(&path), ReadStreamsPolicy::fixed(4), 4)
        .expect_err("a truncated input must fail");
    assert!(err.to_string().contains("truncated BGZF block"), "{err}");
}

/// A `PositionalSource` over a file that records the highest byte requested.
struct RecordingSource {
    file: std::fs::File,
    max_end: Arc<AtomicU64>,
}

impl PositionalSource for RecordingSource {
    fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
        self.max_end.fetch_max(offset + buf.len() as u64, Ordering::SeqCst);
        PositionalSource::read_at(&self.file, buf, offset)
    }
    fn byte_len(&self) -> io::Result<u64> {
        PositionalSource::byte_len(&self.file)
    }
}

/// The lookahead bound, observed outside the planner's ledger: the source
/// records the highest byte requested, the framer's hook the bytes framed;
/// after every framed slice the difference stays within 8 MiB + one fill.
#[test]
fn lookahead_is_bounded() {
    let (_dir, path) = bam_fixture(200_000);
    assert!(std::fs::metadata(&path).unwrap().len() > 32 << 20, "the bound must bind");
    let max_end = Arc::new(AtomicU64::new(0));
    let worst = Arc::new(AtomicU64::new(0));
    let framed = Arc::new(AtomicU64::new(0));
    let src = Arc::new(RecordingSource {
        file: std::fs::File::open(&path).unwrap(),
        max_end: Arc::clone(&max_end),
    });
    let (m2, w2, f2) = (Arc::clone(&max_end), Arc::clone(&worst), Arc::clone(&framed));
    frame_hook.set(move |slice_len| {
        let f = f2.fetch_add(slice_len as u64, Ordering::SeqCst) + slice_len as u64;
        w2.fetch_max(m2.load(Ordering::SeqCst).saturating_sub(f), Ordering::SeqCst);
        // A slow framer lets the planner run ahead if it is not bounded.
        std::thread::sleep(Duration::from_micros(200));
    });
    let result = read_with_native_trio_from(src, ReadStreamsPolicy::fixed(4), 4);
    frame_hook.clear();
    result.expect("native read");
    let worst = worst.load(Ordering::SeqCst);
    assert!(worst > 0, "the hook observed the run");
    assert!(worst <= (8 << 20) + FILL_BYTES as u64, "lookahead exceeded: {worst}");
}

/// A BGZF stream of `bytes` stored blocks (64 KiB payloads) plus an EOF
/// marker.
fn bgzf_stream(bytes: usize) -> Vec<u8> {
    let mut compressor = fgumi_bgzf::InlineBgzfCompressor::new(0);
    let payload: Vec<u8> = (0..bytes).map(|i| u8::try_from(i % 251).unwrap()).collect();
    compressor.write_all(&payload).unwrap();
    compressor.flush().unwrap();
    let mut out = Vec::new();
    compressor.write_blocks_to(&mut out).unwrap();
    out.extend_from_slice(&fgumi_bgzf::BGZF_EOF);
    out
}

/// Ratchet signal wiring smoke test (the deterministic gate is the predicate
/// table in `fgumi_bam_io::pread`): a slow source with an unbounded consumer
/// ratchets; a fast source with a slow framer does not. The delays leave a
/// ≥ 10× margin either way, also on a loaded host: a 40 ms read of a 4 MiB
/// fill (100 MB/s) against framing and inflating 4 MiB of stored blocks (a few
/// ms even when the CPU is oversubscribed — at 4 ms the full suite's load made
/// the consumer, not the device, the bottleneck); a 4 ms frame against µs
/// reads.
#[rstest]
#[case::slow_device(Duration::from_millis(40), Duration::ZERO, true)]
#[case::slow_framer(Duration::ZERO, Duration::from_millis(4), false)]
fn ratchet_is_driven_by_device_starvation_only(
    #[case] read_delay: Duration,
    #[case] frame_delay: Duration,
    #[case] expect_raise: bool,
) {
    let policy = ReadStreamsPolicy::auto();
    let src: Arc<dyn PositionalSource> =
        Arc::new(TimedSource::new(bgzf_stream(64 << 20), read_delay));
    if !frame_delay.is_zero() {
        frame_hook.set(move |_| std::thread::sleep(frame_delay));
    }
    let result = read_with_native_trio_from(src, Arc::clone(&policy), 4);
    frame_hook.clear();
    result.expect("native read");
    assert_eq!(policy.streams() > 1, expect_raise, "history {:?}", policy.history());
}

/// `PreadInputSlices` never runs on the reader worker. The deterministic half
/// is the placement plan; the run checks the thread names a `PreadSlices` test
/// hook recorded.
#[test]
fn input_slices_never_run_on_the_reader_worker() {
    assert_eq!(
        fgumi_pipeline_core::runtime::parallel_hosts(
            PoolPlacement::ExcludeReader,
            Some(0),
            4,
            None
        )
        .workers,
        vec![1, 2, 3]
    );
    let names = Arc::new(parking_lot::Mutex::new(std::collections::BTreeSet::new()));
    let n2 = Arc::clone(&names);
    pread_hook.set(move |()| {
        n2.lock().insert(std::thread::current().name().unwrap_or("<unnamed>").to_owned());
    });
    let (_dir, path) = bam_fixture(50_000);
    let result = read_with_native_trio_from(open_source(&path), ReadStreamsPolicy::fixed(4), 4);
    pread_hook.clear();
    result.expect("native read");
    let names = names.lock();
    assert!(!names.is_empty(), "the hook observed the reads");
    assert!(!names.contains("fgumi-worker-0"), "{names:?}");
}

/// A `PositionalSource` over a file that records the name of every thread a
/// read runs on.
struct ThreadNameSource {
    file: std::fs::File,
    names: Arc<parking_lot::Mutex<std::collections::BTreeSet<String>>>,
}

impl PositionalSource for ThreadNameSource {
    fn read_at(&self, buf: &mut [u8], offset: u64) -> io::Result<usize> {
        let name = std::thread::current().name().unwrap_or("<unnamed>").to_owned();
        self.names.lock().insert(name);
        PositionalSource::read_at(&self.file, buf, offset)
    }
    fn byte_len(&self) -> io::Result<u64> {
        PositionalSource::byte_len(&self.file)
    }
}

/// Every positional read runs on a pipeline thread — never on a thread
/// spawned for the read (the oracle is the thread the source sees, so a
/// scoped per-read thread is caught even though it exits before any sample of
/// the process's thread count).
#[rstest]
fn every_read_runs_on_a_pool_worker(#[values(1usize, 4)] threads: usize) {
    let (_dir, path) = bam_fixture(40_000);
    let names = Arc::new(parking_lot::Mutex::new(std::collections::BTreeSet::new()));
    let src = Arc::new(ThreadNameSource {
        file: std::fs::File::open(&path).unwrap(),
        names: Arc::clone(&names),
    });
    let _ = read_with_native_trio_from(src, ReadStreamsPolicy::fixed(4), threads).unwrap();
    // At one thread this test's three-step chain is fused onto the calling
    // thread; that is a pipeline thread too.
    let caller = std::thread::current().name().unwrap_or("<unnamed>").to_owned();
    let names = names.lock();
    assert!(!names.is_empty());
    assert!(names.iter().all(|n| n.starts_with("fgumi-worker-") || *n == caller), "{names:?}");
}
