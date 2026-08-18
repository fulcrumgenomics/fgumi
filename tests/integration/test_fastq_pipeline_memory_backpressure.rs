//! Queue-memory backpressure from the output side back to the FASTQ Read step.
//!
//! The threaded FASTQ pipeline (`fgumi extract`) advertises a queue memory
//! budget (`--max-memory`, logged as "Queue memory budget: N total"). Every
//! stage between Read and Write is bounded, but the bound on most of them is a
//! *slot count*, not a byte count — so when the output device stalls, each
//! stage saturates with however many bytes its slots happen to hold and the
//! budget is not what decides the ceiling.
//!
//! Read is the pipeline's only source of new bytes, and nothing downstream
//! waits on it to make progress, so gating Read on the accounted in-flight
//! bytes bounds the whole pipeline without risking deadlock. These tests pin
//! that behaviour: with a writer that never completes, the pipeline must stop
//! pulling input rather than draining both FASTQ files into its queues.
//!
//! These are the FASTQ siblings of `test_pipeline_memory_backpressure.rs`.
//!
//! See <https://github.com/fulcrumgenomics/fgumi/issues/766>.

#![allow(clippy::cast_possible_truncation)]

use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record_buf::{QualityScores, Sequence};
use noodles_bgzf::io::Writer as BgzfWriter;
use std::fs::File;
use std::io::{self, BufRead, Cursor, Read, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};
use tempfile::TempDir;

use fgumi_lib::fastq_parse::strip_read_suffix;
use fgumi_lib::grouper::FastqTemplate;
use fgumi_lib::unified_pipeline::{
    FastqPipelineConfig, run_fastq_pipeline, serialize_bam_record_into,
};

/// Queue memory budget used by the tests, in bytes.
const TEST_BUDGET_BYTES: u64 = 2 * 1024 * 1024;

/// A budget generous enough to swallow the whole test input.
const GENEROUS_BUDGET_BYTES: u64 = 64 * 1024 * 1024;

/// Templates in the test input, and bases per read.
///
/// Sized so the input is many times what the pipeline can hold under
/// `TEST_BUDGET_BYTES`; a smaller input would be swallowed whole and the
/// assertions would prove nothing.
const TEMPLATES: usize = 30_000;
const READ_LENGTH: usize = 150;

/// Worker threads for the pipeline under test.
const NUM_THREADS: usize = 4;

/// Bytes the stalled writer accepts before it starts blocking.
///
/// `run_fastq_pipeline` writes the BAM header on the calling thread before it
/// spawns a single worker, so a writer that blocked from the very first byte
/// would stall the header write and the pipeline would never start. This
/// allowance lets the header (a few hundred bytes for a minimal one) through
/// and is negligible beside `TEST_BUDGET_BYTES`.
const WRITER_PRELUDE_BYTES: usize = 64 * 1024;

/// A reader that records how many bytes the pipeline has pulled from a stream.
///
/// `read_n_fastq_records` consumes exclusively through `BufRead::read_until`,
/// so `consume` sees every byte the Read step takes. `Read::read` is counted
/// too for completeness; it cannot double-count because it delegates to the
/// inner cursor rather than to `Self::consume`.
struct CountingReader {
    inner: Cursor<Vec<u8>>,
    consumed: Arc<AtomicU64>,
}

impl Read for CountingReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.consumed.fetch_add(n as u64, Ordering::Relaxed);
        Ok(n)
    }
}

impl BufRead for CountingReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    fn consume(&mut self, amt: usize) {
        self.consumed.fetch_add(amt as u64, Ordering::Relaxed);
        self.inner.consume(amt);
    }
}

/// A writer that blocks every write until it is released, then keeps the bytes.
///
/// Stands in for a saturated output device: the pipeline stays healthy (no
/// error, no panic) but makes no write progress at all. Retaining what it
/// eventually accepts lets the recovery test check record identity rather than
/// just a record count.
struct StalledWriter {
    released: Arc<AtomicBool>,
    sink: Arc<Mutex<Vec<u8>>>,
}

impl Write for StalledWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        // Generous enough never to fire in a passing run; bounded so a failed
        // assertion that unwinds before `released.store(true)` does not leave
        // this thread spinning for the rest of the test binary under
        // `cargo test` (where all tests share one process).
        let deadline = Instant::now() + Duration::from_secs(300);
        loop {
            {
                let sink = self.sink.lock().expect("sink mutex should not be poisoned");
                if self.released.load(Ordering::Acquire) || sink.len() < WRITER_PRELUDE_BYTES {
                    break;
                }
            }
            if Instant::now() >= deadline {
                return Err(io::Error::new(
                    io::ErrorKind::TimedOut,
                    "StalledWriter was never released",
                ));
            }
            std::thread::sleep(Duration::from_millis(5));
        }
        self.sink.lock().expect("sink mutex should not be poisoned").extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Build one FASTQ stream's bytes, and the read names it contains, in order.
///
/// Sequences vary per record so the queues have to hold real bytes rather than
/// a handful of repeated ones.
fn build_fastq_stream(num_records: usize, read_number: u8) -> (Vec<u8>, Vec<String>) {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut data = Vec::with_capacity(num_records * (READ_LENGTH * 2 + 40));
    let mut names = Vec::with_capacity(num_records);
    for i in 0..num_records {
        let name = format!("read{i}");
        let mut lcg = (i as u64).wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        let sequence: Vec<u8> = (0..READ_LENGTH)
            .map(|_| {
                lcg = lcg
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                bases[((lcg >> 33) & 3) as usize]
            })
            .collect();
        data.extend_from_slice(format!("@{name}/{read_number}\n").as_bytes());
        data.extend_from_slice(&sequence);
        data.push(b'\n');
        data.extend_from_slice(b"+\n");
        data.extend(std::iter::repeat_n(b'I', READ_LENGTH));
        data.push(b'\n');
        names.push(name);
    }
    (data, names)
}

/// Read the record names out of a BAM byte stream, in order.
fn record_names(bam_bytes: &[u8]) -> Vec<String> {
    let mut reader = bam::io::Reader::new(Cursor::new(bam_bytes));
    let header = reader.read_header().expect("read header");
    reader
        .record_bufs(&header)
        .map(|r| {
            let record = r.expect("read record");
            String::from_utf8_lossy(record.name().expect("record should have a name").as_ref())
                .into_owned()
        })
        .collect()
}

/// A running pipeline whose writer is stalled, with the handles a test needs to
/// observe it, release it, and inspect what it eventually wrote.
struct StalledRun {
    /// Input bytes the Read step has pulled so far, summed over both streams.
    consumed: Arc<AtomicU64>,
    /// Set to release the writer.
    released: Arc<AtomicBool>,
    /// Output bytes accepted, as a complete BAM stream.
    sink: Arc<Mutex<Vec<u8>>>,
    /// The pipeline itself, yielding the count of records written.
    handle: std::thread::JoinHandle<io::Result<u64>>,
}

/// Spawn a paired-end pass-through pipeline whose writer never completes.
///
/// The streams are handed in already decompressed, which is the plain/gzip
/// input path: Read → Decompress (pass-through) → Pair → Parse → Process →
/// Serialize → Compress → Write.
fn spawn_stalled_pipeline(
    header: &Header,
    streams: Vec<Vec<u8>>,
    queue_memory_limit: u64,
) -> StalledRun {
    let consumed = Arc::new(AtomicU64::new(0));
    let released = Arc::new(AtomicBool::new(false));
    let sink = Arc::new(Mutex::new(Vec::new()));

    let readers: Vec<Box<dyn BufRead + Send>> = streams
        .into_iter()
        .map(|bytes| {
            let reader =
                CountingReader { inner: Cursor::new(bytes), consumed: Arc::clone(&consumed) };
            Box::new(reader) as Box<dyn BufRead + Send>
        })
        .collect();
    let writer = StalledWriter { released: Arc::clone(&released), sink: Arc::clone(&sink) };

    // Production auto-tuning on purpose: the slot counts must stay large enough
    // to swallow the whole test input, so that any bound the test observes comes
    // from the byte budget and not from a queue running out of slots.
    let mut config = FastqPipelineConfig::new(NUM_THREADS, false, 1);
    config.queue_memory_limit = queue_memory_limit;
    // A permanently stalled writer is exactly what the deadlock detector is
    // meant to flag, and flagging it would tear the pipeline down before the
    // assertion. Disable it so the test observes backpressure, not recovery.
    config.deadlock_timeout_secs = 0;

    let owned_header = header.clone();
    let handle = std::thread::spawn(move || {
        run_fastq_pipeline(
            config,
            &[],
            Some(readers),
            &owned_header,
            Box::new(writer),
            |template: FastqTemplate| -> io::Result<RecordBuf> {
                let record = template
                    .records
                    .first()
                    .ok_or_else(|| io::Error::other("template with no records"))?;
                Ok(RecordBuf::builder()
                    .set_name(strip_read_suffix(record.name()).to_vec())
                    .set_sequence(Sequence::from(record.sequence().to_vec()))
                    .set_quality_scores(QualityScores::from(
                        record.quality().iter().map(|q| q - 33).collect::<Vec<_>>(),
                    ))
                    .build())
            },
            |record: RecordBuf, header: &Header, buf: &mut Vec<u8>| -> io::Result<u64> {
                serialize_bam_record_into(&record, header, buf)
            },
        )
    });

    StalledRun { consumed, released, sink, handle }
}

/// Outcome of waiting for the reader to stop pulling input.
struct Consumption {
    /// Input bytes consumed when the wait ended.
    bytes: u64,
    /// Whether consumption actually settled (or reached EOF), as opposed to the
    /// wait timing out with the counter still moving. A timeout means the
    /// caller's premise — that the reader had stopped — never held, so tests
    /// assert on this to fail loudly rather than silently comparing a
    /// still-moving figure.
    settled: bool,
}

/// Wait until input consumption stops growing, or `timeout` elapses.
///
/// Reports the number of input bytes consumed once the counter has held still
/// across `STABLE_SAMPLES` consecutive samples, or as soon as the whole input
/// has been read, with `settled: true`. Several samples rather than one guards
/// against a loaded CI machine descheduling the reader briefly and looking like
/// backpressure. If the deadline is reached with the counter still moving, it
/// reports the last figure with `settled: false`.
fn wait_for_consumption_to_settle(
    consumed: &AtomicU64,
    input_len: u64,
    timeout: Duration,
) -> Consumption {
    /// Consecutive unchanged samples required before calling it settled.
    const STABLE_SAMPLES: u32 = 4;

    let deadline = Instant::now() + timeout;
    let mut previous = u64::MAX;
    let mut stable = 0u32;
    while Instant::now() < deadline {
        std::thread::sleep(Duration::from_millis(250));
        let current = consumed.load(Ordering::Relaxed);
        if current >= input_len {
            return Consumption { bytes: current, settled: true };
        }
        if current == previous {
            stable += 1;
            if stable >= STABLE_SAMPLES {
                return Consumption { bytes: current, settled: true };
            }
        } else {
            stable = 0;
        }
        previous = current;
    }
    Consumption { bytes: consumed.load(Ordering::Relaxed), settled: false }
}

/// Build the paired-end test input: the two stream buffers, the expected
/// template names, and the total input size in bytes.
fn build_paired_input() -> (Vec<Vec<u8>>, Vec<String>, u64) {
    let (r1, names) = build_fastq_stream(TEMPLATES, 1);
    let (r2, _) = build_fastq_stream(TEMPLATES, 2);
    let input_len = (r1.len() + r2.len()) as u64;
    (vec![r1, r2], names, input_len)
}

/// A stalled writer must stop the Read step, not let it drain both FASTQ files
/// into the pipeline's queues.
///
/// Without a byte-aware admission gate on Read, every downstream queue
/// saturates at its slot count and the reader runs to EOF regardless of
/// `--max-memory` — which is how a slow output device turns into an OOM kill on
/// an input far larger than the host's memory (issue #766).
#[test]
fn stalled_writer_stops_the_fastq_reader_within_the_queue_memory_budget() {
    let header = Header::default();
    let (streams, expected_names, input_len) = build_paired_input();
    assert!(
        input_len > TEST_BUDGET_BYTES,
        "test input ({input_len} bytes) must exceed the budget ({TEST_BUDGET_BYTES} bytes) \
         for the assertion to mean anything"
    );

    let StalledRun { consumed, released, sink, handle } =
        spawn_stalled_pipeline(&header, streams, TEST_BUDGET_BYTES);

    let settled = wait_for_consumption_to_settle(&consumed, input_len, Duration::from_secs(30));
    assert!(settled.settled, "reader never stopped under a stalled writer within the timeout");
    assert!(
        settled.bytes < input_len / 2,
        "reader consumed {} of {input_len} input bytes while the writer was stalled; \
         the queue memory budget bounded nothing",
        settled.bytes
    );

    // Throttling the reader is only worth anything if the run still emits every
    // record afterwards, so release and check identity rather than exit status.
    released.store(true, Ordering::Release);
    let written =
        handle.join().expect("pipeline thread should not panic").expect("pipeline should succeed");
    assert_eq!(written, expected_names.len() as u64);
    let output = sink.lock().expect("sink mutex should not be poisoned").clone();
    assert_eq!(record_names(&output), expected_names, "throttling must not change the output");
}

/// How far the reader runs ahead of a stalled writer must depend on the queue
/// memory budget.
///
/// This is the discriminating check. The previous test could in principle be
/// satisfied by the queues simply running out of *slots*; this one holds the
/// input and the slot counts fixed and varies only `--max-memory`, so a budget
/// that bounds nothing shows up as two identical read-ahead figures.
#[test]
fn fastq_read_ahead_under_a_stalled_writer_shrinks_with_the_queue_memory_budget() {
    let header = Header::default();
    let (streams, expected_names, input_len) = build_paired_input();

    let mut consumption = Vec::new();
    for budget in [TEST_BUDGET_BYTES, GENEROUS_BUDGET_BYTES] {
        let StalledRun { consumed, released, sink, handle } =
            spawn_stalled_pipeline(&header, streams.clone(), budget);
        let settled = wait_for_consumption_to_settle(&consumed, input_len, Duration::from_secs(30));
        assert!(
            settled.settled,
            "reader never stopped under a {budget}-byte budget within the timeout"
        );
        released.store(true, Ordering::Release);
        handle.join().expect("pipeline thread should not panic").expect("pipeline should succeed");
        let output = sink.lock().expect("sink mutex should not be poisoned").clone();
        assert_eq!(
            record_names(&output),
            expected_names,
            "budget {budget} changed the output, not just the read-ahead"
        );
        consumption.push(settled.bytes);
    }

    // A margin tied to the budget difference, not bare ordering: a scheduler
    // hiccup that made the tight run read one extra block would still pass
    // `tight < generous`, and the generous budget exceeds the whole input, so
    // that run reads to EOF and a bounded gate would be indistinguishable from
    // none. Requiring the gap to span at least one tight budget pins that the
    // budget — not slot counts or timing — is the bound.
    let (tight, generous) = (consumption[0], consumption[1]);
    assert!(
        generous.saturating_sub(tight) >= TEST_BUDGET_BYTES,
        "read-ahead was {tight} bytes under a {TEST_BUDGET_BYTES}-byte budget and {generous} \
         bytes under a {GENEROUS_BUDGET_BYTES}-byte one; the gap is smaller than one tight \
         budget, so the budget is not bounding the queues"
    );
}

/// The pipeline must still emit every template, in order, once the writer
/// recovers.
///
/// Backpressure that stalls the reader is only safe if releasing the writer
/// lets every record through. The oracle is the input's own read names read
/// back off the output, not a record count: a gate that dropped, duplicated or
/// reordered records while the pipeline was throttled would keep the count
/// intact and change the sequence.
#[test]
fn fastq_reader_resumes_and_completes_after_the_writer_recovers() {
    let header = Header::default();
    let (streams, expected_names, input_len) = build_paired_input();

    let StalledRun { consumed, released, sink, handle } =
        spawn_stalled_pipeline(&header, streams, TEST_BUDGET_BYTES);

    let settled = wait_for_consumption_to_settle(&consumed, input_len, Duration::from_secs(30));
    assert!(settled.settled, "reader never stopped under a stalled writer within the timeout");
    released.store(true, Ordering::Release);

    let written = handle.join().expect("pipeline thread should not panic").expect("pipeline");
    assert_eq!(written, expected_names.len() as u64, "every input template must reach the writer");
    assert_eq!(
        consumed.load(Ordering::Relaxed),
        input_len,
        "the whole input must be consumed once the writer recovers"
    );

    let output = sink.lock().expect("sink mutex should not be poisoned").clone();
    assert_eq!(
        record_names(&output),
        expected_names,
        "output records must match the input exactly, in order"
    );
}

/// A record-count-mismatched FASTQ pair must be rejected under a tight budget
/// too, and rejected the same way it is with no budget at all.
///
/// This is the liveness hazard the Read gate has to design around. Neither the
/// Pair step nor `BlockMerge` can release a batch that is missing a stream
/// until *every* stream has reached EOF, so once the shorter file is exhausted
/// the longer file's surplus piles up unreleasable — and it counts against the
/// budget. A gate that simply stopped reading at the budget would wedge there,
/// because only Read can bring the surviving stream to the EOF that reveals the
/// mismatch. `read_is_required_for_liveness` is what keeps that from happening:
/// Read stays admitted so the pipeline reaches EOF, detects the unequal record
/// counts, and fails with the out-of-sync error (issue `#773`) rather than
/// hanging on the surplus tail.
///
/// What this test pins is that the budget changes neither the outcome (the
/// out-of-sync rejection) nor the ability to reach it within the timeout.
///
/// The writer is deliberately fast here: the stall is not what closes the gate,
/// the unreleasable surplus is.
#[test]
fn mismatched_stream_lengths_are_rejected_under_a_tight_budget() {
    let header = Header::default();
    let short = TEMPLATES / 10;
    let (r1, _) = build_fastq_stream(short, 1);
    let (r2, _) = build_fastq_stream(TEMPLATES, 2);
    let streams = vec![r1, r2];

    let ungated = run_expecting_out_of_sync(&header, streams.clone(), 0);
    let gated = run_expecting_out_of_sync(&header, streams, TEST_BUDGET_BYTES);
    // R1 is the shorter stream, so it must be named as the one that ended first
    // regardless of the budget — the surplus count itself can legitimately vary
    // with which stage detects the mismatch, so it is not part of the oracle.
    for (label, message) in [("ungated", &ungated), ("gated", &gated)] {
        assert!(
            message.contains("R1 ended before R2"),
            "the {label} run must identify R1 as the short stream; got: {message}"
        );
    }
}

/// Run a pipeline with an immediately-released writer and return the terminal
/// error message it produced, failing the test rather than hanging if the Read
/// gate wedges on the surplus tail.
///
/// The message rather than a bare `is_err`: what this pins is that the pipeline
/// reached EOF on the surviving stream and rejected the *mismatch*, not that it
/// stopped for some unrelated reason.
fn run_expecting_out_of_sync(
    header: &Header,
    streams: Vec<Vec<u8>>,
    queue_memory_limit: u64,
) -> String {
    let StalledRun { released, handle, .. } =
        spawn_stalled_pipeline(header, streams, queue_memory_limit);
    // Never stall: this test is about the gate, not about writer backpressure.
    released.store(true, Ordering::Release);

    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(handle.join());
    });
    let joined = rx.recv_timeout(Duration::from_secs(45)).unwrap_or_else(|_| {
        panic!(
            "pipeline did not finish within 45s under a {queue_memory_limit}-byte budget; \
             the Read gate wedged on the surplus tail of the longer stream"
        )
    });
    let error = joined
        .expect("pipeline thread should not panic")
        .expect_err("a record-count-mismatched pair must be rejected");
    assert_eq!(
        error.kind(),
        io::ErrorKind::InvalidData,
        "a mismatched pair must fail with InvalidData, got {error:?}"
    );
    let message = error.to_string();
    assert!(
        message.contains("FASTQ sources out of sync"),
        "expected the out-of-sync rejection under a {queue_memory_limit}-byte budget, got: {message}"
    );
    message
}

// ============================================================================
// BGZF input path
// ============================================================================
//
// The tests above hand the pipeline already-decompressed readers, which is the
// gzip/plain path: Read -> Decompress -> Pair -> Parse. BGZF inputs take a
// different route -- Read -> Decompress -> BlockParseFast -> BlockMerge -- with
// its own byte accounting that nothing above exercises:
// `block_merge_pending_heap_bytes` (the lock-free mirror of
// `BlockMergeState::pending_heap_bytes`), the `q2_block_parsed` charge taken at
// creation rather than at push, and `drain_exhausted_stream`'s refund.
//
// A leak on any of those inflates a counter the Read gate sums, so under a
// tight budget the gate closes on phantom bytes and never reopens -- the run
// wedges instead of finishing. That is what these tests detect, which is why
// they assert completion and record identity rather than a read-ahead figure:
// `run_fastq_pipeline` opens BGZF inputs from paths itself, so there is no
// `CountingReader` to interpose on this path.

/// Write one stream's bytes to a BGZF-compressed FASTQ file.
///
/// Flushes every `RECORDS_PER_BLOCK` records so the input spans many BGZF
/// blocks: `BlockMerge` pairs blocks by index across the two streams, so a
/// single-block input would never exercise its pending maps at all.
fn write_bgzf_fastq(dir: &TempDir, name: &str, num_records: usize, read_number: u8) -> PathBuf {
    /// Records per BGZF block, chosen so the fixture spans hundreds of blocks.
    const RECORDS_PER_BLOCK: usize = 64;

    let path = dir.path().join(name);
    let mut writer = BgzfWriter::new(File::create(&path).expect("create BGZF FASTQ"));
    let (data, _) = build_fastq_stream(num_records, read_number);
    // Read names vary in width, so this is the mean record size rather than an
    // exact one — which is if anything better here, since the resulting block
    // boundaries fall mid-record and exercise the cross-block stitching too.
    let bytes_per_block = (data.len() / num_records).max(1) * RECORDS_PER_BLOCK;
    for chunk in data.chunks(bytes_per_block) {
        writer.write_all(chunk).expect("write BGZF FASTQ");
        writer.flush().expect("flush BGZF block");
    }
    writer.finish().expect("finish BGZF FASTQ");
    path
}

/// Run the BGZF input path to completion under `queue_memory_limit` and return
/// the read names emitted, failing rather than hanging if the Read gate wedges.
fn run_bgzf_to_completion(
    header: &Header,
    paths: &[PathBuf],
    queue_memory_limit: u64,
) -> Vec<String> {
    let sink = Arc::new(Mutex::new(Vec::new()));
    let writer =
        StalledWriter { released: Arc::new(AtomicBool::new(true)), sink: Arc::clone(&sink) };

    let mut config = FastqPipelineConfig::new(NUM_THREADS, true, 1);
    config.queue_memory_limit = queue_memory_limit;
    config.deadlock_timeout_secs = 0;

    let owned_header = header.clone();
    let owned_paths = paths.to_vec();
    let handle = std::thread::spawn(move || {
        run_fastq_pipeline(
            config,
            &owned_paths,
            None,
            &owned_header,
            Box::new(writer),
            |template: FastqTemplate| -> io::Result<RecordBuf> {
                let record = template
                    .records
                    .first()
                    .ok_or_else(|| io::Error::other("template with no records"))?;
                Ok(RecordBuf::builder()
                    .set_name(strip_read_suffix(record.name()).to_vec())
                    .set_sequence(Sequence::from(record.sequence().to_vec()))
                    .set_quality_scores(QualityScores::from(
                        record.quality().iter().map(|q| q - 33).collect::<Vec<_>>(),
                    ))
                    .build())
            },
            |record: RecordBuf, header: &Header, buf: &mut Vec<u8>| -> io::Result<u64> {
                serialize_bam_record_into(&record, header, buf)
            },
        )
    });

    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(handle.join());
    });
    let joined = rx.recv_timeout(Duration::from_secs(60)).unwrap_or_else(|_| {
        panic!(
            "BGZF pipeline did not finish within 60s under a {queue_memory_limit}-byte budget; \
             the Read gate is holding on bytes the BGZF accounting never refunded"
        )
    });
    let written =
        joined.expect("pipeline thread should not panic").expect("pipeline should succeed");
    let output = sink.lock().expect("sink mutex should not be poisoned").clone();
    let names = record_names(&output);
    assert_eq!(names.len() as u64, written, "every reported write must be in the output");
    names
}

/// The BGZF input path must finish under a budget far smaller than its input,
/// and emit exactly what an unbudgeted run of the same input emits.
///
/// This is the BGZF half of the accounting the Read gate depends on. Its
/// charges and refunds are paired across three sites the gzip path never
/// touches — the `q2_block_parsed` charge taken at creation, the
/// `block_merge_pending_heap_bytes` mirror, and `drain_exhausted_stream`'s
/// refund — and an unpaired one there does not merely skew a statistic: the
/// phantom bytes accumulate in a counter `queue_bytes_in_flight` sums, so the
/// gate eventually closes for good and the run never terminates. A tight budget
/// is what makes that reachable inside a test's lifetime; the timeout in
/// `run_bgzf_to_completion` is what reports it.
#[test]
fn bgzf_inputs_complete_under_a_tight_budget_and_match_an_unbudgeted_run() {
    let header = Header::default();
    let dir = TempDir::new().expect("create tempdir");
    let r1 = write_bgzf_fastq(&dir, "r1.fastq.gz", TEMPLATES, 1);
    let r2 = write_bgzf_fastq(&dir, "r2.fastq.gz", TEMPLATES, 2);
    let paths = vec![r1, r2];

    let input_len: u64 =
        paths.iter().map(|p| std::fs::metadata(p).expect("stat BGZF input").len()).sum();
    assert!(
        input_len > TEST_BUDGET_BYTES,
        "compressed input ({input_len} bytes) must exceed the budget ({TEST_BUDGET_BYTES} bytes) \
         for the gate to engage at all"
    );

    let (_, expected_names) = build_fastq_stream(TEMPLATES, 1);
    let ungated = run_bgzf_to_completion(&header, &paths, 0);
    assert_eq!(
        ungated, expected_names,
        "the unbudgeted BGZF run must emit every template in order"
    );

    let gated = run_bgzf_to_completion(&header, &paths, TEST_BUDGET_BYTES);
    assert_eq!(
        gated, ungated,
        "the queue memory budget must throttle the BGZF path without changing its output"
    );
}
