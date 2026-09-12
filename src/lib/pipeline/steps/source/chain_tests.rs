//! End-to-end tests for the read-side source chains.
//!
//! The source steps are `try_run` state machines wrapped around real readers,
//! so their interesting behavior — round-robin stream rotation, chunk
//! alignment to whole 4-line records, per-stream ordinal assignment, drain
//! reporting, and held-item retry under backpressure — only happens when the
//! framework drives them against an actual stream. These tests assemble the
//! the unified K-stream FASTQ topology and the SAM topology and run them
//! through `Pipeline::run`:
//!
//! ```text
//! ReadFastqInputs::new_single(i) ×K → ZipRawFastqK(K) → ParseAndZipFastqN → sink
//! ReadFastqInputs::new_single(0)     → WrapRawFastq1   → ParseAndZipFastqN → sink   (K==1)
//! ReadSamChunks → sink
//! ```
//!
//! One per-stream reader feeds each of the K producer edges of `ZipRawFastqK`
//! (the lockstep aligner), which mints a dense ordinal that the `Parallel`
//! `ParseAndZipFastqN` reorders by. For a lone stream `WrapRawFastq1` is the
//! single-input shim (`ZipRawFastqK` rejects `k == 1`). All K must recover the
//! same templates in the same order from the same input.
//!
//! `zip_raw_fastq_k` and `parse_zip_fastq` carry their own targeted step tests;
//! these cover the plain read-through path those tests stub out with synthetic
//! sources.

use std::io::{self, BufRead};
use std::sync::{Arc, Mutex};

use rstest::rstest;

use crate::pipeline::core::builder::{Pipeline, PipelineBuilder, PipelineConfig};
use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::pipeline::core::step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::core::topology::{BranchIdx, StepIdx};
use crate::pipeline::steps::source::fastq_zip::FastqTemplateBatch;
use crate::pipeline::steps::source::parse_zip_fastq::ParseAndZipFastqN;
use crate::pipeline::steps::source::read_fastq::{FastqOrdinalSequence, ReadFastqInputs};
use crate::pipeline::steps::source::read_sam_chunks::ReadSamChunks;
use crate::pipeline::steps::source::zip_raw_fastq_k::{WrapRawFastq1, ZipRawFastqK};
use crate::pipeline::steps::types::SamChunk;

/// Per-edge byte budget, deliberately small relative to the fixtures so the
/// byte-bounded edges reject pushes mid-run and the sources' held-item retry
/// paths run. `the_edge_budget_binds_on_the_fixture` pins that it has teeth.
const EDGE_LIMIT_BYTES: u64 = 512;

/// FASTQ records per emitted chunk. Small enough that the fixtures below span
/// many chunks, so stream rotation and cross-chunk ordinal assignment are
/// actually exercised rather than fitting in a single chunk.
const BATCH_RECORD_COUNT: usize = 8;

// ============================================================================
// Fixtures
// ============================================================================

/// `n` FASTQ records for stream `stream_idx`, named `read{i}` so all streams
/// agree on the name (which is what `ZipRawFastqK` / `ParseAndZipFastqN` join
/// on) while their sequences differ.
fn fastq_text(n: usize, stream_idx: usize) -> String {
    use std::fmt::Write as _;
    // One distinct base per stream, so a template that mixes streams (or drops
    // one) is visible in the assertion rather than hidden behind identical
    // sequences. Bases cycle every four streams (the 5-stream `Dyn` case reuses
    // stream 0's base for stream 4); `expected_templates` cycles identically, so
    // the oracle still pins each segment to its stream.
    let base = ['A', 'C', 'G', 'T'][stream_idx % 4];
    let mut out = String::new();
    for i in 0..n {
        writeln!(out, "@read{i}\n{base}{base}{base}{base}\n+\nIIII").expect("write to String");
    }
    out
}

fn boxed_reader(text: String) -> Box<dyn BufRead + Send> {
    Box::new(io::Cursor::new(text.into_bytes()))
}

/// SAM record lines (no header) for the `ReadSamChunks` fixture.
fn sam_record_lines(n: usize) -> String {
    let mut out = String::new();
    for i in 0..n {
        use std::fmt::Write as _;
        let pos = 100 + i;
        writeln!(out, "read{i}\t0\tchr1\t{pos}\t60\t4M\t*\t0\t0\tACGT\tIIII")
            .expect("write to String");
    }
    out
}

// ============================================================================
// Collecting sink
// ============================================================================

/// Terminal sink accumulating every item in arrival order. `Exclusive`, so
/// arrival order is the chain's output order with no sink-side interleaving.
struct CollectSink<T: Send + 'static> {
    collected: Arc<Mutex<Vec<T>>>,
}

impl<T: Send + HeapSize + 'static> Step for CollectSink<T> {
    type Input = T;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "CollectSink",
            kind: StepKind::Exclusive,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        match ctx.input.pop() {
            Some(item) => {
                self.collected.lock().expect("sink mutex not poisoned").push(item);
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }
}

/// Flatten collected template batches into `(name, [seq per stream])`, in
/// batch-serial order. Both FASTQ chains must produce the identical result,
/// which is what makes them comparable.
fn templates_in_order(batches: &[FastqTemplateBatch]) -> Vec<(Vec<u8>, Vec<Vec<u8>>)> {
    // Assert arrival order rather than sorting into it. The chain declares
    // `ByItemOrdinal` and terminates in an `Exclusive` sink, so batches must
    // ALREADY arrive in ordinal order — that is the property the reorder stage
    // exists to provide. Sorting here would silently repair a reorder
    // regression and leave the comparison below still passing.
    for pair in batches.windows(2) {
        assert!(
            pair[0].ordinal() < pair[1].ordinal(),
            "batches arrived out of order: ordinal {} before {}",
            pair[0].ordinal(),
            pair[1].ordinal(),
        );
    }
    batches
        .iter()
        .flat_map(|batch| {
            batch.templates.iter().map(|t| {
                (t.name.clone(), t.records.iter().map(|r| r.sequence().to_vec()).collect())
            })
        })
        .collect()
}

/// The templates the fixtures must yield: `read0..read{n-1}`, each pairing an
/// all-`A` R1 sequence with an all-`C` R2 sequence.
fn expected_templates(n: usize, n_streams: usize) -> Vec<(Vec<u8>, Vec<Vec<u8>>)> {
    (0..n)
        .map(|i| {
            let seqs = (0..n_streams)
                .map(|s| {
                    let base = [b'A', b'C', b'G', b'T'][s % 4];
                    vec![base; 4]
                })
                .collect();
            (format!("read{i}").into_bytes(), seqs)
        })
        .collect()
}

// ============================================================================
// Chains
// ============================================================================

/// Append the unified K-stream FASTQ join onto `tails` (one `FastqRawChunk`
/// producer edge per stream) and return the `FastqTemplateBatch`-emitting tail.
///
/// Mirrors `ChainBuilder::append_unified_fastq_join`: `WrapRawFastq1 →
/// ParseAndZipFastqN` for K==1, `ZipRawFastqK → ParseAndZipFastqN` for K>=2.
fn append_unified_join(
    builder: &PipelineBuilder,
    tails: &[(StepIdx, BranchIdx)],
) -> (StepIdx, BranchIdx) {
    let k = tails.len();
    let raw_tail = if k == 1 {
        builder.append_step(WrapRawFastq1::new(EDGE_LIMIT_BYTES), tails[0])
    } else {
        builder.append_step_k(ZipRawFastqK::new(k, EDGE_LIMIT_BYTES), tails)
    };
    builder.append_step(ParseAndZipFastqN::new(EDGE_LIMIT_BYTES), raw_tail)
}

/// The unified K-stream topology: one per-stream `ReadFastqInputs::new_single`
/// reader per producer edge, joined by `ZipRawFastqK` (K>=2) or `WrapRawFastq1`
/// (K==1), then parsed and zipped in the `Parallel` `ParseAndZipFastqN`.
///
/// Swept over `n_streams` including **1** (the `WrapRawFastq1` shim) and **3**
/// (`ZipRawFastqK` handling more than two streams — the whole point of the
/// K-way rewire). Every stream is its own edge with its OWN dense ordinal
/// sequence; cross-stream alignment is by `chunk_serial`.
///
/// Enough threads for K disjoint single-worker reader affinities plus the join
/// and sink: `n_streams + 2`, capped so the small cases still run at a couple
/// of threads.
#[rstest]
fn unified_fastq_chain_recovers_every_template_in_order(
    #[values(1, 8, 100)] n_records: usize,
    // 1 → WrapRawFastq1; 2/3/4 → KInputHandles::Fixed{2,3,4}; 5 → Dyn(Vec).
    #[values(1, 2, 3, 4, 5)] n_streams: usize,
) {
    for threads in [1usize, n_streams + 2] {
        let collected: Arc<Mutex<Vec<FastqTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
        let sink_handle = Arc::clone(&collected);

        let builder = PipelineBuilder::new();
        let mut tails = Vec::with_capacity(n_streams);
        for s in 0..n_streams {
            let worker = (threads - 1).min(s);
            let reader = ReadFastqInputs::new_single(
                boxed_reader(fastq_text(n_records, s)),
                s,
                FastqOrdinalSequence::new(),
                Affinity::Worker(worker),
                BATCH_RECORD_COUNT,
                EDGE_LIMIT_BYTES,
            );
            tails.push(builder.append_source(reader));
        }
        let joined = append_unified_join(&builder, &tails);
        builder.append_step(CollectSink { collected: sink_handle }, joined);

        let pipeline = builder.build().expect("chain builds");
        pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

        let batches = collected.lock().expect("mutex not poisoned");
        let templates = templates_in_order(&batches);
        assert_eq!(templates, expected_templates(n_records, n_streams));
        // Every template must carry exactly one record per stream — a dropped or
        // duplicated stream would otherwise be invisible if the sequences matched.
        assert!(
            templates.iter().all(|(_, seqs)| seqs.len() == n_streams),
            "every template must hold one record per stream",
        );
    }
}

/// Mismatched FASTQ inputs must **fail**, not hang.
///
/// `ZipRawFastqK` is lockstep — it pulls only the streams missing from the
/// front row — so a short stream leaves a permanently-incomplete front row.
/// When every input drains, the aligner reports "FASTQ sources out of sync"
/// rather than stalling the reorder stage. The assertion is that error, not
/// merely "does not hang".
///
/// Stream 0 must outlive stream 1 by several read cycles so the desync is a
/// genuine unequal-length mismatch, not a transient.
#[rstest]
fn mismatched_stream_lengths_error_rather_than_hanging(#[values(1, 4)] threads: usize) {
    let collected: Arc<Mutex<Vec<FastqTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);

    // BATCH_RECORD_COUNT records per chunk per cycle: stream 0 runs for several
    // cycles after stream 1 is exhausted.
    let long_records = BATCH_RECORD_COUNT * 4;
    let short_records = BATCH_RECORD_COUNT;

    let builder = PipelineBuilder::new();
    let r0_worker = 0;
    let r1_worker = (threads - 1).min(1);
    let r0 = builder.append_source(ReadFastqInputs::new_single(
        boxed_reader(fastq_text(long_records, 0)),
        0,
        FastqOrdinalSequence::new(),
        Affinity::Worker(r0_worker),
        BATCH_RECORD_COUNT,
        EDGE_LIMIT_BYTES,
    ));
    let r1 = builder.append_source(ReadFastqInputs::new_single(
        boxed_reader(fastq_text(short_records, 1)),
        1,
        FastqOrdinalSequence::new(),
        Affinity::Worker(r1_worker),
        BATCH_RECORD_COUNT,
        EDGE_LIMIT_BYTES,
    ));
    let joined = append_unified_join(&builder, &[r0, r1]);
    builder.append_step(CollectSink { collected: sink_handle }, joined);
    let pipeline = builder.build().expect("chain builds");

    let err = pipeline
        .run(PipelineConfig { threads, ..Default::default() })
        .expect_err("mismatched FASTQ stream lengths must fail the run");
    let message = format!("{err}");
    assert!(
        message.contains("out of sync"),
        "expected the aligner's out-of-sync diagnostic, got: {message}",
    );
}

/// BGZF-compress `text` into a concatenated block stream terminated by the
/// standard EOF marker — the byte layout `ReadFastqBlocks` reads from a bgzip'd
/// FASTQ file. `write_all` + `flush` + `take_blocks` is the compressor's own
/// documented drive loop; the EOF marker is appended so `read_raw_blocks` sees a
/// clean end of stream (an empty fixture is then just the lone marker).
fn bgzip(text: &str) -> Vec<u8> {
    let mut compressor = fgumi_bgzf::InlineBgzfCompressor::new(6);
    compressor.write_all(text.as_bytes()).expect("buffer FASTQ bytes");
    compressor.flush().expect("flush to BGZF blocks");
    let mut out = Vec::new();
    for block in compressor.take_blocks() {
        out.extend_from_slice(&block.data);
    }
    out.extend_from_slice(&fgumi_bgzf::BGZF_EOF);
    out
}

fn boxed_read(bytes: Vec<u8>) -> Box<dyn io::Read + Send> {
    Box::new(io::Cursor::new(bytes))
}

/// Append one stream's BGZF decode sub-chain (`ReadFastqBlocks →
/// FastqDecompress → FindFastqBoundaries`) and return its `FastqRawChunk` tail.
/// Each stream mints from its OWN `FastqOrdinalSequence` (its edge into the
/// K-way join is `ByItemOrdinal`, so its per-edge ordinals must be dense).
fn append_bgzf_stream(
    builder: &PipelineBuilder,
    bytes: Vec<u8>,
    stream_idx: usize,
    worker: usize,
) -> (StepIdx, BranchIdx) {
    use crate::pipeline::steps::source::fastq_bgzf::{FastqDecompress, ReadFastqBlocks};
    use crate::pipeline::steps::source::find_fastq_boundaries::FindFastqBoundaries;

    let read = builder.append_source(ReadFastqBlocks::new(
        boxed_read(bytes),
        stream_idx,
        Affinity::Worker(worker),
        EDGE_LIMIT_BYTES,
    ));
    let dec = builder.append_step(FastqDecompress::new(EDGE_LIMIT_BYTES, true), read);
    builder.append_step(
        FindFastqBoundaries::new(
            stream_idx,
            BATCH_RECORD_COUNT,
            FastqOrdinalSequence::new(),
            EDGE_LIMIT_BYTES,
        ),
        dec,
    )
}

/// The parallel BGZF decode split (`ReadFastqBlocks → FastqDecompress →
/// FindFastqBoundaries` → unified K-stream join) must recover byte-identical
/// templates, in order, to the fused `ReadFastqInputs` path — for K == 1
/// (`WrapRawFastq1`), K == 2, and K == 3 (`ZipRawFastqK`).
///
/// This is the correctness spine of the split: it separates raw block read
/// from a `Parallel` inflate and re-frames record seams in a distinct `Serial`
/// step, all of which must reassemble the original record stream exactly.
/// Compared against `expected_templates`, the same oracle the fused-path tests
/// use, so a seam/ordinal bug in the split cannot hide behind its own fixture.
#[rstest]
fn bgzf_split_fastq_chain_matches_fused_output(
    #[values(4, 8)] threads: usize,
    #[values(1, 8, 100)] n_records: usize,
    // 1 → WrapRawFastq1; 2/3/4 → KInputHandles::Fixed{2,3,4}; 5 → Dyn(Vec).
    #[values(1, 2, 3, 4, 5)] n_streams: usize,
) {
    let collected: Arc<Mutex<Vec<FastqTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);

    let builder = PipelineBuilder::new();
    let mut tails = Vec::with_capacity(n_streams);
    for s in 0..n_streams {
        let worker = (threads - 1).min(s);
        tails.push(append_bgzf_stream(&builder, bgzip(&fastq_text(n_records, s)), s, worker));
    }
    let joined = append_unified_join(&builder, &tails);
    builder.append_step(CollectSink { collected: sink_handle }, joined);

    let pipeline = builder.build().expect("split chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("split chain runs");

    let batches = collected.lock().expect("mutex not poisoned");
    assert_eq!(templates_in_order(&batches), expected_templates(n_records, n_streams));
}

/// An empty BGZF FASTQ stream must drain cleanly through the split — the raw
/// block reader sees only the EOF marker, so no block, no chunk, no template.
#[rstest]
fn bgzf_split_empty_stream_completes_with_no_templates(#[values(1, 4)] threads: usize) {
    let collected: Arc<Mutex<Vec<FastqTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);

    let builder = PipelineBuilder::new();
    let stream = append_bgzf_stream(&builder, bgzip(""), 0, 0);
    let joined = append_unified_join(&builder, &[stream]);
    builder.append_step(CollectSink { collected: sink_handle }, joined);

    let pipeline = builder.build().expect("split chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("split chain runs");

    let batches = collected.lock().expect("mutex not poisoned");
    let total: usize = batches.iter().map(|b| b.templates.len()).sum();
    assert_eq!(total, 0);
}

/// BGZF-compress `text` with ONE BGZF block per FASTQ record — a deliberately
/// different block layout from [`bgzip`] (which packs records into few blocks).
/// Used to prove the split's chunk cadence is set by record count, not block
/// boundaries: R1 (few blocks) and R2 (one-per-record) must still cut chunks in
/// lockstep for the `chunk_serial` join.
fn bgzip_one_block_per_record(text: &str) -> Vec<u8> {
    let mut out = Vec::new();
    // `text` is whole 4-line records; split on every 4th newline.
    let lines: Vec<&str> = text.split_inclusive('\n').collect();
    for rec_lines in lines.chunks(4) {
        let mut compressor = fgumi_bgzf::InlineBgzfCompressor::new(6);
        for l in rec_lines {
            compressor.write_all(l.as_bytes()).expect("buffer");
        }
        compressor.flush().expect("flush one record to its own block");
        for block in compressor.take_blocks() {
            out.extend_from_slice(&block.data);
        }
    }
    out.extend_from_slice(&fgumi_bgzf::BGZF_EOF);
    out
}

/// The split must stay in lockstep when R1 and R2 have DIFFERENT block layouts —
/// the exact desync hazard the fixed-record-count cut exists to prevent. R1 is
/// packed into few blocks; R2 is one block per record. With a small
/// `batch_record_count`, the drain path must emit multiple full chunks (spanning
/// several `try_run` calls via the held slot) before the short remainder, so the
/// two streams' `chunk_serial` streams line up and the join does not error.
///
/// This is the regression test for the review's B2/B3 finding (drain emitting one
/// oversized final chunk desynced the join). It runs `n_records` that is NOT a
/// multiple of the batch, so a genuine short final chunk exists on both streams.
#[rstest]
fn bgzf_split_mismatched_block_layouts_stay_in_lockstep(
    #[values(1, 4)] threads: usize,
    #[values(7, 50, 101)] n_records: usize,
) {
    let collected: Arc<Mutex<Vec<FastqTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);

    let builder = PipelineBuilder::new();

    // R1: packed blocks. R2: one block per record. Same records, same count.
    let r1_tail = append_bgzf_stream(&builder, bgzip(&fastq_text(n_records, 0)), 0, 0);
    // Clamp R2's worker exactly as the builder does: (threads-1).min(1), so at
    // threads=1 both readers pin to worker 0 (they serialize on the sole worker).
    let r2_worker = (threads - 1).min(1);
    let r2_tail = append_bgzf_stream(
        &builder,
        bgzip_one_block_per_record(&fastq_text(n_records, 1)),
        1,
        r2_worker,
    );

    let joined = append_unified_join(&builder, &[r1_tail, r2_tail]);
    builder.append_step(CollectSink { collected: sink_handle }, joined);

    let pipeline = builder.build().expect("split chain builds");
    pipeline
        .run(PipelineConfig { threads, ..Default::default() })
        .expect("mismatched block layouts must NOT desync the join");

    let batches = collected.lock().expect("mutex not poisoned");
    assert_eq!(templates_in_order(&batches), expected_templates(n_records, 2));
}

/// A bgzip'd FASTQ truncated mid-record (the last record cut short, not merely
/// missing its final newline) must FAIL the run with an error, not silently
/// emit a malformed record or hang. Regression for the review's fail-fast
/// parity finding — the split's `finalize_remainder` surfaces `UnexpectedEof`.
#[rstest]
fn bgzf_split_truncated_final_record_errors(#[values(1, 4)] threads: usize) {
    // Two whole records + a third cut off after its name+sequence lines.
    let mut text = fastq_text(2, 0);
    text.push_str("@r2\nACGT\n"); // truncated: only 2 of 4 lines, no +/qual

    let collected: Arc<Mutex<Vec<FastqTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let builder = PipelineBuilder::new();
    let stream = append_bgzf_stream(&builder, bgzip(&text), 0, 0);
    let joined = append_unified_join(&builder, &[stream]);
    builder.append_step(CollectSink { collected: sink_handle }, joined);

    let pipeline = builder.build().expect("split chain builds");
    let err = pipeline
        .run(PipelineConfig { threads, ..Default::default() })
        .expect_err("a truncated final record must fail the run, not emit garbage");
    let msg = format!("{err}");
    assert!(msg.contains("truncated"), "expected a truncation error, got: {msg}");
}

/// An empty FASTQ stream must drain cleanly and produce nothing — every step
/// sees drain before it ever sees an item.
#[rstest]
fn an_empty_fastq_stream_completes_with_no_templates(#[values(1, 4)] threads: usize) {
    let collected: Arc<Mutex<Vec<FastqTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);

    let builder = PipelineBuilder::new();
    let r0_worker = 0;
    let r1_worker = (threads - 1).min(1);
    let r0 = builder.append_source(ReadFastqInputs::new_single(
        boxed_reader(String::new()),
        0,
        FastqOrdinalSequence::new(),
        Affinity::Worker(r0_worker),
        BATCH_RECORD_COUNT,
        EDGE_LIMIT_BYTES,
    ));
    let r1 = builder.append_source(ReadFastqInputs::new_single(
        boxed_reader(String::new()),
        1,
        FastqOrdinalSequence::new(),
        Affinity::Worker(r1_worker),
        BATCH_RECORD_COUNT,
        EDGE_LIMIT_BYTES,
    ));
    let joined = append_unified_join(&builder, &[r0, r1]);
    builder.append_step(CollectSink { collected: sink_handle }, joined);
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let batches = collected.lock().expect("mutex not poisoned");
    let total: usize = batches.iter().map(|b| b.templates.len()).sum();
    assert_eq!(total, 0);
}

/// `ReadSamChunks` must emit every record line exactly once, split on line
/// boundaries, with a sentinel-form offset table per chunk. A tiny target chunk
/// size forces many chunks so the split path runs repeatedly rather than
/// emitting everything in one go.
#[rstest]
fn sam_chunk_source_emits_every_line_split_on_record_boundaries(
    #[values(1, 4)] threads: usize,
    #[values(1, 200)] n_records: usize,
) {
    let text = sam_record_lines(n_records);
    let collected: Arc<Mutex<Vec<SamChunk>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);

    let builder = Pipeline::builder();
    builder
        .chain(ReadSamChunks::new(
            boxed_reader(text.clone()),
            /* target_chunk_bytes */ 256,
            EDGE_LIMIT_BYTES,
        ))
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let chunks = collected.lock().expect("mutex not poisoned").drain(..).collect::<Vec<_>>();
    // Same reasoning as `templates_in_order`: assert the arrival order the
    // chain promises instead of sorting into it, so an out-of-order emission
    // fails here rather than being repaired before the byte comparison.
    for pair in chunks.windows(2) {
        assert!(
            pair[0].ordinal() < pair[1].ordinal(),
            "chunks arrived out of order: ordinal {} before {}",
            pair[0].ordinal(),
            pair[1].ordinal(),
        );
    }

    // Every chunk's offset table must be sentinel-form and describe complete
    // lines; concatenating them must reproduce the input exactly.
    let mut lines: Vec<String> = Vec::new();
    let mut rebuilt = Vec::new();
    for chunk in &chunks {
        rebuilt.extend_from_slice(&chunk.bytes);
        assert_eq!(
            *chunk.line_offsets.last().expect("sentinel-form table is never empty") as usize,
            chunk.bytes.len(),
            "the final offset must be the sentinel end-of-chunk",
        );
        for w in chunk.line_offsets.windows(2) {
            let line = &chunk.bytes[w[0] as usize..w[1] as usize];
            assert_eq!(line.last(), Some(&b'\n'), "every emitted line must be complete");
            lines.push(String::from_utf8(line.to_vec()).expect("utf8"));
        }
    }
    assert_eq!(rebuilt, text.as_bytes(), "chunks must reproduce the input byte-for-byte");
    assert_eq!(lines.len(), n_records, "one line per input record, no duplicates or drops");
    assert_eq!(lines.concat(), text);
}

/// An empty SAM body drains cleanly with no chunks carrying records.
#[test]
fn an_empty_sam_stream_completes_with_no_records() {
    let collected: Arc<Mutex<Vec<SamChunk>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);

    let builder = Pipeline::builder();
    builder
        .chain(ReadSamChunks::new(boxed_reader(String::new()), 256, EDGE_LIMIT_BYTES))
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
    builder
        .build()
        .expect("chain builds")
        .run(PipelineConfig { threads: 1, ..Default::default() })
        .expect("chain runs");

    let chunks = collected.lock().expect("mutex not poisoned");
    let total: usize = chunks.iter().map(SamChunk::record_count).sum();
    assert_eq!(total, 0);
}

/// `EDGE_LIMIT_BYTES` must be small enough that the multi-record fixtures
/// exceed it, or the byte-bounded edges never reject a push and every chain
/// above silently degrades into an unbounded-queue test.
///
/// This is a **necessary condition, not a proof of backpressure**: a fast
/// -draining downstream can keep occupancy under the limit even when the whole
/// fixture is far larger than one edge. The sibling test below closes the other
/// half — that an edge at this budget really does reject these chunk sizes.
#[test]
fn the_edge_budget_binds_on_the_fixtures() {
    let limit = usize::try_from(EDGE_LIMIT_BYTES).expect("budget fits usize");

    // One FASTQ stream must not fit in a single edge, so the reader is forced
    // to hold a rejected chunk and retry rather than emitting everything in one
    // pass. One chunk must still fit, or the edge would wedge instead.
    let stream_bytes = fastq_text(100, 0).len();
    let chunk_bytes = fastq_text(BATCH_RECORD_COUNT, 0).len();
    assert!(
        stream_bytes > limit,
        "one FASTQ stream ({stream_bytes} B) must exceed the per-edge budget \
         ({limit} B), else the byte-bounded edges never bind and these chains \
         prove nothing about backpressure",
    );
    assert!(
        chunk_bytes < limit,
        "a single chunk ({chunk_bytes} B) must still fit the budget ({limit} B), \
         else the edge admits nothing and the pipeline wedges",
    );

    // Same for the SAM fixture.
    let sam_bytes = sam_record_lines(200).len();
    assert!(sam_bytes > limit, "SAM fixture ({sam_bytes} B) must exceed the budget ({limit} B)");
}

/// A byte-bounded edge at exactly `EDGE_LIMIT_BYTES`, fed the chunk sizes the
/// FASTQ chains put on it, must reject before the fixture is exhausted — so the
/// held-item retry path in `ReadFastqInputs` is genuinely reachable at this
/// budget, not merely reachable in principle.
///
/// What this still does not assert is that a rejection occurred *during the
/// chain runs above*. Observing that needs the per-edge reject counter, which
/// `fgumi-pipeline-core` records (`EdgeMetrics::record_reject`) but keeps
/// `pub(crate)`, so it is unreachable from here without widening that crate's
/// API — out of scope for this PR. The chains cover the retry path indirectly:
/// every record arrives, in order, across an edge this test shows must reject,
/// which could not happen if held items were dropped or never retried. That is
/// an inference from the end state, not a direct assertion, and it is recorded
/// here so the gap is visible rather than assumed closed.
#[test]
fn the_byte_bounded_edge_rejects_at_the_test_budget() {
    use crate::pipeline::core::queues::{ByteBoundedQueue, ItemQueue};
    use crate::pipeline::steps::source::read_fastq::FastqRawChunk;

    let queue = ByteBoundedQueue::<FastqRawChunk>::new(EDGE_LIMIT_BYTES);
    let chunk_data = fastq_text(BATCH_RECORD_COUNT, 0).into_bytes();
    let chunk_bytes = chunk_data.len();

    // Push identical chunks until one is refused. The fixture spans many such
    // chunks, so a queue that never rejects would loop past the whole stream.
    let max_pushes = fastq_text(100, 0).len() / chunk_bytes + 1;
    let mut accepted = 0usize;
    let mut rejected = false;
    for ordinal in 0..max_pushes {
        let chunk = FastqRawChunk {
            ordinal: ordinal as u64,
            stream_idx: 0,
            chunk_serial: ordinal as u64,
            data: chunk_data.clone(),
        };
        if queue.try_push(chunk).is_err() {
            rejected = true;
            break;
        }
        accepted += 1;
    }

    assert!(
        rejected,
        "a {EDGE_LIMIT_BYTES}-byte edge admitted all {max_pushes} chunks of {chunk_bytes} B \
         without rejecting; the byte budget is not binding and the chains above prove \
         nothing about the held-item retry path",
    );
    assert!(
        accepted >= 1,
        "the edge must admit at least one chunk before rejecting, else the chains wedge \
         instead of exercising backpressure",
    );
}
