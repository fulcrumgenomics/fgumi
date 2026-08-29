//! End-to-end tests for the read-side source chains.
//!
//! The source steps are `try_run` state machines wrapped around real readers,
//! so their interesting behavior — round-robin stream rotation, chunk
//! alignment to whole 4-line records, per-stream ordinal assignment, drain
//! reporting, and held-item retry under backpressure — only happens when the
//! framework drives them against an actual stream. These tests assemble the
//! two production FASTQ topologies and the SAM topology and run them through
//! `Pipeline::run`:
//!
//! ```text
//! ReadFastqInputs (all streams)  → ParseFastqChunks → ZipFastqRecords → sink
//! ReadFastqInputs (one per stream) ×2 → PairRawFastq → ParseAndZipFastq → sink
//! ReadSamChunks → sink
//! ```
//!
//! The two FASTQ chains are not redundant: the first is the N>=3 fallback
//! (a single `Affinity::Reader` worker rotating across every stream), the
//! second is the paired-end path (one reader per stream on disjoint workers,
//! joined two-way). They assign ordinals by different rules, so a bug in
//! either is invisible to the other. Both must recover the same templates in
//! the same order from the same input.
//!
//! `pair_fastq` and `zip_fastq` carry their own targeted pipeline tests for
//! backpressure and stream skew; these cover the plain read-through path those
//! tests stub out with synthetic sources.

use std::io::{self, BufRead};
use std::sync::{Arc, Mutex};

use rstest::rstest;

use crate::pipeline::core::builder::{Pipeline, PipelineBuilder, PipelineConfig};
use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::pipeline::core::step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::source::pair_fastq::PairRawFastq;
use crate::pipeline::steps::source::parse_fastq::ParseFastqChunks;
use crate::pipeline::steps::source::parse_zip_fastq::ParseAndZipFastq;
use crate::pipeline::steps::source::read_fastq::{FastqOrdinalSequence, ReadFastqInputs};
use crate::pipeline::steps::source::read_sam_chunks::ReadSamChunks;
use crate::pipeline::steps::source::zip_fastq::{FastqTemplateBatch, ZipFastqRecords};
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

/// `n` FASTQ records for stream `stream_idx`, named `read{i}` so both streams
/// of a pair agree on the name (which is what `ZipFastqRecords` /
/// `ParseAndZipFastq` join on) while their sequences differ.
fn fastq_text(n: usize, stream_idx: usize) -> String {
    use std::fmt::Write as _;
    // One distinct base per stream, so a template that mixes streams (or drops
    // one) is visible in the assertion rather than hidden behind identical
    // sequences. Cycles past four streams, which no fixture here reaches.
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

/// The N-stream fallback topology: one `Affinity::Reader` source rotating across
/// every stream, then parallel parse, then the serial zip join.
///
/// Swept over `n_streams` including **3**, which is the case this topology
/// exists for — `PairRawFastq` handles exactly two, so anything above that
/// routes here. A 2-stream-only test would leave the rotation across a third
/// stream, and its drain, entirely unexercised.
#[rstest]
fn round_robin_fastq_chain_recovers_every_template_in_order(
    #[values(1, 4)] threads: usize,
    #[values(1, 8, 100)] n_records: usize,
    #[values(2, 3)] n_streams: usize,
) {
    let collected: Arc<Mutex<Vec<FastqTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);

    let readers: Vec<Box<dyn BufRead + Send>> =
        (0..n_streams).map(|s| boxed_reader(fastq_text(n_records, s))).collect();

    let builder = Pipeline::builder();
    builder
        .chain(ReadFastqInputs::new(readers, BATCH_RECORD_COUNT, EDGE_LIMIT_BYTES))
        .chain(ParseFastqChunks::new(EDGE_LIMIT_BYTES))
        .chain(ZipFastqRecords::new(n_streams, EDGE_LIMIT_BYTES))
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
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

/// Mismatched FASTQ inputs must **fail**, not hang.
///
/// This is the regression test for the ordinal-density bug. The reader used to
/// derive ordinals as `chunk_serial * n_streams_total + stream_idx`, giving each
/// stream its own residue class. That is unique but only dense when every
/// stream yields the same number of chunks — so a short R2 left a permanently
/// unfilled ordinal, `BranchOrdering::ByItemOrdinal` stalled waiting for it, and
/// the run hung.
///
/// The damage was not just the stall: it *masked the correct diagnostic*.
/// `ZipFastqRecords` already detects this input and reports "FASTQ sources out
/// of sync", but the stalled reorder stage meant those chunks never reached it.
/// So the assertion here is the error, not merely "does not hang" — the fix is
/// what lets the existing check fire.
///
/// Stream 0 must outlive stream 1 by at least two read cycles, since the gap
/// only strands a later ordinal if one is minted *after* it.
#[rstest]
fn mismatched_stream_lengths_error_rather_than_hanging(#[values(1, 4)] threads: usize) {
    let collected: Arc<Mutex<Vec<FastqTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);

    // BATCH_RECORD_COUNT records per chunk per cycle: stream 0 runs for several
    // cycles after stream 1 is exhausted.
    let long_records = BATCH_RECORD_COUNT * 4;
    let short_records = BATCH_RECORD_COUNT;
    let readers: Vec<Box<dyn BufRead + Send>> =
        vec![boxed_reader(fastq_text(long_records, 0)), boxed_reader(fastq_text(short_records, 1))];

    let builder = Pipeline::builder();
    builder
        .chain(ReadFastqInputs::new(readers, BATCH_RECORD_COUNT, EDGE_LIMIT_BYTES))
        .chain(ParseFastqChunks::new(EDGE_LIMIT_BYTES))
        .chain(ZipFastqRecords::new(2, EDGE_LIMIT_BYTES))
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");

    let err = pipeline
        .run(PipelineConfig { threads, ..Default::default() })
        .expect_err("mismatched FASTQ stream lengths must fail the run");
    let message = format!("{err}");
    assert!(
        message.contains("out of sync"),
        "expected the zipper's out-of-sync diagnostic, got: {message}",
    );
}

/// The paired-end topology: one reader per stream pinned to a distinct worker,
/// joined two-way by `PairRawFastq`, then parsed and zipped in one step.
///
/// Needs >= 3 threads: R1 and R2 hold disjoint single-worker affinities and a
/// third worker drives the join and sink.
#[rstest]
fn per_stream_paired_fastq_chain_recovers_every_template_in_order(
    #[values(3, 5)] threads: usize,
    #[values(1, 8, 100)] n_records: usize,
) {
    let collected: Arc<Mutex<Vec<FastqTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);

    // One sequence, cloned into both readers — that shared counter is what
    // keeps the combined ordinal stream gap-free.
    let ordinals = FastqOrdinalSequence::new();
    let r1 = ReadFastqInputs::new_single(
        boxed_reader(fastq_text(n_records, 0)),
        /* global_stream_idx */ 0,
        ordinals.clone(),
        Affinity::Worker(0),
        BATCH_RECORD_COUNT,
        EDGE_LIMIT_BYTES,
    );
    let r2 = ReadFastqInputs::new_single(
        boxed_reader(fastq_text(n_records, 1)),
        1,
        ordinals,
        Affinity::Worker(1),
        BATCH_RECORD_COUNT,
        EDGE_LIMIT_BYTES,
    );

    let builder = PipelineBuilder::new();
    let r1_tail = builder.append_source(r1);
    let r2_tail = builder.append_source(r2);
    let paired = builder.append_step2(PairRawFastq::new(EDGE_LIMIT_BYTES), r1_tail, r2_tail);
    let zipped = builder.append_step(ParseAndZipFastq::new(EDGE_LIMIT_BYTES), paired);
    builder.append_step(CollectSink { collected: sink_handle }, zipped);

    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let batches = collected.lock().expect("mutex not poisoned");
    assert_eq!(templates_in_order(&batches), expected_templates(n_records, 2));
}

/// Both FASTQ topologies draw ordinals from one shared `FastqOrdinalSequence`;
/// they differ only in *who mints* — the round-robin reader is the single
/// instance minting for every stream, while the per-stream readers each hold a
/// clone of the same sequence. Reading the same input through both must
/// therefore recover byte-identical templates in the same order; a wiring
/// mistake that gave one topology its own counter is invisible to that
/// topology's own test, because a per-instance sequence still looks dense from
/// inside a single reader.
#[test]
fn both_fastq_topologies_agree_on_the_same_input() {
    const N: usize = 64;

    let round_robin = {
        let collected: Arc<Mutex<Vec<FastqTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
        let sink_handle = Arc::clone(&collected);
        let builder = Pipeline::builder();
        builder
            .chain(ReadFastqInputs::new(
                vec![boxed_reader(fastq_text(N, 0)), boxed_reader(fastq_text(N, 1))],
                BATCH_RECORD_COUNT,
                EDGE_LIMIT_BYTES,
            ))
            .chain(ParseFastqChunks::new(EDGE_LIMIT_BYTES))
            .chain(ZipFastqRecords::new(2, EDGE_LIMIT_BYTES))
            .chain(CollectSink { collected: sink_handle })
            .into_sink_marker();
        builder
            .build()
            .expect("builds")
            .run(PipelineConfig { threads: 4, ..Default::default() })
            .expect("runs");
        let batches = collected.lock().expect("mutex not poisoned");
        templates_in_order(&batches)
    };

    let per_stream = {
        let collected: Arc<Mutex<Vec<FastqTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
        let sink_handle = Arc::clone(&collected);
        let builder = PipelineBuilder::new();
        let ordinals = FastqOrdinalSequence::new();
        let a = builder.append_source(ReadFastqInputs::new_single(
            boxed_reader(fastq_text(N, 0)),
            0,
            ordinals.clone(),
            Affinity::Worker(0),
            BATCH_RECORD_COUNT,
            EDGE_LIMIT_BYTES,
        ));
        let b = builder.append_source(ReadFastqInputs::new_single(
            boxed_reader(fastq_text(N, 1)),
            1,
            ordinals,
            Affinity::Worker(1),
            BATCH_RECORD_COUNT,
            EDGE_LIMIT_BYTES,
        ));
        let paired = builder.append_step2(PairRawFastq::new(EDGE_LIMIT_BYTES), a, b);
        let zipped = builder.append_step(ParseAndZipFastq::new(EDGE_LIMIT_BYTES), paired);
        builder.append_step(CollectSink { collected: sink_handle }, zipped);
        builder
            .build()
            .expect("builds")
            .run(PipelineConfig { threads: 4, ..Default::default() })
            .expect("runs");
        let batches = collected.lock().expect("mutex not poisoned");
        templates_in_order(&batches)
    };

    assert_eq!(round_robin, expected_templates(N, 2));
    assert_eq!(
        round_robin, per_stream,
        "the round-robin and per-stream reader topologies must agree on identical input",
    );
}

/// An empty FASTQ stream must drain cleanly and produce nothing — every step
/// sees drain before it ever sees an item.
#[rstest]
fn an_empty_fastq_stream_completes_with_no_templates(#[values(1, 4)] threads: usize) {
    let collected: Arc<Mutex<Vec<FastqTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);

    let builder = Pipeline::builder();
    builder
        .chain(ReadFastqInputs::new(
            vec![boxed_reader(String::new()), boxed_reader(String::new())],
            BATCH_RECORD_COUNT,
            EDGE_LIMIT_BYTES,
        ))
        .chain(ParseFastqChunks::new(EDGE_LIMIT_BYTES))
        .chain(ZipFastqRecords::new(2, EDGE_LIMIT_BYTES))
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
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
