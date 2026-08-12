//! Cross-step integration tests for the step subset in this port.
//!
//! Every step in `steps/` is a `try_run` state machine whose interesting
//! behavior — held-item retry, drain detection, `Finished` reporting — only
//! occurs when the framework drives it. Unit tests over a step's helpers
//! therefore leave the step body itself unexercised, so these tests assemble
//! real chains and run them through [`Pipeline::run`]:
//!
//! ```text
//! ReplaySource → BgzfDecompress → FindBamBoundaries → DecodeRecords    → sink
//! ReplaySource → BgzfDecompress → FindBamBoundaries → ParseBamRecords  → sink
//! ReplaySource → BgzfDecompress → FindBamBoundaries → BgzfCompress     → sink
//! ReplaySource → TemplatesToRecordBatch                                → sink
//! ReplaySource → GroupByMi                                             → sink
//! ReplaySource → Process2Ordered                          → sink A + sink B
//! ```
//!
//! `ReplaySource` is a test-local source because the real `source/` steps
//! are not part of this port; the block bytes it replays are built by the same
//! BGZF compressor the production writer uses, so the decode path sees exactly
//! what it would in production.
//!
//! Each chain runs at 1 and 4 threads. That is not redundancy: at 1 thread the
//! `Parallel` steps have a single clone and drain trivially, while at 4 the
//! last-worker barrier (only the final clone may close a shared output) is
//! actually exercised. A drain bug that wedges the pipeline shows up as a
//! hang in the 4-thread case only.
//!
//! The full production chain (through `GroupBam`, serialization, and the file
//! sink) is covered by the integration tests that arrive with the remaining
//! step subtrees.

use std::collections::VecDeque;
use std::io;
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
use std::sync::{Arc, Mutex};

use fgumi_bam_io::GroupKeyConfig;
use fgumi_raw_bam::{RawRecord, SamBuilder};
use rstest::rstest;

use crate::grouper::ProcessedPositionGroup;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::builder::{Pipeline, PipelineConfig};
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::bgzf::compress::BgzfCompress;
use crate::pipeline::steps::bgzf::decompress::BgzfDecompress;
use crate::pipeline::steps::boundaries::bam::FindBamBoundaries;
use crate::pipeline::steps::group::position::BatchedProcessedPositionGroups;
use crate::pipeline::steps::parse::bam::ParseBamRecords;
use crate::pipeline::steps::parse::decode::{DecodeFromRecords, DecodeRecords};
use crate::pipeline::steps::types::{BgzfBlock, DecodedRecordBatch, RecordBatch};
use crate::template::Template;

/// Per-edge byte budget. Deliberately small relative to the multi-record
/// fixtures so the byte-bounded edges bind mid-run rather than swallowing the
/// whole stream in one go. `the_edge_budget_binds_on_the_multi_record_fixtures`
/// pins that this constant actually has teeth — a future edit that raised it
/// past the fixture volume would silently turn every chain below into an
/// unbounded-queue test.
const EDGE_LIMIT_BYTES: u64 = 16 * 1024;

// ============================================================================
// Fixtures
// ============================================================================

/// A minimal binary BAM header: magic, empty text, and `n_ref` references.
/// This is the prefix `FindBamBoundaries` must consume and discard.
fn minimal_binary_bam_header(n_ref: u32) -> Vec<u8> {
    let mut header = Vec::new();
    header.extend_from_slice(fgumi_raw_bam::BAM_MAGIC);
    header.extend_from_slice(&0u32.to_le_bytes()); // l_text = 0
    header.extend_from_slice(&n_ref.to_le_bytes());
    for _ in 0..n_ref {
        header.extend_from_slice(&2u32.to_le_bytes()); // l_name = 2 ("r" + NUL)
        header.extend_from_slice(b"r\0");
        header.extend_from_slice(&1_000_000u32.to_le_bytes()); // l_ref
    }
    header
}

/// `n` mapped, paired records with distinct read names and positions — enough
/// field variety that a group key computed from them is not trivially constant.
fn test_records(n: usize) -> Vec<RawRecord> {
    (0..n)
        .map(|i| {
            let mut b = SamBuilder::new();
            b.read_name(format!("read{i}").as_bytes())
                .flags(fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::FIRST_SEGMENT)
                .ref_id(0)
                .pos(i32::try_from(100 + i * 10).expect("pos fits i32"))
                .cigar_ops(&[4u32 << 4]) // 4M
                .sequence(b"ACGT")
                .qualities(&[30u8; 4]);
            b.build()
        })
        .collect()
}

/// Serialize a header + records into a BAM byte stream, then cut it into BGZF
/// blocks of at most `payload_bytes` uncompressed bytes each.
///
/// Cutting at a fixed byte size (rather than on record boundaries) is the
/// point: records land straddling block edges, which is exactly the carryover
/// case `FindBamBoundaries` exists to handle.
fn bgzf_blocks_for(records: &[RawRecord], payload_bytes: usize) -> Vec<BgzfBlock> {
    let mut stream = minimal_binary_bam_header(1);
    for record in records {
        let bytes = record.as_ref();
        let block_size = u32::try_from(bytes.len()).expect("record fits u32");
        stream.extend_from_slice(&block_size.to_le_bytes());
        stream.extend_from_slice(bytes);
    }
    stream
        .chunks(payload_bytes)
        .enumerate()
        .map(|(i, payload)| {
            let mut compressor = fgumi_bgzf::InlineBgzfCompressor::new(1);
            compressor.write_all(payload).expect("compress payload");
            compressor.flush().expect("flush compressor");
            let mut blocks = compressor.take_blocks();
            assert_eq!(blocks.len(), 1, "payload must fit one BGZF block");
            BgzfBlock {
                batch_serial: i as u64,
                bytes: blocks.remove(0).data,
                uncompressed_size: u32::try_from(payload.len()).expect("payload fits u32"),
            }
        })
        .collect()
}

/// A four-record template for one queryname: primary R1 and R2 plus an R1
/// supplementary and an R2 secondary.
///
/// The split alignments are the point — a flatten written against
/// [`Template::r1`]/[`Template::r2`] would silently drop them, and every record
/// here carries a distinct position, sequence, and length so any dropped,
/// duplicated, or reordered record changes the emitted byte stream.
fn split_alignment_template(qname: &[u8]) -> Template {
    use fgumi_raw_bam::flags::{FIRST_SEGMENT, LAST_SEGMENT, PAIRED, SECONDARY, SUPPLEMENTARY};

    let record = |flags: u16, pos: i32, seq: &[u8], quals: &[u8]| {
        let mut b = SamBuilder::new();
        b.read_name(qname)
            .flags(flags)
            .ref_id(0)
            .pos(pos)
            // (length << 4) | op, where op 0 is `M`.
            .cigar_ops(&[u32::try_from(seq.len()).expect("seq length fits u32") << 4])
            .sequence(seq)
            .qualities(quals);
        b.build()
    };

    Template::from_records(vec![
        record(PAIRED | FIRST_SEGMENT, 100, b"ACGT", &[30u8; 4]),
        record(PAIRED | LAST_SEGMENT, 200, b"TTGCA", &[31u8; 5]),
        record(PAIRED | FIRST_SEGMENT | SUPPLEMENTARY, 300, b"GGGCCC", &[32u8; 6]),
        record(PAIRED | LAST_SEGMENT | SECONDARY, 400, b"TATATATA", &[33u8; 8]),
    ])
    .expect("split-alignment template")
}

/// A minimal record carrying `MI:Z:<mi>` — the only field `GroupByMi` reads.
fn mi_record(mi: &str) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(b"read")
        .flags(0)
        .sequence(b"ACGT")
        .qualities(&[30u8; 4])
        .add_string_tag(*crate::sam::SamTag::MI, mi.as_bytes());
    b.build()
}

/// A minimal record with no `MI` tag at all — dropped by `GroupByMi` and
/// counted against `skipped_no_mi`.
fn record_without_mi() -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(b"read").flags(0).sequence(b"ACGT").qualities(&[30u8; 4]);
    b.build()
}

/// A minimal record carrying `MI:i:<value>` — an integer where a string is
/// required, so it groups nothing and is counted against
/// `skipped_non_string_mi`.
fn record_with_integer_mi(value: i32) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(b"read")
        .flags(0)
        .sequence(b"ACGT")
        .qualities(&[30u8; 4])
        .add_int_tag(*crate::sam::SamTag::MI, value);
    b.build()
}

/// Wrap raw records in a [`DecodedRecordBatch`]. The [`GroupKey`] is left at
/// its default because `GroupByMi` re-reads the MI tag from the record bytes
/// and never consults the pre-computed key.
fn decoded_batch(batch_serial: u64, records: Vec<RawRecord>) -> DecodedRecordBatch {
    DecodedRecordBatch::new(
        batch_serial,
        records
            .into_iter()
            .map(|raw| {
                fgumi_bam_io::DecodedRecord::from_raw_bytes(raw, fgumi_bam_io::GroupKey::default())
            })
            .collect(),
    )
}

// ============================================================================
// Test-local source and sink steps
// ============================================================================

/// Replays pre-built items into the chain. `Exclusive` because it owns mutable
/// cursor state and must not be cloned per worker.
///
/// Holds a rejected item and retries it on a later tick rather than dropping
/// it — dropping would punch a hole in the ordinal sequence, which the
/// downstream `ByItemOrdinal` reorder stages would then wait on forever.
struct ReplaySource<T: Send + HeapSize + Ordered + 'static> {
    items: VecDeque<T>,
    held: HeldSlot<Unpushed<T>>,
}

impl<T: Send + HeapSize + Ordered + 'static> ReplaySource<T> {
    fn new(items: Vec<T>) -> Self {
        Self { items: items.into(), held: HeldSlot::new() }
    }
}

impl<T: Send + HeapSize + Ordered + 'static> Step for ReplaySource<T> {
    type Input = ();
    type Outputs = OrderedBytesSingle<T>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ReplaySource",
            kind: StepKind::Exclusive,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: EDGE_LIMIT_BYTES }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }
        let Some(item) = self.items.pop_front() else {
            return Ok(StepOutcome::Finished);
        };
        if let Err(unpushed) = ctx.outputs.push(item) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }
}

/// Terminal sink that accumulates every item it receives, in arrival order.
/// `Exclusive` so arrival order is the chain's output order with no
/// sink-side interleaving to reason about.
struct CollectSink<T: Send + 'static> {
    collected: Arc<Mutex<Vec<T>>>,
}

impl<T> Step for CollectSink<T>
where
    T: Send + HeapSize + 'static,
{
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

/// Like [`CollectSink`], but ignores its input for the first `remaining_stalls`
/// calls.
///
/// Stalling the terminal step is what lets an upstream byte-bounded output
/// queue actually fill up. With a sink that drains on every tick the queue
/// returns to zero bytes between pushes, so `cur >= limit` never holds at push
/// time and the producer's held-slot retry path is never taken — at one thread
/// it is structurally unreachable, and at four it depends on scheduling. That
/// makes a "this exercises backpressure" claim untrue exactly when the test
/// looks like it passes.
struct StallingCollectSink<T: Send + 'static> {
    collected: Arc<Mutex<Vec<T>>>,
    remaining_stalls: usize,
}

impl<T> Step for StallingCollectSink<T>
where
    T: Send + HeapSize + 'static,
{
    type Input = T;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "StallingCollectSink",
            kind: StepKind::Exclusive,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if self.remaining_stalls > 0 {
            self.remaining_stalls -= 1;
            return Ok(StepOutcome::NoProgress);
        }
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

// ============================================================================
// Chains
// ============================================================================

/// `ReplaySource → BgzfDecompress → FindBamBoundaries → DecodeRecords`.
///
/// Every record in the input must arrive downstream exactly once, in order,
/// with its bytes intact — across BGZF block boundaries, byte backpressure,
/// and (at 4 threads) parallel decompress/decode workers.
#[rstest]
fn decode_chain_round_trips_every_record_in_order(
    #[values(1, 4)] threads: usize,
    #[values(1, 64, 500)] n_records: usize,
) {
    let records = test_records(n_records);
    let blocks = bgzf_blocks_for(&records, 4096);
    let collected: Arc<Mutex<Vec<DecodedRecordBatch>>> = Arc::new(Mutex::new(Vec::new()));

    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(blocks))
        .chain(BgzfDecompress::new(EDGE_LIMIT_BYTES))
        .chain(FindBamBoundaries::new(EDGE_LIMIT_BYTES))
        .chain(DecodeRecords::new(GroupKeyConfig::default(), EDGE_LIMIT_BYTES))
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let batches = collected.lock().expect("mutex not poisoned");
    let seen: Vec<&[u8]> = batches
        .iter()
        .flat_map(|batch| batch.records().iter().map(fgumi_bam_io::DecodedRecord::raw_bytes))
        .collect();
    let expected: Vec<&[u8]> = records.iter().map(AsRef::as_ref).collect();
    assert_eq!(seen, expected, "records must survive the chain intact and in order");
}

/// Same chain but terminating in `ParseBamRecords`, the zero-alloc sibling of
/// `DecodeRecords` that emits shared-buffer `RecordBatch`es instead of owned
/// `DecodedRecord`s. It must yield the same records in the same order.
#[rstest]
fn parse_chain_round_trips_every_record_in_order(
    #[values(1, 4)] threads: usize,
    #[values(1, 64, 500)] n_records: usize,
) {
    let records = test_records(n_records);
    let blocks = bgzf_blocks_for(&records, 4096);
    let collected: Arc<Mutex<Vec<RecordBatch>>> = Arc::new(Mutex::new(Vec::new()));

    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(blocks))
        .chain(BgzfDecompress::new(EDGE_LIMIT_BYTES))
        .chain(FindBamBoundaries::new(EDGE_LIMIT_BYTES))
        .chain(ParseBamRecords::new(EDGE_LIMIT_BYTES))
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let batches = collected.lock().expect("mutex not poisoned");
    let seen: Vec<Vec<u8>> =
        batches.iter().flat_map(|batch| batch.iter_record_bytes().map(<[u8]>::to_vec)).collect();
    let expected: Vec<Vec<u8>> = records.iter().map(|r| r.as_ref().to_vec()).collect();
    assert_eq!(seen, expected, "records must survive the parse chain intact and in order");
}

/// `… → FindBamBoundaries → BgzfCompress` re-compresses the record stream.
/// Decompressing the sink's output must reproduce the record bytes exactly —
/// proving the compress step preserves both content and batch ordering.
#[rstest]
fn compress_chain_output_decompresses_back_to_the_record_stream(#[values(1, 4)] threads: usize) {
    let records = test_records(200);
    let blocks = bgzf_blocks_for(&records, 4096);
    let collected: Arc<Mutex<Vec<BgzfBlock>>> = Arc::new(Mutex::new(Vec::new()));

    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(blocks))
        .chain(BgzfDecompress::new(EDGE_LIMIT_BYTES))
        .chain(FindBamBoundaries::new(EDGE_LIMIT_BYTES))
        .chain(BgzfCompress::new(1, EDGE_LIMIT_BYTES))
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    // Concatenate the emitted BGZF blocks and inflate the whole stream.
    let compressed: Vec<u8> = collected
        .lock()
        .expect("mutex not poisoned")
        .iter()
        .flat_map(|b| b.bytes.clone())
        .collect();
    // Read in bounded batches: `read_raw_blocks` pre-allocates `max_blocks`,
    // so it must be given a sane cap rather than "everything".
    let mut cursor = io::Cursor::new(compressed);
    let mut decompressor = libdeflater::Decompressor::new();
    let mut inflated = Vec::new();
    loop {
        let raw_blocks = fgumi_bgzf::reader::read_raw_blocks(&mut cursor, 64)
            .expect("output is a valid BGZF stream");
        if raw_blocks.is_empty() {
            break;
        }
        for block in &raw_blocks {
            inflated.extend_from_slice(
                &fgumi_bgzf::reader::decompress_block(block, &mut decompressor)
                    .expect("each emitted block inflates"),
            );
        }
    }

    // The compressed stream is the header-stripped record stream: framed
    // `block_size` + body per record, in input order.
    let mut expected = Vec::new();
    for record in &records {
        let bytes = record.as_ref();
        expected.extend_from_slice(&u32::try_from(bytes.len()).expect("fits u32").to_le_bytes());
        expected.extend_from_slice(bytes);
    }
    assert_eq!(inflated, expected, "compress must preserve the record stream byte-for-byte");
}

/// An input with a header but no records must still complete cleanly and
/// produce nothing — the empty-BAM case, where every step sees drain before it
/// ever sees an item.
#[rstest]
fn a_header_only_input_completes_with_no_records(#[values(1, 4)] threads: usize) {
    let blocks = bgzf_blocks_for(&[], 4096);
    let collected: Arc<Mutex<Vec<DecodedRecordBatch>>> = Arc::new(Mutex::new(Vec::new()));

    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(blocks))
        .chain(BgzfDecompress::new(EDGE_LIMIT_BYTES))
        .chain(FindBamBoundaries::new(EDGE_LIMIT_BYTES))
        .chain(DecodeRecords::new(GroupKeyConfig::default(), EDGE_LIMIT_BYTES))
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let batches = collected.lock().expect("mutex not poisoned");
    let total: usize = batches.iter().map(|b| b.records().len()).sum();
    assert_eq!(total, 0, "a header-only BAM yields no records");
}

/// A single BGZF block holding the entire stream exercises the path where
/// `FindBamBoundaries` never carries anything over — the whole header and every
/// record arrive together in one input.
#[rstest]
fn a_single_block_input_needs_no_carryover(#[values(1, 4)] threads: usize) {
    let records = test_records(50);
    // One block big enough for the whole stream.
    let blocks = bgzf_blocks_for(&records, 64 * 1024);
    assert_eq!(blocks.len(), 1, "fixture must fit in a single BGZF block");

    let collected: Arc<Mutex<Vec<DecodedRecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(blocks))
        .chain(BgzfDecompress::new(EDGE_LIMIT_BYTES))
        .chain(FindBamBoundaries::new(EDGE_LIMIT_BYTES))
        .chain(DecodeRecords::new(GroupKeyConfig::default(), EDGE_LIMIT_BYTES))
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let batches = collected.lock().expect("mutex not poisoned");
    let total: usize = batches.iter().map(|b| b.records().len()).sum();
    assert_eq!(total, records.len());
}

/// `FindBamBoundaries::new_no_header` is the runall-spliced entry point: the
/// stream is already past the header, so byte 0 is a record boundary. Feeding
/// it a header-stripped stream must yield every record.
#[rstest]
fn the_no_header_boundary_variant_consumes_a_header_stripped_stream(
    #[values(1, 4)] threads: usize,
) {
    let records = test_records(64);
    // Build blocks over the record stream ONLY — no BAM header prefix.
    let mut stream = Vec::new();
    for record in &records {
        let bytes = record.as_ref();
        stream.extend_from_slice(&u32::try_from(bytes.len()).expect("fits u32").to_le_bytes());
        stream.extend_from_slice(bytes);
    }
    let blocks: Vec<BgzfBlock> = stream
        .chunks(1024)
        .enumerate()
        .map(|(i, payload)| {
            let mut compressor = fgumi_bgzf::InlineBgzfCompressor::new(1);
            compressor.write_all(payload).expect("compress");
            compressor.flush().expect("flush");
            let mut blocks = compressor.take_blocks();
            assert_eq!(blocks.len(), 1);
            BgzfBlock {
                batch_serial: i as u64,
                bytes: blocks.remove(0).data,
                uncompressed_size: u32::try_from(payload.len()).expect("fits u32"),
            }
        })
        .collect();

    let collected: Arc<Mutex<Vec<DecodedRecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(blocks))
        .chain(BgzfDecompress::new(EDGE_LIMIT_BYTES))
        .chain(FindBamBoundaries::new_no_header(EDGE_LIMIT_BYTES))
        .chain(DecodeRecords::new(GroupKeyConfig::default(), EDGE_LIMIT_BYTES))
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let batches = collected.lock().expect("mutex not poisoned");
    let seen: Vec<&[u8]> = batches
        .iter()
        .flat_map(|batch| batch.records().iter().map(fgumi_bam_io::DecodedRecord::raw_bytes))
        .collect();
    let expected: Vec<&[u8]> = records.iter().map(AsRef::as_ref).collect();
    assert_eq!(seen, expected, "no-header mode must treat byte 0 as a record boundary");
}

/// `… → ParseBamRecords → DecodeFromRecords` is the runall fusion shape: parse
/// once into shared-buffer batches, then re-attach group keys downstream
/// instead of re-serializing. It must yield the same records, in the same
/// order, as the fused `DecodeRecords` chain above.
#[rstest]
fn decode_from_records_chain_matches_the_fused_decode_chain(#[values(1, 4)] threads: usize) {
    let records = test_records(120);

    let collect_via = |use_split: bool| -> Vec<Vec<u8>> {
        let blocks = bgzf_blocks_for(&records, 4096);
        let collected: Arc<Mutex<Vec<DecodedRecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
        let sink_handle = Arc::clone(&collected);
        let builder = Pipeline::builder();
        let head = builder
            .chain(ReplaySource::new(blocks))
            .chain(BgzfDecompress::new(EDGE_LIMIT_BYTES))
            .chain(FindBamBoundaries::new(EDGE_LIMIT_BYTES));
        if use_split {
            head.chain(ParseBamRecords::new(EDGE_LIMIT_BYTES))
                .chain(DecodeFromRecords::new(GroupKeyConfig::default(), EDGE_LIMIT_BYTES))
                .chain(CollectSink { collected: sink_handle })
                .into_sink_marker();
        } else {
            head.chain(DecodeRecords::new(GroupKeyConfig::default(), EDGE_LIMIT_BYTES))
                .chain(CollectSink { collected: sink_handle })
                .into_sink_marker();
        }
        let pipeline = builder.build().expect("chain builds");
        pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

        let batches = collected.lock().expect("mutex not poisoned");
        batches
            .iter()
            .flat_map(|batch| batch.records().iter().map(|r| r.raw_bytes().to_vec()))
            .collect()
    };

    let expected: Vec<Vec<u8>> = records.iter().map(|r| r.as_ref().to_vec()).collect();
    assert_eq!(collect_via(false), expected, "fused decode chain must round-trip records");
    assert_eq!(collect_via(true), expected, "split parse+decode chain must agree with it");
}

/// `ReplaySource → ParseSamChunk` covers the SAM ingest path, whose records
/// converge on the same `DecodedRecordBatch` the BAM chain produces. Driving it
/// through the framework (rather than calling `parse_sam_chunk_into_decoded`
/// directly) is what exercises the step's drain and held-item handling.
#[rstest]
fn sam_chain_decodes_every_line_in_order(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::parse::sam::ParseSamChunk;
    use crate::pipeline::steps::types::SamChunk;

    const HEADER_TEXT: &str = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n";
    let header = {
        let mut reader = noodles::sam::io::Reader::new(HEADER_TEXT.as_bytes());
        Arc::new(reader.read_header().expect("header parses"))
    };

    // One chunk per record keeps the ordinal stream long enough that the
    // byte-bounded edge rejects pushes and the retry path is taken.
    let n_records = 200usize;
    let chunks: Vec<SamChunk> = (0..n_records)
        .map(|i| {
            let line = format!("read{i}\t0\tchr1\t{}\t60\t4M\t*\t0\t0\tACGT\tIIII\n", 100 + i);
            let bytes = line.into_bytes();
            let line_offsets = vec![0u32, u32::try_from(bytes.len()).expect("fits u32")];
            SamChunk { batch_serial: i as u64, bytes, line_offsets }
        })
        .collect();

    let key_config =
        GroupKeyConfig::new_raw_no_cell(fgumi_bam_io::LibraryIndex::from_header(&header));
    let collected: Arc<Mutex<Vec<DecodedRecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(chunks))
        .chain(ParseSamChunk::new(Arc::clone(&header), key_config, EDGE_LIMIT_BYTES))
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let batches = collected.lock().expect("mutex not poisoned");
    let names: Vec<Vec<u8>> = batches
        .iter()
        .flat_map(DecodedRecordBatch::records)
        .map(|record| {
            let body = record.raw_bytes();
            let l_read_name = body[8] as usize;
            body[32..32 + l_read_name - 1].to_vec()
        })
        .collect();
    let expected: Vec<Vec<u8>> = (0..n_records).map(|i| format!("read{i}").into_bytes()).collect();
    assert_eq!(names, expected, "SAM lines must decode in order, one record per line");
}

/// A stream that ends mid-record must fail the run, not silently drop the
/// truncated tail. `FindBamBoundaries` carries the partial record forward, and
/// at drain its EOF validation rejects it — so the error surfaces through
/// `Pipeline::run` rather than showing up as quietly missing records.
#[rstest]
fn a_truncated_record_stream_fails_the_run(#[values(1, 4)] threads: usize) {
    let records = test_records(32);
    let mut stream = minimal_binary_bam_header(1);
    for record in &records {
        let bytes = record.as_ref();
        stream.extend_from_slice(&u32::try_from(bytes.len()).expect("fits u32").to_le_bytes());
        stream.extend_from_slice(bytes);
    }
    // Chop the final record in half — a `block_size` prefix promising bytes
    // that never arrive.
    stream.truncate(stream.len() - 20);

    let blocks: Vec<BgzfBlock> = stream
        .chunks(1024)
        .enumerate()
        .map(|(i, payload)| {
            let mut compressor = fgumi_bgzf::InlineBgzfCompressor::new(1);
            compressor.write_all(payload).expect("compress");
            compressor.flush().expect("flush");
            let mut blocks = compressor.take_blocks();
            assert_eq!(blocks.len(), 1);
            BgzfBlock {
                batch_serial: i as u64,
                bytes: blocks.remove(0).data,
                uncompressed_size: u32::try_from(payload.len()).expect("fits u32"),
            }
        })
        .collect();

    let collected: Arc<Mutex<Vec<DecodedRecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(blocks))
        .chain(BgzfDecompress::new(EDGE_LIMIT_BYTES))
        .chain(FindBamBoundaries::new(EDGE_LIMIT_BYTES))
        .chain(DecodeRecords::new(GroupKeyConfig::default(), EDGE_LIMIT_BYTES))
        .chain(CollectSink { collected: Arc::clone(&collected) })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");

    let err = pipeline
        .run(PipelineConfig { threads, ..Default::default() })
        .expect_err("a truncated record stream must fail the run");
    match err {
        crate::pipeline::core::signal::PipelineError::Io { step, source } => {
            assert_eq!(step, "FindBamBoundaries", "the boundary scanner owns EOF validation");
            assert_eq!(source.kind(), io::ErrorKind::UnexpectedEof);
        }
        other => panic!("expected an I/O failure from the boundary step, got {other:?}"),
    }
}

/// `ReplaySource → TemplatesToRecordBatch` is the AAM-fusion adapter: it turns
/// the queryname-template view back into the flat-record view the sort ingest
/// consumes. Driving it through `Pipeline::run` (rather than re-running its
/// flattening loop in a unit test) is what exercises `try_run` itself — the
/// held-slot retry, the drain path, and the last-clone `Finished` gate that a
/// `Parallel` step only reaches under the framework.
///
/// Every record of every template must arrive downstream exactly once,
/// byte-for-byte, in template order, in a batch carrying the input batch's
/// serial — *including* the secondary and supplementary alignments.
#[rstest]
fn templates_to_records_flattens_every_record_including_split_alignments(
    #[values(1, 4)] threads: usize,
) {
    use crate::pipeline::steps::templates_to_records::TemplatesToRecordBatch;
    use crate::pipeline::steps::types::BamTemplateBatch;

    const N_BATCHES: u64 = 200;
    const TEMPLATES_PER_BATCH: usize = 2;
    const RECORDS_PER_TEMPLATE: usize = 4;

    // Record each template's records in `records` order — the order the step is
    // contracted to emit. Reading that off the template rather than hard-coding
    // r1/r2/supplementary/secondary keeps this a test of the flattening step,
    // not of `Template::from_records`' internal layout.
    let mut expected: Vec<Vec<u8>> = Vec::new();
    let mut batches: Vec<BamTemplateBatch> = Vec::new();
    for serial in 0..N_BATCHES {
        let templates: Vec<Template> = (0..TEMPLATES_PER_BATCH)
            .map(|t| split_alignment_template(format!("q{serial}_{t}").as_bytes()))
            .collect();
        for template in &templates {
            expected.extend(template.records.iter().map(|r| r.as_ref().to_vec()));
        }
        batches.push(BamTemplateBatch::new(serial, templates));
    }

    // Same guard as `the_edge_budget_binds_on_the_multi_record_fixtures`: unless
    // the fixture volume exceeds the per-edge budget the output queue never
    // rejects a push, and the step's held-slot retry path — the whole reason for
    // running this through the framework — goes unexercised.
    let fixture_bytes: usize = batches.iter().map(BamTemplateBatch::total_bytes).sum();
    assert!(
        u64::try_from(fixture_bytes).expect("fixture volume fits u64") > EDGE_LIMIT_BYTES,
        "fixture volume ({fixture_bytes} B) must exceed the per-edge budget \
         ({EDGE_LIMIT_BYTES} B), else the held-slot retry path never runs",
    );

    let collected: Arc<Mutex<Vec<RecordBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(batches))
        .chain(TemplatesToRecordBatch::new(EDGE_LIMIT_BYTES))
        .chain(StallingCollectSink {
            collected: sink_handle,
            remaining_stalls: PROCESS_SINK_STALLS,
        })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let emitted = collected.lock().expect("mutex not poisoned");

    // One output batch per input batch, each carrying its input's serial, in
    // ordinal order — the contract `BranchOrdering::ByItemOrdinal` consumers
    // downstream depend on.
    let serials: Vec<u64> = emitted.iter().map(RecordBatch::batch_serial).collect();
    assert_eq!(
        serials,
        (0..N_BATCHES).collect::<Vec<u64>>(),
        "each input batch must yield one output batch carrying its serial, in order",
    );

    // Batch boundaries are preserved 1:1 — the step neither splits an input
    // batch across outputs nor coalesces several into one. Reported as the first
    // offending batch rather than an `assert_eq!` over the whole vector, which
    // would dump 200 elements to name a single bad one.
    let expected_per_batch = TEMPLATES_PER_BATCH * RECORDS_PER_TEMPLATE;
    if let Some((i, n)) =
        emitted.iter().map(RecordBatch::len).enumerate().find(|&(_, n)| n != expected_per_batch)
    {
        panic!(
            "output batch {i} holds {n} records, expected {expected_per_batch} — \
             every output batch must hold exactly its input batch's records",
        );
    }

    let seen: Vec<Vec<u8>> =
        emitted.iter().flat_map(|batch| batch.iter_record_bytes().map(<[u8]>::to_vec)).collect();
    assert_eq!(
        seen.len(),
        expected.len(),
        "every record — primary, supplementary, and secondary — must be emitted exactly once",
    );
    // Likewise pinpointed: a whole-vector compare of ~1,600 records would print
    // the entire fixture to report one differing byte.
    if let Some((i, (got, want))) =
        seen.iter().zip(&expected).enumerate().find(|(_, (got, want))| got != want)
    {
        panic!("record {i} was altered in flight: got {got:02x?}, expected {want:02x?}");
    }
}

/// `ReplaySource → GroupByMi` drives the MI grouper through the framework.
///
/// Its unit tests all call `process_record` / `flush_current_group` directly,
/// which leaves `try_run` — where the ordinals are minted, the held slot is
/// retried, and drain is sequenced — entirely unexercised. Each of those fails
/// by *silently truncating* output rather than erroring, so nothing downstream
/// would notice:
///
/// * **Ordinal contiguity.** `GroupByMi` mints its own batch serials, and the
///   `ByItemOrdinal` reorder stage releases only in contiguous order — a
///   skipped or reused serial stalls it forever.
/// * **Held-slot retry.** `output_byte_limit` is deliberately smaller than one
///   emitted batch, and the sink stalls long enough for the queue to fill, so
///   pushes bounce and `emit_batch` parks the batch in `held`. The retry must
///   not mint a second ordinal for it. Both halves are needed: with a sink that
///   drains every tick the queue empties between pushes and this path is never
///   taken at all (verified by injecting a double-ordinal bug, which a
///   non-stalling sink failed to catch at one thread).
/// * **Drain sequencing.** The molecule count is not a multiple of the target
///   batch count, so the run-in-progress must be flushed and emitted as a final
///   partial batch *before* `Finished`.
#[rstest]
fn group_by_mi_emits_contiguous_batches_ending_in_a_partial(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::group::mi::{BatchedMiGroups, GroupByMi};

    const N_MOLECULES: usize = 337;
    const RECORDS_PER_MOLECULE: usize = 3;
    const TARGET_BATCH_COUNT: usize = 10;
    // Deliberately coprime with `RECORDS_PER_MOLECULE` so molecule runs straddle
    // input-batch boundaries — the run-length state must survive across calls.
    const RECORDS_PER_INPUT_BATCH: usize = 7;
    // Smaller than one emitted batch, so the second push always bounces.
    const OUTPUT_BYTE_LIMIT: u64 = 1024;
    // Long enough for the grouper to emit and then bounce against a full queue,
    // but no longer: each stall costs the driver a back-off interval, so a
    // larger count buys nothing and adds seconds to the multi-threaded case.
    // Calibrated by injecting the double-ordinal bug and finding the smallest
    // count the single-threaded case still catches it at (16), then leaving a
    // margin for scheduling variance. At 8 the bug slips through.
    const SINK_STALLS: usize = 24;

    assert_ne!(
        N_MOLECULES % TARGET_BATCH_COUNT,
        0,
        "the fixture must leave a partial final batch, else the drain arm is untested",
    );

    let records: Vec<RawRecord> = (0..N_MOLECULES)
        .flat_map(|m| (0..RECORDS_PER_MOLECULE).map(move |_| mi_record(&format!("mol{m}"))))
        .collect();
    let batches: Vec<DecodedRecordBatch> = records
        .chunks(RECORDS_PER_INPUT_BATCH)
        .enumerate()
        .map(|(i, chunk)| decoded_batch(i as u64, chunk.to_vec()))
        .collect();

    let collected: Arc<Mutex<Vec<BatchedMiGroups>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(batches))
        .chain(GroupByMi::with_target_batch_count(
            *crate::sam::SamTag::MI,
            OUTPUT_BYTE_LIMIT,
            TARGET_BATCH_COUNT,
        ))
        .chain(StallingCollectSink { collected: sink_handle, remaining_stalls: SINK_STALLS })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let emitted = collected.lock().expect("mutex not poisoned");

    // Several batches, or the contiguity claim below is vacuous.
    assert!(
        emitted.len() > 1,
        "fixture must span several output batches, got {}; the ordinal contiguity and \
         held-slot paths only run when more than one batch is emitted",
        emitted.len(),
    );

    // The load-bearing assertion: serials are dense and start at zero. A held
    // batch that minted a second ordinal on retry, or a dropped emit, shows up
    // here as a gap or a duplicate.
    let serials: Vec<u64> = emitted.iter().map(|b| b.batch_serial).collect();
    assert_eq!(
        serials,
        (0..emitted.len() as u64).collect::<Vec<u64>>(),
        "emitted batch serials must be contiguous from zero",
    );

    // No batch is emitted empty — an empty emit would burn an ordinal for
    // nothing, which the drain arm's `!accumulator.is_empty()` guard prevents.
    assert!(emitted.iter().all(|b| !b.groups.is_empty()), "no batch may be emitted empty");

    // Every molecule arrives exactly once, in order, with all its records.
    let groups: Vec<(String, usize)> =
        emitted.iter().flat_map(|b| &b.groups).map(|g| (g.mi.clone(), g.records.len())).collect();
    let expected: Vec<(String, usize)> =
        (0..N_MOLECULES).map(|m| (format!("mol{m}"), RECORDS_PER_MOLECULE)).collect();
    assert_eq!(
        groups.len(),
        expected.len(),
        "every molecule must be emitted exactly once — a truncated tail means the drain \
         arm dropped the final partial batch",
    );
    if let Some((i, (got, want))) =
        groups.iter().zip(&expected).enumerate().find(|(_, (got, want))| got != want)
    {
        panic!("group {i} differs: got {got:?}, expected {want:?}");
    }

    // The batch cap: no emitted batch may exceed the target, even though one
    // `add_records` call can complete more groups than that. Without this, a
    // regression that let the accumulator grow past the cap before emitting
    // would keep serials dense, keep every molecule present, and still leave a
    // short final batch — passing every other assertion here.
    if let Some((i, n)) =
        emitted.iter().map(|b| b.groups.len()).enumerate().find(|&(_, n)| n > TARGET_BATCH_COUNT)
    {
        panic!("batch {i} holds {n} groups, above the {TARGET_BATCH_COUNT} cap");
    }

    // The last batch is the partial one the drain arm flushed.
    let last = emitted.last().expect("at least one batch");
    assert!(
        last.groups.len() < TARGET_BATCH_COUNT,
        "the final batch must be the short one emitted at drain, got {} groups",
        last.groups.len(),
    );
}

/// Records with no usable MI tag must be dropped without disturbing the
/// grouping of the records around them, and the chain must still reach the
/// drain arm — which is where the skipped-record warnings are reported.
///
/// The counters themselves are unit-tested; what only the framework reaches is
/// `try_run`'s drain arm calling `report_skipped_records` and *still* returning
/// `Finished`. A regression that returned early there would hang the chain
/// rather than truncate it.
#[rstest]
fn group_by_mi_drops_untagged_records_and_still_completes(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::group::mi::{BatchedMiGroups, GroupByMi};

    // Untagged and integer-MI records are interleaved *inside* a molecule run,
    // so a regression that reset the run on a skip would split "mol0" in two.
    let records = [
        record_without_mi(),
        mi_record("mol0"),
        record_with_integer_mi(1),
        mi_record("mol0"),
        record_without_mi(),
        mi_record("mol1"),
        record_with_integer_mi(2),
    ];
    let batches: Vec<DecodedRecordBatch> = records
        .chunks(2)
        .enumerate()
        .map(|(i, chunk)| decoded_batch(i as u64, chunk.to_vec()))
        .collect();

    let collected: Arc<Mutex<Vec<BatchedMiGroups>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(batches))
        .chain(GroupByMi::new(*crate::sam::SamTag::MI, EDGE_LIMIT_BYTES))
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let emitted = collected.lock().expect("mutex not poisoned");
    let groups: Vec<(String, usize)> =
        emitted.iter().flat_map(|b| &b.groups).map(|g| (g.mi.clone(), g.records.len())).collect();
    assert_eq!(
        groups,
        vec![("mol0".to_string(), 2), ("mol1".to_string(), 1)],
        "untagged and wrong-typed records must be dropped without splitting the runs \
         they sit inside",
    );
}

/// A trivial ordered item: `ordinal` is its place in a branch's sequence,
/// `value` identifies which input produced it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Numbered {
    ordinal: u64,
    value: u64,
}

/// A deliberately non-zero weight. `HeapSize`'s default returns 0, and a
/// zero-weight item can never fill a `ByteBounded` queue — `current_bytes`
/// stays at 0, every push is admitted, and the held-slot retry arms of every
/// byte-bounded step below become unreachable. Giving the item a weight is what
/// lets a small `limit_bytes` actually bind.
impl HeapSize for Numbered {
    fn heap_size(&self) -> usize {
        128
    }
}

impl Ordered for Numbered {
    fn ordinal(&self) -> u64 {
        self.ordinal
    }
}

/// Returning `None` for a branch of a [`Process2Ordered`] must not stall the
/// chain **when that branch mints its own dense, branch-local ordinals**.
///
/// This is the safe half of the contract documented on `Process2Output`. The
/// unsafe half — a branch that propagates the *input's* ordinal and then omits
/// an item — leaves a permanent hole, and the `ByItemOrdinal` reorder stage
/// releases only in contiguous order, so it waits for an ordinal that will
/// never arrive. That failure mode is a **hang**, which is precisely why it is
/// not reproduced here: a test that builds the wedged configuration on purpose
/// can only assert "still not finished after N seconds", which is slow, and
/// proves a timeout rather than the contract. Asserting the safe direction
/// positively pins the same rule without ever constructing the stall.
///
/// Note the test is still hang-*sensitive*, just not hang-*dependent*: if a
/// regression made sparse output wedge this chain, the run would never return
/// and nextest would kill it at its slow-test timeout. It fails either way —
/// it just does not need the stall to pass.
#[rstest]
fn sparse_ordered_fan_out_completes_when_each_branch_mints_its_own_ordinals(
    #[values(1, 4)] threads: usize,
) {
    use crate::pipeline::steps::process::{Process2Output, process2_ordered};

    const N_INPUTS: u64 = 300;
    // Every third input goes to B, the rest to A — so *both* branches see
    // `None` for a substantial share of the inputs.
    const B_EVERY: u64 = 3;

    let inputs: Vec<Numbered> = (0..N_INPUTS).map(|i| Numbered { ordinal: i, value: i }).collect();

    // Branch-local ordinal counters. Shared across worker clones (the closure is
    // `Fn` behind an `Arc`), so at four threads the *assignment* of ordinals to
    // values races — but the sequence each branch mints stays dense, which is
    // the property the reorder stage actually requires.
    let next_a = Arc::new(AtomicU64::new(0));
    let next_b = Arc::new(AtomicU64::new(0));
    let (mint_a, mint_b) = (Arc::clone(&next_a), Arc::clone(&next_b));

    let collected_a: Arc<Mutex<Vec<Numbered>>> = Arc::new(Mutex::new(Vec::new()));
    let collected_b: Arc<Mutex<Vec<Numbered>>> = Arc::new(Mutex::new(Vec::new()));
    let (sink_a, sink_b) = (Arc::clone(&collected_a), Arc::clone(&collected_b));

    let builder = Pipeline::builder();
    let multi = builder
        .chain(ReplaySource::new(inputs))
        .chain(process2_ordered(
            "SparseFanOut",
            EDGE_LIMIT_BYTES,
            EDGE_LIMIT_BYTES,
            move |item: Numbered| {
                if item.value.is_multiple_of(B_EVERY) {
                    let ordinal = mint_b.fetch_add(1, AtomicOrdering::Relaxed);
                    Ok(Process2Output::only_b(Numbered { ordinal, value: item.value }))
                } else {
                    let ordinal = mint_a.fetch_add(1, AtomicOrdering::Relaxed);
                    Ok(Process2Output::only_a(Numbered { ordinal, value: item.value }))
                }
            },
        ))
        .into_multi();
    multi.b0.chain(CollectSink { collected: sink_a }).into_sink_marker();
    multi.b1.chain(CollectSink { collected: sink_b }).into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let branch_a = collected_a.lock().expect("mutex not poisoned");
    let branch_b = collected_b.lock().expect("mutex not poisoned");

    // Both branches must have skipped some inputs, or "sparse" is vacuous and
    // this degenerates into an ordinary every-item fan-out test.
    assert!(!branch_a.is_empty() && !branch_b.is_empty(), "both branches must receive items");
    assert_eq!(
        branch_a.len() as u64 + branch_b.len() as u64,
        N_INPUTS,
        "every input must reach exactly one branch",
    );

    // The load-bearing assertion: each branch released a dense sequence from
    // zero, in order. A hole would have stalled the reorder stage instead.
    for (name, branch) in [("a", &*branch_a), ("b", &*branch_b)] {
        let ordinals: Vec<u64> = branch.iter().map(|n| n.ordinal).collect();
        assert_eq!(
            ordinals,
            (0..branch.len() as u64).collect::<Vec<u64>>(),
            "branch {name} must release a dense ordinal sequence from zero",
        );
    }

    // Nothing was routed to the wrong branch, dropped, or duplicated. Sorted
    // because ordinal assignment races across workers at four threads.
    let mut values_a: Vec<u64> = branch_a.iter().map(|n| n.value).collect();
    let mut values_b: Vec<u64> = branch_b.iter().map(|n| n.value).collect();
    values_a.sort_unstable();
    values_b.sort_unstable();
    assert_eq!(
        values_a,
        (0..N_INPUTS).filter(|v| !v.is_multiple_of(B_EVERY)).collect::<Vec<u64>>()
    );
    assert_eq!(values_b, (0..N_INPUTS).filter(|v| v.is_multiple_of(B_EVERY)).collect::<Vec<u64>>());
}

/// A single-end primary record with the given queryname.
fn single_end_record(qname: &[u8]) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(qname)
        .flags(0)
        .ref_id(0)
        .pos(100)
        .cigar_ops(&[4u32 << 4])
        .sequence(b"ACGT")
        .qualities(&[30u8; 4]);
    b.build()
}

/// Wrap raw records in a [`DecodedRecordBatch`], pairing each with an explicit
/// [`GroupKey`]. Position and queryname grouping both key off the pre-computed
/// key rather than re-parsing the record, so the key is the fixture.
fn decoded_batch_with_keys(
    batch_serial: u64,
    records: Vec<(RawRecord, fgumi_bam_io::GroupKey)>,
) -> DecodedRecordBatch {
    DecodedRecordBatch::new(
        batch_serial,
        records
            .into_iter()
            .map(|(raw, key)| fgumi_bam_io::DecodedRecord::from_raw_bytes(raw, key))
            .collect(),
    )
}

/// A single-end position key at `pos` — `strand2` stays `UNKNOWN_STRAND`, so
/// `has_mate_position()` is false and the grouper's MC-tag validation (which
/// only applies to *paired* primaries) is not triggered.
fn position_key_at(pos: i32) -> fgumi_bam_io::GroupKey {
    fgumi_bam_io::GroupKey { ref_id1: 0, pos1: pos, strand1: 0, ..Default::default() }
}

/// Hash a queryname the way the decode step would, so `GroupByQueryname`'s
/// `name_hash` fast-path pre-check is exercised rather than short-circuited by
/// a constant hash.
///
/// This routes through `LibraryIndex::hash_name`, the exact hasher
/// `compute_group_key_from_raw` uses to fill `GroupKey::name_hash` on the decode
/// path (an empty name maps to `hash_name(None)`), so a key built here is
/// byte-for-byte comparable to one the grouper computes rather than merely
/// self-consistent.
fn name_hash_of(qname: &[u8]) -> u64 {
    let name = if qname.is_empty() { None } else { Some(qname) };
    fgumi_bam_io::LibraryIndex::hash_name(name)
}

/// Split a `DecompressedBlock` byte stream — framed `block_size` + body per
/// record — back into record bodies, asserting the framing is well-formed.
///
/// The framing is half of what `SerializeGroups` produces, so parsing it back
/// rather than comparing opaque bytes is what makes a truncated or
/// wrongly-sized prefix fail loudly here instead of much further downstream.
fn unframe_records(bytes: &[u8]) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    let mut pos = 0usize;
    while pos < bytes.len() {
        assert!(pos + 4 <= bytes.len(), "truncated block_size prefix at byte {pos}");
        let len = u32::from_le_bytes(bytes[pos..pos + 4].try_into().expect("4 bytes")) as usize;
        pos += 4;
        assert!(pos + len <= bytes.len(), "block_size {len} at byte {pos} runs past the stream");
        out.push(bytes[pos..pos + len].to_vec());
        pos += len;
    }
    out
}

/// One `ProcessedPositionGroup` holding `n_templates` paired templates, each
/// tagged with the molecule id `mi` (use [`MoleculeId::None`] for the
/// pass-through path, which must emit the record bytes unchanged).
fn processed_group(
    qname_prefix: &str,
    n_templates: usize,
    mi: fgumi_umi::MoleculeId,
) -> ProcessedPositionGroup {
    let templates: Vec<Template> = (0..n_templates)
        .map(|t| {
            let mut template = split_alignment_template(format!("{qname_prefix}_{t}").as_bytes());
            template.mi = mi;
            template
        })
        .collect();
    let input_record_count = templates.iter().map(|t| t.records.len() as u64).sum();
    ProcessedPositionGroup {
        templates,
        family_sizes: ahash::AHashMap::new(),
        filter_metrics: crate::grouper::FilterMetrics::default(),
        input_record_count,
        distinct_mi_count: u64::from(mi.id().is_some()),
    }
}

/// `ReplaySource → SerializeGroups` covers `build_serialize_processed_groups_step`,
/// which had no test at all.
///
/// The step does two things in one pass — rewrite the MI tag into each record's
/// body, and frame each body with its BAM `block_size` prefix — and both are
/// asserted here by parsing the emitted stream back apart.
///
/// The MI half matters beyond line coverage: `emit_templates_raw_with_mi` is
/// deliberately shared with the single-threaded writer in `commands::group`
/// precisely so the two cannot drift in MI-injection semantics. Nothing was
/// pinning that shared behavior, so a change to it would have been invisible
/// on both paths.
///
/// Both branches of the emit are exercised: templates with an assigned
/// `MoleculeId` get `MI:Z:<id>` written into their records, while templates
/// left at `MoleculeId::None` must pass through byte-for-byte untouched.
///
/// Note what the fixture pins about *which* records are emitted. The templates
/// here carry four records each — primary R1 and R2 plus a supplementary and a
/// secondary — and the step emits **only the two primaries**, because
/// `emit_templates_raw_with_mi` iterates `[template.r1(), template.r2()]`. That
/// matches its doc ("primary reads") and is asserted rather than worked around:
/// it is a real output-shaping decision, and a change that started emitting
/// split alignments would otherwise slip through silently.
#[rstest]
fn serialize_processed_groups_frames_every_record_and_injects_mi(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::serialize_processed::build_serialize_processed_groups_step;
    use crate::pipeline::steps::types::DecompressedBlock;
    use fgumi_umi::MoleculeId;

    const N_BATCHES: u64 = 40;
    const TEMPLATES_PER_GROUP: usize = 2;
    // Records per template that the step *emits* (the two primaries). The
    // fixture holds four; the progress counter instead totals every *input*
    // record, which is why it has its own constant.
    const EMITTED_PER_TEMPLATE: usize = 2;
    const INPUT_RECORDS_PER_TEMPLATE: usize = 4;
    // Two groups per batch: one MI-assigned, one left unassigned, so both emit
    // branches run for every batch.
    const GROUPS_PER_BATCH: usize = 2;

    let mut expected: Vec<(Vec<u8>, Option<String>)> = Vec::new();
    let mut batches: Vec<BatchedProcessedPositionGroups> = Vec::new();
    for serial in 0..N_BATCHES {
        let assigned = processed_group(
            &format!("mi{serial}"),
            TEMPLATES_PER_GROUP,
            MoleculeId::Single(serial),
        );
        let untouched =
            processed_group(&format!("raw{serial}"), TEMPLATES_PER_GROUP, MoleculeId::None);
        for group in [&assigned, &untouched] {
            for template in &group.templates {
                let label = template.mi.id().map(|id| id.to_string());
                // ONLY the primary R1/R2 — see the note on the test.
                for raw in [template.r1(), template.r2()].into_iter().flatten() {
                    expected.push((raw.to_vec(), label.clone()));
                }
            }
        }
        batches.push(BatchedProcessedPositionGroups::new(serial, vec![assigned, untouched]));
    }

    let collected: Arc<Mutex<Vec<DecompressedBlock>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let progress = Arc::new(AtomicU64::new(0));
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(batches))
        .chain(build_serialize_processed_groups_step(
            EDGE_LIMIT_BYTES,
            *crate::sam::SamTag::MI,
            Arc::clone(&progress),
        ))
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let emitted = collected.lock().expect("mutex not poisoned");

    let serials: Vec<u64> = emitted.iter().map(|b| b.batch_serial).collect();
    assert_eq!(
        serials,
        (0..N_BATCHES).collect::<Vec<u64>>(),
        "each input batch must yield one block carrying its serial, in order",
    );

    // The progress counter is bumped with each batch's *input record* count.
    let expected_records =
        N_BATCHES * (GROUPS_PER_BATCH * TEMPLATES_PER_GROUP * INPUT_RECORDS_PER_TEMPLATE) as u64;
    assert_eq!(
        progress.load(AtomicOrdering::Relaxed),
        expected_records,
        "the shared progress counter must total every input record",
    );

    let bodies: Vec<Vec<u8>> =
        emitted.iter().flat_map(|block| unframe_records(&block.bytes)).collect();
    assert_eq!(
        bodies.len(),
        expected.len(),
        "every primary record must be framed and emitted exactly once",
    );
    // Stated as its own arithmetic so the primaries-only contract is visible
    // rather than implied by the length of `expected`.
    assert_eq!(
        bodies.len() as u64,
        N_BATCHES * (GROUPS_PER_BATCH * TEMPLATES_PER_GROUP * EMITTED_PER_TEMPLATE) as u64,
        "the step emits primaries only — supplementary and secondary records are excluded",
    );
    assert!(
        bodies.iter().all(|b| {
            let flags = u16::from_le_bytes([b[14], b[15]]);
            flags & (fgumi_raw_bam::flags::SECONDARY | fgumi_raw_bam::flags::SUPPLEMENTARY) == 0
        }),
        "no emitted record may carry the secondary or supplementary flag",
    );

    for (i, (body, (original, want_mi))) in bodies.iter().zip(&expected).enumerate() {
        let got_mi = fgumi_raw_bam::find_string_tag_in_record(body, *crate::sam::SamTag::MI)
            .map(|v| String::from_utf8_lossy(v).into_owned());
        assert_eq!(got_mi.as_deref(), want_mi.as_deref(), "record {i}: MI tag");
        if want_mi.is_none() {
            assert_eq!(
                body, original,
                "record {i}: an unassigned template must pass through byte-for-byte",
            );
        }
    }
}

/// `ReplaySource → GroupByPosition` drives the position grouper through the
/// framework. Like the MI grouper, its unit tests exercise helpers directly and
/// leave `try_run`/`emit_batch` — the ordinal minting, the batch cap, and the
/// `finish()`-once drain path — unrun.
///
/// Positions arrive in runs, so a group closes only when the position changes;
/// the last run has no successor and can only be emitted by the drain arm
/// calling the (non-idempotent) `grouper.finish()`, guarded by `finalized`.
#[rstest]
fn group_by_position_emits_contiguous_batches_ending_in_a_partial(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::group::position::{BatchedRawPositionGroups, GroupByPosition};

    const N_POSITIONS: usize = 205;
    const RECORDS_PER_POSITION: usize = 3;
    const TARGET_BATCH_COUNT: usize = 8;
    const RECORDS_PER_INPUT_BATCH: usize = 7;
    const OUTPUT_BYTE_LIMIT: u64 = 1024;
    const SINK_STALLS: usize = 24;

    assert_ne!(
        N_POSITIONS % TARGET_BATCH_COUNT,
        0,
        "the fixture must leave a partial final batch, else the drain arm is untested",
    );

    let records: Vec<(RawRecord, fgumi_bam_io::GroupKey)> = (0..N_POSITIONS)
        .flat_map(|p| {
            let key = position_key_at(i32::try_from(100 + p * 10).expect("pos fits i32"));
            (0..RECORDS_PER_POSITION)
                .map(move |r| (single_end_record(format!("p{p}_{r}").as_bytes()), key))
        })
        .collect();
    let batches: Vec<DecodedRecordBatch> = records
        .chunks(RECORDS_PER_INPUT_BATCH)
        .enumerate()
        .map(|(i, chunk)| decoded_batch_with_keys(i as u64, chunk.to_vec()))
        .collect();

    let collected: Arc<Mutex<Vec<BatchedRawPositionGroups>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(batches))
        .chain(GroupByPosition::with_target_batch_count(OUTPUT_BYTE_LIMIT, TARGET_BATCH_COUNT))
        .chain(StallingCollectSink { collected: sink_handle, remaining_stalls: SINK_STALLS })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let emitted = collected.lock().expect("mutex not poisoned");
    assert!(emitted.len() > 1, "fixture must span several output batches");

    let serials: Vec<u64> = emitted.iter().map(|b| b.batch_serial).collect();
    assert_eq!(
        serials,
        (0..emitted.len() as u64).collect::<Vec<u64>>(),
        "emitted batch serials must be contiguous from zero",
    );
    assert!(emitted.iter().all(|b| !b.groups.is_empty()), "no batch may be emitted empty");

    // The batch cap: no emitted batch may exceed the target, even though one
    // `add_records` call can complete more groups than that.
    if let Some((i, n)) =
        emitted.iter().map(|b| b.groups.len()).enumerate().find(|&(_, n)| n > TARGET_BATCH_COUNT)
    {
        panic!("batch {i} holds {n} groups, above the {TARGET_BATCH_COUNT} cap");
    }

    let sizes: Vec<usize> =
        emitted.iter().flat_map(|b| &b.groups).map(|g| g.records.len()).collect();
    assert_eq!(
        sizes.len(),
        N_POSITIONS,
        "every position must be emitted exactly once — a short count means the drain arm \
         dropped the final run that only `finish()` can close",
    );
    assert!(
        sizes.iter().all(|&n| n == RECORDS_PER_POSITION),
        "every position group must hold all of its records",
    );
}

/// The default grouper drops secondary/supplementary records — they decode to
/// an `UNKNOWN_REF` position key and have no position of their own — while
/// `GroupByPosition::with_secondary_supplementary` keeps them, which is what
/// `fgumi dedup` needs so the duplicate flag propagates across split
/// alignments.
///
/// Run as one parameterized pair so the two constructors are compared on
/// identical input; asserting only the inclusive case would not show that the
/// default actually excludes.
#[rstest]
#[case::default_grouper_drops_them(false, 1)]
#[case::with_secondary_supplementary_keeps_them(true, 2)]
fn group_by_position_secondary_supplementary_policy(
    #[case] include_secondary_supplementary: bool,
    #[case] expected_groups: usize,
) {
    use crate::pipeline::steps::group::position::{BatchedRawPositionGroups, GroupByPosition};
    use fgumi_bam_io::GroupKey;

    // One positioned run, then a record carrying the UNKNOWN_REF key that a
    // secondary/supplementary alignment decodes to.
    let positioned = position_key_at(100);
    let unknown = GroupKey::default(); // ref_id1 == UNKNOWN_REF
    let records = vec![
        (single_end_record(b"a"), positioned),
        (single_end_record(b"b"), positioned),
        (single_end_record(b"supp"), unknown),
    ];

    let collected: Arc<Mutex<Vec<BatchedRawPositionGroups>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let step = if include_secondary_supplementary {
        GroupByPosition::with_secondary_supplementary(EDGE_LIMIT_BYTES)
    } else {
        GroupByPosition::new(EDGE_LIMIT_BYTES)
    };
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(vec![decoded_batch_with_keys(0, records)]))
        .chain(step)
        .chain(CollectSink { collected: sink_handle })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads: 1, ..Default::default() }).expect("chain runs");

    let emitted = collected.lock().expect("mutex not poisoned");
    let groups: Vec<usize> =
        emitted.iter().flat_map(|b| &b.groups).map(|g| g.records.len()).collect();
    assert_eq!(groups.len(), expected_groups, "position groups emitted: {groups:?}");
    assert_eq!(groups[0], 2, "the positioned run always groups its two records");
}

/// `ReplaySource → GroupByQueryname` drives the queryname grouper through the
/// framework, covering the same `try_run`/`emit_batch` surface: ordinal
/// minting, the held-slot retry, and the drain arm that closes the final
/// run-in-progress before reporting `Finished`.
#[rstest]
fn group_by_queryname_emits_contiguous_batches_ending_in_a_partial(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::group::queryname::GroupByQueryname;
    use crate::pipeline::steps::types::BamTemplateBatch;
    use fgumi_raw_bam::flags::{FIRST_SEGMENT, LAST_SEGMENT, PAIRED};

    const N_TEMPLATES: usize = 203;
    const TARGET_BATCH_COUNT: usize = 8;
    const RECORDS_PER_INPUT_BATCH: usize = 5;
    const OUTPUT_BYTE_LIMIT: u64 = 1024;
    const SINK_STALLS: usize = 24;

    assert_ne!(
        N_TEMPLATES % TARGET_BATCH_COUNT,
        0,
        "the fixture must leave a partial final batch, else the drain arm is untested",
    );

    let paired = |qname: &[u8], flags: u16| {
        let mut b = SamBuilder::new();
        b.read_name(qname)
            .flags(flags)
            .ref_id(0)
            .pos(100)
            .cigar_ops(&[4u32 << 4])
            .sequence(b"ACGT")
            .qualities(&[30u8; 4]);
        b.build()
    };

    let mut records: Vec<(RawRecord, fgumi_bam_io::GroupKey)> = Vec::new();
    for t in 0..N_TEMPLATES {
        let qname = format!("q{t:05}");
        let key = fgumi_bam_io::GroupKey {
            name_hash: name_hash_of(qname.as_bytes()),
            ..Default::default()
        };
        records.push((paired(qname.as_bytes(), PAIRED | FIRST_SEGMENT), key));
        records.push((paired(qname.as_bytes(), PAIRED | LAST_SEGMENT), key));
    }
    let batches: Vec<DecodedRecordBatch> = records
        .chunks(RECORDS_PER_INPUT_BATCH)
        .enumerate()
        .map(|(i, chunk)| decoded_batch_with_keys(i as u64, chunk.to_vec()))
        .collect();

    let collected: Arc<Mutex<Vec<BamTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(batches))
        .chain(GroupByQueryname::with_target_batch_count(OUTPUT_BYTE_LIMIT, TARGET_BATCH_COUNT))
        .chain(StallingCollectSink { collected: sink_handle, remaining_stalls: SINK_STALLS })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let emitted = collected.lock().expect("mutex not poisoned");
    assert!(emitted.len() > 1, "fixture must span several output batches");

    let serials: Vec<u64> = emitted.iter().map(BamTemplateBatch::batch_serial).collect();
    assert_eq!(
        serials,
        (0..emitted.len() as u64).collect::<Vec<u64>>(),
        "emitted batch serials must be contiguous from zero",
    );

    // The batch cap: no emitted batch may exceed the target. Without this, a
    // regression that let the accumulator grow past the cap before emitting
    // would keep serials dense, keep every template present, and still leave a
    // short final batch — passing every other assertion here.
    if let Some((i, n)) = emitted
        .iter()
        .map(|b| b.templates().len())
        .enumerate()
        .find(|&(_, n)| n > TARGET_BATCH_COUNT)
    {
        panic!("batch {i} holds {n} templates, above the {TARGET_BATCH_COUNT} cap");
    }

    // Every queryname becomes exactly one template, in input order, holding
    // both of its records.
    let names: Vec<String> = emitted
        .iter()
        .flat_map(BamTemplateBatch::templates)
        .map(|t| String::from_utf8_lossy(&t.name).into_owned())
        .collect();
    let expected: Vec<String> = (0..N_TEMPLATES).map(|t| format!("q{t:05}")).collect();
    assert_eq!(
        names.len(),
        expected.len(),
        "every queryname must yield one template — a short count means the drain arm \
         dropped the run in progress",
    );
    if let Some((i, (got, want))) =
        names.iter().zip(&expected).enumerate().find(|(_, (got, want))| got != want)
    {
        panic!("template {i} differs: got {got:?}, expected {want:?}");
    }
    assert!(
        emitted.iter().flat_map(BamTemplateBatch::templates).all(|t| t.records.len() == 2),
        "each template must hold both of its records",
    );
}

/// Inputs for the closure-driven mid-step chains below.
fn numbered_inputs(n: u64) -> Vec<Numbered> {
    (0..n).map(|i| Numbered { ordinal: i, value: i }).collect()
}

/// Collect a sink's values, sorted — for the unordered (`BranchOrdering::None`)
/// variants, where arrival order races across workers.
fn sorted_values(sink: &Arc<Mutex<Vec<Numbered>>>) -> Vec<u64> {
    let mut v: Vec<u64> =
        sink.lock().expect("mutex not poisoned").iter().map(|n| n.value).collect();
    v.sort_unstable();
    v
}

const PROCESS_INPUTS: u64 = 200;

/// Small enough that a handful of `Numbered` items exceed it, so the ordered
/// mid-step chains below actually bounce a push and take their held-slot path.
const PROCESS_EDGE_BYTES: u64 = 256;

/// Ditto for the count-bounded variants.
const PROCESS_EDGE_SLOTS: usize = 2;

/// Ticks the terminal sink ignores before it starts draining, so the queue
/// upstream of it fills. See [`StallingCollectSink`].
const PROCESS_SINK_STALLS: usize = 64;

/// The closure-driven mid-steps in `process.rs` are all the same state machine
/// around a different output arity, and every one of their `try_run` bodies —
/// held-slot retry, drain detection, `Finished` on the last clone — was
/// reachable only through the framework. The unit tests in that module cover
/// the `Process*Output` constructors and profiles, not the steps.
///
/// These chains exercise one variant each. They are deliberately small and
/// uniform: the value being pinned is that each variant passes every input
/// through exactly once and terminates, which is what a drain or held-slot
/// regression breaks.
#[rstest]
fn process_chain_passes_every_item_through(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::process::process;

    let collected: Arc<Mutex<Vec<Numbered>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(numbered_inputs(PROCESS_INPUTS)))
        // `Process` is count-bounded and unordered; the doubling makes a
        // pass-through implementation distinguishable from a real one.
        .chain(process("Double", PROCESS_EDGE_SLOTS, |item: Numbered| {
            Ok(Numbered { ordinal: item.ordinal, value: item.value * 2 })
        }))
        .chain(StallingCollectSink {
            collected: sink_handle,
            remaining_stalls: PROCESS_SINK_STALLS,
        })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    assert_eq!(
        sorted_values(&collected),
        (0..PROCESS_INPUTS).map(|v| v * 2).collect::<Vec<u64>>(),
        "every input must be transformed and emitted exactly once",
    );
}

/// `ProcessOrdered` is the byte-bounded, `ByItemOrdinal` sibling of `Process`,
/// so its output must arrive in input order rather than merely intact.
#[rstest]
fn process_ordered_chain_preserves_input_order(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::process::process_ordered;

    let collected: Arc<Mutex<Vec<Numbered>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(numbered_inputs(PROCESS_INPUTS)))
        .chain(process_ordered("DoubleOrdered", PROCESS_EDGE_BYTES, |item: Numbered| {
            Ok(Numbered { ordinal: item.ordinal, value: item.value * 2 })
        }))
        .chain(StallingCollectSink {
            collected: sink_handle,
            remaining_stalls: PROCESS_SINK_STALLS,
        })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let values: Vec<u64> =
        collected.lock().expect("mutex not poisoned").iter().map(|n| n.value).collect();
    assert_eq!(
        values,
        (0..PROCESS_INPUTS).map(|v| v * 2).collect::<Vec<u64>>(),
        "an ordered branch must release in input order, not merely intact",
    );
}

/// `ProcessWithWorkerState` lazily initializes one state value per worker and
/// reuses it across items. The state here counts the items that worker handled,
/// so the per-worker totals must sum to the input count however the framework
/// distributes work.
#[rstest]
fn process_with_worker_state_chain_reuses_state_per_worker(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::process::process_with_worker_state;

    let collected: Arc<Mutex<Vec<Numbered>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let inits = Arc::new(AtomicU64::new(0));
    let init_counter = Arc::clone(&inits);
    let handled = Arc::new(AtomicU64::new(0));
    let handled_counter = Arc::clone(&handled);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(numbered_inputs(PROCESS_INPUTS)))
        .chain(process_with_worker_state(
            "CountPerWorker",
            PROCESS_EDGE_BYTES,
            move || {
                init_counter.fetch_add(1, AtomicOrdering::Relaxed);
                0u64 // per-worker item counter
            },
            move |seen: &mut u64, item: Numbered| {
                *seen += 1;
                handled_counter.fetch_add(1, AtomicOrdering::Relaxed);
                Ok(Numbered { ordinal: item.ordinal, value: item.value })
            },
        ))
        .chain(StallingCollectSink {
            collected: sink_handle,
            remaining_stalls: PROCESS_SINK_STALLS,
        })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let values: Vec<u64> =
        collected.lock().expect("mutex not poisoned").iter().map(|n| n.value).collect();
    assert_eq!(values, (0..PROCESS_INPUTS).collect::<Vec<u64>>(), "ordered output, in order");
    // State is created lazily per worker, so at most one per worker — never one
    // per item, which is the regression this construction exists to avoid.
    let inits = inits.load(AtomicOrdering::Relaxed);
    assert!(
        inits >= 1 && inits <= threads as u64,
        "expected between 1 and {threads} state initializations, got {inits}",
    );
    // The per-worker counters must sum to the input count however the framework
    // distributes work — the claim the doc comment makes, now pinned. Summing
    // through the shared accumulator is what distinguishes "each item handled by
    // some worker" from a drain regression that silently drops items.
    assert_eq!(
        handled.load(AtomicOrdering::Relaxed),
        PROCESS_INPUTS,
        "per-worker counters must total the input count",
    );
}

/// `MiAssign` threads a single running counter through the closure so emitted
/// molecule ids are globally consecutive. It is `Serial`, so the counter needs
/// no synchronization — and the assertion is that the ids come out dense and
/// in order regardless of thread count.
#[rstest]
fn mi_assign_chain_hands_out_consecutive_ids(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::process::mi_assign;

    const IDS_PER_ITEM: u64 = 2;

    let collected: Arc<Mutex<Vec<Numbered>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(numbered_inputs(PROCESS_INPUTS)))
        .chain(mi_assign("AssignIds", PROCESS_EDGE_BYTES, |next_mi: &mut u64, item: Numbered| {
            let base = *next_mi;
            *next_mi += IDS_PER_ITEM;
            Ok(Numbered { ordinal: item.ordinal, value: base })
        }))
        .chain(StallingCollectSink {
            collected: sink_handle,
            remaining_stalls: PROCESS_SINK_STALLS,
        })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let values: Vec<u64> =
        collected.lock().expect("mutex not poisoned").iter().map(|n| n.value).collect();
    assert_eq!(
        values,
        (0..PROCESS_INPUTS).map(|i| i * IDS_PER_ITEM).collect::<Vec<u64>>(),
        "each item must reserve its own block of ids, consecutively and in order",
    );
}

/// `Process2` fans out to two *unordered* branches, each with its own held
/// slot. Splitting on parity sends every input to exactly one branch, so a
/// dropped held item shows up as a missing value on one side.
#[rstest]
fn process2_chain_splits_across_both_branches(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::process::{Process2Output, process2};

    let collected_a: Arc<Mutex<Vec<Numbered>>> = Arc::new(Mutex::new(Vec::new()));
    let collected_b: Arc<Mutex<Vec<Numbered>>> = Arc::new(Mutex::new(Vec::new()));
    let (sink_a, sink_b) = (Arc::clone(&collected_a), Arc::clone(&collected_b));

    let builder = Pipeline::builder();
    let multi = builder
        .chain(ReplaySource::new(numbered_inputs(PROCESS_INPUTS)))
        .chain(process2("SplitParity", PROCESS_EDGE_SLOTS, PROCESS_EDGE_SLOTS, |item: Numbered| {
            Ok(if item.value.is_multiple_of(2) {
                Process2Output::only_a(item)
            } else {
                Process2Output::only_b(item)
            })
        }))
        .into_multi();
    multi
        .b0
        .chain(StallingCollectSink { collected: sink_a, remaining_stalls: PROCESS_SINK_STALLS })
        .into_sink_marker();
    multi
        .b1
        .chain(StallingCollectSink { collected: sink_b, remaining_stalls: PROCESS_SINK_STALLS })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    assert_eq!(
        sorted_values(&collected_a),
        (0..PROCESS_INPUTS).filter(|v| v.is_multiple_of(2)).collect::<Vec<u64>>(),
    );
    assert_eq!(
        sorted_values(&collected_b),
        (0..PROCESS_INPUTS).filter(|v| !v.is_multiple_of(2)).collect::<Vec<u64>>(),
    );
}

/// `Process2WithWorkerState` is the kept/rejects shape with per-worker scratch
/// — `correct`'s cache lives here. Both branches are ordered, and both receive
/// every item, so each must release a dense sequence.
#[rstest]
fn process2_with_worker_state_chain_feeds_both_ordered_branches(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::process::{Process2Output, process2_with_worker_state};

    let collected_a: Arc<Mutex<Vec<Numbered>>> = Arc::new(Mutex::new(Vec::new()));
    let collected_b: Arc<Mutex<Vec<Numbered>>> = Arc::new(Mutex::new(Vec::new()));
    let (sink_a, sink_b) = (Arc::clone(&collected_a), Arc::clone(&collected_b));

    let builder = Pipeline::builder();
    let multi = builder
        .chain(ReplaySource::new(numbered_inputs(PROCESS_INPUTS)))
        .chain(process2_with_worker_state(
            "KeptAndRejects",
            PROCESS_EDGE_BYTES,
            PROCESS_EDGE_BYTES,
            || 0u64,
            |seen: &mut u64, item: Numbered| {
                *seen += 1;
                Ok(Process2Output::both(
                    Numbered { ordinal: item.ordinal, value: item.value },
                    Numbered { ordinal: item.ordinal, value: item.value * 10 },
                ))
            },
        ))
        .into_multi();
    multi
        .b0
        .chain(StallingCollectSink { collected: sink_a, remaining_stalls: PROCESS_SINK_STALLS })
        .into_sink_marker();
    multi
        .b1
        .chain(StallingCollectSink { collected: sink_b, remaining_stalls: PROCESS_SINK_STALLS })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let values_a: Vec<u64> =
        collected_a.lock().expect("mutex not poisoned").iter().map(|n| n.value).collect();
    let values_b: Vec<u64> =
        collected_b.lock().expect("mutex not poisoned").iter().map(|n| n.value).collect();
    assert_eq!(values_a, (0..PROCESS_INPUTS).collect::<Vec<u64>>(), "branch a, in order");
    assert_eq!(
        values_b,
        (0..PROCESS_INPUTS).map(|v| v * 10).collect::<Vec<u64>>(),
        "branch b, in order",
    );
}

/// `Process3WithWorkerState` widens the same shape to three ordered branches —
/// the FASTQ-encode fan-out (R1 / R2 / other).
///
/// This step had no test because it had no way to be *wired*: `OrderedBytesTuple3`
/// carried queue and handle support but no `Chain::into_multi`, unlike its
/// 2-branch sibling, and `append_step` exposes only branch 0. The step
/// type-checked and built and was simply unreachable. `into_multi` for the
/// 3-branch ordered shape is added alongside this test.
#[rstest]
fn process3_with_worker_state_chain_feeds_all_three_branches(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::process::{Process3Output, process3_with_worker_state};

    let sinks: [Arc<Mutex<Vec<Numbered>>>; 3] =
        std::array::from_fn(|_| Arc::new(Mutex::new(Vec::new())));
    let [sink_a, sink_b, sink_c]: [Arc<Mutex<Vec<Numbered>>>; 3] =
        std::array::from_fn(|i| Arc::clone(&sinks[i]));

    let builder = Pipeline::builder();
    let multi = builder
        .chain(ReplaySource::new(numbered_inputs(PROCESS_INPUTS)))
        .chain(process3_with_worker_state(
            "FanOut3",
            PROCESS_EDGE_BYTES,
            PROCESS_EDGE_BYTES,
            PROCESS_EDGE_BYTES,
            || 0u64,
            |seen: &mut u64, item: Numbered| {
                *seen += 1;
                Ok(Process3Output::both3(
                    Numbered { ordinal: item.ordinal, value: item.value },
                    Numbered { ordinal: item.ordinal, value: item.value * 10 },
                    Numbered { ordinal: item.ordinal, value: item.value * 100 },
                ))
            },
        ))
        .into_multi();
    multi
        .b0
        .chain(StallingCollectSink { collected: sink_a, remaining_stalls: PROCESS_SINK_STALLS })
        .into_sink_marker();
    multi
        .b1
        .chain(StallingCollectSink { collected: sink_b, remaining_stalls: PROCESS_SINK_STALLS })
        .into_sink_marker();
    multi
        .b2
        .chain(StallingCollectSink { collected: sink_c, remaining_stalls: PROCESS_SINK_STALLS })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    for (branch, multiplier) in [(0usize, 1u64), (1, 10), (2, 100)] {
        let values: Vec<u64> =
            sinks[branch].lock().expect("mutex not poisoned").iter().map(|n| n.value).collect();
        assert_eq!(
            values,
            (0..PROCESS_INPUTS).map(|v| v * multiplier).collect::<Vec<u64>>(),
            "branch {branch} must release every item, in order",
        );
    }
}

/// `GroupBam` must emit the grouper's batches as formed, never merged.
///
/// [`TemplateGrouper::add_records`] already caps each returned batch at
/// `template_batch_size` and keeps the remainder pending, so the only way to
/// exceed the cap is for the step to `extend` those batches together — which is
/// what it used to do, turning the configured size from a bound into a mere
/// trigger (measured at 624x for a 40k-record input).
///
/// `RECORDS_PER_INPUT_BATCH` is deliberately several times `TARGET_BATCH_COUNT`
/// so one `add_records` call completes several whole batches. That is the case
/// the merge destroyed and the only one that distinguishes the two behaviors:
/// feed it fewer templates per input than the target and both the fixed and the
/// broken step emit one batch per input, and the test proves nothing.
///
/// This drives the step through `Pipeline::run` rather than calling the grouper
/// directly, because the property belongs to `GroupBam`, not to
/// `TemplateGrouper` — a grouper-level assertion passes even if the step
/// reintroduces the merge.
#[rstest]
fn group_bam_emits_grouper_batches_without_merging_them(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::group::bam::GroupBam;
    use crate::pipeline::steps::types::BamTemplateBatch;

    const N_TEMPLATES: usize = 205;
    const TARGET_BATCH_COUNT: usize = 8;
    const RECORDS_PER_INPUT_BATCH: usize = 40;
    const OUTPUT_BYTE_LIMIT: u64 = 1024;
    const SINK_STALLS: usize = 24;

    const _: () = assert!(
        RECORDS_PER_INPUT_BATCH > TARGET_BATCH_COUNT,
        "one input batch must complete more than one output batch or the merge is invisible",
    );
    assert_ne!(
        N_TEMPLATES % TARGET_BATCH_COUNT,
        0,
        "the fixture must leave a partial final batch, else the drain arm is untested",
    );

    // One single-end record per queryname, so templates and records are 1:1 and
    // a short count is unambiguous.
    let records: Vec<(RawRecord, fgumi_bam_io::GroupKey)> = (0..N_TEMPLATES)
        .map(|t| {
            let qname = format!("t{t:04}");
            let key = fgumi_bam_io::GroupKey {
                name_hash: name_hash_of(qname.as_bytes()),
                ..Default::default()
            };
            (single_end_record(qname.as_bytes()), key)
        })
        .collect();
    let batches: Vec<DecodedRecordBatch> = records
        .chunks(RECORDS_PER_INPUT_BATCH)
        .enumerate()
        .map(|(i, chunk)| decoded_batch_with_keys(i as u64, chunk.to_vec()))
        .collect();

    let collected: Arc<Mutex<Vec<BamTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
    let sink_handle = Arc::clone(&collected);
    let builder = Pipeline::builder();
    builder
        .chain(ReplaySource::new(batches))
        .chain(GroupBam::new(TARGET_BATCH_COUNT, OUTPUT_BYTE_LIMIT))
        .chain(StallingCollectSink { collected: sink_handle, remaining_stalls: SINK_STALLS })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    let emitted = collected.lock().expect("mutex not poisoned");
    assert!(emitted.len() > 1, "fixture must span several output batches");

    let serials: Vec<u64> = emitted.iter().map(Ordered::ordinal).collect();
    assert_eq!(
        serials,
        (0..emitted.len() as u64).collect::<Vec<u64>>(),
        "emitted batch serials must be contiguous from zero",
    );
    assert!(emitted.iter().all(|b| !b.templates().is_empty()), "no batch may be emitted empty");

    // The cap, and the whole point of the test: with the merge in place the
    // first batch of each input would hold `RECORDS_PER_INPUT_BATCH` templates.
    if let Some((i, n)) = emitted
        .iter()
        .map(|b| b.templates().len())
        .enumerate()
        .find(|&(_, n)| n > TARGET_BATCH_COUNT)
    {
        panic!("batch {i} holds {n} templates, above the {TARGET_BATCH_COUNT} cap");
    }

    // Nothing was dropped on the way — a short count would mean the `Finished`
    // gate closed the output while batches were still queued.
    let total: usize = emitted.iter().map(|b| b.templates().len()).sum();
    assert_eq!(
        total, N_TEMPLATES,
        "every template must be emitted exactly once; a short count means queued batches \
         were lost when the step reported Finished",
    );
}

/// A `None` on one of `Process3WithWorkerState`'s ordered branches must not
/// stall the run, provided that branch mints its own dense ordinals.
///
/// The sibling test above feeds all three branches on every input, so it never
/// exercises a gap. That matters because the step's named caller is the paired
/// FASTQ encode fan-out (R1 / R2 / **other**), and the "other" branch is `None`
/// for most inputs — the sparse case is the *normal* case in production, not an
/// edge case. All three branches declare [`BranchOrdering::ByItemOrdinal`],
/// whose `ReorderStage` releases only contiguously, so a branch that skipped an
/// input while propagating the *input's* ordinal would leave a permanent hole
/// and wait for it forever.
///
/// This pins the safe half of the rule documented on `Process2Output`: the
/// sparse branch here numbers what it emits, so a skip simply means the next
/// emitted item takes the next ordinal and the sequence stays dense. The
/// failure it guards against is a hang, not a wrong answer, which is why the
/// run completing at all is half the assertion.
#[rstest]
fn process3_tolerates_a_none_branch_that_mints_its_own_ordinals(#[values(1, 4)] threads: usize) {
    use crate::pipeline::steps::process::{Process3Output, process3_with_worker_state};

    // Only every fourth input reaches the third branch.
    const C_EVERY: u64 = 4;

    let sinks: [Arc<Mutex<Vec<Numbered>>>; 3] =
        std::array::from_fn(|_| Arc::new(Mutex::new(Vec::new())));
    let [sink_a, sink_b, sink_c]: [Arc<Mutex<Vec<Numbered>>>; 3] =
        std::array::from_fn(|i| Arc::clone(&sinks[i]));

    // Branch-local counter for the sparse branch. Shared across worker clones,
    // so at four threads which *value* gets which ordinal races — but the
    // sequence stays dense, which is the property the reorder stage requires.
    let next_c = Arc::new(AtomicU64::new(0));
    let mint_c = Arc::clone(&next_c);

    let builder = Pipeline::builder();
    let multi = builder
        .chain(ReplaySource::new(numbered_inputs(PROCESS_INPUTS)))
        .chain(process3_with_worker_state(
            "SparseFanOut3",
            PROCESS_EDGE_BYTES,
            PROCESS_EDGE_BYTES,
            PROCESS_EDGE_BYTES,
            || 0u64,
            move |seen: &mut u64, item: Numbered| {
                *seen += 1;
                // a and b take every input and propagate its ordinal, which
                // stays dense because nothing is skipped on those branches.
                let a = Numbered { ordinal: item.ordinal, value: item.value };
                let b = Numbered { ordinal: item.ordinal, value: item.value * 10 };
                let c = item.value.is_multiple_of(C_EVERY).then(|| Numbered {
                    ordinal: mint_c.fetch_add(1, AtomicOrdering::Relaxed),
                    value: item.value * 100,
                });
                Ok(Process3Output { a: Some(a), b: Some(b), c })
            },
        ))
        .into_multi();
    multi
        .b0
        .chain(StallingCollectSink { collected: sink_a, remaining_stalls: PROCESS_SINK_STALLS })
        .into_sink_marker();
    multi
        .b1
        .chain(StallingCollectSink { collected: sink_b, remaining_stalls: PROCESS_SINK_STALLS })
        .into_sink_marker();
    multi
        .b2
        .chain(StallingCollectSink { collected: sink_c, remaining_stalls: PROCESS_SINK_STALLS })
        .into_sink_marker();
    let pipeline = builder.build().expect("chain builds");
    pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

    // The dense branches are unaffected by the third one's gaps.
    for (branch, multiplier) in [(0usize, 1u64), (1, 10)] {
        let values: Vec<u64> =
            sinks[branch].lock().expect("mutex not poisoned").iter().map(|n| n.value).collect();
        assert_eq!(
            values,
            (0..PROCESS_INPUTS).map(|v| v * multiplier).collect::<Vec<u64>>(),
            "dense branch {branch} must release every item, in order",
        );
    }

    let branch_c = sinks[2].lock().expect("mutex not poisoned");
    let expected_c: Vec<u64> =
        (0..PROCESS_INPUTS).filter(|v| v.is_multiple_of(C_EVERY)).map(|v| v * 100).collect();

    // Sparse, or the test degenerates into the all-branches case above.
    assert!(
        (branch_c.len() as u64) < PROCESS_INPUTS,
        "the third branch must skip inputs or this proves nothing about gaps",
    );
    assert_eq!(branch_c.len(), expected_c.len(), "every routed item reached the sparse branch");

    // The load-bearing assertion: the sparse branch released a dense sequence
    // from zero. A hole would have stalled its reorder stage instead of
    // producing short output.
    let ordinals: Vec<u64> = branch_c.iter().map(|n| n.ordinal).collect();
    assert_eq!(
        ordinals,
        (0..branch_c.len() as u64).collect::<Vec<u64>>(),
        "the sparse branch must release a dense ordinal sequence from zero",
    );

    // Nothing was routed wrongly, dropped, or duplicated. Sorted because
    // ordinal assignment races across workers at four threads.
    let mut values_c: Vec<u64> = branch_c.iter().map(|n| n.value).collect();
    values_c.sort_unstable();
    assert_eq!(values_c, expected_c);
}

/// `EDGE_LIMIT_BYTES` must be small enough that the multi-record fixtures
/// exceed it, or the byte-bounded edges never reject a push and every chain
/// above silently degrades into an unbounded-queue test that proves nothing
/// about backpressure or the steps' held-item retry paths.
#[test]
fn the_edge_budget_binds_on_the_multi_record_fixtures() {
    let decompressed_bytes: usize = bgzf_blocks_for(&test_records(500), 4096)
        .iter()
        .map(|b| b.uncompressed_size as usize)
        .sum();
    assert!(
        u64::try_from(decompressed_bytes).expect("fixture volume fits u64") > EDGE_LIMIT_BYTES,
        "fixture volume ({decompressed_bytes} B) must exceed the per-edge budget \
         ({EDGE_LIMIT_BYTES} B), else the byte-bounded edges never bind",
    );
}
