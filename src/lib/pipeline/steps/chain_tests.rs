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
use std::sync::{Arc, Mutex};

use fgumi_bam_io::GroupKeyConfig;
use fgumi_raw_bam::{RawRecord, SamBuilder};
use rstest::rstest;

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
use crate::pipeline::steps::parse::bam::ParseBamRecords;
use crate::pipeline::steps::parse::decode::{DecodeFromRecords, DecodeRecords};
use crate::pipeline::steps::types::{BgzfBlock, DecodedRecordBatch, RecordBatch};

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
