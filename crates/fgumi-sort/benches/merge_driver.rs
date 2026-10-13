//! Microbenchmark for the slot-backed k-way merge driver
//! (`MergeDriver::from_slots`) as the pipeline's `SortMerge` step drives it:
//! step the driver and append every produced record to a reused output buffer,
//! as the merge's batch builder does.
//!
//! The slots are fed through their own API (`bp_commit_read`, then
//! `bp_insert_drain_finalize` with pre-built decompressed blocks), never by
//! pushing into the decompressed queue directly. Only the first
//! `PHASE2_DECOMP_CAP` blocks reach a slot's queue at once; on a stall the loop
//! drains the stalled slot's next ready blocks (`bp_drain_and_finalize`), as
//! the decompress step would. Records are 450 bytes (a WGS mean) with
//! coordinate keys that interleave across the runs, so the winner switches
//! sources on every record; blocks are cut every 65,280 bytes of framed
//! records, so some records span two blocks.
//!
//!     cargo bench -p fgumi-sort --bench merge_driver

use std::sync::Arc;

use criterion::{BatchSize, Criterion, Throughput, criterion_group, criterion_main};
use fgumi_raw_bam::testutil::make_bam_bytes;
use fgumi_sort::{
    MemorySources, MergeDriver, MergeDriverDyn, MergeStep, RawCoordinateKey, SortMergeSlot,
    SpillCodec, frame_keyed_record_into,
};

/// Records merged per iteration.
const RECORDS: usize = 1_000_000;
/// Bytes per record body.
const RECORD_LEN: usize = 450;
/// Decompressed bytes per block (a full BGZF block's payload).
const BLOCK_BYTES: usize = 65_280;
/// Output buffer size at which the consumer starts a new buffer.
const OUTPUT_BYTES: usize = 64 * 1024;

/// A 450-byte coordinate-sortable BAM record at (tid 0, `pos`).
fn record(i: usize, pos: i32) -> Vec<u8> {
    let name = format!("r{i:09}");
    let fixed = 32 + name.len() + 1;
    // One Z tag pads the body to exactly RECORD_LEN.
    let mut aux = b"XZZ".to_vec();
    aux.extend(std::iter::repeat_n(b'x', RECORD_LEN - fixed - 4));
    aux.push(0);
    let body = make_bam_bytes(0, pos, 0, name.as_bytes(), &[], 0, -1, -1, &aux);
    assert_eq!(body.len(), RECORD_LEN);
    body
}

/// The decompressed blocks of `k` sorted runs: record `i` goes to run `i % k`
/// at position `i`, so the runs interleave record by record.
fn run_blocks(k: usize) -> Vec<Vec<Vec<u8>>> {
    let mut framed: Vec<Vec<u8>> = vec![Vec::new(); k];
    for i in 0..RECORDS {
        let pos = i32::try_from(i).expect("position fits i32");
        frame_keyed_record_into(&mut framed[i % k], &RawCoordinateKey::default(), &record(i, pos))
            .expect("frame record");
    }
    framed.iter().map(|run| run.chunks(BLOCK_BYTES).map(<[u8]>::to_vec).collect()).collect()
}

/// One slot per run, every block committed and inserted through the slot API.
fn slots(blocks: &[Vec<Vec<u8>>]) -> Vec<Arc<SortMergeSlot>> {
    blocks
        .iter()
        .enumerate()
        .map(|(id, run)| {
            let file = tempfile::tempfile().expect("tempfile");
            let id = u32::try_from(id).expect("few runs");
            let slot = Arc::new(SortMergeSlot::new(id, Arc::new(file), 0, 0, SpillCodec::Bgzf));
            slot.bp_commit_read(run.len(), true);
            slot.bp_insert_drain_finalize(0, run.clone(), run.len());
            slot
        })
        .collect()
}

/// Merge every record, appending each to a reused output buffer; returns the
/// record count.
fn merge(slots: Vec<Arc<SortMergeSlot>>, out: &mut Vec<u8>) -> u64 {
    let mut driver = MergeDriver::<RawCoordinateKey>::from_slots(
        slots,
        MemorySources::Shared(Vec::new()),
        RECORDS as u64,
    );
    loop {
        match driver.try_step().expect("merge step") {
            MergeStep::Produced(bytes) => {
                out.extend_from_slice(bytes);
                if out.len() >= OUTPUT_BYTES {
                    out.clear();
                }
            }
            MergeStep::Stalled => {
                let slot = driver.stalled_slot().expect("a slot stalled");
                assert!(slot.bp_drain_and_finalize(), "a stalled slot has blocks to drain");
            }
            MergeStep::Done => return driver.records_merged(),
        }
    }
}

fn bench_merge_driver(c: &mut Criterion) {
    let mut group = c.benchmark_group("merge_driver");
    group.throughput(Throughput::Elements(RECORDS as u64));
    group.sample_size(10);
    for k in [2usize, 27] {
        let blocks = run_blocks(k);
        let mut out = Vec::with_capacity(OUTPUT_BYTES + RECORD_LEN);
        group.bench_function(format!("k{k}_450b_1m"), |b| {
            b.iter_batched(
                || slots(&blocks),
                |slots| assert_eq!(merge(slots, &mut out), RECORDS as u64),
                BatchSize::PerIteration,
            );
        });
    }
    group.finish();
}

criterion_group!(benches, bench_merge_driver);
criterion_main!(benches);
