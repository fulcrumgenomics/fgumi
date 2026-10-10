#![deny(unsafe_code)]

//! The spill gather's per-run cost: planning a sorted chunk into block ranges
//! from its length index (`plan_ranges`, what `SpillGather` does on the
//! coordination thread) against framing every block serially
//! (`frame_blocks`, what it did before the pool framed the ranges).
//!
//! About 5 M records at random offsets over a 2 GiB arena, so the serial
//! framing pays the cold-arena reads a real run does. `--test` shrinks the
//! arena to 256 MiB and the record count with it.
//!
//!     cargo bench -p fgumi-pipeline-io --features test-utils --bench spill_gather

use std::sync::Arc;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use fgumi_bgzf::BGZF_MAX_BLOCK_SIZE;
use fgumi_pipeline_io::sort::protocol::MemoryChunkErased;
use fgumi_pipeline_io::sort::spill_gather::{bench_frame_blocks, bench_plan_ranges};
use fgumi_sort::{PooledSegmentedBuf, RecordRef, SegmentedBuf, coordinate_chunk_from_refs};

/// A coordinate chunk of `records` records of 300-ish bytes at pseudo-random
/// offsets in an arena of `arena_bytes`.
fn chunk(arena_bytes: usize, records: usize) -> MemoryChunkErased {
    let mut arena = SegmentedBuf::new();
    let block = vec![0x5Au8; 1 << 20];
    while arena.len() < arena_bytes {
        arena.extend_from_slice(&block);
    }
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let refs: Vec<RecordRef> = (0..records as u64)
        .map(|i| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let len = 250 + u32::try_from(state % 120).expect("fits");
            let offset = state % (arena_bytes as u64 - u64::from(len));
            RecordRef::new(i, offset, len)
        })
        .collect();
    MemoryChunkErased::Coordinate(coordinate_chunk_from_refs(
        Arc::new(PooledSegmentedBuf::unpooled(arena)),
        refs,
        1,
    ))
}

fn bench_spill_gather(c: &mut Criterion) {
    let quick = std::env::args().any(|a| a == "--test");
    let (arena_bytes, records) =
        if quick { (256 << 20, 625_000) } else { (2usize << 30, 5_000_000) };
    let chunk = chunk(arena_bytes, records);

    let mut group = c.benchmark_group("spill_gather");
    group.throughput(Throughput::Elements(records as u64));
    group.sample_size(10);
    group.bench_function("plan_ranges", |b| {
        b.iter(|| std::hint::black_box(bench_plan_ranges(&chunk, BGZF_MAX_BLOCK_SIZE)));
    });
    group.bench_function("frame_blocks", |b| {
        b.iter(|| {
            std::hint::black_box(bench_frame_blocks(&chunk, BGZF_MAX_BLOCK_SIZE).expect("frames"))
        });
    });
    group.finish();
}

criterion_group!(benches, bench_spill_gather);
criterion_main!(benches);
