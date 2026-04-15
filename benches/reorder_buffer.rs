//! Benchmarks for `ReorderBuffer::insert_with_size` under varying gap patterns.
//!
//! This exists so future data-structure swaps (e.g. `BTreeMap`, `BinaryHeap`,
//! ring buffer) can be measured honestly against the current sparse-`VecDeque`
//! implementation. The sparse `VecDeque` makes inserts O(gap) when a sequence
//! number arrives ahead of the current buffer end (the gap is filled with
//! `None` sentinels), so the cost is sensitive to gap size — the scenarios
//! below sweep a few realistic steady-state patterns.
//!
//! Run with: `cargo bench --bench reorder_buffer`
//! Quick sanity: `cargo bench --bench reorder_buffer -- --quick`
//! View reports in: `target/criterion/report/index.html`
#![deny(unsafe_code)]
#![allow(clippy::cast_possible_truncation)]

use std::hint::black_box;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

use fgumi_lib::reorder_buffer::ReorderBuffer;

/// Number of inserts per bench iteration. Large enough that per-call overhead
/// is amortized, small enough to keep the bench quick.
const BATCH_SIZE: usize = 256;

/// Run a fill+drain cycle where inserts arrive in perfect sequence order.
///
/// With `k == 0`, no gap slots are ever pushed — this is the O(1)-per-insert
/// best case and a baseline for the other scenarios.
fn bench_in_order(buffer: &mut ReorderBuffer<u64>, batch_size: usize) {
    let base = buffer.next_seq();
    for i in 0..batch_size as u64 {
        let seq = base + i;
        buffer.insert_with_size(black_box(seq), black_box(seq), black_box(8));
    }
    // Drain everything so the next iteration starts from a known empty state.
    while let Some((val, size)) = buffer.try_pop_next_with_size() {
        black_box((val, size));
    }
}

/// Insert in chunks of size `k+1` where the far slot is filled first, then the
/// near slots 0..k are filled in order, then drain.
///
/// This models a slow producer that is `k` ahead of the in-order stream:
/// - Chunk layout for k=1: insert [base+1, base+0], drain 2 items
/// - Chunk layout for k=8: insert [base+8, base+0..=7], drain 9 items
/// - Chunk layout for k=64: insert [base+64, base+0..=63], drain 65 items
///
/// The first insert in each chunk walks k slots ahead, pushing k `None`
/// sentinels, so each chunk pays O(k) for the lead insert plus O(1) for the
/// rest — this is exactly the pattern the fixed doc comment is warning about.
fn bench_leading_gap(buffer: &mut ReorderBuffer<u64>, batch_size: usize, k: u64) {
    let chunk = k as usize + 1;
    let mut base = buffer.next_seq();
    let mut remaining = batch_size;
    while remaining > 0 {
        let this_chunk = remaining.min(chunk);
        let lead = this_chunk as u64 - 1;
        // Insert the far slot first (walks `lead` gap slots forward).
        buffer.insert_with_size(black_box(base + lead), black_box(base + lead), black_box(8));
        // Fill the near slots 0..lead in order.
        for i in 0..lead {
            let seq = base + i;
            buffer.insert_with_size(black_box(seq), black_box(seq), black_box(8));
        }
        // Drain this chunk before starting the next so steady-state depth
        // stays bounded and the pattern is repeatable.
        while let Some((val, size)) = buffer.try_pop_next_with_size() {
            black_box((val, size));
        }
        base += this_chunk as u64;
        remaining -= this_chunk;
    }
}

fn bench_reorder_buffer_insert(c: &mut Criterion) {
    let mut group = c.benchmark_group("reorder_buffer_insert");
    group.throughput(Throughput::Elements(BATCH_SIZE as u64));

    // Scenario 1: purely in-order (k=0).
    group.bench_function(BenchmarkId::from_parameter("in_order_k0"), |b| {
        b.iter_batched_ref(
            ReorderBuffer::<u64>::new,
            |buffer| bench_in_order(buffer, BATCH_SIZE),
            criterion::BatchSize::SmallInput,
        );
    });

    // Scenarios 2-4: leading-gap patterns for k in {1, 8, 64}.
    for &k in &[1u64, 8, 64] {
        let id = BenchmarkId::new("leading_gap", k);
        group.bench_with_input(id, &k, |b, &k| {
            b.iter_batched_ref(
                ReorderBuffer::<u64>::new,
                |buffer| bench_leading_gap(buffer, BATCH_SIZE, k),
                criterion::BatchSize::SmallInput,
            );
        });
    }

    group.finish();
}

criterion_group!(benches, bench_reorder_buffer_insert);
criterion_main!(benches);
