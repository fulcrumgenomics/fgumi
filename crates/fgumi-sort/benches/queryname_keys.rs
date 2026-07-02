//! Microbenchmarks isolating the queryname sort's per-record key cost, to
//! decide two things before touching the production queryname path:
//!
//!   1. **Prefix key (solution C)** — does packing the first 8 name bytes into a
//!      `u64` for a cheap first-pass compare actually help? The answer depends
//!      entirely on the name distribution: Illumina read names share a long
//!      common prefix (`instrument:run:flowcell:lane:` is identical for every
//!      read in a lane) and only vary in the `tile:x:y` suffix — so a *front*
//!      prefix key may discriminate nothing and add pure overhead. We therefore
//!      bench with REALISTIC Illumina-style names, not random ones.
//!
//!   2. **Null-key ceilings** — the wall-time floor if the key were free:
//!      - `*_ingest`: owned-name extraction (alloc + memcpy per record) vs a
//!        cheap `u64` extraction (no alloc). The delta is the ceiling on what
//!        solution A (borrow the name from the arena, drop the per-record Vec)
//!        can recover on the ingest side.
//!      - `*_sort`: full-name comparator vs a cheap `u64` comparator. The delta
//!        is the ceiling on what a cheaper comparator (C, or radix on a packed
//!        key) can recover on the sort side.
//!
//!     cargo bench -p fgumi-sort --bench queryname_keys

use std::cmp::Ordering;

use criterion::{BatchSize, Criterion, Throughput, criterion_group, criterion_main};
use fgumi_raw_bam::testutil::make_bam_bytes;
use fgumi_sort::{RawQuerynameKey, RawQuerynameLexKey, RawSortKey};

/// Number of records per bench iteration.
const N: usize = 1_000_000;

/// Build `n` realistic Illumina-style read names: a fixed
/// `instrument:run:flowcell:lane:` prefix shared by every read (22 bytes), then
/// a varying `tile:x:y` suffix. This is the adversarial case for a front-prefix
/// key — the first ~22 bytes are identical across all records.
fn illumina_names(n: usize) -> Vec<Vec<u8>> {
    (0..n)
        .map(|i| {
            // Vary tile (1101..1112), x (0..) and y (0..) like a real lane.
            let tile = 1101 + (i % 12);
            let x = (i.wrapping_mul(2_654_435_761)) % 30_000;
            let y = (i.wrapping_mul(40_503)) % 30_000;
            format!("A00123:45:HGVWXDSXY:1:{tile}:{x}:{y}").into_bytes()
        })
        .collect()
}

/// Pack the first 8 bytes of `name` into a big-endian `u64` so numeric ordering
/// of the u64 matches lexicographic ordering of the first 8 bytes (short names
/// zero-pad on the right).
#[inline]
fn prefix8(name: &[u8]) -> u64 {
    let mut buf = [0u8; 8];
    let take = name.len().min(8);
    buf[..take].copy_from_slice(&name[..take]);
    u64::from_be_bytes(buf)
}

/// Build BAM records carrying `names` (paired R1 flag), so key extraction walks
/// a real record body exactly as production does.
fn bam_records(names: &[Vec<u8>]) -> Vec<Vec<u8>> {
    names.iter().map(|name| make_bam_bytes(0, 0, 0, name, &[], 100, -1, -1, &[])).collect()
}

fn bench_ingest(c: &mut Criterion) {
    let names = illumina_names(N);
    let records = bam_records(&names);

    let mut group = c.benchmark_group("queryname_ingest_1M");
    group.throughput(Throughput::Elements(N as u64));
    group.sample_size(10);

    // Owned natural key: alloc + memcpy of the name per record (current path).
    group.bench_function("natural_owned_extract", |b| {
        b.iter(|| {
            let keys: Vec<RawQuerynameKey> =
                records.iter().map(|r| RawQuerynameKey::extract_from_record(r)).collect();
            std::hint::black_box(keys)
        });
    });

    // Owned lex key: same alloc + memcpy.
    group.bench_function("lex_owned_extract", |b| {
        b.iter(|| {
            let keys: Vec<RawQuerynameLexKey> =
                records.iter().map(|r| RawQuerynameLexKey::extract_from_record(r)).collect();
            std::hint::black_box(keys)
        });
    });

    // Null key: a cheap u64 prefix, NO allocation. The delta vs the owned
    // extracts is the ceiling on solution A (borrow, drop the per-record Vec).
    group.bench_function("null_u64_extract", |b| {
        b.iter(|| {
            let keys: Vec<u64> = records
                .iter()
                .map(|r| {
                    // Mirror the name-offset walk production does, then prefix8.
                    let key = RawQuerynameLexKey::extract_from_record(r);
                    prefix8(key.name())
                })
                .collect();
            std::hint::black_box(keys)
        });
    });

    group.finish();
}

fn bench_sort(c: &mut Criterion) {
    let names = illumina_names(N);
    let records = bam_records(&names);

    // Pre-build the three key representations once, outside the timed loop.
    let nat_keys: Vec<RawQuerynameKey> =
        records.iter().map(|r| RawQuerynameKey::extract_from_record(r)).collect();
    let lex_keys: Vec<RawQuerynameLexKey> =
        records.iter().map(|r| RawQuerynameLexKey::extract_from_record(r)).collect();
    // Prefix + full: (u64 prefix, full lex key). Prefix compared first, full
    // name only on prefix ties.
    let prefix_keys: Vec<(u64, RawQuerynameLexKey)> = records
        .iter()
        .map(|r| {
            let k = RawQuerynameLexKey::extract_from_record(r);
            (prefix8(k.name()), k)
        })
        .collect();
    let u64_keys: Vec<u64> = prefix_keys.iter().map(|(p, _)| *p).collect();

    let mut group = c.benchmark_group("queryname_sort_1M");
    group.throughput(Throughput::Elements(N as u64));
    group.sample_size(10);

    // Full natural comparator (samtools strnum_cmp) — current path.
    group.bench_function("natural_full_cmp", |b| {
        b.iter_batched(
            || nat_keys.clone(),
            |mut keys| {
                keys.sort_unstable();
                std::hint::black_box(keys)
            },
            BatchSize::PerIteration,
        );
    });

    // Full lex comparator (Vec<u8>::cmp) — current lex path.
    group.bench_function("lex_full_cmp", |b| {
        b.iter_batched(
            || lex_keys.clone(),
            |mut keys| {
                keys.sort_unstable();
                std::hint::black_box(keys)
            },
            BatchSize::PerIteration,
        );
    });

    // Prefix-then-full: cheap u64 compare first, full name only on prefix ties.
    // On Illumina names the shared front prefix means this almost always falls
    // through to the full compare — the bench proves whether it helps or hurts.
    group.bench_function("lex_prefix_then_full_cmp", |b| {
        b.iter_batched(
            || prefix_keys.clone(),
            |mut keys| {
                keys.sort_unstable_by(|a, b| match a.0.cmp(&b.0) {
                    Ordering::Equal => a.1.cmp(&b.1),
                    other => other,
                });
                std::hint::black_box(keys)
            },
            BatchSize::PerIteration,
        );
    });

    // Null key: sort the u64 prefix alone. The comparator ceiling — no full
    // name ever compared.
    group.bench_function("null_u64_cmp", |b| {
        b.iter_batched(
            || u64_keys.clone(),
            |mut keys| {
                keys.sort_unstable();
                std::hint::black_box(keys)
            },
            BatchSize::PerIteration,
        );
    });

    group.finish();
}

criterion_group!(benches, bench_ingest, bench_sort);
criterion_main!(benches);
