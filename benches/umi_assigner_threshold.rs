//! Microbenchmark: parallel vs sequential UMI assigner at varying group sizes.
//!
//! Produces the per-strategy crossover points backing the `auto` mode of the
//! `--parallel-group-min-templates` flag in `src/lib/commands/group.rs` (see
//! `parallel_threshold`). For each strategy the gate decides, per position
//! group, whether enough templates are present to justify the parallel
//! assigner; the threshold differs per strategy because the per-template work
//! differs. The hypothesis under test: building a per-group rayon `ThreadPool`
//! (`ParallelXAssigner::new(threads)`) is expensive enough (~120 us) that for
//! small groups the sequential assigner wins.
//!
//! Each bench iteration measures **construct + assign** as a single unit,
//! because the threshold gate's job is to avoid pool construction. Measuring
//! `assign()` alone would hide the cost we're trying to eliminate.
//!
//! Run with: `cargo bench --bench umi_assigner_threshold`
//! Report:   `target/criterion/report/index.html`

#![allow(clippy::cast_possible_truncation)]

use std::hint::black_box;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

use fgumi_lib::umi::assigner::{Strategy, Umi, UmiAssigner};
use fgumi_lib::umi::parallel_assigner::{
    ParallelAdjacencyAssigner, ParallelEditAssigner, ParallelIdentityAssigner,
    ParallelPairedAssigner,
};

/// Group sizes with extra resolution around expected cross-over regions:
///   - Paired: ~128-256 (per earlier coarse bench)
///   - Edit / Adjacency: ~1024-4096
const GROUP_SIZES: &[usize] = &[
    64, 96, 128, 160, 192, 224, 256, 384, 512, 768, 1024, 1536, 2048, 2560, 3072, 4096, 8192, 16384,
];

/// UMI length per single read. Paired UMIs are `LEN-LEN` so the joined string
/// is `2*LEN + 1` chars.
const UMI_LEN: usize = 8;

/// Number of duplicate copies per unique UMI sequence. Real UMI families have
/// 1-20 reads; 8 puts us in the realistic middle for WGS/WES workloads.
const COPIES_PER_UMI: usize = 8;

/// Deterministic seed so all bench inputs are reproducible.
const SEED: u64 = 0xFEED_BEEF;

/// Number of threads the parallel assigner pool spawns. Mirrors what
/// `Command::execute` passes (`num_threads`) on a typical workstation.
fn bench_threads() -> usize {
    num_cpus::get().clamp(2, 8)
}

/// Generate `n_templates` UMIs of length `UMI_LEN` with ~`COPIES_PER_UMI`
/// reads per unique sequence and a small dose of single-base errors. This
/// mimics the input the assigner sees inside one position group.
fn make_umis(n_templates: usize) -> Vec<Umi> {
    let mut rng = StdRng::seed_from_u64(SEED);
    let bases = [b'A', b'C', b'G', b'T'];
    let n_unique = n_templates.div_ceil(COPIES_PER_UMI).max(1);

    let unique: Vec<String> = (0..n_unique)
        .map(|_| (0..UMI_LEN).map(|_| bases[rng.random_range(0..4)] as char).collect::<String>())
        .collect();

    let mut umis = Vec::with_capacity(n_templates);
    while umis.len() < n_templates {
        let pick = &unique[rng.random_range(0..unique.len())];
        // 10% of the time inject a single-base substitution so the
        // edit-distance / adjacency strategies have real clustering work.
        if rng.random_range(0..10) == 0 {
            let mut bytes = pick.clone().into_bytes();
            let pos = rng.random_range(0..bytes.len());
            let cur = bytes[pos];
            let new = *[b'A', b'C', b'G', b'T'].iter().find(|&&b| b != cur).unwrap();
            bytes[pos] = new;
            umis.push(String::from_utf8(bytes).unwrap());
        } else {
            umis.push(pick.clone());
        }
    }
    umis
}

/// Generate paired UMIs in the `AAAAAAAA-CCCCCCCC` format that the Paired
/// strategy expects.
fn make_paired_umis(n_templates: usize) -> Vec<Umi> {
    let left = make_umis(n_templates);
    // Use a different seed offset so the right halves are independent.
    let mut rng = StdRng::seed_from_u64(SEED ^ 0xA5A5_A5A5);
    let bases = [b'A', b'C', b'G', b'T'];
    left.into_iter()
        .map(|l| {
            let r: String = (0..UMI_LEN).map(|_| bases[rng.random_range(0..4)] as char).collect();
            format!("{l}-{r}")
        })
        .collect()
}

/// Index threshold passed to the sequential adjacency/paired assigners. The
/// `group.rs` call sites pass `self.index_threshold` (default 100) so we
/// match that here.
const INDEX_THRESHOLD: usize = 100;

fn bench_identity(c: &mut Criterion) {
    let mut group = c.benchmark_group("umi_assigner/identity");
    let threads = bench_threads();
    for &n in GROUP_SIZES {
        let umis = make_umis(n);
        group.throughput(Throughput::Elements(n as u64));

        group.bench_with_input(BenchmarkId::new("sequential", n), &umis, |b, umis| {
            b.iter(|| {
                let assigner = Strategy::Identity.new_assigner_full(0, 1, INDEX_THRESHOLD);
                black_box(assigner.assign(black_box(umis)))
            });
        });

        group.bench_with_input(BenchmarkId::new("parallel", n), &umis, |b, umis| {
            b.iter(|| {
                let assigner = ParallelIdentityAssigner::new(threads);
                black_box(assigner.assign(black_box(umis)))
            });
        });
    }
    group.finish();
}

fn bench_edit(c: &mut Criterion) {
    let mut group = c.benchmark_group("umi_assigner/edit_1");
    let threads = bench_threads();
    for &n in GROUP_SIZES {
        let umis = make_umis(n);
        group.throughput(Throughput::Elements(n as u64));

        group.bench_with_input(BenchmarkId::new("sequential", n), &umis, |b, umis| {
            b.iter(|| {
                let assigner = Strategy::Edit.new_assigner_full(1, 1, INDEX_THRESHOLD);
                black_box(assigner.assign(black_box(umis)))
            });
        });

        group.bench_with_input(BenchmarkId::new("parallel", n), &umis, |b, umis| {
            b.iter(|| {
                let assigner = ParallelEditAssigner::new(1, threads);
                black_box(assigner.assign(black_box(umis)))
            });
        });
    }
    group.finish();
}

fn bench_adjacency(c: &mut Criterion) {
    let mut group = c.benchmark_group("umi_assigner/adjacency_1");
    let threads = bench_threads();
    for &n in GROUP_SIZES {
        let umis = make_umis(n);
        group.throughput(Throughput::Elements(n as u64));

        group.bench_with_input(BenchmarkId::new("sequential", n), &umis, |b, umis| {
            b.iter(|| {
                let assigner = Strategy::Adjacency.new_assigner_full(1, 1, INDEX_THRESHOLD);
                black_box(assigner.assign(black_box(umis)))
            });
        });

        group.bench_with_input(BenchmarkId::new("parallel", n), &umis, |b, umis| {
            b.iter(|| {
                let assigner = ParallelAdjacencyAssigner::new(1, threads);
                black_box(assigner.assign(black_box(umis)))
            });
        });
    }
    group.finish();
}

fn bench_paired(c: &mut Criterion) {
    let mut group = c.benchmark_group("umi_assigner/paired_1");
    let threads = bench_threads();
    for &n in GROUP_SIZES {
        let umis = make_paired_umis(n);
        group.throughput(Throughput::Elements(n as u64));

        group.bench_with_input(BenchmarkId::new("sequential", n), &umis, |b, umis| {
            b.iter(|| {
                let assigner = Strategy::Paired.new_assigner_full(1, 1, INDEX_THRESHOLD);
                black_box(assigner.assign(black_box(umis)))
            });
        });

        group.bench_with_input(BenchmarkId::new("parallel", n), &umis, |b, umis| {
            b.iter(|| {
                let assigner = ParallelPairedAssigner::new(1, threads);
                black_box(assigner.assign(black_box(umis)))
            });
        });
    }
    group.finish();
}

/// Isolate the cost of `ParallelXAssigner::new(threads)` — i.e. the rayon
/// `ThreadPoolBuilder::build()` call — versus the trivial sequential
/// constructor. This is the cost the threshold gate is designed to avoid.
fn bench_constructor_only(c: &mut Criterion) {
    let mut group = c.benchmark_group("umi_assigner/constructor_only");
    let threads = bench_threads();

    group.bench_function("sequential_identity", |b| {
        b.iter(|| {
            black_box(Strategy::Identity.new_assigner_full(0, 1, INDEX_THRESHOLD));
        });
    });
    group.bench_function("parallel_identity", |b| {
        b.iter(|| {
            black_box(ParallelIdentityAssigner::new(threads));
        });
    });
    group.bench_function("sequential_edit", |b| {
        b.iter(|| {
            black_box(Strategy::Edit.new_assigner_full(1, 1, INDEX_THRESHOLD));
        });
    });
    group.bench_function("parallel_edit", |b| {
        b.iter(|| {
            black_box(ParallelEditAssigner::new(1, threads));
        });
    });
    group.bench_function("sequential_adjacency", |b| {
        b.iter(|| {
            black_box(Strategy::Adjacency.new_assigner_full(1, 1, INDEX_THRESHOLD));
        });
    });
    group.bench_function("parallel_adjacency", |b| {
        b.iter(|| {
            black_box(ParallelAdjacencyAssigner::new(1, threads));
        });
    });
    group.bench_function("sequential_paired", |b| {
        b.iter(|| {
            black_box(Strategy::Paired.new_assigner_full(1, 1, INDEX_THRESHOLD));
        });
    });
    group.bench_function("parallel_paired", |b| {
        b.iter(|| {
            black_box(ParallelPairedAssigner::new(1, threads));
        });
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_constructor_only,
    bench_identity,
    bench_edit,
    bench_adjacency,
    bench_paired,
);
criterion_main!(benches);
