//! Benchmarks for core fgumi functions.
//!
//! Run with: `cargo bench`
//! View reports in: `target/criterion/report/index.html`

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use std::hint::black_box;

use bstr::BString;
use fgumi_lib::bitenc::BitEnc;
use fgumi_lib::consensus::base_builder::ConsensusBaseBuilder;
use fgumi_lib::consensus::caller::ConsensusCaller;
use fgumi_lib::consensus::vanilla_caller::{VanillaUmiConsensusCaller, VanillaUmiConsensusOptions};
use fgumi_lib::dna::{complement_base, reverse_complement};
use fgumi_lib::phred::{
    ln_error_prob_two_trials, ln_sum_exp, phred_to_ln_correct_prob, phred_to_ln_error_prob,
};
use fgumi_lib::sam::record_utils::parse_cigar_string;
use fgumi_lib::umi::assigner::{Strategy, count_mismatches, matches_within_threshold};
use fgumi_lib::vendored::bam_codec::encode_record_buf;
use noodles::core::Position;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::{Cigar, QualityScores, RecordBuf, Sequence};

/// Benchmark PHRED score conversions
fn bench_phred_conversions(c: &mut Criterion) {
    let mut group = c.benchmark_group("phred_conversions");

    // Benchmark direct calculation
    group.bench_function("phred_to_ln_error_direct", |b| {
        b.iter(|| {
            for q in 0..=93u8 {
                black_box(phred_to_ln_error_prob(q));
            }
        });
    });

    group.bench_function("phred_to_ln_correct_direct", |b| {
        b.iter(|| {
            for q in 0..=93u8 {
                black_box(phred_to_ln_correct_prob(q));
            }
        });
    });

    group.finish();
}

/// Benchmark log-space probability operations
fn bench_log_probability_ops(c: &mut Criterion) {
    let mut group = c.benchmark_group("log_probability_ops");

    let ln_p1 = 0.1_f64.ln();
    let ln_p2 = 0.01_f64.ln();

    group.bench_function("ln_sum_exp", |b| {
        b.iter(|| black_box(ln_sum_exp(black_box(ln_p1), black_box(ln_p2))));
    });

    group.bench_function("ln_error_prob_two_trials", |b| {
        b.iter(|| black_box(ln_error_prob_two_trials(black_box(ln_p1), black_box(ln_p2))));
    });

    group.finish();
}

/// Benchmark DNA operations
fn bench_dna_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("dna_operations");

    // Benchmark complement_base
    group.bench_function("complement_base", |b| {
        b.iter(|| {
            for base in [b'A', b'C', b'G', b'T', b'N'] {
                black_box(complement_base(base));
            }
        });
    });

    // Benchmark reverse complement with different sequence lengths
    for len in [50, 100, 150, 200, 300] {
        let seq: Vec<u8> = (0..len).map(|i| b"ACGT"[i % 4]).collect();
        group.throughput(Throughput::Bytes(len as u64));
        group.bench_with_input(BenchmarkId::new("reverse_complement", len), &seq, |b, seq| {
            b.iter(|| black_box(reverse_complement(black_box(seq))));
        });
    }

    group.finish();
}

/// Benchmark CIGAR parsing
fn bench_cigar_parsing(c: &mut Criterion) {
    let mut group = c.benchmark_group("cigar_parsing");

    let cigars = [
        ("simple_50M", "50M"),
        ("with_clips", "5S50M5S"),
        ("with_indels", "25M1I10M2D15M"),
        ("complex", "5H10S30M5I20M3D25M10S5H"),
        ("very_complex", "2H5S10M1I5M2D8M1I6M3D10M2I7M1D15M5S3H"),
    ];

    for (name, cigar) in cigars {
        group.bench_with_input(BenchmarkId::new("parse_cigar_string", name), cigar, |b, cigar| {
            b.iter(|| black_box(parse_cigar_string(black_box(cigar))));
        });
    }

    group.finish();
}

/// Benchmark consensus base building
fn bench_consensus_base_builder(c: &mut Criterion) {
    let mut group = c.benchmark_group("consensus_base_builder");

    // Test with different numbers of observations
    for num_obs in [2_usize, 5, 10, 20, 50] {
        group.bench_with_input(
            BenchmarkId::new("build_consensus", num_obs),
            &num_obs,
            |b, &num_obs| {
                let mut builder = ConsensusBaseBuilder::new(45, 40);
                b.iter(|| {
                    builder.reset();
                    for i in 0..num_obs {
                        // Alternate bases with varying quality
                        let base = if i % 10 == 9 { b'C' } else { b'A' };
                        #[allow(clippy::cast_possible_truncation)]
                        let qual = 30 + (i % 10) as u8;
                        builder.add(base, qual);
                    }
                    black_box(builder.call())
                });
            },
        );
    }

    // Benchmark with unanimous agreement (common case)
    group.bench_function("unanimous_5_obs", |b| {
        let mut builder = ConsensusBaseBuilder::new(45, 40);
        b.iter(|| {
            builder.reset();
            for _ in 0..5 {
                builder.add(b'A', 35);
            }
            black_box(builder.call())
        });
    });

    // Benchmark with disagreement
    group.bench_function("disagreement_5_obs", |b| {
        let mut builder = ConsensusBaseBuilder::new(45, 40);
        b.iter(|| {
            builder.reset();
            builder.add(b'A', 35);
            builder.add(b'A', 35);
            builder.add(b'A', 35);
            builder.add(b'C', 35);
            builder.add(b'G', 35);
            black_box(builder.call())
        });
    });

    group.finish();
}

/// Benchmark UMI mismatch counting (`BitEnc` vs character comparison)
fn bench_umi_mismatch(c: &mut Criterion) {
    let mut group = c.benchmark_group("umi_mismatch");

    // Test various UMI lengths typical in sequencing
    let umi_lengths = [6, 8, 10, 12, 16, 18];

    for len in umi_lengths {
        // Generate two UMIs with ~10% mismatch rate
        let umi_a: String = (0..len).map(|i| if i % 10 == 0 { 'A' } else { 'C' }).collect();
        let umi_b: String = (0..len).map(|i| if i % 10 == 0 { 'T' } else { 'C' }).collect();

        group.bench_with_input(
            BenchmarkId::new("count_mismatches", len),
            &(umi_a.clone(), umi_b.clone()),
            |b, (a, b_str)| {
                b.iter(|| black_box(count_mismatches(black_box(a), black_box(b_str))));
            },
        );

        group.bench_with_input(
            BenchmarkId::new("matches_within_threshold", len),
            &(umi_a.clone(), umi_b.clone()),
            |b, (a, b_str)| {
                b.iter(|| black_box(matches_within_threshold(black_box(a), black_box(b_str), 2)));
            },
        );
    }

    // Benchmark BitEnc directly for comparison
    group.bench_function("bitenc_hamming_12bp", |b| {
        let enc_a = BitEnc::from_bytes(b"ACGTACGTACGT").unwrap();
        let enc_b = BitEnc::from_bytes(b"ACGTACGTACTT").unwrap();
        b.iter(|| black_box(enc_a.hamming_distance(black_box(&enc_b))));
    });

    group.bench_function("bitenc_hamming_18bp", |b| {
        let enc_a = BitEnc::from_bytes(b"ACGTACGTACGTACGTAC").unwrap();
        let enc_b = BitEnc::from_bytes(b"ACGTACGTACGTACGTAT").unwrap();
        b.iter(|| black_box(enc_a.hamming_distance(black_box(&enc_b))));
    });

    group.finish();
}

/// Benchmark many-to-many UMI comparison (simulating UMI clustering)
fn bench_umi_clustering(c: &mut Criterion) {
    let mut group = c.benchmark_group("umi_clustering");

    // Generate a set of UMIs for clustering simulation
    let num_umis = 100;
    let umi_len = 12;
    let umis: Vec<String> = (0..num_umis)
        .map(|i| {
            let bases = ['A', 'C', 'G', 'T'];
            (0..umi_len).map(|j| bases[(i + j) % 4]).collect()
        })
        .collect();

    group.bench_function("pairwise_100_umis_12bp", |b| {
        b.iter(|| {
            let mut count = 0usize;
            for i in 0..umis.len() {
                for j in (i + 1)..umis.len() {
                    if matches_within_threshold(&umis[i], &umis[j], 1) {
                        count += 1;
                    }
                }
            }
            black_box(count)
        });
    });

    group.finish();
}

/// Benchmark adjacency assigner with varying UMI counts to measure MIH improvement
fn bench_adjacency_assigner(c: &mut Criterion) {
    let mut group = c.benchmark_group("adjacency_assigner");
    group.sample_size(20); // Reduce sample size for longer benchmarks

    // Test with different numbers of unique UMIs
    // MIH kicks in at 100 UMIs, so test around that threshold
    for num_unique in [50, 100, 200, 500, 1000] {
        let umi_len = 12;
        let bases = ['A', 'C', 'G', 'T'];

        // Generate unique UMIs with some errors (to create clusters)
        let mut umis = Vec::new();
        for i in 0..num_unique {
            // Base UMI
            let base_umi: String = (0..umi_len).map(|j| bases[(i * 7 + j * 3) % 4]).collect();

            // Add multiple copies (simulating sequencing depth)
            let copies = 5 + (i % 10); // 5-14 copies each
            for _ in 0..copies {
                umis.push(base_umi.clone());
            }

            // Add error variants for ~20% of UMIs
            if i % 5 == 0 {
                let mut error_umi: Vec<char> = base_umi.chars().collect();
                error_umi[0] =
                    bases[(bases.iter().position(|&c| c == error_umi[0]).unwrap() + 1) % 4];
                let error_str: String = error_umi.into_iter().collect();
                umis.push(error_str);
            }
        }

        group.throughput(Throughput::Elements(num_unique as u64));
        group.bench_with_input(
            BenchmarkId::new("adjacency_1_edit", num_unique),
            &umis,
            |b, umis| {
                let assigner = Strategy::Adjacency.new_assigner(1);
                b.iter(|| {
                    let result = assigner.assign(black_box(&umis.clone()));
                    black_box(result)
                });
            },
        );
    }

    group.finish();
}

/// Create a test read for consensus calling benchmarks
fn create_test_read(name: &str, seq: &[u8], qual: u8, start_pos: usize, mi: &str) -> RecordBuf {
    use noodles::sam::alignment::record::cigar::Op;

    let mut rec = RecordBuf::default();
    *rec.name_mut() = Some(BString::from(name.as_bytes().to_vec()));
    *rec.sequence_mut() = Sequence::from(seq.to_vec());
    *rec.quality_scores_mut() = QualityScores::from(vec![qual; seq.len()]);
    *rec.flags_mut() = Flags::SEGMENTED | Flags::PROPERLY_SEGMENTED | Flags::FIRST_SEGMENT;

    if let Ok(pos) = Position::try_from(start_pos) {
        *rec.alignment_start_mut() = Some(pos);
    }

    let mut cigar = Cigar::default();
    cigar.as_mut().push(Op::new(Kind::Match, seq.len()));
    *rec.cigar_mut() = cigar;

    rec.data_mut().insert(Tag::from([b'M', b'I']), Value::from(mi.to_string()));
    rec.data_mut().insert(Tag::from([b'R', b'X']), Value::from("ACGTACGT".to_string()));
    rec
}

/// Benchmark `VanillaUmiConsensusCaller` (simplex consensus calling)
fn bench_vanilla_consensus_caller(c: &mut Criterion) {
    let mut group = c.benchmark_group("vanilla_consensus_caller");
    group.sample_size(50); // Reduce sample size for longer benchmarks

    let read_len = 150;
    let seq: Vec<u8> = (0..read_len).map(|i| b"ACGT"[i % 4]).collect();

    // Benchmark with different numbers of reads per molecule
    for num_reads in [2_usize, 3, 5, 10, 20] {
        // Create reads with slight variation
        let reads: Vec<RecordBuf> = (0..num_reads)
            .map(|i| {
                let mut read_seq = seq.clone();
                // Introduce small errors (~1% error rate)
                if i > 0 && i % 10 == 0 {
                    read_seq[i % read_len] = b"TGCA"[(read_seq[i % read_len] as usize) % 4];
                }
                create_test_read(&format!("read{i}"), &read_seq, 35, 1, "1")
            })
            .collect();

        let header = noodles::sam::Header::default();
        let raw_reads: Vec<Vec<u8>> = reads
            .iter()
            .map(|rec| {
                let mut buf = Vec::new();
                encode_record_buf(&mut buf, &header, rec).unwrap();
                buf
            })
            .collect();

        group.throughput(Throughput::Elements(num_reads as u64));
        group.bench_with_input(
            BenchmarkId::new("call_consensus", num_reads),
            &raw_reads,
            |b, raw_reads| {
                let options = VanillaUmiConsensusOptions::default();
                let mut caller =
                    VanillaUmiConsensusCaller::new("bench".to_string(), "RG1".to_string(), options);
                b.iter(|| {
                    let result = caller.consensus_reads(black_box(raw_reads.clone()));
                    black_box(result)
                });
            },
        );
    }

    // Benchmark with different read lengths
    for read_len in [100_usize, 150, 200, 300] {
        let seq: Vec<u8> = (0..read_len).map(|i| b"ACGT"[i % 4]).collect();
        let reads: Vec<RecordBuf> =
            (0..5).map(|i| create_test_read(&format!("read{i}"), &seq, 35, 1, "1")).collect();
        let header = noodles::sam::Header::default();
        let raw_reads: Vec<Vec<u8>> = reads
            .iter()
            .map(|rec| {
                let mut buf = Vec::new();
                encode_record_buf(&mut buf, &header, rec).unwrap();
                buf
            })
            .collect();

        group.throughput(Throughput::Bytes(read_len as u64 * 5));
        group.bench_with_input(
            BenchmarkId::new("read_length", read_len),
            &raw_reads,
            |b, raw_reads| {
                let options = VanillaUmiConsensusOptions::default();
                let mut caller =
                    VanillaUmiConsensusCaller::new("bench".to_string(), "RG1".to_string(), options);
                b.iter(|| {
                    let result = caller.consensus_reads(black_box(raw_reads.clone()));
                    black_box(result)
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_phred_conversions,
    bench_log_probability_ops,
    bench_dna_operations,
    bench_cigar_parsing,
    bench_consensus_base_builder,
    bench_umi_mismatch,
    bench_umi_clustering,
    bench_adjacency_assigner,
    bench_vanilla_consensus_caller,
);
criterion_main!(benches);
