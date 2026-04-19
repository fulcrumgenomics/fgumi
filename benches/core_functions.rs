//! Benchmarks for core fgumi functions.
//!
//! Run with: `cargo bench`
//! View reports in: `target/criterion/report/index.html`
#![allow(
    clippy::cast_possible_truncation,
    clippy::doc_markdown,
    clippy::needless_range_loop,
    clippy::too_many_lines
)]

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use std::cmp::Ordering;
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
use fgumi_raw_bam::RawRecord;
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
        let raw_reads: Vec<RawRecord> = reads
            .iter()
            .map(|rec| {
                let mut buf = Vec::new();
                encode_record_buf(&mut buf, &header, rec).unwrap();
                RawRecord::from(buf)
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
        let raw_reads: Vec<RawRecord> = reads
            .iter()
            .map(|rec| {
                let mut buf = Vec::new();
                encode_record_buf(&mut buf, &header, rec).unwrap();
                RawRecord::from(buf)
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

// ============================================================================
// Queryname comparison benchmarks
// ============================================================================

/// Compare using only the first 8 bytes as a u64 (O(1), wrong results, upper-bound test).
#[inline]
fn compare_u64_hash(a: &[u8], b: &[u8]) -> Ordering {
    let a_val = if a.len() >= 8 {
        u64::from_le_bytes(a[..8].try_into().unwrap())
    } else {
        let mut buf = [0u8; 8];
        buf[..a.len()].copy_from_slice(a);
        u64::from_le_bytes(buf)
    };
    let b_val = if b.len() >= 8 {
        u64::from_le_bytes(b[..8].try_into().unwrap())
    } else {
        let mut buf = [0u8; 8];
        buf[..b.len()].copy_from_slice(b);
        u64::from_le_bytes(buf)
    };
    a_val.cmp(&b_val)
}

/// Compare with common prefix stripped (caller pre-computes `prefix_len`).
#[inline]
fn natural_compare_prefix_stripped(a: &[u8], b: &[u8], prefix_len: usize) -> Ordering {
    fgumi_lib::sort::natural_compare(&a[prefix_len..], &b[prefix_len..])
}

/// Compute digit-safe common prefix length: back up to the last non-digit byte
/// before the divergence point, so we never split inside a numeric run.
///
/// Example: "x123y" vs "x12ay" diverge at index 3, but index 2 is a digit,
/// so we back up to index 1 (the 'x'). Stripping "x" is safe; stripping "x12" is not.
fn digit_safe_prefix_len(names: &[Vec<u8>]) -> usize {
    if names.len() <= 1 {
        return 0;
    }
    let first = &names[0];
    let mut common = first.len();
    for name in &names[1..] {
        let shared = first.iter().zip(name.iter()).take_while(|(a, b)| a == b).count();
        common = common.min(shared);
        if common == 0 {
            break;
        }
    }
    // Back up past any trailing digits so we don't split a numeric run.
    while common > 0 && first[common - 1].is_ascii_digit() {
        common -= 1;
    }
    common
}

/// Sort with digit-safe common prefix stripping + natural_compare.
fn sort_natural_digit_safe_prefix(indices: &mut [usize], raw: &[Vec<u8>], safe_prefix: usize) {
    indices.sort_by(|&a, &b| {
        fgumi_lib::sort::natural_compare(&raw[a][safe_prefix..], &raw[b][safe_prefix..])
    });
}

/// Faithful port of samtools' `strnum_cmp` using unsafe pointer walking
/// on null-terminated strings — no bounds checks in the hot loop.
///
/// Caller must ensure `a` and `b` are null-terminated.
#[inline]
unsafe fn strnum_cmp_samtools(a: *const u8, b: *const u8) -> i32 {
    unsafe {
        let mut pa = a;
        let mut pb = b;
        while *pa != 0 && *pb != 0 {
            if !(*pa).is_ascii_digit() || !(*pb).is_ascii_digit() {
                if *pa != *pb {
                    return i32::from(*pa) - i32::from(*pb);
                }
                pa = pa.add(1);
                pb = pb.add(1);
            } else {
                while *pa == b'0' {
                    pa = pa.add(1);
                }
                while *pb == b'0' {
                    pb = pb.add(1);
                }
                while (*pa).is_ascii_digit() && *pa == *pb {
                    pa = pa.add(1);
                    pb = pb.add(1);
                }
                let diff = i32::from(*pa) - i32::from(*pb);
                while (*pa).is_ascii_digit() && (*pb).is_ascii_digit() {
                    pa = pa.add(1);
                    pb = pb.add(1);
                }
                if (*pa).is_ascii_digit() {
                    return 1;
                } else if (*pb).is_ascii_digit() {
                    return -1;
                } else if diff != 0 {
                    return diff;
                }
            }
        }
        if *pa != 0 {
            1
        } else if *pb != 0 {
            -1
        } else {
            0
        }
    }
}

/// Wrapper using pre-allocated null-terminated names (avoids alloc in the hot path).
#[inline]
fn compare_strnum_samtools_prealloc(a: &[u8], b: &[u8]) -> Ordering {
    // Assumes a and b are already null-terminated.
    let r = unsafe { strnum_cmp_samtools(a.as_ptr(), b.as_ptr()) };
    r.cmp(&0)
}

/// Wrapper for production `natural_compare_nul` using pre-allocated null-terminated names.
#[inline]
fn compare_natural_nul_prealloc(a: &[u8], b: &[u8]) -> Ordering {
    // Assumes a and b are already null-terminated.
    unsafe { fgumi_lib::sort::keys::natural_compare_nul(a.as_ptr(), b.as_ptr()) }
}

/// Parse colon-delimited Illumina fields and compare as integers.
/// Format: INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y
/// Compares fields 0-3 as bytes (usually equal), then tile/x/y as integers.
#[inline]
fn compare_illumina_structured(a: &[u8], b: &[u8]) -> Ordering {
    let mut a_fields = a.splitn(7, |&c| c == b':');
    let mut b_fields = b.splitn(7, |&c| c == b':');

    // Compare first 4 fields as bytes (instrument:run:flowcell:lane)
    for _ in 0..4 {
        let af = a_fields.next().unwrap_or(b"");
        let bf = b_fields.next().unwrap_or(b"");
        match af.cmp(bf) {
            Ordering::Equal => {}
            ord => return ord,
        }
    }

    // Compare tile, x, y as integers
    for _ in 0..3 {
        let af = a_fields.next().unwrap_or(b"");
        let bf = b_fields.next().unwrap_or(b"");
        let a_num = parse_int_fast(af);
        let b_num = parse_int_fast(bf);
        match a_num.cmp(&b_num) {
            Ordering::Equal => {}
            ord => return ord,
        }
    }

    Ordering::Equal
}

#[inline]
fn parse_int_fast(bytes: &[u8]) -> u64 {
    let mut val: u64 = 0;
    for &b in bytes {
        if b.is_ascii_digit() {
            val = val * 10 + u64::from(b - b'0');
        } else {
            break;
        }
    }
    val
}

/// Generate realistic sorted name pairs for benchmarking.
fn generate_name_pairs(style: &str, count: usize) -> Vec<(Vec<u8>, Vec<u8>)> {
    let mut pairs = Vec::with_capacity(count);
    match style {
        "illumina" => {
            // Illumina-style: A00132:53:HFHJKDSXX:LANE:TILE:X:Y
            // 3 different instrument/flowcell combos (like wgs_chr1)
            let prefixes = ["A00132:53:HFHJKDSXX", "A00132:54:HFH2JDSXX", "A00132:55:HFHKKDSXX"];
            for i in 0..count {
                let p1 = prefixes[i % 3];
                let p2 = prefixes[(i + 1) % 3];
                let lane1 = (i % 4) + 1;
                let lane2 = ((i + 1) % 4) + 1;
                let tile1 = 1100 + (i % 500);
                let tile2 = 1100 + ((i + 7) % 500);
                let x1 = 1000 + (i * 17 % 30000);
                let x2 = 1000 + ((i + 3) * 17 % 30000);
                let y1 = 1000 + (i * 31 % 50000);
                let y2 = 1000 + ((i + 5) * 31 % 50000);
                let a = format!("{p1}:{lane1}:{tile1}:{x1}:{y1}");
                let b = format!("{p2}:{lane2}:{tile2}:{x2}:{y2}");
                pairs.push((a.into_bytes(), b.into_bytes()));
            }
        }
        "illumina_same_prefix" => {
            // Same instrument/flowcell — has common prefix
            for i in 0..count {
                let lane1 = (i % 4) + 1;
                let lane2 = ((i + 1) % 4) + 1;
                let tile1 = 1100 + (i % 500);
                let tile2 = 1100 + ((i + 7) % 500);
                let x1 = 1000 + (i * 17 % 30000);
                let x2 = 1000 + ((i + 3) * 17 % 30000);
                let y1 = 1000 + (i * 31 % 50000);
                let y2 = 1000 + ((i + 5) * 31 % 50000);
                let a = format!("A00132:53:HFHJKDSXX:{lane1}:{tile1}:{x1}:{y1}");
                let b = format!("A00132:53:HFHJKDSXX:{lane2}:{tile2}:{x2}:{y2}");
                pairs.push((a.into_bytes(), b.into_bytes()));
            }
        }
        "srr" => {
            // SRR-style: SRR6109273.NNNNN
            for i in 0..count {
                let a = format!("SRR6109273.{}", 100_000 + i * 7);
                let b = format!("SRR6109273.{}", 100_000 + (i + 1) * 7);
                pairs.push((a.into_bytes(), b.into_bytes()));
            }
        }
        _ => unreachable!(),
    }
    pairs
}

/// Compare using u64 prefix from normalized key + `natural_compare` fallback.
#[inline]
fn compare_normalized_prefix_with_fallback(
    a: &[u8],
    b: &[u8],
    a_norm_prefix: u64,
    b_norm_prefix: u64,
) -> Ordering {
    a_norm_prefix.cmp(&b_norm_prefix).then_with(|| fgumi_lib::sort::natural_compare(a, b))
}

/// Pre-compute normalized key prefixes for benchmark pairs.
fn precompute_normalized_prefixes(pairs: &[(Vec<u8>, Vec<u8>)]) -> Vec<(u64, u64)> {
    pairs
        .iter()
        .map(|(a, b)| {
            let mut norm_a = Vec::with_capacity(a.len() + 8);
            let mut norm_b = Vec::with_capacity(b.len() + 8);
            fgumi_lib::sort::normalize_natural_key(a, &mut norm_a);
            fgumi_lib::sort::normalize_natural_key(b, &mut norm_b);
            let a_prefix = {
                let mut buf = [0u8; 8];
                let n = norm_a.len().min(8);
                buf[..n].copy_from_slice(&norm_a[..n]);
                u64::from_be_bytes(buf)
            };
            let b_prefix = {
                let mut buf = [0u8; 8];
                let n = norm_b.len().min(8);
                buf[..n].copy_from_slice(&norm_b[..n]);
                u64::from_be_bytes(buf)
            };
            (a_prefix, b_prefix)
        })
        .collect()
}

fn bench_queryname_comparators(c: &mut Criterion) {
    let mut group = c.benchmark_group("queryname_compare");
    let n = 10_000;

    for style in ["illumina", "illumina_same_prefix", "srr"] {
        let pairs = generate_name_pairs(style, n);

        // Compute common prefix length for this dataset
        let prefix_len = if let Some((first, _)) = pairs.first() {
            pairs.iter().fold(first.len(), |pl, (a, _)| {
                let common = a.iter().zip(first.iter()).take_while(|(x, y)| x == y).count();
                pl.min(common)
            })
        } else {
            0
        };

        // Pre-compute normalized prefixes for the new benchmark
        let norm_prefixes = precompute_normalized_prefixes(&pairs);

        group.throughput(Throughput::Elements(n as u64));

        group.bench_with_input(BenchmarkId::new("natural_compare", style), &pairs, |b, pairs| {
            b.iter(|| {
                for (a, b_name) in pairs {
                    black_box(fgumi_lib::sort::natural_compare(a, b_name));
                }
            });
        });

        group.bench_with_input(BenchmarkId::new("u64_hash_O1", style), &pairs, |b, pairs| {
            b.iter(|| {
                for (a, b_name) in pairs {
                    black_box(compare_u64_hash(a, b_name));
                }
            });
        });

        // Normalized prefix + natural_compare fallback (the new approach)
        group.bench_with_input(
            BenchmarkId::new("normalized_prefix", style),
            &(pairs.clone(), norm_prefixes.clone()),
            |b, (pairs, prefixes)| {
                b.iter(|| {
                    for ((a, b_name), (ap, bp)) in pairs.iter().zip(prefixes.iter()) {
                        black_box(compare_normalized_prefix_with_fallback(a, b_name, *ap, *bp));
                    }
                });
            },
        );

        // Benchmark normalize_natural_key encoding throughput
        group.bench_with_input(
            BenchmarkId::new("normalize_key_encode", style),
            &pairs,
            |b, pairs| {
                let mut buf = Vec::with_capacity(64);
                b.iter(|| {
                    for (a, _) in pairs {
                        buf.clear();
                        fgumi_lib::sort::normalize_natural_key(black_box(a), &mut buf);
                        black_box(&buf);
                    }
                });
            },
        );

        if prefix_len > 0 {
            group.bench_with_input(
                BenchmarkId::new("prefix_stripped", style),
                &pairs,
                |b, pairs| {
                    b.iter(|| {
                        for (a, b_name) in pairs {
                            black_box(natural_compare_prefix_stripped(a, b_name, prefix_len));
                        }
                    });
                },
            );
        }

        if style.starts_with("illumina") {
            group.bench_with_input(
                BenchmarkId::new("structured_fields", style),
                &pairs,
                |b, pairs| {
                    b.iter(|| {
                        for (a, b_name) in pairs {
                            black_box(compare_illumina_structured(a, b_name));
                        }
                    });
                },
            );
        }

        group.bench_with_input(BenchmarkId::new("byte_cmp_O1", style), &pairs, |b, pairs| {
            b.iter(|| {
                for (a, b_name) in pairs {
                    black_box(a.cmp(b_name));
                }
            });
        });

        // samtools strnum_cmp (natural) — pre-allocate null-terminated names
        let nul_pairs: Vec<(Vec<u8>, Vec<u8>)> = pairs
            .iter()
            .map(|(a, b)| {
                let mut an = Vec::with_capacity(a.len() + 1);
                an.extend_from_slice(a);
                an.push(0);
                let mut bn = Vec::with_capacity(b.len() + 1);
                bn.extend_from_slice(b);
                bn.push(0);
                (an, bn)
            })
            .collect();

        group.bench_with_input(
            BenchmarkId::new("strnum_cmp_samtools", style),
            &nul_pairs,
            |b, nul_pairs| {
                b.iter(|| {
                    for (a, b_name) in nul_pairs {
                        black_box(compare_strnum_samtools_prealloc(a, b_name));
                    }
                });
            },
        );

        // Production natural_compare_nul (null-terminated, no bounds checks)
        group.bench_with_input(
            BenchmarkId::new("natural_compare_nul", style),
            &nul_pairs,
            |b, nul_pairs| {
                b.iter(|| {
                    for (a, b_name) in nul_pairs {
                        black_box(compare_natural_nul_prealloc(a, b_name));
                    }
                });
            },
        );
    }

    group.finish();
}

// ============================================================================
// Queryname sort-strategy benchmarks
// ============================================================================

/// Build a big-endian u64 prefix from the first 8 bytes of a name.
#[inline]
fn name_prefix_u64(name: &[u8]) -> u64 {
    let mut buf = [0u8; 8];
    let n = name.len().min(8);
    buf[..n].copy_from_slice(&name[..n]);
    u64::from_be_bytes(buf)
}

/// Generate a realistic set of querynames for sort benchmarks.
///
/// Returns `(names, u64_prefixes)`.  The names come from a single Illumina
/// flowcell so the first ~20 bytes are identical — the challenging case for
/// prefix-based strategies.
fn generate_sort_names(count: usize) -> (Vec<Vec<u8>>, Vec<u64>) {
    let mut names = Vec::with_capacity(count);
    // Single-run Illumina: common prefix, varying lane:tile:x:y
    for i in 0..count {
        let lane = (i % 4) + 1;
        let tile = 1100 + (i * 7 % 2400);
        let x = 1000 + (i * 17 % 30000);
        let y = 1000 + (i * 31 % 50000);
        names.push(format!("A00558:940:HFHJKDSXX:{lane}:{tile}:{x}:{y}").into_bytes());
    }
    // Shuffle deterministically so the sort actually has work to do.
    // Simple Fisher-Yates with fixed seed.
    let mut seed: u64 = 0xDEAD_BEEF_CAFE_BABE;
    for i in (1..names.len()).rev() {
        // xorshift64
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        let j = (seed as usize) % (i + 1);
        names.swap(i, j);
    }
    let prefixes: Vec<u64> = names.iter().map(|n| name_prefix_u64(n)).collect();
    (names, prefixes)
}

/// (a) Current strategy: comparison sort with u64 prefix fast-reject + byte fallback.
fn sort_comparison_prefix(names: &mut [(u64, usize)], raw: &[Vec<u8>]) {
    names.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| raw[a.1].cmp(&raw[b.1])));
}

/// (b) Radix sort on u64 prefix, then comparison-sort within equal-prefix buckets.
fn sort_radix_prefix_then_comparison(entries: &mut [(u64, usize)], raw: &[Vec<u8>]) {
    // LSD radix sort on the u64 prefix (8 passes, 1 byte each).
    let len = entries.len();
    let mut aux = vec![(0u64, 0usize); len];
    for pass in 0..8u32 {
        let shift = pass * 8;
        let mut counts = [0u32; 256];
        for e in entries.iter() {
            let bucket = ((e.0 >> shift) & 0xFF) as usize;
            counts[bucket] += 1;
        }
        let mut offsets = [0u32; 256];
        for i in 1..256 {
            offsets[i] = offsets[i - 1] + counts[i - 1];
        }
        for e in entries.iter() {
            let bucket = ((e.0 >> shift) & 0xFF) as usize;
            aux[offsets[bucket] as usize] = *e;
            offsets[bucket] += 1;
        }
        entries.copy_from_slice(&aux);
    }
    // Now entries are sorted by prefix.  Find equal-prefix runs and sort them.
    let mut start = 0;
    while start < len {
        let prefix = entries[start].0;
        let mut end = start + 1;
        while end < len && entries[end].0 == prefix {
            end += 1;
        }
        if end - start > 1 {
            entries[start..end].sort_by(|a, b| raw[a.1].cmp(&raw[b.1]));
        }
        start = end;
    }
}

/// (c) Parallel comparison sort (rayon `par_sort_by`).
fn sort_parallel_comparison(names: &mut [(u64, usize)], raw: &[Vec<u8>]) {
    use rayon::prelude::*;
    names.par_sort_by(|a, b| a.0.cmp(&b.0).then_with(|| raw[a.1].cmp(&raw[b.1])));
}

/// (d) Full MSD radix sort on name bytes with common-prefix skip.
///
/// Finds the longest common prefix across all names, then starts the
/// radix sort at the first divergent byte — avoiding wasted passes over
/// identical leading characters (e.g. the ~20-byte Illumina prefix).
fn sort_msd_radix_u64(indices: &mut [usize], raw: &[Vec<u8>]) {
    if indices.len() <= 1 {
        return;
    }
    // Find common prefix length across all names.
    let first = &raw[indices[0]];
    let mut common = first.len();
    for &idx in &indices[1..] {
        let name = &raw[idx];
        let shared = first.iter().zip(name.iter()).take_while(|(a, b)| a == b).count();
        common = common.min(shared);
        if common == 0 {
            break;
        }
    }
    msd_radix_recurse(indices, raw, common);
}

/// Number of buckets for MSD radix: 65536 (u16) + 1 for end-of-string.
const MSD_BUCKETS: usize = 65537;
/// Bytes consumed per radix pass.
const MSD_DIGIT_BYTES: usize = 2;

/// Extract a u16 digit from `name` at byte offset `depth`.
/// Names shorter than `depth + 2` are zero-padded; names ended at `depth`
/// get bucket 0 (sorts before all byte values).
#[inline]
fn msd_digit(name: &[u8], depth: usize) -> usize {
    let remaining = name.len().saturating_sub(depth);
    if remaining == 0 {
        0 // end-of-string sentinel
    } else if remaining == 1 {
        // One byte left: pack as high byte of u16, +1 to leave room for sentinel.
        (usize::from(name[depth]) << 8) + 1
    } else {
        let hi = usize::from(name[depth]);
        let lo = usize::from(name[depth + 1]);
        (hi << 8 | lo) + 1
    }
}

fn msd_radix_recurse(indices: &mut [usize], raw: &[Vec<u8>], depth: usize) {
    if indices.len() <= 32 {
        indices.sort_by(|&a, &b| raw[a][depth..].cmp(&raw[b][depth..]));
        return;
    }

    // Count occurrences per u16 bucket.
    let mut counts = vec![0u32; MSD_BUCKETS];
    for &idx in indices.iter() {
        counts[msd_digit(&raw[idx], depth)] += 1;
    }

    // Build offsets.
    let mut offsets = vec![0u32; MSD_BUCKETS];
    for i in 1..MSD_BUCKETS {
        offsets[i] = offsets[i - 1] + counts[i - 1];
    }

    // Scatter into aux buffer.
    let mut aux = vec![0usize; indices.len()];
    let mut write_pos = offsets.clone();
    for &idx in indices.iter() {
        let bucket = msd_digit(&raw[idx], depth);
        aux[write_pos[bucket] as usize] = idx;
        write_pos[bucket] += 1;
    }
    indices.copy_from_slice(&aux);

    // Recurse into non-trivial buckets (skip bucket 0 = end-of-string).
    let mut start = 0;
    for b in 0..MSD_BUCKETS {
        let end = start + counts[b] as usize;
        if end - start > 1 && b > 0 {
            msd_radix_recurse(&mut indices[start..end], raw, depth + MSD_DIGIT_BYTES);
        }
        start = end;
    }
}

/// (f) Natural sort with normalized prefix (current production approach).
///
/// Pre-computes u64 prefix from `normalize_natural_key`, uses it for fast-reject,
/// falls back to `natural_compare` when prefixes match.
fn sort_natural_with_norm_prefix(entries: &mut [(u64, usize)], raw: &[Vec<u8>]) {
    entries.sort_by(|a, b| {
        a.0.cmp(&b.0).then_with(|| fgumi_lib::sort::natural_compare(&raw[a.1], &raw[b.1]))
    });
}

/// Pre-compute normalized prefixes for sort benchmarks.
fn precompute_norm_prefixes_for_sort(names: &[Vec<u8>]) -> Vec<u64> {
    names
        .iter()
        .map(|n| {
            let mut buf = Vec::with_capacity(n.len() + 8);
            fgumi_lib::sort::normalize_natural_key(n, &mut buf);
            name_prefix_u64(&buf)
        })
        .collect()
}

/// (g) Natural sort, no prefix — pure `natural_compare` on every comparison.
/// This is what samtools effectively does.
fn sort_natural_no_prefix(indices: &mut [usize], raw: &[Vec<u8>]) {
    indices.sort_by(|&a, &b| fgumi_lib::sort::natural_compare(&raw[a], &raw[b]));
}

/// (j) Natural sort with null-terminated pointer-walking comparator.
fn sort_natural_nul(indices: &mut [usize], raw_nul: &[Vec<u8>]) {
    indices.sort_by(|&a, &b| unsafe {
        fgumi_lib::sort::keys::natural_compare_nul(raw_nul[a].as_ptr(), raw_nul[b].as_ptr())
    });
}

/// (k) Natural sort with null-terminated comparator + digit-safe prefix stripping.
fn sort_natural_nul_digit_safe(indices: &mut [usize], raw_nul: &[Vec<u8>], safe_prefix: usize) {
    indices.sort_by(|&a, &b| unsafe {
        fgumi_lib::sort::keys::natural_compare_nul(
            raw_nul[a].as_ptr().add(safe_prefix),
            raw_nul[b].as_ptr().add(safe_prefix),
        )
    });
}

/// (h) Natural sort with allocation-free prefix.
///
/// Computes the u64 prefix by inlining the first 8 bytes of the normalized
/// encoding directly — no Vec allocation needed.
#[inline]
fn normalize_prefix_inline(name: &[u8]) -> u64 {
    let mut buf = [0u8; 8];
    let mut out_pos = 0;
    let mut i = 0;
    while i < name.len() && out_pos < 8 {
        if name[i].is_ascii_digit() {
            // Skip leading zeros
            while i < name.len() && name[i] == b'0' {
                i += 1;
            }
            // Count significant digits
            let sig_start = i;
            while i < name.len() && name[i].is_ascii_digit() {
                i += 1;
            }
            let sig_len = i - sig_start;
            if sig_len == 0 {
                // Was all zeros
                if out_pos < 8 {
                    buf[out_pos] = 1;
                    out_pos += 1;
                }
                if out_pos < 8 {
                    buf[out_pos] = b'0';
                    out_pos += 1;
                }
            } else {
                let len_byte = sig_len.min(255) as u8;
                if out_pos < 8 {
                    buf[out_pos] = len_byte;
                    out_pos += 1;
                }
                let copy_end = sig_start + sig_len.min(255);
                let mut j = sig_start;
                while j < copy_end && out_pos < 8 {
                    buf[out_pos] = name[j];
                    out_pos += 1;
                    j += 1;
                }
            }
        } else {
            buf[out_pos] = name[i];
            out_pos += 1;
            i += 1;
        }
    }
    u64::from_be_bytes(buf)
}

/// (e) Comparison sort with common-prefix skip — simplest possible optimization.
///
/// Find the shared prefix, then compare only the suffix.
fn sort_comparison_prefix_stripped(indices: &mut [usize], raw: &[Vec<u8>]) {
    if indices.len() <= 1 {
        return;
    }
    let first = &raw[indices[0]];
    let mut common = first.len();
    for &idx in &indices[1..] {
        let name = &raw[idx];
        let shared = first.iter().zip(name.iter()).take_while(|(a, b)| a == b).count();
        common = common.min(shared);
        if common == 0 {
            break;
        }
    }
    indices.sort_by(|&a, &b| raw[a][common..].cmp(&raw[b][common..]));
}

fn bench_queryname_sort_strategies(c: &mut Criterion) {
    let mut group = c.benchmark_group("queryname_sort_strategies");
    group.sample_size(10); // Sorts are expensive
    group.measurement_time(std::time::Duration::from_secs(10));

    for count in [100_000, 500_000] {
        let (names, prefixes) = generate_sort_names(count);
        group.throughput(Throughput::Elements(count as u64));

        // (a) Current: comparison sort with u64 prefix fast-reject
        group.bench_with_input(
            BenchmarkId::new("comparison_prefix", count),
            &(&names, &prefixes),
            |b, &(names, prefixes)| {
                b.iter_batched(
                    || prefixes.iter().enumerate().map(|(i, &p)| (p, i)).collect::<Vec<_>>(),
                    |mut entries| {
                        sort_comparison_prefix(&mut entries, names);
                        black_box(entries)
                    },
                    criterion::BatchSize::LargeInput,
                );
            },
        );

        // (b) Radix on u64 prefix, then comparison within equal-prefix buckets
        group.bench_with_input(
            BenchmarkId::new("radix_prefix_then_cmp", count),
            &(&names, &prefixes),
            |b, &(names, prefixes)| {
                b.iter_batched(
                    || prefixes.iter().enumerate().map(|(i, &p)| (p, i)).collect::<Vec<_>>(),
                    |mut entries| {
                        sort_radix_prefix_then_comparison(&mut entries, names);
                        black_box(entries)
                    },
                    criterion::BatchSize::LargeInput,
                );
            },
        );

        // (c) Parallel comparison sort (rayon)
        group.bench_with_input(
            BenchmarkId::new("parallel_comparison", count),
            &(&names, &prefixes),
            |b, &(names, prefixes)| {
                b.iter_batched(
                    || prefixes.iter().enumerate().map(|(i, &p)| (p, i)).collect::<Vec<_>>(),
                    |mut entries| {
                        sort_parallel_comparison(&mut entries, names);
                        black_box(entries)
                    },
                    criterion::BatchSize::LargeInput,
                );
            },
        );

        // (d) MSD radix sort with common-prefix skip
        group.bench_with_input(BenchmarkId::new("msd_radix_full", count), &names, |b, names| {
            b.iter_batched(
                || (0..names.len()).collect::<Vec<_>>(),
                |mut indices| {
                    sort_msd_radix_u64(&mut indices, names);
                    black_box(indices)
                },
                criterion::BatchSize::LargeInput,
            );
        });

        // (e) Comparison sort with common-prefix skip
        group.bench_with_input(
            BenchmarkId::new("cmp_prefix_stripped", count),
            &names,
            |b, names| {
                b.iter_batched(
                    || (0..names.len()).collect::<Vec<_>>(),
                    |mut indices| {
                        sort_comparison_prefix_stripped(&mut indices, names);
                        black_box(indices)
                    },
                    criterion::BatchSize::LargeInput,
                );
            },
        );

        // --- Natural sort variants ---

        // Pre-compute normalized prefixes (with Vec allocation, as production code does)
        let norm_prefixes = precompute_norm_prefixes_for_sort(&names);

        // Pre-compute inline prefixes (no allocation)
        let inline_prefixes: Vec<u64> = names.iter().map(|n| normalize_prefix_inline(n)).collect();

        // (f) Natural sort with normalized prefix (current production approach)
        group.bench_with_input(
            BenchmarkId::new("natural_norm_prefix", count),
            &(&names, &norm_prefixes),
            |b, &(names, prefixes)| {
                b.iter_batched(
                    || prefixes.iter().enumerate().map(|(i, &p)| (p, i)).collect::<Vec<_>>(),
                    |mut entries| {
                        sort_natural_with_norm_prefix(&mut entries, names);
                        black_box(entries)
                    },
                    criterion::BatchSize::LargeInput,
                );
            },
        );

        // (g) Natural sort, no prefix (samtools approach)
        group.bench_with_input(BenchmarkId::new("natural_no_prefix", count), &names, |b, names| {
            b.iter_batched(
                || (0..names.len()).collect::<Vec<_>>(),
                |mut indices| {
                    sort_natural_no_prefix(&mut indices, names);
                    black_box(indices)
                },
                criterion::BatchSize::LargeInput,
            );
        });

        // (h) Natural sort with allocation-free inline prefix
        group.bench_with_input(
            BenchmarkId::new("natural_inline_prefix", count),
            &(&names, &inline_prefixes),
            |b, &(names, prefixes)| {
                b.iter_batched(
                    || prefixes.iter().enumerate().map(|(i, &p)| (p, i)).collect::<Vec<_>>(),
                    |mut entries| {
                        sort_natural_with_norm_prefix(&mut entries, names);
                        black_box(entries)
                    },
                    criterion::BatchSize::LargeInput,
                );
            },
        );

        // (i) Natural sort with digit-safe common prefix stripping
        let safe_prefix = digit_safe_prefix_len(&names);
        group.bench_with_input(
            BenchmarkId::new("natural_digit_safe_prefix", count),
            &(&names, safe_prefix),
            |b, &(names, safe_prefix)| {
                b.iter_batched(
                    || (0..names.len()).collect::<Vec<_>>(),
                    |mut indices| {
                        sort_natural_digit_safe_prefix(&mut indices, names, safe_prefix);
                        black_box(indices)
                    },
                    criterion::BatchSize::LargeInput,
                );
            },
        );

        // Build null-terminated names for pointer-walking benchmarks
        let names_nul: Vec<Vec<u8>> = names
            .iter()
            .map(|n| {
                let mut v = Vec::with_capacity(n.len() + 1);
                v.extend_from_slice(n);
                v.push(0);
                v
            })
            .collect();

        // (j) Natural sort with null-terminated pointer-walking comparator
        group.bench_with_input(
            BenchmarkId::new("natural_nul", count),
            &names_nul,
            |b, names_nul| {
                b.iter_batched(
                    || (0..names_nul.len()).collect::<Vec<_>>(),
                    |mut indices| {
                        sort_natural_nul(&mut indices, names_nul);
                        black_box(indices)
                    },
                    criterion::BatchSize::LargeInput,
                );
            },
        );

        // (k) Natural sort with null-terminated comparator + digit-safe prefix skip
        group.bench_with_input(
            BenchmarkId::new("natural_nul_digit_safe", count),
            &(&names_nul, safe_prefix),
            |b, &(names_nul, safe_prefix)| {
                b.iter_batched(
                    || (0..names_nul.len()).collect::<Vec<_>>(),
                    |mut indices| {
                        sort_natural_nul_digit_safe(&mut indices, names_nul, safe_prefix);
                        black_box(indices)
                    },
                    criterion::BatchSize::LargeInput,
                );
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
    bench_queryname_comparators,
    bench_queryname_sort_strategies,
);
criterion_main!(benches);
