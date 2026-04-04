use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use fgumi_simd_fastq::{SimdFastqReader, find_record_offsets, parse_records};
use seq_io::fastq::{Reader as SeqIoReader, Record};
use std::io::Cursor;

/// Generate synthetic FASTQ data with `n` records, each with a sequence of `seq_len` bases.
fn generate_fastq(n: usize, seq_len: usize) -> Vec<u8> {
    let mut data = Vec::with_capacity(n * (seq_len * 2 + 20));
    let seq: Vec<u8> = (0..seq_len).map(|i| b"ACGT"[i % 4]).collect();
    let qual: Vec<u8> = vec![b'I'; seq_len];

    for i in 0..n {
        data.push(b'@');
        data.extend_from_slice(format!("read{i}").as_bytes());
        data.push(b'\n');
        data.extend_from_slice(&seq);
        data.push(b'\n');
        data.push(b'+');
        data.push(b'\n');
        data.extend_from_slice(&qual);
        data.push(b'\n');
    }
    data
}

/// Previous boundary-finding implementation using memchr (4 calls per record).
fn find_record_offsets_memchr(data: &[u8]) -> Vec<usize> {
    if data.is_empty() {
        return vec![0];
    }

    let mut offsets = vec![0];
    let mut pos = 0;

    while pos < data.len() {
        if data[pos] != b'@' {
            break;
        }
        let mut found = true;
        for _ in 0..4 {
            if let Some(nl) = memchr::memchr(b'\n', &data[pos..]) {
                pos += nl + 1;
            } else {
                found = false;
                break;
            }
        }
        if found {
            offsets.push(pos);
        } else {
            break;
        }
    }

    offsets
}

fn bench_find_record_offsets(c: &mut Criterion) {
    let mut group = c.benchmark_group("find_record_offsets");

    for &seq_len in &[150, 300] {
        let num_records = 100_000;
        let data = generate_fastq(num_records, seq_len);
        let data_size = data.len();

        group.throughput(Throughput::Bytes(data_size as u64));

        group.bench_with_input(
            BenchmarkId::new("simd", format!("{seq_len}bp")),
            &data,
            |b, data| {
                b.iter(|| {
                    let offsets = find_record_offsets(data);
                    assert_eq!(offsets.len(), num_records + 1);
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new("memchr", format!("{seq_len}bp")),
            &data,
            |b, data| {
                b.iter(|| {
                    let offsets = find_record_offsets_memchr(data);
                    assert_eq!(offsets.len(), num_records + 1);
                });
            },
        );
    }

    // Benchmark with real FASTQ data if available
    let real_fastq_path = "/Volumes/scratch-00001/work/rebgzf_test.fastq";
    if std::path::Path::new(real_fastq_path).exists() {
        let data = std::fs::read(real_fastq_path).expect("failed to read real FASTQ file");
        let data_size = data.len();
        group.throughput(Throughput::Bytes(data_size as u64));

        group.bench_with_input(BenchmarkId::new("simd", "real_176MB"), &data, |b, data| {
            b.iter(|| find_record_offsets(data));
        });

        group.bench_with_input(BenchmarkId::new("memchr", "real_176MB"), &data, |b, data| {
            b.iter(|| find_record_offsets_memchr(data));
        });
    }

    group.finish();
}

/// Benchmark full record parsing: SIMD `parse_records` vs `seq_io` reader.
/// Both iterate all records and access name/seq/qual fields.
fn bench_record_parsing(c: &mut Criterion) {
    let mut group = c.benchmark_group("record_parsing");

    for &seq_len in &[150, 300] {
        let num_records = 100_000;
        let data = generate_fastq(num_records, seq_len);
        let data_size = data.len();

        group.throughput(Throughput::Bytes(data_size as u64));

        // SIMD: zero-copy parse from in-memory buffer
        group.bench_with_input(
            BenchmarkId::new("simd_parse_records", format!("{seq_len}bp")),
            &data,
            |b, data| {
                b.iter(|| {
                    let mut count = 0;
                    for rec in parse_records(data) {
                        std::hint::black_box(rec.name);
                        std::hint::black_box(rec.sequence);
                        std::hint::black_box(rec.quality);
                        count += 1;
                    }
                    assert_eq!(count, num_records);
                });
            },
        );

        // SIMD: SimdFastqReader from Cursor (simulates BufRead path)
        group.bench_with_input(
            BenchmarkId::new("simd_reader", format!("{seq_len}bp")),
            &data,
            |b, data| {
                b.iter(|| {
                    let reader = SimdFastqReader::new(Cursor::new(data));
                    let mut count = 0;
                    for result in reader {
                        let rec = result.expect("FASTQ record should parse successfully");
                        std::hint::black_box(&rec.name);
                        std::hint::black_box(&rec.sequence);
                        std::hint::black_box(&rec.quality);
                        count += 1;
                    }
                    assert_eq!(count, num_records);
                });
            },
        );

        // seq_io: the previous FASTQ reader
        group.bench_with_input(
            BenchmarkId::new("seq_io", format!("{seq_len}bp")),
            &data,
            |b, data| {
                b.iter(|| {
                    let reader = SeqIoReader::new(Cursor::new(data));
                    let mut count = 0;
                    for result in reader.into_records() {
                        let rec = result.expect("seq_io FASTQ record should parse successfully");
                        // Access fields to prevent dead-code elimination
                        std::hint::black_box(rec.head());
                        std::hint::black_box(rec.seq());
                        std::hint::black_box(rec.qual());
                        count += 1;
                    }
                    assert_eq!(count, num_records);
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_find_record_offsets, bench_record_parsing);
criterion_main!(benches);
