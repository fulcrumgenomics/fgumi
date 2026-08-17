//! Decode-throughput benchmark for the single-threaded raw-BAM reader unify
//! (#800): the fgumi-bgzf streaming decoder (`FgumiBgzfReader`) — with CRC32
//! verification on and off — versus noodles' single-threaded BGZF reader, all
//! decoding the *same* BGZF byte buffer. This isolates the decode-path change
//! the unify makes on the single-threaded fast paths.

use std::hint::black_box;
use std::io::{Cursor, Read};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use fgumi_bam_io::FgumiBgzfReader;
use fgumi_bgzf::InlineBgzfCompressor;

/// Build a moderately compressible payload of `len` bytes (BAM-like: a mix of a
/// small alphabet so blocks actually compress, deterministic so runs compare).
fn make_payload(len: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(len);
    let mut state: u32 = 0x1234_5678;
    while out.len() < len {
        // xorshift for cheap deterministic pseudo-randomness
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        // Fold into a small-ish alphabet so the deflate stage has real work.
        out.push((state & 0x3F) as u8 + b' ');
    }
    out.truncate(len);
    out
}

/// Compress `payload` into a complete BGZF byte stream (blocks + EOF marker).
fn make_bgzf(payload: &[u8], level: u32) -> Vec<u8> {
    let mut compressor = InlineBgzfCompressor::new(level);
    compressor.write_all(payload).expect("compress payload");
    compressor.flush().expect("flush compressor");
    let mut out = Vec::new();
    compressor.write_blocks_to(&mut out).expect("write blocks");
    out.extend_from_slice(&fgumi_bgzf::BGZF_EOF);
    out
}

fn bench_bgzf_decode(c: &mut Criterion) {
    let payload = make_payload(32 * 1024 * 1024); // 32 MiB uncompressed
    let bgzf = make_bgzf(&payload, 6);

    let mut group = c.benchmark_group("bgzf_decode");
    group.throughput(Throughput::Bytes(payload.len() as u64));

    let mut scratch = Vec::with_capacity(payload.len());

    group.bench_function("fgumi_verify_on", |b| {
        b.iter(|| {
            scratch.clear();
            let mut reader = FgumiBgzfReader::new(Box::new(Cursor::new(bgzf.clone())), true);
            reader.read_to_end(&mut scratch).expect("decode");
            black_box(scratch.len());
        });
    });

    group.bench_function("fgumi_verify_off", |b| {
        b.iter(|| {
            scratch.clear();
            let mut reader = FgumiBgzfReader::new(Box::new(Cursor::new(bgzf.clone())), false);
            reader.read_to_end(&mut scratch).expect("decode");
            black_box(scratch.len());
        });
    });

    group.bench_function("noodles", |b| {
        b.iter(|| {
            scratch.clear();
            let mut reader = noodles_bgzf::io::Reader::new(Cursor::new(bgzf.clone()));
            reader.read_to_end(&mut scratch).expect("decode");
            black_box(scratch.len());
        });
    });

    group.finish();
}

criterion_group!(benches, bench_bgzf_decode);
criterion_main!(benches);
