//! Decode throughput for stored (level-0) vs. deflate BGZF blocks.
//!
//! The stored-block path in `decompress_and_verify` skips libdeflater and
//! memcpy's the payload directly. This bench has three rows:
//!
//! * `stored_bypass` — stored-block input, current code (bypass active).
//!   This is the configuration we ship.
//! * `stored_via_libdeflater` — stored-block input, but we call libdeflater
//!   directly to mimic the pre-bypass behaviour. The delta between this row
//!   and `stored_bypass` is the bypass win.
//! * `deflate_level6` — deflate-compressed input. Regression guard: the
//!   bypass logic must not perturb the compressed path.
//!
//! Producers like `samtools view -u`, htsjdk's level-0 writer, and our own
//! `InlineBgzfCompressor::new(0)` emit stored blocks, so the `stored_bypass`
//! number is what real "uncompressed BAM" pipelines see.
//!
//! Run with: `cargo bench -p fgumi-bgzf --bench stored_block_decode`
//! Quick sanity: `cargo bench -p fgumi-bgzf --bench stored_block_decode -- --quick`

#![deny(unsafe_code)]
#![allow(clippy::cast_possible_truncation)]

use std::hint::black_box;
use std::io::Cursor;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use fgumi_bgzf::reader::{RawBgzfBlock, decompress_block_into, read_raw_blocks};
use fgumi_bgzf::writer::InlineBgzfCompressor;
use libdeflater::Decompressor;

/// Total uncompressed payload size per bench iteration. Sized to:
/// - cross many BGZF blocks (each is up to 64 KiB), so per-block overheads
///   are properly amortised;
/// - stay well under L2 so we're measuring decode work, not RAM bandwidth.
const PAYLOAD_BYTES: usize = 4 * 1024 * 1024;

/// How many blocks to grab per `read_raw_blocks` call. Matches the batch
/// constant used elsewhere in the pipeline.
const BLOCKS_PER_BATCH: usize = 64;

/// Generate a payload that doesn't compress trivially — otherwise the deflate
/// control would collapse to near-zero and not exercise libdeflater realistically.
/// Using a simple xorshift-style PRNG keeps the bench self-contained.
fn make_payload(total_bytes: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(total_bytes);
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    while out.len() + 8 <= total_bytes {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        out.extend_from_slice(&state.to_le_bytes());
    }
    out.resize(total_bytes, 0);
    out
}

/// Encode the payload as a concatenated BGZF stream at the given level.
fn encode_stream(payload: &[u8], level: u32) -> Vec<u8> {
    let mut compressor = InlineBgzfCompressor::new(level);
    compressor.write_all(payload).expect("write");
    compressor.flush().expect("flush");
    let blocks = compressor.take_blocks();
    let total: usize = blocks.iter().map(|b| b.data.len()).sum();
    let mut out = Vec::with_capacity(total);
    for b in &blocks {
        out.extend_from_slice(&b.data);
    }
    out
}

/// Drive the full read + decompress path over a pre-encoded BGZF stream,
/// reusing the decompressor and output buffer between block batches. Uses
/// the public `decompress_block_into`, which transparently picks the stored
/// bypass when applicable.
fn decode_stream(stream: &[u8], scratch: &mut Vec<u8>) {
    scratch.clear();
    let mut reader = Cursor::new(stream);
    let mut decompressor = Decompressor::new();
    loop {
        let blocks = read_raw_blocks(&mut reader, BLOCKS_PER_BATCH).expect("read_raw_blocks");
        if blocks.is_empty() {
            break;
        }
        for block in &blocks {
            decompress_block_into(block, &mut decompressor, scratch).expect("decompress");
        }
    }
    black_box(&*scratch);
}

/// Pre-bypass baseline: call `libdeflater::deflate_decompress` directly on
/// every block, mimicking what `decompress_block_into` did before the
/// stored-block fast path was added. Only used by the bench to size the win.
fn decode_stream_via_libdeflater(stream: &[u8], scratch: &mut Vec<u8>) {
    scratch.clear();
    let mut reader = Cursor::new(stream);
    let mut decompressor = Decompressor::new();
    loop {
        let blocks = read_raw_blocks(&mut reader, BLOCKS_PER_BATCH).expect("read_raw_blocks");
        if blocks.is_empty() {
            break;
        }
        for block in &blocks {
            decompress_block_libdeflater(block, &mut decompressor, scratch);
        }
    }
    black_box(&*scratch);
}

/// Decompress a single block via libdeflater, with no stored-block bypass —
/// the path the reader took before this bench was added. Includes the same
/// CRC32 verification the production code performs, so the comparison to
/// the bypass row stays apples-to-apples (the bypass measures the cost of
/// libdeflater's stored-block memcpy, not whether we verify the CRC).
fn decompress_block_libdeflater(
    block: &RawBgzfBlock,
    decompressor: &mut Decompressor,
    output: &mut Vec<u8>,
) {
    let compressed = block.compressed_data();
    let uncompressed_size = block.uncompressed_size();
    let expected_crc = block.crc32();
    if uncompressed_size == 0 {
        return;
    }
    let start = output.len();
    output.resize(start + uncompressed_size, 0);
    let n = decompressor
        .deflate_decompress(compressed, &mut output[start..])
        .expect("libdeflater decompress");
    assert_eq!(n, uncompressed_size);
    let actual_crc = crc32fast::hash(&output[start..start + uncompressed_size]);
    assert_eq!(actual_crc, expected_crc, "CRC32 mismatch in baseline path");
}

fn bench_decode(c: &mut Criterion) {
    let payload = make_payload(PAYLOAD_BYTES);

    let stored_stream = encode_stream(&payload, 0);
    let deflate_stream = encode_stream(&payload, 6);

    let mut group = c.benchmark_group("bgzf_decode");
    group.throughput(Throughput::Bytes(PAYLOAD_BYTES as u64));

    // Report compressed-stream size in the bench id so a "X% faster" headline
    // can't be misread when the two streams differ a lot in size.
    let stored_id = format!("{}KiB-stream", stored_stream.len() / 1024);
    let deflate_id = format!("{}KiB-stream", deflate_stream.len() / 1024);

    // The "after" — what we ship.
    group.bench_with_input(
        BenchmarkId::new("stored_bypass", &stored_id),
        &stored_stream,
        |b, stream| {
            let mut scratch = Vec::with_capacity(PAYLOAD_BYTES);
            b.iter(|| decode_stream(stream, &mut scratch));
        },
    );

    // The "before" — same stored input, but forced through libdeflater.
    group.bench_with_input(
        BenchmarkId::new("stored_via_libdeflater", &stored_id),
        &stored_stream,
        |b, stream| {
            let mut scratch = Vec::with_capacity(PAYLOAD_BYTES);
            b.iter(|| decode_stream_via_libdeflater(stream, &mut scratch));
        },
    );

    // Regression guard for the compressed path.
    group.bench_with_input(
        BenchmarkId::new("deflate_level6", &deflate_id),
        &deflate_stream,
        |b, stream| {
            let mut scratch = Vec::with_capacity(PAYLOAD_BYTES);
            b.iter(|| decode_stream(stream, &mut scratch));
        },
    );

    group.finish();
}

criterion_group!(benches, bench_decode);
criterion_main!(benches);
