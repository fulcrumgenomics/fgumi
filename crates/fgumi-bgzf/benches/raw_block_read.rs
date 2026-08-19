//! Where the sort pipeline's input reader spends its per-block time.
//!
//! `fgumi sort`'s Phase 1 has one thread doing all input reading (worker 0 owns
//! `ReadInputBlocks` exclusively). On the production cell it occupies that thread
//! for 124.5s across 5,114,408 blocks — **24.3 µs per 8.4 KB block** — while the
//! same volume at the measured single-stream disk ceiling is only ~71s. That
//! leaves ~10.8 µs per block that no disk model explains.
//!
//! This bench bounds how much of that gap userspace framing can possibly account
//! for, by running the production framing path over a BGZF stream already in RAM.
//! Four arms partition the cost additively:
//!
//! * `full` — `read_raw_blocks` exactly as the pipeline calls it. The production path.
//! * `reuse_buf` — the same framing, but reading into one reused `Vec` instead of a
//!   fresh `Vec::with_capacity(block_size)` per block. `full - reuse_buf` is what the
//!   5.11M per-block allocations cost.
//! * `frame_only` — parse and validate each header, then skip the body via
//!   `BufRead::consume` rather than copying it. `reuse_buf - frame_only` is the body
//!   memcpy.
//! * `drain` — walk every byte through the same `BufReader` with no framing at all.
//!   The floor. `frame_only - drain` is header parse, validation, and loop overhead.
//!
//! **What this bench cannot tell you:** whether the production reader actually pays
//! those costs. Reads here come from a `Cursor` over RAM, so the underlying refill is
//! a memcpy rather than a syscall against page cache or EBS. That makes `drain` an
//! optimistic floor by construction. The bench answers "can framing cost 10.8 µs per
//! block?", not "does it, on that host" — in-situ instrumentation answers the latter.
//!
//! Block payloads are pseudo-random filler because this path never inflates them:
//! `read_raw_blocks` frames on `BSIZE` and copies bytes. Only block *size* is
//! load-bearing, so it is swept — 4 KiB, the production mean, 16 KiB, and the 64 KiB
//! maximum — which separates fixed per-block overhead from cost proportional to bytes.
//!
//! Run with: `cargo bench -p fgumi-bgzf --bench raw_block_read`
//! Quick sanity: `cargo bench -p fgumi-bgzf --bench raw_block_read -- --quick`

#![deny(unsafe_code)]

use std::hint::black_box;
use std::io::{self, BufRead, BufReader, Cursor, Read};

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use fgumi_bgzf::header::{BGZF_HEADER_SIZE, MIN_BLOCK_SIZE, block_size_checked, validate};
use fgumi_bgzf::reader::{BGZF_EOF, read_raw_blocks};

/// Input buffer the sort pipeline wraps its BAM in (`SORT_INPUT_BUFFER_SIZE`).
const INPUT_BUFFER_SIZE: usize = 2 * 1024 * 1024;

/// Blocks requested per `read_raw_blocks` call (`INPUT_READ_BATCH_SIZE`).
const BLOCKS_PER_BATCH: usize = 16;

/// Stream size per iteration. Large enough that the 2 MiB buffer refills ~32 times
/// and the per-block allocations cannot all be served from a warm free list, small
/// enough to keep each arm's measurement under a second.
const STREAM_BYTES: usize = 64 * 1024 * 1024;

/// Mean compressed block size on the production cell: 43,131,188,552 bytes across
/// 5,114,408 blocks.
const PRODUCTION_BLOCK_SIZE: usize = 8433;

/// Compressed block sizes to sweep.
const BLOCK_SIZES: [usize; 4] = [4096, PRODUCTION_BLOCK_SIZE, 16 * 1024, 64 * 1024];

/// The fixed 16 bytes of a BGZF header that precede `BSIZE`, per `header::validate`:
/// gzip magic, deflate, `FEXTRA` alone, empty MTIME/XFL, unknown OS, `XLEN` 6, and
/// the `BC` subfield id and length.
const HEADER_PREFIX: [u8; 16] =
    [0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00];

/// Build a BGZF stream of fixed-size blocks, terminated by the EOF marker.
///
/// Bodies are pseudo-random filler: the framing path under test never inflates them.
fn build_stream(block_size: usize, total_bytes: usize) -> Vec<u8> {
    assert!(block_size >= MIN_BLOCK_SIZE, "block size below the BGZF floor");
    let bsize = u16::try_from(block_size - 1).expect("block size must fit in BSIZE");

    let body_len = block_size - BGZF_HEADER_SIZE;
    let mut body = Vec::with_capacity(body_len + 8);
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    while body.len() < body_len {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        body.extend_from_slice(&state.to_le_bytes());
    }
    body.truncate(body_len);

    let mut out = Vec::with_capacity(total_bytes + BGZF_EOF.len());
    while out.len() + block_size <= total_bytes {
        out.extend_from_slice(&HEADER_PREFIX);
        out.extend_from_slice(&bsize.to_le_bytes());
        out.extend_from_slice(&body);
    }
    out.extend_from_slice(&BGZF_EOF);
    out
}

/// Wrap a stream in the same buffered reader the sort pipeline uses.
fn buffered(stream: &[u8]) -> BufReader<Cursor<&[u8]>> {
    BufReader::with_capacity(INPUT_BUFFER_SIZE, Cursor::new(stream))
}

/// `full`: the production path, `read_raw_blocks` in batches of 16.
fn arm_full(stream: &[u8]) -> usize {
    let mut reader = buffered(stream);
    let mut blocks = 0usize;
    loop {
        let batch = read_raw_blocks(&mut reader, BLOCKS_PER_BATCH).expect("framing");
        if batch.is_empty() {
            break;
        }
        blocks += batch.len();
        black_box(&batch);
    }
    blocks
}

/// Frame one block into `buf`, reusing its allocation. Mirrors `read_raw_block`'s
/// header handling exactly so the delta against `full` is only the allocation.
///
/// Returns `false` at a clean EOF.
fn read_block_into(reader: &mut BufReader<Cursor<&[u8]>>, buf: &mut Vec<u8>) -> io::Result<bool> {
    let mut header = [0u8; BGZF_HEADER_SIZE];
    loop {
        match reader.read(&mut header[..1]) {
            Ok(0) => return Ok(false),
            Ok(_) => break,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    reader.read_exact(&mut header[1..])?;
    validate(&header).map_err(|r| io::Error::new(io::ErrorKind::InvalidData, r))?;
    let block_size = block_size_checked(&header)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "BGZF block too small"))?;

    buf.clear();
    buf.extend_from_slice(&header);
    let remaining = u64::try_from(block_size - BGZF_HEADER_SIZE).expect("block size fits in u64");
    reader.take(remaining).read_to_end(buf)?;
    if buf.len() != block_size {
        return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated BGZF block"));
    }
    Ok(true)
}

/// `reuse_buf`: identical framing, one reused block buffer.
fn arm_reuse_buf(stream: &[u8]) -> usize {
    let mut reader = buffered(stream);
    let mut buf: Vec<u8> = Vec::with_capacity(BGZF_MAX_FRAME);
    let mut blocks = 0usize;
    while read_block_into(&mut reader, &mut buf).expect("framing") {
        if buf.as_slice() == BGZF_EOF {
            continue;
        }
        blocks += 1;
        black_box(buf.as_slice());
    }
    blocks
}

/// Largest block a `BSIZE` can describe, so `reuse_buf` never reallocates.
const BGZF_MAX_FRAME: usize = 65536;

/// Advance `reader` by `n` bytes without copying them out of its buffer.
fn skip_exact(reader: &mut BufReader<Cursor<&[u8]>>, mut n: usize) -> io::Result<()> {
    while n > 0 {
        let available = reader.fill_buf()?;
        if available.is_empty() {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated BGZF block"));
        }
        let taken = available.len().min(n);
        reader.consume(taken);
        n -= taken;
    }
    Ok(())
}

/// `frame_only`: parse and validate every header, skip every body.
fn arm_frame_only(stream: &[u8]) -> usize {
    let mut reader = buffered(stream);
    let mut header = [0u8; BGZF_HEADER_SIZE];
    let mut blocks = 0usize;
    loop {
        match reader.read(&mut header[..1]) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => panic!("framing: {e}"),
        }
        reader.read_exact(&mut header[1..]).expect("header");
        validate(&header).expect("valid header");
        let block_size = block_size_checked(&header).expect("block size");
        skip_exact(&mut reader, block_size - BGZF_HEADER_SIZE).expect("body");
        blocks += 1;
        black_box(&header);
    }
    blocks
}

/// `drain`: walk every byte through the same `BufReader`, no framing.
fn arm_drain(stream: &[u8]) -> usize {
    let mut reader = buffered(stream);
    let mut total = 0usize;
    loop {
        let available = reader.fill_buf().expect("fill");
        let n = available.len();
        if n == 0 {
            break;
        }
        black_box(available);
        reader.consume(n);
        total += n;
    }
    total
}

/// Confirm the synthetic stream frames the way every arm assumes before any timing
/// is taken — otherwise an arm that silently bails early looks like a fast one.
fn check_arms_agree(stream: &[u8], block_size: usize) -> usize {
    let expected = STREAM_BYTES / block_size;
    let full = arm_full(stream);
    assert!(full > 0, "synthetic stream framed zero blocks");
    assert_eq!(full, arm_reuse_buf(stream), "reuse_buf framed a different block count");
    // `frame_only` counts the trailing EOF marker, which `read_raw_blocks` drops.
    assert_eq!(full + 1, arm_frame_only(stream), "frame_only framed a different block count");
    assert_eq!(stream.len(), arm_drain(stream), "drain walked a different byte count");
    assert!(
        full.abs_diff(expected) <= 1,
        "expected ~{expected} blocks of {block_size} bytes, framed {full}"
    );
    full
}

fn bench_raw_block_read(c: &mut Criterion) {
    let mut group = c.benchmark_group("raw_block_read");
    for block_size in BLOCK_SIZES {
        let stream = build_stream(block_size, STREAM_BYTES);
        let blocks = check_arms_agree(&stream, block_size);
        group.throughput(Throughput::Elements(u64::try_from(blocks).expect("block count")));

        for (name, arm) in [
            ("full", arm_full as fn(&[u8]) -> usize),
            ("reuse_buf", arm_reuse_buf),
            ("frame_only", arm_frame_only),
            ("drain", arm_drain),
        ] {
            group.bench_with_input(BenchmarkId::new(name, block_size), &stream, |b, stream| {
                b.iter(|| black_box(arm(black_box(stream))));
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_raw_block_read);
criterion_main!(benches);
