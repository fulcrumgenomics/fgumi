//! Arena-lifetime microbenchmark — settles N_arena (in-flight sort arenas) and
//! the gather/compress release point, isolated from the typed-step pipeline.
// Throwaway diagnostic bench (not production): relax pedantic lints.
#![allow(
    clippy::doc_markdown,
    clippy::too_many_lines,
    clippy::map_unwrap_or,
    clippy::cast_possible_truncation
)]
//!
//! Models the proposed Phase-1 dataflow with a **reused arena pool** (no
//! mem-take-and-mint): a producer fills + par-sorts one `RecordBuffer` to a
//! memory limit, hands ownership to a consumer that gathers (frames sorted
//! records into ≤64 KiB blocks) + compresses, then returns the buffer to the
//! pool. The pool depth is `N_arena`; the release point is `split`
//! (return the buffer right after gather, before compress) or `fused`
//! (return only after compress).
//!
//! Uses the REAL kernels: `RecordBuffer` (fill/par_sort/refs/get_record),
//! `frame_keyed_record_into`, `SpillBlockCompressor::compress_block`.
//!
//! Run one config per process (clean peak RSS):
//!   cargo run --release -p fgumi-sort --example arena_bench -- <mode> <N_arena> <chunk_mib> <n_chunks> <workers>
//! e.g. cargo run --release -p fgumi-sort --example arena_bench -- split 1 512 6 4
//!
//! Output (stdout, one TSV line): mode N_arena chunk_mib n_chunks workers wall_s fill_s sort_s gather_s compress_s
//! Peak RSS is sampled externally by the runner wrapper (this process prints its PID first).

use std::sync::mpsc::sync_channel;
use std::time::{Duration, Instant};

use fgumi_sort::{
    RawCoordinateKey, RawSortKey, RecordBuffer, SpillBlockCompressor, SpillCodec,
    frame_keyed_record_into,
};

const BLOCK_SIZE: usize = 65280; // BGZF_MAX_BLOCK_SIZE
const RECORD_BODY: usize = 200; // ~realistic BAM record body bytes

/// Per-record framing overhead `frame_keyed_record_into` adds: a 4-byte `u32`
/// record-length prefix, plus the serialized key ONLY when the key is not
/// already embedded in the record. Derived from `RawCoordinateKey`'s own trait
/// constants so it tracks the framing (and key size) instead of hardcoding `12`,
/// which would silently mis-size blocks if the key encoding ever changed.
const FRAME_OVERHEAD: usize = 4 // u32 record-length prefix
    + if RawCoordinateKey::EMBEDDED_IN_RECORD {
        0
    } else {
        // The coordinate key is fixed-size; a variable-length key (`None`) can't
        // be statically bounded, but this bench only frames `RawCoordinateKey`.
        match RawCoordinateKey::SERIALIZED_SIZE {
            Some(n) => n,
            None => 0,
        }
    };

/// Synthesize one ~200-byte BAM-ish record into `out`: tid (bytes 0-3), pos
/// (4-7), then a semi-compressible payload so zstd does realistic work. `n` is
/// the record ordinal (varies tid/pos so keys differ).
fn synth_record(out: &mut Vec<u8>, n: u64) {
    out.clear();
    let tid: i32 = (n % 25) as i32; // 25 contigs
    let pos: i32 = ((n.wrapping_mul(2_654_435_761) >> 8) % 250_000_000) as i32;
    out.extend_from_slice(&tid.to_le_bytes());
    out.extend_from_slice(&pos.to_le_bytes());
    // bytes 8..16 — remaining fixed BAM header fields (arbitrary)
    out.extend_from_slice(&[0u8; 8]);
    // payload: semi-compressible (a few repeating motifs keyed off n)
    let motif = [(n & 0xFF) as u8, b'A', b'C', b'G', b'T', ((n >> 3) & 0xFF) as u8];
    while out.len() < RECORD_BODY {
        out.extend_from_slice(&motif);
    }
    out.truncate(RECORD_BODY);
}

/// Frame the sorted records of `buf` into ≤64 KiB raw blocks, cut-before-overflow.
/// Reads `buf.refs()` (sorted) + `buf.get_record()` — the real gather access pattern.
fn gather_range(buf: &RecordBuffer, start: usize, end: usize) -> Vec<Vec<u8>> {
    let refs = buf.refs();
    let mut blocks = Vec::new();
    let mut cur = Vec::with_capacity(BLOCK_SIZE + 1024);
    for r in &refs[start..end] {
        let key = RawCoordinateKey { sort_key: r.sort_key };
        let body = buf.get_record(r);
        // would this record overflow a non-empty block? cut first.
        if !cur.is_empty() && cur.len() + FRAME_OVERHEAD + body.len() > BLOCK_SIZE {
            blocks.push(std::mem::replace(&mut cur, Vec::with_capacity(BLOCK_SIZE + 1024)));
        }
        frame_keyed_record_into(&mut cur, &key, body).expect("frame");
    }
    if !cur.is_empty() {
        blocks.push(cur);
    }
    blocks
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mode = args.get(1).map(String::as_str).unwrap_or("split").to_string();
    let n_arena: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1);
    let chunk_mib: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(512);
    let n_chunks: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(6);
    let workers: usize = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(4);
    let codec = SpillCodec::Zstd;
    let level = 1u32;
    let mem_limit = chunk_mib * 1024 * 1024;

    eprintln!("PID {}", std::process::id());
    eprintln!(
        "config: mode={mode} N_arena={n_arena} chunk_mib={chunk_mib} n_chunks={n_chunks} workers={workers} codec=zstd{level}"
    );

    // Per-phase timing accumulators (nanos).
    let fill_ns = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
    let sort_ns = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
    let gather_ns = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
    let compress_ns = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));

    // Arena pool: a free-list channel seeded with N_arena empty buffers, and a
    // work channel of filled+sorted buffers. Bounded so the producer blocks when
    // no arena is free (this IS the N_arena bound).
    let (free_tx, free_rx) = sync_channel::<RecordBuffer>(n_arena);
    let (work_tx, work_rx) = sync_channel::<RecordBuffer>(n_arena);
    // Estimate records/arena for with_capacity pre-sizing.
    let est_records = mem_limit / (RECORD_BODY + 8 + 16);
    for _ in 0..n_arena {
        free_tx.send(RecordBuffer::with_capacity(est_records, mem_limit, 25)).expect("seed pool");
    }

    let wall = Instant::now();

    // Producer: fill + par_sort each chunk, hand the owned buffer to the consumer.
    let prod_fill = fill_ns.clone();
    let prod_sort = sort_ns.clone();
    let producer = std::thread::spawn(move || {
        let mut scratch = Vec::with_capacity(RECORD_BODY + 8);
        let mut ordinal: u64 = 0;
        for _ in 0..n_chunks {
            let mut buf = free_rx.recv().expect("acquire arena");
            let t = Instant::now();
            while buf.memory_usage() < mem_limit {
                synth_record(&mut scratch, ordinal);
                ordinal += 1;
                buf.push_coordinate(&scratch).expect("push");
            }
            prod_fill
                .fetch_add(t.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
            let t = Instant::now();
            buf.par_sort();
            prod_sort
                .fetch_add(t.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
            work_tx.send(buf).expect("hand to consumer");
        }
        // drop work_tx → consumer loop ends
    });

    // Consumer: gather + compress each chunk; release the arena per `mode`.
    let split = mode == "split";
    while let Ok(mut buf) = work_rx.recv() {
        let n = buf.len();
        // contiguous index ranges, one per worker
        let per = n.div_ceil(workers.max(1));
        let ranges: Vec<(usize, usize)> = (0..workers)
            .map(|w| (w * per, ((w + 1) * per).min(n)))
            .filter(|(s, e)| s < e)
            .collect();

        if split {
            // GATHER (parallel) → owned raw blocks; then RELEASE arena; then COMPRESS (parallel).
            let tg = Instant::now();
            let raw: Vec<Vec<u8>> = std::thread::scope(|sc| {
                let bref = &buf;
                let handles: Vec<_> = ranges
                    .iter()
                    .map(|&(s, e)| sc.spawn(move || gather_range(bref, s, e)))
                    .collect();
                handles.into_iter().flat_map(|h| h.join().expect("gather join")).collect()
            });
            gather_ns
                .fetch_add(tg.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
            // RELEASE the arena BEFORE compress (the whole point of "split").
            buf.clear();
            let _ = free_tx.send(buf);
            // COMPRESS the owned raw blocks in parallel (arena already free).
            let tc = Instant::now();
            let chunks: Vec<&[Vec<u8>]> =
                raw.chunks(raw.len().div_ceil(workers.max(1)).max(1)).collect();
            std::thread::scope(|sc| {
                for ch in &chunks {
                    sc.spawn(move || {
                        let mut comp = SpillBlockCompressor::new(codec, level).expect("compressor");
                        for blk in *ch {
                            let out = comp.compress_block(blk).expect("compress");
                            std::hint::black_box(&out);
                        }
                    });
                }
            });
            compress_ns
                .fetch_add(tc.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
        } else {
            // FUSED: each worker gathers+compresses its range from the arena;
            // the arena is released only after ALL workers finish.
            let tf = Instant::now();
            std::thread::scope(|sc| {
                for &(s, e) in &ranges {
                    let bref = &buf;
                    sc.spawn(move || {
                        let mut comp = SpillBlockCompressor::new(codec, level).expect("compressor");
                        for blk in gather_range(bref, s, e) {
                            let out = comp.compress_block(&blk).expect("compress");
                            std::hint::black_box(&out);
                        }
                    });
                }
            });
            // fused gather+compress are interleaved; attribute the whole span to compress.
            compress_ns
                .fetch_add(tf.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
            buf.clear();
            let _ = free_tx.send(buf);
        }
        let _ = n;
    }
    producer.join().expect("producer join");
    let wall_s = wall.elapsed().as_secs_f64();

    let g = |a: &std::sync::atomic::AtomicU64| {
        Duration::from_nanos(a.load(std::sync::atomic::Ordering::Relaxed)).as_secs_f64()
    };
    // headline TSV
    println!(
        "{mode}\t{n_arena}\t{chunk_mib}\t{n_chunks}\t{workers}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}",
        wall_s,
        g(&fill_ns),
        g(&sort_ns),
        g(&gather_ns),
        g(&compress_ns)
    );
}
