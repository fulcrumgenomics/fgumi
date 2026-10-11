//! Soak of the merge slot's raw stash (`stress-tests`): a 10 K-block zstd
//! spill read as randomly cut, out-of-order slices by one ingest thread, four
//! claiming workers, and a consumer that serves itself from the stash whenever
//! its FIFO is empty, under a random FIFO cap. Every block must be delivered
//! exactly once in file order and every slice must return to its pool.

use std::sync::atomic::Ordering;

use fgumi_bam_io::pread::SliceBufferPool;
use fgumi_sort::{SpillBlockDecompressor, SpillCodec};

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

/// One soak run with `seed`.
fn soak_once(dir: &std::path::Path, bytes: &[u8], oracle: &[Vec<u8>], seed: u64) {
    let path = dir.join("run.spill");
    let slot = fgumi_sort::open_spill_slot(&path, 0).unwrap();
    let caps = [1u32, 2, 8];
    slot.fifo_cap.store(caps[usize::try_from(seed % 3).unwrap()], Ordering::Relaxed);
    let pool = SliceBufferPool::new(16);
    let body = &bytes[4..];
    // Random cuts of 1 B .. 64 KiB, ingested in a locally shuffled order.
    let mut cuts = Vec::new();
    let (mut pos, mut r) = (0usize, seed);
    while pos < body.len() {
        r = splitmix64(r);
        let len = (1 + usize::try_from(r % 65_536).unwrap()).min(body.len() - pos);
        cuts.push((pos, len));
        pos += len;
    }
    let mut order: Vec<usize> = (0..cuts.len()).collect();
    for w in order.chunks_mut(4) {
        r = splitmix64(r);
        let k = usize::try_from(r).unwrap_or(0) % w.len();
        w.rotate_left(k);
    }
    let n = cuts.len();
    let got = std::thread::scope(|sc| {
        for _ in 0..4 {
            let slot = &slot;
            sc.spawn(move || {
                let mut dec = SpillBlockDecompressor::new();
                while !slot.queue_eof.load(Ordering::Acquire) {
                    if let Some(b) = slot.bp_claim_raw(1 << 20) {
                        let d = dec.decompress_one(SpillCodec::Zstd, &b.frame).unwrap();
                        slot.bp_insert_drain_finalize(b.seq, vec![d], 1);
                    } else {
                        slot.bp_drain_and_finalize();
                        std::thread::yield_now();
                    }
                }
            });
        }
        let (slot2, pool2) = (&slot, &pool);
        let cuts2 = &cuts;
        let order2 = &order;
        sc.spawn(move || {
            for &i in order2 {
                let (off, len) = cuts2[i];
                slot2.bp_note_issued(len as u64);
                let lease = pool2.lease(body[off..off + len].to_vec());
                slot2.bp_ingest_slice(u32::try_from(i).unwrap(), lease, i == n - 1).unwrap();
            }
        });
        // The consumer: pop; on an empty FIFO serve itself from the stash
        // (claim → decompress → insert), else run the drain-only pass.
        let mut dec = SpillBlockDecompressor::new();
        let mut got = Vec::new();
        loop {
            let popped = slot.pop_decompressed();
            if let Some(b) = popped {
                got.push(b);
                continue;
            }
            if slot.is_drained() {
                break;
            }
            assert!(!slot.has_error());
            if let Some(b) = slot.bp_claim_raw(1 << 20) {
                let d = dec.decompress_one(SpillCodec::Zstd, &b.frame).unwrap();
                slot.bp_insert_drain_finalize(b.seq, vec![d], 1);
            } else {
                slot.bp_drain_and_finalize();
                std::thread::yield_now();
            }
        }
        got
    });
    assert_eq!(got.len(), oracle.len(), "seed {seed}");
    assert!(got == oracle, "seed {seed}: blocks differ");
    drop(slot);
    assert_eq!(pool.resident_bytes(), 0, "seed {seed}: a slice was never returned");
}

#[test]
fn slot_stash_soak_matrix() {
    let dir = tempfile::tempdir().unwrap();
    let mut c = fgumi_sort::SpillBlockCompressor::new(SpillCodec::Zstd, 1).unwrap();
    let mut bytes = fgumi_sort::spill_magic(SpillCodec::Zstd).to_vec();
    let mut raw_blocks = Vec::new();
    for b in 0..10_000u64 {
        let raw: Vec<u8> = (0..200).map(|i| (splitmix64(b * 1000 + i) & 0x3f) as u8).collect();
        bytes.extend_from_slice(&c.compress_block(&raw).unwrap());
        raw_blocks.push(raw);
    }
    std::fs::write(dir.path().join("run.spill"), &bytes).unwrap();
    let runs = if cfg!(coverage) { 2 } else { 100 };
    for seed in 0..runs {
        soak_once(dir.path(), &bytes, &raw_blocks, seed);
    }
}
