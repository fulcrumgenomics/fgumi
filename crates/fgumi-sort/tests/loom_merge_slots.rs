//! Loom model-check of the block-parallel Phase-2 decompress protocol
//! (`file_granularity == false`), driving the **real** `SortMergeSlot`.
//!
//! # What this verifies (and what it does NOT)
//!
//! Under `--cfg loom`, `merge_slots.rs` swaps `std::sync` -> `loom::sync`, so
//! the `SortMergeSlot` constructed here uses loom's atomics and mutexes. loom
//! explores the thread interleavings and the memory reorderings the C11 model
//! permits (preemption-bounded — see "Preemption-bounded exploration" below),
//! running the REAL slot methods each time:
//!
//!   * [`SortMergeSlot::bp_commit_read`] — the publish order (reserve
//!     `in_flight` before setting `reader_eof`) that the original silent-
//!     truncation bug lived in. Production
//!     (`SortSpillDecompress::try_fill_block_parallel_slot`) calls the SAME
//!     method, so the model and the code cannot drift.
//!   * [`SortMergeSlot::bp_insert_drain_finalize`] /
//!     [`SortMergeSlot::bp_drain_and_finalize`] and the private
//!     `drain_locked_and_finalize` (the finalize predicate `!queue_eof &&
//!     reader_eof && in_flight == 0 && reorder.is_empty()`, the lock order
//!     reorder -> decompressed, and the real `in_flight.fetch_sub(AcqRel)`).
//!   * The slot's REAL `reader` mutex serializes reads, and its REAL `reorder`
//!     buffer ([`fgumi_bam_io::reorder::ReorderBuffer`]) reassembles them.
//!
//! Only two things are *not* the production code, both sound and documented:
//!
//!   * **The "read" is simulated.** loom cannot model real file I/O, so a
//!     worker computes `(start_seq, got, hit_eof)` from a test-held total under
//!     the real `reader` lock instead of calling `read_raw`; "decompression"
//!     yields the seq number as the payload. The slot's accounting/finalize
//!     methods that run on the result are 100% production code.
//!   * **Blocking `reader.lock()` instead of production's `try_lock()`.** The
//!     property under test depends only on reads being *serialized* (so
//!     `reader_eof` can never become visible while an unreserved block still
//!     exists); a blocking lock preserves exactly that while collapsing the
//!     `try_lock`-miss-retry fan-out that would otherwise explode loom's state
//!     space (and a spin-retry is a loom anti-pattern). It is a sound
//!     over-approximation of the serialization invariant.
//!
//! ## Modeling choices (documented for honesty)
//!
//!   * **Window/FIFO admission disabled.** `bp_reorder_admits` backpressure is a
//!     *bounded-memory* property, separately covered by the `merge_slots` unit
//!     test `bp_reorder_window_is_bounded_under_straggler`. Workers here always
//!     admit, keeping the model focused on the EOF/truncation invariants.
//!   * **Tiny sizes (2-4 blocks, 2-3 workers).** loom's state space is
//!     super-exponential; these sizes still exercise every out-of-order
//!     insert/finalize interleaving the per-slot protocol can hit (more blocks
//!     add only more of the same kind, not a new kind).
//!   * **One pass per producer, no spin.** `n_workers == reads_needed(N,
//!     batch)`, so every block is read by exactly one producer and the model
//!     always terminates; a protocol that failed to finalize surfaces as the
//!     `queue_eof` assertion below, not as a hang.
//!   * **A concurrent merge consumer** ([`consume_until_drained`]) polls the
//!     FIFO and STOPS at `is_drained()`, mirroring `SortMerge`. This is what
//!     makes a *premature* `queue_eof` observable as truncation (without it the
//!     straggler is appended after the join and the bug hides — verified: the
//!     `bp_commit_read` order swap is caught only with the consumer present).
//!   * **Preemption-bounded exploration.** The consumer's poll loop adds a
//!     scheduling point per iteration; combined with 2-3 producers the fully
//!     exhaustive state space is minutes-long. Every model is therefore explored
//!     under a preemption bound (`Some(k)`) — a recognized technique: essentially
//!     all real concurrency bugs (including the truncation race this guards)
//!     manifest with ≤2-3 preemptions. The bound is verified to still catch the
//!     `bp_commit_read` order swap.
//!
//! Run with:
//! ```text
//! RUSTFLAGS="--cfg loom" cargo test -p fgumi-sort --test loom_merge_slots --release
//! ```
//!
//! # Complementary coverage and its residual
//!
//! This model is one leg of the block-parallel hardening; the others are the
//! `merge_slots` unit tests (bounded-memory window), the
//! `fgumi-pipeline-io` granularity/proptest/soak-matrix tests (the real
//! end-to-end pipeline over real spill files), and a `ThreadSanitizer` pass.
//!
//! **Sanitizer residual (recorded for honesty):** the `ThreadSanitizer` run
//! exercised the real pipeline on **arm64 only**, and the C decompression codecs
//! (`zstd` via the `zstd` crate, `libdeflate` via `libdeflater`) are
//! **uninstrumented** — the sanitizer only sees the Rust side, so a data race
//! *inside* a C codec would be missed. This is acceptable because the codecs are
//! pure per-block transforms with no shared mutable state across threads (each
//! worker decompresses its own block into its own buffer); the cross-thread
//! protocol the sanitizer and loom actually need to clear is the Rust-side slot
//! accounting, which is fully instrumented here.

#![cfg(loom)]
#![deny(unsafe_code)]
// Block counts/seqs in this model are tiny (≤ a handful) and always fit a
// usize; the casts below are between u64 model seqs and usize counts.
#![allow(clippy::cast_possible_truncation)]

use std::fs::File;
use std::io::BufReader;

use fgumi_sort::{SortMergeSlot, SpillCodec};
use loom::sync::Arc;
use loom::sync::atomic::Ordering;

/// A throwaway reader for the slot. The block-parallel slot methods never touch
/// `reader.inner`; the model serializes reads on the `reader` mutex and computes
/// the read result arithmetically, so an empty file is all the struct needs.
fn empty_reader() -> BufReader<File> {
    BufReader::new(tempfile::tempfile().expect("create tempfile"))
}

/// One worker's body: mirrors a single `try_run` of
/// `SortSpillDecompress::try_fill_block_parallel_slot`, but with the file read
/// simulated (see module docs). Reads up to `block_batch` blocks under the REAL
/// `reader` lock, publishes the accounting via the REAL
/// [`SortMergeSlot::bp_commit_read`], "decompresses" outside the lock, then
/// inserts/drains/finalizes via the REAL [`SortMergeSlot::bp_insert_drain_finalize`].
/// A worker that finds the reader already at EOF falls through to the REAL
/// Phase-B drain-only [`SortMergeSlot::bp_drain_and_finalize`].
fn worker_one_pass(slot: &SortMergeSlot, block_batch: u64, total_blocks: u64) {
    if slot.queue_eof.load(Ordering::Acquire) {
        return;
    }
    let mut did_phase_a = false;
    if !slot.reader_eof.load(Ordering::Acquire) {
        // Blocking lock (sound over-approximation of production's `try_lock`;
        // see module docs) — serializes reads on the REAL reader mutex.
        let mut reader = slot.reader.lock().unwrap();
        // Re-check under the lock: another worker may have hit EOF.
        if !slot.reader_eof.load(Ordering::Acquire) {
            let start_seq = reader.next_seq;
            let remaining = total_blocks - start_seq;
            let got = block_batch.min(remaining);
            let hit_eof = got < block_batch;
            // Stamp the read range and commit the accounting BEFORE releasing
            // the lock, via the real publish-order method.
            reader.next_seq += got;
            slot.bp_commit_read(got as usize, hit_eof);
            drop(reader);

            // "Decompress" outside the reader lock: the payload is the seq as
            // 8 little-endian bytes, so the drained FIFO can be checked for
            // in-order, no-loss delivery.
            let blocks: Vec<Vec<u8>> =
                (start_seq..start_seq + got).map(|s| s.to_le_bytes().to_vec()).collect();
            slot.bp_insert_drain_finalize(start_seq, blocks, got as usize);
            did_phase_a = true;
        }
    }
    if !did_phase_a {
        slot.bp_drain_and_finalize();
    }
}

/// Number of reader-lock acquisitions (= worker passes) needed to read every
/// block and then observe the clean EOF: `ceil((N + 1) / batch)`. The `+ 1`
/// accounts for the read that returns fewer than `batch` blocks (possibly
/// empty), which is what sets `reader_eof`.
fn reads_needed(total_blocks: u64, block_batch: u64) -> usize {
    ((total_blocks + 1).div_ceil(block_batch)) as usize
}

/// Decode an 8-byte little-endian seq payload back to its sequence number.
fn seq_of(block: &[u8]) -> u64 {
    let mut buf = [0u8; 8];
    buf.copy_from_slice(&block[..8]);
    u64::from_le_bytes(buf)
}

/// The merge consumer, mirroring `SortMerge`/`slot_try_load_block`: pop every
/// available block, then STOP the instant the slot looks cleanly drained
/// (`is_drained()` == `queue_eof && FIFO empty && !error`). Returns the seqs it
/// collected, in pop (delivery) order.
///
/// This stop condition is what makes a *premature* `queue_eof` observable: if
/// the slot finalizes EOF while a block is still outstanding (the truncation the
/// publish-order protocol prevents), the consumer sees an empty FIFO + EOF and
/// quits early, so `run_model`'s completeness check fails. Without a consumer
/// the straggler would still be appended after the join and the bug would hide.
///
/// `yield_now` between polls is the loom scheduling point; the model is finite
/// (every producer runs exactly one pass), so the wait always terminates.
fn consume_until_drained(slot: &SortMergeSlot) -> Vec<u64> {
    let mut collected = Vec::new();
    loop {
        loop {
            let popped = slot.decompressed.lock().unwrap().pop_front();
            match popped {
                Some(b) => collected.push(seq_of(&b)),
                None => break,
            }
        }
        if slot.is_drained() {
            break;
        }
        loom::thread::yield_now();
    }
    collected
}

/// Drive one worker pass per required read over `total_blocks` blocks with
/// `block_batch` blocks per read PLUS a concurrent merge consumer, under loom,
/// and assert the no-loss / in-order / clean-EOF invariants for every
/// interleaving against the REAL slot state.
fn run_model(total_blocks: u64, block_batch: u64) {
    let slot = Arc::new(SortMergeSlot::new(0, empty_reader(), SpillCodec::Bgzf));
    let n_workers = reads_needed(total_blocks, block_batch);

    let mut handles: Vec<_> = (0..n_workers)
        .map(|_| {
            let slot = Arc::clone(&slot);
            loom::thread::spawn(move || worker_one_pass(&slot, block_batch, total_blocks))
        })
        .collect();
    let consumer = {
        let slot = Arc::clone(&slot);
        loom::thread::spawn(move || consume_until_drained(&slot))
    };
    for h in handles.drain(..) {
        h.join().unwrap();
    }
    let delivered = consumer.join().unwrap();

    // Post-conditions: clean EOF reached, no error, nothing left in flight or
    // buffered, and the consumer collected every block exactly once, in read
    // order, BEFORE it observed the clean EOF (a premature `queue_eof` truncates
    // `delivered`).
    assert!(slot.queue_eof.load(Ordering::Acquire), "slot never reached queue_eof");
    assert!(!slot.has_error(), "spurious decomp_error");
    assert_eq!(slot.in_flight.load(Ordering::Acquire), 0, "blocks left in flight at EOF");
    assert!(slot.reorder.lock().unwrap().is_empty(), "reorder buffer not drained at EOF");
    assert!(slot.decompressed.lock().unwrap().is_empty(), "FIFO not fully consumed at EOF");

    let expected: Vec<u64> = (0..total_blocks).collect();
    assert_eq!(
        delivered, expected,
        "blocks lost, duplicated, reordered, or truncated by early EOF"
    );
}

/// Run `f` under loom with at most `preemption_bound` preemptions. See the
/// module-level "Preemption-bounded exploration" note for why every model is
/// bounded rather than exhaustive.
fn check_model<F: Fn() + Sync + Send + 'static>(preemption_bound: usize, f: F) {
    let mut builder = loom::model::Builder::new();
    builder.preemption_bound = Some(preemption_bound);
    builder.check(f);
}

/// Three blocks, batch 2 ⇒ 2 producer passes: the second read is SHORT (carries
/// seq 2 AND sets `reader_eof`) while the first read's blocks (seq 0,1) may
/// still be in flight — the exact `bp_eof_with_straggler` shape.
#[test]
fn loom_three_blocks_batch2() {
    check_model(3, || run_model(3, 2));
}

/// Two blocks, batch 2 ⇒ 2 producer passes: the second read is the EMPTY
/// EOF-detecting read (got == 0) that sets `reader_eof` while the first read's
/// blocks (seq 0,1) may still be in flight. Exercises the `count == 0` finalize
/// path of `bp_insert_drain_finalize`.
#[test]
fn loom_two_blocks_batch2_empty_eof_read() {
    check_model(3, || run_model(2, 2));
}

/// Four blocks, batch 3 ⇒ 2 producer passes (read seq 0,1,2; short read seq 3
/// sets EOF). A larger in-flight batch straggling behind the EOF read.
#[test]
fn loom_four_blocks_batch3() {
    check_model(3, || run_model(4, 3));
}

/// Two blocks, batch 1 ⇒ 3 producer passes (read seq 0, read seq 1, empty EOF
/// read) plus the consumer — four concurrent threads. Covers the three-way race
/// between two in-flight blocks and the EOF-setter. Bounded tighter (the
/// four-thread × nested-mutex state space is the largest here).
#[test]
fn loom_two_blocks_batch1_three_workers_bounded() {
    check_model(2, || run_model(2, 1));
}

// ── decomp-error-beats-clean-EOF (#399) ──────────────────────────────────────

/// The consumer that observes `queue_eof` under the `decompressed` mutex
/// (`is_drained`) must also observe `decomp_error` whenever the producer took
/// the error path — i.e. a failed slot can never be mistaken for a clean EOF.
///
/// This drives the REAL slot: the producer mirrors
/// `SortSpillDecompress::mark_slot_failed` (store `decomp_error` then
/// `queue_eof`, both under the `decompressed` mutex), and the consumer is the
/// real [`SortMergeSlot::is_drained`] / [`SortMergeSlot::has_error`] pair. The
/// mutex release-acquire makes the two flags jointly visible regardless of which
/// the consumer reads first.
#[test]
fn loom_decomp_error_beats_clean_eof() {
    loom::model(|| {
        let slot = Arc::new(SortMergeSlot::new(0, empty_reader(), SpillCodec::Bgzf));

        // Producer: error path. Mirrors `mark_slot_failed` — set decomp_error
        // THEN queue_eof, both under the `decompressed` mutex.
        let producer = {
            let slot = Arc::clone(&slot);
            loom::thread::spawn(move || {
                let _g = slot.decompressed.lock().unwrap();
                slot.decomp_error.store(true, Ordering::Release);
                slot.queue_eof.store(true, Ordering::Release);
            })
        };

        // Consumer: the real `is_drained()` / `has_error()`. The invariant: if
        // the slot looks drained (clean EOF), it must NOT have an error — but
        // because the producer's only path is the error path, any observation of
        // `queue_eof == true` must come with `decomp_error == true`, so
        // `is_drained()` must return false and `has_error()` true.
        let consumer = {
            let slot = Arc::clone(&slot);
            loom::thread::spawn(move || {
                let drained = slot.is_drained();
                let errored = slot.has_error();
                // A clean drain that hid the error is the bug this guards.
                assert!(!(drained && !errored), "observed clean EOF that hid a decomp error");
            })
        };

        producer.join().unwrap();
        consumer.join().unwrap();
    });
}
