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
//!     the real `reader` lock (`SortMergeSlot::bp_read_batch_for_test`, which
//!     publishes through the real `bp_commit_read`) instead of parsing a
//!     slice; "decompression" yields the seq number as the payload. The slot's
//!     accounting/finalize methods that run on the result are 100% production
//!     code.
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
//!   * **Window/FIFO admission disabled.** The reorder-window backpressure
//!     claims obey is a *bounded-memory* property, separately covered by the
//!     `merge_slots` unit tests `bp_reorder_window_is_bounded_under_straggler`
//!     and `claims_respect_the_reorder_window_under_a_straggler`. Workers here always
//!     admit, keeping the model focused on the EOF/truncation invariants.
//!   * **Tiny sizes (2-4 blocks, 2-3 workers).** loom's state space is
//!     super-exponential; these sizes still exercise every out-of-order
//!     insert/finalize interleaving the per-slot protocol can hit (more blocks
//!     add only more of the same kind, not a new kind).
//!   * **One pass per producer.** `n_workers == reads_needed(N, batch)`, so
//!     every block is read by exactly one producer and every producer
//!     terminates. The *consumer* does spin: [`consume_until_drained`] breaks
//!     only on `is_drained()`, which requires `queue_eof`. So a protocol that
//!     never finalizes does NOT reach [`run_model`]'s assertions — the consumer
//!     never returns from its poll loop, `consumer.join()` never completes, and
//!     the `queue_eof` assertion below is unreachable. It surfaces instead as a
//!     loom `max_branches` panic from the spinning consumer, which names the
//!     harness rather than the protocol. Read that failure as "`queue_eof` was
//!     never finalized", not as a model that needs a bigger branch budget.
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

use fgumi_sort::{SortMergeSlot, SpillCodec};
use loom::sync::Arc;

/// A slot over an empty file. The block-parallel slot methods never read the
/// file; the model serializes reads on the `reader` mutex and computes the read
/// result arithmetically, so an empty file is all the struct needs.
fn empty_slot() -> SortMergeSlot {
    SortMergeSlot::for_test(0, SpillCodec::Bgzf)
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
    if slot.queue_eof() {
        return;
    }
    let mut did_phase_a = false;
    if !slot.reader_eof() {
        // The REAL reader lock serializes the simulated read, which stamps the
        // range and commits the accounting through the real publish order
        // before releasing the lock (`None`: another worker hit EOF first).
        if let Some((start_seq, got)) = slot.bp_read_batch_for_test(total_blocks, block_batch) {
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
/// `yield_now` between polls is the loom scheduling point. Every producer runs
/// exactly one pass, so against a protocol that finalizes `queue_eof` this wait
/// terminates. Against one that does not, it spins forever — see the module
/// header's "One pass per producer" note for why that surfaces as a loom
/// `max_branches` panic rather than as the assertions in [`run_model`].
fn consume_until_drained(slot: &SortMergeSlot) -> Vec<u64> {
    let mut collected = Vec::new();
    loop {
        loop {
            let popped = slot.pop_decompressed();
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
    let slot = Arc::new(empty_slot());
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
    assert!(slot.queue_eof(), "slot never reached queue_eof");
    assert!(!slot.has_error(), "spurious decomp_error");
    assert_eq!(slot.in_flight(), 0, "blocks left in flight at EOF");
    assert_eq!(slot.reorder_len_relaxed(), 0, "reorder buffer not drained at EOF");
    assert_eq!(slot.fifo_len(), 0, "FIFO not fully consumed at EOF");

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

// ── merge wake: no lost wakeup between `await_slot` and a delivery ───────────

/// The merge-wake lost-wakeup model over the REAL consumer sequence
/// ([`fgumi_sort::MergeDemand::await_slot`], which `SortMerge` calls on every
/// stall) and the REAL producer sequence
/// ([`SortMergeSlot::bp_insert_drain_finalize`] then
/// [`fgumi_sort::MergeDemand::notify_delivered`], as
/// `SortSpillDecompress::try_fill_block_parallel_slot` does). The consumer
/// parks only when `await_slot` says it may; if any interleaving lets it park
/// with nobody left to unpark it, loom reports the deadlock. Swapping
/// `await_slot`'s `set_awaited` and its re-check fails this model.
#[test]
fn loom_merge_wake_never_lost() {
    check_model(3, || {
        let slot = Arc::new(empty_slot());
        let demand = Arc::new(fgumi_sort::MergeDemand::new());
        slot.bp_commit_read(1, true); // reserve one block, as the reader would
        let (s2, d2) = (Arc::clone(&slot), Arc::clone(&demand));
        let producer = loom::thread::spawn(move || {
            assert!(s2.bp_insert_drain_finalize(0, vec![0u64.to_le_bytes().to_vec()], 1));
            d2.notify_delivered(0);
        });
        while !demand.await_slot(&slot) {
            loom::thread::park();
        }
        producer.join().unwrap();
        assert_eq!(slot.fifo_len(), 1);
    });
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
        let slot = Arc::new(empty_slot());

        // Producer: error path. Mirrors `mark_slot_failed` — set decomp_error
        // THEN queue_eof, both under the `decompressed` mutex.
        let producer = {
            let slot = Arc::clone(&slot);
            loom::thread::spawn(move || {
                slot.mark_failed();
            })
        };

        // Consumer: the real `is_drained()` / `has_error()`.
        //
        // The assertion is `!drained`, flatly, not `!(drained && !errored)`.
        // The producer's only path is the error path, so there is no
        // interleaving in which a clean drain is legal: before the producer
        // runs, `queue_eof` is unset and `is_drained()` is false; after it,
        // `decomp_error` is set and `is_drained()` must still be false.
        //
        // The composite form is what the old assertion used, and it could not
        // fail: with the `decomp_error` guard present `drained` is always
        // false, and with it deleted `errored` is true, so `drained && !errored`
        // is unsatisfiable either way — the test stayed green against the exact
        // bug it names. Asserting `!drained` alone restores the discrimination;
        // deleting the guard from `is_drained` now fails this model.
        //
        // `has_error()` is still called, before `is_drained()` in one variant
        // and after in the other, so the documented order-independence of the
        // two reads is what the model actually exercises.
        let consumer = {
            let slot = Arc::clone(&slot);
            loom::thread::spawn(move || {
                let errored = slot.has_error();
                let drained = slot.is_drained();
                assert!(!drained, "errored slot reported a clean drain");
                assert!(!drained || errored, "clean EOF hid a decomp error");
            })
        };

        let consumer_reversed = {
            let slot = Arc::clone(&slot);
            loom::thread::spawn(move || {
                let drained = slot.is_drained();
                let _errored = slot.has_error();
                assert!(!drained, "errored slot reported a clean drain (reversed read order)");
            })
        };

        producer.join().unwrap();
        consumer.join().unwrap();
        consumer_reversed.join().unwrap();
    });
}

// ── raw stash: claims, the front escape, EOF with stash outstanding ─────────

/// The frame payload a model stashes for `seq`: the seq as 8 LE bytes, so the
/// consumer can check in-order, exactly-once delivery with [`seq_of`].
fn payload(seq: u64) -> Vec<u8> {
    seq.to_le_bytes().to_vec()
}

/// One claimer pass (bounded, as every model producer is — see the module
/// header): claim the stash head and publish it ("decompress" = copy), or run
/// the drain-only pass when nothing is claimable.
fn claim_once(slot: &SortMergeSlot) {
    match slot.bp_claim_raw(u64::MAX) {
        Some(b) => {
            let d = b.frame.to_vec();
            slot.bp_insert_drain_finalize(b.seq, vec![d], 1);
        }
        None => {
            slot.bp_drain_and_finalize();
        }
    }
}

/// One consumer pass: pop everything queued, then look at the drain state. A
/// premature `queue_eof` shows up as "drained" with blocks still missing.
fn consume_once(slot: &SortMergeSlot, total: u64) -> Vec<u64> {
    let mut got = Vec::new();
    loop {
        let popped = slot.pop_decompressed();
        match popped {
            Some(b) => got.push(seq_of(&b)),
            None => break,
        }
    }
    if slot.is_drained() {
        assert_eq!(got, (0..total).collect::<Vec<_>>(), "clean EOF before every block landed");
    }
    got
}

/// Finish what the bounded passes left (claims, drains) single-threaded and
/// collect the rest of the FIFO.
fn finish(slot: &SortMergeSlot) -> Vec<u64> {
    while let Some(b) = slot.bp_claim_raw(u64::MAX) {
        slot.bp_insert_drain_finalize(b.seq, vec![b.frame.to_vec()], 1);
    }
    slot.bp_drain_and_finalize();
    let mut got = Vec::new();
    while let Some(b) = slot.pop_decompressed() {
        got.push(seq_of(&b));
    }
    got
}

/// One ingest of three frames (last) races two claimer passes and a consumer
/// pass: every seq is delivered exactly once, in order, then `queue_eof` —
/// never a clean EOF with a block missing.
#[test]
fn loom_stash_claims_deliver_each_seq_once() {
    check_model(3, || {
        let slot = Arc::new(empty_slot());
        let s = Arc::clone(&slot);
        let ingest = loom::thread::spawn(move || {
            s.bp_stash_frames_for_test((0..3).map(payload).collect(), true)
        });
        let claimers: Vec<_> = (0..2)
            .map(|_| {
                let s = Arc::clone(&slot);
                loom::thread::spawn(move || claim_once(&s))
            })
            .collect();
        let s = Arc::clone(&slot);
        let consumer = loom::thread::spawn(move || consume_once(&s, 3));
        ingest.join().unwrap();
        for c in claimers {
            c.join().unwrap();
        }
        let mut delivered = consumer.join().unwrap();
        delivered.extend(finish(&slot));
        assert_eq!(delivered, vec![0, 1, 2], "lost, duplicated or reordered");
        assert!(slot.queue_eof() && !slot.has_error());
    });
}

/// The front escape. With the FIFO at its cap (1) and the window budget at 1
/// byte, nothing but the escape admits a claim: the stash head that IS the
/// reorder front must be claimable anyway (a claimer that loops on it
/// terminates without the consumer popping), while a head that is not the
/// front is refused until the front lands.
#[test]
fn loom_front_escape_unsticks_the_window() {
    check_model(3, || {
        let slot = Arc::new(empty_slot());
        slot.set_fifo_cap(1);
        slot.bp_stash_frames_for_test((0..3).map(payload).collect(), false);
        let b0 = slot.bp_claim_raw(1).expect("the front");
        slot.bp_insert_drain_finalize(b0.seq, vec![payload(0)], 1);
        assert_eq!(slot.fifo_len(), 1, "the FIFO is at its cap");
        let s = Arc::clone(&slot);
        let claimer = loom::thread::spawn(move || {
            loop {
                if let Some(b) = s.bp_claim_raw(1) {
                    return b.detach();
                }
                loom::thread::yield_now();
            }
        });
        let b1 = claimer.join().unwrap();
        assert_eq!(b1.seq, 1, "the front is admitted over a full FIFO and window");
        assert!(slot.bp_claim_raw(1).is_none(), "seq 2 is not the front (1 is in flight)");
        // The consumer pops seq 0, so seq 1 drains into the FIFO (refilling it
        // to its cap) and seq 2 becomes the front; only the escape admits it.
        assert_eq!(slot.pop_decompressed(), Some(payload(0)));
        let s = Arc::clone(&slot);
        let inserter = loom::thread::spawn(move || {
            s.bp_insert_drain_finalize(b1.seq, vec![payload(1)], 1);
        });
        let s = Arc::clone(&slot);
        let claimer = loom::thread::spawn(move || {
            loop {
                if let Some(b) = s.bp_claim_raw(1) {
                    return b.seq;
                }
                loom::thread::yield_now();
            }
        });
        inserter.join().unwrap();
        assert_eq!(claimer.join().unwrap(), 2, "seq 2 is admitted once it is the front");
    });
}

/// The last slice's EOF commit races a claim of an earlier block: `queue_eof`
/// is set only after every block (including the one in flight) is inserted.
#[test]
fn loom_eof_slice_with_stash_outstanding() {
    check_model(3, || {
        let slot = Arc::new(empty_slot());
        slot.bp_stash_frames_for_test(vec![payload(0)], false);
        let s = Arc::clone(&slot);
        let claimer = loom::thread::spawn(move || claim_once(&s));
        let s = Arc::clone(&slot);
        let last = loom::thread::spawn(move || s.bp_stash_frames_for_test(vec![payload(1)], true));
        let s = Arc::clone(&slot);
        let consumer = loom::thread::spawn(move || consume_once(&s, 2));
        claimer.join().unwrap();
        last.join().unwrap();
        let mut delivered = consumer.join().unwrap();
        delivered.extend(finish(&slot));
        assert_eq!(delivered, vec![0, 1], "EOF finalized before a claimed block landed");
        assert!(slot.queue_eof());
    });
}

/// A self-serving consumer's decompress failure (`mark_failed`) while a worker
/// claim is in flight: the slot is never reported as a clean drain, and once
/// the consumer has failed its claimed block the slot is visibly failed —
/// `queue_eof` set with `decomp_error` — so a merge waiting on it surfaces the
/// error instead of waiting forever for a block that will never be inserted.
#[test]
fn loom_decomp_error_beats_clean_eof_with_self_serve() {
    check_model(3, || {
        let slot = Arc::new(empty_slot());
        slot.bp_stash_frames_for_test((0..2).map(payload).collect(), true);
        let s = Arc::clone(&slot);
        let worker = loom::thread::spawn(move || {
            if let Some(b) = s.bp_claim_raw(u64::MAX) {
                s.bp_insert_drain_finalize(b.seq, vec![b.frame.to_vec()], 1);
            }
        });
        let s = Arc::clone(&slot);
        let consumer = loom::thread::spawn(move || {
            let claimed = s.bp_claim_raw(u64::MAX).is_some();
            if claimed {
                s.mark_failed();
            }
            claimed
        });
        let s = Arc::clone(&slot);
        let observer = loom::thread::spawn(move || {
            assert!(!s.is_drained(), "a slot that will fail reported a clean drain");
        });
        worker.join().unwrap();
        let consumer_claimed = consumer.join().unwrap();
        observer.join().unwrap();
        assert!(!slot.is_drained());
        if consumer_claimed {
            assert!(
                slot.has_error() && slot.queue_eof(),
                "a failed self-serve must leave the slot failed and at EOF"
            );
        }
    });
}

/// The front escape admits the last block while the FIFO is full; it waits in
/// `reorder`. The consumer pops, then runs the drain-only pass, and reaches
/// `queue_eof` in every interleaving with the worker's insert.
#[test]
fn loom_final_insert_into_full_fifo_drains() {
    check_model(3, || {
        let cap = fgumi_sort::PHASE2_DECOMP_CAP as u64;
        let slot = Arc::new(empty_slot());
        slot.bp_stash_frames_for_test((0..=cap).map(payload).collect(), true);
        for _ in 0..cap {
            let b = slot.bp_claim_raw(u64::MAX).unwrap();
            slot.bp_insert_drain_finalize(b.seq, vec![b.frame.to_vec()], 1);
        }
        let last = slot.bp_claim_raw(u64::MAX).expect("the front over a full FIFO").detach();
        let s = Arc::clone(&slot);
        let worker = loom::thread::spawn(move || {
            s.bp_insert_drain_finalize(last.seq, vec![last.frame.to_vec()], 1);
        });
        let s = Arc::clone(&slot);
        let consumer = loom::thread::spawn(move || {
            let mut got = Vec::new();
            loop {
                loop {
                    let popped = s.pop_decompressed();
                    match popped {
                        Some(b) => got.push(seq_of(&b)),
                        None => break,
                    }
                }
                if s.is_drained() {
                    return got;
                }
                s.bp_drain_and_finalize();
                loom::thread::yield_now();
            }
        });
        worker.join().unwrap();
        let got = consumer.join().unwrap();
        assert_eq!(got, (0..=cap).collect::<Vec<_>>());
    });
}

/// An ingest racing the merge's stall classification: with no claim, a stall
/// sampled at any point of the ingest reads the slot as stashed (or not yet
/// stashed), never as decompressing — `stash_len` is stored before the
/// in-flight reservation is raised, and `awaited_state` loads `in_flight`
/// first.
#[test]
fn loom_awaited_state_never_reads_a_fresh_stash_as_decompressing() {
    check_model(3, || {
        let slot = Arc::new(empty_slot());
        let s = Arc::clone(&slot);
        let ingest = loom::thread::spawn(move || {
            s.bp_stash_frames_for_test((0..2).map(payload).collect(), false);
        });
        let state = slot.awaited_state();
        assert_ne!(state, fgumi_sort::AwaitedSlotState::Decompressing, "no block was claimed");
        ingest.join().unwrap();
        assert_eq!(slot.awaited_state(), fgumi_sort::AwaitedSlotState::Stashed);
    });
}

/// The read planner bounds read-ahead by each slot's charge (`stash_bytes`),
/// so the charge must balance however ingests, claims and the drops that
/// release a claimed frame's charge interleave: a claimer races a second
/// ingest, and once every frame is claimed and dropped the charge is zero —
/// never underflowed, never left behind.
#[test]
fn loom_stash_charge_balances_across_claims_and_ingests() {
    check_model(3, || {
        let slot = Arc::new(empty_slot());
        slot.bp_stash_frames_for_test(vec![payload(0)], false);
        let s = Arc::clone(&slot);
        let claimer = loom::thread::spawn(move || {
            if let Some(b) = s.bp_claim_raw(u64::MAX) {
                let charged = s.stash_bytes();
                drop(b);
                assert!(charged >= 8, "a claimed frame stays charged until dropped");
            }
        });
        let s = Arc::clone(&slot);
        let ingest = loom::thread::spawn(move || {
            s.bp_stash_frames_for_test(vec![payload(1)], true);
        });
        claimer.join().unwrap();
        ingest.join().unwrap();
        while let Some(b) = slot.bp_claim_raw(u64::MAX) {
            drop(b);
        }
        assert_eq!(slot.stash_bytes(), 0, "every charge released exactly once");
    });
}
