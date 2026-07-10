use super::*;
use std::io::BufReader;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use fgumi_sort::{SortMergeSlot, SpillCodec};

/// A resolved-to-zero output budget must not disable the reorder-window byte cap:
/// `bp_reorder_admits` treats `window_budget == 0` as unlimited, so `new()`
/// substitutes the default cap (mirroring the legacy `effective_limit`
/// 0-normalization). A nonzero budget passes through unchanged.
#[test]
fn zero_output_byte_limit_normalizes_reorder_window() {
    let zero = SortSpillDecompress::new(0, SortDecompressTuning::default());
    assert_eq!(zero.window_budget, DEFAULT_REORDER_WINDOW_BYTES);
    assert_ne!(zero.window_budget, 0, "reorder window must stay bounded on a zero budget");

    let nonzero = SortSpillDecompress::new(4 * 1024 * 1024, SortDecompressTuning::default());
    assert_eq!(nonzero.window_budget, 4 * 1024 * 1024, "nonzero budget passes through unchanged");
}

#[test]
fn admission_counter_caps_concurrency() {
    let active = Arc::new(AtomicUsize::new(0));
    let max = Some(2usize);
    // Acquire up to the cap; hold the permits so the count accumulates.
    let p1 = try_acquire(&active, max);
    assert!(p1.is_some());
    let p2 = try_acquire(&active, max);
    assert!(p2.is_some());
    assert!(try_acquire(&active, max).is_none()); // at cap
    drop(p1); // releasing one permit frees a slot
    assert!(try_acquire(&active, max).is_some()); // freed one slot
    drop(p2);
}

#[test]
fn admission_counter_unbounded_when_none() {
    let active = Arc::new(AtomicUsize::new(0));
    for _ in 0..1000 {
        assert!(try_acquire(&active, None).is_some());
    }
}

/// Under real thread contention, the shared counter never exceeds the cap and
/// every acquired permit is released via `DecompressPermit::drop` (the counter
/// returns to zero). Mirrors how `new_worker_copy` clones share one `active`.
#[test]
fn admission_counter_concurrent_never_exceeds_cap() {
    use std::thread;

    let active = Arc::new(AtomicUsize::new(0));
    let cap = 3usize;
    let handles: Vec<_> = (0..16)
        .map(|_| {
            let active = Arc::clone(&active);
            thread::spawn(move || {
                for _ in 0..5000 {
                    // The returned permit IS the ownership token; its Drop at the
                    // end of the block exercises the decrement.
                    if let Some(_permit) = try_acquire(&active, Some(cap)) {
                        // Occupancy observed while holding a permit can never
                        // exceed the cap: `try_acquire` only increments past a
                        // CAS that checks `cur < cap`, and the count only drops
                        // otherwise.
                        let occupancy = active.load(Ordering::Acquire);
                        assert!(occupancy <= cap, "occupancy {occupancy} exceeded cap {cap}");
                    }
                }
            })
        })
        .collect();
    for h in handles {
        h.join().expect("worker thread panicked");
    }
    assert_eq!(active.load(Ordering::Acquire), 0, "every permit must be released on drop");
}

// Most coverage for the decompress step lives in sort/tests.rs (the oracle parity
// suite drives the whole chain). This unit test pins the emptiest-first refill
// ordering in isolation.

#[test]
fn emptiest_first_order_sorts_by_fifo_len_ascending() {
    let mk = |file_id: u32, nblocks: usize| {
        let s = Arc::new(SortMergeSlot::new(
            file_id,
            BufReader::new(tempfile::tempfile().expect("tempfile")),
            SpillCodec::Bgzf,
        ));
        for _ in 0..nblocks {
            s.decompressed.lock().expect("decompressed lock").push_back(vec![0u8]);
        }
        s
    };
    // FIFO depths 5, 1, 3 ⇒ most-starved-first visit order is indices 1, 2, 0.
    let slots = vec![mk(0, 5), mk(1, 1), mk(2, 3)];
    assert_eq!(SortSpillDecompress::emptiest_first_order(&slots), vec![1, 2, 0]);
}

/// A poisoned `reader` mutex (a fill worker panicked while holding the lock)
/// must fail the slot CLOSED — `try_fill_*_slot` returns `Err` and sets
/// `decomp_error`/`queue_eof` — rather than being swallowed as `WouldBlock` and
/// skipped forever. If it were skipped, `queue_eof` would never be set and
/// `SortMerge` would spin on `Contention` and deadlock instead of surfacing the
/// panic. Regression test for the poisoned-vs-would-block conflation.
#[test]
fn poisoned_reader_lock_fails_closed_rather_than_hanging() {
    let make_poisoned_slot = || {
        let slot = Arc::new(SortMergeSlot::new(
            0,
            BufReader::new(tempfile::tempfile().expect("tempfile")),
            SpillCodec::Bgzf,
        ));
        let holder = Arc::clone(&slot);
        // Panic while holding the reader lock; joining the panicked thread leaves
        // the mutex poisoned (mirrors a fill worker dying mid-read).
        let _ = std::thread::spawn(move || {
            let _guard = holder.reader.lock().expect("acquire reader lock");
            panic!("simulated fill-worker panic under the reader lock");
        })
        .join();
        assert!(slot.reader.is_poisoned(), "precondition: reader mutex is poisoned");
        slot
    };

    let mut dec = SortSpillDecompress::new(4 * 1024 * 1024, SortDecompressTuning::default());

    // Inline path.
    let inline_slot = make_poisoned_slot();
    let inline = dec.try_fill_inline_slot(&inline_slot);
    assert!(inline.is_err(), "inline path: poisoned reader must return Err, not Ok(false)");
    assert!(inline_slot.decomp_error.load(Ordering::Acquire), "inline: decomp_error set");
    assert!(inline_slot.queue_eof.load(Ordering::Acquire), "inline: queue_eof set");

    // Block-parallel path.
    let bp_slot = make_poisoned_slot();
    let bp = dec.try_fill_block_parallel_slot(&bp_slot);
    assert!(bp.is_err(), "block-parallel path: poisoned reader must return Err, not skip");
    assert!(bp_slot.decomp_error.load(Ordering::Acquire), "bp: decomp_error set");
    assert!(bp_slot.queue_eof.load(Ordering::Acquire), "bp: queue_eof set");
}
