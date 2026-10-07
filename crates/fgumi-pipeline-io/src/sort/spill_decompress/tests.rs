use super::*;
use std::io::BufReader;
use std::sync::Arc;
use std::sync::atomic::Ordering;

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

/// A fresh step is uncapped (the chain builder decides the phase), and
/// `with_phase_cap` installs the cap on the seed and shares the same counter
/// with every clone.
#[test]
fn with_phase_cap_replaces_the_uncapped_default_on_seed_and_clones() {
    use fgumi_pipeline_core::PhaseCap;
    let plain = SortSpillDecompress::new(4096, SortDecompressTuning::default());
    assert!(plain.phase_cap().is_none(), "uncapped by default");

    let cap = PhaseCap::new("sort-phase2", 1);
    let capped = plain.with_phase_cap(Some(Arc::clone(&cap)));
    assert_eq!(capped.phase_cap().map(PhaseCap::describe).as_deref(), Some("sort-phase2(1)"));
    let clone = capped.clone();
    let _held = cap.try_acquire().expect("the only permit");
    assert!(
        clone.phase_cap().expect("capped").try_acquire().is_none(),
        "the clone shares the exhausted counter"
    );
}

/// The permit gates only the slot fill: with the only permit held elsewhere
/// the step fills nothing and reports `Capped` while a slot is live; once the
/// permit is released it fills (here: reads the empty spill to EOF).
#[test]
fn full_cap_refuses_the_fill_then_fills_after_release() {
    use fgumi_pipeline_core::testing::StepProbe;
    use fgumi_pipeline_core::{PhaseCap, StepOutcome};
    let cap = PhaseCap::new("sort-phase2", 1);
    let mut step = SortSpillDecompress::new(4 * 1024 * 1024, SortDecompressTuning::default())
        .with_phase_cap(Some(Arc::clone(&cap)));
    let slot = Arc::new(SortMergeSlot::new(
        0,
        BufReader::new(tempfile::tempfile().expect("tempfile")),
        SpillCodec::Zstd,
    ));
    step.registry.lock().push(RegisteredSpill { slot: Arc::clone(&slot) });
    let probe = StepProbe::new(&step);

    let held = cap.try_acquire().expect("the only permit");
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::Capped);
    assert!(!slot.queue_eof.load(Ordering::Acquire), "a refused clone must not fill the slot");

    drop(held);
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::Progress);
    assert!(slot.queue_eof.load(Ordering::Acquire), "the admitted clone fills the slot to EOF");
    assert_eq!(cap.active(), 0, "the permit is released when try_run returns");
}

/// An idle scan takes no phase-2 permit: with no live slot and the only permit
/// held elsewhere (by the output compressor, say), the step reports
/// `NoProgress` — not `Capped` — and counts no refusal, so idle decompress
/// clones can never crowd the output compressor out of the shared cap.
#[test]
fn idle_scan_with_no_live_slot_takes_no_permit() {
    use fgumi_pipeline_core::testing::StepProbe;
    use fgumi_pipeline_core::{PhaseCap, StepOutcome};
    let cap = PhaseCap::new("sort-phase2", 1);
    let mut step = SortSpillDecompress::new(4 * 1024 * 1024, SortDecompressTuning::default())
        .with_phase_cap(Some(Arc::clone(&cap)));
    let probe = StepProbe::new(&step);
    let held = cap.try_acquire().expect("the only permit");
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::NoProgress);
    assert_eq!(cap.refused(), 0, "no slot to fill, so no admission attempt");
    probe.close_input();
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::Finished);
    assert_eq!(cap.refused(), 0);
    drop(held);
}

/// A live slot whose FIFO is full cannot be filled, so a capped poll takes no
/// permit, is not a refusal, and reports `NoProgress`; once the merge drains a
/// block the next poll is admitted.
#[test]
fn live_but_unfillable_slot_takes_no_permit() {
    use fgumi_pipeline_core::testing::StepProbe;
    use fgumi_pipeline_core::{PhaseCap, StepOutcome};
    let cap = PhaseCap::new("sort-phase2", 1);
    let mut step = SortSpillDecompress::new(4 * 1024 * 1024, SortDecompressTuning::default())
        .with_phase_cap(Some(Arc::clone(&cap)));
    let slot = Arc::new(SortMergeSlot::new(
        0,
        BufReader::new(tempfile::tempfile().expect("tempfile")),
        SpillCodec::Zstd,
    ));
    for _ in 0..fgumi_sort::PHASE2_DECOMP_CAP {
        slot.decompressed.lock().expect("decompressed lock").push_back(vec![0u8]);
    }
    step.registry.lock().push(RegisteredSpill { slot: Arc::clone(&slot) });
    let probe = StepProbe::new(&step);

    let held = cap.try_acquire().expect("the only permit");
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::NoProgress);
    assert_eq!(cap.refused(), 0, "an unfillable slot is not a refusal");
    drop(held);
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::NoProgress);
    assert_eq!(slot.fifo_len(), fgumi_sort::PHASE2_DECOMP_CAP, "nothing was read");

    slot.decompressed.lock().expect("decompressed lock").pop_front();
    let held = cap.try_acquire().expect("the only permit");
    assert_eq!(
        probe.try_run(&mut step).expect("try_run"),
        StepOutcome::Capped,
        "with room, the poll needs the permit"
    );
    assert_eq!(cap.refused(), 1);
    drop(held);
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::Progress);
}

/// Uncapped, the step keeps its pre-cap behaviour: no registry scan before the
/// fill, and a live slot that could not be filled (here: FIFO full) is
/// contention, not idleness; with no live slot and input drained it finishes.
#[test]
fn uncapped_poll_keeps_the_atomic_only_liveness_check() {
    use fgumi_pipeline_core::StepOutcome;
    use fgumi_pipeline_core::testing::StepProbe;
    let mut step = SortSpillDecompress::new(4 * 1024 * 1024, SortDecompressTuning::default());
    let slot = Arc::new(SortMergeSlot::new(
        0,
        BufReader::new(tempfile::tempfile().expect("tempfile")),
        SpillCodec::Zstd,
    ));
    for _ in 0..fgumi_sort::PHASE2_DECOMP_CAP {
        slot.decompressed.lock().expect("decompressed lock").push_back(vec![0u8]);
    }
    step.registry.lock().push(RegisteredSpill { slot: Arc::clone(&slot) });
    let probe = StepProbe::new(&step);
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::Contention);
    assert_eq!(slot.fifo_len(), fgumi_sort::PHASE2_DECOMP_CAP, "nothing was read");

    slot.queue_eof.store(true, Ordering::Release);
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::NoProgress);
    probe.close_input();
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::Finished);
}

/// Forwarding input events to the merge is bookkeeping, not phase-2 work: it
/// runs with the cap full and takes no permit, so the merge is never starved
/// of announcements behind decompressing clones.
#[test]
fn input_events_are_forwarded_while_the_cap_is_full() {
    use fgumi_pipeline_core::testing::StepProbe;
    use fgumi_pipeline_core::{InputHandle, PhaseCap, StepOutcome};
    let cap = PhaseCap::new("sort-phase2", 1);
    let mut step = SortSpillDecompress::new(4 * 1024 * 1024, SortDecompressTuning::default())
        .with_phase_cap(Some(Arc::clone(&cap)));
    let mut probe = StepProbe::new(&step);
    let out = probe.take_output::<SortPhase2Event>(0);
    probe.push_input(SortPhase1Event::AllAnnounced {
        slot_count: 0,
        memory_chunk_count: 0,
        total_records: 0,
    });
    let held = cap.try_acquire().expect("the only permit");
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::Progress);
    assert!(
        matches!(out.pop(), Some(SortPhase2Event::AllAnnounced { .. })),
        "the event is forwarded while the only permit is held elsewhere"
    );
    assert_eq!(cap.active(), 1, "forwarding took no permit");
    assert_eq!(cap.refused(), 0);
    drop(held);
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

/// Slots that already signalled `queue_eof` are dropped from the refill scan.
///
/// The registry is append-only, so without this filter a drained slot stays in
/// the scan for the rest of the run: every dispatch clones its `Arc`, takes its
/// FIFO lock to sort, and is then rejected immediately by the `queue_eof` guard
/// in `try_fill_*_slot`. With many spill files that drained tail dominates the
/// scan. Filtering is scheduling-only — an EOF slot can never progress — so it
/// changes no output.
#[test]
fn emptiest_first_order_skips_slots_that_reached_eof() {
    use std::sync::atomic::Ordering;

    let mk = |file_id: u32, nblocks: usize, eof: bool| {
        let s = Arc::new(SortMergeSlot::new(
            file_id,
            BufReader::new(tempfile::tempfile().expect("tempfile")),
            SpillCodec::Bgzf,
        ));
        for _ in 0..nblocks {
            s.decompressed.lock().expect("decompressed lock").push_back(vec![0u8]);
        }
        if eof {
            s.queue_eof.store(true, Ordering::Release);
        }
        s
    };

    // Slot 1 is the emptiest but has drained; it must not appear at all.
    let slots = vec![mk(0, 5, false), mk(1, 0, true), mk(2, 3, false)];
    assert_eq!(
        SortSpillDecompress::emptiest_first_order(&slots),
        vec![2, 0],
        "drained slots are skipped; the rest stay most-starved-first"
    );

    // Every slot drained ⇒ nothing to scan.
    let all_done = vec![mk(0, 0, true), mk(1, 0, true)];
    assert!(
        SortSpillDecompress::emptiest_first_order(&all_done).is_empty(),
        "a fully drained registry yields an empty scan order"
    );

    // No slot drained ⇒ unchanged from the pre-filter behaviour.
    let none_done = vec![mk(0, 5, false), mk(1, 1, false), mk(2, 3, false)];
    assert_eq!(SortSpillDecompress::emptiest_first_order(&none_done), vec![1, 2, 0]);
}

/// The `block_batch` clamp is the sibling of the `window_budget` normalization
/// already pinned by `zero_output_byte_limit_normalizes_reorder_window`. A
/// `block_batch` of 0 would declare a phantom EOF after reading nothing on the
/// inline path — silent record loss — so it is clamped to at least one block.
#[test]
fn zero_block_batch_is_clamped_to_one() {
    let tuning = SortDecompressTuning { block_batch: 0, ..Default::default() };
    let clamped = SortSpillDecompress::new(4 * 1024 * 1024, tuning);
    assert_eq!(clamped.tuning.block_batch, 1, "a zero block_batch must be clamped, not honoured");

    // A sane value passes through untouched.
    let tuning = SortDecompressTuning { block_batch: 8, ..Default::default() };
    let passthrough = SortSpillDecompress::new(4 * 1024 * 1024, tuning);
    assert_eq!(passthrough.tuning.block_batch, 8, "a nonzero block_batch is honoured");
}
