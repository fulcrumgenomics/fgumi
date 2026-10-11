use super::*;
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
    let slot = Arc::new(SortMergeSlot::for_test(0, SpillCodec::Zstd));
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
    let slot = Arc::new(SortMergeSlot::for_test(0, SpillCodec::Zstd));
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
    let slot = Arc::new(SortMergeSlot::for_test(0, SpillCodec::Zstd));
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
        let s = Arc::new(SortMergeSlot::for_test(file_id, SpillCodec::Bgzf));
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
        let slot = Arc::new(SortMergeSlot::for_test(0, SpillCodec::Bgzf));
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
        let s = Arc::new(SortMergeSlot::for_test(file_id, SpillCodec::Bgzf));
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

// ── Merge-demand notification ────────────────────────────────────────────────
//
// Every wake below is observed by a real thread parked for 10 s under a 5 s
// watchdog, so only an unpark from the step under test can end the park in time.

/// Write a BGZF spill file holding `blocks` blocks of `per_block` u64-keyed
/// records (prologue, independently compressed blocks, trailer — the layout
/// `SpillWrite` produces) and open it as slot `file_id` with the production
/// opener, `fgumi_sort::open_spill_slot`.
fn spill_slot(
    dir: &std::path::Path,
    file_id: u32,
    blocks: usize,
    per_block: usize,
) -> Arc<SortMergeSlot> {
    use std::io::Write;
    let path = dir.join(format!("run{file_id}.spill"));
    let mut file = std::fs::File::create(&path).expect("create spill file");
    file.write_all(fgumi_sort::spill_magic(SpillCodec::Bgzf)).unwrap();
    let mut compressor = fgumi_sort::SpillBlockCompressor::new(SpillCodec::Bgzf, 1).unwrap();
    for block in 0..blocks {
        let mut raw = Vec::new();
        for i in 0..per_block {
            let key = fgumi_sort::RawCoordinateKey { sort_key: (block * per_block + i) as u64 };
            fgumi_sort::frame_keyed_record_into(&mut raw, &key, &[0u8; 8]).unwrap();
        }
        file.write_all(&compressor.compress_block(&raw).unwrap()).unwrap();
    }
    file.write_all(fgumi_sort::spill_trailer(SpillCodec::Bgzf)).unwrap();
    drop(file);
    fgumi_sort::open_spill_slot(&path, file_id).expect("open spill slot")
}

/// A thread registered as `demand`'s consumer and awaiting `file_id`, then
/// parked for up to 10 s. `assert_woken` fails the test if it has not woken
/// within 5 s — with a 10 s park, only an unpark can end it that early. It
/// registers the way the merge does, through `await_slot`, on an empty open
/// stand-in slot carrying `file_id`.
struct ParkedConsumer {
    done: std::sync::mpsc::Receiver<()>,
    handle: std::thread::JoinHandle<()>,
}

impl ParkedConsumer {
    fn start(demand: &Arc<fgumi_sort::MergeDemand>, file_id: u32) -> Self {
        let d = Arc::clone(demand);
        let (ready_tx, ready_rx) = std::sync::mpsc::channel::<()>();
        let (done_tx, done) = std::sync::mpsc::channel::<()>();
        let handle = std::thread::spawn(move || {
            let stand_in = SortMergeSlot::for_test(file_id, SpillCodec::Bgzf);
            assert!(!d.await_slot(&stand_in), "an empty, open slot never satisfies the re-check");
            ready_tx.send(()).unwrap();
            std::thread::park_timeout(std::time::Duration::from_secs(10));
            let _ = done_tx.send(());
        });
        ready_rx.recv().unwrap();
        Self { done, handle }
    }

    fn assert_woken(self, what: &str) {
        self.done
            .recv_timeout(std::time::Duration::from_secs(5))
            .unwrap_or_else(|_| panic!("{what}: the awaiting consumer was not woken"));
        self.handle.join().unwrap();
    }

    /// Fails if the consumer wakes within 500 ms. A wrong unpark ends the park
    /// at once, so the window only bounds how quickly a woken thread gets to
    /// run; a correct step can never fail this.
    fn assert_still_parked(&self, what: &str) {
        assert!(
            self.done.recv_timeout(std::time::Duration::from_millis(500)).is_err(),
            "{what}: the consumer woke without a matching delivery"
        );
    }
}

/// Register `slots` with the step, as its `SpillReady` arm does.
fn register_slots(step: &mut SortSpillDecompress, slots: &[Arc<SortMergeSlot>]) {
    let mut registry = step.registry.lock();
    for slot in slots {
        registry.push(RegisteredSpill { slot: Arc::clone(slot) });
    }
}

/// The inline (file-granularity) path never raises a slot's `in_flight`, so
/// the merge cannot classify its stalls: wiring the demand there marks the
/// `--sort-stats` awaited-slot line as unclassified; the block-parallel path
/// leaves it classified.
#[rstest::rstest]
#[case::block_parallel(false)]
#[case::inline(true)]
fn inline_decompress_marks_stalls_unclassified(#[case] file_granularity: bool) {
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let _step = SortSpillDecompress::new(
        1 << 20,
        SortDecompressTuning { file_granularity, block_batch: 1 },
    )
    .with_merge_demand(Arc::clone(&demand));
    assert_eq!(demand.snapshot().inline_decompress, file_granularity);
    assert_eq!(
        demand.snapshot().log_lines()[1].contains("not classified"),
        file_granularity,
        "{:?}",
        demand.snapshot().log_lines()
    );
}

/// A delivery to the awaited slot wakes the parked merge; a delivery to another
/// slot does not. Both fill paths. Only slot 1 is registered with the step, so
/// emptiest-first scheduling cannot pick another slot.
#[rstest::rstest]
#[case::block_parallel(false)]
#[case::inline(true)]
fn delivery_to_the_awaited_slot_wakes_the_merge(#[case] file_granularity: bool) {
    let tmp = tempfile::tempdir().unwrap();
    let slot = spill_slot(tmp.path(), 1, 4, 8);
    let demand = Arc::new(fgumi_sort::MergeDemand::new());

    // A delivery to slot 1 while the consumer awaits file 9: no wake. This
    // consumer is never joined; it times out on its own.
    let other = ParkedConsumer::start(&demand, 9);
    let mut step = SortSpillDecompress::new(
        1 << 20,
        SortDecompressTuning { file_granularity, block_batch: 1 },
    )
    .with_merge_demand(Arc::clone(&demand));
    register_slots(&mut step, std::slice::from_ref(&slot));
    assert!(step.try_fill_some_slot().unwrap());
    assert!(slot.fifo_len() > 0, "precondition: the fill delivered a block");
    other.assert_still_parked("delivery to a non-awaited file");
    assert_eq!(demand.awaited(), Some(9));
    demand.clear_awaited();

    // Await slot 1 itself, after draining what was delivered: the next delivery wakes.
    slot.decompressed.lock().unwrap().clear();
    let consumer = ParkedConsumer::start(&demand, 1);
    while slot.fifo_len() == 0 {
        assert!(step.try_fill_some_slot().unwrap());
    }
    consumer.assert_woken("delivery to the awaited slot");
    assert_eq!(demand.awaited(), None);
}

/// End of file with no block to deliver wakes the awaiting merge. The reads
/// are positional, so EOF is known when the cursor reaches the end with
/// nothing carried or pending — normally in the same read as the last blocks.
/// A spill whose body holds no frame (only the BGZF EOF marker) is the case
/// where the EOF-finalizing read delivers nothing (inline: `got == 0`;
/// block-parallel: `bp_commit_read(0, true)` → an empty
/// `bp_insert_drain_finalize`), and that path must still wake the merge.
#[rstest::rstest]
#[case::block_parallel(false)]
#[case::inline(true)]
fn eof_finalize_wakes_the_awaiting_merge(#[case] file_granularity: bool) {
    let tmp = tempfile::tempdir().unwrap();
    let slot = spill_slot(tmp.path(), 0, 0, 8);
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let mut step = SortSpillDecompress::new(
        1 << 20,
        SortDecompressTuning { file_granularity, block_batch: 4 },
    )
    .with_merge_demand(Arc::clone(&demand));
    register_slots(&mut step, std::slice::from_ref(&slot));
    let consumer = ParkedConsumer::start(&demand, 0);
    while !slot.queue_eof.load(Ordering::Acquire) {
        step.try_fill_some_slot().unwrap();
    }
    consumer.assert_woken("EOF finalize");
    assert_eq!(slot.fifo_len(), 0, "EOF with no block");
}

/// A run whose last read delivers its final blocks together with EOF: the
/// first fill (four of eight blocks) must not see EOF — frames remain
/// pending — and the second delivers the rest and finalizes.
#[rstest::rstest]
#[case::block_parallel(false)]
#[case::inline(true)]
fn eof_is_reached_by_position_not_by_count(#[case] file_granularity: bool) {
    let tmp = tempfile::tempdir().unwrap();
    let slot = spill_slot(tmp.path(), 0, 8, 8);
    let mut step = SortSpillDecompress::new(
        1 << 20,
        SortDecompressTuning { file_granularity, block_batch: 4 },
    );
    register_slots(&mut step, std::slice::from_ref(&slot));
    while slot.fifo_len() < 4 {
        assert!(step.try_fill_some_slot().unwrap());
    }
    assert!(!slot.queue_eof.load(Ordering::Acquire), "frames are still pending");
    slot.decompressed.lock().unwrap().clear();
    while !slot.queue_eof.load(Ordering::Acquire) {
        step.try_fill_some_slot().unwrap();
    }
    assert_eq!(slot.fifo_len(), 4, "the last four blocks arrive with EOF");
}

/// The window reader pins no read buffer between fills: after a fill that
/// leaves frames pending, every window buffer is back in the pool (nothing
/// leased), and the pending frames are owned copies. Both fill paths.
#[rstest::rstest]
#[case::block_parallel(false)]
#[case::inline(true)]
fn window_reader_pins_no_read_buffer_between_fills(#[case] file_granularity: bool) {
    let tmp = tempfile::tempdir().unwrap();
    let slot = spill_slot(tmp.path(), 0, 8, 8);
    let mut step = SortSpillDecompress::new(
        1 << 20,
        SortDecompressTuning { file_granularity, block_batch: 1 },
    );
    register_slots(&mut step, std::slice::from_ref(&slot));
    assert!(step.try_fill_some_slot().unwrap());
    let reader = slot.reader.lock().unwrap();
    assert!(!reader.pending.is_empty(), "precondition: the window held more than one frame");
    assert!(
        reader.pending.iter().all(|f| matches!(f, RawFrame::Owned(_))),
        "pending frames are owned"
    );
    drop(reader);
    assert_eq!(step.slices.resident_bytes(), 0, "no window buffer is leased between fills");
}

/// A slot failure (poisoned reader) wakes the awaiting merge, so it surfaces
/// the error instead of sleeping until its timer.
#[rstest::rstest]
#[case::block_parallel(false)]
#[case::inline(true)]
fn failed_slot_wakes_the_awaiting_merge(#[case] file_granularity: bool) {
    let slot = Arc::new(SortMergeSlot::for_test(0, SpillCodec::Bgzf));
    let holder = Arc::clone(&slot);
    let _ = std::thread::spawn(move || {
        let _g = holder.reader.lock().unwrap();
        panic!("simulated fill-worker panic under the reader lock");
    })
    .join();
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let consumer = ParkedConsumer::start(&demand, 0);
    let mut step = SortSpillDecompress::new(
        1 << 20,
        SortDecompressTuning { file_granularity, block_batch: 1 },
    )
    .with_merge_demand(Arc::clone(&demand));
    let r = if file_granularity {
        step.try_fill_inline_slot(&slot)
    } else {
        step.try_fill_block_parallel_slot(&slot)
    };
    assert!(r.is_err());
    assert!(slot.has_error());
    consumer.assert_woken("slot failure");
}

/// Worker copies share the one demand (`Arc::ptr_eq`).
#[test]
fn worker_copies_share_the_merge_demand() {
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let step = SortSpillDecompress::new(1 << 20, SortDecompressTuning::default())
        .with_merge_demand(Arc::clone(&demand));
    let copy = step.new_worker_copy();
    assert!(Arc::ptr_eq(copy.merge_demand_for_test().unwrap(), &demand));
}

/// The block-parallel drain-only pass (Phase B) delivers a block that waited in
/// the reorder buffer behind a full FIFO; that delivery must wake the merge too.
/// The test holds the reader lock so the fill cannot take Phase A.
#[test]
fn phase_b_drain_wakes_the_awaiting_merge() {
    let slot = Arc::new(SortMergeSlot::for_test(0, SpillCodec::Bgzf));
    // A full FIFO, and one in-order block parked in the reorder buffer.
    {
        let mut dec = slot.decompressed.lock().unwrap();
        for _ in 0..fgumi_sort::PHASE2_DECOMP_CAP {
            dec.push_back(vec![0]);
        }
    }
    slot.bp_commit_read(1, false);
    assert!(!slot.bp_insert_drain_finalize(0, vec![vec![7]], 1), "precondition: no FIFO room");
    slot.decompressed.lock().unwrap().clear(); // the merge consumed the FIFO
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let consumer = ParkedConsumer::start(&demand, 0);
    let mut step = SortSpillDecompress::new(
        1 << 20,
        SortDecompressTuning { file_granularity: false, block_batch: 1 },
    )
    .with_merge_demand(Arc::clone(&demand));
    let reader = slot.reader.lock().unwrap();
    assert!(step.try_fill_block_parallel_slot(&slot).unwrap(), "Phase B drained the block");
    drop(reader);
    consumer.assert_woken("Phase-B drain");
    assert_eq!(slot.fifo_len(), 1);
}
