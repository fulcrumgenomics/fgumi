use std::sync::Arc;
use std::time::Duration;

use fgumi_bam_io::pread::SliceBufferPool;
use fgumi_pipeline_core::testing::StepProbe;
use fgumi_pipeline_core::{PhaseCap, StepOutcome};
use fgumi_sort::{SortMergeSlot, SpillCodec};

use super::*;
use crate::pread::{ReadSlice, ReadTarget, SpillClass};

/// Write a BGZF spill file of `blocks` frames of `per_block` keyed records
/// (the format `SpillWrite` produces) and open it as slot `file_id` with the
/// production opener.
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

/// The slot's file body cut at `cuts` (slice lengths; the last slice takes the
/// rest) into leased slices with dense per-slot sequences, each recorded as
/// issued (the planner's half of `issued_bytes`).
fn slices(slot: &Arc<SortMergeSlot>, cuts: &[usize]) -> Vec<ReadSlice> {
    use std::os::unix::fs::FileExt;
    let pool = SliceBufferPool::new(0);
    let start = slot.body_start();
    let mut body = vec![0u8; usize::try_from(slot.len() - start).unwrap()];
    slot.source().read_exact_at(&mut body, start).unwrap();
    let mut sizes = Vec::new();
    let mut left = body.len();
    for &c in cuts {
        let n = c.min(left);
        if n == 0 {
            break;
        }
        sizes.push(n);
        left -= n;
    }
    if left > 0 || sizes.is_empty() {
        sizes.push(left);
    }
    let mut out = Vec::new();
    let mut at = 0usize;
    for (seq, n) in sizes.into_iter().enumerate() {
        let bytes = body[at..at + n].to_vec();
        slot.bp_note_issued(n as u64);
        out.push(ReadSlice {
            ordinal: seq as u64,
            stream: slot.file_id,
            seq: u32::try_from(seq).unwrap(),
            offset: start + at as u64,
            bytes: pool.lease(bytes),
            last: at + n == body.len(),
            target: ReadTarget::Spill { slot: Arc::clone(slot), class: SpillClass::Cold },
        });
        at += n;
    }
    out
}

/// A supply around `demand` with a fresh ledger.
fn supply_with(demand: &Arc<fgumi_sort::MergeDemand>) -> SpillSupply {
    SpillSupply {
        demand: Arc::clone(demand),
        ledger: Arc::new(SupplyLedger::default()),
        retire_epoch: Arc::default(),
    }
}

/// Pop every queued block.
fn drain_fifo(slot: &SortMergeSlot) -> Vec<Vec<u8>> {
    std::iter::from_fn(|| slot.pop_decompressed()).collect()
}

/// Run the step (consuming the FIFO as the merge would) until `slot` is
/// `queue_eof`; returns every block delivered, in order.
fn run_to_eof(
    step: &mut SortSpillDecompress,
    probe: &StepProbe<SortSpillDecompress>,
    slot: &SortMergeSlot,
) -> Vec<Vec<u8>> {
    let mut got = Vec::new();
    for _ in 0..10_000 {
        probe.try_run(step).expect("try_run");
        got.extend(drain_fifo(slot));
        if slot.queue_eof() {
            got.extend(drain_fifo(slot));
            return got;
        }
    }
    panic!("the slot never reached EOF");
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
            std::thread::park_timeout(Duration::from_secs(10));
            let _ = done_tx.send(());
        });
        ready_rx.recv().unwrap();
        Self { done, handle }
    }

    fn assert_woken(self, what: &str) {
        self.done
            .recv_timeout(Duration::from_secs(5))
            .unwrap_or_else(|_| panic!("{what}: the awaiting consumer was not woken"));
        self.handle.join().unwrap();
    }
}

/// The scan visits the merge's demand slots first — awaited, predicted,
/// frontier, deduplicated — then the rotation, and drops a demand slot that
/// is already EOF.
#[test]
fn scan_order_puts_demand_slots_first_and_drops_eof() {
    let slots: Vec<_> =
        (0..5).map(|i| Arc::new(SortMergeSlot::for_test(i, SpillCodec::Bgzf))).collect();
    for s in &slots[..4] {
        s.bp_stash_frames_for_test(vec![vec![1, 2, 3]], false);
    }
    slots[4].bp_stash_frames_for_test(Vec::new(), true);
    assert!(slots[4].queue_eof());
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let step = SortSpillDecompress::new(1 << 20, &supply_with(&demand));
    for s in &slots {
        step.register_for_test(s);
    }
    await_file(&demand, 3);
    demand.set_predicted(Some(1));
    demand.set_frontier(Some(4));
    let mut order = Vec::new();
    let _ = step.scan_order(&mut order);
    assert_eq!(order[..2], [3, 1], "awaited, predicted; the EOF frontier dropped: {order:?}");
    assert_eq!(order.len(), 4, "the rest once each, no EOF slot: {order:?}");
    demand.set_frontier(Some(1));
    demand.set_predicted(Some(0));
    let _ = step.scan_order(&mut order);
    assert_eq!(order[..3], [3, 0, 1], "{order:?}");
    assert_eq!(order.len(), 4, "{order:?}");
}

/// Slices of one slot that land in reverse order are parsed in sequence: the
/// delivered blocks equal those of the same file ingested as one slice.
#[test]
fn out_of_order_slices_for_one_slot_are_reassembled() {
    let tmp = tempfile::tempdir().unwrap();
    let whole = spill_slot(tmp.path(), 0, 8, 16);
    let cut = spill_slot(tmp.path(), 1, 8, 16);
    let mut step = SortSpillDecompress::new(1 << 20, &SpillSupply::new());
    let probe = StepProbe::new(&step);
    for s in slices(&whole, &[]) {
        probe.push_input(s);
    }
    let want = run_to_eof(&mut step, &probe, &whole);
    assert_eq!(want.len(), 8);
    let mut step = SortSpillDecompress::new(1 << 20, &SpillSupply::new());
    let probe = StepProbe::new(&step);
    let mut parts = slices(&cut, &[7, 300, 1, 2000, 18]);
    parts.reverse();
    for s in parts {
        probe.push_input(s);
    }
    assert_eq!(run_to_eof(&mut step, &probe, &cut), want);
}

/// The first slice of a file the step has not seen registers its slot (the
/// planner, not the step, sees `SpillReady`).
#[test]
fn a_slice_for_an_unregistered_file_registers_it() {
    let tmp = tempfile::tempdir().unwrap();
    let slot = spill_slot(tmp.path(), 4, 3, 4);
    let mut step = SortSpillDecompress::new(1 << 20, &SpillSupply::new());
    let probe = StepProbe::new(&step);
    assert!(step.registry.read().slots.is_empty());
    probe.push_input(slices(&slot, &[]).remove(0));
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Progress);
    let r = step.registry.read();
    assert_eq!(r.slots.len(), 1);
    assert!(Arc::ptr_eq(&r.slots[0], &slot));
    assert_eq!(r.by_slot_id.get(&4), Some(&0));
}

/// Declare `file_id` awaited, as the merge does on a stall: through
/// `await_slot` on an empty stand-in slot (`set_awaited` is crate private).
fn await_file(demand: &fgumi_sort::MergeDemand, file_id: u32) {
    let stand_in = SortMergeSlot::for_test(file_id, SpillCodec::Bgzf);
    assert!(!demand.await_slot(&stand_in), "an empty stand-in never resolves");
}

/// The scan reads only lock-free mirrors: with every lock of every slot held
/// by another thread, `scan_order` still returns (a lock would block it past
/// the 5 s watchdog).
#[test]
fn scan_order_takes_zero_slot_locks() {
    let slots: Vec<_> =
        (0..3).map(|i| Arc::new(SortMergeSlot::for_test(i, SpillCodec::Bgzf))).collect();
    for s in &slots {
        s.bp_stash_frames_for_test(vec![vec![1, 2, 3]], false);
    }
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    await_file(&demand, 2);
    let step = SortSpillDecompress::new(1 << 20, &supply_with(&demand));
    for s in &slots {
        step.register_for_test(s);
    }
    let (held_tx, held_rx) = std::sync::mpsc::channel();
    let (release_tx, release_rx) = std::sync::mpsc::channel::<()>();
    let holder_slots = slots.clone();
    let holder = std::thread::spawn(move || {
        let guards: Vec<_> = holder_slots.iter().map(|s| s.lock_all_for_test()).collect();
        held_tx.send(()).unwrap();
        release_rx.recv().unwrap();
        drop(guards);
    });
    held_rx.recv().unwrap();
    let (tx, rx) = std::sync::mpsc::channel();
    let scanner = std::thread::spawn(move || {
        let mut order = Vec::new();
        let scanned = step.scan_order(&mut order);
        tx.send((order, scanned)).unwrap();
    });
    let (order, (_, any_work)) =
        rx.recv_timeout(Duration::from_secs(5)).expect("scan_order blocked on a slot lock");
    assert!(any_work);
    release_tx.send(()).unwrap();
    holder.join().unwrap();
    scanner.join().unwrap();
    assert_eq!(order[0], 2, "the awaited slot is first: {order:?}");
    assert_eq!(order.len(), 3, "every slot with a claimable head: {order:?}");
}

/// A malformed slice fails the slot (the merge sees the error, not EOF), fails
/// the step, and wakes a merge awaiting the slot.
#[test]
fn ingest_error_marks_the_slot_and_fails_the_step() {
    let slot = Arc::new(SortMergeSlot::for_test(0, SpillCodec::Bgzf));
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let consumer = ParkedConsumer::start(&demand, 0);
    let mut step = SortSpillDecompress::new(1 << 20, &supply_with(&demand));
    let probe = StepProbe::new(&step);
    slot.bp_note_issued(64);
    probe.push_input(ReadSlice {
        ordinal: 0,
        stream: 0,
        seq: 0,
        offset: 0,
        bytes: SliceBufferPool::new(0).lease(vec![0xAB; 64]),
        last: true,
        target: ReadTarget::Spill { slot: Arc::clone(&slot), class: SpillClass::Cold },
    });
    assert!(probe.try_run(&mut step).is_err());
    assert!(slot.has_error());
    assert!(!slot.is_drained(), "a failed slot is never a clean drain");
    consumer.assert_woken("ingest failure");
}

/// Each serve claims, decompresses and publishes exactly one block, and a
/// delivery to the awaited slot wakes the parked merge.
#[test]
fn one_block_per_claim_and_notify_on_delivery() {
    let tmp = tempfile::tempdir().unwrap();
    let slot = spill_slot(tmp.path(), 0, 4, 8);
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let mut step = SortSpillDecompress::new(1 << 20, &supply_with(&demand));
    let probe = StepProbe::new(&step);
    probe.push_input(slices(&slot, &[]).remove(0));
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Progress);
    assert_eq!(slot.fifo_len(), 1, "the ingest serves its slot once");
    assert_eq!(slot.stash_len_relaxed(), 3);
    drain_fifo(&slot);
    let consumer = ParkedConsumer::start(&demand, 0);
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Progress);
    consumer.assert_woken("delivery to the awaited slot");
    assert_eq!(slot.fifo_len(), 1, "one block per claim");
    for n in [2, 3] {
        assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Progress);
        assert_eq!(slot.fifo_len(), n);
    }
    assert!(slot.queue_eof(), "the last claim finalizes");
}

/// The ledger's stash tally releases on a claim exactly what the ingest
/// charged — the slice a borrowed frame pins, an owned frame's allocation and
/// the carry — so a slot read to EOF leaves it at zero, and its peak stays
/// within the bytes the slot could hold: the whole file ingested as one slice
/// peaks at that slice; cut into slices with frames straddling the cuts, at
/// most every slice resident plus an owned copy of every frame.
#[test]
fn ledger_stash_returns_to_zero_and_peaks_within_the_resident_bound() {
    let tmp = tempfile::tempdir().unwrap();
    let cuts: [&[usize]; 2] = [&[], &[7, 300, 1, 2000, 18]];
    for (id, cuts) in cuts.into_iter().enumerate() {
        let slot = spill_slot(tmp.path(), u32::try_from(id).unwrap(), 8, 16);
        let supply = SpillSupply::new();
        let mut step = SortSpillDecompress::new(1 << 20, &supply);
        let probe = StepProbe::new(&step);
        let parts = slices(&slot, cuts);
        let body: u64 = parts.iter().map(|s| s.bytes.len() as u64).sum();
        let resident: u64 = parts.iter().map(|s| s.bytes.capacity() as u64).sum();
        let whole = parts.len() == 1;
        for s in parts {
            probe.push_input(s);
        }
        assert_eq!(run_to_eof(&mut step, &probe, &slot).len(), 8);
        let l = &supply.ledger;
        assert_eq!(l.claims(), (8, 0), "cuts {cuts:?}");
        assert_eq!(l.stash_bytes(), 0, "every charge released (cuts {cuts:?})");
        let bound = if whole { resident } else { resident + body };
        assert!(l.stash_peak_bytes() > 0, "cuts {cuts:?}");
        assert!(
            l.stash_peak_bytes() <= bound,
            "peak {} over bound {bound} (cuts {cuts:?})",
            l.stash_peak_bytes()
        );
    }
}

/// The step finishes only once its input is drained and every registered
/// slot is `queue_eof`; until then an idle poll is `NoProgress`.
#[test]
fn finishes_only_when_input_drained_and_every_slot_eof() {
    let a = Arc::new(SortMergeSlot::for_test(0, SpillCodec::Bgzf));
    let b = Arc::new(SortMergeSlot::for_test(1, SpillCodec::Bgzf));
    a.bp_stash_frames_for_test(Vec::new(), true);
    assert!(a.queue_eof());
    let mut step = SortSpillDecompress::new(1 << 20, &SpillSupply::new());
    step.register_for_test(&a);
    step.register_for_test(&b);
    let probe = StepProbe::new(&step);
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::NoProgress);
    probe.close_input();
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::NoProgress, "b is not EOF");
    b.bp_stash_frames_for_test(Vec::new(), true);
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Finished);
}

/// A slot whose last block waits in `reorder` behind a full
/// FIFO, with an empty stash, is still work — the poll takes the permit (it is
/// refused while the only permit is held elsewhere), drains the block once the
/// consumer made room, and the slot reaches `queue_eof`.
#[test]
fn scan_drains_a_slot_with_an_empty_stash() {
    let cap_blocks = u8::try_from(fgumi_sort::PHASE2_DECOMP_CAP).unwrap();
    let slot = Arc::new(SortMergeSlot::for_test(0, SpillCodec::Bgzf));
    slot.bp_stash_frames_for_test((0..=cap_blocks).map(|i| vec![i]).collect(), true);
    while let Some(b) = slot.bp_claim_raw(u64::MAX) {
        slot.bp_insert_drain_finalize(b.seq, vec![b.frame.to_vec()], 1);
    }
    assert_eq!(slot.stash_len_relaxed(), 0);
    assert_eq!(slot.reorder_len_relaxed(), 1, "the last block waits behind the full FIFO");
    assert!(!slot.queue_eof());
    assert_eq!(slot.pop_decompressed(), Some(vec![0]), "the consumer makes room");
    let cap = PhaseCap::new("sort-phase2", 1);
    let mut step = SortSpillDecompress::new(1 << 20, &SpillSupply::new())
        .with_step_cap(fgumi_pipeline_core::StepCap::from(Arc::clone(&cap)));
    step.register_for_test(&slot);
    let probe = StepProbe::new(&step);
    let held = cap.try_acquire().expect("the only permit");
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Capped, "a drainable slot is work");
    drop(held);
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Progress);
    assert!(slot.queue_eof());
    assert_eq!(drain_fifo(&slot).len(), usize::from(cap_blocks));
}

/// A stash head behind a front that waits in `reorder` (the FIFO full) is not
/// work: no claim can admit it and no drain can move it, so the poll takes no
/// permit and reports `NoProgress` rather than spinning on `Capped` or
/// `Contention`.
#[test]
fn a_head_behind_a_waiting_front_is_not_work() {
    let cap_blocks = u8::try_from(fgumi_sort::PHASE2_DECOMP_CAP).unwrap();
    let slot = Arc::new(SortMergeSlot::for_test(0, SpillCodec::Bgzf));
    slot.bp_stash_frames_for_test((0..cap_blocks + 2).map(|i| vec![i]).collect(), true);
    for _ in 0..=cap_blocks {
        let b = slot.bp_claim_raw(u64::MAX).expect("the front, or room");
        slot.bp_insert_drain_finalize(b.seq, vec![b.frame.to_vec()], 1);
    }
    assert_eq!(slot.reorder_len_relaxed(), 1, "the front waits behind the full FIFO");
    assert_eq!(slot.stash_len_relaxed(), 1);
    let cap = PhaseCap::new("sort-phase2", 1);
    let mut step = SortSpillDecompress::new(1 << 20, &SpillSupply::new())
        .with_step_cap(fgumi_pipeline_core::StepCap::from(Arc::clone(&cap)));
    step.register_for_test(&slot);
    let probe = StepProbe::new(&step);
    let held = cap.try_acquire().expect("the only permit");
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::NoProgress);
    assert_eq!(cap.refused(), 0, "no admission attempt");
    drop(held);
    assert_eq!(slot.pop_decompressed(), Some(vec![0]), "the consumer makes room");
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Progress);
}

/// The ingest wake contract: a last slice that carries only
/// the BGZF EOF marker finalizes the slot with no block to deliver, and the
/// step wakes the merge awaiting it.
#[test]
fn eof_finalize_wakes_the_awaiting_merge() {
    let tmp = tempfile::tempdir().unwrap();
    let slot = spill_slot(tmp.path(), 0, 0, 8);
    assert!(!slot.is_empty(), "the body holds the EOF marker");
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let consumer = ParkedConsumer::start(&demand, 0);
    let mut step = SortSpillDecompress::new(1 << 20, &supply_with(&demand));
    let probe = StepProbe::new(&step);
    probe.push_input(slices(&slot, &[]).remove(0));
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Progress);
    assert!(slot.queue_eof());
    assert_eq!(slot.fifo_len(), 0, "EOF with no block");
    consumer.assert_woken("EOF finalize");
}

/// With no slice and no slot with work the step takes no permit: an idle poll
/// is `NoProgress` (not `Capped`) and counts no refusal, so idle clones never
/// crowd the output compressor out of the shared cap.
#[test]
fn idle_poll_takes_no_permit() {
    let cap = PhaseCap::new("sort-phase2", 1);
    let mut step = SortSpillDecompress::new(1 << 20, &SpillSupply::new())
        .with_step_cap(fgumi_pipeline_core::StepCap::from(Arc::clone(&cap)));
    let probe = StepProbe::new(&step);
    let held = cap.try_acquire().expect("the only permit");
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::NoProgress);
    probe.close_input();
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Finished);
    assert_eq!(cap.refused(), 0);
    drop(held);
}

/// A refused admission leaves the slice queued.
#[test]
fn refused_permit_leaves_the_slice_queued() {
    let tmp = tempfile::tempdir().unwrap();
    let slot = spill_slot(tmp.path(), 0, 2, 4);
    let cap = PhaseCap::new("sort-phase2", 1);
    let mut step = SortSpillDecompress::new(1 << 20, &SpillSupply::new())
        .with_step_cap(fgumi_pipeline_core::StepCap::from(Arc::clone(&cap)));
    let probe = StepProbe::new(&step);
    probe.push_input(slices(&slot, &[]).remove(0));
    let held = cap.try_acquire().expect("the only permit");
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Capped);
    assert!(!probe.input_is_empty(), "the slice is still queued");
    drop(held);
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Progress);
    assert_eq!(cap.active(), 0, "the permit is released when try_run returns");
}

/// Worker copies share the demand, the cap, the registry and the rotation.
#[test]
fn worker_copies_share_demand_cap_and_registry() {
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let cap = PhaseCap::new("sort-phase2", 2);
    let step = SortSpillDecompress::new(1 << 20, &supply_with(&demand))
        .with_step_cap(fgumi_pipeline_core::StepCap::from(Arc::clone(&cap)));
    let copy = step.new_worker_copy();
    assert!(Arc::ptr_eq(copy.merge_demand_for_test(), &demand));
    assert!(std::ptr::eq(
        copy.step_cap().and_then(fgumi_pipeline_core::StepCapRef::active).unwrap(),
        &raw const *cap
    ));
    assert!(Arc::ptr_eq(&copy.registry, &step.registry));
    assert!(Arc::ptr_eq(&copy.cursor, &step.cursor));
}

/// A zero output byte limit still bounds the claim window (zero would read as
/// "no bound").
#[test]
fn zero_output_byte_limit_normalizes_reorder_window() {
    assert_eq!(
        SortSpillDecompress::new(0, &SpillSupply::new()).window_budget,
        DEFAULT_REORDER_WINDOW_BYTES
    );
    assert_eq!(SortSpillDecompress::new(4096, &SpillSupply::new()).window_budget, 4096);
}

/// Design issue 8: the consolidation demand is notified only before
/// `AllAnnounced`. `all_announced` is the shared flag the planner sets on
/// `AllAnnounced`. Catches: phase-2 deliveries still paying the consolidation
/// demand's CAS and unpark.
#[test]
fn decompress_notifies_the_consolidation_demand_only_before_all_announced() {
    use std::sync::atomic::{AtomicBool, Ordering};
    let dir = tempfile::TempDir::new().unwrap();
    let flag = Arc::new(AtomicBool::new(false));
    let cons = Arc::new(fgumi_sort::MergeDemand::new());
    let mut step = SortSpillDecompress::new(1 << 20, &SpillSupply::new())
        .with_consolidation_demand(Arc::clone(&cons), Arc::clone(&flag));
    let probe = StepProbe::new(&step);
    let pool = SliceBufferPool::new(4);
    let deliver = |id: u32, step: &mut SortSpillDecompress| {
        let path = dir.path().join(format!("r{id}"));
        crate::sort::spill_fixtures::write_run_file(&path, SpillCodec::Zstd, id, 0..10);
        let slot = fgumi_sort::open_spill_slot(&path, id).unwrap();
        let body =
            std::fs::read(&path).unwrap()[usize::try_from(slot.body_start()).unwrap()..].to_vec();
        slot.bp_note_issued(body.len() as u64);
        // Registers this thread and publishes the slot as awaited (the slot is
        // empty, so the re-check cannot resolve it).
        assert!(!cons.await_slot(&slot));
        probe.push_input(ReadSlice {
            ordinal: 0,
            stream: 0,
            seq: 0,
            offset: slot.body_start(),
            bytes: pool.lease(body),
            last: true,
            target: ReadTarget::Spill {
                slot: Arc::clone(&slot),
                class: SpillClass::ConsolidationCold,
            },
        });
        assert_eq!(probe.try_run(step).unwrap(), StepOutcome::Progress);
        cons.awaited()
    };
    assert_eq!(
        deliver(0, &mut step),
        None,
        "before AllAnnounced: the delivery clears the awaited slot"
    );
    flag.store(true, Ordering::Release);
    let wakes = cons.snapshot().wakes_delivered;
    assert_eq!(deliver(1, &mut step), Some(1), "after AllAnnounced: not notified");
    assert_eq!(cons.snapshot().wakes_delivered, wakes);
}

/// Retired slots leave the registry (it does not grow with every run ever
/// consolidated), and a scan taken before the compaction resolves no slot, so
/// another clone's compaction cannot redirect it to a different slot. Catches:
/// tombstones that are never removed, and an index from a stale generation
/// resolving after the vector shifted.
#[test]
fn retired_slots_are_compacted_and_a_stale_scan_resolves_nothing() {
    let slots: Vec<_> =
        (0..3).map(|i| Arc::new(SortMergeSlot::for_test(i, SpillCodec::Bgzf))).collect();
    for s in &slots {
        s.bp_stash_frames_for_test(vec![vec![1, 2, 3]], false);
    }
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let supply = supply_with(&demand);
    let mut step = SortSpillDecompress::new(1 << 20, &supply);
    for s in &slots {
        step.register_for_test(s);
    }
    let mut order = Vec::new();
    let (generation, _) = step.scan_order(&mut order);
    assert_eq!(order.len(), 3);
    assert!(step.slot_at(generation, 2).is_some());

    slots[0].retire();
    supply.retire_epoch.fetch_add(1, std::sync::atomic::Ordering::Release);
    step.drop_retired();
    assert_eq!(step.registered_slot_ids_for_test(), [1, 2], "the retired slot is removed");
    assert!(step.slot_at(generation, 0).is_none(), "a pre-compaction index resolves nothing");
    let (now, _) = step.scan_order(&mut order);
    assert_ne!(now, generation);
    let ids: Vec<u32> = order.iter().map(|&i| step.slot_at(now, i).unwrap().slot_id()).collect();
    assert_eq!(ids.len(), 2);
    assert!(!ids.contains(&0));
}
