//! Pure planning tests: fake slots over sparse files of a chosen length, one
//! `plan_pass` at a time, no pipeline. Expected sizes are worked from the
//! budget constants, never read back from the planner.

use std::path::PathBuf;
use std::sync::Arc;

use fgumi_bam_io::pread::{ReadStreamsPolicy, SliceBufferPool};
use fgumi_sort::{MergeDemand, SortMergeSlot, SpillCodec};

use super::SpillReadPlanner;
use crate::pread::{ReadRequest, ReadTarget, SpillClass};
use crate::sort::protocol::{SortPhase1Event, SortPhase2Event};
use crate::sort::read_ahead_budget::{
    COLD_INFLIGHT_BYTES, HOT_ALLOWANCE_BYTES, HOT_FILL_BYTES, ReadAheadBudget,
};
use crate::sort::supply_ledger::{SpillSupply, SupplyLedger};

const MIB: u64 = 1 << 20;
/// 12 GiB → R = 512 MiB; at small k the cold allowance is its 2 MiB maximum
/// and the cold fill 1 MiB.
const TOTAL: u64 = 12 << 30;

/// Declare `file_id` awaited, as the merge does on a stall: through
/// `await_slot` on an empty stand-in slot (`set_awaited` is crate private).
fn await_file(demand: &fgumi_sort::MergeDemand, file_id: u32) {
    let stand_in = SortMergeSlot::for_test(file_id, SpillCodec::Bgzf);
    assert!(!demand.await_slot(&stand_in), "an empty stand-in never resolves");
}

/// A slot over a sparse file of `len` bytes (frames from offset 0).
fn slot(id: u32, len: u64) -> Arc<SortMergeSlot> {
    let f = tempfile::tempfile().expect("tempfile");
    f.set_len(len).expect("set_len");
    Arc::new(SortMergeSlot::new(id, Arc::new(f), 0, len, SpillCodec::Zstd))
}

struct Fixture {
    planner: SpillReadPlanner,
    demand: Arc<MergeDemand>,
    ledger: Arc<SupplyLedger>,
    slots: Vec<Arc<SortMergeSlot>>,
}

/// A planner over `lens.len()` registered slots, announced, with `phase2`
/// merge threads and `streams` fixed read streams.
fn fixture(lens: &[u64], phase2: usize, streams: usize) -> Fixture {
    fixture_with(lens, phase2, streams, TOTAL, None)
}

/// [`fixture`] with `total` memory and, when given, the pool replaced by
/// `pool` bytes.
fn fixture_with(
    lens: &[u64],
    phase2: usize,
    streams: usize,
    total: u64,
    pool: Option<u64>,
) -> Fixture {
    let supply = SpillSupply::new();
    let (demand, ledger) = (Arc::clone(&supply.demand), Arc::clone(&supply.ledger));
    let mut planner =
        planner(total, phase2, 16, &supply).with_read_streams(ReadStreamsPolicy::fixed(streams));
    if let Some(pool) = pool {
        planner = planner.with_pool_for_test(pool);
    }
    let slots: Vec<_> =
        lens.iter().enumerate().map(|(i, &l)| slot(u32::try_from(i).unwrap(), l)).collect();
    for s in &slots {
        let _ = planner.on_event_for_test(spill_ready(s));
    }
    let _ = planner.on_event_for_test(all_announced(lens.len()));
    Fixture { planner, demand, ledger, slots }
}

/// A planner over a fresh, empty slice pool.
fn planner(total: u64, phase2: usize, eligible: usize, supply: &SpillSupply) -> SpillReadPlanner {
    SpillReadPlanner::new(total, phase2, eligible, 1 << 20, supply, SliceBufferPool::new(0))
}

/// The read-ahead class a request is charged to.
fn class(r: &ReadRequest) -> SpillClass {
    match &r.target {
        ReadTarget::Spill { class, .. } => *class,
        ReadTarget::Input => panic!("a spill request"),
    }
}

/// Simulate every read of `slots` landing in its stash, unclaimed: the bytes
/// a slot requested become bytes it holds.
fn land_all(slots: &[Arc<SortMergeSlot>]) {
    for s in slots {
        s.set_read_ahead_for_test(0, s.issued_bytes() + s.stash_bytes());
    }
}

/// Bytes every slot holds (requested + stashed).
fn held(slots: &[Arc<SortMergeSlot>]) -> u64 {
    slots.iter().map(|s| s.issued_bytes() + s.stash_bytes()).sum()
}

fn spill_ready(s: &Arc<SortMergeSlot>) -> SortPhase1Event {
    SortPhase1Event::SpillReady {
        slot: Arc::clone(s),
        path: PathBuf::from(format!("spill-{}", s.file_id)),
        records_ingested_so_far: u64::from(s.file_id),
    }
}

fn all_announced(k: usize) -> SortPhase1Event {
    SortPhase1Event::AllAnnounced {
        slot_count: u32::try_from(k).unwrap(),
        memory_chunk_count: 0,
        total_records: 7,
    }
}

fn for_slot(reqs: &[ReadRequest], id: u32) -> Vec<&ReadRequest> {
    reqs.iter().filter(|r| r.stream == id).collect()
}

fn bytes(reqs: &[&ReadRequest]) -> u64 {
    reqs.iter().map(|r| u64::from(r.len)).sum()
}

/// `(stream, seq, offset, len)` of each request, for assertion messages.
fn shape(reqs: &[ReadRequest]) -> Vec<(u32, u32, u64, u32)> {
    reqs.iter().map(|r| (r.stream, r.seq, r.offset, r.len)).collect()
}

/// The hot set is the awaited slot: it is topped up to 16 MiB in 4 MiB fills;
/// every other slot gets the 2 MiB cold allowance in 1 MiB fills.
#[test]
fn hot_set_is_the_awaited_slot() {
    let mut f = fixture(&[64 * MIB, 64 * MIB, 64 * MIB], 16, 1);
    await_file(&f.demand, 1);
    let reqs = f.planner.plan_pass_for_test();
    let hot = for_slot(&reqs, 1);
    assert_eq!(hot.len(), 4, "{:?}", shape(&reqs));
    assert!(hot.iter().all(|r| u64::from(r.len) == 4 * MIB && class(r) == SpillClass::Hot));
    assert_eq!(bytes(&hot), 16 * MIB);
    for id in [0, 2] {
        let cold = for_slot(&reqs, id);
        assert_eq!(cold.len(), 2, "slot {id}: {:?}", shape(&reqs));
        assert!(cold.iter().all(|r| u64::from(r.len) == MIB && class(r) == SpillClass::Cold));
    }
    assert_eq!(f.slots[1].fifo_cap(), 32);
}

/// A cold slot is never issued past its allowance (`issued + stashed + fill ≤
/// allowance`); once a fill's bytes are parsed and claimed, exactly one more
/// fill is issued.
#[test]
fn cold_slots_never_exceed_their_allowance() {
    let mut f = fixture(&[64 * MIB, 64 * MIB], 16, 1);
    let first = f.planner.plan_pass_for_test();
    assert_eq!(bytes(&for_slot(&first, 0)), 2 * MIB);
    assert!(f.planner.plan_pass_for_test().is_empty(), "both slots are at their allowance");
    // One fill of slot 0 parsed and claimed.
    f.slots[0].set_read_ahead_for_test(f.slots[0].issued_bytes() - MIB, 0);
    f.ledger.land(SpillClass::Cold, MIB);
    let next = f.planner.plan_pass_for_test();
    assert_eq!(next.len(), 1, "{:?}", shape(&next));
    assert_eq!((next[0].stream, u64::from(next[0].len)), (0, MIB));
    // Stashed bytes count against the allowance too.
    f.slots[1].set_read_ahead_for_test(0, 2 * MIB);
    assert!(for_slot(&f.planner.plan_pass_for_test(), 1).is_empty());
}

/// Outstanding cold bytes stop at `COLD_INFLIGHT_BYTES` and outstanding slices
/// at twice the phase-2 threads, whatever the slot count.
#[test]
fn cold_inflight_and_slice_caps_are_respected() {
    let lens = vec![64 * MIB; 20];
    let mut f = fixture(&lens, 16, 1);
    let reqs = f.planner.plan_pass_for_test();
    let cold: u64 =
        reqs.iter().filter(|r| class(r) == SpillClass::Cold).map(|r| u64::from(r.len)).sum();
    assert_eq!(cold, COLD_INFLIGHT_BYTES, "20 slots × 2 MiB would be 40 MiB");
    assert_eq!(f.ledger.inflight_cold_bytes(), COLD_INFLIGHT_BYTES);

    let mut f = fixture(&lens, 2, 1);
    assert_eq!(f.planner.max_inflight_slices(), 4);
    let reqs = f.planner.plan_pass_for_test();
    assert_eq!(reqs.len(), 4, "2 phase-2 threads → 4 outstanding slices");
    assert!(f.planner.plan_pass_for_test().is_empty());
    f.ledger.land(SpillClass::Cold, MIB);
    assert_eq!(f.planner.plan_pass_for_test().len(), 1, "a landed slice frees one");
}

/// Per slot, slice sequences are dense from 0, the slices tile the file body
/// exactly, and only the slice covering the last byte is `last`.
#[rstest::rstest]
#[case::one_stream(1)]
#[case::four_streams(4)]
fn seq_is_dense_per_slot_and_last_is_set_once(#[case] streams: usize) {
    let len = 5 * MIB + MIB / 2 + 3;
    let mut f = fixture(&[len], 16, streams);
    await_file(&f.demand, 0);
    let mut all = Vec::new();
    loop {
        let reqs = f.planner.plan_pass_for_test();
        if reqs.is_empty() {
            break;
        }
        // Simulate every issued byte landing and being claimed.
        f.slots[0].set_read_ahead_for_test(0, 0);
        for r in &reqs {
            f.ledger.land(class(r), u64::from(r.len));
        }
        all.extend(reqs);
    }
    let seqs: Vec<u32> = all.iter().map(|r| r.seq).collect();
    assert_eq!(seqs, (0..u32::try_from(all.len()).unwrap()).collect::<Vec<_>>());
    let mut at = 0;
    for r in &all {
        assert_eq!(r.offset, at, "slices tile the body");
        at += u64::from(r.len);
    }
    assert_eq!(at, len);
    assert_eq!(all.iter().filter(|r| r.last).count(), 1);
    assert!(all.last().unwrap().last);
}

/// A spill file with no frame bytes is finalized at
/// `SpillReady` — `queue_eof` set, the merge woken — and never read.
#[test]
fn empty_spill_file_is_finalized_without_a_read() {
    let mut f = fixture(&[0, 4 * MIB], 16, 1);
    await_file(&f.demand, 0);
    assert!(f.slots[0].queue_eof());
    assert!(f.slots[0].is_drained(), "a clean EOF, no error");
    let reqs = f.planner.plan_pass_for_test();
    assert!(for_slot(&reqs, 0).is_empty(), "{:?}", shape(&reqs));
    assert_eq!(f.ledger.awaited_starved_reads(), 0, "an EOF slot is not starved");
}

/// Phase events pass through unchanged and in order (the planner forwards each
/// on its branch 1 as it pops it).
#[test]
fn events_are_forwarded_in_order_on_branch_1() {
    let mut p = planner(TOTAL, 4, 4, &SpillSupply::new());
    let a = slot(3, MIB);
    let b = slot(1, MIB);
    let mut out = Vec::new();
    for e in [spill_ready(&a), spill_ready(&b), all_announced(2)] {
        out.push(p.on_event_for_test(e));
    }
    match &out[..] {
        [
            SortPhase2Event::SpillReady { slot: s0, records_ingested_so_far: 3, .. },
            SortPhase2Event::SpillReady { slot: s1, records_ingested_so_far: 1, .. },
            SortPhase2Event::AllAnnounced {
                slot_count: 2,
                memory_chunk_count: 0,
                total_records: 7,
            },
        ] => {
            assert!(Arc::ptr_eq(s0, &a) && Arc::ptr_eq(s1, &b));
        }
        _ => panic!("events reordered or rewritten"),
    }
}

/// No read is planned before `AllAnnounced`, and the budget it resolves uses
/// the announced slot count.
#[test]
fn budget_resolves_at_all_announced_from_the_actual_k() {
    let mut p = planner(TOTAL, 4, 4, &SpillSupply::new());
    let slots: Vec<_> = (0..5).map(|i| slot(i, 8 * MIB)).collect();
    for s in &slots {
        let _ = p.on_event_for_test(spill_ready(s));
    }
    assert!(p.budget_for_test().is_none());
    assert!(p.plan_pass_for_test().is_empty(), "nothing is read before the budget exists");
    let _ = p.on_event_for_test(all_announced(5));
    assert_eq!(p.budget_for_test(), Some(ReadAheadBudget::resolve(TOTAL, 5)));
    assert!(!p.plan_pass_for_test().is_empty());
}

/// A pass that finds the awaited slot with nothing issued or stashed counts an
/// awaited-starved read; one that finds its reads outstanding does not.
#[test]
fn awaited_starved_reads_are_counted() {
    let mut f = fixture(&[64 * MIB, 64 * MIB], 16, 1);
    await_file(&f.demand, 0);
    let _ = f.planner.plan_pass_for_test();
    assert_eq!(f.ledger.awaited_starved_reads(), 1);
    let _ = f.planner.plan_pass_for_test();
    assert_eq!(f.ledger.awaited_starved_reads(), 1, "its reads are outstanding");
    assert_eq!(f.slots[0].issued_bytes(), HOT_ALLOWANCE_BYTES, "topped up to the hot allowance");
}

/// The spill planner adopts the shared read-stream policy: it splits each
/// fill by `streams()` and never moves the ratchet.
#[test]
fn spill_planner_adopts_and_never_changes_streams() {
    let policy = ReadStreamsPolicy::auto();
    let supply = SpillSupply::new();
    let demand = Arc::clone(&supply.demand);
    let mut p = planner(TOTAL, 16, 16, &supply).with_read_streams(Arc::clone(&policy));
    let s = slot(0, 64 * MIB);
    let _ = p.on_event_for_test(spill_ready(&s));
    let _ = p.on_event_for_test(all_announced(1));
    await_file(&demand, 0);
    for _ in 0..3 {
        let reqs = p.plan_pass_for_test();
        let per_fill = policy.slices_for(usize::try_from(HOT_FILL_BYTES).unwrap(), 16);
        assert!(reqs.len().is_multiple_of(per_fill.max(1)), "{} slices", reqs.len());
        s.set_read_ahead_for_test(0, 0);
    }
    assert_eq!(policy.streams(), 1, "the planner never observes");
    assert!(policy.history().is_empty());
}

/// Moving the hot set across every slot cannot grow read-ahead past
/// `k × cold + pool`: each slot keeps what it read while hot (nothing is
/// claimed here), and later hot slots are topped up only from what the pool
/// has left. Without the pool every slot would keep a full hot allowance —
/// 64 × 16 MiB here.
#[rstest::rstest]
#[case::k64_r64(64, 1 << 30)]
#[case::k1024_r512(1024, 12 << 30)]
fn hot_set_churn_never_exceeds_the_bound(#[case] k: usize, #[case] total: u64) {
    let lens = vec![256 * MIB; k];
    let mut f = fixture_with(&lens, 16, 1, total, None);
    let budget = f.planner.budget_for_test().unwrap();
    let bound = budget.read_ahead_bound();
    let mut hot_max = 0;
    for i in 0..k {
        await_file(&f.demand, u32::try_from(i).unwrap());
        // Several passes per hot slot, so in-flight caps cannot be what holds
        // the bound; every read lands and stays unclaimed.
        for _ in 0..8 {
            let reqs = f.planner.plan_pass_for_test();
            for r in &reqs {
                f.ledger.land(class(r), u64::from(r.len));
            }
            land_all(&f.slots);
        }
        let s = &f.slots[i];
        hot_max = hot_max.max(s.issued_bytes() + s.stash_bytes());
        assert!(held(&f.slots) <= bound, "k={k} slot {i}: {} > {bound}", held(&f.slots));
    }
    assert!(hot_max > budget.cold_allowance, "hot slots read past their cold terms");
}

/// A full pool still lets the awaited slot read within its cold terms: the
/// merge can always make progress on the slot it waits for.
#[test]
fn a_full_pool_still_feeds_the_awaited_slot_its_cold_terms() {
    let mut f = fixture_with(&[64 * MIB; 4], 16, 1, TOTAL, Some(0));
    await_file(&f.demand, 2);
    let reqs = f.planner.plan_pass_for_test();
    let awaited = for_slot(&reqs, 2);
    assert_eq!(bytes(&awaited), 2 * MIB, "the cold allowance: {:?}", shape(&reqs));
    assert!(awaited.iter().all(|r| u64::from(r.len) == MIB), "in cold fills");
}

/// A slot that held more than its cold allowance while hot reads nothing more
/// once it leaves the hot set, until the merge drains it below the allowance;
/// its charge keeps the pool from granting the next hot slot more than is
/// left.
#[test]
fn a_demoted_slot_keeps_its_charge_until_drained() {
    let pool = 20 * MIB;
    let mut f = fixture_with(&[64 * MIB; 3], 16, 1, TOTAL, Some(pool));
    await_file(&f.demand, 0);
    let _ = f.planner.plan_pass_for_test();
    land_all(&f.slots);
    assert_eq!(f.slots[0].stash_bytes(), 16 * MIB, "the hot allowance (14 MiB from the pool)");
    await_file(&f.demand, 1);
    let reqs = f.planner.plan_pass_for_test();
    assert!(for_slot(&reqs, 0).is_empty(), "the demoted slot is over its cold allowance");
    // Slot 1 already holds its 2 MiB cold allowance (read while cold). With
    // 6 MiB of pool left (20 − 14) it gets one hot fill (6 MiB held, 4 above
    // its cold allowance), not a second (8 above, past the pool).
    assert_eq!(bytes(&for_slot(&reqs, 1)), 4 * MIB, "{:?}", shape(&reqs));
    land_all(&f.slots);
    // The merge drains slot 0 to nothing: its charge is released, and slot 1
    // tops up from the freed pool to 14 MiB (a third fill would pass the
    // 16 MiB hot allowance).
    f.slots[0].set_read_ahead_for_test(0, 0);
    let reqs = f.planner.plan_pass_for_test();
    assert_eq!(bytes(&for_slot(&reqs, 1)), 8 * MIB, "{:?}", shape(&reqs));
}

/// Idle slice buffers count against the pool.
#[test]
fn idle_slice_buffers_are_charged_to_the_pool() {
    let supply = SpillSupply::new();
    let slices = SliceBufferPool::new(8);
    for _ in 0..3 {
        drop(slices.lease(vec![0; 4 << 20]));
    }
    assert_eq!(slices.idle_bytes(), 3 * (4 << 20), "the returned buffers are idle");
    let mut p = SpillReadPlanner::new(TOTAL, 16, 16, 1 << 20, &supply, Arc::clone(&slices))
        .with_pool_for_test(12 * MIB);
    let s = slot(0, 64 * MIB);
    let _ = p.on_event_for_test(spill_ready(&s));
    let _ = p.on_event_for_test(all_announced(1));
    await_file(&supply.demand, 0);
    let reqs = p.plan_pass_for_test();
    assert_eq!(bytes(&reqs.iter().collect::<Vec<_>>()), 2 * MIB, "only the cold terms");
}

/// Claim every admissible stash head of `slot`, decompress it (identity) and
/// publish it, until nothing more is admitted. Returns the claims made.
fn serve_until_refused(slot: &SortMergeSlot) -> usize {
    let mut n = 0;
    while let Some(b) = slot.bp_claim_raw(u64::MAX) {
        let seq = b.seq;
        let block = b.frame.to_vec();
        drop(b);
        slot.bp_insert_drain_finalize(seq, vec![block], 1);
        n += 1;
    }
    n
}

/// The planner gives the awaited slot the hot FIFO cap (32) and every
/// other slot the cold one (8). Serving a 40-block stash never fills a FIFO
/// past its cap (the front may be claimed over it, but its block waits in
/// `reorder`), and a slot that leaves the hot set is demoted to the cold cap
/// on the next pass. Every block still reaches the consumer once, in order.
#[test]
fn cold_slot_fifo_stops_at_8_hot_at_32() {
    let mut f = fixture(&[64 * MIB, 64 * MIB], 16, 1);
    await_file(&f.demand, 1);
    let _ = f.planner.plan_pass_for_test();
    let caps = |f: &Fixture| -> Vec<u32> { f.slots.iter().map(|s| s.fifo_cap()).collect() };
    assert_eq!(caps(&f), vec![8, 32]);
    for (s, cap) in f.slots.iter().zip([8usize, 32]) {
        s.bp_stash_frames_for_test((0..40u8).map(|i| vec![i]).collect(), true);
        serve_until_refused(s);
        assert_eq!(
            s.fifo_len(),
            cap,
            "slot {}: the FIFO fills to its cap and no further",
            s.file_id
        );
        let mut popped = Vec::new();
        loop {
            while let Some(b) = s.pop_decompressed() {
                popped.push(b[0]);
            }
            serve_until_refused(s);
            s.bp_drain_and_finalize();
            assert!(s.fifo_len() <= cap, "slot {}: {} > cap {cap}", s.file_id, s.fifo_len());
            if s.is_drained() {
                break;
            }
        }
        assert_eq!(popped, (0..40u8).collect::<Vec<_>>(), "slot {}", s.file_id);
    }
    // Saturate the outstanding-slice cap so the cold rotation visits no slot:
    // only the demotion itself can lower slot 1's cap.
    f.ledger.add_inflight(SpillClass::Hot, 0, f.planner.max_inflight_slices());
    await_file(&f.demand, 0);
    assert!(f.planner.plan_pass_for_test().is_empty(), "no read: the slice cap binds");
    assert_eq!(caps(&f), vec![32, 8], "the former hot slot is demoted on the next pass");
}

/// A slot registered before the budget exists already carries the cold cap,
/// so no slot holds a hot-sized FIFO without being in the hot set.
#[test]
fn registered_slots_start_at_the_cold_fifo_cap() {
    let mut p = planner(TOTAL, 4, 4, &SpillSupply::new());
    let s = slot(0, 8 * MIB);
    assert_eq!(s.fifo_cap(), 32, "the slot's own default");
    let _ = p.on_event_for_test(spill_ready(&s));
    assert_eq!(s.fifo_cap(), 8);
}

/// The hot FIFO cap is charged to the pool: with no room for its 24 blocks
/// above the cold cap, the awaited slot keeps the cold cap (and still reads
/// its cold terms); a demoted slot that still holds more than 8 decoded
/// blocks stays charged for them, which keeps the next hot slot's FIFO cold.
#[test]
fn the_hot_fifo_cap_is_granted_only_from_the_pool() {
    let fifo_reserve = 24 * (64 << 10);
    let mut f = fixture_with(&[64 * MIB; 2], 16, 1, TOTAL, Some(fifo_reserve - 1));
    await_file(&f.demand, 0);
    let reqs = f.planner.plan_pass_for_test();
    assert_eq!(f.slots[0].fifo_cap(), 8, "no room for the hot FIFO");
    assert_eq!(bytes(&for_slot(&reqs, 0)), 2 * MIB, "the cold terms: {:?}", shape(&reqs));

    let mut f = fixture_with(&[64 * MIB; 2], 16, 1, TOTAL, Some(fifo_reserve));
    await_file(&f.demand, 0);
    let _ = f.planner.plan_pass_for_test();
    assert_eq!(f.slots[0].fifo_cap(), 32, "exactly room for the hot FIFO");
    // Slot 0 fills its hot FIFO, then leaves the hot set holding 32 blocks.
    f.slots[0].bp_stash_frames_for_test((0..40u8).map(|i| vec![i]).collect(), false);
    serve_until_refused(&f.slots[0]);
    assert_eq!(f.slots[0].fifo_len(), 32);
    // Its reads landed long ago (within its cold allowance, not the pool).
    f.slots[0].set_read_ahead_for_test(0, f.slots[0].stash_bytes());
    await_file(&f.demand, 1);
    let _ = f.planner.plan_pass_for_test();
    assert_eq!((f.slots[0].fifo_cap(), f.slots[1].fifo_cap()), (8, 8), "still charged");
    // The merge pops slot 0 down to its cold cap: the charge is released.
    for _ in 0..24 {
        f.slots[0].pop_decompressed().unwrap();
    }
    let _ = f.planner.plan_pass_for_test();
    assert_eq!(f.slots[1].fifo_cap(), 32, "granted from the released charge");
}
