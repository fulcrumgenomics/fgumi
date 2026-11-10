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

/// The hot set is the awaited, predicted and frontier slots, in that order,
/// deduplicated: each is topped up to 16 MiB in 4 MiB fills at the hot FIFO
/// cap; every other slot gets the 2 MiB cold allowance in 1 MiB fills.
#[rstest::rstest]
#[case::three_distinct(Some(1), Some(2), Some(3), &[1, 2, 3])]
#[case::predicted_is_awaited(Some(1), Some(1), Some(3), &[1, 3])]
#[case::frontier_is_predicted(Some(1), Some(2), Some(2), &[1, 2])]
#[case::all_one(Some(4), Some(4), Some(4), &[4])]
#[case::no_awaited(None, Some(2), Some(0), &[2, 0])]
fn hot_set_is_awaited_then_predicted_then_frontier(
    #[case] awaited: Option<u32>,
    #[case] predicted: Option<u32>,
    #[case] frontier: Option<u32>,
    #[case] hot: &[u32],
) {
    let mut f = fixture(&[64 * MIB; 5], 16, 1);
    if let Some(a) = awaited {
        await_file(&f.demand, a);
    }
    f.demand.set_predicted(predicted);
    f.demand.set_frontier(frontier);
    assert_eq!(f.planner.hot_set_for_test(), hot, "deduplicated, at most three");
    let reqs = f.planner.plan_pass_for_test();
    let first: Vec<u32> = {
        let mut seen = Vec::new();
        for r in reqs.iter().filter(|r| class(r) == SpillClass::Hot) {
            if !seen.contains(&r.stream) {
                seen.push(r.stream);
            }
        }
        seen
    };
    assert_eq!(first, hot, "hot slots, in issue order: {:?}", shape(&reqs));
    for id in 0..5u32 {
        let mine = for_slot(&reqs, id);
        let cap = f.slots[id as usize].fifo_cap();
        if hot.contains(&id) {
            assert_eq!(bytes(&mine), HOT_ALLOWANCE_BYTES, "slot {id}");
            assert!(mine.iter().all(|r| u64::from(r.len) == HOT_FILL_BYTES), "slot {id}");
            assert_eq!(cap, 32, "slot {id}");
        } else {
            assert_eq!(bytes(&mine), 2 * MIB, "slot {id}: {:?}", shape(&reqs));
            assert!(mine.iter().all(|r| class(r) == SpillClass::Cold), "slot {id}");
            assert_eq!(cap, 8, "slot {id}");
        }
    }
}

/// A slot whose bytes are all issued, or that reached EOF, is not hot: the
/// next demand slot takes no reads from it, and it holds no hot FIFO cap.
#[test]
fn hot_set_skips_a_slot_at_eof() {
    let mut f = fixture(&[0, 64 * MIB, 64 * MIB], 16, 1);
    await_file(&f.demand, 0);
    f.demand.set_predicted(Some(1));
    let reqs = f.planner.plan_pass_for_test();
    assert_eq!(f.slots[0].fifo_cap(), 8, "an EOF slot is not hot");
    assert_eq!(bytes(&for_slot(&reqs, 1)), HOT_ALLOWANCE_BYTES);
}

/// A prediction that moves away does not evict: the former predicted slot
/// keeps what it has read and stashed, drops to the cold class (and FIFO cap),
/// and is issued nothing new until its read-ahead falls below the cold
/// allowance.
#[test]
fn prediction_move_does_not_evict() {
    let mut f = fixture(&[64 * MIB; 3], 16, 1);
    await_file(&f.demand, 0);
    f.demand.set_predicted(Some(1));
    let _ = f.planner.plan_pass_for_test();
    let s1 = &f.slots[1];
    assert_eq!(s1.issued_bytes(), HOT_ALLOWANCE_BYTES);
    // Everything issued lands and is stashed (unclaimed).
    s1.set_read_ahead_for_test(0, HOT_ALLOWANCE_BYTES);
    for _ in 0..4 {
        f.ledger.land(SpillClass::Hot, HOT_FILL_BYTES);
    }
    f.demand.set_predicted(Some(2));
    let reqs = f.planner.plan_pass_for_test();
    assert!(for_slot(&reqs, 1).is_empty(), "16 MiB stashed > 2 MiB cold allowance");
    assert_eq!(s1.stash_bytes(), HOT_ALLOWANCE_BYTES, "nothing evicted");
    assert_eq!(s1.fifo_cap(), 8);
    // The merge claims it down to one cold fill: one more cold fill is issued.
    s1.set_read_ahead_for_test(0, MIB);
    let reqs = f.planner.plan_pass_for_test();
    let again = for_slot(&reqs, 1);
    assert_eq!(again.len(), 1, "{:?}", shape(&reqs));
    assert_eq!((class(again[0]), u64::from(again[0].len)), (SpillClass::Cold, MIB));
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
    assert_eq!(f.planner.merge_max_inflight_slices(), 4);
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

/// A source of `len` zero bytes that records the ranges advised to it.
struct AdviceRecorder {
    len: u64,
    advised: std::sync::Mutex<Vec<(u64, u64)>>,
}

impl fgumi_bam_io::pread::PositionalSource for AdviceRecorder {
    fn read_at(&self, buf: &mut [u8], offset: u64) -> std::io::Result<usize> {
        let n =
            usize::try_from(self.len.saturating_sub(offset)).unwrap_or(usize::MAX).min(buf.len());
        buf[..n].fill(0);
        Ok(n)
    }

    fn byte_len(&self) -> std::io::Result<u64> {
        Ok(self.len)
    }

    fn will_need(&self, offset: u64, len: u64) {
        self.advised.lock().expect("advice log").push((offset, len));
    }
}

/// Each fill is advised to its source as one range when it is planned: the
/// advised ranges tile the file in order, every slice read lies inside one,
/// and there is one hint per fill, not per slice.
#[test]
fn every_planned_fill_is_advised_to_its_source() {
    let len = 5 * MIB + MIB / 2 + 3;
    let source = Arc::new(AdviceRecorder { len, advised: std::sync::Mutex::default() });
    let supply = SpillSupply::new();
    let (demand, ledger) = (Arc::clone(&supply.demand), Arc::clone(&supply.ledger));
    let mut planner =
        planner(TOTAL, 16, 16, &supply).with_read_streams(ReadStreamsPolicy::fixed(4));
    let slot = Arc::new(SortMergeSlot::new(
        0,
        Arc::clone(&source) as Arc<dyn fgumi_bam_io::pread::PositionalSource>,
        0,
        len,
        SpillCodec::Zstd,
    ));
    let _ = planner.on_event_for_test(spill_ready(&slot));
    let _ = planner.on_event_for_test(all_announced(1));
    await_file(&demand, 0);
    let mut all = Vec::new();
    loop {
        let reqs = planner.plan_pass_for_test();
        if reqs.is_empty() {
            break;
        }
        slot.set_read_ahead_for_test(0, 0);
        for r in &reqs {
            ledger.land(class(r), u64::from(r.len));
        }
        all.extend(reqs);
    }
    let advised = source.advised.lock().expect("advice log").clone();
    let mut at = 0;
    for &(offset, n) in &advised {
        assert_eq!(offset, at, "advised ranges tile the file: {advised:?}");
        at += n;
    }
    assert_eq!(at, len, "{advised:?}");
    for r in &all {
        let end = r.offset + u64::from(r.len);
        assert!(
            advised.iter().any(|&(o, n)| r.offset >= o && end <= o + n),
            "slice {:?} outside every advised fill {advised:?}",
            (r.offset, r.len)
        );
    }
    assert!(advised.len() < all.len(), "one hint per fill: {advised:?} vs {}", all.len());
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
    let mut f = fixture(&[64 * MIB, 64 * MIB, 64 * MIB], 16, 1);
    await_file(&f.demand, 1);
    let _ = f.planner.plan_pass_for_test();
    let caps = |f: &Fixture| -> Vec<u32> { f.slots.iter().map(|s| s.fifo_cap()).collect() };
    assert_eq!(caps(&f), vec![8, 32, 8]);
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
    f.ledger.add_inflight(SpillClass::Hot, 0, f.planner.merge_max_inflight_slices());
    await_file(&f.demand, 2);
    assert!(f.planner.plan_pass_for_test().is_empty(), "no read: the slice cap binds");
    assert_eq!(caps(&f), vec![8, 8, 32], "the former hot slot is demoted on the next pass");
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

/// Review r2 I6: a hard error, not a `debug_assert`.
#[test]
fn all_announced_with_a_job_in_flight_is_an_error() {
    use crate::sort::protocol::ConsolidationInput;
    use crate::sort::spill_fixtures::write_run_file;
    use fgumi_pipeline_core::step::StepOutcome;
    use fgumi_pipeline_core::testing::StepProbe;
    let dir = tempfile::TempDir::new().unwrap();
    let inputs: Vec<ConsolidationInput> = (0..2u32)
        .map(|i| {
            let path = dir.path().join(format!("r{i}"));
            let len = write_run_file(&path, SpillCodec::Zstd, i, 0..5);
            ConsolidationInput { file_id: i, path, len }
        })
        .collect();
    let supply = SpillSupply::new();
    let mut pl = planner(64 << 20, 2, 2, &supply).with_consolidation(
        Arc::new(MergeDemand::new()),
        Arc::new(std::sync::atomic::AtomicBool::new(false)),
        None,
    );
    let probe = StepProbe::new(&pl);
    probe.push_input(SortPhase1Event::ConsolidationStart {
        job_id: 9,
        inputs,
        key_kind: fgumi_sort::SpillKeyKind::TemplateK40,
        index_stride: 1,
    });
    probe.push_input(SortPhase1Event::AllAnnounced {
        slot_count: 0,
        memory_chunk_count: 0,
        total_records: 0,
    });
    let err = loop {
        match probe.try_run(&mut pl) {
            Ok(StepOutcome::Finished) => panic!("finished without the error"),
            Ok(_) => {}
            Err(e) => break e,
        }
    };
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    assert_eq!(err.to_string(), "AllAnnounced arrived with consolidation job 9 in flight");
}

/// A second `AllAnnounced` is a protocol violation like the others: an
/// `InvalidData` error, not a panic, and the first announcement's merge mode
/// and budget stand.
#[test]
fn a_second_all_announced_is_an_error() {
    use fgumi_pipeline_core::step::StepOutcome;
    use fgumi_pipeline_core::testing::StepProbe;
    let mut pl = planner(TOTAL, 4, 4, &SpillSupply::new());
    let probe = StepProbe::new(&pl);
    probe.push_input(all_announced(0));
    probe.push_input(all_announced(5));
    let err = loop {
        match probe.try_run(&mut pl) {
            Ok(StepOutcome::Finished) => panic!("finished without the error"),
            Ok(_) => {}
            Err(e) => break e,
        }
    };
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    assert_eq!(err.to_string(), "AllAnnounced arrived twice");
    assert_eq!(pl.budget_for_test(), Some(ReadAheadBudget::resolve(TOTAL, 0)));
}

/// A consolidation input opens at the frame end its writer recorded
/// (`ConsolidationInput::len`): a run with an index footer opens at its
/// frames, and a recorded length that disagrees with the file (one byte
/// short of the frames, or past them into the footer) fails the job's start
/// with `InvalidData`, so a footer is never read as records. Columns: the
/// recorded length's offset from the frames' end (in bytes), a needle of
/// the error (empty: the job starts).
#[rstest::rstest]
#[case::at_the_frames(0, "")]
#[case::one_byte_short(-1, "trailer naming")]
#[case::into_the_footer(1, "trailer naming")]
fn consolidation_input_opens_at_its_written_frame_end(#[case] offset: i64, #[case] needle: &str) {
    use crate::sort::protocol::ConsolidationInput;
    use crate::sort::spill_fixtures::write_run_file;
    let dir = tempfile::TempDir::new().unwrap();
    let path = dir.path().join("r0");
    let frames_end = write_run_file(&path, SpillCodec::Zstd, 0, 0..5);
    let mut index = fgumi_sort::RunIndex::new(40);
    index.push_frame(u32::try_from(frames_end - 4).unwrap(), 1000, 5, &[0u8; 40]);
    index.rebase(4);
    let mut file = std::fs::OpenOptions::new().append(true).open(&path).unwrap();
    index.write_footer(&mut file).unwrap();
    drop(file);
    let len = frames_end.checked_add_signed(offset).unwrap();
    let supply = SpillSupply::new();
    let mut pl = planner(64 << 20, 2, 2, &supply).with_consolidation(
        Arc::new(MergeDemand::new()),
        Arc::new(std::sync::atomic::AtomicBool::new(false)),
        None,
    );
    let inputs = vec![ConsolidationInput { file_id: 0, path, len }];
    match pl.start_job(1, inputs, fgumi_sort::SpillKeyKind::TemplateK40, 1) {
        Ok(job) => {
            assert!(needle.is_empty(), "expected an error naming {needle:?}");
            assert_eq!(job.slots[0].len(), frames_end);
            assert!(job.slots[0].has_index_footer());
        }
        Err(err) => {
            assert!(!needle.is_empty(), "unexpected error: {err}");
            assert_eq!(err.kind(), std::io::ErrorKind::InvalidData, "{err}");
            assert!(err.to_string().contains(needle), "{err}");
        }
    }
}

/// Review I4 (r1): consolidation slots plan under their job budget before the
/// merge budget exists. Catches: `top_up` returning early on `budget == None`.
#[test]
fn planner_plans_consolidation_reads_before_any_merge_budget_exists() {
    use crate::sort::protocol::{ConsolidationInput, ConsolidationJob};
    use crate::sort::spill_fixtures::write_run_file;
    use fgumi_pipeline_core::step::InputHandle;
    use fgumi_pipeline_core::testing::StepProbe;
    let dir = tempfile::TempDir::new().unwrap();
    let inputs: Vec<ConsolidationInput> = (0..2u32)
        .map(|i| {
            let path = dir.path().join(format!("r{i}"));
            let len = write_run_file(&path, SpillCodec::Zstd, i, 0..500);
            ConsolidationInput { file_id: i, path, len }
        })
        .collect();
    let supply = SpillSupply::new();
    let mut pl = planner(64 << 20, 2, 2, &supply).with_consolidation(
        Arc::new(MergeDemand::new()),
        Arc::new(std::sync::atomic::AtomicBool::new(false)),
        None,
    );
    let mut probe = StepProbe::new(&pl);
    let reads = probe.take_output::<ReadRequest>(0);
    let jobs = probe.take_output::<ConsolidationJob>(2);
    probe.push_input(SortPhase1Event::ConsolidationStart {
        job_id: 1,
        inputs,
        key_kind: fgumi_sort::SpillKeyKind::TemplateK40,
        index_stride: 1,
    });
    for _ in 0..4 {
        probe.try_run(&mut pl).unwrap();
    }
    assert!(pl.merge_budget_for_test().is_none());
    let job = jobs.pop().expect("job forwarded on branch 2");
    let mut targets = Vec::new();
    while let Some(r) = reads.pop() {
        assert!(matches!(class(&r), SpillClass::ConsolidationHot | SpillClass::ConsolidationCold));
        let ReadTarget::Spill { slot, .. } = r.target else { panic!() };
        targets.push(slot.slot_id());
    }
    for s in &job.slots {
        assert!(targets.contains(&s.slot_id()), "slot {} got no read", s.slot_id());
    }
}

/// A stream run (a Bgzf-codec spill run: plain BGZF blocks plus the EOF
/// block) of 3,000 frames: its path, bytes and frames.
fn stream_run(dir: &std::path::Path) -> (PathBuf, Vec<u8>, Vec<Vec<u8>>) {
    use crate::sort::spill_fixtures::{keyed_frame, template_key, write_frames_run};
    let path = dir.join("stream");
    let frames: Vec<Vec<u8>> =
        (0..3_000u64).map(|i| keyed_frame(&template_key(i), &[3u8; 90])).collect();
    write_frames_run(&path, SpillCodec::Bgzf, &frames, 60_000);
    let bytes = std::fs::read(&path).unwrap();
    (path, bytes, frames)
}

/// Read the stream at `path` from `src` through the planner, `PreadSlices`
/// and one decompressor, pumped by hand, with the merge's side draining the
/// slot's FIFO so it never caps; assert one read in flight and nothing issued
/// after EOF. Returns the requests issued and the decoded bytes.
fn drive_stream(
    path: &std::path::Path,
    src: Arc<dyn fgumi_bam_io::pread::PositionalSource>,
) -> (u32, Vec<u8>) {
    use crate::pread::{PreadSlices, ReadSlice};
    use crate::sort::spill_decompress::SortSpillDecompress;
    use fgumi_pipeline_core::step::InputHandle;
    use fgumi_pipeline_core::testing::StepProbe;
    let slot = Arc::new(SortMergeSlot::new_sequential(0, src, SpillCodec::Bgzf));
    let supply = SpillSupply::new();
    let slices = SliceBufferPool::new(4);
    let mut pl = SpillReadPlanner::new(64 << 20, 2, 2, 1 << 20, &supply, Arc::clone(&slices));
    let mut pr = PreadSlices::spill(slices, 1 << 20);
    let mut dc = SortSpillDecompress::new(1 << 20, &supply);
    let (mut plp, mut prp, dcp) = (StepProbe::new(&pl), StepProbe::new(&pr), StepProbe::new(&dc));
    let (reads, events, out) = (
        plp.take_output::<ReadRequest>(0),
        plp.take_output::<SortPhase2Event>(1),
        prp.take_output::<ReadSlice>(0),
    );
    plp.push_input(SortPhase1Event::SpillReady {
        slot: Arc::clone(&slot),
        path: path.to_path_buf(),
        records_ingested_so_far: 0,
    });
    plp.push_input(all_announced(1));
    plp.close_input();
    let (mut issued, mut decoded) = (0u32, Vec::new());
    for _ in 0..1_000_000 {
        plp.try_run(&mut pl).unwrap();
        while let Some(r) = reads.pop() {
            assert!(r.sequential && !r.last);
            assert_eq!(
                issued,
                slot.ingested_slices(),
                "a second read issued before the first was ingested"
            );
            issued += 1;
            prp.push_input(r);
        }
        while events.pop().is_some() {}
        prp.try_run(&mut pr).unwrap();
        if let Some(s) = out.pop() {
            dcp.push_input(s);
        }
        dcp.try_run(&mut dc).unwrap();
        while let Some(block) = slot.pop_block_for_test() {
            decoded.extend(block); // the merge's side: drain the FIFO
        }
        if slot.queue_eof() && dcp.input_is_empty() {
            break;
        }
    }
    assert!(slot.eof_seen());
    for _ in 0..4 {
        plp.try_run(&mut pl).unwrap();
    }
    assert!(reads.pop().is_none(), "no request after EOF");
    assert!(!pl.planned_slot_ids_for_test().contains(&slot.slot_id()));
    (issued, decoded)
}

/// Review r2 I1 + Design issue 9: one read in flight per sequential slot, and
/// the planner stops at EOF, over a real `SequentialSource` (its pump fed by a
/// `Trickle`). Catches: two outstanding reads on a stream (a `read_at` at the
/// wrong cursor → `InvalidInput`), and requests issued after EOF.
#[test]
fn planner_keeps_one_read_in_flight_and_stops_after_a_sequential_eof() {
    use fgumi_bam_io::pread::{PositionalSource, SequentialSource, test_sources::Trickle};
    let dir = tempfile::TempDir::new().unwrap();
    let (path, bytes, frames) = stream_run(dir.path());
    let src: Arc<dyn PositionalSource> =
        Arc::new(SequentialSource::new(&path, Box::new(Trickle(std::io::Cursor::new(bytes)))));
    let (issued, decoded) = drive_stream(&path, src);
    // How many reads there are depends on what the pump holds at each; the
    // last is the zero read that is EOF.
    assert!(issued >= 2, "at least one read and the EOF read: {issued}");
    assert_eq!(decoded, frames.concat(), "every byte of the stream, in order");
}

/// A stream's gate in the planner: while its one read is in flight a pass
/// issues nothing, and once the stream's last slice is ingested `top_up`
/// marks it fully issued and never reads it again. (`plan_pass` drops a slot
/// at EOF before topping it up, so the EOF arm is driven through `top_up`.)
#[test]
fn sequential_slot_waits_for_its_read_and_stops_at_eof() {
    use fgumi_bam_io::pread::{PositionalSource, SequentialSource};
    use fgumi_pipeline_core::step::InputHandle;
    use fgumi_pipeline_core::testing::StepProbe;
    let src: Arc<dyn PositionalSource> = Arc::new(SequentialSource::new(
        std::path::Path::new("stream"),
        Box::new(std::io::Cursor::new(Vec::<u8>::new())),
    ));
    // The planner issues a stream's read only once the stream is ready, and
    // this one is ready once its pump has seen its EOF: wait for that, so the
    // passes below are not racing the pump thread.
    crate::pread::tests::wait_ready(&*src, 0);
    let slot = Arc::new(SortMergeSlot::new_sequential(0, src, SpillCodec::Bgzf));
    let supply = SpillSupply::new();
    let slices = SliceBufferPool::new(4);
    let mut pl = SpillReadPlanner::new(64 << 20, 2, 2, 1 << 20, &supply, Arc::clone(&slices));
    let mut plp = StepProbe::new(&pl);
    let reads = plp.take_output::<ReadRequest>(0);
    let _events = plp.take_output::<SortPhase2Event>(1);
    plp.push_input(SortPhase1Event::SpillReady {
        slot: Arc::clone(&slot),
        path: "stream".into(),
        records_ingested_so_far: 0,
    });
    plp.push_input(all_announced(1));
    plp.close_input();
    for _ in 0..4 {
        plp.try_run(&mut pl).unwrap();
    }
    let first = reads.pop().expect("the stream's first read");
    assert_eq!((first.seq, first.offset), (0, 0));
    assert!(reads.pop().is_none(), "one read in flight: nothing more until it is ingested");

    slot.bp_ingest_slice(0, slices.lease(Vec::new()), true).expect("the stream's last slice");
    let budget = ReadAheadBudget::resolve(64 << 20, 1);
    let mut pool_used = 0;
    assert_eq!(pl.top_up(0, SpillClass::Cold, &budget, &mut pool_used), 0, "no read past EOF");
    assert!(pl.slots[0].eof_issued, "the stream is fully issued");
}

/// The most bytes one read of a [`Capped`] stream returns: small and odd, so
/// reads end mid-block and mid-request.
const READ_CAP: usize = 61;

/// A stream whose every read returns at most [`READ_CAP`] bytes, and that
/// fails a read anywhere but its cursor.
struct Capped {
    data: Vec<u8>,
    cursor: parking_lot::Mutex<u64>,
}

impl fgumi_bam_io::pread::PositionalSource for Capped {
    fn read_at(&self, buf: &mut [u8], offset: u64) -> std::io::Result<usize> {
        let mut cursor = self.cursor.lock();
        if offset != *cursor {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("read at {offset}; the stream is at {cursor}"),
            ));
        }
        let at = usize::try_from(offset).unwrap();
        let n = buf.len().min(READ_CAP).min(self.data.len() - at);
        buf[..n].copy_from_slice(&self.data[at..at + n]);
        *cursor += n as u64;
        Ok(n)
    }

    fn byte_len(&self) -> std::io::Result<u64> {
        Err(std::io::Error::from(std::io::ErrorKind::Unsupported))
    }
}

/// Every short read resumes the stream at the bytes it delivered: with reads
/// capped at `READ_CAP` bytes, a stream of `n` bytes takes exactly
/// `ceil(n / READ_CAP)` reads,
/// each at the cursor, plus the zero read that is EOF. Catches: a resume
/// offset that skips (a read past the cursor fails `Capped`) or repeats
/// bytes, and a short read taken as EOF (a truncated stream).
#[test]
fn planner_resumes_a_stream_at_the_bytes_each_short_read_delivered() {
    let dir = tempfile::TempDir::new().unwrap();
    let (path, bytes, frames) = stream_run(dir.path());
    let len = bytes.len();
    let src = Arc::new(Capped { data: bytes, cursor: parking_lot::Mutex::new(0) });
    let (issued, decoded) = drive_stream(&path, src);
    assert_eq!(issued as usize, len.div_ceil(READ_CAP) + 1, "one request per read, plus EOF");
    assert_eq!(decoded, frames.concat(), "every byte of the stream, in order");
}

/// A stream whose `ready_at` the test controls; reads serve `data` from any
/// offset.
struct Gated {
    data: Vec<u8>,
    ready: std::sync::atomic::AtomicBool,
}

impl fgumi_bam_io::pread::PositionalSource for Gated {
    fn read_at(&self, buf: &mut [u8], offset: u64) -> std::io::Result<usize> {
        let at = usize::try_from(offset).unwrap().min(self.data.len());
        let n = buf.len().min(self.data.len() - at);
        buf[..n].copy_from_slice(&self.data[at..at + n]);
        Ok(n)
    }

    fn byte_len(&self) -> std::io::Result<u64> {
        Err(std::io::Error::from(std::io::ErrorKind::Unsupported))
    }

    fn ready_at(&self, _offset: u64) -> bool {
        self.ready.load(std::sync::atomic::Ordering::Relaxed)
    }
}

/// A stream's read is issued only once the stream is ready (has bytes, or has
/// ended), however long the merge has awaited it: the read would otherwise
/// hold a pool worker until the stream's producer writes, and that producer
/// may be waiting on another stream whose read needs the worker. Catches: a
/// stream read issued to the pool before its stream has anything to read.
#[test]
fn planner_issues_a_stream_read_only_once_the_stream_is_ready() {
    use fgumi_bam_io::pread::PositionalSource;
    let gated = Arc::new(Gated { data: vec![1u8; 4096], ready: false.into() });
    let src: Arc<dyn PositionalSource> = Arc::clone(&gated) as Arc<dyn PositionalSource>;
    let slot = Arc::new(SortMergeSlot::new_sequential(0, src, SpillCodec::Bgzf));
    let supply = SpillSupply::new();
    let slices = SliceBufferPool::new(4);
    let mut pl = SpillReadPlanner::new(64 << 20, 2, 2, 1 << 20, &supply, slices);
    let mut plp = fgumi_pipeline_core::testing::StepProbe::new(&pl);
    let reads = plp.take_output::<ReadRequest>(0);
    plp.push_input(SortPhase1Event::SpillReady {
        slot: Arc::clone(&slot),
        path: "stream".into(),
        records_ingested_so_far: 0,
    });
    plp.push_input(all_announced(1));
    for _ in 0..8 {
        plp.try_run(&mut pl).unwrap();
    }
    assert!(
        fgumi_pipeline_core::step::InputHandle::pop(&reads).is_none(),
        "a read issued for a stream with nothing to read"
    );
    gated.ready.store(true, std::sync::atomic::Ordering::Relaxed);
    plp.try_run(&mut pl).unwrap();
    let r = fgumi_pipeline_core::step::InputHandle::pop(&reads).expect("a read once ready");
    assert!(r.sequential);
    assert_eq!(r.offset, 0);
}

/// The range eligibility the planner tests use: indexed runs, phase 2 and
/// the pool `threads` wide, over `merge_memory` with its serial `R`.
fn eligibility(threads: usize, merge_memory: u64) -> crate::sort::RangeEligibility {
    crate::sort::RangeEligibility::for_test(
        threads,
        merge_memory,
        crate::sort::read_ahead_budget::serial_read_ahead(merge_memory),
    )
}

/// A planner with a range eligibility over three registered slots, every one
/// indexed unless `unindexed` names one, announced.
fn range_fixture(unindexed: Option<u32>) -> (SpillReadPlanner, SpillSupply, SortPhase2Event) {
    let supply = SpillSupply::new();
    let mut p = planner(TOTAL, 8, 16, &supply).with_range_eligibility(eligibility(8, TOTAL));
    for id in 0..3 {
        let s = slot(id, 8 * MIB);
        if unindexed != Some(id) {
            s.attach_run_index(Arc::new(fgumi_sort::RunIndex::new(8)));
        }
        let _ = p.on_event_for_test(spill_ready(&s));
    }
    let fwd = p.on_event_for_test(all_announced(3));
    (p, supply, fwd)
}

/// In range mode the planner drops the merge slots at `AllAnnounced`: nothing
/// is planned, and the event still reaches the merge.
#[test]
fn planner_range_mode_reads_nothing() {
    let (mut p, supply, fwd) = range_fixture(None);
    let Some(crate::sort::MergeMode::Range { budget, fit, .. }) = supply.merge_mode.get() else {
        panic!("range mode: {:?}", supply.merge_mode.get());
    };
    // No residual: the whole 12 GiB merge memory is the budget.
    assert_eq!(budget.bytes, TOTAL);
    // 3 runs, 8 workers: two ranges per worker at the 1.05 target.
    assert_eq!((fit.window, fit.target_range_bytes), (16, 20 * 3 * (64 << 10)));
    assert!(p.planned_slot_ids_for_test().is_empty());
    assert!(p.plan_pass_for_test().is_empty());
    assert!(matches!(fwd, SortPhase2Event::AllAnnounced { slot_count: 3, .. }));
}

/// The decision is computed from the merge slots' exact figures: their
/// count, covered frames, index raw bytes and file lengths.
#[test]
fn planner_decides_from_exact_slot_figures() {
    let supply = SpillSupply::new();
    let merge_memory = 256 * MIB;
    let elig = eligibility(8, merge_memory);
    let mut p =
        planner(TOTAL, 8, 16, &supply).with_merge_memory(merge_memory).with_range_eligibility(elig);
    let lens = [100_000, 150_000, 200_000];
    let mut index_bytes = 0u64;
    for (id, &len) in (0u32..).zip(&lens) {
        let mut index = fgumi_sort::RunIndex::new(8);
        for f in 0..=id {
            index.push_frame(1000, 60_000 + f, 10, &[u8::try_from(f).unwrap(); 8]);
        }
        index_bytes += index.heap_bytes() as u64;
        let s = slot(id, len);
        s.attach_run_index(Arc::new(index));
        let _ = p.on_event_for_test(spill_ready(&s));
    }
    let _ = p.on_event_for_test(all_announced(3));
    let inputs = crate::sort::FitInputs {
        k: 3,
        frames: 6,
        key_width: 8,
        raw_bytes: 60_000 + (60_000 + 60_001) + (60_000 + 60_001 + 60_002),
        comp_bytes: 450_000,
        max_frame_raw: 60_002,
    };
    assert_eq!(
        supply.merge_mode.get(),
        Some(&crate::sort::MergeMode::Range {
            budget: elig.budget(),
            inputs,
            fit: crate::sort::range_fit(&elig, &inputs).expect("a fit"),
            index_bytes,
        })
    );
}

/// A slot without a run index keeps the serial merge, read as before.
#[test]
fn planner_serial_mode_unchanged() {
    let (mut p, supply, fwd) = range_fixture(Some(2));
    assert_eq!(
        supply.merge_mode.get(),
        Some(&crate::sort::MergeMode::Serial(crate::sort::SerialReason::SlotNotIndexed {
            file_id: 2
        }))
    );
    assert_eq!(p.planned_slot_ids_for_test(), vec![0, 1, 2]);
    assert!(!p.plan_pass_for_test().is_empty());
    assert!(matches!(fwd, SortPhase2Event::AllAnnounced { slot_count: 3, .. }));
}

/// Without a policy (e.g. `fgumi merge`) the planner still publishes a mode.
#[test]
fn planner_without_policy_publishes_serial_no_policy() {
    let supply = SpillSupply::new();
    let mut p = planner(TOTAL, 8, 16, &supply);
    let s = slot(0, 8 * MIB);
    s.attach_run_index(Arc::new(fgumi_sort::RunIndex::new(8)));
    let _ = p.on_event_for_test(spill_ready(&s));
    let _ = p.on_event_for_test(all_announced(1));
    assert_eq!(
        supply.merge_mode.get(),
        Some(&crate::sort::MergeMode::Serial(crate::sort::SerialReason::NoPolicy))
    );
    assert_eq!(p.planned_slot_ids_for_test(), vec![0]);
}

/// A residual chunk of `records` 64 KiB records (its own arena).
fn residual(records: usize) -> Arc<crate::sort::protocol::MemoryChunkErased> {
    let chunk = fgumi_sort::InMemoryChunk::from_owned_records(
        (0..records)
            .map(|_| (fgumi_sort::RawCoordinateKey::default(), vec![0u8; 64 << 10]))
            .collect(),
    );
    Arc::new(crate::sort::protocol::MemoryChunkErased::Coordinate(chunk))
}

/// A residual chunk in memory fails closed to the serial merge: the range
/// merge needs the final run spilled. 15 runs at 32 workers fit a 4 GiB
/// merge memory without one; with one, the merge stays serial, the merge
/// slots are not dropped, and the serial merge's reads are planned.
/// Columns: residual chunks, the expected serial reason (`None`: range mode).
#[rstest::rstest]
#[case::no_residual_fits(0, None)]
#[case::residual_fails_closed(1, Some(crate::sort::SerialReason::ResidualInMemory { chunks: 1 }))]
fn planner_fails_closed_on_a_residual_chunk(
    #[case] chunks: usize,
    #[case] expected: Option<crate::sort::SerialReason>,
) {
    let supply = SpillSupply::new();
    let merge_memory = 4096 * MIB;
    let elig = eligibility(32, merge_memory);
    let mut p = planner(TOTAL, 32, 16, &supply)
        .with_merge_memory(merge_memory)
        .with_range_eligibility(elig);
    for id in 0..15 {
        let s = slot(id, 8 * MIB);
        s.attach_run_index(Arc::new(fgumi_sort::RunIndex::new(8)));
        let _ = p.on_event_for_test(spill_ready(&s));
    }
    for _ in 0..chunks {
        let _ = p.on_event_for_test(SortPhase1Event::MemoryChunk {
            chunk: residual(16),
            records_ingested_so_far: 0,
        });
    }
    let _ = p.on_event_for_test(SortPhase1Event::AllAnnounced {
        slot_count: 15,
        memory_chunk_count: u32::try_from(chunks).unwrap(),
        total_records: 7,
    });
    match (supply.merge_mode.get().expect("a decided mode"), expected) {
        (crate::sort::MergeMode::Range { budget, .. }, None) => {
            assert_eq!(*budget, elig.budget());
            assert!(p.planned_slot_ids_for_test().is_empty());
        }
        (crate::sort::MergeMode::Serial(reason), Some(expected)) => {
            assert_eq!(*reason, expected);
            assert_eq!(p.planned_slot_ids_for_test(), (0..15).collect::<Vec<u32>>());
            assert!(!p.plan_pass_for_test().is_empty(), "the serial merge is fed");
        }
        (mode, expected) => panic!("mode {mode:?}, expected {expected:?}"),
    }
}

/// An announcement whose residual chunk count differs from the chunks that
/// arrived is a protocol error, not a count to trust. Columns: chunks sent,
/// chunks announced.
#[rstest::rstest]
#[case::announced_more(0, 1)]
#[case::announced_fewer(2, 1)]
fn all_announced_with_a_chunk_count_mismatch_is_an_error(
    #[case] sent: usize,
    #[case] announced: u32,
) {
    let supply = SpillSupply::new();
    let mut p = planner(TOTAL, 32, 16, &supply);
    for _ in 0..sent {
        let _ = p.on_event_for_test(SortPhase1Event::MemoryChunk {
            chunk: residual(16),
            records_ingested_so_far: 0,
        });
    }
    let Err(err) = p.on_event(SortPhase1Event::AllAnnounced {
        slot_count: 0,
        memory_chunk_count: announced,
        total_records: 0,
    }) else {
        panic!("a mismatched announcement is an error");
    };
    assert_eq!(
        err.to_string(),
        format!("AllAnnounced names {announced} residual chunk(s) but {sent} arrived")
    );
}

/// The front spills the final run on an estimate; the planner confirms on
/// the exact figures. When they disagree, the `k + 1` runs (the last one the
/// spilled final run) merge serially with no residual chunk: 15 runs at 32
/// workers do not fit a 64 MiB budget (60 MiB of cursors alone at window 32).
#[test]
fn planner_runs_serial_over_k_plus_one_slots_when_fit_disagrees() {
    let supply = SpillSupply::new();
    let merge_memory = 64 * MIB;
    let elig = eligibility(32, merge_memory);
    let mut p = planner(TOTAL, 32, 16, &supply)
        .with_merge_memory(merge_memory)
        .with_range_eligibility(elig);
    for id in 0..15 {
        let s = slot(id, 8 * MIB);
        s.attach_run_index(Arc::new(fgumi_sort::RunIndex::new(8)));
        let _ = p.on_event_for_test(spill_ready(&s));
    }
    let fwd = p.on_event_for_test(all_announced(15));
    assert_eq!(
        supply.merge_mode.get(),
        Some(&crate::sort::MergeMode::Serial(crate::sort::SerialReason::RangeDoesNotFit {
            k: 15,
            budget: elig.budget(),
        }))
    );
    assert_eq!(p.planned_slot_ids_for_test(), (0..15).collect::<Vec<u32>>());
    assert!(!p.plan_pass_for_test().is_empty(), "the serial merge is fed");
    assert!(matches!(
        fwd,
        SortPhase2Event::AllAnnounced { slot_count: 15, memory_chunk_count: 0, .. }
    ));
}

/// When the phase-1 index policy dropped every index (after the front spilled
/// the final run), the planner keeps the serial merge over the spilled runs
/// alone and feeds it: no residual chunk, every merge slot planned.
#[test]
fn planner_runs_serial_when_the_indexes_were_dropped() {
    let supply = SpillSupply::new();
    let elig = eligibility(8, TOTAL);
    let mut p = planner(TOTAL, 8, 16, &supply).with_range_eligibility(elig);
    for id in 0..3 {
        let _ = p.on_event_for_test(spill_ready(&slot(id, 8 * MIB)));
    }
    let fwd = p.on_event_for_test(all_announced(3));
    assert_eq!(
        supply.merge_mode.get(),
        Some(&crate::sort::MergeMode::Serial(crate::sort::SerialReason::IndexesDropped {
            k: 3,
            budget: elig.budget(),
        }))
    );
    assert_eq!(p.planned_slot_ids_for_test(), vec![0, 1, 2]);
    assert!(!p.plan_pass_for_test().is_empty(), "the serial merge is fed");
    assert!(matches!(
        fwd,
        SortPhase2Event::AllAnnounced { slot_count: 3, memory_chunk_count: 0, .. }
    ));
}
