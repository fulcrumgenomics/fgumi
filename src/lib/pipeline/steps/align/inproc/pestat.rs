//! `CohortPeStatStep` — the `Serial` cohort barrier of the in-process bwa-mem3
//! backend.
//!
//! Accumulates every sub-batch ([`ExtendedWork`]) of a cohort — which arrives
//! in any order off the `Parallel`
//! [`AlignSeedExtendStep`](super::seed_extend::AlignSeedExtendStep) — until the
//! cohort's [`CohortCloser`] (riding on the sub-batch with the highest
//! `index_in_cohort`, only known at cut time) has arrived AND every sub-batch it
//! declares has been received. On completion it runs
//! [`AlignEngine::infer_cohort`] **exactly once**, over the cohort's resident
//! state (every sub-batch has seeded into it by then), then emits one
//! [`PairWork`] per sub-batch
//! (also in `index_in_cohort` order) carrying the shared `Arc<E::PeStat>` (or
//! `None` for a sub-batch with no PE reads) and the sub-batch's computed
//! [`IdBases`](super::engine::IdBases) ([`id_bases`]).
//!
//! **Why an SE-only sub-batch still waits for the closer.** Its id bases need
//! `cohort_n_se` (`id_bases`'s second argument), which — like `cohort_n_pe` —
//! is only known once the closer arrives. So every
//! sub-batch, PE or SE-only, is held in the same per-cohort accumulator; there
//! is no early-emission path for pure-SE sub-batches.
//!
//! **No gather.** The sub-batches' reads and alignment regions live in the
//! cohort's resident state, which [`AlignEngine::infer_cohort`] reads directly,
//! so the barrier only holds each sub-batch's small work item until the cohort
//! completes ([`CohortPeStatStep::finish_cohort`]).
//!
//! This `Serial` cohort barrier is wired between
//! [`AlignSeedExtendStep`](super::seed_extend::AlignSeedExtendStep) and
//! [`AlignPairEmitStep`](super::pair_emit::AlignPairEmitStep) by
//! [`InProcessBwaMem3Backend`](super::InProcessBwaMem3Backend).

use std::collections::{BTreeMap, VecDeque};
use std::io;
use std::sync::Arc;

use super::FlushOutcome;
use super::cohort::{CohortCloser, ExtendedWork, PairWork, id_bases};
use super::engine::{AlignEngine, cohort_of};
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::Single;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};

/// One cohort's in-flight accumulation state.
struct CohortAccum<E: AlignEngine> {
    /// Received sub-batches, keyed by `index_in_cohort`. A `BTreeMap` so that
    /// once the cohort completes, iterating it yields sub-batches in
    /// `index_in_cohort` order for free — independent of arrival order, which
    /// is unordered off the `Parallel` seed/extend step.
    subs: BTreeMap<u32, ExtendedWork<E>>,
    /// Set when the sub-batch carrying the [`CohortCloser`] arrives. May
    /// arrive before or after the cohort's other sub-batches.
    expected: Option<CohortCloser>,
}

impl<E: AlignEngine> CohortAccum<E> {
    fn new() -> Self {
        Self { subs: BTreeMap::new(), expected: None }
    }

    /// A cohort is complete once its closer has arrived and every sub-batch it
    /// declares has been received — independent of arrival order. `subs` holds
    /// only distinct, in-range indices (`accumulate` rejects the rest), so its
    /// length is the received count.
    fn is_complete(&self) -> bool {
        self.expected
            .is_some_and(|closer| usize::try_from(closer.n_sub_batches) == Ok(self.subs.len()))
    }
}

/// The `Serial` cohort barrier.
pub(crate) struct CohortPeStatStep<E: AlignEngine> {
    /// Shared, read-only alignment engine.
    engine: Arc<E>,
    /// Byte limit for the `PairWork` output queue.
    output_byte_limit: u64,
    /// Per-cohort accumulation state, keyed by cohort index.
    cohorts: BTreeMap<u32, CohortAccum<E>>,
    /// `PairWork` of completed cohorts that cannot be released yet because a
    /// lower-indexed cohort is still incomplete, keyed by cohort index. Released
    /// strictly in cohort order (see [`Self::release_in_order`]). Bounded by the
    /// prepare step's two-cohort in-flight gate: at most one completed cohort
    /// can wait here behind the one still accumulating.
    parked: BTreeMap<u32, Vec<PairWork<E>>>,
    /// The next cohort index to release to `ready`. Cohorts are numbered densely
    /// from 0 by the prepare step, so this walks 0, 1, 2, ….
    next_cohort_to_release: u32,
    /// `PairWork` of released cohorts, awaiting push to the output.
    ready: VecDeque<PairWork<E>>,
    /// A single output item bounced by backpressure, retried before new pushes.
    held: HeldSlot<Unpushed<PairWork<E>>>,
}

impl<E: AlignEngine> CohortPeStatStep<E> {
    /// Build the step over a shared `engine`, sizing its output queue at
    /// `output_byte_limit` bytes.
    pub(crate) fn new(engine: Arc<E>, output_byte_limit: u64) -> Self {
        Self {
            engine,
            output_byte_limit,
            cohorts: BTreeMap::new(),
            parked: BTreeMap::new(),
            next_cohort_to_release: 0,
            ready: VecDeque::new(),
            held: HeldSlot::new(),
        }
    }

    /// File one arrived sub-batch into its cohort's accumulator.
    ///
    /// # Errors
    /// A duplicate arrival — a second `CohortCloser` for a cohort, or the same
    /// `index_in_cohort` received twice — is a hard error in every build, not
    /// just a `debug_assert!`: this barrier's whole job is exact-once cohort
    /// completion, and silently overwriting an entry
    /// would corrupt `is_complete()`'s count without any signal,
    /// symmetric with the "lost sub-batch" (too-few) error at drain.
    ///
    /// So is an `index_in_cohort` outside the `0..n_sub_batches` range the
    /// closer declares, whether it arrives before or after the closer: counting
    /// alone would let a cohort declaring 3 sub-batches complete on indices
    /// `{0, 1, 5}` with sub-batch 2 missing. Distinct in-range indices also
    /// bound the received count by `n_sub_batches`, so over-delivery cannot
    /// occur.
    fn accumulate(&mut self, work: ExtendedWork<E>) -> io::Result<()> {
        let cohort = work.id.cohort;
        let index = work.id.index_in_cohort;
        let closer = work.closer;
        let accum = self.cohorts.entry(cohort).or_insert_with(CohortAccum::new);
        if let Some(closer) = closer {
            if accum.expected.is_some() {
                return Err(io::Error::other(format!(
                    "align-and-merge (in-process) pestat: cohort {cohort} received a second \
                     CohortCloser (on sub-batch {index}) — a sub-batch was double-emitted upstream"
                )));
            }
            if let Some((&max_index, _)) = accum.subs.last_key_value()
                && max_index >= closer.n_sub_batches
            {
                return Err(io::Error::other(format!(
                    "align-and-merge (in-process) pestat: cohort {cohort}'s CohortCloser declares \
                     {n} sub-batches but sub-batch {max_index} already arrived — a sub-batch \
                     index was corrupted upstream",
                    n = closer.n_sub_batches,
                )));
            }
            accum.expected = Some(closer);
        }
        if let Some(closer) = accum.expected
            && index >= closer.n_sub_batches
        {
            return Err(io::Error::other(format!(
                "align-and-merge (in-process) pestat: cohort {cohort} received sub-batch {index}, \
                 outside the {n} sub-batches its CohortCloser declares — a sub-batch index was \
                 corrupted upstream",
                n = closer.n_sub_batches,
            )));
        }
        if accum.subs.contains_key(&index) {
            return Err(io::Error::other(format!(
                "align-and-merge (in-process) pestat: cohort {cohort} received sub-batch {index} \
                 twice — a sub-batch was double-emitted upstream"
            )));
        }
        accum.subs.insert(index, work);
        Ok(())
    }

    /// Run `infer_cohort` once (only if the cohort has any PE reads) and
    /// return one `PairWork` per sub-batch, in `index_in_cohort` order (the
    /// order `subs` arrives in).
    ///
    /// The model is computed over the cohort's resident state directly: every
    /// sub-batch has seeded its reads into it by now, so there is nothing to
    /// gather.
    ///
    /// # Errors
    /// Propagates the engine's `infer_cohort` failure.
    fn finish_cohort(
        &self,
        cohort: u32,
        closer: CohortCloser,
        subs: Vec<ExtendedWork<E>>,
    ) -> io::Result<Vec<PairWork<E>>> {
        let pestat: Option<Arc<E::PeStat>> = if closer.cohort_n_pe > 0 {
            // Any sub-batch's lease reaches the cohort's resident state; a
            // cohort with PE reads has at least one sub-batch.
            let lease = &subs.first().expect("a cohort with PE reads has sub-batches").lease;
            let resident = cohort_of(&*self.engine, lease)?;
            let stat = self.engine.infer_cohort(resident).map_err(|e| {
                io::Error::other(format!(
                    "align-and-merge (in-process) pestat: infer_cohort failed for cohort \
                     {cohort}: {e:#}"
                ))
            })?;
            Some(Arc::new(stat))
        } else {
            None
        };

        // Mirrors upstream's `[M::mem_pestat]` stderr line: the FR insert-size
        // summary is only meaningful once `infer_cohort` actually ran, so it's
        // appended only in that branch.
        // Gated on the log level so the summary String is not built when info
        // logging is off.
        if log::log_enabled!(log::Level::Info) {
            let pestat_state = match &pestat {
                Some(stat) => format!("computed ({})", E::pestat_summary(stat)),
                None => "none (no PE reads)".to_string(),
            };
            log::info!(
                "align-and-merge (in-process) pestat: cohort {cohort} complete — \
                 n_sub_batches={n_sub_batches} n_pe={n_pe} n_se={n_se} pestat={pestat_state}",
                n_sub_batches = closer.n_sub_batches,
                n_pe = closer.cohort_n_pe,
                n_se = closer.cohort_n_se,
            );
        }

        // Only sub-batches that hold pairs get the model; an SE-only sub-batch
        // pairs nothing and takes `None`.
        Ok(subs
            .into_iter()
            .map(|work| {
                let ExtendedWork { id, pos, closer: _, layout, unmapped, lease, ranges } = work;
                let pestat = if E::n_pairs(&ranges) > 0 { pestat.clone() } else { None };
                PairWork {
                    id,
                    layout,
                    unmapped,
                    lease,
                    ranges,
                    pestat,
                    ids: id_bases(&pos, closer.cohort_n_se),
                }
            })
            .collect())
    }

    /// Park a completed cohort's `PairWork`, then move every cohort that is now
    /// next in line from `parked` to `ready`, in cohort order.
    ///
    /// Completion order is not cohort order: seed/extend is `Parallel` and
    /// unordered, so a straggling sub-batch of cohort N can outlast all of
    /// N+1's. Releasing N+1 first would let the downstream `ByItemOrdinal`
    /// reorder fill with N+1 while N's `PairWork` waits behind a full output —
    /// a cyclic wait. Releasing in cohort order keeps the reorder's next
    /// expected item always reachable.
    fn release_in_order(&mut self, cohort: u32, pair_works: Vec<PairWork<E>>) {
        self.parked.insert(cohort, pair_works);
        while let Some(works) = self.parked.remove(&self.next_cohort_to_release) {
            self.ready.extend(works);
            self.next_cohort_to_release += 1;
        }
    }

    /// Push staged `ready` items (and any held-back one) to the output.
    fn flush(&mut self, ctx: &mut StepCtx<'_, Self>) -> FlushOutcome {
        if let Some(unpushed) = self.held.take()
            && let Err(again) = ctx.outputs.retry(unpushed)
        {
            self.held.put(again);
            return FlushOutcome::StillHeld;
        }
        while let Some(work) = self.ready.pop_front() {
            if let Err(unpushed) = ctx.outputs.push(work) {
                self.held.put(unpushed);
                return FlushOutcome::NewlyHeld;
            }
        }
        FlushOutcome::Clear
    }
}

impl<E: AlignEngine> Step for CohortPeStatStep<E> {
    type Input = ExtendedWork<E>;
    type Outputs = Single<PairWork<E>>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "CohortPeStat",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        loop {
            // 1. Flush pending output first. A newly held item is forward
            //    progress this call; re-failing an already-held one is not.
            if let Some(outcome) = self.flush(ctx).blocked_outcome() {
                return Ok(outcome);
            }

            // 2. Pull the next sub-batch and file it into its cohort.
            let Some(work) = ctx.input.pop() else {
                if ctx.input.is_drained() {
                    // Every cohort still tracked here is, by construction,
                    // incomplete (a complete one is removed the moment its
                    // last sub-batch arrives) — a lost sub-batch or closer.
                    if let Some((&cohort, accum)) = self.cohorts.iter().next() {
                        return Err(io::Error::other(format!(
                            "align-and-merge (in-process) pestat: input drained with cohort \
                             {cohort} incomplete ({received} of {expected} sub-batches received) \
                             — a sub-batch or its CohortCloser was lost",
                            received = accum.subs.len(),
                            expected = accum
                                .expected
                                .map_or_else(|| "?".to_string(), |c| c.n_sub_batches.to_string()),
                        )));
                    }
                    // A cohort parked here waits on a lower cohort that never
                    // completed and is no longer tracked: its sub-batches were
                    // all lost (no accumulator was ever created).
                    if let Some(&cohort) = self.parked.keys().next() {
                        return Err(io::Error::other(format!(
                            "align-and-merge (in-process) pestat: input drained with cohort \
                             {cohort} complete but cohort {missing} never seen — every \
                             sub-batch of cohort {missing} was lost",
                            missing = self.next_cohort_to_release,
                        )));
                    }
                    return Ok(StepOutcome::Finished);
                }
                return Ok(StepOutcome::NoProgress);
            };

            let cohort = work.id.cohort;
            self.accumulate(work)?;

            // 3. If that completed the cohort, run `infer_cohort` once and
            //    stage its `PairWork`s — in cohort order — for the next flush.
            let complete = self.cohorts.get(&cohort).is_some_and(CohortAccum::is_complete);
            if complete {
                let accum = self.cohorts.remove(&cohort).expect("just checked complete above");
                let closer = accum.expected.expect("a complete cohort carries a closer");
                let subs: Vec<ExtendedWork<E>> = accum.subs.into_values().collect();
                let pair_works = self.finish_cohort(cohort, closer, subs)?;
                self.release_in_order(cohort, pair_works);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::VecDeque;
    use std::sync::{Arc, Mutex};

    use rstest::rstest;

    use super::*;
    use crate::pipeline::core::builder::{Pipeline, PipelineConfig};
    use crate::pipeline::core::item::HeapSize;
    use crate::pipeline::steps::align::inproc::cohort::{CohortPos, Layout, SubBatchId};
    use crate::pipeline::steps::align::inproc::engine::fake::{FakeCohort, FakeEngine, TaggedRead};
    use crate::pipeline::steps::align::inproc::engine::{
        EngineBatch, IdBases, RecordOrigin, RecordSink,
    };
    use crate::pipeline::steps::align::inproc::gate::CohortLease;

    // ----- RecordingEngine: wraps FakeEngine, recording infer_cohort calls -

    /// Wraps [`FakeEngine`], recording each `infer_cohort` call as the number
    /// of pairs its cohort holds, so tests can assert one call per cohort over
    /// that whole cohort without any FFI.
    #[derive(Clone, Default)]
    struct RecordingEngine {
        inner: FakeEngine,
        infer_calls: Arc<Mutex<Vec<usize>>>,
    }

    impl RecordingEngine {
        fn new() -> Self {
            Self::default()
        }

        /// Every recorded `infer_cohort` call (its cohort's pair count), in
        /// call order.
        fn calls(&self) -> Vec<usize> {
            self.infer_calls.lock().expect("mutex not poisoned").clone()
        }
    }

    impl AlignEngine for RecordingEngine {
        type Scratch = <FakeEngine as AlignEngine>::Scratch;
        type Cohort = <FakeEngine as AlignEngine>::Cohort;
        type Ranges = <FakeEngine as AlignEngine>::Ranges;
        type PeStat = <FakeEngine as AlignEngine>::PeStat;

        fn new_scratch(&self) -> anyhow::Result<Self::Scratch> {
            self.inner.new_scratch()
        }

        fn new_cohort(&self) -> anyhow::Result<Self::Cohort> {
            self.inner.new_cohort()
        }

        fn seed_extend(
            &self,
            scratch: &mut Self::Scratch,
            cohort: &Self::Cohort,
            batch: EngineBatch<'_>,
        ) -> anyhow::Result<Self::Ranges> {
            self.inner.seed_extend(scratch, cohort, batch)
        }

        fn infer_cohort(&self, cohort: &Self::Cohort) -> anyhow::Result<Self::PeStat> {
            let pairs = cohort.seeded_pairs.load(std::sync::atomic::Ordering::Relaxed);
            self.infer_calls.lock().expect("mutex not poisoned").push(pairs);
            self.inner.infer_cohort(cohort)
        }

        fn pair_emit(
            &self,
            scratch: &mut Self::Scratch,
            cohort: &Self::Cohort,
            regs: Self::Ranges,
            pestat: Option<&Self::PeStat>,
            ids: IdBases,
            sink: &mut dyn RecordSink,
        ) -> anyhow::Result<()> {
            self.inner.pair_emit(scratch, cohort, regs, pestat, ids, sink)
        }

        fn n_pairs(regs: &Self::Ranges) -> usize {
            FakeEngine::n_pairs(regs)
        }

        fn pestat_summary(pestat: &Self::PeStat) -> String {
            FakeEngine::pestat_summary(pestat)
        }
    }

    // ----- fixtures ------------------------------------------------------

    /// One synthetic sub-batch's regs: `n_pe` pairs (R1, R2 tagged
    /// `Pair(i)`), then `n_se` singles (tagged `Single(j)`) — the same shape
    /// `FakeEngine::seed_extend` produces (engine.rs), built directly since
    /// these tests exercise the barrier, not seed/extend.
    // `n_pe`/`n_se` are the domain terms (also `CohortPos`'s field names); the
    // pair reads more clearly together than renamed apart.
    #[allow(clippy::similar_names)]
    fn regs_for(n_pe: u32, n_se: u32) -> Vec<TaggedRead> {
        let mut regs = Vec::new();
        for i in 0..n_pe as usize {
            regs.push(TaggedRead { origin: RecordOrigin::Pair(i), mate: 0 });
            regs.push(TaggedRead { origin: RecordOrigin::Pair(i), mate: 1 });
        }
        for j in 0..n_se as usize {
            regs.push(TaggedRead { origin: RecordOrigin::Single(j), mate: 0 });
        }
        regs
    }

    /// Build one `ExtendedWork<RecordingEngine>` from its cohort-position
    /// fields, keeping the `n_se`/`n_pe` domain-term pairing (matches
    /// `CohortPos`'s own field names). `lease` is the cohort's shared lease.
    #[allow(clippy::similar_names, clippy::too_many_arguments)]
    fn extended_work(
        lease: &CohortLease,
        serial: u64,
        cohort: u32,
        index_in_cohort: u32,
        se_offset: u64,
        pe_offset: u64,
        n_se: u32,
        n_pe: u32,
        closer: Option<CohortCloser>,
    ) -> ExtendedWork<RecordingEngine> {
        ExtendedWork {
            id: SubBatchId { serial, cohort, index_in_cohort },
            pos: CohortPos { cohort_read_base: 0, se_offset, pe_offset, n_se, n_pe },
            closer,
            layout: Layout::Mixed(Vec::new()),
            unmapped: Vec::new(),
            lease: lease.clone(),
            ranges: regs_for(n_pe, n_se),
        }
    }

    /// Build one cohort's sub-batches from a `(n_se, n_pe)` spec per
    /// sub-batch, in `index_in_cohort` order, stamping the `CohortCloser` on
    /// the last. Returns the sub-batches plus a `(SubBatchId, CohortPos)`
    /// snapshot per sub-batch for computing expected `id_bases` after the
    /// originals are consumed. The sub-batches share one lease, whose resident
    /// cohort is pre-seeded with the cohort's pairs (as the seed/extend step
    /// would leave it), so `infer_cohort` observes the whole cohort.
    #[allow(clippy::similar_names)]
    fn build_cohort(
        cohort: u32,
        serial_start: u64,
        counts: &[(u32, u32)],
    ) -> (Vec<ExtendedWork<RecordingEngine>>, Vec<(SubBatchId, CohortPos)>, CohortCloser) {
        let n_sub_batches = u32::try_from(counts.len()).expect("fits u32");
        let cohort_n_se: u64 = counts.iter().map(|&(se, _)| u64::from(se)).sum();
        let cohort_n_pe: u64 = counts.iter().map(|&(_, pe)| u64::from(pe)).sum();
        let closer = CohortCloser { n_sub_batches, cohort_n_se, cohort_n_pe };
        let lease = CohortLease::for_test();
        let resident: &FakeCohort =
            cohort_of(&RecordingEngine::new(), &lease).expect("fake cohort is infallible");
        resident.seeded_pairs.store(
            usize::try_from(cohort_n_pe).expect("fits usize"),
            std::sync::atomic::Ordering::Relaxed,
        );

        let mut se_offset = 0u64;
        let mut pe_offset = 0u64;
        let mut items = Vec::with_capacity(counts.len());
        let mut snapshot = Vec::with_capacity(counts.len());
        for (i, &(n_se, n_pe)) in counts.iter().enumerate() {
            let idx = u32::try_from(i).expect("fits u32");
            let is_last = idx + 1 == n_sub_batches;
            let work = extended_work(
                &lease,
                serial_start + u64::from(idx),
                cohort,
                idx,
                se_offset,
                pe_offset,
                n_se,
                n_pe,
                is_last.then_some(closer),
            );
            snapshot.push((work.id, work.pos));
            items.push(work);
            se_offset += u64::from(n_se);
            pe_offset += u64::from(n_pe);
        }
        (items, snapshot, closer)
    }

    // ----- test source / sink --------------------------------------------

    /// Replays pre-built items into the chain, in `items`' order. `Exclusive`
    /// (owns cursor) — the queue is a plain FIFO (no reordering declared), so
    /// the consumer sees items in exactly this order.
    struct WorkSource<T: Send + HeapSize + 'static> {
        items: VecDeque<T>,
        held: HeldSlot<Unpushed<T>>,
    }

    impl<T: Send + HeapSize + 'static> Step for WorkSource<T> {
        type Input = ();
        type Outputs = Single<T>;

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "WorkSource",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 1 << 20 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(unpushed) = self.held.take()
                && let Err(again) = ctx.outputs.retry(unpushed)
            {
                self.held.put(again);
                return Ok(StepOutcome::Contention);
            }
            let Some(item) = self.items.pop_front() else {
                return Ok(StepOutcome::Finished);
            };
            if let Err(unpushed) = ctx.outputs.push(item) {
                self.held.put(unpushed);
            }
            Ok(StepOutcome::Progress)
        }
    }

    /// Terminal sink accumulating every input item in arrival order.
    struct CollectSink<T: Send + HeapSize + 'static> {
        collected: Arc<Mutex<Vec<T>>>,
    }

    impl<T: Send + HeapSize + 'static> Step for CollectSink<T> {
        type Input = T;
        type Outputs = ();

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "CollectSink",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(item) => {
                    self.collected.lock().expect("mutex not poisoned").push(item);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    /// Run `items` through a fresh `CohortPeStatStep<RecordingEngine>`,
    /// returning the emitted `PairWork` (arrival order on the unordered
    /// output) or the pipeline error's message.
    fn drive(
        engine: RecordingEngine,
        items: Vec<ExtendedWork<RecordingEngine>>,
    ) -> Result<Vec<PairWork<RecordingEngine>>, String> {
        let collected: Arc<Mutex<Vec<PairWork<RecordingEngine>>>> =
            Arc::new(Mutex::new(Vec::new()));
        let builder = Pipeline::builder();
        builder
            .chain(WorkSource { items: items.into(), held: HeldSlot::new() })
            .chain(CohortPeStatStep::new(Arc::new(engine), 4 * 1024 * 1024))
            .chain(CollectSink { collected: Arc::clone(&collected) })
            .into_sink_marker();
        let pipeline = builder.build().expect("chain builds");
        pipeline
            .run(PipelineConfig { threads: 1, ..Default::default() })
            .map_err(|e| e.to_string())?;
        let out = Arc::try_unwrap(collected)
            .unwrap_or_else(|_| panic!("sole owner after run completes"))
            .into_inner()
            .expect("mutex not poisoned");
        Ok(out)
    }

    /// Assert every `PairWork` in `out` carries the `id_bases` the pure
    /// function computes from its snapshot `CohortPos` and `cohort_n_se`.
    fn assert_ids_match<'a>(
        snapshot: &[(SubBatchId, CohortPos)],
        cohort_n_se: u64,
        out: impl IntoIterator<Item = &'a PairWork<RecordingEngine>>,
    ) {
        for pw in out {
            let (_, pos) = snapshot
                .iter()
                .find(|(id, _)| id.index_in_cohort == pw.id.index_in_cohort)
                .expect("every emitted PairWork's index_in_cohort matches a snapshot sub-batch");
            assert_eq!(
                pw.ids,
                id_bases(pos, cohort_n_se),
                "sub-batch {idx}: ids must match id_bases(pos, cohort_n_se)",
                idx = pw.id.index_in_cohort,
            );
        }
    }

    // ----- tests: single cohort, all PE -----------------------------------

    #[test]
    fn single_cohort_all_pe_infers_once_over_the_whole_cohort() {
        let engine = RecordingEngine::new();
        let (items, snapshot, closer) = build_cohort(0, 0, &[(0, 2), (0, 3), (0, 1)]);

        let out = drive(engine.clone(), items).expect("pipeline runs");
        assert_eq!(out.len(), 3, "one PairWork per sub-batch");
        let in_index_order: Vec<u32> = out.iter().map(|pw| pw.id.index_in_cohort).collect();
        assert_eq!(in_index_order, vec![0, 1, 2], "PairWork emitted in index_in_cohort order");

        assert_eq!(engine.calls(), vec![6], "infer_cohort called once, over all 6 pairs");

        let pestats: Vec<Arc<()>> =
            out.iter().map(|pw| pw.pestat.clone().expect("PE sub-batch has a pestat")).collect();
        assert!(Arc::ptr_eq(&pestats[0], &pestats[1]), "pestat is the same Arc across sub-batches");
        assert!(Arc::ptr_eq(&pestats[1], &pestats[2]), "pestat is the same Arc across sub-batches");

        assert_ids_match(&snapshot, closer.cohort_n_se, &out);
    }

    /// Out-of-order arrival: the closer rides on the LAST sub-batch by
    /// `index_in_cohort`, but sub-batches can arrive in any order off the
    /// unordered `Parallel` seed/extend step. Whatever order they arrive in,
    /// the cohort completes only once all have arrived, `infer_cohort` runs
    /// exactly once over the whole cohort, and the `PairWork`s come out in
    /// `index_in_cohort` order, not arrival order.
    #[rstest]
    #[case::in_order(vec![0, 1, 2])]
    #[case::closer_first(vec![2, 0, 1])]
    #[case::reverse(vec![2, 1, 0])]
    #[case::middle_last(vec![0, 2, 1])]
    fn out_of_order_arrival_completes_with_one_infer_cohort_call(
        #[case] arrival_order: Vec<usize>,
    ) {
        let engine = RecordingEngine::new();
        let (mut items, snapshot, closer) = build_cohort(0, 0, &[(0, 1), (0, 2), (0, 1)]);

        // Reorder `items` to the arrival order under test; `items[i]` still
        // carries its original `index_in_cohort`, so completion detection and
        // the emission order do not depend on this ordering.
        let mut by_index: Vec<Option<ExtendedWork<RecordingEngine>>> =
            items.drain(..).map(Some).collect();
        let ordered: Vec<ExtendedWork<RecordingEngine>> = arrival_order
            .iter()
            .map(|&i| by_index[i].take().expect("each index used once"))
            .collect();

        let out = drive(engine.clone(), ordered).expect("pipeline runs");
        assert_eq!(out.len(), 3, "one PairWork per sub-batch despite out-of-order arrival");
        let in_index_order: Vec<u32> = out.iter().map(|pw| pw.id.index_in_cohort).collect();
        assert_eq!(in_index_order, vec![0, 1, 2], "emitted in index_in_cohort order");
        assert_eq!(engine.calls(), vec![4], "exactly one call, over the cohort's 4 pairs");
        assert_ids_match(&snapshot, closer.cohort_n_se, &out);
    }

    // ----- tests: SE-only sub-batches ---------------------------------------

    /// A pure-SE sub-batch is held until the cohort's closer like a PE one,
    /// and is emitted with `pestat: None`; the PE sub-batches share the one
    /// cohort pestat.
    #[rstest]
    #[case::se_only_first(vec![0, 1, 2])]
    #[case::se_only_before_closer(vec![1, 0, 2])]
    #[case::se_only_after_closer_arrives_last(vec![2, 1, 0])]
    fn se_only_sub_batch_is_held_until_closer_and_emitted_with_no_pestat(
        #[case] arrival_order: Vec<usize>,
    ) {
        let engine = RecordingEngine::new();
        // index 0: SE-only (2 singles, 0 pairs); indices 1, 2: PE. Closer on
        // index 2 (the last by `index_in_cohort`).
        let (mut items, snapshot, closer) = build_cohort(0, 0, &[(2, 0), (0, 2), (0, 1)]);

        let mut by_index: Vec<Option<ExtendedWork<RecordingEngine>>> =
            items.drain(..).map(Some).collect();
        let ordered: Vec<ExtendedWork<RecordingEngine>> = arrival_order
            .iter()
            .map(|&i| by_index[i].take().expect("each index used once"))
            .collect();

        let out = drive(engine.clone(), ordered).expect("pipeline runs");
        assert_eq!(out.len(), 3);
        assert_eq!(engine.calls(), vec![3], "one infer_cohort call over the cohort's 3 pairs");

        let se_only = out.iter().find(|pw| pw.id.index_in_cohort == 0).expect("SE-only sub-batch");
        assert!(se_only.pestat.is_none(), "SE-only sub-batch gets pestat: None");

        let pe_a = out.iter().find(|pw| pw.id.index_in_cohort == 1).expect("PE sub-batch 1");
        let pe_b = out.iter().find(|pw| pw.id.index_in_cohort == 2).expect("PE sub-batch 2");
        let pestat_a = pe_a.pestat.clone().expect("PE sub-batch has a pestat");
        let pestat_b = pe_b.pestat.clone().expect("PE sub-batch has a pestat");
        assert!(Arc::ptr_eq(&pestat_a, &pestat_b), "the two PE sub-batches share one pestat Arc");

        assert_ids_match(&snapshot, closer.cohort_n_se, &out);
    }

    /// An all-SE cohort never calls `infer_cohort` and every sub-batch gets
    /// `pestat: None`.
    #[test]
    fn all_se_cohort_never_calls_infer_cohort() {
        let engine = RecordingEngine::new();
        let (items, snapshot, closer) = build_cohort(0, 0, &[(2, 0), (3, 0)]);

        let out = drive(engine.clone(), items).expect("pipeline runs");
        assert_eq!(out.len(), 2);
        assert!(engine.calls().is_empty(), "no PE reads in the cohort ⇒ infer_cohort never runs");
        assert!(out.iter().all(|pw| pw.pestat.is_none()));
        assert_ids_match(&snapshot, closer.cohort_n_se, &out);
    }

    // ----- tests: multiple independent cohorts ------------------------------

    /// Two cohorts, interleaved on arrival, each get their own `infer_cohort`
    /// call over exactly their own resident cohort.
    #[test]
    fn multiple_cohorts_each_get_their_own_infer_cohort_call() {
        let engine = RecordingEngine::new();
        let (cohort0, snap0, closer0) = build_cohort(0, 0, &[(0, 1), (0, 1)]);
        let (cohort1, snap1, closer1) = build_cohort(1, 100, &[(0, 2), (0, 1)]);

        // Interleave: c0[0], c1[0], c0[1] (closes cohort 0), c1[1] (closes
        // cohort 1) — cohort 0 completes strictly before cohort 1.
        let mut c0 = cohort0.into_iter();
        let mut c1 = cohort1.into_iter();
        let ordered = vec![
            c0.next().expect("c0[0]"),
            c1.next().expect("c1[0]"),
            c0.next().expect("c0[1]"),
            c1.next().expect("c1[1]"),
        ];

        let out = drive(engine.clone(), ordered).expect("pipeline runs");
        assert_eq!(out.len(), 4, "one PairWork per sub-batch across both cohorts");

        assert_eq!(
            engine.calls(),
            vec![2, 3],
            "one call per cohort, over its own pairs: cohort 0 (2) completes first, then 1 (3)"
        );

        let out0: Vec<&PairWork<RecordingEngine>> =
            out.iter().filter(|pw| pw.id.cohort == 0).collect();
        let out1: Vec<&PairWork<RecordingEngine>> =
            out.iter().filter(|pw| pw.id.cohort == 1).collect();
        assert_eq!(out0.len(), 2);
        assert_eq!(out1.len(), 2);
        assert_ids_match(&snap0, closer0.cohort_n_se, out0);
        assert_ids_match(&snap1, closer1.cohort_n_se, out1);
    }

    /// Completion order is not release order: cohort 1 completes first (its
    /// `infer_cohort` runs first), but nothing is emitted until cohort 0
    /// completes, and then cohort 0's `PairWork`s precede cohort 1's.
    #[test]
    fn cohorts_are_released_in_cohort_order_not_completion_order() {
        let engine = RecordingEngine::new();
        let (cohort0, snap0, closer0) = build_cohort(0, 0, &[(0, 1), (0, 1)]);
        let (cohort1, snap1, closer1) = build_cohort(1, 100, &[(0, 2), (0, 1)]);

        // c0[0], then all of cohort 1 (completing it), then c0[1] (completing 0).
        let mut c0 = cohort0.into_iter();
        let mut ordered = vec![c0.next().expect("c0[0]")];
        ordered.extend(cohort1);
        ordered.push(c0.next().expect("c0[1]"));

        let out = drive(engine.clone(), ordered).expect("pipeline runs");
        assert_eq!(engine.calls(), vec![3, 2], "cohort 1 (3 pairs) completes before cohort 0 (2)");
        let released: Vec<(u32, u32)> =
            out.iter().map(|pw| (pw.id.cohort, pw.id.index_in_cohort)).collect();
        assert_eq!(released, vec![(0, 0), (0, 1), (1, 0), (1, 1)], "released in cohort order");

        let out0: Vec<&PairWork<RecordingEngine>> =
            out.iter().filter(|pw| pw.id.cohort == 0).collect();
        let out1: Vec<&PairWork<RecordingEngine>> =
            out.iter().filter(|pw| pw.id.cohort == 1).collect();
        assert_ids_match(&snap0, closer0.cohort_n_se, out0);
        assert_ids_match(&snap1, closer1.cohort_n_se, out1);
    }

    /// A complete cohort parked behind a lower cohort none of whose sub-batches
    /// ever arrived is a hard error at drain, not a silent `Finished` that drops
    /// the parked cohort.
    #[test]
    fn cohort_parked_behind_a_never_seen_cohort_is_an_error_at_drain() {
        let (cohort1, _snapshot, _closer) = build_cohort(1, 100, &[(0, 1)]);
        let Err(err) = drive(RecordingEngine::new(), cohort1) else {
            panic!("a cohort parked behind a lost cohort must error");
        };
        assert!(
            err.contains("cohort 1 complete but cohort 0 never seen"),
            "error names both cohorts: {err}"
        );
    }

    // ----- tests: Finished / error contract ---------------------------------

    /// A cohort missing its closer (and thus incomplete forever) is a hard
    /// error once input drains — not a silent `Finished`.
    #[test]
    fn incomplete_cohort_at_drain_is_an_error() {
        let engine = RecordingEngine::new();
        let (mut items, _snapshot, _closer) = build_cohort(0, 0, &[(0, 1), (0, 1), (0, 1)]);
        items.pop(); // drop the last sub-batch — its CohortCloser never arrives.

        // `PairWork` is not `Debug` (it carries the engine's ranges), so
        // `Result::expect_err` (which requires `T: Debug`) is unavailable —
        // destructure explicitly instead.
        let Err(err) = drive(engine, items) else {
            panic!("an incomplete cohort at drain must error");
        };
        assert!(err.contains("cohort 0 incomplete"), "error names cohort 0: {err}");
        assert!(err.contains("2 of ? sub-batches received"), "error names the counts: {err}");
    }

    /// The upstream mis-deliveries the barrier rejects on arrival, rather than
    /// completing a cohort wrongly or reporting it "lost" at drain.
    #[derive(Clone, Copy, Debug)]
    enum Misdelivery {
        /// Sub-batch 5 arrives after a closer declaring 2 sub-batches.
        IndexOutOfRangeAfterCloser,
        /// Sub-batch 5 arrives before a closer declaring 2 sub-batches.
        IndexOutOfRangeBeforeCloser,
        /// Two sub-batches of one cohort both carry a closer.
        DuplicateCloser,
        /// Sub-batch 0 arrives twice.
        DuplicateIndex,
    }

    /// The arrival sequence for a `Misdelivery`: one cohort (0) of one-pair
    /// sub-batches sharing a lease.
    fn misdelivered_items(case: Misdelivery) -> Vec<ExtendedWork<RecordingEngine>> {
        let lease = CohortLease::for_test();
        let closer = |n_sub_batches: u32| CohortCloser {
            n_sub_batches,
            cohort_n_se: 0,
            cohort_n_pe: u64::from(n_sub_batches),
        };
        let sub = |serial: u64, index: u32, closes: Option<CohortCloser>| {
            extended_work(&lease, serial, 0, index, 0, u64::from(index), 0, 1, closes)
        };
        match case {
            Misdelivery::IndexOutOfRangeAfterCloser => {
                vec![sub(1, 1, Some(closer(2))), sub(5, 5, None), sub(0, 0, None)]
            }
            Misdelivery::IndexOutOfRangeBeforeCloser => {
                vec![sub(5, 5, None), sub(1, 1, Some(closer(2)))]
            }
            Misdelivery::DuplicateCloser => {
                vec![sub(0, 0, None), sub(2, 2, Some(closer(3))), sub(1, 1, Some(closer(3)))]
            }
            Misdelivery::DuplicateIndex => vec![sub(0, 0, None), sub(1, 0, None)],
        }
    }

    #[rstest]
    #[case::index_out_of_range_after_closer(
        Misdelivery::IndexOutOfRangeAfterCloser,
        "received sub-batch 5, outside the 2 sub-batches its CohortCloser declares"
    )]
    #[case::index_out_of_range_before_closer(
        Misdelivery::IndexOutOfRangeBeforeCloser,
        "declares 2 sub-batches but sub-batch 5 already arrived"
    )]
    #[case::duplicate_closer(Misdelivery::DuplicateCloser, "received a second CohortCloser")]
    #[case::duplicate_index(Misdelivery::DuplicateIndex, "received sub-batch 0 twice")]
    fn misdelivered_sub_batch_is_an_error(#[case] case: Misdelivery, #[case] expected: &str) {
        let Err(err) = drive(RecordingEngine::new(), misdelivered_items(case)) else {
            panic!("{case:?} must error");
        };
        assert!(err.contains(expected), "{case:?}: expected {expected:?} in: {err}");
    }

    /// A fully-processed run (every cohort complete) reaches `Finished`
    /// cleanly — the happy-path tests above all rely on this implicitly
    /// (`drive` would otherwise hang or error), but this test pins the
    /// zero-cohort (empty input) degenerate case explicitly.
    #[test]
    fn empty_input_finishes_cleanly_with_no_output() {
        let engine = RecordingEngine::new();
        let out = drive(engine.clone(), Vec::new()).expect("pipeline runs");
        assert!(out.is_empty());
        assert!(engine.calls().is_empty());
    }

    // ----- profile -----------------------------------------------------------

    #[test]
    fn profile_is_serial_bytebounded_unordered() {
        let step = CohortPeStatStep::new(Arc::new(RecordingEngine::new()), 4 * 1024 * 1024);
        let p = step.profile();
        assert_eq!(p.name, "CohortPeStat");
        assert_eq!(p.kind, StepKind::Serial);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::None]);
        assert!(matches!(p.output_queues[0], QueueSpec::ByteBounded { .. }));
    }
}
