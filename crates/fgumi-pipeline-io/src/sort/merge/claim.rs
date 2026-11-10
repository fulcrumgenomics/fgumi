//! The claim table of a range-mode merge: the jobs it posted that nobody has
//! started, oldest first.
//!
//! Both the merge's own thread and the `MergeWork` workers claim from this
//! table, so a posted job the pool has not reached stays within the merge
//! thread's reach: the next result the merge needs in order is always
//! delivered, running on a worker that finishes it within one dispatch, or
//! unstarted and claimable here. The job edge carries one
//! [`MergeTicket`](crate::sort::merge_job::MergeTicket) per posted job; a
//! ticket claims whatever is oldest when a worker runs it, and is spent when
//! the table is empty. The merge thread claims the oldest unstarted assembly
//! (it produces output and releases the window) before the oldest range.
//!
//! Workers claim only while compression has a slot with nothing queued for
//! it (the *gate*): with at least `p` blocks (the phase-2 width) on the edge
//! to `BgzfCompress` and not yet popped ([`OutputDepth`]), compression is
//! the bottleneck, so the pool's threads are needed there and the merge
//! thread runs the posted jobs itself. Without a depth counter (an
//! intermediate sort, which has no compressor) the gate reads the merge's
//! output backlog instead: blocks built that the output edge refused. A
//! refused ticket is re-owed by the merge and pushed again once the gate
//! opens.

use std::collections::VecDeque;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

use std::sync::Arc;

use fgumi_pipeline_core::Job;
use parking_lot::Mutex;

use super::MergeOutput;
use crate::output_compress::OutputDepth;
use crate::sort::merge_job::{JobKind, MergeJob, Runner};

/// Who may claim a range-mode merge's posted jobs, from the hidden
/// diagnostic `FGUMI_SORT_RANGE_CLAIM` (see [`ClaimPolicy::from_env_value`]),
/// which the chain builder resolves and passes to
/// [`super::SortMerge::with_claim_policy`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ClaimPolicy {
    /// Both sides claim (the default).
    #[default]
    Auto,
    /// Only the merge's own thread claims: every ticket is spent and the
    /// merge thread runs every job, range mode's serial-equivalent floor.
    Merge,
    /// Only `MergeWork` workers claim: the merge thread posts, cuts and
    /// emits but runs no job.
    Pool,
}

impl ClaimPolicy {
    /// The policy from a `FGUMI_SORT_RANGE_CLAIM` value: `auto`, `merge` or
    /// `pool`; unset is `auto`, and any other value is `auto` with a warning
    /// that it was ignored.
    #[must_use]
    pub fn from_env_value(value: Option<&str>) -> Self {
        match value {
            None | Some("auto") => Self::Auto,
            Some("merge") => Self::Merge,
            Some("pool") => Self::Pool,
            Some(other) => {
                log::warn!(
                    "FGUMI_SORT_RANGE_CLAIM={other:?} is ignored: set it to auto, merge or pool"
                );
                Self::Auto
            }
        }
    }

    /// The value's name, as `FGUMI_SORT_RANGE_CLAIM` spells it.
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Merge => "merge",
            Self::Pool => "pool",
        }
    }

    /// Whether `runner` may claim under this policy.
    fn allows(self, runner: Runner) -> bool {
        match self {
            Self::Auto => true,
            Self::Merge => runner == Runner::Merge,
            Self::Pool => runner == Runner::Worker,
        }
    }
}

/// Jobs a range-mode merge posted and nobody has started, oldest first.
///
/// The merge posts in ordinal order. Workers claim from the front, so a
/// worker's claim always takes the lowest unstarted ordinal; the merge
/// thread takes the oldest unstarted assembly if there is one
/// ([`Self::claim_prefer_assembly`]).
pub(crate) struct ClaimTable<O: MergeOutput> {
    jobs: Mutex<VecDeque<MergeJob<O>>>,
    /// `jobs.len()`, readable without the lock.
    len: AtomicUsize,
    policy: ClaimPolicy,
    /// Tickets that found nothing to claim.
    spent: AtomicU64,
    /// The most jobs ever unstarted at once.
    high_water: AtomicUsize,
    /// Claims by `[JobKind][Runner]`.
    claims: [[AtomicU64; 2]; 2],
    /// Whether workers' claims are gated: always but under the `pool`
    /// policy, where the merge thread cannot run the jobs workers would
    /// leave.
    gate: bool,
    /// Compression's queue depth and the width `p` it is compared with:
    /// workers claim only while fewer than `p` blocks wait. `None`: the gate
    /// reads `output_backlog` instead.
    depth: Option<(Arc<OutputDepth>, u64)>,
    /// Blocks the merge built that the output edge has not accepted (the held
    /// block plus those queued behind it). Stored by the merge thread after
    /// every push; the gate's signal only without a depth counter.
    output_backlog: AtomicUsize,
    /// Worker tickets refused by the gate, not yet re-owed by the merge.
    refused: AtomicUsize,
    /// Worker tickets refused by the gate, in total (diagnostic).
    refused_total: AtomicU64,
}

impl<O: MergeOutput> ClaimTable<O> {
    /// An empty table under `policy` (gated unless the policy bars the
    /// merge thread from claiming).
    pub(crate) fn new(policy: ClaimPolicy) -> Self {
        Self {
            jobs: Mutex::new(VecDeque::new()),
            len: AtomicUsize::new(0),
            policy,
            spent: AtomicU64::new(0),
            high_water: AtomicUsize::new(0),
            claims: Default::default(),
            gate: policy != ClaimPolicy::Pool,
            depth: None,
            output_backlog: AtomicUsize::new(0),
            refused: AtomicUsize::new(0),
            refused_total: AtomicU64::new(0),
        }
    }

    /// Gate workers on compression's queue depth: they claim only while
    /// fewer than `width` blocks wait in `depth`.
    #[must_use]
    pub(crate) fn with_output_depth(mut self, depth: Arc<OutputDepth>, width: usize) -> Self {
        self.depth = Some((depth, u64::try_from(width.max(1)).unwrap_or(u64::MAX)));
        self
    }

    /// Count `n` blocks the merge pushed onto the output edge (in the depth
    /// counter, when there is one).
    pub(crate) fn note_output_pushed(&self, n: u64) {
        if let Some((depth, _)) = &self.depth
            && n > 0
        {
            depth.note_pushed(n);
        }
    }

    /// Whether compression's queue is at least the gate's width deep (`None`
    /// without a depth counter).
    pub(crate) fn output_depth_full(&self) -> Option<bool> {
        self.depth.as_ref().map(|(depth, width)| depth.depth() >= *width)
    }

    /// The width the depth gate compares with (`None` without one).
    pub(crate) fn depth_width(&self) -> Option<u64> {
        self.depth.as_ref().map(|(_, width)| *width)
    }

    /// Store the merge's output backlog: blocks built and not yet accepted by
    /// the output edge.
    pub(crate) fn set_output_backlog(&self, blocks: usize) {
        self.output_backlog.store(blocks, Ordering::Release);
    }

    /// The merge's output backlog (see [`Self::set_output_backlog`]).
    pub(crate) fn output_backlog(&self) -> usize {
        self.output_backlog.load(Ordering::Acquire)
    }

    /// Whether workers' claims are barred now: the table is gated and
    /// compression's queue holds at least `p` blocks (without a depth
    /// counter: the merge's output is backed up).
    pub(crate) fn workers_gated(&self) -> bool {
        self.gate && self.output_depth_full().unwrap_or_else(|| self.output_backlog() > 0)
    }

    /// Whether workers' claims are gated (for the `--sort-stats` lines).
    pub(crate) fn gate_enabled(&self) -> bool {
        self.gate
    }

    /// Worker tickets the gate refused so far.
    pub(crate) fn refused_total(&self) -> u64 {
        self.refused_total.load(Ordering::Relaxed)
    }

    /// Take the count of worker tickets refused by the gate (each must be
    /// owed again, so every unstarted job keeps a ticket once workers may
    /// claim again).
    pub(crate) fn take_refused(&self) -> usize {
        self.refused.swap(0, Ordering::AcqRel)
    }

    /// Add `job`, the newest posted (ordinals must be posted in order).
    pub(crate) fn post(&self, job: MergeJob<O>) {
        let mut jobs = self.jobs.lock();
        debug_assert!(
            jobs.back().is_none_or(|b| b.ordinal() < job.ordinal()),
            "jobs are posted in ordinal order"
        );
        jobs.push_back(job);
        let len = jobs.len();
        self.len.store(len, Ordering::Release);
        self.high_water.fetch_max(len, Ordering::Relaxed);
    }

    /// Take the oldest unstarted job for `runner`, marking it claimed (its
    /// queue wait ends now), or `None` when the table is empty, the policy
    /// bars `runner`, or `runner` is a worker barred by the gate
    /// (counted as refused, for the merge to re-owe the ticket).
    pub(crate) fn claim(&self, runner: Runner) -> Option<MergeJob<O>> {
        if !self.policy.allows(runner) || self.len.load(Ordering::Acquire) == 0 {
            return None;
        }
        if runner == Runner::Worker && self.workers_gated() {
            self.refused.fetch_add(1, Ordering::AcqRel);
            self.refused_total.fetch_add(1, Ordering::Relaxed);
            return None;
        }
        self.take_at(runner, |_| 0)
    }

    /// The merge thread's claim: the oldest unstarted assembly, else the
    /// oldest unstarted job. An assembly produces output and releases the
    /// window; a lookahead range does neither until much later. `None` when
    /// the table is empty or the policy bars the merge thread.
    pub(crate) fn claim_prefer_assembly(&self) -> Option<MergeJob<O>> {
        if !self.policy.allows(Runner::Merge) || self.len.load(Ordering::Acquire) == 0 {
            return None;
        }
        self.take_at(Runner::Merge, |jobs| {
            jobs.iter().position(|j| j.kind() == Some(JobKind::Assemble)).unwrap_or(0)
        })
    }

    /// Remove the job at the index `at` picks (under the lock) for `runner`,
    /// count the claim and mark it claimed.
    fn take_at(
        &self,
        runner: Runner,
        at: impl FnOnce(&VecDeque<MergeJob<O>>) -> usize,
    ) -> Option<MergeJob<O>> {
        let mut job = {
            let mut jobs = self.jobs.lock();
            let i = at(&jobs);
            let job = jobs.remove(i)?;
            self.len.store(jobs.len(), Ordering::Release);
            job
        };
        if let Some(kind) = job.kind() {
            self.claims[kind as usize][runner as usize].fetch_add(1, Ordering::Relaxed);
        }
        job.claimed(runner);
        Some(job)
    }

    /// Count a ticket that found nothing to claim.
    pub(crate) fn note_spent(&self) {
        self.spent.fetch_add(1, Ordering::Relaxed);
    }

    /// Unstarted jobs now (lock-free).
    pub(crate) fn len(&self) -> usize {
        self.len.load(Ordering::Acquire)
    }

    /// The oldest unstarted job's ordinal, if any.
    pub(crate) fn front_ordinal(&self) -> Option<u64> {
        self.jobs.lock().front().map(Job::ordinal)
    }

    /// Claims of `kind` by `runner` so far.
    pub(crate) fn claims(&self, kind: JobKind, runner: Runner) -> u64 {
        self.claims[kind as usize][runner as usize].load(Ordering::Relaxed)
    }

    /// The `--sort-stats` `Range claims:` line.
    pub(crate) fn log(&self) {
        log::info!(
            "Range claims (policy {}): merge thread {} ranges, {} assemblies; workers {} ranges, \
             {} assemblies; spent tickets {} ({} refused while compression was behind, gate \
             {}); table high-water {}",
            self.policy.name(),
            self.claims(JobKind::Range, Runner::Merge),
            self.claims(JobKind::Assemble, Runner::Merge),
            self.claims(JobKind::Range, Runner::Worker),
            self.claims(JobKind::Assemble, Runner::Worker),
            self.spent.load(Ordering::Relaxed),
            self.refused_total.load(Ordering::Relaxed),
            if self.gate { "on" } else { "off" },
            self.high_water.load(Ordering::Relaxed),
        );
    }
}

#[cfg(test)]
mod tests {
    use std::marker::PhantomData;
    use std::sync::Arc;

    use rstest::rstest;

    use super::*;
    use crate::sort::merge::BlockOutput;

    /// An assembly job with no batches at inbox ordinal `job`.
    fn job(job: u64) -> MergeJob<BlockOutput> {
        MergeJob::Assemble {
            batches: Vec::new(),
            records: 0,
            framed_bytes: 0,
            job,
            target_batch_count: 1,
            lead: None,
            timing: None,
            recycle: None,
            _output: PhantomData,
        }
    }

    /// A range job that is never run, at inbox ordinal `job`.
    fn range(job: u64) -> MergeJob<BlockOutput> {
        struct Never;
        impl fgumi_sort::RangeTask for Never {
            fn read(&mut self) -> std::io::Result<()> {
                unreachable!("a claim test never runs its jobs")
            }
            fn merge(&mut self, _: u64) -> std::io::Result<Option<fgumi_sort::RangeOut>> {
                unreachable!("a claim test never runs its jobs")
            }
        }
        MergeJob::Range {
            task: Box::new(Never),
            range: usize::try_from(job).unwrap(),
            job,
            read: false,
            timing: None,
            buffers: None,
        }
    }

    /// While the merge's output backlog is positive the gate refuses
    /// workers' claims (counting each refusal once, for the merge to re-owe)
    /// and never the merge thread's; with the backlog back at zero, or under
    /// the `pool` policy (no gate), workers claim. This is the gate's signal
    /// without a depth counter. Columns: policy, backlog, whether a worker
    /// claims, refusals.
    #[rstest]
    #[case::backlog_bars_workers(ClaimPolicy::Auto, 3, false, 1)]
    #[case::no_backlog(ClaimPolicy::Auto, 0, true, 0)]
    #[case::pool_has_no_gate(ClaimPolicy::Pool, 3, true, 0)]
    fn backlog_gate_bars_workers_only(
        #[case] policy: ClaimPolicy,
        #[case] backlog: usize,
        #[case] worker_claims: bool,
        #[case] refused: usize,
    ) {
        let table = ClaimTable::<BlockOutput>::new(policy);
        table.post(job(0));
        table.post(job(1));
        table.set_output_backlog(backlog);
        assert_eq!(table.output_backlog(), backlog);
        assert_eq!(table.claim(Runner::Worker).map(|j| j.ordinal()), worker_claims.then_some(0));
        assert_eq!(table.take_refused(), refused);
        assert_eq!(table.take_refused(), 0, "a refusal is re-owed once");
        if policy.allows(Runner::Merge) {
            let expected = u64::from(worker_claims);
            assert_eq!(table.claim_prefer_assembly().map(|j| j.ordinal()), Some(expected));
        }
        // With no backlog, a worker claims whatever is left.
        table.set_output_backlog(0);
        let left = table.len();
        assert_eq!(table.claim(Runner::Worker).is_some(), left > 0);
    }

    /// With a depth counter the gate refuses a worker iff at least `p`
    /// blocks wait for compression, whatever the backlog; the merge thread
    /// is never refused. Columns: width `p`, blocks pushed, blocks popped,
    /// whether a worker claims.
    #[rstest]
    #[case::p1_empty(1, 0, 0, true)]
    #[case::p1_one_waiting(1, 1, 0, false)]
    #[case::p1_drained(1, 3, 3, true)]
    #[case::p2_one_waiting(2, 1, 0, true)]
    #[case::p2_two_waiting(2, 2, 0, false)]
    #[case::p4_three_waiting(4, 5, 2, true)]
    #[case::p4_four_waiting(4, 6, 2, false)]
    fn depth_gate_bars_workers_at_width(
        #[case] width: usize,
        #[case] pushed: u64,
        #[case] popped: u64,
        #[case] worker_claims: bool,
    ) {
        let depth = Arc::new(OutputDepth::new());
        let table = ClaimTable::<BlockOutput>::new(ClaimPolicy::Auto)
            .with_output_depth(Arc::clone(&depth), width);
        table.post(job(0));
        table.post(job(1));
        // The backlog no longer gates: a positive one is ignored.
        table.set_output_backlog(7);
        table.note_output_pushed(pushed);
        for _ in 0..popped {
            depth.note_popped();
        }
        assert_eq!(table.claim(Runner::Worker).is_some(), worker_claims);
        assert_eq!(table.take_refused(), usize::from(!worker_claims));
        assert!(table.claim_prefer_assembly().is_some(), "the merge thread is never refused");
    }

    /// The merge thread takes the oldest unstarted assembly ahead of older
    /// ranges, else the oldest job; workers stay oldest-first.
    /// Columns: posted jobs (true: an assembly), claims (true: the merge
    /// thread), expected ordinals claimed.
    #[rstest]
    #[case::assembly_first(
        vec![false, false, true, true],
        vec![true, true, true, true],
        vec![2, 3, 0, 1]
    )]
    #[case::ranges_only(vec![false, false], vec![true, true], vec![0, 1])]
    #[case::workers_oldest_first(
        vec![false, true, false, true],
        vec![false, true, false, true],
        vec![0, 1, 2, 3]
    )]
    #[case::mixed(
        vec![false, false, true, false, true],
        vec![true, false, true, false, true],
        vec![2, 0, 4, 1, 3]
    )]
    fn merge_thread_prefers_assemblies(
        #[case] posted: Vec<bool>,
        #[case] order: Vec<bool>,
        #[case] expected: Vec<u64>,
    ) {
        let table = ClaimTable::<BlockOutput>::new(ClaimPolicy::Auto);
        for (j, &assembly) in (0u64..).zip(&posted) {
            table.post(if assembly { job(j) } else { range(j) });
        }
        let got: Vec<u64> = order
            .iter()
            .map(|&merge| {
                let job =
                    if merge { table.claim_prefer_assembly() } else { table.claim(Runner::Worker) };
                job.expect("a job remains").ordinal()
            })
            .collect();
        assert_eq!(got, expected);
        assert_eq!(table.len(), 0);
    }

    #[rstest]
    #[case::unset(None, ClaimPolicy::Auto)]
    #[case::auto(Some("auto"), ClaimPolicy::Auto)]
    #[case::merge(Some("merge"), ClaimPolicy::Merge)]
    #[case::pool(Some("pool"), ClaimPolicy::Pool)]
    #[case::other_is_auto(Some("workers"), ClaimPolicy::Auto)]
    #[case::empty_is_auto(Some(""), ClaimPolicy::Auto)]
    fn claim_policy_parses(#[case] value: Option<&str>, #[case] expected: ClaimPolicy) {
        assert_eq!(ClaimPolicy::from_env_value(value), expected);
    }

    /// Claims from either side interleave and always take the lowest
    /// unstarted ordinal; `len` tracks the table.
    /// Columns: jobs posted, claim order (true: the merge thread), expected
    /// ordinals claimed.
    #[rstest]
    #[case::merge_only(3, vec![true, true, true, true], vec![Some(0), Some(1), Some(2), None])]
    #[case::workers_only(2, vec![false, false, false], vec![Some(0), Some(1), None])]
    #[case::interleaved(
        4,
        vec![false, true, false, true],
        vec![Some(0), Some(1), Some(2), Some(3)]
    )]
    fn claims_are_oldest_first(
        #[case] posted: u64,
        #[case] order: Vec<bool>,
        #[case] expected: Vec<Option<u64>>,
    ) {
        let table = ClaimTable::<BlockOutput>::new(ClaimPolicy::Auto);
        for j in 0..posted {
            table.post(job(j));
        }
        assert_eq!(table.len(), usize::try_from(posted).unwrap());
        let got: Vec<Option<u64>> = order
            .iter()
            .map(|&merge| {
                let runner = if merge { Runner::Merge } else { Runner::Worker };
                table.claim(runner).map(|j| j.ordinal())
            })
            .collect();
        assert_eq!(got, expected);
        assert_eq!(table.len(), 0);
    }

    /// The policy bars one side from claiming. With one job posted per side,
    /// each side claims its own job or is barred, and only the claims made
    /// are counted, against their runner.
    /// Columns: policy, whether the merge thread / a worker gets a job.
    #[rstest]
    #[case::auto(ClaimPolicy::Auto, true, true)]
    #[case::merge(ClaimPolicy::Merge, true, false)]
    #[case::pool(ClaimPolicy::Pool, false, true)]
    fn policy_bars_a_side(
        #[case] policy: ClaimPolicy,
        #[case] merge_claims: bool,
        #[case] worker_claims: bool,
    ) {
        let table = ClaimTable::<BlockOutput>::new(policy);
        table.post(job(0));
        table.post(job(1));
        assert_eq!(table.claim(Runner::Worker).map(|j| j.ordinal()), worker_claims.then_some(0));
        let merge_job = u64::from(worker_claims);
        assert_eq!(
            table.claim(Runner::Merge).map(|j| j.ordinal()),
            merge_claims.then_some(merge_job)
        );
        assert_eq!(table.claims(JobKind::Assemble, Runner::Worker), u64::from(worker_claims));
        assert_eq!(table.claims(JobKind::Assemble, Runner::Merge), u64::from(merge_claims));
        assert_eq!(table.len(), 2 - usize::from(worker_claims) - usize::from(merge_claims));
    }

    /// Concurrent claims never hand out a job twice and never lose one.
    #[test]
    fn concurrent_claims_never_duplicate() {
        let table = Arc::new(ClaimTable::<BlockOutput>::new(ClaimPolicy::Auto));
        let n = 4000u64;
        let poster = {
            let table = Arc::clone(&table);
            std::thread::spawn(move || {
                for j in 0..n {
                    table.post(job(j));
                }
            })
        };
        let claimers: Vec<_> = (0..4)
            .map(|i| {
                let table = Arc::clone(&table);
                std::thread::spawn(move || {
                    let runner = if i == 0 { Runner::Merge } else { Runner::Worker };
                    let mut got = Vec::new();
                    let mut idle = 0;
                    while idle < 100_000 {
                        match table.claim(runner) {
                            Some(j) => {
                                got.push(j.ordinal());
                                idle = 0;
                            }
                            None => idle += 1,
                        }
                    }
                    assert!(got.windows(2).all(|w| w[0] < w[1]), "oldest first: {got:?}");
                    got
                })
            })
            .collect();
        poster.join().unwrap();
        let mut all: Vec<u64> = claimers.into_iter().flat_map(|c| c.join().unwrap()).collect();
        while let Some(j) = table.claim(Runner::Worker) {
            all.push(j.ordinal());
        }
        all.sort_unstable();
        assert_eq!(all, (0..n).collect::<Vec<_>>());
        assert_eq!(table.len(), 0);
    }
}
