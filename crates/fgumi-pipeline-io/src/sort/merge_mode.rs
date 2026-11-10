//! How the sort's final merge runs: serially (one k-way merge fed by the
//! [`super::SpillReadPlanner`]) or as key ranges merged in parallel from the
//! spill runs' in-memory indexes.
//!
//! The planner decides once, at `AllAnnounced`, before any merge read is
//! planned, and publishes the decision through [`super::SpillSupply`].
//!
//! Range mode spills the final run, so the merge reads only spill runs and
//! holds no sort buffer. `range_fit` sizes it in one place: it chooses the
//! index stride `N`, the range size `T` and the window of live ranges per
//! worker. The target is decode amplification `A = 1 + kNF / T` of 1.05
//! (`A_TARGET_RECIP`), bounded at 1.25 (`A_MAX_RECIP`); a live range holds
//! `rho = 2 + S_c / S` bytes per raw byte; and the memory used (run indexes,
//! per-run cursors, live ranges) must fit `max(R, merge_memory)`, where `R` is
//! the serial read-ahead floor, with the run indexes at most a quarter of the
//! merge memory (`INDEX_BUDGET_DIVISOR`). The phase-1 front owns the run-index
//! policy ([`RunIndexPolicy`]): at each run boundary it picks the stride from
//! a conservative estimate of this fit. The indexes are merge-phase state:
//! each run's index is written to its spill file's footer and, when the final
//! run was spilled, read back at `AllAnnounced`, after the phase-1 arena is
//! retired, and the merge budget pays for them. Phase 1 holds only the open
//! run's index and an in-flight consolidation's, each at its stride; a
//! resident final run leaves the footers unread (the merge is serial). When
//! no sizing fits, the merge is serial ([`SerialReason`]).

use std::fmt;
use std::sync::Arc;

use fgumi_bgzf::BGZF_MAX_BLOCK_SIZE;
use fgumi_sort::{RunIndex, SortMergeSlot, SortOrder, SpillCodec, SpillKeyKind, spill_key_width};

use super::read_ahead_budget::serial_read_ahead;

/// The fewest phase-2 threads a range merge runs with; below this (a
/// [`RangeEligibility`] built with none) the serial merge is used.
///
/// Range mode runs on the merge's own thread plus the `p` phase-2 slots, as
/// the serial merge does; the slots compress first and the merge thread runs
/// every job they have not reached, and compresses itself once compression
/// is behind. It is measured at least as fast as the serial merge on both
/// orders from one slot up, so one slot is enough.
pub const RANGE_MIN_PHASE2_THREADS: usize = 1;

/// Key ranges live per worker in range mode at most: posted, merged but not
/// yet cut into batches, or still referenced by a batch not yet emitted.
/// `range_fit` tries this many per worker (at least [`RANGE_MIN_WINDOW`])
/// and one per worker.
pub(crate) const RANGE_WINDOW_PER_WORKER: usize = 2;

/// The smallest window `range_fit` prefers, whatever the phase-2 width: at
/// one or two phase-2 threads, `2p` live ranges leave the merge thread one
/// or three ranges of lookahead, and this floor gives it slack against a
/// slow worker (the memory model charges every live range, so the fit pays
/// for it like any other window). The one-per-worker fallback stays below
/// it when memory is too tight for the floor.
pub(crate) const RANGE_MIN_WINDOW: usize = 4;

/// The raw bytes of one spill frame the fit assumes (`F`): the spill path cuts
/// frames at [`SPILL_FRAME_BYTES`] (one BGZF block's payload, 65,280 bytes),
/// rounded up to 64 KiB. A record longer than that is a frame of its own;
/// [`FitInputs::max_frame_raw`] carries the largest such frame, and the fit
/// charges `max(F, max_frame_raw)` per frame ([`FitInputs::frame_bytes`]).
pub(crate) const FRAME_RAW_BYTES: u64 = 64 << 10;

/// The framed bytes at which the spill path cuts a frame (`B`): `SpillGather`,
/// `SpillBlockCompress` and `ConsolidateRuns` cut a frame before the record
/// that would take it past one BGZF block's payload (65,280 bytes). A record
/// longer than that is a frame of its own. So for any frame but a run's last,
/// it plus the next frame's first record exceed `B`: two consecutive frames
/// hold more than `B` bytes, and a frame holds more than `B - m` when no
/// record exceeds `m` bytes framed.
pub(crate) const SPILL_FRAME_BYTES: u64 = BGZF_MAX_BLOCK_SIZE as u64;

const _: () = assert!(FRAME_RAW_BYTES >= SPILL_FRAME_BYTES);

/// Frames a running range holds per run whatever its size: the run's
/// compressed boundary frame (at most one frame) and its cursor's decode
/// buffer (one frame raw), each charged [`FitInputs::frame_bytes`]. Only
/// started ranges hold these: a posted range holds nothing until its job
/// starts, and a merged range holds only its output.
pub(crate) const RANGE_CURSOR_FRAMES: u64 = 2;

/// Smallest target range size (`T_min`): a range job costs about one pool
/// round trip and `k` preads against about 5 ms of merge per MiB per worker,
/// so 1 MiB keeps per-range overhead at about 2 %.
pub(crate) const MIN_RANGE_BYTES: u64 = 1 << 20;

/// `1 / (A_target - 1)`: the solver aims for decode amplification
/// `A_target = 1.05` when memory allows, and spends no memory below it.
pub(crate) const A_TARGET_RECIP: u64 = 20;

/// `1 / (A_max - 1)`: decode amplification is bounded at `A_max = 1.25`;
/// a range merge that cannot reach it within the budget runs serially.
pub(crate) const A_MAX_RECIP: u64 = 4;

/// The sparsest run index `range_fit` considers: one entry per this many
/// frames. A sparsified entry sums its frames' raw lengths in a `u32`; at
/// `MAX_INDEX_STRIDE * FRAME_RAW_BYTES` = 1024 x 64 KiB = 64 MiB that fits
/// with room to spare for frames of at most `F`. Frames of very long records
/// (up to 64 MiB each) can overflow it; [`RunIndex::sparsify`] then fails,
/// leaving the index unchanged, and `SpillWrite` drops that index, so the
/// merge falls back to serial.
pub(crate) const MAX_INDEX_STRIDE: u32 = 1 << 10;

/// The run indexes a fit may hold, as a share of the merge memory: a feasible
/// fit needs `I(N) <= merge_memory / INDEX_BUDGET_DIVISOR`. A merge-phase
/// bound only (the indexes are read back from the runs' footers once phase 1
/// is over); the phase-1 front's estimate uses the same cap, so a front fit
/// implies a planner fit. For `-m` of 1 GiB and more it never binds (the
/// cursor term fails the fit first). It binds only where the `R` floor lifts
/// the budget far above `-m`, and keeps a tiny sort from holding an index up
/// to `R`.
pub(crate) const INDEX_BUDGET_DIVISOR: u64 = 4;

/// What the range merge needs to know to decide whether it may run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RangeEligibility {
    /// The sort writes indexed spill runs (coordinate or template-coordinate
    /// order with the zstd spill codec).
    pub indexable: bool,
    /// The resolved phase-2 thread count (`p`).
    pub phase2_threads: usize,
    /// `MergeWork` worker copies (one per pool worker): each runs at most
    /// one claimed range at a time, so this (plus the merge's own thread)
    /// bounds how many ranges have read their spans at once. Taken as at
    /// least `phase2_threads`.
    pub pool_threads: usize,
    /// Force the serial merge (`FGUMI_SORT_SERIAL_MERGE=1`).
    pub force_serial: bool,
    /// The memory the final merge may use (bytes).
    pub merge_memory: u64,
    /// The serial merge's read-ahead `R` (bytes) at `merge_memory`, the range
    /// budget's floor.
    pub serial_r: u64,
}

impl RangeEligibility {
    /// The eligibility of a sort in `sort_order` that spills with `codec`:
    /// its runs are indexable when the order is coordinate or
    /// template-coordinate and the codec zstd (the one policy every indexing
    /// step follows, through [`Self::can_range`]); the budget floor is the
    /// serial merge's read-ahead `R` at `merge_memory`.
    #[must_use]
    pub fn for_sort(
        sort_order: SortOrder,
        codec: SpillCodec,
        phase2_threads: usize,
        pool_threads: usize,
        force_serial: bool,
        merge_memory: u64,
    ) -> Self {
        Self {
            indexable: matches!(sort_order, SortOrder::Coordinate | SortOrder::TemplateCoordinate)
                && codec == SpillCodec::Zstd,
            phase2_threads,
            pool_threads,
            force_serial,
            merge_memory,
            serial_r: serial_read_ahead(merge_memory),
        }
    }

    /// The range budget once nothing but the merge holds memory:
    /// `max(R, merge_memory)`.
    #[must_use]
    pub fn budget(&self) -> RangeBudget {
        RangeBudget::new(self.merge_memory, self.serial_r)
    }

    /// The run indexes a fit may hold: `merge_memory / INDEX_BUDGET_DIVISOR`,
    /// one cap for the planner and the phase-1 front alike, so a front fit
    /// implies a planner fit.
    #[must_use]
    pub fn index_cap(&self) -> u64 {
        self.merge_memory / INDEX_BUDGET_DIVISOR
    }

    /// An indexable, not forced-serial eligibility with `threads` phase-2 and
    /// pool threads (tests).
    #[cfg(test)]
    pub(crate) fn for_test(threads: usize, merge_memory: u64, serial_r: u64) -> Self {
        Self {
            indexable: true,
            phase2_threads: threads,
            pool_threads: threads,
            force_serial: false,
            merge_memory,
            serial_r,
        }
    }

    /// Whether the merge could ever run by key range under this eligibility:
    /// indexable runs, not forced serial, and enough phase-2 threads. Only
    /// then are spill runs indexed at all (a sort that cannot range builds no
    /// frames and writes no index footers); `indexable` itself is
    /// kept as configured so the merge decision still reports
    /// [`SerialReason::ForcedByEnv`] / [`SerialReason::TooFewThreads`].
    #[must_use]
    pub fn can_range(&self) -> bool {
        self.indexable && !self.force_serial && self.phase2_threads >= RANGE_MIN_PHASE2_THREADS
    }

    /// Whether a merge over `inputs` would run as ranges under this
    /// eligibility: [`Self::can_range`], at least two runs, and a fit within
    /// [`Self::budget`] (tests).
    #[cfg(test)]
    #[must_use]
    pub fn spills_final_run(&self, inputs: &FitInputs) -> bool {
        self.can_range() && inputs.k >= 2 && range_fit(self, inputs).is_some()
    }
}

/// The phase-1 front's run-index policy, shared by both fronts (`ReadBlocks`
/// and `SortBuffer`): the single owner of each run's index stride. The spill
/// steps only follow the stride the front stamps on the run's events (and on
/// `AllAnnounced`). It charges nothing to the sort buffer: the indexes live in
/// the runs' footers until the merge, so the front's run boundaries are the
/// serial sort's.
///
/// At every run boundary the front calls [`Self::at_run`] with the runs
/// sealed so far plus the run starting (estimated at the last run's size):
/// a fit over [`FitInputs::estimated`] at the current stride or sparser, with
/// the index at most the planner's cap ([`INDEX_BUDGET_DIVISOR`]), sets the
/// stride to the fit's; no fit turns indexing off for good (stride `0` from
/// then on). Every input is a function of the input and the settings, so the
/// strides are too.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct RunIndexPolicy {
    /// The sort's eligibility; `None` indexes nothing.
    elig: Option<RangeEligibility>,
    /// The current run's index stride (`0` = no index).
    stride: u32,
}

impl RunIndexPolicy {
    /// The policy for a sort under `elig`: stride 1 (dense) from the first
    /// run when the sort can merge by range, else 0 (no index, exactly the
    /// unindexed configuration).
    pub(crate) fn new(elig: Option<RangeEligibility>) -> Self {
        let stride = u32::from(elig.is_some_and(|e| e.can_range()));
        Self { elig, stride }
    }

    /// The current run's index stride (`0` = no index).
    pub(crate) fn stride(&self) -> u32 {
        self.stride
    }

    /// Decide the stride of a run that makes `k` runs with the runs sealed so
    /// far, of `kind`, over `spill` (the sealed runs' framed bytes plus the
    /// starting run's estimate). The first run (`k < 2`) keeps the stride:
    /// there is nothing to merge it with yet. An unknown kind (no record
    /// provisioned the key context yet) keeps it for this run only; a kind
    /// with no index turns indexing off.
    pub(crate) fn at_run(&mut self, k: usize, kind: Option<SpillKeyKind>, spill: SpillEstimate) {
        if self.stride == 0 || k < 2 {
            return;
        }
        let Some(kind) = kind else { return };
        self.stride = self.fit(k, kind, spill).map_or(0, |fit| fit.stride);
    }

    /// The final run's stride when it spills, `None` when it stays the
    /// residual: the sort still indexes, the final run makes `k >= 2` runs
    /// over `spill`, and a range merge over them fits at the current stride
    /// or sparser. The stride is at least the current one.
    pub(crate) fn final_stride(
        &self,
        k: usize,
        kind: Option<SpillKeyKind>,
        spill: SpillEstimate,
    ) -> Option<u32> {
        if self.stride == 0 || k < 2 {
            return None;
        }
        self.fit(k, kind?, spill).map(|fit| fit.stride)
    }

    /// The fit over the estimate at the current stride or sparser, under the
    /// planner's index cap; `None` for a sort that cannot range or a kind
    /// with no index.
    fn fit(&self, k: usize, kind: SpillKeyKind, spill: SpillEstimate) -> Option<RangeFit> {
        let elig = self.elig.as_ref().filter(|e| e.can_range())?;
        spill_key_width(kind)?;
        range_fit_from(elig, &FitInputs::estimated(k, kind, spill), self.stride)
    }
}

/// What a phase-1 front knows of the spill bytes a merge would read, in the
/// spill framing (`[key prefix][u32 len][body]` per record).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SpillEstimate {
    /// Framed bytes of runs whose largest framed record is known.
    pub framed: u64,
    /// The largest framed record of those runs (`None`: unknown, and they
    /// are bounded as [`Self::unbounded`] is).
    pub max_record: Option<u64>,
    /// Framed bytes of runs whose largest record is not known (a run still
    /// being read, estimated at the previous run's size).
    pub unbounded: u64,
}

impl SpillEstimate {
    /// Framed bytes of `records` records of `kind` over `len_prefixed` bytes
    /// (each record's `u32` length prefix included, as the arena holds them):
    /// the template key prefix is added per record, the length prefix is not
    /// counted twice.
    #[must_use]
    pub fn framed_bytes(kind: Option<SpillKeyKind>, len_prefixed: u64, records: u64) -> u64 {
        let prefix = kind.map_or(0, |k| k.serialized_key_bytes() as u64);
        len_prefixed.saturating_add(records.saturating_mul(prefix))
    }

    /// An upper bound on the frames `k` runs over these bytes are cut into.
    /// Within one run every frame but the last, plus the next frame, holds
    /// more than `B` = [`SPILL_FRAME_BYTES`], so a run of `S` framed bytes has
    /// fewer than `2S / B + 1` frames; when no record exceeds `m < B` bytes,
    /// every frame but the last holds more than `B - m`, so fewer than
    /// `S / (B - m) + 1`. Summed over runs: `frames <= bound(known) +
    /// ceil(2 x unbounded / B) + k`.
    #[must_use]
    pub fn frame_bound(&self, k: usize) -> u64 {
        let b = SPILL_FRAME_BYTES;
        let loose = |framed: u64| framed.saturating_mul(2).div_ceil(b);
        let known = match self.max_record {
            Some(m) if m < b => loose(self.framed).min(self.framed.div_ceil(b - m)),
            _ => loose(self.framed),
        };
        let k = u64::try_from(k).unwrap_or(u64::MAX);
        known.saturating_add(loose(self.unbounded)).saturating_add(k)
    }

    /// Every framed byte, known and unbounded.
    #[must_use]
    pub fn total(&self) -> u64 {
        self.framed.saturating_add(self.unbounded)
    }
}

/// The spill runs' figures the range fit is computed from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FitInputs {
    /// Spill runs (`k`).
    pub k: usize,
    /// Frames across every run (the dense, stride-1 entry count).
    pub frames: u64,
    /// Width in bytes of each index entry's key.
    pub key_width: u16,
    /// Raw (decompressed) bytes across every run (`S`).
    pub raw_bytes: u64,
    /// Compressed (file) bytes across every run (`S_c`).
    pub comp_bytes: u64,
    /// The largest raw bytes of any one frame, when a record longer than one
    /// frame made a frame above `F`; `0` (or anything at most `F`) when every
    /// frame is at most `F`.
    pub max_frame_raw: u64,
}

impl FitInputs {
    /// The fit inputs a phase-1 front estimates at a run boundary, from
    /// deterministic figures only: `k` runs of `key_kind` over `spill`. Every
    /// figure is an upper bound on what the planner later reads from the
    /// runs' indexes:
    ///
    /// - frames: `SpillEstimate::frame_bound`, which holds for any record
    ///   length (the tight form only where the largest record is known);
    /// - raw bytes: the framed bytes, which is what the indexes' raw lengths
    ///   sum to;
    /// - compressed bytes: the worst case `S_c = S` (`rho = 3`), as
    ///   compression only shrinks what a live range holds.
    ///
    /// - the largest frame: `max_record` of the runs whose largest record is
    ///   known; a run still being read may hold a longer record, which only
    ///   the planner's exact figure ([`RunIndex::max_frame_raw`]) sees.
    ///
    /// So when the estimate fits, the planner's exact fit at the same stride
    /// fits too, except where a record above one frame in a run still being
    /// read is not yet known (the BAM front bounds the final run's record
    /// count at EOF by the smallest BAM record, so its key bytes are an upper
    /// bound too). Then the planner, which decides from the runs' indexes, may
    /// still fall back to serial; the output and the `-m` bound are the same.
    ///
    /// # Panics
    ///
    /// Panics if `key_kind` is not an indexed kind ([`spill_key_width`]):
    /// only coordinate and template-coordinate runs carry an index.
    #[must_use]
    pub fn estimated(k: usize, key_kind: SpillKeyKind, spill: SpillEstimate) -> Self {
        let raw_bytes = spill.total();
        let key_width = spill_key_width(key_kind).expect("an indexed spill key kind");
        Self {
            k,
            frames: spill.frame_bound(k),
            key_width,
            raw_bytes,
            comp_bytes: raw_bytes,
            max_frame_raw: spill.max_record.unwrap_or(0),
        }
    }

    /// The raw bytes the fit charges for one frame: `F`, or the largest frame
    /// when a long record made one larger.
    #[must_use]
    pub fn frame_bytes(&self) -> u64 {
        FRAME_RAW_BYTES.max(self.max_frame_raw)
    }

    /// Heap bytes the runs' indexes hold at one entry per `stride` frames
    /// (`I(N)`): `(k + (frames - k) / N) x e`, rounded up. Every run keeps
    /// its first entry, so sparsifying cannot shrink a one-frame run; fewer
    /// frames than runs count as one frame per run.
    #[must_use]
    pub fn index_bytes(&self, stride: u32) -> u64 {
        let k = u64::try_from(self.k).unwrap_or(u64::MAX);
        let rest = self.frames.saturating_sub(k).div_ceil(u64::from(stride.max(1)));
        k.saturating_add(rest).saturating_mul(RunIndex::entry_bytes(self.key_width))
    }

    /// Bytes a live range holds per raw byte (`rho = S_c / S + 2`): its
    /// compressed span, its merged output and its assembled batches. 3 (the
    /// `S_c = S` worst case) when no raw byte is known.
    #[must_use]
    #[allow(clippy::cast_precision_loss, reason = "a ratio for display")]
    pub fn rho(&self) -> f64 {
        if self.raw_bytes == 0 { 3.0 } else { self.comp_bytes as f64 / self.raw_bytes as f64 + 2.0 }
    }
}

/// The range merge's sizing, as `range_fit` chose it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RangeFit {
    /// Frames per run-index entry (`N`, a power of two).
    pub stride: u32,
    /// Live ranges the merge keeps.
    pub window: usize,
    /// Raw bytes per range the plan aims for (`T`).
    pub target_range_bytes: u64,
    /// Bytes decoded per range beyond its own: `k x stride x F`, one boundary
    /// entry per run (`F` = [`FitInputs::frame_bytes`]).
    pub boundary_bytes: u64,
    /// Bytes range mode holds at once at this sizing.
    pub mem: u64,
}

impl RangeFit {
    /// Planned decode amplification: `1 + boundary_bytes / target_range_bytes`.
    #[must_use]
    #[allow(clippy::cast_precision_loss, reason = "a ratio for display and tests")]
    pub fn amplification(&self) -> f64 {
        1.0 + self.boundary_bytes as f64 / self.target_range_bytes as f64
    }
}

/// Bytes range data may hold in range mode, and where that figure came from.
///
/// `bytes = max(R, merge_memory)`: in range mode the sort buffer is empty, so
/// the merge draws on the buffer's own allowance; the `R` floor keeps tiny
/// `-m` sorts in range mode on the grounds that the serial merge spends `R`
/// above the buffer there.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RangeBudget {
    /// The budget range data is charged to (bytes).
    pub bytes: u64,
    /// The merge memory (bytes).
    pub merge_memory: u64,
    /// The serial merge's read-ahead `R` (bytes), the budget's floor.
    pub serial_r: u64,
}

impl RangeBudget {
    /// The range budget for a merge with `merge_memory` and serial read-ahead
    /// `serial_r`.
    #[must_use]
    pub fn new(merge_memory: u64, serial_r: u64) -> Self {
        Self { bytes: merge_memory.max(serial_r), merge_memory, serial_r }
    }
}

impl fmt::Display for RangeBudget {
    /// `range budget X MiB (merge memory Y MiB, floor R Z MiB)`, every figure
    /// rounded down.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "range budget {} MiB (merge memory {} MiB, floor R {} MiB)",
            self.bytes >> 20,
            self.merge_memory >> 20,
            self.serial_r >> 20,
        )
    }
}

/// Size the range merge over `inputs` within `elig`'s budget
/// ([`RangeEligibility::budget`]) at any stride; `None` when no sizing keeps
/// decode amplification within `A_max` (1.25) inside it (tests: production
/// fits through [`range_fit_at`] and [`RunIndexPolicy`]).
///
/// See `range_fit_within` for the model.
#[cfg(test)]
#[must_use]
pub(crate) fn range_fit(elig: &RangeEligibility, inputs: &FitInputs) -> Option<RangeFit> {
    range_fit_from(elig, inputs, 1)
}

/// The range fit at index strides of at least `min_stride` (a power of two)
/// and at most [`MAX_INDEX_STRIDE`]: how [`RunIndexPolicy`] fits.
#[must_use]
pub(crate) fn range_fit_from(
    elig: &RangeEligibility,
    inputs: &FitInputs,
    min_stride: u32,
) -> Option<RangeFit> {
    range_fit_within(elig, inputs, min_stride.max(1), MAX_INDEX_STRIDE)
}

/// The range fit at exactly index stride `stride`: the merge decision fits
/// the runs' indexes as they are, so the published `fit.stride` (and with it
/// the boundary bytes and amplification) is the indexes' own.
#[must_use]
pub(crate) fn range_fit_at(
    elig: &RangeEligibility,
    inputs: &FitInputs,
    stride: u32,
) -> Option<RangeFit> {
    let stride = stride.max(1);
    range_fit_within(elig, inputs, stride, stride)
}

/// Size the range merge over `inputs` within `elig`'s budget
/// ([`RangeEligibility::budget`]), the run indexes within its
/// [`RangeEligibility::index_cap`].
///
/// With `k` runs, `p` phase-2 threads, `w` pool threads and `F` =
/// [`FitInputs::frame_bytes`] ([`FRAME_RAW_BYTES`], or the largest frame
/// when a long record made one larger), every `stride` `N` (powers of two from `min_stride`
/// to `max_stride`, stopping once the index `I(N)` is under 1 % of the
/// budget) and `window` in `[max(2p, 4), p]` is tried, with `started = min(window,
/// max(w, p) + 1)` the ranges that hold per-run cursors (each `MergeWork`
/// worker runs at most one claimed range at a time, and the merge's own
/// thread one more):
///
/// ```text
/// B        = k x N x F                               (boundary bytes)
/// cursors  = started x k x 2F
/// avail    = budget - I(N) - cursors                  (skipped if <= 0)
/// T_target = 20 x B                                   (A = 1.05)
/// T_avail  = avail x S / (window x (2S + S_c)) - B    (S = 0: as if S_c = S)
/// T        = max(1 MiB, min(T_target, T_avail))
/// mem      = I(N) + cursors + window x (2 + S_c / S) x (T + B)
/// feasible iff 4 x B <= T (A <= 1.25), mem <= budget and
///              I(N) <= index_cap
///   (index_cap: merge_memory / INDEX_BUDGET_DIVISOR)
/// ```
///
/// A live range spans up to `T + B` raw bytes (its own `T` plus one boundary
/// entry per run): it holds that span compressed (`(T + B) x S_c / S`), its
/// merged output and the batches cut from it (`T + B` each). A stride whose
/// index is over `index_cap` is skipped and sparser ones
/// are still tried. Of the feasible sizings the one
/// with the smallest amplification `1 + kNF / T` wins, then the larger
/// window, then the smaller `mem`. More memory gives every sizing more
/// room, so it should never raise the amplification;
/// `range_fit_is_monotone_in_memory` checks that over sampled inputs, but it
/// is not proved for every input.
fn range_fit_within(
    elig: &RangeEligibility,
    inputs: &FitInputs,
    min_stride: u32,
    max_stride: u32,
) -> Option<RangeFit> {
    let (budget, index_cap) = (elig.budget(), elig.index_cap());
    let phase2 = elig.phase2_threads;
    let pool = elig.pool_threads.max(phase2);
    let k = u128::from(u64::try_from(inputs.k).unwrap_or(u64::MAX));
    let (raw, comp) = if inputs.raw_bytes == 0 {
        (1u128, 1u128)
    } else {
        (u128::from(inputs.raw_bytes), u128::from(inputs.comp_bytes))
    };
    let bytes = u128::from(budget.bytes);
    let frame = u128::from(inputs.frame_bytes());
    let mut best: Option<RangeFit> = None;
    let mut stride = min_stride;
    while stride <= max_stride {
        let index = inputs.index_bytes(stride);
        let boundary = k * u128::from(stride) * frame;
        let index_fits = index <= index_cap;
        let preferred = phase2.saturating_mul(RANGE_WINDOW_PER_WORKER).max(RANGE_MIN_WINDOW);
        for window in [preferred, phase2] {
            if window == 0 || !index_fits {
                continue;
            }
            let w = window as u128;
            // One started range per `MergeWork` worker, plus the one the
            // merge's own thread claims and runs beside them.
            let started = window.min(pool.saturating_add(1)) as u128;
            let cursors = started * k * u128::from(RANGE_CURSOR_FRAMES) * frame;
            let Some(avail) = bytes.checked_sub(u128::from(index) + cursors).filter(|&a| a > 0)
            else {
                continue;
            };
            let t_target = u128::from(A_TARGET_RECIP) * boundary;
            let t_avail = (avail * raw / (w * (2 * raw + comp))).saturating_sub(boundary);
            let t = t_target.min(t_avail).max(u128::from(MIN_RANGE_BYTES));
            let span = t + boundary;
            let mem = u128::from(index) + cursors + w * (2 * span + span * comp / raw);
            if u128::from(A_MAX_RECIP) * boundary > t || mem > bytes {
                continue;
            }
            let candidate = RangeFit {
                stride,
                window,
                target_range_bytes: u64::try_from(t).unwrap_or(u64::MAX),
                boundary_bytes: u64::try_from(boundary).unwrap_or(u64::MAX),
                mem: u64::try_from(mem).unwrap_or(u64::MAX),
            };
            if best.is_none_or(|b| fits_better(&candidate, &b)) {
                best = Some(candidate);
            }
        }
        if index_fits && index < budget.bytes / 100 {
            break;
        }
        let Some(next) = stride.checked_mul(2) else { break };
        stride = next;
    }
    best
}

/// Whether `a` beats `b`: smaller amplification (`boundary / target`,
/// compared exactly), then the larger window, then the smaller `mem`.
fn fits_better(a: &RangeFit, b: &RangeFit) -> bool {
    let lhs = u128::from(a.boundary_bytes) * u128::from(b.target_range_bytes);
    let rhs = u128::from(b.boundary_bytes) * u128::from(a.target_range_bytes);
    lhs.cmp(&rhs).reverse().then(a.window.cmp(&b.window)).then(b.mem.cmp(&a.mem)).is_gt()
}

/// The final merge's mode.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MergeMode {
    /// Today's serial k-way merge, and why.
    Serial(SerialReason),
    /// A range-parallel merge: range workers read the spill files themselves.
    Range {
        /// The budget range data is charged to (in place of the serial
        /// read-ahead `R`, which range mode does not spend).
        budget: RangeBudget,
        /// The spill runs' figures the fit was computed from. The fit's
        /// memory model charges the index by the estimate
        /// `inputs.index_bytes(fit.stride)`; the actual figure is
        /// `index_bytes`.
        inputs: FitInputs,
        /// The sizing the merge uses; `fit.stride` is every run index's
        /// stride.
        fit: RangeFit,
        /// Heap bytes the runs' indexes hold (the sum of their
        /// [`RunIndex::heap_bytes`]).
        index_bytes: u64,
    },
}

/// Why the final merge runs serially.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SerialReason {
    /// No range eligibility was configured (e.g. `fgumi merge`).
    NoPolicy,
    /// `FGUMI_SORT_SERIAL_MERGE=1`.
    ForcedByEnv,
    /// Every record stayed in memory.
    NoSpill,
    /// Fewer phase-2 threads than `RANGE_MIN_PHASE2_THREADS`.
    TooFewThreads {
        /// The phase-2 thread count.
        have: usize,
        /// The minimum.
        need: usize,
    },
    /// A spill run carries no run index (a BGZF spill, a queryname sort, or a
    /// consolidated run one of whose inputs had none).
    SlotNotIndexed {
        /// The lowest such run's file id.
        file_id: u32,
    },
    /// No sizing keeps decode amplification within 1.25 over this many runs
    /// inside the range budget (see `range_fit`).
    RangeDoesNotFit {
        /// Spill runs.
        k: usize,
        /// The range budget and its origin.
        budget: RangeBudget,
    },
    /// The final run was spilled, but no spill run carries an index although
    /// the sort indexes its runs: the phase-1 front's run-index policy found
    /// no range fit over the runs spilled so far plus the next and turned
    /// indexing off, so `SpillWrite` dropped them.
    IndexesDropped {
        /// Spill runs.
        k: usize,
        /// The range budget and its origin.
        budget: RangeBudget,
    },
    /// The spill runs' indexes do not share one stride (the phase-1 policy
    /// sparsifies every run alike, so this is a defect; fail closed).
    IndexStridesDiffer {
        /// The smallest stride.
        min: u32,
        /// The largest stride.
        max: u32,
    },
    /// Residual chunks are in memory: the range merge needs the final run
    /// spilled, and the phase-1 front kept it resident (its estimate did not
    /// fit, or indexing was off), so the merge is serial. `SpillWrite` does
    /// not read the run indexes back then, so no index-shaped reason is
    /// checked.
    ResidualInMemory {
        /// Residual chunks carried into the merge.
        chunks: u32,
    },
}

impl fmt::Display for SerialReason {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoPolicy => f.write_str("no policy"),
            Self::ForcedByEnv => f.write_str("forced by FGUMI_SORT_SERIAL_MERGE"),
            Self::NoSpill => f.write_str("no spill runs"),
            Self::TooFewThreads { have, need } => {
                write!(f, "too few phase-2 threads ({have} < {need})")
            }
            Self::SlotNotIndexed { file_id } => write!(f, "spill run {file_id} has no run index"),
            Self::RangeDoesNotFit { k, budget } => {
                write!(f, "no range window fits {k} runs within {budget}")
            }
            Self::IndexesDropped { k, budget } => write!(
                f,
                "no range window fits {k} runs within {budget} (run indexes dropped while \
                 spilling; final run spilled)"
            ),
            Self::IndexStridesDiffer { min, max } => {
                write!(f, "spill run indexes have different strides ({min} and {max})")
            }
            Self::ResidualInMemory { chunks } => write!(
                f,
                "range merge needs the final run spilled but {chunks} residual chunk(s) are in \
                 memory"
            ),
        }
    }
}

/// Decide the final merge's mode from `elig`, the merge's spill `slots` and
/// the number of residual chunks carried into the merge in memory. Checks, in
/// order: no eligibility, the environment override, no spill runs, too few
/// phase-2 threads, then, for a sort that indexes its runs, residual chunks
/// ([`SerialReason::ResidualInMemory`]: the range merge needs the final run
/// spilled, its budget assumes the sort buffer is empty, and `SpillWrite`
/// leaves the run indexes unread beside a resident final run), then no run
/// indexed at all although the sort indexes its runs
/// ([`SerialReason::IndexesDropped`]: the phase-1 policy found no fit), then
/// the lowest `file_id` without a run index (every slot's when the sort
/// writes no indexed runs), then indexes of different strides, then
/// [`range_fit_at`] the indexes' common stride over the slots' exact figures
/// within [`RangeEligibility::budget`]. Only coordinate and template zstd
/// runs carry an index, so its presence implies an order the range merge
/// supports.
///
/// The fit is computed from `FitInputs { k: slots, frames: sum of covered
/// frames, raw: sum of index raw bytes, comp: sum of slot lengths,
/// max_frame_raw: the indexes' largest frame }` at the
/// indexes' own stride (the decision never re-sparsifies, so `fit.stride`
/// equals every index's [`RunIndex::stride`]) and published in
/// [`MergeMode::Range`] with the indexes' exact heap bytes, so the merge uses
/// exactly what was checked.
#[must_use]
pub fn decide_merge_mode(
    elig: Option<&RangeEligibility>,
    slots: &[Arc<SortMergeSlot>],
    memory_chunk_count: u32,
) -> MergeMode {
    let Some(elig) = elig else {
        return MergeMode::Serial(SerialReason::NoPolicy);
    };
    if elig.force_serial {
        return MergeMode::Serial(SerialReason::ForcedByEnv);
    }
    if slots.is_empty() {
        return MergeMode::Serial(SerialReason::NoSpill);
    }
    if elig.phase2_threads < RANGE_MIN_PHASE2_THREADS {
        return MergeMode::Serial(SerialReason::TooFewThreads {
            have: elig.phase2_threads,
            need: RANGE_MIN_PHASE2_THREADS,
        });
    }
    if elig.indexable && memory_chunk_count > 0 {
        return MergeMode::Serial(SerialReason::ResidualInMemory { chunks: memory_chunk_count });
    }
    let budget = elig.budget();
    if elig.indexable && slots.iter().all(|s| s.run_index().is_none()) {
        return MergeMode::Serial(SerialReason::IndexesDropped { k: slots.len(), budget });
    }
    if let Some(file_id) =
        slots.iter().filter(|s| !elig.indexable || s.run_index().is_none()).map(|s| s.file_id).min()
    {
        return MergeMode::Serial(SerialReason::SlotNotIndexed { file_id });
    }
    let indexes = || slots.iter().filter_map(|s| s.run_index());
    let (min, max) =
        indexes().map(|ix| ix.stride()).fold((u32::MAX, 0), |(lo, hi), s| (lo.min(s), hi.max(s)));
    if min != max {
        return MergeMode::Serial(SerialReason::IndexStridesDiffer { min, max });
    }
    let inputs = FitInputs {
        k: slots.len(),
        frames: indexes().map(|ix| ix.covered_frames()).sum(),
        key_width: indexes().next().map_or(0, |ix| ix.key_width()),
        raw_bytes: indexes().map(|ix| ix.total_raw_bytes()).sum(),
        comp_bytes: slots.iter().map(|s| s.len()).sum(),
        max_frame_raw: indexes().map(|ix| u64::from(ix.max_frame_raw())).max().unwrap_or(0),
    };
    match range_fit_at(elig, &inputs, min) {
        None => MergeMode::Serial(SerialReason::RangeDoesNotFit { k: inputs.k, budget }),
        Some(fit) => MergeMode::Range {
            budget,
            inputs,
            fit,
            index_bytes: indexes().map(|ix| ix.heap_bytes() as u64).sum(),
        },
    }
}

#[cfg(test)]
mod tests {
    use fgumi_sort::QuerynameComparator;
    use proptest::prelude::*;
    use rstest::rstest;

    use super::*;

    const MIB: u64 = 1 << 20;

    /// An eligibility with indexed runs, not forced serial.
    fn el(phase2: usize, pool: usize, merge_memory: u64, serial_r: u64) -> RangeEligibility {
        RangeEligibility {
            pool_threads: pool,
            ..RangeEligibility::for_test(phase2, merge_memory, serial_r)
        }
    }

    // Columns: (k, frames, key width, raw bytes, compressed bytes), (merge
    // memory, R, phase-2 threads, pool threads), expected fit.
    #[rstest]
    #[case::wgs_coordinate((18, 5_700_000, 8, 330_000_000_000, 115_500_000_000), (20_971_520_000, 512 << 20, 32, 32), Some(RangeFit { stride: 1, window: 64, target_range_bytes: 23_592_960, boundary_bytes: 1_179_648, mem: 3_940_456_960 }))]
    #[case::wgs_template((19, 5_900_000, 40, 330_000_000_000, 115_500_000_000), (20_971_520_000, 512 << 20, 32, 32), Some(RangeFit { stride: 1, window: 64, target_range_bytes: 24_903_680, boundary_bytes: 1_245_184, mem: 4_345_371_264 }))]
    // Window 64 would need A = 1.28 > A_max.
    #[case::low_memory_k90_takes_window_p((90, 28_000_000, 8, 1_650_000_000_000, 577_500_000_000), (4_194_304_000, 256 << 20, 32, 32), Some(RangeFit { stride: 1, window: 32, target_range_bytes: 35_921_130, boundary_bytes: 5_898_240, mem: 4_194_303_968 }))]
    // The `R` floor is the budget; window 16 only reaches A = 1.077.
    #[case::tiny_sort_reaches_target_at_window_p((2, 40, 8, 2_500_000, 1_000_000), (4 << 20, 64 << 20, 8, 8), Some(RangeFit { stride: 1, window: 8, target_range_bytes: 2_621_440, boundary_bytes: 131_072, mem: 54_946_336 }))]
    // Strides 1 and 2 hold an index (683.6 / 341.8 MB) over a quarter of the
    // 800 MiB merge memory and are skipped; stride 4 reaches A = 1.138 at
    // window 32: mem = I(4) 170_898_616 + cursors 16_777_216 + 32 x 2.35 x
    // (T + B).
    #[case::index_dominated_sparsifies((4, 12_207_031, 40, 800_000_000_000, 280_000_000_000), (800 << 20, 64 << 20, 32, 32), Some(RangeFit { stride: 4, window: 32, target_range_bytes: 7_610_798, boundary_bytes: 1_048_576, mem: 838_860_728 }))]
    #[case::k200_infeasible((200, 60_000_000, 40, 3_900_000_000_000, 1_365_000_000_000), (4_194_304_000, 256 << 20, 32, 32), None)]
    #[case::pool_wider_than_phase2((15, 1_000_000, 8, 60_000_000_000, 21_000_000_000), (512 << 20, 512 << 20, 4, 32), Some(RangeFit { stride: 1, window: 8, target_range_bytes: 19_660_800, boundary_bytes: 983_040, mem: 427_832_832 }))]
    #[case::pool_equals_phase2((15, 1_000_000, 8, 60_000_000_000, 21_000_000_000), (512 << 20, 512 << 20, 4, 4), Some(RangeFit { stride: 1, window: 8, target_range_bytes: 19_660_800, boundary_bytes: 983_040, mem: 421_934_592 }))]
    // At one and two phase-2 threads the window floor (4) replaces 2p; the
    // memory model charges the four live ranges: mem = I(1) 24_000_000 +
    // cursors (started = min(4, p + 1) ranges) + 4 x 2.35 x (T + B).
    #[case::window_floor_at_p1((15, 1_000_000, 8, 60_000_000_000, 21_000_000_000), (512 << 20, 512 << 20, 1, 1), Some(RangeFit { stride: 1, window: 4, target_range_bytes: 19_660_800, boundary_bytes: 983_040, mem: 221_984_256 }))]
    #[case::window_floor_at_p2((15, 1_000_000, 8, 60_000_000_000, 21_000_000_000), (512 << 20, 512 << 20, 2, 2), Some(RangeFit { stride: 1, window: 4, target_range_bytes: 19_660_800, boundary_bytes: 983_040, mem: 223_950_336 }))]
    // Too little memory for four live ranges at p = 1 within A <= 1.25: the
    // one-per-worker fallback (window 1, stride 2) still fits.
    #[case::window_floor_falls_back_to_p((15, 1_000_000, 8, 60_000_000_000, 21_000_000_000), (64 << 20, 64 << 20, 1, 1), Some(RangeFit { stride: 2, window: 1, target_range_bytes: 20_647_788, boundary_bytes: 1_966_080, mem: 67_108_861 }))]
    // Stride 1's index (200_001 x 56 = 11.2 MB) is over a quarter of the
    // 32 MiB merge memory and is skipped although the 64 MiB `R` floor would
    // hold it; stride 2's (100_001 x 56) is under it.
    #[case::index_over_a_quarter_of_memory_sparsifies((2, 200_000, 40, 10_000_000, 10_000_000), (32 << 20, 64 << 20, 8, 8), Some(RangeFit { stride: 2, window: 8, target_range_bytes: 2_213_341, boundary_bytes: 262_144, mem: 67_108_848 }))]
    // At WGS sizes the index is a few percent of `-m`: the bound never binds.
    #[case::quarter_bound_never_binds_at_wgs((18, 5_700_000, 8, 330_000_000_000, 115_500_000_000), (20_971_520_000, 512 << 20, 32, 32), Some(RangeFit { stride: 1, window: 64, target_range_bytes: 23_592_960, boundary_bytes: 1_179_648, mem: 3_940_456_960 }))]
    fn range_fit_table(
        #[case] (k, frames, key_width, raw_bytes, comp_bytes): (usize, u64, u16, u64, u64),
        #[case] (merge_memory, serial_r, phase2, pool): (u64, u64, usize, usize),
        #[case] expected: Option<RangeFit>,
    ) {
        let inputs = FitInputs { k, frames, key_width, raw_bytes, comp_bytes, max_frame_raw: 0 };
        let elig = el(phase2, pool, merge_memory, serial_r);
        assert_eq!(range_fit(&elig, &inputs), expected);
    }

    proptest! {
        /// More memory never makes the fit worse or infeasible, and every fit
        /// stays within its budget and within `A_max`.
        #[test]
        fn range_fit_is_monotone_in_memory(
            k in 2usize..64,
            frames in 10u64..10_000_000,
            key_width in prop::sample::select(vec![8u16, 24, 32, 40]),
            m1 in MIB..(64u64 << 30),
            extra in 0u64..(64u64 << 30),
        ) {
            let m2 = m1 + extra;
            // About 60 KB raw per frame, compressed to 35 %.
            let raw_bytes = frames * 60_000;
            let inputs = FitInputs {
                k, frames, key_width, raw_bytes, comp_bytes: raw_bytes * 35 / 100, max_frame_raw: 0,
            };
            // `R = clamp(M / 16, 64 MiB, 512 MiB)`, monotone in `M`.
            let r = |m: u64| (m / 16).clamp(64 * MIB, 512 * MIB);
            let f1 = range_fit(&el(32, 32, m1, r(m1)), &inputs);
            let f2 = range_fit(&el(32, 32, m2, r(m2)), &inputs);
            if f1.is_some() {
                prop_assert!(f2.is_some(), "fit at {m1} but none at {m2}");
            }
            if let (Some(a), Some(b)) = (f1, f2) {
                prop_assert!(a.amplification() >= b.amplification(), "{a:?} vs {b:?}");
            }
            for (m, fit) in [(m1, f1), (m2, f2)] {
                if let Some(fit) = fit {
                    prop_assert!(fit.mem <= m.max(r(m)), "{fit:?} over {m}");
                    prop_assert!(
                        inputs.index_bytes(fit.stride) <= m / INDEX_BUDGET_DIVISOR,
                        "{fit:?}: index over a quarter of {m}"
                    );
                    prop_assert!(fit.amplification() <= 1.25 + 1e-9, "{fit:?}");
                }
            }
        }
    }

    /// A frame above `F` (a record longer than one frame) is charged for every
    /// boundary entry and cursor: the boundary bytes follow it, and a frame
    /// large enough fails the fit that `F` passes. Columns: the largest
    /// frame, the expected boundary bytes (`None`: no fit).
    #[rstest]
    #[case::within_f(60_000, Some(18 * 65_536))]
    #[case::at_f(65_536, Some(18 * 65_536))]
    #[case::long_record(1 << 20, Some(18 << 20))]
    #[case::too_long_to_fit(64 << 20, None)]
    fn the_fit_charges_the_largest_frame(
        #[case] max_frame_raw: u64,
        #[case] boundary: Option<u64>,
    ) {
        let inputs = FitInputs { max_frame_raw, ..wgs_inputs(18) };
        let fit = range_fit_at(&el(32, 32, 20_971_520_000, 512 * MIB), &inputs, 1);
        assert_eq!(fit.map(|f| f.boundary_bytes), boundary, "{fit:?}");
    }

    /// Every run keeps its first index entry, so sparsifying cannot shrink a
    /// one-frame run; a dense index is one entry per frame. Columns: k,
    /// frames, key width, stride, expected bytes.
    #[rstest]
    #[case::one_frame_runs(3, 3, 40, 4, 3 * 56)]
    #[case::runs_keep_their_first_entry(3, 11, 40, 4, 5 * 56)]
    #[case::dense_wgs(18, 5_700_000, 8, 1, 5_700_000 * 24)]
    #[case::fewer_frames_than_runs(3, 0, 8, 1, 3 * 24)]
    fn index_bytes_keeps_one_entry_per_run(
        #[case] k: usize,
        #[case] frames: u64,
        #[case] key_width: u16,
        #[case] stride: u32,
        #[case] expected: u64,
    ) {
        let inputs =
            FitInputs { k, frames, key_width, raw_bytes: 0, comp_bytes: 0, max_frame_raw: 0 };
        assert_eq!(inputs.index_bytes(stride), expected);
    }

    /// The WGS coordinate inputs scaled to `k` runs.
    fn wgs_inputs(k: usize) -> FitInputs {
        let k64 = k as u64;
        FitInputs {
            k,
            frames: 5_700_000 * k64 / 18,
            key_width: 8,
            raw_bytes: 330_000_000_000 * k64 / 18,
            comp_bytes: 115_500_000_000 * k64 / 18,
            max_frame_raw: 0,
        }
    }

    /// Only an indexable, not forced-serial sort with enough phase-2 threads
    /// can range (and so indexes its runs). Columns: indexable, forced
    /// serial, phase-2 threads, expected.
    #[rstest]
    #[case::eligible(true, false, 4, true)]
    #[case::one_thread(true, false, 1, true)]
    #[case::not_indexable(false, false, 32, false)]
    #[case::forced(true, true, 32, false)]
    #[case::no_phase2_threads(true, false, 0, false)]
    fn can_range_table(
        #[case] indexable: bool,
        #[case] force_serial: bool,
        #[case] phase2: usize,
        #[case] expected: bool,
    ) {
        let elig = RangeEligibility { indexable, force_serial, ..el(phase2, phase2, MIB, MIB) };
        assert_eq!(elig.can_range(), expected);
    }

    /// A sort's eligibility indexes coordinate and template-coordinate zstd
    /// runs only, and floors its budget at the serial read-ahead `R`.
    /// Columns: order, codec, merge memory, expected (indexable, R).
    #[rstest]
    #[case::coordinate_zstd(SortOrder::Coordinate, SpillCodec::Zstd, 20_000 * MIB, (true, 512 * MIB))]
    #[case::template_zstd(SortOrder::TemplateCoordinate, SpillCodec::Zstd, 4 * MIB, (true, 64 * MIB))]
    #[case::coordinate_bgzf(SortOrder::Coordinate, SpillCodec::Bgzf, 2048 * MIB, (false, 128 * MIB))]
    #[case::queryname(SortOrder::Queryname(QuerynameComparator::Lexicographic), SpillCodec::Zstd, 2048 * MIB, (false, 128 * MIB))]
    fn for_sort_table(
        #[case] order: SortOrder,
        #[case] codec: SpillCodec,
        #[case] merge_memory: u64,
        #[case] (indexable, serial_r): (bool, u64),
    ) {
        let elig = RangeEligibility::for_sort(order, codec, 8, 8, false, merge_memory);
        assert_eq!(
            elig,
            RangeEligibility {
                indexable,
                phase2_threads: 8,
                pool_threads: 8,
                force_serial: false,
                merge_memory,
                serial_r
            }
        );
    }

    // Columns: indexable, forced serial, phase-2 threads, k, expected.
    #[rstest]
    #[case::eligible(true, false, 32, 2, true)]
    #[case::one_run_only(true, false, 32, 1, false)]
    #[case::not_indexable(false, false, 32, 2, false)]
    #[case::forced(true, true, 32, 2, false)]
    #[case::two_threads(true, false, 2, 2, true)]
    #[case::no_phase2_threads(true, false, 0, 2, false)]
    fn spills_final_run_table(
        #[case] indexable: bool,
        #[case] force_serial: bool,
        #[case] phase2: usize,
        #[case] k: usize,
        #[case] expected: bool,
    ) {
        let elig = RangeEligibility {
            indexable,
            force_serial,
            ..el(phase2, phase2, 20_971_520_000, 512 * MIB)
        };
        assert_eq!(elig.spills_final_run(&wgs_inputs(k)), expected);
    }

    /// `n` slots in descending `file_id` order; when `all_indexed` is false,
    /// every slot but file 0 lacks an index, so the first unindexed slot in
    /// the slice (the highest id) is not the one reported.
    fn slots(n: u32, all_indexed: bool, index: &RunIndex) -> Vec<Arc<SortMergeSlot>> {
        (0..n)
            .rev()
            .map(|id| {
                let s = SortMergeSlot::for_test(id, SpillCodec::Zstd);
                if all_indexed || id == 0 {
                    s.attach_run_index(Arc::new(index.clone()));
                }
                Arc::new(s)
            })
            .collect()
    }

    /// The fit inputs of `k` empty indexed runs.
    fn empty_inputs(k: usize) -> FitInputs {
        FitInputs { k, frames: 0, key_width: 8, raw_bytes: 0, comp_bytes: 0, max_frame_raw: 0 }
    }

    /// Range mode over `k` empty runs within `budget`.
    fn range(
        k: usize,
        budget: RangeBudget,
        window: usize,
        target_range_bytes: u64,
        mem: u64,
    ) -> MergeMode {
        MergeMode::Range {
            budget,
            inputs: empty_inputs(k),
            fit: RangeFit {
                stride: 1,
                window,
                target_range_bytes,
                boundary_bytes: k as u64 * FRAME_RAW_BYTES,
                mem,
            },
            // An empty index holds its one trailing offset.
            index_bytes: 8 * k as u64,
        }
    }

    const BUDGET: u64 = 512 * MIB;

    /// An eligibility over `BUDGET` with `R = BUDGET` and a pool as wide as
    /// phase 2 (the default config, no `--merge-threads`).
    fn p(phase2: usize) -> RangeEligibility {
        el(phase2, phase2, BUDGET, BUDGET)
    }

    // Columns: eligibility, slot count, all indexed, residual chunks,
    // expected. The runs are empty (no frames, no bytes), so the fit
    // reaches `A = 1.05` (`T = 20 x k x 64 KiB`) unless the budget binds,
    // and a live range holds `3 T` (the `S_c = S` worst case).
    #[rstest]
    #[case::no_policy(None, 3, true, 0, MergeMode::Serial(SerialReason::NoPolicy))]
    #[case::forced(
        Some(RangeEligibility { force_serial: true, ..p(8) }),
        3,
        true,
        0,
        MergeMode::Serial(SerialReason::ForcedByEnv)
    )]
    #[case::no_spill(Some(p(8)), 0, true, 0, MergeMode::Serial(SerialReason::NoSpill))]
    #[case::no_phase2_threads(
        Some(p(0)),
        3,
        true,
        0,
        MergeMode::Serial(SerialReason::TooFewThreads { have: 0, need: 1 })
    )]
    #[case::not_indexed(
        Some(p(8)),
        3,
        false,
        0,
        MergeMode::Serial(SerialReason::SlotNotIndexed { file_id: 1 })
    )]
    #[case::not_indexable(
        Some(RangeEligibility { indexable: false, ..p(8) }),
        3,
        true,
        0,
        MergeMode::Serial(SerialReason::SlotNotIndexed { file_id: 0 })
    )]
    // A sort that indexes nothing reports its unindexed runs, residual or
    // not: the residual reason is for a sort that could range.
    #[case::not_indexable_beside_a_residual(
        Some(RangeEligibility { indexable: false, ..p(8) }),
        3,
        true,
        1,
        MergeMode::Serial(SerialReason::SlotNotIndexed { file_id: 0 })
    )]
    // 3 runs, 8 workers: window 16 at the target, 3 MiB of cursors.
    #[case::fits_at_two_per_worker(
        Some(p(8)),
        3,
        true,
        0,
        range(3, RangeBudget::new(BUDGET, BUDGET), 16, 3_932_160, 201_719_880)
    )]
    #[case::fits_k15_four_workers(
        Some(p(4)),
        15,
        true,
        0,
        range(15, RangeBudget::new(BUDGET, BUDGET), 8, 19_660_800, 505_282_920)
    )]
    // 15 runs, 32 workers: window 64 would need `A > 1.25`; window 32 spends
    // the budget at `A = 1.249`.
    #[case::fits_only_at_one_per_worker(
        Some(p(32)),
        15,
        true,
        0,
        range(15, RangeBudget::new(BUDGET, BUDGET), 32, 3_954_001, 536_870_856)
    )]
    // 200 runs, 32 workers: the cursors alone (800 MiB at window 32) exceed
    // 128 MiB.
    #[case::does_not_fit_k200(
        Some(el(32, 32, 128 * MIB, 128 * MIB)),
        200,
        true,
        0,
        MergeMode::Serial(SerialReason::RangeDoesNotFit {
            k: 200,
            budget: RangeBudget::new(128 * MIB, 128 * MIB),
        })
    )]
    // 100 MiB, 15 runs, 4 of 32 workers: window 8 would need `A > 1.25`.
    #[case::merge_threads_below_pool_shrinks(
        Some(el(4, 32, 100 * MIB, 100 * MIB)),
        15,
        true,
        0,
        range(15, RangeBudget::new(100 * MIB, 100 * MIB), 4, 7_099_703, 104_857_596)
    )]
    // 15 runs at 32 workers fit a 4 GiB merge memory with no residual ...
    #[case::fits_without_a_residual(
        Some(el(32, 32, 4096 * MIB, 64 * MIB)),
        15,
        true,
        0,
        range(15, RangeBudget::new(4096 * MIB, 64 * MIB), 64, 19_660_800, 4_028_498_280)
    )]
    // ... but residual chunks in memory fail closed: the range merge needs
    // the final run spilled.
    #[case::residual_in_memory(
        Some(el(32, 32, 4096 * MIB, 64 * MIB)),
        15,
        true,
        3,
        MergeMode::Serial(SerialReason::ResidualInMemory { chunks: 3 })
    )]
    // The residual is reported before the fit, which a resident final run
    // makes moot.
    #[case::residual_before_the_fit(
        Some(el(32, 32, 128 * MIB, 128 * MIB)),
        200,
        true,
        1,
        MergeMode::Serial(SerialReason::ResidualInMemory { chunks: 1 })
    )]
    // ... and before the indexes, which `SpillWrite` leaves unread beside a
    // resident final run.
    #[case::residual_with_unread_indexes(
        Some(el(32, 32, 4096 * MIB, 64 * MIB)),
        15,
        false,
        1,
        MergeMode::Serial(SerialReason::ResidualInMemory { chunks: 1 })
    )]
    fn decide_merge_mode_table(
        #[case] elig: Option<RangeEligibility>,
        #[case] slot_count: u32,
        #[case] all_indexed: bool,
        #[case] memory_chunk_count: u32,
        #[case] expected: MergeMode,
    ) {
        let slots = slots(slot_count, all_indexed, &RunIndex::new(8));
        assert_eq!(decide_merge_mode(elig.as_ref(), &slots, memory_chunk_count), expected);
    }

    /// Slots with indexes of `strides` (`None`: no index), `file_id` = position.
    fn strided_slots(strides: &[Option<u32>]) -> Vec<Arc<SortMergeSlot>> {
        (0u32..)
            .zip(strides)
            .map(|(id, stride)| {
                let s = SortMergeSlot::for_test(id, SpillCodec::Zstd);
                if let Some(stride) = stride {
                    let mut index = RunIndex::new(8);
                    for f in 0..8u8 {
                        index.push_frame(100, 1000, 1, &[f; 8]);
                    }
                    index.sparsify(*stride).unwrap();
                    s.attach_run_index(Arc::new(index));
                }
                Arc::new(s)
            })
            .collect()
    }

    /// The index-shaped serial reasons: no index at all although the sort
    /// indexes its runs (the phase-1 policy dropped them, reported as a
    /// failed fit), some run unindexed, or indexes of different strides.
    /// Columns: indexable, slot strides, expected reason.
    #[rstest]
    #[case::all_dropped(true, vec![None, None, None], SerialReason::IndexesDropped { k: 3, budget: RangeBudget::new(BUDGET, BUDGET) })]
    #[case::not_indexable(false, vec![None, None, None], SerialReason::SlotNotIndexed { file_id: 0 })]
    #[case::one_unindexed(true, vec![Some(1), None, Some(1)], SerialReason::SlotNotIndexed { file_id: 1 })]
    #[case::strides_differ(true, vec![Some(2), Some(1), Some(4)], SerialReason::IndexStridesDiffer { min: 1, max: 4 })]
    fn decide_merge_mode_index_reasons(
        #[case] indexable: bool,
        #[case] strides: Vec<Option<u32>>,
        #[case] expected: SerialReason,
    ) {
        let elig = RangeEligibility { indexable, ..p(8) };
        let slots = strided_slots(&strides);
        assert_eq!(decide_merge_mode(Some(&elig), &slots, 0), MergeMode::Serial(expected));
    }

    /// The decision fits at the indexes' common stride and publishes their
    /// exact heap bytes, whichever stride a free fit would pick.
    #[rstest]
    #[case::dense(1)]
    #[case::stride_two(2)]
    #[case::stride_four(4)]
    fn decide_merge_mode_fits_at_the_index_stride(#[case] stride: u32) {
        let slots = strided_slots(&[Some(stride); 3]);
        let elig = p(8);
        let MergeMode::Range { inputs, fit, index_bytes, .. } =
            decide_merge_mode(Some(&elig), &slots, 0)
        else {
            panic!("range mode");
        };
        assert_eq!(fit.stride, stride);
        assert_eq!(Some(fit), range_fit_at(&elig, &inputs, stride));
        let heap: u64 = slots.iter().map(|s| s.run_index().unwrap().heap_bytes() as u64).sum();
        assert_eq!(index_bytes, heap);
    }

    /// `range_fit_from` never goes below its minimum stride, and at stride 1
    /// it is `range_fit`.
    #[test]
    fn range_fit_from_respects_the_minimum_stride() {
        let inputs = wgs_inputs(18);
        let elig = el(32, 32, 20_971_520_000, 512 * MIB);
        assert_eq!(range_fit_from(&elig, &inputs, 1), range_fit(&elig, &inputs));
        let fit = range_fit_from(&elig, &inputs, 8).expect("a sparser fit");
        assert_eq!(fit.stride, 8, "WGS coordinate is dense at best; the floor binds");
    }

    /// The decision sums every run's covered frames and raw bytes and every
    /// slot's length, and its fit is `range_fit` on those figures.
    #[test]
    fn decide_merge_mode_sums_the_slot_figures() {
        let mut index = RunIndex::new(8);
        index.push_frame(100, 200, 3, &[7; 8]);
        index.push_frame(50, 300, 2, &[8; 8]);
        // Slots of 1,000, 2,000 and 3,000 file bytes: the compressed bytes are
        // their lengths summed, not the indexes'.
        let file = Arc::new(tempfile::tempfile().expect("tempfile"));
        let slots: Vec<Arc<SortMergeSlot>> = (0..3u32)
            .map(|id| {
                let len = 1_000 * u64::from(id + 1);
                let s = SortMergeSlot::new(id, Arc::clone(&file) as _, 0, len, SpillCodec::Zstd);
                s.attach_run_index(Arc::new(index.clone()));
                Arc::new(s)
            })
            .collect();
        let index_bytes = slots.iter().map(|s| s.run_index().unwrap().heap_bytes() as u64).sum();
        let elig = p(8);
        let got = decide_merge_mode(Some(&elig), &slots, 0);
        let expected_inputs = FitInputs {
            k: 3,
            frames: 6,
            key_width: 8,
            raw_bytes: 1500,
            comp_bytes: 6_000,
            max_frame_raw: 300,
        };
        assert_eq!(
            got,
            MergeMode::Range {
                budget: elig.budget(),
                inputs: expected_inputs,
                fit: range_fit(&elig, &expected_inputs).expect("a fit"),
                index_bytes,
            }
        );
    }

    #[rstest]
    #[case::no_policy(SerialReason::NoPolicy, "no policy")]
    #[case::forced(SerialReason::ForcedByEnv, "forced by FGUMI_SORT_SERIAL_MERGE")]
    #[case::no_spill(SerialReason::NoSpill, "no spill runs")]
    #[case::too_few_threads(
        SerialReason::TooFewThreads { have: 0, need: 1 },
        "too few phase-2 threads (0 < 1)"
    )]
    #[case::not_indexed(SerialReason::SlotNotIndexed { file_id: 5 }, "spill run 5 has no run index")]
    #[case::does_not_fit(
        SerialReason::RangeDoesNotFit { k: 200, budget: RangeBudget::new(300 * MIB, BUDGET) },
        "no range window fits 200 runs within range budget 512 MiB (merge memory 300 MiB, floor \
         R 512 MiB)"
    )]
    #[case::indexes_dropped(
        SerialReason::IndexesDropped { k: 9, budget: RangeBudget::new(64 * MIB, 64 * MIB) },
        "no range window fits 9 runs within range budget 64 MiB (merge memory 64 MiB, floor R \
         64 MiB) (run indexes dropped while spilling; final run spilled)"
    )]
    #[case::strides_differ(
        SerialReason::IndexStridesDiffer { min: 1, max: 2 },
        "spill run indexes have different strides (1 and 2)"
    )]
    #[case::residual_in_memory(
        SerialReason::ResidualInMemory { chunks: 2 },
        "range merge needs the final run spilled but 2 residual chunk(s) are in memory"
    )]
    fn serial_reason_display(#[case] reason: SerialReason, #[case] expected: &str) {
        assert_eq!(reason.to_string(), expected);
    }

    /// The range budget: `max(R, merge_memory)`. Columns: merge memory, R,
    /// expected budget.
    #[rstest]
    #[case::merge_memory_wins(20_000 << 20, 512 << 20, 20_000 << 20)]
    #[case::floor_r_wins(4 << 20, 64 << 20, 64 << 20)]
    fn range_budget_table(#[case] merge_memory: u64, #[case] serial_r: u64, #[case] expected: u64) {
        assert_eq!(
            RangeBudget::new(merge_memory, serial_r),
            RangeBudget { bytes: expected, merge_memory, serial_r }
        );
    }

    #[test]
    fn range_budget_display_names_its_origin() {
        assert_eq!(
            RangeBudget::new(4096 * MIB, 512 * MIB).to_string(),
            "range budget 4096 MiB (merge memory 4096 MiB, floor R 512 MiB)"
        );
    }

    /// The front's frame bound, `B` = 65,280: `ceil(2S / B)` for bytes whose
    /// largest record is unknown, the smaller of that and `ceil(S / (B - m))`
    /// when it is `m < B`, plus one frame per run; raw and compressed bytes
    /// are every framed byte. Columns: k, kind, known framed bytes, their
    /// largest record, unbounded framed bytes, and the expected (frames, key
    /// width, largest frame).
    #[rstest]
    // ceil(1_010_000 / 64_876) = 16 (tighter than 31), + 2.
    #[case::short_records_tight(2, SpillKeyKind::Coordinate, 1_010_000, Some(404), 0, (18, 8, 404))]
    // ceil(2_020_000 / 65_280) = 31, + 2.
    #[case::largest_unknown(2, SpillKeyKind::Coordinate, 1_010_000, None, 0, (33, 8, 0))]
    #[case::unbounded(2, SpillKeyKind::TemplateK40, 0, None, 1_010_000, (33, 40, 0))]
    // 16 known + 31 unbounded + 2.
    #[case::known_and_unbounded(2, SpillKeyKind::Coordinate, 1_010_000, Some(404), 1_010_000, (49, 8, 404))]
    // 40 KB records: ceil(1_010_000 / 25_280) = 40 is looser than 31.
    #[case::long_records(2, SpillKeyKind::Coordinate, 1_010_000, Some(40_000), 0, (33, 8, 40_000))]
    // A record over one block is a frame of its own: the 2x form holds, and
    // the fit charges that frame's size for every frame.
    #[case::oversized_record(2, SpillKeyKind::Coordinate, 1_010_000, Some(70_000), 0, (33, 8, 70_000))]
    #[case::empty(1, SpillKeyKind::Coordinate, 0, None, 0, (1, 8, 0))]
    fn estimated_inputs_table(
        #[case] k: usize,
        #[case] kind: SpillKeyKind,
        #[case] framed: u64,
        #[case] max_record: Option<u64>,
        #[case] unbounded: u64,
        #[case] (frames, key_width, max_frame_raw): (u64, u16, u64),
    ) {
        let raw_bytes = framed + unbounded;
        assert_eq!(
            FitInputs::estimated(k, kind, SpillEstimate { framed, max_record, unbounded }),
            FitInputs { k, frames, key_width, raw_bytes, comp_bytes: raw_bytes, max_frame_raw }
        );
    }

    /// Bytes framed from records with their length prefix and a key prefix
    /// per record (template), the length prefix counted once.
    #[rstest]
    #[case::coordinate(Some(SpillKeyKind::Coordinate), 1_000, 10, 1_000)]
    #[case::template_k24(Some(SpillKeyKind::TemplateK24), 1_000, 10, 1_240)]
    #[case::unknown_kind(None, 1_000, 10, 1_000)]
    fn framed_bytes_add_the_key_prefix(
        #[case] kind: Option<SpillKeyKind>,
        #[case] len_prefixed: u64,
        #[case] records: u64,
        #[case] expected: u64,
    ) {
        assert_eq!(SpillEstimate::framed_bytes(kind, len_prefixed, records), expected);
    }

    /// The bound holds against the spill path's own cut rule for records of
    /// every length, long reads included: frame a run of `n` records of
    /// `len` bytes with the gather's rule and compare.
    #[rstest]
    fn frame_bound_holds_for_long_records(
        #[values(100, 4_000, 8_000, 20_000, 33_000, 40_000, 64_000, 70_000)] len: usize,
        #[values(1, 7, 50)] n: usize,
    ) {
        use fgumi_sort::{InMemoryChunk, RawCoordinateKey};
        let chunk = crate::sort::protocol::MemoryChunkErased::Coordinate(
            InMemoryChunk::from_owned_records(
                (0..n).map(|i| (RawCoordinateKey { sort_key: i as u64 }, vec![0u8; len])).collect(),
            ),
        );
        let frames = crate::sort::spill_gather::frame_chunk_into_blocks(
            &chunk,
            usize::try_from(SPILL_FRAME_BYTES).unwrap(),
        )
        .unwrap()
        .len() as u64;
        let framed = (n * (len + 4)) as u64;
        let known = SpillEstimate { framed, max_record: Some(len as u64 + 4), unbounded: 0 };
        let unknown = SpillEstimate { framed: 0, max_record: None, unbounded: framed };
        assert!(known.frame_bound(1) >= frames, "{len} x {n}: {frames} frames");
        assert!(unknown.frame_bound(1) >= frames, "{len} x {n}: {frames} frames");
    }

    /// The policy only ever raises the stride or turns indexing off: once
    /// off it stays off even when a later run's figures would fit, and once
    /// sparsified it never goes back to dense.
    #[test]
    fn index_stride_is_monotone() {
        let unbounded = |b| SpillEstimate { framed: 0, max_record: None, unbounded: b };
        // Off: 16 runs at 32 workers need 64 MiB of cursors, the whole budget.
        let mut policy = RunIndexPolicy::new(Some(el(32, 32, 1_000_000, 64 * MIB)));
        assert_eq!(policy.stride(), 1);
        policy.at_run(16, Some(SpillKeyKind::Coordinate), unbounded(14_400_000));
        assert_eq!(policy.stride(), 0);
        policy.at_run(2, Some(SpillKeyKind::Coordinate), unbounded(1_800_000));
        assert_eq!(policy.stride(), 0, "off for good");
        assert_eq!(
            policy.final_stride(2, Some(SpillKeyKind::Coordinate), unbounded(1_800_000)),
            None
        );
        // Sparsified: two K40 runs over 9 GB hold a dense index (275_738
        // entries x 56 B) over a quarter of the 32 MiB merge memory, so the
        // policy sparsifies to stride 2; a later small run that would fit
        // densely keeps stride 2.
        let mut policy = RunIndexPolicy::new(Some(el(8, 8, 32 * MIB, 64 * MIB)));
        policy.at_run(2, Some(SpillKeyKind::TemplateK40), unbounded(9_000_000_000));
        assert_eq!(policy.stride(), 2);
        policy.at_run(2, Some(SpillKeyKind::TemplateK40), unbounded(1_000_000));
        assert_eq!(policy.stride(), 2);
        assert_eq!(
            policy.final_stride(2, Some(SpillKeyKind::TemplateK40), unbounded(1_000_000)),
            Some(2)
        );
    }

    /// The policy indexes nothing for a sort that cannot merge by range
    /// (forced serial, too few threads, no policy), or for an unindexed key
    /// kind; the first run and a run whose key kind is not yet known keep the
    /// stride. Columns: the eligibility, runs, kind, expected stride.
    #[rstest]
    // One spilled run of 900_000 framed bytes whose largest record is 300
    // (ceil(900_000 / 64_980) = 14 frames) and the next estimated at its
    // size (ceil(1_800_000 / 65_280) = 28), + 2 = 44 entries x 24 B, well
    // under a quarter of the 1 MB merge memory: dense.
    #[case::indexed(Some(el(8, 8, 1_000_000, 64 * MIB)), 2, Some(SpillKeyKind::Coordinate), 1)]
    #[case::first_run(Some(el(8, 8, 1_000_000, 64 * MIB)), 1, Some(SpillKeyKind::Coordinate), 1)]
    #[case::forced_serial(Some(RangeEligibility { force_serial: true, ..el(8, 8, 1_000_000, 64 * MIB) }), 2, Some(SpillKeyKind::Coordinate), 0)]
    #[case::no_phase2_threads(Some(el(0, 0, 1_000_000, 64 * MIB)), 2, Some(SpillKeyKind::Coordinate), 0)]
    #[case::no_policy(None, 2, Some(SpillKeyKind::Coordinate), 0)]
    #[case::queryname(Some(el(8, 8, 1_000_000, 64 * MIB)), 2, Some(SpillKeyKind::QuerynameLex), 0)]
    #[case::kind_unknown(Some(el(8, 8, 1_000_000, 64 * MIB)), 2, None, 1)]
    fn run_index_policy_table(
        #[case] elig: Option<RangeEligibility>,
        #[case] k: usize,
        #[case] kind: Option<SpillKeyKind>,
        #[case] expected: u32,
    ) {
        let mut policy = RunIndexPolicy::new(elig);
        let spill = SpillEstimate { framed: 900_000, max_record: Some(300), unbounded: 900_000 };
        policy.at_run(k, kind, spill);
        assert_eq!(policy.stride(), expected);
    }

    /// The estimate is conservative: 10 runs over 64 MiB at four workers fit
    /// a budget of 31 MiB plus the index when nothing is known compressed
    /// (`S_c = 0` understates `rho`: window 4 at `T` = 2.6 MiB, `A` = 1.24),
    /// but not at `S_c = S` (`A` would exceed 1.25), so the front keeps the
    /// residual where a comp-less estimate would have spilled it.
    #[test]
    fn estimated_fit_inputs_err_toward_the_residual() {
        let raw = 64 * MIB;
        let spill = SpillEstimate { framed: raw, max_record: None, unbounded: 0 };
        let estimated = FitInputs::estimated(10, SpillKeyKind::Coordinate, spill);
        let elig = el(4, 4, 31 * MIB + estimated.index_bytes(1), MIB);
        let comp_less = FitInputs { comp_bytes: 0, ..estimated };
        assert!(range_fit(&elig, &comp_less).is_some(), "S_c = 0 fits");
        assert!(elig.spills_final_run(&comp_less));
        assert_eq!(range_fit(&elig, &estimated), None);
        assert!(!elig.spills_final_run(&estimated));
    }
}
