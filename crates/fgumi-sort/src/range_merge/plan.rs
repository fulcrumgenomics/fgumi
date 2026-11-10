//! The key-range plan of a range-parallel merge.
//!
//! [`RangePlan::new`] takes the merge's leaves (indexed spill slots; range
//! mode runs only when the final run spilled, so there is no in-memory
//! residual), validates them, and chooses bounds `b_0 < b_1 < ... < b_{R-2}`
//! that cut the merge's record order into `R` ranges of roughly
//! `target_range_bytes` raw bytes each. Range `j` is `[b_{j-1}, b_j)`; range 0
//! is unbounded below and the last range unbounded above. The plan then
//! answers, per source and range, which run entries a range merge must read.
//! Bound choice affects only how evenly work is spread, never the merged
//! output.
//!
//! The merge orders records by key, then source (leaf) index, then position
//! in the source. A bound is a [`RangeBound`] in that same order, placed at
//! the start of an index entry, so a bound may fall inside a group of equal
//! keys (every unmapped read of a coordinate sort, a deep pileup) and each
//! range stays near the target however many records share a key.
//!
//! Each range's sampled entries (those whose start lies in the range) hold at
//! most `ceil(total / R) + W` raw bytes, where `total` is the runs' raw bytes,
//! `R = clamp(ceil(total / target), 1, max_ranges)` and `W` the largest index
//! entry. A range also reads one boundary entry of every other run whose
//! records straddle its lower bound, so its raw span is at most
//! `ceil(total / R) + k x W` for `k` runs (at most `target + k x W` when
//! `max_ranges` does not bind).

use std::io;
use std::ops::Range;
use std::sync::Arc;

use crate::keys::RawSortKey;
use crate::merge_slots::SortMergeSlot;
use crate::run_index::RunIndex;

/// Upper bound on the number of ranges a plan may cut.
pub const MAX_RANGES: usize = 1 << 20;

/// How many ranges a plan cuts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RangeSizing {
    /// Raw bytes per range to aim for (0 is treated as 1).
    pub target_range_bytes: u64,
    /// Most ranges to cut (0 is treated as 1; never more than [`MAX_RANGES`]).
    /// Production cuts up to [`MAX_RANGES`] ([`RangeSizing::new`]); tests
    /// lower it to pin a plan's range count.
    pub max_ranges: usize,
}

impl RangeSizing {
    /// Ranges of about `target_range_bytes` raw bytes, at most [`MAX_RANGES`].
    #[must_use]
    pub fn new(target_range_bytes: u64) -> Self {
        Self { target_range_bytes, max_ranges: MAX_RANGES }
    }

    /// The number of ranges `R = clamp(ceil(total / target), 1, max_ranges)`
    /// cut over `total` raw bytes.
    fn range_count(&self, total: u128) -> usize {
        let target = u128::from(self.target_range_bytes.max(1));
        let max_ranges = self.max_ranges.clamp(1, MAX_RANGES);
        let wanted = total.div_ceil(target).max(1);
        usize::try_from(wanted).map_or(max_ranges, |r| r.min(max_ranges))
    }
}

/// One leaf of the merge (an indexed spill run), in the serial merge's leaf
/// order.
pub(crate) struct RangeSource {
    /// The run's slot.
    pub(crate) slot: Arc<SortMergeSlot>,
    /// The run's index (absolute file offsets), read in place.
    pub(crate) index: Arc<RunIndex>,
}

/// A point in the merge's record order: a key, a source (leaf) index and an
/// index entry of that source. Bounds order by key, then source, then entry,
/// as the merge orders records.
///
/// A record with key `k`, in entry `p` of source `r`, sorts below the bound
/// `(key, run, entry)` when `k < key`, or `k == key` and either `r < run` or
/// `r == run` and `p < entry`. A plan only places a bound at the start of an
/// entry of `run` whose first key is `key`, so every record of `run` before
/// that entry has a key of at most `key`, and none after it a smaller one.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct RangeBound<K> {
    /// The key at the bound.
    pub key: K,
    /// The source (leaf) index, in the plan's source order.
    pub run: usize,
    /// The entry of `run`'s index the bound starts.
    pub entry: usize,
}

impl<K: Ord> RangeBound<K> {
    /// Whether a record with `key`, in entry `entry` of source `run`, sorts
    /// below this bound (so a range `[lo, hi)` keeps it when it is not below
    /// `lo` and is below `hi`).
    #[must_use]
    pub fn has_below(&self, key: &K, run: usize, entry: usize) -> bool {
        match key.cmp(&self.key) {
            std::cmp::Ordering::Less => true,
            std::cmp::Ordering::Greater => false,
            std::cmp::Ordering::Equal => (run, entry) < (self.run, self.entry),
        }
    }
}

/// The key ranges of a range-parallel merge and the sources they read.
pub struct RangePlan<K> {
    sources: Vec<RangeSource>,
    splitters: Vec<RangeBound<K>>,
    total_records: u64,
    /// Raw bytes of the run frames: the denominator of
    /// [`RangePlan::planned_decode_amplification`].
    run_raw_bytes: u64,
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// Entries of source `run`'s `index` whose start sorts below `bound` (entry
/// `i` starts at `(key_at(i), run, i)`, which increases with `i`).
fn entries_below<K: RawSortKey>(
    index: &RunIndex,
    run: usize,
    bound: &RangeBound<K>,
) -> io::Result<usize> {
    let (mut lo, mut hi) = (0usize, index.len());
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if bound.has_below(&index.key_at::<K>(mid)?, run, mid) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    Ok(lo)
}

/// Entries of source `run` that can hold a record in `[lo, hi)`: every entry
/// starting in the range, plus the entry before the first of them when `run`
/// is not `lo`'s source (a record at or above `lo` may sit in its tail). On
/// `lo`'s own source the range starts exactly at `lo`'s entry.
///
/// # Errors
///
/// If an entry's first key does not decode as `K`.
fn frame_span_of<K: RawSortKey>(
    index: &RunIndex,
    run: usize,
    lo: Option<&RangeBound<K>>,
    hi: Option<&RangeBound<K>>,
) -> io::Result<Range<usize>> {
    let start = match lo {
        Some(lo) if lo.run == run => entries_below(index, run, lo)?,
        Some(lo) => entries_below(index, run, lo)?.saturating_sub(1),
        None => 0,
    };
    let end = match hi {
        Some(hi) => entries_below(index, run, hi)?,
        None => index.len(),
    };
    Ok(if start >= end { start..start } else { start..end })
}

/// Byte-weighted quantile bounds over `samples` (one per index entry, its
/// start and raw bytes; consumed and sorted).
///
/// Cuts `R = clamp(ceil(total / target), 1, max_ranges)` ranges, where `total`
/// is the summed sample weight: for `r` in `1..R`, the first sample after the
/// previous cut (and never the first sample, which would leave range 0 empty)
/// at which the cumulative weight reaches `total * r / R`. Samples are
/// distinct, so the bounds strictly increase; equal keys may repeat across
/// consecutive bounds. Each range's samples then weigh less than
/// `ceil(total / R)` plus its first sample's weight.
fn splitters_of<K: Ord + Copy>(
    mut samples: Vec<(RangeBound<K>, u64)>,
    sizing: RangeSizing,
) -> Vec<RangeBound<K>> {
    if samples.len() < 2 {
        return Vec::new();
    }
    samples.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    let total: u128 = samples.iter().map(|&(_, w)| u128::from(w)).sum();
    let ranges = sizing.range_count(total);

    let mut splitters = Vec::with_capacity(ranges - 1);
    // The first sample is never a cut: count it before the first threshold.
    let mut cumulative = u128::from(samples[0].1);
    let mut next = samples.iter().skip(1);
    'cuts: for r in 1..ranges {
        let threshold = total * r as u128 / ranges as u128;
        // Advance until the cumulative weight reaches the threshold; the
        // sample that reaches it is the cut.
        let bound = loop {
            match next.next() {
                Some(&(b, w)) => {
                    cumulative += u128::from(w);
                    if cumulative >= threshold {
                        break b;
                    }
                }
                None => break 'cuts,
            }
        };
        splitters.push(bound);
    }
    // The capacity was sized for the wanted ranges; fewer may be cut.
    splitters.shrink_to_fit();
    splitters
}

impl<K: RawSortKey + Default + Copy + Send + Sync + 'static> RangePlan<K> {
    /// Plan ranges over `slots`.
    ///
    /// Sources are ordered as the serial merge orders its leaves: slots sorted
    /// by `file_id`.
    ///
    /// # Errors
    ///
    /// Returns [`io::ErrorKind::InvalidData`] when a slot has no run index,
    /// is a stream or skips decompressed bytes before its first record (a
    /// range task preads frames at indexed offsets from record 0), an index's
    /// key width is not `K`'s serialized size, an index does not start at the
    /// slot's body or runs past the slot's file, a run's entry first keys
    /// decrease, or the index record counts do not sum to `total_records`.
    ///
    /// # Memory
    ///
    /// Choosing the bounds sorts one sample per index entry (its bound and
    /// raw bytes, `size_of::<K>() + 24` bytes), a transient freed before this
    /// returns. The merge decision caps the indexes at a quarter of the merge
    /// memory, and an entry costs `16 + key_width` bytes, so the samples are
    /// at most about 1.5 times the indexes (coordinate: 32 B against 24 B),
    /// under 3/8 of the budget, while no range is live yet and the rest of
    /// the budget is free.
    pub fn new(
        slots: Vec<Arc<SortMergeSlot>>,
        total_records: u64,
        sizing: RangeSizing,
    ) -> io::Result<Self> {
        let mut slots = slots;
        slots.sort_by_key(|s| s.file_id);

        let mut sources = Vec::with_capacity(slots.len());
        // Sized up front: one sample per index entry. Grown by doubling, the
        // vector overshoots by up to its own size while the merge's slots are
        // being built.
        let mut samples: Vec<(RangeBound<K>, u64)> =
            Vec::with_capacity(slots.iter().filter_map(|s| s.run_index()).map(|ix| ix.len()).sum());
        let mut counted_records = 0u64;
        let mut run_raw_bytes = 0u64;

        for slot in slots {
            let index = Arc::clone(slot.run_index().ok_or_else(|| {
                invalid(format!("spill slot {} has no run index", slot.slot_id()))
            })?);
            if slot.is_sequential() || slot.skip_bytes() != 0 {
                return Err(invalid(format!(
                    "spill slot {} cannot be range-read (sequential: {}, skip bytes: {})",
                    slot.slot_id(),
                    slot.is_sequential(),
                    slot.skip_bytes()
                )));
            }
            if K::SERIALIZED_SIZE != Some(usize::from(index.key_width())) {
                return Err(invalid(format!(
                    "spill slot {} run index key width {} does not match the merge key ({:?})",
                    slot.slot_id(),
                    index.key_width(),
                    K::SERIALIZED_SIZE
                )));
            }
            if !index.is_empty() {
                if index.offset(0) != slot.body_start() {
                    return Err(invalid(format!(
                        "spill slot {} run index starts at byte {}, not the body at {}",
                        slot.slot_id(),
                        index.offset(0),
                        slot.body_start()
                    )));
                }
                let end = index.end(index.len() - 1);
                if end > slot.len() {
                    return Err(invalid(format!(
                        "spill slot {} run index ends at byte {end}, past the file's {}",
                        slot.slot_id(),
                        slot.len()
                    )));
                }
            }
            let mut prev: Option<K> = None;
            for i in 0..index.len() {
                let key: K = index.key_at(i)?;
                if prev.is_some_and(|prev| prev > key) {
                    return Err(invalid(format!(
                        "spill slot {} run index frame {i} first key decreases",
                        slot.slot_id()
                    )));
                }
                let start = RangeBound { key, run: sources.len(), entry: i };
                samples.push((start, u64::from(index.raw_len(i))));
                prev = Some(key);
            }
            counted_records += index.total_records();
            run_raw_bytes += index.total_raw_bytes();
            sources.push(RangeSource { slot, index });
        }

        if counted_records != total_records {
            return Err(invalid(format!(
                "run indexes hold {counted_records} records, expected {total_records}"
            )));
        }

        let splitters = splitters_of(samples, sizing);
        Ok(Self { sources, splitters, total_records, run_raw_bytes })
    }

    /// Number of ranges (one more than the bounds between them).
    #[must_use]
    pub fn range_count(&self) -> usize {
        self.splitters.len() + 1
    }

    /// The merge's leaves, in leaf order.
    #[must_use]
    pub(crate) fn sources(&self) -> &[RangeSource] {
        &self.sources
    }

    /// Range `j`'s bounds `[lo, hi)`; `None` is unbounded.
    ///
    /// # Panics
    ///
    /// If `j >= self.range_count()`.
    #[must_use]
    pub fn bounds(&self, j: usize) -> (Option<&RangeBound<K>>, Option<&RangeBound<K>>) {
        assert!(j < self.range_count(), "range {j} out of {}", self.range_count());
        let lo = j.checked_sub(1).map(|p| &self.splitters[p]);
        (lo, self.splitters.get(j))
    }

    /// Index entries of run source `source` that range `j` reads.
    ///
    /// # Panics
    ///
    /// If `source` or `j` is out of bounds. Entry keys were decoded when the
    /// plan was built, so decoding them again cannot fail.
    #[must_use]
    pub fn frame_span(&self, source: usize, j: usize) -> Range<usize> {
        let (lo, hi) = self.bounds(j);
        frame_span_of(&self.sources[source].index, source, lo, hi)
            .expect("entry keys were decoded when planned")
    }

    /// Raw bytes and largest frame of run source `s`'s frame span for range
    /// `j`, plus the span's file bytes; `None` for an empty span.
    fn run_span_bytes(&self, s: usize, j: usize) -> Option<RunSpanBytes> {
        let index = &self.sources[s].index;
        let span = self.frame_span(s, j);
        if span.is_empty() {
            return None;
        }
        let (raw, max_frame) = span
            .clone()
            .map(|f| u64::from(index.raw_len(f)))
            .fold((0u64, 0u64), |(sum, max), n| (sum + n, max.max(n)));
        Some(RunSpanBytes {
            comp: index.end(span.end - 1) - index.offset(span.start),
            raw,
            max_frame,
        })
    }

    /// Bytes range `j` holds while it runs: the file bytes of its frame
    /// spans (read up front), the raw bytes of those spans (bounding the
    /// merged output it keeps from them) and one decoded frame per run with a
    /// non-empty span (each cursor's decode buffer, the span's largest
    /// frame).
    ///
    /// Every run whose key span overlaps the range contributes at least its
    /// boundary frame, which sits in the span, so the charge grows with the
    /// run count even for a small range.
    ///
    /// This is the exact per-range figure of the plan as cut; the merge
    /// decision's fit (`range_fit` in `fgumi-pipeline-io`) models the same
    /// range before any plan exists, as at most `T + kNF` raw bytes held
    /// compressed, merged and batched. The module docs' bound
    /// (`target + k x W`) ties the two: they agree whenever no entry exceeds
    /// `N x F`. The merge charges this figure to its inbox as a gauge; the
    /// fit is what decided the memory.
    #[must_use]
    pub fn estimated_bytes(&self, j: usize) -> u64 {
        (0..self.sources.len())
            .filter_map(|s| self.run_span_bytes(s, j))
            .map(|b| b.comp + b.raw + b.max_frame)
            .sum()
    }

    /// Raw run bytes every range decodes, summed over ranges, divided by the
    /// runs' raw bytes: 1.0 when no frame straddles a range edge, and more as
    /// boundary frames are decoded by both ranges they straddle. 1.0 when the
    /// runs hold no bytes.
    #[must_use]
    #[allow(clippy::cast_precision_loss, reason = "a diagnostic ratio")]
    pub fn planned_decode_amplification(&self) -> f64 {
        if self.run_raw_bytes == 0 {
            return 1.0;
        }
        let decoded: u64 = (0..self.range_count())
            .flat_map(|j| (0..self.sources.len()).map(move |s| (s, j)))
            .filter_map(|(s, j)| self.run_span_bytes(s, j))
            .map(|b| b.raw)
            .sum();
        decoded as f64 / self.run_raw_bytes as f64
    }

    /// Records across all sources.
    #[must_use]
    pub fn total_records(&self) -> u64 {
        self.total_records
    }

    /// Raw bytes of range `j`'s entry spans over every source: what the
    /// module docs bound by `ceil(total / R) + k x W`.
    #[cfg(test)]
    pub(crate) fn span_raw_bytes(&self, j: usize) -> u64 {
        (0..self.sources.len()).filter_map(|s| self.run_span_bytes(s, j)).map(|b| b.raw).sum()
    }

    /// The largest index entry's raw bytes over every source (`W`).
    #[cfg(test)]
    pub(crate) fn max_entry_raw(&self) -> u64 {
        self.sources.iter().map(|s| u64::from(s.index.max_entry_raw())).max().unwrap_or(0)
    }

    /// The range-size bound of the module docs, `ceil(total / R) + k x W`,
    /// for the sizing the plan was built with: `R` from the summed entry
    /// weights, exactly as [`splitters_of`] computes it.
    #[cfg(test)]
    pub(crate) fn range_bound(&self, sizing: RangeSizing) -> u64 {
        let total: u64 = self
            .sources
            .iter()
            .flat_map(|s| (0..s.index.len()).map(|i| u64::from(s.index.raw_len(i))))
            .sum();
        assert_eq!(total, self.run_raw_bytes, "entry weights sum to the runs' raw bytes");
        let ranges = sizing.range_count(u128::from(total)) as u64;
        total.div_ceil(ranges) + self.sources.len() as u64 * self.max_entry_raw()
    }
}

/// File and raw bytes of one run's frame span for one range.
struct RunSpanBytes {
    /// File bytes (`[u32 len][zstd frame]` units) of the span.
    comp: u64,
    /// Decompressed bytes of the span.
    raw: u64,
    /// Decompressed bytes of the span's largest frame.
    max_frame: u64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SpillCodec;
    use crate::keys::RawCoordinateKey;
    use proptest::prelude::*;
    use rstest::rstest;

    const MIB: u64 = 1 << 20;
    const GIB: u64 = 1 << 30;

    fn key(n: u64) -> RawCoordinateKey {
        RawCoordinateKey { sort_key: n }
    }

    /// The bound at entry `entry` of source `run`, whose first key is `k`.
    fn bound(k: u64, run: usize, entry: usize) -> RangeBound<RawCoordinateKey> {
        RangeBound { key: key(k), run, entry }
    }

    /// A slot over `frames` of `(first_key, raw_len, records)`, each frame
    /// spanning `comp_len` file bytes, with the index attached.
    fn indexed_slot(
        file_id: u32,
        body_start: u64,
        comp_len: u32,
        frames: &[(u64, u32, u32)],
    ) -> Arc<SortMergeSlot> {
        let mut index = RunIndex::new(8);
        for &(k, raw_len, records) in frames {
            index.push_frame(comp_len, raw_len, records, &k.to_le_bytes());
        }
        index.rebase(body_start);
        slot_with(file_id, body_start, comp_len, Some(index))
    }

    fn slot_with(
        file_id: u32,
        body_start: u64,
        comp_len: u32,
        index: Option<RunIndex>,
    ) -> Arc<SortMergeSlot> {
        let frames = index.as_ref().map_or(0, RunIndex::covered_frames);
        let len = body_start + frames * u64::from(comp_len);
        let file = tempfile::tempfile().expect("tempfile");
        let slot = SortMergeSlot::new(file_id, Arc::new(file), body_start, len, SpillCodec::Zstd);
        if let Some(index) = index {
            slot.attach_run_index(Arc::new(index));
        }
        Arc::new(slot)
    }

    fn sizing(target_range_bytes: u64, max_ranges: usize) -> RangeSizing {
        RangeSizing { target_range_bytes, max_ranges }
    }

    #[rstest]
    #[case::one_range(10 * MIB, 64 * MIB, MAX_RANGES, false, 1)]
    #[case::ten(100 * MIB, 10 * MIB, MAX_RANGES, false, 10)]
    #[case::capped(10 * GIB, MIB, 16, false, 16)]
    // Equal keys no longer pin one range: bounds fall between their entries.
    #[case::all_equal_keys(100 * MIB, MIB, MAX_RANGES, true, 100)]
    fn range_count_follows_bytes(
        #[case] total: u64,
        #[case] target: u64,
        #[case] max_ranges: usize,
        #[case] equal_keys: bool,
        #[case] expected: usize,
    ) {
        const FRAMES: u64 = 1024;
        let raw_len = u32::try_from(total / FRAMES).unwrap();
        let frames: Vec<_> =
            (0..FRAMES).map(|f| (if equal_keys { 7 } else { f }, raw_len, 1)).collect();
        let plan: RangePlan<RawCoordinateKey> = RangePlan::new(
            vec![indexed_slot(0, 0, 100, &frames)],
            FRAMES,
            sizing(target, max_ranges),
        )
        .unwrap();
        assert_eq!(plan.range_count(), expected);
        assert_eq!(plan.run_raw_bytes, total);
    }

    /// Bounds are `(key, source, entry)`; the index (frames keyed
    /// `[10, 20, 20, 30]`) is source `run`'s.
    #[rstest]
    #[case::first(0, None, Some((20, 0, 1)), 0..1)]
    #[case::own_lo_starts_at_its_entry(0, Some((20, 0, 1)), Some((30, 0, 3)), 1..3)]
    #[case::last(0, Some((30, 0, 3)), None, 3..4)]
    // `lo` on an earlier source: this source's 20s sort above it.
    #[case::lo_on_an_earlier_source(1, Some((20, 0, 1)), None, 0..4)]
    // `lo` on a later source: this source's 20s sort below it; the last entry
    // starting at or below 20 may hold larger keys.
    #[case::lo_on_a_later_source(1, Some((20, 2, 1)), None, 2..4)]
    #[case::hi_on_a_later_source(1, None, Some((20, 2, 1)), 0..3)]
    #[case::hi_on_an_earlier_source(1, None, Some((20, 0, 1)), 0..1)]
    #[case::below_all(1, None, Some((5, 0, 0)), 0..0)]
    #[case::above_all(1, Some((40, 0, 0)), None, 3..4)]
    fn frame_span_rule(
        #[case] run: usize,
        #[case] lo: Option<(u64, usize, usize)>,
        #[case] hi: Option<(u64, usize, usize)>,
        #[case] expected: Range<usize>,
    ) {
        let mut index = RunIndex::new(8);
        for k in [10u64, 20, 20, 30] {
            index.push_frame(10, 100, 1, &k.to_le_bytes());
        }
        let lo = lo.map(|(k, s, e)| bound(k, s, e));
        let hi = hi.map(|(k, s, e)| bound(k, s, e));
        let got = frame_span_of(&index, run, lo.as_ref(), hi.as_ref()).unwrap();
        assert_eq!(got, expected);
    }

    /// A run of five frames (100 file bytes, 1000 raw bytes, 2 records each)
    /// keyed `[10, 20, 30, 40, 50]`, sparsified by 2 into entries keyed
    /// `[10, 30, 50]` (the last covering one frame), cut at entries 1 and 2.
    #[rstest]
    #[case::below_thirty(0, 0..1, 200)]
    #[case::thirty_to_fifty(1, 1..2, 200)]
    #[case::fifty_up(2, 2..3, 100)]
    fn frame_span_on_a_sparsified_index(
        #[case] j: usize,
        #[case] entries: Range<usize>,
        #[case] comp_bytes: u64,
    ) {
        let mut index = RunIndex::new(8);
        for k in [10u64, 20, 30, 40, 50] {
            index.push_frame(100, 1000, 2, &k.to_le_bytes());
        }
        index.rebase(16);
        index.sparsify(2).unwrap();
        assert_eq!((index.len(), index.stride(), index.covered_frames()), (3, 2, 5));
        let mut plan =
            RangePlan::new(vec![slot_with(0, 16, 100, Some(index))], 10, sizing(MIB, MAX_RANGES))
                .unwrap();
        plan.splitters = vec![bound(30, 0, 1), bound(50, 0, 2)];
        assert_eq!(plan.frame_span(0, j), entries);
        assert_eq!(plan.run_span_bytes(0, j).unwrap().comp, comp_bytes);
    }

    /// `runs` runs, each with frames keyed `[10, 20, 20, 30]` (1000 raw
    /// bytes, 100 file bytes, 2 records each), cut at `splitters`.
    fn run_plan(runs: u32, splitters: &[(u64, usize, usize)]) -> RangePlan<RawCoordinateKey> {
        let frames = [(10, 1000, 2), (20, 1000, 2), (20, 1000, 2), (30, 1000, 2)];
        let slots = (0..runs).map(|id| indexed_slot(id, 16, 100, &frames)).collect();
        let mut plan = RangePlan::new(slots, 8 * u64::from(runs), sizing(1, MAX_RANGES)).unwrap();
        plan.splitters = splitters.iter().map(|&(k, s, e)| bound(k, s, e)).collect();
        plan
    }

    #[rstest]
    // Run span file + raw bytes, then one decoded frame (1000) for the run's
    // cursor.
    #[case::below_twenty(0, (None, Some(20)), 0..1, 1100 + 1000)]
    #[case::twenty_to_thirty(1, (Some(20), Some(30)), 1..3, 2 * 1100 + 1000)]
    #[case::thirty_up(2, (Some(30), None), 3..4, 1100 + 1000)]
    fn spans_and_estimated_bytes(
        #[case] j: usize,
        #[case] bounds: (Option<u64>, Option<u64>),
        #[case] frames: Range<usize>,
        #[case] bytes: u64,
    ) {
        let plan = run_plan(1, &[(20, 0, 1), (30, 0, 3)]);
        let (lo, hi) = plan.bounds(j);
        assert_eq!((lo.map(|b| b.key), hi.map(|b| b.key)), (bounds.0.map(key), bounds.1.map(key)));
        assert_eq!(plan.frame_span(0, j), frames);
        assert_eq!(plan.estimated_bytes(j), bytes);
    }

    /// The planned decode amplification of [`run_plan`]'s runs (1000 raw
    /// bytes a frame) under `splitters`.
    #[rstest]
    #[case::one_range(1, vec![], 1.0)]
    // Spans 0..1, 1..3 and 3..4: a bound on a run's own entry decodes nothing
    // twice.
    #[case::three_ranges(1, vec![(20, 0, 1), (30, 0, 3)], 1.0)]
    #[case::between_equal_keys(1, vec![(20, 0, 2)], 1.0)]
    // A bound inside the 20s at run 1's entry 1: run 0 (an earlier source)
    // reads 0..3 then 2..4, run 1 reads 0..1 then 1..4. 9000 of 8000 bytes.
    #[case::inside_a_tie_across_runs(2, vec![(20, 1, 1)], 1.125)]
    fn planned_decode_amplification_counts_boundary_frames(
        #[case] runs: u32,
        #[case] splitters: Vec<(u64, usize, usize)>,
        #[case] expected: f64,
    ) {
        let plan = run_plan(runs, &splitters);
        assert!((plan.planned_decode_amplification() - expected).abs() < 1e-9);
    }

    #[test]
    fn sources_follow_the_merge_leaf_order() {
        let plan = RangePlan::<RawCoordinateKey>::new(
            vec![indexed_slot(9, 0, 10, &[(1, 10, 1)]), indexed_slot(2, 0, 10, &[(1, 10, 1)])],
            2,
            sizing(MIB, MAX_RANGES),
        )
        .unwrap();
        let order: Vec<String> = plan
            .sources()
            .iter()
            .map(|RangeSource { slot, index }| {
                format!("run{}@[{}, {}]", slot.file_id, index.offset(0), index.end(0))
            })
            .collect();
        assert_eq!(order, ["run2@[0, 10]", "run9@[0, 10]"]);
        assert_eq!(plan.total_records(), 2);
    }

    /// Per-run frames `(first_key, raw_len)`, target range bytes, the range
    /// cap, and the index stride.
    type PlanInputs = (Vec<Vec<(u64, u32)>>, u64, usize, u32);

    /// Runs of frame keys from a small alphabet (heavy duplicates); [`build`]
    /// sorts them.
    fn plan_inputs() -> impl Strategy<Value = PlanInputs> {
        (
            prop::collection::vec(prop::collection::vec((0u64..4, 1u32..50_000), 0..40), 0..6),
            1u64..400_000,
            1usize..64,
            1u32..=4,
        )
    }

    fn build(
        runs: &[Vec<(u64, u32)>],
        target: u64,
        max_ranges: usize,
        stride: u32,
    ) -> RangePlan<RawCoordinateKey> {
        let mut slots = Vec::new();
        let mut records = 0u64;
        for (id, run) in runs.iter().enumerate() {
            let mut run = run.clone();
            run.sort_unstable();
            let mut index = RunIndex::new(8);
            for &(k, raw) in &run {
                index.push_frame(64, raw, 3, &k.to_le_bytes());
            }
            index.sparsify(stride).unwrap();
            records += 3 * run.len() as u64;
            slots.push(slot_with(u32::try_from(id).unwrap(), 0, 64, Some(index)));
        }
        RangePlan::new(slots, records, sizing(target, max_ranges)).unwrap()
    }

    proptest! {
        #[test]
        fn splitters_are_strictly_increasing_and_cover(
            (runs, target, max_ranges, stride) in plan_inputs()
        ) {
            let plan = build(&runs, target, max_ranges, stride);
            prop_assert!(plan.splitters.windows(2).all(|w| w[0] < w[1]));
            // Every bound starts an entry of its source whose first key is
            // the bound's key.
            for b in &plan.splitters {
                let index = &plan.sources[b.run].index;
                prop_assert!(b.entry < index.len());
                prop_assert_eq!(index.key_at::<RawCoordinateKey>(b.entry).unwrap(), b.key);
            }
            prop_assert_eq!(plan.range_count(), plan.splitters.len() + 1);
            prop_assert!(plan.range_count() <= max_ranges);
            let last = plan.range_count() - 1;
            prop_assert!(plan.bounds(0).0.is_none());
            prop_assert!(plan.bounds(last).1.is_none());
            // Consecutive ranges share their boundary.
            for j in 1..plan.range_count() {
                prop_assert_eq!(plan.bounds(j - 1).1, plan.bounds(j).0);
            }
        }

        /// However many entries share a key, every range's span stays within
        /// `ceil(total / R) + k x W`, and the ranges together read every entry.
        #[test]
        fn ranges_are_bounded_under_heavy_duplication(
            (runs, target, max_ranges, stride) in plan_inputs()
        ) {
            let s = sizing(target, max_ranges);
            let plan = build(&runs, target, max_ranges, stride);
            let limit = plan.range_bound(s);
            for j in 0..plan.range_count() {
                let raw = plan.span_raw_bytes(j);
                prop_assert!(raw <= limit, "range {} spans {} raw bytes > {}", j, raw, limit);
            }
            for (src, source) in plan.sources.iter().enumerate() {
                let mut read = vec![false; source.index.len()];
                for j in 0..plan.range_count() {
                    for e in plan.frame_span(src, j) {
                        read[e] = true;
                    }
                }
                prop_assert!(read.iter().all(|&r| r), "source {} has an unread entry", src);
            }
        }
    }

    #[test]
    fn plan_is_deterministic() {
        let runs = vec![
            vec![(1, 40_000), (3, 9_000), (3, 70_000), (8, 1_000)],
            vec![(0, 20_000), (3, 20_000), (5, 60_000), (9, 5_000)],
        ];
        let a = build(&runs, 30_000, MAX_RANGES, 1);
        let b = build(&runs, 30_000, MAX_RANGES, 1);
        assert!(a.range_count() > 2, "the inputs should cut several ranges");
        assert_eq!(a.splitters, b.splitters);
    }

    #[derive(Debug, Clone, Copy)]
    enum Reject {
        SlotWithoutIndex,
        IndexWidthMismatch,
        IndexRecordsDisagree,
        IndexDoesNotStartAtBody,
        IndexPastFile,
        SequentialSlot,
        SkipBytes,
    }

    /// Every rejection is `InvalidData`, and its message names the cause.
    /// Columns: the defect, a needle of its message.
    #[rstest]
    #[case::slot_without_index(Reject::SlotWithoutIndex, "spill slot 1 has no run index")]
    #[case::index_width_mismatch(Reject::IndexWidthMismatch, "run index key width 24")]
    #[case::index_records_disagree(Reject::IndexRecordsDisagree, "hold 8 records, expected 9")]
    #[case::index_does_not_start_at_body(
        Reject::IndexDoesNotStartAtBody,
        "starts at byte 0, not the body at 4"
    )]
    #[case::index_past_file(Reject::IndexPastFile, "ends at byte 15, past the file's 14")]
    #[case::sequential_slot(Reject::SequentialSlot, "sequential: true")]
    #[case::skip_bytes(Reject::SkipBytes, "skip bytes: 4")]
    fn plan_rejects(#[case] reject: Reject, #[case] needle: &str) {
        let good = indexed_slot(0, 0, 10, &[(1, 100, 4)]);
        let (bad, total) = match reject {
            Reject::SlotWithoutIndex => (slot_with(1, 0, 10, None), 4),
            Reject::IndexWidthMismatch => {
                let mut index = RunIndex::new(24);
                index.push_frame(10, 100, 4, &[0u8; 24]);
                (slot_with(1, 0, 10, Some(index)), 4 + 4)
            }
            // Indexes hold 8 records, not 9.
            Reject::IndexRecordsDisagree => (indexed_slot(1, 0, 10, &[(2, 100, 4)]), 9),
            // The slot's body starts at 4 but its index starts at 0.
            Reject::IndexDoesNotStartAtBody => {
                let mut index = RunIndex::new(8);
                index.push_frame(10, 100, 4, &2u64.to_le_bytes());
                (slot_with(1, 4, 10, Some(index)), 4 + 4)
            }
            // The index ends at byte 14; the slot's file holds 4 + 10.
            Reject::IndexPastFile => {
                let mut index = RunIndex::new(8);
                index.push_frame(11, 100, 4, &2u64.to_le_bytes());
                index.rebase(4);
                (slot_with(1, 4, 10, Some(index)), 4 + 4)
            }
            Reject::SequentialSlot | Reject::SkipBytes => {
                let file = Arc::new(tempfile::tempfile().expect("tempfile"));
                let slot = if matches!(reject, Reject::SequentialSlot) {
                    SortMergeSlot::new_sequential(1, file, SpillCodec::Zstd)
                } else {
                    SortMergeSlot::new(1, file, 0, 10, SpillCodec::Zstd).with_skip_bytes(4)
                };
                let mut index = RunIndex::new(8);
                index.push_frame(10, 100, 4, &2u64.to_le_bytes());
                slot.attach_run_index(Arc::new(index));
                (Arc::new(slot), 4 + 4)
            }
        };
        let err =
            RangePlan::<RawCoordinateKey>::new(vec![good, bad], total, sizing(MIB, MAX_RANGES))
                .err()
                .expect("plan should be rejected");
        assert_eq!(err.kind(), io::ErrorKind::InvalidData, "{err}");
        assert!(err.to_string().contains(needle), "{err}");
    }
}
