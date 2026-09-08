//! Inline consensus-metrics accumulator: the per-worker state fed by the
//! fused (T1, Task 9) and standalone (T2, Task 8) adapters to compute the
//! same CS/SS/DS family-size, UMI, and downsampling-yield metrics the
//! separate-pass simplex-metrics/duplex-metrics commands compute, inline,
//! during the existing streaming pass. See
//! docs/superpowers/specs/2026-09-08-inline-consensus-metrics-design.md §3-§6.
//!
//! **Ported ahead of its caller.** Every item here is `pub(crate)`. Until
//! Task 8's `Serial` collector step and Task 9's fused T1 accumulator wire
//! this module in, nothing outside its own tests constructs a
//! `ConsensusMetricsAccumulator`, so the whole module is dead code to the
//! compiler. The `allow` below is scoped to this module and should be
//! deleted once a caller lands (Task 8 or 9).
//!
//! **Fix round 1 (post-Task-7 review):** the original per-slot design called
//! `record_simplex_coordinate_group`/`record_duplex_coordinate_group` once
//! *per fraction slot*, each time with a 1-element `fractions` window. Both
//! reducers gate their once-per-group work (UMI decoding for Simplex;
//! `duplex_umi_counts`'s UMI-string building for Duplex) on
//! `idx == fractions.len() - 1` — with a 1-element window that's always
//! `0 == 0`, so the gate fired on *every* fraction instead of only the 100%
//! one (spec §5.6). Calling per-slot also recomputed
//! `compute_template_metadata(group)` 20x per coordinate group instead of
//! once. The fix restructures this type from a per-slot `[InlineCollector; 20]`
//! array into two homogeneous-array variants (`Simplex`/`Duplex`) so each
//! reducer is called exactly ONCE per group with the full 20-length
//! `DOWNSAMPLING_FRACTIONS`/`collectors`/`fraction_template_counts` arrays —
//! matching how the separate-pass `simplex_metrics`/`duplex_metrics` commands
//! call them (`simplex_metrics.rs:131-147`, `duplex_metrics.rs:176`).

#![allow(dead_code)]

use std::collections::HashMap;

use crate::commands::shared_metrics::{
    DOWNSAMPLING_FRACTIONS, Interval, ReadInfoKey, TemplateInfo, build_template_info,
    overlaps_intervals, record_duplex_coordinate_group, record_simplex_coordinate_group,
};
use crate::mi_group::MiGroup;
use crate::simple_umi_consensus::SimpleUmiConsensusCaller;
use anyhow::Result;
use fgumi_bam_io::LibraryIndex;
use fgumi_metrics::duplex::DuplexMetricsCollector;
use fgumi_metrics::simplex::SimplexMetricsCollector;
use fgumi_pipeline_core::{HeapSize, MetricsReducer, Ordered};
use fgumi_raw_bam::{RawRecord, flags as raw_flags};

/// Accumulator state for one consensus stage's inline metrics. Shared by T1
/// (wrapped in `PerThreadAccumulator`, N sharded instances folded via
/// `merge()`) and T2 (a single un-sharded instance owned directly by the
/// `Serial` collector step, Task 8 — `merge()` is never called on that path).
///
/// One variant per mode, each holding a homogeneous 20-element array (one
/// slot per `DOWNSAMPLING_FRACTIONS` entry) so the shared reducer functions
/// can be called once per coordinate group with the FULL array — required
/// for their internal `idx == fractions.len() - 1` ("only at 100%") gates to
/// see the true last index, and to amortize `compute_template_metadata`
/// across all 20 fractions instead of recomputing it per fraction.
///
/// `fraction_template_counts` mirrors the separate-pass commands'
/// `process_templates_from_bam` accumulator (`shared_metrics.rs:893`): a
/// running per-fraction template count across every coordinate group this
/// accumulator ever sees, not per-call scratch. It is persisted here (not
/// reallocated per `record_coordinate_group` call) so a future finalize step
/// can build yield metrics from it the same way the separate-pass commands
/// do (`collector.into_yield_metric(fraction, read_pairs, min_reads)` with
/// `read_pairs` sourced from this array).
pub(crate) enum ConsensusMetricsAccumulator {
    Simplex {
        // Boxed: `SimplexMetricsCollector` and `DuplexMetricsCollector` are
        // different sizes, so an unboxed `[Collector; 20]` field on each
        // variant makes this enum's stack footprint the size of its LARGEST
        // variant (clippy::large_enum_variant) — every `Simplex` value would
        // carry the unused space for a `Duplex`-sized payload. Boxing shrinks
        // each variant to a pointer.
        collectors: Box<[SimplexMetricsCollector; 20]>,
        umi_caller: SimpleUmiConsensusCaller,
        fraction_template_counts: [usize; 20],
    },
    Duplex {
        collectors: Box<[DuplexMetricsCollector; 20]>,
        umi_caller: SimpleUmiConsensusCaller,
        fraction_template_counts: [usize; 20],
        /// Mirrors `duplex-metrics --duplex-umi-counts` (`self.duplex_umi_counts`
        /// on `DuplexMetrics`, `duplex_metrics.rs:104-105`), threaded into
        /// `record_duplex_coordinate_group`'s `duplex_umi_counts` parameter
        /// (Task 5).
        duplex_umi_counts: bool,
    },
}

impl ConsensusMetricsAccumulator {
    pub(crate) fn new_simplex() -> Self {
        Self::Simplex {
            collectors: Box::new(std::array::from_fn(|_| SimplexMetricsCollector::new())),
            umi_caller: SimpleUmiConsensusCaller::default(),
            fraction_template_counts: [0; 20],
        }
    }

    pub(crate) fn new_duplex(collect_duplex_umi_counts: bool) -> Self {
        Self::Duplex {
            collectors: Box::new(std::array::from_fn(|_| {
                DuplexMetricsCollector::new(collect_duplex_umi_counts)
            })),
            umi_caller: SimpleUmiConsensusCaller::default(),
            fraction_template_counts: [0; 20],
            duplex_umi_counts: collect_duplex_umi_counts,
        }
    }

    /// Records one already-assembled coordinate/strand group, applying the
    /// per-template interval filter itself (matching
    /// `process_templates_from_bam`'s ordering: interval filter strictly
    /// before any fraction bucketing) — callers pass the RAW group here.
    ///
    /// Calls the shared reducer exactly ONCE with the full 20-length
    /// `DOWNSAMPLING_FRACTIONS` array (see the type doc for why this matters).
    pub(crate) fn record_coordinate_group(
        &mut self,
        group: &[TemplateInfo],
        intervals: &[Interval],
    ) -> Result<()> {
        let filtered: Vec<TemplateInfo> =
            group.iter().filter(|t| overlaps_intervals(t, intervals)).cloned().collect();
        if filtered.is_empty() {
            return Ok(());
        }
        match self {
            Self::Simplex { collectors, umi_caller, fraction_template_counts } => {
                record_simplex_coordinate_group(
                    &filtered,
                    &DOWNSAMPLING_FRACTIONS,
                    &mut collectors[..],
                    umi_caller,
                    fraction_template_counts,
                )
            }
            Self::Duplex {
                collectors,
                umi_caller,
                fraction_template_counts,
                duplex_umi_counts,
            } => record_duplex_coordinate_group(
                &filtered,
                &DOWNSAMPLING_FRACTIONS,
                &mut collectors[..],
                umi_caller,
                fraction_template_counts,
                *duplex_umi_counts,
            ),
        }
    }

    /// Folds `other` into `self` by zipping and merging the 20 per-fraction
    /// collectors (and per-fraction template counts) pairwise. **T1 only** —
    /// T2's `Serial` collector step (Task 8) owns exactly one instance and
    /// never calls this.
    ///
    /// The `duplex_umi_counts` flag is left untouched — like
    /// `DuplexMetricsCollector::merge`'s `collect_duplex_umi_counts` field
    /// (Task 3), it is a per-stage config flag, not accumulated state; every
    /// accumulator for one stage is always constructed with the same value.
    pub(crate) fn merge(&mut self, other: Self) -> Result<()> {
        match (self, other) {
            (
                Self::Simplex { collectors: mine, fraction_template_counts: mine_counts, .. },
                Self::Simplex {
                    collectors: theirs, fraction_template_counts: theirs_counts, ..
                },
            ) => {
                for (m, t) in mine.iter_mut().zip(*theirs) {
                    m.merge(t);
                }
                for (m, t) in mine_counts.iter_mut().zip(theirs_counts) {
                    *m += t;
                }
            }
            (
                Self::Duplex { collectors: mine, fraction_template_counts: mine_counts, .. },
                Self::Duplex {
                    collectors: theirs, fraction_template_counts: theirs_counts, ..
                },
            ) => {
                for (m, t) in mine.iter_mut().zip(*theirs) {
                    m.merge(t);
                }
                for (m, t) in mine_counts.iter_mut().zip(theirs_counts) {
                    *m += t;
                }
            }
            _ => anyhow::bail!(
                "internal error: attempted to merge a Simplex accumulator with a Duplex \
                 one — every ConsensusMetricsAccumulator for one stage must be constructed \
                 via the same new_simplex()/new_duplex() call, never mixed"
            ),
        }
        Ok(())
    }
}

/// A `CoordinateGroupFragment` carries one **batch**'s worth of already-paired
/// `(TemplateInfo, ReadInfoKey)` entries — matching the verified precedent the
/// existing rejects branch already establishes (X5-001, spec §3: exactly one
/// item pushed onto the ordered third branch per input batch, always, even
/// when empty) — so `Task 11`'s wiring pushes one of these per `MiGroup`
/// batch, tagged with that batch's `batch_serial`.
///
/// `Ordered::ordinal` returns `batch_serial` so the framework's
/// `ReorderStage` (Task 11) can restore true cross-batch emission order in
/// front of [`CoordinateGroupCollector`] — the whole reason this collector
/// never needs to reconstruct a group split across worker threads: by the
/// time it sees fragments, they arrive strictly in order, adjacent.
pub(crate) struct CoordinateGroupFragment {
    pub(crate) batch_serial: u64,
    pub(crate) entries: Vec<(TemplateInfo, ReadInfoKey)>,
}

impl Ordered for CoordinateGroupFragment {
    fn ordinal(&self) -> u64 {
        self.batch_serial
    }
}

impl HeapSize for CoordinateGroupFragment {
    /// Approximate heap footprint: the entries `Vec`'s allocated capacity.
    /// `TemplateInfo`/`ReadInfoKey` themselves hold their own heap
    /// allocations (`String`s, an optional `Box<[u8]>`), but this collector
    /// is never routed through a `ByteBoundedQueue` in this task (Task 8 unit
    /// tests it directly, with no queue at all) — Task 11's wiring is what
    /// actually exercises this bound, and a capacity-only estimate matches
    /// the level of precision `BatchedMiGroups::heap_size` already accepts
    /// for its own `Vec<MiGroup>` (`crate::pipeline::steps::group::mi`).
    fn heap_size(&self) -> usize {
        self.entries.capacity() * std::mem::size_of::<(TemplateInfo, ReadInfoKey)>()
    }
}

/// Re-pairs one `MiGroup`'s flat record list into R1/R2 pairs by read name,
/// applying the same paired/mapped/primary filter
/// `process_templates_from_bam` uses. Records within one `MiGroup` all
/// belong to the same UMI family (already clustered at essentially one
/// physical position by upstream grouping), so pairing via an unordered
/// `HashMap` does not risk interleaving templates from genuinely different
/// `ReadInfoKey`s within a single `MiGroup`.
fn pair_records_by_read_name(records: &[RawRecord]) -> Vec<(RawRecord, RawRecord)> {
    let mut by_name: HashMap<Vec<u8>, (Option<RawRecord>, Option<RawRecord>)> = HashMap::new();
    for record in records {
        let flags = record.flags();
        let qualifies = (flags & raw_flags::PAIRED) != 0
            && (flags & raw_flags::UNMAPPED) == 0
            && (flags & raw_flags::MATE_UNMAPPED) == 0
            && (flags & raw_flags::SECONDARY) == 0
            && (flags & raw_flags::SUPPLEMENTARY) == 0;
        if !qualifies {
            continue;
        }
        let name = fgumi_raw_bam::read_name(record.as_ref()).to_vec();
        let entry = by_name.entry(name).or_default();
        if (flags & raw_flags::FIRST_SEGMENT) != 0 {
            entry.0 = Some(record.clone());
        } else if (flags & raw_flags::LAST_SEGMENT) != 0 {
            entry.1 = Some(record.clone());
        }
    }
    by_name
        .into_values()
        .filter_map(|(r1, r2)| match (r1, r2) {
            (Some(r1), Some(r2)) => Some((r1, r2)),
            _ => None,
        })
        .collect()
}

/// Appends one `MiGroup`'s entries onto a caller-owned, per-**batch**
/// `entries` accumulator. **One `CoordinateGroupFragment` is emitted per
/// BATCH, not per `MiGroup`** — Task 11's wiring calls this once per
/// `MiGroup` inside a batch's loop, into ONE shared `Vec`, then constructs
/// exactly one `CoordinateGroupFragment { batch_serial, entries }` at the end
/// of that batch's processing. (Within one batch, `MiGroup`s are already
/// processed strictly in order by one worker thread — only *cross-batch*
/// order needs `ReorderStage` to restore, which is exactly what tagging the
/// whole batch's fragment with that batch's `batch_serial` achieves.) Every
/// paired, successfully-converted template in the group becomes one
/// `(TemplateInfo, ReadInfoKey)` entry; a template that fails to produce one
/// (missing R1/R2, unmapped, no CIGAR) is silently omitted, matching
/// `process_templates_from_bam`'s own behavior.
///
/// # Errors
///
/// Returns an error if a qualifying R1/R2 pair is missing a required `MI`/`RX`
/// tag (propagated from [`build_template_info`]).
pub(crate) fn push_mi_group_entries(
    mi_group: &MiGroup,
    header: &noodles::sam::Header,
    library_index: &LibraryIndex,
    entries: &mut Vec<(TemplateInfo, ReadInfoKey)>,
) -> Result<()> {
    for (r1, r2) in pair_records_by_read_name(&mi_group.records) {
        if let Some((info, key)) = build_template_info(&r1, &r2, header, library_index)? {
            entries.push((info, key));
        }
    }
    Ok(())
}

/// `T2` standalone-consensus reducer: owns ONE un-sharded
/// `ConsensusMetricsAccumulator`, buffers consecutive same-`ReadInfoKey`
/// entries across fragments, and flushes a coordinate group into the
/// accumulator when the key changes. Fed by an ordered stream of
/// `CoordinateGroupFragment`s (Task 11's `ReorderStage` restores true
/// emission order upstream), so — unlike the retired
/// `MiGroupCoordinateBuffer`/`open_groups` carry-state design — there is
/// never more than one instance and never any question of which worker saw
/// which fragment first.
pub(crate) struct CoordinateGroupCollector {
    current_key: Option<ReadInfoKey>,
    current_group: Vec<TemplateInfo>,
    accumulator: ConsensusMetricsAccumulator,
    intervals: Vec<Interval>,
    on_finish: Box<dyn FnOnce(ConsensusMetricsAccumulator) -> Result<()> + Send>,
}

impl CoordinateGroupCollector {
    pub(crate) fn new(
        accumulator: ConsensusMetricsAccumulator,
        intervals: Vec<Interval>,
        on_finish: Box<dyn FnOnce(ConsensusMetricsAccumulator) -> Result<()> + Send>,
    ) -> Self {
        Self { current_key: None, current_group: Vec::new(), accumulator, intervals, on_finish }
    }
}

impl MetricsReducer for CoordinateGroupCollector {
    type Item = CoordinateGroupFragment;

    /// A near-verbatim port of `process_templates_from_bam`'s inner loop
    /// (`shared_metrics.rs`) and of `GroupByMi::process_record`'s own
    /// buffer-then-flush shape (`mi.rs`) — the only difference is that this
    /// reads from an already-ordered stream of fragments (each fragment
    /// holding one batch's worth of entries) rather than raw BAM records one
    /// at a time.
    fn record(&mut self, fragment: CoordinateGroupFragment) -> Result<()> {
        for (info, key) in fragment.entries {
            match &self.current_key {
                Some(k) if *k == key => self.current_group.push(info),
                Some(_) => {
                    let finished = std::mem::take(&mut self.current_group);
                    self.accumulator.record_coordinate_group(&finished, &self.intervals)?;
                    self.current_group = vec![info];
                    self.current_key = Some(key);
                }
                None => {
                    self.current_group.push(info);
                    self.current_key = Some(key);
                }
            }
        }
        Ok(())
    }

    /// End-of-stream flush, mirroring `process_templates_from_bam`'s own
    /// end-of-loop flush, then hands the finished accumulator to
    /// `on_finish` — real file-writing is supplied by the caller (Task 11),
    /// not this type; Task 8 tests `on_finish` with a plain observer closure.
    fn finish(mut self) -> Result<()> {
        if !self.current_group.is_empty() {
            self.accumulator.record_coordinate_group(&self.current_group, &self.intervals)?;
        }
        (self.on_finish)(self.accumulator)
    }
}

/// Shared `#[cfg(test)]` fixtures, used by both `mod tests` (Task 7's
/// `ConsensusMetricsAccumulator` tests) and `mod coordinate_group_collector_tests`
/// below.
#[cfg(test)]
pub(crate) fn template(
    mi: &str,
    ref_name: &str,
    position: i32,
    hash_fraction: f64,
) -> TemplateInfo {
    TemplateInfo {
        mi: mi.to_string(),
        rx: "AAA".to_string(),
        ref_name: Some(ref_name.to_string()),
        position: Some(position),
        end_position: Some(position + 50),
        r1_positive: true,
        hash_fraction,
    }
}

/// Builds a minimal `ReadInfoKey` distinguished only by `ref_index`/`start1`
/// — the two fields the ordered-stream tests below vary — with every other
/// field pinned to a fixed, arbitrary value so two keys built from the same
/// `(ref_index, start)` always compare equal.
#[cfg(test)]
pub(crate) fn key(ref_index: usize, start: i32) -> ReadInfoKey {
    ReadInfoKey {
        ref_index1: ref_index,
        start1: start,
        strand1: false,
        ref_index2: ref_index,
        start2: start + 50,
        strand2: true,
        library: 0,
        cell_barcode: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn record_coordinate_group_drops_templates_outside_intervals() {
        let mut acc = ConsensusMetricsAccumulator::new_simplex();
        let group = vec![
            template("1", "chr1", 100, 1.0),
            template("2", "chr2", 100, 1.0), // different ref, outside the interval
        ];
        let intervals = vec![Interval { ref_name: "chr1".to_string(), start: 0, end: 1000 }];

        acc.record_coordinate_group(&group, &intervals).expect("records");

        let ConsensusMetricsAccumulator::Simplex { collectors, .. } = &acc else {
            panic!("expected Simplex variant");
        };
        let metrics = collectors[19].family_size_metrics();
        let cs_count_at_size_1: usize =
            metrics.iter().filter(|m| m.family_size == 1).map(|m| m.cs_count).sum();
        assert_eq!(cs_count_at_size_1, 1);
    }

    #[test]
    fn merge_folds_two_accumulators_fraction_collectors() {
        let mut a = ConsensusMetricsAccumulator::new_simplex();
        a.record_coordinate_group(&[template("1", "chr1", 100, 1.0)], &[]).expect("records");

        let mut b = ConsensusMetricsAccumulator::new_simplex();
        b.record_coordinate_group(&[template("2", "chr1", 200, 1.0)], &[]).expect("records");

        a.merge(b).expect("merge succeeds");

        let ConsensusMetricsAccumulator::Simplex { collectors, .. } = &a else {
            panic!("expected Simplex variant");
        };
        let metrics = collectors[19].family_size_metrics();
        let cs_count_at_size_1: usize =
            metrics.iter().filter(|m| m.family_size == 1).map(|m| m.cs_count).sum();
        assert_eq!(
            cs_count_at_size_1, 2,
            "both accumulators' size-1 families must be present after merge"
        );
    }

    /// Guards against the fix-round-1 regression: the original per-slot
    /// implementation called `record_simplex_coordinate_group` once per
    /// fraction slot with a 1-element `fractions` window, so the reducer's
    /// internal `idx == fractions.len() - 1` ("only build UMI metrics at the
    /// 100% fraction") gate saw `0 == 0` on every call and fired at every
    /// fraction, not just the 100% one. This test uses `hash_fraction = 0.05`
    /// so the group is included in EVERY downsampling fraction's subset
    /// (0.05 <= every fraction from 0.05 to 1.00) — under the old per-slot
    /// code, `collectors[0].umi_metrics()` would therefore be non-empty here,
    /// failing the first assertion below. With the fix (the reducer called
    /// once with the FULL 20-length fractions array), only `collectors[19]`
    /// (the 100% slot, where `idx == last_fraction_idx` is genuinely true)
    /// gets UMI metrics.
    #[test]
    fn record_coordinate_group_only_calls_record_umi_at_full_fraction() {
        let mut acc = ConsensusMetricsAccumulator::new_simplex();
        let group = vec![template("1", "chr1", 100, 0.05), template("1", "chr1", 100, 0.05)];
        acc.record_coordinate_group(&group, &[]).expect("records");

        let ConsensusMetricsAccumulator::Simplex { collectors, .. } = &acc else {
            panic!("expected Simplex variant");
        };

        assert!(
            collectors[0].umi_metrics().is_empty(),
            "sub-1.0 fraction (5%) must not record UMI metrics"
        );
        assert!(!collectors[19].umi_metrics().is_empty(), "100% fraction must record UMI metrics");
    }
}

#[cfg(test)]
mod coordinate_group_fragment_tests {
    use super::*;

    #[test]
    fn coordinate_group_fragment_ordinal_is_its_batch_serial() {
        let fragment = CoordinateGroupFragment { batch_serial: 7, entries: vec![] };
        assert_eq!(fragment.ordinal(), 7);
    }
}

#[cfg(test)]
mod coordinate_group_collector_tests {
    use super::*;

    /// Fresh `CoordinateGroupCollector` over a simplex accumulator, plus a
    /// handle to observe the accumulator `finish` hands back.
    fn observing_collector() -> (
        CoordinateGroupCollector,
        std::sync::Arc<std::sync::Mutex<Option<ConsensusMetricsAccumulator>>>,
    ) {
        let observed = std::sync::Arc::new(std::sync::Mutex::new(None));
        let observed_clone = std::sync::Arc::clone(&observed);
        let collector = CoordinateGroupCollector::new(
            ConsensusMetricsAccumulator::new_simplex(),
            vec![],
            Box::new(move |acc| {
                *observed_clone.lock().unwrap() = Some(acc);
                Ok(())
            }),
        );
        (collector, observed)
    }

    fn size_1_family_count(acc: &ConsensusMetricsAccumulator) -> usize {
        let ConsensusMetricsAccumulator::Simplex { collectors, .. } = acc else {
            panic!("expected Simplex accumulator");
        };
        collectors[19]
            .family_size_metrics()
            .iter()
            .filter(|m| m.family_size == 1)
            .map(|m| m.cs_count)
            .sum()
    }

    fn size_2_family_count(acc: &ConsensusMetricsAccumulator) -> usize {
        let ConsensusMetricsAccumulator::Simplex { collectors, .. } = acc else {
            panic!("expected Simplex accumulator");
        };
        collectors[19]
            .family_size_metrics()
            .iter()
            .filter(|m| m.family_size == 2)
            .map(|m| m.cs_count)
            .sum()
    }

    #[test]
    fn one_fragment_per_group_is_the_common_case() {
        let (mut collector, observed) = observing_collector();
        collector
            .record(CoordinateGroupFragment {
                batch_serial: 0,
                entries: vec![(template("0", "chr1", 100, 1.0), key(0, 100))],
            })
            .expect("records");
        collector.finish().expect("finishes");
        assert_eq!(size_1_family_count(observed.lock().unwrap().as_ref().unwrap()), 1);
    }

    /// The `GroupByMi`-batch-boundary-split case that motivated the whole
    /// redesign — now just "the buffer doesn't flush until the key changes,"
    /// because the `ReorderStage` upstream (Task 11) guarantees these two
    /// fragments arrive in true emission order, adjacent. This is the core
    /// ordered-stream correctness claim: a coordinate group split across
    /// MULTIPLE fragments/batches must aggregate to ONE group when delivered
    /// in order.
    #[test]
    fn a_group_split_across_consecutive_fragments_with_the_same_key_stays_one_family() {
        let (mut collector, observed) = observing_collector();
        collector
            .record(CoordinateGroupFragment {
                batch_serial: 0,
                entries: vec![(template("0", "chr1", 100, 1.0), key(0, 100))],
            })
            .expect("records");
        collector
            .record(CoordinateGroupFragment {
                batch_serial: 1,
                entries: vec![(template("1", "chr1", 100, 1.0), key(0, 100))],
            })
            .expect("records");
        collector.finish().expect("finishes");
        let acc = observed.lock().unwrap();
        let acc = acc.as_ref().unwrap();
        assert_eq!(
            size_2_family_count(acc),
            1,
            "both fragments' templates must land in one CS family of size 2"
        );
        assert_eq!(size_1_family_count(acc), 0);
    }

    /// The X5-001 placeholder case (spec §3): every input batch must push
    /// SOMETHING on the third branch, even an empty fragment, so the reorder
    /// stage never stalls. Assert it doesn't spuriously flush.
    #[test]
    fn an_empty_fragment_in_the_middle_of_the_stream_is_a_pure_no_op() {
        let (mut collector, observed) = observing_collector();
        collector
            .record(CoordinateGroupFragment {
                batch_serial: 0,
                entries: vec![(template("0", "chr1", 100, 1.0), key(0, 100))],
            })
            .expect("records");
        collector
            .record(CoordinateGroupFragment { batch_serial: 1, entries: vec![] })
            .expect("records");
        collector
            .record(CoordinateGroupFragment {
                batch_serial: 2,
                entries: vec![(template("1", "chr1", 100, 1.0), key(0, 100))],
            })
            .expect("records");
        collector.finish().expect("finishes");
        assert_eq!(
            size_2_family_count(observed.lock().unwrap().as_ref().unwrap()),
            1,
            "the empty fragment must not break up the group spanning around it"
        );
    }

    #[test]
    fn finish_flushes_whatever_is_still_buffered_at_end_of_stream() {
        let (mut collector, observed) = observing_collector();
        collector
            .record(CoordinateGroupFragment {
                batch_serial: 0,
                entries: vec![(template("0", "chr1", 100, 1.0), key(0, 100))],
            })
            .expect("records");
        // No key change ever occurs — the group is still open when finish() runs.
        collector.finish().expect("finishes");
        assert_eq!(
            size_1_family_count(observed.lock().unwrap().as_ref().unwrap()),
            1,
            "finish() must flush the still-open group, not drop it"
        );
    }

    /// A key change across TWO different `ReadInfoKey`s (rather than the
    /// same one repeated) must close the first group and open a second —
    /// pins that the buffer-then-flush-on-key-change branch itself (not just
    /// the same-key and end-of-stream branches above) works correctly across
    /// fragment boundaries.
    #[test]
    fn a_genuine_key_change_across_fragments_closes_one_group_and_opens_another() {
        let (mut collector, observed) = observing_collector();
        collector
            .record(CoordinateGroupFragment {
                batch_serial: 0,
                entries: vec![(template("0", "chr1", 100, 1.0), key(0, 100))],
            })
            .expect("records");
        collector
            .record(CoordinateGroupFragment {
                batch_serial: 1,
                entries: vec![(template("1", "chr2", 200, 1.0), key(1, 200))],
            })
            .expect("records");
        collector.finish().expect("finishes");
        let acc = observed.lock().unwrap();
        let acc = acc.as_ref().unwrap();
        assert_eq!(size_1_family_count(acc), 2, "two distinct size-1 families, one per key");
        assert_eq!(size_2_family_count(acc), 0);
    }
}
