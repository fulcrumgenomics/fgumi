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

use crate::commands::shared_metrics::{
    DOWNSAMPLING_FRACTIONS, Interval, TemplateInfo, overlaps_intervals,
    record_duplex_coordinate_group, record_simplex_coordinate_group,
};
use crate::simple_umi_consensus::SimpleUmiConsensusCaller;
use anyhow::Result;
use fgumi_metrics::duplex::DuplexMetricsCollector;
use fgumi_metrics::simplex::SimplexMetricsCollector;

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

#[cfg(test)]
mod tests {
    use super::*;

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
