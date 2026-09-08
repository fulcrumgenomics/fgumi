//! Inline consensus-metrics accumulator: the per-worker state fed by the
//! fused (T1, Task 9) and standalone (T2, Task 8) adapters to compute the
//! same CS/SS/DS family-size, UMI, and downsampling-yield metrics the
//! separate-pass simplex-metrics/duplex-metrics commands compute, inline,
//! during the existing streaming pass. See
//! docs/superpowers/specs/2026-09-08-inline-consensus-metrics-design.md §3-§6.
//!
//! **Ported ahead of its caller.** Every item here is `pub(crate)`. Until
//! Task 8's `Serial` collector step and Task 9's fused T1 accumulator wire
//! this module in, nothing outside its own tests constructs an
//! `InlineCollector`/`ConsensusMetricsAccumulator`, so the whole module is
//! dead code to the compiler. The `allow` below is scoped to this module and
//! should be deleted once a caller lands (Task 8 or 9).

#![allow(dead_code)]

use crate::commands::shared_metrics::{
    DOWNSAMPLING_FRACTIONS, Interval, TemplateInfo, overlaps_intervals,
    record_duplex_coordinate_group, record_simplex_coordinate_group,
};
use crate::simple_umi_consensus::SimpleUmiConsensusCaller;
use anyhow::Result;
use fgumi_metrics::duplex::DuplexMetricsCollector;
use fgumi_metrics::simplex::SimplexMetricsCollector;

/// One mode's per-fraction collector. `Duplex` covers both `duplex` and
/// `codec` modes — both are DS-capable and write the `duplex_*` file set.
pub(crate) enum InlineCollector {
    Simplex(SimplexMetricsCollector),
    Duplex(DuplexMetricsCollector),
}

impl InlineCollector {
    fn record(
        &mut self,
        group: &[TemplateInfo],
        fractions: &[f64],
        caller: &mut SimpleUmiConsensusCaller,
        fraction_template_counts: &mut [usize],
        idx: usize,
        duplex_umi_counts: bool,
    ) -> Result<()> {
        match self {
            InlineCollector::Simplex(c) => record_simplex_coordinate_group(
                group,
                &fractions[idx..=idx],
                std::slice::from_mut(c),
                caller,
                &mut fraction_template_counts[idx..=idx],
            ),
            InlineCollector::Duplex(c) => record_duplex_coordinate_group(
                group,
                &fractions[idx..=idx],
                std::slice::from_mut(c),
                caller,
                &mut fraction_template_counts[idx..=idx],
                duplex_umi_counts,
            ),
        }
    }
}

/// Accumulator state for one consensus stage's inline metrics. Shared by T1
/// (wrapped in `PerThreadAccumulator`, N sharded instances folded via
/// `merge()`) and T2 (a single un-sharded instance owned directly by the
/// `Serial` collector step, Task 8 — `merge()` is never called on that path).
pub(crate) struct ConsensusMetricsAccumulator {
    fraction_collectors: [InlineCollector; 20],
    umi_caller: SimpleUmiConsensusCaller,
    /// Mirrors `duplex-metrics --duplex-umi-counts` (`self.duplex_umi_counts`
    /// on `DuplexMetrics`, `duplex_metrics.rs:104-105`), threaded into
    /// `record_duplex_coordinate_group`'s `duplex_umi_counts` parameter
    /// (Task 5). Ignored (unused) when `fraction_collectors` holds
    /// `InlineCollector::Simplex`.
    duplex_umi_counts: bool,
}

impl ConsensusMetricsAccumulator {
    pub(crate) fn new_simplex() -> Self {
        Self {
            fraction_collectors: std::array::from_fn(|_| {
                InlineCollector::Simplex(SimplexMetricsCollector::new())
            }),
            duplex_umi_counts: false, // unused by the Simplex variant
            umi_caller: SimpleUmiConsensusCaller::default(),
        }
    }

    pub(crate) fn new_duplex(collect_duplex_umi_counts: bool) -> Self {
        Self {
            fraction_collectors: std::array::from_fn(|_| {
                InlineCollector::Duplex(DuplexMetricsCollector::new(collect_duplex_umi_counts))
            }),
            duplex_umi_counts: collect_duplex_umi_counts,
            umi_caller: SimpleUmiConsensusCaller::default(),
        }
    }

    /// Records one already-assembled coordinate/strand group, applying the
    /// per-template interval filter itself (matching
    /// `process_templates_from_bam`'s ordering: interval filter strictly
    /// before any fraction bucketing) — callers pass the RAW group here.
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
        let mut fraction_template_counts = vec![0usize; DOWNSAMPLING_FRACTIONS.len()];
        for (idx, collector) in self.fraction_collectors.iter_mut().enumerate() {
            collector.record(
                &filtered,
                &DOWNSAMPLING_FRACTIONS,
                &mut self.umi_caller,
                &mut fraction_template_counts,
                idx,
                self.duplex_umi_counts,
            )?;
        }
        Ok(())
    }

    /// Folds `other` into `self` by zipping and merging the 20 per-fraction
    /// collectors pairwise. **T1 only** — T2's `Serial` collector step
    /// (Task 8) owns exactly one instance and never calls this.
    ///
    /// `self.duplex_umi_counts` is left untouched — like
    /// `DuplexMetricsCollector::merge`'s `collect_duplex_umi_counts` field
    /// (Task 3), it is a per-stage config flag, not accumulated state; every
    /// accumulator for one stage is always constructed with the same value.
    pub(crate) fn merge(&mut self, other: Self) -> Result<()> {
        for (mine, theirs) in self.fraction_collectors.iter_mut().zip(other.fraction_collectors) {
            match (mine, theirs) {
                (InlineCollector::Simplex(m), InlineCollector::Simplex(t)) => m.merge(t),
                (InlineCollector::Duplex(m), InlineCollector::Duplex(t)) => m.merge(t),
                _ => anyhow::bail!(
                    "internal error: attempted to merge a Simplex accumulator with a Duplex \
                     one — every ConsensusMetricsAccumulator for one stage must be constructed \
                     via the same new_simplex()/new_duplex() call, never mixed"
                ),
            }
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

        let InlineCollector::Simplex(collector) = &acc.fraction_collectors[19] else {
            panic!("expected Simplex variant");
        };
        let metrics = collector.family_size_metrics();
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

        let InlineCollector::Simplex(collector) = &a.fraction_collectors[19] else {
            panic!("expected Simplex variant");
        };
        let metrics = collector.family_size_metrics();
        let cs_count_at_size_1: usize =
            metrics.iter().filter(|m| m.family_size == 1).map(|m| m.cs_count).sum();
        assert_eq!(
            cs_count_at_size_1, 2,
            "both accumulators' size-1 families must be present after merge"
        );
    }
}
