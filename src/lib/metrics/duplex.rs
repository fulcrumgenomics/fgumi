//! Metrics for duplex sequencing QC.
//!
//! This module provides comprehensive QC metrics for duplex sequencing experiments,
//! including family size distributions, UMI frequencies, and duplex yield at multiple
//! sampling levels.
//!
//! Used by both `duplex` and `duplex-metrics` commands.

use crate::simple_umi_consensus::SimpleUmiConsensusCaller;
use crate::umi::extract_mi_base;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use super::Metric;

/// Metrics quantifying the distribution of different kinds of read family sizes.
///
/// Three kinds of families are described:
/// - **CS** (Coordinate & Strand): families grouped by unclipped 5' genomic positions and strands
/// - **SS** (Single Strand): single-strand families using UMIs, not linking opposing strands
/// - **DS** (Double Strand): families combining single-strand families from opposite strands
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FamilySizeMetric {
    /// The family size (number of read pairs grouped together)
    pub family_size: usize,
    /// Count of CS families with this size
    pub cs_count: usize,
    /// Fraction of all CS families with this size
    pub cs_fraction: f64,
    /// Fraction of CS families with size >= `family_size`
    pub cs_fraction_gt_or_eq_size: f64,
    /// Count of SS families with this size
    pub ss_count: usize,
    /// Fraction of all SS families with this size
    pub ss_fraction: f64,
    /// Fraction of SS families with size >= `family_size`
    pub ss_fraction_gt_or_eq_size: f64,
    /// Count of DS families with this size
    pub ds_count: usize,
    /// Fraction of all DS families with this size
    pub ds_fraction: f64,
    /// Fraction of DS families with size >= `family_size`
    pub ds_fraction_gt_or_eq_size: f64,
}

impl FamilySizeMetric {
    /// Creates a new family size metric
    #[must_use]
    pub fn new(family_size: usize) -> Self {
        Self {
            family_size,
            cs_count: 0,
            cs_fraction: 0.0,
            cs_fraction_gt_or_eq_size: 0.0,
            ss_count: 0,
            ss_fraction: 0.0,
            ss_fraction_gt_or_eq_size: 0.0,
            ds_count: 0,
            ds_fraction: 0.0,
            ds_fraction_gt_or_eq_size: 0.0,
        }
    }
}

impl Default for FamilySizeMetric {
    fn default() -> Self {
        Self::new(0)
    }
}

impl Metric for FamilySizeMetric {
    fn metric_name() -> &'static str {
        "duplex family size"
    }
}

/// Metrics describing double-stranded (duplex) tag families by AB and BA strand sizes.
///
/// For a given tag family, `ab` is the larger sub-family and `ba` is the smaller one.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DuplexFamilySizeMetric {
    /// Number of reads in the AB sub-family (larger)
    pub ab_size: usize,
    /// Number of reads in the BA sub-family (smaller)
    pub ba_size: usize,
    /// Count of families with these AB/BA sizes
    pub count: usize,
    /// Fraction of all duplex families with these sizes
    pub fraction: f64,
    /// Fraction of duplex families with AB >= `ab_size` and BA >= `ba_size`
    pub fraction_gt_or_eq_size: f64,
}

impl DuplexFamilySizeMetric {
    /// Creates a new duplex family size metric
    #[must_use]
    pub fn new(ab_size: usize, ba_size: usize) -> Self {
        Self { ab_size, ba_size, count: 0, fraction: 0.0, fraction_gt_or_eq_size: 0.0 }
    }
}

impl Default for DuplexFamilySizeMetric {
    fn default() -> Self {
        Self::new(0, 0)
    }
}

impl Metric for DuplexFamilySizeMetric {
    fn metric_name() -> &'static str {
        "duplex AB/BA family size"
    }
}

impl Ord for DuplexFamilySizeMetric {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.ab_size.cmp(&other.ab_size).then_with(|| self.ba_size.cmp(&other.ba_size))
    }
}

impl PartialOrd for DuplexFamilySizeMetric {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for DuplexFamilySizeMetric {}

impl PartialEq for DuplexFamilySizeMetric {
    fn eq(&self, other: &Self) -> bool {
        self.ab_size == other.ab_size && self.ba_size == other.ba_size
    }
}

/// Metrics sampled at various levels of coverage via random downsampling.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DuplexYieldMetric {
    /// Approximate fraction of full dataset used
    pub fraction: f64,
    /// Number of read pairs upon which metrics are based
    pub read_pairs: usize,
    /// Number of CS (Coordinate & Strand) families
    pub cs_families: usize,
    /// Number of SS (Single-Strand by UMI) families
    pub ss_families: usize,
    /// Number of DS (Double-Strand by UMI) families
    pub ds_families: usize,
    /// Number of DS families that are duplexes (min reads on both strands)
    pub ds_duplexes: usize,
    /// Fraction of DS families that are duplexes
    pub ds_fraction_duplexes: f64,
    /// Expected fraction of DS families that should be duplexes under ideal model
    pub ds_fraction_duplexes_ideal: f64,
}

impl DuplexYieldMetric {
    /// Creates a new yield metric
    #[must_use]
    pub fn new(fraction: f64) -> Self {
        Self {
            fraction,
            read_pairs: 0,
            cs_families: 0,
            ss_families: 0,
            ds_families: 0,
            ds_duplexes: 0,
            ds_fraction_duplexes: 0.0,
            ds_fraction_duplexes_ideal: 0.0,
        }
    }
}

impl Default for DuplexYieldMetric {
    fn default() -> Self {
        Self::new(0.0)
    }
}

impl Metric for DuplexYieldMetric {
    fn metric_name() -> &'static str {
        "duplex yield"
    }
}

/// Metrics describing observed UMI sequences and their observation frequencies.
///
/// UMI sequences may be corrected using information within a double-stranded tag family.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UmiMetric {
    /// The UMI sequence (possibly corrected)
    pub umi: String,
    /// Number of read pairs observing this UMI (after correction)
    pub raw_observations: usize,
    /// Subset of raw observations that underwent correction
    pub raw_observations_with_errors: usize,
    /// Number of double-stranded tag families observing this UMI
    pub unique_observations: usize,
    /// Fraction of all raw observations
    pub fraction_raw_observations: f64,
    /// Fraction of all unique observations
    pub fraction_unique_observations: f64,
}

impl UmiMetric {
    /// Creates a new UMI metric
    #[must_use]
    pub fn new(umi: String) -> Self {
        Self {
            umi,
            raw_observations: 0,
            raw_observations_with_errors: 0,
            unique_observations: 0,
            fraction_raw_observations: 0.0,
            fraction_unique_observations: 0.0,
        }
    }
}

impl Default for UmiMetric {
    fn default() -> Self {
        Self::new(String::new())
    }
}

impl Metric for UmiMetric {
    fn metric_name() -> &'static str {
        "UMI"
    }
}

/// Metrics describing observed duplex UMI sequences and their frequencies.
///
/// Duplex UMIs are normalized to F1R2 orientation (positive strand first).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DuplexUmiMetric {
    /// The duplex UMI sequence (possibly corrected, F1R2 normalized)
    pub umi: String,
    /// Number of read pairs observing this duplex UMI
    pub raw_observations: usize,
    /// Subset of raw observations that underwent correction
    pub raw_observations_with_errors: usize,
    /// Number of double-stranded tag families observing this duplex UMI
    pub unique_observations: usize,
    /// Fraction of all raw observations
    pub fraction_raw_observations: f64,
    /// Fraction of all unique observations
    pub fraction_unique_observations: f64,
    /// Expected fraction based on individual UMI frequencies
    pub fraction_unique_observations_expected: f64,
}

impl DuplexUmiMetric {
    /// Creates a new duplex UMI metric
    #[must_use]
    pub fn new(umi: String) -> Self {
        Self {
            umi,
            raw_observations: 0,
            raw_observations_with_errors: 0,
            unique_observations: 0,
            fraction_raw_observations: 0.0,
            fraction_unique_observations: 0.0,
            fraction_unique_observations_expected: 0.0,
        }
    }
}

impl Default for DuplexUmiMetric {
    fn default() -> Self {
        Self::new(String::new())
    }
}

impl Metric for DuplexUmiMetric {
    fn metric_name() -> &'static str {
        "duplex UMI"
    }
}

/// Collector for duplex sequencing metrics.
///
/// Tracks family sizes, UMI frequencies, and yields at multiple sampling levels.
pub struct DuplexMetricsCollector {
    /// Minimum AB reads for duplex
    _min_ab_reads: usize,
    /// Minimum BA reads for duplex
    _min_ba_reads: usize,
    /// Whether to collect duplex UMI counts (memory intensive)
    collect_duplex_umi_counts: bool,

    // Family size tracking
    cs_family_sizes: HashMap<usize, usize>,
    ss_family_sizes: HashMap<usize, usize>,
    ds_family_sizes: HashMap<usize, usize>,
    duplex_family_sizes: HashMap<(usize, usize), usize>,

    // UMI tracking
    umi_raw_counts: HashMap<String, usize>,
    umi_raw_error_counts: HashMap<String, usize>,
    umi_unique_counts: HashMap<String, usize>,
    duplex_umi_raw_counts: HashMap<String, usize>,
    duplex_umi_raw_error_counts: HashMap<String, usize>,
    duplex_umi_unique_counts: HashMap<String, usize>,

    // Yield metrics at different sampling fractions
    yield_metrics: Vec<DuplexYieldMetric>,
}

#[allow(clippy::cast_precision_loss)]
impl DuplexMetricsCollector {
    /// Creates a new metrics collector
    #[must_use]
    pub fn new(min_ab_reads: usize, min_ba_reads: usize, collect_duplex_umi_counts: bool) -> Self {
        Self {
            _min_ab_reads: min_ab_reads,
            _min_ba_reads: min_ba_reads,
            collect_duplex_umi_counts,
            cs_family_sizes: HashMap::new(),
            ss_family_sizes: HashMap::new(),
            ds_family_sizes: HashMap::new(),
            duplex_family_sizes: HashMap::new(),
            umi_raw_counts: HashMap::new(),
            umi_raw_error_counts: HashMap::new(),
            umi_unique_counts: HashMap::new(),
            duplex_umi_raw_counts: HashMap::new(),
            duplex_umi_raw_error_counts: HashMap::new(),
            duplex_umi_unique_counts: HashMap::new(),
            yield_metrics: Vec::new(),
        }
    }

    /// Records a CS (Coordinate+Strand) family
    pub fn record_cs_family(&mut self, size: usize) {
        *self.cs_family_sizes.entry(size).or_insert(0) += 1;
    }

    /// Records an SS (Single-Strand) family
    pub fn record_ss_family(&mut self, size: usize) {
        *self.ss_family_sizes.entry(size).or_insert(0) += 1;
    }

    /// Records a DS (Double-Strand) family
    pub fn record_ds_family(&mut self, size: usize) {
        *self.ds_family_sizes.entry(size).or_insert(0) += 1;
    }

    /// Records a duplex family with AB and BA sizes
    pub fn record_duplex_family(&mut self, ab_size: usize, ba_size: usize) {
        // Ensure ab >= ba by convention
        let (ab, ba) = if ab_size >= ba_size { (ab_size, ba_size) } else { (ba_size, ab_size) };
        *self.duplex_family_sizes.entry((ab, ba)).or_insert(0) += 1;
    }

    /// Records a UMI observation
    ///
    /// # Arguments
    /// * `umi` - The consensus UMI sequence
    /// * `raw_count` - Number of raw observations of this UMI
    /// * `error_count` - Number of raw observations that had errors (differed from consensus)
    /// * `is_unique` - Whether this is a unique family observation
    pub fn record_umi(&mut self, umi: &str, raw_count: usize, error_count: usize, is_unique: bool) {
        *self.umi_raw_counts.entry(umi.to_string()).or_insert(0) += raw_count;
        *self.umi_raw_error_counts.entry(umi.to_string()).or_insert(0) += error_count;
        if is_unique {
            *self.umi_unique_counts.entry(umi.to_string()).or_insert(0) += 1;
        }
    }

    /// Merge another collector's metrics into this one.
    ///
    /// Used for aggregating metrics from parallel batch processing.
    pub fn merge(&mut self, other: &Self) {
        for (size, count) in &other.cs_family_sizes {
            *self.cs_family_sizes.entry(*size).or_insert(0) += count;
        }
        for (size, count) in &other.ss_family_sizes {
            *self.ss_family_sizes.entry(*size).or_insert(0) += count;
        }
        for (size, count) in &other.ds_family_sizes {
            *self.ds_family_sizes.entry(*size).or_insert(0) += count;
        }
        for (sizes, count) in &other.duplex_family_sizes {
            *self.duplex_family_sizes.entry(*sizes).or_insert(0) += count;
        }
        for (umi, count) in &other.umi_raw_counts {
            *self.umi_raw_counts.entry(umi.clone()).or_insert(0) += count;
        }
        for (umi, count) in &other.umi_raw_error_counts {
            *self.umi_raw_error_counts.entry(umi.clone()).or_insert(0) += count;
        }
        for (umi, count) in &other.umi_unique_counts {
            *self.umi_unique_counts.entry(umi.clone()).or_insert(0) += count;
        }
        for (umi, count) in &other.duplex_umi_raw_counts {
            *self.duplex_umi_raw_counts.entry(umi.clone()).or_insert(0) += count;
        }
        for (umi, count) in &other.duplex_umi_raw_error_counts {
            *self.duplex_umi_raw_error_counts.entry(umi.clone()).or_insert(0) += count;
        }
        for (umi, count) in &other.duplex_umi_unique_counts {
            *self.duplex_umi_unique_counts.entry(umi.clone()).or_insert(0) += count;
        }
    }

    /// Records a duplex UMI observation
    ///
    /// # Arguments
    /// * `umi` - The duplex consensus UMI sequence
    /// * `raw_count` - Number of raw observations of this duplex UMI
    /// * `error_count` - Number of raw observations that had errors
    /// * `is_unique` - Whether this is a unique family observation
    pub fn record_duplex_umi(
        &mut self,
        umi: &str,
        raw_count: usize,
        error_count: usize,
        is_unique: bool,
    ) {
        if !self.collect_duplex_umi_counts {
            return;
        }
        *self.duplex_umi_raw_counts.entry(umi.to_string()).or_insert(0) += raw_count;
        *self.duplex_umi_raw_error_counts.entry(umi.to_string()).or_insert(0) += error_count;
        if is_unique {
            *self.duplex_umi_unique_counts.entry(umi.to_string()).or_insert(0) += 1;
        }
    }

    /// Updates UMI metrics for a duplex family.
    ///
    /// This method:
    /// 1. Uses RX tags (raw UMI sequences) to extract individual UMI observations
    /// 2. Separates by strand (/A and /B suffixes in MI tags), swapping UMI parts for B strand
    /// 3. Calls consensus for each UMI position
    /// 4. Records raw observations, errors, and unique observations for each individual UMI
    /// 5. Records duplex UMI metrics if enabled
    ///
    /// # Arguments
    /// * `umi_consensus_caller` - Caller for computing UMI consensus sequences
    /// * `mi_rx_pairs` - Pairs of (MI tag value, RX tag value) for each read in the family
    /// * `base_mi` - The base MI (without /A or /B suffix) for this family
    pub fn update_umi_metrics_for_family(
        &mut self,
        umi_consensus_caller: &mut SimpleUmiConsensusCaller,
        mi_rx_pairs: &[(String, String)],
        base_mi: &str,
    ) {
        // Collect individual UMI parts from each strand's RX tag
        // For /A strand: split RX "AAA-TTT" → umi1s += "AAA", umi2s += "TTT"
        // For /B strand: split RX "TTT-AAA" → umi1s += "AAA", umi2s += "TTT" (swapped)
        let mut umi1s = Vec::new();
        let mut umi2s = Vec::new();

        for (mi, rx) in mi_rx_pairs {
            // Check if this MI tag belongs to the current base_mi family
            let mi_base = extract_mi_base(mi);

            if mi_base != base_mi {
                continue;
            }

            // Split the RX tag to get individual UMI parts
            let parts: Vec<&str> = rx.split('-').collect();
            if parts.len() != 2 {
                // Not a valid duplex UMI, skip
                continue;
            }

            // Check that both components are non-empty
            if parts[0].is_empty() || parts[1].is_empty() {
                // Empty component, skip
                continue;
            }

            // Add UMI parts based on strand
            if mi.ends_with("/A") {
                // For /A strand: u1 goes to umi1s, u2 goes to umi2s
                umi1s.push(parts[0].to_string());
                umi2s.push(parts[1].to_string());
            } else if mi.ends_with("/B") {
                // For /B strand: u2 goes to umi1s, u1 goes to umi2s (swapped)
                umi1s.push(parts[1].to_string());
                umi2s.push(parts[0].to_string());
            }
        }

        // Call consensus for each UMI position and record metrics
        let mut consensus_umis = Vec::new();

        if !umi1s.is_empty() {
            let (consensus, _had_errors) = umi_consensus_caller.consensus(&umi1s);
            let raw_count = umi1s.len();
            let error_count = umi1s.iter().filter(|u| **u != consensus).count();
            self.record_umi(&consensus, raw_count, error_count, true);
            consensus_umis.push(consensus);
        }

        if !umi2s.is_empty() {
            let (consensus, _had_errors) = umi_consensus_caller.consensus(&umi2s);
            let raw_count = umi2s.len();
            let error_count = umi2s.iter().filter(|u| **u != consensus).count();
            self.record_umi(&consensus, raw_count, error_count, true);
            consensus_umis.push(consensus);
        }

        // Record duplex UMI metrics if enabled
        if self.collect_duplex_umi_counts && consensus_umis.len() == 2 {
            let duplex_umi = format!("{}-{}", consensus_umis[0], consensus_umis[1]);
            // Each read pair contributes one observation to the duplex UMI
            let total_raw = umi1s.len();

            // Count how many raw RX tags had errors (don't match either duplex orientation)
            let expected_duplex1 = format!("{}-{}", consensus_umis[0], consensus_umis[1]);
            let expected_duplex2 = format!("{}-{}", consensus_umis[1], consensus_umis[0]);
            let error_count = mi_rx_pairs
                .iter()
                .filter(|(mi, rx)| {
                    let mi_base = extract_mi_base(mi);
                    mi_base == base_mi && *rx != expected_duplex1 && *rx != expected_duplex2
                })
                .count();

            self.record_duplex_umi(&duplex_umi, total_raw, error_count, true);
        }
    }

    /// Updates UMI metrics using the consensus UMI from the main caller's output.
    ///
    /// This method takes the consensus UMI that was already computed by `DuplexConsensusCaller`
    /// and compares it against the raw UMI observations to count errors.
    ///
    /// # Arguments
    /// * `consensus_rx` - The consensus UMI from the output read's RX tag (e.g., "AAT-CCG")
    /// * `strand_rx_pairs` - Raw UMI observations: (is_a_strand, rx_bytes) for each read
    pub fn update_umi_metrics_from_consensus(
        &mut self,
        consensus_rx: &[u8],
        strand_rx_pairs: &[(bool, Vec<u8>)],
    ) {
        // Parse consensus UMI into two parts
        let Some(consensus_sep) = consensus_rx.iter().position(|&b| b == b'-') else {
            return;
        };
        let consensus_umi1 = &consensus_rx[..consensus_sep];
        let consensus_umi2 = &consensus_rx[consensus_sep + 1..];

        if consensus_umi1.is_empty() || consensus_umi2.is_empty() {
            return;
        }

        // Count raw observations and errors for each UMI position
        let mut umi1_count = 0usize;
        let mut umi1_errors = 0usize;
        let mut umi2_count = 0usize;
        let mut umi2_errors = 0usize;

        for (is_a_strand, rx) in strand_rx_pairs {
            // Split raw RX tag
            let Some(sep_pos) = rx.iter().position(|&b| b == b'-') else {
                continue;
            };
            let part1 = &rx[..sep_pos];
            let part2 = &rx[sep_pos + 1..];

            if part1.is_empty() || part2.is_empty() {
                continue;
            }

            // Normalize by strand (B strand is swapped)
            let (raw_umi1, raw_umi2) = if *is_a_strand {
                (part1, part2)
            } else {
                (part2, part1)
            };

            // Count UMI1
            umi1_count += 1;
            if raw_umi1 != consensus_umi1 {
                umi1_errors += 1;
            }

            // Count UMI2
            umi2_count += 1;
            if raw_umi2 != consensus_umi2 {
                umi2_errors += 1;
            }
        }

        // Record metrics for each UMI position
        if umi1_count > 0 {
            let umi1_str = String::from_utf8_lossy(consensus_umi1);
            self.record_umi(&umi1_str, umi1_count, umi1_errors, true);
        }
        if umi2_count > 0 {
            let umi2_str = String::from_utf8_lossy(consensus_umi2);
            self.record_umi(&umi2_str, umi2_count, umi2_errors, true);
        }

        // Record duplex UMI metrics if enabled
        if self.collect_duplex_umi_counts && umi1_count > 0 && umi2_count > 0 {
            let duplex_umi = String::from_utf8_lossy(consensus_rx);
            // Count errors: raw RX doesn't match consensus in either orientation
            let reversed_consensus: Vec<u8> = consensus_umi2
                .iter()
                .chain(std::iter::once(&b'-'))
                .chain(consensus_umi1.iter())
                .copied()
                .collect();

            let duplex_errors = strand_rx_pairs
                .iter()
                .filter(|(_, rx)| rx.as_slice() != consensus_rx && rx.as_slice() != reversed_consensus.as_slice())
                .count();

            self.record_duplex_umi(&duplex_umi, umi1_count, duplex_errors, true);
        }
    }

    /// Generates family size metrics
    #[must_use]
    #[allow(clippy::similar_names)]
    pub fn family_size_metrics(&self) -> Vec<FamilySizeMetric> {
        // Find max family size across all types
        let max_size = *[
            self.cs_family_sizes.keys().max().unwrap_or(&0),
            self.ss_family_sizes.keys().max().unwrap_or(&0),
            self.ds_family_sizes.keys().max().unwrap_or(&0),
        ]
        .iter()
        .max()
        .unwrap();

        let cs_total: usize = self.cs_family_sizes.values().sum();
        let ss_total: usize = self.ss_family_sizes.values().sum();
        let ds_total: usize = self.ds_family_sizes.values().sum();

        let mut metrics = Vec::new();
        for size in 1..=*max_size {
            let mut metric = FamilySizeMetric::new(size);

            metric.cs_count = *self.cs_family_sizes.get(&size).unwrap_or(&0);
            metric.cs_fraction =
                if cs_total > 0 { metric.cs_count as f64 / cs_total as f64 } else { 0.0 };

            metric.ss_count = *self.ss_family_sizes.get(&size).unwrap_or(&0);
            metric.ss_fraction =
                if ss_total > 0 { metric.ss_count as f64 / ss_total as f64 } else { 0.0 };

            metric.ds_count = *self.ds_family_sizes.get(&size).unwrap_or(&0);
            metric.ds_fraction =
                if ds_total > 0 { metric.ds_count as f64 / ds_total as f64 } else { 0.0 };

            metrics.push(metric);
        }

        // Calculate cumulative fractions (>= size)
        for i in (0..metrics.len()).rev() {
            let next_cs =
                if i + 1 < metrics.len() { metrics[i + 1].cs_fraction_gt_or_eq_size } else { 0.0 };
            let next_ss =
                if i + 1 < metrics.len() { metrics[i + 1].ss_fraction_gt_or_eq_size } else { 0.0 };
            let next_ds =
                if i + 1 < metrics.len() { metrics[i + 1].ds_fraction_gt_or_eq_size } else { 0.0 };

            metrics[i].cs_fraction_gt_or_eq_size = metrics[i].cs_fraction + next_cs;
            metrics[i].ss_fraction_gt_or_eq_size = metrics[i].ss_fraction + next_ss;
            metrics[i].ds_fraction_gt_or_eq_size = metrics[i].ds_fraction + next_ds;
        }

        metrics
    }

    /// Generates duplex family size metrics
    #[must_use]
    pub fn duplex_family_size_metrics(&self) -> Vec<DuplexFamilySizeMetric> {
        let total: usize = self.duplex_family_sizes.values().sum();

        let mut metrics: Vec<_> = self
            .duplex_family_sizes
            .iter()
            .map(|((ab, ba), count)| {
                let mut metric = DuplexFamilySizeMetric::new(*ab, *ba);
                metric.count = *count;
                metric.fraction = if total > 0 { *count as f64 / total as f64 } else { 0.0 };
                metric
            })
            .collect();

        metrics.sort();

        // Calculate 2D cumulative fractions: fraction of families with AB >= ab AND BA >= ba
        // This matches fgbio's definition in DuplexFamilySizeMetric
        for metric in &mut metrics {
            let cumulative_count: usize = self
                .duplex_family_sizes
                .iter()
                .filter(|((a, b), _)| *a >= metric.ab_size && *b >= metric.ba_size)
                .map(|(_, count)| count)
                .sum();
            metric.fraction_gt_or_eq_size =
                if total > 0 { cumulative_count as f64 / total as f64 } else { 0.0 };
        }

        metrics
    }

    /// Generates UMI metrics
    #[must_use]
    pub fn umi_metrics(&self) -> Vec<UmiMetric> {
        let total_raw: usize = self.umi_raw_counts.values().sum();
        let total_unique: usize = self.umi_unique_counts.values().sum();

        let mut metrics: Vec<_> = self
            .umi_raw_counts
            .keys()
            .map(|umi| {
                let mut metric = UmiMetric::new(umi.clone());
                metric.raw_observations = *self.umi_raw_counts.get(umi).unwrap_or(&0);
                metric.raw_observations_with_errors =
                    *self.umi_raw_error_counts.get(umi).unwrap_or(&0);
                metric.unique_observations = *self.umi_unique_counts.get(umi).unwrap_or(&0);
                metric.fraction_raw_observations = if total_raw > 0 {
                    metric.raw_observations as f64 / total_raw as f64
                } else {
                    0.0
                };
                metric.fraction_unique_observations = if total_unique > 0 {
                    metric.unique_observations as f64 / total_unique as f64
                } else {
                    0.0
                };
                metric
            })
            .collect();

        // Sort by UMI string (matching fgbio's sortBy(_.umi))
        metrics.sort_by(|a, b| a.umi.cmp(&b.umi));
        metrics
    }

    /// Generates duplex UMI metrics
    ///
    /// # Arguments
    /// * `umi_metrics` - Individual UMI metrics used to calculate expected frequencies
    #[must_use]
    pub fn duplex_umi_metrics(&self, umi_metrics: &[UmiMetric]) -> Vec<DuplexUmiMetric> {
        if !self.collect_duplex_umi_counts {
            return Vec::new();
        }

        // Build a map of individual UMI -> fraction_unique_observations for lookup
        let single_umi_fractions: HashMap<&str, f64> =
            umi_metrics.iter().map(|m| (m.umi.as_str(), m.fraction_unique_observations)).collect();

        let total_raw: usize = self.duplex_umi_raw_counts.values().sum();
        let total_unique: usize = self.duplex_umi_unique_counts.values().sum();

        let mut metrics: Vec<_> = self
            .duplex_umi_raw_counts
            .keys()
            .map(|umi| {
                let mut metric = DuplexUmiMetric::new(umi.clone());
                metric.raw_observations = *self.duplex_umi_raw_counts.get(umi).unwrap_or(&0);
                metric.raw_observations_with_errors =
                    *self.duplex_umi_raw_error_counts.get(umi).unwrap_or(&0);
                metric.unique_observations = *self.duplex_umi_unique_counts.get(umi).unwrap_or(&0);
                metric.fraction_raw_observations = if total_raw > 0 {
                    metric.raw_observations as f64 / total_raw as f64
                } else {
                    0.0
                };
                metric.fraction_unique_observations = if total_unique > 0 {
                    metric.unique_observations as f64 / total_unique as f64
                } else {
                    0.0
                };

                // Calculate expected fraction based on individual UMI frequencies
                // Expected frequency = freq(umi1) * freq(umi2) assuming independence
                metric.fraction_unique_observations_expected =
                    if let Some((umi1, umi2)) = umi.split_once('-') {
                        let freq1 = single_umi_fractions.get(umi1).copied().unwrap_or(0.0);
                        let freq2 = single_umi_fractions.get(umi2).copied().unwrap_or(0.0);
                        freq1 * freq2
                    } else {
                        0.0
                    };

                metric
            })
            .collect();

        // Sort by unique observations descending (matching Scala's sort order)
        metrics.sort_by(|a, b| b.unique_observations.cmp(&a.unique_observations));
        metrics
    }

    /// Returns the yield metrics collected at different sampling fractions
    #[must_use]
    pub fn yield_metrics(&self) -> &[DuplexYieldMetric] {
        &self.yield_metrics
    }
}

// ============================================================================
// Ideal Duplex Fraction Calculation - Shared binomial model
// ============================================================================

use statrs::distribution::{Binomial, DiscreteCDF};

/// Calculates the ideal duplex fraction using binomial model.
///
/// For each family of size N, calculates the probability that both strands
/// have sufficient reads (A >= `min_ab` AND B >= `min_ba` where A + B = N).
/// Assumes each read has 0.5 probability of being on each strand.
///
/// # Arguments
///
/// * `family_sizes` - Slice of DS family sizes
/// * `min_ab` - Minimum reads required on AB strand
/// * `min_ba` - Minimum reads required on BA strand
///
/// # Returns
///
/// Expected fraction of DS families that should be duplexes under ideal model.
#[must_use]
pub fn calculate_ideal_duplex_fraction(
    family_sizes: &[usize],
    min_ab: usize,
    min_ba: usize,
) -> f64 {
    if family_sizes.is_empty() {
        return 0.0;
    }

    let mut ideal_duplexes = 0.0;
    let total_families = family_sizes.len() as f64;

    for &size in family_sizes {
        if size < min_ab + min_ba {
            // Impossible to form a duplex with this family size
            continue;
        }

        // Calculate P(A >= min_ab AND B >= min_ba) where A ~ Binomial(n=size, p=0.5)
        // and B = size - A
        //
        // This is equivalent to: P(min_ba <= A <= size - min_ab)
        // = P(A <= size - min_ab) - P(A < min_ba)
        // = CDF(size - min_ab) - CDF(min_ba - 1)

        let binomial = match Binomial::new(0.5, size as u64) {
            Ok(b) => b,
            Err(_) => continue, // Skip if binomial creation fails
        };

        let upper_bound = size - min_ba;
        let lower_bound = min_ab;

        // P(A >= lower_bound AND A <= upper_bound)
        let prob = if upper_bound >= lower_bound {
            let p_upper = binomial.cdf(upper_bound as u64);
            let p_lower =
                if lower_bound > 0 { binomial.cdf((lower_bound - 1) as u64) } else { 0.0 };
            p_upper - p_lower
        } else {
            0.0
        };

        ideal_duplexes += prob;
    }

    ideal_duplexes / total_families
}

// ============================================================================
// Metrics Writer - Shared logic for writing duplex metrics files
// ============================================================================

use std::sync::OnceLock;

/// Embedded R script for PDF plot generation (bundled with binary)
const R_SCRIPT: &str = include_str!("../../../resources/CollectDuplexSeqMetrics.R");

/// Cached R availability check (computed once per process)
static R_AVAILABLE: OnceLock<bool> = OnceLock::new();

/// Writer for duplex metrics files.
///
/// This struct provides shared logic for writing metrics files from both
/// the `duplex` and `duplex-metrics` commands.
pub struct DuplexMetricsWriter<'a> {
    /// Output prefix for metrics files
    prefix: &'a std::path::Path,
    /// Optional description for PDF plots
    description: Option<&'a str>,
}

impl<'a> DuplexMetricsWriter<'a> {
    /// Create a new metrics writer with the given output prefix.
    #[must_use]
    pub fn new(prefix: &'a std::path::Path) -> Self {
        Self { prefix, description: None }
    }

    /// Set the description for PDF plot titles.
    #[must_use]
    pub fn with_description(mut self, description: Option<&'a str>) -> Self {
        self.description = description;
        self
    }

    /// Write all metrics files from a collector.
    ///
    /// Writes:
    /// - `PREFIX.family_sizes.txt`
    /// - `PREFIX.duplex_family_sizes.txt`
    /// - `PREFIX.duplex_yield_metrics.txt`
    /// - `PREFIX.umi_counts.txt`
    /// - `PREFIX.duplex_qc.pdf` (if R is available)
    ///
    /// # Errors
    ///
    /// Returns an error if any file cannot be written.
    pub fn write_all(
        &self,
        collector: &DuplexMetricsCollector,
        yield_metrics: &[DuplexYieldMetric],
    ) -> anyhow::Result<()> {
        use anyhow::Context;
        use fgoxide::io::DelimFile;
        use log::info;

        // Write family_sizes.txt
        let family_size_metrics = collector.family_size_metrics();
        let family_size_path = format!("{}.family_sizes.txt", self.prefix.display());
        DelimFile::default()
            .write_tsv(&family_size_path, family_size_metrics)
            .with_context(|| format!("Failed to write family size metrics: {family_size_path}"))?;
        info!("Wrote family size metrics to {family_size_path}");

        // Write duplex_family_sizes.txt
        let duplex_family_size_metrics = collector.duplex_family_size_metrics();
        let duplex_family_size_path = format!("{}.duplex_family_sizes.txt", self.prefix.display());
        DelimFile::default()
            .write_tsv(&duplex_family_size_path, duplex_family_size_metrics)
            .with_context(|| {
                format!("Failed to write duplex family size metrics: {duplex_family_size_path}")
            })?;
        info!("Wrote duplex family size metrics to {duplex_family_size_path}");

        // Write umi_counts.txt
        let umi_metrics = collector.umi_metrics();
        let umi_path = format!("{}.umi_counts.txt", self.prefix.display());
        DelimFile::default()
            .write_tsv(&umi_path, umi_metrics.clone())
            .with_context(|| format!("Failed to write UMI metrics: {umi_path}"))?;
        info!("Wrote UMI metrics to {umi_path}");

        // Write duplex_yield_metrics.txt
        let yield_path = format!("{}.duplex_yield_metrics.txt", self.prefix.display());
        DelimFile::default()
            .write_tsv(&yield_path, yield_metrics.to_vec())
            .with_context(|| format!("Failed to write yield metrics: {yield_path}"))?;
        info!("Wrote yield metrics to {yield_path}");

        // Generate PDF plots using R script (optional)
        let pdf_path = format!("{}.duplex_qc.pdf", self.prefix.display());
        if Self::is_r_available() {
            let description = self.description.unwrap_or("Sample");
            match Self::execute_r_script(
                &family_size_path,
                &duplex_family_size_path,
                &yield_path,
                &umi_path,
                &pdf_path,
                description,
            ) {
                Ok(()) => info!("Generated PDF plots: {pdf_path}"),
                Err(e) => {
                    log::warn!("Failed to generate PDF plots: {e}. Continuing without plots.");
                    log::warn!(
                        "To enable PDF generation, ensure R is installed with ggplot2 and scales packages:"
                    );
                    log::warn!("  install.packages(c(\"ggplot2\", \"scales\"))");
                }
            }
        } else {
            log::warn!(
                "R or required packages (ggplot2, scales) not available. Skipping PDF generation."
            );
            log::warn!("To enable PDF generation, install R and required packages:");
            log::warn!("  install.packages(c(\"ggplot2\", \"scales\"))");
        }

        Ok(())
    }

    /// Checks if R and required packages (ggplot2, scales) are available.
    /// Result is cached for the lifetime of the process.
    fn is_r_available() -> bool {
        use std::process::Command;

        *R_AVAILABLE.get_or_init(|| {
            Command::new("Rscript")
                .args(["-e", "stopifnot(require(ggplot2)); stopifnot(require(scales))"])
                .output()
                .map(|output| output.status.success())
                .unwrap_or(false)
        })
    }

    /// Executes the R script to generate PDF plots.
    fn execute_r_script(
        family_size_path: &str,
        duplex_family_size_path: &str,
        yield_path: &str,
        umi_path: &str,
        pdf_path: &str,
        description: &str,
    ) -> anyhow::Result<()> {
        use anyhow::Context;
        use log::info;
        use std::process::Command;

        // Write embedded R script to temp file
        let temp_dir = std::env::temp_dir();
        let r_script_path = temp_dir.join("fgumi_CollectDuplexSeqMetrics.R");
        std::fs::write(&r_script_path, R_SCRIPT)
            .context("Failed to write embedded R script to temp file")?;

        info!("Executing R script to generate PDF plots...");

        let output = Command::new("Rscript")
            .arg(&r_script_path)
            .arg(family_size_path)
            .arg(duplex_family_size_path)
            .arg(yield_path)
            .arg(umi_path)
            .arg(pdf_path)
            .arg(description)
            .output()
            .context("Failed to execute Rscript command")?;

        // Clean up temp file (ignore errors)
        let _ = std::fs::remove_file(&r_script_path);

        if output.status.success() {
            Ok(())
        } else {
            let stderr = String::from_utf8_lossy(&output.stderr);
            anyhow::bail!(
                "R script execution failed with exit code {:?}. Error: {}",
                output.status.code(),
                stderr
            )
        }
    }
}

#[cfg(test)]
#[allow(clippy::cast_precision_loss)]
mod tests {
    use super::*;

    // =========================================================================
    // FamilySizeMetric tests
    // =========================================================================

    #[test]
    fn test_family_size_metric_new() {
        let metric = FamilySizeMetric::new(5);
        assert_eq!(metric.family_size, 5);
        assert_eq!(metric.cs_count, 0);
        assert!(metric.cs_fraction.abs() < f64::EPSILON);
        assert_eq!(metric.ss_count, 0);
        assert_eq!(metric.ds_count, 0);
    }

    // =========================================================================
    // DuplexFamilySizeMetric tests
    // =========================================================================

    #[test]
    fn test_duplex_family_size_metric_new() {
        let metric = DuplexFamilySizeMetric::new(10, 5);
        assert_eq!(metric.ab_size, 10);
        assert_eq!(metric.ba_size, 5);
        assert_eq!(metric.count, 0);
        assert!(metric.fraction.abs() < f64::EPSILON);
    }

    #[test]
    fn test_duplex_family_size_metric_ordering() {
        let m1 = DuplexFamilySizeMetric::new(5, 3);
        let m2 = DuplexFamilySizeMetric::new(5, 4);
        let m3 = DuplexFamilySizeMetric::new(6, 2);

        // m1 < m2 (same ab, smaller ba)
        assert!(m1 < m2);
        // m1 < m3 (smaller ab)
        assert!(m1 < m3);
        // m2 < m3 (smaller ab)
        assert!(m2 < m3);
    }

    #[test]
    fn test_duplex_family_size_metric_equality() {
        let m1 = DuplexFamilySizeMetric::new(5, 3);
        let m2 = DuplexFamilySizeMetric::new(5, 3);
        let m3 = DuplexFamilySizeMetric::new(5, 4);

        assert_eq!(m1, m2);
        assert_ne!(m1, m3);
    }

    // =========================================================================
    // DuplexYieldMetric tests
    // =========================================================================

    #[test]
    fn test_duplex_yield_metric_new() {
        let metric = DuplexYieldMetric::new(0.5);
        assert!((metric.fraction - 0.5).abs() < f64::EPSILON);
        assert_eq!(metric.read_pairs, 0);
        assert_eq!(metric.cs_families, 0);
        assert_eq!(metric.ds_duplexes, 0);
    }

    // =========================================================================
    // UmiMetric tests
    // =========================================================================

    #[test]
    fn test_umi_metric_new() {
        let metric = UmiMetric::new("ACGT".to_string());
        assert_eq!(metric.umi, "ACGT");
        assert_eq!(metric.raw_observations, 0);
        assert_eq!(metric.unique_observations, 0);
    }

    // =========================================================================
    // DuplexUmiMetric tests
    // =========================================================================

    #[test]
    fn test_duplex_umi_metric_new() {
        let metric = DuplexUmiMetric::new("ACGT-TGCA".to_string());
        assert_eq!(metric.umi, "ACGT-TGCA");
        assert_eq!(metric.raw_observations, 0);
        assert!(metric.fraction_unique_observations_expected.abs() < f64::EPSILON);
    }

    // =========================================================================
    // DuplexMetricsCollector tests
    // =========================================================================

    #[test]
    fn test_collector_new() {
        let collector = DuplexMetricsCollector::new(1, 1, true);
        assert!(collector.yield_metrics().is_empty());
    }

    #[test]
    fn test_record_cs_family() {
        let mut collector = DuplexMetricsCollector::new(1, 1, false);
        collector.record_cs_family(5);
        collector.record_cs_family(5);
        collector.record_cs_family(10);

        let metrics = collector.family_size_metrics();
        // Size 5 should have count 2
        let size_5 = metrics.iter().find(|m| m.family_size == 5).unwrap();
        assert_eq!(size_5.cs_count, 2);
        // Size 10 should have count 1
        let size_10 = metrics.iter().find(|m| m.family_size == 10).unwrap();
        assert_eq!(size_10.cs_count, 1);
    }

    #[test]
    fn test_record_ss_family() {
        let mut collector = DuplexMetricsCollector::new(1, 1, false);
        collector.record_ss_family(3);
        collector.record_ss_family(3);
        collector.record_ss_family(3);

        let metrics = collector.family_size_metrics();
        let size_3 = metrics.iter().find(|m| m.family_size == 3).unwrap();
        assert_eq!(size_3.ss_count, 3);
    }

    #[test]
    fn test_record_ds_family() {
        let mut collector = DuplexMetricsCollector::new(1, 1, false);
        collector.record_ds_family(2);

        let metrics = collector.family_size_metrics();
        let size_2 = metrics.iter().find(|m| m.family_size == 2).unwrap();
        assert_eq!(size_2.ds_count, 1);
    }

    #[test]
    fn test_record_duplex_family_normalization() {
        let mut collector = DuplexMetricsCollector::new(1, 1, false);
        // Record with ab < ba - should be normalized
        collector.record_duplex_family(3, 5);
        // Record with ab > ba - already normalized
        collector.record_duplex_family(5, 3);

        let metrics = collector.duplex_family_size_metrics();
        // Both should map to (5, 3)
        assert_eq!(metrics.len(), 1);
        assert_eq!(metrics[0].ab_size, 5);
        assert_eq!(metrics[0].ba_size, 3);
        assert_eq!(metrics[0].count, 2);
    }

    #[test]
    fn test_record_umi() {
        let mut collector = DuplexMetricsCollector::new(1, 1, false);
        collector.record_umi("AAAA", 10, 2, true);
        collector.record_umi("AAAA", 5, 1, false); // Not unique
        collector.record_umi("CCCC", 8, 0, true);

        let metrics = collector.umi_metrics();
        assert_eq!(metrics.len(), 2);

        let aaaa = metrics.iter().find(|m| m.umi == "AAAA").unwrap();
        assert_eq!(aaaa.raw_observations, 15); // 10 + 5
        assert_eq!(aaaa.raw_observations_with_errors, 3); // 2 + 1
        assert_eq!(aaaa.unique_observations, 1); // Only one unique

        let cccc = metrics.iter().find(|m| m.umi == "CCCC").unwrap();
        assert_eq!(cccc.raw_observations, 8);
        assert_eq!(cccc.unique_observations, 1);
    }

    #[test]
    fn test_record_duplex_umi_disabled() {
        let mut collector = DuplexMetricsCollector::new(1, 1, false); // Disabled
        collector.record_duplex_umi("AAAA-TTTT", 10, 0, true);

        let umi_metrics = collector.umi_metrics();
        let duplex_metrics = collector.duplex_umi_metrics(&umi_metrics);
        assert!(duplex_metrics.is_empty());
    }

    #[test]
    fn test_record_duplex_umi_enabled() {
        let mut collector = DuplexMetricsCollector::new(1, 1, true); // Enabled
        collector.record_duplex_umi("AAAA-TTTT", 10, 2, true);

        // Need to also record individual UMIs for expected calculation
        collector.record_umi("AAAA", 5, 0, true);
        collector.record_umi("TTTT", 5, 0, true);

        let umi_metrics = collector.umi_metrics();
        let duplex_metrics = collector.duplex_umi_metrics(&umi_metrics);
        assert_eq!(duplex_metrics.len(), 1);
        assert_eq!(duplex_metrics[0].umi, "AAAA-TTTT");
        assert_eq!(duplex_metrics[0].raw_observations, 10);
    }

    #[test]
    fn test_family_size_metrics_fractions() {
        let mut collector = DuplexMetricsCollector::new(1, 1, false);
        // 4 families total: 2 size-1, 1 size-2, 1 size-3
        collector.record_cs_family(1);
        collector.record_cs_family(1);
        collector.record_cs_family(2);
        collector.record_cs_family(3);

        let metrics = collector.family_size_metrics();

        let size_1 = metrics.iter().find(|m| m.family_size == 1).unwrap();
        assert_eq!(size_1.cs_count, 2);
        assert!((size_1.cs_fraction - 0.5).abs() < 0.001); // 2/4 = 0.5
        assert!((size_1.cs_fraction_gt_or_eq_size - 1.0).abs() < 0.001); // All >= 1

        let size_3 = metrics.iter().find(|m| m.family_size == 3).unwrap();
        assert_eq!(size_3.cs_count, 1);
        assert!((size_3.cs_fraction - 0.25).abs() < 0.001); // 1/4 = 0.25
        assert!((size_3.cs_fraction_gt_or_eq_size - 0.25).abs() < 0.001); // Only size 3
    }

    #[test]
    fn test_duplex_family_size_metrics_sorting() {
        let mut collector = DuplexMetricsCollector::new(1, 1, false);
        collector.record_duplex_family(5, 3);
        collector.record_duplex_family(2, 1);
        collector.record_duplex_family(5, 2);

        let metrics = collector.duplex_family_size_metrics();
        // Should be sorted by (ab, ba)
        assert_eq!(metrics[0].ab_size, 2);
        assert_eq!(metrics[0].ba_size, 1);
        assert_eq!(metrics[1].ab_size, 5);
        assert_eq!(metrics[1].ba_size, 2);
        assert_eq!(metrics[2].ab_size, 5);
        assert_eq!(metrics[2].ba_size, 3);
    }

    #[test]
    fn test_umi_metrics_sorting() {
        let mut collector = DuplexMetricsCollector::new(1, 1, false);
        collector.record_umi("ZZZZ", 1, 0, true);
        collector.record_umi("AAAA", 1, 0, true);
        collector.record_umi("MMMM", 1, 0, true);

        let metrics = collector.umi_metrics();
        // Should be sorted alphabetically
        assert_eq!(metrics[0].umi, "AAAA");
        assert_eq!(metrics[1].umi, "MMMM");
        assert_eq!(metrics[2].umi, "ZZZZ");
    }

    #[test]
    fn test_duplex_umi_expected_frequency() {
        let mut collector = DuplexMetricsCollector::new(1, 1, true);

        // Individual UMIs with known frequencies
        collector.record_umi("AAAA", 10, 0, true);
        collector.record_umi("TTTT", 10, 0, true);
        // Total unique = 2, so each has fraction 0.5

        // Duplex UMI
        collector.record_duplex_umi("AAAA-TTTT", 5, 0, true);

        let umi_metrics = collector.umi_metrics();
        let duplex_metrics = collector.duplex_umi_metrics(&umi_metrics);

        assert_eq!(duplex_metrics.len(), 1);
        // Expected = 0.5 * 0.5 = 0.25
        assert!((duplex_metrics[0].fraction_unique_observations_expected - 0.25).abs() < 0.001);
    }

    #[test]
    fn test_empty_collector() {
        let collector = DuplexMetricsCollector::new(1, 1, false);

        let family_metrics = collector.family_size_metrics();
        assert!(family_metrics.is_empty());

        let duplex_metrics = collector.duplex_family_size_metrics();
        assert!(duplex_metrics.is_empty());

        let umi_metrics = collector.umi_metrics();
        assert!(umi_metrics.is_empty());
    }

    #[test]
    fn test_yield_metrics_empty() {
        let collector = DuplexMetricsCollector::new(1, 1, false);
        assert!(collector.yield_metrics().is_empty());
    }

    #[test]
    fn test_metric_trait_impl() {
        assert_eq!(FamilySizeMetric::metric_name(), "duplex family size");
        assert_eq!(DuplexFamilySizeMetric::metric_name(), "duplex AB/BA family size");
        assert_eq!(DuplexYieldMetric::metric_name(), "duplex yield");
        assert_eq!(UmiMetric::metric_name(), "UMI");
        assert_eq!(DuplexUmiMetric::metric_name(), "duplex UMI");
    }
}
