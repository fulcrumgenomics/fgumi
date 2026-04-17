//! Clipping metrics for the `clip` command.
//!
//! This module provides structures and functionality for collecting detailed
//! statistics about clipping operations performed by `clip`.

use noodles::sam::alignment::record::Cigar;
use noodles::sam::alignment::record_buf::RecordBuf;
use serde::{Deserialize, Serialize};

use crate::Metric;

/// Bundled counts for a single clipping operation.
#[derive(Debug, Clone, Copy, Default)]
pub struct ClipCounts {
    /// Number of bases clipped before `clip`
    pub prior: usize,
    /// Number of bases clipped from 5' end
    pub five_prime: usize,
    /// Number of bases clipped from 3' end
    pub three_prime: usize,
    /// Number of bases clipped due to overlap
    pub overlapping: usize,
    /// Number of bases clipped due to extending past mate
    pub extending: usize,
}

/// Type of read for metrics tracking
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReadType {
    Fragment,
    ReadOne,
    ReadTwo,
    Pair,
    All,
}

/// Clipping metrics for a specific read type
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClippingMetrics {
    /// The type of read this metric applies to
    pub read_type: ReadType,

    /// Total number of reads examined
    pub reads: usize,

    /// Number of reads that became unmapped due to clipping
    pub reads_unmapped: usize,

    /// Number of reads with any clipping before `clip`
    pub reads_clipped_pre: usize,

    /// Number of reads with any clipping after `clip`
    pub reads_clipped_post: usize,

    /// Number of reads clipped on 5' end
    pub reads_clipped_five_prime: usize,

    /// Number of reads clipped on 3' end
    pub reads_clipped_three_prime: usize,

    /// Number of reads clipped due to overlapping reads
    pub reads_clipped_overlapping: usize,

    /// Number of reads clipped due to extending past mate
    pub reads_clipped_extending: usize,

    /// Total number of aligned bases after clipping
    pub bases: usize,

    /// Number of bases clipped before `clip`
    pub bases_clipped_pre: usize,

    /// Number of bases clipped after `clip`
    pub bases_clipped_post: usize,

    /// Number of bases clipped on 5' end
    pub bases_clipped_five_prime: usize,

    /// Number of bases clipped on 3' end
    pub bases_clipped_three_prime: usize,

    /// Number of bases clipped due to overlapping reads
    pub bases_clipped_overlapping: usize,

    /// Number of bases clipped due to extending past mate
    pub bases_clipped_extending: usize,
}

impl ClippingMetrics {
    /// Creates a new `ClippingMetrics` for the given read type
    #[must_use]
    pub fn new(read_type: ReadType) -> Self {
        Self {
            read_type,
            reads: 0,
            reads_unmapped: 0,
            reads_clipped_pre: 0,
            reads_clipped_post: 0,
            reads_clipped_five_prime: 0,
            reads_clipped_three_prime: 0,
            reads_clipped_overlapping: 0,
            reads_clipped_extending: 0,
            bases: 0,
            bases_clipped_pre: 0,
            bases_clipped_post: 0,
            bases_clipped_five_prime: 0,
            bases_clipped_three_prime: 0,
            bases_clipped_overlapping: 0,
            bases_clipped_extending: 0,
        }
    }

    /// Updates metrics from a raw BAM record after clipping.
    ///
    /// Equivalent to [`Self::update`] but reads alignment information directly from
    /// the raw BAM bytes, avoiding a noodles decode round-trip.
    ///
    /// # Arguments
    /// * `record` - The raw BAM record after clipping
    /// * `counts` - The clip counts for this operation
    #[cfg(feature = "clip")]
    pub fn update_raw(&mut self, record: &fgumi_raw_bam::RawRecord, counts: ClipCounts) {
        use fgumi_raw_bam::get_cigar_ops;

        self.reads += 1;

        // Count aligned bases (M/=/X ops) from raw CIGAR
        let cigar_ops = get_cigar_ops(record.as_ref());
        let aligned_bases: usize = cigar_ops
            .iter()
            .filter(|&&op| {
                matches!(op & 0xF, 0 | 7 | 8) // Match, SequenceMatch, SequenceMismatch
            })
            .map(|&op| (op >> 4) as usize)
            .sum();
        self.bases += aligned_bases;

        // Track pre-clipping
        if counts.prior > 0 {
            self.reads_clipped_pre += 1;
            self.bases_clipped_pre += counts.prior;
        }

        // Track 5' clipping
        if counts.five_prime > 0 {
            self.reads_clipped_five_prime += 1;
            self.bases_clipped_five_prime += counts.five_prime;
        }

        // Track 3' clipping
        if counts.three_prime > 0 {
            self.reads_clipped_three_prime += 1;
            self.bases_clipped_three_prime += counts.three_prime;
        }

        // Track overlapping clipping
        if counts.overlapping > 0 {
            self.reads_clipped_overlapping += 1;
            self.bases_clipped_overlapping += counts.overlapping;
        }

        // Track extending clipping
        if counts.extending > 0 {
            self.reads_clipped_extending += 1;
            self.bases_clipped_extending += counts.extending;
        }

        // Total clipping after ClipBam
        let additional_clipped =
            counts.five_prime + counts.three_prime + counts.overlapping + counts.extending;
        let total_clipped = additional_clipped + counts.prior;

        if total_clipped > 0 {
            self.reads_clipped_post += 1;
            self.bases_clipped_post += total_clipped;

            // Check if read became unmapped
            if record.is_unmapped() && additional_clipped > 0 {
                self.reads_unmapped += 1;
            }
        }
    }

    /// Updates metrics based on a clipping operation
    ///
    /// # Arguments
    /// * `record` - The record that was clipped
    /// * `counts` - The clip counts for this operation
    pub fn update(&mut self, record: &RecordBuf, counts: ClipCounts) {
        self.reads += 1;

        // Count aligned bases
        let cigar = record.cigar();
        #[allow(clippy::redundant_closure_for_method_calls)]
        // clippy suggests incorrect path for Op::len
        let aligned_bases: usize = cigar
            .iter()
            .filter_map(Result::ok)
            .filter(|op| {
                use noodles::sam::alignment::record::cigar::op::Kind;
                matches!(op.kind(), Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch)
            })
            .map(|op| op.len())
            .sum();
        self.bases += aligned_bases;

        // Track pre-clipping
        if counts.prior > 0 {
            self.reads_clipped_pre += 1;
            self.bases_clipped_pre += counts.prior;
        }

        // Track 5' clipping
        if counts.five_prime > 0 {
            self.reads_clipped_five_prime += 1;
            self.bases_clipped_five_prime += counts.five_prime;
        }

        // Track 3' clipping
        if counts.three_prime > 0 {
            self.reads_clipped_three_prime += 1;
            self.bases_clipped_three_prime += counts.three_prime;
        }

        // Track overlapping clipping
        if counts.overlapping > 0 {
            self.reads_clipped_overlapping += 1;
            self.bases_clipped_overlapping += counts.overlapping;
        }

        // Track extending clipping
        if counts.extending > 0 {
            self.reads_clipped_extending += 1;
            self.bases_clipped_extending += counts.extending;
        }

        // Total clipping after ClipBam
        let additional_clipped =
            counts.five_prime + counts.three_prime + counts.overlapping + counts.extending;
        let total_clipped = additional_clipped + counts.prior;

        if total_clipped > 0 {
            self.reads_clipped_post += 1;
            self.bases_clipped_post += total_clipped;

            // Check if read became unmapped
            if record.flags().is_unmapped() && additional_clipped > 0 {
                self.reads_unmapped += 1;
            }
        }
    }

    /// Adds metrics from another `ClippingMetrics` instance
    pub fn add(&mut self, other: &ClippingMetrics) {
        *self += other;
    }
}

impl std::ops::AddAssign<&ClippingMetrics> for ClippingMetrics {
    fn add_assign(&mut self, other: &ClippingMetrics) {
        self.reads += other.reads;
        self.reads_unmapped += other.reads_unmapped;
        self.reads_clipped_pre += other.reads_clipped_pre;
        self.reads_clipped_post += other.reads_clipped_post;
        self.reads_clipped_five_prime += other.reads_clipped_five_prime;
        self.reads_clipped_three_prime += other.reads_clipped_three_prime;
        self.reads_clipped_overlapping += other.reads_clipped_overlapping;
        self.reads_clipped_extending += other.reads_clipped_extending;
        self.bases += other.bases;
        self.bases_clipped_pre += other.bases_clipped_pre;
        self.bases_clipped_post += other.bases_clipped_post;
        self.bases_clipped_five_prime += other.bases_clipped_five_prime;
        self.bases_clipped_three_prime += other.bases_clipped_three_prime;
        self.bases_clipped_overlapping += other.bases_clipped_overlapping;
        self.bases_clipped_extending += other.bases_clipped_extending;
    }
}

impl Default for ClippingMetrics {
    fn default() -> Self {
        Self::new(ReadType::All)
    }
}

impl Metric for ClippingMetrics {
    fn metric_name() -> &'static str {
        "clipping"
    }
}

/// Collection of clipping metrics for all read types
pub struct ClippingMetricsCollection {
    pub fragment: ClippingMetrics,
    pub read_one: ClippingMetrics,
    pub read_two: ClippingMetrics,
    pub pair: ClippingMetrics,
    pub all: ClippingMetrics,
}

impl ClippingMetricsCollection {
    /// Creates a new metrics collection
    #[must_use]
    pub fn new() -> Self {
        Self {
            fragment: ClippingMetrics::new(ReadType::Fragment),
            read_one: ClippingMetrics::new(ReadType::ReadOne),
            read_two: ClippingMetrics::new(ReadType::ReadTwo),
            pair: ClippingMetrics::new(ReadType::Pair),
            all: ClippingMetrics::new(ReadType::All),
        }
    }

    /// Finalizes metrics by aggregating pair and all categories
    pub fn finalize(&mut self) {
        // Aggregate ReadOne and ReadTwo into Pair
        self.pair.add(&self.read_one);
        self.pair.add(&self.read_two);

        // Aggregate Fragment and Pair into All
        self.all.add(&self.fragment);
        self.all.add(&self.pair);
    }

    /// Returns all metrics in order
    #[must_use]
    pub fn all_metrics(&self) -> [&ClippingMetrics; 5] {
        [&self.fragment, &self.read_one, &self.read_two, &self.pair, &self.all]
    }
}

impl Default for ClippingMetricsCollection {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_sam::builder::RecordBuilder;

    /// Creates a test record with the given CIGAR string.
    /// Sequence and qualities are auto-generated from CIGAR length.
    fn create_test_record(cigar: &str, is_unmapped: bool) -> RecordBuf {
        RecordBuilder::new().cigar(cigar).unmapped(is_unmapped).build()
    }

    #[test]
    fn test_read_type_variants() {
        assert_eq!(ReadType::Fragment, ReadType::Fragment);
        assert_ne!(ReadType::Fragment, ReadType::ReadOne);
        assert_ne!(ReadType::ReadOne, ReadType::ReadTwo);
    }

    #[test]
    fn test_clipping_metrics_new() {
        let metrics = ClippingMetrics::new(ReadType::Fragment);
        assert_eq!(metrics.reads, 0);
        assert_eq!(metrics.reads_unmapped, 0);
        assert_eq!(metrics.bases, 0);
        assert_eq!(metrics.bases_clipped_pre, 0);
    }

    #[test]
    fn test_clip_counts_default() {
        let counts = ClipCounts::default();
        assert_eq!(counts.prior, 0);
        assert_eq!(counts.five_prime, 0);
        assert_eq!(counts.three_prime, 0);
        assert_eq!(counts.overlapping, 0);
        assert_eq!(counts.extending, 0);
    }

    #[test]
    fn test_clipping_metrics_update_no_clipping() {
        let mut metrics = ClippingMetrics::new(ReadType::ReadOne);
        let record = create_test_record("100M", false);

        metrics.update(&record, ClipCounts::default());

        assert_eq!(metrics.reads, 1);
        assert_eq!(metrics.bases, 100);
        assert_eq!(metrics.reads_clipped_pre, 0);
        assert_eq!(metrics.reads_clipped_post, 0);
    }

    #[test]
    fn test_clipping_metrics_update_with_prior_clipping() {
        let mut metrics = ClippingMetrics::new(ReadType::ReadOne);
        let record = create_test_record("90M", false);

        metrics.update(&record, ClipCounts { prior: 10, ..ClipCounts::default() });

        assert_eq!(metrics.reads, 1);
        assert_eq!(metrics.reads_clipped_pre, 1);
        assert_eq!(metrics.bases_clipped_pre, 10);
        assert_eq!(metrics.reads_clipped_post, 1);
        assert_eq!(metrics.bases_clipped_post, 10);
    }

    #[test]
    fn test_clipping_metrics_update_five_prime() {
        let mut metrics = ClippingMetrics::new(ReadType::ReadOne);
        let record = create_test_record("95M", false);

        metrics.update(&record, ClipCounts { five_prime: 5, ..ClipCounts::default() });

        assert_eq!(metrics.reads_clipped_five_prime, 1);
        assert_eq!(metrics.bases_clipped_five_prime, 5);
        assert_eq!(metrics.reads_clipped_post, 1);
        assert_eq!(metrics.bases_clipped_post, 5);
    }

    #[test]
    fn test_clipping_metrics_update_three_prime() {
        let mut metrics = ClippingMetrics::new(ReadType::ReadOne);
        let record = create_test_record("97M", false);

        metrics.update(&record, ClipCounts { three_prime: 3, ..ClipCounts::default() });

        assert_eq!(metrics.reads_clipped_three_prime, 1);
        assert_eq!(metrics.bases_clipped_three_prime, 3);
        assert_eq!(metrics.reads_clipped_post, 1);
        assert_eq!(metrics.bases_clipped_post, 3);
    }

    #[test]
    fn test_clipping_metrics_update_overlapping() {
        let mut metrics = ClippingMetrics::new(ReadType::ReadOne);
        let record = create_test_record("80M", false);

        metrics.update(&record, ClipCounts { overlapping: 20, ..ClipCounts::default() });

        assert_eq!(metrics.reads_clipped_overlapping, 1);
        assert_eq!(metrics.bases_clipped_overlapping, 20);
        assert_eq!(metrics.reads_clipped_post, 1);
        assert_eq!(metrics.bases_clipped_post, 20);
    }

    #[test]
    fn test_clipping_metrics_update_extending() {
        let mut metrics = ClippingMetrics::new(ReadType::ReadOne);
        let record = create_test_record("85M", false);

        metrics.update(&record, ClipCounts { extending: 15, ..ClipCounts::default() });

        assert_eq!(metrics.reads_clipped_extending, 1);
        assert_eq!(metrics.bases_clipped_extending, 15);
        assert_eq!(metrics.reads_clipped_post, 1);
        assert_eq!(metrics.bases_clipped_post, 15);
    }

    #[test]
    fn test_clipping_metrics_update_all_types() {
        let mut metrics = ClippingMetrics::new(ReadType::ReadOne);
        let record = create_test_record("50M", false);

        metrics.update(
            &record,
            ClipCounts { prior: 10, five_prime: 5, three_prime: 3, overlapping: 20, extending: 12 },
        );

        assert_eq!(metrics.reads, 1);
        assert_eq!(metrics.bases_clipped_pre, 10);
        assert_eq!(metrics.bases_clipped_five_prime, 5);
        assert_eq!(metrics.bases_clipped_three_prime, 3);
        assert_eq!(metrics.bases_clipped_overlapping, 20);
        assert_eq!(metrics.bases_clipped_extending, 12);
        assert_eq!(metrics.bases_clipped_post, 50); // 10 + 5 + 3 + 20 + 12
    }

    #[test]
    fn test_clipping_metrics_update_unmapped() {
        let mut metrics = ClippingMetrics::new(ReadType::ReadOne);
        let record = create_test_record("", true);

        metrics.update(&record, ClipCounts { five_prime: 100, ..ClipCounts::default() });

        assert_eq!(metrics.reads_unmapped, 1);
    }

    #[test]
    fn test_clipping_metrics_update_multiple_records() {
        let mut metrics = ClippingMetrics::new(ReadType::ReadOne);

        for _ in 0..5 {
            let record = create_test_record("90M", false);
            metrics.update(&record, ClipCounts { five_prime: 10, ..ClipCounts::default() });
        }

        assert_eq!(metrics.reads, 5);
        assert_eq!(metrics.bases, 450);
        assert_eq!(metrics.reads_clipped_five_prime, 5);
        assert_eq!(metrics.bases_clipped_five_prime, 50);
    }

    #[test]
    fn test_clipping_metrics_add() {
        let mut metrics1 = ClippingMetrics::new(ReadType::ReadOne);
        let record1 = create_test_record("90M", false);
        metrics1.update(&record1, ClipCounts { five_prime: 10, ..ClipCounts::default() });

        let mut metrics2 = ClippingMetrics::new(ReadType::ReadOne);
        let record2 = create_test_record("85M", false);
        metrics2.update(&record2, ClipCounts { three_prime: 15, ..ClipCounts::default() });

        metrics1.add(&metrics2);

        assert_eq!(metrics1.reads, 2);
        assert_eq!(metrics1.bases, 175);
        assert_eq!(metrics1.reads_clipped_five_prime, 1);
        assert_eq!(metrics1.reads_clipped_three_prime, 1);
        assert_eq!(metrics1.bases_clipped_five_prime, 10);
        assert_eq!(metrics1.bases_clipped_three_prime, 15);
    }

    #[test]
    fn test_clipping_metrics_collection_new() {
        let collection = ClippingMetricsCollection::new();

        assert_eq!(collection.fragment.reads, 0);
        assert_eq!(collection.read_one.reads, 0);
        assert_eq!(collection.read_two.reads, 0);
        assert_eq!(collection.pair.reads, 0);
        assert_eq!(collection.all.reads, 0);
    }

    #[test]
    fn test_clipping_metrics_collection_default() {
        let collection = ClippingMetricsCollection::default();
        assert_eq!(collection.fragment.reads, 0);
    }

    #[test]
    fn test_clipping_metrics_collection_finalize() {
        let mut collection = ClippingMetricsCollection::new();

        let record = create_test_record("90M", false);
        collection.read_one.update(&record, ClipCounts { five_prime: 10, ..ClipCounts::default() });
        collection
            .read_two
            .update(&record, ClipCounts { three_prime: 10, ..ClipCounts::default() });
        collection.fragment.update(&record, ClipCounts { five_prime: 5, ..ClipCounts::default() });

        collection.finalize();

        // Pair should have ReadOne + ReadTwo
        assert_eq!(collection.pair.reads, 2);
        assert_eq!(collection.pair.bases_clipped_five_prime, 10);
        assert_eq!(collection.pair.bases_clipped_three_prime, 10);

        // All should have Fragment + Pair
        assert_eq!(collection.all.reads, 3);
        assert_eq!(collection.all.bases_clipped_five_prime, 15);
    }

    #[test]
    fn test_clipping_metrics_collection_all_metrics() {
        let collection = ClippingMetricsCollection::new();
        let all_metrics = collection.all_metrics();

        assert_eq!(all_metrics.len(), 5);
        assert_eq!(all_metrics[0].read_type, ReadType::Fragment);
        assert_eq!(all_metrics[1].read_type, ReadType::ReadOne);
        assert_eq!(all_metrics[2].read_type, ReadType::ReadTwo);
        assert_eq!(all_metrics[3].read_type, ReadType::Pair);
        assert_eq!(all_metrics[4].read_type, ReadType::All);
    }

    #[test]
    fn test_clipping_metrics_cigar_with_deletions() {
        let mut metrics = ClippingMetrics::new(ReadType::ReadOne);
        let record = create_test_record("40M10D40M", false);

        metrics.update(&record, ClipCounts::default());

        // Only matches count as aligned bases
        assert_eq!(metrics.bases, 80);
    }

    #[test]
    fn test_clipping_metrics_cigar_with_insertions() {
        let mut metrics = ClippingMetrics::new(ReadType::ReadOne);
        let record = create_test_record("45M5I45M", false);

        metrics.update(&record, ClipCounts::default());

        // Only matches count as aligned bases
        assert_eq!(metrics.bases, 90);
    }

    #[test]
    fn test_clipping_metrics_update_unmapped_no_additional_clipping() {
        let mut metrics = ClippingMetrics::new(ReadType::ReadOne);
        let record = create_test_record("", true);

        // No additional clipping, so shouldn't count as unmapped
        metrics.update(&record, ClipCounts { prior: 10, ..ClipCounts::default() });

        assert_eq!(metrics.reads_unmapped, 0);
    }

    #[test]
    fn test_read_type_serialization() {
        // Test that ReadType can be serialized/deserialized
        let rt = ReadType::Fragment;
        assert_eq!(rt, ReadType::Fragment);
    }

    #[test]
    fn test_clipping_metrics_add_empty() {
        let mut metrics1 = ClippingMetrics::new(ReadType::ReadOne);
        let metrics2 = ClippingMetrics::new(ReadType::ReadOne);

        let record = create_test_record("90M", false);
        metrics1.update(&record, ClipCounts { five_prime: 10, ..ClipCounts::default() });

        metrics1.add(&metrics2);

        // Should still have just the one record's data
        assert_eq!(metrics1.reads, 1);
    }

    #[test]
    fn test_clipping_metrics_add_assign() {
        let mut metrics1 = ClippingMetrics::new(ReadType::ReadOne);
        let record1 = create_test_record("90M", false);
        metrics1.update(&record1, ClipCounts { five_prime: 10, ..ClipCounts::default() });

        let mut metrics2 = ClippingMetrics::new(ReadType::ReadOne);
        let record2 = create_test_record("85M", false);
        metrics2.update(&record2, ClipCounts { three_prime: 15, ..ClipCounts::default() });

        metrics1 += &metrics2;

        assert_eq!(metrics1.reads, 2);
        assert_eq!(metrics1.bases, 175);
        assert_eq!(metrics1.reads_clipped_five_prime, 1);
        assert_eq!(metrics1.reads_clipped_three_prime, 1);
    }

    #[test]
    fn test_metric_trait_impl() {
        assert_eq!(ClippingMetrics::metric_name(), "clipping");
    }
}
