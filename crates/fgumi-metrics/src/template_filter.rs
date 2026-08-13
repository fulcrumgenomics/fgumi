//! Template-level filter accounting shared by the `group` and `dedup` commands.
//!
//! Both commands run the same pre-grouping template filter but report it
//! differently: `group` emits fgbio's five-column `UmiGroupingMetric` in
//! **primary-record** units, `dedup` emits fgumi-native `filtered_*` columns in
//! **template** units. [`TemplateFilterCounts`] is the single measurement those
//! two views render from, so it records every decision in both units and keys
//! rejections by [`TemplateFilterReason`] — one variant per rejection site in the
//! filter, independent of any output file's column vocabulary.

/// Why the pre-grouping template filter rejected a template.
///
/// One variant per rejection site, not per output column: fgbio's
/// `UmiGroupingMetric` collapses six of these into a single
/// `discarded_poor_alignment` column, which is a property of that file format
/// rather than of the filter. Adding a rejection site requires adding a variant
/// here, which breaks the exhaustive matches in every renderer until each decides
/// what to do with it.
///
/// `#[repr(usize)]` is load-bearing: [`TemplateFilterCounts`] indexes its arrays
/// by `reason as usize`, and [`Self::ALL`] is ordered to match.
// `EnumIter` (test-only) auto-generates `iter()` so coverage assertions enumerate
// every variant without a hand-maintained list that can drift.
#[cfg_attr(test, derive(strum::EnumIter))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(usize)]
pub enum TemplateFilterReason {
    /// A record was shorter than the minimum BAM record length (corrupt input).
    MalformedRecord,
    /// The template had neither an R1 nor an R2 primary record.
    NoPrimaryReads,
    /// Every primary read in the template was unmapped.
    Unmapped,
    /// A primary read carried the QC-fail (`0x200`) flag.
    NotPassingFilter,
    /// A mapped primary read's mapping quality was below the threshold.
    LowMappingQuality,
    /// A primary read's mate mapping quality (`MQ` tag) was below the threshold.
    LowMateMappingQuality,
    /// A primary read had no UMI tag.
    MissingUmi,
    /// The UMI contained at least one `N` base.
    NsInUmi,
    /// The UMI was shorter than the configured minimum length.
    UmiTooShort,
}

impl TemplateFilterReason {
    /// Number of variants.
    pub const COUNT: usize = 9;

    /// Every variant, in discriminant order.
    pub const ALL: [Self; Self::COUNT] = [
        Self::MalformedRecord,
        Self::NoPrimaryReads,
        Self::Unmapped,
        Self::NotPassingFilter,
        Self::LowMappingQuality,
        Self::LowMateMappingQuality,
        Self::MissingUmi,
        Self::NsInUmi,
        Self::UmiTooShort,
    ];

    /// Snake-case identifier for this reason, used as the suffix of the
    /// `filtered_*` column in `dedup`'s metrics file and in log summaries.
    #[must_use]
    pub fn column_suffix(self) -> &'static str {
        match self {
            Self::MalformedRecord => "malformed_record",
            Self::NoPrimaryReads => "no_primary_reads",
            Self::Unmapped => "unmapped",
            Self::NotPassingFilter => "not_passing_filter",
            Self::LowMappingQuality => "low_mapping_quality",
            Self::LowMateMappingQuality => "low_mate_mapping_quality",
            Self::MissingUmi => "missing_umi",
            Self::NsInUmi => "ns_in_umi",
            Self::UmiTooShort => "umi_too_short",
        }
    }
}

/// Per-reason accounting for the pre-grouping template filter.
///
/// The filter is the only writer: each template that enters is recorded exactly
/// once, via [`Self::record_accepted`] or [`Self::record_rejected`]. That makes
/// `total == accepted + sum(rejected)` true by construction in both units.
///
/// Counts are tracked in two units simultaneously:
/// - **templates**: one per template, regardless of how many records it holds
/// - **primary reads**: the number of primary (non-secondary, non-supplementary)
///   records in the template, which is fgbio's unit for `UmiGroupingMetric`
///
/// A template with neither an R1 nor an R2 contributes `1` template and `0`
/// primary reads, so the units are not proportional and neither is derivable
/// from the other.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct TemplateFilterCounts {
    accepted_templates: u64,
    accepted_primary_reads: u64,
    rejected_templates: [u64; TemplateFilterReason::COUNT],
    rejected_primary_reads: [u64; TemplateFilterReason::COUNT],
}

impl TemplateFilterCounts {
    /// Creates an empty set of counts.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Records a template that passed the filter, carrying `primary_reads`
    /// primary records.
    pub fn record_accepted(&mut self, primary_reads: u64) {
        self.accepted_templates += 1;
        self.accepted_primary_reads += primary_reads;
    }

    /// Records a template rejected for `reason`, carrying `primary_reads` primary
    /// records.
    pub fn record_rejected(&mut self, reason: TemplateFilterReason, primary_reads: u64) {
        self.rejected_templates[reason as usize] += 1;
        self.rejected_primary_reads[reason as usize] += primary_reads;
    }

    /// Adds every count in `other` into `self`.
    pub fn merge(&mut self, other: &Self) {
        self.accepted_templates += other.accepted_templates;
        self.accepted_primary_reads += other.accepted_primary_reads;
        for index in 0..TemplateFilterReason::COUNT {
            self.rejected_templates[index] += other.rejected_templates[index];
            self.rejected_primary_reads[index] += other.rejected_primary_reads[index];
        }
    }

    /// Templates that passed the filter.
    #[must_use]
    pub fn accepted_templates(&self) -> u64 {
        self.accepted_templates
    }

    /// Primary reads in templates that passed the filter.
    #[must_use]
    pub fn accepted_primary_reads(&self) -> u64 {
        self.accepted_primary_reads
    }

    /// Templates rejected for `reason`.
    #[must_use]
    pub fn rejected_templates(&self, reason: TemplateFilterReason) -> u64 {
        self.rejected_templates[reason as usize]
    }

    /// Primary reads in templates rejected for `reason`.
    #[must_use]
    pub fn rejected_primary_reads(&self, reason: TemplateFilterReason) -> u64 {
        self.rejected_primary_reads[reason as usize]
    }

    /// Templates rejected for any reason.
    #[must_use]
    pub fn total_rejected_templates(&self) -> u64 {
        self.rejected_templates.iter().sum()
    }

    /// Primary reads in templates rejected for any reason.
    #[must_use]
    pub fn total_rejected_primary_reads(&self) -> u64 {
        self.rejected_primary_reads.iter().sum()
    }

    /// Templates that entered the filter (accepted plus rejected).
    #[must_use]
    pub fn total_templates(&self) -> u64 {
        self.accepted_templates + self.total_rejected_templates()
    }

    /// Primary reads that entered the filter (accepted plus rejected).
    #[must_use]
    pub fn total_primary_reads(&self) -> u64 {
        self.accepted_primary_reads + self.total_rejected_primary_reads()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use strum::IntoEnumIterator;

    /// `TemplateFilterCounts` indexes its arrays by `reason as usize`, so `ALL`
    /// must list every variant in discriminant order or counts land in the wrong
    /// bucket.
    #[test]
    fn all_covers_every_variant_in_discriminant_order() {
        assert_eq!(TemplateFilterReason::ALL.len(), TemplateFilterReason::COUNT);
        assert_eq!(TemplateFilterReason::iter().count(), TemplateFilterReason::COUNT);
        for (position, &reason) in TemplateFilterReason::ALL.iter().enumerate() {
            assert_eq!(reason as usize, position, "{reason:?} is out of order in ALL");
        }
    }

    /// Column suffixes become metric column names: unique and `snake_case`.
    #[test]
    fn column_suffixes_are_unique_snake_case() {
        let mut seen = std::collections::HashSet::new();
        for reason in TemplateFilterReason::ALL {
            let suffix = reason.column_suffix();
            assert!(seen.insert(suffix), "duplicate column suffix {suffix}");
            assert!(
                suffix.bytes().all(|b| b.is_ascii_lowercase() || b == b'_'),
                "{suffix} is not snake_case"
            );
        }
    }

    /// The central invariant: every template entering the filter is recorded
    /// exactly once, so totals reconcile in both units.
    #[rstest::rstest]
    #[case::only_accepts(&[], 3, 2)]
    #[case::only_rejects(&[(TemplateFilterReason::NsInUmi, 2)], 0, 0)]
    #[case::mixed(&[(TemplateFilterReason::NotPassingFilter, 2), (TemplateFilterReason::MalformedRecord, 1)], 4, 2)]
    #[case::zero_primary_reads(&[(TemplateFilterReason::NoPrimaryReads, 0)], 1, 2)]
    fn totals_reconcile(
        #[case] rejects: &[(TemplateFilterReason, u64)],
        #[case] accepts: u64,
        #[case] primary_reads_per_accept: u64,
    ) {
        let mut counts = TemplateFilterCounts::new();
        for _ in 0..accepts {
            counts.record_accepted(primary_reads_per_accept);
        }
        for &(reason, primary_reads) in rejects {
            counts.record_rejected(reason, primary_reads);
        }

        assert_eq!(
            counts.total_templates(),
            counts.accepted_templates() + counts.total_rejected_templates()
        );
        assert_eq!(
            counts.total_primary_reads(),
            counts.accepted_primary_reads() + counts.total_rejected_primary_reads()
        );
        assert_eq!(counts.accepted_templates(), accepts);
        assert_eq!(counts.accepted_primary_reads(), accepts * primary_reads_per_accept);
        assert_eq!(counts.total_rejected_templates(), rejects.len() as u64);
        assert_eq!(
            counts.total_rejected_primary_reads(),
            rejects.iter().map(|&(_, n)| n).sum::<u64>()
        );
    }

    /// A rejection lands in its own bucket and nowhere else.
    #[rstest::rstest]
    #[case::malformed(TemplateFilterReason::MalformedRecord)]
    #[case::no_primary(TemplateFilterReason::NoPrimaryReads)]
    #[case::unmapped(TemplateFilterReason::Unmapped)]
    #[case::non_pf(TemplateFilterReason::NotPassingFilter)]
    #[case::low_mapq(TemplateFilterReason::LowMappingQuality)]
    #[case::low_mate_mapq(TemplateFilterReason::LowMateMappingQuality)]
    #[case::missing_umi(TemplateFilterReason::MissingUmi)]
    #[case::ns_in_umi(TemplateFilterReason::NsInUmi)]
    #[case::umi_too_short(TemplateFilterReason::UmiTooShort)]
    fn rejection_lands_in_exactly_one_bucket(#[case] recorded: TemplateFilterReason) {
        let mut counts = TemplateFilterCounts::new();
        counts.record_rejected(recorded, 2);

        for reason in TemplateFilterReason::ALL {
            let expected_templates = u64::from(reason == recorded);
            let expected_reads = if reason == recorded { 2 } else { 0 };
            assert_eq!(counts.rejected_templates(reason), expected_templates, "{reason:?}");
            assert_eq!(counts.rejected_primary_reads(reason), expected_reads, "{reason:?}");
        }
    }

    /// Merge aggregates per-worker counters; it must be additive on every field.
    #[test]
    fn merge_is_additive() {
        let mut left = TemplateFilterCounts::new();
        left.record_accepted(2);
        left.record_rejected(TemplateFilterReason::NsInUmi, 2);

        let mut right = TemplateFilterCounts::new();
        right.record_accepted(1);
        right.record_rejected(TemplateFilterReason::NsInUmi, 2);
        right.record_rejected(TemplateFilterReason::Unmapped, 2);

        left.merge(&right);

        assert_eq!(left.accepted_templates(), 2);
        assert_eq!(left.accepted_primary_reads(), 3);
        assert_eq!(left.rejected_templates(TemplateFilterReason::NsInUmi), 2);
        assert_eq!(left.rejected_primary_reads(TemplateFilterReason::NsInUmi), 4);
        assert_eq!(left.rejected_templates(TemplateFilterReason::Unmapped), 1);
        assert_eq!(left.total_templates(), 5);
        assert_eq!(left.total_primary_reads(), 9);
    }
}
