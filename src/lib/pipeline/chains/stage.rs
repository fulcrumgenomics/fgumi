//! Stage labels for the unified chain-builder.

/// One pipeline stage. Each variant is a label, not data — the canonical
/// multi-step sequence each variant expands to lives in
/// [`chains::build_for`](crate::pipeline::chains::build_for)'s
/// match arms.
///
/// **Mutual exclusions** (enforced by
/// [`crate::pipeline::chains::validate::validate_stage_progression`]):
///
/// - [`Stage::Align`] and [`Stage::Zipper`] both produce
///   `BamTemplateBatch` from different sources — at most one per chain.
/// - [`Stage::Simplex`], [`Stage::Duplex`], [`Stage::Codec`] are
///   terminal consensus stages — at most one per chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Stage {
    Extract,
    Correct,
    Align,
    Zipper,
    Sort,
    Group,
    Simplex,
    Duplex,
    Codec,
    Clip,
    Filter,
    Dedup,
    Downsample,
    /// Terminal BAM → FASTQ encode (interleaved or paired split output).
    Fastq,
}

impl Stage {
    /// Returns `true` if this stage is a terminal consensus stage.
    /// Used by validators that reject post-consensus stages.
    #[must_use]
    pub fn is_consensus(self) -> bool {
        matches!(self, Stage::Simplex | Stage::Duplex | Stage::Codec)
    }

    /// Returns `true` if this stage produces `BamTemplateBatch` from
    /// an external source (AAM subprocess or zipper merge).
    #[must_use]
    pub fn is_template_producer(self) -> bool {
        matches!(self, Stage::Align | Stage::Zipper)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every `Stage` variant, so the classifier tests below cover the whole
    /// enum rather than a hand-picked subset. A new variant that is not added
    /// here is still forced through the exhaustive `match`es in the
    /// `expected_*` helpers, which fail to compile until it is classified.
    const ALL_STAGES: [Stage; 14] = [
        Stage::Extract,
        Stage::Correct,
        Stage::Align,
        Stage::Zipper,
        Stage::Sort,
        Stage::Group,
        Stage::Simplex,
        Stage::Duplex,
        Stage::Codec,
        Stage::Clip,
        Stage::Filter,
        Stage::Dedup,
        Stage::Downsample,
        Stage::Fastq,
    ];

    /// Expected `is_consensus` truth per variant, spelled out with an
    /// exhaustive `match` so adding a `Stage` variant is a compile error here
    /// until it is deliberately classified.
    fn expected_is_consensus(stage: Stage) -> bool {
        match stage {
            Stage::Simplex | Stage::Duplex | Stage::Codec => true,
            Stage::Extract
            | Stage::Correct
            | Stage::Align
            | Stage::Zipper
            | Stage::Sort
            | Stage::Group
            | Stage::Clip
            | Stage::Filter
            | Stage::Dedup
            | Stage::Downsample
            | Stage::Fastq => false,
        }
    }

    /// Expected `is_template_producer` truth per variant, exhaustive for the
    /// same reason as [`expected_is_consensus`].
    fn expected_is_template_producer(stage: Stage) -> bool {
        match stage {
            Stage::Align | Stage::Zipper => true,
            Stage::Extract
            | Stage::Correct
            | Stage::Sort
            | Stage::Group
            | Stage::Simplex
            | Stage::Duplex
            | Stage::Codec
            | Stage::Clip
            | Stage::Filter
            | Stage::Dedup
            | Stage::Downsample
            | Stage::Fastq => false,
        }
    }

    #[test]
    fn consensus_classifier_matches_every_variant() {
        for s in ALL_STAGES {
            assert_eq!(s.is_consensus(), expected_is_consensus(s), "is_consensus({s:?})");
        }
    }

    #[test]
    fn template_producer_classifier_matches_every_variant() {
        for s in ALL_STAGES {
            assert_eq!(
                s.is_template_producer(),
                expected_is_template_producer(s),
                "is_template_producer({s:?})"
            );
        }
    }
}
