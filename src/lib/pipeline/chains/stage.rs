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

    #[test]
    fn consensus_classifier_is_exhaustive() {
        for s in [Stage::Simplex, Stage::Duplex, Stage::Codec] {
            assert!(s.is_consensus(), "{s:?} should be consensus");
        }
        for s in [Stage::Correct, Stage::Group, Stage::Sort, Stage::Clip] {
            assert!(!s.is_consensus(), "{s:?} should not be consensus");
        }
    }

    #[test]
    fn template_producer_classifier_is_exhaustive() {
        assert!(Stage::Align.is_template_producer());
        assert!(Stage::Zipper.is_template_producer());
        for s in [Stage::Correct, Stage::Group, Stage::Simplex] {
            assert!(!s.is_template_producer());
        }
    }
}
