//! What feeds the chain.

use std::path::PathBuf;

use anyhow::{Result, ensure};
use read_structure::ReadStructure;

/// The source feeding a chain. Variants encode the source shape; the
/// chain builder resolves each variant to the appropriate source step(s).
#[derive(Debug, Clone)]
pub enum SourceSpec {
    /// Single BAM file (or `-` for stdin). Used by standalone correct,
    /// sort, group, consensus, clip, filter, dedup, downsample, and by
    /// runall when `--start-from` lands at any of those stages.
    Bam(PathBuf),

    /// Paired BAMs joined by `ZipperMergeStep`. Used by standalone
    /// zipper and by `runall --start-from zipper`.
    PairedBams { unmapped: PathBuf, mapped: PathBuf, reference: PathBuf },

    /// FASTQ inputs feeding `fgumi extract`. Used by standalone extract
    /// and by `runall --start-from extract` (the latter unlocks in
    /// Phase 5).
    ///
    /// `paths` and `read_structures` are positional pairs and **must** be
    /// the same length: `read_structures[i]` describes `paths[i]`. Build
    /// this variant via [`SourceSpec::fastqs`] (which enforces the
    /// invariant) rather than a struct literal; [`SourceSpec::validate`]
    /// re-checks it for any spec that was constructed directly.
    Fastqs { paths: Vec<PathBuf>, read_structures: Vec<ReadStructure> },

    /// SAM text file (or `-` for stdin). Standalone-only — runall does
    /// not accept SAM input (it always reads from a BAM intermediate).
    Sam(PathBuf),
}

impl SourceSpec {
    /// Constructs a [`SourceSpec::Fastqs`] source, enforcing that `paths`
    /// and `read_structures` are positional pairs of equal length.
    ///
    /// Prefer this over a `SourceSpec::Fastqs { .. }` struct literal so a
    /// mis-sized pairing fails here — at construction — rather than being
    /// pushed downstream into the chain builder.
    ///
    /// # Errors
    ///
    /// Returns an error if `paths.len() != read_structures.len()`.
    pub fn fastqs(paths: Vec<PathBuf>, read_structures: Vec<ReadStructure>) -> Result<Self> {
        ensure!(
            paths.len() == read_structures.len(),
            "FASTQ source is mis-sized: {} path(s) but {} read structure(s); \
             each path needs exactly one read structure",
            paths.len(),
            read_structures.len(),
        );
        Ok(Self::Fastqs { paths, read_structures })
    }

    /// Re-checks structural invariants that the variant types cannot
    /// encode. Currently this validates that a [`SourceSpec::Fastqs`]
    /// source has one read structure per path. Call this when accepting a
    /// `SourceSpec` that may have been built via a struct literal rather
    /// than [`SourceSpec::fastqs`].
    ///
    /// # Errors
    ///
    /// Returns an error if a `Fastqs` source has diverging `paths` and
    /// `read_structures` lengths.
    pub fn validate(&self) -> Result<()> {
        if let Self::Fastqs { paths, read_structures } = self {
            ensure!(
                paths.len() == read_structures.len(),
                "FASTQ source is mis-sized: {} path(s) but {} read structure(s); \
                 each path needs exactly one read structure",
                paths.len(),
                read_structures.len(),
            );
        }
        Ok(())
    }

    /// Returns `true` if this source is permissible inside runall.
    /// SAM input is standalone-only; runall always reads from one of
    /// `Bam`, `PairedBams`, or `Fastqs`.
    #[must_use]
    pub fn is_runall_compatible(&self) -> bool {
        match self {
            Self::Bam(_) | Self::PairedBams { .. } | Self::Fastqs { .. } => true,
            Self::Sam(_) => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runall_compatibility() {
        assert!(SourceSpec::Bam(PathBuf::from("a.bam")).is_runall_compatible());
        assert!(!SourceSpec::Sam(PathBuf::from("a.sam")).is_runall_compatible());
    }

    fn read_structure() -> ReadStructure {
        use std::str::FromStr;
        ReadStructure::from_str("8M+T").expect("valid read structure")
    }

    #[test]
    fn fastqs_constructor_accepts_matched_lengths() {
        let spec = SourceSpec::fastqs(
            vec![PathBuf::from("r1.fq"), PathBuf::from("r2.fq")],
            vec![read_structure(), read_structure()],
        )
        .expect("matched lengths should construct");
        assert!(spec.validate().is_ok());
    }

    #[test]
    fn fastqs_constructor_rejects_mismatched_lengths() {
        let err = SourceSpec::fastqs(
            vec![PathBuf::from("r1.fq"), PathBuf::from("r2.fq")],
            vec![read_structure()],
        )
        .expect_err("mismatched lengths should fail");
        assert!(err.to_string().contains("mis-sized"), "got: {err}");
    }

    #[test]
    fn validate_catches_mismatched_struct_literal() {
        // A directly-constructed (struct literal) mis-sized Fastqs source
        // is caught by validate(), even though the constructor was bypassed.
        let spec = SourceSpec::Fastqs {
            paths: vec![PathBuf::from("r1.fq")],
            read_structures: vec![read_structure(), read_structure()],
        };
        assert!(spec.validate().is_err());
    }

    #[test]
    fn validate_is_ok_for_non_fastq_sources() {
        assert!(SourceSpec::Bam(PathBuf::from("a.bam")).validate().is_ok());
        assert!(SourceSpec::Sam(PathBuf::from("a.sam")).validate().is_ok());
    }
}
