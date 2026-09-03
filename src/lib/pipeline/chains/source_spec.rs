//! What feeds the chain.

use std::path::PathBuf;

use crate::read_structure::ReadStructure;
use anyhow::{Result, ensure};

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

    /// A single interleaved FASTQ stream (`R1, R2, R1, R2, …`, or `-` for
    /// stdin) that the chain de-interleaves into two logical reads. Used by
    /// `fgumi extract --interleaved`.
    ///
    /// Distinct from [`SourceSpec::Fastqs`] because interleaved input is one
    /// physical `path` carrying **two** `read_structures` — the R1/R2 pair —
    /// which would violate the `Fastqs` "one read structure per path"
    /// invariant. Build this via [`SourceSpec::interleaved_fastq`] (which
    /// enforces exactly two read structures); [`SourceSpec::validate`]
    /// re-checks it.
    InterleavedFastq { path: PathBuf, read_structures: Vec<ReadStructure> },

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
        let spec = Self::Fastqs { paths, read_structures };
        // Route through validate() so the length invariant lives in exactly one
        // place; a struct-literal Fastqs is re-checked by the same code.
        spec.validate()?;
        Ok(spec)
    }

    /// Constructs a [`SourceSpec::InterleavedFastq`] source, enforcing that
    /// the single interleaved `path` carries exactly two read structures
    /// (the R1/R2 pair the stream de-interleaves into).
    ///
    /// Prefer this over an `InterleavedFastq { .. }` struct literal so a
    /// wrong read-structure count fails here rather than downstream.
    ///
    /// # Errors
    ///
    /// Returns an error if `read_structures.len() != 2`.
    pub fn interleaved_fastq(path: PathBuf, read_structures: Vec<ReadStructure>) -> Result<Self> {
        let spec = Self::InterleavedFastq { path, read_structures };
        spec.validate()?;
        Ok(spec)
    }

    /// Re-checks structural invariants that the variant types cannot
    /// encode: a [`SourceSpec::Fastqs`] source has one read structure per
    /// path, and a [`SourceSpec::InterleavedFastq`] source has exactly two.
    /// Call this when accepting a `SourceSpec` that may have been built via a
    /// struct literal rather than [`SourceSpec::fastqs`] /
    /// [`SourceSpec::interleaved_fastq`].
    ///
    /// # Errors
    ///
    /// Returns an error if a `Fastqs` source has diverging `paths` and
    /// `read_structures` lengths, or an `InterleavedFastq` source does not
    /// carry exactly two read structures.
    pub fn validate(&self) -> Result<()> {
        match self {
            Self::Fastqs { paths, read_structures } => {
                ensure!(
                    paths.len() == read_structures.len(),
                    "FASTQ source is mis-sized: {} path(s) but {} read structure(s); \
                     each path needs exactly one read structure",
                    paths.len(),
                    read_structures.len(),
                );
            }
            Self::InterleavedFastq { read_structures, .. } => {
                ensure!(
                    read_structures.len() == 2,
                    "interleaved FASTQ source needs exactly two read structures \
                     (R1/R2); got {}",
                    read_structures.len(),
                );
            }
            Self::Bam(_) | Self::PairedBams { .. } | Self::Sam(_) => {}
        }
        Ok(())
    }

    /// Returns `true` if this source is permissible inside runall.
    /// SAM input is standalone-only; runall always reads from one of
    /// `Bam`, `PairedBams`, or `Fastqs`.
    #[must_use]
    pub fn is_runall_compatible(&self) -> bool {
        match self {
            Self::Bam(_)
            | Self::PairedBams { .. }
            | Self::Fastqs { .. }
            | Self::InterleavedFastq { .. } => true,
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

    #[test]
    fn interleaved_constructor_accepts_exactly_two_read_structures() {
        let spec = SourceSpec::interleaved_fastq(
            PathBuf::from("interleaved.fq"),
            vec![read_structure(), read_structure()],
        )
        .expect("two read structures should construct");
        assert!(spec.validate().is_ok());
        assert!(spec.is_runall_compatible());
    }

    #[rstest::rstest]
    #[case::one(1)]
    #[case::three(3)]
    fn interleaved_constructor_rejects_non_two_read_structures(#[case] n: usize) {
        let err = SourceSpec::interleaved_fastq(
            PathBuf::from("interleaved.fq"),
            std::iter::repeat_with(read_structure).take(n).collect(),
        )
        .expect_err("only two read structures are valid");
        assert!(err.to_string().contains("exactly two read structures"), "got: {err}");
    }

    #[test]
    fn validate_catches_mismatched_interleaved_struct_literal() {
        let spec = SourceSpec::InterleavedFastq {
            path: PathBuf::from("interleaved.fq"),
            read_structures: vec![read_structure()],
        };
        assert!(spec.validate().is_err());
    }
}
