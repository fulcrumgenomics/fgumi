//! Where the chain writes its output.

use std::path::PathBuf;

/// The sink for a chain.
///
/// **Note for programmatic callers constructing a `ChainSpec` directly:**
/// `Bam` is the default — picking it on a sort-terminal chain produces
/// a coordinate-sorted BAM with **no** companion `.bai`. If you want a
/// sidecar BAI, use `BamWithIndex`; the chain-builder will register an
/// `IndexBamFinalizeHook` that reads the finished BAM and emits
/// `<output>.bam.bai`. There is no "auto-index" inference from sort
/// order; the caller must opt in via this variant.
///
/// `BamWithIndex` is sort-only (other stages can't produce a coordinate-
/// sorted output that warrants a `.bai` write). The sort command's
/// post-pipeline `--write-index` flag selects this variant; standalone
/// commands that don't sort use `Bam`. Cross-stage validator Rule 3
/// (in `chains::validate`) rejects `BamWithIndex` on any chain whose
/// terminal stage isn't `Stage::Sort`.
///
/// Future CRAM output should be parallel `Cram(PathBuf)` /
/// `CramWithIndex(PathBuf)` variants — not a `format: IndexFormat`
/// field grafted onto `BamWithIndex`, because BAI vs. CRAI is a
/// file-format-physics distinction (different on-disk encoding for the
/// same logical chunk-offset index), not a runtime choice.
#[derive(Debug, Clone)]
pub enum SinkSpec {
    /// BAM file write (or `-` for stdout). No sidecar index.
    Bam(PathBuf),
    /// BAM file write + `.bai` index write after the BAM finishes.
    /// Only valid when the final stage emits coordinate-sorted output.
    BamWithIndex(PathBuf),
}

impl SinkSpec {
    /// The output path, regardless of variant.
    #[must_use]
    pub fn path(&self) -> &PathBuf {
        match self {
            SinkSpec::Bam(p) | SinkSpec::BamWithIndex(p) => p,
        }
    }
}
