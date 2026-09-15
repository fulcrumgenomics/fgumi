//! Where the chain writes its output.

use std::path::PathBuf;

/// The sink for a chain.
///
/// **Note for programmatic callers constructing a `ChainSpec` directly:**
/// `Bam` is the default — picking it on a sort-terminal chain produces
/// a coordinate-sorted BAM with **no** companion `.bai`. If you want a
/// sidecar BAI, use `BamWithIndex`; the chain builder's `add_sink` attaches an
/// inline indexer directly to the output sink, which builds the `.bai`
/// alongside the BAM as it writes (no post-pipeline re-read) and emits
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
    /// Interleaved FASTQ file write (or `-` for stdout). A `.gz`/`.bgz`
    /// path is written as BGZF. Only valid on a `Stage::Fastq`-terminal chain.
    ///
    /// Paired split output uses the [`SinkSpec::FastqPaired`] variant.
    Fastq(PathBuf),
    /// Paired split FASTQ output: R1 → `out1`, R2 → `out2`, and "other"
    /// (single-end / ambiguous) reads → `out0` if given, else stdout. Each
    /// `.gz`/`.bgz` path is written as BGZF. Only valid on a `Stage::Fastq`-
    /// terminal chain. The three streams are fanned out by the paired FASTQ
    /// encode step (a 3-output `Process3WithWorkerState`).
    FastqPaired { out1: PathBuf, out2: PathBuf, out0: Option<PathBuf> },
    /// No output file. The terminal stage produces no BAM — it terminates the
    /// chain itself with a `type Outputs = ()` sink step (e.g. the metrics
    /// stage records into a per-thread accumulator and writes its own TSVs via
    /// a finalize hook). `add_sink` is a no-op for this variant. Only valid on
    /// a `Stage::Metrics`-terminal chain (enforced by a cross-stage validator
    /// biconditional).
    None,
}

impl SinkSpec {
    /// A representative output path for the sink, used for logging (and
    /// available to a caller's input-clobber guard). For [`SinkSpec::FastqPaired`]
    /// this is `out1`; the chains layer itself does not currently guard against
    /// clobbering the input. Returns `None` for [`SinkSpec::None`] (no output
    /// file at all).
    #[must_use]
    pub fn path(&self) -> Option<&PathBuf> {
        match self {
            SinkSpec::Bam(p) | SinkSpec::BamWithIndex(p) | SinkSpec::Fastq(p) => Some(p),
            SinkSpec::FastqPaired { out1, .. } => Some(out1),
            SinkSpec::None => None,
        }
    }

    /// The output path for a file-writing sink, for the BAM/FASTQ stage methods
    /// that always run against a real output. Every stage that calls this pairs
    /// with a file sink (`Bam`/`BamWithIndex`/`Fastq`/`FastqPaired`); only the
    /// metrics stage uses [`SinkSpec::None`], and it never calls this. The
    /// biconditional in `validate_cross_stage_constraints` guarantees a
    /// `None` sink only ever appears with a `Stage::Metrics` terminal, so this
    /// is unreachable for those callers.
    ///
    /// # Panics
    ///
    /// Panics if the sink is [`SinkSpec::None`] — a framework bug (a
    /// file-output stage wired against a no-output sink).
    #[must_use]
    pub fn output_path(&self) -> &PathBuf {
        self.path().expect(
            "SinkSpec::output_path called on a SinkSpec::None sink; only the metrics stage \
             uses a None sink and it never writes an output file",
        )
    }
}
