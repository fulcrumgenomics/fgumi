//! Where the chain writes its output.

use std::path::{Path, PathBuf};

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
/// CRAM output has its own parallel `Cram` / `CramWithIndex` variants
/// rather than a `format: IndexFormat` field grafted onto `BamWithIndex`,
/// because BAI vs. CRAI is a file-format-physics distinction (different
/// on-disk encoding for the same logical chunk-offset index), not a runtime
/// choice.
#[derive(Debug, Clone)]
pub enum SinkSpec {
    /// BAM file write (or `-` for stdout). No sidecar index.
    Bam(PathBuf),
    /// BAM file write + `.bai` index write after the BAM finishes.
    /// Only valid when the final stage emits coordinate-sorted output.
    BamWithIndex(PathBuf),
    /// CRAM file write (feature `cram`; without it the chain builder rejects
    /// it). Chosen by [`SinkSpec::alignment`] when the user-facing output name
    /// ends in `.cram`. No sidecar index. `path` is where the bytes go (a
    /// staged temporary for a sort); `name` is the destination the user named,
    /// which becomes the CRAM file-definition ID.
    Cram { path: PathBuf, name: PathBuf },
    /// CRAM file write + `.crai` (see `fgumi_bam_io::crai_sidecar_path`) built
    /// as the file is written. The same rules as [`SinkSpec::BamWithIndex`] apply.
    CramWithIndex { path: PathBuf, name: PathBuf },
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

pub(crate) use fgumi_bam_io::is_cram_path;

/// The kind of coordinate index a sink writes beside its output.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexKind {
    /// A `.bai` beside a BAM.
    Bai,
    /// A `.crai` beside a CRAM.
    Crai,
}

impl IndexKind {
    /// The index an output named `output` gets: CRAI for `*.cram`, else BAI.
    #[must_use]
    pub fn for_path(output: &Path) -> Self {
        if is_cram_path(output) { IndexKind::Crai } else { IndexKind::Bai }
    }

    /// The index's name for messages (`"BAI"` or `"CRAI"`).
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            IndexKind::Bai => "BAI",
            IndexKind::Crai => "CRAI",
        }
    }

    /// How an output-collision error names the sidecar
    /// (`--write-index (.crai sidecar)`).
    #[must_use]
    pub fn sidecar_label(self) -> &'static str {
        match self {
            IndexKind::Bai => "--write-index (.bai sidecar)",
            IndexKind::Crai => "--write-index (.crai sidecar)",
        }
    }

    /// The sidecar written beside `output` (`<output>.bai` / `<output>.crai`).
    #[must_use]
    pub fn sidecar(self, output: &Path) -> PathBuf {
        match self {
            IndexKind::Bai => fgumi_bam_io::bai_sidecar_path(output),
            IndexKind::Crai => fgumi_bam_io::crai_sidecar_path(output),
        }
    }
}

impl SinkSpec {
    /// An alignment-file sink writing to `write_path`, in the format the
    /// user-facing output name `named` asks for: CRAM for `*.cram`, else BAM.
    /// The two paths differ when the output is staged under a temporary name
    /// and renamed into place.
    #[must_use]
    pub fn alignment(write_path: PathBuf, named: &Path, index: bool) -> Self {
        match (is_cram_path(named), index) {
            (true, true) => SinkSpec::CramWithIndex { path: write_path, name: named.to_path_buf() },
            (true, false) => SinkSpec::Cram { path: write_path, name: named.to_path_buf() },
            (false, true) => SinkSpec::BamWithIndex(write_path),
            (false, false) => SinkSpec::Bam(write_path),
        }
    }

    /// Whether this sink writes CRAM.
    #[must_use]
    pub fn is_cram(&self) -> bool {
        matches!(self, SinkSpec::Cram { .. } | SinkSpec::CramWithIndex { .. })
    }

    /// Whether this sink also writes a coordinate index (`.bai` or `.crai`).
    #[must_use]
    pub fn has_index(&self) -> bool {
        self.index().is_some()
    }

    /// The coordinate index this sink writes, if any: its kind and the
    /// sidecar path (beside the path the bytes go to).
    #[must_use]
    pub fn index(&self) -> Option<(IndexKind, PathBuf)> {
        match self {
            SinkSpec::BamWithIndex(path) => Some((IndexKind::Bai, IndexKind::Bai.sidecar(path))),
            SinkSpec::CramWithIndex { path, .. } => {
                Some((IndexKind::Crai, IndexKind::Crai.sidecar(path)))
            }
            _ => None,
        }
    }

    /// For a CRAM sink, the destination name the user gave (see
    /// [`SinkSpec::Cram`]).
    #[must_use]
    pub fn cram_name(&self) -> Option<&Path> {
        match self {
            SinkSpec::Cram { name, .. } | SinkSpec::CramWithIndex { name, .. } => Some(name),
            _ => None,
        }
    }

    /// The variant's name, for error messages (`SinkSpec::CramWithIndex`).
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            SinkSpec::Bam(_) => "SinkSpec::Bam",
            SinkSpec::BamWithIndex(_) => "SinkSpec::BamWithIndex",
            SinkSpec::Cram { .. } => "SinkSpec::Cram",
            SinkSpec::CramWithIndex { .. } => "SinkSpec::CramWithIndex",
            SinkSpec::Fastq(_) => "SinkSpec::Fastq",
            SinkSpec::FastqPaired { .. } => "SinkSpec::FastqPaired",
            SinkSpec::None => "SinkSpec::None",
        }
    }

    /// The BGZF compression settings a terminal sort's merge may compress its
    /// output with for this sink: `None` for CRAM (the CRAM sink takes raw
    /// blocks), else BAM at `level`, indexing when the sink writes a `.bai`.
    #[must_use]
    pub fn merge_output_compression(
        &self,
        level: u32,
    ) -> Option<fgumi_pipeline_io::sort::OutputCompression> {
        (!self.is_cram()).then(|| fgumi_pipeline_io::sort::OutputCompression {
            level,
            index_bam: self.has_index(),
        })
    }

    /// A representative output path for the sink, used for logging (and
    /// available to a caller's input-clobber guard). For [`SinkSpec::FastqPaired`]
    /// this is `out1`; the chains layer itself does not currently guard against
    /// clobbering the input. Returns `None` for [`SinkSpec::None`] (no output
    /// file at all).
    #[must_use]
    pub fn path(&self) -> Option<&PathBuf> {
        match self {
            SinkSpec::Bam(p)
            | SinkSpec::BamWithIndex(p)
            | SinkSpec::Cram { path: p, .. }
            | SinkSpec::CramWithIndex { path: p, .. }
            | SinkSpec::Fastq(p) => Some(p),
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
