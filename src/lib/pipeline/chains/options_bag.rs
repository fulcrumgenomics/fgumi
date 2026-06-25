//! Typed option bag for [`super::ChainSpec`].

use std::path::PathBuf;

// Re-export each command's options type so spec construction reads
// naturally at the call site. Each is imported from its existing
// home, not redefined.
pub use crate::aligner::AlignerOptions;
pub use crate::commands::clip::Clip;
#[cfg(feature = "consensus")]
pub use crate::commands::codec::CodecOptions;
pub use crate::commands::correct::CorrectOptions;
pub use crate::commands::dedup::MarkDuplicates;
#[cfg(feature = "consensus")]
pub use crate::commands::duplex::DuplexOptions;
pub use crate::commands::extract::ExtractOptions;
pub use crate::commands::filter::FilterOptions;
pub use crate::commands::group::GroupOptions;
#[cfg(feature = "consensus")]
pub use crate::commands::simplex::SimplexOptions;
pub use crate::commands::sort::SortOptions;
pub use crate::commands::zipper::ZipperOptions;

/// Per-stage options for [`crate::pipeline::chains::Stage::Align`].
///
/// Bundles the aligner tuning knobs ([`AlignerOptions`]) with the two
/// additional parameters that [`AlignerOptions::resolve`] needs but
/// that `AlignerOptions` does not carry directly:
///
/// - `reference` — the reference FASTA path (required; used to locate
///   the accompanying `.dict` file and to validate aligner index files
///   in preset mode).
/// - `aligner_bin` — optional binary path override (preset mode only;
///   `None` in command mode).
///
/// Constructed by runall's `--start-from align` dispatch (T3b.3) and
/// by integration tests that exercise the chain builder directly.
pub struct AlignOptions {
    /// Per-stage aligner tuning: `--aligner::preset`, `--aligner::command`,
    /// `--aligner::threads`, `--aligner::chunk-size`.
    pub aligner: AlignerOptions,
    /// Reference FASTA path. Must have a `.dict` sibling (or
    /// `<ref>.dict` fallback). Passed to [`AlignerOptions::resolve`].
    pub reference: PathBuf,
    /// Optional aligner binary override for preset mode (`--aligner-bin`
    /// in runall). Ignored in command mode.
    pub aligner_bin: Option<PathBuf>,
}

/// Per-stage option bag for [`crate::pipeline::chains::ChainSpec`].
///
/// Each field is `Option<T>` so callers fill in only the stages they
/// need. [`crate::pipeline::chains::build_for`] validates
/// that every stage in `spec.stages` has its matching options
/// present before constructing any chain steps.
///
/// Twelve fields are wired: Correct, Zipper, Sort, Group, Duplex, Codec,
/// Extract (free-standing `<Command>Options` structs) and Dedup, Filter,
/// Clip, and Simplex (hold the command struct directly — no separate
/// `DedupOptions`/`FilterOptions`/`ClipOptions`/`SimplexOptions` is
/// extracted because test code constructs those structs via
/// `clap::Parser::parse_from` and extracting a sub-struct would
/// destabilise that surface). Remaining slot (Downsample) is added
/// incrementally in T2.19–T2.22.
#[derive(Default)]
pub struct StageOptionsBag {
    pub correct: Option<CorrectOptions>,
    pub zipper: Option<ZipperOptions>,
    pub sort: Option<SortOptions>,
    pub group: Option<GroupOptions>,
    #[cfg(feature = "consensus")]
    pub duplex: Option<DuplexOptions>,
    #[cfg(feature = "consensus")]
    pub codec: Option<CodecOptions>,
    /// Dedup options. Holds `MarkDuplicates` directly (approach 2 from
    /// the T2.12 design note) to keep the test surface stable.
    pub dedup: Option<MarkDuplicates>,
    /// Filter options. Holds the free-standing `FilterOptions` struct
    /// (same approach as Sort/Group/Correct), so runall can re-expose
    /// the tuning knobs via `--filter::*` through `MultiFilterOptions`.
    pub filter: Option<FilterOptions>,
    /// Clip options. Holds `Clip` directly (same approach as Dedup/Filter).
    pub clip: Option<Clip>,
    /// Simplex options. Holds the free-standing `SimplexOptions` struct
    /// (same approach as Duplex/Codec), so runall can re-expose the tuning
    /// knobs via `--simplex::*` through `MultiSimplexOptions`.
    #[cfg(feature = "consensus")]
    pub simplex: Option<SimplexOptions>,
    /// Align (`AlignAndMerge`) options. Carries [`AlignerOptions`] +
    /// reference path + optional aligner-binary override.
    pub aligner: Option<AlignOptions>,
    /// Extract options. Carries header-synthesis fields and behavior knobs
    /// for the FASTQ→unmapped-BAM extract stage (Phase 5 T5.3).
    pub extract: Option<ExtractOptions>,
    // Additional fields added during command migrations:
    //   - downsample (T2.13).
    // Cross-stage fields (read_group, methylation_mode, reference, …)
    // resolve in Phase 3 when runall lands on the façade.
}
