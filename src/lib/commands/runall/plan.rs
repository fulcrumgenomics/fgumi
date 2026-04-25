//! Runall pipeline plan.
//!
//! A [`Plan`] is a pure-data description of a pipeline: one [`SourceSpec`], a
//! `Vec<StageSpec>` of middle stages, and one [`SinkSpec`]. The planner
//! ([`super::planner`]) builds a `Plan` from `Runall` CLI arguments; the
//! runner ([`super::runner`]) instantiates a `Plan` into concrete types and
//! dispatches to [`run_multi_stage`](crate::runall::engine::driver::run_multi_stage).
//!
//! This split decouples "what pipeline to run" (planning) from "how to run it"
//! (execution), removes duplicated per-entry-point builder code, and makes the
//! pipeline observable as data.

use std::path::PathBuf;
use std::sync::Arc;

use fgumi_consensus::filter::FilterConfig;
use fgumi_umi::{Strategy, TagInfo};
use noodles::sam::Header;
use read_structure::ReadStructure;

use crate::commands::filter::CollectedFilterMetrics;
use crate::correct::EncodedUmiSet;
use crate::extract::QualityEncoding;
use crate::grouper::GroupFilterConfig;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::runall::engine::driver::PipelineConfig;
use crate::runall::engine::stages::consensus::CallerFactory;
use crate::runall::engine::stages::correct::CorrectMetricsAccumulator;

// ---------------------------------------------------------------------------
// Source specs
// ---------------------------------------------------------------------------

/// Source stage configuration — the first stage in a pipeline.
#[derive(Debug, Clone)]
pub enum SourceSpec {
    /// Read FASTQ files (gzip/bgzf/plain). Emits `PerStreamChunk`.
    FastqFileRead { paths: Vec<PathBuf> },
    /// Read a BAM file and emit batches of length-prefixed records.
    /// Emits `SerializedBatch`.
    BamBatchedSource { path: PathBuf, threads: usize },
}

// ---------------------------------------------------------------------------
// Sink specs
// ---------------------------------------------------------------------------

/// Sink stage configuration — the last stage in a pipeline.
#[derive(Debug, Clone)]
#[allow(
    clippy::large_enum_variant,
    reason = "SinkSpec is constructed once per pipeline run; size is not hot"
)]
pub enum SinkSpec {
    /// Write a BAM file via the compressed byte path. Consumes
    /// `CompressedBatch` (typically produced by `BgzfCompress`).
    BamFileWrite {
        path: PathBuf,
        header: Header,
        secondary_path: Option<PathBuf>,
        secondary_header: Option<Header>,
        queue_capacity: usize,
    },
    /// Write interleaved FASTQ. Consumes `SerializedFastqBatch`.
    FastqFileSink { path: PathBuf, queue_capacity: usize },
    /// Consensus-output BAM sink (internal compression). Consumes
    /// `ConsensusOutput`.
    ///
    /// When `deferred_header` is `Some`, the sink blocks until the slot is
    /// populated by `StageSpec::AlignAndMerge` and uses that value as the
    /// BAM writer header — otherwise it uses `header` (the plan-time
    /// approximation built from the dict + unmapped header, missing the
    /// aligner's @PG lines).
    BamSink {
        path: PathBuf,
        header: Header,
        threads: usize,
        compression_level: u32,
        deferred_header: Option<crate::runall::engine::sink::DeferredHeader>,
    },
}

// ---------------------------------------------------------------------------
// Middle-stage specs
// ---------------------------------------------------------------------------

/// Middle-stage configuration — zero or more of these live between a
/// [`SourceSpec`] and a [`SinkSpec`] in a [`Plan`].
#[derive(Clone)]
pub enum StageSpec {
    // BGZF I/O
    BgzfDecompress,
    BgzfCompress {
        level: u32,
    },

    // FASTQ parsing (record-aligned path)
    FastqParse,
    FastqPair {
        num_streams: usize,
    },

    // FASTQ parsing (raw-block path)
    FastqBlockParse,
    FastqBlockMerge {
        num_streams: usize,
    },

    // Extract: FASTQ templates → unmapped BAM bytes
    Extract {
        params: Arc<ExtractParamsPlanData>,
        read_structures: Arc<Vec<ReadStructure>>,
        encoding: QualityEncoding,
    },

    // UMI correction
    Correct {
        state: CorrectionPlanState,
    },

    // BAM → interleaved FASTQ bytes
    ToFastq,

    // Subprocess aligner + zipper merge
    AlignAndMerge {
        aligner_command: String,
        reference: PathBuf,
        dict_path: PathBuf,
        unmapped_header: Header,
        tag_info: TagInfo,
        skip_pa_tags: bool,
        no_read_suffix: bool,
        /// Shared slot the stage populates at runtime with the merged output
        /// header (built from unmapped + aligner-produced mapped header +
        /// dict). When paired with a `SinkSpec::BamSink` that has the same
        /// slot, the sink uses this header instead of its plan-time value.
        deferred_header: Option<crate::runall::engine::sink::DeferredHeader>,
    },

    // Template-coordinate sort (SpecialStage, barrier)
    Sort {
        memory_limit: usize,
        temp_dir: Option<PathBuf>,
        threads: usize,
        header: Header,
    },

    // Group records by position key, emit PositionGroupBatch
    PositionBatch {
        header: Header,
    },

    // Group-assign: PositionGroupBatch → MiGroupBatch (local MoleculeIds,
    // MI tag NOT yet injected — MiAssignStage does that downstream)
    GroupAssign {
        umi_tag: [u8; 2],
        assign_tag: [u8; 2],
        strategy: Strategy,
        max_edits: u32,
        filter_config: GroupFilterConfig,
    },

    // MI-assign: serial-ordered cumulative offset on local MoleculeIds.
    // Rewrites every MiGroup's local id into a global id and injects the
    // final MI tag bytes into each record. Mirrors the `mi_assign_fn` zone
    // PR #319 added to the unified BAM pipeline.
    MiAssign,

    // Consensus calling: MiGroupBatch → ConsensusOutput
    Consensus {
        factory: CallerFactory,
        /// When true, apply overlapping-bases correction (matches standalone
        /// `fgumi simplex --consensus-call-overlapping-bases`, default true).
        call_overlapping_bases: bool,
    },

    // Filter consensus output
    Filter {
        config: FilterConfig,
        filter_by_template: bool,
        /// Optional per-thread accumulator for filter stats. When set, each
        /// worker merges into its own slot; the planner folds slots and
        /// writes the `--filter::stats` TSV at pipeline end.
        metrics: Option<Arc<PerThreadAccumulator<CollectedFilterMetrics>>>,
    },

    // Reorder + coalesce SerializedBatch → SerializedBatch
    CoalesceSerialized {
        target_chunk_size: usize,
    },

    // Reorder + coalesce MiGroupBatch → SerializedBatch
    CoalesceMiGroup {
        target_chunk_size: usize,
    },
}

impl std::fmt::Debug for StageSpec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

impl StageSpec {
    /// Short debug name. Matches the stage's own `name()` where applicable.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            StageSpec::BgzfDecompress => "BgzfDecompress",
            StageSpec::BgzfCompress { .. } => "BgzfCompress",
            StageSpec::FastqParse => "FastqParse",
            StageSpec::FastqPair { .. } => "FastqPair",
            StageSpec::FastqBlockParse => "FastqBlockParse",
            StageSpec::FastqBlockMerge { .. } => "FastqBlockMerge",
            StageSpec::Extract { .. } => "Extract",
            StageSpec::Correct { .. } => "Correct",
            StageSpec::ToFastq => "ToFastq",
            StageSpec::AlignAndMerge { .. } => "AlignAndMerge",
            StageSpec::Sort { .. } => "Sort",
            StageSpec::PositionBatch { .. } => "PositionBatch",
            StageSpec::GroupAssign { .. } => "GroupAssign",
            StageSpec::MiAssign => "MiAssign",
            StageSpec::Consensus { .. } => "Consensus",
            StageSpec::Filter { .. } => "Filter",
            StageSpec::CoalesceSerialized { .. } => "CoalesceSerialized",
            StageSpec::CoalesceMiGroup { .. } => "CoalesceMiGroup",
        }
    }
}

// ---------------------------------------------------------------------------
// Correction plan state (needed because CorrectStage isn't Clone)
// ---------------------------------------------------------------------------

/// Shared correction state, cloned by reference into the per-worker factory
/// that [`StageSpec::Correct`] produces.
#[derive(Clone)]
pub struct CorrectionPlanState {
    pub encoded: Arc<EncodedUmiSet>,
    pub umi_length: usize,
    pub max_mismatches: usize,
    pub min_distance_diff: usize,
    pub revcomp: bool,
    pub dont_store_original_umis: bool,
    pub track_rejects: bool,
    pub cache_size: usize,
    pub metrics: Arc<PerThreadAccumulator<CorrectMetricsAccumulator>>,
}

// ---------------------------------------------------------------------------
// Extract params (needed because ExtractParams is constructed from CLI args)
// ---------------------------------------------------------------------------

/// Plan-side wrapper for `ExtractStage` parameters — this is just a re-alias
/// of [`crate::extract::ExtractParams`] to keep the `StageSpec` import
/// surface narrow.
pub type ExtractParamsPlanData = crate::extract::ExtractParams;

// ---------------------------------------------------------------------------
// Plan
// ---------------------------------------------------------------------------

/// A complete pipeline description.
///
/// Build with [`crate::commands::runall::Runall::build_plan`], run with
/// [`super::runner::run_plan`] (or [`Plan::run`]).
pub struct Plan {
    pub source: SourceSpec,
    pub stages: Vec<StageSpec>,
    pub sink: SinkSpec,
    pub config: PipelineConfig,
    /// Optional path to write per-stage metrics TSV after the pipeline
    /// completes successfully. Columns: `stage`, `wall_time_secs`,
    /// `records_in`, `records_out`.
    pub metrics_tsv: Option<PathBuf>,
}

impl Plan {
    /// Human-readable plan summary for diagnostics (`INFO` log).
    pub fn summary(&self) -> String {
        let source_name = match &self.source {
            SourceSpec::FastqFileRead { .. } => "FastqFileRead",
            SourceSpec::BamBatchedSource { .. } => "BamBatchedSource",
        };
        let sink_name = match &self.sink {
            SinkSpec::BamFileWrite { .. } => "BamFileWrite",
            SinkSpec::FastqFileSink { .. } => "FastqFileSink",
            SinkSpec::BamSink { .. } => "BamSink",
        };
        let stages: Vec<&str> = self.stages.iter().map(StageSpec::name).collect();
        format!("{source_name} → {} → {sink_name}", stages.join(" → "))
    }

    /// Format this plan as a human-readable multi-line report.
    ///
    /// Shows the resolved source, each stage with its configurable parameters,
    /// the sink, and runtime configuration. Intended for the `--explain` flag
    /// on `fgumi runall`, and structured to be roughly stable across runs so
    /// it could be snapshot-tested in the future.
    ///
    /// `command_line` is displayed in the header so the report identifies the
    /// exact invocation that produced it.
    pub fn explain(&self, command_line: &str) -> String {
        use std::fmt::Write as _;

        let mut out = String::new();
        let rule: String = "─".repeat(70);

        let _ = writeln!(out, "Pipeline plan for: {command_line}");
        let _ = writeln!(out, "{rule}");

        // --- Source ---
        let _ = writeln!(out, "Source:");
        match &self.source {
            SourceSpec::FastqFileRead { paths } => {
                let _ = writeln!(out, "    type: FASTQ batched");
                let joined =
                    paths.iter().map(|p| p.display().to_string()).collect::<Vec<_>>().join(", ");
                let _ = writeln!(out, "    inputs: {joined}");
            }
            SourceSpec::BamBatchedSource { path, threads } => {
                let _ = writeln!(out, "    type: BAM batched");
                let _ = writeln!(out, "    input: {}", path.display());
                let _ = writeln!(out, "    threads: {threads}");
            }
        }
        let _ = writeln!(out);

        // --- Stages ---
        let _ = writeln!(out, "Stages:");
        for (idx, stage) in self.stages.iter().enumerate() {
            let n = idx + 1;
            let name = stage.name();
            match stage {
                StageSpec::BgzfDecompress | StageSpec::FastqParse | StageSpec::FastqBlockParse => {
                    let _ = writeln!(out, "  {n:>2}. {name}");
                }
                StageSpec::ToFastq => {
                    let _ = writeln!(out, "  {n:>2}. {name}");
                }
                StageSpec::BgzfCompress { level } => {
                    let _ = writeln!(out, "  {n:>2}. {name:<20} (level: {level})");
                }
                StageSpec::FastqPair { num_streams }
                | StageSpec::FastqBlockMerge { num_streams } => {
                    let _ = writeln!(out, "  {n:>2}. {name:<20} (streams: {num_streams})");
                }
                StageSpec::Extract { params, read_structures, encoding } => {
                    let _ = writeln!(out, "  {n:>2}. {name}");
                    let rs = read_structures
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join(", ");
                    let _ = writeln!(out, "      read structures: {rs}");
                    let _ = writeln!(out, "      quality encoding: {encoding:?}");
                    let _ = writeln!(out, "      read group id: {}", params.read_group_id);
                    let _ = writeln!(
                        out,
                        "      umi tag: {}   cell tag: {}",
                        bytes_to_tag(params.umi_tag),
                        bytes_to_tag(params.cell_tag),
                    );
                    let _ = writeln!(
                        out,
                        "      annotate read names: {}   extract umis from read names: {}",
                        params.annotate_read_names, params.extract_umis_from_read_names
                    );
                }
                StageSpec::Correct { state } => {
                    let _ = writeln!(out, "  {n:>2}. {name}");
                    let _ = writeln!(out, "      umi length: {}", state.umi_length);
                    let _ = writeln!(out, "      max mismatches: {}", state.max_mismatches);
                    let _ = writeln!(out, "      min distance diff: {}", state.min_distance_diff);
                    let _ = writeln!(out, "      revcomp: {}", state.revcomp);
                    let _ = writeln!(out, "      cache size: {}", state.cache_size);
                    let _ = writeln!(
                        out,
                        "      track rejects: {}   store original umis: {}",
                        state.track_rejects, !state.dont_store_original_umis
                    );
                }
                StageSpec::AlignAndMerge {
                    aligner_command,
                    reference,
                    dict_path,
                    tag_info: _,
                    skip_pa_tags,
                    no_read_suffix,
                    unmapped_header: _,
                    deferred_header: _,
                } => {
                    let _ = writeln!(out, "  {n:>2}. {name}");
                    let _ = writeln!(out, "      aligner command: {aligner_command}");
                    let _ = writeln!(out, "      reference: {}", reference.display());
                    let _ = writeln!(out, "      dict: {}", dict_path.display());
                    let _ = writeln!(
                        out,
                        "      skip PA tags: {skip_pa_tags}   no read suffix: {no_read_suffix}"
                    );
                }
                StageSpec::Sort { memory_limit, temp_dir, threads, header: _ } => {
                    let _ = writeln!(out, "  {n:>2}. {name}");
                    let _ = writeln!(out, "      memory limit: {}", fmt_bytes(*memory_limit));
                    let _ = writeln!(out, "      threads: {threads}");
                    let _ = writeln!(
                        out,
                        "      temp dir: {}",
                        temp_dir
                            .as_ref()
                            .map(|p| p.display().to_string())
                            .unwrap_or_else(|| "(system temp)".to_string())
                    );
                }
                StageSpec::PositionBatch { header: _ } => {
                    let _ = writeln!(out, "  {n:>2}. {name}");
                }
                StageSpec::GroupAssign {
                    umi_tag,
                    assign_tag,
                    strategy,
                    max_edits,
                    filter_config,
                } => {
                    let _ = writeln!(out, "  {n:>2}. {name}");
                    let _ = writeln!(
                        out,
                        "      umi tag: {}   assign tag: {}",
                        bytes_to_tag(*umi_tag),
                        bytes_to_tag(*assign_tag)
                    );
                    let _ = writeln!(out, "      strategy: {strategy:?}");
                    let _ = writeln!(out, "      max edits: {max_edits}");
                    let _ = writeln!(
                        out,
                        "      min mapq: {}   include non-PF: {}",
                        filter_config.min_mapq, filter_config.include_non_pf
                    );
                    let _ = writeln!(
                        out,
                        "      allow unmapped: {}   no umi: {}",
                        filter_config.allow_unmapped, filter_config.no_umi
                    );
                    if let Some(min_umi) = filter_config.min_umi_length {
                        let _ = writeln!(out, "      min umi length: {min_umi}");
                    } else {
                        let _ = writeln!(out, "      min umi length: (not set)");
                    }
                }
                StageSpec::MiAssign => {
                    let _ = writeln!(out, "  {n:>2}. {name}");
                }
                StageSpec::Consensus { factory: _, call_overlapping_bases } => {
                    let _ = writeln!(out, "  {n:>2}. {name}");
                    let _ = writeln!(out, "      call overlapping bases: {call_overlapping_bases}");
                }
                StageSpec::Filter { config, filter_by_template, metrics: _ } => {
                    let _ = writeln!(out, "  {n:>2}. {name}");
                    let _ = writeln!(out, "      filter by template: {filter_by_template}");
                    let _ = writeln!(
                        out,
                        "      max no-call fraction: {}",
                        config.max_no_call_fraction
                    );
                    let _ = writeln!(
                        out,
                        "      min base quality: {}",
                        opt_display(config.min_base_quality)
                    );
                    let _ = writeln!(
                        out,
                        "      min mean base quality: {}",
                        opt_display(config.min_mean_base_quality)
                    );
                    if let Some(t) = config.single_strand_thresholds.as_ref() {
                        let _ = writeln!(
                            out,
                            "      single-strand thresholds: min reads: {}, max read err: {}, max base err: {}",
                            t.min_reads, t.max_read_error_rate, t.max_base_error_rate
                        );
                    }
                    if let Some(t) = config.duplex_thresholds.as_ref() {
                        let _ = writeln!(
                            out,
                            "      duplex thresholds: min reads: {}, max read err: {}, max base err: {}",
                            t.min_reads, t.max_read_error_rate, t.max_base_error_rate
                        );
                    }
                    if let Some(t) = config.ab_thresholds.as_ref() {
                        let _ = writeln!(
                            out,
                            "      AB thresholds: min reads: {}, max read err: {}, max base err: {}",
                            t.min_reads, t.max_read_error_rate, t.max_base_error_rate
                        );
                    }
                    if let Some(t) = config.ba_thresholds.as_ref() {
                        let _ = writeln!(
                            out,
                            "      BA thresholds: min reads: {}, max read err: {}, max base err: {}",
                            t.min_reads, t.max_read_error_rate, t.max_base_error_rate
                        );
                    }
                }
                StageSpec::CoalesceSerialized { target_chunk_size }
                | StageSpec::CoalesceMiGroup { target_chunk_size } => {
                    let _ = writeln!(
                        out,
                        "  {n:>2}. {name:<20} (target chunk size: {})",
                        fmt_bytes(*target_chunk_size)
                    );
                }
            }
        }
        let _ = writeln!(out);

        // --- Sink ---
        let _ = writeln!(out, "Sink:");
        match &self.sink {
            SinkSpec::BamFileWrite {
                path,
                secondary_path,
                queue_capacity,
                header: _,
                secondary_header: _,
            } => {
                let _ = writeln!(out, "    type: BAM file write (compressed-byte path)");
                let _ = writeln!(out, "    path: {}", path.display());
                let _ = writeln!(
                    out,
                    "    secondary path: {}",
                    secondary_path
                        .as_ref()
                        .map(|p| p.display().to_string())
                        .unwrap_or_else(|| "(not set)".to_string())
                );
                let _ = writeln!(out, "    queue capacity: {queue_capacity}");
            }
            SinkSpec::FastqFileSink { path, queue_capacity } => {
                let _ = writeln!(out, "    type: interleaved FASTQ");
                let _ = writeln!(out, "    path: {}", path.display());
                let _ = writeln!(out, "    queue capacity: {queue_capacity}");
            }
            SinkSpec::BamSink {
                path,
                threads,
                compression_level,
                header: _,
                deferred_header: _,
            } => {
                let _ = writeln!(out, "    type: BAM sink (internal compression)");
                let _ = writeln!(out, "    path: {}", path.display());
                let _ = writeln!(out, "    threads: {threads}");
                let _ = writeln!(out, "    compression level: {compression_level}");
            }
        }
        let _ = writeln!(out);

        // --- Runtime ---
        let _ = writeln!(out, "Runtime:");
        let _ = writeln!(out, "    worker threads: {}", self.config.worker_threads);
        let _ = writeln!(out, "    queue capacity: {} items", self.config.queue_capacity);
        let _ =
            writeln!(out, "    queue memory limit: {}", fmt_bytes(self.config.queue_memory_limit));
        let _ = writeln!(
            out,
            "    global memory limit: {}",
            fmt_bytes(self.config.global_memory_limit)
        );
        let _ = writeln!(
            out,
            "    metrics tsv: {}",
            self.metrics_tsv
                .as_ref()
                .map(|p| p.display().to_string())
                .unwrap_or_else(|| "(not set)".to_string())
        );

        out
    }
}

// ---------------------------------------------------------------------------
// Explain-report helpers
// ---------------------------------------------------------------------------

fn bytes_to_tag(b: [u8; 2]) -> String {
    String::from_utf8(b.to_vec()).unwrap_or_else(|_| format!("{b:?}"))
}

fn opt_display<T: std::fmt::Display>(v: Option<T>) -> String {
    match v {
        Some(x) => format!("{x}"),
        None => "(not set)".to_string(),
    }
}

/// Format a byte count as KiB / MiB / GiB for readability.
fn fmt_bytes(bytes: usize) -> String {
    const KIB: usize = 1024;
    const MIB: usize = 1024 * KIB;
    const GIB: usize = 1024 * MIB;
    if bytes >= GIB && bytes.is_multiple_of(GIB) {
        format!("{} GiB", bytes / GIB)
    } else if bytes >= MIB && bytes.is_multiple_of(MIB) {
        format!("{} MiB", bytes / MIB)
    } else if bytes >= KIB && bytes.is_multiple_of(KIB) {
        format!("{} KiB", bytes / KIB)
    } else {
        format!("{bytes} B")
    }
}
