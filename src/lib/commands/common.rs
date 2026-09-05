//! Common CLI options shared across commands.
//!
//! This module provides shared argument structures that can be composed into
//! command structs using `#[command(flatten)]`.

use std::collections::HashMap;
use std::ffi::OsString;
use std::path::{Path, PathBuf};
#[cfg(feature = "simplex")]
use std::sync::Arc;

use crate::assigner::Strategy;
#[cfg(feature = "simplex")]
use crate::logging::OperationTimer;
use crate::pipeline::backpressure::{
    BACKPRESSURE_THRESHOLD_BYTES, Q5_BACKPRESSURE_THRESHOLD_BYTES, stage_high_water_mark,
};
use crate::scheduler_strategy::SchedulerStrategy;
use crate::unified_pipeline::BamPipelineConfig;
use crate::validation::validate_input_exists;
use bytesize::ByteSize;
use clap::Args;
use fgumi_bam_io::is_stdout_path;
use fgumi_consensus::methylation::RefBaseProvider;
use fgumi_umi::IndexThreshold;
#[cfg(feature = "simplex")]
use log::{info, warn};
use noodles::sam::Header;

/// CLI argument value for `--methylation-mode`.
///
/// Maps to [`fgumi_consensus::MethylationMode`] variants (excluding `Disabled`,
/// which is represented by the absence of the flag).
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum MethylationModeArg {
    /// EM-Seq (enzymatic methyl-seq): unmethylated C is converted to T.
    /// C in read at ref-C = methylated (protected); T = unmethylated (converted).
    #[value(name = "em-seq")]
    EmSeq,
    /// TAPs/Illumina 5-base: methylated C is converted to T.
    /// C in read at ref-C = unmethylated (not a target); T = methylated (converted).
    #[value(name = "taps")]
    Taps,
}

impl From<MethylationModeArg> for fgumi_consensus::MethylationMode {
    fn from(arg: MethylationModeArg) -> Self {
        match arg {
            MethylationModeArg::EmSeq => Self::EmSeq,
            MethylationModeArg::Taps => Self::Taps,
        }
    }
}

/// CLI argument value for `--tie-rule`.
///
/// Maps to [`fgumi_consensus::TieRule`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, clap::ValueEnum)]
pub enum TieRuleArg {
    /// No-call a maximum separated from the runner-up by only a few ULPs.
    ///
    /// A separation that small is summation-order noise rather than evidence, and the rule is
    /// independent of both read order and family depth — the better rule on the merits, but
    /// not the fgbio-compatible one, which is why it is not the default.
    #[value(name = "ulp-relative")]
    UlpRelative,
    /// Reproduce fgbio's tie rule exactly, defects included.
    ///
    /// The default: fgumi is a drop-in replacement for fgbio, so matching its output is the
    /// contract. For byte-parity with fgbio output — chiefly cross-tool equivalency testing.
    /// Inherits fgbio's order dependence and scale dependence, so it will call a base off one
    /// ULP of accumulation noise.
    #[default]
    #[value(name = "fgbio-compat")]
    FgbioCompat,
}

impl From<TieRuleArg> for fgumi_consensus::TieRule {
    fn from(arg: TieRuleArg) -> Self {
        match arg {
            TieRuleArg::UlpRelative => Self::UlpRelative,
            TieRuleArg::FgbioCompat => Self::FgbioCompat,
        }
    }
}

impl From<fgumi_consensus::TieRule> for TieRuleArg {
    /// Reverse of the `TieRuleArg → TieRule` resolution. The two enums are a
    /// 1:1 pairing, so the `Options` projections that store the resolved
    /// [`fgumi_consensus::TieRule`] can reconstruct the CLI-facing
    /// [`ConsensusCallingOptions`] without losing information.
    fn from(rule: fgumi_consensus::TieRule) -> Self {
        match rule {
            fgumi_consensus::TieRule::UlpRelative => Self::UlpRelative,
            fgumi_consensus::TieRule::FgbioCompat => Self::FgbioCompat,
        }
    }
}

/// Resolves an optional `--methylation-mode` CLI arg to a [`MethylationMode`].
///
/// Returns `Disabled` when `None` (flag not provided).
///
/// [`MethylationMode`]: fgumi_consensus::MethylationMode
pub fn resolve_methylation_mode(
    arg: Option<MethylationModeArg>,
) -> fgumi_consensus::MethylationMode {
    arg.map_or(fgumi_consensus::MethylationMode::Disabled, Into::into)
}

/// Methylation reference pair: reference base provider + contig name mapping.
#[cfg(feature = "simplex")]
pub type MethylationRef = Option<(
    Arc<dyn fgumi_consensus::methylation::RefBaseProvider + Send + Sync>,
    Arc<Vec<String>>,
)>;

/// Loads the reference FASTA and builds contig name mapping for methylation-aware modes.
///
/// Returns `None` if methylation mode is disabled. Errors if enabled but `reference` is `None`.
#[cfg(feature = "simplex")]
pub fn load_methylation_reference(
    methylation_mode: fgumi_consensus::MethylationMode,
    reference: &Option<PathBuf>,
    header: &Header,
) -> anyhow::Result<MethylationRef> {
    if !methylation_mode.is_enabled() {
        return Ok(None);
    }
    let mode_name = match methylation_mode {
        fgumi_consensus::MethylationMode::EmSeq => "EM-Seq",
        fgumi_consensus::MethylationMode::Taps => "TAPs",
        fgumi_consensus::MethylationMode::Disabled => unreachable!(),
    };
    let ref_path = reference
        .as_ref()
        .ok_or_else(|| anyhow::anyhow!("--ref is required when --methylation-mode is set"))?;
    let ref_timer = OperationTimer::new("Loading reference FASTA");
    let reference = Arc::new(crate::reference::ReferenceReader::new(ref_path)?);
    ref_timer.log_completion(0);

    let ref_names: Vec<String> =
        header.reference_sequences().keys().map(|name| name.to_string()).collect();

    // Fail fast if any BAM header contigs are missing from the reference FASTA
    let missing_contigs: Vec<&String> =
        ref_names.iter().filter(|name| reference.sequence_for(name).is_none()).collect();
    if !missing_contigs.is_empty() {
        anyhow::bail!(
            "Reference FASTA is missing {} contig(s) from the BAM header: {}",
            missing_contigs.len(),
            missing_contigs.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(", ")
        );
    }

    info!("{mode_name} mode enabled with {} reference contigs", ref_names.len());
    Ok(Some((reference, Arc::new(ref_names))))
}

/// Add a @PG record to an existing header, using the current fgumi version.
///
/// Wraps [`fgumi_bam_io::header::add_pg_record`] with the binary's version string.
pub fn add_pg_record(header: Header, command_line: &str) -> anyhow::Result<Header> {
    fgumi_bam_io::header::add_pg_record(header, crate::version::VERSION.as_str(), command_line)
}

/// Ensure the header carries an `@HD` line, synthesizing `@HD VN:1.6 SO:unsorted`
/// when it is absent (matching fgbio).
///
/// Wraps [`fgumi_bam_io::header::ensure_hd_record`].
pub fn ensure_hd_record(header: Header) -> anyhow::Result<Header> {
    fgumi_bam_io::header::ensure_hd_record(header)
}

/// Add a @PG record to a header builder, using the current fgumi version.
///
/// Wraps [`fgumi_bam_io::header::add_pg_to_builder`] with the binary's version string.
pub fn add_pg_to_builder(
    builder: noodles::sam::header::Builder,
    command_line: &str,
) -> anyhow::Result<noodles::sam::header::Builder> {
    fgumi_bam_io::header::add_pg_to_builder(builder, crate::version::VERSION.as_str(), command_line)
}

/// Builds the `Index threshold:` startup line for a UMI-assignment command.
///
/// `group` and `dedup` both report the N-gram index threshold that is actually in
/// effect, which is not always the raw `--index-threshold` value:
///
/// - [`Strategy::Edit`] indexes only at one mismatch, and floors a numeric flag at its
///   own measured crossover ([`fgumi_umi::EDIT_INDEX_THRESHOLD`], higher than the shared
///   default), so it reports `max(flag, EDIT_INDEX_THRESHOLD)` at one mismatch and
///   "not used" otherwise.
/// - [`Strategy::Adjacency`] also indexes only at one mismatch, so it reports "not used"
///   at any other, and the flag verbatim at one.
/// - [`Strategy::Paired`] uses the flag verbatim at every mismatch count.
/// - [`Strategy::Identity`] never indexes, so it reports nothing.
///
/// # Arguments
///
/// * `strategy` - the strategy actually in effect (after any `--no-umi` override)
/// * `effective_edits` - the mismatch count actually in effect
/// * `index_threshold` - the raw `--index-threshold` flag value
///
/// # Returns
///
/// The line to log, or `None` when the strategy reports no threshold.
pub fn index_threshold_log_message(
    strategy: Strategy,
    effective_edits: u32,
    index_threshold: IndexThreshold,
) -> Option<String> {
    match strategy {
        Strategy::Edit if effective_edits != 1 => {
            Some("Index threshold: not used (edit indexes only at --edits 1)".to_string())
        }
        Strategy::Edit => Some(format!(
            "Index threshold: {} (edit)",
            index_threshold.floored_at(fgumi_umi::EDIT_INDEX_THRESHOLD)
        )),
        Strategy::Adjacency if effective_edits != 1 => {
            Some("Index threshold: not used (adjacency indexes only at --edits 1)".to_string())
        }
        Strategy::Adjacency | Strategy::Paired => {
            Some(format!("Index threshold: {index_threshold}"))
        }
        Strategy::Identity => None,
    }
}

/// The `--index-threshold` doc-comments on `MarkDuplicates` (`dedup.rs`) and
/// `GroupReadsByUmi` (`group.rs`) both recite edit's floor as a literal `200`,
/// because clap's `help` must be a string literal and cannot interpolate a
/// constant. Breaking the build here is what keeps those two help texts honest
/// when the constant moves.
const _: () = assert!(
    fgumi_umi::EDIT_INDEX_THRESHOLD == 200,
    "EDIT_INDEX_THRESHOLD changed: update the `--index-threshold` help text in \
     src/lib/commands/dedup.rs and src/lib/commands/group.rs, which hardcode 200."
);

/// Logs the `Index threshold:` startup line built by [`index_threshold_log_message`].
///
/// Shared by `group` and `dedup` so the two cannot drift apart when a threshold
/// changes. See [`index_threshold_log_message`] for the arguments and the
/// per-strategy behavior.
pub fn log_index_threshold(
    strategy: Strategy,
    effective_edits: u32,
    index_threshold: IndexThreshold,
) {
    if let Some(message) = index_threshold_log_message(strategy, effective_edits, index_threshold) {
        log::info!("{message}");
    }
}

/// Verifies that a consensus-calling input BAM is sorted appropriately.
///
/// Mirrors fgbio's `UmiConsensusCaller.checkSortOrder`
/// (`UmiConsensusCaller.scala:165-173`), which the `CallMolecularConsensusReads`,
/// `CallDuplexConsensusReads`, and `CallCodecConsensusReads` tools all invoke on
/// their input header. Consensus calling groups reads by molecule as they stream,
/// so it requires template-coordinate order; coordinate-sorted (or otherwise
/// mis-ordered) input silently interleaves molecules and every molecule is split,
/// corrupting the entire output.
///
/// Behavior, matching fgbio exactly:
/// - Template-coordinate (`SO:unsorted`, `GO:query`, `SS:template-coordinate`):
///   accepted silently.
/// - Query-grouped but unsorted without the `SS` sub-sort (`SO:unsorted`,
///   `GO:query`): accepted with a warning (probably compatible, e.g. an fgbio
///   `GroupReadsByUmi` output that omits the sub-sort).
/// - Anything else (coordinate, queryname, ungrouped, …): rejected with an error.
///
/// This reads only the header and never touches or re-opens the record stream, so
/// it is safe to call on a stdin-backed reader after the header has been read.
///
/// # Arguments
///
/// * `header` - the already-read input BAM header
/// * `source` - a human-readable description of the input (path or `<stdin>`)
///
/// # Errors
///
/// Returns an error if the header advertises an order that is definitely
/// incompatible with consensus calling.
pub fn check_consensus_sort_order(header: &Header, source: &str) -> anyhow::Result<()> {
    if fgumi_sam::is_template_coordinate_sorted(header) {
        return Ok(());
    }
    if fgumi_sam::is_query_grouped_unsorted(header) {
        warn!(
            "File {source} may not be sorted correctly for consensus read generation \
             (query-grouped but missing the template-coordinate sub-sort). Continuing, \
             but sort with `fgumi sort --order template-coordinate` if the output looks wrong."
        );
        return Ok(());
    }
    anyhow::bail!(
        "File {source} is not sorted correctly for consensus calling. The input must be \
         template-coordinate sorted (header must advertise SO:unsorted, GO:query, and \
         SS:template-coordinate).\n\nSort it first, e.g.:\n  \
         fgumi sort -i input.bam -o sorted.bam --order template-coordinate"
    );
}

/// Requires that the input header is queryname sorted or query grouped, matching
/// fgbio's `Bams.requireQueryGrouped` (used by `FilterConsensusReads` and `ClipBam`).
///
/// Template-based commands group a template's reads by *adjacency*. On coordinate-
/// sorted input the mates scatter, every read degrades to its own single-read
/// "template", and pair-level logic (overlap clipping, mate-info fix,
/// both-primaries-pass) is silently wrong with a success exit. fgbio hard-fails
/// here; fgumi must too.
///
/// This is the query-grouped guard, weaker than [`check_consensus_sort_order`]:
/// a plain `SO:queryname` file (no template-coordinate sub-sort) is accepted here
/// but rejected there. Use this for filter/clip; use `check_consensus_sort_order`
/// for consensus callers.
///
/// Reads only the header, so it is safe to call on a stdin-backed reader after the
/// header has been read.
///
/// # Arguments
///
/// * `header` - the already-read input BAM header
/// * `source` - a human-readable description of the input (path or `<stdin>`)
///
/// # Errors
///
/// Returns an error if the header is neither queryname sorted (`SO:queryname`) nor
/// query grouped (`GO:query`).
pub fn require_query_grouped(header: &Header, source: &str) -> anyhow::Result<()> {
    if fgumi_sam::is_query_grouped(header) {
        return Ok(());
    }
    let (so, go, ss) = header_sort_and_group_order(header);
    // Mirror fgbio's requireQueryGrouped: append " SS:{ss}" only when present.
    let ss_suffix = ss.map_or_else(String::new, |ss| format!(" SS:{ss}"));
    anyhow::bail!(
        "File {source} was not queryname sorted or query grouped, found: SO:{so} GO:{go}{ss_suffix}. \
         A template's reads must be adjacent, so the input must advertise SO:queryname \
         or GO:query.\n\nSort it first, e.g.:\n  \
         fgumi sort -i input.bam -o sorted.bam --order queryname"
    );
}

/// Returns the header's declared sort order (`SO`, default `unsorted` per htsjdk),
/// group order (`GO`, default `none`), and sub-sort (`SS`, `None` when absent) as
/// display strings for diagnostics — mirroring fgbio's `requireQueryGrouped`.
fn header_sort_and_group_order(header: &Header) -> (String, String, Option<String>) {
    let Some(hdr_map) = header.header() else {
        return ("unsorted".to_string(), "none".to_string(), None);
    };
    let other = hdr_map.other_fields();
    let read = |tag: &[u8; 2]| -> Option<String> {
        other.get(tag).map(|v| String::from_utf8_lossy(<_ as AsRef<[u8]>>::as_ref(v)).into_owned())
    };
    (
        read(b"SO").unwrap_or_else(|| "unsorted".to_string()),
        read(b"GO").unwrap_or_else(|| "none".to_string()),
        read(b"SS"),
    )
}

/// The fgbio `ConsensusCallingIterator` pre-group filter
/// (`ConsensusCallingIterator.scala:56-58`): drop secondary/supplementary
/// alignments and any record that is unmapped.
///
/// fgbio applies this to *every* consensus caller (simplex, duplex, codec),
/// for both the single- and multi-threaded execution paths, before grouping
/// reads by molecular identifier. fgumi replicates it in the consensus
/// commands.
///
/// An unmapped record is dropped even when its mate is mapped. It has no cigar,
/// and an empty cigar is a prefix of every cigar, so it would match every
/// alignment group in `select_most_common_alignment_group` rather than forming
/// one of its own, and would contribute bases and depth to the consensus with no
/// alignment supporting it. The *mapped* end of a half-mapped pair is still kept.
/// See fulcrumgenomics/fgbio#1168.
///
/// `--allow-unmapped` relaxes **only** the mapped-record rule, admitting
/// unmapped primary records into grouping (mirroring `fgumi group
/// --allow-unmapped`). Whether they then produce a consensus is the *caller's*
/// decision, and the callers differ: `VanillaUmiConsensusCaller` (simplex and
/// duplex) calls fully-unmapped input, while `CodecConsensusCaller` requires a
/// mapped primary FR pair and rejects anything else as
/// `CallerRejectionReason::NotPrimaryFrPair`, so `codec` emits no consensus for
/// an unmapped pair however the flag is set. Where unmapped input *is* called,
/// unmapped reads never displace mapped ones: within any set being consensus
/// called, `VanillaUmiConsensusCaller::drop_unmapped_if_any_mapped` prefers the
/// mapped reads. Secondary/supplementary alignments are **always** dropped,
/// matching fgbio — `--allow-unmapped` never lets non-primary alignments into
/// grouping.
///
/// Returns `true` if the record should be kept.
#[must_use]
pub fn consensus_pregroup_keep_flags(flags: u16, allow_unmapped: bool) -> bool {
    use fgumi_raw_bam::flags;
    // Secondary/supplementary alignments are always excluded, regardless of
    // --allow-unmapped: fgbio never groups non-primary alignments.
    if flags & flags::SECONDARY != 0 || flags & flags::SUPPLEMENTARY != 0 {
        return false;
    }
    // --allow-unmapped relaxes only the mapped-record rule. Without it a primary
    // record is kept only when it is itself mapped; a mapped mate is not an exception.
    if allow_unmapped {
        return true;
    }
    flags & flags::UNMAPPED == 0
}

/// Raw-BAM-bytes wrapper around [`consensus_pregroup_keep_flags`], suitable as
/// a `MiGrouper` / raw-record-iterator record filter.
#[must_use]
pub fn consensus_pregroup_keep_raw(raw: &[u8], allow_unmapped: bool) -> bool {
    consensus_pregroup_keep_flags(fgumi_raw_bam::RawRecordView::new(raw).flags(), allow_unmapped)
}

/// EM-Seq methylation-aware consensus calling options.
#[derive(Debug, Clone, Default, Args)]
pub struct EmSeqOptions {
    /// Enable EM-Seq (enzymatic methyl-seq) methylation-aware consensus calling.
    /// Requires --ref. C→T conversions at reference cytosine positions are treated
    /// as bisulfite/enzymatic conversion, and cu/ct per-base count tags
    /// and MM/ML methylation tags are emitted on consensus reads.
    #[arg(long = "em-seq", default_value_t = false, requires = "reference")]
    pub em_seq: bool,

    /// Path to the reference FASTA file (required when --em-seq is enabled)
    #[arg(long = "ref")]
    pub reference: Option<PathBuf>,
}

/// Common input/output options for commands that read a BAM and write a BAM.
#[derive(Debug, Clone, Args)]
pub struct BamIoOptions {
    /// Input BAM file
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Output BAM file
    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,

    /// Enable async userspace prefetch on the input BAM.
    ///
    /// Spawns a dedicated I/O thread that reads raw bytes into a bounded
    /// queue ahead of the decompression step, so processing threads do
    /// not block on disk. Prototype flag; defaults to off.
    #[arg(long = "async-reader", default_value_t = false, hide = true)]
    pub async_reader: bool,

    /// Verify each BGZF block's CRC32 checksum while decoding the input.
    ///
    /// Without either `--check-crc` or `--no-check-crc`, fgumi verifies for
    /// file input and skips verification for stdin input: a freshly-piped
    /// aligner stream (e.g. `bwa-mem3 ... | fgumi ... -i /dev/stdin`) is
    /// trusted, since any corruption there is a bug in the upstream process
    /// rather than data at rest, while a file may have been archived,
    /// transferred, or copied since it was written, where a flipped bit is
    /// exactly what CRC32 exists to catch. Pass `--check-crc` to force
    /// verification on (e.g. for stdin input you don't trust). Mutually
    /// exclusive with `--no-check-crc`. Honored on every command's input decode,
    /// in both single- and multi-threaded modes: single-threaded decodes route
    /// through fgumi-bgzf, and `--threads N` runs (e.g. `clip --threads N`) decode
    /// through the unified pipeline, which takes its CRC policy from the same
    /// flag. Every run logs a `CRC verify:` line at startup stating what actually
    /// happened.
    #[arg(long = "check-crc", default_value_t = false, conflicts_with = "no_check_crc")]
    pub check_crc: bool,

    /// Skip CRC32 verification while decoding the input.
    ///
    /// Trades the CRC32 integrity check for faster decode. See `--check-crc`
    /// for the default policy this overrides and which commands/modes this
    /// has no effect on. Mutually exclusive with `--check-crc`.
    #[arg(long = "no-check-crc", default_value_t = false, conflicts_with = "check_crc")]
    pub no_check_crc: bool,
}

impl Default for BamIoOptions {
    fn default() -> Self {
        Self {
            input: PathBuf::new(),
            output: PathBuf::new(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
        }
    }
}

impl BamIoOptions {
    /// Construct a `BamIoOptions` from input and output paths. Leaves
    /// opt-in tuning flags (e.g. `async_reader`) and the CRC-verification
    /// override flags at their default (unset) values, so
    /// [`effective_check_crc`](Self::effective_check_crc) falls back to the
    /// file-vs-stdin policy.
    pub fn new(input: impl Into<PathBuf>, output: impl Into<PathBuf>) -> Self {
        Self {
            input: input.into(),
            output: output.into(),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
        }
    }

    /// Resolve the effective CRC-verification policy from `--check-crc` /
    /// `--no-check-crc` and the input's stdin-vs-file status.
    ///
    /// Policy (matches dupblaster): `--check-crc` forces verification on;
    /// `--no-check-crc` forces it off; with neither given, verification
    /// defaults on for file input and off for stdin, since a fresh
    /// aligner-piped stream is trusted while a file may have been archived
    /// or transferred since it was written. `check_crc` and `no_check_crc`
    /// are mutually exclusive at the CLI layer (`conflicts_with`), so at most
    /// one is ever true.
    #[must_use]
    pub fn effective_check_crc(&self) -> bool {
        if self.check_crc {
            true
        } else if self.no_check_crc {
            false
        } else {
            !fgumi_bam_io::is_stdin_path(&self.input)
        }
    }

    /// Log the effective CRC-verification setting at info level, once, at
    /// run start. Makes the `effective_check_crc` policy — a default-behavior
    /// change from fgumi's previous always-verify default — visible in every
    /// run's log rather than a silent decision.
    pub fn log_effective_check_crc(&self) {
        let effective = self.effective_check_crc();
        let reason = if self.check_crc {
            " (--check-crc)"
        } else if self.no_check_crc {
            " (--no-check-crc)"
        } else if fgumi_bam_io::is_stdin_path(&self.input) {
            " (trusted stdin)"
        } else {
            ""
        };
        log::info!("CRC verify: {}{reason}", if effective { "on" } else { "off" });
    }

    /// Build [`fgumi_bam_io::PipelineReaderOpts`] from the async-reader flag
    /// and the [`effective_check_crc`](Self::effective_check_crc) policy.
    pub fn pipeline_reader_opts(&self) -> fgumi_bam_io::PipelineReaderOpts {
        fgumi_bam_io::PipelineReaderOpts {
            async_reader: self.async_reader,
            verify_crc: self.effective_check_crc(),
            // This helper serves non-sort command paths, which do not expose a
            // read-stream knob; keep the plain sequential/async reader.
            read_streams: fgumi_bam_io::ReadStreams::Fixed(1),
        }
    }

    /// Validates that the input file exists (skipped for stdin paths).
    ///
    /// # Errors
    ///
    /// Returns an error if the input file does not exist.
    pub fn validate(&self) -> anyhow::Result<()> {
        validate_input_exists(&self.input, "Input BAM")?;
        Ok(())
    }
}

/// Options for writing rejected reads to a separate file.
#[derive(Debug, Clone, Default, Args)]
pub struct RejectsOptions {
    /// Optional output BAM file for rejected reads
    #[arg(short = 'r', long = "rejects")]
    pub rejects: Option<PathBuf>,
}

impl RejectsOptions {
    /// Returns true if rejects tracking is enabled.
    #[must_use]
    pub fn is_enabled(&self) -> bool {
        self.rejects.is_some()
    }
}

/// Refuse a command whose output paths do not all resolve to distinct
/// destinations.
///
/// Given every output a command will open — `--output`, `--rejects`, `--stats`,
/// `--metrics` (including a prefix's expanded files), histograms — as
/// `(path, flag-label)` pairs, this rejects any two that resolve to the same
/// destination. A command's writers are opened up front and written
/// independently, so pointing two at one destination does not make them take
/// turns — it interleaves or overwrites, corrupting the output (a metrics/stats
/// path equal to `--output` truncates a finished BAM to a TSV). How that fails
/// depends on the destination, and both ways are silent:
///
/// - **stdout.** `-` and `/dev/stdout` both resolve to fd 1, where the two
///   writers share one file description. Every byte lands, but the stream
///   carries two headers and two EOF blocks and no reader can parse it.
/// - **A regular file.** Two `File::create` calls give two file *descriptions*
///   with independent offsets, both starting at zero, so the streams overwrite
///   each other byte for byte — the worse of the two, since bytes are destroyed
///   rather than merely reordered. `fgumi correct -o same.bam --rejects same.bam`
///   exited zero having logged "kept 400 and rejected 0", after which
///   `samtools view` reported `Invalid BGZF header at offset 392`. The secondary
///   writer is created eagerly and emits its header at offset 0 whether or not
///   anything is ever written to it, so this needs no rejected reads to happen.
///
/// Both exit zero today, so the combination is rejected before any writer is
/// opened — `File::create` truncates, and a guard that fired later would already
/// have clobbered the primary output.
///
/// File identity is approximated deliberately. Paths are compared as
/// `(canonicalized parent, file name)`, which catches `out.bam` against
/// `./out.bam` or `sub/../out.bam` but not two symlinks onto one target, two
/// hard links, or one name in two cases on a case-insensitive filesystem. A true
/// `dev`+`ino` comparison needs both files to exist, and by the time they do both
/// `File::create` calls have already truncated; closing that gap means moving the
/// check into the writer layer — a follow-up to #715 (item 1) tracked separately.
///
/// The null device is the one exempt destination: it discards every byte, so two
/// writers on it cannot corrupt each other and `-o /dev/null --rejects /dev/null`
/// is a legitimate compute-only run. Nothing else is exempt for merely not being
/// a regular file — a FIFO or `/dev/fd/N` is a *shared* byte stream, where the
/// two writers append through one file description and interleave two headers,
/// two block sequences, and two EOF markers into a stream no reader can parse.
/// That is the same failure as `-o - --rejects -`, so it is rejected the same way.
///
/// Each target carries the flag label that named it, so the error can point at
/// the offending options. Absent (`None`) outputs are simply left out of the
/// slice by the caller.
///
/// # Errors
///
/// Returns an error if more than one target names stdout, or if two targets
/// resolve to the same non-null destination.
pub fn reject_output_collisions(targets: &[(&Path, &str)]) -> anyhow::Result<()> {
    // stdout is a single shared stream, so at most one target may name it —
    // regardless of spelling (`-`, `/dev/stdout`). Two writers on it would
    // interleave two headers and two block sequences into one unreadable stream.
    let stdout_flags: Vec<&str> =
        targets.iter().filter(|(path, _)| is_stdout_path(path)).map(|&(_, flag)| flag).collect();
    if stdout_flags.len() > 1 {
        anyhow::bail!(
            "{} cannot all write to stdout: the streams would interleave into one unreadable \
             output; give all but one an explicit path",
            join_flags(&stdout_flags)
        );
    }

    // Group the remaining file targets by the identity that can be established
    // before the file exists. The null device is the one destination multiple
    // writers may share (it discards every byte), so it is skipped rather than
    // grouped. The first flag to claim an identity wins the error's "prior" slot.
    let mut seen: HashMap<(PathBuf, Option<OsString>), &str> = HashMap::new();
    for &(path, flag) in targets {
        if is_stdout_path(path) || is_null_device(path) {
            continue;
        }
        let identity = resolve_output_identity_owned(path);
        if let Some(&prior_flag) = seen.get(&identity) {
            anyhow::bail!(
                "{prior_flag} and {flag} both write to {}: two outputs on one destination would \
                 overwrite each other byte for byte, or interleave on a shared stream, and either \
                 way nothing can read the result; give one a different path",
                path.display()
            );
        }
        seen.insert(identity, flag);
    }
    Ok(())
}

/// Join flag labels into an English list for an error message: `"--output and
/// --stats"`, or `"--output, --stats, and --metrics"`.
fn join_flags(flags: &[&str]) -> String {
    match flags {
        [] => String::new(),
        [only] => (*only).to_string(),
        [a, b] => format!("{a} and {b}"),
        [rest @ .., last] => format!("{}, and {last}", rest.join(", ")),
    }
}

/// [`resolve_output_identity`] with the file name owned, so the identity can key
/// a map across targets rather than borrowing each path.
fn resolve_output_identity_owned(path: &Path) -> (PathBuf, Option<OsString>) {
    let (parent, name) = resolve_output_identity(path);
    (parent, name.map(std::ffi::OsStr::to_os_string))
}

/// Whether `path` is the null device, the one destination two writers can share.
///
/// Identity comes from `stat`, not from spelling: a path is the null device when
/// it resolves to the same `(device, inode)` as `/dev/null`, which also accepts
/// `/dev/./null` and a symlink onto it. A path that does not exist, or that
/// cannot be stat'd, is not the null device — the usual case for an output — so
/// the collision check below still runs rather than being skipped on a failed
/// `stat`.
///
/// Off Unix there is no null device to compare against, so nothing is exempt.
#[cfg(unix)]
fn is_null_device(path: &Path) -> bool {
    use std::os::unix::fs::MetadataExt;

    let (Ok(candidate), Ok(null)) = (std::fs::metadata(path), std::fs::metadata(NULL_DEVICE_PATH))
    else {
        return false;
    };
    candidate.dev() == null.dev() && candidate.ino() == null.ino()
}

#[cfg(not(unix))]
fn is_null_device(_path: &Path) -> bool {
    false
}

/// The discard sink two output writers are allowed to share.
#[cfg(unix)]
const NULL_DEVICE_PATH: &str = "/dev/null";

/// A path reduced to the identity that can be established before the file exists.
///
/// [`std::fs::canonicalize`] resolves `.`, `..`, and symlinks, but only for a
/// path that already exists — which an output usually does not. Canonicalizing
/// the *parent* and keeping the file name verbatim resolves as much as can be
/// resolved up front.
///
/// Two fallbacks keep this from inventing failures. `Path::parent` is `Some("")`
/// for a bare file name, which `canonicalize` rejects with `ENOENT`; the intent
/// there is the current directory. And a parent that cannot be canonicalized at
/// all falls back to the path as written rather than raising, because a run whose
/// output directory does not exist should get the writer's message about that,
/// not a collision error from a guard that only meant to compare two names.
fn resolve_output_identity(path: &Path) -> (PathBuf, Option<&std::ffi::OsStr>) {
    let parent = match path.parent() {
        Some(parent) if parent.as_os_str().is_empty() => Path::new("."),
        Some(parent) => parent,
        None => return (path.to_path_buf(), None),
    };
    let resolved = std::fs::canonicalize(parent).unwrap_or_else(|_| parent.to_path_buf());
    (resolved, path.file_name())
}

/// Options for writing statistics to a file.
#[derive(Debug, Clone, Default, Args)]
pub struct StatsOptions {
    /// Optional output file for statistics
    #[arg(short = 's', long = "stats")]
    pub stats: Option<PathBuf>,
}

impl StatsOptions {
    /// Returns true if stats output is enabled.
    #[must_use]
    pub fn is_enabled(&self) -> bool {
        self.stats.is_some()
    }
}

/// Common options for consensus calling (simplex, duplex, codec).
#[derive(Debug, Clone, Args)]
pub struct ConsensusCallingOptions {
    /// Phred-scaled error rate prior to UMI integration
    #[arg(short = '1', long = "error-rate-pre-umi", default_value = "45")]
    pub error_rate_pre_umi: u8,

    /// Phred-scaled error rate post UMI integration
    #[arg(short = '2', long = "error-rate-post-umi", default_value = "40")]
    pub error_rate_post_umi: u8,

    /// Minimum base quality in raw reads to use for consensus
    #[arg(short = 'm', long = "min-input-base-quality", default_value = "10")]
    pub min_input_base_quality: u8,

    /// Produce per-base tags (cd, ce) in addition to per-read tags. The default (true) matches
    /// fgbio, which always writes per-base tags; setting this to false drops per-base tags that
    /// fgbio emits unconditionally.
    #[arg(short = 'B', long = "output-per-base-tags", value_name = "true|false", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub output_per_base_tags: bool,

    /// Quality-trim reads before consensus calling (removes low-quality bases from ends)
    #[arg(long = "trim", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub trim: bool,

    /// Minimum consensus base quality (output consensus bases below this are masked to N). The
    /// default (2) matches fgbio, which hardcodes the consensus base-quality minimum to `MIN_PHRED`
    /// (2) and defers further masking to `fgumi filter` / fgbio `FilterConsensusReads`; exposing
    /// this as a flag is an fgumi superset.
    #[arg(long = "min-consensus-base-quality", default_value = "2")]
    pub min_consensus_base_quality: u8,

    /// How to resolve a near-tie between the two most likely consensus bases.
    ///
    /// `fgbio-compat` (default) reproduces fgbio's rule, calling a base off a one-ULP
    /// separation, because matching fgbio is the contract. `ulp-relative` treats a separation
    /// of a few ULPs as summation-order noise and emits a no-call instead — the better rule on
    /// the merits, and independent of read order and family depth.
    ///
    /// Applies to consensus base calling. Consensus UMI calling (the `RX` tag) always uses the
    /// fgbio-compatible rule.
    ///
    /// Hidden: this is a cross-tool equivalency-testing knob, not a routine analysis option.
    #[arg(long = "tie-rule", default_value = "fgbio-compat", value_enum, hide = true)]
    pub tie_rule: TieRuleArg,
}

impl Default for ConsensusCallingOptions {
    fn default() -> Self {
        Self {
            error_rate_pre_umi: 45,
            error_rate_post_umi: 40,
            min_input_base_quality: 10,
            output_per_base_tags: true,
            trim: false,
            min_consensus_base_quality: 2,
            tie_rule: TieRuleArg::FgbioCompat,
        }
    }
}

impl ConsensusCallingOptions {
    /// Maximum valid Phred score (Illumina 1.8+ uses 0-41, but we allow up to 93).
    ///
    /// 93 is the largest score representable in the SAM `!`..`~` QUAL alphabet;
    /// anything above it would serialize to a byte outside that range.
    pub(crate) const MAX_PHRED: u8 = 93;

    /// Validates the consensus calling options.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Any Phred quality value exceeds `MAX_PHRED` (93)
    /// - `min_consensus_base_quality` is less than 2 (`MIN_PHRED`)
    pub fn validate(&self) -> anyhow::Result<()> {
        use anyhow::bail;

        if self.error_rate_pre_umi > Self::MAX_PHRED {
            bail!(
                "error-rate-pre-umi ({}) exceeds maximum Phred score ({})",
                self.error_rate_pre_umi,
                Self::MAX_PHRED
            );
        }
        if self.error_rate_post_umi > Self::MAX_PHRED {
            bail!(
                "error-rate-post-umi ({}) exceeds maximum Phred score ({})",
                self.error_rate_post_umi,
                Self::MAX_PHRED
            );
        }
        if self.min_input_base_quality > Self::MAX_PHRED {
            bail!(
                "min-input-base-quality ({}) exceeds maximum Phred score ({})",
                self.min_input_base_quality,
                Self::MAX_PHRED
            );
        }
        if self.min_consensus_base_quality < 2 {
            bail!(
                "min-consensus-base-quality ({}) must be at least 2 (MIN_PHRED)",
                self.min_consensus_base_quality
            );
        }
        if self.min_consensus_base_quality > Self::MAX_PHRED {
            bail!(
                "min-consensus-base-quality ({}) exceeds maximum Phred score ({})",
                self.min_consensus_base_quality,
                Self::MAX_PHRED
            );
        }

        Ok(())
    }
}

/// Options for read group and read name prefix in consensus output.
#[derive(Debug, Clone, Args)]
pub struct ReadGroupOptions {
    /// Prefix for consensus read names
    #[arg(short = 'p', long = "read-name-prefix")]
    pub read_name_prefix: Option<String>,

    /// Read group ID for consensus reads
    #[arg(short = 'R', long = "read-group-id", default_value = "A")]
    pub read_group_id: String,
}

impl ReadGroupOptions {
    /// Returns the configured read name prefix, or derives one from the header.
    #[must_use]
    pub fn prefix_or_from_header(&self, header: &noodles::sam::Header) -> String {
        self.read_name_prefix
            .clone()
            .unwrap_or_else(|| crate::consensus_caller::make_prefix_from_header(header))
    }
}

impl Default for ReadGroupOptions {
    fn default() -> Self {
        Self { read_name_prefix: None, read_group_id: "A".to_string() }
    }
}

/// Options for overlapping bases consensus calling.
#[derive(Debug, Clone, Args)]
pub struct OverlappingConsensusOptions {
    /// Consensus call overlapping bases in read pairs before UMI consensus calling
    #[arg(long = "consensus-call-overlapping-bases", value_name = "true|false", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub consensus_call_overlapping_bases: bool,
}

impl Default for OverlappingConsensusOptions {
    fn default() -> Self {
        Self { consensus_call_overlapping_bases: true }
    }
}

impl OverlappingConsensusOptions {
    /// Returns true if overlapping consensus calling is enabled.
    #[must_use]
    pub fn is_enabled(&self) -> bool {
        self.consensus_call_overlapping_bases
    }
}

/// Threading mode for parallel processing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThreadingMode {
    /// Single-threaded mode (no parallelism).
    SingleThreaded,
    /// Thread cap mode: never exceed N total threads.
    /// Suitable for HPC/shared systems where predictable resource usage matters.
    Threads(usize),
}

impl ThreadingMode {
    /// Returns true if this mode enables parallel processing.
    #[must_use]
    pub fn is_parallel(&self) -> bool {
        !matches!(self, Self::SingleThreaded)
    }

    /// Returns the target thread count for this mode.
    #[must_use]
    pub fn num_threads(&self) -> usize {
        match self {
            Self::SingleThreaded => 1,
            Self::Threads(n) => *n,
        }
    }
}

/// Common threading options for parallel BAM processing.
///
/// The `--threads N` option caps total threads at N.
/// Use this on HPC/shared systems where you need predictable resource usage.
///
/// # Examples
///
/// ```bash
/// fgumi group --threads 8 ...
/// # Uses up to 8 threads with work-stealing scheduler
/// ```
#[derive(Debug, Clone, Args)]
pub struct ThreadingOptions {
    /// Number of threads for the multi-threaded pipeline.
    ///
    /// If not specified, uses a single-threaded fast path optimized for
    /// simple streaming. When specified (even with --threads 1), uses the
    /// 7-step parallel pipeline with work-stealing scheduler.
    #[arg(long = "threads")]
    pub threads: Option<usize>,
}

/// Options for output compression.
///
/// Controls BGZF compression level for BAM output files.
#[derive(Debug, Clone, Args)]
pub struct CompressionOptions {
    /// Compression level for output BAM (0-12).
    ///
    /// Level 0 writes uncompressed BGZF (valid BAM, no DEFLATE) — useful for
    /// piped/intermediate outputs where downstream tools will recompress.
    /// Level 1 is fastest DEFLATE with larger files.
    /// Level 12 produces smallest files but is slowest.
    #[arg(long, default_value_t = 1, value_parser = clap::value_parser!(u32).range(0..=12))]
    pub compression_level: u32,
}

impl Default for CompressionOptions {
    /// Mirrors the clap `default_value_t = 1` so programmatic and default-constructed callers
    /// emit level-1 compression rather than the `u32` default of `0`, which would silently
    /// write uncompressed BGZF.
    fn default() -> Self {
        Self { compression_level: 1 }
    }
}

/// Option controlling whether unmapped reads are processed by a consensus caller.
///
/// Shared by the `simplex`, `duplex`, and `codec` commands so the `--allow-unmapped` flag
/// (name, default, and boolean parsing) stays identical across all three.
#[derive(Debug, Clone, Args)]
pub struct AllowUnmappedOptions {
    /// Process unmapped reads. By default (fgbio parity,
    /// `ConsensusCallingIterator.scala:56-58`) an unmapped read is dropped before
    /// consensus calling whether or not its mate is mapped, as are all
    /// secondary/supplementary alignments. Enable for consensus on unmapped input
    /// (e.g. ribosome/protein display), mirroring `fgumi group --allow-unmapped`.
    /// Note that `codec` requires a mapped FR pair, so it emits no consensus for
    /// an unmapped pair even with this flag.
    #[arg(long = "allow-unmapped", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub enabled: bool,
}

/// Options for pipeline scheduler configuration.
///
/// Controls which scheduling strategy is used for thread work assignment
/// in the unified pipeline. Also controls pipeline statistics output.
#[derive(Debug, Clone, Default, Args)]
pub struct SchedulerOptions {
    /// Scheduler strategy for thread work assignment.
    ///
    /// - `chase-bottleneck` (default): Threads dynamically follow work through
    ///   the pipeline, moving downstream when output is blocked and upstream
    ///   when input is empty. Shows ~10% improvement at medium thread counts.
    ///
    /// - `fixed-priority`: Assigns fixed thread roles (reader, writer, workers).
    ///   Thread 0 prioritizes reading, Thread N-1 prioritizes writing.
    #[arg(long = "scheduler", value_enum, default_value_t = SchedulerStrategy::default(), hide = true)]
    pub scheduler: SchedulerStrategy,

    /// Print detailed pipeline statistics at completion.
    ///
    /// Shows per-step timing, throughput, contention metrics, and
    /// per-thread work distribution.
    #[arg(long = "pipeline-stats", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true, hide = true)]
    pub pipeline_stats: bool,

    /// Timeout in seconds for deadlock detection (default: 10, 0 = disabled).
    ///
    /// When no progress is made for this duration, a warning is logged with
    /// diagnostic info (queue depths, memory usage, per-queue timestamps).
    #[arg(long = "deadlock-timeout", default_value_t = 10, hide = true)]
    pub deadlock_timeout: u64,

    /// Enable automatic deadlock recovery (default: false, detection only).
    ///
    /// Uses progressive doubling: 2x -> 4x -> unbind, with restoration
    /// after 30s of sustained progress.
    #[arg(long = "deadlock-recover", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true, hide = true)]
    pub deadlock_recover: bool,
}

impl SchedulerOptions {
    /// Returns the scheduler strategy.
    #[must_use]
    pub fn strategy(&self) -> SchedulerStrategy {
        self.scheduler
    }

    /// Returns true if pipeline stats should be collected and printed.
    #[must_use]
    pub fn collect_stats(&self) -> bool {
        self.pipeline_stats
    }

    /// Returns the deadlock detection timeout in seconds (0 = disabled).
    #[must_use]
    pub fn deadlock_timeout_secs(&self) -> u64 {
        self.deadlock_timeout
    }

    /// Returns true if automatic deadlock recovery is enabled.
    #[must_use]
    pub fn deadlock_recover_enabled(&self) -> bool {
        self.deadlock_recover
    }
}

impl ThreadingOptions {
    /// Default batch size for MI group processing.
    ///
    /// This determines how many MI groups are processed together in a single batch
    /// when using parallel processing. Smaller values reduce latency but increase
    /// synchronization overhead; larger values improve throughput but may cause
    /// uneven work distribution.
    pub const DEFAULT_BATCH_SIZE: usize = 100;

    /// Creates threading options with N threads (uses pipeline).
    #[must_use]
    pub fn new(threads: usize) -> Self {
        Self { threads: Some(threads) }
    }

    /// Creates threading options with no threads specified (uses single-threaded fast path).
    #[must_use]
    pub fn none() -> Self {
        Self { threads: None }
    }

    /// Returns the threading mode based on CLI options.
    ///
    /// - `None` -> `SingleThreaded` (fast path, no pipeline)
    /// - `Some(n)` -> `Threads(n)` (uses pipeline, even when n=1)
    #[must_use]
    pub fn mode(&self) -> ThreadingMode {
        match self.threads {
            None => ThreadingMode::SingleThreaded,
            Some(n) => ThreadingMode::Threads(n),
        }
    }

    /// Returns the number of threads.
    #[must_use]
    pub fn num_threads(&self) -> usize {
        self.mode().num_threads()
    }

    /// Returns true if parallel processing should be used.
    #[must_use]
    pub fn is_parallel(&self) -> bool {
        self.mode().is_parallel()
    }

    /// Returns true if running in single-threaded mode.
    #[must_use]
    pub fn is_single_threaded(&self) -> bool {
        matches!(self.mode(), ThreadingMode::SingleThreaded)
    }

    /// Returns the queue length for parallel processing channels.
    ///
    /// The queue length determines how many batches can be buffered between
    /// the reader and worker threads. A value of `2 * num_threads` provides
    /// good overlap between I/O and compute while limiting memory usage.
    #[must_use]
    pub fn queue_len(&self) -> usize {
        self.num_threads() * 2
    }

    /// Returns a log message describing the threading configuration.
    #[must_use]
    pub fn log_message(&self) -> String {
        match self.mode() {
            ThreadingMode::SingleThreaded => "Single-threaded mode".to_string(),
            ThreadingMode::Threads(n) => format!("Using {n} threads"),
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
// Shared memory-budget types
//
// `--max-memory` / `--memory-reserve` are used both by `fgumi sort` (to bound
// the in-memory sort buffer before spilling to disk) and by every pipeline
// command via [`QueueMemoryOptions`] (to bound the inter-stage pipeline queue).
// The types, parsers, and resolution logic live here so the two surfaces stay
// in lockstep.
//////////////////////////////////////////////////////////////////////////////

/// A memory limit, either auto-detected from the host or a fixed byte count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MemoryLimit {
    /// Detect the (cgroup-aware) host memory and subtract the reserve.
    Auto,
    /// Use a fixed memory limit in bytes.
    Fixed(usize),
}

/// How much memory to reserve for other processes (OS, aligners, etc.) when a
/// memory limit is set to [`MemoryLimit::Auto`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MemoryReserve {
    /// Automatic: `min(10 GiB, 50% of host memory)`.
    Auto,
    /// Reserve a fixed number of bytes.
    Fixed(usize),
}

/// How many spilled runs may accumulate before the oldest are consolidated.
///
/// Mirrors [`MemoryLimit`]: both are host-derived resource budgets, so both are
/// spelled with an explicit `auto` rather than inferred from the flag's absence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaxTempFiles {
    /// Size the limit to the process's open-file budget (`ulimit -n`).
    Auto,
    /// Use a fixed limit.
    Fixed(usize),
}

/// The smallest temp-file limit a merge can act on: it needs at least two
/// inputs to merge.
pub(crate) const MIN_MAX_TEMP_FILES: usize = 2;

/// The minimum per-thread memory budget (256 MiB).
pub(crate) const MIN_MEMORY_PER_THREAD: usize = 256 * 1024 * 1024;

/// Default auto-reserve cap: 10 GiB.
pub(crate) const AUTO_RESERVE_CAP: usize = 10 * 1024 * 1024 * 1024;

/// Fraction of host memory (percent) above which a budget that still fits is
/// warned about anyway. A budget that does not exceed total RAM can still leave
/// too little headroom for a co-resident process — the common case being an
/// aligner holding a resident genome index while its output pipes into this
/// command.
const BUDGET_HIGH_FRACTION_PCT: usize = 80;

/// Whether a resolved `budget` should trigger the co-residency warning: an
/// explicit (`Fixed`) budget that fits host RAM but takes more than
/// [`BUDGET_HIGH_FRACTION_PCT`] of it, leaving little headroom for a co-resident
/// process. `Auto` budgets are excluded — they delegate headroom to
/// `--memory-reserve`, so the warning's "use auto" advice would never apply to a
/// caller already on auto. Extracted as a pure predicate so the threshold and
/// the gate are directly testable (the warning itself has no return-value
/// effect to assert on).
fn budget_crowds_host(is_fixed: bool, budget: usize, total: usize) -> bool {
    is_fixed
        && budget <= total
        && budget.saturating_mul(100) > total.saturating_mul(BUDGET_HIGH_FRACTION_PCT)
}

/// Parse a memory size string into `usize` bytes, suitable for use in clap
/// value parsers.
///
/// Delegates to [`parse_memory_size`] for numeric parsing. Plain numbers are
/// interpreted as MiB (e.g. "768" = 768 MiB). Supports human-readable formats
/// like "2GB", "1GiB", "512MiB". See [`parse_memory_size`] for full details.
fn parse_memory_bytes(s: &str, label: &str) -> Result<usize, String> {
    let bytes = parse_memory_size(s).map_err(|e| e.to_string())?;
    usize::try_from(bytes).map_err(|_| format!("{label} too large: {bytes}"))
}

/// Parse a memory-limit string (e.g. "512M", "1G", "768", "auto").
pub(crate) fn parse_memory(s: &str) -> Result<MemoryLimit, String> {
    let s = s.trim();
    if s.eq_ignore_ascii_case("auto") {
        return Ok(MemoryLimit::Auto);
    }
    Ok(MemoryLimit::Fixed(parse_memory_bytes(s, "Memory size")?))
}

/// Parse a max-temp-files string (e.g. "64", "auto").
///
/// Rejects values below [`MIN_MAX_TEMP_FILES`] at parse time rather than
/// clamping them, so a caller who asks for something the merge cannot do is told
/// so instead of silently getting a different limit.
pub(crate) fn parse_max_temp_files(s: &str) -> Result<MaxTempFiles, String> {
    let s = s.trim();
    if s.eq_ignore_ascii_case("auto") {
        return Ok(MaxTempFiles::Auto);
    }
    let limit: usize = s
        .parse()
        .map_err(|_| format!("invalid value '{s}': expected `auto` or a positive integer"))?;
    if limit < MIN_MAX_TEMP_FILES {
        return Err(format!("invalid value '{s}': must be at least {MIN_MAX_TEMP_FILES}"));
    }
    Ok(MaxTempFiles::Fixed(limit))
}

/// Parse a memory-reserve string (e.g. "10G", "auto").
pub(crate) fn parse_memory_reserve(s: &str) -> Result<MemoryReserve, String> {
    let s = s.trim();
    if s.eq_ignore_ascii_case("auto") {
        return Ok(MemoryReserve::Auto);
    }
    Ok(MemoryReserve::Fixed(parse_memory_bytes(s, "Memory reserve")?))
}

/// Resolve a [`MemoryReserve`] to a concrete byte count given total host memory.
pub(crate) fn resolve_reserve(reserve: MemoryReserve, total_memory: usize) -> usize {
    match reserve {
        MemoryReserve::Fixed(bytes) => bytes,
        // min(10 GiB, 50% of host memory)
        MemoryReserve::Auto => AUTO_RESERVE_CAP.min(total_memory / 2),
    }
}

/// Resolve a memory budget to a concrete byte count.
///
/// For [`MemoryLimit::Auto`]: detects total host memory (cgroup-aware via
/// [`detect_total_memory`]), subtracts the reserve, and—when `per_thread` is
/// set—targets each thread's share with a 256 MiB floor. The result is then
/// **capped to the available (post-reserve) memory**, so the floor can never
/// push the total past what the host has. The reserve makes the budget shrink
/// to fit the host, which is what lets pipeline commands self-throttle instead
/// of OOM-ing.
///
/// For [`MemoryLimit::Fixed`]: multiplies by `threads` when `per_thread` is set;
/// the reserve and host size are ignored.
///
/// The thread count that sizes a sort's per-thread memory budget: `--threads`
/// is the floor, raised to `--sort-threads` when that override is larger, but
/// **0 when `--threads` is 0** so `resolve_memory_budget` still rejects a
/// zero-thread run rather than being clamped up here.
///
/// Single source of truth for both the owned-engine banner path
/// (`Sort::memory_budget_threads`) and the chain builder's `add_sort`
/// (`sort_budget_threads`), so the two cannot drift.
#[must_use]
pub(crate) fn sort_memory_budget_threads(num_threads: usize, sort_threads: Option<usize>) -> usize {
    if num_threads == 0 { 0 } else { num_threads.max(sort_threads.unwrap_or(0)) }
}

/// Calls [`detect_total_memory`] exactly once (it invokes `sysinfo`, which is
/// not free).
pub(crate) fn resolve_memory_budget(
    limit: MemoryLimit,
    reserve: MemoryReserve,
    threads: usize,
    per_thread: bool,
) -> anyhow::Result<usize> {
    // Call once — detect_total_memory() invokes sysinfo, which is not free.
    resolve_memory_budget_with_total(limit, reserve, threads, per_thread, detect_total_memory())
}

/// Pure resolver behind [`resolve_memory_budget`], with `total` (host memory)
/// injected so the `Auto` math is unit-testable on simulated small hosts.
fn resolve_memory_budget_with_total(
    limit: MemoryLimit,
    reserve: MemoryReserve,
    threads: usize,
    per_thread: bool,
    total: usize,
) -> anyhow::Result<usize> {
    if threads == 0 {
        anyhow::bail!("--threads must be at least 1");
    }

    // The co-residency warning below applies only to an explicit budget: `auto`
    // delegates headroom to `--memory-reserve` (so this warning's "use auto"
    // advice would never apply to a caller already on auto), which is why
    // `budget_crowds_host` gates on it.
    let is_fixed = matches!(limit, MemoryLimit::Fixed(_));

    let budget = match limit {
        MemoryLimit::Fixed(bytes) => {
            if per_thread {
                bytes
                    .checked_mul(threads)
                    .ok_or_else(|| anyhow::anyhow!("memory limit × {threads} threads overflowed"))?
            } else {
                bytes
            }
        }
        MemoryLimit::Auto => {
            let margin = resolve_reserve(reserve, total);
            let available = total.saturating_sub(margin);
            // The per-thread floor is a *target*, not a guarantee. On a small
            // host (or high thread count) the floor-based budget can exceed what
            // is actually available; cap it to `available` so `auto` truly
            // self-throttles instead of multiplying the floor past physical
            // memory — the exact OOM this feature exists to prevent (#380).
            let target = if per_thread {
                (available / threads)
                    .max(MIN_MEMORY_PER_THREAD)
                    .checked_mul(threads)
                    .ok_or_else(|| anyhow::anyhow!("auto memory budget overflowed"))?
            } else {
                available.max(MIN_MEMORY_PER_THREAD)
            };
            let budget = target.min(available);
            if budget < target {
                log::warn!(
                    "Auto memory: capping budget to host-available {} ({}/thread target × {} threads \
                     exceeds it after reserve {}); throughput may drop but the run stays within memory",
                    ByteSize(budget as u64),
                    ByteSize(MIN_MEMORY_PER_THREAD as u64),
                    threads,
                    ByteSize(margin as u64),
                );
            }
            // Only the per-thread arm allocates `budget` as `threads` independent slices;
            // otherwise it is one shared pool and `budget / threads` would misdescribe it.
            // Mirrors the framing `log_memory_config` already uses below.
            if per_thread {
                log::debug!(
                    "Auto memory: {} of {} ({}/thread × {} threads, reserve {})",
                    ByteSize(budget as u64),
                    ByteSize(total as u64),
                    ByteSize((budget / threads) as u64),
                    threads,
                    ByteSize(margin as u64),
                );
            } else {
                log::debug!(
                    "Auto memory: {} of {} (shared across {} threads, reserve {})",
                    ByteSize(budget as u64),
                    ByteSize(total as u64),
                    threads,
                    ByteSize(margin as u64),
                );
            }
            budget
        }
    };

    if budget > total {
        log::warn!(
            "Memory budget {} exceeds total host memory {}; this may cause OOM (or, for sort, earlier spill-to-disk)",
            ByteSize(budget as u64),
            ByteSize(total as u64),
        );
    } else if budget_crowds_host(is_fixed, budget, total) {
        // Fits, but only just. In a pipeline the budget is rarely the only
        // resident consumer: an aligner feeding this command holds its genome
        // index (many GiB), so a budget at this fraction of RAM can still OOM
        // even though it is under `total` on its own. The resolver cannot see
        // the co-resident process, so it flags the risk and leaves the choice to
        // the operator rather than silently shrinking a budget they set
        // explicitly.
        let pct = budget.saturating_mul(100) / total.max(1);
        log::warn!(
            "Memory budget {} is {}% of total host memory {}; in a pipeline this may leave too \
             little headroom for co-resident processes (e.g. an aligner holding a genome index) \
             and risk OOM. Consider a smaller --max-memory, or --max-memory auto (which subtracts \
             --memory-reserve to self-throttle).",
            ByteSize(budget as u64),
            pct,
            ByteSize(total as u64),
        );
    }

    Ok(budget)
}

/// Options for pipeline queue memory limits.
///
/// Bounds the memory held in the pipeline's inter-stage queues so a command
/// self-throttles instead of OOM-ing. Shared (via `#[command(flatten)]`) by
/// every multi-threaded command, so the `--max-memory` surface matches
/// `fgumi sort`.
///
/// `--max-memory` is the dominant *controllable* consumer; it is a budget, not
/// a hard RSS cap — a single pathological position group is still processed
/// whole, and each worker has transient working-set memory on top of the queue.
///
/// Both pipelines enforce the budget at the Read step: they stop admitting
/// input once the bytes queued between their stages reach it, so a slow or
/// contended output device backs pressure up to the reader instead of filling
/// every queue to its slot count. See `BamPipelineState::read_admission_allowed`
/// and `FastqPipelineState::read_admission_allowed`.
///
/// Memory held *outside* those aggregates is not covered: a worker's
/// in-progress batch is bounded by the thread count rather than by this budget.
/// The write reorder buffers, by contrast, *are* summed into the aggregates
/// (the FASTQ aggregate counts its write reorder state and the BAM aggregate its
/// input and write reorder states), so this budget bounds them through the Read
/// gate; they additionally apply their own threshold, capped at
/// [`BACKPRESSURE_THRESHOLD_BYTES`].
///
/// The budget is a **total**, and it is not the same thing as what any one stage
/// may hold. Stages back off at their own high-water marks —
/// [`BACKPRESSURE_THRESHOLD_BYTES`]
/// (512 MiB) for the queues and reorder buffers,
/// [`Q5_BACKPRESSURE_THRESHOLD_BYTES`]
/// (256 MiB) for the processed queue. A budget below one of those marks tightens
/// it; a budget above it leaves it alone. So raising `--max-memory` raises how
/// much the pipeline may hold in total, not how much any one stage may hold —
/// see [`stage_high_water_mark`]
/// for why a per-stage trigger should not scale with the total (issue #765).
#[derive(Debug, Clone, Args)]
pub struct QueueMemoryOptions {
    /// Maximum total memory for the pipeline queues.
    ///
    /// Default is "768" (768 MiB) per thread. Pass "auto" to detect host memory
    /// (cgroup-aware) and subtract --memory-reserve, so the budget shrinks to
    /// fit the host and the command self-throttles instead of OOM-ing. Explicit
    /// values like "512M", "1G", "4GiB" are per-thread when --memory-per-thread
    /// is enabled (default); plain numbers are MiB. Mirrors `fgumi sort`'s
    /// --max-memory.
    ///
    /// This bounds the bytes held in the queues *between* pipeline stages, in
    /// total. It is not a hard RSS cap: per-thread working memory and a single
    /// oversized group sit outside it. Individual stages also back off at their
    /// own lower high-water marks (512 MiB, and 256 MiB for the processed
    /// queue); a smaller budget tightens those marks, a larger one raises the
    /// total but leaves them where they are.
    #[arg(long = "max-memory", default_value = "768", value_parser = parse_memory)]
    pub max_memory: MemoryLimit,

    /// Memory to reserve for other processes when --max-memory=auto.
    ///
    /// "auto" (default) reserves min(10 GiB, 50% of host memory). Explicit
    /// values like "10G", "8GiB" set a fixed reservation. Ignored when
    /// --max-memory is an explicit value.
    #[arg(long = "memory-reserve", default_value = "auto", value_parser = parse_memory_reserve)]
    pub memory_reserve: MemoryReserve,

    /// Scale the memory limit by thread count.
    ///
    /// When enabled (default), --max-memory is per thread, so total memory =
    /// `max_memory` × threads. Disable for a fixed total budget regardless of
    /// thread count (recommended on fixed-RAM hosts).
    #[arg(long = "memory-per-thread", value_name = "true|false", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub memory_per_thread: bool,
}

impl Default for QueueMemoryOptions {
    fn default() -> Self {
        Self {
            // 768 MiB per thread — byte-identical to the historical default.
            max_memory: MemoryLimit::Fixed(768 * 1024 * 1024),
            memory_reserve: MemoryReserve::Auto,
            memory_per_thread: true,
        }
    }
}

impl QueueMemoryOptions {
    /// Resolve the total queue memory budget in bytes for `num_threads`.
    ///
    /// Under `--max-memory auto` the budget is detected from (cgroup-aware) host
    /// memory minus `--memory-reserve`, so it shrinks to fit the host. Under an
    /// explicit value it is `max_memory` (× threads when per-thread).
    ///
    /// # Errors
    /// Returns an error if `num_threads` is 0 or the multiplication overflows.
    pub fn calculate_memory_limit(&self, num_threads: usize) -> anyhow::Result<u64> {
        let bytes = resolve_memory_budget(
            self.max_memory,
            self.memory_reserve,
            num_threads,
            self.memory_per_thread,
        )?;
        Ok(bytes as u64)
    }

    /// Renders the resolved memory configuration as a single log line.
    ///
    /// The budget is reported as what it is — a **total** — alongside the
    /// per-stage high-water marks actually in force, so the line cannot be read
    /// as a promise that raising `--max-memory` raises what any one stage holds.
    /// It does not: a budget above a stage's default mark leaves that mark
    /// alone, and only a budget below it pulls it down (issue #765).
    ///
    /// # Arguments
    /// * `num_threads` - Number of threads used for the calculation.
    /// * `total_memory` - The resolved total budget from `calculate_memory_limit`.
    #[must_use]
    pub fn memory_config_summary(&self, num_threads: usize, total_memory: u64) -> String {
        let budget = if self.memory_per_thread && num_threads > 1 {
            format!(
                "{} total ({}/thread × {} threads)",
                ByteSize(total_memory),
                ByteSize(total_memory / num_threads as u64),
                num_threads,
            )
        } else {
            format!("{} total", ByteSize(total_memory))
        };
        format!(
            "Queue memory budget: {budget}; per-stage high-water marks {} (reorder-buffered stages), {} (processed queue)",
            ByteSize(stage_high_water_mark(total_memory, BACKPRESSURE_THRESHOLD_BYTES)),
            ByteSize(stage_high_water_mark(total_memory, Q5_BACKPRESSURE_THRESHOLD_BYTES)),
        )
    }

    /// Logs the resolved memory configuration.
    ///
    /// # Arguments
    /// * `num_threads` - Number of threads used for the calculation.
    /// * `total_memory` - The resolved total budget from `calculate_memory_limit`.
    pub fn log_memory_config(&self, num_threads: usize, total_memory: u64) {
        log::info!("{}", self.memory_config_summary(num_threads, total_memory));
    }
}

/// Frame raw BAM record bytes into the pipeline's serialize buffer.
///
/// Each element of `records` is treated as the body of a BAM record (no
/// `block_size` prefix). This function writes `<u32 LE block_size><body>`
/// for each record, matching BAM's on-disk record framing. The pipeline's
/// downstream BGZF compression stage operates on the framed bytes verbatim.
///
/// Used by command-level `serialize_fn` and `secondary_serialize_fn`
/// closures that produce raw record bytes (e.g. `commands::filter`,
/// the `--rejects` paths in `commands::correct`/`simplex`/`duplex`/`codec`).
///
/// Generic over `R: AsRef<[u8]>` so callers can pass either
/// `&[fgumi_raw_bam::RawRecord]` or `&[Vec<u8>]` without copying.
///
/// # Errors
///
/// Returns [`std::io::ErrorKind::InvalidData`] if:
/// - a record body exceeds `u32::MAX` bytes (BAM records cannot be larger
///   than ~4 GiB by spec), or
/// - the summed framed size of `records` overflows `usize` (only reachable
///   on pathological inputs; surfaces as a hard error instead of wrapping
///   silently in release builds).
///
/// Both error paths are checked in a single up-front validation pass, so
/// `output` is never partially appended on error: either every record is
/// written or none of `output`'s bytes are touched.
pub(crate) fn serialize_raw_bam_records<R: AsRef<[u8]>>(
    records: &[R],
    output: &mut Vec<u8>,
) -> std::io::Result<u64> {
    // Reserve total framed size up front so the per-record extend_from_slice
    // loop doesn't repeatedly grow `output`. Validates each record's body fits
    // a `u32` `block_size` *and* that the summed framed size fits `usize` in
    // this same pass: by the time we start writing, every record is known
    // good, so the write loop cannot fail and leave `output` partially
    // appended. Either no bytes are written on error or all records are.
    let additional = records.iter().try_fold(0usize, |acc, record| {
        let len = record.as_ref().len();
        u32::try_from(len).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("BAM record too large ({len} bytes) for u32 block_size"),
            )
        })?;
        len.checked_add(4).and_then(|frame| acc.checked_add(frame)).ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "serialized BAM batch size overflowed usize",
            )
        })
    })?;
    output.reserve(additional);

    for record in records {
        let body = record.as_ref();
        // Pre-validated in the first pass above; `body.len()` is known to fit
        // a `u32`, so this conversion cannot fail.
        let block_size = u32::try_from(body.len()).expect("body length pre-validated to fit u32");
        output.extend_from_slice(&block_size.to_le_bytes());
        output.extend_from_slice(body);
    }
    Ok(records.len() as u64)
}

/// Checks whether R1 is genomically earlier than R2 from their raw BAM bytes.
///
/// Shared by `group` and `dedup` so their paired-UMI assignment can't drift.
/// Uses zero-allocation CIGAR iteration; unmapped reads return position 0
/// (matching noodles `unwrap_or(0)` behavior). Ordering is by reference id,
/// then unclipped 5' coordinate, then strand.
///
/// The strand tie-break mirrors fgbio `GroupReadsByUmi.umiForRead`'s
/// `pos1 == pos2 && r1.positiveStrand`: when both mates share an unclipped 5'
/// coordinate (fully-overlapping / short-insert pairs), R1 is "earlier" iff it
/// is on the forward strand. A bare `r1_pos <= r2_pos` returns `true`
/// unconditionally on a tie, which assigns the lower/higher paired-UMI prefix
/// inconsistently between the two duplex strands — so the strands fail to
/// reverse-match and one duplex molecule is incorrectly split into two.
pub(crate) fn is_r1_genomically_earlier_raw(r1: &[u8], r2: &[u8]) -> bool {
    use fgumi_raw_bam::RawRecordView;

    let ref1 = fgumi_raw_bam::ref_id(r1);
    let ref2 = fgumi_raw_bam::ref_id(r2);
    if ref1 != ref2 {
        return ref1 < ref2;
    }
    let r1_pos = fgumi_raw_bam::unclipped_5prime_from_raw_bam(r1);
    let r2_pos = fgumi_raw_bam::unclipped_5prime_from_raw_bam(r2);
    if r1_pos != r2_pos {
        return r1_pos < r2_pos;
    }
    (RawRecordView::new(r1).flags() & fgumi_raw_bam::flags::REVERSE) == 0
}

// Re-export from the library crate for backward compatibility.
pub(crate) use crate::system::detect_total_memory;
pub use crate::validation::parse_memory_size;

/// Builds a [`BamPipelineConfig`] from the common CLI option structs.
///
/// This consolidates the pipeline configuration boilerplate that is repeated
/// across all multi-threaded commands: auto-tuning, scheduler strategy,
/// stats collection, deadlock settings, queue memory limits, and the
/// `--check-crc`/`--no-check-crc` CRC-verification policy (`io`).
/// Centralizing `verify_crc` here (rather than each command setting
/// `config.pipeline.verify_crc` individually after calling this) means every
/// pipeline-backed command is guaranteed to populate it from
/// [`BamIoOptions::effective_check_crc`] — a command that builds a
/// `PipelineConfig` without going through here can't silently leave
/// `verify_crc` at its default and log a setting it doesn't actually apply.
///
/// After calling this, commands can further customize the returned config
/// (e.g. setting `group_key_config` for raw-byte mode).
pub fn build_pipeline_config(
    scheduler_opts: &SchedulerOptions,
    compression: &CompressionOptions,
    queue_memory: &QueueMemoryOptions,
    io: &BamIoOptions,
    num_threads: usize,
) -> anyhow::Result<BamPipelineConfig> {
    let mut config = BamPipelineConfig::auto_tuned(num_threads, compression.compression_level);
    config.pipeline.scheduler_strategy = scheduler_opts.strategy();
    if scheduler_opts.collect_stats() {
        config.pipeline = config.pipeline.with_stats(true);
    }
    config.pipeline.deadlock_timeout_secs = scheduler_opts.deadlock_timeout_secs();
    config.pipeline.deadlock_recover_enabled = scheduler_opts.deadlock_recover_enabled();
    config.pipeline.verify_crc = io.effective_check_crc();

    let queue_memory_limit_bytes = queue_memory.calculate_memory_limit(num_threads)?;
    config.pipeline.queue_memory_limit = queue_memory_limit_bytes;
    queue_memory.log_memory_config(num_threads, queue_memory_limit_bytes);
    Ok(config)
}

/// Reject an `--index-threshold` that demands indexing under a configuration that can
/// never index.
///
/// `always` is an assertion that the UMI index will be used, so a strategy/edits
/// combination with no index code path has to be an error rather than a flag that is
/// quietly ignored. Pass the EFFECTIVE strategy and edits, so `--no-umi` (which forces
/// identity) is caught too.
///
/// A bare integer is deliberately not checked. It is a tuning knob allowed to end up
/// inert — otherwise the default `100` could not coexist with `--strategy identity` —
/// and that includes `0`, even though `0` admits every group exactly as `always` does.
///
/// # Errors
///
/// Returns an error naming the conflict and how to resolve it.
pub fn validate_index_threshold(
    index_threshold: IndexThreshold,
    strategy: Strategy,
    edits: u32,
) -> anyhow::Result<()> {
    if !index_threshold.demands_indexing() || strategy.can_use_index(edits) {
        return Ok(());
    }

    let name = format!("{strategy:?}").to_lowercase();
    let (context, reason) = match strategy {
        // Edit and adjacency both have an index, just not at this edit distance, so the
        // reason has to name the distance rather than claim they never index.
        Strategy::Edit | Strategy::Adjacency => {
            (format!(" --edits {edits}"), format!("the {name} strategy only indexes at --edits 1"))
        }
        _ => (String::new(), format!("the {name} strategy never uses the UMI index")),
    };

    anyhow::bail!(
        "--index-threshold always cannot be honoured with --strategy {name}{context}: \
         {reason}. Drop --index-threshold to leave indexing to the default threshold, \
         or pass --index-threshold never to state that a linear scan is intended."
    )
}

// ==== ported from feat-runall for the chain builder (R2) ====
/// Log warnings for `SchedulerOptions` / `QueueMemoryOptions` flags that the
/// typed-step pipeline doesn't honor. Called from the chain builder's
/// `add_*` stage methods (which back every command's multi-threaded chain
/// dispatch) so users see why their flags might appear to be ignored.
/// `--pipeline-stats` is honored separately via `attach_new_pipeline_stats`;
/// this only warns about the others.
pub(crate) fn warn_unwired_pipeline_flags(scheduler_opts: &SchedulerOptions) {
    // --scheduler selects a legacy unified-pipeline scheduler strategy that the
    // typed-step chain engine does not consume, so setting it (a hidden dev flag)
    // has no effect on any chain-backed command. Mirror the --deadlock-recover
    // diagnostic below: surface it when set to a non-default value so a developer
    // sees the flag was inert rather than silently ignored. Use `warn!` (not
    // `info!`) to match the sibling "flag has no effect on the chain path"
    // diagnostics in `group::execute_chain` (--debug-memory, FGUMI_SHORT_CIRCUIT).
    let requested_scheduler = scheduler_opts.strategy();
    if requested_scheduler != SchedulerStrategy::default() {
        log::warn!(
            "--scheduler has no effect in the typed-step pipeline: the chain engine \
             does not use a pluggable scheduler strategy, so the requested strategy \
             ({requested_scheduler:?}) is ignored"
        );
    }
    // --deadlock-recover: the legacy progressive-doubling recovery
    // addressed a failure mode (a worker pinned on a stuck step under
    // legacy's static scheduler) that doesn't exist in the typed-step
    // dispatch model — every worker round-robins through every step
    // each iteration, so there's nothing to "recover" from.
    if scheduler_opts.deadlock_recover_enabled() {
        log::warn!(
            "--deadlock-recover has no effect in the typed-step pipeline: \
             the dispatch model round-robins all workers across all steps, \
             so the failure mode legacy's progressive recovery addressed \
             does not occur"
        );
    }
    // `--queue-memory` needs no warning here: it is honored by the chain
    // engine — the total bytes flow into `PipelineConfig::queue_memory_total`,
    // which seeds the initial per-queue budget AND enables the runtime
    // rebalancer that shifts budget between consistently-full / empty queues.
}

// ==== ported from feat-runall for the chain builder (R2) ====
/// Print the new-pipeline `PipelineStats` snapshot to the log if any
/// were collected. Pairs with `attach_new_pipeline_stats`.
pub(crate) fn log_new_pipeline_stats(
    stats: Option<std::sync::Arc<crate::pipeline::core::runtime::stats::PipelineStats>>,
) {
    if let Some(stats) = stats {
        let snapshot = stats.snapshot();
        log::info!("=== Pipeline statistics ===");
        for line in format!("{snapshot}").lines() {
            log::info!("{line}");
        }
    }
}

/// How an input header satisfies `group`'s record-ordering requirement.
///
/// The three cases are mutually exclusive and are what the diagnostics and the
/// rejection message key off. They must be decided in priority order: a
/// template-coordinate header also advertises `GO:query`, so it satisfies
/// [`crate::sam::is_query_grouped`] too and would otherwise be misreported as
/// merely query-grouped whenever `--allow-unmapped` is set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum InputOrdering {
    /// Template-coordinate sorted: records sharing a position key are adjacent,
    /// so mapped templates group by position as intended.
    TemplateCoordinate,
    /// Query-grouped but not template-coordinate sorted. Accepted only under
    /// `--allow-unmapped`, which groups unmapped reads by UMI alone.
    QueryGroupedUnmapped,
    /// Neither ordering holds, so a template's records may not be adjacent.
    Unusable,
}

/// Classify how an input header orders its records.
///
/// Template-coordinate is the ordering `group` is built around. The streaming
/// `RecordPositionGrouper` emits a group whenever the position key changes
/// between consecutive records, so any input where records sharing a position
/// key aren't already adjacent (queryname-sorted, coordinate-sorted,
/// FASTQ-order, etc.) would be split into many small groups and assigned
/// distinct `MoleculeIds` -- silently wrong. Template-coordinate sort is what
/// makes adjacency match the grouping key.
///
/// `--allow-unmapped` places every unmapped read in a single position group, so
/// the position adjacency that template-coordinate sort exists to provide is not
/// needed: every unmapped record lands in that one group regardless of order.
/// Requiring a template-coordinate sort there would tell a user with unaligned
/// data to sort by coordinates that do not exist, which is not actionable -- it
/// blocked the unaligned-only workflow (ribosome display, for example)
/// outright.
///
/// Query grouping is still required, and that is a correctness requirement
/// rather than a policy choice: templates are assembled by comparing each
/// record's name to the previous one, so if a template's records are not
/// adjacent an R1/R2 pair splits into two partial templates.
///
/// Mapped reads in a query-grouped (rather than template-coordinate) input are
/// not position-adjacent, so each mapped template forms its own position group
/// and therefore its own family. That under-groups rather than merging unrelated
/// molecules, and `--allow-unmapped` already warns that it groups by UMI alone;
/// the caller's warning names the case explicitly so it is not a surprise.
pub(crate) fn classify_input_ordering(header: &Header, allow_unmapped: bool) -> InputOrdering {
    if crate::sam::is_template_coordinate_sorted(header) {
        InputOrdering::TemplateCoordinate
    } else if allow_unmapped && crate::sam::is_query_grouped(header) {
        InputOrdering::QueryGroupedUnmapped
    } else {
        InputOrdering::Unusable
    }
}

/// Validate that `header` advertises a record ordering the group stage accepts,
/// emitting the accompanying info/warn logging and bailing with the
/// order-specific remediation hint on rejection.
///
/// This is the single implementation shared by `Group::execute` and the chain
/// builder's `add_group`, so the two orchestrations of the same stage cannot
/// drift on the accepted orders, the error text, or the info/warn logging. It
/// replaces an earlier `require_group_input_sort` whose doc claimed the same
/// sharing but which `Group::execute` never actually called — the two paths had
/// silently diverged (a stricter `SO:queryname` predicate and different
/// messages), which is exactly the failure this consolidation removes.
pub(crate) fn require_group_input_ordering(
    header: &Header,
    allow_unmapped: bool,
) -> anyhow::Result<()> {
    use anyhow::bail;
    use log::{info, warn};

    match classify_input_ordering(header, allow_unmapped) {
        InputOrdering::TemplateCoordinate => {
            info!("Input is template-coordinate sorted");
        }
        InputOrdering::QueryGroupedUnmapped => {
            info!("Input is query-grouped; grouping unmapped reads by UMI only");
            if !header.reference_sequences().is_empty() {
                warn!(
                    "Input declares reference sequences but is not template-coordinate \
                     sorted. Any mapped templates will not be position-grouped, so each \
                     will form its own family. Sort with --order template-coordinate if \
                     the input contains mapped reads you want grouped by position."
                );
            }
        }
        InputOrdering::Unusable => {
            if allow_unmapped {
                bail!(
                    "With --allow-unmapped the input must be either template-coordinate \
                     sorted or query-grouped, so that a template's records are adjacent \
                     (header must advertise SO:queryname or GO:query).\n\n\
                     To sort your BAM file, run:\n  \
                     fgumi sort -i input.bam -o sorted.bam --order queryname"
                );
            }
            bail!(
                "Input BAM must be template-coordinate sorted (header must advertise \
                 SO:unsorted, GO:query, and SS:template-coordinate).\n\n\
                 To sort your BAM file, run:\n  \
                 fgumi sort -i input.bam -o sorted.bam --order template-coordinate"
            );
        }
    }
    Ok(())
}

/// Shared, process-global capturing logger for tests across the crate.
///
/// `log::set_logger` may be called only once per process. `nextest` isolates
/// each test in its own process, but the plain `cargo test` / `cargo t`
/// workflow shares one process across every test, so two modules each
/// installing their own logger panic on the second install. Routing every
/// capture-using test through this one shared, `Once`-guarded installer keeps
/// `cargo t` working. `pub(crate)` so any test module can share it (e.g. the
/// standalone-sort summary test in `pipeline::chains::commands::sort`).
#[cfg(test)]
pub(crate) mod test_log_capture {
    use std::sync::{Mutex, MutexGuard, Once};

    static CAPTURED_LOGS: Mutex<Vec<String>> = Mutex::new(Vec::new());

    /// Serializes whole capture sessions. `CAPTURED_LOGS` is process-global, so
    /// under plain `cargo test` (one process, tests running concurrently) two
    /// capture-using tests would otherwise interleave: one test's `capture_logs`
    /// clear wipes another's records mid-assertion, or its emitted lines pollute
    /// the other's snapshot. Holding this session lock from the clear through the
    /// final `captured()` makes each capture session exclusive.
    ///
    /// Deliberately a *separate* lock from `CAPTURED_LOGS`: `CaptureLogger::log`
    /// takes `CAPTURED_LOGS` on every emitted record, so reusing that buffer's
    /// lock to serialize sessions would deadlock the instant a held session
    /// logged anything.
    static CAPTURE_SESSION: Mutex<()> = Mutex::new(());

    struct CaptureLogger;
    static CAPTURE_LOGGER: CaptureLogger = CaptureLogger;

    impl log::Log for CaptureLogger {
        fn enabled(&self, _metadata: &log::Metadata) -> bool {
            true
        }
        fn log(&self, record: &log::Record) {
            CAPTURED_LOGS
                .lock()
                .expect("captured-log lock poisoned")
                .push(record.args().to_string());
        }
        fn flush(&self) {}
    }

    /// Exclusive log-capture session. Bind it for the whole test and keep it
    /// alive through the final `captured()` call; dropping it releases the
    /// session so another capture-using test may run. The held `MutexGuard` is
    /// never read — it exists purely for its `Drop`.
    #[must_use = "bind the capture session; dropping it immediately ends the exclusive capture window"]
    pub(crate) struct CaptureSession(#[allow(dead_code)] MutexGuard<'static, ()>);

    /// Install the shared capturing logger at most once per process.
    fn install_capture_logger() {
        static INSTALL: Once = Once::new();
        INSTALL.call_once(|| {
            log::set_logger(&CAPTURE_LOGGER).expect("no logger installed yet in this test process");
            log::set_max_level(log::LevelFilter::Trace);
        });
    }

    /// Begin an exclusive capture session and enable the logger (so `log` macros
    /// evaluate their arguments) *without* clearing prior records.
    pub(crate) fn enable_logging() -> CaptureSession {
        // Recover from a poisoned session lock: a prior test that panicked
        // mid-session must not wedge every later capture-using test.
        let guard = CAPTURE_SESSION.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
        install_capture_logger();
        CaptureSession(guard)
    }

    /// Begin an exclusive capture session, install the logger (idempotently),
    /// and clear prior records, so the caller asserts only on what its own
    /// operation emits.
    pub(crate) fn capture_logs() -> CaptureSession {
        let session = enable_logging();
        CAPTURED_LOGS.lock().expect("captured-log lock poisoned").clear();
        session
    }

    /// Snapshot the captured log lines so far. Call while still holding the
    /// `CaptureSession` returned by `enable_logging` / `capture_logs`.
    pub(crate) fn captured() -> Vec<String> {
        CAPTURED_LOGS.lock().expect("captured-log lock poisoned").clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `effective_check_crc` truth table (dupblaster policy): an explicit flag
    /// always wins; with neither flag given, the policy falls back to
    /// file-vs-stdin (verify for files, skip for stdin).
    #[rstest::rstest]
    #[case::check_crc_flag_forces_on_even_for_stdin(true, false, "-", true)]
    #[case::no_check_crc_flag_forces_off_even_for_a_file(false, true, "input.bam", false)]
    #[case::neither_flag_file_input_defaults_on(false, false, "input.bam", true)]
    #[case::neither_flag_stdin_defaults_off(false, false, "-", false)]
    #[case::neither_flag_dev_stdin_defaults_off(false, false, "/dev/stdin", false)]
    fn effective_check_crc_truth_table(
        #[case] check_crc: bool,
        #[case] no_check_crc: bool,
        #[case] input: &str,
        #[case] expected: bool,
    ) {
        let io = BamIoOptions {
            input: PathBuf::from(input),
            output: PathBuf::from("output.bam"),
            async_reader: false,
            check_crc,
            no_check_crc,
        };
        assert_eq!(io.effective_check_crc(), expected);
    }

    /// `reject_output_collisions` compares destinations, not the strings naming
    /// them: `out.bam` and `./out.bam` are one file under two names, and two
    /// writers on one file overwrite each other byte for byte.
    ///
    /// Each case names its two paths relative to a scratch directory the test
    /// creates. `sub` exists and `missing` does not, which puts both sides of the
    /// canonicalization fallback under test — the fallback must still catch the
    /// spelling users actually type, and must not turn a run whose output
    /// directory is absent into a collision error, since that run's real problem
    /// is the missing directory and the writer reports it better.
    #[rstest]
    #[case::identical("out.bam", "out.bam", true)]
    #[case::dot_prefixed("out.bam", "./out.bam", true)]
    #[case::through_an_existing_parent("out.bam", "sub/../out.bam", true)]
    #[case::distinct_names("out.bam", "rejects.bam", false)]
    #[case::identical_under_a_missing_parent("missing/out.bam", "missing/out.bam", true)]
    #[case::distinct_under_a_missing_parent("missing/out.bam", "missing/rejects.bam", false)]
    fn reject_output_collisions_compares_resolved_files(
        #[case] output: &str,
        #[case] secondary: &str,
        #[case] collides: bool,
    ) {
        let dir = tempfile::tempdir().expect("create temp dir");
        std::fs::create_dir_all(dir.path().join("sub")).expect("create subdirectory");
        let output = dir.path().join(output);
        let secondary = dir.path().join(secondary);

        let result = reject_output_collisions(&[
            (output.as_path(), "--output"),
            (secondary.as_path(), "--rejects"),
        ]);

        assert_eq!(
            result.is_err(),
            collides,
            "`-o {}` against `--rejects {}`: got {result:?}",
            output.display(),
            secondary.display()
        );
        if let Err(error) = result {
            let message = error.to_string();
            assert!(
                message.contains("--rejects") && message.contains("out.bam"),
                "the error must name the colliding option and the path, got: {message}"
            );
        }
    }

    /// Stdout collides with itself under either spelling, and the null device is
    /// the one destination two writers may share.
    ///
    /// `-` and `/dev/stdout` are the same fd, so the two must be caught in any
    /// combination — comparing the paths as written would let the mixed pairs
    /// through. `/dev/null` is the opposite case: it discards both streams, so
    /// `-o /dev/null --rejects /dev/null` is a legitimate compute-only run that a
    /// name-based check would wrongly reject.
    #[rstest]
    #[case::both_dash("-", "-", true)]
    #[case::both_dev_stdout("/dev/stdout", "/dev/stdout", true)]
    #[case::dash_then_dev_stdout("-", "/dev/stdout", true)]
    #[case::dev_stdout_then_dash("/dev/stdout", "-", true)]
    #[case::stdout_then_a_file("-", "rejects.bam", false)]
    #[case::a_file_then_stdout("out.bam", "-", false)]
    #[case::both_dev_null("/dev/null", "/dev/null", false)]
    #[case::dev_null_then_a_file("/dev/null", "rejects.bam", false)]
    fn reject_output_collisions_handles_stdout_and_the_null_device(
        #[case] output: &str,
        #[case] secondary: &str,
        #[case] collides: bool,
    ) {
        let secondary = PathBuf::from(secondary);

        let result = reject_output_collisions(&[
            (Path::new(output), "--output"),
            (secondary.as_path(), "--rejects"),
        ]);

        assert_eq!(
            result.is_err(),
            collides,
            "`-o {output}` against `--rejects {}`: got {result:?}",
            secondary.display()
        );
        if let Err(error) = result {
            let message = error.to_string();
            assert!(
                message.contains("--rejects") && message.contains("stdout"),
                "the error must name the colliding option and the stream, got: {message}"
            );
        }
    }

    /// A command that was given no secondary output has nothing to collide with.
    ///
    /// The commands call this guard unconditionally, so the `None` path is the
    /// one taken by nearly every real run — including `-o -`, which is the whole
    /// point of the surrounding change and must not be caught by its own guard.
    #[rstest]
    #[case::stdout("-")]
    #[case::dev_stdout("/dev/stdout")]
    #[case::a_file("out.bam")]
    fn reject_output_collisions_allows_an_absent_secondary(#[case] output: &str) {
        let result = reject_output_collisions(&[(Path::new(output), "--output")]);
        assert!(result.is_ok(), "`-o {output}` with no --rejects must be allowed: {result:?}");
    }

    /// A FIFO named by both outputs is a collision, not an exemption.
    ///
    /// The guard used to exempt every destination that merely was not a regular
    /// file, which let `-o pipe --rejects pipe` through. Both writers then append
    /// through one file description, so the FIFO carries two headers, two block
    /// sequences, and two EOF markers — the same unreadable stream that
    /// `-o - --rejects -` produces, and which that pair is already rejected for.
    /// Two *distinct* FIFOs are two streams and stay legal.
    #[cfg(unix)]
    #[rstest]
    #[case::one_fifo_named_twice("primary", "primary", true)]
    #[case::two_distinct_fifos("primary", "secondary", false)]
    fn reject_output_collisions_rejects_a_shared_fifo(
        #[case] output: &str,
        #[case] secondary: &str,
        #[case] collides: bool,
    ) {
        let dir = tempfile::tempdir().expect("create temp dir");
        let output = dir.path().join(output);
        let secondary = dir.path().join(secondary);
        make_fifo(&output);
        if secondary != output {
            make_fifo(&secondary);
        }

        let result = reject_output_collisions(&[
            (output.as_path(), "--output"),
            (secondary.as_path(), "--rejects"),
        ]);

        assert_eq!(
            result.is_err(),
            collides,
            "`-o {}` against `--rejects {}`: got {result:?}",
            output.display(),
            secondary.display()
        );
        if let Err(error) = result {
            let message = error.to_string();
            assert!(
                message.contains("interleave"),
                "two writers on one stream interleave rather than overwrite, got: {message}"
            );
        }
    }

    /// The guard checks every pair among a command's outputs, not just
    /// output-vs-rejects, and a distinct set of paths is always allowed.
    #[rstest]
    #[case::stats_hits_output(&[("out.bam", "--output"), ("r.bam", "--rejects"), ("out.bam", "--stats")], true)]
    #[case::metrics_hits_rejects(&[("out.bam", "--output"), ("m.tsv", "--rejects"), ("m.tsv", "--metrics")], true)]
    #[case::all_distinct(&[("out.bam", "--output"), ("r.bam", "--rejects"), ("s.tsv", "--stats")], false)]
    #[case::single_output(&[("out.bam", "--output")], false)]
    #[case::no_outputs(&[], false)]
    fn reject_output_collisions_checks_every_pair(
        #[case] specs: &[(&str, &str)],
        #[case] collides: bool,
    ) {
        let dir = tempfile::tempdir().expect("create temp dir");
        let owned: Vec<(PathBuf, &str)> =
            specs.iter().map(|(p, flag)| (dir.path().join(p), *flag)).collect();
        let targets: Vec<(&Path, &str)> =
            owned.iter().map(|(p, flag)| (p.as_path(), *flag)).collect();

        let result = reject_output_collisions(&targets);
        assert_eq!(result.is_err(), collides, "{specs:?}: got {result:?}");
    }

    /// stdout is one shared stream, so at most one output may name it — and the
    /// error names every flag that tried to.
    #[test]
    fn reject_output_collisions_rejects_more_than_one_stdout() {
        let result = reject_output_collisions(&[
            (Path::new("-"), "--output"),
            (Path::new("out.bam"), "--rejects"),
            (Path::new("/dev/stdout"), "--stats"),
        ]);
        let message = result.expect_err("two stdout writers must be rejected").to_string();
        assert!(
            message.contains("stdout")
                && message.contains("--output")
                && message.contains("--stats"),
            "the error must name stdout and both offending flags, got: {message}"
        );
    }

    /// Create a FIFO at `path`.
    ///
    /// The standard library has no `mkfifo(3)` binding, and the crate's `nix`
    /// dependency is gated to Linux, so the coreutils front end stands in — it is
    /// present wherever this test compiles.
    #[cfg(unix)]
    fn make_fifo(path: &Path) {
        let status = std::process::Command::new("mkfifo").arg(path).status().expect("run mkfifo");
        assert!(status.success(), "mkfifo {} failed: {status}", path.display());
    }

    /// A default-constructed `CompressionOptions` must agree with the clap
    /// `default_value_t = 1`. The derived `Default` would yield `0`, which is a valid but
    /// *uncompressed* BGZF level, so programmatic callers would silently write uncompressed
    /// BAM while the same command run from the CLI wrote level 1.
    #[test]
    fn compression_options_default_matches_the_cli_default() {
        assert_eq!(CompressionOptions::default().compression_level, 1);
    }

    /// `index_threshold_log_message` reports the threshold that is actually in effect,
    /// which is not always the raw `--index-threshold` value. `Edit` floors a numeric
    /// flag at its own measured crossover and indexes only at one mismatch; `Adjacency`
    /// uses the flag verbatim and also only indexes at one mismatch; `Paired` uses the
    /// flag verbatim always; `Identity` never indexes and so reports nothing.
    #[rstest]
    // Edit floors a numeric flag at EDIT_INDEX_THRESHOLD (200)...
    #[case::edit_flag_below_floor(
        Strategy::Edit,
        1,
        IndexThreshold::MinUmis(100),
        Some("Index threshold: 200 (edit)")
    )]
    // ...a flag above the floor wins, so the user can still raise it...
    #[case::edit_flag_above_floor(
        Strategy::Edit,
        1,
        IndexThreshold::MinUmis(500),
        Some("Index threshold: 500 (edit)")
    )]
    // ...and `0` does not disable indexing (the gate is `distinct >= threshold`), so
    // the floor still applies to it.
    #[case::edit_flag_zero(
        Strategy::Edit,
        1,
        IndexThreshold::MinUmis(0),
        Some("Index threshold: 200 (edit)")
    )]
    // The keywords state an intent rather than a group size, so the floor leaves them be.
    #[case::edit_always(
        Strategy::Edit,
        1,
        IndexThreshold::Always,
        Some("Index threshold: always (edit)")
    )]
    #[case::edit_never(
        Strategy::Edit,
        1,
        IndexThreshold::Never,
        Some("Index threshold: never (edit)")
    )]
    // Edit indexes only at one mismatch, so every other distance reports "not used".
    #[case::edit_two_mismatches(
        Strategy::Edit,
        2,
        IndexThreshold::MinUmis(100),
        Some("Index threshold: not used (edit indexes only at --edits 1)")
    )]
    #[case::edit_zero_mismatches(
        Strategy::Edit,
        0,
        IndexThreshold::MinUmis(100),
        Some("Index threshold: not used (edit indexes only at --edits 1)")
    )]
    // The edit-distance condition outranks the keyword: `never` at two mismatches is
    // "not used" for the stronger reason, not "never". (`always` cannot reach here at
    // all -- `validate_index_threshold` rejects it before anything is logged.)
    #[case::edit_never_two_mismatches(
        Strategy::Edit,
        2,
        IndexThreshold::Never,
        Some("Index threshold: not used (edit indexes only at --edits 1)")
    )]
    // Adjacency reports the raw flag, with no floor, but only indexes at one mismatch.
    #[case::adjacency(
        Strategy::Adjacency,
        1,
        IndexThreshold::MinUmis(100),
        Some("Index threshold: 100")
    )]
    #[case::adjacency_two_mismatches(
        Strategy::Adjacency,
        2,
        IndexThreshold::MinUmis(100),
        Some("Index threshold: not used (adjacency indexes only at --edits 1)")
    )]
    // Paired reports the raw flag at every mismatch count.
    #[case::paired(Strategy::Paired, 1, IndexThreshold::MinUmis(37), Some("Index threshold: 37"))]
    #[case::paired_two_mismatches(
        Strategy::Paired,
        2,
        IndexThreshold::Never,
        Some("Index threshold: never")
    )]
    // Identity never indexes, so it emits no line at all.
    #[case::identity(Strategy::Identity, 0, IndexThreshold::MinUmis(100), None)]
    fn test_index_threshold_log_message(
        #[case] strategy: Strategy,
        #[case] effective_edits: u32,
        #[case] index_threshold: IndexThreshold,
        #[case] expected: Option<&str>,
    ) {
        assert_eq!(
            index_threshold_log_message(strategy, effective_edits, index_threshold).as_deref(),
            expected
        );
    }

    /// The floor the Edit branch applies is the crate constant, not a copy of it — if
    /// `EDIT_INDEX_THRESHOLD` moves, the reported value moves with it.
    #[test]
    fn test_index_threshold_log_message_uses_crate_floor() {
        let floor = fgumi_umi::EDIT_INDEX_THRESHOLD;
        assert_eq!(
            index_threshold_log_message(Strategy::Edit, 1, IndexThreshold::MinUmis(0)),
            Some(format!("Index threshold: {floor} (edit)"))
        );
    }

    /// `log_index_threshold` is the logging wrapper; it must not panic for any strategy,
    /// including the `Identity` case where there is no message to emit.
    #[rstest]
    fn test_log_index_threshold_emits_without_panicking(
        #[values(Strategy::Identity, Strategy::Edit, Strategy::Adjacency, Strategy::Paired)]
        strategy: Strategy,
        #[values(0, 1, 2)] effective_edits: u32,
        #[values(IndexThreshold::Always, IndexThreshold::Never, IndexThreshold::MinUmis(100))]
        index_threshold: IndexThreshold,
    ) {
        let _session = enable_logging();
        log_index_threshold(strategy, effective_edits, index_threshold);
    }

    // Capture-logging helpers are shared crate-wide (see `super::test_log_capture`)
    // so the plain `cargo t` process, which runs every test in one process,
    // installs exactly one global logger instead of panicking on a second install.
    use super::test_log_capture::{capture_logs, captured, enable_logging};

    /// CONS-01: `check_consensus_sort_order` mirrors fgbio's `UmiConsensusCaller.checkSortOrder`
    /// (the cases pinned by `UmiConsensusCallerTest.scala:34-60`): template-coordinate is
    /// accepted silently; query-grouped-unsorted (no `SS`) is accepted with a warning;
    /// coordinate-sorted (silently splits every molecule) and ungrouped headers both error
    /// with a message naming the source and suggesting template-coordinate.
    #[rstest]
    #[case::template_coordinate(
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\tSS:template-coordinate\n",
        true
    )]
    #[case::unsorted_query_warns("@HD\tVN:1.6\tSO:unsorted\tGO:query\n", true)]
    #[case::coordinate("@HD\tVN:1.6\tSO:coordinate\n", false)]
    #[case::ungrouped("@HD\tVN:1.6\n", false)]
    fn test_check_consensus_sort_order(#[case] header_str: &str, #[case] should_accept: bool) {
        let _session = enable_logging();
        let header: Header = header_str.parse().expect("parse");
        let result = check_consensus_sort_order(&header, "foo.bam");
        if should_accept {
            assert!(result.is_ok(), "header must be accepted: {header_str:?}");
        } else {
            let err = result.expect_err("header must be rejected");
            assert!(err.to_string().contains("foo.bam"), "error must name the source: {err}");
            assert!(
                err.to_string().contains("template-coordinate"),
                "error must suggest template-coordinate: {err}"
            );
        }
    }

    /// FILT3-02 / CLIP3-05: `require_query_grouped` mirrors fgbio's
    /// `Bams.requireQueryGrouped` (`isQueryGrouped = SO:queryname || GO:query`).
    /// Error cases assert the found-order echo and the fgbio-style message.
    #[rstest]
    // A plain queryname sort qualifies (unlike the stricter consensus guard).
    #[case::queryname("@HD\tVN:1.6\tSO:queryname\n", true, &[])]
    // GO:query alone qualifies.
    #[case::go_query("@HD\tVN:1.6\tSO:unsorted\tGO:query\n", true, &[])]
    // Template-coordinate consensus output qualifies via GO:query.
    #[case::template_coordinate(
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\tSS:template-coordinate\n",
        true,
        &[]
    )]
    // The FILT3-02/CLIP3-05 footgun: coordinate-sorted input scatters mates.
    #[case::coordinate(
        "@HD\tVN:1.6\tSO:coordinate\n",
        false,
        &["foo.bam", "queryname sorted or query grouped", "SO:coordinate", "GO:none"]
    )]
    // A rejected header with an SS sub-sort echoes it like fgbio (" SS:{ss}").
    #[case::coordinate_with_ss(
        "@HD\tVN:1.6\tSO:coordinate\tSS:coordinate:natural\n",
        false,
        &["SO:coordinate", "GO:none", "SS:coordinate:natural"]
    )]
    // No SO/GO at all (htsjdk default unsorted, GO none) → rejected.
    #[case::bare_unsorted("@HD\tVN:1.6\n", false, &["queryname sorted or query grouped"])]
    fn test_require_query_grouped(
        #[case] header_str: &str,
        #[case] expect_ok: bool,
        #[case] expected_substrings: &[&str],
    ) {
        let header: Header = header_str.parse().expect("parse");
        let result = require_query_grouped(&header, "foo.bam");
        assert_eq!(result.is_ok(), expect_ok, "for header {header_str:?}");
        if let Err(err) = result {
            let msg = err.to_string();
            for sub in expected_substrings {
                assert!(msg.contains(sub), "message {msg:?} missing {sub:?}");
            }
        }
    }

    #[test]
    fn test_none_is_single_threaded() {
        let opts = ThreadingOptions::none();
        assert!(opts.is_single_threaded());
        assert!(!opts.is_parallel());
        assert_eq!(opts.mode(), ThreadingMode::SingleThreaded);
        assert_eq!(opts.threads, None);
    }

    #[test]
    fn test_mode_detection() {
        // None -> SingleThreaded (fast path)
        assert_eq!(ThreadingOptions::none().mode(), ThreadingMode::SingleThreaded);
        // Some(1) -> Threads(1) (uses pipeline, even with 1 thread)
        assert_eq!(ThreadingOptions::new(1).mode(), ThreadingMode::Threads(1));
        // Some(8) -> Threads(8)
        assert_eq!(ThreadingOptions::new(8).mode(), ThreadingMode::Threads(8));
    }

    #[test]
    fn test_num_threads() {
        assert_eq!(ThreadingOptions::none().num_threads(), 1);
        assert_eq!(ThreadingOptions::new(1).num_threads(), 1);
        assert_eq!(ThreadingOptions::new(8).num_threads(), 8);
    }

    #[test]
    fn test_queue_len() {
        assert_eq!(ThreadingOptions::new(1).queue_len(), 2);
        assert_eq!(ThreadingOptions::new(8).queue_len(), 16);
    }

    #[test]
    fn test_log_message() {
        let opts = ThreadingOptions::new(8);
        let msg = opts.log_message();
        assert!(msg.contains("8 threads"));

        let opts = ThreadingOptions::none();
        let msg = opts.log_message();
        assert!(msg.contains("Single-threaded"));
    }

    #[test]
    fn test_new_uses_pipeline() {
        // new(1) should use pipeline (Threads mode), not fast path
        let opts = ThreadingOptions::new(1);
        assert!(!opts.is_single_threaded());
        assert!(opts.is_parallel());
        assert_eq!(opts.threads, Some(1));
    }

    // ========== Tests for option struct validation ==========

    #[test]
    fn test_consensus_calling_options_validate_defaults() {
        let opts = ConsensusCallingOptions::default();
        assert!(opts.validate().is_ok());
    }

    #[test]
    fn test_consensus_calling_options_validate_valid() {
        let opts = ConsensusCallingOptions {
            error_rate_pre_umi: 45,
            error_rate_post_umi: 40,
            min_input_base_quality: 10,
            output_per_base_tags: true,
            trim: false,
            min_consensus_base_quality: 13,
            tie_rule: TieRuleArg::default(),
        };
        assert!(opts.validate().is_ok());
    }

    #[test]
    fn test_consensus_calling_options_validate_error_rate_pre_umi_too_high() {
        let opts = ConsensusCallingOptions {
            error_rate_pre_umi: 94, // Exceeds MAX_PHRED
            ..ConsensusCallingOptions::default()
        };
        let err = opts.validate().unwrap_err();
        assert!(err.to_string().contains("error-rate-pre-umi"));
    }

    #[test]
    fn test_consensus_calling_options_validate_error_rate_post_umi_too_high() {
        let opts = ConsensusCallingOptions {
            error_rate_post_umi: 94, // Exceeds MAX_PHRED
            ..ConsensusCallingOptions::default()
        };
        let err = opts.validate().unwrap_err();
        assert!(err.to_string().contains("error-rate-post-umi"));
    }

    #[test]
    fn test_consensus_calling_options_validate_min_consensus_too_low() {
        let opts = ConsensusCallingOptions {
            min_consensus_base_quality: 1, // Below MIN_PHRED
            ..ConsensusCallingOptions::default()
        };
        let err = opts.validate().unwrap_err();
        assert!(err.to_string().contains("min-consensus-base-quality"));
    }

    #[test]
    fn test_consensus_calling_options_validate_min_consensus_at_min() {
        let opts = ConsensusCallingOptions {
            min_consensus_base_quality: 2, // Exactly MIN_PHRED
            ..ConsensusCallingOptions::default()
        };
        assert!(opts.validate().is_ok());
    }

    // ========== Tests for SchedulerOptions ==========

    #[test]
    fn test_scheduler_options_default() {
        let opts = SchedulerOptions::default();
        assert_eq!(opts.strategy(), SchedulerStrategy::BalancedChaseDrain);
        assert!(!opts.collect_stats());
    }

    #[test]
    fn test_scheduler_options_strategy() {
        let opts = SchedulerOptions {
            scheduler: SchedulerStrategy::FixedPriority,
            pipeline_stats: false,
            deadlock_timeout: 10,
            deadlock_recover: false,
        };
        assert_eq!(opts.strategy(), SchedulerStrategy::FixedPriority);
    }

    #[test]
    fn test_scheduler_options_collect_stats() {
        let opts = SchedulerOptions {
            scheduler: SchedulerStrategy::default(),
            pipeline_stats: true,
            deadlock_timeout: 10,
            deadlock_recover: false,
        };
        assert!(opts.collect_stats());
    }

    #[test]
    fn test_scheduler_options_deadlock_timeout() {
        let opts = SchedulerOptions {
            scheduler: SchedulerStrategy::default(),
            pipeline_stats: false,
            deadlock_timeout: 30,
            deadlock_recover: false,
        };
        assert_eq!(opts.deadlock_timeout_secs(), 30);
    }

    #[test]
    fn test_scheduler_options_deadlock_recover() {
        let opts = SchedulerOptions {
            scheduler: SchedulerStrategy::default(),
            pipeline_stats: false,
            deadlock_timeout: 10,
            deadlock_recover: true,
        };
        assert!(opts.deadlock_recover_enabled());
    }

    // ========== Tests for QueueMemoryOptions ==========

    #[test]
    fn test_queue_memory_default_is_768_mib_per_thread() {
        // The default must stay byte-identical to the historical default
        // (768 MiB per thread).
        let opts = QueueMemoryOptions::default();
        let result = opts
            .calculate_memory_limit(4)
            .expect("calculate_memory_limit should succeed for the default");
        assert_eq!(result, 768 * 1024 * 1024 * 4);
    }

    #[test]
    fn test_queue_memory_max_memory_per_thread_scaling() {
        // 100 MiB × 4 threads = 400 MiB
        let opts = QueueMemoryOptions {
            max_memory: MemoryLimit::Fixed(100 * 1024 * 1024),
            ..QueueMemoryOptions::default()
        };
        let result = opts.calculate_memory_limit(4).expect("should succeed");
        assert_eq!(result, 100 * 1024 * 1024 * 4);
    }

    #[test]
    fn test_queue_memory_fixed_total_does_not_scale() {
        let opts = QueueMemoryOptions {
            max_memory: MemoryLimit::Fixed(200 * 1024 * 1024),
            memory_per_thread: false,
            ..QueueMemoryOptions::default()
        };
        let result = opts.calculate_memory_limit(8).expect("should succeed");
        assert_eq!(result, 200 * 1024 * 1024);
    }

    #[test]
    fn test_queue_memory_zero_threads_is_error() {
        assert!(QueueMemoryOptions::default().calculate_memory_limit(0).is_err());
    }

    #[test]
    fn test_queue_memory_human_readable_via_parse() {
        // "2GB" parses as decimal gigabytes (matching `fgumi sort` / ByteSize).
        let opts = QueueMemoryOptions {
            max_memory: parse_memory("2GB").expect("parse 2GB"),
            memory_per_thread: false,
            ..QueueMemoryOptions::default()
        };
        let result = opts.calculate_memory_limit(4).expect("should succeed");
        assert_eq!(result, 2 * 1000 * 1000 * 1000);
    }

    #[test]
    fn test_queue_memory_auto_is_bounded_by_host() {
        // `auto` self-throttles: the budget must never exceed (cgroup-aware)
        // host memory, regardless of thread count.
        let total = detect_total_memory() as u64;
        let opts =
            QueueMemoryOptions { max_memory: MemoryLimit::Auto, ..QueueMemoryOptions::default() };
        let result = opts.calculate_memory_limit(4).expect("auto should resolve");
        assert!(result > 0, "auto resolved to a zero budget");
        assert!(result <= total, "auto budget {result} exceeded host total {total}");
    }

    #[test]
    fn test_auto_never_oversubscribes_small_host() {
        let _session = enable_logging(); // exercise the cap-warning and auto-debug log branches
        // Simulated 4 GiB host, 16 threads: the 256 MiB/thread floor would want
        // 4 GiB before reserve, which cannot fit after the auto reserve. The
        // budget must be capped to `available`, never `floor × threads`.
        let total = 4 * 1024 * 1024 * 1024; // 4 GiB
        let margin = resolve_reserve(MemoryReserve::Auto, total); // min(10 GiB, 2 GiB) = 2 GiB
        let available = total - margin;
        let budget = resolve_memory_budget_with_total(
            MemoryLimit::Auto,
            MemoryReserve::Auto,
            16,
            true,
            total,
        )
        .expect("should resolve");
        assert!(budget <= available, "budget {budget} oversubscribed available {available}");
        assert!(budget <= total, "budget {budget} oversubscribed host {total}");
    }

    #[test]
    fn test_auto_uses_floor_when_host_is_ample() {
        // Simulated 256 GiB host, 4 threads: plenty of room, so the budget is
        // the per-thread share and stays under available.
        let total = 256 * 1024 * 1024 * 1024;
        let margin = resolve_reserve(MemoryReserve::Auto, total); // 10 GiB cap
        let available = total - margin;
        let budget = resolve_memory_budget_with_total(
            MemoryLimit::Auto,
            MemoryReserve::Auto,
            4,
            true,
            total,
        )
        .expect("should resolve");
        assert!(budget >= MIN_MEMORY_PER_THREAD * 4, "budget {budget} fell below the floor");
        assert!(budget <= available, "budget {budget} exceeded available {available}");
    }

    #[test]
    fn test_fixed_budget_independent_of_host() {
        let _session = enable_logging(); // exercise the "budget exceeds host total" warn branch
        // Fixed limits ignore host size entirely (reserve is irrelevant).
        let tiny_host = 512 * 1024 * 1024;
        let budget = resolve_memory_budget_with_total(
            MemoryLimit::Fixed(2 * 1024 * 1024 * 1024),
            MemoryReserve::Auto,
            4,
            false,
            tiny_host,
        )
        .expect("should resolve");
        assert_eq!(budget, 2 * 1024 * 1024 * 1024);
    }

    #[test]
    fn budget_crowds_host_fires_only_for_tight_fixed_budgets() {
        let g = 1024 * 1024 * 1024;
        // Fixed budget above the 80% line warns; below it does not.
        assert!(budget_crowds_host(true, 14 * g, 16 * g), "87.5% fixed should warn");
        assert!(budget_crowds_host(true, 9 * g, 10 * g), "90% fixed should warn");
        assert!(!budget_crowds_host(true, 12 * g, 16 * g), "75% fixed should not warn");
        // Strict `>`: exactly 80% must NOT warn (off-by-one guard).
        assert!(!budget_crowds_host(true, 8 * g, 10 * g), "exactly 80% must not warn");
        // The `is_fixed` gate: an `auto` budget never takes this branch, however
        // high its fraction of host RAM.
        assert!(!budget_crowds_host(false, 15 * g, 16 * g), "auto must never warn");
        // A budget that exceeds total is the other branch's job, not this one.
        assert!(
            !budget_crowds_host(true, 20 * g, 16 * g),
            "over-total is the exceeds-total branch"
        );
    }

    #[test]
    fn test_high_fraction_fixed_budget_is_returned_unshrunk() {
        // A 14 GiB budget on a 16 GiB host (87.5%) takes the co-residency warn
        // branch (see `budget_crowds_host_fires_only_for_tight_fixed_budgets`
        // for the branch logic). This asserts the complementary facts: the branch
        // (a) emits the co-residency warning and (b) only warns — it does not
        // shrink the explicit budget. `capture_logs` records emitted records so
        // the warning is asserted rather than merely evaluated.
        let _session = capture_logs();
        let host = 16 * 1024 * 1024 * 1024;
        let budget = resolve_memory_budget_with_total(
            MemoryLimit::Fixed(14 * 1024 * 1024 * 1024),
            MemoryReserve::Auto,
            1,
            false,
            host,
        )
        .expect("should resolve");
        assert_eq!(budget, 14 * 1024 * 1024 * 1024, "explicit budget must not be shrunk");

        let logs = captured();
        assert!(
            logs.iter().any(|line| line.contains("of total host memory")
                && line.contains("co-resident processes")),
            "co-residency warning must be emitted; captured logs: {logs:?}"
        );
    }

    #[test]
    fn test_queue_memory_auto_reserve_shrinks_budget() {
        let large_reserve = QueueMemoryOptions {
            max_memory: MemoryLimit::Auto,
            memory_reserve: MemoryReserve::Fixed(512 * 1024 * 1024),
            ..QueueMemoryOptions::default()
        }
        .calculate_memory_limit(4)
        .expect("should succeed");
        let small_reserve = QueueMemoryOptions {
            max_memory: MemoryLimit::Auto,
            memory_reserve: MemoryReserve::Fixed(128 * 1024 * 1024),
            ..QueueMemoryOptions::default()
        }
        .calculate_memory_limit(4)
        .expect("should succeed");
        assert!(large_reserve <= small_reserve);
    }

    #[test]
    fn test_fixed_per_thread_overflow_is_error() {
        let opts = QueueMemoryOptions {
            max_memory: MemoryLimit::Fixed(usize::MAX),
            ..QueueMemoryOptions::default()
        };
        assert!(opts.calculate_memory_limit(4).is_err());
    }

    #[test]
    fn test_auto_per_thread_overflow_is_error() {
        // A pathological thread count makes the per-thread floor × threads
        // overflow; this must surface as an error, not wrap.
        let result = resolve_memory_budget_with_total(
            MemoryLimit::Auto,
            MemoryReserve::Auto,
            usize::MAX,
            true,
            1024 * 1024 * 1024,
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_log_memory_config_exercises_both_branches() {
        // Logging only (no return value); call both paths to keep them covered
        // and to guard against a panic in the formatting/division.
        let per_thread = QueueMemoryOptions::default();
        per_thread.log_memory_config(8, 8 * 768 * 1024 * 1024); // per-thread, threads > 1
        per_thread.log_memory_config(1, 768 * 1024 * 1024); // threads == 1 → total branch

        let fixed_total =
            QueueMemoryOptions { memory_per_thread: false, ..QueueMemoryOptions::default() };
        fixed_total.log_memory_config(8, 1024 * 1024 * 1024); // fixed total branch
    }

    /// The startup line must name the budget as a total *and* report the
    /// per-stage marks in force, so it cannot be read as promising that a
    /// larger `--max-memory` gives any one stage more room (issue #765).
    #[test]
    fn test_memory_config_summary_reports_the_per_stage_marks() {
        let per_thread = QueueMemoryOptions::default();
        let line = per_thread.memory_config_summary(8, 8 * 768 * 1024 * 1024);
        assert!(line.contains("6.0 GiB total"), "{line}");
        assert!(line.contains("768.0 MiB/thread × 8 threads"), "{line}");
        assert!(line.contains("512.0 MiB (reorder-buffered stages)"), "{line}");
        assert!(line.contains("256.0 MiB (processed queue)"), "{line}");
    }

    /// A budget larger than every per-stage mark must still report the marks at
    /// their defaults — that is precisely the fact the old line hid.
    #[test]
    fn test_memory_config_summary_does_not_inflate_the_marks_for_a_huge_budget() {
        let fixed_total =
            QueueMemoryOptions { memory_per_thread: false, ..QueueMemoryOptions::default() };
        let line = fixed_total.memory_config_summary(8, 32 * 1024 * 1024 * 1024);
        assert!(line.contains("32.0 GiB total"), "{line}");
        assert!(line.contains("512.0 MiB (reorder-buffered stages)"), "{line}");
        assert!(line.contains("256.0 MiB (processed queue)"), "{line}");
    }

    /// A budget below the marks pulls them down, and the line says so.
    #[test]
    fn test_memory_config_summary_reports_tightened_marks_for_a_small_budget() {
        let fixed_total =
            QueueMemoryOptions { memory_per_thread: false, ..QueueMemoryOptions::default() };
        let line = fixed_total.memory_config_summary(8, 64 * 1024 * 1024);
        assert!(line.contains("64.0 MiB total"), "{line}");
        assert!(line.contains("64.0 MiB (reorder-buffered stages)"), "{line}");
        assert!(line.contains("64.0 MiB (processed queue)"), "{line}");
    }

    #[test]
    #[cfg(feature = "memory-debug")]
    fn test_sysinfo_returns_reasonable_values() {
        use sysinfo::System;
        let mut system = System::new();
        system.refresh_memory();

        let total = system.total_memory();
        let available = system.available_memory();

        assert!(total > 100_000_000); // > 100MB
        assert!(available > 0);
        assert!(available <= total);
    }

    use clap::Parser;

    /// Test-only wrapper to exercise clap parsing of flattened Args structs.
    #[derive(Debug, Parser)]
    #[command(name = "test")]
    struct TestBoolFlags {
        #[command(flatten)]
        consensus: ConsensusCallingOptions,
        #[command(flatten)]
        overlapping: OverlappingConsensusOptions,
        #[command(flatten)]
        queue_memory: QueueMemoryOptions,
    }

    use rstest::rstest;

    /// The `ConsensusCallingIterator` pre-group filter (fgbio parity,
    /// `ConsensusCallingIterator.scala:56-58`): the default keeps only mapped primary
    /// records. Secondary and supplementary alignments are always dropped, and an
    /// unmapped read is dropped whether or not its mate is mapped — a mapped mate is
    /// not an exception, because the unmapped end has no cigar to group on.
    ///
    /// `--allow-unmapped` (the `allow_unmapped` column) relaxes **only** the
    /// mapped-record rule: it makes primary alignments pass regardless of mapping, but
    /// secondary/supplementary alignments are still dropped. The
    /// `*_dropped_even_with_allow_unmapped` cases pin that the flag never lets a non-primary
    /// alignment into grouping (the regression this exists to prevent).
    #[rstest]
    // allow_unmapped = false (fgbio parity default)
    #[case::mapped_primary_kept(0, false, true)]
    #[case::secondary_dropped(fgumi_raw_bam::flags::SECONDARY, false, false)]
    #[case::supplementary_dropped(fgumi_raw_bam::flags::SUPPLEMENTARY, false, false)]
    #[case::unmapped_unpaired_dropped(fgumi_raw_bam::flags::UNMAPPED, false, false)]
    #[case::unmapped_paired_mate_unmapped_dropped(
        fgumi_raw_bam::flags::UNMAPPED | fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::MATE_UNMAPPED,
        false,
        false
    )]
    // The unmapped end of a half-mapped pair is dropped along with any other unmapped record: it
    // has no cigar to group on. The mapped end is kept, since it is mapped in its own right.
    #[case::unmapped_paired_mate_mapped_dropped(
        fgumi_raw_bam::flags::UNMAPPED | fgumi_raw_bam::flags::PAIRED,
        false,
        false
    )]
    #[case::mapped_paired_mate_unmapped_kept(
        fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::MATE_UNMAPPED,
        false,
        true
    )]
    // allow_unmapped = true relaxes only the unmapped rule
    #[case::mapped_primary_kept_with_allow_unmapped(0, true, true)]
    #[case::unmapped_unpaired_kept_with_allow_unmapped(fgumi_raw_bam::flags::UNMAPPED, true, true)]
    #[case::unmapped_paired_mate_unmapped_kept_with_allow_unmapped(
        fgumi_raw_bam::flags::UNMAPPED | fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::MATE_UNMAPPED,
        true,
        true
    )]
    #[case::secondary_dropped_even_with_allow_unmapped(
        fgumi_raw_bam::flags::SECONDARY,
        true,
        false
    )]
    #[case::supplementary_dropped_even_with_allow_unmapped(
        fgumi_raw_bam::flags::SUPPLEMENTARY,
        true,
        false
    )]
    #[case::unmapped_secondary_dropped_even_with_allow_unmapped(
        fgumi_raw_bam::flags::UNMAPPED | fgumi_raw_bam::flags::SECONDARY,
        true,
        false
    )]
    fn test_consensus_pregroup_keep_flags(
        #[case] flags_bits: u16,
        #[case] allow_unmapped: bool,
        #[case] keep: bool,
    ) {
        assert_eq!(consensus_pregroup_keep_flags(flags_bits, allow_unmapped), keep);
    }

    #[rstest]
    // --output-per-base-tags (default true)
    #[case(&["test"], true)]
    #[case(&["test", "--output-per-base-tags"], true)]
    #[case(&["test", "--output-per-base-tags", "true"], true)]
    #[case(&["test", "--output-per-base-tags", "false"], false)]
    #[case(&["test", "--output-per-base-tags=true"], true)]
    #[case(&["test", "--output-per-base-tags=false"], false)]
    fn test_output_per_base_tags_parsing(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.consensus.output_per_base_tags, expected);
    }

    #[rstest]
    // --trim (default false)
    #[case(&["test"], false)]
    #[case(&["test", "--trim"], true)]
    #[case(&["test", "--trim", "true"], true)]
    #[case(&["test", "--trim", "false"], false)]
    #[case(&["test", "--trim=true"], true)]
    #[case(&["test", "--trim=false"], false)]
    fn test_trim_parsing(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.consensus.trim, expected);
    }

    #[rstest]
    // --consensus-call-overlapping-bases (default true)
    #[case(&["test"], true)]
    #[case(&["test", "--consensus-call-overlapping-bases"], true)]
    #[case(&["test", "--consensus-call-overlapping-bases", "true"], true)]
    #[case(&["test", "--consensus-call-overlapping-bases", "false"], false)]
    #[case(&["test", "--consensus-call-overlapping-bases=true"], true)]
    #[case(&["test", "--consensus-call-overlapping-bases=false"], false)]
    fn test_overlapping_bases_parsing(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.overlapping.consensus_call_overlapping_bases, expected);
    }

    /// The tie-rule default is `fgbio-compat`, and all three declarations of it agree.
    ///
    /// It is declared three times, independently: `#[default]` on the [`TieRuleArg`]
    /// variant, `default_value = "fgbio-compat"` on the `--tie-rule` arg, and the
    /// hand-written [`ConsensusCallingOptions::default`]. Nothing tied them together,
    /// which is how the variant doc came to claim `ulp-relative` was the default long
    /// after #666 switched it — a wrong doc on the flag that decides whether fgumi
    /// reproduces fgbio's base calls.
    ///
    /// Pinning the *value* matters as much as pinning the agreement: `fgbio-compat` is
    /// a compatibility contract, so a silent flip would change consensus output on
    /// near-ties rather than fail anything.
    #[test]
    fn tie_rule_defaults_to_fgbio_compat_everywhere() {
        assert_eq!(
            TieRuleArg::default(),
            TieRuleArg::FgbioCompat,
            "the #[default] variant must stay fgbio-compat",
        );
        assert_eq!(
            TestBoolFlags::try_parse_from(["test"]).expect("no args parses").consensus.tie_rule,
            TieRuleArg::FgbioCompat,
            "clap's default_value must agree with the #[default] variant",
        );
        assert_eq!(
            ConsensusCallingOptions::default().tie_rule,
            TieRuleArg::FgbioCompat,
            "the hand-written Default impl must agree too",
        );
    }

    /// The two `TieRuleArg` variants must map onto the matching
    /// [`fgumi_consensus::TieRule`] variants.
    ///
    /// The CLI enum mirrors the domain enum, so the `From` impl is the seam where a
    /// mirrored pair could be crossed — and crossing them would silently select the
    /// opposite tie rule rather than fail to compile.
    #[rstest]
    #[case::ulp_relative(TieRuleArg::UlpRelative, fgumi_consensus::TieRule::UlpRelative)]
    #[case::fgbio_compat(TieRuleArg::FgbioCompat, fgumi_consensus::TieRule::FgbioCompat)]
    fn tie_rule_arg_maps_onto_the_matching_domain_variant(
        #[case] arg: TieRuleArg,
        #[case] expected: fgumi_consensus::TieRule,
    ) {
        assert_eq!(fgumi_consensus::TieRule::from(arg), expected);
    }

    #[rstest]
    // --memory-per-thread (default true)
    #[case(&["test"], true)]
    #[case(&["test", "--memory-per-thread"], true)]
    #[case(&["test", "--memory-per-thread", "true"], true)]
    #[case(&["test", "--memory-per-thread", "false"], false)]
    #[case(&["test", "--memory-per-thread=true"], true)]
    #[case(&["test", "--memory-per-thread=false"], false)]
    fn test_memory_per_thread_parsing(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.queue_memory.memory_per_thread, expected);
    }

    #[rstest]
    // --max-memory parses to a MemoryLimit; "auto" and explicit values both work.
    #[case(&["test"], MemoryLimit::Fixed(768 * 1024 * 1024))]
    #[case(&["test", "--max-memory", "auto"], MemoryLimit::Auto)]
    #[case(&["test", "--max-memory", "2GiB"], MemoryLimit::Fixed(2 * 1024 * 1024 * 1024))]
    #[case(&["test", "--max-memory=512M"], MemoryLimit::Fixed(512 * 1000 * 1000))]
    fn test_max_memory_parsing(#[case] args: &[&str], #[case] expected: MemoryLimit) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.queue_memory.max_memory, expected);
    }

    /// The accepted boolean spellings, asserted through the CLI rather than a
    /// parser function, so the test survives a change of parser.
    ///
    /// `on`/`off`/`1`/`0` are accepted since the switch to clap's
    /// `BoolishValueParser`; they were rejected by the previous hand-rolled
    /// `parse_bool`, and the cases below moved up from the rejected table.
    #[rstest]
    #[case::lower_true("true", true)]
    #[case::lower_false("false", false)]
    #[case::yes("yes", true)]
    #[case::no("no", false)]
    #[case::t("t", true)]
    #[case::f("f", false)]
    #[case::y("y", true)]
    #[case::n("n", false)]
    #[case::on("on", true)]
    #[case::off("off", false)]
    #[case::one("1", true)]
    #[case::zero("0", false)]
    #[case::mixed_case_true("TrUe", true)]
    #[case::upper_false("FALSE", false)]
    #[case::upper_yes("YES", true)]
    #[case::mixed_case_no("No", false)]
    #[case::upper_on("ON", true)]
    #[case::mixed_case_off("oFf", false)]
    fn bool_flag_accepts_spelling(#[case] value: &str, #[case] expected: bool) {
        let cmd = TestBoolFlags::try_parse_from(["test", "--memory-per-thread", value])
            .expect("boolean spelling should be accepted");
        assert_eq!(cmd.queue_memory.memory_per_thread, expected);
    }

    /// Anything outside that set is still rejected -- the parser widens what is
    /// accepted, it does not fall back to `true` for unrecognized input.
    #[rstest]
    #[case::empty("")]
    #[case::truncated_true("tru")]
    #[case::truncated_false("fals")]
    #[case::overlong_true("truee")]
    #[case::overlong_no("noo")]
    #[case::overlong_yes("yess")]
    #[case::not_a_bool("maybe")]
    #[case::leading_space(" true")]
    #[case::trailing_space("true ")]
    #[case::numeric_two("2")]
    #[case::negative_one("-1")]
    #[case::size_suffix("10G")]
    fn bool_flag_rejects_spelling(#[case] value: &str) {
        assert!(
            TestBoolFlags::try_parse_from(["test", "--memory-per-thread", value]).is_err(),
            "expected rejection for input: {value:?}"
        );
    }

    #[rstest]
    #[case(&["test", "--trim", "yes"], true)]
    #[case(&["test", "--trim", "no"], false)]
    #[case(&["test", "--trim", "y"], true)]
    #[case(&["test", "--trim", "n"], false)]
    #[case(&["test", "--trim", "t"], true)]
    #[case(&["test", "--trim", "f"], false)]
    #[case(&["test", "--trim", "YES"], true)]
    #[case(&["test", "--trim", "NO"], false)]
    #[case(&["test", "--trim=yes"], true)]
    #[case(&["test", "--trim=no"], false)]
    // Accepted since the switch to clap's `BoolishValueParser`; these four
    // moved up from the rejected table below.
    #[case(&["test", "--trim", "on"], true)]
    #[case(&["test", "--trim", "off"], false)]
    #[case(&["test", "--trim", "1"], true)]
    #[case(&["test", "--trim", "0"], false)]
    fn test_extended_bool_values_in_cli(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.consensus.trim, expected);
    }

    #[rstest]
    #[case(&["test", "--trim", "maybe"])]
    #[case(&["test", "--trim", "2"])]
    #[case(&["test", "--trim", "-1"])]
    #[case(&["test", "--trim", "tru"])]
    #[case(&["test", "--trim", "10G"])]
    fn test_extended_bool_values_in_cli_invalid(#[case] args: &[&str]) {
        assert!(TestBoolFlags::try_parse_from(args).is_err());
    }

    // -------------------------------------------------------------------------
    // cgroup-aware memory detection
    // -------------------------------------------------------------------------

    #[test]
    fn test_detect_total_memory_returns_nonzero() {
        let total = detect_total_memory();
        assert!(total > 0, "detect_total_memory returned 0");
        // Sanity: at least 64 MiB (even the smallest CI runner has more than this).
        assert!(total >= 64 * 1024 * 1024, "detect_total_memory returned < 64 MiB: {total}");
    }

    #[test]
    fn test_detect_total_memory_bounded_by_sysinfo() {
        // cgroup_limits().total_memory is min(cgroup_max, physical_ram), so
        // detect_total_memory() can never exceed what sysinfo reports.
        let total = detect_total_memory();
        let mut system = sysinfo::System::new();
        system.refresh_memory();
        let sysinfo_total = usize::try_from(system.total_memory()).unwrap_or(usize::MAX);
        assert!(
            total <= sysinfo_total,
            "cgroup-limited total {total} exceeded sysinfo total {sysinfo_total}"
        );
    }

    #[test]
    fn test_serialize_raw_bam_records_empty() {
        let records: Vec<Vec<u8>> = Vec::new();
        let mut output = Vec::new();
        let count = serialize_raw_bam_records(&records, &mut output).unwrap();
        assert_eq!(count, 0);
        assert!(output.is_empty());
    }

    #[test]
    fn test_serialize_raw_bam_records_single_frames_correctly() {
        let records = vec![vec![0xDEu8, 0xAD, 0xBE, 0xEF]];
        let mut output = Vec::new();
        let count = serialize_raw_bam_records(&records, &mut output).unwrap();
        assert_eq!(count, 1);
        // <u32 LE block_size = 4><body = DE AD BE EF>
        assert_eq!(output, vec![0x04, 0x00, 0x00, 0x00, 0xDE, 0xAD, 0xBE, 0xEF]);
    }

    #[test]
    fn test_serialize_raw_bam_records_multiple_frames_concatenated() {
        let records = vec![vec![0x11u8, 0x22], vec![0x33u8, 0x44, 0x55]];
        let mut output = Vec::new();
        let count = serialize_raw_bam_records(&records, &mut output).unwrap();
        assert_eq!(count, 2);
        assert_eq!(
            output,
            vec![0x02, 0x00, 0x00, 0x00, 0x11, 0x22, 0x03, 0x00, 0x00, 0x00, 0x33, 0x44, 0x55],
        );
    }

    #[test]
    fn test_serialize_raw_bam_records_reserves_capacity_upfront() {
        // Build a batch large enough that a naive per-record extend would
        // trigger at least one realloc growth of `output`. We then assert
        // that after one call `output.capacity()` is at least the exact
        // serialized size — proof that the upfront reserve happened.
        let records: Vec<Vec<u8>> = (0..32).map(|i| vec![i as u8; 64]).collect();
        let expected_size: usize = records.iter().map(|r| 4 + r.len()).sum();

        let mut output = Vec::new();
        serialize_raw_bam_records(&records, &mut output).unwrap();
        assert_eq!(output.len(), expected_size);
        assert!(
            output.capacity() >= expected_size,
            "capacity {} should be >= expected serialized size {expected_size}",
            output.capacity(),
        );
    }

    #[test]
    fn test_serialize_raw_bam_records_preserves_existing_output_content() {
        // Reserve must extend, not clear: writing into a non-empty buffer
        // (as the pipeline can do when reusing serialization scratch
        // buffers) must leave the existing prefix untouched.
        let mut output = vec![0xAAu8, 0xBB];
        let records = vec![vec![0x01u8]];
        serialize_raw_bam_records(&records, &mut output).unwrap();
        assert_eq!(output, vec![0xAA, 0xBB, 0x01, 0x00, 0x00, 0x00, 0x01]);
    }
}
