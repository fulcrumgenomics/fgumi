//! Common CLI options shared across commands.
//!
//! This module provides shared argument structures that can be composed into
//! command structs using `#[command(flatten)]`.

use std::path::{Path, PathBuf};
#[cfg(feature = "simplex")]
use std::sync::Arc;

use crate::assigner::Strategy;
#[cfg(feature = "simplex")]
use crate::logging::OperationTimer;
use crate::unified_pipeline::{BamPipelineConfig, SchedulerStrategy};
use crate::validation::validate_input_exists;
use bytesize::ByteSize;
use clap::Args;
#[cfg(feature = "simplex")]
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
    /// The default. A separation that small is summation-order noise rather than evidence,
    /// and the rule is independent of both read order and family depth.
    #[value(name = "ulp-relative")]
    UlpRelative,
    /// Reproduce fgbio's tie rule exactly, defects included.
    ///
    /// For byte-parity with fgbio output — chiefly cross-tool equivalency testing. Inherits
    /// fgbio's order dependence and scale dependence, so it will call a base off one ULP of
    /// accumulation noise.
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
/// alignments and any record that is unmapped **and** lacks a mapped mate.
///
/// fgbio applies this to *every* consensus caller (simplex, duplex, codec),
/// for both the single- and multi-threaded execution paths, before grouping
/// reads by molecular identifier. fgumi replicates it in the consensus
/// commands.
///
/// `--allow-unmapped` relaxes **only** the unmapped-without-mapped-mate rule so
/// fully-unmapped input can be consensus-called (mirroring
/// `fgumi group --allow-unmapped`). Secondary/supplementary alignments are
/// **always** dropped, matching fgbio — `--allow-unmapped` never lets
/// non-primary alignments into grouping.
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
    // --allow-unmapped relaxes only the mapped-or-mate-mapped eligibility check.
    if allow_unmapped {
        return true;
    }
    let is_mapped = flags & flags::UNMAPPED == 0;
    let has_mapped_mate = flags & flags::PAIRED != 0 && flags & flags::MATE_UNMAPPED == 0;
    is_mapped || has_mapped_mate
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
}

impl Default for BamIoOptions {
    fn default() -> Self {
        Self { input: PathBuf::new(), output: PathBuf::new(), async_reader: false }
    }
}

impl BamIoOptions {
    /// Construct a `BamIoOptions` from input and output paths. Leaves
    /// opt-in tuning flags (e.g. `async_reader`) at their default values.
    pub fn new(input: impl Into<PathBuf>, output: impl Into<PathBuf>) -> Self {
        Self { input: input.into(), output: output.into(), async_reader: false }
    }

    /// Build [`fgumi_bam_io::PipelineReaderOpts`] from the async-reader flag.
    pub fn pipeline_reader_opts(&self) -> fgumi_bam_io::PipelineReaderOpts {
        fgumi_bam_io::PipelineReaderOpts { async_reader: self.async_reader }
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

/// Refuse a secondary output that resolves to the same destination as `--output`.
///
/// A command's primary and secondary writers are both opened up front and
/// written from the same loop, so pointing them at one destination does not make
/// them take turns — it interleaves two BAMs. How that fails depends on the
/// destination, and both ways are silent:
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
/// check into the writer layer. See #715, which also covers `--stats`/`--metrics`
/// sharing a path with `--output`.
///
/// The null device is the one exempt destination: it discards every byte, so two
/// writers on it cannot corrupt each other and `-o /dev/null --rejects /dev/null`
/// is a legitimate compute-only run. Nothing else is exempt for merely not being
/// a regular file — a FIFO or `/dev/fd/N` is a *shared* byte stream, where the
/// two writers append through one file description and interleave two headers,
/// two block sequences, and two EOF markers into a stream no reader can parse.
/// That is the same failure as `-o - --rejects -`, so it is rejected the same way.
///
/// `secondary_flag` names the offending option in the error, since a command may
/// have more than one secondary output.
///
/// # Errors
///
/// Returns an error if `output` and `secondary` both name stdout, or both resolve
/// to the same destination and that destination is not the null device.
pub fn reject_colliding_outputs(
    output: &Path,
    secondary: Option<&PathBuf>,
    secondary_flag: &str,
) -> anyhow::Result<()> {
    let Some(secondary) = secondary else { return Ok(()) };

    if is_stdout_path(output) || is_stdout_path(secondary) {
        anyhow::ensure!(
            !(is_stdout_path(output) && is_stdout_path(secondary)),
            "--output and {secondary_flag} cannot both write to stdout: the two BAM streams \
             would interleave into one unreadable stream; give {secondary_flag} a path"
        );
        return Ok(());
    }

    if is_null_device(output) && is_null_device(secondary) {
        return Ok(());
    }

    if resolve_output_identity(output) == resolve_output_identity(secondary) {
        anyhow::bail!(
            "--output and {secondary_flag} both write to {}: the two BAM streams would \
             interleave on a shared stream, or overwrite each other byte for byte on a \
             regular file, and either way nothing can read the result; give \
             {secondary_flag} a different path",
            secondary.display()
        );
    }
    Ok(())
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
    #[arg(short = 'B', long = "output-per-base-tags", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub output_per_base_tags: bool,

    /// Quality-trim reads before consensus calling (removes low-quality bases from ends)
    #[arg(long = "trim", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
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
            tie_rule: TieRuleArg::default(),
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
    #[arg(long = "consensus-call-overlapping-bases", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
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
#[derive(Debug, Clone, Default, Args)]
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

/// Option controlling whether unmapped reads are processed by a consensus caller.
///
/// Shared by the `simplex`, `duplex`, and `codec` commands so the `--allow-unmapped` flag
/// (name, default, and boolean parsing) stays identical across all three.
#[derive(Debug, Clone, Args)]
pub struct AllowUnmappedOptions {
    /// Process reads that are unmapped and lack a mapped mate. By default
    /// (fgbio parity, `ConsensusCallingIterator.scala:56-58`) such reads — and
    /// all secondary/supplementary alignments — are dropped before consensus
    /// calling. Enable for consensus on unmapped input (e.g. ribosome/protein
    /// display), mirroring `fgumi group --allow-unmapped`.
    #[arg(long = "allow-unmapped", value_name = "ALLOW_UNMAPPED", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
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
    #[arg(long = "pipeline-stats", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool, hide = true)]
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
    #[arg(long = "deadlock-recover", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool, hide = true)]
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

/// The minimum per-thread memory budget (256 MiB).
pub(crate) const MIN_MEMORY_PER_THREAD: usize = 256 * 1024 * 1024;

/// Default auto-reserve cap: 10 GiB.
pub(crate) const AUTO_RESERVE_CAP: usize = 10 * 1024 * 1024 * 1024;

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
            log::debug!(
                "Auto memory: {} of {} ({}/thread × {} threads, reserve {})",
                ByteSize(budget as u64),
                ByteSize(total as u64),
                ByteSize((budget / threads) as u64),
                threads,
                ByteSize(margin as u64),
            );
            budget
        }
    };

    if budget > total {
        log::warn!(
            "Memory budget {} exceeds total host memory {}; this may cause OOM (or, for sort, earlier spill-to-disk)",
            ByteSize(budget as u64),
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
#[derive(Debug, Clone, Args)]
pub struct QueueMemoryOptions {
    /// Maximum memory for the pipeline queues.
    ///
    /// Default is "768" (768 MiB) per thread. Pass "auto" to detect host memory
    /// (cgroup-aware) and subtract --memory-reserve, so the budget shrinks to
    /// fit the host and the command self-throttles instead of OOM-ing. Explicit
    /// values like "512M", "1G", "4GiB" are per-thread when --memory-per-thread
    /// is enabled (default); plain numbers are MiB. Mirrors `fgumi sort`'s
    /// --max-memory.
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
    #[arg(long = "memory-per-thread", default_value = "true", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
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

    /// Logs the resolved memory configuration.
    ///
    /// # Arguments
    /// * `num_threads` - Number of threads used for the calculation.
    /// * `total_memory` - The resolved total budget from `calculate_memory_limit`.
    pub fn log_memory_config(&self, num_threads: usize, total_memory: u64) {
        let per_thread = self.memory_per_thread;
        if per_thread && num_threads > 1 {
            log::info!(
                "Queue memory budget: {} total ({}/thread × {} threads)",
                ByteSize(total_memory),
                ByteSize(total_memory / num_threads as u64),
                num_threads
            );
        } else {
            log::info!("Queue memory budget: {} total", ByteSize(total_memory));
        }
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

/// Parses a boolean value from a string, accepting: true/false, yes/no, y/n, t/f
/// (case-insensitive). Matches sopt/fgbio behavior.
pub(crate) fn parse_bool(s: &str) -> Result<bool, String> {
    match s.to_ascii_lowercase().as_str() {
        "true" | "t" | "yes" | "y" => Ok(true),
        "false" | "f" | "no" | "n" => Ok(false),
        _ => Err(format!("Invalid boolean value '{s}'. Expected: true|false|yes|no|y|n|t|f")),
    }
}

// Re-export from the library crate for backward compatibility.
pub(crate) use crate::system::detect_total_memory;
pub use crate::validation::parse_memory_size;

/// Builds a [`BamPipelineConfig`] from the common CLI option structs.
///
/// This consolidates the pipeline configuration boilerplate that is repeated
/// across all multi-threaded commands: auto-tuning, scheduler strategy,
/// stats collection, deadlock settings, and queue memory limits.
///
/// After calling this, commands can further customize the returned config
/// (e.g. setting `group_key_config` for raw-byte mode).
pub fn build_pipeline_config(
    scheduler_opts: &SchedulerOptions,
    compression: &CompressionOptions,
    queue_memory: &QueueMemoryOptions,
    num_threads: usize,
) -> anyhow::Result<BamPipelineConfig> {
    let mut config = BamPipelineConfig::auto_tuned(num_threads, compression.compression_level);
    config.pipeline.scheduler_strategy = scheduler_opts.strategy();
    if scheduler_opts.collect_stats() {
        config.pipeline = config.pipeline.with_stats(true);
    }
    config.pipeline.deadlock_timeout_secs = scheduler_opts.deadlock_timeout_secs();
    config.pipeline.deadlock_recover_enabled = scheduler_opts.deadlock_recover_enabled();

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

#[cfg(test)]
mod tests {
    use super::*;

    /// `reject_colliding_outputs` compares destinations, not the strings naming
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
    fn reject_colliding_outputs_compares_resolved_files(
        #[case] output: &str,
        #[case] secondary: &str,
        #[case] collides: bool,
    ) {
        let dir = tempfile::tempdir().expect("create temp dir");
        std::fs::create_dir_all(dir.path().join("sub")).expect("create subdirectory");
        let output = dir.path().join(output);
        let secondary = dir.path().join(secondary);

        let result = reject_colliding_outputs(&output, Some(&secondary), "--rejects");

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
    fn reject_colliding_outputs_handles_stdout_and_the_null_device(
        #[case] output: &str,
        #[case] secondary: &str,
        #[case] collides: bool,
    ) {
        let secondary = PathBuf::from(secondary);

        let result = reject_colliding_outputs(Path::new(output), Some(&secondary), "--rejects");

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
    fn reject_colliding_outputs_allows_an_absent_secondary(#[case] output: &str) {
        let result = reject_colliding_outputs(Path::new(output), None, "--rejects");
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
    fn reject_colliding_outputs_rejects_a_shared_fifo(
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

        let result = reject_colliding_outputs(&output, Some(&secondary), "--rejects");

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
        enable_logging();
        log_index_threshold(strategy, effective_edits, index_threshold);
    }

    /// Enable an at-Trace logger so `log::warn!`/`debug!` macros evaluate their
    /// arguments — without an enabled logger the `log` crate skips argument
    /// evaluation, leaving the formatting expressions inside the memory-budget
    /// warn/debug branches unexecuted under test. nextest runs each test in its
    /// own process, so `try_init` is local and idempotent.
    fn enable_logging() {
        let _ =
            env_logger::builder().is_test(true).filter_level(log::LevelFilter::Trace).try_init();
    }

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
        enable_logging();
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
        enable_logging(); // exercise the cap-warning and auto-debug log branches
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
        enable_logging(); // exercise the "budget exceeds host total" warn branch
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
    /// `ConsensusCallingIterator.scala:56-58`): secondary and supplementary alignments are
    /// always dropped, and an unmapped read is kept only when it has a mapped mate.
    ///
    /// `--allow-unmapped` (the `allow_unmapped` column) relaxes **only** the
    /// unmapped-without-mapped-mate rule: it makes otherwise-eligible primary alignments pass
    /// regardless of mapping, but secondary/supplementary alignments are still dropped. The
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
    #[case::unmapped_paired_mate_mapped_kept(
        fgumi_raw_bam::flags::UNMAPPED | fgumi_raw_bam::flags::PAIRED,
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

    #[rstest]
    #[case("true", true)]
    #[case("false", false)]
    #[case("yes", true)]
    #[case("no", false)]
    #[case("t", true)]
    #[case("f", false)]
    #[case("y", true)]
    #[case("n", false)]
    #[case("True", true)]
    #[case("TRUE", true)]
    #[case("False", false)]
    #[case("FALSE", false)]
    #[case("Yes", true)]
    #[case("YES", true)]
    #[case("No", false)]
    #[case("NO", false)]
    #[case("T", true)]
    #[case("F", false)]
    #[case("Y", true)]
    #[case("N", false)]
    #[case("tRuE", true)]
    #[case("fAlSe", false)]
    #[case("yEs", true)]
    fn test_parse_bool_valid(#[case] input: &str, #[case] expected: bool) {
        assert_eq!(parse_bool(input).expect("should parse"), expected);
    }

    #[rstest]
    #[case("")]
    #[case("tru")]
    #[case("fals")]
    #[case("truee")]
    #[case("noo")]
    #[case("yess")]
    #[case("maybe")]
    #[case("0")]
    #[case("1")]
    #[case("on")]
    #[case("off")]
    #[case(" true")]
    #[case("true ")]
    fn test_parse_bool_invalid(#[case] input: &str) {
        assert!(parse_bool(input).is_err(), "expected error for input: {input:?}");
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
    fn test_extended_bool_values_in_cli(#[case] args: &[&str], #[case] expected: bool) {
        let cmd = TestBoolFlags::try_parse_from(args).expect("valid CLI args should parse");
        assert_eq!(cmd.consensus.trim, expected);
    }

    #[rstest]
    #[case(&["test", "--trim", "maybe"])]
    #[case(&["test", "--trim", "0"])]
    #[case(&["test", "--trim", "1"])]
    #[case(&["test", "--trim", "on"])]
    #[case(&["test", "--trim", "off"])]
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
