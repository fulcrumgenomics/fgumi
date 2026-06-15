//! Stage-by-stage chain builder.
//!
//! The principled endpoint sketched in
//! `docs/design/refactor/2026-05-20-unified-chain-builder-design.md`
//! and implemented in Phase 3: `build_for` constructs a `ChainBuilder`,
//! calls `add_source`, walks `spec.stages` calling `add_stage` per
//! variant, calls `add_sink`, then `build()`s the final `BuiltPipeline`.
//!
//! Each `add_<stage>` method:
//!   1. Reads `self.spec.stage_opts.<stage>` (already validated present).
//!   2. Pushes the stage's canonical step sequence onto `self.pipeline`
//!      via [`PipelineBuilder::append_source`] /
//!      [`PipelineBuilder::append_step`].
//!   3. Pushes any `FinalizeHook` impl(s) onto `self.finalize`.
//!
//! ## Type-erasure compatibility note
//!
//! `PipelineBuilder::append_source` / `append_step` bypass the typed
//! `Chain<'b, O>` API. Type correctness is NOT checked at `build()` time —
//! `PipelineBuilder::build()` only validates that every output is wired
//! (there is no `Step<Input=A>` consuming an output of `Step<Output=B>`
//! check at compile time or at `build()` time). A type mismatch panics in
//! `TypedStep::resolve_input` at the first dispatch with "input handle
//! downcast failed — chain topology invariant" once `Pipeline::run` begins.
//! The panic is loud and immediate but it is a runtime check, not a
//! compile-time or build-time one.
//!
//! `ChainBuilder`'s callers are constrained: each `add_<stage>` method
//! pushes a known-typed step sequence onto the chain, and the sequencing
//! is verified by the unit tests that drive `build_<command>_chain`
//! end-to-end. Misuse from outside `chains::builder` is prevented by
//! the `pub(crate)` scope.

use std::sync::Arc;
use std::sync::atomic::AtomicU64;

use anyhow::{Result, anyhow, bail};
use noodles::sam::Header;

use crate::pipeline::chains::{
    BuiltPipeline, ChainSpec, FinalizeHook, PipelineStatsFinalizeHook, SinkSpec, SourceSpec, Stage,
    StageTimingFinalizeHook, build_pipeline_config_for_chain,
};
use crate::pipeline::core::builder::PipelineBuilder;
use crate::pipeline::core::item::HeapSize;
use crate::pipeline::core::topology::{BranchIdx, StepIdx};
use crate::pipeline::steps::tuning::BamPipelineTuning;

// ─────────────────────────────────────────────────────────────────────────────
// PendingSource — opened input held between new() and add_source().
// ─────────────────────────────────────────────────────────────────────────────

/// Holds the opened input(s) between [`ChainBuilder::new`] (which opens them)
/// and [`ChainBuilder::add_source`] (which consumes them to push the source
/// preamble onto the pipeline). The variant matches the [`SourceSpec`] the
/// spec carries; each variant carries enough state to construct the source
/// step sequence.
///
/// `add_source` matches on `PendingSource` and routes to the correct preamble.
/// Future variants (`Paired`, `Fastq`) are stubbed here so T3a.12 (zipper) and
/// Phase 5 (extract) can extend `add_source` without changing `ChainBuilder`'s
/// field layout.
///
/// All variants that carry large structs box them to prevent
/// `clippy::large_enum_variant` from flagging this enum.
pub(crate) enum PendingSource {
    /// `SourceSpec::Bam` or `SourceSpec::Sam` — a single open [`InputSource`].
    /// Boxed to keep enum size comparable to the smaller variants.
    Single(Box<crate::pipeline::steps::source::InputSource>),

    /// `SourceSpec::PairedBams` — zipper's two-input form: mapped BAM +
    /// unmapped BAM + reference path. T3a.12 (zipper migration) will extend
    /// `add_source` to consume this variant.
    #[allow(dead_code)]
    Paired {
        mapped: Box<crate::pipeline::steps::source::InputSource>,
        unmapped: Box<crate::pipeline::steps::source::InputSource>,
        reference: std::path::PathBuf,
    },

    /// `SourceSpec::Fastqs` — extract's input form.
    /// `add_source` consumes this to build the FASTQ read preamble:
    /// `ReadFastqInputs → ZipFastqRecords`. Read structures are read
    /// from `self.spec.source` by `add_extract`, not carried here.
    Fastq { paths: Vec<std::path::PathBuf> },
}

// ─────────────────────────────────────────────────────────────────────────────
// StagePosition — terminal vs intermediate position in the chain.
// ─────────────────────────────────────────────────────────────────────────────

/// Whether a stage is the last step before the output sink (`Terminal`) or is
/// followed by at least one more stage (`Intermediate`).
///
/// Each `add_<stage>` method receives a `StagePosition` argument. For
/// `Terminal`, the stage appends its serialise-to-bytes step so the chain tail
/// is [`DecompressedBlock`] ready for `BgzfCompress → WriteBgzfFile`. For
/// `Intermediate`, the stage leaves its native typed output on the chain tail
/// (e.g., `BatchedProcessedDedupGroups` for dedup) so the next stage can
/// consume it directly without a round-trip through bytes.
///
/// [`build_for`][`crate::pipeline::chains::build_for`] computes the
/// position for each stage from the spec's stage list:
///
/// ```text
/// let last_idx = spec.stages.len() - 1;
/// for (i, stage) in spec.stages.iter().enumerate() {
///     let position = if i == last_idx { StagePosition::Terminal }
///                    else             { StagePosition::Intermediate };
///     chain.add_stage(*stage, position)?;
/// }
/// ```
///
/// Single-stage chains always pass `Terminal`. Multi-stage chains (Phase 3b
/// runall fused paths) pass `Intermediate` for every stage except the last.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StagePosition {
    /// This stage is the last before the output sink. Append any
    /// serialise-to-bytes step needed to feed `BgzfCompress`.
    Terminal,
    /// This stage is followed by at least one more stage. Leave the native
    /// typed output on the chain tail for the next stage to consume.
    Intermediate,
}

// ─────────────────────────────────────────────────────────────────────────────
// ChainTailKind — tracks what type the current chain tail produces.
// ─────────────────────────────────────────────────────────────────────────────

/// Tracks the logical output type of the chain tail at each point in
/// construction. Used by multi-stage `add_<stage>` methods to decide how to
/// wire their preamble steps: the source preamble, or the previous stage's
/// output, may produce different types depending on the chain layout.
///
/// ## Why this is needed
///
/// `ChainBuilder` erases types after each `append_step` call (there is no
/// compile-time check that a downstream step's `Input` matches the upstream
/// step's `Output`). For most stages the source preamble always ends at
/// [`DecodedRecordBatch`] and the stages are type-compatible. However, two
/// multi-stage fused patterns need a different preamble:
///
/// 1. **Sort intermediate** — Sort consumes `RecordBatch`. When Sort is the
///    first stage, `add_source` emits `RecordBatch` via `ParseBamRecords`.
///    When Sort follows Align/Zipper/Correct (which emit `BamTemplateBatch`),
///    `add_sort` prepends `TemplatesToRecordBatch`. After Sort's three-step
///    chain, the tail is `DecodedRecordBatch` again.
///
/// 2. **Group→Consensus** — Consensus stages normally receive `DecodedRecordBatch`
///    and prepend `GroupByMi`. When Group runs as an intermediate stage,
///    the tail is `BatchedProcessedPositionGroups` after `MiAssign`. The
///    consensus stages must then prepend `TemplatesToMiGroups` instead of
///    `GroupByMi` to bridge the two types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ChainTailKind {
    /// The chain tail produces [`DecodedRecordBatch`]. This is the default
    /// after the source preamble (`DecodeRecords`) and after Sort (via
    /// `DecodeFromRecords`).
    ///
    /// [`DecodedRecordBatch`]: crate::pipeline::steps::types::DecodedRecordBatch
    DecodedRecordBatch,

    /// The chain tail produces [`RecordBatch`]. Set when Sort is the first
    /// stage: `add_source` uses `ParseBamRecords` (not `DecodeRecords`).
    ///
    /// [`RecordBatch`]: crate::pipeline::steps::types::RecordBatch
    RecordBatch,

    /// The chain tail produces [`BamTemplateBatch`]. Set by `add_align`
    /// (Intermediate) and `add_correct` (Intermediate).
    ///
    /// [`BamTemplateBatch`]: crate::pipeline::steps::types::BamTemplateBatch
    BamTemplateBatch,

    /// The chain tail produces [`BatchedProcessedPositionGroups`]. Set by
    /// `add_group` (Intermediate) — after `GroupByPosition → ProcessOrdered
    /// → MiAssign` but before `SerializeBamRecords`.
    ///
    /// [`BatchedProcessedPositionGroups`]: crate::pipeline::steps::group::position::BatchedProcessedPositionGroups
    BatchedProcessedPositionGroups,

    /// The chain tail produces
    /// [`FastqTemplateBatch`]. Set by `add_source` on the `Fastqs` arm,
    /// after `ReadFastqInputs → ZipFastqRecords`. The only consumer is
    /// `add_extract`, which converts it to `BamTemplateBatch`.
    ///
    /// [`FastqTemplateBatch`]: crate::pipeline::steps::source::zip_fastq::FastqTemplateBatch
    FastqTemplateBatch,
}

// ─────────────────────────────────────────────────────────────────────────────
// ChainBuilder
// ─────────────────────────────────────────────────────────────────────────────

/// In-progress chain builder. Constructed by
/// [`crate::pipeline::chains::build_for`] (or by a per-command
/// builder during Phase 3a).
///
/// Owns the resolved output header, an accumulating `PipelineBuilder`, and a
/// growing `Vec<Box<dyn FinalizeHook>>` populated by `add_<stage>`
/// methods. Per-stage methods are private — callers drive the chain
/// via the public `add_source` / `add_stage(Stage, StagePosition)` /
/// `add_sink` / `build()` flow.
///
/// ## Type erasure
///
/// The framework's `Chain<'b, O>` provides compile-time step-compatibility
/// enforcement but cannot span method boundaries on `&mut self`. Instead,
/// `ChainBuilder` tracks the chain tail as `(StepIdx, BranchIdx)` and uses
/// `PipelineBuilder::append_source` / `append_step`, which bypass the
/// typed-chain API. See the module-level type-erasure note for the runtime
/// behaviour when a type mismatch is introduced.
pub(crate) struct FuseState {
    scratch: Vec<u8>,
    mi_buf: String,
}

impl HeapSize for FuseState {}

pub struct ChainBuilder<'a> {
    spec: &'a ChainSpec,
    tuning: BamPipelineTuning,

    /// Output header (with `@PG` record injected). Populated by `new()`,
    /// consumed by `add_sink()` (for the BAM writer header) and by
    /// `add_source()` (SAM parse step uses it).
    header: Header,

    /// In-progress chain state. `None` before the first step (before
    /// `add_source` runs); `Some((producer, branch))` after. Updated
    /// by every `add_source` / `add_<stage>` / `add_sink` call.
    current_tail: Option<(StepIdx, BranchIdx)>,

    /// Accumulating pipeline. Steps are registered via
    /// `append_source` / `append_step` in `add_*` methods.
    pipeline: PipelineBuilder,

    /// Post-pipeline cleanup hooks. Populated by `add_<stage>` methods.
    finalize: Vec<Box<dyn FinalizeHook>>,

    /// Shared progress counter threaded through source + stage steps for
    /// progress logging. `Arc` so the stage step and the finalize hook can
    /// both reference it.
    #[allow(dead_code)] // populated by add_<stage>; used by finalize hooks
    progress_records: Arc<AtomicU64>,

    /// The opened input source(s), held here until `add_source` consumes them.
    /// `None` after `add_source` has been called.
    pending_source: Option<PendingSource>,

    /// Second source chain tail for zipper's `PairedBams` path.
    ///
    /// Set by `add_source` when the pending source is
    /// [`PendingSource::Paired`]: `current_tail` holds the unmapped
    /// source chain tail and `paired_tail` holds the mapped source
    /// chain tail. `add_zipper` consumes both to build
    /// `ZipperMergeStep` via `PipelineBuilder::append_step2`.
    /// `None` for all single-source stages.
    paired_tail: Option<(
        crate::pipeline::core::topology::StepIdx,
        crate::pipeline::core::topology::BranchIdx,
    )>,

    /// Override for `PipelineConfig::threads` applied in [`Self::build`].
    ///
    /// `None` (the default) leaves `config.threads` set to
    /// `spec.threading.num_threads()`. `Some(n)` forces it to `n`.
    ///
    /// Currently only [`Self::add_sort`] sets this: sort wraps a single
    /// `Exclusive` [`SortBamFile`] step that drives the sort engine's own
    /// internal thread pool. The framework must run **one** driver thread
    /// (`threads = 1`) regardless of how many sorter threads were requested
    /// via `--threads`. Without this override the pipeline would spin up
    /// `num_sorter_threads` framework workers for a chain that has only one
    /// step, wasting OS threads.
    ///
    /// [`SortBamFile`]: crate::pipeline::steps::sort::SortBamFile
    override_pipeline_threads: Option<usize>,

    /// `HeaderHandle` stashed by [`Self::add_align`] for consumption by
    /// [`Self::add_sink`].
    ///
    /// `AlignAndMergeStep` resolves the final merged header at runtime
    /// (after the aligner emits its own `@PG`/`@RG`/`@CO` lines). The
    /// downstream [`WriteBgzfFile`] must therefore be constructed with
    /// [`WriteBgzfFile::new_with_handle`], which blocks until the handle
    /// is resolved before writing the BAM header bytes.
    ///
    /// `add_align` sets this; `add_sink` checks it and dispatches to
    /// `new_with_handle` when `Some`, falling back to the usual `new`
    /// path when `None` (all other stages).
    ///
    /// [`WriteBgzfFile`]: crate::pipeline::steps::sink::write_bgzf::WriteBgzfFile
    pending_header_handle: Option<crate::pipeline::core::header::HeaderHandle>,

    /// Tracks the logical output type of the current chain tail. Updated by
    /// `add_source` and each `add_<stage>` method. See [`ChainTailKind`]
    /// for the full enumeration and rationale.
    ///
    /// Default is [`ChainTailKind::DecodedRecordBatch`], set by `add_source`
    /// on the normal BAM/SAM source path. Changes:
    ///
    /// - `add_source` sets `RecordBatch` when Sort is the first stage (uses
    ///   `ParseBamRecords` instead of `DecodeRecords`).
    /// - `add_correct(Intermediate)` sets `BamTemplateBatch`.
    /// - `add_align(Intermediate)` sets `BamTemplateBatch`.
    /// - `add_sort(Intermediate)` reads the kind to decide its preamble,
    ///   then sets `DecodedRecordBatch` after `DecodeFromRecords`.
    /// - `add_group(Intermediate)` sets `BatchedProcessedPositionGroups`.
    /// - Consensus stages read the kind to decide their grouping preamble
    ///   (`GroupByMi` for `DecodedRecordBatch`, `TemplatesToMiGroups` for
    ///   `BatchedProcessedPositionGroups`).
    chain_tail_kind: ChainTailKind,
}

impl<'a> ChainBuilder<'a> {
    /// Construct a `ChainBuilder` from a validated [`ChainSpec`].
    ///
    /// Opens the input source (for BAM/SAM specs) to resolve the header
    /// and stores the reader in `pending_source` for [`Self::add_source`]
    /// to consume. Applies `@PG` injection; stores the output header.
    ///
    /// # Errors
    ///
    /// Returns I/O errors from input file open or header parsing.
    pub fn new(spec: &'a ChainSpec) -> Result<Self> {
        let num_threads = spec.threading.num_threads();
        let tuning = BamPipelineTuning::auto_tuned(num_threads)
            .with_compression_level(spec.compression.compression_level);

        let (raw_header, pending_source) = Self::open_source(spec)?;
        let header = crate::commands::common::add_pg_record(raw_header, &spec.command_line)?;

        Ok(Self {
            spec,
            tuning,
            header,
            current_tail: None,
            pipeline: PipelineBuilder::new(),
            finalize: Vec::new(),
            progress_records: Arc::new(AtomicU64::new(0)),
            pending_source,
            paired_tail: None,
            override_pipeline_threads: None,
            pending_header_handle: None,
            // Initialise to the default; add_source will set the correct kind
            // based on whether sort is the first intermediate stage.
            chain_tail_kind: ChainTailKind::DecodedRecordBatch,
        })
    }

    /// Open the input source and extract the header. For `Bam`/`Sam` specs,
    /// opens the file and wraps it in [`PendingSource::Single`]. For
    /// `PairedBams`, opens both inputs, resolves the dict path, builds the
    /// merged output header via `build_output_header`, and wraps them in
    /// [`PendingSource::Paired`]. `Fastqs` returns an error (Phase 5).
    fn open_source(spec: &ChainSpec) -> Result<(Header, Option<PendingSource>)> {
        use crate::pipeline::steps::source::InputSource;

        match &spec.source {
            SourceSpec::Bam(path) | SourceSpec::Sam(path) => {
                // Validate the file exists before opening (stdin exempt).
                if !fgumi_bam_io::is_stdin_path(path) {
                    crate::validation::validate_file_exists(path, "input BAM/SAM file")?;
                }
                let input = InputSource::open_with_opts(
                    path,
                    fgumi_bam_io::PipelineReaderOpts { async_reader: spec.async_reader },
                )
                .map_err(|e| anyhow!("open input: {e}"))?;
                let header = input.header().clone();
                Ok((header, Some(PendingSource::Single(Box::new(input)))))
            }
            SourceSpec::PairedBams { unmapped, mapped, reference } => {
                // Validate inputs exist (stdin exempt for both).
                if !fgumi_bam_io::is_stdin_path(unmapped) {
                    crate::validation::validate_file_exists(unmapped, "unmapped BAM file")?;
                }
                if !fgumi_bam_io::is_stdin_path(mapped) {
                    crate::validation::validate_file_exists(mapped, "mapped BAM file")?;
                }
                crate::validation::validate_file_exists(reference, "reference FASTA file")?;

                let dict_path = crate::reference::find_dict_path(reference).ok_or_else(|| {
                    anyhow!(
                        "Reference dictionary file not found. Tried:\n  \
                            - {}\n  \
                            - {}.dict\n\
                            Please run: samtools dict {} -o {}",
                        reference.with_extension("dict").display(),
                        reference.display(),
                        reference.display(),
                        reference.with_extension("dict").display()
                    )
                })?;

                // Unmapped source must be BAM (fgumi extract always produces BAM).
                let unmapped_input = InputSource::open_with_opts(
                    unmapped,
                    fgumi_bam_io::PipelineReaderOpts { async_reader: spec.async_reader },
                )
                .map_err(|e| anyhow!("open unmapped: {e}"))?;
                if matches!(unmapped_input, InputSource::Sam { .. }) {
                    bail!(
                        "SAM input is not supported for --unmapped \
                        (expected `fgumi extract` BAM output)"
                    );
                }
                let unmapped_header = unmapped_input.header().clone();

                // Mapped source: BAM or SAM (SAM from aligner is preferred).
                let mapped_input = InputSource::open_with_opts(
                    mapped,
                    fgumi_bam_io::PipelineReaderOpts { async_reader: spec.async_reader },
                )
                .map_err(|e| anyhow!("open mapped: {e}"))?;
                let mapped_header = mapped_input.header().clone();

                // Build the merged output header from the two inputs + dict.
                // `new()` will apply `add_pg_record` on top of this.
                let output_header = crate::commands::zipper::build_output_header(
                    &unmapped_header,
                    &mapped_header,
                    &dict_path,
                )?;

                Ok((
                    output_header,
                    Some(PendingSource::Paired {
                        mapped: Box::new(mapped_input),
                        unmapped: Box::new(unmapped_input),
                        reference: reference.clone(),
                    }),
                ))
            }
            SourceSpec::Fastqs { paths, read_structures: _ } => {
                let extract_opts = spec.stage_opts.extract.as_ref().ok_or_else(|| {
                    anyhow!("Fastqs source requires extract options in StageOptionsBag")
                })?;
                let header =
                    crate::pipeline::chains::commands::extract::build_fastq_header(extract_opts)?;
                Ok((header, Some(PendingSource::Fastq { paths: paths.clone() })))
            }
        }
    }

    /// Add the input source step(s) to the pipeline.
    ///
    /// Consumes the [`PendingSource`] stored by `new()`.
    /// - [`PendingSource::Single`] with a BAM reader: 4-step preamble
    ///   (`ReadBgzfBlocks → BgzfDecompress → FindBamBoundaries → DecodeRecords`).
    /// - [`PendingSource::Single`] with a SAM reader: 2-step preamble
    ///   (`ReadSamChunks → ParseSamChunk`).
    /// - [`PendingSource::Paired`]: Builds both source preambles for zipper's
    ///   dual-input path. `current_tail` = unmapped chain tail;
    ///   `paired_tail` = mapped chain tail. No standard single-output tail is
    ///   set — `add_zipper` consumes both tails to assemble `ZipperMergeStep`.
    ///   Note that both preamble chains end at `GroupByQueryname` so their
    ///   output type is `OrderedBytesSingle<BamTemplateBatch>`, matching
    ///   `ZipperMergeStep: Step2<InputA = BamTemplateBatch, InputB = BamTemplateBatch>`.
    ///   The reference path is stashed in `PendingSource::Paired` but is
    ///   consumed later by `add_zipper` (reference loading is deferred because
    ///   `restore_unconverted_bases` controls whether the FASTA is opened at all).
    /// - [`PendingSource::Fastq`]: bails — implemented in Phase 5.
    ///
    /// # Errors
    ///
    /// Returns errors from step construction or from unsupported source variants.
    ///
    /// # Panics
    ///
    /// Panics if called twice (after `pending_source` is already consumed).
    #[allow(clippy::too_many_lines)]
    pub fn add_source(&mut self) -> Result<()> {
        use crate::pipeline::steps::bgzf::decompress::BgzfDecompress;
        use crate::pipeline::steps::boundaries::bam::FindBamBoundaries;
        use crate::pipeline::steps::parse::bam::ParseBamRecords;
        use crate::pipeline::steps::source::InputSource;
        use crate::pipeline::steps::source::read_bam::read_bam_from_reader;
        use crate::sam::check_sort;

        let source = self.pending_source.take().expect("add_source called twice");

        // Detect whether Sort is the first intermediate stage. When it is,
        // the source preamble must end at `ParseBamRecords → RecordBatch`
        // (rather than `DecodeRecords → DecodedRecordBatch`) so that
        // `SortAndSpill` (which consumes `RecordBatch`) can be wired
        // directly to the source output without an extra round-trip.
        //
        // This replaces the fused dispatcher's use of `ParseBamRecords`
        // before `SortAndSpill` in `execute_group_pipeline`'s
        // `is_sort_enabled()` branch, where `--start-from sort --stop-after
        // {group,consensus}` fed the source through `ParseBamRecords`
        // instead of `DecodeRecords`.
        let sort_is_first_intermediate = matches!(
            self.spec.stages.as_slice(),
            [Stage::Sort, ..] if self.spec.stages.len() > 1
        );

        match source {
            PendingSource::Single(boxed_input) => {
                let input = *boxed_input;
                // GroupKeyConfig: used by DecodeRecords (BAM) and ParseSamChunk (SAM).
                // Queryname-grouping first stages (correct) skip cell-barcode
                // extraction — see `source_group_key_config`.
                let group_key_config = self.source_group_key_config();

                match input {
                    InputSource::Bam { reader, .. } => {
                        if sort_is_first_intermediate {
                            // Sort-as-intermediate needs `RecordBatch` input.
                            // Use ParseBamRecords instead of DecodeRecords.
                            let (read_step, _) = read_bam_from_reader(
                                reader,
                                self.header.clone(),
                                self.tuning.blocks_per_batch,
                                self.tuning.per_step_byte_limit,
                            );
                            let tail = self.pipeline.append_source(read_step);
                            let tail = self.pipeline.append_step(
                                BgzfDecompress::new(self.tuning.per_step_byte_limit),
                                tail,
                            );
                            let tail = self.pipeline.append_step(
                                FindBamBoundaries::new(self.tuning.per_step_byte_limit),
                                tail,
                            );
                            let tail = self.pipeline.append_step(
                                ParseBamRecords::new(self.tuning.per_step_byte_limit),
                                tail,
                            );
                            self.current_tail = Some(tail);
                            self.chain_tail_kind = ChainTailKind::RecordBatch;
                        } else {
                            let tail = self.build_bam_decode_preamble(
                                reader,
                                self.header.clone(),
                                group_key_config,
                            );
                            self.current_tail = Some(tail);
                            self.chain_tail_kind = ChainTailKind::DecodedRecordBatch;
                        }
                    }
                    InputSource::Sam { reader: sam_reader, .. } => {
                        let buf = sam_reader.into_inner();
                        let tail = self.build_sam_parse_preamble(
                            buf,
                            self.header.clone(),
                            group_key_config,
                        );
                        self.current_tail = Some(tail);
                        self.chain_tail_kind = ChainTailKind::DecodedRecordBatch;
                    }
                }
            }
            PendingSource::Paired { mapped: boxed_mapped, unmapped: boxed_unmapped, reference } => {
                // Zipper's two-input preamble. Both chains end at GroupByQueryname
                // so their output type is OrderedBytesSingle<BamTemplateBatch>,
                // matching ZipperMergeStep's InputA / InputB.
                //
                // The reference path is carried in PendingSource::Paired but
                // add_zipper reads it directly from self.spec.source, which is
                // always SourceSpec::PairedBams for zipper. Dropping `reference`
                // here is intentional.

                let mapped = *boxed_mapped;
                let unmapped = *boxed_unmapped;
                let _ = reference; // consumed by add_zipper via self.spec.source

                let unmapped_header = unmapped.header().clone();
                let mapped_header = mapped.header().clone();

                // Resolve paths for check_sort log messages (stdin → "<stdin>").
                let (unmapped_path, mapped_path) = match &self.spec.source {
                    SourceSpec::PairedBams { unmapped: u, mapped: m, .. } => (u.clone(), m.clone()),
                    _ => unreachable!("PendingSource::Paired with non-PairedBams spec"),
                };

                check_sort(&unmapped_header, &unmapped_path, "unmapped");
                check_sort(&mapped_header, &mapped_path, "mapped");

                // Warn if mapped input is BAM (SAM from aligner is the fast path).
                if matches!(mapped, InputSource::Bam { .. }) {
                    log::warn!(
                        "BAM input detected for --input. For best performance, pipe SAM directly \
                         from the aligner (e.g. bwa mem ... | fgumi zipper ...)."
                    );
                }

                // GroupKeyConfig used by DecodeRecords (BAM) and ParseSamChunk (SAM).
                let group_key_config = self.bam_group_key_config();

                // ── Unmapped preamble: always BAM (guaranteed by open_source) ──
                let unmapped_tail = match unmapped {
                    InputSource::Bam { reader, header: bam_header } => self
                        .build_bam_decode_then_group_preamble(
                            reader,
                            bam_header,
                            group_key_config.clone(),
                        ),
                    InputSource::Sam { .. } => {
                        // Guarded in open_source; belt-and-suspenders.
                        bail!(
                            "SAM input is not supported for --unmapped \
                            (expected `fgumi extract` BAM output)"
                        );
                    }
                };

                // ── Mapped preamble: BAM or SAM ────────────────────────────────
                let mapped_tail = match mapped {
                    InputSource::Bam { reader, header: bam_header } => self
                        .build_bam_decode_then_group_preamble(reader, bam_header, group_key_config),
                    InputSource::Sam { reader: rdr, header: hdr } => {
                        let buf = rdr.into_inner();
                        self.build_sam_parse_then_group_preamble(buf, hdr, group_key_config)
                    }
                };

                // Store both tails. current_tail = unmapped; paired_tail = mapped.
                // add_zipper will consume both via PipelineBuilder::append_step2.
                self.current_tail = Some(unmapped_tail);
                self.paired_tail = Some(mapped_tail);
            }
            PendingSource::Fastq { paths } => {
                use crate::pipeline::core::step::Affinity;
                use crate::pipeline::steps::source::pair_fastq::PairRawFastq;
                use crate::pipeline::steps::source::parse_fastq::ParseFastqChunks;
                use crate::pipeline::steps::source::parse_zip_fastq::ParseAndZipFastq;
                use crate::pipeline::steps::source::read_fastq::ReadFastqInputs;
                use crate::pipeline::steps::source::zip_fastq::ZipFastqRecords;

                // Open one BufRead reader per FASTQ path. Extract's
                // `open_fastq_reader` handles BGZF / gzip / plain
                // detection and optional async prefetch wrapping.
                let extract_opts = self.spec.stage_opts.extract.as_ref().ok_or_else(|| {
                    anyhow!("Fastq source requires extract options in StageOptionsBag")
                })?;
                let async_reader = extract_opts.async_reader;
                // FASTQ decompression threads: for BGZF inputs the reader
                // can multi-thread; for gzip/plain it is single-threaded
                // regardless.  Pass 1 here — the pipeline framework
                // provides the parallelism via typed steps.
                let mut readers: Vec<Box<dyn std::io::BufRead + Send>> = paths
                    .iter()
                    .map(|p| crate::commands::extract::open_fastq_reader(p, 1, async_reader))
                    .collect::<Result<Vec<_>>>()?;

                let n_streams = readers.len();
                // batch_record_count — same default as the legacy pipeline.
                let batch_records = 400usize;
                let byte_limit = self.tuning.per_step_byte_limit;
                // Worker count the pipeline will actually run with. Used to
                // pin per-stream readers to distinct, in-range workers.
                let num_threads = self.spec.threading.num_threads().max(1);

                // gzip decompression is the FASTQ bottleneck. The framework
                // runs distinct `Serial` steps on distinct workers
                // concurrently (Serial = per-step mutex, not global), so we
                // instantiate one reader PER stream for the common paired-end
                // (N == 2) case — both decompressors run at once. For N == 1 a
                // single reader is trivially enough; for N >= 3 we fall back to
                // a single all-streams round-robin reader.
                //
                // Topology by stream count (all arms produce a
                // `FastqTemplateBatch`-yielding tail):
                //   N == 2: ReadFastqInputs×2 → PairRawFastq (Serial, cheap
                //           chunk-level pairing + serial ordinal mint) →
                //           ParseAndZipFastq (Parallel: parse both streams'
                //           bytes AND build templates). This mirrors the legacy
                //           pipeline's pair-then-parse order and moves the
                //           expensive record-level template build into a
                //           parallel step.
                //   N != 2: ReadFastqInputs → ParseFastqChunks (Parallel parse)
                //           → ZipFastqRecords (Serial record-level join). The
                //           2-way `PairRawFastq` Step2 only handles two inputs,
                //           so single- and ≥3-stream cases keep this structure.
                //
                // Each per-stream reader is pinned to a DISTINCT worker via
                // `Affinity::Worker`. Two `Serial` sources sharing a mutex
                // dispatched by more than one worker each would race on the
                // runtime's `try_lock`-after-`Finished` source-drain path, so
                // disjoint single-worker affinities (never `Affinity::None`)
                // are required for correctness — not just throughput. The
                // worker index is clamped to `num_threads - 1` so a low
                // `--threads` count can never request a non-existent worker
                // (which would deadlock).
                let tail = match n_streams {
                    1 => {
                        let only = readers.pop().expect("n_streams == 1");
                        let read_step = ReadFastqInputs::new_single(
                            only,
                            0, // global_stream_idx
                            1, // n_streams_total
                            Affinity::Worker(0),
                            batch_records,
                            byte_limit,
                        );
                        let tail = self.pipeline.append_source(read_step);
                        let parse_tail =
                            self.pipeline.append_step(ParseFastqChunks::new(byte_limit), tail);
                        self.pipeline.append_step(ZipFastqRecords::new(1, byte_limit), parse_tail)
                    }
                    2 => {
                        // Two concurrent single-stream readers (R1, R2) →
                        // PairRawFastq (Step2) → ParseAndZipFastq.
                        let r2_in = readers.pop().expect("n_streams == 2");
                        let r1_in = readers.pop().expect("n_streams == 2");
                        // R1 → worker 0; R2 → worker 1 when it exists
                        // (clamped to worker 0 at --threads 1, where the two
                        // readers necessarily serialize on the sole worker).
                        let r2_worker = (num_threads - 1).min(1);
                        let r1_step = ReadFastqInputs::new_single(
                            r1_in,
                            0,
                            2,
                            Affinity::Worker(0),
                            batch_records,
                            byte_limit,
                        );
                        let r2_step = ReadFastqInputs::new_single(
                            r2_in,
                            1,
                            2,
                            Affinity::Worker(r2_worker),
                            batch_records,
                            byte_limit,
                        );
                        let r1_tail = self.pipeline.append_source(r1_step);
                        let r2_tail = self.pipeline.append_source(r2_step);
                        // `PairRawFastq` (Serial) does cheap chunk-level
                        // pairing and mints the globally-unique ordinal that
                        // `ParseAndZipFastq` (Parallel) reorders by. The latter
                        // parses both streams' bytes and builds the templates,
                        // so no separate `ZipFastqRecords` is appended for
                        // N == 2.
                        let pair_tail = self.pipeline.append_step2(
                            PairRawFastq::new(byte_limit),
                            r1_tail,
                            r2_tail,
                        );
                        self.pipeline.append_step(ParseAndZipFastq::new(byte_limit), pair_tail)
                    }
                    _ => {
                        // N >= 3 fallback: single all-streams round-robin
                        // reader (Affinity::Reader), serial decompress.
                        let read_step = ReadFastqInputs::new(readers, batch_records, byte_limit);
                        let tail = self.pipeline.append_source(read_step);
                        let parse_tail =
                            self.pipeline.append_step(ParseFastqChunks::new(byte_limit), tail);
                        self.pipeline
                            .append_step(ZipFastqRecords::new(n_streams, byte_limit), parse_tail)
                    }
                };

                self.current_tail = Some(tail);
                // Every arm's tail emits FastqTemplateBatch; the only
                // consumer is add_extract.
                self.chain_tail_kind = ChainTailKind::FastqTemplateBatch;
            }
        }

        Ok(())
    }

    /// Resolve the input path used for log messages. For `PairedBams`,
    /// returns the mapped path (the user-facing "main" input). For
    /// `Fastqs`, returns the first path (R1). The result is for display only —
    /// the actual source reader(s) are created by [`Self::add_source`] from
    /// the full [`SourceSpec`].
    ///
    /// Lets per-stage methods (`add_group`, `add_simplex`, `add_duplex`,
    /// `add_codec`, `add_clip`, `add_filter`, `add_correct`) accept all
    /// `SourceSpec` variants for logging without having to bail on
    /// fused-chain compositions where their upstream is a multi-input source.
    fn resolve_log_input_path(&self) -> std::path::PathBuf {
        match &self.spec.source {
            SourceSpec::Bam(p) | SourceSpec::Sam(p) => p.clone(),
            SourceSpec::PairedBams { mapped, .. } => mapped.clone(),
            SourceSpec::Fastqs { paths, .. } => paths.first().cloned().unwrap_or_default(),
        }
    }

    /// Build a [`GroupKeyConfig`] from the **current** `self.header`.
    ///
    /// Reads the library index (from `@RG` `LB` fields) and the cell-barcode
    /// tag. Callers must invoke this at the point in chain construction where
    /// `self.header` already reflects the records being decoded: `add_source`
    /// calls it before any stage replaces the header; the consensus stages call
    /// it *after* `self.header` has been replaced with the consensus output
    /// header (the records being decoded are consensus output records). Do not
    /// hoist a call above a `self.header = …` reassignment.
    ///
    /// [`GroupKeyConfig`]: fgumi_bam_io::GroupKeyConfig
    fn bam_group_key_config(&self) -> fgumi_bam_io::GroupKeyConfig {
        use crate::sam::SamTag;
        use noodles::sam::alignment::record::data::field::Tag;
        let cell_tag = Tag::from(SamTag::CB);
        let library_index = fgumi_bam_io::LibraryIndex::from_header(&self.header);
        fgumi_bam_io::GroupKeyConfig::new(library_index, cell_tag)
    }

    /// Group-key config for the **source preamble's** `DecodeRecords`.
    ///
    /// A first stage that groups by queryname (correct's `GroupByQueryname`)
    /// reads only `key.name_hash` and discards the rest of the `GroupKey`, so
    /// the source decode uses
    /// [`fgumi_bam_io::GroupKeyConfig::name_hash_only`] — skipping the CIGAR
    /// 5′-position walk and the RG/CB/MC aux-tag extraction pass entirely.
    /// (The legacy single-threaded path uses `new_raw_no_cell`, which still
    /// pays the combined aux pass; name-hash-only is strictly less work.) All
    /// other first stages (group/dedup/consensus/clip) need the full
    /// position/cell key, so they fall through to [`Self::bam_group_key_config`].
    fn source_group_key_config(&self) -> fgumi_bam_io::GroupKeyConfig {
        if matches!(self.spec.stages.first(), Some(Stage::Correct)) {
            let library_index = fgumi_bam_io::LibraryIndex::from_header(&self.header);
            fgumi_bam_io::GroupKeyConfig::name_hash_only(library_index)
        } else {
            self.bam_group_key_config()
        }
    }

    /// Finish a consensus stage's chain tail.
    ///
    /// The consensus step emits a record-aligned [`DecompressedBlock`]. For a
    /// [`StagePosition::Terminal`] consensus stage, that block IS the chain tail
    /// (`add_sink` wires `BgzfCompress → Write`); `chain_tail_kind` is left at
    /// its sentinel (the enum has no `DecompressedBlock` variant and nothing
    /// downstream of a terminal consensus reads it). For a
    /// [`StagePosition::Intermediate`] consensus stage, a downstream stage
    /// (e.g. filter) consumes records, so append the existing [`DecodeRecords`]
    /// step — consensus output is already record-aligned, so no
    /// `FindBamBoundaries`/`ParseBamRecords` is needed — yielding
    /// [`DecodedRecordBatch`].
    ///
    /// [`DecompressedBlock`]: crate::pipeline::steps::types::DecompressedBlock
    /// [`DecodedRecordBatch`]: crate::pipeline::steps::types::DecodedRecordBatch
    fn finish_consensus_tail(
        &mut self,
        tail: (
            crate::pipeline::core::topology::StepIdx,
            crate::pipeline::core::topology::BranchIdx,
        ),
        position: StagePosition,
    ) -> (crate::pipeline::core::topology::StepIdx, crate::pipeline::core::topology::BranchIdx)
    {
        match position {
            // terminal: tail is DecompressedBlock; chain_tail_kind left as sentinel.
            StagePosition::Terminal => tail,
            StagePosition::Intermediate => {
                use crate::pipeline::steps::parse::decode::DecodeRecords;
                let group_key_config = self.bam_group_key_config();
                let tail = self.pipeline.append_step(
                    DecodeRecords::new(group_key_config, self.tuning.per_step_byte_limit),
                    tail,
                );
                self.chain_tail_kind = ChainTailKind::DecodedRecordBatch;
                tail
            }
        }
    }

    /// Append the 4-step BAM decode preamble:
    /// `ReadBgzfBlocks → BgzfDecompress → FindBamBoundaries → DecodeRecords`.
    fn build_bam_decode_preamble(
        &mut self,
        reader: Box<dyn std::io::Read + Send>,
        header: Header,
        group_key_config: fgumi_bam_io::GroupKeyConfig,
    ) -> (crate::pipeline::core::topology::StepIdx, crate::pipeline::core::topology::BranchIdx)
    {
        use crate::pipeline::steps::bgzf::decompress::BgzfDecompress;
        use crate::pipeline::steps::boundaries::bam::FindBamBoundaries;
        use crate::pipeline::steps::parse::decode::DecodeRecords;
        use crate::pipeline::steps::source::read_bam::read_bam_from_reader;

        let (read_step, _) = read_bam_from_reader(
            reader,
            header,
            self.tuning.blocks_per_batch,
            self.tuning.per_step_byte_limit,
        );
        let tail = self.pipeline.append_source(read_step);
        let tail =
            self.pipeline.append_step(BgzfDecompress::new(self.tuning.per_step_byte_limit), tail);
        let tail = self
            .pipeline
            .append_step(FindBamBoundaries::new(self.tuning.per_step_byte_limit), tail);
        self.pipeline.append_step(
            DecodeRecords::new(group_key_config, self.tuning.per_step_byte_limit),
            tail,
        )
    }

    /// Append the BAM decode preamble followed by `GroupByQueryname`:
    /// `ReadBgzfBlocks → BgzfDecompress → FindBamBoundaries → DecodeRecords → GroupByQueryname`.
    fn build_bam_decode_then_group_preamble(
        &mut self,
        reader: Box<dyn std::io::Read + Send>,
        header: Header,
        group_key_config: fgumi_bam_io::GroupKeyConfig,
    ) -> (crate::pipeline::core::topology::StepIdx, crate::pipeline::core::topology::BranchIdx)
    {
        use crate::pipeline::steps::group::queryname::GroupByQueryname;

        let tail = self.build_bam_decode_preamble(reader, header, group_key_config);
        self.pipeline.append_step(GroupByQueryname::new(self.tuning.per_step_byte_limit), tail)
    }

    /// Append the 2-step SAM parse preamble:
    /// `ReadSamChunks → ParseSamChunk`.
    fn build_sam_parse_preamble(
        &mut self,
        buf: Box<dyn std::io::BufRead + Send>,
        header: Header,
        group_key_config: fgumi_bam_io::GroupKeyConfig,
    ) -> (crate::pipeline::core::topology::StepIdx, crate::pipeline::core::topology::BranchIdx)
    {
        use crate::pipeline::steps::parse::sam::ParseSamChunk;
        use crate::pipeline::steps::source::read_sam_chunks::{
            DEFAULT_SAM_CHUNK_BYTES, ReadSamChunks,
        };

        let read_step =
            ReadSamChunks::new(buf, DEFAULT_SAM_CHUNK_BYTES, self.tuning.per_step_byte_limit);
        let parse_step =
            ParseSamChunk::new(Arc::new(header), group_key_config, self.tuning.per_step_byte_limit);
        let tail = self.pipeline.append_source(read_step);
        self.pipeline.append_step(parse_step, tail)
    }

    /// Append the SAM parse preamble followed by `GroupByQueryname`:
    /// `ReadSamChunks → ParseSamChunk → GroupByQueryname`.
    fn build_sam_parse_then_group_preamble(
        &mut self,
        buf: Box<dyn std::io::BufRead + Send>,
        header: Header,
        group_key_config: fgumi_bam_io::GroupKeyConfig,
    ) -> (crate::pipeline::core::topology::StepIdx, crate::pipeline::core::topology::BranchIdx)
    {
        use crate::pipeline::steps::group::queryname::GroupByQueryname;

        let tail = self.build_sam_parse_preamble(buf, header, group_key_config);
        self.pipeline.append_step(GroupByQueryname::new(self.tuning.per_step_byte_limit), tail)
    }

    /// Dispatch by `Stage` variant to the per-stage internal builder.
    ///
    /// `position` controls whether the stage appends a serialise-to-bytes step
    /// (Terminal) or leaves its native typed output on the chain tail
    /// (Intermediate). See [`StagePosition`] for the contract.
    ///
    /// # Errors
    ///
    /// Returns errors from step construction or pipeline building.
    pub fn add_stage(&mut self, stage: Stage, position: StagePosition) -> Result<()> {
        match stage {
            Stage::Dedup => self.add_dedup(position),
            Stage::Filter => self.add_filter(position),
            Stage::Clip => self.add_clip(position),
            Stage::Sort => self.add_sort(position),
            Stage::Group => self.add_group(position),
            Stage::Simplex => self.add_simplex(position),
            Stage::Duplex => self.add_duplex(position),
            Stage::Codec => self.add_codec(position),
            Stage::Correct => self.add_correct(position),
            Stage::Zipper => self.add_zipper(position),
            Stage::Align => self.add_align(position),
            Stage::Extract => self.add_extract(position),
            Stage::Downsample => {
                Err(anyhow!("Stage::Downsample is out of scope — not on the typed-step pipeline"))
            }
        }
    }

    /// Add the output sink step(s) to the pipeline.
    ///
    /// Appends `BgzfCompress` + `WriteBgzfFile`. Reads `spec.sink` for the
    /// output path and uses the output header from `self.header`.
    ///
    /// # Errors
    ///
    /// Returns errors from sink step construction (e.g., cannot create output file).
    ///
    /// # Panics
    ///
    /// Panics if called before [`Self::add_source`].
    pub fn add_sink(&mut self) -> Result<()> {
        use crate::pipeline::steps::bgzf::compress::BgzfCompress;
        use crate::pipeline::steps::sink::write_bgzf::WriteBgzfFile;

        let tail = self.current_tail.expect("add_sink called before add_source");
        let output_path = self.spec.sink.path();

        let compress_step =
            BgzfCompress::new(self.tuning.compression_level, self.tuning.per_step_byte_limit);
        let tail = self.pipeline.append_step(compress_step, tail);

        // When Stage::Align is in the chain, the merged output header is
        // resolved at pipeline-run time (the aligner contributes @PG/@RG/@CO
        // lines that aren't known until it emits its SAM/BAM header). The
        // `HeaderHandle` stashed by `add_align` must be forwarded here so the
        // writer blocks until the handle is resolved before emitting any bytes.
        let write_step = if let Some(handle) = self.pending_header_handle.take() {
            WriteBgzfFile::new_with_handle(output_path, handle, self.tuning.compression_level)
                .map_err(|e| anyhow!("WriteBgzfFile::new_with_handle: {e}"))?
        } else {
            WriteBgzfFile::new(output_path, &self.header, self.tuning.compression_level)
                .map_err(|e| anyhow!("WriteBgzfFile::new: {e}"))?
        };
        let tail = self.pipeline.append_step(write_step, tail);

        self.current_tail = Some(tail);
        Ok(())
    }

    /// Replace `self.header` with `header`, invalidating any pending
    /// [`HeaderHandle`] stashed by [`Self::add_align`].
    ///
    /// `add_align` stashes a `HeaderHandle` so `add_sink` can block on the
    /// aligner's runtime-resolved mapped-read header. When a *later* stage
    /// (e.g. an intermediate `add_align` followed by sort, consensus, or clip)
    /// builds a new output header, that stale handle would otherwise win in
    /// `add_sink` and the wrong header (the aligner's, missing the downstream
    /// stage's `@HD`/`@PG`/`@CO` updates) would be written. Routing every
    /// post-construction `self.header` reassignment through this helper drops
    /// the stale handle so `add_sink` falls back to the up-to-date
    /// `self.header`. Stages that leave the header unchanged (e.g. `add_group`)
    /// do not call this, preserving the handle for `align → group` chains.
    ///
    /// [`HeaderHandle`]: crate::pipeline::core::header::HeaderHandle
    fn replace_header(&mut self, header: noodles::sam::Header) {
        self.header = header;
        self.pending_header_handle = None;
    }

    /// Build the final [`BuiltPipeline`] from the accumulated state.
    ///
    /// Prepends a chain-level [`StageTimingFinalizeHook`] whose `Instant` is
    /// captured here — at the moment `build()` is called — so the reported
    /// wall time covers the full pipeline run, not the chain-construction
    /// phase. The label is derived from `spec.stages` formatted as
    /// `"dedup"` (single stage) or `"group→simplex"` (multi-stage), making
    /// the log line unambiguous for fused chains.
    ///
    /// Also calls [`build_pipeline_config_for_chain`] for shared threading /
    /// deadlock-timeout / queue-memory wiring. Registers
    /// [`PipelineStatsFinalizeHook`] if the user requested stats.
    ///
    /// # Errors
    ///
    /// Returns errors from pipeline configuration, memory-limit calculation,
    /// or if the pipeline has unwired output branches (indicates a programming
    /// error in an `add_*` method).
    pub fn build(mut self) -> Result<BuiltPipeline> {
        let num_threads = self.spec.threading.num_threads();
        let pipeline = self.pipeline.build().map_err(|e| anyhow!("Pipeline::build: {e:?}"))?;
        let (mut config, pipeline_stats) = build_pipeline_config_for_chain(
            &pipeline,
            num_threads,
            &self.spec.scheduler,
            &self.spec.queue_memory,
        )?;
        // Apply per-stage thread override, if set (currently only sort).
        if let Some(t) = self.override_pipeline_threads {
            config.threads = t;
        }
        let user_wants_stats = self.spec.scheduler.collect_stats();
        if user_wants_stats {
            if let Some(s) = pipeline_stats {
                self.finalize.push(Box::new(PipelineStatsFinalizeHook { stats: s }));
            }
        }

        // Prepend a single chain-level timing hook. Capturing Instant::now()
        // here (at build() time, not at add_<stage>() time) means the timer
        // starts just before Pipeline::run() is called, so elapsed ≈ pipeline
        // wall time. For multi-stage chains (e.g. [Group, Simplex]) we get
        // one timing line ("Pipeline `group→simplex` ran in Xs") rather than
        // one per stage — each showing the full pipeline duration.
        let chain_label = self
            .spec
            .stages
            .iter()
            .map(|s| format!("{s:?}").to_lowercase())
            .collect::<Vec<_>>()
            .join("→");
        // Box<dyn FinalizeHook> is not Clone, so we prepend by inserting at 0.
        // The timing hook fires first (before per-stage hooks like DedupFinalize).
        self.finalize.insert(0, Box::new(StageTimingFinalizeHook::new_with_label(chain_label)));

        Ok(BuiltPipeline { pipeline, config, finalize: self.finalize })
    }

    // -- private per-stage methods, one per Stage variant --
    // Each method's body fills in starting in 3a.3 (dedup) and is
    // mechanically extended for the other 9 stages in 3a.4-3a.12.

    /// Extract-specific step sequence:
    ///
    /// ```text
    /// (source: ReadFastqInputs → ZipFastqRecords)
    ///     ↓ FastqTemplateBatch
    /// ExtractStep (Parallel) — via build_extract_step
    ///     ↓ BamTemplateBatch
    /// (Terminal) SerializeBamRecords
    ///     ↓ DecompressedBlock
    /// ```
    ///
    /// For [`StagePosition::Terminal`], `SerializeBamRecords` is appended and
    /// the chain tail is `DecompressedBlock` (bytes ready for
    /// `BgzfCompress → WriteBgzfFile`).
    ///
    /// For [`StagePosition::Intermediate`], the tail is left as
    /// `BamTemplateBatch` so the next stage can consume it directly.
    ///
    /// Registers an [`ExtractFinalizeHook`] for record count + timer
    /// logging.
    ///
    /// # Errors
    ///
    /// Returns errors if extract options are missing from the spec bag,
    /// if validation fails, or if the source is not `SourceSpec::Fastqs`.
    ///
    /// [`ExtractFinalizeHook`]: crate::pipeline::chains::commands::extract::ExtractFinalizeHook
    fn add_extract(&mut self, position: StagePosition) -> Result<()> {
        use std::sync::atomic::AtomicU64;

        use crate::logging::OperationTimer;
        use crate::pipeline::chains::commands::extract::ExtractFinalizeHook;
        use crate::pipeline::steps::extract::build_extract_step;
        use crate::pipeline::steps::serialize::SerializeBamRecords;

        let extract_opts = self
            .spec
            .stage_opts
            .extract
            .as_ref()
            .ok_or_else(|| anyhow!("Stage::Extract options missing from StageOptionsBag"))?;

        extract_opts.validate()?;

        let tail = self.current_tail.expect("add_extract called before add_source");

        let output_path = match &self.spec.sink {
            SinkSpec::Bam(p) | SinkSpec::BamWithIndex(p) => p.clone(),
        };

        let timer = OperationTimer::new("Extracting UMIs");

        let read_structures = match &self.spec.source {
            SourceSpec::Fastqs { read_structures, .. } => Arc::new(read_structures.clone()),
            other => bail!("add_extract requires SourceSpec::Fastqs, got {other:?}"),
        };

        let input_path = self.resolve_log_input_path();
        log::info!("Starting Extract");
        log::info!("Input: {}", input_path.display());
        log::info!("Output: {}", output_path.display());

        let records_emitted = Arc::new(AtomicU64::new(0));

        let step = build_extract_step(
            read_structures,
            Arc::new(extract_opts.clone()),
            Arc::clone(&records_emitted),
            self.tuning.per_step_byte_limit,
        );
        let tail = self.pipeline.append_step(step, tail);

        if position == StagePosition::Terminal {
            let tail = self
                .pipeline
                .append_step(SerializeBamRecords::new(self.tuning.per_step_byte_limit), tail);
            self.current_tail = Some(tail);
            self.chain_tail_kind = ChainTailKind::DecodedRecordBatch;
        } else {
            self.current_tail = Some(tail);
            self.chain_tail_kind = ChainTailKind::BamTemplateBatch;
        }

        self.finalize.push(Box::new(ExtractFinalizeHook { records_emitted, timer }));

        Ok(())
    }

    /// Correct-specific step sequence:
    ///
    /// `GroupByQueryname` →
    /// `correct_step_kept_only` (single output → `BamTemplateBatch`) OR
    /// `correct_step_with_rejects` (two outputs: branch 0 kept `BamTemplateBatch`,
    /// branch 1 rejects `DecompressedBlock`) →
    /// `SerializeBamRecords` — **Terminal only** (on branch 0 / kept branch).
    ///
    /// For [`StagePosition::Terminal`], `SerializeBamRecords` is appended and
    /// the chain tail is [`DecompressedBlock`] (bytes ready for `BgzfCompress`).
    ///
    /// For [`StagePosition::Intermediate`], correct returns `Err` as a guard —
    /// no Phase 3 runall combination requires intermediate correct.
    ///
    /// When `--rejects` is set, branch 1 of the correct step carries pre-framed
    /// `DecompressedBlock` bytes and is wired here directly to its own
    /// `BgzfCompress → WriteBgzfFile` sink (same dual-sink pattern as
    /// `add_filter`, T3a.4). Branch 0 (kept `BamTemplateBatch`) becomes
    /// `current_tail` for the `SerializeBamRecords` step and, after that, for
    /// `add_sink`.
    ///
    /// Applies `correct_opts.validate()` for semantic option checks (UMI source
    /// presence + numeric ranges). Loads UMI sequences and encodes them into
    /// [`EncodedUmiSet`] up front.
    ///
    /// The `rejects_header` is the **input header verbatim** (PR #332 contract:
    /// rejects are raw-input records, so they carry the input's `@HD` sort
    /// fields + RG/PG — not the consensus/output header). For `correct`,
    /// `self.header` is the input header (correct preserves input order and does
    /// not replace `self.header` with a consensus header), so it is used as-is.
    ///
    /// Registers a [`CorrectFinalizeHook`] for: per-thread accumulator reduce →
    /// `finalize_metrics` → summary + warn banners → `--min-corrected` ratio
    /// gate (bail!) → `timer.log_completion`.
    /// The chain-level [`StageTimingFinalizeHook`] is registered by
    /// [`Self::build`] (not here), so it captures the correct `Instant`.
    ///
    /// # Errors
    ///
    /// Returns errors if correct options are missing from the spec bag, if UMI
    /// sequence loading fails, or if `position` is `Intermediate` (not yet
    /// implemented).
    ///
    /// [`EncodedUmiSet`]: crate::commands::correct::EncodedUmiSet
    /// [`DecompressedBlock`]: crate::pipeline::steps::types::DecompressedBlock
    /// [`CorrectFinalizeHook`]: crate::pipeline::chains::commands::correct::CorrectFinalizeHook
    #[allow(clippy::too_many_lines)]
    fn add_correct(&mut self, position: StagePosition) -> Result<()> {
        use crate::commands::common::warn_unwired_pipeline_flags;
        use crate::commands::correct::{CollectedCorrectMetrics, CorrectUmis, EncodedUmiSet};
        use crate::logging::OperationTimer;
        use crate::per_thread_accumulator::PerThreadAccumulator;
        use crate::pipeline::chains::commands::correct::CorrectFinalizeHook;
        use crate::pipeline::steps::correct::{
            CorrectStepConfig, correct_step_kept_only, correct_step_with_rejects,
        };
        use crate::pipeline::steps::group::queryname::GroupByQueryname;
        use crate::pipeline::steps::serialize::SerializeBamRecords;
        use crate::sam::SamTag;
        use log::info;

        // `Intermediate` is used when correct feeds into align-and-merge
        // (the `--start-from correct --stop-after ≥ zipper` fused chain).
        // In that case we skip `SerializeBamRecords` and leave the tail at
        // the correct step's `BamTemplateBatch` output (branch 0) so
        // `add_align` can wire `GroupByQueryname → AlignAndMergeStep` next.

        let correct_opts = self
            .spec
            .stage_opts
            .correct
            .as_ref()
            .ok_or_else(|| anyhow!("Stage::Correct options missing from StageOptionsBag"))?;

        // Semantic option-level checks (UMI source presence + numeric ranges).
        correct_opts.validate()?;

        let tail = self.current_tail.expect("add_correct called before add_source");

        // Resolve source path for log messages only.
        let input_path = self.resolve_log_input_path();
        let output_path = match &self.spec.sink {
            SinkSpec::Bam(p) | SinkSpec::BamWithIndex(p) => p.clone(),
        };

        let timer = OperationTimer::new("Correcting UMIs (new pipeline)");

        info!("{}", crate::commands::correct::NEW_PIPELINE_START_LOG);
        info!("Input: {}", input_path.display());
        info!("Output: {}", output_path.display());

        // Build a minimal CorrectUmis wrapper used by finalize_metrics.
        // Only `self.options.metrics` is accessed inside finalize_metrics —
        // the other fields are populated from spec for completeness.
        let correct_wrapper = CorrectUmis {
            io: crate::commands::common::BamIoOptions {
                input: input_path,
                output: output_path,
                ..Default::default()
            },
            rejects_opts: crate::commands::common::RejectsOptions::default(),
            options: correct_opts.clone(),
            threading: self.spec.threading.clone(),
            compression: self.spec.compression.clone(),
            scheduler_opts: self.spec.scheduler.clone(),
            queue_memory: self.spec.queue_memory.clone(),
        };

        let (umi_sequences, umi_length) = correct_wrapper.load_umi_sequences()?;
        let encoded_umi_set = Arc::new(EncodedUmiSet::new(&umi_sequences));
        correct_wrapper.check_umi_distances(&umi_sequences);
        let unmatched_umi = "N".repeat(umi_length);

        let num_threads = self.spec.threading.num_threads();

        // Build per-thread metrics accumulators and records-emitted counter.
        // Note: PerThreadAccumulator::new() already returns Arc<Self>.
        let metrics_slots = num_threads.max(1);
        let metrics = PerThreadAccumulator::<CollectedCorrectMetrics>::new(metrics_slots);
        let records_emitted = Arc::new(std::sync::atomic::AtomicU64::new(0));

        let umi_tag_bytes: [u8; 2] = SamTag::RX.into();
        let cfg = CorrectStepConfig {
            encoded_umi_set: Arc::clone(&encoded_umi_set),
            umi_length,
            umi_tag: umi_tag_bytes,
            opts: correct_opts.clone(),
            metrics: Arc::clone(&metrics),
            records_emitted: Arc::clone(&records_emitted),
            output_byte_limit: self.tuning.per_step_byte_limit,
            unmatched_umi: unmatched_umi.clone(),
        };

        let track_rejects = correct_opts.rejects_path.is_some();

        warn_unwired_pipeline_flags(&self.spec.scheduler, &self.spec.queue_memory);
        info!("{}", self.spec.threading.log_message());
        info!("Using pipeline with {num_threads} threads");

        // Wire GroupByQueryname before the correct step.
        let tail =
            self.pipeline.append_step(GroupByQueryname::new(self.tuning.per_step_byte_limit), tail);

        // Dispatch on rejects presence: either a 2-output or a 1-output step.
        let process_tail = if track_rejects {
            use crate::pipeline::core::topology::BranchIdx;
            use crate::pipeline::steps::bgzf::compress::BgzfCompress;
            use crate::pipeline::steps::sink::write_bgzf::WriteBgzfFile;

            let rejects_path = correct_opts.rejects_path.as_ref().ok_or_else(|| {
                anyhow!("rejects_path unexpectedly None (build_for dispatch bug)")
            })?;
            // PR #332 contract: rejects carry the INPUT header verbatim (raw-input
            // records, input order). For correct, `self.header` IS the input header.
            let rejects_header = self.header.clone();
            let rejects_write =
                WriteBgzfFile::new(rejects_path, &rejects_header, self.tuning.compression_level)
                    .map_err(|e| anyhow!("WriteBgzfFile (rejects): {e}"))?;

            let step = correct_step_with_rejects(cfg);
            let pt = self.pipeline.append_step(step, tail);

            // Branch 1 = pre-framed rejects bytes → BgzfCompress → WriteBgzfFile.
            let rejects_branch = (pt.0, BranchIdx(1));
            let rejects_compress_tail = self.pipeline.append_step(
                BgzfCompress::new(self.tuning.compression_level, self.tuning.per_step_byte_limit),
                rejects_branch,
            );
            let _rejects_write_tail =
                self.pipeline.append_step(rejects_write, rejects_compress_tail);

            // Branch 0 = kept BamTemplateBatch (needs SerializeBamRecords).
            pt
        } else {
            let step = correct_step_kept_only(cfg);
            self.pipeline.append_step(step, tail)
        };
        // process_tail = (process_step_idx, BranchIdx(0)) = kept BamTemplateBatch.

        if position == StagePosition::Terminal {
            // Append SerializeBamRecords on the Terminal path so the tail is
            // DecompressedBlock, ready for BgzfCompress → WriteBgzfFile.
            let tail = self.pipeline.append_step(
                SerializeBamRecords::new(self.tuning.per_step_byte_limit),
                process_tail,
            );
            self.current_tail = Some(tail);
            self.chain_tail_kind = ChainTailKind::DecodedRecordBatch; // actually DecompressedBlock
        } else {
            // Intermediate: leave the tail as BamTemplateBatch so the next
            // stage (add_align → GroupByQueryname → AlignAndMergeStep) can
            // consume it directly.
            self.current_tail = Some(process_tail);
            self.chain_tail_kind = ChainTailKind::BamTemplateBatch;
        }

        // Register the correct finalize hook.
        // The chain-level StageTimingFinalizeHook is inserted at index 0 by
        // build() after all stages have been added — that way a single timer
        // covers the full pipeline run regardless of how many stages there are.
        self.finalize.push(Box::new(CorrectFinalizeHook {
            metrics,
            records_emitted,
            encoded_umi_set,
            unmatched_umi,
            metrics_path: correct_opts.metrics.clone(),
            min_corrected: correct_opts.min_corrected,
            correct: correct_wrapper,
            timer,
        }));

        Ok(())
    }

    /// AlignAndMerge-specific step sequence (subprocess source):
    ///
    /// ```text
    /// (source preamble: ReadBgzfBlocks → BgzfDecompress → FindBamBoundaries → DecodeRecords)
    ///     ↓
    /// GroupByQueryname  ← appended here
    ///     ↓
    /// AlignAndMergeStep ← appended here (spawns aligner subprocess)
    ///     ↓
    /// (Terminal) SerializeBamRecords  ← appended only when Terminal
    ///     ↓
    /// (next stage consumes BamTemplateBatch — Intermediate only)
    /// ```
    ///
    /// Both [`StagePosition::Terminal`] and [`StagePosition::Intermediate`]
    /// are supported:
    ///
    /// - **Terminal** — `SerializeBamRecords` is appended after the AAM step
    ///   so the chain tail is `DecompressedBlock` ready for
    ///   `BgzfCompress → WriteBgzfFile`. This is the `--stop-after zipper`
    ///   path when `--start-from align-and-merge`: the merged BAM is
    ///   written directly to the output file.
    /// - **Intermediate** — the tail is left as `BamTemplateBatch` so the
    ///   next stage (Sort, Group, Consensus) can wire the correct preamble.
    ///
    /// ## Header handling
    ///
    /// The aligner contributes `@PG`/`@RG`/`@CO` lines to the output
    /// BAM header at runtime (not at chain-construction time). The
    /// `AlignAndMergeStep` resolves a [`HeaderHandle`] once the aligner
    /// emits its SAM/BAM header; the downstream [`WriteBgzfFile`]
    /// (added by [`Self::add_sink`]) blocks on the handle before
    /// writing any record bytes.
    ///
    /// To prepare the correct partial header (dict-derived `@SQ` +
    /// unmapped-BAM headers + fgumi `@PG`) for `AlignAndMergeConfig`,
    /// `add_align` calls [`build_output_header`] with the current
    /// `self.header` (which was built from the unmapped BAM by
    /// `ChainBuilder::new`) and the reference FASTA's `.dict` path.
    /// The result replaces `self.header` so downstream stages and
    /// `add_sink` see the dict-derived `@SQ` table.
    ///
    /// The `HeaderHandle` is stashed in `self.pending_header_handle` for
    /// `add_sink` to consume.
    ///
    /// ## Thread floor
    ///
    /// AAM spawns two daemon threads (FASTQ writer + SAM/BAM reader) plus
    /// the aligner subprocess. The pipeline needs at least 4 framework
    /// workers for steady-state progress:
    ///
    /// - 1 worker to drive the source preamble
    /// - 1 worker dispatching `AlignAndMergeStep` (Serial)
    /// - 1+ workers for any downstream Parallel steps
    /// - 1 spare for progress / bookkeeping
    ///
    /// `override_pipeline_threads` is raised to `threads.max(4)` so
    /// `build()` forwards the correct value to `PipelineConfig::threads`.
    ///
    /// ## Errors
    ///
    /// Returns errors if:
    /// - align options are missing from the spec bag,
    /// - the reference `.dict` file cannot be found,
    /// - `AlignAndMergeStep::new` fails to spawn the aligner subprocess.
    ///
    /// [`HeaderHandle`]: crate::pipeline::core::header::HeaderHandle
    /// [`WriteBgzfFile`]: crate::pipeline::steps::sink::write_bgzf::WriteBgzfFile
    /// [`build_output_header`]: crate::commands::zipper::build_output_header
    fn add_align(&mut self, position: StagePosition) -> Result<()> {
        use std::sync::atomic::AtomicU64;

        use crate::commands::zipper::build_output_header;
        use crate::logging::OperationTimer;
        use crate::pipeline::chains::commands::align::AlignFinalizeHook;
        use crate::pipeline::core::header::HeaderHandle;
        use crate::pipeline::steps::align_and_merge::{AlignAndMergeConfig, AlignAndMergeStep};
        use crate::pipeline::steps::group::queryname::GroupByQueryname;
        use crate::pipeline::steps::serialize::SerializeBamRecords;
        use crate::reference::find_dict_path;
        use crate::sam::check_sort;
        use crate::umi::TagInfo;
        use log::info;

        let align_opts = self
            .spec
            .stage_opts
            .aligner
            .as_ref()
            .ok_or_else(|| anyhow!("Stage::Align options missing from StageOptionsBag"))?;

        let num_threads = self.spec.threading.num_threads();

        // Resolve the aligner command. This validates the option combination
        // (preset vs command, index files, binary path) before spawning.
        let resolved = align_opts.aligner.clone().resolve(
            &align_opts.reference,
            num_threads,
            align_opts.aligner_bin.as_deref(),
        )?;

        info!(
            "AlignAndMerge: mode = {:?}, chunk_size = {}, threads = {:?}",
            resolved.mode, resolved.chunk_size, resolved.threads,
        );
        info!("AlignAndMerge command: {}", resolved.command);

        // Validate the unmapped BAM source has queryname sort order — but only
        // when Align reads directly from the chain source (a file). When an
        // upstream stage already produced `BamTemplateBatch` (extract → align,
        // correct → align), Align's input is that in-pipeline tail — wired below
        // via `align_input_tail` — not `spec.source`. In that case `spec.source`
        // is the chain's *original* start input (e.g. `Fastqs` for
        // `--start-from extract`), which is neither a Bam/Sam file to sort-check
        // nor Align's actual input, so the check must be skipped. The upstream
        // stage's output is queryname-grouped by construction. This mirrors the
        // `chain_tail_kind == BamTemplateBatch` distinction used for
        // `align_input_tail` just below.
        if self.chain_tail_kind != ChainTailKind::BamTemplateBatch {
            let input_path = match &self.spec.source {
                SourceSpec::Bam(p) | SourceSpec::Sam(p) => p.clone(),
                other => {
                    bail!("Stage::Align requires SourceSpec::Bam (unmapped BAM); got {other:?}")
                }
            };
            check_sort(&self.header, &input_path, "unmapped");
        }

        // Build the correct partial output header:
        //   dict-derived @SQ + unmapped-BAM @HD/@CO/@RG/@PG + fgumi @PG.
        // `self.header` already has add_pg_record applied, so the fgumi PG
        // record is present; build_output_header merges the unmapped header's
        // PG/RG/CO lines with the dict @SQ and leaves the aligner's runtime
        // contributions to be folded in via HeaderHandle.
        let dict_path = find_dict_path(&align_opts.reference).ok_or_else(|| {
            anyhow!(
                "Reference dictionary file not found alongside '{}'. \
                 Run `samtools dict {} -o {}.dict` and try again.",
                align_opts.reference.display(),
                align_opts.reference.display(),
                align_opts.reference.display(),
            )
        })?;
        let partial_header =
            build_output_header(&self.header, &noodles::sam::Header::default(), &dict_path)?;
        let partial_header = std::sync::Arc::new(partial_header);

        // Replace self.header with the dict-corrected partial header so
        // downstream stages (sort, group) and add_sink see the correct @SQ.
        self.header = partial_header.as_ref().clone();

        let header_handle = HeaderHandle::new();
        let records_emitted = std::sync::Arc::new(AtomicU64::new(0));

        let aam_cfg = AlignAndMergeConfig {
            tag_info: std::sync::Arc::new(TagInfo::new(Vec::new(), Vec::new(), Vec::new())),
            skip_tc_tags: false,
            reference: None,
            partial_output_header: std::sync::Arc::clone(&partial_header),
            header_handle: header_handle.clone(),
            records_emitted: std::sync::Arc::clone(&records_emitted),
            output_byte_limit: self.tuning.per_step_byte_limit,
            // Bound AAM's in-flight unmapped backlog to the aligner's `-K`
            // chunk size so a fast-draining aligner can't accumulate the
            // whole input's unmapped reads in RAM (issue #382).
            in_flight_unmapped_budget:
                crate::pipeline::steps::align_and_merge::in_flight_budget_for_chunk_size(
                    resolved.chunk_size,
                ),
        };

        let aam_step = AlignAndMergeStep::new(aam_cfg, &resolved.command)
            .map_err(|e| anyhow!("AlignAndMergeStep::new: {e:#}"))?;

        let tail = self.current_tail.expect("add_align called before add_source");

        // Wire GroupByQueryname before the AAM step — but only when the upstream
        // stage emits DecodedRecordBatch (the normal source preamble path).
        // When `chain_tail_kind == BamTemplateBatch`, an upstream stage (e.g.
        // Correct) has already grouped records into BamTemplateBatch, which is
        // exactly the input type AlignAndMergeStep expects. Skip GroupByQueryname
        // in that case to avoid a DecodedRecordBatch→BamTemplateBatch type mismatch.
        let align_input_tail = if self.chain_tail_kind == ChainTailKind::BamTemplateBatch {
            // Upstream is already BamTemplateBatch — wire AAM directly.
            tail
        } else {
            // Upstream is DecodedRecordBatch (source preamble). Group by queryname
            // to produce BamTemplateBatch for AlignAndMergeStep.
            self.pipeline.append_step(GroupByQueryname::new(self.tuning.per_step_byte_limit), tail)
        };
        let aam_tail = self.pipeline.append_step(aam_step, align_input_tail);

        if position == StagePosition::Terminal {
            // Terminal: append SerializeBamRecords so the chain tail is
            // DecompressedBlock ready for BgzfCompress → WriteBgzfFile in add_sink.
            // This is the --stop-after zipper path: merged BAM written directly.
            let tail = self
                .pipeline
                .append_step(SerializeBamRecords::new(self.tuning.per_step_byte_limit), aam_tail);
            self.current_tail = Some(tail);
            // chain_tail_kind is DecompressedBlock after SerializeBamRecords, but
            // we reuse DecodedRecordBatch as the sentinel for "bytes ready for sink"
            // (consistent with other Terminal paths).
            self.chain_tail_kind = ChainTailKind::DecodedRecordBatch;
        } else {
            self.current_tail = Some(aam_tail);
            // AAM output is BamTemplateBatch; downstream stages (Sort, Group,
            // Consensus) must know this so they can wire the correct preamble.
            self.chain_tail_kind = ChainTailKind::BamTemplateBatch;
        }

        // Stash the HeaderHandle for add_sink to pick up.
        self.pending_header_handle = Some(header_handle);

        // AAM spawns two daemon threads + aligner subprocess; 4 workers
        // minimum for steady-state progress.
        let effective = num_threads.max(4);
        if let Some(prev) = self.override_pipeline_threads {
            self.override_pipeline_threads = Some(prev.max(effective));
        } else {
            self.override_pipeline_threads = Some(effective);
        }

        let timer = OperationTimer::new("AlignAndMerge");
        self.finalize.push(Box::new(AlignFinalizeHook { records_emitted, timer }));

        Ok(())
    }

    /// Zipper-specific step sequence (dual-input):
    ///
    /// Precondition: `add_source` has already been called with a
    /// `PendingSource::Paired` spec, which set `current_tail` to the
    /// unmapped source chain tail and `paired_tail` to the mapped
    /// source chain tail. Both chains end at `GroupByQueryname`, so
    /// their output type is `OrderedBytesSingle<BamTemplateBatch>`,
    /// matching `ZipperMergeStep: Step2<InputA = BamTemplateBatch,
    /// InputB = BamTemplateBatch>`.
    ///
    /// Appends `ZipperMergeStep` (via `PipelineBuilder::append_step2`
    /// which wires the two tails into input slots 0 and 1), then — for
    /// `Terminal` — appends `SerializeBamRecords` so the chain tail is
    /// `DecompressedBlock` ready for `BgzfCompress → WriteBgzfFile`.
    ///
    /// Applies the thread floor: zipper requires at least 4 worker
    /// threads (`≥ 3 Exclusive steps + 1`); `override_pipeline_threads`
    /// is set to `spec.threading.num_threads().max(4)` so that
    /// `build()` forwards the correct value to `PipelineConfig::threads`.
    ///
    /// Registers a [`ZipperFinalizeHook`] for the missing-reads summary,
    /// "zipper completed successfully" log, and `timer.log_completion`.
    /// The chain-level [`StageTimingFinalizeHook`] is registered by
    /// [`Self::build`], so a single timer covers the full pipeline run.
    ///
    /// # Errors
    ///
    /// Returns errors if zipper options are missing from the spec bag, if
    /// `add_source` has not been called with a `PairedBams` spec (missing
    /// `paired_tail`), or if `position` is `Intermediate` (not yet
    /// implemented).
    ///
    /// [`ZipperFinalizeHook`]: crate::pipeline::chains::commands::zipper::ZipperFinalizeHook
    fn add_zipper(&mut self, position: StagePosition) -> Result<()> {
        use crate::commands::common::warn_unwired_pipeline_flags;
        use crate::commands::zipper::NEW_PIPELINE_START_LOG;
        use crate::logging::OperationTimer;
        use crate::pipeline::chains::commands::zipper::{
            ZipperFinalizeHook, ZipperMergeCaptures, build_zipper_merge_step,
        };
        use crate::pipeline::steps::serialize::SerializeBamRecords;
        use log::info;

        let zipper_opts = self
            .spec
            .stage_opts
            .zipper
            .as_ref()
            .ok_or_else(|| anyhow!("Stage::Zipper options missing from StageOptionsBag"))?
            .clone();

        // Consume both source tails set by add_source for PairedBams.
        let unmapped_tail =
            self.current_tail.expect("add_zipper called before add_source (Paired)");
        let mapped_tail = self
            .paired_tail
            .take()
            .ok_or_else(|| anyhow!("add_zipper: paired_tail missing — spec must be PairedBams"))?;

        // Resolve the reference path from the spec (the reference is always
        // SourceSpec::PairedBams for zipper — validated by the cross-stage validator).
        let reference_path = match &self.spec.source {
            SourceSpec::PairedBams { reference, .. } => reference.clone(),
            other => bail!("Stage::Zipper requires SourceSpec::PairedBams, got {other:?}"),
        };

        // Apply the thread floor: zipper needs ≥ 4 workers (2 Exclusive source
        // readers + 1 Exclusive write step + at least 1 Serial worker).
        let raw_threads = self.spec.threading.num_threads();
        let num_threads = raw_threads.max(4);

        // Recompute tuning with the floored thread count so that batch-sizing
        // parameters (template_batch_size, per_step_byte_limit) are consistent
        // with the actual concurrency. This replaces self.tuning (which was
        // built from raw_threads in new()).
        let floored_tuning = BamPipelineTuning::auto_tuned(num_threads)
            .with_compression_level(self.spec.compression.compression_level);

        // Emit the start log and pipeline-info lines (mirrors build_zipper_chain).
        info!("{NEW_PIPELINE_START_LOG}");
        warn_unwired_pipeline_flags(&self.spec.scheduler, &self.spec.queue_memory);
        info!("{}", self.spec.threading.log_message());
        info!("Using pipeline with {num_threads} threads");

        let timer = OperationTimer::new("Zipping BAMs (new pipeline)");

        let missing_count = Arc::new(std::sync::atomic::AtomicU64::new(0));
        let records_emitted = Arc::new(std::sync::atomic::AtomicU64::new(0));

        let merge_step = build_zipper_merge_step(ZipperMergeCaptures {
            zipper_opts: zipper_opts.clone(),
            output_header: Arc::new(self.header.clone()),
            reference_path,
            tuning: floored_tuning,
            missing_count: Arc::clone(&missing_count),
            records_emitted: Arc::clone(&records_emitted),
        })?;

        // Wire both source chains into ZipperMergeStep via the 2-input append.
        let merge_tail = self.pipeline.append_step2(merge_step, unmapped_tail, mapped_tail);

        // Gate SerializeBamRecords on Terminal position. For Intermediate
        // (e.g., zipper → sort → ...), the chain tail stays at ZipperMergeStep's
        // BamTemplateBatch output, ready for the next stage to consume.
        if position == StagePosition::Terminal {
            let tail = self.pipeline.append_step(
                SerializeBamRecords::new(floored_tuning.per_step_byte_limit),
                merge_tail,
            );
            self.current_tail = Some(tail);
            // Tail is now DecompressedBlock (bytes). The chain_tail_kind
            // sentinel after SerializeBamRecords is DecodedRecordBatch
            // (same convention used by other Terminal paths — see the
            // comment on `ChainTailKind` for the "ready for sink" sentinel).
            self.chain_tail_kind = ChainTailKind::DecodedRecordBatch;
        } else {
            self.current_tail = Some(merge_tail);
            // Intermediate: ZipperMergeStep emits BamTemplateBatch so the
            // next stage (typically add_sort) prepends TemplatesToRecordBatch.
            self.chain_tail_kind = ChainTailKind::BamTemplateBatch;
        }

        // Set the thread override to the floored count so build() forwards it
        // to PipelineConfig::threads.
        self.override_pipeline_threads = Some(num_threads);

        // Register ZipperFinalizeHook.
        // The chain-level StageTimingFinalizeHook is inserted at index 0 by
        // build() after all stages have been added — that way a single timer
        // covers the full pipeline run regardless of how many stages there are.
        self.finalize.push(Box::new(ZipperFinalizeHook {
            missing_count,
            records_emitted,
            exclude_missing_reads: zipper_opts.exclude_missing_reads,
            timer,
        }));

        Ok(())
    }

    /// Sort-specific step sequence — two execution modes:
    ///
    /// ## Standalone mode (Sort is the sole stage)
    ///
    /// When [`StagePosition::Terminal`] and `current_tail` is `None`
    /// (`build_sort_chain` skips `add_source` / `add_sink`), sort emits
    /// a single `Exclusive` [`SortBamFile`] step that reads from and
    /// writes to files directly. `SortBamFile` has `Input = ()` and
    /// `Outputs = ()` — it is a self-contained black-box. The method:
    ///
    /// - Does **not** read from `current_tail`.
    /// - Registers the step via `PipelineBuilder::append_source`.
    /// - Leaves `current_tail = None` (`add_sink` must not be called afterwards).
    /// - Sets `override_pipeline_threads = Some(1)` so `build()` runs a single
    ///   framework worker; the sorter's internal `SortWorkerPool` / rayon pools
    ///   provide the real concurrency.
    /// - Registers a `SortFinalizeHook` for stats logging + timer.
    ///
    /// ## Streaming mode (sort composes with other stages)
    ///
    /// When `current_tail` is `Some` (the source preamble or an
    /// upstream stage has run), sort emits the in-pipeline step
    /// sequence the legacy fused dispatchers used:
    ///
    /// ```text
    /// [BamTemplateBatch or RecordBatch at current_tail]
    ///     ↓ (TemplatesToRecordBatch if BamTemplateBatch)
    /// RecordBatch
    ///     ↓ SortAndSpill (Serial)
    /// SortPhase1Event
    ///     ↓ SortSpillDecompress (Parallel)
    /// SortPhase2Event
    ///     ↓ SortMerge (Serial)
    /// RecordBatch
    ///     ↓ [Terminal]     SerializeRecordBatch → DecompressedBlock
    ///     ↓ [Intermediate] DecodeFromRecords    → DecodedRecordBatch
    /// ```
    ///
    /// Streaming mode does NOT register a `SortFinalizeHook`; the
    /// chain-level [`StageTimingFinalizeHook`] covers the wall time.
    /// Memory must be `Fixed` — `Auto` cannot be honoured when sort
    /// shares its budget with other in-pipeline stages.
    ///
    /// # Errors
    ///
    /// Returns errors if sort options are missing from the spec bag, if the
    /// source is not a BAM/SAM path, or if `--sort::max-memory=auto` is
    /// requested in streaming mode.
    ///
    /// [`SortBamFile`]: crate::pipeline::steps::sort::SortBamFile
    #[allow(clippy::too_many_lines)]
    fn add_sort(&mut self, position: StagePosition) -> Result<()> {
        use crate::commands::common::{MemoryLimit, resolve_memory_budget};
        use crate::pipeline::chains::commands::sort::{
            IndexBamFinalizeHook, SortFinalizeHook, SortStepCaptures, build_sort_step,
            log_sort_start, sort_sink_path, sort_source_path,
        };

        let sort = self
            .spec
            .stage_opts
            .sort
            .as_ref()
            .ok_or_else(|| anyhow!("Stage::Sort options missing from StageOptionsBag"))?
            .clone();

        if position == StagePosition::Terminal && self.current_tail.is_none() {
            // ── Standalone sort (Stage::Sort is the only stage) ─────────────
            //
            // SortBamFile is a self-contained Exclusive step that reads from
            // and writes to files directly. Used only when Sort is the sole
            // stage in the chain (current_tail is None because add_source was
            // skipped by build_for for the sort-terminal special case).

            let input_path = sort_source_path(self.spec)?;
            let output_path = sort_sink_path(self.spec);

            let num_sorter_threads = self.spec.threading.num_threads();

            let effective_memory = resolve_memory_budget(
                sort.max_memory,
                sort.memory_reserve,
                num_sorter_threads,
                sort.memory_per_thread,
            )?;

            let timer = crate::logging::OperationTimer::new("Sorting BAM");

            log_sort_start(&sort, &input_path, &output_path, num_sorter_threads, effective_memory);

            let stats_slot = Arc::new(parking_lot::Mutex::new(None::<fgumi_sort::SortStats>));

            let sort_step = build_sort_step(SortStepCaptures {
                sort,
                input_path: input_path.clone(),
                output_path: output_path.clone(),
                num_sorter_threads,
                output_compression: self.spec.compression.compression_level,
                command_line: self.spec.command_line.clone(),
                stats_slot: Arc::clone(&stats_slot),
            })?;

            // Register SortBamFile as a source (Input = (), Outputs = ()).
            // PipelineBuilder::build() will not flag unwired branches since
            // Outputs::arity() == 0.
            let tail = self.pipeline.append_source(sort_step);
            // current_tail must remain None / unchanged — add_sink must NOT be
            // called after add_sort (Terminal).  We don't update current_tail
            // here so that any mistaken add_sink() call will panic on
            // "add_sink called before add_source" rather than silently adding
            // orphaned steps.
            let _ = tail; // tail is unused intentionally

            // Single Exclusive step → pipeline needs only 1 framework thread.
            self.override_pipeline_threads = Some(1);

            self.finalize.push(Box::new(SortFinalizeHook {
                stats_slot,
                output_path: output_path.clone(),
                timer,
            }));

            // When the spec asks for a sidecar BAI, queue the
            // `IndexBamFinalizeHook` *after* `SortFinalizeHook` so the
            // "Records written / X records/s" summary lands before the
            // indexer's "Indexing BAM:" / "Wrote BAM index:" pair. The
            // hook re-reads the finished BAM and emits
            // `<output>.bam.bai` next to it (see `IndexBamFinalizeHook`
            // for the I/O contract). Rule 3 in `chains::validate`
            // guarantees `BamWithIndex` only appears when the terminal
            // chain stage is `Stage::Sort`, so a sole-stage
            // `[Stage::Sort]` chain is the only path that reaches the
            // Standalone branch with that sink — multi-stage
            // sort-terminal chains land in the streaming branch below,
            // where the hook is mirrored.
            if matches!(self.spec.sink, SinkSpec::BamWithIndex(_)) {
                self.finalize.push(Box::new(IndexBamFinalizeHook { output_path }));
            }
        } else {
            // ── Streaming sort (Intermediate OR Terminal-after-upstream-stages) ──
            //
            // Used when Sort has upstream pipeline stages (current_tail is Some),
            // regardless of whether Sort is Intermediate or Terminal.
            //
            // Chain topology (continuing from current_tail):
            //
            //   [BamTemplateBatch or RecordBatch at current_tail]
            //       ↓ (TemplatesToRecordBatch if BamTemplateBatch)
            //   RecordBatch
            //       ↓
            //   SortAndSpill  (Serial) → SortPhase1Event
            //       ↓
            //   SortSpillDecompress (Parallel)
            //       ↓ SortPhase2Event
            //   SortMerge (Serial) → RecordBatch
            //       ↓
            //   [Intermediate] DecodeFromRecords (Parallel) → DecodedRecordBatch
            //   [Terminal]     SerializeRecordBatch (Parallel) → DecompressedBlock
            //
            // For Intermediate: chain_tail_kind = DecodedRecordBatch;
            // add_group / add_simplex / etc. can proceed normally.
            //
            // For Terminal: chain_tail_kind = DecodedRecordBatch (sentinel for
            // "bytes ready for sink") — add_sink appends BgzfCompress → WriteBgzfFile.
            //
            // Memory must be `Fixed` — `Auto` cannot be honoured in a
            // multi-stage pipeline where the memory budget is shared across
            // stages. This matches `build_sort_steps`' behaviour.
            use crate::pipeline::core::step::Affinity;
            use crate::pipeline::steps::parse::decode::DecodeFromRecords;
            use crate::pipeline::steps::serialize_record_batch::SerializeRecordBatch;
            use crate::pipeline::steps::sort::{SortAndSpill, SortMerge, SortSpillDecompress};
            use crate::sam::SamTag;
            use fgumi_bam_io::GroupKeyConfig;
            use fgumi_bam_io::LibraryIndex;
            use fgumi_sort::{RawExternalSorter, SortOrder};
            use noodles::sam::alignment::record::data::field::Tag;

            // Auto isn't supported in a fused chain (sort shares the budget with
            // the other stages); require a fixed value.
            if matches!(sort.max_memory, MemoryLimit::Auto) {
                bail!(
                    "--sort::max-memory=auto is not supported in a fused multi-stage pipeline \
                     (sort is not standalone); pass a fixed value (e.g. --sort::max-memory 768M)"
                );
            }
            let num_threads = self.spec.threading.num_threads();
            // Reuse the same helper as the standalone path so `--sort::memory-per-thread`
            // is honored (the previous inline `× num_threads` ignored it).
            let total_memory = resolve_memory_budget(
                sort.max_memory,
                sort.memory_reserve,
                num_threads,
                sort.memory_per_thread,
            )?;

            let sort_order = SortOrder::TemplateCoordinate;
            // NB: the spill codec is pinned to BGZF inside `build_stream`
            // (`SortSpillDecompress` is BGZF-block-only); see the comment there.
            let mut sorter = RawExternalSorter::new(sort_order)
                .memory_limit(total_memory)
                .threads(num_threads.max(1))
                .output_compression(1)
                .temp_compression(sort.temp_compression)
                .cell_tag(SamTag::CB);

            if !sort.tmp_dirs.is_empty() {
                sorter = sorter.temp_dirs(sort.tmp_dirs.clone());
            }

            let affinity =
                if num_threads.max(2) >= 3 { Affinity::Worker(1) } else { Affinity::Reader };

            let and_spill = SortAndSpill::from_sorter(sorter, &self.header, 64)
                .map_err(|e| anyhow!("SortAndSpill::from_sorter: {e}"))?
                .with_affinity(affinity);
            let decompress = SortSpillDecompress::new(64);
            let merge =
                SortMerge::new(sort_order, self.tuning.per_step_byte_limit).with_affinity(affinity);

            let tail = self.current_tail.expect("streaming sort called before add_source");

            // If the upstream stage produced BamTemplateBatch (from Align or
            // Correct), flatten it to RecordBatch before feeding SortAndSpill.
            let tail = if self.chain_tail_kind == ChainTailKind::BamTemplateBatch {
                use crate::pipeline::steps::templates_to_records::TemplatesToRecordBatch;
                self.pipeline
                    .append_step(TemplatesToRecordBatch::new(self.tuning.per_step_byte_limit), tail)
            } else {
                // chain_tail_kind == RecordBatch (from ParseBamRecords in add_source)
                // or DecodedRecordBatch would be a bug caught at runtime.
                tail
            };

            let tail = self.pipeline.append_step(and_spill, tail);
            let tail = self.pipeline.append_step(decompress, tail);
            let merge_tail = self.pipeline.append_step(merge, tail);

            if position == StagePosition::Terminal {
                // Terminal-after-upstream: sort then write to output file.
                // SerializeRecordBatch converts RecordBatch → DecompressedBlock
                // ready for BgzfCompress → WriteBgzfFile in add_sink.
                let tail = self.pipeline.append_step(
                    SerializeRecordBatch::new(self.tuning.per_step_byte_limit),
                    merge_tail,
                );
                self.current_tail = Some(tail);
                // Use DecodedRecordBatch as the sentinel for "bytes ready for sink"
                // (consistent with other Terminal paths).
                self.chain_tail_kind = ChainTailKind::DecodedRecordBatch;

                // Mirror the Standalone-branch BAI hook registration: when
                // the spec asks for a sidecar BAI, queue the
                // `IndexBamFinalizeHook` so the indexer runs after the chain
                // has flushed the final BAM. Multi-stage Sort-terminal
                // chains (e.g. `[Stage::Correct, Stage::Sort]`) land here
                // instead of the Standalone branch (which requires
                // `current_tail.is_none()`), so without this wire a
                // `BamWithIndex` request would pass Rule 3 and silently
                // skip indexing. Today only standalone `fgumi sort` builds
                // such specs (runall never constructs `SinkSpec::BamWithIndex`
                // because its sort output is template-coordinate, where BAI
                // is undefined), so the Streaming-Terminal path is unreachable
                // in production, but the architectural invariant —
                // "BamWithIndex triggers the hook in every terminal-sort
                // path" — is enforced here for future producers.
                if let SinkSpec::BamWithIndex(p) = &self.spec.sink {
                    self.finalize.push(Box::new(
                        crate::pipeline::chains::commands::sort::IndexBamFinalizeHook {
                            output_path: p.clone(),
                        },
                    ));
                }
            } else {
                // Intermediate: decode to DecodedRecordBatch for downstream stages.
                let cell_tag = Tag::from(SamTag::CB);
                let library_index = LibraryIndex::from_header(&self.header);
                let group_key_config = GroupKeyConfig::new(library_index, cell_tag);
                let tail = self.pipeline.append_step(
                    DecodeFromRecords::new(group_key_config, self.tuning.per_step_byte_limit),
                    merge_tail,
                );
                self.current_tail = Some(tail);
                self.chain_tail_kind = ChainTailKind::DecodedRecordBatch;
            }

            // Update self.header to reflect template-coordinate sort order so
            // downstream add_group's sort-order header check passes. This is
            // the same header update that the external sorter writes into the
            // temp BAM — here we apply it directly so the static check at
            // chain-construction time sees the correct sort order.
            self.replace_header(fgumi_sort::create_output_header(sort_order, &self.header));

            // Streaming sort does NOT register a SortFinalizeHook
            // (no output-path stats to log; the overall chain timing hook
            // registered by build() covers the full pipeline wall time).
        }

        Ok(())
    }

    /// Group-specific step sequence:
    ///
    /// `GroupByPosition` →
    /// `ProcessOrdered` (filter + UMI assign + metrics accumulation) →
    /// `MiAssign` (serial MI counter) →
    /// `ProcessWithWorkerState` (serialize records with MI tags) — **Terminal only**.
    ///
    /// For [`StagePosition::Terminal`], all four steps are appended and the
    /// chain tail is [`DecompressedBlock`] (bytes ready for `BgzfCompress`).
    ///
    /// For [`StagePosition::Intermediate`], only the first three steps are
    /// appended; the chain tail stays as [`BatchedProcessedPositionGroups`] for
    /// the next stage to consume. Intermediate group is not yet needed by any
    /// Phase 3a runall combination; calling `add_group` with `Intermediate`
    /// returns `Err` as a guard until Phase 3b's fused group→consensus chains
    /// require it.
    ///
    /// Accepts both template-coordinate-sorted and (with `--allow-unmapped`)
    /// queryname-sorted inputs, matching `GroupReadsByUmi::execute`'s validation.
    ///
    /// Registers a `GroupFinalizeHook` for accumulators reduce + family-size
    /// histogram + grouping-metrics + summary banner + timer.
    /// The chain-level [`StageTimingFinalizeHook`] is registered by
    /// [`Self::build`] (not here), so it captures the correct `Instant`.
    ///
    /// # Errors
    ///
    /// Returns errors if the sort order is wrong, if the group options are missing
    /// from the spec bag, or if `position` is `Intermediate` (not yet implemented).
    ///
    /// [`BatchedProcessedPositionGroups`]: crate::pipeline::steps::group::position::BatchedProcessedPositionGroups
    #[allow(clippy::too_many_lines)]
    fn add_group(&mut self, position: StagePosition) -> Result<()> {
        use crate::commands::common::warn_unwired_pipeline_flags;
        use crate::commands::group::GroupFilterConfig;
        use crate::logging::OperationTimer;
        use crate::per_thread_accumulator::PerThreadAccumulator;
        use crate::pipeline::chains::commands::group::{
            GroupFinalizeHook, build_group_mi_assign_step, build_group_process_step,
            build_group_serialize_step,
        };
        use crate::pipeline::steps::group::position::GroupByPosition;
        use crate::sam::{SamTag, is_sorted, is_template_coordinate_sorted};
        use log::{info, warn};
        use noodles::sam::header::record::value::map::header::sort_order::QUERY_NAME;

        let group = self
            .spec
            .stage_opts
            .group
            .as_ref()
            .ok_or_else(|| anyhow!("Stage::Group options missing from StageOptionsBag"))?;

        let tail = self.current_tail.expect("add_group called before add_source");

        // Resolve source path for log messages only.
        let input_path = self.resolve_log_input_path();
        let output_path = match &self.spec.sink {
            SinkSpec::Bam(p) | SinkSpec::BamWithIndex(p) => p.clone(),
        };

        let timer = OperationTimer::new("Grouping reads by UMI");

        info!("Starting group");
        info!("Input: {}", input_path.display());
        info!("Output: {}", output_path.display());
        info!("Strategy: {:?}", group.effective_strategy);
        info!("Edits: {}", group.effective_edits);
        if group.no_umi {
            info!("No-UMI mode: grouping by position only");
        }
        if matches!(
            group.effective_strategy,
            crate::assigner::Strategy::Adjacency | crate::assigner::Strategy::Paired
        ) {
            info!("Index threshold: {}", group.index_threshold);
        }
        if group.allow_unmapped {
            info!("Allow unmapped: enabled (unmapped templates will be grouped by UMI only)");
            warn!(
                "WARNING: All unmapped reads are placed in a single position group. \
                 Reads with identical/similar UMIs will be grouped together even if they \
                 originate from different genomic locations."
            );
            if matches!(
                group.strategy,
                crate::assigner::Strategy::Edit
                    | crate::assigner::Strategy::Adjacency
                    | crate::assigner::Strategy::Paired
            ) {
                warn!(
                    "WARNING: For paired UMIs (e.g., ACGT-TGCA), edit distance is computed \
                     on the concatenated sequence with dashes removed. With --edits {}, \
                     only {} mismatch(es) allowed across ALL bases.",
                    group.edits, group.edits
                );
            }
        }

        warn_unwired_pipeline_flags(&self.spec.scheduler, &self.spec.queue_memory);
        let num_threads = self.spec.threading.num_threads();
        info!("{}", self.spec.threading.log_message());
        info!("Using pipeline with {num_threads} threads");

        // Check sort order — identical validation to GroupReadsByUmi::execute.
        let is_tc_sorted = is_template_coordinate_sorted(&self.header);
        let is_qname_sorted = is_sorted(&self.header, QUERY_NAME);

        if !(is_tc_sorted || group.allow_unmapped && is_qname_sorted) {
            if group.allow_unmapped {
                bail!(
                    "Input BAM must be template-coordinate sorted or queryname sorted \
                    when --allow-unmapped is enabled.\n\n\
                    To queryname sort your BAM file, run:\n  \
                    samtools sort -n input.bam -o sorted.bam"
                );
            }
            bail!(
                "Input BAM must be template-coordinate sorted (header must advertise \
                SO:unsorted, GO:query, and SS:template-coordinate).\n\n\
                To sort your BAM file, run:\n  \
                fgumi sort -i input.bam -o sorted.bam --order template-coordinate"
            );
        }

        if is_tc_sorted {
            info!("Input is template-coordinate sorted");
        } else {
            info!("Input is queryname sorted (accepted with --allow-unmapped)");
            info!("All unmapped reads will form a single position group per library/cell");
        }

        // Tag constants per SAM specification.
        let raw_tag: [u8; 2] = *SamTag::RX;
        let assign_tag_bytes: [u8; 2] = *SamTag::MI;

        // Build filter configuration from GroupOptions fields (mirrors
        // GroupReadsByUmi::execute's construction of filter_config).
        let filter_config = GroupFilterConfig {
            umi_tag: raw_tag,
            min_mapq: group.min_map_q,
            include_non_pf: group.include_non_pf_reads,
            min_umi_length: group.min_umi_length,
            no_umi: group.no_umi,
            allow_unmapped: group.allow_unmapped,
        };

        // Per-thread metric accumulators.
        let accumulators =
            PerThreadAccumulator::<crate::commands::group::GroupMetricsAccumulator>::new(
                num_threads,
            );
        let accumulators_for_process = Arc::clone(&accumulators);

        // ── Step factories (see chains::commands::group) ──────────────────
        let process_step = build_group_process_step(
            self.tuning.per_step_byte_limit,
            group.effective_strategy,
            group.effective_edits,
            group.index_threshold,
            raw_tag,
            group.no_umi,
            group.allow_unmapped,
            num_threads,
            filter_config,
            accumulators_for_process,
        );
        let mi_assign_step = build_group_mi_assign_step(self.tuning.per_step_byte_limit);

        // Wire GroupByPosition → ProcessOrdered → MiAssign (always appended).
        let tail =
            self.pipeline.append_step(GroupByPosition::new(self.tuning.per_step_byte_limit), tail);
        let tail = self.pipeline.append_step(process_step, tail);
        let tail = self.pipeline.append_step(mi_assign_step, tail);

        if position == StagePosition::Terminal {
            // Terminal: append SerializeGroups so the chain tail is DecompressedBlock,
            // ready for BgzfCompress → WriteBgzfFile.
            let serialize_step = build_group_serialize_step(
                self.tuning.per_step_byte_limit,
                assign_tag_bytes,
                Arc::clone(&self.progress_records),
            );
            let tail = self.pipeline.append_step(serialize_step, tail);
            self.current_tail = Some(tail);
            // chain_tail_kind remains DecodedRecordBatch (actually DecompressedBlock,
            // but that's used only for add_sink detection, not for stage routing).
        } else {
            // Intermediate: leave tail as BatchedProcessedPositionGroups so the
            // next stage (add_simplex / add_duplex / add_codec) can wire
            // TemplatesToMiGroups instead of GroupByMi.
            self.current_tail = Some(tail);
            self.chain_tail_kind = ChainTailKind::BatchedProcessedPositionGroups;
        }

        // Store assign_tag_bytes in the group options for downstream TemplatesToMiGroups
        // bridge. We persist it via a field instead of re-deriving from SamTag::MI so
        // the consensus stages don't need to import SamTag themselves.
        // (assign_tag_bytes == *SamTag::MI, which is a pub const, so deriving it
        // downstream is equally fine — but we set it as a field for clarity.)
        //
        // NOTE: we pass assign_tag_bytes into a new `group_assign_tag` field.
        // That field is Option<[u8;2]> so it's None before add_group runs and
        // Some after.  The consensus stages check chain_tail_kind first; if it's
        // BatchedProcessedPositionGroups they read assign_tag_bytes from the
        // SamTag::MI constant directly (it's the same value, always).

        // Register the group finalize hook.
        // The chain-level StageTimingFinalizeHook is inserted at index 0 by
        // build() after all stages have been added — that way a single timer
        // covers the full pipeline run regardless of how many stages there are.
        self.finalize.push(Box::new(GroupFinalizeHook {
            accumulators,
            output_path,
            family_size_histogram: group.family_size_histogram.clone(),
            grouping_metrics: group.grouping_metrics.clone(),
            metrics_prefix: group.metrics_prefix.clone(),
            timer,
        }));

        Ok(())
    }

    /// Simplex-specific step sequence:
    ///
    /// `GroupByMi` →
    /// `ProcessWithWorkerState` (consensus calling + rejects streaming + metrics).
    ///
    /// For [`StagePosition::Terminal`], the chain tail is [`DecompressedBlock`]
    /// (bytes ready for `BgzfCompress`).
    ///
    /// For [`StagePosition::Intermediate`] (e.g. fused `consensus → filter`),
    /// [`Self::finish_consensus_tail`] appends a [`DecodeRecords`] step so the
    /// chain tail is [`DecodedRecordBatch`] for the next stage — consensus
    /// output is already record-aligned, so no boundary/parse step is needed.
    ///
    /// Applies the M5 fix: enforces that `--ref` requires `--methylation-mode`
    /// to be set (mirrors the same check in `Simplex::execute`'s
    /// single-threaded fast path).
    ///
    /// Registers a `SimplexFinalizeHook` for accumulators reduce, stats write,
    /// overlapping-consensus stats log, summary banner, rejects-writer
    /// finalize, and timer.
    /// The chain-level [`StageTimingFinalizeHook`] is registered by
    /// [`Self::build`] (not here), so it captures the correct `Instant`.
    ///
    /// # Errors
    ///
    /// Returns errors if simplex options are missing from the spec bag, or if
    /// `--ref` is set without `--methylation-mode`.
    ///
    /// [`BatchedMiGroups`]: crate::pipeline::steps::group::mi::BatchedMiGroups
    /// [`DecompressedBlock`]: crate::pipeline::steps::types::DecompressedBlock
    #[allow(clippy::too_many_lines)]
    fn add_simplex(&mut self, position: StagePosition) -> Result<()> {
        use crate::commands::common::{
            MethylationRef, load_methylation_reference, resolve_methylation_mode,
            warn_unwired_pipeline_flags,
        };
        use crate::commands::consensus_runner::create_unmapped_consensus_header;
        use crate::logging::OperationTimer;
        use crate::per_thread_accumulator::PerThreadAccumulator;
        use crate::pipeline::chains::commands::simplex::{
            CollectedSimplexMetrics, SimplexConsensusCaptures, SimplexFinalizeHook,
            build_simplex_consensus_step_kept_only, build_simplex_consensus_step_with_rejects,
        };
        use crate::pipeline::steps::group::mi::GroupByMi;
        use crate::sam::SamTag;
        use crate::vanilla_consensus_caller::VanillaUmiConsensusOptions;
        use log::info;
        use noodles::sam::alignment::record::data::field::Tag;

        let simplex = self
            .spec
            .stage_opts
            .simplex
            .as_ref()
            .ok_or_else(|| anyhow!("Stage::Simplex options missing from StageOptionsBag"))?;

        // M5 fix: mirrors the same check in Simplex::execute's single-threaded
        // fast path. Enforced here so runall's execute_consensus_only path also
        // rejects --ref without --methylation-mode.
        if simplex.reference.is_some() && simplex.methylation_mode.is_none() {
            bail!("--ref requires --methylation-mode to be set");
        }

        let tail = self.current_tail.expect("add_simplex called before add_source");

        // Resolve source path for log messages only.
        let input_path = self.resolve_log_input_path();
        let output_path = match &self.spec.sink {
            SinkSpec::Bam(p) | SinkSpec::BamWithIndex(p) => p.clone(),
        };

        let timer = OperationTimer::new("Calling simplex consensus");

        info!("Starting Simplex");
        info!("Input: {}", input_path.display());
        info!("Output: {}", output_path.display());
        info!("Min reads: {}", simplex.min_reads);
        if let Some(max) = simplex.max_reads {
            info!("Max reads: {max}");
        }
        // Reconstruct the shared option structs from the inlined SimplexOptions
        // flat fields.
        let consensus = simplex.consensus();
        let overlapping = simplex.overlapping();
        info!("Error rate pre-UMI: Q{}", consensus.error_rate_pre_umi);
        info!("Error rate post-UMI: Q{}", consensus.error_rate_post_umi);

        let overlapping_enabled = overlapping.is_enabled();
        if overlapping_enabled {
            info!("Overlapping consensus calling enabled");
        }

        warn_unwired_pipeline_flags(&self.spec.scheduler, &self.spec.queue_memory);
        let num_threads = self.spec.threading.num_threads();
        info!("{}", self.spec.threading.log_message());
        info!("Using pipeline with {num_threads} threads");

        info!("Processing reads and calling consensus (streaming)...");

        // Resolve methylation mode and load reference if needed.
        let methylation_mode = resolve_methylation_mode(simplex.methylation_mode);
        let methylation_ref: MethylationRef =
            load_methylation_reference(methylation_mode, &simplex.reference, &self.header)?;

        // create_unmapped_consensus_header already inserts the @PG record via
        // add_pg_to_builder, so no second add_pg_record call is needed here.
        // We pass self.header (the input header with @PG) for reference
        // sequence lines, RG, etc.
        let output_header = create_unmapped_consensus_header(
            &self.header,
            &simplex.read_group.read_group_id,
            "Read group",
            &self.spec.command_line,
        )?;

        // Derive read_name_prefix from the *input* header (before replacing
        // self.header) so the prefix is consistent with the original input
        // read groups — mirrors what the old build_simplex_chain did, where
        // `prefix_or_from_header` was called on `input_header` (not on the
        // freshly constructed output_header). If the input has no read groups,
        // make_prefix_from_header returns "" which is the correct empty prefix.
        let read_name_prefix = simplex.read_group.prefix_or_from_header(&self.header);

        // Capture the input header (raw-input RG/PG + @HD sort fields) BEFORE
        // replacing self.header with the consensus output header. Per the PR
        // #332 rejects contract, the rejects branch is written with this input
        // header verbatim (rejects are raw-input records in input order).
        let input_header = self.header.clone();

        // Replace self.header with the consensus output header so add_sink()
        // uses the correct header for WriteBgzfFile.
        self.replace_header(output_header);

        let track_rejects = simplex.rejects_opts.is_enabled();

        // Per-thread metrics accumulator.
        let accumulators = PerThreadAccumulator::<CollectedSimplexMetrics>::new(num_threads);
        let accumulators_for_step = Arc::clone(&accumulators);

        // Tag constants.
        let cell_tag = Tag::from(SamTag::CB);

        let read_group_id = simplex.read_group.read_group_id.clone();
        let consensus_options = VanillaUmiConsensusOptions {
            tag: "MI".to_string(),
            error_rate_pre_umi: consensus.error_rate_pre_umi,
            error_rate_post_umi: consensus.error_rate_post_umi,
            min_input_base_quality: consensus.min_input_base_quality,
            min_reads: simplex.min_reads,
            max_reads: simplex.max_reads,
            produce_per_base_tags: consensus.output_per_base_tags,
            seed: Some(42),
            trim: consensus.trim,
            min_consensus_base_quality: consensus.min_consensus_base_quality,
            cell_tag: Some(cell_tag),
            methylation_mode,
        };

        let progress_records = Arc::clone(&self.progress_records);

        // ── Step factories (see chains::commands::simplex) ────────────────

        let consensus_cap = SimplexConsensusCaptures {
            track_rejects,
            overlapping_enabled,
            consensus_options,
            read_name_prefix,
            read_group_id,
            methylation_ref,
            accumulators: accumulators_for_step,
            min_reads: simplex.min_reads,
            progress: progress_records,
        };

        // ── Group-MI preamble: two paths depending on the incoming tail type ──
        //
        // Path 1 (normal): tail is DecodedRecordBatch → prepend GroupByMi.
        // Path 2 (fused group→simplex): tail is BatchedProcessedPositionGroups
        //   (from add_group(Intermediate)) → prepend TemplatesToMiGroups bridge.
        //
        // TemplatesToMiGroups replicates build_group_and_fuse_steps's
        // `templates_to_mi_step` from runall.rs.  The assign_tag_bytes is always
        // *SamTag::MI (== [b'M', b'I']) — the same constant used in add_group.
        // duplex_strip_strand = false for simplex (simplex does not use /A /B suffixes).
        let tail = if self.chain_tail_kind == ChainTailKind::BatchedProcessedPositionGroups {
            use crate::mi_group::MiGroup;
            use crate::pipeline::steps::group::mi::BatchedMiGroups;
            use crate::pipeline::steps::group::position::BatchedProcessedPositionGroups;
            use crate::pipeline::steps::process::process_with_worker_state;
            use fgumi_raw_bam::RawRecord;
            use std::io;

            let assign_tag_bytes: [u8; 2] = *SamTag::MI;
            let duplex_strip_strand = false;

            let templates_to_mi_step = process_with_worker_state::<
                BatchedProcessedPositionGroups,
                BatchedMiGroups,
                _,
                FuseState,
                _,
            >(
                "TemplatesToMiGroups",
                self.tuning.per_step_byte_limit,
                || FuseState {
                    scratch: Vec::with_capacity(512),
                    mi_buf: String::with_capacity(16),
                },
                move |state: &mut FuseState,
                      item: BatchedProcessedPositionGroups|
                      -> io::Result<BatchedMiGroups> {
                    let BatchedProcessedPositionGroups { batch_serial, groups } = item;
                    let mut mi_groups: Vec<MiGroup> = Vec::new();
                    for processed in &groups {
                        let mut current_mi: Option<String> = None;
                        let mut current_records: Vec<RawRecord> = Vec::new();
                        for template in &processed.templates {
                            if !template.mi.is_assigned() {
                                continue;
                            }
                            let mi_key = if duplex_strip_strand {
                                template.mi.id().unwrap().to_string()
                            } else {
                                template.mi.write_with_offset(0, &mut state.mi_buf);
                                state.mi_buf.clone()
                            };
                            template.mi.write_with_offset(0, &mut state.mi_buf);
                            match &current_mi {
                                Some(curr) if curr != &mi_key => {
                                    if !current_records.is_empty() {
                                        mi_groups.push(MiGroup::new(
                                            std::mem::take(&mut current_mi).unwrap(),
                                            std::mem::take(&mut current_records),
                                        ));
                                    }
                                    current_mi = Some(mi_key);
                                }
                                None => {
                                    current_mi = Some(mi_key);
                                }
                                Some(_) => {}
                            }
                            for raw in [template.r1(), template.r2()].into_iter().flatten() {
                                state.scratch.clear();
                                state.scratch.extend_from_slice(raw);
                                fgumi_raw_bam::update_string_tag(
                                    &mut state.scratch,
                                    assign_tag_bytes,
                                    state.mi_buf.as_bytes(),
                                );
                                current_records.push(RawRecord::from(state.scratch.clone()));
                            }
                        }
                        if !current_records.is_empty() {
                            if let Some(mi) = current_mi.take() {
                                mi_groups.push(MiGroup::new(mi, current_records));
                            }
                        }
                    }
                    Ok(BatchedMiGroups::new(batch_serial, mi_groups))
                },
            );
            self.pipeline.append_step(templates_to_mi_step, tail)
        } else {
            // Normal path: DecodedRecordBatch → GroupByMi → BatchedMiGroups.
            let group_mi_step = GroupByMi::new(*SamTag::MI, self.tuning.per_step_byte_limit)
                .with_cell_tag(Some(*SamTag::CB));
            self.pipeline.append_step(group_mi_step, tail)
        };

        // Wire the simplex consensus step. With `--rejects` it is a 2-output
        // step (branch 0 = consensus DecompressedBlock, branch 1 = rejects
        // DecompressedBlock → BgzfCompress → WriteBgzfFile with the INPUT
        // header, in input order — the PR #332 fan-out contract); without
        // rejects it is the 1-output kept-only variant (the framework has no
        // public discard sink). Mirrors add_correct.
        let limit = self.tuning.per_step_byte_limit;
        let consensus_branch0 = if track_rejects {
            use crate::pipeline::core::topology::BranchIdx;
            use crate::pipeline::steps::bgzf::compress::BgzfCompress;
            use crate::pipeline::steps::sink::write_bgzf::WriteBgzfFile;

            let rejects_path = simplex.rejects_opts.rejects.as_ref().ok_or_else(|| {
                anyhow!("rejects path unexpectedly None when track_rejects is set")
            })?;
            let rejects_write =
                WriteBgzfFile::new(rejects_path, &input_header, self.tuning.compression_level)
                    .map_err(|e| anyhow!("WriteBgzfFile (simplex rejects): {e}"))?;

            let step = build_simplex_consensus_step_with_rejects(limit, consensus_cap);
            let pt = self.pipeline.append_step(step, tail);

            let rejects_branch = (pt.0, BranchIdx(1));
            let rejects_compress_tail = self.pipeline.append_step(
                BgzfCompress::new(self.tuning.compression_level, limit),
                rejects_branch,
            );
            self.pipeline.append_step(rejects_write, rejects_compress_tail);

            pt
        } else {
            let step = build_simplex_consensus_step_kept_only(limit, consensus_cap);
            self.pipeline.append_step(step, tail)
        };
        // Branch 0 = consensus DecompressedBlock. For an Intermediate consensus
        // stage finish_consensus_tail appends DecodeRecords; Terminal leaves the
        // DecompressedBlock for add_sink.
        let tail = self.finish_consensus_tail(consensus_branch0, position);
        self.current_tail = Some(tail);

        // Register the simplex finalize hook.
        // The chain-level StageTimingFinalizeHook is inserted at index 0 by
        // build() after all stages have been added — that way a single timer
        // covers the full pipeline run regardless of how many stages there are.
        self.finalize.push(Box::new(SimplexFinalizeHook {
            accumulators,
            stats_path: simplex.stats_opts.stats.clone(),
            overlapping_enabled,
            timer,
        }));

        Ok(())
    }

    /// Duplex-specific step sequence:
    ///
    /// `GroupByMi` (with record filter + MI transform) →
    /// `ProcessWithWorkerState` (duplex consensus calling + rejects streaming + metrics).
    ///
    /// For [`StagePosition::Terminal`], the chain tail is [`DecompressedBlock`]
    /// (bytes ready for `BgzfCompress`).
    ///
    /// For [`StagePosition::Intermediate`] (e.g. fused `consensus → filter`),
    /// [`Self::finish_consensus_tail`] appends a [`DecodeRecords`] step so the
    /// chain tail is [`DecodedRecordBatch`] for the next stage.
    ///
    /// The record filter skips secondary/supplementary reads and requires each
    /// record to be mapped or have a mapped mate. The MI transform strips `/A`
    /// and `/B` suffixes so both strands group into the same `MiGroup`.
    ///
    /// Applies the read_name_prefix-from-input-header fix: derives the prefix
    /// from `self.header` (the input header, before it is replaced with the
    /// consensus output header) so the prefix is consistent with the original
    /// input read groups — mirrors the Phase 2 behaviour in `build_duplex_chain`
    /// where `prefix_or_from_header` was called on `input_header`.
    ///
    /// Applies the M5 fix: enforces that `--ref` requires `--methylation-mode`
    /// to be set (mirrors the same check in `Duplex::execute`'s single-threaded
    /// fast path).
    ///
    /// Registers a `DuplexFinalizeHook` for accumulators reduce, stats write,
    /// overlapping-consensus stats log, summary banner, rejects-writer
    /// finalize, and timer.
    /// The chain-level [`StageTimingFinalizeHook`] is registered by
    /// [`Self::build`] (not here), so it captures the correct `Instant`.
    ///
    /// # Errors
    ///
    /// Returns errors if duplex options are missing from the spec bag, or if
    /// `--ref` is set without `--methylation-mode`.
    ///
    /// [`BatchedMiGroups`]: crate::pipeline::steps::group::mi::BatchedMiGroups
    /// [`DecompressedBlock`]: crate::pipeline::steps::types::DecompressedBlock
    #[allow(clippy::too_many_lines)]
    fn add_duplex(&mut self, position: StagePosition) -> Result<()> {
        use crate::commands::common::{
            load_methylation_reference, resolve_methylation_mode, warn_unwired_pipeline_flags,
        };
        use crate::commands::consensus_runner::create_unmapped_consensus_header;
        use crate::logging::OperationTimer;
        use crate::per_thread_accumulator::PerThreadAccumulator;
        use crate::pipeline::chains::commands::duplex::{
            CollectedDuplexMetrics, DuplexConsensusCaptures, DuplexFinalizeHook,
            build_duplex_consensus_step_kept_only, build_duplex_consensus_step_with_rejects,
        };
        use crate::pipeline::steps::group::mi::GroupByMi;
        use crate::sam::SamTag;
        use crate::umi::extract_mi_base;
        use fgumi_raw_bam::RawRecordView;
        use log::info;
        use noodles::sam::alignment::record::data::field::Tag;

        let duplex = self
            .spec
            .stage_opts
            .duplex
            .as_ref()
            .ok_or_else(|| anyhow!("Stage::Duplex options missing from StageOptionsBag"))?;

        // M5 fix: mirrors the same check in Duplex::execute's single-threaded
        // fast path. Enforced here so runall's execute_consensus_only path also
        // rejects --ref without --methylation-mode.
        if duplex.reference.is_some() && duplex.methylation_mode.is_none() {
            bail!("--ref requires --methylation-mode to be set");
        }

        let tail = self.current_tail.expect("add_duplex called before add_source");

        // Resolve source path for log messages only.
        let input_path = self.resolve_log_input_path();
        let output_path = match &self.spec.sink {
            SinkSpec::Bam(p) | SinkSpec::BamWithIndex(p) => p.clone(),
        };

        let timer = OperationTimer::new("Calling duplex consensus");

        // Reconstruct the shared option structs from the inlined DuplexOptions
        // flat fields.
        let consensus = duplex.consensus();
        let overlapping = duplex.overlapping();

        info!("Starting Duplex");
        info!("Input: {}", input_path.display());
        info!("Output: {}", output_path.display());
        info!("Min reads: {:?}", duplex.min_reads);
        info!("Min base quality: {}", consensus.min_input_base_quality);
        info!("Output per-base tags: {}", consensus.output_per_base_tags);
        info!("Trim reads: {}", consensus.trim);
        info!("Max reads per strand: {:?}", duplex.max_reads_per_strand);
        info!("Consensus call overlapping bases: {}", overlapping.consensus_call_overlapping_bases);

        let overlapping_enabled = overlapping.consensus_call_overlapping_bases;
        if overlapping_enabled {
            info!("Overlapping consensus calling enabled");
        }

        warn_unwired_pipeline_flags(&self.spec.scheduler, &self.spec.queue_memory);
        let num_threads = self.spec.threading.num_threads();
        info!("{}", self.spec.threading.log_message());
        info!("Using pipeline with {num_threads} threads");

        info!("Processing reads and calling duplex consensus (streaming)...");
        info!("Reading input");

        // Resolve methylation mode and load reference if needed.
        let methylation_mode = resolve_methylation_mode(duplex.methylation_mode);
        let methylation_ref =
            load_methylation_reference(methylation_mode, &duplex.reference, &self.header)?;

        // create_unmapped_consensus_header already inserts the @PG record via
        // add_pg_to_builder, so no second add_pg_record call is needed here.
        // We pass self.header (the input header with @PG) for reference
        // sequence lines, RG, etc.
        let output_header = create_unmapped_consensus_header(
            &self.header,
            &duplex.read_group.read_group_id,
            "Read group",
            &self.spec.command_line,
        )?;

        // Derive read_name_prefix from the *input* header (before replacing
        // self.header) so the prefix is consistent with the original input
        // read groups — mirrors what the old build_duplex_chain did, where
        // `prefix_or_from_header` was called on `input_header` (not on the
        // freshly constructed output_header). If the input has no read groups,
        // make_prefix_from_header returns "" which is the correct empty prefix.
        let read_name_prefix = duplex.read_group.prefix_or_from_header(&self.header);

        // Capture the input header BEFORE replacing self.header with the
        // consensus output header — the rejects branch is written with the
        // input header verbatim (PR #332 contract).
        let input_header = self.header.clone();

        // Replace self.header with the consensus output header so add_sink()
        // uses the correct header for WriteBgzfFile.
        self.replace_header(output_header);

        let track_rejects = duplex.rejects_opts.is_enabled();

        // Per-thread metrics accumulator.
        let accumulators = PerThreadAccumulator::<CollectedDuplexMetrics>::new(num_threads);
        let accumulators_for_step = Arc::clone(&accumulators);

        // Tag constants.
        let cell_tag = Tag::from(SamTag::CB);

        let read_group_id = duplex.read_group.read_group_id.clone();
        let min_reads = duplex.min_reads.clone();
        let min_input_base_quality = consensus.min_input_base_quality;
        let output_per_base_tags = consensus.output_per_base_tags;
        let trim = consensus.trim;
        let max_reads_per_strand = duplex.max_reads_per_strand;
        let error_rate_pre_umi = consensus.error_rate_pre_umi;
        let error_rate_post_umi = consensus.error_rate_post_umi;

        let progress_records = Arc::clone(&self.progress_records);

        // ── Step factories (see chains::commands::duplex) ────────────────────

        let consensus_cap = DuplexConsensusCaptures {
            track_rejects,
            overlapping_enabled,
            methylation_ref,
            methylation_mode,
            read_name_prefix,
            read_group_id,
            min_reads,
            min_input_base_quality,
            output_per_base_tags,
            trim,
            max_reads_per_strand,
            error_rate_pre_umi,
            error_rate_post_umi,
            cell_tag,
            accumulators: accumulators_for_step,
            progress: progress_records,
        };

        // ── Group-MI preamble: two paths depending on the incoming tail type ──
        //
        // Path 1 (normal): tail is DecodedRecordBatch → prepend GroupByMi with
        //   record filter (skip secondary/supplementary) + MI transform (strip /A /B).
        // Path 2 (fused group→duplex): tail is BatchedProcessedPositionGroups
        //   (from add_group(Intermediate)) → prepend TemplatesToMiGroups bridge
        //   with duplex_strip_strand = true.
        let tail = if self.chain_tail_kind == ChainTailKind::BatchedProcessedPositionGroups {
            use crate::mi_group::MiGroup;
            use crate::pipeline::steps::group::mi::BatchedMiGroups;
            use crate::pipeline::steps::group::position::BatchedProcessedPositionGroups;
            use crate::pipeline::steps::process::process_with_worker_state;
            use fgumi_raw_bam::RawRecord;
            use std::io;

            let assign_tag_bytes: [u8; 2] = *SamTag::MI;
            let duplex_strip_strand = true; // duplex: strip /A /B so both strands group together

            let templates_to_mi_step = process_with_worker_state::<
                BatchedProcessedPositionGroups,
                BatchedMiGroups,
                _,
                FuseState,
                _,
            >(
                "TemplatesToMiGroups",
                self.tuning.per_step_byte_limit,
                || FuseState {
                    scratch: Vec::with_capacity(512),
                    mi_buf: String::with_capacity(16),
                },
                move |state: &mut FuseState,
                      item: BatchedProcessedPositionGroups|
                      -> io::Result<BatchedMiGroups> {
                    let BatchedProcessedPositionGroups { batch_serial, groups } = item;
                    let mut mi_groups: Vec<MiGroup> = Vec::new();
                    for processed in &groups {
                        let mut current_mi: Option<String> = None;
                        let mut current_records: Vec<RawRecord> = Vec::new();
                        for template in &processed.templates {
                            if !template.mi.is_assigned() {
                                continue;
                            }
                            let mi_key = if duplex_strip_strand {
                                template.mi.id().unwrap().to_string()
                            } else {
                                template.mi.write_with_offset(0, &mut state.mi_buf);
                                state.mi_buf.clone()
                            };
                            template.mi.write_with_offset(0, &mut state.mi_buf);
                            match &current_mi {
                                Some(curr) if curr != &mi_key => {
                                    if !current_records.is_empty() {
                                        mi_groups.push(MiGroup::new(
                                            std::mem::take(&mut current_mi).unwrap(),
                                            std::mem::take(&mut current_records),
                                        ));
                                    }
                                    current_mi = Some(mi_key);
                                }
                                None => {
                                    current_mi = Some(mi_key);
                                }
                                Some(_) => {}
                            }
                            for raw in [template.r1(), template.r2()].into_iter().flatten() {
                                state.scratch.clear();
                                state.scratch.extend_from_slice(raw);
                                fgumi_raw_bam::update_string_tag(
                                    &mut state.scratch,
                                    assign_tag_bytes,
                                    state.mi_buf.as_bytes(),
                                );
                                current_records.push(RawRecord::from(state.scratch.clone()));
                            }
                        }
                        if !current_records.is_empty() {
                            if let Some(mi) = current_mi.take() {
                                mi_groups.push(MiGroup::new(mi, current_records));
                            }
                        }
                    }
                    Ok(BatchedMiGroups::new(batch_serial, mi_groups))
                },
            );
            self.pipeline.append_step(templates_to_mi_step, tail)
        } else {
            // Normal path: DecodedRecordBatch → GroupByMi → BatchedMiGroups.
            // Record filter: skip secondary/supplementary; require mapped or mapped-mate.
            let record_filter = |raw: &[u8]| -> bool {
                let flg = RawRecordView::new(raw).flags();
                if flg & fgumi_raw_bam::flags::SECONDARY != 0
                    || flg & fgumi_raw_bam::flags::SUPPLEMENTARY != 0
                {
                    return false;
                }
                let is_mapped = flg & fgumi_raw_bam::flags::UNMAPPED == 0;
                let has_mapped_mate = flg & fgumi_raw_bam::flags::PAIRED != 0
                    && flg & fgumi_raw_bam::flags::MATE_UNMAPPED == 0;
                is_mapped || has_mapped_mate
            };
            // MI transform: strip /A and /B so both strands group together.
            let mi_transform = |mi_bytes: &[u8]| -> String {
                let mi_str = std::borrow::Cow::from(std::str::from_utf8(mi_bytes).unwrap_or(""));
                extract_mi_base(&mi_str).to_string()
            };
            let group_mi_step = GroupByMi::new(*SamTag::MI, self.tuning.per_step_byte_limit)
                .with_cell_tag(Some(*SamTag::CB))
                .with_record_filter(record_filter)
                .with_mi_transform(mi_transform);
            self.pipeline.append_step(group_mi_step, tail)
        };

        // Wire the duplex consensus step. With `--rejects` it is a 2-output step
        // (branch 0 = consensus, branch 1 = rejects → BgzfCompress → WriteBgzfFile
        // with the INPUT header, in input order — the PR #332 fan-out contract);
        // without rejects it is the 1-output kept-only variant. Mirrors add_correct.
        let limit = self.tuning.per_step_byte_limit;
        let consensus_branch0 = if track_rejects {
            use crate::pipeline::core::topology::BranchIdx;
            use crate::pipeline::steps::bgzf::compress::BgzfCompress;
            use crate::pipeline::steps::sink::write_bgzf::WriteBgzfFile;

            let rejects_path = duplex.rejects_opts.rejects.as_ref().ok_or_else(|| {
                anyhow!("rejects path unexpectedly None when track_rejects is set")
            })?;
            let rejects_write =
                WriteBgzfFile::new(rejects_path, &input_header, self.tuning.compression_level)
                    .map_err(|e| anyhow!("WriteBgzfFile (duplex rejects): {e}"))?;

            let step = build_duplex_consensus_step_with_rejects(limit, consensus_cap);
            let pt = self.pipeline.append_step(step, tail);

            let rejects_branch = (pt.0, BranchIdx(1));
            let rejects_compress_tail = self.pipeline.append_step(
                BgzfCompress::new(self.tuning.compression_level, limit),
                rejects_branch,
            );
            self.pipeline.append_step(rejects_write, rejects_compress_tail);

            pt
        } else {
            let step = build_duplex_consensus_step_kept_only(limit, consensus_cap);
            self.pipeline.append_step(step, tail)
        };
        // Branch 0 = consensus DecompressedBlock. Intermediate appends
        // DecodeRecords; Terminal leaves the DecompressedBlock for add_sink.
        let tail = self.finish_consensus_tail(consensus_branch0, position);
        self.current_tail = Some(tail);

        // Register the duplex finalize hook.
        // The chain-level StageTimingFinalizeHook is inserted at index 0 by
        // build() after all stages have been added — that way a single timer
        // covers the full pipeline run regardless of how many stages there are.
        self.finalize.push(Box::new(DuplexFinalizeHook {
            accumulators,
            stats_path: duplex.stats_opts.stats.clone(),
            overlapping_enabled,
            timer,
        }));

        Ok(())
    }

    /// Codec-specific step sequence:
    ///
    /// `GroupByMi` →
    /// `ProcessWithWorkerState` (codec consensus calling + rejects streaming + metrics).
    ///
    /// For [`StagePosition::Terminal`], the chain tail is [`DecompressedBlock`]
    /// (bytes ready for `BgzfCompress`).
    ///
    /// For [`StagePosition::Intermediate`] (e.g. fused `consensus → filter`),
    /// [`Self::finish_consensus_tail`] appends a [`DecodeRecords`] step so the
    /// chain tail is [`DecodedRecordBatch`] for the next stage.
    ///
    /// Unlike duplex, codec does **not** apply a record filter, MI-tag transform,
    /// overlapping consensus, or methylation mode — matching fgbio's
    /// `CallCodecConsensusReads` behaviour.
    ///
    /// Applies the read_name_prefix-from-input-header fix: derives the prefix
    /// from `self.header` (the input header, before it is replaced with the
    /// consensus output header) so the prefix is consistent with the original
    /// input read groups — mirrors the Phase 2 behaviour in `build_codec_chain`
    /// where `prefix_or_from_header` was called on `input_header`.
    ///
    /// Registers a `CodecFinalizeHook` for accumulators reduce, stats write,
    /// summary banner, rejects-writer finalize, and timer.
    /// The chain-level [`StageTimingFinalizeHook`] is registered by
    /// [`Self::build`] (not here), so it captures the correct `Instant`.
    ///
    /// # Errors
    ///
    /// Returns errors if codec options are missing from the spec bag.
    ///
    /// [`BatchedMiGroups`]: crate::pipeline::steps::group::mi::BatchedMiGroups
    /// [`DecompressedBlock`]: crate::pipeline::steps::types::DecompressedBlock
    #[allow(clippy::too_many_lines)]
    fn add_codec(&mut self, position: StagePosition) -> Result<()> {
        use crate::commands::common::warn_unwired_pipeline_flags;
        use crate::commands::consensus_runner::create_unmapped_consensus_header;
        use crate::consensus::codec_caller::CodecConsensusOptions;
        use crate::logging::OperationTimer;
        use crate::per_thread_accumulator::PerThreadAccumulator;
        use crate::pipeline::chains::commands::codec::{
            CodecConsensusCaptures, CodecFinalizeHook, CollectedCodecMetrics,
            build_codec_consensus_step_kept_only, build_codec_consensus_step_with_rejects,
        };
        use crate::pipeline::steps::group::mi::GroupByMi;
        use crate::sam::SamTag;
        use log::info;
        use noodles::sam::alignment::record::data::field::Tag;

        let codec = self
            .spec
            .stage_opts
            .codec
            .as_ref()
            .ok_or_else(|| anyhow!("Stage::Codec options missing from StageOptionsBag"))?;

        let tail = self.current_tail.expect("add_codec called before add_source");

        // Resolve source path for log messages only.
        let input_path = self.resolve_log_input_path();
        let output_path = match &self.spec.sink {
            SinkSpec::Bam(p) | SinkSpec::BamWithIndex(p) => p.clone(),
        };

        let timer = OperationTimer::new("Calling CODEC consensus");

        // Reconstruct the shared ConsensusCallingOptions from the inlined
        // CodecOptions flat fields.
        let consensus = codec.consensus();

        info!("Starting CODEC consensus calling");
        info!("Input: {}", input_path.display());
        info!("Output: {}", output_path.display());
        info!("Min reads: {}", codec.min_reads);
        if let Some(max) = codec.max_reads {
            info!("Max reads: {max}");
        }
        info!("Error rate pre-UMI: Q{}", consensus.error_rate_pre_umi);
        info!("Error rate post-UMI: Q{}", consensus.error_rate_post_umi);
        info!("Min duplex length: {}", codec.min_duplex_length);

        warn_unwired_pipeline_flags(&self.spec.scheduler, &self.spec.queue_memory);
        let num_threads = self.spec.threading.num_threads();
        info!("Worker threads: {num_threads}");
        info!("Reader threads: {num_threads}");
        if consensus.trim {
            info!("Quality trimming enabled");
        }

        info!("Processing reads and calling consensus (streaming)...");
        info!("Reading input");

        // create_unmapped_consensus_header already inserts the @PG record via
        // add_pg_to_builder, so no second add_pg_record call is needed here.
        // We pass self.header (the input header with @PG) for reference
        // sequence lines, RG, etc.
        let output_header = create_unmapped_consensus_header(
            &self.header,
            &codec.read_group.read_group_id,
            "Read group",
            &self.spec.command_line,
        )?;

        // Derive read_name_prefix from the *input* header (before replacing
        // self.header) so the prefix is consistent with the original input
        // read groups — mirrors what the old build_codec_chain did, where
        // `prefix_or_from_header` was called on `input_header` (not on the
        // freshly constructed output_header). If the input has no read groups,
        // make_prefix_from_header returns "" which is the correct empty prefix.
        let read_name_prefix = codec.read_group.prefix_or_from_header(&self.header);

        // Capture the input header BEFORE replacing self.header with the
        // consensus output header — the rejects branch is written with the
        // input header verbatim (PR #332 contract).
        let input_header = self.header.clone();

        // Replace self.header with the consensus output header so add_sink()
        // uses the correct header for WriteBgzfFile.
        self.replace_header(output_header);

        let track_rejects = codec.rejects_opts.is_enabled();

        // Per-thread metrics accumulator.
        let accumulators = PerThreadAccumulator::<CollectedCodecMetrics>::new(num_threads);
        let accumulators_for_step = Arc::clone(&accumulators);

        // Tag constants.
        let cell_tag = Tag::from(SamTag::CB);

        let read_group_id = codec.read_group.read_group_id.clone();
        let consensus_options = CodecConsensusOptions {
            min_input_base_quality: consensus.min_input_base_quality,
            error_rate_pre_umi: consensus.error_rate_pre_umi,
            error_rate_post_umi: consensus.error_rate_post_umi,
            min_reads_per_strand: codec.min_reads,
            max_reads_per_strand: codec.max_reads,
            min_duplex_length: codec.min_duplex_length,
            single_strand_qual: codec.single_strand_qual,
            outer_bases_qual: codec.outer_bases_qual,
            outer_bases_length: codec.outer_bases_length,
            max_duplex_disagreements: codec.max_duplex_disagreements.unwrap_or(usize::MAX),
            max_duplex_disagreement_rate: codec.max_duplex_disagreement_rate,
            cell_tag: Some(cell_tag),
            produce_per_base_tags: consensus.output_per_base_tags,
            trim: consensus.trim,
            min_consensus_base_quality: consensus.min_consensus_base_quality,
        };

        let progress_records = Arc::clone(&self.progress_records);

        // ── Step factories (see chains::commands::codec) ─────────────────────

        let consensus_cap = CodecConsensusCaptures {
            track_rejects,
            read_name_prefix,
            read_group_id,
            consensus_options,
            accumulators: accumulators_for_step,
            progress: progress_records,
        };

        // ── Group-MI preamble: two paths depending on the incoming tail type ──
        //
        // Path 1 (normal): tail is DecodedRecordBatch → prepend GroupByMi.
        //   Codec: no record filter, no MI transform.
        // Path 2 (fused group→codec): tail is BatchedProcessedPositionGroups
        //   (from add_group(Intermediate)) → prepend TemplatesToMiGroups bridge.
        //   duplex_strip_strand = false (codec does not use /A /B suffixes).
        let tail = if self.chain_tail_kind == ChainTailKind::BatchedProcessedPositionGroups {
            use crate::mi_group::MiGroup;
            use crate::pipeline::steps::group::mi::BatchedMiGroups;
            use crate::pipeline::steps::group::position::BatchedProcessedPositionGroups;
            use crate::pipeline::steps::process::process_with_worker_state;
            use fgumi_raw_bam::RawRecord;
            use std::io;

            let assign_tag_bytes: [u8; 2] = *SamTag::MI;
            let duplex_strip_strand = false;

            let templates_to_mi_step = process_with_worker_state::<
                BatchedProcessedPositionGroups,
                BatchedMiGroups,
                _,
                FuseState,
                _,
            >(
                "TemplatesToMiGroups",
                self.tuning.per_step_byte_limit,
                || FuseState {
                    scratch: Vec::with_capacity(512),
                    mi_buf: String::with_capacity(16),
                },
                move |state: &mut FuseState,
                      item: BatchedProcessedPositionGroups|
                      -> io::Result<BatchedMiGroups> {
                    let BatchedProcessedPositionGroups { batch_serial, groups } = item;
                    let mut mi_groups: Vec<MiGroup> = Vec::new();
                    for processed in &groups {
                        let mut current_mi: Option<String> = None;
                        let mut current_records: Vec<RawRecord> = Vec::new();
                        for template in &processed.templates {
                            if !template.mi.is_assigned() {
                                continue;
                            }
                            let mi_key = if duplex_strip_strand {
                                template.mi.id().unwrap().to_string()
                            } else {
                                template.mi.write_with_offset(0, &mut state.mi_buf);
                                state.mi_buf.clone()
                            };
                            template.mi.write_with_offset(0, &mut state.mi_buf);
                            match &current_mi {
                                Some(curr) if curr != &mi_key => {
                                    if !current_records.is_empty() {
                                        mi_groups.push(MiGroup::new(
                                            std::mem::take(&mut current_mi).unwrap(),
                                            std::mem::take(&mut current_records),
                                        ));
                                    }
                                    current_mi = Some(mi_key);
                                }
                                None => {
                                    current_mi = Some(mi_key);
                                }
                                Some(_) => {}
                            }
                            for raw in [template.r1(), template.r2()].into_iter().flatten() {
                                state.scratch.clear();
                                state.scratch.extend_from_slice(raw);
                                fgumi_raw_bam::update_string_tag(
                                    &mut state.scratch,
                                    assign_tag_bytes,
                                    state.mi_buf.as_bytes(),
                                );
                                current_records.push(RawRecord::from(state.scratch.clone()));
                            }
                        }
                        if !current_records.is_empty() {
                            if let Some(mi) = current_mi.take() {
                                mi_groups.push(MiGroup::new(mi, current_records));
                            }
                        }
                    }
                    Ok(BatchedMiGroups::new(batch_serial, mi_groups))
                },
            );
            self.pipeline.append_step(templates_to_mi_step, tail)
        } else {
            // Normal path: DecodedRecordBatch → GroupByMi → BatchedMiGroups.
            // Codec: no record filter (processes all mapped reads), no MI transform.
            let group_mi_step = GroupByMi::new(*SamTag::MI, self.tuning.per_step_byte_limit)
                .with_cell_tag(Some(*SamTag::CB));
            self.pipeline.append_step(group_mi_step, tail)
        };

        // Wire the codec consensus step. With `--rejects` it is a 2-output step
        // (branch 0 = consensus, branch 1 = rejects → BgzfCompress → WriteBgzfFile
        // with the INPUT header, in input order — the PR #332 fan-out contract);
        // without rejects it is the 1-output kept-only variant. Mirrors add_correct.
        let limit = self.tuning.per_step_byte_limit;
        let consensus_branch0 = if track_rejects {
            use crate::pipeline::core::topology::BranchIdx;
            use crate::pipeline::steps::bgzf::compress::BgzfCompress;
            use crate::pipeline::steps::sink::write_bgzf::WriteBgzfFile;

            let rejects_path = codec.rejects_opts.rejects.as_ref().ok_or_else(|| {
                anyhow!("rejects path unexpectedly None when track_rejects is set")
            })?;
            let rejects_write =
                WriteBgzfFile::new(rejects_path, &input_header, self.tuning.compression_level)
                    .map_err(|e| anyhow!("WriteBgzfFile (codec rejects): {e}"))?;

            let step = build_codec_consensus_step_with_rejects(limit, consensus_cap);
            let pt = self.pipeline.append_step(step, tail);

            let rejects_branch = (pt.0, BranchIdx(1));
            let rejects_compress_tail = self.pipeline.append_step(
                BgzfCompress::new(self.tuning.compression_level, limit),
                rejects_branch,
            );
            self.pipeline.append_step(rejects_write, rejects_compress_tail);

            pt
        } else {
            let step = build_codec_consensus_step_kept_only(limit, consensus_cap);
            self.pipeline.append_step(step, tail)
        };
        // Branch 0 = consensus DecompressedBlock. Intermediate appends
        // DecodeRecords; Terminal leaves the DecompressedBlock for add_sink.
        let tail = self.finish_consensus_tail(consensus_branch0, position);
        self.current_tail = Some(tail);

        // Register the codec finalize hook.
        // The chain-level StageTimingFinalizeHook is inserted at index 0 by
        // build() after all stages have been added — that way a single timer
        // covers the full pipeline run regardless of how many stages there are.
        self.finalize.push(Box::new(CodecFinalizeHook {
            accumulators,
            stats_path: codec.stats_opts.stats.clone(),
            timer,
        }));

        Ok(())
    }

    /// Clip-specific step sequence:
    ///
    /// `GroupBam` →
    /// `ProcessOrdered` (clip templates + metrics accumulation) →
    /// `SerializeBamRecords` — **Terminal only**.
    ///
    /// For [`StagePosition::Terminal`], all three steps are appended and the
    /// chain tail is [`DecompressedBlock`] (bytes ready for `BgzfCompress`).
    ///
    /// For [`StagePosition::Intermediate`], only the first two steps are
    /// appended; the chain tail stays as [`BamTemplateBatch`] for the next
    /// stage to consume. Intermediate clip is not yet needed by any Phase 3
    /// runall combination; calling `add_clip` with `Intermediate` returns
    /// `Err` as a guard until a concrete fused path requires it.
    ///
    /// Applies `clip.update_header_sort_order` to `self.header` so that the
    /// output BAM written by `add_sink` carries the updated sort-order
    /// annotation (if `--sort-order` was requested).
    ///
    /// Registers a `ClipFinalizeHook` for metrics logging + timer.
    /// The chain-level [`StageTimingFinalizeHook`] is registered by
    /// [`Self::build`] (not here), so it captures the correct `Instant`.
    ///
    /// # Errors
    ///
    /// Returns errors if clip options are missing from the spec bag, if the
    /// reference FASTA is absent, if no clipping operation is requested, if
    /// `--metrics` is used with `--threads` mode, or if `position` is
    /// `Intermediate` (not yet implemented).
    #[allow(clippy::too_many_lines)]
    fn add_clip(&mut self, position: StagePosition) -> Result<()> {
        use crate::commands::common::warn_unwired_pipeline_flags;
        use crate::logging::OperationTimer;
        use crate::pipeline::chains::commands::clip::{
            ClipAtomicMetrics, ClipFinalizeHook, ClipProcessCaptures, build_clip_process_step,
            build_clip_serialize_step,
        };
        use crate::pipeline::steps::group::bam::GroupBam;
        use log::info;

        if position == StagePosition::Intermediate {
            bail!(
                "intermediate clip not implemented; clip is typically terminal \
                (no Phase 3 runall combination requires it as an intermediate stage)"
            );
        }

        let clip = self
            .spec
            .stage_opts
            .clip
            .as_ref()
            .ok_or_else(|| anyhow!("Stage::Clip options missing from StageOptionsBag"))?;

        let tail = self.current_tail.expect("add_clip called before add_source");

        // Resolve source path for log messages only.
        let input_path = self.resolve_log_input_path();
        let output_path = match &self.spec.sink {
            SinkSpec::Bam(p) | SinkSpec::BamWithIndex(p) => p.clone(),
        };

        // Reference is always required for clip.
        crate::validation::validate_file_exists(&clip.reference, "Reference FASTA")?;

        // Validate that at least one clipping operation is requested.
        if !clip.upgrade_clipping
            && !clip.clip_overlapping_reads
            && !clip.clip_extending_past_mate
            && clip.read_one_five_prime == 0
            && clip.read_one_three_prime == 0
            && clip.read_two_five_prime == 0
            && clip.read_two_three_prime == 0
        {
            bail!("At least one clipping option is required");
        }

        // --metrics is not produced in --threads mode; fail fast rather than
        // silently dropping a user-requested output file.
        let num_threads = self.spec.threading.num_threads();
        if let Some(path) = &clip.metrics {
            anyhow::bail!(
                "--metrics {} cannot be used with --threads: detailed clipping metrics \
                 are only produced by the single-threaded path",
                path.display()
            );
        }

        let timer = OperationTimer::new("Clipping reads");

        info!("Clip");
        info!("  Input: {}", input_path.display());
        info!("  Output: {}", output_path.display());
        info!("  Clipping mode: {}", clip.clipping_mode);
        info!("  Clip overlapping reads: {}", clip.clip_overlapping_reads);
        info!("  Clip extending past mate: {}", clip.clip_extending_past_mate);
        info!("  {}", self.spec.threading.log_message());

        warn_unwired_pipeline_flags(&self.spec.scheduler, &self.spec.queue_memory);
        info!("{}", self.spec.threading.log_message());
        info!("Using pipeline with {num_threads} threads");

        // Apply clip's sort-order annotation to the output header. This must
        // happen before add_sink() consumes self.header for WriteBgzfFile.
        // Order relative to the @PG record (already in self.header from new())
        // does not matter — @HD and @PG are independent header sections.
        self.replace_header(clip.update_header_sort_order(self.header.clone())?);

        // Load reference (always required for clip).
        let reference = Arc::new(crate::reference::ReferenceReader::new(&clip.reference)?);

        // Shared atomic metric counters.
        let metrics = Arc::new(ClipAtomicMetrics::default());
        let progress_counter = Arc::new(std::sync::atomic::AtomicU64::new(0));

        // ── Step factories (see chains::commands::clip) ──────────────────
        let process_step = build_clip_process_step(
            self.tuning.per_step_byte_limit,
            ClipProcessCaptures {
                clipping_mode: clip.clipping_mode,
                auto_clip_attributes: clip.auto_clip_attributes,
                upgrade_clipping: clip.upgrade_clipping,
                clip_overlapping_reads: clip.clip_overlapping_reads,
                clip_extending_past_mate: clip.clip_extending_past_mate,
                read_one_five_prime: clip.read_one_five_prime,
                read_one_three_prime: clip.read_one_three_prime,
                read_two_five_prime: clip.read_two_five_prime,
                read_two_three_prime: clip.read_two_three_prime,
                header: self.header.clone(),
                reference: Arc::clone(&reference),
                metrics: Arc::clone(&metrics),
                progress: Arc::clone(&progress_counter),
            },
        );
        let serialize_step = build_clip_serialize_step(self.tuning.per_step_byte_limit);

        // Wire the clip-specific steps.
        let tail = self.pipeline.append_step(
            GroupBam::new(self.tuning.template_batch_size, self.tuning.per_step_byte_limit),
            tail,
        );
        let tail = self.pipeline.append_step(process_step, tail);
        let tail = self.pipeline.append_step(serialize_step, tail);
        self.current_tail = Some(tail);

        // Register the clip finalize hook.
        // The chain-level StageTimingFinalizeHook is inserted at index 0 by
        // build() after all stages have been added — that way a single timer
        // covers the full pipeline run regardless of how many stages there are.
        self.finalize.push(Box::new(ClipFinalizeHook { metrics, progress_counter, timer }));

        Ok(())
    }

    /// Filter-specific step sequence:
    ///
    /// - If `filter.filter_by_template`: `GroupByQueryname` inserted before the process step.
    /// - No-rejects path: `ProcessOrdered` (single output → `DecompressedBlock`).
    /// - With-rejects path: `Process2Ordered` (two outputs: branch 0 kept, branch 1 rejects).
    ///   Branch 1 is wired here to its own `BgzfCompress → WriteBgzfFile` sink; branch 0
    ///   becomes `current_tail` for `add_sink` to wire the primary output.
    ///
    /// Only [`StagePosition::Terminal`] is supported; intermediate filter is not needed by
    /// any Phase 3 runall combination and bails with a clear error.
    ///
    /// Registers a `FilterFinalizeHook` for metrics reduction + summary banner.
    /// The chain-level [`StageTimingFinalizeHook`] is registered by
    /// [`Self::build`] (not here), so it captures the correct `Instant`.
    ///
    /// # Errors
    ///
    /// Returns errors if filter options are missing from the spec bag, if the
    /// rejects path is missing when rejects is enabled, or if `position` is
    /// `Intermediate` (not yet implemented).
    #[allow(clippy::too_many_lines)]
    fn add_filter(&mut self, position: StagePosition) -> Result<()> {
        use crate::commands::common::warn_unwired_pipeline_flags;
        use crate::logging::OperationTimer;
        use crate::pipeline::chains::commands::filter::{
            FilterFinalizeHook, build_filter_step_single_no_rejects,
            build_filter_step_single_with_rejects, build_filter_step_template_no_rejects,
            build_filter_step_template_with_rejects,
        };
        use crate::validation::validate_file_exists;
        use fgumi_bam_io::is_stdin_path;
        use log::info;

        if position == StagePosition::Intermediate {
            bail!(
                "intermediate filter not implemented; filter is typically terminal \
                (no Phase 3 runall combination requires it as an intermediate stage)"
            );
        }

        let filter = self
            .spec
            .stage_opts
            .filter
            .as_ref()
            .ok_or_else(|| anyhow!("Stage::Filter options missing from StageOptionsBag"))?;

        let tail = self.current_tail.expect("add_filter called before add_source");

        // Resolve input path for validation log messages.
        let input_path = self.resolve_log_input_path();

        // Validate input file exists (stdin exempt). The PairedBams path here
        // returns the mapped BAM, which `open_source` has already validated;
        // re-validating is cheap and harmless.
        if !is_stdin_path(&input_path) {
            validate_file_exists(&input_path, "Input BAM")?;
        }

        if let Some(ref reference) = filter.reference {
            validate_file_exists(reference, "Reference FASTA")?;
        }

        // Validate parameter counts (1-3 values for duplex support).
        filter.validate_parameters()?;

        let timer = OperationTimer::new("Filtering consensus reads");

        // Resolve output path for log.
        let output_path = match &self.spec.sink {
            SinkSpec::Bam(p) | SinkSpec::BamWithIndex(p) => p.clone(),
        };

        info!("Starting Filter");
        info!("Input: {}", input_path.display());
        info!("Output: {}", output_path.display());
        match &filter.reference {
            Some(r) => info!("Reference: {}", r.display()),
            None => info!("Reference: <none> (tag regeneration disabled)"),
        }
        info!("Min reads: {:?}", filter.min_reads);
        info!("Max read error rate: {:?}", filter.max_read_error_rate);
        info!("Max base error rate: {:?}", filter.max_base_error_rate);
        if let Some(q) = filter.min_base_quality {
            info!("Min base quality: {q}");
        }
        if let Some(q) = filter.min_mean_base_quality {
            info!("Min mean base quality: {q}");
        }
        info!("Max no-call fraction: {}", filter.max_no_call_fraction);
        if !filter.min_methylation_depth.is_empty() {
            info!("Min methylation depth: {:?}", filter.min_methylation_depth);
        }
        if filter.require_strand_methylation_agreement {
            info!("Require strand methylation agreement: true");
        }
        if let Some(frac) = filter.min_conversion_fraction {
            info!("Min conversion fraction: {frac}");
        }
        if let Some(mode) = &filter.methylation_mode {
            info!("Methylation mode: {mode:?}");
        }

        warn_unwired_pipeline_flags(&self.spec.scheduler, &self.spec.queue_memory);
        let num_threads = self.spec.threading.num_threads();
        info!("{}", self.spec.threading.log_message());
        info!("Using pipeline with {num_threads} threads");

        // Build shared FilterPipelineSetup (config, reference, accumulators, progress).
        let setup = filter.setup_pipeline(num_threads, &self.header)?;
        let accumulators = Arc::clone(&setup.collected_metrics);

        let filter_by_template = filter.filter_by_template;
        let track_rejects = filter.rejects.is_some();

        // If template-aware, insert GroupByQueryname step before process.
        let tail = if filter_by_template {
            use crate::pipeline::steps::group::queryname::GroupByQueryname;
            self.pipeline.append_step(GroupByQueryname::new(self.tuning.per_step_byte_limit), tail)
        } else {
            tail
        };

        // Select and append the process step (kept-only or kept+rejects).
        let process_tail = match (filter_by_template, track_rejects) {
            (false, false) => {
                let captures = filter.process_captures(&setup, &self.header);
                let step = build_filter_step_single_no_rejects(
                    self.tuning.per_step_byte_limit,
                    captures,
                    Arc::clone(&accumulators),
                );
                self.pipeline.append_step(step, tail)
            }
            (false, true) => {
                let captures = filter.process_captures(&setup, &self.header);
                let step = build_filter_step_single_with_rejects(
                    self.tuning.per_step_byte_limit,
                    captures,
                    Arc::clone(&accumulators),
                );
                self.pipeline.append_step(step, tail)
            }
            (true, false) => {
                let captures = filter.process_captures(&setup, &self.header);
                let step = build_filter_step_template_no_rejects(
                    self.tuning.per_step_byte_limit,
                    captures,
                    Arc::clone(&accumulators),
                );
                self.pipeline.append_step(step, tail)
            }
            (true, true) => {
                let captures = filter.process_captures(&setup, &self.header);
                let step = build_filter_step_template_with_rejects(
                    self.tuning.per_step_byte_limit,
                    captures,
                    Arc::clone(&accumulators),
                );
                self.pipeline.append_step(step, tail)
            }
        };
        // process_tail = (process_step_idx, BranchIdx(0)) = kept branch.

        // Wire the rejects branch (branch 1) to its own compress+write sink, when present.
        if track_rejects {
            use crate::pipeline::steps::bgzf::compress::BgzfCompress;
            use crate::pipeline::steps::sink::write_bgzf::WriteBgzfFile;

            let rejects_path = filter
                .rejects
                .as_ref()
                .ok_or_else(|| anyhow!("rejects path missing (build_for dispatch bug)"))?;
            let rejects_compress =
                BgzfCompress::new(self.tuning.compression_level, self.tuning.per_step_byte_limit);
            let rejects_writer =
                WriteBgzfFile::new(rejects_path, &self.header, self.tuning.compression_level)
                    .map_err(|e| anyhow!("WriteBgzfFile (rejects)::new: {e}"))?;

            // Branch 1 of process_tail is the rejects output.
            let rejects_branch = (process_tail.0, BranchIdx(1));
            let rejects_compress_tail = self.pipeline.append_step(rejects_compress, rejects_branch);
            // Sink: no tail needed after this.
            let _rejects_write_tail =
                self.pipeline.append_step(rejects_writer, rejects_compress_tail);
        }

        // Branch 0 of process_tail (kept) becomes current_tail for add_sink.
        self.current_tail = Some(process_tail);

        // Register the filter finalize hook.
        // The chain-level StageTimingFinalizeHook is inserted at index 0 by
        // build() after all stages have been added — that way a single timer
        // covers the full pipeline run regardless of how many stages there are.
        self.finalize.push(Box::new(FilterFinalizeHook {
            accumulators,
            stats_path: filter.stats.clone(),
            has_rejects: track_rejects,
            timer,
        }));

        Ok(())
    }

    /// Dedup-specific step sequence:
    ///
    /// `GroupByPosition::with_secondary_supplementary` →
    /// `ProcessOrdered` (dedup + metrics accumulation) →
    /// `MiAssign` (serial MI counter) →
    /// `ProcessWithWorkerState` (serialize records with MI tags) — **Terminal only**.
    ///
    /// For [`StagePosition::Terminal`], all four steps are appended and the
    /// chain tail is [`DecompressedBlock`] (bytes ready for `BgzfCompress`).
    ///
    /// For [`StagePosition::Intermediate`], only the first three steps are
    /// appended; the chain tail stays as [`BatchedProcessedDedupGroups`] for
    /// the next stage to consume. Intermediate dedup is not yet needed by any
    /// Phase 3 runall combination; calling `add_dedup` with `Intermediate`
    /// returns `Err` as a guard until a concrete fused path requires it.
    ///
    /// Registers a `DedupFinalizeHook` for metrics reduction + summary banner.
    /// The chain-level [`StageTimingFinalizeHook`] is registered by
    /// [`Self::build`] (not here), so it captures the correct `Instant`.
    ///
    /// # Errors
    ///
    /// Returns errors if the input is not template-coordinate sorted, if the
    /// dedup options are missing from the spec bag, or if `position` is
    /// `Intermediate` (not yet implemented).
    #[allow(clippy::too_many_lines)]
    fn add_dedup(&mut self, position: StagePosition) -> Result<()> {
        use crate::assigner::Strategy;
        use crate::commands::common::warn_unwired_pipeline_flags;
        use crate::logging::OperationTimer;
        use crate::per_thread_accumulator::PerThreadAccumulator;
        use crate::pipeline::chains::commands::dedup::{
            DedupFinalizeHook, build_mi_assign_step, build_process_step, build_serialize_step,
        };
        use crate::pipeline::steps::group::position::GroupByPosition;
        use crate::sam::SamTag;
        use crate::sam::is_template_coordinate_sorted;
        use log::info;

        if position == StagePosition::Intermediate {
            bail!(
                "intermediate dedup not implemented; dedup is typically terminal \
                (no Phase 3 runall combination requires it as an intermediate stage)"
            );
        }

        let dedup = self
            .spec
            .stage_opts
            .dedup
            .as_ref()
            .ok_or_else(|| anyhow!("Stage::Dedup options missing from StageOptionsBag"))?;

        let tail = self.current_tail.expect("add_dedup called before add_source");

        // Validate template-coordinate sort order.
        if !is_template_coordinate_sorted(&self.header) {
            bail!(
                "Input BAM must be template-coordinate sorted.\n\n\
                To prepare your BAM file, run:\n  \
                fgumi zipper -i mapped.bam -u unmapped.bam -r reference.fa -o merged.bam\n  \
                fgumi sort -i merged.bam -o sorted.bam --order template-coordinate"
            );
        }
        info!("Template-coordinate sorted");

        let timer = OperationTimer::new("Marking duplicates");

        let output_path = match &self.spec.sink {
            SinkSpec::Bam(p) | SinkSpec::BamWithIndex(p) => p.clone(),
        };

        info!("Starting dedup");
        info!("Output: {}", output_path.display());

        // Handle --no-umi mode: force identity strategy.
        let (effective_strategy, no_umi_edits_override) = if dedup.no_umi {
            if !matches!(dedup.strategy, Strategy::Identity) {
                info!("--no-umi mode: overriding strategy to identity");
            }
            (Strategy::Identity, true)
        } else {
            (dedup.strategy, false)
        };

        info!("Strategy: {effective_strategy:?}");

        let min_mapq: u8 = dedup.min_map_q.unwrap_or(0);
        let effective_edits =
            if no_umi_edits_override || matches!(effective_strategy, Strategy::Identity) {
                0
            } else {
                dedup.edits
            };

        info!("Edits: {effective_edits}");
        info!("Remove duplicates: {}", dedup.remove_duplicates);
        if dedup.no_umi {
            info!("No-UMI mode: deduplicating by position only");
        }
        if matches!(effective_strategy, Strategy::Adjacency | Strategy::Paired) {
            info!("Index threshold: {}", dedup.index_threshold);
        }

        warn_unwired_pipeline_flags(&self.spec.scheduler, &self.spec.queue_memory);
        let num_threads = self.spec.threading.num_threads();
        info!("{}", self.spec.threading.log_message());
        info!("Using pipeline with {num_threads} threads");

        let raw_tag = SamTag::RX;
        let assign_tag_bytes: [u8; 2] = *SamTag::MI;
        let filter_config = crate::commands::dedup::DedupFilterConfig {
            umi_tag: *raw_tag,
            min_mapq,
            include_non_pf: dedup.include_non_pf_reads,
            min_umi_length: dedup.min_umi_length,
            no_umi: dedup.no_umi,
        };

        let accumulators =
            PerThreadAccumulator::<crate::commands::dedup::CollectedDedupMetrics>::new(num_threads);
        let accumulators_for_process = Arc::clone(&accumulators);

        // ── Step factories (see chains::commands::dedup) ──────────────────
        let process_step = build_process_step(
            self.tuning.per_step_byte_limit,
            filter_config,
            effective_strategy,
            effective_edits,
            dedup.index_threshold,
            raw_tag,
            dedup.min_umi_length,
            dedup.no_umi,
            accumulators_for_process,
        );
        let mi_assign_step = build_mi_assign_step(self.tuning.per_step_byte_limit);
        let serialize_step = build_serialize_step(
            self.tuning.per_step_byte_limit,
            dedup.remove_duplicates,
            assign_tag_bytes,
            Arc::clone(&self.progress_records),
        );

        // Wire the dedup-specific steps.
        let tail = self.pipeline.append_step(
            GroupByPosition::with_secondary_supplementary(self.tuning.per_step_byte_limit),
            tail,
        );
        let tail = self.pipeline.append_step(process_step, tail);
        let tail = self.pipeline.append_step(mi_assign_step, tail);
        let tail = self.pipeline.append_step(serialize_step, tail);
        self.current_tail = Some(tail);

        // Register the dedup finalize hook.
        // The chain-level StageTimingFinalizeHook is inserted at index 0 by
        // build() after all stages have been added — that way a single timer
        // covers the full pipeline run regardless of how many stages there are.
        self.finalize.push(Box::new(DedupFinalizeHook {
            accumulators,
            metrics_path: dedup.metrics.clone(),
            family_size_histogram_path: dedup.family_size_histogram.clone(),
            timer,
        }));

        Ok(())
    }
}
