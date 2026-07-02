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
// Only the `consensus`-gated `FuseState` implements `HeapSize` in this module.
#[cfg(feature = "consensus")]
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

    /// The chain tail produces serialized BGZF-ready bytes
    /// ([`DecompressedBlock`]) after a Terminal serialize step
    /// (`SerializeBamRecords` / `SerializeGroups`, the terminal
    /// `SortMerge<BlockOutput>` for sort, or the record-aligned block emitted by
    /// a terminal consensus stage). No
    /// stage downstream of a Terminal stage reads `chain_tail_kind`, so this is
    /// a terminal-only marker; it exists so the kind is never a lie. `add_sink`
    /// wires `BgzfCompress → WriteBgzfFile` regardless of this value.
    ///
    /// [`DecompressedBlock`]: crate::pipeline::steps::types::DecompressedBlock
    SerializedBytes,

    /// The chain tail produces [`BamTemplateBatch`]. Set by `add_align`
    /// (Intermediate) and `add_correct` (Intermediate).
    ///
    /// [`BamTemplateBatch`]: crate::pipeline::steps::types::BamTemplateBatch
    BamTemplateBatch,

    /// The chain tail produces `BgzfBlock` (raw compressed BGZF blocks, NOT yet
    /// decompressed), with an arena-backed front waiting for the parallel-inflate
    /// sort ingest steps. Set by `add_source` when the first stage is a sort over
    /// a BAM source (all four orders): `ReadBgzfBlocks` is emitted but
    /// `BgzfDecompress` is NOT — the new front (`ReadBlocks` → `InflateToArena`
    /// → `FindBoundariesAndSort`) handles decompression + boundary-find + sort
    /// inline.
    ///
    /// `add_sort` consumes this tail for every BAM sort order.
    BgzfBlockArena,

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

/// Per-worker scratch state for the `templates_to_mi_step` bridge: a reusable
/// byte buffer and an MI-key string buffer, reset per template.
#[cfg(feature = "consensus")]
pub(crate) struct FuseState {
    scratch: Vec<u8>,
    mi_buf: String,
}

#[cfg(feature = "consensus")]
impl HeapSize for FuseState {}

/// The per-record filter applied by the duplex consensus stage before
/// MI-grouping: reject secondary/supplementary alignments and any record that
/// is unmapped without a mapped mate (the duplex caller works on primary,
/// alignment-anchored reads). Shared by both the non-fused
/// `GroupByMi::with_record_filter` path and the fused `templates_to_mi_step`
/// bridge so the two paths cannot drift (NEW-002).
#[cfg(feature = "consensus")]
fn duplex_record_filter(raw: &[u8]) -> bool {
    use fgumi_raw_bam::RawRecordView;
    let flg = RawRecordView::new(raw).flags();
    if flg & fgumi_raw_bam::flags::SECONDARY != 0 || flg & fgumi_raw_bam::flags::SUPPLEMENTARY != 0
    {
        return false;
    }
    let is_mapped = flg & fgumi_raw_bam::flags::UNMAPPED == 0;
    let has_mapped_mate =
        flg & fgumi_raw_bam::flags::PAIRED != 0 && flg & fgumi_raw_bam::flags::MATE_UNMAPPED == 0;
    is_mapped || has_mapped_mate
}

/// Build the `TemplatesToMiGroups` bridge step that fuses an intermediate
/// `add_group` (which emits `BatchedProcessedPositionGroups`) into a consensus
/// stage (which consumes `BatchedMiGroups`). For each template it splices the
/// assigned molecular identifier into the `MI` tag and runs the records into
/// per-MI groups in a single pass. This is the in-builder successor of the
/// former `runall.rs` `templates_to_mi_step` (folded into the chain builder).
///
/// Grouping on the MI alone (without a cell-barcode partition) is exactly
/// equivalent to the non-fused `GroupByMi::with_cell_tag(Some(CB))` path
/// (S5b1-003): MI values are globally unique across cell barcodes. The cell
/// barcode is part of the upstream `GroupByPosition` key (so reads from
/// different cells land in different position groups), and `assign_mi_offsets`
/// shifts every group's local MI ids by a single global monotonic counter, so
/// every distinct molecule receives a unique integer MI regardless of cell.
/// Two templates with the same MI but different CB therefore cannot exist, so
/// the bridge needs no cell-tag partition.
///
/// `duplex_strip_strand` controls whether the `/A` `/B` strand suffix is
/// stripped from the MI key so both strands of a duplex molecule group
/// together (`true` for duplex; `false` for simplex and codec).
///
/// `record_filter` is an optional per-record predicate applied to each emitted
/// `r1()`/`r2()` raw record; records for which it returns `false` are dropped.
/// Duplex passes [`duplex_record_filter`] here so the fused
/// group→duplex path applies the **same** record filter that the non-fused
/// `GroupByMi::with_record_filter` applies (NEW-002); simplex and codec pass
/// `None` (they install no record filter on either path).
#[cfg(feature = "consensus")]
fn templates_to_mi_step<F>(
    limit_bytes: u64,
    duplex_strip_strand: bool,
    record_filter: Option<F>,
) -> impl crate::pipeline::core::step::Step<
    Input = crate::pipeline::steps::group::position::BatchedProcessedPositionGroups,
    Outputs = crate::pipeline::core::outputs::OrderedBytesSingle<
        crate::pipeline::steps::group::mi::BatchedMiGroups,
    >,
>
where
    F: Fn(&[u8]) -> bool + Send + Sync + 'static,
{
    use crate::mi_group::MiGroup;
    use crate::pipeline::steps::group::mi::BatchedMiGroups;
    use crate::pipeline::steps::group::position::BatchedProcessedPositionGroups;
    use crate::pipeline::steps::process::process_with_worker_state;
    use crate::sam::SamTag;
    use fgumi_raw_bam::RawRecord;
    use std::io;

    // The assign tag is always `MI` (== [b'M', b'I']), matching `add_group`.
    let assign_tag_bytes: [u8; 2] = *SamTag::MI;

    process_with_worker_state::<BatchedProcessedPositionGroups, BatchedMiGroups, _, FuseState, _>(
        "TemplatesToMiGroups",
        limit_bytes,
        || FuseState { scratch: Vec::with_capacity(512), mi_buf: String::with_capacity(16) },
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
                    // Write the full strand-suffixed MI into `mi_buf` for the
                    // per-record `MI` tag exactly once. On the strip path the
                    // grouping key is the strand-stripped base id (so both
                    // strands group together); on the non-strip path the full
                    // MI is also the grouping key (reuse `mi_buf`).
                    template.mi.write_with_offset(0, &mut state.mi_buf);
                    let mi_key = if duplex_strip_strand {
                        template.mi.id().unwrap().to_string()
                    } else {
                        state.mi_buf.clone()
                    };
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
                        // Apply the optional per-record filter (duplex) so the
                        // fused bridge drops the same records the non-fused
                        // `GroupByMi::with_record_filter` would (NEW-002).
                        if let Some(filter) = &record_filter {
                            if !filter(raw) {
                                continue;
                            }
                        }
                        state.scratch.clear();
                        state.scratch.extend_from_slice(raw);
                        fgumi_raw_bam::update_string_tag(
                            &mut state.scratch,
                            assign_tag_bytes,
                            state.mi_buf.as_bytes(),
                        );
                        // Move the tagged bytes into the record instead of cloning
                        // them; reinstall a fresh scratch buffer pre-sized to the
                        // same capacity so the next record reuses it without a
                        // realloc (the per-record `clear()` above keeps it warm).
                        let cap = state.scratch.capacity();
                        let record = std::mem::replace(&mut state.scratch, Vec::with_capacity(cap));
                        current_records.push(RawRecord::from(record));
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
    )
}

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
    /// Always drained, even when the run fails (see [`BuiltPipeline::run`]).
    finalize: Vec<Box<dyn FinalizeHook>>,

    /// Post-pipeline hooks that run **only on a fully successful run**.
    /// Populated by `add_<stage>` methods for actions that publish a derived
    /// artifact from the run's output (currently only `add_sort`'s
    /// `IndexBamFinalizeHook`, which writes the `.bai` sidecar by re-reading
    /// the finished BAM). Gating these mirrors the standalone `fgumi sort`
    /// flow, where the index write is behind `run_result?`: a failed run
    /// leaves a partial BAM, so publishing its index would be stale.
    finalize_on_success: Vec<Box<dyn FinalizeHook>>,

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
    /// Set by stages whose framework thread count differs from the requested
    /// `--threads` (e.g. [`Self::add_zipper`] / [`Self::add_align`], which
    /// reserve a driver thread for an external aligner process). Sort no longer
    /// sets it: standalone sort streams through the normal multi-step pipeline
    /// and uses the requested thread count like every other stage.
    override_pipeline_threads: Option<usize>,

    /// Set when a BAM sort source is wired through the parallel-inflate arena
    /// front (`ReadBlocks → InflateToArena → FindBoundariesAndSort`), for any
    /// sort order. [`Self::build`] then opts the ENTIRE sort-first pipeline —
    /// including fused `Sort → Group/Simplex/Duplex/Codec` chains — into the
    /// downstream-first [`DrainFirstScheduler`], so the sort's serial
    /// boundary/key scan overlaps the parallel inflate (drained ahead of
    /// production) instead of starving behind it on the shared pool. Chains
    /// without a sort source keep the default upstream-first scheduler.
    use_drain_first_scheduler: bool,

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
    /// - `add_source` sets `BgzfBlockArena` when Sort is the first stage over a
    ///   BAM source (leaves raw BGZF blocks for the parallel-inflate arena front
    ///   instead of `DecodeRecords`; SAM sort-first input is rejected).
    /// - `add_correct(Intermediate)` sets `BamTemplateBatch`.
    /// - `add_align(Intermediate)` sets `BamTemplateBatch`.
    /// - `add_sort(Intermediate)` reads the kind to decide its preamble,
    ///   then sets `DecodedRecordBatch` after `DecodeFromRecords`.
    /// - `add_group(Intermediate)` sets `BatchedProcessedPositionGroups`.
    /// - Consensus stages read the kind to decide their grouping preamble
    ///   (`GroupByMi` for `DecodedRecordBatch`, `TemplatesToMiGroups` for
    ///   `BatchedProcessedPositionGroups`).
    chain_tail_kind: ChainTailKind,

    /// Lever 2: run the terminal `WriteBgzfFile` on its own dedicated
    /// `StepKind::Detached` thread instead of as a pool-scheduled
    /// `Serial + Affinity::Writer` step. Set ONLY by `add_sort` on the
    /// standalone-sort terminal (`spec.is_sort_terminal()`), so every other
    /// chain (runall / group / consensus / zipper / correct / dedup / filter /
    /// clip) keeps the pool-scheduled writer. `add_sink` reads this flag.
    /// Default `false`.
    detached_writer: bool,
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

        // Every chain — including the sole-`[Stage::Sort]` chain — opens its
        // source exactly once here (header read; the reader continues in
        // `add_source` for the body, which is the single-open contract stdin
        // relies on). The former sort-terminal skip (where `SortBamFile` opened
        // the input itself) is gone: standalone sort now streams through the
        // normal source path, so `@PG` injection applies to it too.
        let (raw_header, pending_source) = Self::open_source(spec)?;
        let header = crate::commands::common::add_pg_record(raw_header, &spec.command_line)?;

        Ok(Self {
            spec,
            tuning,
            header,
            current_tail: None,
            pipeline: PipelineBuilder::new(),
            finalize: Vec::new(),
            finalize_on_success: Vec::new(),
            progress_records: Arc::new(AtomicU64::new(0)),
            pending_source,
            paired_tail: None,
            override_pipeline_threads: None,
            use_drain_first_scheduler: false,
            pending_header_handle: None,
            // Initialise to the default; add_source will set the correct kind
            // based on whether sort is the first intermediate stage.
            chain_tail_kind: ChainTailKind::DecodedRecordBatch,
            // Pool-scheduled writer by default; add_sort opts the standalone
            // sort terminal into a Detached writer (lever 2).
            detached_writer: false,
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
        use crate::pipeline::steps::source::InputSource;
        use crate::pipeline::steps::source::read_bam::read_bam_from_reader;
        use crate::sam::check_sort;

        let source = self.pending_source.take().expect("add_source called twice");

        // Detect whether Sort is the FIRST stage (whether or not it is the
        // only one). When it is over a BAM source, the preamble stops at the
        // raw `ReadBgzfBlocks` step (tail = `ChainTailKind::BgzfBlockArena`)
        // instead of `DecodeRecords → DecodedRecordBatch`, so `add_sort`'s
        // parallel-inflate arena front (`ReadBlocks → InflateToArena →
        // FindBoundariesAndSort`) can consume it directly without an extra
        // round-trip. SAM sort-first is rejected outright (see the `Sam` arm
        // below); `SortBuffer` is only reached by the non-arena path in
        // `add_sort`, where sort is *not* first (a fused earlier stage feeds it
        // `DecodedRecordBatch` / `BamTemplateBatch`). This covers both the
        // sole-`[Sort]` chain (standalone `fgumi sort`) and
        // sort-as-first-intermediate (`--start-from sort --stop-after
        // {group,consensus}`); the former used to be special-cased out of
        // `add_source` entirely.
        let sort_is_first_intermediate = matches!(self.spec.stages.as_slice(), [Stage::Sort, ..]);

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
                            // All four sort orders use the arena-backed
                            // parallel-inflate front (ReadBgzfBlocks only — no
                            // BgzfDecompress); ReadBlocks → InflateToArena →
                            // FindBoundariesAndSort (with the per-order strategy)
                            // are wired in add_sort. Leave the BgzfBlock tail for
                            // that front. SAM / fused inputs still use SortBuffer.
                            let (read_step, _) = read_bam_from_reader(
                                reader,
                                self.header.clone(),
                                self.tuning.blocks_per_batch,
                                self.tuning.per_step_byte_limit,
                            );
                            let tail = self.pipeline.append_source(read_step);
                            self.current_tail = Some(tail);
                            self.chain_tail_kind = ChainTailKind::BgzfBlockArena;
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
                        // Sort-as-first-stage needs a `RecordBatch` source, but the
                        // SAM preamble only produces `DecodedRecordBatch`; feeding
                        // that into the sort ingest (`SortBuffer`) is a handle-type
                        // mismatch (see `add_sort`). SAM input to a sort-first chain was never
                        // supported (the legacy file→file sorter read BAM only), so
                        // reject it with a clear message instead of a downstream
                        // typed-step panic.
                        if sort_is_first_intermediate {
                            anyhow::bail!(
                                "sort requires BAM input; SAM input is not supported \
                                 (convert with `samtools view -b` first)"
                            );
                        }
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
        let config = fgumi_bam_io::GroupKeyConfig::new(library_index, cell_tag);
        // When a Group stage is in the chain, cache the UMI (RX) value position
        // during decode so `fgumi group`'s per-template UMI-assignment lookup can
        // slice it without re-scanning aux data (#334). Group's UMI tag is fixed
        // to RX (see `add_group`'s `raw_tag`). Chains without a Group stage skip
        // the cache to avoid the extra per-record aux scan.
        if self.spec.stages.contains(&Stage::Group) {
            config.with_umi_tag(*SamTag::RX)
        } else {
            config
        }
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
    /// (`add_sink` wires `BgzfCompress → Write`); `chain_tail_kind` is set to
    /// [`ChainTailKind::SerializedBytes`] (the honest marker for "serialized
    /// bytes ready for sink"; nothing downstream of a terminal consensus reads
    /// it). For a
    /// [`StagePosition::Intermediate`] consensus stage, a downstream stage
    /// (e.g. filter) consumes records, so append the existing [`DecodeRecords`]
    /// step — consensus output is already record-aligned, so no
    /// `FindBamBoundaries`/`ParseBamRecords` is needed — yielding
    /// [`DecodedRecordBatch`].
    ///
    /// [`DecompressedBlock`]: crate::pipeline::steps::types::DecompressedBlock
    /// [`DecodedRecordBatch`]: crate::pipeline::steps::types::DecodedRecordBatch
    #[cfg(feature = "consensus")]
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
            // terminal: tail is DecompressedBlock (serialized bytes) → SerializedBytes.
            StagePosition::Terminal => {
                self.chain_tail_kind = ChainTailKind::SerializedBytes;
                tail
            }
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

    /// Wire the rejects fan-out branch (branch 1) of a 2-output consensus step:
    /// `BgzfCompress` → `WriteBgzfFile` with the **input** header, in input
    /// order — the PR #332 fan-out contract (rejects flow through the ordered
    /// serialize/compress stages, not a mutex side-channel).
    ///
    /// `consensus_pt` is the just-appended consensus step (branch 0 is its kept
    /// output, consumed by the caller); `label` names the command for error
    /// context. Shared by `add_simplex`/`add_duplex`/`add_codec`, whose
    /// rejects wiring is otherwise identical — only the consensus step's type
    /// (built inline by the caller) differs.
    #[cfg(feature = "consensus")]
    fn wire_consensus_rejects_branch(
        &self,
        consensus_pt: (StepIdx, BranchIdx),
        rejects_path: Option<&std::path::Path>,
        input_header: &Header,
        label: &str,
    ) -> Result<()> {
        use crate::pipeline::steps::bgzf::compress::BgzfCompress;
        use crate::pipeline::steps::sink::write_bgzf::WriteBgzfFile;

        let rejects_path = rejects_path
            .ok_or_else(|| anyhow!("rejects path unexpectedly None when track_rejects is set"))?;
        let rejects_write =
            WriteBgzfFile::new(rejects_path, input_header, self.tuning.compression_level)
                .map_err(|e| anyhow!("WriteBgzfFile ({label} rejects): {e}"))?;

        let rejects_branch = (consensus_pt.0, BranchIdx(1));
        let rejects_compress_tail = self.pipeline.append_step(
            BgzfCompress::new(self.tuning.compression_level, self.tuning.per_step_byte_limit),
            rejects_branch,
        );
        self.pipeline.append_step(rejects_write, rejects_compress_tail);
        Ok(())
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
        // Lever 2: on the standalone-sort terminal only, run the writer on the
        // sort's I/O driver thread (legacy "N + 2"), sharing it with the phase-1
        // `SpillWrite`. Every other chain keeps the pool-scheduled Serial +
        // Affinity::Writer writer.
        let write_step = if self.detached_writer {
            use crate::pipeline::core::step::DetachedGroup;
            use crate::pipeline::steps::sort::SORT_IO_GROUP;
            write_step.with_detached(DetachedGroup::Shared(SORT_IO_GROUP))
        } else {
            write_step
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
        // Opt the whole sort-first pipeline into downstream-first dispatch so the
        // arena front's serial scan overlaps the parallel inflate (flag set in
        // add_sort for any sort order). Chains without a sort source keep the
        // default upstream-first scheduler.
        if self.use_drain_first_scheduler {
            config = config.with_scheduler(std::sync::Arc::new(
                crate::pipeline::core::runtime::DrainFirstScheduler,
            ));
        }
        // DIAGNOSTIC: the stats `Arc` already exists whenever the deadlock
        // monitor is on (default), so allow `FGUMI_PIPELINE_STATS=1` to dump the
        // per-step timing report even on paths (e.g. standalone `fgumi sort`)
        // that do not flatten `--pipeline-stats`.
        // `--pipeline-trace` (any instrumentation level above Off) also produces a
        // stats `Arc` in `build_pipeline_config_for_chain`; without this clause a
        // trace-only run would collect edge metrics but never emit the end-of-run
        // report, since neither `collect_stats()` nor the env var is set.
        let user_wants_stats = self.spec.scheduler.collect_stats()
            || super::build_helpers::env_flag_enabled("FGUMI_PIPELINE_STATS")
            || self.spec.scheduler.instrumentation_level().is_on();
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

        Ok(BuiltPipeline {
            pipeline,
            config,
            finalize: self.finalize,
            finalize_on_success: self.finalize_on_success,
        })
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
            // SerializeBamRecords emits DecompressedBlock (serialized bytes), not
            // DecodedRecordBatch. Mark the tail honestly as SerializedBytes, matching
            // every other terminal serialize path (e.g. add_dedup); nothing downstream
            // of a Terminal stage reads this, but the kind must never lie.
            self.chain_tail_kind = ChainTailKind::SerializedBytes;
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
    /// For [`StagePosition::Intermediate`], `SerializeBamRecords` is **not**
    /// appended: the chain tail is left as `BamTemplateBatch` so the next stage
    /// (`add_align` → `AlignAndMergeStep`) can consume the correct step's kept
    /// output (branch 0) directly; `add_align` skips `GroupByQueryname` on an
    /// incoming `BamTemplateBatch`. This is the correct→align fused path.
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
    /// Returns errors if correct options are missing from the spec bag or if UMI
    /// sequence loading fails.
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
            // tail is DecompressedBlock (serialized bytes) → SerializedBytes.
            self.chain_tail_kind = ChainTailKind::SerializedBytes;
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
            // tail is DecompressedBlock (serialized bytes) after SerializeBamRecords
            // → SerializedBytes.
            self.chain_tail_kind = ChainTailKind::SerializedBytes;
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
            // tail is DecompressedBlock (serialized bytes) after SerializeBamRecords
            // → SerializedBytes.
            self.chain_tail_kind = ChainTailKind::SerializedBytes;
        } else {
            self.current_tail = Some(merge_tail);
            // Intermediate: ZipperMergeStep emits BamTemplateBatch so the
            // next stage (typically add_sort) prepends TemplatesToRecordBatch.
            self.chain_tail_kind = ChainTailKind::BamTemplateBatch;
        }

        // Set the thread override to the floored count so build() forwards it
        // to PipelineConfig::threads. Take the max with any prior override (e.g.
        // a preceding align stage) so we never lower an earlier stage's floor
        // (S5b1-004).
        if let Some(prev) = self.override_pipeline_threads {
            self.override_pipeline_threads = Some(prev.max(num_threads));
        } else {
            self.override_pipeline_threads = Some(num_threads);
        }

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

    /// Sort-specific step sequence — a single streaming pipeline used for
    /// every sort, whether standalone (`[Stage::Sort]`), intermediate, or
    /// terminal-after-upstream-stages.
    ///
    /// `add_source` always runs before this method (even for a sole-`[Sort]`
    /// chain), so `current_tail` is always `Some`. Sort emits the in-pipeline
    /// step sequence:
    ///
    /// ```text
    /// [BamTemplateBatch or RecordBatch at current_tail]
    ///     ↓ (TemplatesToRecordBatch if BamTemplateBatch)
    /// RecordBatch (SAM) │ BgzfBlockArena (BAM)
    ///     ↓ arena sort ingest:
    ///        SAM → SortBuffer (Serial)
    ///        BAM → ReadBlocks → InflateToArena → FindBoundariesAndSort
    ///     ↓ SpillGather (Serial) → SpillCompress (Parallel) → SpillWrite (Serial)
    /// SortPhase1Event
    ///     ↓ SortSpillDecompress (Parallel)
    /// SortPhase2Event
    ///     ↓ [Terminal]     SortMerge<BlockOutput> (Serial) → DecompressedBlock
    ///     ↓ [Intermediate] SortMerge (Serial) → RecordBatch → DecodeFromRecords
    /// ```
    ///
    /// For a `SinkSpec::BamWithIndex` (terminal sort only), an
    /// [`IndexBamFinalizeHook`] is registered on the success-gated finalize
    /// list. `--sort::max-memory=auto` is only honoured for a sole-`[Sort]`
    /// chain (it owns the whole budget); a fused chain must pass a fixed value.
    ///
    /// The chain-level [`StageTimingFinalizeHook`] (registered by [`Self::build`])
    /// covers the wall time; sort does not register a per-stage summary hook.
    ///
    /// # Errors
    ///
    /// Returns errors if sort options are missing from the spec bag, or if
    /// `--sort::max-memory=auto` is requested in a fused multi-stage pipeline.
    ///
    /// [`IndexBamFinalizeHook`]: crate::pipeline::chains::commands::sort::IndexBamFinalizeHook
    #[allow(clippy::too_many_lines)]
    fn add_sort(&mut self, position: StagePosition) -> Result<()> {
        use crate::commands::common::{MemoryLimit, resolve_memory_budget};
        use crate::pipeline::chains::commands::sort::IndexBamFinalizeHook;

        let sort = self
            .spec
            .stage_opts
            .sort
            .as_ref()
            .ok_or_else(|| anyhow!("Stage::Sort options missing from StageOptionsBag"))?
            .clone();

        {
            // ── Streaming sort (standalone, Intermediate, OR Terminal-after-upstream) ──
            //
            // Every sort — including a sole-`[Stage::Sort]` standalone sort —
            // runs through the streaming source → arena sort ingest →
            // SpillGather → SpillCompress → SpillWrite → SortSpillDecompress →
            // SortMerge → sink pipeline. `add_source` always runs first (even
            // for `[Sort]`), so `current_tail` is always `Some` here.
            //
            // Chain topology (continuing from current_tail):
            //
            //   [BamTemplateBatch or RecordBatch at current_tail]
            //       ↓ (TemplatesToRecordBatch if BamTemplateBatch)
            //   RecordBatch (SAM) │ BgzfBlockArena (BAM)
            //       ↓ arena sort ingest (SortBuffer, or ReadBlocks →
            //         InflateToArena → FindBoundariesAndSort for BAM)
            //       ↓
            //   SpillGather (Serial) → SpillCompress (Parallel) → SpillWrite
            //       (Serial) → SortPhase1Event
            //       ↓
            //   SortSpillDecompress (Parallel)
            //       ↓ SortPhase2Event
            //   SortMerge<RecordBatchOutput> (Serial) → RecordBatch      [Intermediate]
            //   SortMerge<BlockOutput>       (Serial) → DecompressedBlock [Terminal]
            //       ↓
            //   [Intermediate] DecodeFromRecords (Parallel) → DecodedRecordBatch
            //   [Terminal]     (none — SortMerge<BlockOutput> emits DecompressedBlock directly)
            //
            // For Intermediate: chain_tail_kind = DecodedRecordBatch;
            // add_group / add_simplex / etc. can proceed normally.
            //
            // For Terminal: chain_tail_kind = SerializedBytes (DecompressedBlock
            // bytes ready for sink) — add_sink appends BgzfCompress → WriteBgzfFile.
            //
            // `Auto` memory is honored only for a sole-`[Sort]` chain (see
            // `is_standalone_sort` below), where sort owns the whole budget. In a
            // multi-stage pipeline the budget is shared across stages, so `Auto`
            // cannot be resolved here and a fixed value is required.
            use crate::pipeline::core::step::Affinity;
            use crate::pipeline::steps::parse::decode::DecodeFromRecords;
            use crate::pipeline::steps::sort::{
                BlockOutput, RecordBatchOutput, SortBuffer, SortDecompressTuning, SortMerge,
                SortSpillDecompress, SpillCompress, SpillGather, SpillWrite,
            };
            use fgumi_sort::{RawExternalSorter, SortOrder};

            // `auto` is honorable only for a sole-`[Sort]` chain (standalone
            // sort): it owns the whole budget. In a fused chain sort shares the
            // budget with the other stages, so `auto` can't be split — require
            // a fixed value there.
            // Single-sourced from the spec classifier (sole-`[Sort]` layout).
            let is_standalone_sort = self.spec.is_sort_terminal();

            // Lever 2: on the standalone-sort terminal, run the output writer on
            // its own dedicated thread (legacy "N + 2"). Gated to the standalone
            // sort terminal ONLY — any intermediate sort, or a sort inside a
            // fused runall chain, keeps the pool-scheduled writer so no streaming
            // command's terminal is affected. `add_sink` reads this flag.
            if is_standalone_sort && matches!(position, StagePosition::Terminal) {
                self.detached_writer = true;
            }
            if matches!(sort.max_memory, MemoryLimit::Auto) && !is_standalone_sort {
                bail!(
                    "--sort::max-memory=auto is not supported in a fused multi-stage pipeline \
                     (sort is not standalone); pass a fixed value (e.g. --sort::max-memory 768M)"
                );
            }
            let num_threads = self.spec.threading.num_threads();
            // Per-phase worker counts default to --threads; runall forwards
            // --sort::sort-threads / --sort::merge-threads via SortOptions.
            // Phase 1 caps the streaming sort/spill worker pool; Phase 2 caps
            // concurrent spill-decompression (the merge step itself is serial).
            let num_phase1_threads = resolve_phase_threads(sort.sort_threads, num_threads);
            let num_phase2_threads = resolve_phase_threads(sort.merge_threads, num_threads);
            if sort.sort_threads.is_some() || sort.merge_threads.is_some() {
                log::debug!(
                    "streaming sort: per-phase thread split (phase1={num_phase1_threads}, \
                     phase2={num_phase2_threads}); phase2 caps decompress concurrency, merge is serial"
                );
            }
            // Reuse the same helper as the standalone path so `--sort::memory-per-thread`
            // is honored (the previous inline `× num_threads` ignored it).
            let total_memory = resolve_memory_budget(
                sort.max_memory,
                sort.memory_reserve,
                num_threads,
                sort.memory_per_thread,
            )?;

            // Honor the requested order (standalone sort can request any of the
            // four; runall always passes `TemplateCoordinate`). The cell tag is
            // only meaningful for cell-grouped orders — `parse_cell_tag` returns
            // `Some(CB)` for template-coordinate and `None` otherwise — so it is
            // applied conditionally rather than hardcoded.
            let sort_order: SortOrder = sort.order.into();
            let cell_tag = crate::commands::sort::parse_cell_tag(sort.order)?;
            // Honor the requested spill codec (`SortSpillDecompress` handles both
            // BGZF and zstd). Standalone sort may request `bgzf` to pair with
            // `--temp-compression 0` (uncompressed spill); runall leaves it at
            // the zstd default. Omitting this previously forced the sorter's
            // default codec, breaking the uncompressed-bgzf spill path.
            let mut sorter = RawExternalSorter::new(sort_order)
                .memory_limit(total_memory)
                .threads(num_threads.max(1))
                .sort_threads(num_phase1_threads)
                .merge_threads(num_phase2_threads)
                .output_compression(1)
                .temp_compression(sort.temp_compression)
                .spill_codec(sort.temp_codec);
            if let Some(ct) = cell_tag {
                sorter = sorter.cell_tag(ct);
            }
            // Honor `--key-types` (template-coordinate only): the arena template
            // accumulator (`TemplateArenaAccumulator`) provisions the dropped-lane
            // variant from this and validates that dropped lanes are constant.
            // Previously the streaming production path ignored it entirely.
            if let Some(kt) = sort.key_types {
                sorter = sorter.key_types(kt);
            }

            if !sort.tmp_dirs.is_empty() {
                sorter = sorter.temp_dirs(sort.tmp_dirs.clone());
            }

            // NOTE: the legacy `build_sort_step` clamped the sorter's *initial*
            // in-memory capacity (768 MiB/thread) for standalone `--max-memory
            // auto`, relying on the old `RawExternalSorter`'s lazy buffer growth.
            // That knob does not carry over to the arena front used below
            // (`ReadBlocks`): it sizes its arena segment to the FULL budget
            // (`run_cap = memory_limit`) so a run that fits `--max-memory` sorts
            // in memory with zero spills (see `ReadBlocks::new`). Capping the
            // segment would cap `run_cap` and force premature spills — the exact
            // wall-clock regression that doc warns against — so there is no
            // separate initial-capacity clamp here.

            let affinity =
                if num_threads.max(2) >= 3 { Affinity::Worker(1) } else { Affinity::Reader };

            // Phase-2 decompression granularity knobs (`--sort::file-granularity`,
            // `--sort::block-batch`). The default is the block-parallel path
            // (P5c, hardened by loom/TSan + the soak matrix); its reorder window
            // is bounded by the per-step byte limit (see SortSpillDecompress).
            // `block_batch` defaults to 4 (the original MAX_BATCH_PER_CALL),
            // pending a fleet decompress-throughput bench.
            let decompress_tuning = SortDecompressTuning {
                file_granularity: sort.file_granularity,
                block_batch: sort.block_batch,
            };
            // `--merge-threads` (Phase-2) caps how many decompress workers run
            // concurrently; the block-parallel path still does the work, but the
            // admission counter bounds concurrency to `num_phase2_threads`.
            let decompress =
                SortSpillDecompress::new(self.tuning.per_step_byte_limit, decompress_tuning)
                    .with_max_concurrency(Some(num_phase2_threads));
            // Standalone sort gets an end-of-run summary (records processed /
            // written / temporary chunks); the fused runall path does not (the
            // chain-level timing hook covers it). The slot is filled by
            // `SortMerge` on completion and read by `SortSummaryFinalizeHook`.
            let sort_stats_slot = is_standalone_sort
                .then(|| Arc::new(parking_lot::Mutex::new(None::<fgumi_sort::SortStats>)));
            // `SortMerge` is constructed inside the Terminal/Intermediate branch
            // below, because its output-framing strategy differs by position:
            //
            // - Terminal: `SortMerge<BlockOutput>` frames each merged record as
            //   `[u32 LE block_size][body]` directly into a `DecompressedBlock`,
            //   wired straight to `BgzfCompress` — the former
            //   `SortMerge → SerializeRecordBatch → BgzfCompress` triple folded to
            //   `SortMerge → BgzfCompress` (lever 1: one fewer pool step, one
            //   fewer reorder stage, one fewer memcpy per record). The framing is
            //   byte-for-byte identical to the removed `SerializeRecordBatch`.
            // - Intermediate: `SortMerge<RecordBatchOutput>` (the default) emits
            //   `RecordBatch` for the downstream `DecodeFromRecords`.

            let tail = self.current_tail.expect("streaming sort called before add_source");

            // If the upstream stage produced BamTemplateBatch (from Align or
            // Correct), flatten it to RecordBatch before feeding the sort ingest.
            let tail = if self.chain_tail_kind == ChainTailKind::BamTemplateBatch {
                use crate::pipeline::steps::templates_to_records::TemplatesToRecordBatch;
                self.pipeline
                    .append_step(TemplatesToRecordBatch::new(self.tuning.per_step_byte_limit), tail)
            } else {
                // chain_tail_kind == RecordBatch (from ParseBamRecords in add_source)
                // or DecodedRecordBatch would be a bug caught at runtime.
                tail
            };

            // Phase-1 head — the block-parallel spill-write split for ALL sort
            // orders: `SortBuffer` (Serial: ingest + par-sort + emit sorted
            // chunks) → `SpillGather` (Serial: fan each chunk into raw blocks,
            // mint a dense ordinal) → `SpillCompress` (Parallel + ByItemOrdinal:
            // compress blocks across the work-stealing pool) → `SpillWrite`
            // (Serial + Writer: demux blocks to per-file spill files, emit
            // `SortPhase1Event`). This mirrors the output BAM tail
            // (Serialize → BgzfCompress → WriteBgzfFile) so a single spill's
            // compression saturates cores instead of pinning one worker, while
            // keeping `SortBuffer` Serial (no N² rayon oversubscription). The head
            // emits `SortPhase1Event`, so `SortSpillDecompress`/`SortMerge` are
            // identical downstream regardless of order.
            let phase1_tail = {
                let (temp_dirs, alloc) = sorter
                    .create_spill_dirs()
                    .map_err(|e| anyhow!("create spill temp dirs: {e}"))?;
                let temp_codec = sort.temp_codec;
                let temp_compression = sort.temp_compression;
                let temp_dirs = Arc::new(temp_dirs);
                let alloc = Arc::new(parking_lot::Mutex::new(alloc));
                let byte_limit = self.tuning.per_step_byte_limit;
                // Sort head: two cases depending on the source tail kind.
                // BgzfBlockArena (BAM): use the parallel-inflate front
                // (ReadBlocks → InflateToArena → FindBoundariesAndSort).
                // Else (RecordBatch from SAM or upstream stage): use SortBuffer.
                let head_tail = if self.chain_tail_kind == ChainTailKind::BgzfBlockArena {
                    // BAM sort through the parallel-inflate arena front
                    // (ReadBlocks → InflateToArena → FindBoundariesAndSort).
                    // Each order selects its per-order strategy.
                    use fgumi_pipeline_io::sort::protocol::MemoryChunkErased;
                    use fgumi_pipeline_io::sort::{
                        CoordinateStrategy, FindBoundariesAndSort, InflateToArena,
                        QuerynameStrategy, ReadBlocks, TemplateStrategy,
                    };
                    use fgumi_sort::{QuerynameComparator, RawQuerynameKey, RawQuerynameLexKey};
                    let n_ref = u32::try_from(self.header.reference_sequences().len())
                        .map_err(|_| anyhow!("reference sequence count overflows u32"))?;
                    // sorter is not consumed by this branch; the spill dirs/alloc it
                    // owns were already taken above (create_spill_dirs). Drop it to
                    // silence the unused-variable warning.
                    drop(sorter);
                    // `total_memory` is the full in-memory budget (--max-memory ×
                    // threads). ReadBlocks sizes its arena segment to hold one
                    // full-budget run, so data that fits the budget sorts entirely
                    // in memory with zero spills — matching legacy.
                    let read_blocks = ReadBlocks::new(total_memory, byte_limit);
                    let inflate = InflateToArena::new(byte_limit);
                    // Hand the resolved Phase-1 count (honoring `--sort-threads`,
                    // else `--threads`) to the per-chunk sort so large chunks use
                    // the parallel radix. Using the raw global `num_threads()` here
                    // silently dropped the `--sort-threads` override.
                    let sort_threads = num_phase1_threads;
                    let t = self.pipeline.append_step(read_blocks, tail);
                    let t = self.pipeline.append_step(inflate, t);
                    let out = match sort_order {
                        SortOrder::TemplateCoordinate => {
                            // Template-coordinate: the accumulator carries the
                            // library / cell-barcode / MI and `--key-types`
                            // narrowed-lane machinery the template key needs, built
                            // from the header exactly as the legacy
                            // TemplateChunkSorter.
                            let key_types = sort.key_types.unwrap_or_default();
                            let acc = fgumi_sort::TemplateArenaAccumulator::from_header(
                                &self.header,
                                cell_tag,
                                key_types,
                            );
                            self.pipeline.append_step(
                                FindBoundariesAndSort::new(
                                    TemplateStrategy::new(acc),
                                    sort_threads,
                                    byte_limit,
                                ),
                                t,
                            )
                        }
                        SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
                            // Queryname (lexicographic): the read name IS the key,
                            // embedded in the record; the strategy comparator-sorts
                            // arena refs into an InMemoryChunk<RawQuerynameLexKey>.
                            self.pipeline.append_step(
                                FindBoundariesAndSort::new(
                                    QuerynameStrategy::<RawQuerynameLexKey>::new(
                                        MemoryChunkErased::QuerynameLex,
                                    ),
                                    sort_threads,
                                    byte_limit,
                                ),
                                t,
                            )
                        }
                        SortOrder::Queryname(QuerynameComparator::Natural) => {
                            self.pipeline.append_step(
                                FindBoundariesAndSort::new(
                                    QuerynameStrategy::<RawQuerynameKey>::new(
                                        MemoryChunkErased::QuerynameNatural,
                                    ),
                                    sort_threads,
                                    byte_limit,
                                ),
                                t,
                            )
                        }
                        SortOrder::Coordinate => self.pipeline.append_step(
                            FindBoundariesAndSort::new(
                                CoordinateStrategy::new(n_ref),
                                sort_threads,
                                byte_limit,
                            ),
                            t,
                        ),
                    };
                    // Opt this chain into downstream-first dispatch: the serial
                    // FindBoundariesAndSort scan then overlaps the 8-way parallel
                    // InflateToArena instead of starving behind it on the pool.
                    self.use_drain_first_scheduler = true;
                    out
                } else {
                    let sort_buffer = SortBuffer::from_sorter(sorter, &self.header, byte_limit)
                        .map_err(|e| anyhow!("SortBuffer::from_sorter: {e}"))?
                        .with_affinity(affinity);
                    self.pipeline.append_step(sort_buffer, tail)
                };
                // `SpillGather` runs on the off-pool coordination driver
                // (`StepKind::Detached`, group `sort-coord`) alongside the sort
                // head, so it never contends a pool worker — the driver frames the
                // just-sorted chunk into blocks while the pool inflates/compresses.
                let serialize = SpillGather::new(byte_limit);
                let compress = SpillCompress::new(temp_codec, temp_compression, byte_limit);
                let write = SpillWrite::new(alloc, temp_codec, byte_limit, temp_dirs);
                // Lever 2, Phase-1 analogue: on the standalone-sort terminal,
                // detach the spill writer onto its own dedicated thread as well,
                // so the single serial spill-write stream stops occupying a pool
                // worker that could be running the compression-bound
                // `SpillCompress`. One persistent writer thread for the whole run
                // (not one per spill chunk). Same gate as the terminal output
                // writer set above; every other chain keeps the pool-scheduled
                // `Serial + Affinity::Writer` spill writer.
                let write = if self.detached_writer { write.with_detached() } else { write };
                let tail = self.pipeline.append_step(serialize, head_tail);
                let tail = self.pipeline.append_step(compress, tail);
                self.pipeline.append_step(write, tail)
            };
            let decompress_tail = self.pipeline.append_step(decompress, phase1_tail);

            if position == StagePosition::Terminal {
                // Terminal sort: `SortMerge<BlockOutput>` frames merged records
                // directly into `DecompressedBlock`s (lever 1), wired straight to
                // `BgzfCompress → WriteBgzfFile` in add_sink. No
                // `SerializeRecordBatch` step (and so no extra reorder stage) —
                // the framing is byte-for-byte identical to what that step
                // produced, but emitted on the already-in-order `Detached` merge
                // thread, dropping one pool step and one memcpy per record.
                let mut merge =
                    SortMerge::<BlockOutput>::new(sort_order, self.tuning.per_step_byte_limit);
                if let Some(slot) = &sort_stats_slot {
                    merge = merge.with_stats_slot(Arc::clone(slot));
                }
                let merge_tail = self.pipeline.append_step(merge, decompress_tail);
                self.current_tail = Some(merge_tail);
                // tail is DecompressedBlock (serialized bytes) directly from
                // SortMerge → SerializedBytes.
                self.chain_tail_kind = ChainTailKind::SerializedBytes;

                // Standalone sort: log the `=== Summary ===` block from the
                // SortMerge stats slot. Success-gated (`finalize_on_success`), NOT
                // always-run: the always-run `finalize` hooks drain even on the
                // error path, and the stats slot is unset (default zeros) whenever
                // SortMerge never reached `Done`, so an always-run summary would
                // report a bogus "Records processed: 0 ... completed" line on a
                // failed run. Registered before the index hook so on success the
                // summary still lands before any "Indexing BAM:" line, mirroring
                // the legacy SortFinalizeHook ordering.
                if let Some(slot) = sort_stats_slot {
                    let output_path = match &self.spec.sink {
                        SinkSpec::Bam(p) | SinkSpec::BamWithIndex(p) => p.clone(),
                    };
                    self.finalize_on_success.push(Box::new(
                        crate::pipeline::chains::commands::sort::SortSummaryFinalizeHook {
                            stats_slot: slot,
                            output_path,
                            timer: fgumi_cli_common::OperationTimer::new("Sorting BAM"),
                        },
                    ));
                }

                // When the spec asks for a sidecar BAI, queue the
                // `IndexBamFinalizeHook` on the *success-gated* finalize list
                // so the indexer runs after the chain has flushed the final
                // BAM on a successful run only (a failed run leaves a partial
                // BAM whose index would be stale). Every Sort-terminal chain —
                // standalone `[Stage::Sort]` or fused (e.g.
                // `[Stage::Correct, Stage::Sort]`) — lands here in the unified
                // `add_sort` flow, so without this wire a `BamWithIndex` request
                // would pass Rule 3 and silently skip indexing. Today only
                // standalone `fgumi sort` builds such specs (runall never
                // constructs `SinkSpec::BamWithIndex` because its sort output is
                // template-coordinate, where BAI is undefined), but the
                // architectural invariant — "BamWithIndex triggers the hook in
                // every terminal-sort path" — is enforced here for future
                // producers.
                if let SinkSpec::BamWithIndex(p) = &self.spec.sink {
                    self.finalize_on_success
                        .push(Box::new(IndexBamFinalizeHook { output_path: p.clone() }));
                }
            } else {
                // Intermediate: `SortMerge<RecordBatchOutput>` (the default)
                // emits `RecordBatch`, decoded to `DecodedRecordBatch` for the
                // downstream stages. Route through `bam_group_key_config` so a
                // `[Sort, Group]` chain still caches the UMI (RX) value position
                // during this decode — otherwise group would re-scan aux per
                // template and lose the #334/#343 optimization on the sort→group
                // path. The intermediate sort is never standalone, so the stats
                // slot is `None`; no `with_stats_slot` call is needed.
                let merge = SortMerge::<RecordBatchOutput>::new(
                    sort_order,
                    self.tuning.per_step_byte_limit,
                );
                let merge_tail = self.pipeline.append_step(merge, decompress_tail);
                let group_key_config = self.bam_group_key_config();
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
    /// appended; the chain tail stays as [`BatchedProcessedPositionGroups`] so
    /// the consensus stages (`add_simplex`/`add_duplex`/`add_codec`) can prepend
    /// `templates_to_mi_step` and consume it. This is the fused group→consensus
    /// path.
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
    /// Returns errors if the sort order is wrong or if the group options are
    /// missing from the spec bag.
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

        // Reject the same incompatible flag combinations that standalone
        // `fgumi group` rejects up front, so the runall path enforces them too
        // (S5c2-003). These read the user-supplied `strategy`/`no_umi`/
        // `min_umi_length` (not the derived `effective_*`), matching
        // `GroupReadsByUmi::execute`.
        if group.min_umi_length.is_some()
            && matches!(group.strategy, crate::assigner::Strategy::Paired)
        {
            bail!("Paired strategy cannot be used with --min-umi-length");
        }
        if group.no_umi && matches!(group.strategy, crate::assigner::Strategy::Paired) {
            bail!("--no-umi cannot be used with --strategy paired");
        }

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
            group.parallel_group_min_templates.clone(),
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
            // tail is DecompressedBlock (serialized bytes) after SerializeGroups
            // → SerializedBytes.
            self.chain_tail_kind = ChainTailKind::SerializedBytes;
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

    // ----------------------------------------------------------------------
    // Consensus stages (simplex / duplex / codec).
    //
    // The full method bodies below are compiled only with the `consensus`
    // feature, which gates the consensus command modules and their option
    // types. When the crate is built without `consensus`, the consensus CLI
    // commands and option-bag slots are compiled out, so these stages are
    // unreachable — but `add_stage`'s match must still resolve all three
    // `Stage` variants. These stubs provide that resolution and fail loudly
    // if ever reached.
    // ----------------------------------------------------------------------
    // `&mut self` matches the consensus-enabled signatures so the call sites are
    // identical; the stub bails before touching any field, hence `unused_self`.
    #[cfg(not(feature = "consensus"))]
    #[allow(clippy::unused_self)]
    fn add_simplex(&mut self, _position: StagePosition) -> Result<()> {
        bail!("fgumi was built without consensus support; rebuild with `--features consensus`")
    }

    #[cfg(not(feature = "consensus"))]
    #[allow(clippy::unused_self)]
    fn add_duplex(&mut self, _position: StagePosition) -> Result<()> {
        bail!("fgumi was built without consensus support; rebuild with `--features consensus`")
    }

    #[cfg(not(feature = "consensus"))]
    #[allow(clippy::unused_self)]
    fn add_codec(&mut self, _position: StagePosition) -> Result<()> {
        bail!("fgumi was built without consensus support; rebuild with `--features consensus`")
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
    #[cfg(feature = "consensus")]
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
        //   (from add_group(Intermediate)) → prepend the `templates_to_mi_step`
        //   bridge. Simplex passes `duplex_strip_strand = false` (no /A /B
        //   suffixes).
        let tail = if self.chain_tail_kind == ChainTailKind::BatchedProcessedPositionGroups {
            // simplex: no `/A` `/B` strand stripping, no record filter.
            self.pipeline.append_step(
                templates_to_mi_step(
                    self.tuning.per_step_byte_limit,
                    false,
                    None::<fn(&[u8]) -> bool>,
                ),
                tail,
            )
        } else {
            // Normal path: DecodedRecordBatch → GroupByMi → BatchedMiGroups.
            // The cell-tag partition is defensive/redundant given MI uniqueness
            // (see `templates_to_mi_step` doc, S5b1-003): MI is globally unique
            // across cell barcodes, so MI-only grouping is already cell-correct.
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
            let step = build_simplex_consensus_step_with_rejects(limit, consensus_cap);
            let pt = self.pipeline.append_step(step, tail);
            self.wire_consensus_rejects_branch(
                pt,
                simplex.rejects_opts.rejects.as_deref(),
                &input_header,
                "simplex",
            )?;
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
    #[cfg(feature = "consensus")]
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
            // duplex: strip `/A` `/B` so both strands group together, and apply
            // the same per-record filter the non-fused path applies (NEW-002).
            self.pipeline.append_step(
                templates_to_mi_step(
                    self.tuning.per_step_byte_limit,
                    true,
                    Some(duplex_record_filter),
                ),
                tail,
            )
        } else {
            // Normal path: DecodedRecordBatch → GroupByMi → BatchedMiGroups.
            // Record filter: skip secondary/supplementary; require mapped or
            // mapped-mate (shared with the fused bridge above).
            // MI transform: strip /A and /B so both strands group together.
            let mi_transform = |mi_bytes: &[u8]| -> String {
                let mi_str = std::borrow::Cow::from(std::str::from_utf8(mi_bytes).unwrap_or(""));
                extract_mi_base(&mi_str).to_string()
            };
            // `with_cell_tag` is defensive/redundant given MI uniqueness
            // (S5b1-003); `with_record_filter` is the load-bearing duplex filter.
            let group_mi_step = GroupByMi::new(*SamTag::MI, self.tuning.per_step_byte_limit)
                .with_cell_tag(Some(*SamTag::CB))
                .with_record_filter(duplex_record_filter)
                .with_mi_transform(mi_transform);
            self.pipeline.append_step(group_mi_step, tail)
        };

        // Wire the duplex consensus step. With `--rejects` it is a 2-output step
        // (branch 0 = consensus, branch 1 = rejects → BgzfCompress → WriteBgzfFile
        // with the INPUT header, in input order — the PR #332 fan-out contract);
        // without rejects it is the 1-output kept-only variant. Mirrors add_correct.
        let limit = self.tuning.per_step_byte_limit;
        let consensus_branch0 = if track_rejects {
            let step = build_duplex_consensus_step_with_rejects(limit, consensus_cap);
            let pt = self.pipeline.append_step(step, tail);
            self.wire_consensus_rejects_branch(
                pt,
                duplex.rejects_opts.rejects.as_deref(),
                &input_header,
                "duplex",
            )?;
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
    #[cfg(feature = "consensus")]
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
            // codec: no `/A` `/B` strand stripping, no record filter.
            self.pipeline.append_step(
                templates_to_mi_step(
                    self.tuning.per_step_byte_limit,
                    false,
                    None::<fn(&[u8]) -> bool>,
                ),
                tail,
            )
        } else {
            // Normal path: DecodedRecordBatch → GroupByMi → BatchedMiGroups.
            // Codec: no record filter (processes all mapped reads), no MI
            // transform. The cell-tag partition is defensive/redundant given
            // MI uniqueness (S5b1-003).
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
            let step = build_codec_consensus_step_with_rejects(limit, consensus_cap);
            let pt = self.pipeline.append_step(step, tail);
            self.wire_consensus_rejects_branch(
                pt,
                codec.rejects_opts.rejects.as_deref(),
                &input_header,
                "codec",
            )?;
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

/// Resolve a sort phase's worker-thread count.
///
/// A per-phase override (`--sort::sort-threads` / `--sort::merge-threads`, or
/// the standalone `--sort-threads` / `--merge-threads`) is used as-is when
/// present; otherwise the phase falls back to the chain's base sorter-thread
/// count (`--threads`). The result is clamped to at least 1 so a `0` override
/// never disables a phase. This is the single source of truth for both the
/// sole-stage (`SortStepCaptures`) and streaming (`RawExternalSorter`) sort
/// paths in `add_sort`.
fn resolve_phase_threads(override_threads: Option<usize>, num_sorter_threads: usize) -> usize {
    override_threads.unwrap_or(num_sorter_threads).max(1)
}

#[cfg(test)]
mod tests {
    use super::resolve_phase_threads;

    /// Pin the per-phase thread resolution contract shared by the standalone and
    /// streaming sort branches in `add_sort`: each phase falls back to the base
    /// sorter-thread count when no explicit override is supplied, and is clamped
    /// to at least 1. Exercises the production `resolve_phase_threads` helper
    /// directly so it cannot drift from the code `add_sort` actually runs.
    #[test]
    fn sole_stage_sort_resolves_per_phase_threads() {
        let num_sorter_threads = 8usize;

        // Explicit overrides are used as-is (clamped ≥ 1).
        assert_eq!(resolve_phase_threads(Some(1), num_sorter_threads), 1);
        assert_eq!(resolve_phase_threads(Some(2), num_sorter_threads), 2);
        // Zero is clamped to 1.
        assert_eq!(resolve_phase_threads(Some(0), num_sorter_threads), 1);
        // None falls back to num_sorter_threads.
        assert_eq!(resolve_phase_threads(None, num_sorter_threads), 8);
    }
}
