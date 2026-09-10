//! Inline consensus-metrics accumulator: the per-worker state fed by the
//! fused (T1) and standalone (T2) adapters to compute the same CS/SS/DS
//! family-size, UMI, and downsampling-yield metrics the separate-pass
//! simplex-metrics/duplex-metrics commands compute, inline, during the
//! existing streaming pass. See
//! docs/superpowers/specs/2026-09-08-inline-consensus-metrics-design.md §3-§6.
//!
//! Every item here is `pub(crate)`, consumed by the fused (T1) accumulator in
//! `add_group` for all three modes and by the standalone (T2) simplex,
//! duplex, and codec producers (`run_{simplex,duplex,codec}_consensus_batch_with_metrics`
//! in `pipeline/chains/commands/{simplex,duplex,codec}.rs`). Both paths now
//! record metrics fully in parallel, sharded across worker threads via
//! `PerThreadAccumulator<ConsensusMetricsSlot>`: T1's fused tap calls
//! `record_coordinate_group` directly on each worker, since its coordinate
//! groups are always complete and never split across a batch boundary; T2's
//! producers split each batch's entries into same-key runs
//! (`split_into_runs`/`classify_batch_runs`), record interior (fully-within-
//! one-batch) runs directly, and defer the first/last run of each batch as a
//! `BoundaryRun` for `reassemble_boundary` to fold back into stream order at
//! finalize. The earlier `CoordinateGroupFragment`/`CoordinateGroupCollector`/
//! `MetricsCollectorStep` serial-collector design — a third output branch
//! reordered and reduced, single-threaded, by a dedicated pipeline step after
//! the fact — has been fully retired: every mode now records metrics inline
//! in the consensus worker body, in parallel with consensus calling itself.
//!
//! **Fix round 1 (post-Task-7 review):** the original per-slot design called
//! `record_simplex_coordinate_group`/`record_duplex_coordinate_group` once
//! *per fraction slot*, each time with a 1-element `fractions` window. Both
//! reducers gate their once-per-group work (UMI decoding for Simplex;
//! `duplex_umi_counts`'s UMI-string building for Duplex) on
//! `idx == fractions.len() - 1` — with a 1-element window that's always
//! `0 == 0`, so the gate fired on *every* fraction instead of only the 100%
//! one (spec §5.6). Calling per-slot also recomputed
//! `compute_template_metadata(group)` 20x per coordinate group instead of
//! once. The fix restructures this type from a per-slot `[InlineCollector; 20]`
//! array into two homogeneous-array variants (`Simplex`/`Duplex`) so each
//! reducer is called exactly ONCE per group with the full 20-length
//! `DOWNSAMPLING_FRACTIONS`/`collectors`/`fraction_template_counts` arrays —
//! matching how the separate-pass `simplex_metrics`/`duplex_metrics` commands
//! call them (`simplex_metrics.rs:131-147`, `duplex_metrics.rs:176`).
//!
//! **Parallel T2 consensus metrics primitives:** `RunKind`, `BoundaryRun`,
//! `ConsensusMetricsSlot`, `split_into_runs`, `classify_batch_runs`, and
//! `reassemble_boundary` are pure, fully-unit-tested order-free building
//! blocks for the parallel-worker T2 wiring. `split_into_runs` and
//! `classify_batch_runs` are called directly by the standalone simplex,
//! duplex, and codec producers' per-batch bodies; `reassemble_boundary` is
//! called by `ConsensusMetricsFinalizeHook::finalize` for every mode.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::commands::group::with_extension;
use crate::commands::shared_metrics::{
    DOWNSAMPLING_FRACTIONS, Interval, ReadInfoKey, TemplateInfo, build_template_info,
    build_template_info_with_mi, overlaps_intervals, record_duplex_coordinate_group,
    record_simplex_coordinate_group,
};
use crate::mi_group::MiGroup;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::FinalizeHook;
use crate::simple_umi_consensus::SimpleUmiConsensusCaller;
use crate::template::Template;
use anyhow::Result;
use fgumi_bam_io::LibraryIndex;
use fgumi_metrics::duplex::DuplexMetricsCollector;
use fgumi_metrics::simplex::SimplexMetricsCollector;
use fgumi_raw_bam::{RawRecord, flags as raw_flags};
use indexmap::IndexMap;

/// Accumulator state for one consensus stage's inline metrics. Shared by T1
/// and T2 alike: both wrap it (inside a [`ConsensusMetricsSlot`]) in a
/// `PerThreadAccumulator`, N sharded instances (one per worker) folded via
/// `merge()` in `ConsensusMetricsFinalizeHook::finalize`.
///
/// One variant per mode, each holding a homogeneous 20-element array (one
/// slot per `DOWNSAMPLING_FRACTIONS` entry) so the shared reducer functions
/// can be called once per coordinate group with the FULL array — required
/// for their internal `idx == fractions.len() - 1` ("only at 100%") gates to
/// see the true last index, and to amortize `compute_template_metadata`
/// across all 20 fractions instead of recomputing it per fraction.
///
/// `fraction_template_counts` mirrors the separate-pass commands'
/// `process_templates_from_bam` accumulator (`shared_metrics.rs:893`): a
/// running per-fraction template count across every coordinate group this
/// accumulator ever sees, not per-call scratch. It is persisted here (not
/// reallocated per `record_coordinate_group` call) so a future finalize step
/// can build yield metrics from it the same way the separate-pass commands
/// do (`collector.to_yield_metric(fraction, read_pairs, min_reads)` with
/// `read_pairs` sourced from this array).
pub(crate) enum ConsensusMetricsAccumulator {
    Simplex {
        // Boxed: `SimplexMetricsCollector` and `DuplexMetricsCollector` are
        // different sizes, so an unboxed `[Collector; 20]` field on each
        // variant makes this enum's stack footprint the size of its LARGEST
        // variant (clippy::large_enum_variant) — every `Simplex` value would
        // carry the unused space for a `Duplex`-sized payload. Boxing shrinks
        // each variant to a pointer.
        collectors: Box<[SimplexMetricsCollector; 20]>,
        umi_caller: SimpleUmiConsensusCaller,
        fraction_template_counts: [usize; 20],
    },
    Duplex {
        collectors: Box<[DuplexMetricsCollector; 20]>,
        umi_caller: SimpleUmiConsensusCaller,
        fraction_template_counts: [usize; 20],
        /// Mirrors `duplex-metrics --duplex-umi-counts` (`self.duplex_umi_counts`
        /// on `DuplexMetrics`, `duplex_metrics.rs:104-105`), threaded into
        /// `record_duplex_coordinate_group`'s `duplex_umi_counts` parameter
        /// (Task 5).
        duplex_umi_counts: bool,
    },
}

impl ConsensusMetricsAccumulator {
    pub(crate) fn new_simplex() -> Self {
        Self::Simplex {
            collectors: Box::new(std::array::from_fn(|_| SimplexMetricsCollector::new())),
            umi_caller: SimpleUmiConsensusCaller::default(),
            fraction_template_counts: [0; 20],
        }
    }

    pub(crate) fn new_duplex(collect_duplex_umi_counts: bool) -> Self {
        Self::Duplex {
            collectors: Box::new(std::array::from_fn(|_| {
                DuplexMetricsCollector::new(collect_duplex_umi_counts)
            })),
            umi_caller: SimpleUmiConsensusCaller::default(),
            fraction_template_counts: [0; 20],
            duplex_umi_counts: collect_duplex_umi_counts,
        }
    }

    /// Records one already-assembled coordinate/strand group, applying the
    /// per-template interval filter itself (matching
    /// `process_templates_from_bam`'s ordering: interval filter strictly
    /// before any fraction bucketing) — callers pass the RAW group here.
    ///
    /// Calls the shared reducer exactly ONCE with the full 20-length
    /// `DOWNSAMPLING_FRACTIONS` array (see the type doc for why this matters).
    pub(crate) fn record_coordinate_group(
        &mut self,
        group: &[TemplateInfo],
        intervals: &[Interval],
    ) -> Result<()> {
        // In the common case (no `--intervals`) `overlaps_intervals` is an
        // unconditional pass-through, so borrow the caller's group directly
        // and only materialize a filtered `Vec` when an interval filter is
        // actually set — avoiding a full clone of every `TemplateInfo` on the
        // hot path.
        let owned_filtered;
        let filtered: &[TemplateInfo] = if intervals.is_empty() {
            group
        } else {
            owned_filtered = group
                .iter()
                .filter(|t| overlaps_intervals(t, intervals))
                .cloned()
                .collect::<Vec<_>>();
            &owned_filtered
        };
        if filtered.is_empty() {
            return Ok(());
        }
        match self {
            Self::Simplex { collectors, umi_caller, fraction_template_counts } => {
                record_simplex_coordinate_group(
                    filtered,
                    &DOWNSAMPLING_FRACTIONS,
                    &mut collectors[..],
                    umi_caller,
                    fraction_template_counts,
                )
            }
            Self::Duplex {
                collectors,
                umi_caller,
                fraction_template_counts,
                duplex_umi_counts,
            } => record_duplex_coordinate_group(
                filtered,
                &DOWNSAMPLING_FRACTIONS,
                &mut collectors[..],
                umi_caller,
                fraction_template_counts,
                *duplex_umi_counts,
            ),
        }
    }

    /// Folds `other` into `self` by zipping and merging the 20 per-fraction
    /// collectors (and per-fraction template counts) pairwise. Called by both
    /// T1 and T2's `ConsensusMetricsFinalizeHook::finalize` to fold every
    /// worker's per-thread accumulator into one before writing the metrics
    /// files.
    ///
    /// The `duplex_umi_counts` flag is left untouched — like
    /// `DuplexMetricsCollector::merge`'s `collect_duplex_umi_counts` field
    /// (Task 3), it is a per-stage config flag, not accumulated state; every
    /// accumulator for one stage is always constructed with the same value.
    pub(crate) fn merge(&mut self, other: Self) -> Result<()> {
        match (self, other) {
            (
                Self::Simplex { collectors: mine, fraction_template_counts: mine_counts, .. },
                Self::Simplex {
                    collectors: theirs, fraction_template_counts: theirs_counts, ..
                },
            ) => {
                for (m, t) in mine.iter_mut().zip(*theirs) {
                    m.merge(t);
                }
                for (m, t) in mine_counts.iter_mut().zip(theirs_counts) {
                    *m += t;
                }
            }
            (
                Self::Duplex { collectors: mine, fraction_template_counts: mine_counts, .. },
                Self::Duplex {
                    collectors: theirs, fraction_template_counts: theirs_counts, ..
                },
            ) => {
                for (m, t) in mine.iter_mut().zip(*theirs) {
                    m.merge(t);
                }
                for (m, t) in mine_counts.iter_mut().zip(theirs_counts) {
                    *m += t;
                }
            }
            _ => anyhow::bail!(
                "internal error: attempted to merge a Simplex accumulator with a Duplex \
                 one — every ConsensusMetricsAccumulator for one stage must be constructed \
                 via the same new_simplex()/new_duplex() call, never mixed"
            ),
        }
        Ok(())
    }
}

/// Bundles everything a consensus stage's inline-metrics accumulator needs
/// across its lifetime. Shared by both paths: `add_group` constructs one for
/// the fused (T1) case whenever a downstream consensus stage's `metrics`
/// field is set, so `add_simplex`/`add_duplex`/`add_codec` reuse the SAME
/// accumulator `Arc` the T1 tap closure already writes into; `add_simplex`/
/// `add_duplex`/`add_codec` construct their own (via
/// `build_consensus_metrics_captures`) for the standalone (T2) case, whose
/// per-batch producers record into it directly via
/// `split_into_runs`/`classify_batch_runs`.
pub(crate) struct ConsensusMetricsCaptures {
    pub(crate) accumulator: Arc<PerThreadAccumulator<ConsensusMetricsSlot>>,
    pub(crate) intervals: Vec<Interval>,
    pub(crate) output_prefix: PathBuf,
}

/// Simplex needs one min-reads threshold; duplex and codec need two
/// (AB/BA) — codec's caller has a single symmetric `min_reads_per_strand`
/// (verified: `crates/fgumi-consensus/src/codec_caller.rs:155-156`, checked
/// identically against both strands), so codec constructs this variant with
/// `min_ab_reads == min_ba_reads == codec's min_reads` (Task 11).
pub(crate) enum MetricsThresholds {
    Simplex { min_reads: usize },
    Duplex { min_ab_reads: usize, min_ba_reads: usize },
}

/// Finalize hook shared by T1 and T2 alike: merges every `PerThreadAccumulator`
/// slot's `ConsensusMetricsSlot` into one — folding each worker's
/// `ConsensusMetricsAccumulator` via `merge` and reassembling each worker's
/// deferred `BoundaryRun`s via [`reassemble_boundary`] — then writes the same
/// file set the separate-pass simplex-metrics/duplex-metrics commands write.
/// T1's fused tap always calls `record_coordinate_group` directly (its
/// coordinate groups are never split across a batch boundary), so its
/// `boundary` list is always empty and the reassembly step is a no-op for
/// that path — the same hook handles both cases uniformly.
pub(crate) struct ConsensusMetricsFinalizeHook {
    pub(crate) accumulators: Arc<PerThreadAccumulator<ConsensusMetricsSlot>>,
    pub(crate) output_prefix: PathBuf,
    pub(crate) intervals: Vec<Interval>,
    pub(crate) thresholds: MetricsThresholds,
}

impl FinalizeHook for ConsensusMetricsFinalizeHook {
    fn finalize(self: Box<Self>) -> anyhow::Result<()> {
        let ConsensusMetricsFinalizeHook { accumulators, output_prefix, intervals, thresholds } =
            *self;

        // `ConsensusMetricsSlot` has no mode-less `Default` (its inner
        // `ConsensusMetricsAccumulator` has none — Task 7), so draining the
        // sharded slots goes through `into_slots_with` (Task 6's
        // Default-free constructor's sibling) rather than `into_slots`. The
        // `init` closure only matters for the lossy fallback path
        // (outstanding `Arc` holders, a caller bug) — pick the constructor
        // matching this hook's own mode so the fallback's type checks out;
        // its *values* are never read since that path replaces, rather than
        // merges, an in-progress slot.
        let init: fn() -> ConsensusMetricsSlot = match &thresholds {
            MetricsThresholds::Simplex { .. } => ConsensusMetricsSlot::new_simplex,
            MetricsThresholds::Duplex { .. } => duplex_fallback_init,
        };

        let slots = accumulators.into_slots_with(init);
        let mut slots_iter = slots.into_iter();
        let Some(mut merged) = slots_iter.next() else {
            return Ok(()); // zero worker threads — nothing to merge
        };
        for slot in slots_iter {
            merged.merge(slot)?;
        }

        // Fold each worker's deferred boundary runs back into stream order
        // and close them into coordinate groups (`reassemble_boundary`,
        // Task 1), then record each reassembled group into the merged
        // accumulator. T1's fused path never defers any runs to the
        // boundary (its tap always calls `record_coordinate_group`
        // directly — see `group.rs`), so `boundary` is always empty here
        // and this loop is a no-op, keeping T1's behavior unchanged.
        let ConsensusMetricsSlot { mut acc, boundary } = merged;
        for group in reassemble_boundary(boundary) {
            acc.record_coordinate_group(&group, &intervals)?;
        }

        match thresholds {
            MetricsThresholds::Simplex { min_reads } => {
                write_simplex_metrics_files(&acc, &output_prefix, min_reads)
            }
            MetricsThresholds::Duplex { min_ab_reads, min_ba_reads } => {
                write_duplex_metrics_files(&acc, &output_prefix, min_ab_reads, min_ba_reads)
            }
        }
    }
}

/// Non-capturing fallback constructor for `ConsensusMetricsFinalizeHook`'s
/// `Duplex` `into_slots_with` `init` argument — see the SAFETY-style comment
/// at that call site for why the `duplex_umi_counts` value here (`false`)
/// is inconsequential.
fn duplex_fallback_init() -> ConsensusMetricsSlot {
    ConsensusMetricsSlot::new_duplex(false)
}

/// Writes the same three files the separate-pass `simplex-metrics` command
/// writes (`simplex_metrics.rs`'s `execute()` tail): `<prefix>.family_sizes.txt`,
/// `<prefix>.umi_counts.txt`, `<prefix>.simplex_yield_metrics.txt`. The 100%
/// fraction collector (the last slot) supplies the family-size and UMI
/// metrics; the yield curve is built across all 20 fraction slots, each paired
/// with its running per-fraction template count.
pub(crate) fn write_simplex_metrics_files(
    merged: &ConsensusMetricsAccumulator,
    output_prefix: &Path,
    min_reads: usize,
) -> anyhow::Result<()> {
    let ConsensusMetricsAccumulator::Simplex { collectors, fraction_template_counts, .. } = merged
    else {
        anyhow::bail!("internal error: Simplex thresholds passed to a Duplex-mode accumulator");
    };
    // The 100% fraction is the last slot — matching `simplex_metrics.rs`'s
    // `collectors.pop()` (the final entry is always the 1.0 fraction).
    let main_collector = &collectors[DOWNSAMPLING_FRACTIONS.len() - 1];
    let family_size_metrics = main_collector.family_size_metrics();
    let umi_metrics = main_collector.umi_metrics();
    let yield_metrics: Vec<_> = collectors
        .iter()
        .zip(DOWNSAMPLING_FRACTIONS.iter())
        .zip(fraction_template_counts.iter())
        .map(|((collector, &fraction), &read_pairs)| {
            collector.to_yield_metric(fraction, read_pairs, min_reads)
        })
        .collect();

    crate::metrics::writer::write_metrics_auto(
        with_extension(output_prefix, "family_sizes.txt"),
        &family_size_metrics,
    )?;
    crate::metrics::writer::write_metrics_auto(
        with_extension(output_prefix, "umi_counts.txt"),
        &umi_metrics,
    )?;
    crate::metrics::writer::write_metrics_auto(
        with_extension(output_prefix, "simplex_yield_metrics.txt"),
        &yield_metrics,
    )?;
    Ok(())
}

/// Writes the same file set the separate-pass `duplex-metrics` command writes
/// (`duplex_metrics.rs`'s `execute()` tail): `<prefix>.family_sizes.txt`,
/// `<prefix>.duplex_family_sizes.txt`, `<prefix>.umi_counts.txt`, optionally
/// `<prefix>.duplex_umi_counts.txt` (only when `duplex_umi_counts` is on), and
/// `<prefix>.duplex_yield_metrics.txt`.
pub(crate) fn write_duplex_metrics_files(
    merged: &ConsensusMetricsAccumulator,
    output_prefix: &Path,
    min_ab_reads: usize,
    min_ba_reads: usize,
) -> anyhow::Result<()> {
    let ConsensusMetricsAccumulator::Duplex {
        collectors,
        fraction_template_counts,
        duplex_umi_counts,
        ..
    } = merged
    else {
        anyhow::bail!("internal error: Duplex thresholds passed to a Simplex-mode accumulator");
    };
    let main_collector = &collectors[DOWNSAMPLING_FRACTIONS.len() - 1];
    let family_size_metrics = main_collector.family_size_metrics();
    let duplex_family_size_metrics = main_collector.duplex_family_size_metrics();
    let umi_metrics = main_collector.umi_metrics();
    let yield_metrics: Vec<_> = collectors
        .iter()
        .zip(DOWNSAMPLING_FRACTIONS.iter())
        .zip(fraction_template_counts.iter())
        .map(|((collector, &fraction), &read_pairs)| {
            collector.to_yield_metric(fraction, read_pairs, min_ab_reads, min_ba_reads)
        })
        .collect();

    crate::metrics::writer::write_metrics_auto(
        with_extension(output_prefix, "family_sizes.txt"),
        &family_size_metrics,
    )?;
    crate::metrics::writer::write_metrics_auto(
        with_extension(output_prefix, "duplex_family_sizes.txt"),
        &duplex_family_size_metrics,
    )?;
    crate::metrics::writer::write_metrics_auto(
        with_extension(output_prefix, "umi_counts.txt"),
        &umi_metrics,
    )?;
    if *duplex_umi_counts {
        let duplex_umi_metrics = main_collector.duplex_umi_metrics(&umi_metrics);
        crate::metrics::writer::write_metrics_auto(
            with_extension(output_prefix, "duplex_umi_counts.txt"),
            &duplex_umi_metrics,
        )?;
    }
    crate::metrics::writer::write_metrics_auto(
        with_extension(output_prefix, "duplex_yield_metrics.txt"),
        &yield_metrics,
    )?;
    Ok(())
}

/// Re-pairs one `MiGroup`'s flat record list into R1/R2 pairs by read name,
/// applying the same paired/mapped/primary filter
/// `process_templates_from_bam` uses. Uses an insertion-ordered `IndexMap`
/// keyed by read name, so the returned pairs come back in first-appearance
/// order of each read name within `records` (i.e. BAM order) rather than
/// hash-bucket order — matching the standalone `simplex-metrics`/
/// `duplex-metrics` oracle's within-family read order, which is required for
/// `umi_counts` parity at coordinate positions where UMI tags tie.
fn pair_records_by_read_name(records: &[RawRecord]) -> Vec<(&RawRecord, &RawRecord)> {
    let mut by_name: IndexMap<Vec<u8>, (Option<&RawRecord>, Option<&RawRecord>)> = IndexMap::new();
    for record in records {
        let flags = record.flags();
        let qualifies = (flags & raw_flags::PAIRED) != 0
            && (flags & raw_flags::UNMAPPED) == 0
            && (flags & raw_flags::MATE_UNMAPPED) == 0
            && (flags & raw_flags::SECONDARY) == 0
            && (flags & raw_flags::SUPPLEMENTARY) == 0;
        if !qualifies {
            continue;
        }
        let name = fgumi_raw_bam::read_name(record.as_ref()).to_vec();
        let entry = by_name.entry(name).or_default();
        if (flags & raw_flags::FIRST_SEGMENT) != 0 {
            entry.0 = Some(record);
        } else if (flags & raw_flags::LAST_SEGMENT) != 0 {
            entry.1 = Some(record);
        }
    }
    by_name
        .into_values()
        .filter_map(|(r1, r2)| match (r1, r2) {
            (Some(r1), Some(r2)) => Some((r1, r2)),
            _ => None,
        })
        .collect()
}

/// Appends one `MiGroup`'s entries onto a caller-owned, per-**batch**
/// `entries` accumulator. The per-batch worker body (see
/// `run_{simplex,duplex,codec}_consensus_batch_with_metrics`) calls this once
/// per `MiGroup` inside a batch's loop, into one shared `Vec`, then splits
/// that batch's entries into same-key runs via `split_into_runs` /
/// `classify_batch_runs`. Every paired, successfully-converted template in
/// the group becomes one
/// `(TemplateInfo, ReadInfoKey)` entry; a template that fails to produce one
/// (missing R1/R2, unmapped, no CIGAR) is silently omitted, matching
/// `process_templates_from_bam`'s own behavior.
///
/// # Errors
///
/// Returns an error if a qualifying R1/R2 pair is missing a required `MI`/`RX`
/// tag (propagated from [`build_template_info`]).
pub(crate) fn push_mi_group_entries(
    mi_group: &MiGroup,
    header: &noodles::sam::Header,
    library_index: &LibraryIndex,
    entries: &mut Vec<(TemplateInfo, ReadInfoKey)>,
) -> Result<()> {
    for (r1, r2) in pair_records_by_read_name(&mi_group.records) {
        if let Some((info, key)) = build_template_info(r1, r2, header, library_index)? {
            entries.push((info, key));
        }
    }
    Ok(())
}

/// Converts a fused position group's `Template`s into `TemplateInfo`s, the
/// shape `ConsensusMetricsAccumulator::record_coordinate_group` consumes.
/// This is the T1 (fused/runall) adapter — zero cross-batch risk, because
/// `chains/commands/group.rs`'s per-position closure only ever sees a whole,
/// already-complete position group (`GroupByPosition` only batches whole
/// `RawPositionGroup`s, `src/lib/pipeline/steps/group/position.rs`).
///
/// A template that fails to produce a `TemplateInfo` (missing R1 or R2, an
/// unmapped reference, or no CIGAR) is silently omitted, matching
/// `process_templates_from_bam`'s own behavior for the same cases.
///
/// Sources `mi` from the `Template`'s own `mi` FIELD (via
/// [`build_template_info_with_mi`]) rather than reading the `MI` aux tag —
/// at this point in the fused pipeline, `group` has already assigned
/// `template.mi`, but the tag is not written onto the record until BAM
/// serialization runs later, so reading it from the aux data here would
/// fail with "missing the required MI tag" on every fused
/// `runall --start-from group --consensus <mode> --<mode>::metrics=<prefix>`
/// run. `template.mi` is the LOCAL (pre-`MiAssignGroups`-offset) id, which is
/// correct here: the family partition within one coordinate group is
/// identical whether local or globally offset, and the local id still
/// carries any `/A`,`/B` duplex suffix.
///
/// # Errors
///
/// Returns an error if a qualifying R1/R2 pair is missing a required `RX`
/// tag (propagated from [`build_template_info_with_mi`]).
pub(crate) fn coordinate_group_from_processed_position(
    templates: &[Template],
    header: &noodles::sam::Header,
    library_index: &LibraryIndex,
) -> Result<Vec<TemplateInfo>> {
    let mut infos = Vec::with_capacity(templates.len());
    for template in templates {
        let (Some(r1), Some(r2)) = (template.r1(), template.r2()) else {
            continue;
        };
        if let Some((info, _key)) =
            build_template_info_with_mi(r1, r2, header, library_index, template.mi.to_string())?
        {
            infos.push(info);
        }
    }
    Ok(infos)
}

/// Which end (if either) of a batch's entry stream a boundary run sits on.
/// Declared in this order so `as u8` gives `Head=0 < Whole=1 < Tail=2` — the
/// ordering [`reassemble_boundary`] sorts on within one `batch_serial` so a
/// `Tail` run always sorts after a `Head`/`Whole` run from the same batch.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum RunKind {
    Head,
    Whole,
    Tail,
}

/// One batch's first or last (or, if the batch produced a single run, only)
/// run of contiguous same-`ReadInfoKey` entries, deferred for cross-batch
/// reassembly because — unlike an interior run — it might continue into an
/// adjacent batch's boundary run sharing the same key.
pub(crate) struct BoundaryRun {
    pub(crate) batch_serial: u64,
    pub(crate) kind: RunKind,
    pub(crate) key: ReadInfoKey,
    pub(crate) templates: Vec<TemplateInfo>,
}

/// Splits a batch's entries (coordinate-sorted stream order) into maximal
/// contiguous runs of equal `ReadInfoKey`, preserving order.
pub(crate) fn split_into_runs(
    entries: Vec<(TemplateInfo, ReadInfoKey)>,
) -> Vec<(ReadInfoKey, Vec<TemplateInfo>)> {
    let mut runs: Vec<(ReadInfoKey, Vec<TemplateInfo>)> = Vec::new();
    for (info, key) in entries {
        match runs.last_mut() {
            Some((k, templates)) if *k == key => templates.push(info),
            _ => runs.push((key, vec![info])),
        }
    }
    runs
}

/// Classifies a batch's runs: the first is Head, the last is Tail, a single run
/// is Whole (open on both ends); everything strictly between is interior (a
/// complete coordinate group). Returns (interior runs, boundary runs).
pub(crate) fn classify_batch_runs(
    batch_serial: u64,
    runs: Vec<(ReadInfoKey, Vec<TemplateInfo>)>,
) -> (Vec<(ReadInfoKey, Vec<TemplateInfo>)>, Vec<BoundaryRun>) {
    let n = runs.len();
    let mut interior = Vec::new();
    let mut boundary = Vec::new();
    for (idx, (key, templates)) in runs.into_iter().enumerate() {
        let kind = match idx {
            _ if n == 1 => Some(RunKind::Whole),
            0 => Some(RunKind::Head),
            i if i == n - 1 => Some(RunKind::Tail),
            _ => None,
        };
        match kind {
            Some(kind) => boundary.push(BoundaryRun { batch_serial, kind, key, templates }),
            None => interior.push((key, templates)),
        }
    }
    (interior, boundary)
}

/// Restores global stream order over the deferred boundary runs and closes
/// coordinate groups with the SAME rule the retired serial collector used:
/// a `Tail` always closes the open group and opens a new one; a
/// `Head`/`Whole` extends the open group iff its key matches, else flushes and
/// opens a new one. Returns the reassembled groups in order. Interior runs were
/// already recorded on the workers and are not present here.
pub(crate) fn reassemble_boundary(mut boundary: Vec<BoundaryRun>) -> Vec<Vec<TemplateInfo>> {
    boundary.sort_by_key(|r| (r.batch_serial, r.kind as u8));
    let mut groups: Vec<Vec<TemplateInfo>> = Vec::new();
    let mut open: Option<(ReadInfoKey, Vec<TemplateInfo>)> = None;
    for run in boundary {
        let extends =
            !matches!(run.kind, RunKind::Tail) && open.as_ref().is_some_and(|(k, _)| *k == run.key);
        if extends {
            open.as_mut().expect("checked Some").1.extend(run.templates);
        } else {
            if let Some((_, g)) = open.take() {
                groups.push(g);
            }
            open = Some((run.key, run.templates));
        }
    }
    if let Some((_, g)) = open.take() {
        groups.push(g);
    }
    groups
}

/// Per-worker slot for the parallel T2 consensus-metrics reduction: pairs an
/// order-free `ConsensusMetricsAccumulator` (fed by each worker's interior
/// runs, which need no cross-batch reassembly) with the worker's deferred
/// boundary runs (which do). `merge` folds two slots' accumulators and
/// concatenates their boundary runs — global ordering and group-closing is
/// deferred to [`reassemble_boundary`] at finalize time, not done here.
pub(crate) struct ConsensusMetricsSlot {
    pub(crate) acc: ConsensusMetricsAccumulator,
    pub(crate) boundary: Vec<BoundaryRun>,
}

impl ConsensusMetricsSlot {
    pub(crate) fn new_simplex() -> Self {
        Self { acc: ConsensusMetricsAccumulator::new_simplex(), boundary: Vec::new() }
    }

    pub(crate) fn new_duplex(collect_duplex_umi_counts: bool) -> Self {
        Self {
            acc: ConsensusMetricsAccumulator::new_duplex(collect_duplex_umi_counts),
            boundary: Vec::new(),
        }
    }

    pub(crate) fn merge(&mut self, mut other: Self) -> anyhow::Result<()> {
        self.acc.merge(other.acc)?;
        self.boundary.append(&mut other.boundary);
        Ok(())
    }
}

/// Shared `#[cfg(test)]` fixtures, used by `mod tests` (Task 7's
/// `ConsensusMetricsAccumulator` tests) and the parallel-T2 primitive test
/// modules below (`run_split_tests`, `reassemble_boundary_tests`,
/// `consensus_metrics_slot_tests`).
#[cfg(test)]
pub(crate) fn template(
    mi: &str,
    ref_name: &str,
    position: i32,
    hash_fraction: f64,
) -> TemplateInfo {
    TemplateInfo {
        mi: mi.to_string(),
        rx: "AAA".to_string(),
        ref_name: Some(ref_name.to_string()),
        position: Some(position),
        end_position: Some(position + 50),
        r1_positive: true,
        hash_fraction,
    }
}

/// Builds a minimal `ReadInfoKey` distinguished only by `ref_index`/`start1`
/// — the two fields the ordered-stream tests below vary — with every other
/// field pinned to a fixed, arbitrary value so two keys built from the same
/// `(ref_index, start)` always compare equal.
#[cfg(test)]
pub(crate) fn key(ref_index: usize, start: i32) -> ReadInfoKey {
    ReadInfoKey {
        ref_index1: ref_index,
        start1: start,
        strand1: false,
        ref_index2: ref_index,
        start2: start + 50,
        strand2: true,
        library: 0,
        cell_barcode: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn record_coordinate_group_drops_templates_outside_intervals() {
        let mut acc = ConsensusMetricsAccumulator::new_simplex();
        let group = vec![
            template("1", "chr1", 100, 1.0),
            template("2", "chr2", 100, 1.0), // different ref, outside the interval
        ];
        let intervals = vec![Interval { ref_name: "chr1".to_string(), start: 0, end: 1000 }];

        acc.record_coordinate_group(&group, &intervals).expect("records");

        let ConsensusMetricsAccumulator::Simplex { collectors, .. } = &acc else {
            panic!("expected Simplex variant");
        };
        let metrics = collectors[19].family_size_metrics();
        let cs_count_at_size_1: usize =
            metrics.iter().filter(|m| m.family_size == 1).map(|m| m.cs_count).sum();
        assert_eq!(cs_count_at_size_1, 1);
    }

    #[test]
    fn merge_folds_two_accumulators_fraction_collectors() {
        let mut a = ConsensusMetricsAccumulator::new_simplex();
        a.record_coordinate_group(&[template("1", "chr1", 100, 1.0)], &[]).expect("records");

        let mut b = ConsensusMetricsAccumulator::new_simplex();
        b.record_coordinate_group(&[template("2", "chr1", 200, 1.0)], &[]).expect("records");

        a.merge(b).expect("merge succeeds");

        let ConsensusMetricsAccumulator::Simplex { collectors, .. } = &a else {
            panic!("expected Simplex variant");
        };
        let metrics = collectors[19].family_size_metrics();
        let cs_count_at_size_1: usize =
            metrics.iter().filter(|m| m.family_size == 1).map(|m| m.cs_count).sum();
        assert_eq!(
            cs_count_at_size_1, 2,
            "both accumulators' size-1 families must be present after merge"
        );
    }

    /// Guards against the fix-round-1 regression: the original per-slot
    /// implementation called `record_simplex_coordinate_group` once per
    /// fraction slot with a 1-element `fractions` window, so the reducer's
    /// internal `idx == fractions.len() - 1` ("only build UMI metrics at the
    /// 100% fraction") gate saw `0 == 0` on every call and fired at every
    /// fraction, not just the 100% one. This test uses `hash_fraction = 0.05`
    /// so the group is included in EVERY downsampling fraction's subset
    /// (0.05 <= every fraction from 0.05 to 1.00) — under the old per-slot
    /// code, `collectors[0].umi_metrics()` would therefore be non-empty here,
    /// failing the first assertion below. With the fix (the reducer called
    /// once with the FULL 20-length fractions array), only `collectors[19]`
    /// (the 100% slot, where `idx == last_fraction_idx` is genuinely true)
    /// gets UMI metrics.
    #[test]
    fn record_coordinate_group_only_calls_record_umi_at_full_fraction() {
        let mut acc = ConsensusMetricsAccumulator::new_simplex();
        let group = vec![template("1", "chr1", 100, 0.05), template("1", "chr1", 100, 0.05)];
        acc.record_coordinate_group(&group, &[]).expect("records");

        let ConsensusMetricsAccumulator::Simplex { collectors, .. } = &acc else {
            panic!("expected Simplex variant");
        };

        assert!(
            collectors[0].umi_metrics().is_empty(),
            "sub-1.0 fraction (5%) must not record UMI metrics"
        );
        assert!(!collectors[19].umi_metrics().is_empty(), "100% fraction must record UMI metrics");
    }
}

#[cfg(test)]
mod consensus_metrics_finalize_hook_tests {
    use super::*;

    #[test]
    fn consensus_metrics_finalize_hook_merges_slots_and_writes() {
        let accumulators = PerThreadAccumulator::new_with(1, ConsensusMetricsSlot::new_simplex);
        accumulators.with_slot(|slot| {
            slot.acc.record_coordinate_group(&[template("0", "chr1", 100, 1.0)], &[]).unwrap();
        });

        let hook = ConsensusMetricsFinalizeHook {
            accumulators,
            output_prefix: std::env::temp_dir().join("consensus_metrics_finalize_hook_test"),
            intervals: Vec::new(),
            thresholds: MetricsThresholds::Simplex { min_reads: 1 },
        };
        Box::new(hook).finalize().expect("finalize succeeds");
    }
}

#[cfg(test)]
mod coordinate_group_from_processed_position_tests {
    use super::*;

    /// Builds a two-record `Template` (R1 + R2) from a `shared_metrics::tests::build_pair`
    /// pair, going through the same raw-record encode path production code uses
    /// (`encode_record_buf_to_raw`), then `Template::from_records` — the real
    /// constructor (`template.rs` has no `Builder` type).
    /// `mi` is a plain non-negative integer string (e.g. `"0"`, `"1"`) — set on
    /// both the record's `MI` aux tag (via `build_pair`, unused by production
    /// once the adapter sources `mi` from the field, but kept so a stray
    /// tag-read regression would still be caught) AND, since the fix, on the
    /// `Template.mi` FIELD as a `MoleculeId::Single` — matching what `group`
    /// has actually assigned by the time the fused tap runs.
    fn template_from_pair(name: &str, header: &noodles::sam::Header, mi: &str) -> Template {
        let (r1_buf, r2_buf) =
            crate::commands::shared_metrics::tests::build_pair(name, 0, 100, 0, 150, mi);
        let r1 = fgumi_raw_bam::encode_record_buf_to_raw(&r1_buf, header).expect("encode r1");
        let r2 = fgumi_raw_bam::encode_record_buf_to_raw(&r2_buf, header).expect("encode r2");
        let mut template = Template::from_records(vec![r1, r2]).expect("builds template");
        template.mi = fgumi_umi::MoleculeId::Single(mi.parse().expect("mi is a plain integer"));
        template
    }

    #[test]
    fn coordinate_group_from_processed_position_converts_every_paired_template() {
        let header = crate::commands::shared_metrics::tests::test_header();
        let template1 = template_from_pair("t1", &header, "0");
        let template2 = template_from_pair("t2", &header, "1");

        let library_index = LibraryIndex::from_header(&header);
        let infos = coordinate_group_from_processed_position(
            &[template1, template2],
            &header,
            &library_index,
        )
        .expect("converts");

        assert_eq!(infos.len(), 2);
        let mis: Vec<&str> = infos.iter().map(|info| info.mi.as_str()).collect();
        assert_eq!(mis, vec!["0", "1"], "each Template's MI tag must survive the conversion");
    }

    #[test]
    fn coordinate_group_from_processed_position_omits_a_template_missing_r2() {
        let header = crate::commands::shared_metrics::tests::test_header();
        let (r1_buf, _r2_buf) =
            crate::commands::shared_metrics::tests::build_pair("unpaired", 0, 100, 0, 150, "0");
        let r1 = fgumi_raw_bam::encode_record_buf_to_raw(&r1_buf, &header).expect("encode r1");
        let template = Template::from_records(vec![r1]).expect("builds R1-only template");

        let library_index = LibraryIndex::from_header(&header);
        let infos = coordinate_group_from_processed_position(&[template], &header, &library_index)
            .expect("converts");

        assert!(infos.is_empty(), "a template missing R2 must be silently omitted");
    }

    #[test]
    fn coordinate_group_from_processed_position_handles_an_empty_group() {
        let header = crate::commands::shared_metrics::tests::test_header();
        let library_index = LibraryIndex::from_header(&header);
        let infos = coordinate_group_from_processed_position(&[], &header, &library_index)
            .expect("converts");
        assert!(infos.is_empty());
    }

    /// Regression test for the verified T1 bug: at the fused inline-metrics
    /// tap in `chains/commands/group.rs`, `group` has already assigned the
    /// `Template.mi` FIELD, but the `MI` aux tag is not written onto the
    /// records until BAM serialization runs later. Build an R1/R2 pair
    /// carrying an `RX` tag but deliberately **no `MI` tag** (mirroring the
    /// fused-tap reality), set `Template.mi` directly, and assert the
    /// adapter still produces the correct `TemplateInfo.mi` — sourced from
    /// the field, not a tag read. Against the pre-fix
    /// `build_template_info` (which unconditionally read the `MI` aux tag)
    /// this failed with "missing the required MI tag".
    #[test]
    fn coordinate_group_from_processed_position_sources_mi_from_the_template_field_when_the_mi_tag_is_absent()
     {
        use crate::sam::SamTag;
        use fgumi_raw_bam::{SamBuilder as RawSamBuilder, testutil::encode_op};

        let header = crate::commands::shared_metrics::tests::test_header();
        let seq = vec![b'A'; 100];
        let quals = vec![30u8; 100];
        let cigar = encode_op(0, 100); // 100M

        let mut b1 = RawSamBuilder::new();
        b1.read_name(b"no-mi-tag")
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT | raw_flags::MATE_REVERSE)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[cigar])
            .sequence(&seq)
            .qualities(&quals)
            .mate_ref_id(0)
            .mate_pos(149);
        b1.add_string_tag(SamTag::RX, b"ACGT-TGCA");
        // Deliberately no MI tag — this is what distinguishes the fused-tap
        // reality from the existing `build_pair`-based tests above, whose
        // records always carry one.
        let r1_buf =
            fgumi_raw_bam::raw_record_to_record_buf(&b1.build(), &noodles::sam::Header::default())
                .expect("decode r1");

        let mut b2 = RawSamBuilder::new();
        b2.read_name(b"no-mi-tag")
            .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE)
            .ref_id(0)
            .pos(149)
            .mapq(60)
            .cigar_ops(&[cigar])
            .sequence(&seq)
            .qualities(&quals)
            .mate_ref_id(0)
            .mate_pos(99);
        b2.add_string_tag(SamTag::RX, b"ACGT-TGCA");
        let r2_buf =
            fgumi_raw_bam::raw_record_to_record_buf(&b2.build(), &noodles::sam::Header::default())
                .expect("decode r2");

        let r1 = fgumi_raw_bam::encode_record_buf_to_raw(&r1_buf, &header).expect("encode r1");
        let r2 = fgumi_raw_bam::encode_record_buf_to_raw(&r2_buf, &header).expect("encode r2");

        let mut template = Template::from_records(vec![r1, r2]).expect("builds template");
        // Mirrors what `group` has done by the time the fused tap runs: the
        // in-memory field is assigned, independent of the (not-yet-written)
        // aux tag.
        template.mi = fgumi_umi::MoleculeId::PairedA(7);

        let library_index = LibraryIndex::from_header(&header);
        let infos = coordinate_group_from_processed_position(&[template], &header, &library_index)
            .expect("must succeed by reading Template.mi, not error on the absent MI aux tag");

        assert_eq!(infos.len(), 1);
        assert_eq!(
            infos[0].mi, "7/A",
            "TemplateInfo.mi must come from the Template.mi field's Display format"
        );
    }
}

#[cfg(test)]
mod pair_and_push_tests {
    use super::*;

    /// Encode a mapped, paired, primary R1/R2 pair through the same raw-record
    /// path production code uses (`build_pair` → `encode_record_buf_to_raw`).
    fn raw_pair(name: &str, mi: &str, header: &noodles::sam::Header) -> (RawRecord, RawRecord) {
        let (r1_buf, r2_buf) =
            crate::commands::shared_metrics::tests::build_pair(name, 0, 100, 0, 150, mi);
        let r1 = fgumi_raw_bam::encode_record_buf_to_raw(&r1_buf, header).expect("encode r1");
        let r2 = fgumi_raw_bam::encode_record_buf_to_raw(&r2_buf, header).expect("encode r2");
        (r1, r2)
    }

    #[test]
    fn pair_records_by_read_name_pairs_two_complete_pairs() {
        let header = crate::commands::shared_metrics::tests::test_header();
        let (a1, a2) = raw_pair("a", "0", &header);
        let (b1, b2) = raw_pair("b", "1", &header);
        let records = [a1, a2, b1, b2];
        let pairs = pair_records_by_read_name(&records);
        assert_eq!(pairs.len(), 2, "two complete R1/R2 pairs must produce two tuples");
    }

    #[test]
    fn pair_records_by_read_name_drops_a_record_missing_its_mate() {
        let header = crate::commands::shared_metrics::tests::test_header();
        let (a1, _a2) = raw_pair("a", "0", &header);
        let records = [a1];
        let pairs = pair_records_by_read_name(&records);
        assert!(pairs.is_empty(), "an R1 with no matching R2 must be dropped");
    }

    #[test]
    fn pair_records_by_read_name_preserves_first_appearance_order() {
        let header = crate::commands::shared_metrics::tests::test_header();
        // R1s appear in order c, a, b; R2s interleaved.
        let (c1, c2) = raw_pair("c", "0", &header);
        let (a1, a2) = raw_pair("a", "1", &header);
        let (b1, b2) = raw_pair("b", "2", &header);
        let records = [c1, a1, b1, c2, a2, b2]; // all R1s before their R2s → first-appearance = c,a,b
        let pairs = pair_records_by_read_name(&records);
        let names: Vec<&[u8]> =
            pairs.iter().map(|(r1, _)| fgumi_raw_bam::read_name(r1.as_ref())).collect();
        assert_eq!(
            names,
            vec![&b"c"[..], &b"a"[..], &b"b"[..]],
            "pairs must come back in first-appearance order of R1, not HashMap order"
        );
    }

    #[test]
    fn pair_records_by_read_name_excludes_a_secondary_alignment() {
        use fgumi_raw_bam::{SamBuilder, testutil::encode_op};
        let cigar = encode_op(0, 100); // 100M
        let seq = vec![b'A'; 100];
        let quals = vec![30u8; 100];
        let mut builder = SamBuilder::new();
        builder
            .read_name(b"sec")
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT | raw_flags::SECONDARY)
            .ref_id(0)
            .pos(99)
            .mapq(60)
            .cigar_ops(&[cigar])
            .sequence(&seq)
            .qualities(&quals)
            .mate_ref_id(0)
            .mate_pos(149);
        let secondary = builder.build();
        let records = [secondary];
        let pairs = pair_records_by_read_name(&records);
        assert!(pairs.is_empty(), "a SECONDARY record must never qualify for pairing");
    }

    #[test]
    fn push_mi_group_entries_converts_every_qualifying_pair() {
        let header = crate::commands::shared_metrics::tests::test_header();
        let library_index = LibraryIndex::from_header(&header);
        let (a1, a2) = raw_pair("a", "0", &header);
        let (b1, b2) = raw_pair("b", "1", &header);
        let group = MiGroup::new("0".to_string(), vec![a1, a2, b1, b2]);

        let mut entries = Vec::new();
        push_mi_group_entries(&group, &header, &library_index, &mut entries).expect("pushes");

        assert_eq!(entries.len(), 2, "each complete pair becomes one entry");
        let mut mis: Vec<&str> = entries.iter().map(|(info, _key)| info.mi.as_str()).collect();
        mis.sort_unstable();
        assert_eq!(mis, vec!["0", "1"], "each pair's MI tag must survive the conversion");
    }

    #[test]
    fn push_mi_group_entries_appends_onto_a_non_empty_accumulator() {
        let header = crate::commands::shared_metrics::tests::test_header();
        let library_index = LibraryIndex::from_header(&header);
        let (a1, a2) = raw_pair("a", "0", &header);
        let group = MiGroup::new("0".to_string(), vec![a1, a2]);

        // Pre-seed with one unrelated entry to prove the function appends
        // rather than replacing the caller-owned accumulator (batch semantics).
        let mut entries = vec![(template("seed", "chr1", 10, 1.0), key(0, 10))];
        push_mi_group_entries(&group, &header, &library_index, &mut entries).expect("pushes");

        assert_eq!(entries.len(), 2, "the seed entry must be retained and the new pair appended");
    }
}

#[cfg(test)]
mod run_split_tests {
    use super::*;

    fn entry(mi: &str, ref_index: usize, start: i32) -> (TemplateInfo, ReadInfoKey) {
        (template(mi, "chr1", start, 1.0), key(ref_index, start))
    }

    #[test]
    fn split_empty_is_empty() {
        assert!(split_into_runs(vec![]).is_empty());
    }

    #[test]
    fn split_groups_contiguous_keys_and_splits_on_change() {
        // K1 K1 K2 K1 -> [K1 K1], [K2], [K1]
        let runs = split_into_runs(vec![
            entry("0", 0, 100),
            entry("1", 0, 100),
            entry("2", 1, 200),
            entry("3", 0, 100),
        ]);
        assert_eq!(runs.iter().map(|(_, t)| t.len()).collect::<Vec<_>>(), vec![2, 1, 1]);
        // `ReadInfoKey` has no `Debug` impl, so compare with `assert!` rather
        // than `assert_eq!` (which requires `Debug` for its failure message).
        assert!(runs[0].0 == key(0, 100));
        assert!(runs[1].0 == key(1, 200));
        assert!(runs[2].0 == key(0, 100));
    }

    #[test]
    fn classify_single_run_is_whole_no_interior() {
        let runs = split_into_runs(vec![entry("0", 0, 100), entry("1", 0, 100)]);
        let (interior, boundary) = classify_batch_runs(7, runs);
        assert!(interior.is_empty());
        assert_eq!(boundary.len(), 1);
        assert!(matches!(boundary[0].kind, RunKind::Whole));
        assert_eq!(boundary[0].batch_serial, 7);
    }

    #[test]
    fn classify_three_runs_head_interior_tail() {
        let runs = split_into_runs(vec![
            entry("0", 0, 100), // head (K@ref0,100)
            entry("1", 1, 200), // interior (K@ref1,200)
            entry("2", 2, 300), // tail (K@ref2,300)
        ]);
        let (interior, boundary) = classify_batch_runs(3, runs);
        assert_eq!(interior.len(), 1);
        // The middle entry ("1", ref_index=1, start=200) builds key(1, 200) via
        // the `entry` fixture above — not key(2, 200); the brief's snippet used
        // an illustrative placeholder here. `ReadInfoKey` has no `Debug` impl,
        // so compare with `assert!` rather than `assert_eq!`.
        assert!(interior[0].0 == key(1, 200)); // interior is the middle key
        assert_eq!(boundary.len(), 2);
        assert!(matches!(boundary[0].kind, RunKind::Head));
        assert!(matches!(boundary[1].kind, RunKind::Tail));
    }
}

#[cfg(test)]
mod reassemble_boundary_tests {
    use super::*;

    fn run(batch: u64, kind: RunKind, k: ReadInfoKey, mis: &[&str]) -> BoundaryRun {
        BoundaryRun {
            batch_serial: batch,
            kind,
            key: k,
            templates: mis.iter().map(|m| template(m, "chr1", 100, 1.0)).collect(),
        }
    }

    #[test]
    fn k1_k2_k1_in_one_batch_stays_two_groups() {
        // one batch [Head K1, (interior K2 recorded elsewhere), Tail K1]
        let groups = reassemble_boundary(vec![
            run(0, RunKind::Head, key(0, 100), &["a"]),
            run(0, RunKind::Tail, key(0, 100), &["b"]),
        ]);
        assert_eq!(
            groups.iter().map(Vec::len).collect::<Vec<_>>(),
            vec![1, 1],
            "Tail must open a new group; the two K1 pieces are distinct groups"
        );
    }

    #[test]
    fn genuine_cross_batch_continuation_merges() {
        // b0 tail K1 + b1 head K1 = ONE group
        let groups = reassemble_boundary(vec![
            run(0, RunKind::Tail, key(0, 100), &["a"]),
            run(1, RunKind::Head, key(0, 100), &["b"]),
        ]);
        assert_eq!(groups.iter().map(Vec::len).collect::<Vec<_>>(), vec![2]);
    }

    #[test]
    fn a_key_spanning_three_batches_is_one_group() {
        let groups = reassemble_boundary(vec![
            run(0, RunKind::Tail, key(0, 100), &["a"]),
            run(1, RunKind::Whole, key(0, 100), &["b", "c"]),
            run(2, RunKind::Head, key(0, 100), &["d"]),
        ]);
        assert_eq!(groups.iter().map(Vec::len).collect::<Vec<_>>(), vec![4]);
    }

    #[test]
    fn distinct_keys_close_separately_and_order_is_restored_from_shuffle() {
        // Deliberately out of order to prove the (batch_serial, kind) sort.
        let groups = reassemble_boundary(vec![
            run(1, RunKind::Head, key(1, 200), &["y"]),
            run(0, RunKind::Tail, key(0, 100), &["x"]),
        ]);
        assert_eq!(groups.iter().map(Vec::len).collect::<Vec<_>>(), vec![1, 1]);
    }

    /// Extracts each group's template MIs, in order, for content comparison
    /// (`TemplateInfo` has no `PartialEq`/`Debug` impl, so groups can't be
    /// `assert_eq!`'d directly — the `mi` field is a unique, human-readable
    /// stand-in for full template identity).
    fn group_mis(groups: &[Vec<TemplateInfo>]) -> Vec<Vec<&str>> {
        groups.iter().map(|g| g.iter().map(|t| t.mi.as_str()).collect()).collect()
    }

    /// Sharper complement to
    /// `distinct_keys_close_separately_and_order_is_restored_from_shuffle`:
    /// that test only asserts group *sizes* match across a shuffle of two
    /// non-merging runs. This test builds a single richer multiset of six
    /// `BoundaryRun`s spanning four batches — one immediate-close `Whole`, an
    /// unmatched `Head` that closes the prior group, a `Head` that DOES merge
    /// into a same-key `Tail` from an earlier batch, and a `Whole` that
    /// merges into a same-key `Tail` — in two completely different `Vec`
    /// orders (forward-built vs. reversed), and asserts `reassemble_boundary`
    /// returns byte-identical output (same group count, same per-group size,
    /// same per-group membership in the same order) from both, not merely
    /// same sizes. Because every run here has a distinct `(batch_serial,
    /// kind)` pair, the `(batch_serial, kind)` sort key fully determines
    /// processing order regardless of input order, which is exactly the
    /// determinism `reassemble_boundary` must guarantee against
    /// slot-merge/thread-scheduling-dependent collection order.
    #[test]
    fn two_full_orderings_of_the_same_boundary_runs_produce_identical_output() {
        let build = || {
            vec![
                run(0, RunKind::Whole, key(0, 100), &["a"]),
                run(1, RunKind::Head, key(1, 200), &["b1"]),
                run(1, RunKind::Tail, key(2, 300), &["c1"]),
                run(2, RunKind::Head, key(2, 300), &["c2"]),
                run(2, RunKind::Tail, key(3, 400), &["d1"]),
                run(3, RunKind::Whole, key(3, 400), &["d2"]),
            ]
        };

        let forward = build();
        let mut reversed = build();
        reversed.reverse();
        // Confirm the two input orderings are genuinely different, so a pass
        // below is not vacuous.
        assert_ne!(
            forward.iter().map(|r| (r.batch_serial, r.kind as u8)).collect::<Vec<_>>(),
            reversed.iter().map(|r| (r.batch_serial, r.kind as u8)).collect::<Vec<_>>(),
        );

        let groups_forward = reassemble_boundary(forward);
        let groups_reversed = reassemble_boundary(reversed);

        let expected = vec![vec!["a"], vec!["b1"], vec!["c1", "c2"], vec!["d1", "d2"]];
        assert_eq!(
            group_mis(&groups_forward),
            expected,
            "forward-ordered input must reassemble to the expected groups"
        );
        assert_eq!(
            group_mis(&groups_reversed),
            expected,
            "reversed-ordered input must reassemble to the SAME expected groups"
        );
        assert_eq!(
            group_mis(&groups_forward),
            group_mis(&groups_reversed),
            "reassembly must be independent of the collection order of the same boundary runs"
        );
    }
}

#[cfg(test)]
mod consensus_metrics_slot_tests {
    use super::*;

    fn size1_cs(acc: &ConsensusMetricsAccumulator) -> usize {
        let ConsensusMetricsAccumulator::Simplex { collectors, .. } = acc else { panic!() };
        collectors[19]
            .family_size_metrics()
            .iter()
            .filter(|m| m.family_size == 1)
            .map(|m| m.cs_count)
            .sum()
    }

    #[test]
    fn merge_folds_acc_and_concatenates_boundary() {
        let mut a = ConsensusMetricsSlot::new_simplex();
        a.acc.record_coordinate_group(&[template("0", "chr1", 100, 1.0)], &[]).unwrap();
        a.boundary.push(BoundaryRun {
            batch_serial: 0,
            kind: RunKind::Tail,
            key: key(0, 200),
            templates: vec![template("1", "chr1", 200, 1.0)],
        });
        let mut b = ConsensusMetricsSlot::new_simplex();
        b.acc.record_coordinate_group(&[template("2", "chr1", 300, 1.0)], &[]).unwrap();
        b.boundary.push(BoundaryRun {
            batch_serial: 1,
            kind: RunKind::Head,
            key: key(0, 200),
            templates: vec![template("3", "chr1", 200, 1.0)],
        });
        a.merge(b).unwrap();
        assert_eq!(size1_cs(&a.acc), 2, "both recorded size-1 families present after merge");
        assert_eq!(a.boundary.len(), 2, "boundary runs concatenated (ordered at finalize)");
    }
}
