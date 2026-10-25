//! Chain builder support for [`crate::pipeline::chains::Stage::Realign`].
//!
//! `add_realign` (on [`crate::pipeline::chains::builder::ChainBuilder`]) wires
//!
//! ```text
//! queryname grouper → RealignWindowProbe → <align backend> → RealignSelectPrimary
//! ```
//!
//! This module holds the two realign-specific steps and the finalize hook:
//!
//! - [`WindowProbeStep`] rejects an input that was never window-grouped (no
//!   `rw` on any primary of its first [`WINDOW_PROBE_TEMPLATES`] templates, or
//!   anywhere in a shorter input) before any of it reaches the aligner; see
//!   [`WindowProbe`]. Ends without `rw` in a window-grouped input are expected
//!   and pass. It also gives every lone mate a placeholder mate
//!   ([`add_lone_mate_placeholder`]), so the aligner is fed whole pairs.
//! - [`build_select_primary_step`] runs
//!   [`select_primary`] on every merged
//!   template and folds its outcomes into per-thread [`RealignCounts`].
//! - [`RealignFinalizeHook`] reduces those counts and logs them.

use std::collections::VecDeque;
use std::io;
use std::sync::Arc;

use anyhow::{Context, Result};
use log::{info, warn};
use noodles::sam::Header;

use fgumi_raw_bam::SamTag;

use crate::aligner::{BwaScoring, ResolvedAligner};
use crate::commands::realign::RealignOptions;
use crate::logging::OperationTimer;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::pipeline::chains::FinalizeHook;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::align::is_primary_for_alignment;
use crate::pipeline::steps::align::merge::add_lone_mate_placeholder;
use crate::pipeline::steps::process::{ProcessOrdered, process_ordered};
use crate::pipeline::steps::types::BamTemplateBatch;
use crate::realign::select::{RealignCounts, select_primary};
use crate::reference::ReferenceSource;
use crate::template::Template;

/// How many templates (in input order) [`WindowProbe`] reads before deciding the
/// input was never window-grouped: if none of the first this many templates has
/// a primary record carrying `rw`, realign fails.
pub(crate) const WINDOW_PROBE_TEMPLATES: u64 = 1000;

/// The error for an input whose leading templates carry no `rw` at all.
fn no_window_input_error(first_name: &[u8], templates: u64) -> io::Error {
    io::Error::other(format!(
        "realign: none of the first {templates} consensus templates carries an rw tag (the \
         molecule's raw-read window), starting with consensus read '{}'. realign needs \
         consensus reads called from a BAM grouped with `fgumi group --realign-window` (or \
         `fgumi runall --realign`); an input grouped without it (or by another tool) has \
         no windows to choose alignments with.",
        String::from_utf8_lossy(first_name),
    ))
}

/// Whether any primary record of `template` carries `rw`.
///
/// Secondary and supplementary records are ignored, matching the records the
/// align backends send to the aligner ([`is_primary_for_alignment`]).
fn has_window(template: &Template) -> bool {
    template.records.iter().any(|record| {
        is_primary_for_alignment(record.flags())
            && fgumi_raw_bam::find_string_tag(fgumi_raw_bam::aux_data_slice(record), SamTag::RW)
                .is_some()
    })
}

/// Decides, in input order, whether an input was ever window-grouped: rejects
/// an input that was not (an fgbio- or older fgumi-grouped consensus BAM),
/// without rejecting the ordinary consensus ends that legitimately lack `rw`.
///
/// `group --realign-window` stamps only mapped primaries, so a consensus end
/// whose raw reads never mapped (a mate-unmapped molecule, `--allow-unmapped`,
/// a duplex strand with no mapped reads) has no `rw`; such an end is aligned and
/// left as bwa placed it (counted as without a window). What is an error is an
/// input with no windows at all: if none of the first [`WINDOW_PROBE_TEMPLATES`]
/// templates in input order has a primary carrying `rw`, or a shorter input has
/// none anywhere, it fails naming the first read. An empty input passes.
///
/// [`WindowProbeStep`] feeds it batches in input order and withholds them until
/// it has passed, so an input it rejects never reaches the aligner or the output.
#[derive(Debug, Default)]
pub(crate) struct WindowProbe {
    /// Templates observed so far.
    templates: u64,
    /// The first template's name.
    first_name: Option<Vec<u8>>,
    /// A window was found within the first [`WINDOW_PROBE_TEMPLATES`] templates.
    passed: bool,
}

impl WindowProbe {
    /// Folds the next batch (in input order) into the probe; returns whether the
    /// input has passed (a windowed template among the first
    /// [`WINDOW_PROBE_TEMPLATES`]). Once passed, later batches are not scanned.
    ///
    /// # Errors
    ///
    /// Returns an error once [`WINDOW_PROBE_TEMPLATES`] templates have been
    /// observed without a window.
    pub(crate) fn observe(&mut self, templates: &[Template]) -> io::Result<bool> {
        if self.passed {
            return Ok(true);
        }
        if self.first_name.is_none() {
            self.first_name = templates.first().map(|t| t.name().to_vec());
        }
        let before = self.templates;
        self.templates += templates.len() as u64;
        if let Some(i) = templates.iter().position(has_window)
            && before + (i as u64) < WINDOW_PROBE_TEMPLATES
        {
            self.passed = true;
            return Ok(true);
        }
        if self.templates >= WINDOW_PROBE_TEMPLATES {
            let name = self.first_name.as_deref().unwrap_or_default();
            return Err(no_window_input_error(name, WINDOW_PROBE_TEMPLATES));
        }
        Ok(false)
    }

    /// The end-of-input decision for an input shorter than
    /// [`WINDOW_PROBE_TEMPLATES`] templates that has not passed: an empty input
    /// passes; a non-empty one had no windowed primary and fails.
    ///
    /// # Errors
    ///
    /// Returns an error naming the first read if the input had templates but no
    /// `rw` on any primary.
    pub(crate) fn finish(&mut self) -> io::Result<()> {
        if self.passed {
            return Ok(());
        }
        match &self.first_name {
            Some(name) if self.templates > 0 => Err(no_window_input_error(name, self.templates)),
            _ => {
                self.passed = true;
                Ok(())
            }
        }
    }
}

/// A `Serial` pass-through step that runs [`WindowProbe`] on the input in
/// order and holds every batch back until the probe has passed: an input that
/// was never window-grouped fails before any of it reaches the aligner or the
/// output. The input arrives in input order (its producer's branch is ordered),
/// so no reordering is needed here. At most [`WINDOW_PROBE_TEMPLATES`]
/// templates (plus the rest of the batch that reaches that count) are held.
/// It also gives every lone mate a placeholder mate ([`pair_lone_mates`]).
pub(crate) struct WindowProbeStep {
    probe: WindowProbe,
    /// Batches held back while the probe is undecided, then released in order.
    pending: VecDeque<BamTemplateBatch>,
    held: HeldSlot<Unpushed<BamTemplateBatch>>,
    limit_bytes: u64,
}

impl WindowProbeStep {
    /// A probe step whose output queue holds up to `limit_bytes`.
    pub(crate) fn new(limit_bytes: u64) -> Self {
        Self {
            probe: WindowProbe::default(),
            pending: VecDeque::new(),
            held: HeldSlot::new(),
            limit_bytes,
        }
    }

    /// Pushes `batch`, parking it in `held` if the output is full.
    fn push(&mut self, ctx: &mut StepCtx<'_, Self>, batch: BamTemplateBatch) {
        if let Err(unpushed) = ctx.outputs.push(batch) {
            self.held.put(unpushed);
        }
    }
}

impl Step for WindowProbeStep {
    type Input = BamTemplateBatch;
    type Outputs = OrderedBytesSingle<BamTemplateBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "RealignWindowProbe",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.limit_bytes }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take()
            && let Err(again) = ctx.outputs.retry(unpushed)
        {
            self.held.put(again);
            return Ok(StepOutcome::Contention);
        }
        // Once passed, release the held-back batches, oldest first.
        if self.probe.passed
            && let Some(batch) = self.pending.pop_front()
        {
            self.push(ctx, batch);
            return Ok(StepOutcome::Progress);
        }
        let Some(batch) = ctx.input.pop() else {
            if !ctx.input.is_drained() {
                return Ok(StepOutcome::NoProgress);
            }
            if !self.probe.passed {
                // A short input: decide now, before anything is released.
                self.probe.finish()?;
                return Ok(StepOutcome::Progress);
            }
            return Ok(if self.pending.is_empty() {
                StepOutcome::Finished
            } else {
                StepOutcome::NoProgress
            });
        };
        let batch = pair_lone_mates(batch)?;
        if self.probe.passed {
            self.push(ctx, batch);
        } else {
            self.probe.observe(batch.templates())?;
            self.pending.push_back(batch);
        }
        Ok(StepOutcome::Progress)
    }
}

/// `batch` with a placeholder mate added to every lone mate
/// ([`add_lone_mate_placeholder`]), so the aligner is fed only whole pairs and
/// bwa's even-parity `-K` cut cannot fall between the two ends of a real pair.
/// A batch without lone mates is returned as is.
fn pair_lone_mates(batch: BamTemplateBatch) -> io::Result<BamTemplateBatch> {
    if !batch.templates().iter().any(|t| t.lone_mate_segment().is_some()) {
        return Ok(batch);
    }
    let (serial, mut templates) = batch.into_parts();
    for template in &mut templates {
        add_lone_mate_placeholder(template)?;
    }
    // The placeholders add records, so the batch's heap size is recomputed.
    Ok(BamTemplateBatch::new(serial, templates))
}

/// The realign stage's aligner: `opts`' `--aligner::*` selection resolved as the
/// align stage resolves it, with bwa's `-h` applied only when the user set
/// `--xa-max-hits` / `--realign::xa-max-hits`.
///
/// # Errors
///
/// Returns an error if the aligner selection does not resolve.
pub(crate) fn resolve_realign_aligner(
    opts: &RealignOptions,
    num_threads: usize,
) -> Result<ResolvedAligner> {
    let align = &opts.align;
    let resolved = align
        .aligner
        .clone()
        .resolve(&align.reference, num_threads, align.aligner_bin.as_deref())
        .context("realign aligns the consensus reads with the --aligner::* selection")?;
    Ok(match opts.xa_max_hits {
        Some(xa) => resolved.with_xa_max_hits(xa.primary, xa.alt, opts.xa_max_hits_option),
        None => resolved,
    })
}

/// What [`build_select_primary_step`] needs, shared by every worker.
pub(crate) struct SelectPrimaryCaptures {
    /// The reference the promoted alignments are scored against.
    pub(crate) reference: ReferenceSource,
    /// The output header: names the contigs records refer to by id.
    pub(crate) header: Arc<Header>,
    /// The scoring `AS` is recomputed with for a promoted alignment.
    pub(crate) scoring: BwaScoring,
    /// Per-thread selection outcome counts.
    pub(crate) counts: Arc<PerThreadAccumulator<RealignCounts>>,
}

/// A `Parallel` step choosing each merged template's primary alignments with
/// [`select_primary`] (which also strips `rw`), keeping input order.
pub(crate) fn build_select_primary_step(
    limit_bytes: u64,
    captures: SelectPrimaryCaptures,
) -> ProcessOrdered<
    BamTemplateBatch,
    BamTemplateBatch,
    impl Fn(BamTemplateBatch) -> io::Result<BamTemplateBatch> + Send + Sync + 'static,
> {
    process_ordered("RealignSelectPrimary", limit_bytes, move |batch: BamTemplateBatch| {
        let SelectPrimaryCaptures { reference, header, scoring, counts } = &captures;
        let (serial, mut templates) = batch.into_parts();
        let mut outcomes = Vec::with_capacity(templates.len());
        for template in &mut templates {
            let selected = select_primary(template, reference, header, scoring).map_err(|e| {
                io::Error::other(format!(
                    "realign: selecting the primary of '{}': {e:#}",
                    String::from_utf8_lossy(template.name()),
                ))
            })?;
            outcomes.push(selected);
        }
        counts.with_slot(|c| outcomes.iter().for_each(|o| c.tally(o)));
        // Selection edits record bytes, so the batch's heap size is recomputed.
        Ok(BamTemplateBatch::new(serial, templates))
    })
}

/// Reduce the per-thread counts. The pipeline has drained before any finalize
/// hook runs, so every slot holds its final state.
fn reduce(counts: &PerThreadAccumulator<RealignCounts>) -> RealignCounts {
    let mut totals = RealignCounts::default();
    for slot in counts.slots() {
        totals.merge(&slot.lock());
    }
    totals
}

/// Success-only finalize hook: logs the realign summary and completion time.
pub(crate) struct RealignFinalizeHook {
    pub(crate) counts: Arc<PerThreadAccumulator<RealignCounts>>,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for RealignFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let RealignFinalizeHook { counts, timer } = *self;
        let totals = reduce(&counts);
        info!("=== Realign Summary ===");
        info!("Consensus reads realigned: {}", totals.reads);
        info!("Primary already in window: {}", totals.in_window);
        info!("Primary swapped to an in-window XA hit: {}", totals.swapped);
        info!(
            "Primary swapped to an in-window supplementary piece: {}",
            totals.swapped_supplementary
        );
        info!(
            "Primary kept outside window (proper pair with an in-window mate): {}",
            totals.kept_proper_pair
        );
        info!(
            "Primary kept outside window (proper pair whose mate cannot move into its window): {}",
            totals.kept_proper_pair_mate_stays
        );
        info!(
            "Primary outside window (no in-window XA hit or supplementary piece): {}",
            totals.outside_window
        );
        info!("Unmapped: {}", totals.unmapped);
        if totals.no_window > 0 {
            info!("Without a window: {}", totals.no_window);
        }
        if totals.malformed_xa > 0 {
            warn!("XA tags that did not parse (treated as absent): {}", totals.malformed_xa);
        }
        if totals.skipped_xa_hits > 0 {
            warn!("XA hits skipped for an invalid CIGAR: {}", totals.skipped_xa_hits);
        }
        timer.log_completion(totals.reads);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use fgumi_raw_bam::flags::{FIRST_SEGMENT, LAST_SEGMENT, PAIRED, SECONDARY, UNMAPPED};
    use fgumi_raw_bam::{RawRecord, SamBuilder};
    use rstest::rstest;

    use super::*;
    use crate::commands::zipper::merge_raw;
    use crate::umi::TagInfo;

    /// bwa's `-h` reaches the realign aligner of either backend only when set;
    /// unset, the command is passed through byte for byte.
    #[rstest]
    #[case::command_unset(None, None)]
    #[case::command_set(
        Some(crate::commands::realign::XaMaxHits { primary: 7, alt: 9 }),
        Some("bwa-mem3 mem -h 7,9 -p REF /dev/stdin")
    )]
    fn resolve_realign_aligner_applies_xa_max_hits_to_command(
        #[case] xa: Option<crate::commands::realign::XaMaxHits>,
        #[case] expected: Option<&str>,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let reference = tmp.path().join("ref.fa");
        let opts = realign_options(
            crate::aligner::AlignerOptions {
                command: Some("bwa-mem3 mem -p {ref} /dev/stdin".to_string()),
                ..crate::aligner::AlignerOptions::default()
            },
            &reference,
            xa,
        );
        let resolved = resolve_realign_aligner(&opts, 2).expect("resolve");
        let crate::aligner::ResolvedBackend::Subprocess { command, .. } = resolved.backend else {
            panic!("command mode resolves to the subprocess backend");
        };
        let unset = format!("bwa-mem3 mem -p {} /dev/stdin", reference.display());
        let want = expected.map_or(unset, |e| e.replace("REF", &reference.display().to_string()));
        assert_eq!(command, want);
    }

    /// The in-process backend records `-h` for its `MemOpts` only when set.
    #[cfg(feature = "aligner-bwa-mem3")]
    #[rstest]
    #[case::unset(None, None)]
    #[case::set(Some(crate::commands::realign::XaMaxHits { primary: 7, alt: 9 }), Some((7, 9)))]
    fn resolve_realign_aligner_applies_xa_max_hits_inproc(
        #[case] xa: Option<crate::commands::realign::XaMaxHits>,
        #[case] expected: Option<(u32, u32)>,
    ) {
        let tmp = tempfile::tempdir().unwrap();
        let reference = tmp.path().join("ref.fa");
        std::fs::write(&reference, ">chr1\nACGT\n").unwrap();
        for suffix in [".amb", ".ann", ".bwt.2bit.64", ".pac"] {
            std::fs::write(tmp.path().join(format!("ref.fa{suffix}")), b"").unwrap();
        }
        let opts = realign_options(
            crate::aligner::AlignerOptions {
                preset: Some(crate::aligner::AlignerPreset::BwaMem3InProc),
                ..crate::aligner::AlignerOptions::default()
            },
            &reference,
            xa,
        );
        let resolved = resolve_realign_aligner(&opts, 2).expect("resolve");
        let crate::aligner::ResolvedBackend::InProcessBwaMem3 { xa_max_hits, .. } =
            resolved.backend
        else {
            panic!("the in-process preset resolves to the in-process backend");
        };
        assert_eq!(xa_max_hits, expected);
    }

    fn realign_options(
        aligner: crate::aligner::AlignerOptions,
        reference: &std::path::Path,
        xa_max_hits: Option<crate::commands::realign::XaMaxHits>,
    ) -> RealignOptions {
        RealignOptions {
            align: crate::pipeline::chains::options_bag::AlignOptions {
                aligner,
                reference: reference.to_path_buf(),
                aligner_bin: None,
            },
            xa_max_hits,
            xa_max_hits_option: "--xa-max-hits",
            in_memory_reference: false,
        }
    }

    const RW: &str = "chr1:101-150:+";

    fn record(name: &str, flags: u16, rw: Option<&str>) -> RawRecord {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes()).flags(flags).sequence(b"ACGTACGT").qualities(&[30; 8]);
        if let Some(rw) = rw {
            b.add_string_tag(SamTag::RW, rw.as_bytes());
        }
        b.add_string_tag(SamTag::MI, b"7");
        b.build()
    }

    fn pair(rw1: Option<&str>, rw2: Option<&str>, extra: Option<RawRecord>) -> Template {
        let mut records = vec![
            record("q1", PAIRED | FIRST_SEGMENT | UNMAPPED, rw1),
            record("q1", PAIRED | LAST_SEGMENT | UNMAPPED, rw2),
        ];
        records.extend(extra);
        Template::from_records(records).expect("template")
    }

    /// `has_window` looks at primaries only: one windowed end is enough, and a
    /// secondary record's `rw` does not count.
    #[rstest]
    #[case::both_ends(pair(Some(RW), Some(RW), None), true)]
    #[case::r1_only(pair(Some(RW), None, None), true)]
    #[case::r2_only(pair(None, Some(RW), None), true)]
    #[case::neither(pair(None, None, None), false)]
    #[case::secondary_only(
        pair(None, None, Some(record("q1", PAIRED | FIRST_SEGMENT | SECONDARY, Some(RW)))),
        false
    )]
    fn has_window_checks_primaries(#[case] template: Template, #[case] expected: bool) {
        assert_eq!(has_window(&template), expected);
    }

    /// One batch for the probe: `n` templates named `b<serial>t<i>`, of which
    /// only the one at `window_at` (if any) has a windowed R1.
    fn batch(serial: u64, n: usize, window_at: Option<usize>) -> Vec<Template> {
        (0..n)
            .map(|i| {
                let name = format!("b{serial}t{i}");
                let rw = (window_at == Some(i)).then_some(RW);
                Template::from_records(vec![
                    record(&name, PAIRED | FIRST_SEGMENT | UNMAPPED, rw),
                    record(&name, PAIRED | LAST_SEGMENT | UNMAPPED, None),
                ])
                .expect("template")
            })
            .collect()
    }

    /// `WindowProbe` decides on the first `WINDOW_PROBE_TEMPLATES` (1000)
    /// templates in input order. `observe_err` is the index of the `observe`
    /// call expected to fail; `finish_ok` is the end-of-input check (run only
    /// when nothing failed); `named` is the read the error must name.
    #[rstest]
    #[case::window_first(vec![batch(0, 3, Some(0))], None, true, "")]
    #[case::window_late_in_small_input(vec![batch(0, 3, Some(2))], None, true, "")]
    #[case::no_window_small_input(vec![batch(0, 3, None), batch(1, 2, None)], None, false, "b0t0")]
    #[case::empty_input(vec![], None, true, "")]
    #[case::no_window_in_first_1000(vec![batch(0, 1000, None)], Some(0), false, "b0t0")]
    #[case::window_at_template_999(vec![batch(0, 999, None), batch(1, 5, Some(0))], None, true, "")]
    #[case::window_at_template_1000(
        vec![batch(0, 1000, None), batch(1, 5, Some(0))], Some(0), false, "b0t0"
    )]
    #[case::window_past_1000_in_one_batch(vec![batch(0, 1500, Some(1000))], Some(0), false, "b0t0")]
    #[case::settles_on_a_window(vec![batch(0, 2, Some(1)), batch(1, 2000, None)], None, true, "")]
    fn window_probe_decides_on_input_order(
        #[case] batches: Vec<Vec<Template>>,
        #[case] observe_err: Option<usize>,
        #[case] finish_ok: bool,
        #[case] named: &str,
    ) {
        let mut probe = WindowProbe::default();
        let mut failed = None;
        for (i, templates) in batches.iter().enumerate() {
            if let Err(e) = probe.observe(templates) {
                failed = Some((i, e));
                break;
            }
        }
        assert_eq!(failed.as_ref().map(|(i, _)| *i), observe_err, "{failed:?}");
        let error = if let Some((_, e)) = failed {
            Some(e)
        } else {
            let finished = probe.finish();
            assert_eq!(finished.is_ok(), finish_ok, "{finished:?}");
            finished.err()
        };
        if let Some(e) = error {
            let msg = e.to_string();
            assert!(msg.contains(&format!("'{named}'")), "error must name {named}: {msg}");
            assert!(msg.contains("--realign-window"), "error must say how to fix it: {msg}");
        }
    }

    /// Drives a [`WindowProbeStep`] over `batches` (in order, then end of
    /// input) and returns the serials it released before it finished or
    /// failed, and whether it failed.
    fn run_probe_step(batches: Vec<Vec<Template>>) -> (Vec<u64>, bool) {
        use crate::pipeline::core::Ordered;
        use crate::pipeline::core::step::InputHandle;
        use crate::pipeline::core::testing::StepProbe;

        let mut step = WindowProbeStep::new(u64::MAX);
        let mut harness = StepProbe::new(&step);
        let output = harness.take_output::<BamTemplateBatch>(0);
        for (serial, templates) in batches.into_iter().enumerate() {
            harness.push_input(BamTemplateBatch::new(serial as u64, templates));
        }
        harness.close_input();
        let mut released = Vec::new();
        let failed = loop {
            match harness.try_run(&mut step) {
                Err(_) => break true,
                Ok(StepOutcome::Finished) => break false,
                Ok(_) => {}
            }
            while let Some(batch) = output.pop() {
                released.push(batch.ordinal());
            }
        };
        while let Some(batch) = output.pop() {
            released.push(batch.ordinal());
        }
        (released, failed)
    }

    /// The step releases nothing until the probe has decided: a rejected input
    /// fails with nothing released (so nothing reaches the aligner or the
    /// output), and an accepted one is released whole, in order.
    #[rstest]
    #[case::short_input_without_rw(vec![batch(0, 3, None), batch(1, 2, None)], vec![], true)]
    #[case::long_input_without_rw(
        vec![batch(0, 600, None), batch(1, 600, None), batch(2, 5, None)], vec![], true
    )]
    #[case::window_in_second_batch(
        vec![batch(0, 3, None), batch(1, 3, Some(1)), batch(2, 3, None)], vec![0, 1, 2], false
    )]
    #[case::short_input_with_window(vec![batch(0, 3, Some(2))], vec![0], false)]
    #[case::empty_input(vec![], vec![], false)]
    fn window_probe_step_withholds_until_decided(
        #[case] batches: Vec<Vec<Template>>,
        #[case] expected_released: Vec<u64>,
        #[case] expected_failed: bool,
    ) {
        let (released, failed) = run_probe_step(batches);
        assert_eq!(failed, expected_failed);
        assert_eq!(released, expected_released);
    }

    /// R2: the align stage's zipper merge copies `rw` from the unmapped
    /// consensus record onto the mapped record, so `select_primary` sees it.
    #[test]
    fn zipper_merge_carries_rw_onto_mapped_records() {
        let unmapped = pair(Some(RW), Some("chr1:201-250:-"), None);
        let mapped_record = |flags: u16, pos: i32| {
            let mut b = SamBuilder::new();
            b.read_name(b"q1")
                .flags(flags)
                .ref_id(0)
                .pos(pos)
                .mapq(60)
                .cigar_ops(&[8 << 4])
                .sequence(b"ACGTACGT")
                .qualities(&[30; 8]);
            b.build()
        };
        let mut mapped = Template::from_records(vec![
            mapped_record(PAIRED | FIRST_SEGMENT, 100),
            mapped_record(PAIRED | LAST_SEGMENT | fgumi_raw_bam::flags::REVERSE, 200),
        ])
        .expect("template");
        let tag_info =
            TagInfo::new(vec![], vec!["Consensus".to_string()], vec!["Consensus".to_string()]);
        merge_raw(&unmapped, &mut mapped, &tag_info, false).expect("merge");
        let windows: Vec<Option<String>> = mapped
            .records
            .iter()
            .map(|r| {
                fgumi_raw_bam::find_string_tag(fgumi_raw_bam::aux_data_slice(r), SamTag::RW)
                    .map(|v| String::from_utf8_lossy(v).into_owned())
            })
            .collect();
        assert_eq!(windows, [Some(RW.to_string()), Some("chr1:201-250:-".to_string())]);
    }
}
