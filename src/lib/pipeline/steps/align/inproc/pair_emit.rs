//! `AlignPairEmitStep` — the `Parallel` pair + emit stage of the in-process
//! bwa-mem3 backend, and its last step: it also runs the zipper merge.
//!
//! Each pool thread borrows its [`AlignEngine::Scratch`] from the
//! [`ScratchPool`] it shares with the seed/extend step, and each worker keeps a
//! reusable record-collecting sink, so the per-item work allocates little after
//! warm-up.
//! Per [`PairWork`] the step:
//!
//! 1. runs the engine's `pair_emit`, which classifies pairs, rescues mates,
//!    marks primaries, and emits one packed BAM record **body** per output
//!    record (no `u32 block_size` prefix) through the sink, tagged with its
//!    [`RecordOrigin`] in bwa's emission order — pairs (R1 side then R2 side,
//!    primary then supplementary) then singles;
//! 2. groups those record bodies back into `Vec<Template>` in **input order** by
//!    mapping each origin through the sub-batch's [`Layout`] to its template
//!    index (for [`Layout::AllPairs`], `Pair(i)` is template `i`; for
//!    [`Layout::Mixed`], each origin resolves through the stored slots — this is
//!    where a pair Prepare split across a cohort boundary reappears as its two
//!    mid-pair-split `Single` templates), building each with
//!    [`Template::from_records`] exactly as the subprocess reader does; and
//! 3. pairs those `mapped` templates with the original `unmapped` half in a
//!    [`ZipperBatch`], then drops the sub-batch's clone of the cohort
//!    [`CohortLease`](super::gate::CohortLease) — moved through every stage
//!    rather than cloned-and-dropped early, so the cohort's
//!    [`CohortGate`](super::gate::CohortGate) reservation is released only once
//!    its last sub-batch has been emitted; and
//! 4. merges it with [`merge_zipper_batch`] — the same merge the shared
//!    `MergeAlignedStep` runs for the subprocess backend — while the records are
//!    still hot, emitting the `BamTemplateBatch` in input (serial) order.
//!
//! Three always-on guards mirror the subprocess reader's semantics
//! ([`super::super::subprocess`]): an input template with zero emitted records
//! (b), a record whose origin indexes beyond the sub-batch's templates (c), and
//! a per-template queryname mismatch between the mapped and unmapped halves (d).
//! They are error paths that never fire on valid input, but reproduce the
//! subprocess reader's protection against a backend that drops, adds, or
//! reorders records.
//!
//! A `RecordSink` body maps to a [`RawRecord`] with **no** `block_size` surgery:
//! `RawRecord` is a BAM record body *without* the 4-byte length prefix
//! (`crates/fgumi-raw-bam/src/raw_bam_record.rs`), which is exactly what the sink
//! receives, so the sink does `RawRecord::from(body.to_vec())` and nothing else.
//!
//! This `Parallel` step is the in-process backend's last step, wired after
//! [`CohortPeStatStep`](super::pestat::CohortPeStatStep) by
//! [`InProcessBwaMem3Backend`](super::InProcessBwaMem3Backend).

use std::io;
use std::sync::Arc;

use fgumi_raw_bam::RawRecord;

use crate::commands::zipper::ZipperTags;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::align::merge::{MergeConfig, merge_zipper_batch};
use crate::pipeline::steps::align::{ZipperBatch, is_primary_for_alignment};
use crate::pipeline::steps::types::BamTemplateBatch;
use crate::template::Template;

use super::cohort::{Layout, PairWork, TemplateSlot};
use super::engine::{AlignEngine, RecordOrigin, RecordSink, cohort_of};
use super::scratch::ScratchPool;

// ---------------------------------------------------------------------------
// Record-collecting sink
// ---------------------------------------------------------------------------

/// A reusable [`RecordSink`] that collects each emitted body as an owned
/// [`RawRecord`], tagged with its [`RecordOrigin`], in emission order.
///
/// The body a `RecordSink` receives is a packed BAM record with **no**
/// `block_size` prefix, which is exactly a [`RawRecord`]'s representation — so
/// `RawRecord::from(body.to_vec())` is the whole conversion (no prefix surgery).
/// [`Self::clear`] resets the length but retains the `Vec`'s capacity, so after
/// warm-up the sink allocates only the per-record body copies.
#[derive(Default)]
struct CollectingSink {
    records: Vec<(RecordOrigin, RawRecord)>,
}

impl CollectingSink {
    /// Reset for a new sub-batch, retaining the outer `Vec`'s capacity.
    fn clear(&mut self) {
        self.records.clear();
    }
}

impl RecordSink for CollectingSink {
    fn emit(&mut self, origin: RecordOrigin, body: &[u8]) {
        // No `block_size` surgery: the body IS a `RawRecord` (a BAM record body
        // without the 4-byte length prefix). The per-record copy is inherent:
        // each output `Template` owns its `RawRecord`s, and the zipper merge
        // then edits them in place (tags, flags), exactly as on the subprocess
        // route, whose reader also allocates one `RawRecord` per record.
        self.records.push((origin, RawRecord::from(body.to_vec())));
    }
}

// ---------------------------------------------------------------------------
// Origin -> template-index mapping (input order)
// ---------------------------------------------------------------------------

/// Dense inverse maps from a `-p` group ordinal to the template's input-order
/// index: `pair_to_template[i]` is the template index of pair `i`, and
/// `single_to_template[j]` that of single `j`.
///
/// For [`Layout::AllPairs`] every template is a pair, so pair `i` maps to
/// template `i` (the identity) and there are no singles — represented without
/// materializing the identity map. For [`Layout::Mixed`] the maps are built by
/// walking the stored slots in template order; because [`Layout::classify`]
/// assigns pair/single ordinals densely and in increasing order, pushing each
/// slot's template index in order yields the inverse map directly.
enum InverseLayout {
    /// Pair `i` -> template `i`; no singles. The identity map is not allocated —
    /// [`Self::template_index`] computes it (bounds-checked against
    /// `n_templates`) so the common all-pairs sub-batch does no heap work here.
    AllPairs { n_templates: usize },
    /// Explicit dense maps, in template order.
    Mixed { pair_to_template: Vec<usize>, single_to_template: Vec<usize> },
}

impl InverseLayout {
    fn build(layout: &Layout, n_templates: usize) -> Self {
        match layout {
            Layout::AllPairs => Self::AllPairs { n_templates },
            Layout::Mixed(slots) => {
                let mut pair_to_template = Vec::new();
                let mut single_to_template = Vec::new();
                for (t_idx, slot) in slots.iter().enumerate() {
                    match slot {
                        TemplateSlot::Pair { .. } => pair_to_template.push(t_idx),
                        TemplateSlot::Single { .. } => single_to_template.push(t_idx),
                    }
                }
                Self::Mixed { pair_to_template, single_to_template }
            }
        }
    }

    /// The input-order template index a record with this `origin` belongs to, or
    /// `None` if the origin indexes beyond the sub-batch's classification (the
    /// caller turns `None` into guard (c)'s trailing/extra-record error).
    fn template_index(&self, origin: RecordOrigin) -> Option<usize> {
        match self {
            // Identity for in-range pairs; an AllPairs sub-batch has no singles,
            // so any `Single` origin (or an out-of-range pair) resolves to `None`
            // exactly as an empty `single_to_template` would have.
            Self::AllPairs { n_templates } => match origin {
                RecordOrigin::Pair(i) if i < *n_templates => Some(i),
                _ => None,
            },
            Self::Mixed { pair_to_template, single_to_template } => match origin {
                RecordOrigin::Pair(i) => pair_to_template.get(i).copied(),
                RecordOrigin::Single(j) => single_to_template.get(j).copied(),
            },
        }
    }
}

/// Group the sink's emitted record bodies into one [`Template`] per input
/// template, in input order, running the three always-on guards.
///
/// Drains `sink.records` (retaining its capacity for the next item). The
/// per-template record lists preserve the sink's emission order, so each
/// template's records stay in the R1/R2, primary/supplementary order bwa emitted
/// them.
///
/// The three guards mirror the subprocess reader's protections
/// ([`super::super::subprocess`]) against a backend that drops, adds, or
/// reorders records — error paths that never fire on valid input:
///
/// - **(c) trailing/extra record** (`subprocess.rs` EOF trailing-record guard):
///   a record whose origin indexes beyond the sub-batch's `-p` classification
///   (or a resolved template index `>= unmapped.len()`) has no home template.
/// - **(b) fewer alignments** (`subprocess.rs` "fewer alignments" guard): an
///   input template that received zero emitted records, or whose emitted
///   primary count differs from its unmapped half's primary (input read) count
///   — bwa emits exactly one primary per input read, so the engine dropped (or
///   duplicated) a mate.
/// - **(d) queryname match** (`subprocess.rs` per-template queryname guard): the
///   assembled mapped template's name must equal the paired unmapped template's,
///   catching any residual out-of-order emission that origin indexing missed.
fn assemble_templates(
    sink: &mut CollectingSink,
    layout: &Layout,
    unmapped: &[Template],
) -> io::Result<Vec<Template>> {
    let n_templates = unmapped.len();
    let inverse = InverseLayout::build(layout, n_templates);

    let mut per_template: Vec<Vec<RawRecord>> = (0..n_templates).map(|_| Vec::new()).collect();
    for (origin, record) in sink.records.drain(..) {
        // Guard (c): an origin that resolves to no template (out-of-range group
        // ordinal) is an extra/trailing record — the engine emitted more than
        // the input reads. `InverseLayout::template_index` returns `None` for any
        // origin outside the layout's classification; the `< n_templates` filter
        // also bounds a resolved index, so a `Mixed` slot list longer than
        // `unmapped` errors here instead of panicking at `per_template[t]`.
        let t = inverse.template_index(origin).filter(|&t| t < n_templates).ok_or_else(|| {
            io::Error::other(format!(
                "align-and-merge (in-process) pair+emit: record origin {origin:?} indexes beyond \
                 the sub-batch's {n_templates} templates — the engine emitted an extra/trailing \
                 record (more alignments than input reads)"
            ))
        })?;
        per_template[t].push(record);
    }

    let mut mapped = Vec::with_capacity(n_templates);
    for (i, records) in per_template.into_iter().enumerate() {
        // Guard (b): a template with no emitted records means the engine dropped
        // one — fewer alignments than input reads.
        if records.is_empty() {
            return Err(io::Error::other(format!(
                "align-and-merge (in-process) pair+emit: template '{name}' produced no emitted \
                 records — the engine emitted fewer alignments than input reads",
                name = String::from_utf8_lossy(unmapped[i].name()),
            )));
        }
        // Guard (b), per read: bwa emits exactly one primary record per input
        // read (an unmapped one if nothing aligns), so the emitted primaries must
        // match the unmapped half's primaries — a pair whose R2 was dropped would
        // otherwise pass the emptiness check above as a one-mate template.
        let n_emitted = primary_count(&records);
        let n_input = primary_count(unmapped[i].records());
        if n_emitted != n_input {
            let relation = if n_emitted < n_input { "fewer" } else { "more" };
            return Err(io::Error::other(format!(
                "align-and-merge (in-process) pair+emit: template '{name}' produced \
                 {n_emitted} primary records for {n_input} input reads — the engine emitted \
                 {relation} alignments than input reads",
                name = String::from_utf8_lossy(unmapped[i].name()),
            )));
        }
        let template = Template::from_records(records).map_err(|e| {
            io::Error::other(format!(
                "align-and-merge (in-process) pair+emit: Template::from_records for template \
                 '{name}' failed: {e}",
                name = String::from_utf8_lossy(unmapped[i].name()),
            ))
        })?;
        // Guard (d): the mapped and unmapped halves are paired by position, so
        // their querynames must match — otherwise the engine emitted templates
        // out of input order.
        if template.name() != unmapped[i].name() {
            return Err(io::Error::other(format!(
                "align-and-merge (in-process) pair+emit: queryname mismatch — unmapped[{i}]='{u}' \
                 but mapped='{m}'; the engine emitted templates out of input order",
                u = String::from_utf8_lossy(unmapped[i].name()),
                m = String::from_utf8_lossy(template.name()),
            )));
        }
        mapped.push(template);
    }
    Ok(mapped)
}

/// The number of primary (non-secondary, non-supplementary) records in
/// `records`: one per input read on the unmapped side, and one per aligned read
/// on the emitted side.
fn primary_count(records: &[RawRecord]) -> usize {
    records.iter().filter(|record| is_primary_for_alignment(record.flags())).count()
}

// ---------------------------------------------------------------------------
// The step
// ---------------------------------------------------------------------------

/// The `Parallel` pair + emit step. Each worker copy carries
/// its own engine scratch and record-collecting sink; the engine itself is
/// shared by reference (`Arc<E>`, `E: Sync`).
pub(crate) struct AlignPairEmitStep<E: AlignEngine> {
    /// Shared, read-only alignment engine (index + options behind `Arc`s).
    engine: Arc<E>,
    /// Byte limit for the `ZipperBatch` output queue.
    output_byte_limit: u64,
    /// Engine scratches keyed by pool thread, shared with the seed/extend step
    /// so a thread keeps one warm scratch whichever step it runs.
    scratch: Arc<ScratchPool<E::Scratch>>,
    /// Per-worker, reused record-collecting sink.
    sink: CollectingSink,
    /// The zipper merge each sub-batch goes through right after it is emitted,
    /// and its tag bitsets, built once.
    merge: Arc<MergeConfig>,
    tags: Arc<ZipperTags>,
    /// A single output item bounced by backpressure, retried before new pushes.
    held: HeldSlot<Unpushed<BamTemplateBatch>>,
}

impl<E: AlignEngine> AlignPairEmitStep<E> {
    /// Build the step over a shared `engine`, scratch pool and merge
    /// configuration, with the given output byte limit.
    pub(crate) fn new(
        engine: Arc<E>,
        scratch: Arc<ScratchPool<E::Scratch>>,
        merge: Arc<MergeConfig>,
        output_byte_limit: u64,
    ) -> Self {
        let tags = Arc::new(ZipperTags::from_tag_info(&merge.tag_info));
        Self {
            engine,
            output_byte_limit,
            merge,
            tags,
            scratch,
            sink: CollectingSink::default(),
            held: HeldSlot::new(),
        }
    }

    /// Pair + emit one sub-batch and assemble its [`ZipperBatch`].
    ///
    /// Consumes `work`: its `ranges` are moved into the engine call (the reads
    /// they name stay in the cohort's resident state, freed with the lease's
    /// last clone), and `unmapped` is moved into the resulting `ZipperBatch`.
    /// The sub-batch's lease clone is dropped once its records are emitted and
    /// assembled — before the merge, as the merged batch no longer needs the
    /// cohort's resident state — so the cohort's gate reservation (and resident
    /// state) is freed when its last sub-batch gets here.
    fn process_item(&mut self, work: PairWork<E>) -> io::Result<ZipperBatch> {
        let PairWork { id, layout, unmapped, lease, ranges, pestat, ids } = work;

        self.sink.clear();
        let cohort = cohort_of(&*self.engine, &lease)?;
        let engine = &*self.engine;
        let sink = &mut self.sink;
        self.scratch.with(
            || {
                engine.new_scratch().map_err(|e| {
                    io::Error::other(format!(
                        "align-and-merge (in-process) pair+emit: allocating engine scratch: {e:#}"
                    ))
                })
            },
            |scratch| {
                engine.pair_emit(scratch, cohort, ranges, pestat.as_deref(), ids, sink).map_err(
                    |e| {
                        io::Error::other(format!(
                            "align-and-merge (in-process) pair+emit: engine pair_emit failed: \
                             {e:#}"
                        ))
                    },
                )
            },
        )?;

        let mapped = assemble_templates(&mut self.sink, &layout, &unmapped)?;
        drop(lease);

        let serial = id.serial;
        Ok(ZipperBatch { serial, mapped, unmapped: BamTemplateBatch::new(serial, unmapped) })
    }
}

impl<E: AlignEngine> Step for AlignPairEmitStep<E> {
    type Input = PairWork<E>;
    type Outputs = OrderedBytesSingle<BamTemplateBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "AlignPairEmit",
            kind: StepKind::Parallel,
            sticky: false,
            // Each sub-batch is merged here, so this step restores input order
            // itself from each batch's serial — what the shared
            // MergeAlignedStep does for the subprocess backend. `ByteBounded` is
            // mandatory under the armed deadlock monitor.
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // Drain the held output slot first, before popping more input.
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    // `Contention` (not `NoProgress`): a Parallel step returning
                    // `NoProgress` on a held slot risks the framework marking
                    // this worker `Skip` and dropping the held item. Mirrors
                    // `AlignSeedExtendStep` / `MergeAlignedStep`.
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(work) = ctx.input.pop() else {
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };

        // Merge while the sub-batch's records are still hot, instead of handing
        // them to a separate merge step through a queue.
        let merged = merge_zipper_batch(self.process_item(work)?, &self.merge, &self.tags)?;
        match ctx.outputs.push(merged) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }

    fn new_worker_copy(&self) -> Self {
        Self {
            engine: Arc::clone(&self.engine),
            output_byte_limit: self.output_byte_limit,
            // The shared scratch pool and merge config, and a fresh per-worker
            // sink.
            scratch: Arc::clone(&self.scratch),
            merge: Arc::clone(&self.merge),
            tags: Arc::clone(&self.tags),
            sink: CollectingSink::default(),
            held: HeldSlot::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::VecDeque;
    use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
    use std::sync::{Arc, Mutex};

    use fgumi_raw_bam::{SamBuilder, flags};
    use rstest::rstest;

    use super::*;
    use crate::pipeline::core::builder::{Pipeline, PipelineConfig};
    use crate::pipeline::core::item::Ordered;
    use crate::pipeline::core::outputs::Single;
    use crate::pipeline::steps::align::inproc::cohort::SubBatchId;
    use crate::pipeline::steps::align::inproc::engine::IdBases;
    use crate::pipeline::steps::align::inproc::engine::fake::{
        FakeEngine, TaggedRead, fake_record_name,
    };
    use crate::pipeline::steps::align::inproc::gate::{CohortGate, CohortLease};
    use crate::pipeline::steps::align::merge::MergeConfig;
    use crate::umi::TagInfo;

    // ----- fixtures --------------------------------------------------------

    /// A paired-end unmapped template (R1 + R2) with the given queryname.
    fn pair_template(name: &[u8]) -> Template {
        let mut r1 = SamBuilder::new();
        r1.read_name(name)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::UNMAPPED)
            .sequence(b"ACGT")
            .qualities(&[30u8, 30, 30, 30]);
        let mut r2 = SamBuilder::new();
        r2.read_name(name)
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::UNMAPPED)
            .sequence(b"ACGT")
            .qualities(&[30u8, 30, 30, 30]);
        Template::from_records(vec![r1.build(), r2.build()]).expect("pair template")
    }

    /// A single-end unmapped template with the given queryname.
    fn single_template(name: &[u8]) -> Template {
        let mut r = SamBuilder::new();
        r.read_name(name).flags(flags::UNMAPPED).sequence(b"ACGT").qualities(&[30u8, 30, 30, 30]);
        Template::from_records(vec![r.build()]).expect("single template")
    }

    /// The two `TaggedRead`s (R1, R2) `FakeEngine::seed_extend` would record for
    /// pair `i`.
    fn regs_pair(i: usize) -> [TaggedRead; 2] {
        [
            TaggedRead { origin: RecordOrigin::Pair(i), mate: 0 },
            TaggedRead { origin: RecordOrigin::Pair(i), mate: 1 },
        ]
    }

    /// The one `TaggedRead` `FakeEngine::seed_extend` would record for single `j`.
    fn regs_single(j: usize) -> TaggedRead {
        TaggedRead { origin: RecordOrigin::Single(j), mate: 0 }
    }

    /// Build a `PairWork<FakeEngine>` from its parts. `pestat` is `Some(())` here
    /// since these sub-batches classify pairs; `FakeEngine` never reads its value.
    fn make_pair_work(
        serial: u64,
        layout: Layout,
        unmapped: Vec<Template>,
        regs: Vec<TaggedRead>,
    ) -> PairWork<FakeEngine> {
        PairWork {
            id: SubBatchId { serial, cohort: 0, index_in_cohort: 0 },
            layout,
            unmapped,
            lease: CohortLease::for_test(),
            ranges: regs,
            pestat: Some(Arc::new(())),
            ids: IdBases::default(),
        }
    }

    // ----- (a) grouping by origin into templates in input order ------------

    /// The two input shapes exercised for grouping: the production `AllPairs`
    /// case and a `Mixed` case whose `Single` templates model a pair Prepare
    /// split across a cohort boundary into two mid-pair-split singles.
    #[derive(Clone, Copy)]
    enum Shape {
        AllPairs,
        MixedWithSplitSingles,
    }

    /// Build (layout, unmapped templates, regs, expected mapped names in input
    /// order) for a `Shape`. Every unmapped template's name equals the name
    /// `FakeEngine` will stamp on the records of its origin, so the name-match
    /// guard passes on this valid input.
    fn inputs_for(shape: Shape) -> (Layout, Vec<Template>, Vec<TaggedRead>, Vec<Vec<u8>>) {
        match shape {
            Shape::AllPairs => {
                let names: Vec<Vec<u8>> =
                    (0..3).map(|i| fake_record_name(RecordOrigin::Pair(i))).collect();
                let unmapped = names.iter().map(|n| pair_template(n)).collect();
                let mut regs = Vec::new();
                for i in 0..3 {
                    regs.extend_from_slice(&regs_pair(i));
                }
                (Layout::AllPairs, unmapped, regs, names)
            }
            Shape::MixedWithSplitSingles => {
                // Input order: single, pair, single. Emission order (pairs then
                // singles): Pair(0) R1/R2, then Single(0), Single(1).
                let layout = Layout::Mixed(vec![
                    TemplateSlot::Single { idx: 0 },
                    TemplateSlot::Pair { idx: 0 },
                    TemplateSlot::Single { idx: 1 },
                ]);
                let names = vec![
                    fake_record_name(RecordOrigin::Single(0)),
                    fake_record_name(RecordOrigin::Pair(0)),
                    fake_record_name(RecordOrigin::Single(1)),
                ];
                let unmapped = vec![
                    single_template(&names[0]),
                    pair_template(&names[1]),
                    single_template(&names[2]),
                ];
                let mut regs = Vec::new();
                regs.extend_from_slice(&regs_pair(0));
                regs.push(regs_single(0));
                regs.push(regs_single(1));
                (layout, unmapped, regs, names)
            }
        }
    }

    #[rstest]
    #[case::all_pairs(Shape::AllPairs)]
    #[case::mixed_with_split_singles(Shape::MixedWithSplitSingles)]
    fn records_grouped_by_origin_into_templates_in_input_order(#[case] shape: Shape) {
        let (layout, unmapped, regs, expected_names) = inputs_for(shape);
        let n_templates = unmapped.len();
        let mut step = AlignPairEmitStep::new(
            Arc::new(FakeEngine),
            Arc::default(),
            test_merge_config(),
            4 * 1024 * 1024,
        );

        let zb = step
            .process_item(make_pair_work(7, layout, unmapped, regs))
            .expect("valid input assembles a ZipperBatch");

        assert_eq!(zb.serial, 7, "the ZipperBatch keeps the sub-batch serial");
        assert_eq!(zb.mapped.len(), n_templates, "one mapped template per input template");
        assert_eq!(
            zb.unmapped.templates().len(),
            n_templates,
            "the unmapped half carries every input template"
        );

        let got_names: Vec<Vec<u8>> = zb.mapped.iter().map(|t| t.name().to_vec()).collect();
        assert_eq!(
            got_names, expected_names,
            "mapped templates are grouped by origin into input order"
        );
    }

    /// `process_item` drops the sub-batch's lease clone once the sub-batch is
    /// emitted and assembled, so the last sub-batch of a cohort frees the
    /// cohort's gate reservation here, and a clone still held elsewhere (another
    /// sub-batch of the same cohort) keeps it reserved.
    #[rstest]
    #[case::last_clone(false)]
    #[case::another_clone_alive(true)]
    fn process_item_drops_the_sub_batch_lease(#[case] keep_another_clone: bool) {
        let (layout, unmapped, regs, _) = inputs_for(Shape::AllPairs);
        let bound = 100;
        let gate = Arc::new(CohortGate::new(2 * bound));
        let lease = gate.try_acquire_cohort(bound).expect("an empty gate admits a cohort");
        let other = keep_another_clone.then(|| lease.clone());
        let mut work = make_pair_work(0, layout, unmapped, regs);
        work.lease = lease;
        assert_eq!(gate.in_flight_bytes(), bound, "the cohort is reserved before pair+emit");

        let mut step = AlignPairEmitStep::new(
            Arc::new(FakeEngine),
            Arc::default(),
            test_merge_config(),
            4 * 1024 * 1024,
        );
        let zb = step.process_item(work).expect("valid input assembles a ZipperBatch");

        let expected = if keep_another_clone { bound } else { 0 };
        assert_eq!(
            gate.in_flight_bytes(),
            expected,
            "the reservation is released only by the cohort's last lease clone"
        );
        drop((zb, other));
        assert_eq!(gate.in_flight_bytes(), 0, "every clone dropped frees the reservation");
    }

    /// Under the real `Parallel` runtime at several thread counts, every input
    /// `PairWork` yields exactly one merged `BamTemplateBatch`, delivered in
    /// serial order, with every emitted record counted.
    #[rstest]
    fn parallel_runtime_emits_one_merged_batch_per_pair_work_in_order(
        #[values(1, 2, 4)] threads: usize,
    ) {
        let inputs: Vec<PairWork<FakeEngine>> = (0..6u64)
            .map(|serial| {
                let names: Vec<Vec<u8>> =
                    (0..2).map(|i| fake_record_name(RecordOrigin::Pair(i))).collect();
                let unmapped = names.iter().map(|n| pair_template(n)).collect();
                let mut regs = Vec::new();
                for i in 0..2 {
                    regs.extend_from_slice(&regs_pair(i));
                }
                make_pair_work(serial, Layout::AllPairs, unmapped, regs)
            })
            .collect();
        let n_inputs = inputs.len();

        let collected: Arc<Mutex<Vec<BamTemplateBatch>>> = Arc::new(Mutex::new(Vec::new()));
        let merge = test_merge_config();
        let builder = Pipeline::builder();
        builder
            .chain(WorkSource { items: inputs.into(), held: HeldSlot::new() })
            .chain(AlignPairEmitStep::new(
                Arc::new(FakeEngine),
                Arc::default(),
                Arc::clone(&merge),
                4 * 1024 * 1024,
            ))
            .chain(CollectSink { collected: Arc::clone(&collected) })
            .into_sink_marker();
        let pipeline = builder.build().expect("chain builds");
        pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

        let out = collected.lock().expect("mutex not poisoned");
        assert_eq!(out.len(), n_inputs, "one merged batch per input PairWork");
        let serials: Vec<u64> = out.iter().map(Ordered::ordinal).collect();
        assert_eq!(serials, (0..n_inputs as u64).collect::<Vec<_>>(), "delivered in serial order");
        let n_records: usize = out
            .iter()
            .map(|b| {
                assert_eq!(b.templates().len(), 2, "two merged templates per sub-batch");
                b.templates().iter().map(|t| t.records.len()).sum::<usize>()
            })
            .sum();
        assert_eq!(merge.records_emitted.load(AtomicOrdering::Relaxed), n_records as u64);
    }

    // ----- (b) zero-emitted-records guard ----------------------------------

    /// An input template with zero emitted records is a hard error (the engine
    /// emitted fewer alignments than input reads). Here `unmapped` has three
    /// pair templates but `regs` omits pair 2, so template 2 gets no records.
    #[test]
    fn template_with_zero_emitted_records_is_a_hard_error() {
        let names: Vec<Vec<u8>> = (0..3).map(|i| fake_record_name(RecordOrigin::Pair(i))).collect();
        let unmapped: Vec<Template> = names.iter().map(|n| pair_template(n)).collect();
        // regs cover only pairs 0 and 1 — pair 2's template gets no records.
        let mut regs = Vec::new();
        regs.extend_from_slice(&regs_pair(0));
        regs.extend_from_slice(&regs_pair(1));

        let mut step = AlignPairEmitStep::new(
            Arc::new(FakeEngine),
            Arc::default(),
            test_merge_config(),
            4 * 1024 * 1024,
        );
        let err = step
            .process_item(make_pair_work(0, Layout::AllPairs, unmapped, regs))
            .expect_err("a template with no emitted records must error");
        let msg = err.to_string();
        assert!(
            msg.contains("fewer alignments"),
            "error mirrors the subprocess 'fewer alignments' guard: {msg}"
        );
    }

    /// A pair template that received only one mate's records is a hard error:
    /// the emitted primary count (1) differs from the unmapped half's (2), so
    /// the engine dropped R2. Without the per-read count this template would
    /// pass the emptiness check and merge as a one-mate template.
    #[test]
    fn pair_template_missing_a_mate_is_a_hard_error() {
        let name = fake_record_name(RecordOrigin::Pair(0));
        let unmapped = vec![pair_template(&name)];
        // Only mate 0 (R1) of Pair(0) is emitted.
        let regs = vec![regs_pair(0)[0]];

        let mut step = AlignPairEmitStep::new(
            Arc::new(FakeEngine),
            Arc::default(),
            test_merge_config(),
            4 * 1024 * 1024,
        );
        let err = step
            .process_item(make_pair_work(0, Layout::AllPairs, unmapped, regs))
            .expect_err("a pair with one emitted mate must error");
        let msg = err.to_string();
        assert!(
            msg.contains("produced 1 primary records for 2 input reads"),
            "error names both counts: {msg}"
        );
        assert!(msg.contains("fewer alignments"), "error is the 'fewer alignments' guard: {msg}");
    }

    // ----- (c) origin index beyond the templates ---------------------------

    /// A record whose origin indexes beyond the sub-batch's templates is a
    /// trailing/extra-record error. Here `unmapped` has two pair templates but
    /// `regs` includes a pair 2, whose `Pair(2)` origin has no template.
    #[test]
    fn origin_index_beyond_templates_is_a_trailing_record_error() {
        let names: Vec<Vec<u8>> = (0..2).map(|i| fake_record_name(RecordOrigin::Pair(i))).collect();
        let unmapped: Vec<Template> = names.iter().map(|n| pair_template(n)).collect();
        let mut regs = Vec::new();
        regs.extend_from_slice(&regs_pair(0));
        regs.extend_from_slice(&regs_pair(1));
        regs.extend_from_slice(&regs_pair(2)); // no template for Pair(2)

        let mut step = AlignPairEmitStep::new(
            Arc::new(FakeEngine),
            Arc::default(),
            test_merge_config(),
            4 * 1024 * 1024,
        );
        let err = step
            .process_item(make_pair_work(0, Layout::AllPairs, unmapped, regs))
            .expect_err("an out-of-range origin must error");
        let msg = err.to_string();
        assert!(
            msg.contains("trailing") || msg.contains("extra"),
            "error mirrors the subprocess trailing/extra-record guard: {msg}"
        );
    }

    /// A `Mixed` layout with more slots than `unmapped` templates must not index
    /// past `unmapped`: an origin resolving to the extra slot is the
    /// trailing/extra-record error, not an out-of-bounds panic. Here the layout
    /// classifies two pairs but only one unmapped template exists.
    #[test]
    fn mixed_layout_slot_beyond_unmapped_is_a_trailing_record_error() {
        let unmapped = vec![pair_template(&fake_record_name(RecordOrigin::Pair(0)))];
        let layout =
            Layout::Mixed(vec![TemplateSlot::Pair { idx: 0 }, TemplateSlot::Pair { idx: 1 }]);
        let mut regs = regs_pair(0).to_vec();
        regs.extend_from_slice(&regs_pair(1)); // resolves to slot 1: no unmapped template

        let mut step = AlignPairEmitStep::new(
            Arc::new(FakeEngine),
            Arc::default(),
            test_merge_config(),
            4 * 1024 * 1024,
        );
        let err = step
            .process_item(make_pair_work(0, layout, unmapped, regs))
            .expect_err("an origin past `unmapped` must error, not panic");
        let msg = err.to_string();
        assert!(msg.contains("extra/trailing"), "error is the trailing/extra-record guard: {msg}");
    }

    // ----- (d) per-template queryname mismatch -----------------------------

    /// If a mapped template's name does not equal the paired unmapped template's
    /// name, that is a queryname mismatch — the engine emitted records out of
    /// input order. Here the unmapped template is named "not-the-fake-name",
    /// while `FakeEngine` stamps `fake_record_name(Pair(0))`.
    #[test]
    fn queryname_mismatch_between_mapped_and_unmapped_is_an_error() {
        let unmapped = vec![pair_template(b"not-the-fake-name")];
        let regs = regs_pair(0).to_vec();

        let mut step = AlignPairEmitStep::new(
            Arc::new(FakeEngine),
            Arc::default(),
            test_merge_config(),
            4 * 1024 * 1024,
        );
        let err = step
            .process_item(make_pair_work(0, Layout::AllPairs, unmapped, regs))
            .expect_err("a queryname mismatch must error");
        let msg = err.to_string();
        assert!(
            msg.contains("queryname mismatch"),
            "error mirrors the subprocess queryname-match guard: {msg}"
        );
    }

    // ----- (e) RecordSink body -> RawRecord, no block_size surgery ----------

    /// A `RecordSink` body is a BAM record body with no `block_size` prefix, so
    /// the sink's `RawRecord::from(body.to_vec())` round-trips: the collected
    /// `RawRecord` equals the record whose bytes were emitted, byte-for-byte.
    #[test]
    fn record_sink_body_maps_to_raw_record_without_prefix_surgery() {
        let mut builder = SamBuilder::new();
        builder
            .read_name(b"roundtrip")
            .flags(flags::UNMAPPED)
            .sequence(b"ACGTAC")
            .qualities(&[30u8, 31, 32, 33, 34, 35]);
        let original: RawRecord = builder.build();

        let mut sink = CollectingSink::default();
        sink.emit(RecordOrigin::Single(0), original.as_ref());

        assert_eq!(sink.records.len(), 1, "one body collected");
        let (origin, collected) = &sink.records[0];
        assert_eq!(*origin, RecordOrigin::Single(0), "the origin tag is preserved");
        assert_eq!(
            collected, &original,
            "RawRecord::from(body) equals the emitted record (no block_size surgery)"
        );
        assert_eq!(collected.read_name(), b"roundtrip", "the parsed name round-trips");
    }

    // ----- profile / worker-copy -------------------------------------------

    #[test]
    fn profile_is_parallel_bytebounded_ordered_by_serial() {
        let step = AlignPairEmitStep::new(
            Arc::new(FakeEngine),
            Arc::default(),
            test_merge_config(),
            4 * 1024 * 1024,
        );
        let p = step.profile();
        assert_eq!(p.name, "AlignPairEmit");
        assert_eq!(p.kind, StepKind::Parallel);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
        assert!(matches!(p.output_queues[0], QueueSpec::ByteBounded { .. }));
    }

    #[test]
    fn new_worker_copy_shares_engine_and_resets_per_worker_state() {
        let step = AlignPairEmitStep::new(
            Arc::new(FakeEngine),
            Arc::default(),
            test_merge_config(),
            4 * 1024 * 1024,
        );
        let copy = step.new_worker_copy();
        assert!(Arc::ptr_eq(&step.engine, &copy.engine), "engine Arc is shared");
        assert!(
            Arc::ptr_eq(&step.scratch, &copy.scratch),
            "worker copies share the step's per-thread scratch pool"
        );
        assert!(copy.sink.records.is_empty(), "worker copy gets a fresh sink");
        assert_eq!(copy.output_byte_limit, step.output_byte_limit);
    }

    // ----- test source / sink (real-runtime tests) -------------------------

    /// Replays pre-built `PairWork`s into the chain. `Exclusive` (owns cursor).
    struct WorkSource {
        items: VecDeque<PairWork<FakeEngine>>,
        held: HeldSlot<Unpushed<PairWork<FakeEngine>>>,
    }

    impl Step for WorkSource {
        type Input = ();
        type Outputs = Single<PairWork<FakeEngine>>;

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "WorkSource",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 1 << 20 }],
                branch_ordering: vec![BranchOrdering::None],
            }
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            if let Some(unpushed) = self.held.take()
                && let Err(again) = ctx.outputs.retry(unpushed)
            {
                self.held.put(again);
                return Ok(StepOutcome::Contention);
            }
            let Some(item) = self.items.pop_front() else {
                return Ok(StepOutcome::Finished);
            };
            if let Err(unpushed) = ctx.outputs.push(item) {
                self.held.put(unpushed);
            }
            Ok(StepOutcome::Progress)
        }
    }

    /// Terminal sink accumulating every merged batch in arrival order.
    struct CollectSink {
        collected: Arc<Mutex<Vec<BamTemplateBatch>>>,
    }

    /// A merge configuration with no tag rules, the one the chain builder uses.
    fn test_merge_config() -> Arc<MergeConfig> {
        Arc::new(MergeConfig {
            tag_info: Arc::new(TagInfo::new(vec![], vec![], vec![])),
            skip_tc_tags: false,
            reference: None,
            partial_output_header: Arc::new(noodles::sam::Header::default()),
            records_emitted: Arc::new(AtomicU64::new(0)),
            output_byte_limit: 4 * 1024 * 1024,
        })
    }

    impl Step for CollectSink {
        type Input = BamTemplateBatch;
        type Outputs = ();

        fn profile(&self) -> StepProfile {
            StepProfile {
                name: "CollectSink",
                kind: StepKind::Exclusive,
                sticky: false,
                output_queues: vec![],
                branch_ordering: vec![],
            }
        }

        fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
            match ctx.input.pop() {
                Some(item) => {
                    self.collected.lock().expect("mutex not poisoned").push(item);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }
}
