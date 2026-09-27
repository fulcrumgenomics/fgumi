//! `AlignSeedExtendStep` — the `Parallel` seed + SE-extension stage of the
//! in-process bwa-mem3 backend.
//!
//! Each pool thread borrows its [`AlignEngine::Scratch`] from the [`ScratchPool`]
//! it shares with the pair/emit step, and each worker keeps a reusable
//! [`ReadArena`], so the per-item work — unpack every surviving read's SEQ/QUAL,
//! run the engine's `seed_extend`, and forward the result — allocates nothing after
//! warm-up. The step consumes an [`AlignWork`] (from
//! [`AlignPrepareStep`](super::prepare::AlignPrepareStep)) and produces an
//! [`ExtendedWork`], moving the [`CohortLease`](super::gate::CohortLease) forward
//! untouched (dropping it early would release the cohort gate before the cohort
//! finishes and over-admit).
//!
//! **Read selection is byte-identical to Prepare and to the subprocess FASTQ
//! writer.** Both backends must feed bwa-mem3 the same reads, or their outputs
//! diverge, so this step reconstructs, from `AlignWork` alone, exactly the reads
//! Prepare classified:
//!
//! - it walks each template's records in record order, keeping only the
//!   **primary** reads (the shared
//!   [`is_primary_for_alignment`]), and
//! - it groups them into bwa's `-p` SE/PE batch via the sub-batch's
//!   [`Layout`]: a `Pair` template contributes its two primaries as an
//!   [`EnginePair`] (R1 then R2), a `Single` template its one primary as an
//!   [`EngineRead`]. A pair that Prepare split across a mid-pair cohort boundary
//!   is already two single-record templates in `unmapped`, each classified as a
//!   `Single` — so it reconstructs with no special case here.
//!
//! Each read's SEQ and QUAL are decoded by the same `decode_fastq_seq_qual`
//! (`commands/fastq.rs`) that `write_fastq_record` uses, **including the
//! reverse-complement of a `REVERSE`-flagged record**, so the bytes bwa sees
//! match the subprocess path by construction.
//!
//! This `Parallel` step is wired between
//! [`AlignPrepareStep`](super::prepare::AlignPrepareStep) and
//! [`CohortPeStatStep`](super::pestat::CohortPeStatStep) by
//! [`InProcessBwaMem3Backend`](super::InProcessBwaMem3Backend).

use std::io;
use std::sync::Arc;

use fgumi_raw_bam::RawRecord;

use crate::commands::fastq::decode_fastq_seq_qual;
use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::Single;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::align::is_primary_for_alignment;
use crate::template::Template;

use super::cohort::{AlignWork, ExtendedWork, Layout, TemplateSlot};
use super::engine::{AlignEngine, EngineBatch, EnginePair, EngineRead, cohort_of};
use super::scratch::ScratchPool;

/// Estimated resident bytes of one [`ExtendedWork`] at the default sub-batch
/// size: about 1 MiB, originally ≈ 256 pairs times the
/// per-pair unmapped (≈ 900 B), seq (≈ 700 B), and regs (≈ 600 B) footprint.
/// The seq and regs now live in the cohort's resident state, so an item is
/// charged only its templates (≈ 460 KiB at the 512-pair aarch64 default); kept
/// as the conservative sizing until re-tuned against measurements.
/// Used only to size the step's output queue so that `T` worker clones can each
/// hold one pushed item without the queue's byte limit being the bottleneck.
const EST_EXTENDED_WORK_BYTES: u64 = 1024 * 1024;

/// Byte limit for the `SeedExtend → PeStat` queue:
/// `max(per_step_byte_limit, 2·T·EST_EXTENDED_WORK_BYTES)`, so that `T` clones
/// pushing at once do not all hold and stall the pool. `threads` is the pool
/// worker count `T`; it is clamped to `>= 1`.
pub(crate) fn seed_extend_output_byte_limit(per_step_byte_limit: u64, threads: usize) -> u64 {
    let t = u64::try_from(threads.max(1)).unwrap_or(u64::MAX);
    let clones = 2u64.saturating_mul(t).saturating_mul(EST_EXTENDED_WORK_BYTES);
    per_step_byte_limit.max(clones)
}

// ---------------------------------------------------------------------------
// Per-worker read arena
// ---------------------------------------------------------------------------

/// Byte spans of one unpacked read inside a [`ReadArena`]'s `bytes` buffer.
/// SEQ and QUAL are appended contiguously; the engine batch borrows the two
/// slices these spans describe.
#[derive(Clone, Copy, Debug)]
struct ReadSpan {
    seq_off: usize,
    seq_len: usize,
    qual_off: usize,
    qual_len: usize,
}

/// A per-worker, reused scratch arena for one sub-batch's unpacked reads.
///
/// `bytes` holds every read's decoded SEQ + QUAL back-to-back; `spans` records
/// where each read's two slices sit. `seq_scratch`/`qual_scratch` are the
/// decode buffers `decode_fastq_seq_qual` writes into
/// (both clear their destination), whose contents are then appended into
/// `bytes`. [`ReadArena::clear`] resets the lengths but retains every buffer's
/// capacity, so after warm-up a sub-batch allocates nothing.
struct ReadArena {
    bytes: Vec<u8>,
    spans: Vec<ReadSpan>,
    seq_scratch: Vec<u8>,
    qual_scratch: Vec<u8>,
}

impl ReadArena {
    fn new() -> Self {
        Self {
            bytes: Vec::new(),
            spans: Vec::new(),
            seq_scratch: Vec::new(),
            qual_scratch: Vec::new(),
        }
    }

    /// Reset the arena for a new sub-batch, retaining every buffer's capacity.
    fn clear(&mut self) {
        self.bytes.clear();
        self.spans.clear();
        // seq/qual scratch are cleared by the decode helpers on next use.
    }

    fn seq_of(&self, span: &ReadSpan) -> &[u8] {
        &self.bytes[span.seq_off..span.seq_off + span.seq_len]
    }

    fn qual_of(&self, span: &ReadSpan) -> &[u8] {
        &self.bytes[span.qual_off..span.qual_off + span.qual_len]
    }
}

// ---------------------------------------------------------------------------
// Read reconstruction (shared with Prepare's selection, verified byte-parity)
// ---------------------------------------------------------------------------

/// A template's primary reads in record order — the exact reads Prepare fed to
/// the cutter and the subprocess FASTQ writer emits.
fn primary_records(template: &Template) -> impl Iterator<Item = &RawRecord> {
    template.records().iter().filter(|record| is_primary_for_alignment(record.flags()))
}

/// Whether the `t_idx`-th template of a sub-batch is a PE pair, per its
/// [`Layout`] (`AllPairs` ⇒ always a pair; `Mixed` ⇒ its stored slot).
fn template_is_pair(layout: &Layout, t_idx: usize) -> io::Result<bool> {
    match layout {
        Layout::AllPairs => Ok(true),
        Layout::Mixed(slots) => match slots.get(t_idx) {
            Some(TemplateSlot::Pair { .. }) => Ok(true),
            Some(TemplateSlot::Single { .. }) => Ok(false),
            None => Err(io::Error::other(format!(
                "align-and-merge (in-process) seed/extend: layout carries no slot for template \
                 index {t_idx}"
            ))),
        },
    }
}

/// Decode every surviving read's SEQ + QUAL into `arena`, one [`ReadSpan`] per
/// primary read, in template then record order.
///
/// The SEQ/QUAL bytes come from the same [`decode_fastq_seq_qual`] the
/// subprocess `write_fastq_record` uses, so both backends feed bwa-mem3
/// identical bytes.
fn fill_arena(arena: &mut ReadArena, work: &AlignWork) {
    arena.clear();
    for template in &work.unmapped {
        for record in primary_records(template) {
            decode_fastq_seq_qual(
                record,
                record.flags(),
                &mut arena.seq_scratch,
                &mut arena.qual_scratch,
            );

            let seq_off = arena.bytes.len();
            arena.bytes.extend_from_slice(&arena.seq_scratch);
            let seq_len = arena.seq_scratch.len();

            let qual_off = arena.bytes.len();
            arena.bytes.extend_from_slice(&arena.qual_scratch);
            let qual_len = arena.qual_scratch.len();

            arena.spans.push(ReadSpan { seq_off, seq_len, qual_off, qual_len });
        }
    }
}

/// Nth span, or a descriptive error if the arena and the layout disagree on the
/// read count (an internal inconsistency, not user input).
fn span_at(arena: &ReadArena, idx: usize) -> io::Result<&ReadSpan> {
    arena.spans.get(idx).ok_or_else(|| {
        io::Error::other(format!(
            "align-and-merge (in-process) seed/extend: read span {idx} missing — arena and layout \
             disagree on the sub-batch's read count"
        ))
    })
}

/// Build the engine's SE/PE batch, borrowing the filled `arena` for SEQ/QUAL and
/// the `work` templates for read names, in the exact order Prepare classified
/// them: pairs (R1, R2) in pair-index order, singles in single-index order.
fn build_engine_reads<'a>(
    arena: &'a ReadArena,
    work: &'a AlignWork,
) -> io::Result<(Vec<EnginePair<'a>>, Vec<EngineRead<'a>>)> {
    let mut pairs: Vec<EnginePair<'a>> = Vec::new();
    let mut singles: Vec<EngineRead<'a>> = Vec::new();
    let mut next_span = 0usize;

    for (t_idx, template) in work.unmapped.iter().enumerate() {
        let mut prims = primary_records(template);
        if template_is_pair(&work.layout, t_idx)? {
            let (Some(r1), Some(r2)) = (prims.next(), prims.next()) else {
                return Err(io::Error::other(format!(
                    "align-and-merge (in-process) seed/extend: template '{name}' is classified as a \
                     pair but has fewer than two primary reads",
                    name = String::from_utf8_lossy(template.name()),
                )));
            };
            let s1 = span_at(arena, next_span)?;
            let s2 = span_at(arena, next_span + 1)?;
            next_span += 2;
            pairs.push(EnginePair {
                r1: EngineRead::new(r1.read_name(), arena.seq_of(s1), Some(arena.qual_of(s1))),
                r2: EngineRead::new(r2.read_name(), arena.seq_of(s2), Some(arena.qual_of(s2))),
            });
        } else {
            let Some(read) = prims.next() else {
                return Err(io::Error::other(format!(
                    "align-and-merge (in-process) seed/extend: template '{name}' is classified as a \
                     single but has no primary read",
                    name = String::from_utf8_lossy(template.name()),
                )));
            };
            let s = span_at(arena, next_span)?;
            next_span += 1;
            singles.push(EngineRead::new(
                read.read_name(),
                arena.seq_of(s),
                Some(arena.qual_of(s)),
            ));
        }
    }

    if next_span != arena.spans.len() {
        return Err(io::Error::other(format!(
            "align-and-merge (in-process) seed/extend: consumed {next_span} of {total} unpacked \
             reads — arena and layout disagree on the sub-batch's read count",
            total = arena.spans.len(),
        )));
    }

    Ok((pairs, singles))
}

// ---------------------------------------------------------------------------
// The step
// ---------------------------------------------------------------------------

/// The `Parallel` seed + SE-extension step. Each worker copy
/// carries its own read arena and borrows its pool thread's engine scratch from
/// the shared [`ScratchPool`]; the engine itself is shared by reference
/// (`Arc<E>`, `E: Sync`).
pub(crate) struct AlignSeedExtendStep<E: AlignEngine> {
    /// Shared, read-only alignment engine (index + options behind `Arc`s).
    engine: Arc<E>,
    /// Byte limit for the `ExtendedWork` output queue.
    output_byte_limit: u64,
    /// Engine scratches keyed by pool thread, shared with the pair/emit step so
    /// a thread keeps one warm scratch whichever step it runs.
    scratch: Arc<ScratchPool<E::Scratch>>,
    /// Per-worker, reused read arena.
    arena: ReadArena,
    /// A single output item bounced by backpressure, retried before new pushes.
    held: HeldSlot<Unpushed<ExtendedWork<E>>>,
}

impl<E: AlignEngine> AlignSeedExtendStep<E> {
    /// Build the step over a shared `engine` and scratch pool with the given
    /// output byte limit (size it via [`seed_extend_output_byte_limit`]).
    pub(crate) fn new(
        engine: Arc<E>,
        scratch: Arc<ScratchPool<E::Scratch>>,
        output_byte_limit: u64,
    ) -> Self {
        Self { engine, output_byte_limit, scratch, arena: ReadArena::new(), held: HeldSlot::new() }
    }

    /// Seed + SE-extend one sub-batch: unpack its reads into the arena, run the
    /// engine against the cohort's resident state (created through the lease by
    /// whichever sub-batch of the cohort arrives first), and move the
    /// carried-forward fields (id/pos/closer/layout/unmapped/**lease**) plus the
    /// sub-batch's ranges into an [`ExtendedWork`].
    fn process_item(&mut self, work: AlignWork) -> io::Result<ExtendedWork<E>> {
        fill_arena(&mut self.arena, &work);
        let (pairs, singles) = build_engine_reads(&self.arena, &work)?;
        let batch = EngineBatch { pairs: &pairs, singles: &singles };

        let cohort = cohort_of(&*self.engine, &work.lease)?;
        let engine = &*self.engine;
        let ranges = self.scratch.with(
            || {
                engine.new_scratch().map_err(|e| {
                    io::Error::other(format!(
                        "align-and-merge (in-process) seed/extend: allocating engine scratch: \
                         {e:#}"
                    ))
                })
            },
            |scratch| {
                engine.seed_extend(scratch, cohort, batch).map_err(|e| {
                    io::Error::other(format!(
                        "align-and-merge (in-process) seed/extend: engine seed_extend failed: \
                         {e:#}"
                    ))
                })
            },
        )?;

        // End the arena/record borrows before moving `work`'s fields out, then
        // reset the arena (capacity retained) for the next item.
        drop((pairs, singles));
        self.arena.clear();

        let AlignWork { id, pos, closer, layout, unmapped, lease } = work;
        Ok(ExtendedWork { id, pos, closer, layout, unmapped, lease, ranges })
    }
}

impl<E: AlignEngine> Step for AlignSeedExtendStep<E> {
    type Input = AlignWork;
    type Outputs = Single<ExtendedWork<E>>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "AlignSeedExtend",
            kind: StepKind::Parallel,
            sticky: false,
            // Unordered: PeStat reassembles by (cohort, index_in_cohort), so
            // the downstream consumer does not depend on seed/extend output
            // order. `ByteBounded` is mandatory under the armed deadlock monitor.
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
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
                    // `MergeAlignedStep`.
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

        let extended = self.process_item(work)?;
        match ctx.outputs.push(extended) {
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
            // The shared scratch pool, and a fresh per-worker arena.
            scratch: Arc::clone(&self.scratch),
            arena: ReadArena::new(),
            held: HeldSlot::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::VecDeque;
    use std::sync::{Arc, Mutex};

    use fgumi_raw_bam::{SamBuilder, flags};
    use rstest::rstest;

    use super::*;
    use crate::pipeline::core::builder::{Pipeline, PipelineConfig};
    use crate::pipeline::core::item::HeapSize;
    use crate::pipeline::steps::align::inproc::cohort::{CohortPos, SubBatchId};
    use crate::pipeline::steps::align::inproc::engine::RecordOrigin;
    use crate::pipeline::steps::align::inproc::engine::fake::FakeEngine;
    use crate::pipeline::steps::align::inproc::gate::CohortLease;

    // ----- fixtures --------------------------------------------------------

    /// Cyclic ASCII sequence of `len` bases (never a bare 2-byte literal).
    fn seq_of(len: u32) -> Vec<u8> {
        b"ACGTN".iter().copied().cycle().take(len as usize).collect()
    }

    /// A single-end template: one primary record of `len` bases.
    fn single_template(name: &[u8], len: u32) -> Template {
        let seq = seq_of(len);
        let mut b = SamBuilder::new();
        b.read_name(name)
            .flags(flags::UNMAPPED)
            .sequence(&seq)
            .qualities(&vec![30u8; len as usize]);
        Template::from_records(vec![b.build()]).expect("single template")
    }

    /// A single-end template flagged REVERSE (to exercise the RC unpack path).
    fn reverse_single_template(name: &[u8], seq: &[u8], quals: &[u8]) -> Template {
        let mut b = SamBuilder::new();
        b.read_name(name).flags(flags::UNMAPPED | flags::REVERSE).sequence(seq).qualities(quals);
        Template::from_records(vec![b.build()]).expect("reverse single template")
    }

    /// A paired-end template: R1 + R2, each `len` bases.
    fn pair_template(name: &[u8], len: u32) -> Template {
        let seq = seq_of(len);
        let quals = vec![30u8; len as usize];
        let mut r1 = SamBuilder::new();
        r1.read_name(name)
            .flags(flags::PAIRED | flags::FIRST_SEGMENT | flags::UNMAPPED)
            .sequence(&seq)
            .qualities(&quals);
        let mut r2 = SamBuilder::new();
        r2.read_name(name)
            .flags(flags::PAIRED | flags::LAST_SEGMENT | flags::UNMAPPED)
            .sequence(&seq)
            .qualities(&quals);
        Template::from_records(vec![r1.build(), r2.build()]).expect("pair template")
    }

    /// Build an `AlignWork` from a sub-batch's templates, classifying its layout
    /// from the primary-read counts exactly as Prepare does.
    // `n_se`/`n_pe` are the domain terms (the `CohortPos` field names); the pair
    // reads more clearly together than renamed apart.
    #[allow(clippy::similar_names)]
    fn make_work(serial: u64, templates: Vec<Template>) -> AlignWork {
        let counts: Vec<u8> = templates
            .iter()
            .map(|t| {
                u8::try_from(primary_records(t).count()).expect("<= 2 primary reads per template")
            })
            .collect();
        let layout = Layout::classify(&counts);
        let mut n_pe = 0u32;
        let mut n_se = 0u32;
        for &count in &counts {
            if count == 2 {
                n_pe += 1;
            } else {
                n_se += 1;
            }
        }
        AlignWork {
            id: SubBatchId { serial, cohort: 0, index_in_cohort: u32::try_from(serial).unwrap() },
            pos: CohortPos { cohort_read_base: 0, se_offset: 0, pe_offset: 0, n_se, n_pe },
            closer: None,
            layout,
            unmapped: templates,
            lease: CohortLease::for_test(),
        }
    }

    /// The record origins `FakeEngine` tags for a work item: every pair (R1, R2)
    /// in pair-index order, then every single in single-index order (bwa's
    /// emission order).
    #[allow(clippy::similar_names)] // `n_se`/`n_pe` are the domain terms.
    fn expected_origins(work: &AlignWork) -> Vec<RecordOrigin> {
        let mut n_pe = 0usize;
        let mut n_se = 0usize;
        for (t_idx, _) in work.unmapped.iter().enumerate() {
            if template_is_pair(&work.layout, t_idx).unwrap() {
                n_pe += 1;
            } else {
                n_se += 1;
            }
        }
        let mut origins = Vec::new();
        for i in 0..n_pe {
            origins.push(RecordOrigin::Pair(i));
            origins.push(RecordOrigin::Pair(i));
        }
        for j in 0..n_se {
            origins.push(RecordOrigin::Single(j));
        }
        origins
    }

    // ----- test source / sink ---------------------------------------------

    /// Replays pre-built `AlignWork`s into the chain. `Exclusive` (owns cursor).
    struct WorkSource {
        items: VecDeque<AlignWork>,
        held: HeldSlot<Unpushed<AlignWork>>,
    }

    impl Step for WorkSource {
        type Input = ();
        type Outputs = Single<AlignWork>;

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

    /// Terminal sink accumulating every `ExtendedWork` in arrival order.
    struct CollectSink {
        collected: Arc<Mutex<Vec<ExtendedWork<FakeEngine>>>>,
    }

    impl Step for CollectSink {
        type Input = ExtendedWork<FakeEngine>;
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
                    self.collected.lock().expect("sink mutex not poisoned").push(item);
                    Ok(StepOutcome::Progress)
                }
                None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
                None => Ok(StepOutcome::NoProgress),
            }
        }
    }

    // ----- tests -----------------------------------------------------------

    /// End-to-end on the real pipeline at 1–8 threads: every `AlignWork` yields
    /// exactly one `ExtendedWork`, and its regs tag the same input read indices,
    /// in emission order (pairs then singles) — including a mixed sub-batch
    /// whose singles model mid-pair-split reads.
    #[rstest]
    fn each_work_yields_one_extended_work_tagging_its_input_reads(
        #[values(1, 2, 4, 8)] threads: usize,
    ) {
        // A mix of shapes: all-pairs, single-only, and mixed sub-batches. The
        // mixed sub-batch's `Single` templates are single-record templates,
        // exactly how Prepare represents a mid-pair-split read.
        let inputs: Vec<AlignWork> = vec![
            make_work(0, vec![pair_template(b"p0a", 12), pair_template(b"p0b", 12)]),
            make_work(
                1,
                vec![
                    single_template(b"split_left", 20),
                    pair_template(b"mid_pair", 15),
                    single_template(b"split_right", 20),
                ],
            ),
            make_work(2, vec![pair_template(b"solo_pair", 30)]),
            make_work(3, vec![single_template(b"sng0", 10), single_template(b"sng1", 10)]),
            make_work(
                4,
                (0..5).map(|i| pair_template(format!("bulk{i}").as_bytes(), 40)).collect(),
            ),
        ];
        // Per input serial: the expected regs origins and the template count.
        let expected: Vec<(u64, Vec<RecordOrigin>, usize)> =
            inputs.iter().map(|w| (w.id.serial, expected_origins(w), w.unmapped.len())).collect();
        let n_inputs = inputs.len();

        let collected: Arc<Mutex<Vec<ExtendedWork<FakeEngine>>>> = Arc::new(Mutex::new(Vec::new()));
        let limit = seed_extend_output_byte_limit(4 * 1024 * 1024, threads);

        let builder = Pipeline::builder();
        builder
            .chain(WorkSource { items: inputs.into(), held: HeldSlot::new() })
            .chain(AlignSeedExtendStep::new(Arc::new(FakeEngine), Arc::default(), limit))
            .chain(CollectSink { collected: Arc::clone(&collected) })
            .into_sink_marker();
        let pipeline = builder.build().expect("chain builds");
        pipeline.run(PipelineConfig { threads, ..Default::default() }).expect("chain runs");

        let out = collected.lock().expect("mutex not poisoned");
        assert_eq!(out.len(), n_inputs, "one ExtendedWork per input AlignWork");

        for (serial, want_origins, n_templates) in expected {
            let ew = out
                .iter()
                .find(|w| w.id.serial == serial)
                .unwrap_or_else(|| panic!("no ExtendedWork for serial {serial}"));
            let got: Vec<RecordOrigin> = ew.ranges.iter().map(|tagged| tagged.origin).collect();
            assert_eq!(got, want_origins, "serial {serial}: regs tag the input reads in order");
            // The carried-forward templates survive unchanged.
            assert_eq!(
                ew.unmapped.len(),
                n_templates,
                "serial {serial}: templates carried forward"
            );
        }
    }

    /// Driving two items through one step reuses the arena: after each item the
    /// arena is emptied (len 0) but its capacity is retained across items.
    #[test]
    fn arena_capacity_is_retained_across_items() {
        let mut step = AlignSeedExtendStep::new(Arc::new(FakeEngine), Arc::default(), 1 << 20);

        // A large first item grows the arena buffers.
        let big = make_work(
            0,
            (0..8).map(|i| pair_template(format!("big{i}").as_bytes(), 120)).collect(),
        );
        let ew0 = step.process_item(big).expect("process big");
        assert_eq!(FakeEngine::n_pairs(&ew0.ranges), 8);
        let bytes_cap = step.arena.bytes.capacity();
        let spans_cap = step.arena.spans.capacity();
        assert!(bytes_cap > 0, "arena grew for the first item");
        assert_eq!(step.arena.bytes.len(), 0, "arena is emptied after an item");
        assert_eq!(step.arena.spans.len(), 0, "spans are emptied after an item");

        // A small second item must not shrink the retained capacity.
        let small = make_work(1, vec![single_template(b"tiny", 4)]);
        let ew1 = step.process_item(small).expect("process small");
        assert_eq!(ew1.ranges.len(), 1, "one read for the single-template item");
        assert!(
            step.arena.bytes.capacity() >= bytes_cap,
            "arena byte capacity is retained (reused), not reallocated smaller"
        );
        assert!(step.arena.spans.capacity() >= spans_cap, "arena span capacity is retained");
        assert_eq!(step.arena.bytes.len(), 0, "arena emptied after the small item too");
    }

    /// Sub-batches of one cohort (clones of one lease) seed into a single
    /// resident cohort, created by whichever arrives first, and each item is
    /// charged only its moved templates: the resident reads live in the cohort.
    #[test]
    fn sub_batches_of_a_cohort_share_one_resident_cohort() {
        let mut step = AlignSeedExtendStep::new(Arc::new(FakeEngine), Arc::default(), 1 << 20);
        let lease = CohortLease::for_test();
        let mut first = make_work(0, vec![pair_template(b"hp0", 20), pair_template(b"hp1", 20)]);
        let mut second =
            make_work(1, vec![pair_template(b"hp2", 20), single_template(b"sng0", 20)]);
        first.lease = lease.clone();
        second.lease = lease.clone();

        let ew0 = step.process_item(first).expect("process first");
        let ew1 = step.process_item(second).expect("process second");

        let cohort = cohort_of(&FakeEngine, &lease).expect("cohort exists");
        assert_eq!(
            cohort.seeded_pairs.load(std::sync::atomic::Ordering::Relaxed),
            3,
            "both sub-batches seeded into the one cohort",
        );
        for ew in [&ew0, &ew1] {
            let templates_bytes: usize = ew.unmapped.iter().map(Template::heap_size).sum();
            assert_eq!(ew.heap_size(), templates_bytes, "heap_size = templates only");
        }
    }

    /// The `REVERSE` unpack path matches `write_fastq_record`: SEQ is
    /// reverse-complemented and QUAL is reversed into the arena.
    #[test]
    fn reverse_flagged_read_is_reverse_complemented_in_the_arena() {
        // Not a palindrome: identity, reverse-only, complement-only and RC
        // all differ, so dropping the REVERSE branch fails the SEQ assert.
        let seq: &[u8] = b"AAACGGTC";
        let quals: [u8; 8] = [10, 20, 30, 40, 50, 60, 70, 80];
        let work = make_work(0, vec![reverse_single_template(b"rev0", seq, &quals)]);

        let mut arena = ReadArena::new();
        fill_arena(&mut arena, &work);
        assert_eq!(arena.spans.len(), 1, "one primary read unpacked");
        let span = arena.spans[0];

        assert_eq!(arena.seq_of(&span), b"GACCGTTT", "SEQ is reverse-complemented");

        // QUAL is the Phred+33 encoding, reversed.
        let mut expected_qual: Vec<u8> = quals.iter().map(|&q| q + 33).collect();
        expected_qual.reverse();
        assert_eq!(arena.qual_of(&span), expected_qual.as_slice(), "QUAL is reversed (Phred+33)");
    }

    /// The arena holds exactly the SEQ/QUAL bytes the subprocess route writes to
    /// FASTQ, forward and `REVERSE`-flagged alike — including a non-IUPAC `=`
    /// base, which the FASTQ complement folds to `N`.
    #[rstest]
    #[case::forward(flags::UNMAPPED, b"AAACGGTC".as_slice())]
    #[case::reverse(flags::UNMAPPED | flags::REVERSE, b"AAACGGTC".as_slice())]
    #[case::reverse_with_eq_base(flags::UNMAPPED | flags::REVERSE, b"AAC=GGTN".as_slice())]
    fn arena_bytes_match_the_fastq_writer(#[case] record_flags: u16, #[case] seq: &[u8]) {
        let quals: Vec<u8> = (0u8..).step_by(5).take(seq.len()).collect();
        let mut b = SamBuilder::new();
        b.read_name(b"r0").flags(record_flags).sequence(seq).qualities(&quals);
        let record = b.build();

        let mut fastq = Vec::new();
        let mut buffers = crate::commands::fastq::FastqRecordBuffers::with_capacity(16);
        crate::commands::fastq::write_fastq_record(
            &mut fastq,
            &record,
            record_flags,
            true,
            &mut buffers,
            None,
        )
        .expect("write FASTQ");
        let lines: Vec<&[u8]> = fastq.split(|&c| c == b'\n').collect();

        let work = make_work(0, vec![Template::from_records(vec![record]).expect("template")]);
        let mut arena = ReadArena::new();
        fill_arena(&mut arena, &work);
        let span = arena.spans[0];
        assert_eq!(arena.seq_of(&span), lines[1], "SEQ matches the FASTQ writer");
        assert_eq!(arena.qual_of(&span), lines[3], "QUAL matches the FASTQ writer");
    }

    /// The output queue is sized `max(per_step, 2·T·EST_EXTENDED_WORK_BYTES)`.
    #[rstest]
    #[case::per_step_wins(64 * 1024 * 1024, 2, 64 * 1024 * 1024)]
    #[case::formula_wins(1024 * 1024, 8, 16 * 1024 * 1024)]
    #[case::threads_clamped_to_one(0, 0, 2 * 1024 * 1024)]
    fn output_byte_limit_is_the_max_of_per_step_and_formula(
        #[case] per_step: u64,
        #[case] threads: usize,
        #[case] expected: u64,
    ) {
        assert_eq!(seed_extend_output_byte_limit(per_step, threads), expected);
    }

    /// Profile: `Parallel`, single `ByteBounded` unordered output.
    #[test]
    fn profile_is_parallel_bytebounded_unordered() {
        let step = AlignSeedExtendStep::new(Arc::new(FakeEngine), Arc::default(), 4 * 1024 * 1024);
        let p = step.profile();
        assert_eq!(p.name, "AlignSeedExtend");
        assert_eq!(p.kind, StepKind::Parallel);
        assert!(!p.sticky);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::None]);
        let QueueSpec::ByteBounded { limit_bytes } = p.output_queues[0] else {
            panic!("output queue must be ByteBounded, got {:?}", p.output_queues[0]);
        };
        assert_eq!(
            limit_bytes,
            4 * 1024 * 1024,
            "output queue is sized at the constructor's output_byte_limit"
        );
    }

    /// `new_worker_copy` shares the engine `Arc` and the per-thread scratch
    /// pool, and gives each worker its own fresh arena.
    #[test]
    fn new_worker_copy_shares_engine_and_resets_per_worker_state() {
        let step = AlignSeedExtendStep::new(Arc::new(FakeEngine), Arc::default(), 4 * 1024 * 1024);
        let copy = step.new_worker_copy();
        assert!(Arc::ptr_eq(&step.engine, &copy.engine), "engine Arc is shared");
        assert!(
            Arc::ptr_eq(&step.scratch, &copy.scratch),
            "worker copies share the step's per-thread scratch pool"
        );
        assert_eq!(copy.arena.bytes.capacity(), 0, "worker copy gets a fresh arena");
        assert_eq!(copy.output_byte_limit, step.output_byte_limit);
    }
}
