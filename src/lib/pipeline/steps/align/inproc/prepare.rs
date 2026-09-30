//! `AlignPrepareStep` — the `Serial` front of the in-process bwa-mem3 backend.
//!
//! This step turns the queryname-grouped `BamTemplateBatch` stream (from
//! `GroupByQueryname`) into cohort-cut, sub-batch-sliced [`AlignWork`] items,
//! admitted through the [`CohortGate`] at **cohort granularity**. It is
//! engine-agnostic: it touches no FFI and no
//! [`AlignEngine`](super::engine::AlignEngine); `AlignWork` carries only the
//! moved-out unmapped templates plus the identity/position/layout the downstream
//! parallel steps need.
//!
//! What it does, per template (in input order), reproducing the subprocess
//! FASTQ-writer selection exactly so the two backends feed bwa byte-identical
//! reads:
//!
//! 1. Reject a read name containing whitespace in C `isspace`'s sense (space,
//!    `\t`, `\n`, `\v`, `\f`, `\r`). The subprocess path would truncate it at
//!    the first such byte — kseq splits the FASTQ name there — and then fail its
//!    reader's queryname-mismatch guard, so the in-process path rejects it up
//!    front rather than diverge; the byte-parity contract only covers valid
//!    unmapped-BAM names.
//!    Likewise reject a name ending in `/<digit>` (longer than two bytes): bwa's
//!    `trim_readno` strips that suffix on read, so the subprocess path's output
//!    names would not match the input and its reader fails on the mismatch.
//! 2. Select the template's **primary** records in record order (skip
//!    `SECONDARY | SUPPLEMENTARY` via the shared
//!    [`is_primary_for_alignment`]); a
//!    template with none is a hard error worded by the shared
//!    [`no_primary_records_message`].
//! 3. Feed each primary read's `l_seq` to the [`CohortCutter`] one read at a
//!    time, classify the template as a pair (2 primaries) or a single, and slice
//!    it into `sub_batch_templates`-bounded sub-batches that never span a cohort.
//!
//! **Mid-pair cut.** With mixed SE/PE input a `-K` cohort boundary can fall
//! *between* a pair's two reads (an odd running read count shifts parity). bwa
//! then classifies the two reads as two separate singles in adjacent cohorts.
//! The cutter reports this via its per-read
//! [`CohortCutter::push_read`], and this step splits such a pair into two
//! single-record templates, one placed in each cohort.
//!
//! **Cohort admission.** When a cohort opens the step reserves its whole byte
//! budget with a single non-blocking [`CohortGate::try_acquire_cohort`]; on
//! refusal it holds the popped input and reports `NoProgress` (never blocking a
//! worker). The reserved [`CohortLease`] is cloned into every `AlignWork` of the
//! cohort and released when the last one drains downstream. The final partial
//! cohort is closed at end-of-input, mirroring bwa closing its last `-K` chunk
//! at EOF.
//!
//! This step is the `Serial` front of the in-process backend wired by
//! [`InProcessBwaMem3Backend`](super::InProcessBwaMem3Backend).

use std::collections::VecDeque;
use std::io;
use std::sync::Arc;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::Single;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::align::{
    is_primary_for_alignment, no_primary_records_message, split_pair_into_singles,
};
use crate::pipeline::steps::types::BamTemplateBatch;
use crate::template::Template;

use super::FlushOutcome;
use super::cohort::{AlignWork, CohortCloser, CohortCutter, CohortPos, Layout, SubBatchId};
use super::gate::{CohortGate, CohortLease, cohort_bound_for_chunk_size, in_process_gate_budget};

/// One template (or one half of a mid-pair-split template) staged for placement
/// into a cohort's sub-batches. Computed by `template_to_placements`, which has
/// already driven the [`CohortCutter`] over the read(s); the `opens`/`closes`
/// flags record where cohort boundaries fall so placement needs no further
/// cutter interaction.
struct Placement {
    /// The template to move into a sub-batch. For a mid-pair split this is a
    /// single-record template holding one of the pair's reads.
    template: Template,
    /// Primary-read count for SE/PE classification: `2` = pair, else single. A
    /// split half is always `1`.
    count: u8,
    /// A fresh cohort must be opened (lease acquired) before placing this — true
    /// for the first read overall and for the read that opens a cohort after a
    /// cut.
    opens: bool,
    /// This placement closes its cohort (its last read triggered the `-K` cut).
    closes: bool,
}

/// The current sub-batch being accumulated within the open cohort.
#[derive(Default)]
struct SubBatchAccum {
    /// Templates moved into this sub-batch, in input order.
    templates: Vec<Template>,
    /// Per-template primary-read count (`1` or `2`), the input to
    /// [`Layout::classify`].
    counts: Vec<u8>,
    /// Single-end reads placed in this sub-batch so far.
    n_se: u32,
    /// Paired templates placed in this sub-batch so far.
    n_pe: u32,
    /// SE reads committed in earlier sub-batches of this cohort (this sub-batch's
    /// `CohortPos::se_offset`), captured when the sub-batch started.
    se_offset: u64,
    /// PE pairs committed in earlier sub-batches of this cohort (this sub-batch's
    /// `CohortPos::pe_offset`).
    pe_offset: u64,
}

/// The `Serial` step: cohort cutting, sub-batch slicing, dense serials, and
/// cohort-granularity gate admission.
pub(crate) struct AlignPrepareStep {
    /// Max templates per sub-batch; a sub-batch never exceeds this and never
    /// spans a cohort.
    sub_batch_templates: usize,
    /// Byte limit for the `AlignWork` output queue.
    output_byte_limit: u64,
    /// Cohort-granularity byte gate; the step reserves one cohort's `bound` on
    /// each cohort open.
    gate: Arc<CohortGate>,
    /// Per-cohort byte reservation passed to
    /// [`CohortGate::try_acquire_cohort`].
    cohort_bound: u64,

    /// Templates popped from input batches, not yet converted to placements.
    pending_input: VecDeque<Template>,
    /// Converted placements awaiting placement (front is next); a mid-pair split
    /// or a gate stall can leave items here across `try_run` calls.
    pending_placements: VecDeque<Placement>,

    /// Built `AlignWork` items awaiting push to the output.
    ready: VecDeque<AlignWork>,
    /// A single output item bounced by backpressure, retried before new pushes.
    held: HeldSlot<Unpushed<AlignWork>>,

    /// The `-K` even-parity cut state machine.
    cutter: CohortCutter,
    /// `true` when the next read opens a fresh cohort: initially, and after any
    /// cut.
    at_cohort_boundary: bool,

    /// Next dense serial to assign (dense from 0 across the whole run).
    next_serial: u64,
    /// Current cohort index.
    cohort: u32,
    /// Global reads before the current cohort — bwa's `n_processed` at cohort
    /// start.
    cohort_read_base: u64,
    /// The open cohort's lease; `None` between cohorts. Cloned into every
    /// `AlignWork` of the cohort.
    lease: Option<CohortLease>,
    /// Index of the next sub-batch within the current cohort.
    index_in_cohort: u32,
    /// SE reads committed in earlier (already-emitted) sub-batches of this cohort.
    cohort_se_committed: u64,
    /// PE pairs committed in earlier (already-emitted) sub-batches of this cohort.
    cohort_pe_committed: u64,
    /// The sub-batch currently being accumulated.
    sub: SubBatchAccum,

    /// Whether the final partial cohort has been closed at end-of-input.
    finalized: bool,
}

impl AlignPrepareStep {
    /// Build the step for a `-K` chunk size, sub-batch cap, and output byte
    /// limit, deriving the two-cohort in-flight gate from the chunk size.
    pub(crate) fn new(chunk_size: u64, sub_batch_templates: usize, output_byte_limit: u64) -> Self {
        let cohort_bound = cohort_bound_for_chunk_size(chunk_size);
        let gate = Arc::new(CohortGate::with_refill_signal(
            in_process_gate_budget(cohort_bound),
            cohort_bound,
        ));
        Self::with_gate(chunk_size, sub_batch_templates, output_byte_limit, gate, cohort_bound)
    }

    /// The cohort gate's pool-scheduler refill signal (see
    /// [`CohortGate::with_refill_signal`]).
    pub(crate) fn refill_signal(&self) -> Option<Arc<std::sync::atomic::AtomicBool>> {
        self.gate.refill_signal()
    }

    /// Build the step around a caller-provided gate + per-cohort bound. Lets
    /// tests size the gate independently of `chunk_size` to exercise admission
    /// stalls.
    fn with_gate(
        chunk_size: u64,
        sub_batch_templates: usize,
        output_byte_limit: u64,
        gate: Arc<CohortGate>,
        cohort_bound: u64,
    ) -> Self {
        assert!(sub_batch_templates >= 1, "sub_batch_templates must be >= 1");
        assert!(chunk_size >= 1, "chunk_size must be >= 1");
        Self {
            sub_batch_templates,
            output_byte_limit,
            gate,
            cohort_bound,
            pending_input: VecDeque::new(),
            pending_placements: VecDeque::new(),
            ready: VecDeque::new(),
            held: HeldSlot::new(),
            cutter: CohortCutter::new(chunk_size),
            at_cohort_boundary: true,
            next_serial: 0,
            cohort: 0,
            cohort_read_base: 0,
            lease: None,
            index_in_cohort: 0,
            cohort_se_committed: 0,
            cohort_pe_committed: 0,
            sub: SubBatchAccum::default(),
            finalized: false,
        }
    }

    /// Convert one input template into 1–2 [`Placement`]s (pushed onto
    /// `pending_placements`), driving the [`CohortCutter`] over its primary reads
    /// in record order. Splits a pair whose two reads straddle a cohort boundary
    /// into two singles.
    ///
    /// # Errors
    /// Whitespace read name, a name ending in `/<digit>`, no primary records, or
    /// (out of contract) more than two primary reads.
    fn template_to_placements(&mut self, template: Template) -> io::Result<()> {
        let name = template.name();
        // kseq's KS_SEP_SPACE splits on C `isspace`, which (unlike
        // `u8::is_ascii_whitespace`) also includes vertical tab (0x0B).
        if name.iter().any(|&b| b.is_ascii_whitespace() || b == 0x0B) {
            return Err(io::Error::other(format!(
                "align-and-merge (in-process): read name '{name}' contains ASCII whitespace, which \
                 the aligner truncates at the first whitespace (kseq splits the FASTQ name there). \
                 The in-process aligner rejects it rather than silently diverge from the subprocess \
                 path; byte-parity is contracted only for valid unmapped-BAM read names.",
                name = String::from_utf8_lossy(name),
            )));
        }
        if has_trimmed_readno_suffix(name) {
            return Err(io::Error::other(format!(
                "align-and-merge (in-process): read name '{name}' ends in '/<digit>', which the \
                 aligner strips from every read name (bwa's trim_readno), so its output names \
                 would not match the input. The in-process aligner rejects it rather than \
                 silently diverge from the subprocess path, which fails on the mismatched name; \
                 byte-parity is contracted only for valid unmapped-BAM read names.",
                name = String::from_utf8_lossy(name),
            )));
        }

        // Primary reads in record order, exactly the reads the subprocess writer
        // would emit to FASTQ. A template is single-end
        // (1) or paired-end (2); more is an out-of-contract error. So capture
        // only the first two lengths and the total count on the stack — no
        // per-template Vec — then classify on the count below.
        let mut primary_lens = [0u32; 2];
        let mut n_primary = 0usize;
        for record in template.records().iter().filter(|r| is_primary_for_alignment(r.flags())) {
            if n_primary < primary_lens.len() {
                primary_lens[n_primary] = record.l_seq();
            }
            n_primary += 1;
        }

        match n_primary {
            0 => Err(io::Error::other(format!(
                "align-and-merge (in-process): {}",
                no_primary_records_message(name)
            ))),
            1 => {
                let opens = self.at_cohort_boundary;
                let closes = self.cutter.push_read(primary_lens[0]);
                self.at_cohort_boundary = closes;
                self.pending_placements.push_back(Placement { template, count: 1, opens, closes });
                Ok(())
            }
            2 => {
                let opens0 = self.at_cohort_boundary;
                let cut0 = self.cutter.push_read(primary_lens[0]);
                if cut0 {
                    // Mid-pair cut: read 0 closes this cohort as a single; read 1
                    // opens the next cohort as a single. Split into two
                    // single-record templates.
                    let cut1 = self.cutter.push_read(primary_lens[1]);
                    // Read 1 is alone in the fresh cohort — an odd read count —
                    // so the even-parity cut rule cannot fire on it.
                    debug_assert!(!cut1, "a lone read cannot close a cohort");
                    self.at_cohort_boundary = false;
                    let (first, second) = split_pair_into_singles(template)?;
                    self.pending_placements.push_back(Placement {
                        template: first,
                        count: 1,
                        opens: opens0,
                        closes: true,
                    });
                    self.pending_placements.push_back(Placement {
                        template: second,
                        count: 1,
                        opens: true,
                        closes: false,
                    });
                } else {
                    let closes = self.cutter.push_read(primary_lens[1]);
                    self.at_cohort_boundary = closes;
                    self.pending_placements.push_back(Placement {
                        template,
                        count: 2,
                        opens: opens0,
                        closes,
                    });
                }
                Ok(())
            }
            n => Err(io::Error::other(format!(
                "align-and-merge (in-process): template '{name}' has {n} primary records; the \
                 in-process aligner supports single-end (1) or paired-end (2) unmapped templates.",
                name = String::from_utf8_lossy(name),
            ))),
        }
    }

    /// Process staged input into `ready` until it is exhausted or the gate
    /// refuses a cohort open. Returns `true` iff it stalled on the gate (the
    /// offending placement stays at the front of `pending_placements`).
    ///
    /// # Errors
    /// Propagates a `template_to_placements` error.
    fn advance(&mut self) -> io::Result<bool> {
        loop {
            if self.pending_placements.is_empty() {
                let Some(template) = self.pending_input.pop_front() else {
                    return Ok(false); // nothing staged
                };
                self.template_to_placements(template)?;
                continue;
            }

            // Peek the front placement: if it opens a cohort, admit it first.
            let opens = self.pending_placements.front().expect("non-empty checked above").opens;
            if opens && self.lease.is_none() {
                match self.gate.try_acquire_cohort(self.cohort_bound) {
                    Some(lease) => self.begin_cohort(lease),
                    None => return Ok(true), // stalled: leave the placement queued
                }
            }

            let placement = self.pending_placements.pop_front().expect("non-empty checked above");
            self.place(placement);
        }
    }

    /// Open a fresh cohort with the acquired `lease`, starting its first
    /// sub-batch.
    fn begin_cohort(&mut self, lease: CohortLease) {
        debug_assert!(self.lease.is_none(), "opening a cohort while one is already open");
        self.lease = Some(lease);
        self.index_in_cohort = 0;
        self.cohort_se_committed = 0;
        self.cohort_pe_committed = 0;
        self.start_sub_batch();
    }

    /// Reset the sub-batch accumulator in place, capturing the cohort offsets at
    /// its start.
    ///
    /// `templates` was just handed downstream by [`Self::emit_sub_batch`]'s
    /// `mem::take` (or is the empty initial buffer), so it always starts fresh —
    /// pre-size it to the sub-batch cap so it never reallocates while filling.
    /// `counts` is only borrowed by [`Layout::classify`] and then reset here, so
    /// its buffer is pooled: `clear()` retains its capacity across sub-batches.
    fn start_sub_batch(&mut self) {
        self.sub.templates = Vec::with_capacity(self.sub_batch_templates);
        self.sub.counts.clear();
        self.sub.n_se = 0;
        self.sub.n_pe = 0;
        self.sub.se_offset = self.cohort_se_committed;
        self.sub.pe_offset = self.cohort_pe_committed;
    }

    /// Place one placement into the open cohort: flush the current sub-batch
    /// first if adding this template would exceed the cap, add the template, then
    /// close the cohort if this placement carries the cut.
    fn place(&mut self, placement: Placement) {
        debug_assert!(self.lease.is_some(), "placing a template with no cohort open");

        // Cap the sub-batch: flush the (full) current one as non-last before
        // adding, so a sub-batch never exceeds the cap and the flushed one is
        // known not to be the cohort's last.
        if self.sub.templates.len() >= self.sub_batch_templates {
            self.emit_sub_batch(false);
            self.index_in_cohort += 1;
            self.start_sub_batch();
        }

        self.sub.templates.push(placement.template);
        self.sub.counts.push(placement.count);
        if placement.count == 2 {
            self.sub.n_pe += 1;
        } else {
            self.sub.n_se += 1;
        }

        if placement.closes {
            self.emit_sub_batch(true);
            self.close_cohort();
        }
    }

    /// Emit the current sub-batch as an [`AlignWork`], stamping a
    /// [`CohortCloser`] iff it is the cohort's last, and commit its SE/PE counts
    /// to the cohort running totals. Does not advance `index_in_cohort` or start
    /// a new sub-batch — the caller owns those transitions.
    // `n_se`/`n_pe` are the domain terms (they are also the `CohortPos`/`CohortCloser`
    // field names); the pair reads more clearly together than renamed apart.
    #[allow(clippy::similar_names)]
    fn emit_sub_batch(&mut self, is_last: bool) {
        let n_se = self.sub.n_se;
        let n_pe = self.sub.n_pe;
        let closer = is_last.then(|| CohortCloser {
            n_sub_batches: self.index_in_cohort + 1,
            cohort_n_se: self.cohort_se_committed + u64::from(n_se),
            cohort_n_pe: self.cohort_pe_committed + u64::from(n_pe),
        });
        let layout = Layout::classify(&self.sub.counts);
        let pos = CohortPos {
            cohort_read_base: self.cohort_read_base,
            se_offset: self.sub.se_offset,
            pe_offset: self.sub.pe_offset,
            n_se,
            n_pe,
        };
        let id = SubBatchId {
            serial: self.next_serial,
            cohort: self.cohort,
            index_in_cohort: self.index_in_cohort,
        };
        let unmapped = std::mem::take(&mut self.sub.templates);
        let lease = self.lease.clone().expect("cohort lease held while emitting a sub-batch");
        self.ready.push_back(AlignWork { id, pos, closer, layout, unmapped, lease });

        self.next_serial += 1;
        self.cohort_se_committed += u64::from(n_se);
        self.cohort_pe_committed += u64::from(n_pe);
    }

    /// Close the open cohort: advance the global read base by this cohort's total
    /// reads, bump the cohort index, and drop the step's lease clone (downstream
    /// clones release the gate reservation as they drain).
    fn close_cohort(&mut self) {
        self.cohort_read_base += self.cohort_se_committed + 2 * self.cohort_pe_committed;
        self.cohort += 1;
        self.lease = None;
        self.gate.cohort_filled();
    }

    /// Close the final partial cohort at end-of-input (bwa closes its last `-K`
    /// chunk at EOF too). The current sub-batch is the cohort's last.
    fn finalize(&mut self) {
        if self.lease.is_some() {
            debug_assert!(
                !self.sub.templates.is_empty(),
                "an open cohort always has a non-empty current sub-batch"
            );
            self.emit_sub_batch(true);
            self.close_cohort();
        }
    }

    /// Push staged `ready` items (and any held-back one) to the output.
    fn flush(&mut self, ctx: &mut StepCtx<'_, Self>) -> FlushOutcome {
        if let Some(unpushed) = self.held.take()
            && let Err(again) = ctx.outputs.retry(unpushed)
        {
            self.held.put(again);
            return FlushOutcome::StillHeld;
        }
        while let Some(work) = self.ready.pop_front() {
            if let Err(unpushed) = ctx.outputs.push(work) {
                self.held.put(unpushed);
                return FlushOutcome::NewlyHeld;
            }
        }
        FlushOutcome::Clear
    }
}

/// Whether bwa would strip a read-number suffix from `name`: an exact port of
/// bwa-mem3's `trim_readno` / `fr_trim_readno_len` (`bwa.cpp`,
/// `fast_reader_bseq.c`), which drop a trailing `/` + one ASCII digit from any
/// name longer than two bytes.
fn has_trimmed_readno_suffix(name: &[u8]) -> bool {
    let l = name.len();
    l > 2 && name[l - 2] == b'/' && name[l - 1].is_ascii_digit()
}

impl Step for AlignPrepareStep {
    type Input = BamTemplateBatch;
    type Outputs = Single<AlignWork>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "AlignPrepare",
            kind: StepKind::Serial,
            sticky: false,
            // `AlignWork` carries no cross-batch ordinal that a consumer must see
            // in order (the merge restores order via each `BamTemplateBatch`'s
            // serial), so the output queue is byte-bounded, unordered.
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        loop {
            // 1. Flush pending output first. A newly held item changed state
            //    (`Progress`); re-failing an already-held one did not
            //    (`Contention`).
            if let Some(outcome) = self.flush(ctx).blocked_outcome() {
                return Ok(outcome);
            }

            // 2. Convert staged input into `AlignWork`.
            if !self.pending_placements.is_empty() || !self.pending_input.is_empty() {
                let stalled = self.advance()?;
                if !self.ready.is_empty() {
                    continue; // flush what we produced
                }
                if stalled {
                    // Gate full: hold the staged input and wait for a downstream
                    // release. Do not drop data.
                    return Ok(StepOutcome::NoProgress);
                }
                continue; // input absorbed into the open sub-batch; pull more
            }

            // 3. Pull the next input batch, moving its templates out wholesale.
            let Some(batch) = ctx.input.pop() else {
                if ctx.input.is_drained() {
                    if !self.finalized {
                        self.finalize();
                        self.finalized = true;
                        continue; // flush the final cohort's last sub-batch
                    }
                    return Ok(StepOutcome::Finished);
                }
                return Ok(StepOutcome::NoProgress);
            };
            let (_serial, templates) = batch.into_parts();
            self.pending_input.extend(templates);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::{SamBuilder, flags};
    use rstest::rstest;

    // Distinct read names (never 2 bytes: `ci-tag-literals` rejects bare 2-byte
    // byte-string literals).
    fn seq_of(len: u32) -> Vec<u8> {
        b"ACGTN".iter().copied().cycle().take(len as usize).collect()
    }

    /// A single-end template with one primary record of `len` bases.
    fn single_template(name: &[u8], len: u32) -> Template {
        let seq = seq_of(len);
        let mut b = SamBuilder::new();
        b.read_name(name)
            .flags(flags::UNMAPPED)
            .sequence(&seq)
            .qualities(&vec![30u8; len as usize]);
        Template::from_records(vec![b.build()]).expect("single template")
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

    /// A template whose single record carries only the SECONDARY flag.
    fn secondary_only_template(name: &[u8]) -> Template {
        let mut b = SamBuilder::new();
        b.read_name(name).flags(flags::SECONDARY).sequence(b"ACGT").qualities(&[30u8; 4]);
        // `from_records` categorizes by flags; a lone secondary yields a template
        // with zero primaries, which is exactly the fixture we want.
        Template::from_records(vec![b.build()]).expect("secondary-only template")
    }

    /// A step with an effectively unbounded gate (so admission never stalls) for
    /// the correctness tests. Distinct-name pairs cannot collide.
    fn step_unbounded(chunk_size: u64, sub_batch_templates: usize) -> AlignPrepareStep {
        let gate = Arc::new(CohortGate::new(u64::MAX));
        AlignPrepareStep::with_gate(chunk_size, sub_batch_templates, 1 << 20, gate, 1)
    }

    /// Drive a set of input batches through the step's internal machinery,
    /// draining every produced `AlignWork` (so leases release and the gate never
    /// fills), then finalize. Returns all produced work in emission order.
    fn drive(step: &mut AlignPrepareStep, batches: Vec<Vec<Template>>) -> Vec<AlignWork> {
        let mut out = Vec::new();
        for templates in batches {
            step.pending_input.extend(templates);
            loop {
                let stalled = step.advance().expect("advance");
                out.extend(step.ready.drain(..));
                assert!(!stalled, "drive() uses an unbounded gate; it must never stall");
                if step.pending_input.is_empty() && step.pending_placements.is_empty() {
                    break;
                }
            }
        }
        step.finalize();
        out.extend(step.ready.drain(..));
        out
    }

    // ─────────────────────────────────────────────────────────────────────
    // (a) FASTQ-writer selection semantics
    // ─────────────────────────────────────────────────────────────────────

    #[test]
    fn all_secondary_template_is_a_hard_error() {
        let mut step = step_unbounded(1_000_000, 8);
        step.pending_input.push_back(secondary_only_template(b"secondary_only"));
        let err = step.advance().expect_err("all-secondary template must be a hard error");
        let msg = err.to_string();
        assert!(msg.contains("no primary records"), "message: {msg}");
        assert!(msg.contains("secondary_only"), "message names the queryname: {msg}");
    }

    #[rstest]
    #[case::space(b"read one".as_slice())]
    #[case::tab(b"read\tone".as_slice())]
    #[case::newline(b"read\none".as_slice())]
    #[case::vertical_tab(b"read\x0Bone".as_slice())]
    #[case::carriage_return(b"read\rone".as_slice())]
    #[case::form_feed(b"read\x0Cone".as_slice())]
    fn whitespace_read_name_is_rejected(#[case] name: &[u8]) {
        let mut step = step_unbounded(1_000_000, 8);
        step.pending_input.push_back(single_template(name, 10));
        let err = step.advance().expect_err("a whitespace read name must be rejected");
        assert!(err.to_string().contains("whitespace"), "message: {err}");
    }

    /// bwa's `trim_readno` strips `/` + one digit from a name longer than two
    /// bytes, so those names are rejected; anything the rule leaves alone is
    /// accepted.
    #[rstest]
    #[case::slash_one(b"read/1".as_slice(), true)]
    #[case::slash_two(b"read/2".as_slice(), true)]
    #[case::slash_other_digit(b"read/3".as_slice(), true)]
    #[case::slash_digit_mid_name(b"a/1b".as_slice(), false)]
    #[case::two_digits_after_slash(b"read/12".as_slice(), false)]
    #[case::too_short_to_trim(b"/1".as_slice(), false)]
    #[case::slash_letter(b"read/a".as_slice(), false)]
    fn readno_suffix_name_is_rejected(#[case] name: &[u8], #[case] rejected: bool) {
        let mut step = step_unbounded(1_000_000, 8);
        step.pending_input.push_back(single_template(name, 10));
        let result = step.advance();
        if rejected {
            let err = result.expect_err("a '/<digit>' read name must be rejected");
            assert!(err.to_string().contains("ends in '/<digit>'"), "message: {err}");
        } else {
            result.expect("a name bwa does not trim is accepted");
        }
    }

    // ─────────────────────────────────────────────────────────────────────
    // (b) dense serials from 0 across multiple input batches
    // ─────────────────────────────────────────────────────────────────────

    #[test]
    fn serials_are_dense_from_zero_across_batches() {
        // Large chunk → one cohort; small sub-batch cap → many sub-batches that
        // straddle the two input batches.
        let mut step = step_unbounded(1_000_000, 2);
        let batch_a: Vec<Template> =
            (0..5).map(|i| pair_template(format!("a{i}").as_bytes(), 10)).collect();
        let batch_b: Vec<Template> =
            (0..4).map(|i| pair_template(format!("b{i}").as_bytes(), 10)).collect();
        let work = drive(&mut step, vec![batch_a, batch_b]);

        let serials: Vec<u64> = work.iter().map(|w| w.id.serial).collect();
        assert_eq!(serials, (0..serials.len() as u64).collect::<Vec<_>>(), "dense from 0, no gaps");
        // 9 pairs, cap 2 → ceil(9/2) = 5 sub-batches, all one cohort.
        assert_eq!(work.len(), 5);
        assert!(work.iter().all(|w| w.id.cohort == 0), "one cohort for a large -K");
    }

    // ─────────────────────────────────────────────────────────────────────
    // (c) sub-batches never span a cohort and are <= sub_batch_templates
    // ─────────────────────────────────────────────────────────────────────

    #[test]
    fn sub_batches_respect_cap_and_never_span_a_cohort() {
        let cap = 3;
        let mut step = step_unbounded(1_000_000, cap);
        let batch: Vec<Template> =
            (0..10).map(|i| pair_template(format!("t{i}").as_bytes(), 10)).collect();
        let work = drive(&mut step, vec![batch]);

        // 10 pairs, cap 3 → sub-batches of 3,3,3,1.
        assert_eq!(work.iter().map(|w| w.unmapped.len()).collect::<Vec<_>>(), vec![3, 3, 3, 1]);
        assert!(work.iter().all(|w| w.unmapped.len() <= cap), "no sub-batch exceeds the cap");
        assert!(work.iter().all(|w| w.id.cohort == 0), "all in one cohort");
        // index_in_cohort is dense 0..n within the cohort.
        assert_eq!(work.iter().map(|w| w.id.index_in_cohort).collect::<Vec<_>>(), vec![0, 1, 2, 3]);
    }

    // ─────────────────────────────────────────────────────────────────────
    // (d) CohortCloser on exactly the last sub-batch, with correct counts
    // ─────────────────────────────────────────────────────────────────────

    #[test]
    fn closer_stamped_on_last_sub_batch_only_with_correct_counts() {
        let mut step = step_unbounded(1_000_000, 3);
        let batch: Vec<Template> =
            (0..10).map(|i| pair_template(format!("p{i}").as_bytes(), 10)).collect();
        let work = drive(&mut step, vec![batch]);

        let with_closer: Vec<&AlignWork> = work.iter().filter(|w| w.closer.is_some()).collect();
        assert_eq!(with_closer.len(), 1, "exactly one sub-batch carries the closer");
        let last = with_closer[0];
        assert_eq!(last.id.index_in_cohort, 3, "closer is on the cohort's last sub-batch");
        let closer = last.closer.expect("closer present");
        assert_eq!(closer.n_sub_batches, 4);
        assert_eq!(closer.cohort_n_pe, 10, "all 10 templates are pairs");
        assert_eq!(closer.cohort_n_se, 0);
    }

    #[test]
    fn multiple_cohorts_each_get_a_closer_and_advancing_read_base() {
        // chunk 20, reads of 10 bases → each all-pairs template (2 reads = 20
        // bases, even parity) is its own cohort.
        let mut step = step_unbounded(20, 8);
        let batch: Vec<Template> =
            (0..3).map(|i| pair_template(format!("c{i}").as_bytes(), 10)).collect();
        let work = drive(&mut step, vec![batch]);

        assert_eq!(work.len(), 3, "one sub-batch per single-pair cohort");
        assert_eq!(work.iter().map(|w| w.id.cohort).collect::<Vec<_>>(), vec![0, 1, 2]);
        // Every cohort is complete → every sub-batch is its own cohort's last.
        assert!(work.iter().all(|w| w.closer.is_some()));
        // cohort_read_base advances by 2 reads per (single-pair) cohort.
        assert_eq!(work.iter().map(|w| w.pos.cohort_read_base).collect::<Vec<_>>(), vec![0, 2, 4],);
        for w in &work {
            let closer = w.closer.expect("closer");
            assert_eq!((closer.cohort_n_pe, closer.cohort_n_se, closer.n_sub_batches), (1, 0, 1));
        }
    }

    // ─────────────────────────────────────────────────────────────────────
    // Mid-pair split (reviewer-requested)
    // ─────────────────────────────────────────────────────────────────────

    #[test]
    fn mid_pair_cut_splits_a_pair_into_two_singles_across_cohorts() {
        // chunk 10, reads of 10 bases. Sequence: one single, then a pair.
        //   single S: read → size 10, n=1 (odd) → no cut.
        //   pair P read0: size 20, n=2 (even) ≥ 10 → CUT between P's two reads.
        // So P's read0 closes cohort 0 (as a single, with S); read1 opens
        // cohort 1 (as a single).
        let mut step = step_unbounded(10, 8);
        let work =
            drive(&mut step, vec![vec![single_template(b"solo", 10), pair_template(b"duo", 10)]]);

        assert_eq!(work.len(), 2, "two cohorts: {{solo, duo/1}} then {{duo/2}}");

        let c0 = &work[0];
        assert_eq!(c0.id.cohort, 0);
        let c0_names: Vec<&[u8]> = c0.unmapped.iter().map(Template::name).collect();
        assert_eq!(
            c0_names,
            vec![b"solo".as_slice(), b"duo".as_slice()],
            "cohort 0 holds solo + duo's first read"
        );
        // Both are singles: layout is Mixed, and the closer counts two singles.
        assert!(matches!(c0.layout, Layout::Mixed(_)), "cohort 0 is not all-pairs");
        let c0_closer = c0.closer.expect("cohort 0 closes at the mid-pair cut");
        assert_eq!((c0_closer.cohort_n_se, c0_closer.cohort_n_pe), (2, 0));

        let c1 = &work[1];
        assert_eq!(c1.id.cohort, 1);
        assert_eq!(
            c1.unmapped.iter().map(Template::name).collect::<Vec<_>>(),
            vec![b"duo".as_slice()]
        );
        assert_eq!(c1.unmapped.len(), 1, "duo's second read is a lone single in cohort 1");
        let c1_closer = c1.closer.expect("cohort 1 closed at finalize");
        assert_eq!((c1_closer.cohort_n_se, c1_closer.cohort_n_pe), (1, 0));

        // The split preserved both of the pair's records (one per half).
        assert_eq!(records_named(c0, b"duo") + records_named(c1, b"duo"), 2);
    }

    /// Count records across a sub-batch's templates whose queryname is `name`.
    fn records_named(work: &AlignWork, name: &[u8]) -> usize {
        work.unmapped.iter().filter(|t| t.name() == name).map(|t| t.records().len()).sum()
    }

    // ─────────────────────────────────────────────────────────────────────
    // (e) cohort-granularity gate admission at the step level
    // ─────────────────────────────────────────────────────────────────────

    #[test]
    fn advance_stalls_when_the_gate_is_full_then_resumes_after_release() {
        // Two-cohort budget; chunk 20 / 10-base reads → one cohort per pair.
        let bound = 100;
        let gate = Arc::new(CohortGate::new(2 * bound));
        let mut step = AlignPrepareStep::with_gate(20, 8, 1 << 20, Arc::clone(&gate), bound);
        // Three single-pair cohorts staged.
        step.pending_input.extend((0..3).map(|i| pair_template(format!("g{i}").as_bytes(), 10)));

        // First advance admits two cohorts, then stalls opening the third.
        let stalled = step.advance().expect("advance");
        assert!(stalled, "the third cohort must stall the two-cohort gate");
        let first_two: Vec<AlignWork> = step.ready.drain(..).collect();
        assert_eq!(first_two.len(), 2, "two cohorts admitted before the stall");
        assert_eq!(gate.in_flight_bytes(), 2 * bound, "both leases still reserve the gate");

        // Dropping the drained work releases their leases; the gate frees.
        drop(first_two);
        assert_eq!(gate.in_flight_bytes(), 0, "released after the admitted cohorts drop");

        // Re-advancing now admits the held third cohort.
        let stalled = step.advance().expect("advance");
        assert!(!stalled, "the third cohort admits once the gate freed");
        let third: Vec<AlignWork> = step.ready.drain(..).collect();
        assert_eq!(third.len(), 1);
        assert_eq!(third[0].id.cohort, 2, "serial cohort numbering continues");
    }

    /// A gate stall landing *between* the two halves of a mid-pair split keeps
    /// the second half queued as a cohort opener and admits it once the gate
    /// frees. chunk 10 / 10-base reads: `solo` + `duo`'s first read close cohort
    /// 0; `duo`'s second read must open cohort 1, which a one-cohort gate refuses
    /// while cohort 0's work is still alive.
    #[test]
    fn gate_stall_between_mid_pair_halves_keeps_second_half_queued() {
        let bound = 100;
        let gate = Arc::new(CohortGate::new(bound));
        let mut step = AlignPrepareStep::with_gate(10, 8, 1 << 20, Arc::clone(&gate), bound);
        step.pending_input.extend([single_template(b"solo", 10), pair_template(b"duo", 10)]);

        let stalled = step.advance().expect("advance");
        assert!(stalled, "the second half's cohort must stall the one-cohort gate");
        let cohort0: Vec<AlignWork> = step.ready.drain(..).collect();
        assert_eq!(cohort0.len(), 1, "cohort 0 (solo + duo's first read) was emitted");
        assert_eq!(records_named(&cohort0[0], b"duo"), 1, "cohort 0 holds one duo read");
        assert_eq!(step.pending_placements.len(), 1, "the second half stays queued");
        let queued = step.pending_placements.front().expect("queued second half");
        assert!(queued.opens, "the queued second half still opens a cohort");
        assert_eq!(queued.count, 1, "the queued second half is a single");

        drop(cohort0);
        let stalled = step.advance().expect("advance");
        assert!(!stalled, "the second half admits once cohort 0 released the gate");
        step.finalize();
        let cohort1: Vec<AlignWork> = step.ready.drain(..).collect();
        assert_eq!(cohort1.len(), 1);
        assert_eq!(cohort1[0].id.cohort, 1, "the second half opens cohort 1");
        assert_eq!(records_named(&cohort1[0], b"duo"), 1, "cohort 1 holds duo's second read");
    }

    // ─────────────────────────────────────────────────────────────────────
    // profile
    // ─────────────────────────────────────────────────────────────────────

    #[test]
    fn profile_is_serial_bytebounded_unordered() {
        let step = step_unbounded(1_000_000, 8);
        let p = step.profile();
        assert_eq!(p.name, "AlignPrepare");
        assert_eq!(p.kind, StepKind::Serial);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::None]);
        assert!(matches!(p.output_queues[0], QueueSpec::ByteBounded { .. }));
    }

    /// A single-end-only cohort classifies its templates as singles and closes at
    /// finalize with the right SE count. Guards that the SE path (not just PE)
    /// carries through layout, counts, and the closer.
    #[test]
    fn single_end_only_cohort_counts_singles() {
        let mut step = step_unbounded(1_000_000, 8);
        let work = drive(
            &mut step,
            vec![(0..4).map(|i| single_template(format!("s{i}").as_bytes(), 10)).collect()],
        );
        assert_eq!(work.len(), 1);
        let w = &work[0];
        assert!(matches!(w.layout, Layout::Mixed(_)), "all-singles is Mixed, not AllPairs");
        assert_eq!((w.pos.n_se, w.pos.n_pe), (4, 0));
        let closer = w.closer.expect("closer at finalize");
        assert_eq!((closer.cohort_n_se, closer.cohort_n_pe), (4, 0));
    }
}
