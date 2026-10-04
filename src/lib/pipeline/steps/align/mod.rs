//! Align-stage backend abstraction and shared pieces.
//!
//! `runall`'s align stage pairs `N` aligner-emitted (mapped) templates with the
//! `N` original (unmapped) templates of each batch (`ZipperBatch`, tagged with
//! a dense serial from 0) and zipper-merges them into `BamTemplateBatch`es with
//! `merge::merge_zipper_batch`.
//!
//! What produces the aligned half is a *backend*, behind the `AlignBackend`
//! trait: the subprocess backend (a `Serial` `subprocess::SubprocessAlignStep`
//! wrapping an aligner subprocess, followed by the `Parallel`
//! `merge::MergeAlignedStep` that merges its `ZipperBatch` stream), or the
//! in-process bwa-mem3 subgraph, whose pair/emit step merges each sub-batch
//! itself as soon as it is emitted. The trait's `wire` appends the backend's
//! steps after a queryname-grouped `BamTemplateBatch` tail and returns its
//! merged `BamTemplateBatch` tail plus the scheduling facts the chain builder
//! needs (`AlignWired`). The merge itself is common code. `backend_for`
//! builds the backend a resolved `--aligner::*` selection names.
//!
//! This module holds the trait, its wiring context/result types, the
//! backend-agnostic `ZipperBatch`, the subprocess backend's `InFlightGate`
//! byte-budget gate, the mimalloc purge-delay setter both backends call at
//! wiring, and the header helpers (`validate_sq_consistency`,
//! `merge_aligner_header`) shared by both backends.

#[cfg(feature = "aligner-bwa-mem3")]
pub(crate) mod inproc;
// `bwa-mem3-rs` is a target-gated dependency, so on any other architecture the
// feature would enable `inproc` without the crate it links; reject it clearly.
#[cfg(all(
    feature = "aligner-bwa-mem3",
    not(any(target_arch = "x86_64", target_arch = "aarch64"))
))]
compile_error!("feature `aligner-bwa-mem3` is supported only on x86_64 and aarch64");
pub(crate) mod merge;
pub(crate) mod subprocess;

use std::io;
use std::sync::Arc;

use noodles::sam::Header;
use parking_lot::Mutex;

use crate::aligner::{ResolvedAligner, ResolvedBackend};
use crate::pipeline::core::builder::PipelineBuilder;
use crate::pipeline::core::header::HeaderHandle;
use crate::pipeline::core::item::{HeapSize, Ordered};
use crate::pipeline::core::topology::{BranchIdx, StepIdx};
use crate::pipeline::steps::types::{BamTemplateBatch, container_bytes};
use crate::template::Template;

/// Flag bits excluded when selecting a template's primary reads to align.
/// `SECONDARY | SUPPLEMENTARY` (`0x900`): if a previously-aligned BAM is passed
/// by mistake, dropping these prevents duplicate reads for the same template.
///
/// This is the single source of truth shared by BOTH align backends — the
/// subprocess FASTQ writer ([`subprocess`]) and the in-process prepare step
/// (`inproc::prepare`) — so the two select byte-identical reads and cannot
/// drift.
pub(crate) const FASTQ_WRITER_EXCLUDE_FLAGS: u16 =
    fgumi_raw_bam::flags::SECONDARY | fgumi_raw_bam::flags::SUPPLEMENTARY;

/// Whether a record with these `flags` is a *primary* read to align — i.e. not
/// SECONDARY/SUPPLEMENTARY. Both backends walk a template's records in record
/// order and keep exactly the records this returns `true` for, so the reads fed
/// to the aligner (and the `-K` byte accounting over them) are identical.
#[must_use]
pub(crate) fn is_primary_for_alignment(flags: u16) -> bool {
    (flags & FASTQ_WRITER_EXCLUDE_FLAGS) == 0
}

/// The shared body of the "template has no primary records" hard error, so the
/// subprocess writer and the in-process prepare step word it identically and
/// cannot drift. A template whose every record is SECONDARY/SUPPLEMENTARY
/// produces zero aligner input; align-and-merge expects unmapped BAM input
/// (from `fgumi extract`), which has no secondaries, so this means a re-aligned
/// BAM was passed by mistake (the subprocess writer's `!wrote_any` guard in
/// `writer_loop`).
#[must_use]
pub(crate) fn no_primary_records_message(name: &[u8]) -> String {
    format!(
        "template '{name}' has no primary records (all flagged SECONDARY/SUPPLEMENTARY). \
         align-and-merge expects unmapped BAM input — did you pass a re-aligned BAM by mistake? \
         Run `fgumi extract` first or pre-filter the input.",
        name = String::from_utf8_lossy(name),
    )
}

/// Estimated in-flight bytes per aligner-input base, used to derive the
/// in-flight-unmapped byte budget from the aligner's `-K` chunk size.
/// One base costs ~1 byte of sequence + ~1 byte of quality held in the
/// unmapped `RawRecord`, plus read-name / tag overhead — 3 is a
/// deliberately generous estimate so the budget never under-shoots the
/// aligner's real buffering (under-shooting risks a writer stall).
const IN_FLIGHT_BYTES_PER_BASE: u64 = 3;

/// Multiplier over a single `-K` chunk for the in-flight budget. bwa's
/// `kt_pipeline` double-buffers (`p_nt=2`: one chunk aligning, one
/// queued), so the aligner can hold ~2 chunks before emitting; 4× adds
/// slack for the FASTQ-write buffer and pipeline jitter so the writer
/// never blocks before the aligner has a full chunk to emit (which
/// would deadlock — see `InFlightGate`).
const IN_FLIGHT_CHUNK_MULTIPLIER: u64 = 4;

/// Floor for the in-flight-unmapped budget. Keeps the budget workable
/// for small `-K` values (custom aligner commands, tests) where the
/// `-K`-derived value would be tiny, and covers aligners whose internal
/// buffering exceeds their nominal `-K`.
const IN_FLIGHT_MIN_BUDGET: u64 = 512 * 1024 * 1024;

/// Derive the in-flight-unmapped byte budget for the subprocess align path from
/// the aligner's `-K` chunk size (bases per batch). The budget bounds the
/// otherwise unbounded writer→reader token backlog: the writer blocks once the
/// in-flight unmapped bytes reach this, so a fast-draining aligner can't
/// accumulate the whole input's unmapped reads in RAM (issue #382).
///
/// Sized `≥` the aligner's real buffering so a **streaming** aligner
/// (one that emits output after ≤ `-K` bases — every real aligner)
/// always has a full chunk to emit before the writer blocks, so the
/// block is transient, not a deadlock. A non-streaming aligner (one that
/// reads all stdin before emitting) is out of contract: it would stall
/// the writer here rather than OOM. For the default `-K` of 150M bases
/// this is ~1.8 GiB.
#[must_use]
pub(crate) fn in_flight_budget_for_chunk_size(chunk_size_bases: u64) -> u64 {
    chunk_size_bases
        .saturating_mul(IN_FLIGHT_BYTES_PER_BASE)
        .saturating_mul(IN_FLIGHT_CHUNK_MULTIPLIER)
        .max(IN_FLIGHT_MIN_BUDGET)
}

/// Byte-budget gate bounding the subprocess align path's in-flight unmapped
/// reads (the writer→reader token backlog). The writer
/// [`acquire`](InFlightGate::acquire)s before feeding a batch to the aligner and
/// the reader [`release`](InFlightGate::release)s when it consumes the matching
/// token, so resident in-flight unmapped bytes stay near `budget`.
///
/// Deadlock-safety: the gate blocks the writer only while in-flight is at
/// budget AND non-empty — a single oversized batch always passes when the
/// gate is empty. Because the budget is sized `≥` a streaming aligner's
/// `-K` buffering (see [`in_flight_budget_for_chunk_size`]), the aligner
/// always has a full chunk to emit before the writer blocks, so the
/// reader drains and the block lifts. If the reader exits (EOF or error),
/// it latches `consumer_gone` so a blocked writer wakes and bails rather
/// than hanging.
pub(crate) struct InFlightGate {
    budget: u64,
    inner: Mutex<GateInner>,
    cond: parking_lot::Condvar,
}

struct GateInner {
    in_flight: u64,
    consumer_gone: bool,
}

impl InFlightGate {
    pub(crate) fn new(budget: u64) -> Self {
        Self {
            budget,
            inner: Mutex::new(GateInner { in_flight: 0, consumer_gone: false }),
            cond: parking_lot::Condvar::new(),
        }
    }

    /// Reserve `n` in-flight bytes, blocking while the gate is full and
    /// non-empty. Returns `false` if the consumer (reader) has gone — the
    /// caller (writer) should then stop.
    pub(crate) fn acquire(&self, n: u64) -> bool {
        let mut g = self.inner.lock();
        // Block while non-empty AND this reservation would exceed budget.
        // The `in_flight != 0` guard lets a single oversized batch through
        // when the gate is empty (it can't be split, so holding it is
        // unavoidable) — without it the writer would deadlock on a batch
        // larger than the whole budget.
        while !g.consumer_gone && g.in_flight != 0 && g.in_flight.saturating_add(n) > self.budget {
            self.cond.wait(&mut g);
        }
        if g.consumer_gone {
            return false;
        }
        g.in_flight = g.in_flight.saturating_add(n);
        true
    }

    /// Release `n` in-flight bytes (reader consumed a token) and wake the
    /// writer if it is blocked.
    pub(crate) fn release(&self, n: u64) {
        let mut g = self.inner.lock();
        g.in_flight = g.in_flight.saturating_sub(n);
        drop(g);
        self.cond.notify_all();
    }

    /// Latch that the consumer (reader) has exited, so a blocked writer
    /// wakes and bails instead of hanging forever.
    pub(crate) fn mark_consumer_gone(&self) {
        let mut g = self.inner.lock();
        g.consumer_gone = true;
        drop(g);
        self.cond.notify_all();
    }
}

/// RAII guard that latches consumer-gone on drop — on normal return, on an
/// error return, AND on a panic unwinding out of `reader_loop`. Without the
/// guard the panic path would skip `mark_consumer_gone`, leaving a writer
/// parked in `gate.acquire` blocked forever and `Drop`'s `writer_thread.join()`
/// hung with it.
pub(crate) struct ConsumerGoneGuard<'a>(pub(crate) &'a InFlightGate);

impl Drop for ConsumerGoneGuard<'_> {
    fn drop(&mut self) {
        self.0.mark_consumer_gone();
    }
}

/// Aligner-parsed but pre-merge templates paired with their original unmapped
/// halves — the backend-agnostic edge every align backend produces, and the
/// input to [`merge::merge_zipper_batch`] (`merge_raw` + optional bisulfite
/// restore). The subprocess backend hands it to the shared
/// [`merge::MergeAlignedStep`]; the in-process pair/emit step merges it itself.
#[derive(Debug)]
pub(crate) struct ZipperBatch {
    pub(crate) serial: u64,
    pub(crate) mapped: Vec<Template>,
    pub(crate) unmapped: BamTemplateBatch,
}

impl HeapSize for ZipperBatch {
    fn heap_size(&self) -> usize {
        // Matches `BamTemplateBatch::total_bytes` semantics: summed
        // `Template::heap_size` plus the `Vec`'s reserved backing store, so
        // both halves are counted the same way. `Template` carries only an
        // inherent `heap_size` method (not a `HeapSize` trait impl) to keep
        // `template.rs` framework-agnostic, so the sum is hand-rolled here.
        let mapped_heap: usize = self.mapped.iter().map(Template::heap_size).sum::<usize>()
            + container_bytes::<Template>(self.mapped.capacity());
        mapped_heap + self.unmapped.heap_size()
    }
}

impl Ordered for ZipperBatch {
    fn ordinal(&self) -> u64 {
        self.serial
    }
}

/// What `add_align` needs to wire any align backend. Built once by the chain
/// builder and passed to [`AlignBackend::wire`].
///
/// Carries the output header and its handle, the per-step queue byte limit, the
/// pool's thread budget, and the shared zipper-merge config. The thread budget
/// is read only by the in-process backend, so it is `allow(dead_code)` in
/// feature-off builds.
pub(crate) struct AlignWiringCtx {
    /// Partial output header (dict `@SQ` + unmapped `@HD`/`@CO`/`@RG`/`@PG` +
    /// fgumi `@PG`) the backend validates the aligner's `@SQ` against and merges
    /// the aligner's runtime `@PG`/`@RG`/`@CO` into.
    pub(crate) partial_output_header: Arc<Header>,
    /// One-shot handle the backend resolves (set or poison) exactly once with
    /// the merged header; the downstream `WriteBgzfFile` polls it.
    pub(crate) header_handle: HeaderHandle,
    /// Byte limit for the backend's output queue (`ChainTuning::per_step_byte_limit`).
    pub(crate) per_step_byte_limit: u64,
    /// Pool worker count (`spec.threading.num_threads()`) — the single thread
    /// budget. The in-process backend loads the index with this many threads and
    /// sizes its seed/extend output queue for this many in-flight worker clones;
    /// the subprocess backend ignores it (its parallelism is the aligner
    /// subprocess's own `-t`), so it is unread in a feature-off build.
    #[cfg_attr(not(feature = "aligner-bwa-mem3"), allow(dead_code))]
    pub(crate) num_threads: usize,
    /// The zipper merge every aligned template goes through: the subprocess
    /// backend runs it in the shared [`merge::MergeAlignedStep`] it appends, the
    /// in-process pair/emit step inline.
    pub(crate) merge: Arc<merge::MergeConfig>,
}

/// The result of wiring a backend: its tail plus the scheduling facts the chain
/// builder folds into pool sizing and scheduler selection.
pub(crate) struct AlignWired {
    /// The backend's tail: merged `BamTemplateBatch`es in input order (dense
    /// serials from 0).
    pub(crate) tail: (StepIdx, BranchIdx),
    /// Minimum pool workers this backend needs for steady-state progress.
    pub(crate) min_workers: usize,
    /// Whether the chain builder should prefer drain-first dispatch.
    pub(crate) prefers_drain_first: bool,
    /// When set and drain-first is chosen automatically, the chain uses
    /// [`RefillDrainScheduler`](crate::pipeline::core::runtime::RefillDrainScheduler)
    /// on this hint.
    pub(crate) refill: Option<RefillHint>,
}

/// A backend's input-refill hint: while `signal` is raised and the queue
/// `feed` (a producer step and output branch) holds less than `cap_bytes`, the
/// steps up to and including the producer should be walked upstream-first.
#[derive(Clone)]
pub(crate) struct RefillHint {
    pub(crate) signal: Arc<std::sync::atomic::AtomicBool>,
    pub(crate) feed: (StepIdx, BranchIdx),
    pub(crate) cap_bytes: u64,
}

/// A pluggable align backend: appends its steps after a queryname-grouped
/// `BamTemplateBatch` tail — aligning each batch and zipper-merging the aligned
/// half with the unmapped half — and returns the merged `BamTemplateBatch` tail.
pub(crate) trait AlignBackend: Send + 'static {
    /// Identity for logs and for the synthesized/merged `@PG`.
    fn describe(&self) -> String;

    /// Append this backend's steps after `input` (a `BamTemplateBatch` tail,
    /// queryname-grouped) and return its merged tail. Must resolve or
    /// arrange to resolve `ctx.header_handle` (set or poison) exactly once.
    ///
    /// # Errors
    ///
    /// Returns an error if the backend fails to construct its step(s) (e.g. the
    /// subprocess aligner fails to spawn).
    fn wire(
        self: Box<Self>,
        pipeline: &PipelineBuilder,
        input: (StepIdx, BranchIdx),
        ctx: &AlignWiringCtx,
    ) -> anyhow::Result<AlignWired>;
}

/// Construct the boxed align backend a [`ResolvedAligner`] selected.
///
/// Called once by the chain builder (`chains/builder.rs::add_align`) after
/// [`AlignerOptions::resolve`](crate::aligner::AlignerOptions::resolve) has
/// validated the option combination. Both arms only package already-validated
/// fields, so this is infallible: the in-process backend loads its index later,
/// inside [`AlignBackend::wire`], which is where that failure surfaces.
///
/// Lives here rather than on `ResolvedAligner` so the framework-agnostic
/// `aligner` module does not depend on the pipeline's align stage.
pub(crate) fn backend_for(resolved: ResolvedAligner) -> Box<dyn AlignBackend> {
    match resolved.backend {
        ResolvedBackend::Subprocess { command, accept_mid_pair_split } => {
            Box::new(subprocess::SubprocessBackend {
                command,
                chunk_size: resolved.chunk_size,
                accept_mid_pair_split,
            })
        }
        #[cfg(feature = "aligner-bwa-mem3")]
        ResolvedBackend::InProcessBwaMem3 { reference, sub_batch_templates, dedup_reads } => {
            // The index is loaded (once) and the header synthesized/resolved
            // inside `AlignBackend::wire`, not here — this only packages the
            // resolved parameters. `chunk_size` drives the cohort cutter/gate and
            // the synthesized `@PG CL`.
            Box::new(inproc::InProcessBwaMem3Backend {
                reference,
                sub_batch_templates,
                dedup_reads,
                chunk_size: resolved.chunk_size,
            })
        }
        #[cfg(not(feature = "aligner-bwa-mem3"))]
        ResolvedBackend::InProcessBwaMem3 { .. } => {
            unreachable!(
                "AlignerOptions::resolve bails for AlignerPreset::BwaMem3InProc \
                 when the `aligner-bwa-mem3` feature is off, so this variant \
                 cannot be constructed in a feature-off build"
            )
        }
    }
}

/// Compare the aligner's emitted `@SQ` table to the partial output
/// header's. The partial header's `@SQ` came from the reference dict
/// at construction time; the aligner's `@SQ` comes from whatever
/// FASTA the aligner was indexed against. If they don't match, the
/// aligner's per-record `tid` integers index into a different
/// `@SQ` ordering and the merged BAM would silently have corrupt
/// reference IDs.
///
/// Comparison is name + length only (M5/UR/AS/SP are dict-specific
/// fields the aligner doesn't propagate, so we don't require them).
pub(crate) fn validate_sq_consistency(partial: &Header, aligner: &Header) -> io::Result<()> {
    let partial_refs = partial.reference_sequences();
    let aligner_refs = aligner.reference_sequences();

    if partial_refs.len() != aligner_refs.len() {
        return Err(io::Error::other(format!(
            "align-and-merge: aligner @SQ count ({}) does not match reference dict @SQ count ({}). \
             The aligner was indexed against a different FASTA than the supplied --ref. \
             Re-index or fix the --ref path.",
            aligner_refs.len(),
            partial_refs.len(),
        )));
    }

    for (idx, ((p_name, p_map), (a_name, a_map))) in
        partial_refs.iter().zip(aligner_refs.iter()).enumerate()
    {
        if p_name != a_name {
            return Err(io::Error::other(format!(
                "align-and-merge: aligner @SQ name mismatch at position {idx}: dict='{dict_name}' \
                 aligner='{aln_name}'. The aligner was indexed against a different \
                 FASTA than the supplied --ref.",
                dict_name = String::from_utf8_lossy(p_name),
                aln_name = String::from_utf8_lossy(a_name),
            )));
        }
        if p_map.length() != a_map.length() {
            return Err(io::Error::other(format!(
                "align-and-merge: aligner @SQ length mismatch for '{name}': dict={dict_len} \
                 aligner={aln_len}. The aligner was indexed against a different \
                 FASTA than the supplied --ref.",
                name = String::from_utf8_lossy(p_name),
                dict_len = p_map.length(),
                aln_len = a_map.length(),
            )));
        }
    }

    Ok(())
}

/// The pairing bits [`split_pair_into_singles`] clears on each half's unmapped
/// record: everything that marks it as one segment of a pair. `QC_FAIL`,
/// `REVERSE` and the rest are kept.
const SPLIT_HALF_CLEARED_FLAGS: u16 = fgumi_raw_bam::flags::PAIRED
    | fgumi_raw_bam::flags::PROPER_PAIR
    | fgumi_raw_bam::flags::MATE_UNMAPPED
    | fgumi_raw_bam::flags::MATE_REVERSE
    | fgumi_raw_bam::flags::FIRST_SEGMENT
    | fgumi_raw_bam::flags::LAST_SEGMENT;

/// Split a two-primary-read template into two single-record templates, in record
/// order: bwa's mid-pair split. With mixed SE/PE input a `-K`
/// chunk boundary can fall between a pair's two reads, and bwa then aligns them as
/// two unpaired reads. Both backends split the pair's unmapped template the same
/// way so each half zips with its own read's alignment. Errors if the template
/// carries records beyond its two primaries (secondaries make the split ambiguous —
/// out of contract for the valid unmapped input the parity claim covers).
///
/// Each half's record has its pairing bits ([`SPLIT_HALF_CLEARED_FLAGS`]) cleared,
/// so it is an unpaired read like the aligner's record for it. The zipper merge
/// picks the mapped segment to copy tags and the QC-fail flag onto from the
/// *unmapped* record's `PAIRED`/`FIRST_SEGMENT` bits; left set, the second half
/// (still `PAIRED | LAST_SEGMENT`) would look for a mapped R2, find none (the
/// aligner emitted it unpaired, which files as R1), and silently copy nothing.
/// The merged record's own flags come from the aligner, so clearing these bits
/// on the unmapped side changes nothing else in the output.
pub(crate) fn split_pair_into_singles(template: Template) -> io::Result<(Template, Template)> {
    // Reject anything but exactly two primaries before consuming the template,
    // so the error path reads `template.name()` by reference and the success
    // path never allocates a name copy.
    if template.records().len() != 2 {
        return Err(io::Error::other(format!(
            "align-and-merge: cannot split template '{name}' across a mid-pair -K chunk \
             boundary: it carries {n} records, not exactly two primaries (secondary/supplementary \
             records are out of contract for mixed single/paired input).",
            name = String::from_utf8_lossy(template.name()),
            n = template.records().len(),
        )));
    }
    let mut records = template.into_records();
    let second = records.pop().expect("len checked == 2");
    let first = records.pop().expect("len checked == 2");
    let to_template = |mut record: fgumi_raw_bam::RawRecord| {
        record.set_flags(record.flags() & !SPLIT_HALF_CLEARED_FLAGS);
        Template::from_records(vec![record])
            .map_err(|e| io::Error::other(format!("align-and-merge: split-half template: {e:#}")))
    };
    Ok((to_template(first)?, to_template(second)?))
}

/// Delay mimalloc returning freed pages to the OS (by
/// [`fgumi_sort::RETAINED_PURGE_DELAY_MS`]), unless the user set the purge
/// delay. The align stage frees and reallocates large per-batch buffers on
/// every pool thread, and mimalloc's default purge (decommit after 1 s) turns
/// that reuse into millions of page faults.
///
/// The setting is process-wide and is never restored, so it covers every stage
/// of the `runall` process, not just align. The delay is bounded rather than
/// never-purge on purpose: with `-1`, memory a fused sort frees at each spill
/// is never returned, and peak RSS grows by roughly the sort's memory budget
/// (+5.5 GB in-process on 30M pairs, enough to exceed a 30 GB limit). Measured
/// on 3M pairs at 32 threads (extract through simplex consensus, subprocess and
/// in-process aligners, with and without a spilling sort), the 60 s delay
/// matches never-purging on wall time, CPU and page faults (~30k vs ~750k at
/// mimalloc's default, which is ~3% slower); on 30M pairs it gives back most of
/// never-purging's extra memory (fgumi 13.4 GB vs 17.6 GB, 11.6 GB at the
/// default). Set `MIMALLOC_PURGE_DELAY` to choose otherwise.
pub(crate) fn retain_freed_memory_unless_user_set() {
    if !crate::aligner::user_set_mimalloc_purge() {
        fgumi_sort::retain_freed_memory();
    }
    log::debug!("align: mimalloc purge delay {} ms", fgumi_sort::mi_purge_delay_ms());
}

/// Merge aligner-emitted header lines into the partial output header.
///
/// The aligner contributes:
/// - `@PG` lines (with bwa version, command line, etc.) — appended.
///   On duplicate ID, the partial's existing PG wins (the aligner's
///   PG with that ID is dropped). **Limitation:** this does not
///   maintain the BAM spec's `@PG.PP` chain when the aligner adds
///   its own PG to a chain that already has one with the same ID
///   (rare in practice — bwa's PG ID is "bwa", which won't be in
///   `fgumi extract`'s unmapped BAM output). Maintaining a proper
///   PP chain requires constructing fresh IDs and threading the
///   PP pointer; deferred to a follow-up commit if real-world data
///   surfaces a collision.
/// - `@RG` lines (if `-R` was passed to the aligner) — appended;
///   partial's existing RG wins on duplicate ID.
/// - `@CO` comment lines — concatenated (partial first, then
///   aligner).
///
/// The aligner's `@SQ` lines are deliberately discarded:
/// [`validate_sq_consistency`] runs first to ensure the two `@SQ`
/// tables agree, after which the partial's (dict-derived) version
/// is authoritative.
pub(crate) fn merge_aligner_header(partial: &Header, aligner: &Header) -> Header {
    use bstr::BString;

    let mut builder = Header::builder();

    if let Some(hd) = partial.header() {
        builder = builder.set_header(hd.clone());
    }

    for (name, map) in partial.reference_sequences() {
        builder = builder.add_reference_sequence(name.clone(), map.clone());
    }

    let mut rg_seen: std::collections::HashSet<BString> = std::collections::HashSet::new();
    for (id, rg) in partial.read_groups() {
        builder = builder.add_read_group(id.clone(), rg.clone());
        rg_seen.insert(id.clone());
    }
    for (id, rg) in aligner.read_groups() {
        if !rg_seen.contains(id) {
            builder = builder.add_read_group(id.clone(), rg.clone());
        }
    }

    let mut pg_seen: std::collections::HashSet<BString> = std::collections::HashSet::new();
    for (id, pg) in partial.programs().as_ref() {
        builder = builder.add_program(id.clone(), pg.clone());
        pg_seen.insert(id.clone());
    }
    for (id, pg) in aligner.programs().as_ref() {
        if !pg_seen.contains(id) {
            builder = builder.add_program(id.clone(), pg.clone());
        }
    }

    for c in partial.comments() {
        builder = builder.add_comment(c.clone());
    }
    for c in aligner.comments() {
        builder = builder.add_comment(c.clone());
    }

    builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Only a template of exactly two records is split: any other count is out
    /// of contract and errors, naming the template and its record count.
    #[rstest::rstest]
    #[case::one_record(&[fgumi_raw_bam::flags::UNMAPPED], 1)]
    #[case::pair_with_supplementary(
        &[
            fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::FIRST_SEGMENT,
            fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::LAST_SEGMENT,
            fgumi_raw_bam::flags::PAIRED
                | fgumi_raw_bam::flags::FIRST_SEGMENT
                | fgumi_raw_bam::flags::SUPPLEMENTARY,
        ],
        3
    )]
    fn split_pair_into_singles_rejects_other_than_two_records(
        #[case] flags: &[u16],
        #[case] n: usize,
    ) {
        let records = flags
            .iter()
            .map(|&f| {
                let mut b = fgumi_raw_bam::SamBuilder::new();
                b.read_name(b"q1").flags(f).sequence(b"ACGT").qualities(b"IIII");
                b.build()
            })
            .collect();
        let template = Template::from_records(records).expect("template");
        let err = split_pair_into_singles(template).expect_err("must not split").to_string();
        assert!(
            err.contains("cannot split template 'q1'")
                && err.contains(&format!("carries {n} records")),
            "got: {err}"
        );
    }

    /// A `ZipperBatch` sits in a byte-bounded queue, so its `heap_size` must
    /// count both halves' `Vec` backing store, as `BamTemplateBatch` already
    /// does for the unmapped half. Empty `Vec`s with reserved capacity isolate
    /// the container term: no `Template` contributes any heap of its own.
    #[test]
    fn zipper_batch_heap_size_counts_both_halves_container_capacity() {
        let unmapped = BamTemplateBatch::new(0, Vec::with_capacity(8));
        let unmapped_bytes = unmapped.total_bytes();
        assert_eq!(unmapped_bytes, 8 * std::mem::size_of::<Template>());
        let batch = ZipperBatch { serial: 0, mapped: Vec::with_capacity(16), unmapped };
        assert_eq!(batch.heap_size(), 16 * std::mem::size_of::<Template>() + unmapped_bytes);
    }

    #[test]
    fn in_flight_budget_derivation() {
        // Default -K (150M bases) → 4 × 150M × 3 = 1.8 GiB.
        assert_eq!(in_flight_budget_for_chunk_size(150_000_000), 4 * 150_000_000 * 3);
        // Small / zero -K floors at IN_FLIGHT_MIN_BUDGET.
        assert_eq!(in_flight_budget_for_chunk_size(0), IN_FLIGHT_MIN_BUDGET);
        assert_eq!(in_flight_budget_for_chunk_size(1_000_000), IN_FLIGHT_MIN_BUDGET);
    }

    #[test]
    fn consumer_gone_latched_even_when_reader_scope_panics() {
        // A panic unwinding out of the reader scope must still latch
        // consumer-gone via `ConsumerGoneGuard`, so a writer parked in
        // `acquire` bails (returns false) instead of blocking forever. Without
        // the guard the panic would skip `mark_consumer_gone`, deadlocking the
        // writer and the `writer_thread.join()` in `Drop`.
        let gate = InFlightGate::new(1024);
        assert!(gate.acquire(10), "acquire succeeds before the consumer is gone");

        let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _consumer_gone = ConsumerGoneGuard(&gate);
            panic!("simulated reader_loop panic");
        }));
        assert!(panicked.is_err(), "the panic must propagate out of the guarded scope");

        // The guard's `Drop` latched consumer-gone during unwind, so a
        // subsequent writer reservation bails instead of blocking.
        assert!(!gate.acquire(10), "a writer must bail once the reader scope has panicked");
    }

    #[test]
    fn in_flight_gate_oversized_single_batch_passes_when_empty() {
        // A batch larger than the whole budget must still pass when the gate
        // is empty — it can't be split, so holding it is unavoidable, and
        // blocking would deadlock.
        let gate = InFlightGate::new(100);
        assert!(gate.acquire(1000), "oversized batch must pass when gate empty");
    }

    #[test]
    fn in_flight_gate_blocks_until_release() {
        use std::time::Duration;
        let gate = Arc::new(InFlightGate::new(100));
        assert!(gate.acquire(100), "first acquire fills the budget");

        let g2 = Arc::clone(&gate);
        let (tx, rx) = std::sync::mpsc::channel();
        let h = std::thread::spawn(move || {
            // in_flight is 100 (full, non-empty) → this blocks until release.
            assert!(g2.acquire(50));
            tx.send(()).unwrap();
        });
        // Still blocked.
        assert!(
            rx.recv_timeout(Duration::from_millis(150)).is_err(),
            "second acquire must block while the gate is full"
        );
        gate.release(100);
        rx.recv_timeout(Duration::from_secs(5)).expect("acquire must unblock after release");
        h.join().unwrap();
    }

    #[test]
    fn in_flight_gate_consumer_gone_unblocks_and_bails() {
        use std::time::Duration;
        let gate = Arc::new(InFlightGate::new(100));
        assert!(gate.acquire(100), "fill the budget");

        let g2 = Arc::clone(&gate);
        let (tx, rx) = std::sync::mpsc::channel();
        let h = std::thread::spawn(move || {
            let ok = g2.acquire(50); // blocks (full) until consumer-gone
            tx.send(ok).unwrap();
        });
        assert!(rx.recv_timeout(Duration::from_millis(150)).is_err(), "must be blocked");
        gate.mark_consumer_gone();
        let ok = rx.recv_timeout(Duration::from_secs(5)).expect("must wake on consumer-gone");
        assert!(!ok, "acquire must return false once the consumer is gone");
        h.join().unwrap();
    }

    fn make_sq_header(refs: &[(&str, usize)]) -> Header {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReferenceSequence;
        let mut b = Header::builder();
        for (name, length) in refs {
            let map: Map<ReferenceSequence> = Map::<ReferenceSequence>::new(
                std::num::NonZeroUsize::new(*length).expect("nonzero ref length"),
            );
            b = b.add_reference_sequence(bstr::BString::from(*name), map);
        }
        b.build()
    }

    #[test]
    fn validate_sq_consistency_accepts_matching_refs() {
        let partial = make_sq_header(&[("chr1", 1000), ("chr2", 2000)]);
        let aligner = make_sq_header(&[("chr1", 1000), ("chr2", 2000)]);
        validate_sq_consistency(&partial, &aligner).expect("matching refs are ok");
    }

    #[test]
    fn validate_sq_consistency_rejects_count_mismatch() {
        let partial = make_sq_header(&[("chr1", 1000), ("chr2", 2000)]);
        let aligner = make_sq_header(&[("chr1", 1000)]);
        let err =
            validate_sq_consistency(&partial, &aligner).expect_err("count mismatch must reject");
        let msg = err.to_string();
        assert!(msg.contains("aligner @SQ count (1) does not match reference dict @SQ count (2)"));
    }

    #[test]
    fn validate_sq_consistency_rejects_length_mismatch() {
        let partial = make_sq_header(&[("chr1", 1000)]);
        let aligner = make_sq_header(&[("chr1", 999)]);
        let err =
            validate_sq_consistency(&partial, &aligner).expect_err("length mismatch must reject");
        let msg = err.to_string();
        assert!(msg.contains("aligner @SQ length mismatch for 'chr1'"), "got: {msg}");
    }

    #[test]
    fn validate_sq_consistency_rejects_name_mismatch() {
        // Two refs so the position-index reporting can be exercised
        // (round-1 had a format-args bug that put the dict name in
        // the position slot; a single-ref test wouldn't have caught
        // it).
        let partial = make_sq_header(&[("chr1", 1000), ("chr2", 1000)]);
        let aligner = make_sq_header(&[("chr1", 1000), ("chrZ", 1000)]);
        let err =
            validate_sq_consistency(&partial, &aligner).expect_err("name mismatch must reject");
        let msg = err.to_string();
        assert!(msg.contains("aligner @SQ name mismatch at position 1"), "got: {msg}");
        assert!(msg.contains("dict='chr2'"), "dict name in message: {msg}");
        assert!(msg.contains("aligner='chrZ'"), "aligner name in message: {msg}");
    }

    /// `merge_aligner_header` must keep the *partial* header's `@PG` on a
    /// duplicate ID (the aligner's PG with that ID is dropped, not merged)
    /// while still appending aligner PGs whose IDs are unique to the
    /// aligner. This pins the dedup behavior documented on the function.
    #[test]
    fn merge_aligner_header_keeps_partial_pg_on_duplicate_id() {
        use bstr::BString;
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::Program;
        use noodles::sam::header::record::value::map::program::tag as pg_tag;

        // Helper: build a @PG map carrying a distinguishing PN field so we
        // can tell the partial's PG apart from the aligner's.
        let program_with_name = |name: &str| -> Map<Program> {
            Map::<Program>::builder().insert(pg_tag::NAME, name).build().expect("valid @PG")
        };

        // Partial header: a "dup" PG (PN=partial) plus a partial-only PG.
        let partial = Header::builder()
            .add_program(BString::from("dup"), program_with_name("partial"))
            .add_program(BString::from("onlypartial"), program_with_name("partial"))
            .build();

        // Aligner header: a "dup" PG (PN=aligner, must be dropped) plus a
        // unique "bwa" PG (must be appended).
        let aligner = Header::builder()
            .add_program(BString::from("dup"), program_with_name("aligner"))
            .add_program(BString::from("bwa"), program_with_name("aligner"))
            .build();

        let merged = merge_aligner_header(&partial, &aligner);
        let programs = merged.programs();
        let programs = programs.as_ref();

        // Exactly three PGs survive: dup (partial's), onlypartial, bwa.
        assert_eq!(programs.len(), 3, "expected dup + onlypartial + bwa");

        // Read back the PN field for a given @PG ID.
        let pn_of = |id: &str| -> Option<String> {
            programs
                .get(&BString::from(id))
                .and_then(|pg| pg.other_fields().get(&pg_tag::NAME).map(ToString::to_string))
        };

        // The duplicate-ID PG is the PARTIAL's (PN=partial), not the
        // aligner's — the aligner's same-ID PG was dropped.
        assert_eq!(
            pn_of("dup").as_deref(),
            Some("partial"),
            "duplicate @PG ID must retain the partial's PG, not the aligner's"
        );
        // The partial-only PG is preserved.
        assert_eq!(pn_of("onlypartial").as_deref(), Some("partial"));
        // The aligner-unique PG is appended.
        assert_eq!(
            pn_of("bwa").as_deref(),
            Some("aligner"),
            "aligner @PG with a non-duplicate ID must be appended"
        );
    }
}
