//! In-process bwa-mem3 aligner backend (feature `aligner-bwa-mem3`).
//!
//! This subgraph replaces the subprocess `bwa-mem3` invocation with bwa-mem3
//! linked into the fgumi binary through [`bwa_mem3_rs`], run as work items on
//! the pipeline's shared work-stealing pool. The whole module
//! tree is compiled only under `aligner-bwa-mem3`; default builds stay pure
//! Rust with no C toolchain (the module is gated at its declaration in the
//! parent `align` module).
//!
//! [`InProcessBwaMem3Backend`] is the [`AlignBackend`] that ties the four
//! in-process steps together: it loads the `Arc<BwaIndex>` once, builds the
//! shared `Arc<MemOpts>`, synthesizes + merges + resolves the output header
//! **synchronously at wire time** (unlike the subprocess backend, which resolves
//! it at runtime off the aligner's stdout), and appends
//! `AlignPrepare → AlignSeedExtend → CohortPeStat → AlignPairEmit`, returning the
//! `AlignPairEmit` tail, which already carries merged `BamTemplateBatch`es:
//! pair/emit runs the zipper merge itself.
//!
//! The parity-critical *pure-logic foundation* lives in [`cohort`]: the `-K`
//! even-parity cohort cut, the SE/PE template classification, and the
//! per-sub-batch id-base arithmetic — each verified by a proptest against a
//! literal Rust port of the corresponding upstream bwa-mem3 rule.

pub(crate) mod cohort;
pub(crate) mod engine;
pub(crate) mod gate;
pub(crate) mod header;
pub(crate) mod pair_emit;
pub(crate) mod pestat;
pub(crate) mod prepare;
pub(crate) mod scratch;
pub(crate) mod seed_extend;

use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use engine::{BwaIndex, MemOpts};

use crate::pipeline::core::builder::PipelineBuilder;
use crate::pipeline::core::step::StepOutcome;
use crate::pipeline::core::topology::{BranchIdx, StepIdx};

use super::{
    AlignBackend, AlignWired, AlignWiringCtx, RefillHint, merge_aligner_header,
    validate_sq_consistency,
};
use engine::BwaMem3Engine;
use gate::cohort_bound_for_chunk_size;
use header::{load_index_header_sidecar, synthesize_aligner_header};
use pair_emit::AlignPairEmitStep;
use pestat::CohortPeStatStep;
use prepare::AlignPrepareStep;
use scratch::ScratchPool;
use seed_extend::{AlignSeedExtendStep, seed_extend_output_byte_limit};

/// The in-process bwa-mem3 align backend.
///
/// Constructed by [`backend_for`](super::backend_for)
/// for the `bwa-mem3-inproc` preset. [`AlignBackend::wire`] loads the index once,
/// synthesizes + resolves the output header at wire time, and appends the four
/// in-process steps.
pub(crate) struct InProcessBwaMem3Backend {
    /// Reference FASTA path — the bwa-mem3 index prefix
    /// ([`BwaIndex::load_with_threads`] appends the index suffixes).
    pub(crate) reference: PathBuf,
    /// Templates per sub-batch (`default_sub_batch_templates()` default, or the hidden
    /// `--aligner::sub-batch-templates` override), passed to
    /// [`AlignPrepareStep`].
    pub(crate) sub_batch_templates: usize,
    /// The aligner's `-K` chunk size in bases, driving [`AlignPrepareStep`]'s
    /// cohort cutter and gate, and the synthesized `@PG CL`.
    pub(crate) chunk_size: u64,
}

/// What one `flush` of a `Serial` step's staged output did, so `try_run` can
/// report liveness honestly: only a call that actually moved or newly held an
/// item is `Progress`; re-failing to push an already-held item changed nothing
/// and is `Contention`, the convention every held-retry step in the tree follows
/// (reporting `Progress` there would make a wedged downstream look alive to the
/// deadlock monitor, and spin the worker instead of letting it park).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum FlushOutcome {
    /// Everything staged (and any held item) reached the output.
    Clear,
    /// This call pushed or popped work, and backpressure then held an item.
    NewlyHeld,
    /// The item already held from an earlier call was refused again; nothing
    /// changed.
    StillHeld,
}

impl FlushOutcome {
    /// The outcome `try_run` must return when this flush left the output
    /// blocked, or `None` when the output is clear and the step should carry on.
    #[must_use]
    pub(crate) fn blocked_outcome(self) -> Option<StepOutcome> {
        match self {
            Self::Clear => None,
            Self::NewlyHeld => Some(StepOutcome::Progress),
            Self::StillHeld => Some(StepOutcome::Contention),
        }
    }
}

impl InProcessBwaMem3Backend {
    /// Minimum pool workers this backend needs for steady-state progress:
    /// the in-process backend runs entirely on the
    /// shared pool with no dedicated aligner/daemon threads, unlike the
    /// subprocess backend's floor of 4
    /// ([`SubprocessAlignStep`](super::subprocess::SubprocessAlignStep)). This
    /// is also what keeps `--threads 1` eligible for the fused/inline path
    /// (`should_fuse_single_thread`'s `n_threads == 1` precondition,
    /// `crate::pipeline::core::runtime::fused`) — the subprocess backend's
    /// floor of 4 forecloses that path unconditionally.
    pub(crate) const MIN_WORKERS: usize = 1;
    /// Whether this backend prefers the chain builder's drain-first scheduler:
    /// downstream work (merge/serialize/compress/sort
    /// ingest) is ~10% of the total per-pair cost and cheap to keep drained,
    /// so keeping it drained is a throughput-and-latency win here (unlike the
    /// subprocess backend, whose `AlignWired::prefers_drain_first` stays
    /// `false`, preserving its pre-existing scheduling behavior).
    pub(crate) const PREFERS_DRAIN_FIRST: bool = true;
}

impl AlignBackend for InProcessBwaMem3Backend {
    fn describe(&self) -> String {
        format!(
            "in-process bwa-mem3 (bwa-mem3-rs {version}, sub-batch {sub}, -K {chunk})",
            version = engine::version(),
            sub = self.sub_batch_templates,
            chunk = self.chunk_size,
        )
    }

    fn wire(
        self: Box<Self>,
        pipeline: &PipelineBuilder,
        input: (StepIdx, BranchIdx),
        ctx: &AlignWiringCtx,
    ) -> anyhow::Result<AlignWired> {
        // No poison-on-error guard is needed for the header handle: if `wire`
        // returns `Err`, `ChainBuilder::add_align` propagates it before stashing
        // the handle for the sink, so no consumer can be blocked on it.
        super::retain_freed_memory_unless_user_set();

        // Load the index once, synchronously, on the shared thread budget. Log
        // wall time and whether a staged shm segment was attached. A staged
        // segment means bwa-mem3 attached to it instead of reading from disk; the
        // probe is best-effort (failure => assume disk).
        let shm_attached = engine::shm::is_staged(&self.reference).unwrap_or(false);
        let load_start = Instant::now();
        let idx = BwaIndex::load_with_threads(&self.reference, ctx.num_threads)?;
        let idx = Arc::new(idx);
        log::info!(
            "in-process bwa-mem3: loaded index '{reference}' in {elapsed:.2?} \
             ({source}, {n_contigs} contigs, {threads} load threads)",
            reference = self.reference.display(),
            elapsed = load_start.elapsed(),
            source = if shm_attached { "attached shared-memory segment" } else { "read from disk" },
            n_contigs = idx.n_contigs(),
            threads = ctx.num_threads,
        );

        // Build the shared, read-only options once. `set_pe(true)` documents the
        // paired-end intent (the engine sets/clears `MEM_F_PE` per group).
        let mut opts = MemOpts::new()?;
        opts.set_pe(true);
        let opts = Arc::new(opts);

        // Synthesize the aligner header from the index contigs and resolve the
        // shared output header at WIRE time, following the subprocess order:
        // validate @SQ, merge, set. An @SQ mismatch is a wire-time `Err` (the
        // subprocess path catches it later, on the reader thread).
        // bwa-mem3's per-index header sidecar (`<prefix>.hdr`/`<baseprefix>.dict`),
        // whose @RG/@PG/@CO the CLI copies into its header, so the subprocess
        // preset's output carries them too.
        let sidecar = load_index_header_sidecar(&self.reference)?;
        let synth = synthesize_aligner_header(
            idx.contigs(),
            sidecar.as_ref(),
            self.chunk_size,
            &self.reference,
        );
        validate_sq_consistency(&ctx.partial_output_header, &synth)?;
        let merged = merge_aligner_header(&ctx.partial_output_header, &synth);
        ctx.header_handle.set(merged).map_err(|_| {
            anyhow::anyhow!(
                "in-process bwa-mem3 backend: output HeaderHandle was already resolved before \
                 wiring — a wiring bug (this backend is the sole producer)"
            )
        })?;

        let engine = Arc::new(BwaMem3Engine::new(idx, opts));

        // Append the four in-process steps: Prepare (Serial) → SeedExtend
        // (Parallel) → PeStat (Serial) → PairEmit (Parallel). The Prepare step
        // owns the cohort cutter + in-flight gate (both derived from -K); the two
        // parallel steps and the barrier share the engine by `Arc`. SeedExtend's
        // output is sized for `T` in-flight worker clones.
        let prepare = AlignPrepareStep::new(
            self.chunk_size,
            self.sub_batch_templates,
            ctx.per_step_byte_limit,
        );
        let refill_signal = prepare.refill_signal();
        let prepare_tail = pipeline.append_step(prepare, input);
        // While the cohort gate wants input, read and decode ahead up to (not
        // including) AlignPrepare, about one cohort's worth. Stopping short of
        // AlignPrepare keeps the next cohort out of seed/extend until the pool
        // drains the current cohort's pair/emit: running the two together costs
        // ~6% CPU at 8 threads.
        let refill = refill_signal.map(|signal| RefillHint {
            signal,
            feed: input,
            cap_bytes: cohort_bound_for_chunk_size(self.chunk_size),
        });

        let seed_limit = seed_extend_output_byte_limit(ctx.per_step_byte_limit, ctx.num_threads);
        // One scratch per pool thread, shared by the two Parallel aligner steps.
        let scratch = Arc::new(ScratchPool::new());
        let seed = AlignSeedExtendStep::new(Arc::clone(&engine), Arc::clone(&scratch), seed_limit);
        let seed_tail = pipeline.append_step(seed, prepare_tail);

        let pestat = CohortPeStatStep::new(Arc::clone(&engine), ctx.per_step_byte_limit);
        let pestat_tail = pipeline.append_step(pestat, seed_tail);

        // Pair+emit also runs the zipper merge on each sub-batch while its
        // records are hot, so no separate merge step follows this backend.
        let pair_emit = AlignPairEmitStep::new(
            engine,
            scratch,
            Arc::clone(&ctx.merge),
            ctx.per_step_byte_limit,
        );
        let pair_tail = pipeline.append_step(pair_emit, pestat_tail);

        Ok(AlignWired {
            tail: pair_tail,
            min_workers: Self::MIN_WORKERS,
            prefers_drain_first: Self::PREFERS_DRAIN_FIRST,
            refill,
        })
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::{FlushOutcome, InProcessBwaMem3Backend};
    use crate::pipeline::core::step::StepOutcome;

    /// Only a flush that moved or newly held work reports `Progress`; re-failing
    /// an already-held item is `Contention`, so a wedged downstream does not look
    /// alive to the deadlock monitor.
    #[rstest]
    #[case::clear(FlushOutcome::Clear, None)]
    #[case::newly_held(FlushOutcome::NewlyHeld, Some(StepOutcome::Progress))]
    #[case::still_held(FlushOutcome::StillHeld, Some(StepOutcome::Contention))]
    fn flush_outcome_maps_to_step_outcome(
        #[case] flush: FlushOutcome,
        #[case] expected: Option<StepOutcome>,
    ) {
        assert_eq!(flush.blocked_outcome(), expected);
    }

    /// Pin the in-process backend's scheduling hints directly, without going
    /// through `wire()` — which needs a real bwa-mem3 index
    /// (`BwaIndex::load_with_threads`) and so cannot run as a plain unit test.
    /// `ChainBuilder::add_align` folds these two constants into the chain-wide
    /// thread floor and drain-first flag via `fold_align_wired_scheduling`
    /// (`chains/builder.rs`), whose own test pins the fold's behavior against
    /// these same values (`1`/`true`) plus the subprocess backend's `4`/`false`.
    /// A real `wire()`-level smoke test, and the full fused-chain `--threads 1`
    /// run, are covered by the env-gated real-index tests (see
    /// `engine::tests::bwa_mem3_engine_smoke_pair_emit` for the established
    /// `FGUMI_BWA_MEM3_TEST_REF`-gated pattern this module tree follows).
    #[test]
    fn scheduling_hints_are_one_worker_and_drain_first() {
        assert_eq!(InProcessBwaMem3Backend::MIN_WORKERS, 1);
        // Not `assert!(PREFERS_DRAIN_FIRST)`: clippy's `assertions_on_constants`
        // flags an `assert!` whose argument is itself a `const`, since the
        // compiler can prove it statically. Route the constant through a
        // non-`const` binding first so the check still runs as real test
        // assertions (this constant IS the fact under test — see the doc
        // comment above).
        let prefers_drain_first: bool =
            std::hint::black_box(InProcessBwaMem3Backend::PREFERS_DRAIN_FIRST);
        assert!(prefers_drain_first);
    }
}
