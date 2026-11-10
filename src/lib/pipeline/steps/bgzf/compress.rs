//! `BgzfCompress` mid-step. `Parallel + ByItemOrdinal`. Compresses each
//! `DecompressedBlock` into one logical `BgzfBlock` whose `bytes` field
//! contains the concatenation of however many physical 64-KiB BGZF blocks
//! `libdeflater` produces internally.
//!
//! ## 1:1 input-to-output mapping (Phase 3 design)
//!
//! Each `try_run` consumes exactly one input and emits exactly one output
//! with `batch_serial = input.batch_serial`. This preserves the
//! consecutive-ordinal invariant required by the downstream
//! `ReorderStage` without any sub-serial encoding.
//!
//! Internally, `InlineBgzfCompressor` may produce multiple physical BGZF
//! blocks when the input exceeds `BGZF_MAX_BLOCK_SIZE = 64 KiB`. We
//! concatenate those physical blocks into a single `BgzfBlock`'s bytes —
//! a valid BGZF stream is just a concatenation of independent blocks, so
//! the downstream writer can emit them verbatim.
//!
//! The compression itself is [`OutputCompressor`]'s, shared with the
//! range-mode sort merge, which compresses the blocks its output edge
//! refuses: `BgzfCompress<OutputBlock>` (the standalone sort's terminal)
//! passes those [`OutputBlock::Compressed`] items through unchanged. Every
//! other chain uses the default `BgzfCompress<DecompressedBlock>`.

use std::io;
use std::marker::PhantomData;
use std::sync::Arc;

use fgumi_pipeline_core::CounterSpec;
use fgumi_pipeline_io::output_compress::{
    CompressInput, OutputBlock, OutputCompressor, OutputDepth,
};

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::core::{PhaseCap, admit_input};
use crate::pipeline::steps::types::{BgzfBlock, DecompressedBlock};

/// Counter slot: blocks the producer had already compressed, passed through.
const PASSED_THROUGH: usize = 0;

/// `Parallel + ByItemOrdinal` BGZF compressor over input `I`: a
/// [`DecompressedBlock`] (the default) or an [`OutputBlock`], whose
/// already-compressed blocks are passed through.
pub struct BgzfCompress<I: CompressInput = DecompressedBlock> {
    /// Per-worker compressor (each worker holds its own fresh copy).
    compressor: OutputCompressor,
    held: HeldSlot<Unpushed<BgzfBlock>>,
    output_byte_limit: u64,
    /// Admission cap on how many pool workers compress at once. `None`
    /// (uncapped) by default. `ChainBuilder::add_sink` sets the standalone
    /// sort's phase-2 cap (`--merge-threads`) on its terminal, and
    /// `build_merge` caps `fgumi merge`'s terminal at `--threads` so the
    /// merge's extra pool workers only read and decompress; every other
    /// chain's output compressor is untouched.
    cap: Option<Arc<PhaseCap>>,
    /// Counts this step's pops in the producer's queue-depth counter (the
    /// standalone sort's terminal: the range merge gates its pool claims on
    /// that depth). `None` for every other chain.
    depth: Option<Arc<OutputDepth>>,
    _input: PhantomData<fn(I)>,
}

impl BgzfCompress {
    /// A compressor of [`DecompressedBlock`]s at `compression_level`; with
    /// `index_bam` it emits each block's `BamIndexManifest` (built from the
    /// uncompressed framed records) for the inline BAI indexer.
    #[must_use]
    pub fn new(compression_level: u32, output_byte_limit: u64, index_bam: bool) -> Self {
        Self::for_input(compression_level, output_byte_limit, index_bam)
    }
}

impl<I: CompressInput> BgzfCompress<I> {
    /// [`BgzfCompress::new`] over input `I`.
    #[must_use]
    pub fn for_input(compression_level: u32, output_byte_limit: u64, index_bam: bool) -> Self {
        Self {
            compressor: OutputCompressor::new(compression_level, index_bam),
            held: HeldSlot::new(),
            output_byte_limit,
            cap: None,
            depth: None,
            _input: PhantomData,
        }
    }

    /// Share a phase admission cap with the other pool steps of that phase.
    #[must_use]
    pub fn with_phase_cap(mut self, cap: Option<Arc<PhaseCap>>) -> Self {
        self.cap = cap;
        self
    }

    /// Count every block this step pops in `depth` (shared with the
    /// producer, which counts its pushes).
    #[must_use]
    pub fn with_output_depth(mut self, depth: Option<Arc<OutputDepth>>) -> Self {
        self.depth = depth;
        self
    }

    /// One `try_run` with the phase cap moved out of `self` (see `try_run`).
    fn run_once(
        &mut self,
        ctx: &mut StepCtx<'_, Self>,
        cap: Option<&PhaseCap>,
    ) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        // Admission (see the pipeline-core `admission` module): a capped step
        // takes its phase permit BEFORE popping (an empty poll takes none and
        // is recorded as an empty pop), so a refused clone leaves the item for
        // a clone that may run and reports `Capped`. Uncapped, this is a no-op
        // and the pop below runs exactly as without a cap. The permit is held
        // to the end of this call.
        let _permit = match admit_input(ctx.input, cap) {
            Ok(permit) => permit,
            Err(outcome) => return Ok(outcome),
        };

        let Some(parent) = ctx.input.pop() else {
            // No input this call. If upstream is drained, every item has been
            // processed (held output was flushed by the Contention preamble
            // above) and this step will never push again — report Finished.
            // For a Parallel step only the last clone to finish closes the
            // shared output (gated by the StepDrainCounter in the driver).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };
        if let Some(depth) = &self.depth {
            depth.note_popped();
        }
        let out = match parent.into_output_block() {
            OutputBlock::Raw(raw) => self.compressor.compress(raw)?,
            // Compressed by the producer with the same settings: its bytes
            // and manifest are what this step would have built.
            OutputBlock::Compressed(block) => {
                ctx.counters.add(PASSED_THROUGH, 1);
                block
            }
        };
        match ctx.outputs.push(out) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }
}

impl<I: CompressInput> Clone for BgzfCompress<I> {
    fn clone(&self) -> Self {
        Self {
            compressor: self.compressor.clone(),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
            cap: self.cap.clone(),
            depth: self.depth.clone(),
            _input: PhantomData,
        }
    }
}

impl<I: CompressInput> Step for BgzfCompress<I> {
    type Input = I;
    type Outputs = OrderedBytesSingle<BgzfBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "BgzfCompress",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn counters(&self) -> &'static [CounterSpec] {
        const SPECS: &[CounterSpec] = &[CounterSpec::new("passed_through", "blocks")];
        SPECS
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        // Moved out for the call: the permit borrows the cap while the block assembly
        // needs `&mut self`. A move, so no reference-count traffic.
        let cap = self.cap.take();
        let outcome = self.run_once(ctx, cap.as_deref());
        self.cap = cap;
        outcome
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }

    fn step_cap(&self) -> Option<fgumi_pipeline_core::StepCapRef<'_>> {
        self.cap.as_deref().map(fgumi_pipeline_core::StepCapRef::Fixed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The output compressor is phase-neutral: uncapped by default (every
    /// chain's terminal but the standalone sort's and the merge's shows no
    /// `cap=` token), capped only when a builder hands it a cap (the standalone
    /// sort's phase-2 cap, or the merge's `--threads` compression cap) — and
    /// worker copies share that same counter, not an equal-looking new one.
    #[test]
    fn phase_cap_is_none_by_default_and_shared_by_worker_copies() {
        use crate::pipeline::core::{PhaseCap, Step};
        let plain = BgzfCompress::new(1, 1 << 20, false);
        assert!(plain.step_cap().is_none(), "uncapped by default");

        let cap = PhaseCap::new("sort-phase2", 2);
        let capped = plain.with_phase_cap(Some(std::sync::Arc::clone(&cap)));
        assert_eq!(
            capped.step_cap().map(fgumi_pipeline_core::StepCapRef::label).as_deref(),
            Some("sort-phase2(2)")
        );
        let copy = capped.new_worker_copy();
        let _p1 = cap.try_acquire().expect("1 of 2");
        let _p2 = cap.try_acquire().expect("2 of 2");
        assert!(
            copy.step_cap()
                .and_then(fgumi_pipeline_core::StepCapRef::active)
                .expect("capped")
                .try_acquire()
                .is_none(),
            "the worker copy shares the exhausted counter"
        );
    }

    /// Output compression under `--merge-threads`: the phase-2 admission
    /// contract, for both inputs.
    #[test]
    fn capped_step_follows_the_admission_contract() {
        crate::pipeline::core::testing::assert_admission_contract(
            "sort-phase2",
            |cap| BgzfCompress::new(1, 1 << 20, false).with_phase_cap(cap),
            DecompressedBlock { batch_serial: 0, bytes: b"ACGT".repeat(16) },
        );
        crate::pipeline::core::testing::assert_admission_contract(
            "sort-phase2",
            |cap| BgzfCompress::<OutputBlock>::for_input(1, 1 << 20, false).with_phase_cap(cap),
            OutputBlock::Raw(DecompressedBlock { batch_serial: 0, bytes: b"ACGT".repeat(16) }),
        );
    }

    #[test]
    fn profile_advertises_parallel_byordinal() {
        let s = BgzfCompress::new(1, 1024, false);
        let p = s.profile();
        assert_eq!(p.name, "BgzfCompress");
        assert_eq!(p.kind, StepKind::Parallel);
        assert_eq!(p.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
        assert_eq!(s.counters()[PASSED_THROUGH].name, "passed_through");
    }

    // Compression itself (byte identity with a fresh compressor, the
    // single/multi-block assembly, manifests) is pinned in
    // `fgumi_pipeline_io::output_compress`; the driven `try_run` contract
    // (manifest or none, pass-through of pre-compressed blocks) in
    // `steps/chain_tests.rs`.
}
