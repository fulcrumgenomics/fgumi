//! `SortAndSpill` — first step of the runall-sort three-step chain.
//!
//! ```text
//! [RecordBatch] → SortAndSpill ──SortPhase1Event──> SortSpillDecompress ──SortPhase2Event──> SortMerge → [RecordBatch]
//!                   (Serial)        (Parallel)                                                    (Serial)
//! ```
//!
//! Drives the streaming-sort engine's ingestion + spill phase
//! (`*SortStream::push_records` and internal spill management). On input
//! drain, finalizes via `*SortStream::into_slot_setup`, converts the
//! resulting `SlotSetup` into `SortPhase1Event`s, and emits them
//! downstream:
//!
//! * One `SpillReady` event per spill file (carries the shared
//!   `Arc<SortMergeSlot>` constructed by the slot-setup conversion).
//! * One `MemoryChunk` event per residual in-memory chunk (after the
//!   par-sort prologue).
//!
//! See `docs/design/sort-step-split.md` for the locked design.
//!
//! ## State machine
//!
//! 1. **`Ingesting`** — the streaming-sort handle is collecting input
//!    `RecordBatch`es and managing internal spills via `push_records`.
//!    Transitions to `Emitting` when `try_run` observes
//!    `ctx.input.is_drained()`.
//! 2. **`Emitting`** — the conversion to `SlotSetup` has completed and
//!    we hold a `VecDeque<SortPhase1Event>` of events to push downstream
//!    plus the spill-directory `Vec<TempDir>` (kept alive for the
//!    duration of the merge — slot readers reference files inside).
//!    Transitions to `Done` when the queue is empty.
//! 3. **`Done`** — terminal; `try_run` returns `Finished` (every event has
//!    been emitted and the step will never push again).

use std::collections::VecDeque;
use std::io;
use std::sync::Arc;

use anyhow::Result;
use fgumi_sort::{
    CoordinateSortStream, QuerynameSortStream, QuerynameSortStreamGeneric, RawExternalSorter,
    SortOrder, TemplateCoordinateSortStream,
};
use noodles::sam::Header;
use tempfile::TempDir;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::Single;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::sort::protocol::{MemoryChunkErased, SortPhase1Event};
use crate::pipeline::steps::types::RecordBatch;

/// Max input batches consumed per `try_run` invocation in the `Ingesting`
/// state. Bounds how long this worker holds the Sort lock before yielding
/// to the round-robin so other steps make progress.
const MAX_INGEST_BATCHES_PER_LOCK: usize = 8;

/// Max events emitted per `try_run` invocation in the `Emitting` state.
/// Keeps the per-call work bounded; events are typically a handful, so
/// this is rarely the bottleneck.
const MAX_EVENTS_PER_LOCK: usize = 8;

/// Return type of `SortStream::finalize` and the per-queryname-variant
/// helpers. Factored out because the inner tuple is too wide for clippy's
/// `type_complexity` threshold.
type FinalizeResult =
    (Vec<Arc<fgumi_sort::SortMergeSlot>>, Vec<MemoryChunkErased>, Vec<TempDir>, u64);

/// Wraps the three `*SortStream` variants so `SortAndSpill` can hold one
/// concrete state regardless of sort order. Boxed (see `SortAndSpillState`)
/// because the variants are ~hundreds of bytes.
enum SortStream {
    Coordinate(CoordinateSortStream),
    Queryname(QuerynameSortStream),
    TemplateCoordinate(TemplateCoordinateSortStream),
}

impl SortStream {
    fn push_records<'a, I: IntoIterator<Item = &'a [u8]>>(&mut self, iter: I) -> Result<()> {
        match self {
            Self::Coordinate(s) => s.push_records(iter),
            Self::Queryname(s) => s.push_records(iter),
            Self::TemplateCoordinate(s) => s.push_records(iter),
        }
    }

    /// Finalize the stream: drain pending spill, par-sort residual buffer
    /// (matching legacy tie-break), open `Arc<SortMergeSlot>`s for each
    /// spill file. Returns the slots, K-erased memory chunks (as
    /// `Vec<MemoryChunkErased>`), and the temp-dir RAII handles that
    /// must outlive the slots' readers.
    fn finalize(self) -> Result<FinalizeResult> {
        match self {
            Self::Coordinate(s) => {
                let setup = s.into_slot_setup()?;
                let erased: Vec<MemoryChunkErased> =
                    setup.memory_chunks.into_iter().map(MemoryChunkErased::Coordinate).collect();
                Ok((setup.slots, erased, setup.temp_dirs, setup.total_records))
            }
            Self::Queryname(QuerynameSortStream::Lex(s)) => {
                let (slots, erased, dirs, total) = finalize_queryname_lex(s)?;
                Ok((slots, erased, dirs, total))
            }
            Self::Queryname(QuerynameSortStream::Natural(s)) => {
                let (slots, erased, dirs, total) = finalize_queryname_natural(s)?;
                Ok((slots, erased, dirs, total))
            }
            Self::TemplateCoordinate(s) => {
                let setup = s.into_slot_setup()?;
                let erased: Vec<MemoryChunkErased> = setup
                    .memory_chunks
                    .into_iter()
                    .map(MemoryChunkErased::TemplateCoordinate)
                    .collect();
                Ok((setup.slots, erased, setup.temp_dirs, setup.total_records))
            }
        }
    }
}

fn finalize_queryname_lex(
    s: QuerynameSortStreamGeneric<fgumi_sort::RawQuerynameLexKey>,
) -> Result<FinalizeResult> {
    let setup = s.into_slot_setup()?;
    let erased: Vec<MemoryChunkErased> =
        setup.memory_chunks.into_iter().map(MemoryChunkErased::QuerynameLex).collect();
    Ok((setup.slots, erased, setup.temp_dirs, setup.total_records))
}

fn finalize_queryname_natural(
    s: QuerynameSortStreamGeneric<fgumi_sort::RawQuerynameKey>,
) -> Result<FinalizeResult> {
    let setup = s.into_slot_setup()?;
    let erased: Vec<MemoryChunkErased> =
        setup.memory_chunks.into_iter().map(MemoryChunkErased::QuerynameNatural).collect();
    Ok((setup.slots, erased, setup.temp_dirs, setup.total_records))
}

fn build_stream(sorter: RawExternalSorter, header: &Header) -> Result<SortStream> {
    match sorter.sort_order() {
        SortOrder::Coordinate => Ok(SortStream::Coordinate(sorter.into_coordinate_stream(header)?)),
        SortOrder::Queryname(_) => Ok(SortStream::Queryname(sorter.into_queryname_stream(header)?)),
        SortOrder::TemplateCoordinate => {
            Ok(SortStream::TemplateCoordinate(sorter.into_template_coordinate_stream(header)?))
        }
    }
}

enum SortAndSpillState {
    Ingesting(Box<SortStream>),
    Emitting {
        pending_events: VecDeque<SortPhase1Event>,
        /// RAII for the spill directories. Slots' `BufReader<File>`
        /// references point inside these directories; the temp dirs
        /// must live until the slots are dropped (after `SortMerge`
        /// finishes consuming them).
        ///
        /// Held inside `SortAndSpill` rather than embedded in each
        /// event because the `SortAndSpill` step instance is alive for
        /// the entire pipeline run.
        _temp_dirs: Vec<TempDir>,
    },
    Done,
}

/// `Serial` step that ingests `RecordBatch`es into a streaming sort
/// engine, manages internal spills, and on input drain emits
/// `SortPhase1Event`s describing the spill files + residual in-memory
/// chunks to the downstream `SortSpillDecompress` step.
pub struct SortAndSpill {
    state: SortAndSpillState,
    held: HeldSlot<Unpushed<SortPhase1Event>>,
    output_capacity: usize,
    /// Running `records_ingested` counter — snapshotted into each
    /// emitted event. Stamped onto every spill ready by the streaming
    /// sort engine's `total_records` after `into_slot_setup`; for
    /// memory chunks we use the same value (records ingested at
    /// finalization time equals the final total).
    ///
    /// V1 (this commit): all events emitted on input drain, so the
    /// snapshot is always the final value. V2 (Phase 1/Phase 2 overlap
    /// follow-up) will move `SpillReady` emission into `try_run` per
    /// spill close and snapshot the counter at each spill.
    affinity: Affinity,
}

impl SortAndSpill {
    /// Build a `SortAndSpill` from a configured `RawExternalSorter` and
    /// the output `Header`.
    ///
    /// # Errors
    ///
    /// Propagates any error from the chosen `into_*_stream` entry.
    pub fn from_sorter(
        sorter: RawExternalSorter,
        header: &Header,
        output_capacity: usize,
    ) -> Result<Self> {
        let stream = build_stream(sorter, header)?;
        Ok(Self {
            state: SortAndSpillState::Ingesting(Box::new(stream)),
            held: HeldSlot::new(),
            output_capacity,
            affinity: Affinity::None,
        })
    }

    /// Override the affinity hint. Not load-bearing for correctness.
    #[must_use]
    pub fn with_affinity(mut self, affinity: Affinity) -> Self {
        self.affinity = affinity;
        self
    }

    fn flush_held(&mut self, ctx: &mut StepCtx<'_, Self>) -> bool {
        let Some(unpushed) = self.held.take() else {
            return true;
        };
        match ctx.outputs.retry(unpushed) {
            Ok(()) => true,
            Err(again) => {
                self.held.put(again);
                false
            }
        }
    }

    /// Consume the `Ingesting` state's `SortStream`, finalize it, and
    /// build the `pending_events` queue (plus its companion
    /// `_temp_dirs` RAII bundle) that the `Emitting` state will push
    /// to `out_chan`.
    ///
    /// Returns `Ok(None)` for empty inputs — no slots, no memory
    /// chunks. Callers leave `self.state == Done` (already set by the
    /// caller's `mem::replace`) and return early.
    ///
    /// Called from `try_run` once `ctx.input.is_drained()` is observed. The
    /// block was previously inlined at two call sites (the drain-detected path
    /// and a now-deleted `on_input_drained` race fallback); factoring it out
    /// kept the two in sync.
    ///
    /// `caller_label` is interpolated into the error message so
    /// future post-mortems can attribute a finalize failure to the
    /// correct call site without changing the error surface.
    ///
    /// # Panics
    ///
    /// Panics if the number of spill slots exceeds `u32::MAX`. Each slot
    /// corresponds to one on-disk spill file, so reaching `u32::MAX`
    /// (~4 billion) slots is physically unreachable in any real run — the
    /// process would exhaust file descriptors and disk long before. The
    /// `u32` cast is required because `slot_count` is carried in the
    /// `AllAnnounced` sentinel as a `u32`.
    fn finalize_into_pending(
        stream: SortStream,
        caller_label: &str,
    ) -> io::Result<Option<(VecDeque<SortPhase1Event>, Vec<TempDir>)>> {
        let (slots, memory_chunks, temp_dirs, total_records) = stream.finalize().map_err(|e| {
            io::Error::other(format!("SortAndSpill: {caller_label} stream finalize failed: {e:#}"))
        })?;

        let mut pending_events: VecDeque<SortPhase1Event> = VecDeque::new();
        for slot in slots {
            pending_events.push_back(SortPhase1Event::SpillReady {
                slot,
                path: std::path::PathBuf::new(),
                records_ingested_so_far: total_records,
            });
        }
        let slot_count = u32::try_from(pending_events.len()).expect("spill count fits in u32");
        let mut memory_chunk_count: u32 = 0;
        for chunk in memory_chunks {
            if chunk.is_empty() {
                continue;
            }
            memory_chunk_count += 1;
            pending_events.push_back(SortPhase1Event::MemoryChunk {
                chunk: Arc::new(chunk),
                records_ingested_so_far: total_records,
            });
        }
        if pending_events.is_empty() {
            // Empty input — no slots, no memory chunks. Don't emit
            // even AllAnnounced; SortMerge falls through to its
            // `is_drained` fallback with zero sources.
            return Ok(None);
        }
        // AllAnnounced sentinel — always emitted LAST so SortMerge
        // can transition WaitingForSetup → Merging as soon as the
        // count predicate matches. See
        // `docs/design/sort-step-split-parity-fix.md` Change 1.
        pending_events.push_back(SortPhase1Event::AllAnnounced {
            slot_count,
            memory_chunk_count,
            total_records,
        });
        Ok(Some((pending_events, temp_dirs)))
    }

    /// `try_run` Phase 1: pop up to `MAX_INGEST_BATCHES_PER_LOCK`
    /// batches from upstream and push records into the `SortStream`.
    /// Returns `true` if any records were ingested (caller returns
    /// `Progress`), `false` if upstream was empty this call.
    ///
    /// # Panics
    ///
    /// Panics if `self.state` is not `Ingesting`. Callers must check
    /// the state before invoking.
    fn ingest_batches(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<bool> {
        let SortAndSpillState::Ingesting(stream) = &mut self.state else {
            unreachable!("ingest_batches called outside Ingesting state");
        };
        let mut did_work = false;
        for _ in 0..MAX_INGEST_BATCHES_PER_LOCK {
            let Some(batch) = ctx.input.pop() else { break };
            did_work = true;
            stream.push_records(batch.iter_record_bytes()).map_err(|e| {
                io::Error::other(format!("SortAndSpill: push_records failed: {e:#}"))
            })?;
        }
        Ok(did_work)
    }

    /// If currently in `Ingesting`, finalize the stream and transition
    /// to `Emitting`. On empty input, leave state at `Done` (set by
    /// the internal `mem::replace`). No-op when already in `Emitting`
    /// or `Done`. Called from `try_run` once the input is drained.
    fn transition_to_emitting(&mut self, caller_label: &str) -> io::Result<()> {
        if !matches!(&self.state, SortAndSpillState::Ingesting(_)) {
            return Ok(());
        }
        let SortAndSpillState::Ingesting(stream) =
            std::mem::replace(&mut self.state, SortAndSpillState::Done)
        else {
            unreachable!("just matched Ingesting")
        };
        if let Some((pending_events, temp_dirs)) =
            Self::finalize_into_pending(*stream, caller_label)?
        {
            self.state = SortAndSpillState::Emitting { pending_events, _temp_dirs: temp_dirs };
        }
        // else: empty input — state stays Done.
        Ok(())
    }

    /// Cooperative event emit. Push up to `MAX_EVENTS_PER_LOCK`
    /// events from the `Emitting` state's `pending_events` queue,
    /// stashing the first rejected push in `held`, returning the
    /// appropriate `StepOutcome` so the framework can interleave
    /// other steps' dispatches. Transitions state to `Done` when
    /// the queue is empty; `try_run` then reports `Finished`.
    ///
    /// # Panics
    ///
    /// Panics if `self.state` is not `Emitting`. Callers must check
    /// the state before invoking.
    fn emit_pending_cooperative(&mut self, ctx: &mut StepCtx<'_, Self>) -> StepOutcome {
        let mut emitted = 0usize;
        let drained;
        {
            let SortAndSpillState::Emitting { pending_events, .. } = &mut self.state else {
                unreachable!("emit_pending_cooperative called outside Emitting state");
            };
            while emitted < MAX_EVENTS_PER_LOCK {
                let Some(event) = pending_events.pop_front() else { break };
                if let Err(unpushed) = ctx.outputs.push(event) {
                    self.held.put(unpushed);
                    return StepOutcome::Progress;
                }
                emitted += 1;
            }
            drained = pending_events.is_empty();
        }
        if drained {
            self.state = SortAndSpillState::Done;
        }
        if emitted > 0 || !drained { StepOutcome::Progress } else { StepOutcome::NoProgress }
    }
}

impl Step for SortAndSpill {
    type Input = RecordBatch;
    type Outputs = Single<SortPhase1Event>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "SortAndSpill",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::CountBounded { capacity: self.output_capacity }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn affinity(&self) -> Affinity {
        self.affinity
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if !self.flush_held(ctx) {
            return Ok(StepOutcome::Contention);
        }

        // Phase 1: ingest from upstream while we're in Ingesting. If
        // upstream is drained AND we ingested nothing this call,
        // transition Ingesting → Emitting (or directly to Done on
        // empty input) and fall through.
        if matches!(&self.state, SortAndSpillState::Ingesting(_)) {
            if self.ingest_batches(ctx)? {
                return Ok(StepOutcome::Progress);
            }
            if !ctx.input.is_drained() {
                return Ok(StepOutcome::NoProgress);
            }
            self.transition_to_emitting("try_run")?;
        }

        // Phase 2: cooperative emit if we're in Emitting. Done reports
        // `Finished` — the step has emitted every event and will never push
        // again (it only transitions to Emitting on input drain, so reaching
        // Done means input is drained and all events are out). Ingesting is
        // unreachable — Phase 1 either returned early or transitioned us out.
        match &self.state {
            SortAndSpillState::Emitting { .. } => Ok(self.emit_pending_cooperative(ctx)),
            SortAndSpillState::Done => Ok(StepOutcome::Finished),
            SortAndSpillState::Ingesting(_) => {
                unreachable!("Phase 1 must have left state non-Ingesting")
            }
        }
    }
}

#[cfg(test)]
mod tests;
