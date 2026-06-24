//! `SortAndSpill` — first step of the runall-sort three-step chain.

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

use crate::sort::protocol::{MemoryChunkErased, SortPhase1Event};
use crate::types::RecordBatch;
use fgumi_pipeline_core::{
    Unpushed,
    held::HeldSlot,
    outputs::Single,
    queues::QueueSpec,
    reorder::BranchOrdering,
    step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// Max input batches consumed per `try_run` invocation in the `Ingesting` state.
const MAX_INGEST_BATCHES_PER_LOCK: usize = 8;

/// Max events emitted per `try_run` invocation in the `Emitting` state.
const MAX_EVENTS_PER_LOCK: usize = 8;

/// Return type of `SortStream::finalize`.
type FinalizeResult =
    (Vec<Arc<fgumi_sort::SortMergeSlot>>, Vec<MemoryChunkErased>, Vec<TempDir>, u64);

/// Wraps the three `*SortStream` variants.
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
    Emitting { pending_events: VecDeque<SortPhase1Event>, _temp_dirs: Vec<TempDir> },
    Done,
}

/// `Serial` step that ingests `RecordBatch`es into a streaming sort engine.
pub struct SortAndSpill {
    state: SortAndSpillState,
    held: HeldSlot<Unpushed<SortPhase1Event>>,
    output_byte_limit: u64,
    affinity: Affinity,
}

impl SortAndSpill {
    /// Build a `SortAndSpill` from a configured `RawExternalSorter` and the output `Header`.
    ///
    /// `output_byte_limit` byte-bounds the output event queue. The emitted
    /// `SortPhase1Event::MemoryChunk` variant retains sorted record chunks, so
    /// this queue must budget on bytes (`HeapSize`), not event count, to keep
    /// retained memory a function of configuration.
    ///
    /// # Errors
    ///
    /// Propagates any error from the chosen `into_*_stream` entry.
    pub fn from_sorter(
        sorter: RawExternalSorter,
        header: &Header,
        output_byte_limit: u64,
    ) -> Result<Self> {
        let stream = build_stream(sorter, header)?;
        Ok(Self {
            state: SortAndSpillState::Ingesting(Box::new(stream)),
            held: HeldSlot::new(),
            output_byte_limit,
            affinity: Affinity::None,
        })
    }

    /// Override the affinity hint.
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

    /// # Panics
    ///
    /// Panics if the number of spill slots exceeds `u32::MAX`, or if the number
    /// of in-memory chunks exceeds `u32::MAX`.
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
        // Derive both announced counts the same way — `u32::try_from` on the
        // number of pushed events — so the two overflow guards are symmetric and
        // can't drift (a `slot_count` cast that panics differently from a
        // `memory_chunk_count` accumulator would mask a protocol bug).
        let slot_count = u32::try_from(pending_events.len()).expect("spill slot count fits in u32");
        let memory_events_start = pending_events.len();
        for chunk in memory_chunks {
            if chunk.is_empty() {
                continue;
            }
            pending_events.push_back(SortPhase1Event::MemoryChunk {
                chunk: Arc::new(chunk),
                records_ingested_so_far: total_records,
            });
        }
        let memory_chunk_count = u32::try_from(pending_events.len() - memory_events_start)
            .expect("memory chunk count fits in u32");
        if pending_events.is_empty() {
            return Ok(None);
        }
        pending_events.push_back(SortPhase1Event::AllAnnounced {
            slot_count,
            memory_chunk_count,
            total_records,
        });
        Ok(Some((pending_events, temp_dirs)))
    }

    /// # Panics
    ///
    /// Panics if `self.state` is not `Ingesting`.
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
        Ok(())
    }

    /// # Panics
    ///
    /// Panics if `self.state` is not `Emitting`.
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
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
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

        if matches!(&self.state, SortAndSpillState::Ingesting(_)) {
            if self.ingest_batches(ctx)? {
                return Ok(StepOutcome::Progress);
            }
            if !ctx.input.is_drained() {
                return Ok(StepOutcome::NoProgress);
            }
            self.transition_to_emitting("try_run")?;
        }

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
