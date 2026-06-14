//! Typed-event protocol between the three sort steps in the runall-sort
//! fused chain.
//!
//! ```text
//! [RecordBatch] → SortAndSpill ──SortPhase1Event──> SortSpillDecompress ──SortPhase2Event──> SortMerge → [RecordBatch]
//!                   (Serial)        (Parallel)                                                    (Serial)
//! ```
//!
//! v4 design: Decompressed BGZF block bytes do NOT flow on the typed
//! event chain. They live exclusively in `slot.decompressed` (per-slot
//! bounded queue on `SortMergeSlot`); `SortSpillDecompress` pushes,
//! `SortMerge` pops. The events on the chain carry only:
//!
//! * `SpillReady` — registers a slot (with `Arc<SortMergeSlot>`) for
//!   `SortSpillDecompress` to operate on AND for `SortMerge` to install
//!   in its slot table.
//! * `MemoryChunk` — residual in-memory chunk emitted by `SortAndSpill`'s
//!   drained-completion path; nothing for `SortSpillDecompress` to do
//!   (forwarded verbatim).
//! * `AllAnnounced` — sentinel emitted LAST by `SortAndSpill`. Carries
//!   the final counts so `SortMerge` can transition to `Merging` as
//!   soon as it has received the announced slot/chunk count, without
//!   waiting for `ctx.input.is_drained()` (which would require upstream
//!   workers to Skip themselves first — see
//!   `docs/design/sort-step-split-parity-fix.md`).
//!
//! See `docs/design/sort-step-split-parity-fix.md` for the full
//! locked design and rationale.

use std::path::PathBuf;
use std::sync::Arc;

use fgumi_raw_bam::RawRecord;
use fgumi_sort::{
    RawCoordinateKey, RawQuerynameKey, RawQuerynameLexKey, SortMergeSlot, TemplateKey,
};

use crate::pipeline::core::item::HeapSize;

/// Approximate fixed overhead of a `Vec<(K, RawRecord)>` entry (key + record
/// header). Used by `MemoryChunkErased::heap_size` to estimate buffer-bounded
/// queue pressure when a memory chunk flows through the typed chain.
///
/// Conservatively chosen to match the per-record memory budget used by
/// `RawExternalSorter` for spill triggering (`bytes_per_record = 354` in the
/// sort engine). Slightly over-counting is harmless — under-counting risks
/// queue overshoot.
const PER_MEMORY_RECORD_OVERHEAD: usize = 354;

/// In-memory sorted residual chunk produced by `SortAndSpill`, type-erased
/// over the sort-key variant `K` so it can flow through the typed-event
/// chain without `K` leaking into step signatures.
///
/// `SortMerge` is constructed with `sort_order` known statically; it
/// `match`es on the variant to recover the typed `Vec<(K, RawRecord)>` for
/// `MergeDriver::from_slots`'s `memory_chunks` argument.
pub enum MemoryChunkErased {
    /// Coordinate-sort residual. `K = RawCoordinateKey`.
    Coordinate(Vec<(RawCoordinateKey, RawRecord)>),
    /// Queryname-sort residual, lexicographic comparator. `K = RawQuerynameLexKey`.
    QuerynameLex(Vec<(RawQuerynameLexKey, RawRecord)>),
    /// Queryname-sort residual, natural comparator. `K = RawQuerynameKey`.
    QuerynameNatural(Vec<(RawQuerynameKey, RawRecord)>),
    /// Template-coordinate-sort residual. `K = TemplateKey`.
    TemplateCoordinate(Vec<(TemplateKey, RawRecord)>),
}

impl MemoryChunkErased {
    /// Number of records in this chunk.
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            Self::Coordinate(v) => v.len(),
            Self::QuerynameLex(v) => v.len(),
            Self::QuerynameNatural(v) => v.len(),
            Self::TemplateCoordinate(v) => v.len(),
        }
    }

    /// `true` iff the chunk has zero records.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Approximate heap footprint in bytes. Sum of each record's payload
    /// length plus a fixed per-entry overhead approximating `K + RawRecord`
    /// inline-struct cost. Used by `HeapSize` for byte-bounded queue
    /// accounting.
    #[must_use]
    pub fn approx_heap_bytes(&self) -> usize {
        let (count, payload): (usize, usize) = match self {
            Self::Coordinate(v) => {
                (v.len(), v.iter().map(|(_, r)| r.as_ref().len()).sum::<usize>())
            }
            Self::QuerynameLex(v) => {
                (v.len(), v.iter().map(|(_, r)| r.as_ref().len()).sum::<usize>())
            }
            Self::QuerynameNatural(v) => {
                (v.len(), v.iter().map(|(_, r)| r.as_ref().len()).sum::<usize>())
            }
            Self::TemplateCoordinate(v) => {
                (v.len(), v.iter().map(|(_, r)| r.as_ref().len()).sum::<usize>())
            }
        };
        count * PER_MEMORY_RECORD_OVERHEAD + payload
    }
}

/// Events from `SortAndSpill` → `SortSpillDecompress`.
pub enum SortPhase1Event {
    /// A spill chunk file has been closed and is ready for Phase 2
    /// decompression. Carries the per-file shared state, the on-disk
    /// path (informational), and a snapshot of
    /// `records_ingested_so_far` at spill-close.
    SpillReady { slot: Arc<SortMergeSlot>, path: PathBuf, records_ingested_so_far: u64 },
    /// A par-sorted residual in-memory chunk produced by `SortAndSpill`'s
    /// drained-completion path. Forwarded by `SortSpillDecompress` to
    /// `SortMerge` (no decompression needed).
    MemoryChunk { chunk: Arc<MemoryChunkErased>, records_ingested_so_far: u64 },
    /// Sentinel emitted as the LAST event by `SortAndSpill`'s drained-completion
    /// path (after all `SpillReady` and `MemoryChunk` events). Carries the final
    /// counts so `SortMerge`
    /// can transition to `Merging` early — without waiting for
    /// `ctx.input.is_drained()`. See
    /// `docs/design/sort-step-split-parity-fix.md` Change 1.
    AllAnnounced { slot_count: u32, memory_chunk_count: u32, total_records: u64 },
}

impl HeapSize for SortPhase1Event {
    fn heap_size(&self) -> usize {
        match self {
            Self::SpillReady { path, .. } => path.as_os_str().len(),
            Self::MemoryChunk { chunk, .. } => chunk.approx_heap_bytes(),
            // Three small numeric fields plus the discriminant.
            Self::AllAnnounced { .. } => 0,
        }
    }
}

/// Events from `SortSpillDecompress` → `SortMerge`.
///
/// Same shape as `SortPhase1Event`; `SortSpillDecompress` forwards all
/// three variants verbatim. (Decompressed Block bytes do NOT flow on
/// this chain — they live in `slot.decompressed`.)
pub enum SortPhase2Event {
    /// Forwarded `SortPhase1Event::SpillReady` — registers the slot in
    /// `SortMerge`'s slot table.
    SpillReady { slot: Arc<SortMergeSlot>, path: PathBuf, records_ingested_so_far: u64 },
    /// Forwarded `SortPhase1Event::MemoryChunk`. `SortMerge` consumes
    /// these directly as `ChunkSource::Memory` sources in
    /// `MergeDriver::from_slots`.
    MemoryChunk { chunk: Arc<MemoryChunkErased>, records_ingested_so_far: u64 },
    /// Forwarded `SortPhase1Event::AllAnnounced`. Gates `SortMerge`'s
    /// early transition from `WaitingForSetup` to `Merging`.
    AllAnnounced { slot_count: u32, memory_chunk_count: u32, total_records: u64 },
}

impl HeapSize for SortPhase2Event {
    fn heap_size(&self) -> usize {
        match self {
            Self::SpillReady { path, .. } => path.as_os_str().len(),
            Self::MemoryChunk { chunk, .. } => chunk.approx_heap_bytes(),
            Self::AllAnnounced { .. } => 0,
        }
    }
}

impl SortPhase1Event {
    /// Running snapshot of records ingested at the moment this event was
    /// emitted. `SortMerge` takes `max()` over all popped events to recover
    /// the final `total_records` count.
    #[must_use]
    pub fn records_ingested_so_far(&self) -> u64 {
        match self {
            Self::SpillReady { records_ingested_so_far, .. }
            | Self::MemoryChunk { records_ingested_so_far, .. } => *records_ingested_so_far,
            Self::AllAnnounced { total_records, .. } => *total_records,
        }
    }
}

impl SortPhase2Event {
    /// See [`SortPhase1Event::records_ingested_so_far`].
    #[must_use]
    pub fn records_ingested_so_far(&self) -> u64 {
        match self {
            Self::SpillReady { records_ingested_so_far, .. }
            | Self::MemoryChunk { records_ingested_so_far, .. } => *records_ingested_so_far,
            Self::AllAnnounced { total_records, .. } => *total_records,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn memory_chunk_len_and_is_empty() {
        let chunk: MemoryChunkErased = MemoryChunkErased::Coordinate(Vec::new());
        assert_eq!(chunk.len(), 0);
        assert!(chunk.is_empty());

        let chunk = MemoryChunkErased::Coordinate(vec![
            (RawCoordinateKey::default(), RawRecord::from(vec![0xAA; 16])),
            (RawCoordinateKey::default(), RawRecord::from(vec![0xBB; 32])),
        ]);
        assert_eq!(chunk.len(), 2);
        assert!(!chunk.is_empty());
    }

    #[test]
    fn memory_chunk_approx_heap_bytes_counts_overhead_plus_payload() {
        let chunk = MemoryChunkErased::Coordinate(vec![
            (RawCoordinateKey::default(), RawRecord::from(vec![0xAA; 100])),
            (RawCoordinateKey::default(), RawRecord::from(vec![0xBB; 200])),
        ]);
        assert_eq!(chunk.approx_heap_bytes(), 2 * PER_MEMORY_RECORD_OVERHEAD + 300);
    }

    #[test]
    fn records_ingested_so_far_accessors() {
        let dummy_path = std::path::PathBuf::from("/tmp/x");
        let slot =
            Arc::new(SortMergeSlot::new(0, std::io::BufReader::new(tempfile::tempfile().unwrap())));
        let ev1 = SortPhase1Event::SpillReady {
            slot: Arc::clone(&slot),
            path: dummy_path.clone(),
            records_ingested_so_far: 100,
        };
        assert_eq!(ev1.records_ingested_so_far(), 100);

        let chunk = Arc::new(MemoryChunkErased::Coordinate(Vec::new()));
        let ev2 = SortPhase1Event::MemoryChunk {
            chunk: Arc::clone(&chunk),
            records_ingested_so_far: 250,
        };
        assert_eq!(ev2.records_ingested_so_far(), 250);

        let ev3 =
            SortPhase2Event::SpillReady { slot, path: dummy_path, records_ingested_so_far: 300 };
        assert_eq!(ev3.records_ingested_so_far(), 300);

        let ev4 = SortPhase2Event::MemoryChunk { chunk, records_ingested_so_far: 400 };
        assert_eq!(ev4.records_ingested_so_far(), 400);

        let ev5 = SortPhase1Event::AllAnnounced {
            slot_count: 4,
            memory_chunk_count: 1,
            total_records: 500,
        };
        assert_eq!(ev5.records_ingested_so_far(), 500);
        let ev6 = SortPhase2Event::AllAnnounced {
            slot_count: 4,
            memory_chunk_count: 1,
            total_records: 600,
        };
        assert_eq!(ev6.records_ingested_so_far(), 600);
    }

    #[test]
    fn all_announced_heap_size_is_zero() {
        let ev1 = SortPhase1Event::AllAnnounced {
            slot_count: 16,
            memory_chunk_count: 4,
            total_records: 1_000_000,
        };
        assert_eq!(ev1.heap_size(), 0);
        let ev2 = SortPhase2Event::AllAnnounced {
            slot_count: 16,
            memory_chunk_count: 4,
            total_records: 1_000_000,
        };
        assert_eq!(ev2.heap_size(), 0);
    }
}
