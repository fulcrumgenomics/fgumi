//! Typed-event protocol between the three sort steps in the runall-sort
//! fused chain.

use std::path::PathBuf;
use std::sync::Arc;

use fgumi_raw_bam::RawRecord;
use fgumi_sort::{
    RawCoordinateKey, RawQuerynameKey, RawQuerynameLexKey, SortMergeSlot, TemplateKey,
};

use fgumi_pipeline_core::item::HeapSize;

/// Approximate fixed per-entry overhead of a `Vec<(K, RawRecord)>`, in bytes,
/// added once per record in [`MemoryChunkErased::approx_heap_bytes`] on top of
/// the variable record payload (`RawRecord::as_ref().len()`).
///
/// This is the size of one `(K, RawRecord)` tuple slot — the sort key `K`
/// (largest variant: `TemplateKey`) inline plus the `RawRecord` handle — together
/// with the per-`RawRecord` heap-allocation and allocator-bucket slack that the
/// payload byte count does not capture. It is an intentionally conservative
/// constant, not a `size_of` expression, so the queue accounting over-counts
/// rather than under-counts memory. If the key or record types grow materially,
/// re-derive it from `size_of::<(TemplateKey, RawRecord)>()` plus allocator slack.
const PER_MEMORY_RECORD_OVERHEAD: usize = 354;

/// In-memory sorted residual chunk produced by `SortAndSpill`, type-erased
/// over the sort-key variant `K`.
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

    /// Approximate heap footprint in bytes.
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
    /// A spill chunk file has been closed and is ready for Phase 2 decompression.
    SpillReady { slot: Arc<SortMergeSlot>, path: PathBuf, records_ingested_so_far: u64 },
    /// A par-sorted residual in-memory chunk produced by `SortAndSpill`'s
    /// drained-completion path.
    MemoryChunk { chunk: Arc<MemoryChunkErased>, records_ingested_so_far: u64 },
    /// Sentinel emitted as the LAST event by `SortAndSpill`'s drained-completion path.
    AllAnnounced { slot_count: u32, memory_chunk_count: u32, total_records: u64 },
}

impl HeapSize for SortPhase1Event {
    fn heap_size(&self) -> usize {
        // Charge a fixed per-event base so the byte-bounded transport queues
        // cannot accept an unbounded count of near-zero-cost control events
        // (`SpillReady` with an empty path, `AllAnnounced`). Memory stays a
        // function of configuration rather than event count.
        let base = std::mem::size_of::<Self>();
        match self {
            Self::SpillReady { path, .. } => base + path.as_os_str().len(),
            Self::MemoryChunk { chunk, .. } => base + chunk.approx_heap_bytes(),
            Self::AllAnnounced { .. } => base,
        }
    }
}

/// Events from `SortSpillDecompress` → `SortMerge`.
pub enum SortPhase2Event {
    /// Forwarded `SortPhase1Event::SpillReady`.
    SpillReady { slot: Arc<SortMergeSlot>, path: PathBuf, records_ingested_so_far: u64 },
    /// Forwarded `SortPhase1Event::MemoryChunk`.
    MemoryChunk { chunk: Arc<MemoryChunkErased>, records_ingested_so_far: u64 },
    /// Forwarded `SortPhase1Event::AllAnnounced`.
    AllAnnounced { slot_count: u32, memory_chunk_count: u32, total_records: u64 },
}

impl HeapSize for SortPhase2Event {
    fn heap_size(&self) -> usize {
        // See `SortPhase1Event::heap_size`: a fixed per-event base keeps the
        // byte-bounded queues from absorbing unbounded control-event counts.
        let base = std::mem::size_of::<Self>();
        match self {
            Self::SpillReady { path, .. } => base + path.as_os_str().len(),
            Self::MemoryChunk { chunk, .. } => base + chunk.approx_heap_bytes(),
            Self::AllAnnounced { .. } => base,
        }
    }
}

impl SortPhase1Event {
    /// Running snapshot of records ingested at the moment this event was emitted.
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
        let slot = Arc::new(SortMergeSlot::new(
            0,
            std::io::BufReader::new(tempfile::tempfile().unwrap()),
            fgumi_sort::SpillCodec::Bgzf,
        ));
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
    fn all_announced_heap_size_charges_base_cost() {
        // Control events carry no heap payload but must still cost a fixed,
        // non-zero amount so the byte-bounded queues cannot accept an unbounded
        // count of them.
        let ev1 = SortPhase1Event::AllAnnounced {
            slot_count: 16,
            memory_chunk_count: 4,
            total_records: 1_000_000,
        };
        assert_eq!(ev1.heap_size(), std::mem::size_of::<SortPhase1Event>());
        let ev2 = SortPhase2Event::AllAnnounced {
            slot_count: 16,
            memory_chunk_count: 4,
            total_records: 1_000_000,
        };
        assert_eq!(ev2.heap_size(), std::mem::size_of::<SortPhase2Event>());
    }
}
