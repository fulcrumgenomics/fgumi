//! Typed-event protocol between the three sort steps in the runall-sort
//! fused chain.

use std::path::PathBuf;
use std::sync::Arc;

use fgumi_sort::{
    InMemoryChunk, RawCoordinateKey, RawQuerynameKey, RawQuerynameLexKey, SortMergeSlot,
    TemplateMemChunk,
};

use fgumi_pipeline_core::item::HeapSize;

/// Approximate fixed per-record index overhead of an arena-backed
/// [`InMemoryChunk<K>`], in bytes, added once per record in
/// [`MemoryChunkErased::approx_heap_bytes`] on top of the variable record payload
/// (`InMemoryChunk::payload_bytes`).
///
/// This covers one `(K, offset, len)` index slot — the sort key `K` (largest
/// variant: `TemplateKey`) plus the offset/len into the shared buffer — together
/// with allocator-bucket slack the payload byte count does not capture. It is an
/// intentionally conservative constant, not a `size_of` expression, so the queue
/// accounting over-counts rather than under-counts memory.
const PER_MEMORY_RECORD_OVERHEAD: usize = 354;

/// In-memory sorted residual chunk produced by the Phase-1 sort head
/// (`SortBuffer` or `FindBoundariesAndSort`), type-erased over the sort-key
/// variant `K`.
pub enum MemoryChunkErased {
    /// Coordinate-sort residual. `K = RawCoordinateKey`. Zero-copy arena-backed
    /// chunk (shares the sort buffer's `Arc<SegmentedBuf>`).
    Coordinate(InMemoryChunk<RawCoordinateKey>),
    /// Queryname-sort residual, lexicographic comparator. `K = RawQuerynameLexKey`.
    /// Arena-backed like [`Coordinate`](Self::Coordinate): the record bodies live
    /// in the chunk's shared buffer; the key owns its (small) name bytes.
    QuerynameLex(InMemoryChunk<RawQuerynameLexKey>),
    /// Queryname-sort residual, natural comparator. `K = RawQuerynameKey`.
    QuerynameNatural(InMemoryChunk<RawQuerynameKey>),
    /// Template-coordinate-sort residual, carried as a variant-tagged
    /// [`TemplateMemChunk`] so the chosen `--key-types` narrow lane rides through
    /// merge and spill like every other order's `K`. Zero-copy arena-backed chunk
    /// (shares the sort buffer's `Arc<SegmentedBuf>`), like
    /// [`Coordinate`](Self::Coordinate).
    TemplateCoordinate(TemplateMemChunk),
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

    /// Borrow the `i`th record's raw BAM body bytes, in this chunk's sorted order.
    ///
    /// Zero-copy for every variant — each is arena-backed and slices its shared
    /// `SegmentedBuf`. Used by `SortMerge`'s single-source fast path to gather a
    /// fully-sorted in-memory chunk into output blocks without a (k = 1)
    /// loser-tree merge.
    ///
    /// # Panics
    ///
    /// Panics if `i >= self.len()`.
    #[must_use]
    pub fn record_bytes(&self, i: usize) -> &[u8] {
        match self {
            Self::Coordinate(v) => v.record_bytes(i),
            Self::QuerynameLex(v) => v.record_bytes(i),
            Self::QuerynameNatural(v) => v.record_bytes(i),
            Self::TemplateCoordinate(v) => v.record_bytes(i),
        }
    }

    /// Approximate heap footprint in bytes.
    #[must_use]
    pub fn approx_heap_bytes(&self) -> usize {
        let (count, payload): (usize, usize) = match self {
            Self::Coordinate(v) => (v.len(), v.payload_bytes()),
            Self::QuerynameLex(v) => (v.len(), v.payload_bytes()),
            Self::QuerynameNatural(v) => (v.len(), v.payload_bytes()),
            Self::TemplateCoordinate(v) => (v.len(), v.payload_bytes()),
        };
        count * PER_MEMORY_RECORD_OVERHEAD + payload
    }
}

/// Events from `SortBuffer` → `CompressSpill` (the P6 Phase-1 split).
///
/// `SortBuffer` (Serial) sorts each filled buffer and emits the sorted records
/// as a chunk; `CompressSpill` (Parallel) then compresses spill chunks to disk
/// or passes the in-memory residual through. The seam carries already-sorted
/// `MemoryChunkErased`s, so it is byte-bounded by [`HeapSize`] just like the
/// downstream `SortPhase1Event` queue.
pub enum SortChunkEvent {
    /// A sorted chunk to be compressed and written to disk by `CompressSpill`.
    ///
    /// `seq` is the chunk's logical spill index, assigned monotonically by
    /// `SortBuffer`. `CompressSpill` uses it as the opened slot's `file_id` so
    /// the merge tie-break order matches the legacy spill order regardless of
    /// which Parallel worker writes the file.
    Spill { seq: u32, chunk: MemoryChunkErased, records_ingested_so_far: u64 },
    /// A sorted in-memory residual chunk to pass straight through as a
    /// `SortPhase1Event::MemoryChunk` — no disk round-trip (the fast path).
    Residual { chunk: MemoryChunkErased, records_ingested_so_far: u64 },
    /// Terminal sentinel carrying the final counts (number of `Spill` chunks =
    /// `slot_count`, number of `Residual` chunks = `memory_chunk_count`).
    /// `CompressSpill` forwards it verbatim as `SortPhase1Event::AllAnnounced`.
    AllAnnounced { slot_count: u32, memory_chunk_count: u32, total_records: u64 },
}

impl HeapSize for SortChunkEvent {
    fn heap_size(&self) -> usize {
        // Mirror `SortPhase1Event::heap_size`: a fixed per-event base so the
        // byte-bounded transport queue cannot absorb an unbounded count of
        // near-zero-cost control events (`AllAnnounced`).
        let base = std::mem::size_of::<Self>();
        match self {
            Self::Spill { chunk, .. } | Self::Residual { chunk, .. } => {
                base + chunk.approx_heap_bytes()
            }
            Self::AllAnnounced { .. } => base,
        }
    }
}

impl SortChunkEvent {
    /// Running snapshot of records ingested at the moment this event was emitted.
    #[must_use]
    pub fn records_ingested_so_far(&self) -> u64 {
        match self {
            Self::Spill { records_ingested_so_far, .. }
            | Self::Residual { records_ingested_so_far, .. } => *records_ingested_so_far,
            Self::AllAnnounced { total_records, .. } => *total_records,
        }
    }
}

/// Events from `SpillGather` → `SpillBlockCompress` → `SpillWrite` (the
/// block-parallel spill-write split that replaces the monolithic single-worker
/// `CompressSpill`).
///
/// `SpillGather` (Serial) fans each `SortChunkEvent::Spill` chunk into
/// record-aligned raw [`Block`](Self::Block)s and forwards `Residual` /
/// `AllAnnounced`. `SpillBlockCompress` (Parallel) compresses each `Block`'s `bytes`
/// in place. `SpillWrite` (Serial) demultiplexes blocks back to per-`file_id`
/// spill files and emits the existing [`SortPhase1Event`].
///
/// **Ordinal contract:** `SpillGather` mints `ordinal` monotonically across
/// **every** emitted item (every block of every file, plus the `Residual` /
/// `AllAnnounced` passthroughs), so the stream is dense and gap-free for the
/// framework's single-cursor `ByItemOrdinal` reorder. Because `SortBuffer`
/// (Serial) emits `Spill` events one-at-a-time in `seq` order and
/// `SpillGather` (Serial) drains them in order, each file's blocks are
/// **contiguous** in the ordinal stream — so `SpillWrite` only ever has one
/// spill file open at a time.
pub enum SpillBlockEvent {
    /// One raw (pre-compression) or compressed (post-`SpillBlockCompress`) block of a
    /// spill file. `file_id` is the spill `seq` (the eventual slot `file_id`),
    /// `is_last_in_file` marks the final block so `SpillWrite` can finalize the
    /// file and emit `SpillReady`. `SpillWrite` (Serial) owns the path allocation,
    /// so the block carries no path — only the `file_id` that names the file.
    Block {
        ordinal: u64,
        file_id: u32,
        is_last_in_file: bool,
        records_ingested_so_far: u64,
        bytes: Vec<u8>,
    },
    /// A sorted in-memory residual chunk, passed straight through to
    /// `SortPhase1Event::MemoryChunk` (no disk round-trip).
    Residual { ordinal: u64, chunk: MemoryChunkErased, records_ingested_so_far: u64 },
    /// Terminal sentinel forwarded verbatim as `SortPhase1Event::AllAnnounced`.
    AllAnnounced { ordinal: u64, slot_count: u32, memory_chunk_count: u32, total_records: u64 },
}

impl SpillBlockEvent {
    /// The dense ordinal that drives `ByItemOrdinal` reordering into `SpillWrite`.
    #[must_use]
    pub fn ordinal(&self) -> u64 {
        match self {
            Self::Block { ordinal, .. }
            | Self::Residual { ordinal, .. }
            | Self::AllAnnounced { ordinal, .. } => *ordinal,
        }
    }
}

impl HeapSize for SpillBlockEvent {
    fn heap_size(&self) -> usize {
        // Fixed per-event base (so the byte-bounded queue can't absorb an
        // unbounded count of near-empty control events) plus the variable
        // payload: a block's bytes, or a residual chunk's records.
        let base = std::mem::size_of::<Self>();
        match self {
            Self::Block { bytes, .. } => base + bytes.capacity(),
            Self::Residual { chunk, .. } => base + chunk.approx_heap_bytes(),
            Self::AllAnnounced { .. } => base,
        }
    }
}

impl fgumi_pipeline_core::item::Ordered for SpillBlockEvent {
    fn ordinal(&self) -> u64 {
        SpillBlockEvent::ordinal(self)
    }
}

/// Events from the Phase-1 producers (`CompressSpill` / `SpillWrite`) →
/// `SortSpillDecompress`.
pub enum SortPhase1Event {
    /// A spill chunk file has been closed and is ready for Phase 2 decompression.
    SpillReady { slot: Arc<SortMergeSlot>, path: PathBuf, records_ingested_so_far: u64 },
    /// A par-sorted residual in-memory chunk passed through by the Phase-1
    /// producer (`CompressSpill` / `SpillWrite`) on its drained-completion path.
    MemoryChunk { chunk: Arc<MemoryChunkErased>, records_ingested_so_far: u64 },
    /// Sentinel emitted as the LAST event by the Phase-1 producer
    /// (`CompressSpill` / `SpillWrite`) on its drained-completion path.
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

    /// Build the arena-backed coordinate chunk the `Coordinate` variant carries,
    /// from raw record payloads (keys are `default()`; irrelevant to these tests).
    fn coord(payloads: Vec<Vec<u8>>) -> InMemoryChunk<RawCoordinateKey> {
        InMemoryChunk::from_owned_records(
            payloads.into_iter().map(|b| (RawCoordinateKey::default(), b)).collect(),
        )
    }

    #[test]
    fn memory_chunk_len_and_is_empty() {
        let chunk: MemoryChunkErased = MemoryChunkErased::Coordinate(coord(Vec::new()));
        assert_eq!(chunk.len(), 0);
        assert!(chunk.is_empty());

        let chunk = MemoryChunkErased::Coordinate(coord(vec![vec![0xAA; 16], vec![0xBB; 32]]));
        assert_eq!(chunk.len(), 2);
        assert!(!chunk.is_empty());
    }

    #[test]
    fn memory_chunk_approx_heap_bytes_counts_overhead_plus_payload() {
        let chunk = MemoryChunkErased::Coordinate(coord(vec![vec![0xAA; 100], vec![0xBB; 200]]));
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

        let chunk = Arc::new(MemoryChunkErased::Coordinate(coord(Vec::new())));
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
    fn sort_chunk_event_heap_size_and_accessors() {
        let chunk = MemoryChunkErased::Coordinate(coord(vec![vec![0xAA; 100], vec![0xBB; 200]]));
        let payload = chunk.approx_heap_bytes();

        let spill = SortChunkEvent::Spill { seq: 3, chunk, records_ingested_so_far: 42 };
        assert_eq!(spill.heap_size(), std::mem::size_of::<SortChunkEvent>() + payload);
        assert_eq!(spill.records_ingested_so_far(), 42);

        let residual = SortChunkEvent::Residual {
            chunk: MemoryChunkErased::Coordinate(coord(Vec::new())),
            records_ingested_so_far: 7,
        };
        assert_eq!(residual.records_ingested_so_far(), 7);

        // Control events carry no heap payload but still cost a fixed base so the
        // byte-bounded queue cannot accept an unbounded count of them.
        let announced = SortChunkEvent::AllAnnounced {
            slot_count: 4,
            memory_chunk_count: 1,
            total_records: 500,
        };
        assert_eq!(announced.heap_size(), std::mem::size_of::<SortChunkEvent>());
        assert_eq!(announced.records_ingested_so_far(), 500);
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
