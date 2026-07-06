//! Lean in-memory buffering sorter for the P6 `SortBuffer` step.
//!
//! [`CoordinateChunkSorter`] is the buffer-management slice of the former
//! `CoordinateSortStream` (the streaming sort engine retired in P7) with **all
//! the disk and pool machinery removed**: it ingests records into a
//! [`RecordBuffer`], and on
//! demand par-sorts the buffer and materializes it into an owned sorted chunk
//! (`Vec<(RawCoordinateKey, RawRecord)>`). It never writes to disk, owns no
//! `SortWorkerPool`, no temp dirs, and no spill files — compressing and writing
//! a chunk is the `CompressSpill` step's job (P6 Phase-1 split).
//!
//! The driver (`SortBuffer`) calls [`push`](CoordinateChunkSorter::push) per
//! record; when it returns `true` (buffer at the memory limit), the driver takes
//! a sorted chunk via [`take_sorted_chunk`](CoordinateChunkSorter::take_sorted_chunk)
//! and emits it as a spill chunk. At end of input the driver takes one final
//! residual chunk the same way.
//!
//! # Parity
//!
//! `take_sorted_chunk` uses the global stable [`RecordBuffer::par_sort`] plus a
//! parallel `par_iter` materialization — exactly the no-spill residual path in
//! `CoordinateSortStream::into_slot_setup`. The legacy with-spill multi-thread
//! path used `par_sort_into_chunks` to emit *k* sub-runs; emitting one
//! globally-sorted chunk instead is record-order identical (the sub-runs are
//! merged back by key with a lower-source-index tie-break, which for a stable
//! global sort is the same original-insertion order), and it keeps the residual
//! a single merge source so a `Parallel` `CompressSpill` reordering chunks can
//! never perturb the tie-break.

use std::sync::Arc;

use anyhow::Result;
#[cfg(test)]
use fgumi_raw_bam::RawRecord;
use fgumi_raw_bam::SamTag;
use rayon::ThreadPool;

use crate::arena_pool::ArenaPool;
use crate::external::{
    KeyTypesSpec, LibraryLookup, TemplateKeyVariant, dropped_lane_error,
    extract_template_key_inline, select_template_variant, verify_dropped_lanes,
};
use crate::inline::{
    CbKey32, InMemoryChunk, RecordBuffer, TemplateKey, TemplateKey24, TemplateKey40,
    TemplateLaneKey, TemplateRecordBuffer, TertKey32,
};
use crate::keys::RawCoordinateKey;
use crate::memory_probe::force_mi_collect;

/// In-memory coordinate buffering sorter that emits owned sorted chunks.
pub struct CoordinateChunkSorter {
    buffer: RecordBuffer,
    /// Private rayon pool sized to `--threads`, used for `par_sort` and the
    /// parallel chunk materialization.
    rayon_pool: ThreadPool,
    /// Memory-usage threshold (bytes) at which `push` signals a chunk is due.
    memory_limit: usize,
    total_records: u64,
    /// Bounded pool of reusable sort arenas (lever-1 RSS fix). The buffer's
    /// backing `SegmentedBuf` is acquired here and returns (reset-for-reuse)
    /// when the emitted chunk's `Arc` drops, so at most `N_arena` arenas are
    /// ever live. Capacity 1 (default) matches legacy's one-arena-at-a-time
    /// `drain_pending_spill` model.
    arena_pool: Arc<ArenaPool>,
    /// `true` once an arena is installed in `buffer` and not yet drained.
    has_arena: bool,
}

impl CoordinateChunkSorter {
    /// Construct from a pre-sized buffer + rayon pool + arena pool. Built via
    /// [`RawExternalSorter::into_coordinate_chunk_sorter`](crate::RawExternalSorter::into_coordinate_chunk_sorter).
    /// `buffer` starts drained; the first [`ensure_arena`](Self::ensure_arena)
    /// installs one from `arena_pool`.
    pub(crate) fn new(
        buffer: RecordBuffer,
        rayon_pool: ThreadPool,
        memory_limit: usize,
        arena_pool: Arc<ArenaPool>,
    ) -> Self {
        Self { buffer, rayon_pool, memory_limit, total_records: 0, arena_pool, has_arena: false }
    }

    /// Ensure the buffer has a backing arena before a fill. Returns `false`
    /// when all `N_arena` pool arenas are in flight — the caller (`SortBuffer`)
    /// must backpressure (`NoProgress`) and retry once an in-flight chunk is
    /// spilled and its arena returns to the pool. Idempotent while installed.
    #[must_use]
    pub fn ensure_arena(&mut self) -> bool {
        if self.has_arena {
            return true;
        }
        // Only acquire at the START of a fill (drained buffer). Never swap a
        // partially-filled (default) buffer for a pool arena mid-fill — the
        // arena install replaces the backing store, which would orphan the
        // records already buffered. A caller that began filling the default
        // (because the bounded pool was momentarily empty) finishes into it and
        // drains it unpooled.
        if !self.buffer.is_empty() {
            return false;
        }
        if let Some(arena) = self.arena_pool.try_acquire() {
            self.buffer.install_arena(arena);
            self.has_arena = true;
            true
        } else {
            false
        }
    }

    /// Push one raw BAM record into the buffer.
    ///
    /// Returns `true` when the buffer's memory usage has reached the configured
    /// limit and the caller should take a sorted (spill) chunk via
    /// [`Self::take_sorted_chunk`] before pushing more. The caller must have
    /// called [`ensure_arena`](Self::ensure_arena) (returning `true`) first.
    ///
    /// # Errors
    ///
    /// Returns an error if the record cannot be parsed for its coordinate key
    /// (e.g. truncated record).
    pub fn push(&mut self, bam_bytes: &[u8]) -> Result<bool> {
        // Best-effort: acquire a pool arena if one is free. A bounded, RSS-
        // critical driver would call `ensure_arena` explicitly first and
        // backpressure on `false`, so on that path the arena is always installed
        // here. Callers (the `SortBuffer` RecordBatch path, unit tests) that
        // don't pre-check simply fill the buffer's current (possibly pool,
        // possibly default) storage — `take_sorted_chunk` returns it to the pool
        // only when it is a pool arena (`has_arena`).
        let _ = self.ensure_arena();
        // Count only after the record is successfully buffered, so a truncated
        // record (push_coordinate error) never bumps `total_records` for a
        // record that was never stored.
        self.buffer.push_coordinate(bam_bytes)?;
        self.total_records += 1;
        Ok(self.buffer.memory_usage() >= self.memory_limit)
    }

    /// Par-sort the current buffer, materialize it into an owned sorted chunk
    /// (parallel `par_iter`), and clear the buffer for reuse.
    ///
    /// Returns an empty chunk if the buffer is empty (so the caller can skip
    /// emitting an empty chunk).
    ///
    /// **Zero-copy materialisation.** The sorted buffer is drained into an
    /// [`InMemoryChunk`] that *moves* the arena (`Arc<SegmentedBuf>`) and keeps
    /// per-record `(key, offset, len)` references into it — no per-record byte
    /// copy or allocation. (The pre-fix path materialised an owned
    /// `Vec<(K, RawRecord)>`, copying + `malloc`-ing every record and doubling
    /// peak RSS; the A/B smoke campaign flagged that regression vs the legacy
    /// path, which already used this primitive.)
    #[must_use]
    pub fn take_sorted_chunk(&mut self) -> InMemoryChunk<RawCoordinateKey> {
        if self.buffer.is_empty() {
            return InMemoryChunk::empty();
        }
        {
            let buffer = &mut self.buffer;
            let rayon_pool = &self.rayon_pool;
            rayon_pool.install(|| buffer.par_sort());
        }
        // If the buffer is a pool arena, drain into a POOL-MANAGED chunk so the
        // arena returns (reset-for-reuse) when the chunk's last `Arc` drops
        // (after the gather frames it) — this is the reuse + cap-1 bound on the
        // standalone coordinate path. If it is the default storage (a caller
        // that bypassed the pool because it was momentarily empty), drain
        // unpooled so a non-pool buffer never enters the pool. `has_arena` is
        // cleared either way so the next fill re-acquires.
        let chunk = if self.has_arena {
            self.buffer.drain_into_pooled_chunk(&self.arena_pool)
        } else {
            self.buffer.drain_into_single_chunk()
        };
        self.has_arena = false;
        chunk
    }

    /// Total records pushed so far.
    #[must_use]
    pub fn total_records(&self) -> u64 {
        self.total_records
    }

    /// `true` iff the buffer currently holds no records.
    #[must_use]
    pub fn buffer_is_empty(&self) -> bool {
        self.buffer.is_empty()
    }

    /// Construct a minimal `CoordinateChunkSorter` for unit tests.
    ///
    /// Mirrors the construction in
    /// [`RawExternalSorter::into_coordinate_chunk_sorter`](crate::external::RawExternalSorter::into_coordinate_chunk_sorter):
    /// a zero-capacity `RecordBuffer` (the pool hands out the first arena on the
    /// initial `ensure_arena`/`push`), a single-thread rayon pool, and a
    /// capacity-1 `ArenaPool` with the standard segment size.
    ///
    /// # Panics
    ///
    /// Panics if rayon cannot build a single-thread pool (which never occurs in
    /// practice but is the contract of `ThreadPoolBuilder::build`).
    #[cfg(any(test, feature = "test-utils"))]
    #[must_use]
    pub fn for_test(memory_limit: usize, n_ref: u32) -> Self {
        let buffer = crate::inline::RecordBuffer::with_capacity(0, 0, n_ref);
        let rayon_pool =
            rayon::ThreadPoolBuilder::new().num_threads(1).build().expect("rayon pool for test");
        let arena_pool = crate::arena_pool::ArenaPool::new(1, crate::inline::SORT_SEGMENT_SIZE);
        Self::new(buffer, rayon_pool, memory_limit, arena_pool)
    }
}

// ============================================================================
// TemplateChunkSorter — template-coordinate (stable radix, like coordinate)
// ============================================================================

/// The narrowed lane-key buffer the [`TemplateChunkSorter`] sorts on. The
/// variant is chosen on the first `push` from `--key-types` (via
/// [`select_template_variant`]) and matches the batch `RawExternalSorter::sort`
/// dispatch one-for-one: `(cb, tertiary)` → `(false,false)`=`TemplateKey24`
/// (3 lanes), `(true,false)`=`CbKey32`, `(false,true)`=`TertKey32` (4 lanes),
/// `(true,true)`=`TemplateKey40` (5 lanes, the full key). The radix sort runs
/// on the narrower key for speed; the emitted chunk re-widens to the full
/// [`TemplateKey`] (see [`TemplateChunkSorter::take_sorted_chunk`]).
enum TemplateBuffer {
    K24(TemplateRecordBuffer<TemplateKey24>),
    Cb32(TemplateRecordBuffer<CbKey32>),
    Tert32(TemplateRecordBuffer<TertKey32>),
    K40(TemplateRecordBuffer<TemplateKey40>),
}

/// Run `$body` against the inner `TemplateRecordBuffer<K>` of whichever variant
/// `$self` holds, with the buffer bound to `$buf`. `$body` is monomorphized per
/// arm, so type inference resolves `K` (e.g. `TemplateLaneKey::from_full`,
/// `b.push`) from the matched buffer type.
macro_rules! with_template_buffer {
    ($self:expr, $buf:ident => $body:expr) => {
        match $self {
            TemplateBuffer::K24($buf) => $body,
            TemplateBuffer::Cb32($buf) => $body,
            TemplateBuffer::Tert32($buf) => $body,
            TemplateBuffer::K40($buf) => $body,
        }
    };
}

impl TemplateBuffer {
    /// Build the buffer for the chosen `variant`, sized from the shared
    /// (key-width-agnostic) capacity hints. The slight per-`K` ref-size
    /// difference is immaterial for a capacity hint — the buffer grows as
    /// needed — so a single estimate is used across variants.
    fn from_variant(
        variant: TemplateKeyVariant,
        estimated_records: usize,
        estimated_data_bytes: usize,
    ) -> Self {
        match (variant.cb, variant.tertiary) {
            (false, false) => Self::K24(TemplateRecordBuffer::with_capacity(
                estimated_records,
                estimated_data_bytes,
            )),
            (true, false) => Self::Cb32(TemplateRecordBuffer::with_capacity(
                estimated_records,
                estimated_data_bytes,
            )),
            (false, true) => Self::Tert32(TemplateRecordBuffer::with_capacity(
                estimated_records,
                estimated_data_bytes,
            )),
            (true, true) => Self::K40(TemplateRecordBuffer::with_capacity(
                estimated_records,
                estimated_data_bytes,
            )),
        }
    }

    /// Push a record under the narrowed key derived from its full key.
    fn push_full(&mut self, bam_bytes: &[u8], full: &TemplateKey) -> Result<()> {
        with_template_buffer!(self, b => b.push(bam_bytes, TemplateLaneKey::from_full(full)))
    }

    fn memory_usage(&self) -> usize {
        with_template_buffer!(self, b => b.memory_usage())
    }

    fn is_empty(&self) -> bool {
        with_template_buffer!(self, b => b.is_empty())
    }

    fn par_sort(&mut self) {
        with_template_buffer!(self, b => b.par_sort());
    }

    // Test-only: only the owned oracle path (`take_sorted_chunk_owned`) clears
    // the buffer explicitly; the production arena drain empties it in place.
    #[cfg(test)]
    fn clear(&mut self) {
        with_template_buffer!(self, b => b.clear());
    }

    /// Materialize the (already narrow-key-sorted) buffer into an owned chunk
    /// keyed by the **full** [`TemplateKey`], re-extracted per record. Narrow-K
    /// order equals full-key order because every dropped lane is verified
    /// constant, so the chunk is correctly ordered for the downstream merge.
    /// Must be invoked inside the sorter's `rayon_pool.install`.
    ///
    /// Test-only: the owned materialisation is retained purely as the byte-parity
    /// oracle for [`drain_into_full_key_chunk`](Self::drain_into_full_key_chunk),
    /// which is the production path.
    #[cfg(test)]
    fn materialize_full(
        &self,
        lib_lookup: &LibraryLookup,
        cell_tag: Option<SamTag>,
        cb_hasher: &ahash::RandomState,
    ) -> Vec<(TemplateKey, RawRecord)> {
        use rayon::prelude::*;
        with_template_buffer!(self, b => b
            .refs()
            .par_iter()
            .map(|r| {
                let bam_bytes = b.get_record(r);
                let full = extract_template_key_inline(bam_bytes, lib_lookup, cell_tag, cb_hasher);
                (full, RawRecord::from(bam_bytes.to_vec()))
            })
            .collect::<Vec<_>>())
    }

    /// Zero-copy analogue of [`materialize_full`](Self::materialize_full): drain
    /// the (already narrow-key-sorted) buffer into an arena-backed
    /// `InMemoryChunk<TemplateKey>` keyed by the re-extracted full key, WITHOUT
    /// copying record bodies (the arena is moved into the chunk). Produces the
    /// same sorted order, keys, and body bytes as `materialize_full` — the only
    /// difference is the records reference the moved arena instead of owned
    /// `RawRecord`s. Must be invoked inside the sorter's `rayon_pool.install`.
    fn drain_into_full_key_chunk(
        &mut self,
        lib_lookup: &LibraryLookup,
        cell_tag: Option<SamTag>,
        cb_hasher: &ahash::RandomState,
    ) -> InMemoryChunk<TemplateKey> {
        with_template_buffer!(self, b => b.drain_into_full_key_chunk(|body| {
            extract_template_key_inline(body, lib_lookup, cell_tag, cb_hasher)
        }))
    }
}

/// State that exists only once the first record has been seen: the chosen
/// narrowed-key variant, the full key of the first record (the dropped-lane
/// verify baseline), and the matching buffer. Bundling the three behind one
/// `Option` captures their shared "set together on the first push" invariant in
/// the type, so the hot path never unwraps.
struct Provisioned {
    /// Full key of the first record; later records' dropped lanes are verified
    /// against it.
    first_key: TemplateKey,
    /// The chosen narrowed-key variant (`--key-types` + first record).
    variant: TemplateKeyVariant,
    /// The narrowed-key buffer matching `variant`.
    buffer: TemplateBuffer,
}

/// In-memory template-coordinate buffering sorter that emits owned sorted
/// chunks. The template analogue of [`CoordinateChunkSorter`]: a
/// [`TemplateRecordBuffer`] + private rayon pool, with the library / cell-barcode
/// / MI key-extraction state the template key needs. Uses the same global stable
/// radix [`TemplateRecordBuffer::par_sort`], so the single-residual-chunk parity
/// argument is identical to coordinate's.
///
/// **Narrowed radix key (`--key-types`).** Like the batch `RawExternalSorter::sort`
/// path, the in-memory radix sort runs on the `KeyTypesSpec`-narrowed key
/// (`TemplateBuffer` — `TemplateKey24`/`CbKey32`/`TertKey32`/`TemplateKey40`),
/// chosen lazily on the first `push` once the first record is available to seed
/// `select_template_variant`. The narrower key is a **speed** optimization
/// (fewer radix passes). The emitted chunk re-widens to the full [`TemplateKey`]
/// in [`take_sorted_chunk`](Self::take_sorted_chunk) so the buffer-chain protocol
/// (`CompressSpill`/`SortMerge`) stays on the full key.
///
/// Narrowing is order-preserving and validated: `--key-types` may only drop a
/// lane the [dropped-lane validation](verify_dropped_lanes) proves constant, so
/// the narrow-key sort yields the identical order to a full-key sort. This path
/// therefore produces byte-for-byte the same output as the batch path AND rejects
/// the same dropped-lane violations. (The pre-P6 streaming path skipped this
/// validation entirely; restoring it makes `--key-types` actually function in
/// production.)
pub struct TemplateChunkSorter {
    /// Variant + first-key + buffer, provisioned on the first `push` (`None`
    /// until then; an empty stream leaves it `None`).
    state: Option<Provisioned>,
    rayon_pool: ThreadPool,
    memory_limit: usize,
    /// Capacity hints for the lazily-built `TemplateBuffer` (key-width-agnostic).
    estimated_records: usize,
    estimated_data_bytes: usize,
    lib_lookup: LibraryLookup,
    cb_hasher: ahash::RandomState,
    cell_tag: Option<SamTag>,
    /// `--key-types` spec; selects which optional lanes may be dropped.
    key_types: KeyTypesSpec,
    /// Whether the header realizes >1 library ordinal (informs `Auto` selection).
    header_library_varies: bool,
    total_records: u64,
}

impl TemplateChunkSorter {
    /// Built via
    /// [`RawExternalSorter::into_template_chunk_sorter`](crate::RawExternalSorter::into_template_chunk_sorter).
    /// The narrowed-key buffer is deferred to the first `push` (the variant is
    /// chosen from the first record), so this takes capacity hints rather than a
    /// pre-built buffer.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        rayon_pool: ThreadPool,
        memory_limit: usize,
        estimated_records: usize,
        estimated_data_bytes: usize,
        lib_lookup: LibraryLookup,
        cb_hasher: ahash::RandomState,
        cell_tag: Option<SamTag>,
        key_types: KeyTypesSpec,
        header_library_varies: bool,
    ) -> Self {
        Self {
            state: None,
            rayon_pool,
            memory_limit,
            estimated_records,
            estimated_data_bytes,
            lib_lookup,
            cb_hasher,
            cell_tag,
            key_types,
            header_library_varies,
            total_records: 0,
        }
    }

    /// Push one record (extracting the template key); returns `true` when the
    /// buffer has reached the memory limit.
    ///
    /// On the first record the dropped-lane variant is provisioned from
    /// `--key-types`; every record then has its dropped lanes verified against
    /// the first record's, matching the batch `RawExternalSorter::sort` path.
    ///
    /// # Errors
    ///
    /// Returns an error if the record cannot be buffered (e.g. truncated), or if
    /// a record carries a dropped-lane value (CB / MI / library) absent from the
    /// input's first record.
    pub fn push(&mut self, bam_bytes: &[u8]) -> Result<bool> {
        let key = extract_template_key_inline(
            bam_bytes,
            &self.lib_lookup,
            self.cell_tag,
            &self.cb_hasher,
        );

        let memory_usage = if let Some(state) = self.state.as_mut() {
            // Subsequent records: verify the lanes the chosen variant drops are
            // constant relative to the first record, then buffer.
            if let Some(violation) = verify_dropped_lanes(&state.first_key, &key, state.variant) {
                let name = fgumi_raw_bam::RawRecordView::new(bam_bytes).read_name();
                return Err(dropped_lane_error(&String::from_utf8_lossy(name), violation));
            }
            state.buffer.push_full(bam_bytes, &key)?;
            state.buffer.memory_usage()
        } else {
            // First record: provision the variant from --key-types and build the
            // matching narrowed-key buffer. Buffer the record BEFORE committing
            // `self.state`, so a `push_full` error leaves the sorter unprovisioned
            // rather than holding a first_key/variant from a never-buffered record.
            let variant =
                select_template_variant(Some(&key), self.key_types, self.header_library_varies);
            let mut buffer = TemplateBuffer::from_variant(
                variant,
                self.estimated_records,
                self.estimated_data_bytes,
            );
            buffer.push_full(bam_bytes, &key)?;
            let memory_usage = buffer.memory_usage();
            self.state = Some(Provisioned { first_key: key, variant, buffer });
            memory_usage
        };

        // Count only after the record is successfully buffered (past the
        // dropped-lane check and push_full), so an error never bumps
        // `total_records` for a record that was never stored.
        self.total_records += 1;
        Ok(memory_usage >= self.memory_limit)
    }

    /// Par-sort (on the narrowed key) + drain the buffer into one **arena-backed**
    /// `InMemoryChunk<TemplateKey>` keyed by the full [`TemplateKey`] — zero body
    /// copies (the sort arena is moved into the chunk, the full key re-extracted
    /// per record). Symmetric with [`CoordinateChunkSorter::take_sorted_chunk`].
    /// Both `data` and `refs` are cleared by the drain.
    #[must_use]
    pub fn take_sorted_chunk(&mut self) -> InMemoryChunk<TemplateKey> {
        let Some(state) = self.state.as_mut() else {
            return InMemoryChunk::default();
        };
        if state.buffer.is_empty() {
            return InMemoryChunk::default();
        }
        self.rayon_pool.install(|| state.buffer.par_sort());
        let lib_lookup = &self.lib_lookup;
        let cell_tag = self.cell_tag;
        let cb_hasher = &self.cb_hasher;
        let chunk = self
            .rayon_pool
            .install(|| state.buffer.drain_into_full_key_chunk(lib_lookup, cell_tag, cb_hasher));
        force_mi_collect();
        chunk
    }

    /// Test-only owned materialisation — the byte-parity oracle reference for
    /// [`take_sorted_chunk`](Self::take_sorted_chunk). Par-sorts then copies each
    /// record into an owned `Vec<(TemplateKey, RawRecord)>` via `materialize_full`
    /// (the pre-arena behaviour), so a test can assert the arena chunk is
    /// byte-identical to what the owned path produced.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn take_sorted_chunk_owned(&mut self) -> Vec<(TemplateKey, RawRecord)> {
        let Some(state) = self.state.as_mut() else {
            return Vec::new();
        };
        if state.buffer.is_empty() {
            return Vec::new();
        }
        self.rayon_pool.install(|| state.buffer.par_sort());
        let buffer_ref = &state.buffer;
        let lib_lookup = &self.lib_lookup;
        let cell_tag = self.cell_tag;
        let cb_hasher = &self.cb_hasher;
        let chunk = self
            .rayon_pool
            .install(|| buffer_ref.materialize_full(lib_lookup, cell_tag, cb_hasher));
        state.buffer.clear();
        force_mi_collect();
        chunk
    }

    /// Total records pushed so far.
    #[must_use]
    pub fn total_records(&self) -> u64 {
        self.total_records
    }

    /// The narrowed-key variant chosen on the first push (`None` before any
    /// push). Exposed for tests asserting `--key-types` narrowing selects the
    /// expected lane width.
    #[cfg(test)]
    pub(crate) fn chosen_variant(&self) -> Option<TemplateKeyVariant> {
        self.state.as_ref().map(|s| s.variant)
    }
}

#[cfg(test)]
mod tests {
    use crate::RawExternalSorter;
    use crate::inline::InMemoryChunk;
    use crate::keys::{RawCoordinateKey, SortOrder};
    use fgumi_raw_bam::testutil::make_bam_bytes;
    use noodles::sam::Header;

    /// Build `n` coordinate records at descending positions on tid 0, so a
    /// correct sort must reorder them to ascending position.
    fn descending_pos_records(n: usize) -> Vec<Vec<u8>> {
        (0..n)
            .map(|i| {
                let pos = i32::try_from(n - 1 - i).expect("pos fits i32");
                let name = format!("r{i:05}");
                make_bam_bytes(0, pos, 0, name.as_bytes(), &[], 80, -1, -1, &[])
            })
            .collect()
    }

    /// Assert a chunk's keys are non-decreasing (each chunk is a sorted run).
    fn assert_chunk_sorted(chunk: &InMemoryChunk<RawCoordinateKey>) {
        for i in 1..chunk.len() {
            assert!(
                chunk.key_at(i - 1).cmp(chunk.key_at(i)).is_le(),
                "chunk not sorted by coordinate key"
            );
        }
    }

    /// A generous memory limit keeps everything in one residual chunk; `push`
    /// never signals full, and the single chunk is fully sorted.
    #[test]
    fn no_spill_single_sorted_residual_chunk() {
        let records = descending_pos_records(500);
        let mut sorter = RawExternalSorter::new(SortOrder::Coordinate)
            .memory_limit(256 * 1024 * 1024)
            .threads(2)
            .into_coordinate_chunk_sorter(&Header::default())
            .expect("build chunk sorter");

        for r in &records {
            assert!(!sorter.push(r).expect("push"), "no spill expected under a large memory limit");
        }
        let residual = sorter.take_sorted_chunk();
        assert_eq!(residual.len(), records.len(), "residual must hold every record");
        assert_chunk_sorted(&residual);
        assert_eq!(sorter.total_records(), records.len() as u64);
        assert!(sorter.buffer_is_empty(), "buffer cleared after take");
    }

    /// A tiny memory limit forces mid-stream spill chunks; every chunk is a
    /// sorted run and the union covers all records exactly once.
    #[test]
    fn tiny_limit_emits_multiple_sorted_chunks_covering_all_records() {
        let records = descending_pos_records(2_000);
        let mut sorter = RawExternalSorter::new(SortOrder::Coordinate)
            .memory_limit(16 * 1024) // tiny → frequent spills
            .threads(2)
            .into_coordinate_chunk_sorter(&Header::default())
            .expect("build chunk sorter");

        let mut chunks = Vec::new();
        for r in &records {
            if sorter.push(r).expect("push") {
                let chunk = sorter.take_sorted_chunk();
                assert!(!chunk.is_empty());
                assert_chunk_sorted(&chunk);
                chunks.push(chunk);
            }
        }
        let residual = sorter.take_sorted_chunk();
        if !residual.is_empty() {
            assert_chunk_sorted(&residual);
            chunks.push(residual);
        }

        assert!(chunks.len() >= 2, "tiny limit should force at least one spill plus residual");
        // Assert identity, not just count: a spill path that drops record A and
        // duplicates B keeps the total unchanged, so compare the multiset of
        // record bytes across all chunks against the input.
        let mut seen: Vec<Vec<u8>> = chunks
            .iter()
            .flat_map(|c| (0..c.len()).map(move |i| c.record_bytes(i).to_vec()))
            .collect();
        seen.sort();
        let mut expected: Vec<Vec<u8>> = records.clone();
        expected.sort();
        assert_eq!(seen, expected, "chunk union must equal the input set (no drops/dups)");
    }

    /// Equal-key stability: when every record shares one coordinate, the global
    /// stable `par_sort` must emit them in input order. This is the
    /// output-identity-critical tie-break the single-residual design relies on
    /// (vs samtools sort), so pin it directly — not just sortedness.
    #[test]
    fn equal_keys_preserve_input_order_within_chunk() {
        // All at tid 0, pos 0 → identical coordinate key; distinct names make
        // each record's bytes unique so input order is observable.
        let records: Vec<Vec<u8>> = (0..400)
            .map(|i| {
                let name = format!("r{i:05}");
                make_bam_bytes(0, 0, 0, name.as_bytes(), &[], 60, -1, -1, &[])
            })
            .collect();
        let mut sorter = RawExternalSorter::new(SortOrder::Coordinate)
            .memory_limit(256 * 1024 * 1024)
            .threads(2)
            .into_coordinate_chunk_sorter(&Header::default())
            .expect("build chunk sorter");
        for r in &records {
            assert!(!sorter.push(r).expect("push"));
        }
        let chunk = sorter.take_sorted_chunk();
        assert_eq!(chunk.len(), records.len());
        for (i, rec) in records.iter().enumerate() {
            assert_eq!(
                chunk.record_bytes(i),
                rec.as_slice(),
                "equal-key record {i} emerged out of input order"
            );
        }
    }

    /// `--key-types` is honored in the buffer path: dropping the tertiary (MI)
    /// lane and then feeding a record whose MI differs from the first record's is
    /// rejected with the same actionable error the batch `sort()` path produces.
    /// (Pre-P6 the streaming production path skipped this validation entirely.)
    #[test]
    fn template_chunk_sorter_rejects_dropped_lane_violation() {
        use crate::external::KeyTypesSpec;

        // Raw BAM aux for `MI:i:<v>` (tag bytes + 'i' type + i32 LE value).
        fn mi_aux(v: i32) -> Vec<u8> {
            let mut a = vec![b'M', b'I', b'i'];
            a.extend_from_slice(&v.to_le_bytes());
            a
        }

        let mut sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .memory_limit(256 * 1024 * 1024)
            .threads(1)
            .key_types(KeyTypesSpec::None) // drop all optional lanes incl. tertiary/MI
            .into_template_chunk_sorter(&Header::default())
            .expect("build template chunk sorter");

        let r1 = make_bam_bytes(0, 10, 0, b"r1", &[], 40, -1, -1, &mi_aux(1));
        let r2 = make_bam_bytes(0, 10, 0, b"r2", &[], 40, -1, -1, &mi_aux(2));
        assert!(!sorter.push(&r1).expect("first push provisions the variant"));
        let err =
            sorter.push(&r2).expect_err("a differing MI under --key-types none must be rejected");
        let msg = err.to_string();
        assert!(
            msg.contains("MI") && msg.contains("--key-types mi"),
            "unexpected dropped-lane error: {msg}"
        );
    }

    /// With a consistent input (the dropped lanes are constant), `--key-types`
    /// validation passes and the sort proceeds — the common case.
    #[test]
    fn template_chunk_sorter_accepts_constant_dropped_lanes() {
        use crate::external::KeyTypesSpec;
        let mut sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .memory_limit(256 * 1024 * 1024)
            .threads(1)
            .key_types(KeyTypesSpec::None)
            .into_template_chunk_sorter(&Header::default())
            .expect("build");
        // No optional lanes present at all → nothing to violate.
        for i in 0..200i32 {
            let pos = 200 - i;
            let rec = make_bam_bytes(0, pos, 0, format!("r{i}").as_bytes(), &[], 40, -1, -1, &[]);
            sorter.push(&rec).expect("constant (absent) dropped lanes accepted");
        }
        assert_eq!(sorter.take_sorted_chunk().len(), 200);
    }

    /// Plain template records (no cell barcode, no MI, single library) under
    /// `Auto` select the narrowest 3-lane [`TemplateKey24`] variant — confirming
    /// the radix sort runs on the narrowed key, not the full key — while still
    /// producing a correctly ordered chunk keyed by the full [`TemplateKey`].
    #[test]
    fn template_chunk_sorter_selects_narrow_variant_without_cb_or_mi() {
        use crate::external::KeyTypesSpec;
        let mut sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .memory_limit(256 * 1024 * 1024)
            .threads(2)
            .key_types(KeyTypesSpec::Auto)
            .into_template_chunk_sorter(&Header::default())
            .expect("build template chunk sorter");
        assert!(sorter.chosen_variant().is_none(), "no variant before the first push");

        // Descending positions on tid 0 → a correct sort reorders to ascending.
        for i in 0..300i32 {
            let pos = 300 - i;
            let rec =
                make_bam_bytes(0, pos, 0, format!("r{i:05}").as_bytes(), &[], 40, -1, -1, &[]);
            sorter.push(&rec).expect("plain record accepted");
        }

        let variant = sorter.chosen_variant().expect("variant provisioned on first push");
        assert_eq!(
            variant.lanes(),
            3,
            "no cb / no MI / single library must narrow to the 3-lane TemplateKey24"
        );

        let chunk = sorter.take_sorted_chunk_owned();
        assert_eq!(chunk.len(), 300, "every record retained");
        for w in chunk.windows(2) {
            assert!(
                w[0].0.cmp(&w[1].0).is_le(),
                "narrow-key sort must yield a full-key-ordered chunk"
            );
        }
    }

    /// Output-identity: the narrowed-key (3-lane `TemplateKey24`) sort must emit
    /// byte-identical record order to the full-key (5-lane `TemplateKey40`) sort,
    /// **including equal-key ties** — template-coordinate order is
    /// output-identity-critical, so pin identity, not just sortedness. The input
    /// deliberately contains records sharing a full template key (same tid / pos /
    /// name) but differing in a non-key field (mapq), so a stable sort must keep
    /// each tie group in input order under both key widths.
    #[test]
    fn template_chunk_sorter_narrow_matches_full_key_order() {
        use crate::external::KeyTypesSpec;

        let mut records: Vec<Vec<u8>> = Vec::new();
        for i in 0..120i32 {
            let pos = (i % 8) + 1; // 8 distinct positions → coordinate-level ties
            let name = format!("t{:03}", i % 30); // 30 names → full-key ties within a pos
            for seq_len in [20usize, 40usize] {
                // Same (tid, pos, name) but distinct seq length → identical
                // template key, distinct bytes: an observable equal-key tie.
                records.push(make_bam_bytes(0, pos, 0, name.as_bytes(), &[], seq_len, -1, -1, &[]));
            }
        }

        let sort_with = |key_types: KeyTypesSpec| -> (usize, Vec<Vec<u8>>) {
            let mut sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
                .memory_limit(256 * 1024 * 1024)
                .threads(2)
                .key_types(key_types)
                .into_template_chunk_sorter(&Header::default())
                .expect("build template chunk sorter");
            for r in &records {
                sorter.push(r).expect("push");
            }
            let lanes = sorter.chosen_variant().expect("variant provisioned").lanes();
            let out = sorter
                .take_sorted_chunk_owned()
                .into_iter()
                .map(|(_, rec)| rec.as_ref().to_vec())
                .collect();
            (lanes, out)
        };

        let (narrow_lanes, narrow) = sort_with(KeyTypesSpec::Auto);
        let (full_lanes, full) = sort_with(KeyTypesSpec::Full);
        assert_eq!(narrow_lanes, 3, "Auto over cb/MI-free data must narrow to TemplateKey24");
        assert_eq!(full_lanes, 5, "Full must force the 5-lane TemplateKey40");
        assert_eq!(narrow.len(), records.len(), "every record retained");
        assert_eq!(
            narrow, full,
            "narrowed-key sort must emit byte-identical record order (incl. equal-key ties) \
             to the full-key sort"
        );
    }

    /// Byte-parity oracle: the arena-backed production drain
    /// (`take_sorted_chunk`, zero body copies) must produce byte-for-byte the
    /// same sorted (full-key, body) sequence as the owned `materialize_full`
    /// path (`take_sorted_chunk_owned`). Same records + the same deterministic stable
    /// radix on two identical sorters ⇒ identical output, so any divergence in
    /// the arena offset math, full-key re-extraction, or ordering fails here.
    /// Exercises equal-key ties (same tid/pos/name, distinct seq length) under
    /// both the narrowed (`Auto` → `TemplateKey24`) and full (`Full` →
    /// `TemplateKey40`) key widths.
    #[test]
    fn template_arena_drain_matches_owned_materialize() {
        use crate::external::KeyTypesSpec;

        let mut records: Vec<Vec<u8>> = Vec::new();
        for i in 0..120i32 {
            let pos = (i % 8) + 1;
            let name = format!("t{:03}", i % 30);
            for seq_len in [20usize, 40usize] {
                records.push(make_bam_bytes(0, pos, 0, name.as_bytes(), &[], seq_len, -1, -1, &[]));
            }
        }

        let build = |key_types: KeyTypesSpec| {
            RawExternalSorter::new(SortOrder::TemplateCoordinate)
                .memory_limit(256 * 1024 * 1024)
                .threads(2)
                .key_types(key_types)
                .into_template_chunk_sorter(&Header::default())
                .expect("build template chunk sorter")
        };

        for key_types in [KeyTypesSpec::Auto, KeyTypesSpec::Full] {
            let mut owned_sorter = build(key_types);
            let mut arena_sorter = build(key_types);
            for r in &records {
                owned_sorter.push(r).expect("push owned");
                arena_sorter.push(r).expect("push arena");
            }
            let owned = owned_sorter.take_sorted_chunk_owned();
            let arena = arena_sorter.take_sorted_chunk();

            assert_eq!(arena.len(), owned.len(), "record count must match ({key_types:?})");
            assert_eq!(owned.len(), records.len(), "every record retained ({key_types:?})");
            for (i, (key, rec)) in owned.iter().enumerate() {
                assert_eq!(
                    arena.key_at(i),
                    key,
                    "full key at {i} must match owned path ({key_types:?})"
                );
                assert_eq!(
                    arena.record_bytes(i),
                    rec.as_ref(),
                    "record bytes at {i} must be byte-identical to owned path ({key_types:?})"
                );
            }
        }
    }

    /// Byte-parity: the arena-front `TemplateArenaAccumulator` must produce the
    /// same sorted (full-key, body) sequence as the owned `TemplateChunkSorter`
    /// — same library/CB provisioning, same `--key-types` narrowed-lane radix,
    /// same full-key re-extraction — for a mix of tied and distinct keys, under
    /// both the narrowed (`Auto` → `TemplateKey24`) and full (`Full` →
    /// `TemplateKey40`) widths. This pins the arena front (Inc 3) to the legacy
    /// path before it is wired into the pipeline.
    #[test]
    #[allow(unsafe_code)]
    fn template_arena_accumulator_matches_owned_sorter() {
        use crate::arena_pool::PooledSegmentedBuf;
        use crate::external::KeyTypesSpec;
        use crate::segmented_buf::SegmentedBuf;
        use crate::template_arena::TemplateArenaAccumulator;
        use std::sync::Arc;

        let mut records: Vec<Vec<u8>> = Vec::new();
        for i in 0..120i32 {
            let pos = (i % 8) + 1;
            let name = format!("t{:03}", i % 30);
            for seq_len in [20usize, 40usize] {
                records.push(make_bam_bytes(0, pos, 0, name.as_bytes(), &[], seq_len, -1, -1, &[]));
            }
        }

        for key_types in [KeyTypesSpec::Auto, KeyTypesSpec::Full] {
            // ---- Owned oracle ----
            let mut owned = RawExternalSorter::new(SortOrder::TemplateCoordinate)
                .memory_limit(256 * 1024 * 1024)
                .threads(2)
                .key_types(key_types)
                .into_template_chunk_sorter(&Header::default())
                .expect("build owned template sorter");
            for r in &records {
                owned.push(r).expect("owned push");
            }
            let owned_chunk = owned.take_sorted_chunk_owned();

            // ---- Arena accumulator (bodies resident in a shared arena) ----
            let mut arena = SegmentedBuf::with_capacity(0, 1 << 20);
            arena.reserve_full_capacity();
            let mut acc =
                TemplateArenaAccumulator::from_header(&Header::default(), None, key_types);
            for r in &records {
                // SAFETY: each slot is fully written (copy_from_slice) before any read.
                let off = unsafe { arena.grow_uninit(r.len()) };
                unsafe { arena.slice_mut(off, r.len()) }.copy_from_slice(r);
                acc.push(r, off as u64, u32::try_from(r.len()).unwrap()).expect("arena push");
            }
            let arena_arc = Arc::new(PooledSegmentedBuf::unpooled(arena));
            let arena_chunk = acc.seal(arena_arc, 2);

            assert_eq!(arena_chunk.len(), owned_chunk.len(), "record count ({key_types:?})");
            assert_eq!(owned_chunk.len(), records.len(), "every record retained ({key_types:?})");
            // Byte-parity is the invariant: the arena chunk's records, in sorted
            // order, must be byte-identical to the owned oracle's. The arena chunk
            // now carries the chosen narrow lane key (variant-agnostic here), and
            // narrow-lane order equals full-key order, so record order matches the
            // owned full-key sort exactly.
            for (i, (_key, rec)) in owned_chunk.iter().enumerate() {
                assert_eq!(
                    arena_chunk.record_bytes(i),
                    rec.as_ref(),
                    "record body at {i} must be byte-identical to owned path ({key_types:?})"
                );
            }
        }
    }

    /// An empty template stream never builds a buffer and yields an empty
    /// residual (the lazy `Option<TemplateBuffer>` stays `None`).
    #[test]
    fn template_chunk_sorter_empty_stream_yields_empty_residual() {
        use crate::external::KeyTypesSpec;
        let mut sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .memory_limit(256 * 1024 * 1024)
            .threads(1)
            .key_types(KeyTypesSpec::Auto)
            .into_template_chunk_sorter(&Header::default())
            .expect("build template chunk sorter");
        assert!(sorter.take_sorted_chunk().is_empty());
        assert_eq!(sorter.total_records(), 0);
        assert!(sorter.chosen_variant().is_none(), "no push → no variant");
    }

    /// An empty stream yields an empty residual (the caller emits no chunk).
    #[test]
    fn empty_stream_yields_empty_residual() {
        let mut sorter = RawExternalSorter::new(SortOrder::Coordinate)
            .memory_limit(256 * 1024 * 1024)
            .threads(1)
            .into_coordinate_chunk_sorter(&Header::default())
            .expect("build chunk sorter");
        assert!(sorter.take_sorted_chunk().is_empty());
        assert_eq!(sorter.total_records(), 0);
    }
}
