//! Template-coordinate arena-front sort: build a sorted
//! [`InMemoryChunk<TemplateKey>`] from record bodies already resident in a shared
//! arena, WITHOUT copying record bytes — the template analogue of the coordinate
//! [`coordinate_chunk_from_refs`](crate::ref_sort::coordinate_chunk_from_refs).
//!
//! [`TemplateArenaAccumulator`] encapsulates everything the template key needs
//! that the arena-front pipeline step (`FindBoundariesAndSort` in
//! `fgumi-pipeline-io`) does not have: the [`LibraryLookup`], cell-barcode tag +
//! hasher, the `--key-types` narrowed-lane selection, and the per-record
//! dropped-lane verification. The pipeline step calls [`push`](TemplateArenaAccumulator::push)
//! per record (arena offset + len) and [`seal`](TemplateArenaAccumulator::seal)
//! per run, exactly mirroring the owned [`TemplateChunkSorter`](crate::TemplateChunkSorter)
//! so the output is byte-for-byte identical — only the record bodies stay in the
//! shared arena instead of being copied into owned `RawRecord`s.

use std::io;
use std::path::Path;
use std::sync::Arc;

use anyhow::Result;
use noodles::sam::Header;

use crate::arena_pool::PooledSegmentedBuf;
use crate::external::{
    KeyTypesSpec, LibraryLookup, TemplateKeyVariant, cb_hasher, dropped_lane_error,
    extract_template_key_inline, select_template_variant, verify_dropped_lanes,
};
use crate::inline::{
    CbKey32, InMemoryChunk, TemplateKey, TemplateKey24, TemplateKey40, TemplateLaneKey,
    TemplateRecordRef, TertKey32, parallel_radix_sort_template_refs, radix_sort_template_refs,
};
use crate::ref_sort::PARALLEL_SORT_THRESHOLD;
use crate::{SpillCodec, frame_keyed_record_into, write_sorted_chunk_inmem};
use fgumi_raw_bam::{RawRecordView, SamTag};

/// Variant-carrying erased template residual chunk: one arm per `--key-types`
/// narrowed lane. Lets `MemoryChunkErased::TemplateCoordinate` hold whichever
/// lane variant the sort chose, so template-coordinate rides its natural narrow
/// key end-to-end through merge and spill — exactly like every other sort order
/// rides its own `K` — instead of being pinned to the full 40-byte
/// [`TemplateKey`]. Narrow-lane order equals full-key order (every dropped lane
/// is verified constant on the ingest path), so merging and spilling the narrow
/// key is byte-identical to the full key. The [`K40`](Self::K40) arm is the full
/// key (all lanes), used by the legacy owned path and the full variant.
pub enum TemplateMemChunk {
    /// 24-byte core-only lane (neither cb nor tertiary optional word present).
    K24(InMemoryChunk<TemplateKey24>),
    /// 32-byte lane whose optional word carries `cb_hash`.
    Cb32(InMemoryChunk<CbKey32>),
    /// 32-byte lane whose optional word carries the tertiary (library<<48 | mi).
    Tert32(InMemoryChunk<TertKey32>),
    /// Full 40-byte key (all lanes) — the legacy owned path and full variant.
    K40(InMemoryChunk<TemplateKey>),
}

/// Run `$body` against the inner `InMemoryChunk<K>` of whichever variant
/// `$self` holds, with the chunk bound to `$chunk`. `$body` is monomorphized per
/// arm, so type inference resolves `K` (e.g. `c.key_at`) from the matched chunk
/// type. Mirrors `chunk_sorter::with_template_buffer!`.
///
/// Every per-variant dispatch on this enum goes through here: writing the
/// four-arm match out by hand in each method makes a fifth lane variant a
/// six-site edit, and a copy-paste that dispatches the wrong arm still compiles.
macro_rules! with_template_chunk {
    ($self:expr, $chunk:ident => $body:expr) => {
        match $self {
            TemplateMemChunk::K24($chunk) => $body,
            TemplateMemChunk::Cb32($chunk) => $body,
            TemplateMemChunk::Tert32($chunk) => $body,
            TemplateMemChunk::K40($chunk) => $body,
        }
    };
}

impl TemplateMemChunk {
    /// Number of records in the chunk.
    #[must_use]
    pub fn len(&self) -> usize {
        with_template_chunk!(self, c => c.len())
    }

    /// `true` iff the chunk holds zero records.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Total record-payload bytes (sum of record lengths; excludes keys and
    /// index overhead). Used for byte-budget accounting at the chunk boundary.
    #[must_use]
    pub fn payload_bytes(&self) -> usize {
        with_template_chunk!(self, c => c.payload_bytes())
    }

    /// Borrow the `i`th record's raw BAM body bytes, in this chunk's sorted order.
    ///
    /// # Panics
    ///
    /// Panics if `i >= self.len()`.
    #[must_use]
    pub fn record_bytes(&self, i: usize) -> &[u8] {
        with_template_chunk!(self, c => c.record_bytes(i))
    }

    /// Frame the `i`th record into `out` in the spill layout
    /// `[key][u32 LE len][record]`, using this chunk's narrow-lane key.
    ///
    /// Encapsulates the per-variant key dispatch so the pipeline-io spill
    /// serializer stays variant-agnostic.
    ///
    /// # Errors
    ///
    /// Returns an error if writing to `out` fails.
    pub fn frame_record_into(&self, i: usize, out: &mut Vec<u8>) -> io::Result<()> {
        with_template_chunk!(self, c => frame_keyed_record_into(out, c.key_at(i), c.record_bytes(i)))
    }

    /// Write the whole chunk to a spill file at `path` via
    /// [`write_sorted_chunk_inmem`], keyed by this chunk's narrow lane.
    ///
    /// Encapsulates the per-variant key dispatch so the pipeline-io compress
    /// step stays variant-agnostic.
    ///
    /// # Errors
    ///
    /// Returns an error if the spill write fails.
    pub fn write_spill(&self, path: &Path, codec: SpillCodec, compression: u32) -> Result<()> {
        with_template_chunk!(self, c => write_sorted_chunk_inmem(path, codec, compression, c))
    }
}

/// Build a sorted template chunk from narrow-lane refs pointing into `arena`.
///
/// Featherweight seal (the template analogue of the coordinate
/// [`coordinate_chunk_from_refs`](crate::ref_sort::coordinate_chunk_from_refs)):
/// sorts `refs` on their cached narrow lane key (stable radix — parallel for
/// large multi-threaded runs, matching the owned path's `par_sort`) and copies
/// that narrow key straight into the returned [`InMemoryChunk<K>`]. NO full-key
/// re-extraction, NO arena body access — the key is already resident in each
/// ref (computed once at `push`). No record bytes are copied: the records
/// reference their bodies in `arena` at `(offset, len)`. Narrow-lane order
/// equals full-key order (every dropped lane is verified constant on the ingest
/// path), so the chunk is correctly ordered for the downstream `MergeDriver<K>`,
/// and spilling/merging the narrow key is byte-identical to the full key.
///
/// The parallel-vs-serial decision mirrors the coordinate front
/// ([`sort_coordinate_refs`](crate::ref_sort)): it uses the parallel radix only
/// when `sort_threads > 1` AND the run exceeds the shared
/// `PARALLEL_SORT_THRESHOLD`, so a small run does not pay the parallel radix's
/// partition/coordination overhead. Both paths produce byte-identical output
/// (the parallel template radix is stability-tested against the serial one).
#[must_use]
pub fn template_chunk_from_arena_refs<K: TemplateLaneKey>(
    arena: Arc<PooledSegmentedBuf>,
    mut refs: Vec<TemplateRecordRef<K>>,
    sort_threads: usize,
) -> InMemoryChunk<K> {
    if sort_threads > 1 && refs.len() >= PARALLEL_SORT_THRESHOLD {
        parallel_radix_sort_template_refs(&mut refs);
    } else {
        radix_sort_template_refs(&mut refs);
    }
    let records: Vec<(K, u64, u32)> = refs.into_iter().map(|r| (r.key, r.offset, r.len)).collect();
    InMemoryChunk::from_parts(arena, records)
}

/// Narrow-lane ref accumulator, one arm per `--key-types` variant (mirrors the
/// owned `TemplateBuffer`'s variant set). Each holds `TemplateRecordRef<K>`s
/// pointing at record bodies in the shared inflate arena — no owned byte storage.
enum ArenaRefs {
    K24(Vec<TemplateRecordRef<TemplateKey24>>),
    Cb32(Vec<TemplateRecordRef<CbKey32>>),
    Tert32(Vec<TemplateRecordRef<TertKey32>>),
    K40(Vec<TemplateRecordRef<TemplateKey40>>),
}

/// [`with_template_chunk!`]'s analogue for [`ArenaRefs`]: run `$body` against
/// the inner `Vec<TemplateRecordRef<K>>` of whichever variant `$self` holds.
macro_rules! with_arena_refs {
    ($self:expr, $refs:ident => $body:expr) => {
        match $self {
            ArenaRefs::K24($refs) => $body,
            ArenaRefs::Cb32($refs) => $body,
            ArenaRefs::Tert32($refs) => $body,
            ArenaRefs::K40($refs) => $body,
        }
    };
}

impl ArenaRefs {
    fn for_variant(v: TemplateKeyVariant) -> Self {
        match (v.cb, v.tertiary) {
            (false, false) => Self::K24(Vec::new()),
            (true, false) => Self::Cb32(Vec::new()),
            (false, true) => Self::Tert32(Vec::new()),
            (true, true) => Self::K40(Vec::new()),
        }
    }

    #[inline]
    fn push(&mut self, full: &TemplateKey, offset: u64, len: u32) {
        // The `key` field type resolves `K` per arm, so `from_full` narrows to the
        // matching lane width (identical to the owned `TemplateRecordBuffer::push`).
        with_arena_refs!(self, v => v.push(TemplateRecordRef {
            key: TemplateLaneKey::from_full(full),
            offset,
            len,
            padding: 0,
        }));
    }

    fn reserve(&mut self, n: usize) {
        with_arena_refs!(self, v => v.reserve(n));
    }
}

/// State that exists only once the first record has been seen: the chosen
/// narrowed-key variant, the first record's full key (the dropped-lane verify
/// baseline), and the accumulated refs. Set together on the first push and
/// persisted across runs (spills) — matching the owned `TemplateChunkSorter`,
/// whose variant/baseline are chosen once and never re-selected.
struct AccState {
    first_key: TemplateKey,
    variant: TemplateKeyVariant,
    refs: ArenaRefs,
}

/// Arena-front template-coordinate accumulator: the template analogue of the
/// owned [`TemplateChunkSorter`](crate::TemplateChunkSorter), accumulating
/// arena-pointing refs instead of copying records. Produces byte-identical
/// output (same provisioning, same dropped-lane rejection, same sorted order).
pub struct TemplateArenaAccumulator {
    lib_lookup: LibraryLookup,
    cell_tag: Option<SamTag>,
    cb_hasher: ahash::RandomState,
    key_types: KeyTypesSpec,
    header_library_varies: bool,
    /// Variant + baseline + refs; `None` until the first record provisions it.
    state: Option<AccState>,
    /// Reserve hint received before the first record (variant unknown), applied
    /// when the ref buffer is provisioned.
    pending_reserve: usize,
    /// Bounded rayon pool (sized to `sort_threads`) that [`seal`](Self::seal)
    /// installs the per-run radix + full-key gather into, so they run on exactly
    /// `sort_threads` threads instead of the GLOBAL rayon pool. Without this the
    /// parallel radix and `par_iter` gather fan out over every core and
    /// oversubscribe the pipeline's own worker pool during a spill (the pipeline
    /// runs the sort on one worker while the others inflate/compress). Mirrors
    /// the owned `TemplateChunkSorter`'s `rayon_pool.install(...)`.
    ///
    /// **Shared across worker copies, and it has to be.** `sort_threads` is only
    /// known at the first [`seal`](Self::seal), which happens *after*
    /// [`fresh`](Self::fresh) has already produced the copies — so a plain
    /// `Option<ThreadPool>`, or even an `Option<Arc<ThreadPool>>`, gives every
    /// worker its own pool and puts `workers × sort_threads` threads on the box.
    /// That is precisely the oversubscription this field exists to prevent, just
    /// reintroduced one level up. `Arc<OnceLock<_>>` lets the copies be made
    /// first and the single pool be built once, by whichever worker seals first.
    sort_pool: Arc<std::sync::OnceLock<rayon::ThreadPool>>,
}

impl TemplateArenaAccumulator {
    /// Build an accumulator from the BAM `header`, the sort's cell-barcode tag,
    /// and the `--key-types` spec, deriving the library lookup and CB hasher
    /// exactly as
    /// [`RawExternalSorter::into_template_chunk_sorter`](crate::RawExternalSorter::into_template_chunk_sorter).
    #[must_use]
    pub fn from_header(header: &Header, cell_tag: Option<SamTag>, key_types: KeyTypesSpec) -> Self {
        let lib_lookup = LibraryLookup::from_header(header);
        let header_library_varies = lib_lookup.distinct_header_ordinals() > 1;
        Self {
            lib_lookup,
            cell_tag,
            cb_hasher: cb_hasher(),
            key_types,
            header_library_varies,
            state: None,
            pending_reserve: 0,
            sort_pool: Arc::new(std::sync::OnceLock::new()),
        }
    }

    /// Reserve capacity for approximately `est_records` refs for the current run.
    pub fn reserve(&mut self, est_records: usize) {
        if let Some(state) = self.state.as_mut() {
            state.refs.reserve(est_records);
        } else {
            self.pending_reserve = self.pending_reserve.max(est_records);
        }
    }

    /// Extract the template key from `body` (the record's BAM body, `block_size`
    /// prefix excluded, at arena offset `body_off`, length `len`), provision the
    /// narrowed-lane variant on the first record, verify the dropped lanes on
    /// every subsequent record, and accumulate a ref into the arena.
    ///
    /// # Errors
    ///
    /// Returns an error if a record carries a dropped-lane value (CB / MI /
    /// library) absent from the first record — the same rejection the owned
    /// `TemplateChunkSorter::push` performs.
    pub fn push(&mut self, body: &[u8], body_off: u64, len: u32) -> Result<()> {
        let key =
            extract_template_key_inline(body, &self.lib_lookup, self.cell_tag, &self.cb_hasher);
        if let Some(state) = self.state.as_mut() {
            if let Some(violation) = verify_dropped_lanes(&state.first_key, &key, state.variant) {
                let name = RawRecordView::new(body).read_name();
                return Err(dropped_lane_error(&String::from_utf8_lossy(name), violation));
            }
            state.refs.push(&key, body_off, len);
        } else {
            let variant =
                select_template_variant(Some(&key), self.key_types, self.header_library_varies);
            let mut refs = ArenaRefs::for_variant(variant);
            if self.pending_reserve > 0 {
                refs.reserve(self.pending_reserve);
            }
            refs.push(&key, body_off, len);
            self.state = Some(AccState { first_key: key, variant, refs });
        }
        Ok(())
    }

    /// Sort the refs accumulated for the current run and drain them into an
    /// arena-backed [`TemplateMemChunk`] (zero body copies). The chosen variant +
    /// baseline are RETAINED for subsequent runs (spills), matching the owned
    /// sorter. Empty if nothing was pushed.
    ///
    /// Featherweight seal: each arm emits the variant-matching narrow
    /// [`TemplateMemChunk`] lane using the key already resident in the ref — no
    /// per-record re-extraction and no body access — so the downstream
    /// merge/spill ride the narrow key lane chosen for this run.
    ///
    /// # `sort_threads` is per-sort, not per-call
    ///
    /// The signature takes it per call, which reads as though each call sizes its
    /// own execution. It does not: the value builds the shared bounded pool on
    /// the FIRST seal across the whole worker fan-out, and every later call —
    /// this accumulator's next run, or another worker copy — runs on that pool at
    /// its original width, silently ignoring a different value. Pass the same
    /// number every time; see `sort_pool`.
    ///
    /// # Panics
    ///
    /// Panics if the bounded sort rayon pool cannot be built (an infrastructure
    /// failure, e.g. the OS refuses the `sort_threads` worker threads).
    #[must_use]
    pub fn seal(
        &mut self,
        arena: Arc<PooledSegmentedBuf>,
        sort_threads: usize,
    ) -> TemplateMemChunk {
        // Drain this run's refs, leaving a fresh empty accumulator of the SAME
        // variant so the next run re-uses the once-chosen variant + baseline.
        // (Scoped so the `self.state` borrow ends before we touch `sort_pool`.)
        let refs = {
            let Some(state) = self.state.as_mut() else {
                return TemplateMemChunk::K40(InMemoryChunk::default());
            };
            std::mem::replace(&mut state.refs, ArenaRefs::for_variant(state.variant))
        };
        let pool = self.bounded_sort_pool(sort_threads);
        // Featherweight seal: each arm produces the variant-matching narrow chunk
        // (the key is already in the ref — no re-extraction, no body access), and
        // tags it with the chosen variant so the downstream merge/spill ride the
        // narrow lane. The variant is global, so every run's chunk shares one arm.
        pool.install(move || match refs {
            ArenaRefs::K24(r) => {
                TemplateMemChunk::K24(template_chunk_from_arena_refs(arena, r, sort_threads))
            }
            ArenaRefs::Cb32(r) => {
                TemplateMemChunk::Cb32(template_chunk_from_arena_refs(arena, r, sort_threads))
            }
            ArenaRefs::Tert32(r) => {
                TemplateMemChunk::Tert32(template_chunk_from_arena_refs(arena, r, sort_threads))
            }
            ArenaRefs::K40(r) => {
                TemplateMemChunk::K40(template_chunk_from_arena_refs(arena, r, sort_threads))
            }
        })
    }

    /// The shared bounded pool, built on first use at `sort_threads` threads.
    ///
    /// Running the radix + gather inside `pool.install` bounds BOTH to
    /// `sort_threads` (via the pool's `current_num_threads`), so they do not
    /// oversubscribe the pipeline's worker pool on a spill. Built once for the
    /// whole fan-out — see the `sort_pool` field.
    ///
    /// **First call wins.** `get_or_init` builds the pool once, so a later call
    /// with a different `sort_threads` silently gets the existing pool at the
    /// original width. That is deliberate — one pool for the fan-out is the whole
    /// point — but it means `sort_threads` is a property of the *sort*, not of an
    /// individual call.
    fn bounded_sort_pool(&self, sort_threads: usize) -> &rayon::ThreadPool {
        self.sort_pool.get_or_init(|| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(sort_threads.max(1))
                .thread_name(|i| format!("tmpl-sort-{i}"))
                .build()
                .expect("build bounded template-sort rayon pool")
        })
    }

    /// The shared bounded pool, for tests that assert copies share one (the
    /// sharing is otherwise unobservable — `seal` returns early before touching
    /// the pool when the accumulator was never provisioned).
    #[cfg(test)]
    pub fn sort_pool_for_test(&mut self, sort_threads: usize) -> &rayon::ThreadPool {
        self.bounded_sort_pool(sort_threads)
    }

    /// A worker copy: same configuration and provisioning, empty refs.
    ///
    /// **Provisioning is per-sort, not per-worker, and this carries it across.**
    /// The chosen lane variant and the dropped-lane baseline (`first_key`) are
    /// selected once from the first record and must then be identical for every
    /// worker in the sort. Resetting them here instead would let each `Auto`
    /// worker select a lane from *its own* first record, so one sort could emit
    /// chunks of two different key widths; would give each worker a different
    /// baseline, so a lane constant within every worker but varying across them
    /// would pass `verify_dropped_lanes`; and would seal a worker that received
    /// no records as [`TemplateMemChunk::K40`] no matter what the others chose.
    ///
    /// The refs are NOT carried — each worker accumulates its own — so a copy
    /// starts empty but already knows which lane it is filling.
    ///
    /// # Preconditions
    ///
    /// Call this only **after** the parent has been provisioned (i.e. after its
    /// first [`push`](Self::push)). Cloning a worker from an unprovisioned
    /// parent yields `state: None`, and that worker will provision itself
    /// independently — the exact divergence above. There is no production
    /// caller yet; the pipeline step that fans workers out arrives with
    /// `fgumi-pipeline-io`, and it owns honouring this.
    #[must_use]
    pub fn fresh(&self) -> Self {
        Self {
            lib_lookup: self.lib_lookup.clone(),
            cell_tag: self.cell_tag,
            cb_hasher: self.cb_hasher.clone(),
            key_types: self.key_types,
            header_library_varies: self.header_library_varies,
            // Carry variant + baseline, not refs: same lane, own accumulation.
            state: self.state.as_ref().map(|s| AccState {
                first_key: s.first_key,
                variant: s.variant,
                refs: ArenaRefs::for_variant(s.variant),
            }),
            pending_reserve: 0,
            // Share the holder, not a new one: see the field's doc.
            sort_pool: Arc::clone(&self.sort_pool),
        }
    }
}
