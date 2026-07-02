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

impl TemplateMemChunk {
    /// Number of records in the chunk.
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            Self::K24(c) => c.len(),
            Self::Cb32(c) => c.len(),
            Self::Tert32(c) => c.len(),
            Self::K40(c) => c.len(),
        }
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
        match self {
            Self::K24(c) => c.payload_bytes(),
            Self::Cb32(c) => c.payload_bytes(),
            Self::Tert32(c) => c.payload_bytes(),
            Self::K40(c) => c.payload_bytes(),
        }
    }

    /// Borrow the `i`th record's raw BAM body bytes, in this chunk's sorted order.
    ///
    /// # Panics
    ///
    /// Panics if `i >= self.len()`.
    #[must_use]
    pub fn record_bytes(&self, i: usize) -> &[u8] {
        match self {
            Self::K24(c) => c.record_bytes(i),
            Self::Cb32(c) => c.record_bytes(i),
            Self::Tert32(c) => c.record_bytes(i),
            Self::K40(c) => c.record_bytes(i),
        }
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
        match self {
            Self::K24(c) => frame_keyed_record_into(out, c.key_at(i), c.record_bytes(i)),
            Self::Cb32(c) => frame_keyed_record_into(out, c.key_at(i), c.record_bytes(i)),
            Self::Tert32(c) => frame_keyed_record_into(out, c.key_at(i), c.record_bytes(i)),
            Self::K40(c) => frame_keyed_record_into(out, c.key_at(i), c.record_bytes(i)),
        }
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
        match self {
            Self::K24(c) => write_sorted_chunk_inmem(path, codec, compression, c),
            Self::Cb32(c) => write_sorted_chunk_inmem(path, codec, compression, c),
            Self::Tert32(c) => write_sorted_chunk_inmem(path, codec, compression, c),
            Self::K40(c) => write_sorted_chunk_inmem(path, codec, compression, c),
        }
    }
}

/// Build a sorted template chunk from narrow-lane refs pointing into `arena`.
///
/// Featherweight seal (the template analogue of the coordinate
/// [`coordinate_chunk_from_refs`](crate::ref_sort::coordinate_chunk_from_refs)):
/// sorts `refs` on their cached narrow lane key (stable radix — parallel for
/// `sort_threads > 1`, matching the owned path's `par_sort`) and copies that
/// narrow key straight into the returned [`InMemoryChunk<K>`]. NO full-key
/// re-extraction, NO arena body access — the key is already resident in each
/// ref (computed once at `push`). No record bytes are copied: the records
/// reference their bodies in `arena` at `(offset, len)`. Narrow-lane order
/// equals full-key order (every dropped lane is verified constant on the ingest
/// path), so the chunk is correctly ordered for the downstream `MergeDriver<K>`,
/// and spilling/merging the narrow key is byte-identical to the full key.
#[must_use]
pub fn template_chunk_from_arena_refs<K: TemplateLaneKey>(
    arena: Arc<PooledSegmentedBuf>,
    mut refs: Vec<TemplateRecordRef<K>>,
    sort_threads: usize,
) -> InMemoryChunk<K> {
    if sort_threads > 1 {
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
        match self {
            Self::K24(v) => {
                v.push(TemplateRecordRef {
                    key: TemplateLaneKey::from_full(full),
                    offset,
                    len,
                    padding: 0,
                });
            }
            Self::Cb32(v) => {
                v.push(TemplateRecordRef {
                    key: TemplateLaneKey::from_full(full),
                    offset,
                    len,
                    padding: 0,
                });
            }
            Self::Tert32(v) => {
                v.push(TemplateRecordRef {
                    key: TemplateLaneKey::from_full(full),
                    offset,
                    len,
                    padding: 0,
                });
            }
            Self::K40(v) => {
                v.push(TemplateRecordRef {
                    key: TemplateLaneKey::from_full(full),
                    offset,
                    len,
                    padding: 0,
                });
            }
        }
    }

    fn reserve(&mut self, n: usize) {
        match self {
            Self::K24(v) => v.reserve(n),
            Self::Cb32(v) => v.reserve(n),
            Self::Tert32(v) => v.reserve(n),
            Self::K40(v) => v.reserve(n),
        }
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
    /// runs the sort on one worker while the others inflate/compress). Built once
    /// on the first seal (`sort_threads` is only known then) and reused; `None`
    /// on a fresh worker copy. Mirrors the owned `TemplateChunkSorter`'s
    /// `rayon_pool.install(...)`.
    sort_pool: Option<rayon::ThreadPool>,
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
            sort_pool: None,
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
    /// In this increment the chunk is still the full-key [`K40`](TemplateMemChunk::K40)
    /// arm (the full [`TemplateKey`] re-extracted per record), so behavior is
    /// byte-identical to the pre-`TemplateMemChunk` path; the featherweight
    /// narrow-lane arm is a follow-up.
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
        // Bounded pool sized to `sort_threads`, built once and reused. Running the
        // radix + gather inside `pool.install` bounds BOTH to `sort_threads` (via
        // the pool's `current_num_threads`), so they no longer oversubscribe the
        // pipeline's worker pool on a spill.
        let pool = self.sort_pool.get_or_insert_with(|| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(sort_threads.max(1))
                .thread_name(|i| format!("tmpl-sort-{i}"))
                .build()
                .expect("build bounded template-sort rayon pool")
        });
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

    /// A fresh accumulator with the same configuration and empty state — used to
    /// build a worker copy of the arena-front step.
    #[must_use]
    pub fn fresh(&self) -> Self {
        Self {
            lib_lookup: self.lib_lookup.clone(),
            cell_tag: self.cell_tag,
            cb_hasher: self.cb_hasher.clone(),
            key_types: self.key_types,
            header_library_varies: self.header_library_varies,
            state: None,
            pending_reserve: 0,
            sort_pool: None,
        }
    }
}
