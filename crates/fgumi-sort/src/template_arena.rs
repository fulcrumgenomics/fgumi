//! Template-coordinate arena-front sort: build a sorted
//! [`InMemoryChunk<TemplateKey>`] from record bodies already resident in a shared
//! arena, WITHOUT copying record bytes — the template analogue of the coordinate
//! [`coordinate_chunk_from_refs`](crate::ref_sort::coordinate_chunk_from_refs).
//!
//! Two fronts share one extraction ([`TemplateKeyContext::extract_one`]), one
//! ref container ([`TemplateArenaRefs`]) and one seal ([`seal_template_refs`]):
//!
//! - [`TemplateArenaAccumulator`] is the per-record front. The arena-front
//!   pipeline step (`FindBoundariesAndSort` in `fgumi-pipeline-io`) calls
//!   [`push`](TemplateArenaAccumulator::push) per record (arena offset + len)
//!   and [`seal`](TemplateArenaAccumulator::seal) per run; the accumulator
//!   provisions the [`TemplateKeyContext`] from the first record and calls
//!   `extract_one` for every later one.
//! - [`TemplateKeyContext::extract_batch_erased`] is the batched front: it
//!   re-walks a range of `[u32 len][body]*` frames, and its
//!   [`ExtractedRefs`] are appended to a [`TemplateArenaRefs`] that
//!   [`seal_template_refs`] sorts.
//!
//! Both mirror the owned [`TemplateChunkSorter`](crate::TemplateChunkSorter) —
//! same library/CB provisioning, same `--key-types` lane selection, same
//! dropped-lane rejection — so the output is byte-for-byte identical; only the
//! record bodies stay in the shared arena instead of being copied into owned
//! `RawRecord`s.
//!
//! [`TemplateKeyContext`]: crate::template_key::TemplateKeyContext
//! [`TemplateKeyContext::extract_one`]: crate::template_key::TemplateKeyContext::extract_one
//! [`TemplateKeyContext::extract_batch_erased`]: crate::template_key::TemplateKeyContext::extract_batch_erased

use std::io;
use std::path::Path;
use std::sync::Arc;

use anyhow::Result;
use noodles::sam::Header;

use crate::arena_pool::PooledSegmentedBuf;
use crate::external::{KeyTypesSpec, TemplateKeyVariant, dropped_lane_error, violation_read_name};
use crate::inline::{
    CbKey32, InMemoryChunk, TemplateKey, TemplateKey24, TemplateKey40, TemplateLaneKey,
    TemplateRecordRef, TertKey32, parallel_radix_sort_template_refs, radix_sort_template_refs,
};
use crate::run_bound::RunBound;
use crate::sort_pool::BoundedSortPool;
use crate::template_key::{
    ExtractedRefs, TemplateKeyContext, TemplateKeySeed, lane_variant, with_lane_refs,
};
use crate::{SpillCodec, frame_keyed_record_into, write_sorted_chunk_inmem};
use fgumi_raw_bam::{MIN_BAM_RECORD_LEN, SamTag};

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

    /// The `i`th record's body length in bytes, WITHOUT touching the shared data
    /// buffer (reads only the per-record `len` index). See
    /// [`InMemoryChunk::record_len`](crate::InMemoryChunk::record_len).
    ///
    /// # Panics
    ///
    /// Panics if `i >= self.len()`.
    #[must_use]
    pub fn record_len(&self, i: usize) -> u32 {
        with_template_chunk!(self, c => c.record_len(i))
    }

    /// The chunk's minimum sort key (`key_at(0)`) as an owned [`RunBound`], or
    /// `None` if the chunk is empty. The variant matches this chunk's narrowed
    /// template lane. Records are pre-sorted, so this is `O(1)`.
    #[must_use]
    pub fn min_bound(&self) -> Option<RunBound> {
        if self.is_empty() {
            return None;
        }
        Some(match self {
            TemplateMemChunk::K24(c) => RunBound::TemplateK24(*c.key_at(0)),
            TemplateMemChunk::Cb32(c) => RunBound::TemplateCb32(*c.key_at(0)),
            TemplateMemChunk::Tert32(c) => RunBound::TemplateTert32(*c.key_at(0)),
            TemplateMemChunk::K40(c) => RunBound::TemplateK40(*c.key_at(0)),
        })
    }

    /// The chunk's maximum sort key (`key_at(len - 1)`) as an owned [`RunBound`],
    /// or `None` if the chunk is empty. The variant matches this chunk's
    /// narrowed template lane. Records are pre-sorted, so this is `O(1)`.
    #[must_use]
    pub fn max_bound(&self) -> Option<RunBound> {
        let last = self.len().checked_sub(1)?;
        Some(match self {
            TemplateMemChunk::K24(c) => RunBound::TemplateK24(*c.key_at(last)),
            TemplateMemChunk::Cb32(c) => RunBound::TemplateCb32(*c.key_at(last)),
            TemplateMemChunk::Tert32(c) => RunBound::TemplateTert32(*c.key_at(last)),
            TemplateMemChunk::K40(c) => RunBound::TemplateK40(*c.key_at(last)),
        })
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
    if crate::ref_sort::radix_sorts_in_parallel(refs.len(), sort_threads) {
        parallel_radix_sort_template_refs(&mut refs);
    } else {
        radix_sort_template_refs(&mut refs);
    }
    let records: Vec<(K, u64, u32)> = refs.into_iter().map(|r| (r.key, r.offset, r.len)).collect();
    InMemoryChunk::from_parts(arena, records)
}

/// Narrow-lane ref accumulator, one arm per `--key-types` variant. Each arm holds
/// `TemplateRecordRef<K>`s pointing at record bodies in the shared inflate arena.
/// Public: the per-record [`TemplateArenaAccumulator`] and any batched caller of
/// [`TemplateKeyContext::extract_batch_erased`](crate::template_key::TemplateKeyContext::extract_batch_erased)
/// share one container and one seal ([`seal_template_refs`]).
#[derive(Debug)]
pub enum TemplateArenaRefs {
    /// 24-byte lane.
    K24(Vec<TemplateRecordRef<TemplateKey24>>),
    /// 32-byte lane carrying `cb_hash`.
    Cb32(Vec<TemplateRecordRef<CbKey32>>),
    /// 32-byte lane carrying the tertiary word.
    Tert32(Vec<TemplateRecordRef<TertKey32>>),
    /// Full 40-byte key.
    K40(Vec<TemplateRecordRef<TemplateKey40>>),
}

impl TemplateArenaRefs {
    /// An empty container of the arm matching `v`.
    #[must_use]
    pub fn for_variant(v: TemplateKeyVariant) -> Self {
        match (v.cb, v.tertiary) {
            (false, false) => Self::K24(Vec::new()),
            (true, false) => Self::Cb32(Vec::new()),
            (false, true) => Self::Tert32(Vec::new()),
            (true, true) => Self::K40(Vec::new()),
        }
    }

    /// The variant this container's arm encodes.
    #[must_use]
    pub fn variant(&self) -> TemplateKeyVariant {
        lane_variant!(Self, self)
    }

    /// Reserve room for `n` more refs.
    pub fn reserve(&mut self, n: usize) {
        with_lane_refs!(Self, self, v => v.reserve(n));
    }

    /// Number of refs held.
    #[must_use]
    pub fn len(&self) -> usize {
        with_lane_refs!(Self, self, v => v.len())
    }

    /// `true` iff no refs are held.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Append one extracted batch, in order. Append only a batch that
    /// reported no violation: one that did fails the run instead (see
    /// [`TemplateKeyContext::extract_batch`](crate::template_key::TemplateKeyContext::extract_batch)).
    ///
    /// # Panics
    ///
    /// Panics if `batch`'s arm differs from this container's — the variant is
    /// chosen once per sort, so a mismatch is a wiring bug.
    pub fn append(&mut self, batch: ExtractedRefs) {
        match (self, batch) {
            (Self::K24(v), ExtractedRefs::K24(b)) => v.extend_from_slice(&b),
            (Self::Cb32(v), ExtractedRefs::Cb32(b)) => v.extend_from_slice(&b),
            (Self::Tert32(v), ExtractedRefs::Tert32(b)) => v.extend_from_slice(&b),
            (Self::K40(v), ExtractedRefs::K40(b)) => v.extend_from_slice(&b),
            (me, b) => panic!(
                "TemplateArenaRefs::append: container arm {} cannot take a {} batch \
                 (the lane variant is per-sort; this is a wiring bug)",
                arm_name(me.variant()),
                arm_name(b.variant()),
            ),
        }
    }

    /// Push one ref from a FULL key, narrowed to this arm's lane (the `key`
    /// field's type resolves `K` per arm, identical to the owned
    /// `TemplateRecordBuffer::push`).
    #[inline]
    pub fn push_full(&mut self, full: &TemplateKey, body_off: u64, len: u32) {
        with_lane_refs!(Self, self, v => v.push(TemplateRecordRef {
            key: TemplateLaneKey::from_full(full),
            offset: body_off,
            len,
            padding: 0,
        }));
    }
}

/// The arm name for a variant, for wiring-bug panics.
fn arm_name(v: TemplateKeyVariant) -> &'static str {
    match (v.cb, v.tertiary) {
        (false, false) => "K24",
        (true, false) => "Cb32",
        (false, true) => "Tert32",
        (true, true) => "K40",
    }
}

/// Sort the run's refs into the matching narrow-lane [`TemplateMemChunk`] arm
/// (featherweight: the key is already in each ref; no body access) inside
/// `pool.install`. `sort_threads` is the phase-1 width; `pool` is the shared
/// bounded pool built at that width. When `--sort-threads` caps phase 1, the
/// caller holds the whole cap for the duration of this call, so pool steps and
/// this sort never run together beyond the cap. `refs` is drained, leaving an
/// empty container of the SAME arm for the next run (the variant is chosen once
/// per sort).
#[must_use]
pub fn seal_template_refs(
    refs: &mut TemplateArenaRefs,
    arena: Arc<PooledSegmentedBuf>,
    sort_threads: usize,
    pool: &rayon::ThreadPool,
) -> TemplateMemChunk {
    let drained = std::mem::replace(refs, TemplateArenaRefs::for_variant(refs.variant()));
    pool.install(move || match drained {
        TemplateArenaRefs::K24(r) => {
            TemplateMemChunk::K24(template_chunk_from_arena_refs(arena, r, sort_threads))
        }
        TemplateArenaRefs::Cb32(r) => {
            TemplateMemChunk::Cb32(template_chunk_from_arena_refs(arena, r, sort_threads))
        }
        TemplateArenaRefs::Tert32(r) => {
            TemplateMemChunk::Tert32(template_chunk_from_arena_refs(arena, r, sort_threads))
        }
        TemplateArenaRefs::K40(r) => {
            TemplateMemChunk::K40(template_chunk_from_arena_refs(arena, r, sort_threads))
        }
    })
}

/// Arena-front template-coordinate accumulator: the template analogue of the
/// owned [`TemplateChunkSorter`](crate::TemplateChunkSorter), accumulating
/// arena-pointing refs instead of copying records. Produces byte-identical
/// output (same provisioning, same dropped-lane rejection, same sorted order).
pub struct TemplateArenaAccumulator {
    seed: TemplateKeySeed,
    /// Variant + baseline and this run's refs; `None` until the first record
    /// provisions the context, then carried (context only) by `fresh` —
    /// provisioning is per sort, not per worker. One `Option` so "refs exist
    /// exactly when the context does" holds by construction.
    state: Option<(Arc<TemplateKeyContext>, TemplateArenaRefs)>,
    /// Reserve hint received before the first record (variant unknown), applied
    /// when the ref buffer is provisioned.
    pending_reserve: usize,
    /// See [`BoundedSortPool`]: shared across `fresh` copies.
    sort_pool: BoundedSortPool,
}

impl TemplateArenaAccumulator {
    /// Build an accumulator from the BAM `header`, the sort's cell-barcode tag,
    /// and the `--key-types` spec, deriving the library lookup and CB hasher
    /// exactly as
    /// [`RawExternalSorter::into_template_chunk_sorter`](crate::RawExternalSorter::into_template_chunk_sorter).
    #[must_use]
    pub fn from_header(header: &Header, cell_tag: Option<SamTag>, key_types: KeyTypesSpec) -> Self {
        Self {
            seed: TemplateKeySeed::from_header(header, cell_tag, key_types),
            state: None,
            pending_reserve: 0,
            sort_pool: BoundedSortPool::new("tmpl-sort"),
        }
    }

    /// Refs accumulated for the current run (what the next seal sorts).
    #[must_use]
    pub fn pending_len(&self) -> usize {
        self.state.as_ref().map_or(0, |(_, refs)| refs.len())
    }

    /// Reserve capacity for approximately `est_records` refs for the current run.
    pub fn reserve(&mut self, est_records: usize) {
        match self.state.as_mut() {
            Some((_, refs)) => refs.reserve(est_records),
            None => self.pending_reserve = self.pending_reserve.max(est_records),
        }
    }

    /// Extract the template key from `body` (the record's BAM body, `block_size`
    /// prefix excluded, at arena offset `body_off`, length `len`), provision the
    /// narrowed-lane variant on the first record, verify the dropped lanes on
    /// every subsequent record, and accumulate a ref into the arena.
    ///
    /// # Errors
    ///
    /// Returns an error if `body` is shorter than the 32 fixed BAM bytes (it
    /// cannot be a record, and the key extractor indexes those fields), or if a
    /// record carries a dropped-lane value (CB / MI / library) absent from the
    /// first record — the same rejection the owned `TemplateChunkSorter::push`
    /// performs.
    pub fn push(&mut self, body: &[u8], body_off: u64, len: u32) -> Result<()> {
        if body.len() < MIN_BAM_RECORD_LEN {
            anyhow::bail!(
                "record at arena offset {body_off}: body of {} bytes is shorter than a BAM \
                 record ({MIN_BAM_RECORD_LEN} bytes)",
                body.len()
            );
        }
        let Some((ctx, refs)) = self.state.as_mut() else {
            let ctx = self.seed.provision(body);
            let mut refs = TemplateArenaRefs::for_variant(ctx.variant);
            if self.pending_reserve > 0 {
                refs.reserve(self.pending_reserve);
            }
            refs.push_full(&ctx.first_key, body_off, len);
            self.state = Some((ctx, refs));
            return Ok(());
        };
        let (key, violation) = ctx.extract_one(body);
        if let Some(v) = violation {
            return Err(dropped_lane_error(&violation_read_name(body), v));
        }
        refs.push_full(&key, body_off, len);
        Ok(())
    }

    /// Sort the refs accumulated for the current run and drain them into an
    /// arena-backed [`TemplateMemChunk`] (zero body copies) via
    /// [`seal_template_refs`]. The chosen variant + baseline are RETAINED for
    /// subsequent runs (spills), matching the owned sorter. Empty if nothing was
    /// pushed.
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
    /// number every time; see [`BoundedSortPool`].
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
        let Some((_, refs)) = self.state.as_mut() else {
            return TemplateMemChunk::K40(InMemoryChunk::default());
        };
        let pool = self.sort_pool.bounded(sort_threads);
        seal_template_refs(refs, arena, sort_threads, pool)
    }

    /// The shared bounded pool, for tests that assert copies share one (the
    /// sharing is otherwise unobservable — `seal` returns early before touching
    /// the pool when the accumulator was never provisioned).
    #[cfg(test)]
    pub fn sort_pool_for_test(&mut self, sort_threads: usize) -> &rayon::ThreadPool {
        self.sort_pool.bounded(sort_threads)
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
    /// independently — the exact divergence above. The caller that fans workers
    /// out owns honouring this.
    #[must_use]
    pub fn fresh(&self) -> Self {
        Self {
            seed: self.seed.clone(),
            // Carry the context, not the refs: same lane, own accumulation.
            state: self
                .state
                .as_ref()
                .map(|(ctx, _)| (Arc::clone(ctx), TemplateArenaRefs::for_variant(ctx.variant))),
            pending_reserve: 0,
            // Share the holder, not a new one: see `BoundedSortPool`.
            sort_pool: self.sort_pool.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A `TemplateKey24` with the given `primary` lane (the leading comparison
    /// field); other lanes zeroed. Enough to order bounds in these tests.
    fn k24(primary: u64) -> TemplateKey24 {
        TemplateKey24 { primary, secondary: 0, name_hash_upper: 0 }
    }

    /// Build a K24 lane chunk from pre-sorted `(primary, body)` records (the seal
    /// path always produces sorted chunks; `from_owned_records` preserves order).
    fn k24_chunk(records: Vec<(u64, Vec<u8>)>) -> TemplateMemChunk {
        TemplateMemChunk::K24(InMemoryChunk::from_owned_records(
            records.into_iter().map(|(p, b)| (k24(p), b)).collect(),
        ))
    }

    #[test]
    fn empty_template_chunk_has_no_bounds() {
        let chunk = k24_chunk(Vec::new());
        assert_eq!(chunk.min_bound(), None);
        assert_eq!(chunk.max_bound(), None);
    }

    #[test]
    fn template_bounds_are_first_and_last_keys() {
        let chunk = k24_chunk(vec![(10, vec![0xAA; 4]), (20, vec![0xBB; 4]), (30, vec![0xCC; 4])]);
        assert_eq!(chunk.min_bound(), Some(RunBound::TemplateK24(k24(10))));
        assert_eq!(chunk.max_bound(), Some(RunBound::TemplateK24(k24(30))));
    }

    #[test]
    fn contiguous_template_chunks_extend_overlapping_do_not() {
        let first = k24_chunk(vec![(10, vec![1; 4]), (20, vec![1; 4])]);
        let contiguous = k24_chunk(vec![(20, vec![1; 4]), (30, vec![1; 4])]);
        let overlapping = k24_chunk(vec![(15, vec![1; 4]), (25, vec![1; 4])]);

        let open_max = first.max_bound().unwrap();
        // second chunk's min (20) >= open run max (20) → extend (content tie).
        assert!(open_max.precedes_or_equal(&contiguous.min_bound().unwrap()));
        // overlapping chunk's min (15) < open run max (20) → new run.
        assert!(!open_max.precedes_or_equal(&overlapping.min_bound().unwrap()));
    }

    #[test]
    fn append_joins_batches_in_order_and_tracks_len() {
        let v = TemplateKeyVariant { cb: false, tertiary: false };
        let mut refs = TemplateArenaRefs::for_variant(v);
        let mk = |p: u64| TemplateRecordRef { key: k24(p), offset: p, len: 1, padding: 0 };
        refs.append(ExtractedRefs::K24(vec![mk(1), mk(2)].into_boxed_slice()));
        refs.append(ExtractedRefs::K24(vec![mk(3)].into_boxed_slice()));
        assert_eq!(refs.len(), 3);
        assert_eq!(refs.variant(), v);
        let TemplateArenaRefs::K24(inner) = &refs else { panic!("arm") };
        assert_eq!(inner.iter().map(|r| r.offset).collect::<Vec<_>>(), vec![1, 2, 3]);
    }

    #[test]
    #[should_panic(expected = "container arm K24 cannot take a K40 batch")]
    fn append_of_a_different_arm_panics_naming_both() {
        let mut refs =
            TemplateArenaRefs::for_variant(TemplateKeyVariant { cb: false, tertiary: false });
        refs.append(ExtractedRefs::K40(Box::new([])));
    }

    /// A later record whose `l_read_name` overruns its body is rejected with a
    /// placeholder name, not a slice panic in the violation branch.
    #[test]
    fn accumulator_reports_a_violation_whose_read_name_overruns_its_body() {
        use fgumi_raw_bam::testutil::make_bam_bytes;
        let mut cb = b"CBZ".to_vec();
        cb.extend_from_slice(b"AAAA\0");
        let first = make_bam_bytes(0, 100, 0x3, b"r0a", &[], 20, 0, 150, &cb);
        let mut bad = make_bam_bytes(0, 101, 0x3, b"r1a", &[], 20, 0, 151, &[]);
        bad[8] = 200;
        let mut acc = TemplateArenaAccumulator::from_header(
            &Header::default(),
            Some(SamTag::CB),
            KeyTypesSpec::None,
        );
        acc.push(&first, 4, u32::try_from(first.len()).unwrap()).unwrap();
        let err = acc.push(&bad, 200, u32::try_from(bad.len()).unwrap()).unwrap_err();
        assert!(
            format!("{err:#}").contains("record <malformed read name> carries a CB"),
            "{err:#}"
        );
    }

    /// A body shorter than the 32 fixed BAM bytes is an error, not an
    /// out-of-bounds panic in the key extractor — as the first record (which
    /// would provision the context) and as a later one.
    #[test]
    fn accumulator_rejects_a_body_shorter_than_a_bam_record() {
        use fgumi_raw_bam::testutil::make_bam_bytes;
        let short = vec![0u8; 12];
        let mut acc =
            TemplateArenaAccumulator::from_header(&Header::default(), None, KeyTypesSpec::Full);
        let err = acc.push(&short, 4, 12).unwrap_err();
        assert!(format!("{err:#}").contains("shorter than a BAM record"), "{err:#}");
        assert_eq!(acc.pending_len(), 0, "a rejected first record provisions nothing");

        let first = make_bam_bytes(0, 100, 0x3, b"r0a", &[], 20, 0, 150, &[]);
        acc.push(&first, 4, u32::try_from(first.len()).unwrap()).unwrap();
        let err = acc.push(&short, 200, 12).unwrap_err();
        assert!(format!("{err:#}").contains("shorter than a BAM record"), "{err:#}");
        assert_eq!(acc.pending_len(), 1);
    }

    /// Three ways to sort one run must agree record for record: the owned
    /// `TemplateChunkSorter` (which shares the key extractor with `extract_one`
    /// but none of the ref/seal plumbing), the per-record accumulator, and the
    /// batched path (provision → `extract_batch_erased` per 3-record batch →
    /// append → `seal_template_refs`). The fixture is 30 names × 8 positions ×
    /// 2 lane values × 2 sequence lengths, written in an order that is not
    /// sorted. Each (position, name) appears with both lane values, so where a
    /// case varies CB and/or MI those lanes decide the order; the two sequence
    /// lengths of one record share a full key, so tie order and batch boundaries
    /// (3 does not divide the tie groups) are both exercised. Cases cover all
    /// four lane widths, `auto` selecting each of K24, Cb32 and K40 from the
    /// first record, and the CB tag is passed to all three paths. Below the
    /// parallel radix threshold (`PARALLEL_SORT_THRESHOLD`, 256 Ki records), so
    /// this pins the serial radix.
    #[rstest::rstest]
    #[case::auto_plain(KeyTypesSpec::Auto, &[], &[], (false, false))]
    #[case::none_plain(KeyTypesSpec::None, &[], &[], (false, false))]
    #[case::full_plain(KeyTypesSpec::Full, &[], &[], (true, true))]
    #[case::cb32_varying_cb(
        KeyTypesSpec::Explicit { cb: true, tertiary: false },
        &["AAAA", "CCCC"],
        &[],
        (true, false)
    )]
    #[case::tert32_varying_mi(
        KeyTypesSpec::Explicit { cb: false, tertiary: true },
        &[],
        &["7", "3"],
        (false, true)
    )]
    #[case::full_varying_cb_and_mi(KeyTypesSpec::Full, &["GGGG", "AAAA"], &["2", "9"], (true, true))]
    #[case::auto_selects_cb32(KeyTypesSpec::Auto, &["CCCC", "AAAA"], &[], (true, false))]
    #[case::auto_selects_k40(KeyTypesSpec::Auto, &["AAAA", "TTTT"], &["5", "1"], (true, true))]
    fn accumulator_batched_path_and_owned_sorter_agree(
        #[case] key_types: KeyTypesSpec,
        #[case] cbs: &[&str],
        #[case] mis: &[&str],
        #[case] want: (bool, bool),
    ) {
        use crate::segmented_buf::SegmentedBuf;
        use crate::{RawExternalSorter, SortOrder};
        use fgumi_raw_bam::testutil::make_bam_bytes;

        let z_tag = |tag: SamTag, value: &str| {
            let mut aux = (*tag).to_vec();
            aux.push(b'Z');
            aux.extend_from_slice(value.as_bytes());
            aux.push(0);
            aux
        };
        let mut records: Vec<Vec<u8>> = Vec::new();
        for i in 0..120i32 {
            let pos = (i % 8) + 1;
            let name = format!("t{:03}", i % 30);
            // Alternate which lane value comes first, so input order is not the
            // lane order.
            for lane in [usize::try_from(i % 2).unwrap(), usize::try_from(1 - i % 2).unwrap()] {
                let mut aux = Vec::new();
                if !cbs.is_empty() {
                    aux.extend_from_slice(&z_tag(SamTag::CB, cbs[lane]));
                }
                if !mis.is_empty() {
                    aux.extend_from_slice(&z_tag(SamTag::MI, mis[lane]));
                }
                for seq_len in [20usize, 40usize] {
                    records.push(make_bam_bytes(
                        0,
                        pos,
                        0,
                        name.as_bytes(),
                        &[],
                        seq_len,
                        -1,
                        -1,
                        &aux,
                    ));
                }
            }
        }
        let cell_tag = Some(SamTag::CB);
        let want = TemplateKeyVariant { cb: want.0, tertiary: want.1 };

        let mut frames = Vec::new();
        let mut bodies = Vec::new();
        for r in &records {
            let at = frames.len();
            frames.extend_from_slice(&u32::try_from(r.len()).unwrap().to_le_bytes());
            frames.extend_from_slice(r);
            bodies.push((at as u64 + 4, u32::try_from(r.len()).unwrap()));
        }
        let mut buf = SegmentedBuf::with_capacity(0, frames.len() + 64);
        assert_eq!(buf.extend_from_slice(&frames), 0);
        let arena = Arc::new(PooledSegmentedBuf::unpooled(buf));

        // Owned oracle.
        let mut owned = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .memory_limit(256 * 1024 * 1024)
            .threads(2)
            .cell_tag(SamTag::CB)
            .key_types(key_types)
            .into_template_chunk_sorter(&Header::default())
            .expect("build owned template sorter");
        for r in &records {
            owned.push(r).expect("owned push");
        }
        let owned_chunk = owned.take_sorted_chunk_owned();
        let mut owned_keys: Vec<TemplateKey> = owned_chunk.iter().map(|(k, _)| *k).collect();
        owned_keys.dedup();
        assert!(owned_keys.len() < records.len(), "the fixture must contain full-key ties");
        // Where the case varies a lane, that lane must decide order: some
        // (position, name) group's first output record is not its first input
        // record (input order alternates which lane value comes first).
        let group = |r: &[u8]| (r[4..8].to_vec(), r[32..36].to_vec());
        let mut first_in: std::collections::HashMap<_, &[u8]> = std::collections::HashMap::new();
        for r in &records {
            first_in.entry(group(r)).or_insert(&r[..]);
        }
        let mut seen = std::collections::HashSet::new();
        let reordered = owned_chunk
            .iter()
            .filter(|(_, r)| seen.insert(group(r)) && first_in[&group(r)] != r.as_ref())
            .count();
        assert_eq!(reordered > 0, !cbs.is_empty() || !mis.is_empty(), "{reordered} groups");

        // Per-record accumulator.
        let mut acc =
            TemplateArenaAccumulator::from_header(&Header::default(), cell_tag, key_types);
        for (r, &(off, len)) in records.iter().zip(&bodies) {
            acc.push(r, off, len).unwrap();
        }
        let a = acc.seal(Arc::clone(&arena), 2);

        // Batched path.
        let ctx = TemplateKeySeed::from_header(&Header::default(), cell_tag, key_types)
            .provision(&records[0]);
        assert_eq!(ctx.variant, want, "the case selects the lane width it names");
        let mut refs = TemplateArenaRefs::for_variant(ctx.variant);
        let mut idx = 0u64;
        for chunk in bodies.chunks(3) {
            let start = chunk[0].0 - 4;
            let end = chunk.last().map(|&(o, l)| o + u64::from(l)).unwrap();
            let n = u32::try_from(chunk.len()).unwrap();
            let (batch, v) = ctx.extract_batch_erased(&arena, start, end, n, idx).unwrap();
            assert!(v.is_none());
            refs.append(batch);
            idx += u64::from(n);
        }
        let pool = BoundedSortPool::new("tmpl-sort");
        let b = seal_template_refs(&mut refs, Arc::clone(&arena), 2, pool.bounded(2));
        assert!(refs.is_empty(), "seal drains the container");
        assert_eq!(refs.variant(), ctx.variant, "and leaves the same arm");

        assert_eq!(owned_chunk.len(), records.len());
        assert_eq!((a.len(), b.len()), (records.len(), records.len()));
        assert_eq!(std::mem::discriminant(&a), std::mem::discriminant(&b));
        assert_eq!((a.min_bound(), a.max_bound()), (b.min_bound(), b.max_bound()));
        for (i, (_key, rec)) in owned_chunk.iter().enumerate() {
            assert_eq!(a.record_bytes(i), rec.as_ref(), "accumulator record {i} vs owned");
            assert_eq!(b.record_bytes(i), rec.as_ref(), "batched record {i} vs owned");
        }
    }
}
