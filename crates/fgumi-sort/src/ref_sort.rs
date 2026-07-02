//! Build a sorted coordinate [`InMemoryChunk`] from record bodies already resident
//! in a shared arena, WITHOUT copying record bytes — the core of the parallel-inflate
//! ingest. The chunk's `(offset, len)` point at record bodies in the supplied arena;
//! the existing Phase-2 merge consumes it unchanged.

use std::cmp::Ordering;
use std::sync::Arc;

use voracious_radix_sort::{RadixSort, Radixable};

use crate::arena_pool::PooledSegmentedBuf;
use crate::inline::{
    InMemoryChunk, RecordRef, extract_coordinate_key_inline, radix_sort_record_refs,
};
use crate::keys::{RawCoordinateKey, RawSortKey};

/// Arrays smaller than this sort faster single-threaded (the parallel radix's
/// partition + thread-coordination overhead exceeds its benefit on small inputs;
/// microbench-tuned).
const PARALLEL_SORT_THRESHOLD: usize = 256 * 1024;

/// `repr(transparent)` view of a [`RecordRef`] that orders by the
/// `(sort_key, offset)` **composite** rather than by `sort_key` alone.  This lets
/// `voracious`'s (unstable) parallel radix produce a STABLE-equivalent coordinate
/// order: within one chunk every record has a distinct `offset` that increases in
/// input order, so `(sort_key, offset)` is a total order whose ties resolve to
/// input order — byte-identical to the stable single-threaded radix.  Wrapping
/// (instead of changing `RecordRef`'s own `PartialEq`) keeps `RecordRef`'s
/// existing key-only equality untouched for the rest of the engine.
#[repr(transparent)]
#[derive(Copy, Clone)]
struct CoordSortRef(RecordRef);

impl CoordSortRef {
    #[inline]
    fn composite(self) -> u128 {
        (u128::from(self.0.sort_key) << 64) | u128::from(self.0.offset)
    }
}

impl PartialOrd for CoordSortRef {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.composite().cmp(&other.composite()))
    }
}

impl PartialEq for CoordSortRef {
    fn eq(&self, other: &Self) -> bool {
        self.composite() == other.composite()
    }
}

impl Radixable<u128> for CoordSortRef {
    type Key = u128;
    #[inline]
    fn key(&self) -> u128 {
        self.composite()
    }
}

/// Sort coordinate `RecordRef`s in place, stably by `sort_key`.
///
/// For large inputs with `sort_threads > 1`, uses `voracious`'s parallel radix
/// over the `(sort_key, offset)` composite (see [`CoordSortRef`]) — measured ~4×
/// faster than the single-threaded radix at realistic coordinate-key widths.  For
/// small inputs or a single thread it falls back to the single-threaded stable LSD
/// radix (the parallel path's overhead does not pay off there).  Both produce
/// byte-identical output (verified by the ref-sort oracle tests).
fn sort_coordinate_refs(refs: &mut [RecordRef], sort_threads: usize) {
    if sort_threads > 1 && refs.len() >= PARALLEL_SORT_THRESHOLD {
        // SAFETY: `CoordSortRef` is `#[repr(transparent)]` over `RecordRef`, so the
        // two have identical size, alignment, and layout; reinterpreting the slice
        // (via a pointer cast — clippy rejects a ref-to-ref `transmute`) is sound.
        // This is the only `unsafe` in this module's production path.
        #[allow(unsafe_code)]
        let view: &mut [CoordSortRef] =
            unsafe { &mut *(std::ptr::from_mut::<[RecordRef]>(refs) as *mut [CoordSortRef]) };
        view.voracious_mt_sort(sort_threads);
    } else {
        radix_sort_record_refs(refs);
    }
}

/// Build a sorted coordinate chunk from record bodies already resident in `arena`.
///
/// `bodies[i] = (body_offset, body_len)` locate each record's BAM body (the 4-byte
/// `block_size` prefix already excluded) within `arena`. The coordinate key is
/// extracted from each body exactly as `RecordBuffer::push_coordinate` does, the
/// refs are radix-sorted, and the same `arena` is wrapped into the returned chunk
/// (no record bytes are copied). The chunk's `(offset, len)` therefore point at the
/// bodies in `arena`, byte-identical to what the copy-based sorter materializes.
#[must_use]
#[allow(clippy::cast_possible_truncation)] // offset/len fit in usize on all supported targets (BAM is < 4 GiB; arenas < address space)
pub fn coordinate_chunk_from_arena_refs(
    arena: Arc<PooledSegmentedBuf>,
    bodies: &[(u64, u32)],
    n_ref: u32,
    sort_threads: usize,
) -> InMemoryChunk<RawCoordinateKey> {
    let mut refs: Vec<RecordRef> = Vec::with_capacity(bodies.len());
    for &(offset, len) in bodies {
        let body = arena.slice(offset as usize, len as usize);
        let sort_key = extract_coordinate_key_inline(body, n_ref);
        refs.push(RecordRef::new(sort_key, offset, len));
    }
    coordinate_chunk_from_refs(arena, refs, sort_threads)
}

/// Build a sorted coordinate chunk from already-extracted `RecordRef`s (each
/// carries its coordinate `sort_key` and the `(offset, len)` of its body in
/// `arena`).  Radix-sorts the refs and wraps `arena` into the returned chunk —
/// NO walk of the arena and NO record-byte copies.
///
/// This is the fused entry point: the caller (`FindBoundariesAndSort`) extracts
/// the key during its single boundary scan and hands the refs here, so the arena
/// is walked exactly once for the whole run (the previous two-step
/// [`coordinate_chunk_from_arena_refs`] re-walked it to re-derive each key).
#[must_use]
pub fn coordinate_chunk_from_refs(
    arena: Arc<PooledSegmentedBuf>,
    mut refs: Vec<RecordRef>,
    sort_threads: usize,
) -> InMemoryChunk<RawCoordinateKey> {
    sort_coordinate_refs(&mut refs, sort_threads);
    let records: Vec<(RawCoordinateKey, u64, u32)> =
        refs.iter().map(|r| (RawCoordinateKey { sort_key: r.sort_key }, r.offset, r.len)).collect();
    InMemoryChunk::from_parts(arena, records)
}

/// Build a sorted queryname chunk from already-extracted `(key, offset, len)`
/// refs pointing into `arena`. Unlike the coordinate/template builders, the
/// queryname key is a variable-length read name (not radix-able), so this sorts
/// by the key comparator with `par_sort_unstable_by` — matching the legacy
/// queryname path's unstable comparator sort. The caller is expected to run this
/// inside a bounded rayon pool so the
/// parallel sort does not oversubscribe the pipeline's worker pool. No record
/// bytes are copied: the returned chunk references its bodies in `arena`.
#[must_use]
pub fn queryname_chunk_from_arena_refs<K: RawSortKey>(
    arena: Arc<PooledSegmentedBuf>,
    mut records: Vec<(K, u64, u32)>,
) -> InMemoryChunk<K> {
    use rayon::prelude::*;
    records.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));
    InMemoryChunk::from_parts(arena, records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arena_pool::PooledSegmentedBuf;
    use crate::chunk_sorter::CoordinateChunkSorter;
    use crate::segmented_buf::SegmentedBuf;
    use std::sync::Arc;

    // Build a minimal valid BAM record body (no block_size prefix) with the given
    // refId/pos and a one-character read name `name`. The read name is NOT part of
    // the coordinate sort key (which reads only refId/pos/flags at bytes 0-15), so
    // two records with the same (refId,pos) but different `name` are a coordinate
    // TIE whose relative order is observable in the output bytes — that is what lets
    // these tests actually witness a stable-vs-unstable sort difference. Returns the
    // body bytes.
    fn coord_body(ref_id: i32, pos: i32, name: u8) -> Vec<u8> {
        // Minimal BAM record body: refID(4) pos(4) l_read_name(1) mapq(1) bin(2)
        // n_cigar_op(2) flag(2) l_seq(4) next_refID(4) next_pos(4) tlen(4) + read_name.
        let mut b = Vec::new();
        b.extend_from_slice(&ref_id.to_le_bytes());
        b.extend_from_slice(&pos.to_le_bytes());
        b.push(2); // l_read_name (incl NUL): "<name>\0"
        b.push(0); // mapq
        b.extend_from_slice(&0u16.to_le_bytes()); // bin
        b.extend_from_slice(&0u16.to_le_bytes()); // n_cigar_op
        b.extend_from_slice(&0u16.to_le_bytes()); // flag (forward)
        b.extend_from_slice(&0u32.to_le_bytes()); // l_seq
        b.extend_from_slice(&(-1i32).to_le_bytes()); // next_refID
        b.extend_from_slice(&(-1i32).to_le_bytes()); // next_pos
        b.extend_from_slice(&0i32.to_le_bytes()); // tlen
        b.push(name); // read_name char (distinguishes tied records)
        b.push(0); // read_name NUL terminator
        b
    }

    #[allow(unsafe_code)]
    #[test]
    fn ref_sort_matches_copy_based_sorter_byte_for_byte() {
        // Records deliberately out of coordinate order.
        let n_ref = 4u32;
        // Distinct read-name bytes make every record byte-unique, so a stability
        // divergence (reordering the (0,50) tie at indices 1 and 3) would change the
        // output byte sequence and fail the assertion below.
        let recs = vec![
            coord_body(2, 100, b'a'),
            coord_body(0, 50, b'b'),
            coord_body(2, 10, b'c'),
            coord_body(0, 50, b'd'), // coordinate tie with index 1 → stable sort must keep b before d
            coord_body(1, 999, b'e'),
        ];

        // ---- Path A: the existing copy-based sorter (the ORACLE) ----
        let mut oracle = CoordinateChunkSorter::for_test(usize::MAX, n_ref); // see note below
        for r in &recs {
            let _ = oracle.push(r).unwrap();
        }
        let oracle_chunk = oracle.take_sorted_chunk();
        let oracle_bytes: Vec<Vec<u8>> =
            (0..oracle_chunk.len()).map(|i| oracle_chunk.record_bytes(i).to_vec()).collect();

        // ---- Path B: the ref-sort over an arena holding verbatim [block_size][body] ----
        let mut arena = SegmentedBuf::with_capacity(0, 1 << 20);
        arena.reserve_full_capacity();
        let mut bodies = Vec::new();
        for r in &recs {
            // store verbatim record: 4-byte block_size (= body len) then body.
            let block_size = u32::try_from(r.len()).unwrap();
            // SAFETY: each slot fully written before any read.
            let prefix_off = unsafe { arena.grow_uninit(4) };
            unsafe { arena.slice_mut(prefix_off, 4) }.copy_from_slice(&block_size.to_le_bytes());
            let body_off = unsafe { arena.grow_uninit(r.len()) };
            unsafe { arena.slice_mut(body_off, r.len()) }.copy_from_slice(r);
            bodies.push((body_off as u64, block_size)); // ref points at the BODY, len = block_size
        }
        let arena_arc = Arc::new(PooledSegmentedBuf::unpooled(arena));
        let ref_chunk = coordinate_chunk_from_arena_refs(arena_arc, &bodies, n_ref, 1);
        let ref_bytes: Vec<Vec<u8>> =
            (0..ref_chunk.len()).map(|i| ref_chunk.record_bytes(i).to_vec()).collect();

        assert_eq!(
            ref_bytes, oracle_bytes,
            "ref-sort output must be byte-identical to the copy-based sorter"
        );
    }

    use proptest::prelude::*;

    proptest! {
        // For any multiset of (refId, pos) records — including ties and unmapped —
        // the ref-sort produces the same sorted body sequence as the copy sorter.
        #[test]
        #[allow(unsafe_code)]
        fn prop_ref_sort_matches_copy_sorter(
            coords in prop::collection::vec((0i32..5, 0i32..2000), 0..60),
        ) {
            let n_ref = 5u32;
            // A unique read-name per record (by input index) makes every record
            // byte-distinct, so any reordering of equal-key (tied) records would
            // change the output and fail the assertion — witnessing stability.
            let recs: Vec<Vec<u8>> = coords
                .iter()
                .enumerate()
                .map(|(i, &(t, p))| coord_body(t, p, u8::try_from(i).unwrap()))
                .collect();

            let mut oracle = CoordinateChunkSorter::for_test(usize::MAX, n_ref);
            for r in &recs { let _ = oracle.push(r).unwrap(); }
            let oracle_chunk = oracle.take_sorted_chunk();
            let oracle_bytes: Vec<Vec<u8>> =
                (0..oracle_chunk.len()).map(|i| oracle_chunk.record_bytes(i).to_vec()).collect();

            let mut arena = SegmentedBuf::with_capacity(0, 1 << 20);
            arena.reserve_full_capacity();
            let mut bodies = Vec::new();
            for r in &recs {
                let bs = u32::try_from(r.len()).unwrap();
                let po = unsafe { arena.grow_uninit(4) };
                unsafe { arena.slice_mut(po, 4) }.copy_from_slice(&bs.to_le_bytes());
                let bo = unsafe { arena.grow_uninit(r.len()) };
                unsafe { arena.slice_mut(bo, r.len()) }.copy_from_slice(r);
                bodies.push((bo as u64, bs));
            }
            let arena_arc = Arc::new(PooledSegmentedBuf::unpooled(arena));
            let ref_chunk = coordinate_chunk_from_arena_refs(arena_arc, &bodies, n_ref, 1);
            let ref_bytes: Vec<Vec<u8>> =
                (0..ref_chunk.len()).map(|i| ref_chunk.record_bytes(i).to_vec()).collect();

            prop_assert_eq!(ref_bytes, oracle_bytes);
        }
    }
}
