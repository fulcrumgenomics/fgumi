//! Arena-borrowed queryname sort refs.
//!
//! A queryname run sorts 40-byte refs `(QnKey, offset, len)` over the run's
//! arena instead of owned name keys. A [`QnKey`] holds a 24-byte window of the
//! read name taken *after* an elided common prefix (ECP) — the
//! `instrument:run:flowcell:lane:` part every read of a lane shares — plus the
//! name's flag class, the tail length and the prefix's class id in a per-sort
//! [`QnClassTable`]. The prefixes are learned per batch of
//! [`QN_CLASS_BATCH_RECORDS`] records by [`QnRefBuilder`].
//!
//! [`cmp_qn_refs`] decides almost every comparison from the refs alone:
//! - two refs of one class compare their windows;
//! - refs of different classes compare `ECP ++ window` on the stack;
//! - only when a truncated window was read to its end does it read the two
//!   names from the arena.
//!
//! Then it compares flags, then arena offset. That is exactly the order of the
//! owned keys ([`QnComparator::Key`]) followed by the offset, so the sorted run
//! is byte-identical to sorting owned keys: equal keys keep input order.
//!
//! **The class rule.** An ECP holds no NUL and is empty or ends on a byte that
//! is not an ASCII digit (`trim_to_class_rule`). The natural order ends a name
//! at its first NUL, so a prefix holding one would decide nothing. The natural comparator enters its digit-run
//! state only when both current bytes are digits; after an identical prefix
//! that ends on a non-digit, both walkers stand at the tail start in the
//! plain-byte state, so walking the tails is walking the names from there.
//! Without the rule a digit run could be split between prefix and window
//! (`r1` + `0` vs `r1` + `a`).

use std::cmp::Ordering;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering as AtomicOrdering};
use std::sync::{Arc, OnceLock};

use crossbeam_utils::CachePadded;
use fgumi_raw_bam::sort::natural_cmp_bounded;

use crate::arena_pool::PooledSegmentedBuf;
use crate::inline::{ChunkRef, InMemoryChunk};
use crate::keys::{RawQuerynameKey, RawQuerynameLexKey, RawSortKey, qn_name_and_flags};
use crate::run_consolidate::SpillKeyKind;

/// Name bytes a [`QnKey`] carries after its class's elided common prefix.
pub const QN_WINDOW: usize = 24;
/// Most classes a [`QnClassTable`] holds; class ids are `u8`, and class 0 (the
/// empty prefix) is always valid.
pub const QN_MAX_CLASSES: usize = 256;
/// Records per classification batch: [`QnRefBuilder::reclassify`] runs once
/// per batch this size.
pub const QN_CLASS_BATCH_RECORDS: usize = 4096;
/// Longest elided common prefix a batch learns.
const QN_MAX_ECP: usize = 255;

/// Per-record queryname sort key stored in the ref array.
#[derive(Clone, Copy, Debug, Default)]
pub struct QnKey {
    /// `name[ecp_len..][..QN_WINDOW]`, zero-padded. Compared only against a key
    /// of the same class, or behind its class's prefix.
    pub window: [u8; QN_WINDOW],
    /// `queryname_flag_order(raw flags)`.
    pub flags: u16,
    /// The name's length after the class prefix (`l_read_name` is a `u8`, so
    /// at most 254). `> QN_WINDOW` means the window is truncated.
    pub tail_len: u8,
    /// The class id in the sort's [`QnClassTable`].
    pub class: u8,
}

/// The 40-byte queryname ref: `(key, body offset, body length)` — the record
/// tuple of an `InMemoryChunk<QnKey>`.
pub type QnRef = (QnKey, u64, u32);

const _: () = assert!(size_of::<QnKey>() == 28 && size_of::<QnRef>() == 40);

impl QnKey {
    /// A key for a name whose tail after the class prefix is `tail`.
    #[inline]
    fn new(tail: &[u8], flags: u16, class: u8) -> Self {
        let mut window = [0u8; QN_WINDOW];
        let n = tail.len().min(QN_WINDOW);
        window[..n].copy_from_slice(&tail[..n]);
        // A name is at most 254 bytes (`l_read_name` counts the NUL in a `u8`).
        let tail_len = u8::try_from(tail.len()).expect("a read name is at most 254 bytes");
        Self { window, flags, tail_len, class }
    }

    /// The known window bytes and whether the name continues past them.
    #[inline]
    fn known(&self) -> (&[u8], bool) {
        let t = usize::from(self.tail_len);
        (&self.window[..t.min(QN_WINDOW)], t > QN_WINDOW)
    }
}

/// The two queryname orders as zero-sized comparator types.
pub trait QnComparator: Copy + Default + Send + Sync + 'static {
    /// The canonical owned key: the definition of the order.
    type Key: RawSortKey + Default;
    /// The spill key kind of this order.
    const KIND: SpillKeyKind;
    /// Compare two full names (the slices `qn_name_and_flags` returns).
    fn cmp_names(a: &[u8], b: &[u8]) -> Ordering;
    /// Compare two known name prefixes; `a_trunc` says `a`'s name continues
    /// past `a`. `None` when the answer depends on bytes past a truncated
    /// prefix.
    fn cmp_window(a: &[u8], a_trunc: bool, b: &[u8], b_trunc: bool) -> Option<Ordering>;
    /// [`Self::cmp_window`] of two keys of one class. An order may answer
    /// common cases from the fixed, zero-padded windows directly.
    #[inline]
    #[must_use]
    fn cmp_same_class(a: &QnKey, b: &QnKey) -> Option<Ordering> {
        let ((wa, ta), (wb, tb)) = (a.known(), b.known());
        Self::cmp_window(wa, ta, wb, tb)
    }
}

/// A window as three big-endian words: comparing them compares the
/// zero-padded windows byte by byte.
#[inline]
fn window_words(w: &[u8; QN_WINDOW]) -> [u64; 3] {
    let word = |i: usize| u64::from_be_bytes(w[i..i + 8].try_into().expect("8 bytes"));
    [word(0), word(8), word(16)]
}

/// The first index where two windows differ, if they differ.
#[inline]
fn first_difference(a: &[u8; QN_WINDOW], b: &[u8; QN_WINDOW]) -> Option<usize> {
    let (wa, wb) = (window_words(a), window_words(b));
    (0..3).find(|&i| wa[i] != wb[i]).map(|i| {
        // Big-endian: the first differing byte holds the highest differing bit.
        8 * i + (wa[i] ^ wb[i]).leading_zeros() as usize / 8
    })
}

/// Lexicographic queryname order ([`RawQuerynameLexKey`]).
#[derive(Clone, Copy, Debug, Default)]
pub struct QnLex;

/// Natural queryname order ([`RawQuerynameKey`]).
#[derive(Clone, Copy, Debug, Default)]
pub struct QnNatural;

impl QnComparator for QnLex {
    type Key = RawQuerynameLexKey;
    const KIND: SpillKeyKind = SpillKeyKind::QuerynameLex;

    /// Compares the zero-padded windows as three words. Where the padded
    /// windows first differ, either both bytes are real (the slice compare
    /// decides there too) or one side is padding against a non-zero real byte —
    /// that side is a proper prefix, which the slice compare also orders first.
    /// Equal padded windows leave only the lengths, as in [`Self::cmp_window`].
    #[inline]
    fn cmp_same_class(a: &QnKey, b: &QnKey) -> Option<Ordering> {
        match window_words(&a.window).cmp(&window_words(&b.window)) {
            Ordering::Equal => {}
            decided => return Some(decided),
        }
        let ((wa, ta), (wb, tb)) = (a.known(), b.known());
        match wa.len().cmp(&wb.len()) {
            Ordering::Less => (!ta).then_some(Ordering::Less),
            Ordering::Greater => (!tb).then_some(Ordering::Greater),
            Ordering::Equal => match (ta, tb) {
                (false, false) => Some(Ordering::Equal),
                (true, false) => Some(Ordering::Greater),
                (false, true) => Some(Ordering::Less),
                (true, true) => None,
            },
        }
    }

    #[inline]
    fn cmp_names(a: &[u8], b: &[u8]) -> Ordering {
        a.cmp(b)
    }

    #[inline]
    fn cmp_window(a: &[u8], a_trunc: bool, b: &[u8], b_trunc: bool) -> Option<Ordering> {
        let n = a.len().min(b.len());
        match a[..n].cmp(&b[..n]) {
            Ordering::Equal => {}
            decided => return Some(decided),
        }
        // One known prefix is a prefix of the other: the shorter side's next
        // byte decides, unless that side's name is complete.
        match a.len().cmp(&b.len()) {
            Ordering::Less => (!a_trunc).then_some(Ordering::Less),
            Ordering::Greater => (!b_trunc).then_some(Ordering::Greater),
            Ordering::Equal => match (a_trunc, b_trunc) {
                (false, false) => Some(Ordering::Equal),
                (true, false) => Some(Ordering::Greater),
                (false, true) => Some(Ordering::Less),
                (true, true) => None,
            },
        }
    }
}

impl QnComparator for QnNatural {
    type Key = RawQuerynameKey;
    const KIND: SpillKeyKind = SpillKeyKind::QuerynameNatural;

    /// The windows are identical up to their first difference `p`, so the walk
    /// reaches the start `r` of the digit run holding `p` (or `p` itself) in
    /// lockstep, in the plain-byte state — unless a NUL ended both names
    /// before it. When the bytes at `p` are real and neither is a digit, the
    /// walk decides there by the byte values; otherwise it resumes at `r`.
    #[inline]
    fn cmp_same_class(a: &QnKey, b: &QnKey) -> Option<Ordering> {
        let ((wa, ta), (wb, tb)) = (a.known(), b.known());
        if let Some(p) = first_difference(&a.window, &b.window)
            && p < wa.len()
            && p < wb.len()
            && !wa[..p].contains(&0)
        {
            if !wa[p].is_ascii_digit() && !wb[p].is_ascii_digit() {
                return Some(wa[p].cmp(&wb[p]));
            }
            let mut r = p;
            while r > 0 && wa[r - 1].is_ascii_digit() {
                r -= 1;
            }
            return Self::cmp_window(&wa[r..], ta, &wb[r..], tb);
        }
        Self::cmp_window(wa, ta, wb, tb)
    }

    #[inline]
    fn cmp_names(a: &[u8], b: &[u8]) -> Ordering {
        natural_cmp_bounded(a, b).0
    }

    #[inline]
    fn cmp_window(a: &[u8], a_trunc: bool, b: &[u8], b_trunc: bool) -> Option<Ordering> {
        // The walk decides from the bytes it has read; if it never read a
        // truncated side's end, the full walk reads the same bytes and decides
        // the same way.
        let (ord, ia, ib) = natural_cmp_bounded(a, b);
        if (a_trunc && ia >= a.len()) || (b_trunc && ib >= b.len()) { None } else { Some(ord) }
    }
}

/// Append-only table of elided common prefixes, one per sort, shared by every
/// extractor and seal of that sort. Reads are lock-free; entry 0 is the empty
/// prefix. Class ids are a performance hint only: [`cmp_qn_refs`] is correct
/// for any assignment, so concurrent registration order does not matter.
pub struct QnClassTable {
    ecps: Box<[OnceLock<Box<[u8]>>]>,
    /// Entries `[0, registered)` are initialised.
    registered: AtomicUsize,
    /// Serialises appends; readers never take it.
    register: parking_lot::Mutex<()>,
}

impl QnClassTable {
    /// A table holding only class 0, the empty prefix.
    #[must_use]
    pub fn new() -> Arc<Self> {
        let ecps: Box<[OnceLock<Box<[u8]>>]> =
            (0..QN_MAX_CLASSES).map(|_| OnceLock::new()).collect();
        let _ = ecps[0].set(Box::default());
        Arc::new(Self {
            ecps,
            registered: AtomicUsize::new(1),
            register: parking_lot::Mutex::new(()),
        })
    }

    /// The prefix of `class`. Lock-free.
    #[inline]
    #[must_use]
    pub fn ecp(&self, class: u8) -> &[u8] {
        self.ecps[usize::from(class)].get().map_or(&[], |b| b)
    }

    /// Classes registered, class 0 included.
    #[must_use]
    pub fn len(&self) -> usize {
        self.registered.load(AtomicOrdering::Acquire)
    }

    /// Always `false`: class 0 is always registered.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        false
    }

    /// The class of `ecp`, registering it if new. `ecp` must satisfy the class
    /// rule (no NUL; empty, or not ending on an ASCII digit). Returns 0 for the empty
    /// prefix and when the table is full.
    ///
    /// # Panics
    ///
    /// Never in practice: class ids are below [`QN_MAX_CLASSES`], so they fit `u8`.
    pub fn class_for(&self, ecp: &[u8]) -> u8 {
        debug_assert!(
            ecp.last().is_none_or(|b| !b.is_ascii_digit()) && !ecp.contains(&0),
            "a class prefix must hold no NUL and not end on a digit: {ecp:?}"
        );
        if ecp.is_empty() {
            return 0;
        }
        let seen = self.len();
        if let Some(class) = self.find(ecp, 1, seen) {
            return class;
        }
        let _append = self.register.lock();
        let now = self.len();
        if let Some(class) = self.find(ecp, seen, now) {
            return class;
        }
        if now == QN_MAX_CLASSES {
            return 0;
        }
        let _ = self.ecps[now].set(ecp.into());
        self.registered.store(now + 1, AtomicOrdering::Release);
        u8::try_from(now).expect("class ids fit u8")
    }

    fn find(&self, ecp: &[u8], from: usize, to: usize) -> Option<u8> {
        (from..to)
            .find(|&i| self.ecps[i].get().is_some_and(|e| &e[..] == ecp))
            .map(|i| u8::try_from(i).expect("class ids fit u8"))
    }
}

/// The queryname key state one sort shares across its fronts, seals and
/// merge: the class table, the counters, and whether compares are counted
/// (`--sort-stats`).
#[derive(Clone)]
pub struct QnSortState {
    /// The sort's class table.
    pub table: Arc<QnClassTable>,
    /// The sort's counters.
    pub stats: Arc<QnCmpStats>,
    /// Whether the run sorts count their compares into `stats`.
    pub count_compares: bool,
}

impl QnSortState {
    /// A fresh table and zeroed counters.
    #[must_use]
    pub fn new(count_compares: bool) -> Self {
        Self { table: QnClassTable::new(), stats: Arc::new(QnCmpStats::new()), count_compares }
    }

    /// The counters to pass a run sort: `Some` only when compares are counted.
    #[must_use]
    pub fn compare_stats(&self) -> Option<&QnCmpStats> {
        self.count_compares.then_some(&*self.stats)
    }
}

/// `lcp` cut to at most 255 bytes, before its first NUL, and then to its last
/// non-digit byte, so it satisfies the class rule (empty when nothing remains).
///
/// The NUL cut matters to the natural order, which ends a name at its first
/// NUL: a prefix holding one would make every name of its class compare equal
/// on the name, whatever its window says.
#[must_use]
pub(crate) fn trim_to_class_rule(lcp: &[u8]) -> &[u8] {
    let mut n = lcp.len().min(QN_MAX_ECP);
    if let Some(nul) = lcp[..n].iter().position(|&b| b == 0) {
        n = nul;
    }
    while n > 0 && lcp[n - 1].is_ascii_digit() {
        n -= 1;
    }
    &lcp[..n]
}

/// Running longest common prefix of one batch's names, bounded at 255 bytes.
#[derive(Clone, Debug)]
pub struct QnBatchLcp {
    buf: [u8; QN_MAX_ECP],
    len: usize,
    seen: bool,
}

impl Default for QnBatchLcp {
    fn default() -> Self {
        Self { buf: [0; QN_MAX_ECP], len: 0, seen: false }
    }
}

impl QnBatchLcp {
    /// Shrink the prefix to the part `name` shares.
    #[inline]
    pub fn fold(&mut self, name: &[u8]) {
        if self.seen {
            // The common case — the name shares the whole prefix — is one memcmp.
            if !name.starts_with(&self.buf[..self.len]) {
                let n = self.len.min(name.len());
                self.len = self.buf[..n].iter().zip(name).take_while(|(a, b)| a == b).count();
            }
        } else {
            self.len = name.len().min(QN_MAX_ECP);
            self.buf[..self.len].copy_from_slice(&name[..self.len]);
            self.seen = true;
        }
    }

    /// The prefix every folded name shares.
    #[must_use]
    pub fn prefix(&self) -> &[u8] {
        &self.buf[..self.len]
    }
}

/// Builds refs, carrying the class guess across consecutive batches. One per
/// extractor (or per `SortBuffer` strategy).
#[derive(Clone, Debug, Default)]
pub struct QnRefBuilder {
    guess: u8,
}

impl QnRefBuilder {
    /// A builder whose first guess is `class` (tests).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn with_guess(class: u8) -> Self {
        Self { guess: class }
    }

    /// The class the next batch is windowed against (tests).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn guess(&self) -> u8 {
        self.guess
    }

    /// One ref from one record body: windowed behind the guessed class's prefix
    /// when the name starts with it, else in class 0. Folds the name into `lcp`.
    #[inline]
    #[must_use]
    pub fn build(
        &self,
        table: &QnClassTable,
        lcp: &mut QnBatchLcp,
        body: &[u8],
        offset: u64,
        len: u32,
    ) -> QnRef {
        let (name, flags) = qn_name_and_flags(body);
        lcp.fold(name);
        (key_in_class(table, name, flags, self.guess), offset, len)
    }

    /// After a batch: when more than half its refs fell to class 0 and its names
    /// share a non-empty class-rule prefix, register that prefix, make it the
    /// guess and return it. The caller then re-windows the batch
    /// ([`Self::rewindow`]).
    pub fn adopt_batch_class(
        &mut self,
        table: &QnClassTable,
        refs: &[QnRef],
        lcp: &QnBatchLcp,
    ) -> Option<u8> {
        let class0 = refs.iter().filter(|r| r.0.class == 0).count();
        if class0 * 2 <= refs.len() {
            return None;
        }
        let ecp = trim_to_class_rule(lcp.prefix());
        if ecp.is_empty() {
            return None;
        }
        let class = table.class_for(ecp);
        (class != 0).then(|| {
            self.guess = class;
            class
        })
    }

    /// Re-window `refs` against `class`, reading each name from `arena`.
    pub fn rewindow(
        table: &QnClassTable,
        arena: &PooledSegmentedBuf,
        refs: &mut [QnRef],
        class: u8,
    ) {
        for r in refs {
            let body = arena_body(arena, r);
            let (name, flags) = qn_name_and_flags(body);
            r.0 = key_in_class(table, name, flags, class);
        }
    }

    /// [`Self::adopt_batch_class`] then [`Self::rewindow`]; `true` if it
    /// re-windowed.
    pub fn reclassify(
        &mut self,
        table: &QnClassTable,
        arena: &PooledSegmentedBuf,
        refs: &mut [QnRef],
        lcp: &QnBatchLcp,
    ) -> bool {
        match self.adopt_batch_class(table, refs, lcp) {
            Some(class) => {
                Self::rewindow(table, arena, refs, class);
                true
            }
            None => false,
        }
    }
}

/// The key of `name` in `class` when the name starts with that class's
/// prefix, else in class 0.
#[inline]
fn key_in_class(table: &QnClassTable, name: &[u8], flags: u16, class: u8) -> QnKey {
    let ecp = table.ecp(class);
    if class != 0 && name.starts_with(ecp) {
        QnKey::new(&name[ecp.len()..], flags, class)
    } else {
        QnKey::new(name, flags, 0)
    }
}

#[inline]
#[allow(clippy::cast_possible_truncation)] // arena offsets fit usize on supported targets
fn arena_body<'a>(arena: &'a PooledSegmentedBuf, r: &QnRef) -> &'a [u8] {
    arena.slice(r.1 as usize, r.2 as usize)
}

/// Longest `ECP ++ window` a cross-class compare builds on the stack.
const QN_STACK: usize = QN_MAX_ECP + QN_WINDOW;

/// The comparator over refs: name (by window, by prefix + window, or from the
/// arena), then flags, then arena offset. Counts its paths into `stats` when
/// `STATS`.
#[inline]
pub fn cmp_qn_refs<C: QnComparator, const STATS: bool>(
    table: &QnClassTable,
    arena: &PooledSegmentedBuf,
    stats: &QnCmpStats,
    a: &QnRef,
    b: &QnRef,
) -> Ordering {
    let (ka, kb) = (&a.0, &b.0);
    // The hot path — two refs of one class decided by their windows — stays
    // small enough to inline into the sort; the rest is out of line.
    let names = if ka.class == kb.class {
        if STATS {
            stats.window();
        }
        match C::cmp_same_class(ka, kb) {
            Some(decided) => decided,
            None => cmp_names_in_arena::<C, STATS>(arena, stats, a, b),
        }
    } else {
        cmp_across_classes::<C, STATS>(table, arena, stats, a, b)
    };
    match names {
        Ordering::Equal => ka.flags.cmp(&kb.flags).then_with(|| a.1.cmp(&b.1)),
        decided => decided,
    }
}

/// Compare two refs of different classes: each side's class prefix plus its
/// window, concatenated on the stack.
#[inline(never)]
fn cmp_across_classes<C: QnComparator, const STATS: bool>(
    table: &QnClassTable,
    arena: &PooledSegmentedBuf,
    stats: &QnCmpStats,
    a: &QnRef,
    b: &QnRef,
) -> Ordering {
    let (ka, kb) = (&a.0, &b.0);
    if STATS {
        stats.cross_class();
    }
    let ((wa, trunc_a), (wb, trunc_b)) = (ka.known(), kb.known());
    let (ea, eb) = (table.ecp(ka.class), table.ecp(kb.class));
    let mut la = [0u8; QN_STACK];
    let mut lb = [0u8; QN_STACK];
    la[..ea.len()].copy_from_slice(ea);
    la[ea.len()..ea.len() + wa.len()].copy_from_slice(wa);
    lb[..eb.len()].copy_from_slice(eb);
    lb[eb.len()..eb.len() + wb.len()].copy_from_slice(wb);
    match C::cmp_window(&la[..ea.len() + wa.len()], trunc_a, &lb[..eb.len() + wb.len()], trunc_b) {
        Some(decided) => decided,
        None => cmp_names_in_arena::<C, STATS>(arena, stats, a, b),
    }
}

/// Compare the two full names, read from the arena: the windows could not
/// decide.
#[cold]
#[inline(never)]
fn cmp_names_in_arena<C: QnComparator, const STATS: bool>(
    arena: &PooledSegmentedBuf,
    stats: &QnCmpStats,
    a: &QnRef,
    b: &QnRef,
) -> Ordering {
    if STATS {
        stats.fallback();
    }
    C::cmp_names(
        qn_name_and_flags(arena_body(arena, a)).0,
        qn_name_and_flags(arena_body(arena, b)).0,
    )
}

/// Sort a run's refs by [`cmp_qn_refs`] into a chunk over `arena`. The refs
/// move into the chunk; nothing is copied. Runs on the calling rayon pool.
#[must_use]
pub fn queryname_chunk_from_qn_refs<C: QnComparator>(
    table: &QnClassTable,
    arena: Arc<PooledSegmentedBuf>,
    mut refs: Vec<QnRef>,
    stats: Option<&QnCmpStats>,
) -> InMemoryChunk<QnKey> {
    use rayon::prelude::*;
    match stats {
        Some(stats) => {
            refs.par_sort_unstable_by(|a, b| cmp_qn_refs::<C, true>(table, &arena, stats, a, b));
        }
        None => refs
            .par_sort_unstable_by(|a, b| cmp_qn_refs::<C, false>(table, &arena, &UNCOUNTED, a, b)),
    }
    InMemoryChunk::from_parts(arena, refs)
}

/// The final run of a spilled sort as `pieces` pre-sorted chunks: `refs` split
/// into contiguous input-order ranges, each sorted by [`cmp_qn_refs`] (a total
/// order — ties break by arena offset, so each range's sort is stable), the
/// ranges in parallel on the current rayon pool. Each range is copied into its
/// own exactly-sized chunk over `arena` and `refs` is freed before returning,
/// so the pieces hold the run's only ref bytes (splitting the Vec instead
/// would leave the first piece owning the whole run's allocation). An empty
/// `refs` yields no chunks.
#[must_use]
pub fn queryname_chunks_from_qn_refs<C: QnComparator>(
    table: &QnClassTable,
    arena: &Arc<PooledSegmentedBuf>,
    mut refs: Vec<QnRef>,
    stats: Option<&QnCmpStats>,
    pieces: usize,
) -> Vec<InMemoryChunk<QnKey>> {
    use rayon::prelude::*;
    let size = crate::ref_sort::piece_len(refs.len(), pieces);
    let piece = |range: &mut [QnRef]| {
        match stats {
            Some(stats) => {
                range.sort_unstable_by(|a, b| cmp_qn_refs::<C, true>(table, arena, stats, a, b));
            }
            None => range
                .sort_unstable_by(|a, b| cmp_qn_refs::<C, false>(table, arena, &UNCOUNTED, a, b)),
        }
        InMemoryChunk::from_parts(Arc::clone(arena), range.to_vec())
    };
    refs.par_chunks_mut(size).map(piece).collect()
}

impl ChunkRef<RawQuerynameLexKey> for QnKey {
    #[inline]
    fn merge_key(&self, body: &[u8]) -> RawQuerynameLexKey {
        RawQuerynameLexKey::extract_from_record(body)
    }
}

impl ChunkRef<RawQuerynameKey> for QnKey {
    #[inline]
    fn merge_key(&self, body: &[u8]) -> RawQuerynameKey {
        RawQuerynameKey::extract_from_record(body)
    }
}

/// Number of per-thread counter cells in a [`QnCmpStats`].
const QN_STAT_CELLS: usize = 64;

/// One thread's compare counters (each cell sits in a [`CachePadded`]).
struct QnStatCell {
    window: AtomicU64,
    cross_class: AtomicU64,
    fallback: AtomicU64,
}

impl QnStatCell {
    const fn new() -> Self {
        Self {
            window: AtomicU64::new(0),
            cross_class: AtomicU64::new(0),
            fallback: AtomicU64::new(0),
        }
    }
}

/// The stats a sort that does not count compares passes to the comparator.
static UNCOUNTED: QnCmpStats = QnCmpStats::new();

/// Per-sort queryname key counters. The batch counters are always kept; the
/// compare counters only when the sort counts compares (`--sort-stats`), via
/// the `STATS` comparator.
pub struct QnCmpStats {
    cells: [CachePadded<QnStatCell>; QN_STAT_CELLS],
    batches: AtomicU64,
    reclassified: AtomicU64,
    refs: AtomicU64,
    class0_refs: AtomicU64,
}

impl Default for QnCmpStats {
    fn default() -> Self {
        Self::new()
    }
}

/// A [`QnCmpStats`] reading.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct QnStatsSnapshot {
    /// Same-class compares (decided by the windows unless they fell back).
    pub window: u64,
    /// Cross-class compares (decided on the stack unless they fell back).
    pub cross_class: u64,
    /// Compares that read both names from the arena.
    pub arena_fallbacks: u64,
    /// Batches classified.
    pub batches: u64,
    /// Batches re-windowed into a newly learned class.
    pub reclassified: u64,
    /// Refs built.
    pub refs: u64,
    /// Refs left in class 0 after classification.
    pub class0_refs: u64,
}

impl QnCmpStats {
    /// Zeroed counters.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            cells: [const { CachePadded::new(QnStatCell::new()) }; QN_STAT_CELLS],
            batches: AtomicU64::new(0),
            reclassified: AtomicU64::new(0),
            refs: AtomicU64::new(0),
            class0_refs: AtomicU64::new(0),
        }
    }

    #[inline]
    fn cell(&self) -> &QnStatCell {
        &self.cells[rayon::current_thread_index().unwrap_or(0) % QN_STAT_CELLS]
    }

    #[inline]
    fn window(&self) {
        self.cell().window.fetch_add(1, AtomicOrdering::Relaxed);
    }

    #[inline]
    fn cross_class(&self) {
        self.cell().cross_class.fetch_add(1, AtomicOrdering::Relaxed);
    }

    #[inline]
    fn fallback(&self) {
        self.cell().fallback.fetch_add(1, AtomicOrdering::Relaxed);
    }

    /// Count one classified batch of `refs` refs, `class0` of them left in
    /// class 0, and whether it was re-windowed.
    pub fn note_batch(&self, refs: usize, class0: usize, reclassified: bool) {
        self.batches.fetch_add(1, AtomicOrdering::Relaxed);
        self.reclassified.fetch_add(u64::from(reclassified), AtomicOrdering::Relaxed);
        self.refs.fetch_add(refs as u64, AtomicOrdering::Relaxed);
        self.class0_refs.fetch_add(class0 as u64, AtomicOrdering::Relaxed);
    }

    /// The counters now.
    #[must_use]
    pub fn snapshot(&self) -> QnStatsSnapshot {
        let sum = |f: fn(&QnStatCell) -> &AtomicU64| {
            self.cells.iter().map(|c| f(c).load(AtomicOrdering::Relaxed)).sum()
        };
        QnStatsSnapshot {
            window: sum(|c| &c.window),
            cross_class: sum(|c| &c.cross_class),
            arena_fallbacks: sum(|c| &c.fallback),
            batches: self.batches.load(AtomicOrdering::Relaxed),
            reclassified: self.reclassified.load(AtomicOrdering::Relaxed),
            refs: self.refs.load(AtomicOrdering::Relaxed),
            class0_refs: self.class0_refs.load(AtomicOrdering::Relaxed),
        }
    }

    /// The `Queryname keys:` summary line. The compare shares appear only when
    /// compares were counted.
    #[must_use]
    #[allow(clippy::cast_precision_loss)] // display percentages
    pub fn summary_line(&self, table: &QnClassTable) -> String {
        use std::fmt::Write as _;
        let s = self.snapshot();
        let mut line = format!(
            "Queryname keys: {} classes; {} of {} batches re-windowed; {} of {} refs in class 0",
            table.len(),
            s.reclassified,
            s.batches,
            s.class0_refs,
            s.refs
        );
        // Decided-without-the-arena and fallbacks partition the compares; the
        // cross-class share cuts across both, so it is stated apart.
        let compares = s.window + s.cross_class;
        if compares > 0 {
            let pct = |n: u64| 100.0 * n as f64 / compares as f64;
            let _ = write!(
                line,
                "; {compares} compares ({:.2}% cross-class): decided without the arena \
                 {:.2}%, arena fallbacks {:.2}%",
                pct(s.cross_class),
                pct(compares - s.arena_fallbacks),
                pct(s.arena_fallbacks)
            );
        }
        line
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::segmented_buf::SegmentedBuf;
    use fgumi_raw_bam::testutil::make_bam_bytes;
    use proptest::prelude::*;
    use rstest::rstest;

    /// A raw record body with name `name` (any bytes; at most 254) and flags.
    fn body_named(name: &[u8], flags: u16) -> Vec<u8> {
        make_bam_bytes(0, 0, flags, name, &[], 4, -1, -1, &[])
    }

    /// One arena holding `bodies`, and each body's `(offset, len)`.
    fn arena_of(bodies: &[Vec<u8>]) -> (Arc<PooledSegmentedBuf>, Vec<(u64, u32)>) {
        let mut buf = SegmentedBuf::new();
        let spans = bodies
            .iter()
            .map(|b| (buf.extend_from_slice(b) as u64, u32::try_from(b.len()).unwrap()))
            .collect();
        (Arc::new(PooledSegmentedBuf::unpooled(buf)), spans)
    }

    /// Records named `names`, each prefixed with `ecp` when a seeded draw falls
    /// under `prefix_p`, alternating R1/R2 flags; capped at 254 name bytes.
    fn arena_with_named_records(
        names: &[Vec<u8>],
        ecp: &[u8],
        prefix_p: f64,
    ) -> (Arc<PooledSegmentedBuf>, Vec<(u64, u32)>) {
        let bodies: Vec<Vec<u8>> = names
            .iter()
            .enumerate()
            .map(|(i, n)| {
                #[allow(clippy::cast_precision_loss)]
                let draw = ((i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 11) as f64
                    / (1u64 << 53) as f64;
                let mut name = if draw < prefix_p { [ecp, n].concat() } else { n.clone() };
                name.truncate(254);
                body_named(&name, if i % 2 == 0 { 0x41 } else { 0x81 })
            })
            .collect();
        arena_of(&bodies)
    }

    /// Refs over `spans`, built in batches of 7 by a builder whose first guess
    /// is `guess` (0: class 0; 1: the trimmed ECP; 2: a wrong prefix), each
    /// batch reclassified as production does.
    fn build_refs_with_guess(
        table: &QnClassTable,
        arena: &PooledSegmentedBuf,
        spans: &[(u64, u32)],
        ecp: &[u8],
        guess: u8,
    ) -> Vec<QnRef> {
        let start = match guess {
            1 => table.class_for(trim_to_class_rule(ecp)),
            2 => table.class_for(b"zz:"),
            _ => 0,
        };
        let mut builder = QnRefBuilder::with_guess(start);
        let mut refs = Vec::new();
        for batch in spans.chunks(7) {
            let mut lcp = QnBatchLcp::default();
            let mut out: Vec<QnRef> = batch
                .iter()
                .map(|&(o, l)| {
                    builder.build(
                        table,
                        &mut lcp,
                        arena_body(arena, &(QnKey::default(), o, l)),
                        o,
                        l,
                    )
                })
                .collect();
            builder.reclassify(table, arena, &mut out, &lcp);
            refs.extend(out);
        }
        refs
    }

    /// The owned-key order: `(C::Key::extract_from_record(body), offset)`.
    fn owned_order<C: QnComparator>(arena: &PooledSegmentedBuf, a: &QnRef, b: &QnRef) -> Ordering
    where
        C::Key: Ord,
    {
        let ka = C::Key::extract_from_record(arena_body(arena, a));
        let kb = C::Key::extract_from_record(arena_body(arena, b));
        ka.cmp(&kb).then(a.1.cmp(&b.1))
    }

    fn cmp<C: QnComparator>(
        table: &QnClassTable,
        arena: &PooledSegmentedBuf,
        a: &QnRef,
        b: &QnRef,
    ) -> Ordering {
        cmp_qn_refs::<C, false>(table, arena, &UNCOUNTED, a, b)
    }

    /// Name bytes for the comparator properties: letters, digits, `:`, `.`,
    /// `_`, high bytes and, rarely, NUL.
    fn qn_name_strategy() -> impl Strategy<Value = Vec<u8>> {
        proptest::collection::vec(
            prop_oneof![
                6 => proptest::sample::select(b"0123456789".to_vec()),
                4 => proptest::sample::select(b"ab:._".to_vec()),
                1 => 0x80u8..=0xFF,
                1 => Just(0u8),
            ],
            0..60,
        )
    }

    #[test]
    fn qn_key_is_28_bytes_and_the_ref_40() {
        assert_eq!(std::mem::size_of::<QnKey>(), 28);
        assert_eq!(std::mem::size_of::<QnRef>(), 40);
    }

    proptest! {
        /// `cmp_qn_refs` equals the `(owned key, offset)` order for every pair,
        /// under every class guess, for both comparators.
        #[test]
        fn window_compare_agrees_with_the_full_compare(
            names in proptest::collection::vec(qn_name_strategy(), 2..48),
            ecp in qn_name_strategy(),
            prefix_p in 0.0f64..1.0,
            guess in 0u8..3,
        ) {
            let (arena, spans) = arena_with_named_records(&names, &ecp, prefix_p);
            let table = QnClassTable::new();
            let refs = build_refs_with_guess(&table, &arena, &spans, &ecp, guess);
            for a in &refs {
                for b in &refs {
                    prop_assert_eq!(cmp::<QnLex>(&table, &arena, a, b), owned_order::<QnLex>(&arena, a, b));
                    prop_assert_eq!(cmp::<QnNatural>(&table, &arena, a, b), owned_order::<QnNatural>(&arena, a, b));
                }
            }
        }
    }

    /// Named edge cases, each compared both ways in both orders, in every
    /// class assignment: both class 0, both behind a shared prefix, and one
    /// behind it (cross-class).
    #[rstest]
    #[case::digit_run_cut_one_side(b"A00132:53:HFHJ:1:1646:33332222222222222222".as_slice(), b"A00132:53:HFHJ:1:1646:3334".as_slice())]
    #[case::digit_run_cut_both_sides(b"P:1111111111111111111111111111".as_slice(), b"P:1111111111111111111111111112".as_slice())]
    #[case::leading_zeros_equal(b"x:01", b"x:1")]
    #[case::leading_zeros_inside(b"x:a001", b"x:a01")]
    #[case::leading_zeros_at_the_window_start(b"P:007", b"P:7")]
    #[case::equal_through_the_window(
        b"P:aaaaaaaaaaaaaaaaaaaaaaaaZ",
        b"P:aaaaaaaaaaaaaaaaaaaaaaaaY"
    )]
    #[case::longer_names_behind_equal_windows(
        b"P:aaaaaaaaaaaaaaaaaaaaaaaaaaaa",
        b"P:aaaaaaaaaaaaaaaaaaaaaaaaa"
    )]
    #[case::name_equals_the_prefix(b"P:", b"P:7")]
    #[case::empty_name(b"", b"a")]
    #[case::non_ascii(b"P:\xc3\xa91", b"P:\xc3\xa910")]
    #[case::raw_high_bytes(b"P:\xff9", b"P:\x809")]
    #[case::embedded_nul(b"P:ab\0c", b"P:ab\0b")]
    #[case::truncated_prefix_of_a_longer_cross_class_name(
        b"P:abcdefghijklmnopqrstuvwxyzABCD",
        b"P:abcdefghijklmnopqrstuvwx"
    )]
    fn natural_edge_cases_match_the_owned_key(#[case] a: &[u8], #[case] b: &[u8]) {
        for (x, y) in [(a, b), (b, a)] {
            let (arena, spans) = arena_of(&[body_named(x, 0x41), body_named(y, 0x41)]);
            let table = QnClassTable::new();
            let class = table.class_for(b"P:");
            let mut refs: Vec<QnRef> =
                spans.iter().map(|&(o, l)| (QnKey::default(), o, l)).collect();
            for (ca, cb) in [(0, 0), (class, class), (0, class), (class, 0)] {
                QnRefBuilder::rewindow(&table, &arena, &mut refs[..1], ca);
                QnRefBuilder::rewindow(&table, &arena, &mut refs[1..], cb);
                let (ra, rb) = (&refs[0], &refs[1]);
                assert_eq!(
                    cmp::<QnLex>(&table, &arena, ra, rb),
                    owned_order::<QnLex>(&arena, ra, rb),
                    "lex {x:?} vs {y:?} classes {ca}/{cb}"
                );
                assert_eq!(
                    cmp::<QnNatural>(&table, &arena, ra, rb),
                    owned_order::<QnNatural>(&arena, ra, rb),
                    "natural {x:?} vs {y:?} classes {ca}/{cb}"
                );
            }
        }
    }

    /// 254-byte names (the longest a BAM name can be) and a declared name length
    /// past the end of the body both compare exactly as the owned keys do.
    #[test]
    fn longest_and_overlong_names_match_the_owned_key() {
        let long_a = [b"P:".as_slice(), &[b'7'; 252]].concat();
        let long_b = [b"P:".as_slice(), &[b'7'; 251], b"8"].concat();
        let mut overlong = body_named(b"P:abc", 0x41);
        overlong[8] = 200; // l_read_name past the body end
        let bodies = vec![body_named(&long_a, 0x41), body_named(&long_b, 0x81), overlong];
        let (arena, spans) = arena_of(&bodies);
        let table = QnClassTable::new();
        let refs = build_refs_with_guess(&table, &arena, &spans, b"P:", 1);
        for a in &refs {
            for b in &refs {
                assert_eq!(cmp::<QnLex>(&table, &arena, a, b), owned_order::<QnLex>(&arena, a, b));
                assert_eq!(
                    cmp::<QnNatural>(&table, &arena, a, b),
                    owned_order::<QnNatural>(&arena, a, b)
                );
            }
        }
    }

    #[test]
    fn class_rule_prevents_the_digit_run_split() {
        // An ECP ending in a digit could cut a digit run between ECP and window; the rule trims it.
        assert_eq!(trim_to_class_rule(b"r1"), b"r");
        assert_eq!(trim_to_class_rule(b"A00132:53:HFHJ:1:"), b"A00132:53:HFHJ:1:");
        // "r10" vs "r1a": the full names give Greater; behind the trimmed prefix
        // `r` the windows "10" vs "1a" give the same.
        let (arena, spans) = arena_of(&[body_named(b"r10", 0x41), body_named(b"r1a", 0x41)]);
        let table = QnClassTable::new();
        let mut lcp = QnBatchLcp::default();
        lcp.fold(b"r10");
        lcp.fold(b"r1a");
        assert_eq!(lcp.prefix(), b"r1");
        let class = table.class_for(trim_to_class_rule(lcp.prefix()));
        let mut refs: Vec<QnRef> = spans.iter().map(|&(o, l)| (QnKey::default(), o, l)).collect();
        QnRefBuilder::rewindow(&table, &arena, &mut refs, class);
        assert_eq!(cmp::<QnNatural>(&table, &arena, &refs[0], &refs[1]), Ordering::Greater);
    }

    #[rstest]
    #[case::all_digits(b"0123".as_slice(), b"".as_slice())]
    #[case::trailing_colon_kept(b"A00132:53:", b"A00132:53:")]
    #[case::trailing_digits_popped(b"SRR36097899.12", b"SRR36097899.")]
    #[case::empty(b"", b"")]
    #[case::cut_before_a_nul(b"abc\0cd:", b"abc")]
    fn trim_to_class_rule_cases(#[case] lcp: &[u8], #[case] expected: &[u8]) {
        assert_eq!(trim_to_class_rule(lcp), expected);
    }

    #[test]
    fn trim_to_class_rule_caps_at_255_bytes() {
        let long = vec![b'a'; 300];
        assert_eq!(trim_to_class_rule(&long).len(), 255);
    }

    #[test]
    fn class_table_dedupes_exact_prefixes() {
        let t = QnClassTable::new();
        let a = t.class_for(b"L1:");
        assert_eq!(t.class_for(b"L1:"), a);
        assert_ne!(t.class_for(b"L2:"), a);
        assert_eq!(t.class_for(b""), 0);
        assert_eq!(t.len(), 3);
        assert_eq!(t.ecp(a), b"L1:");
    }

    #[test]
    fn class_table_overflow_falls_back_to_class_0() {
        let t = QnClassTable::new();
        for i in 0..300u32 {
            let _ = t.class_for(format!("P{i:04}:").as_bytes());
        }
        assert_eq!(t.len(), QN_MAX_CLASSES);
        assert_eq!(
            t.class_for(b"Q9999:"),
            0,
            "the 257th distinct prefix maps to class 0 (empty ECP), no panic"
        );
    }

    #[test]
    fn concurrent_class_for_registers_each_prefix_once() {
        let t = QnClassTable::new();
        let got: Vec<Vec<u8>> = std::thread::scope(|s| {
            let hs: Vec<_> = (0..8)
                .map(|_| {
                    let t = &t;
                    s.spawn(move || {
                        (0..64u32)
                            .map(|i| t.class_for(format!("L{i}:").as_bytes()))
                            .collect::<Vec<u8>>()
                    })
                })
                .collect();
            hs.into_iter().map(|h| h.join().unwrap()).collect()
        });
        assert!(
            got.windows(2).all(|w| w[0] == w[1]),
            "every thread sees the same class per prefix"
        );
        assert_eq!(t.len(), 64 + 1 /* class 0 */);
    }

    /// A batch whose guess fails is re-windowed into its learned prefix, which
    /// is registered and becomes the guess for the next batch.
    #[test]
    fn a_failed_guess_rewindows_the_batch_and_registers_its_prefix() {
        let names: Vec<Vec<u8>> =
            (0..10).map(|i| format!("A00132:53:HFHJ:1:{}:{}", 1100 + i, i).into_bytes()).collect();
        let (arena, spans) =
            arena_of(&names.iter().map(|n| body_named(n, 0x41)).collect::<Vec<_>>());
        let table = QnClassTable::new();
        let mut builder = QnRefBuilder::default();
        let mut lcp = QnBatchLcp::default();
        let mut refs: Vec<QnRef> = spans
            .iter()
            .map(|&(o, l)| {
                builder.build(&table, &mut lcp, arena_body(&arena, &(QnKey::default(), o, l)), o, l)
            })
            .collect();
        assert!(refs.iter().all(|r| r.0.class == 0));
        assert!(builder.reclassify(&table, &arena, &mut refs, &lcp));
        let class = builder.guess();
        assert_eq!(table.ecp(class), b"A00132:53:HFHJ:1:");
        assert!(refs.iter().all(|r| r.0.class == class && usize::from(r.0.tail_len) < QN_WINDOW));
        // The guess carries: the next batch is built in the class directly.
        let mut lcp2 = QnBatchLcp::default();
        let next =
            builder.build(&table, &mut lcp2, arena_body(&arena, &refs[0]), refs[0].1, refs[0].2);
        assert_eq!(next.0.class, class);
        assert!(
            !builder.reclassify(&table, &arena, &mut [next], &lcp2),
            "a matching guess is kept"
        );
    }

    /// A batch straddling two lanes keeps its minority in class 0: the majority
    /// matched the guess, so nothing is re-windowed.
    #[test]
    fn a_lane_straddling_batch_keeps_its_minority_in_class_0() {
        let table = QnClassTable::new();
        let lane1 = table.class_for(b"F:1:");
        let names: Vec<Vec<u8>> = (0..10)
            .map(|i| if i < 7 { format!("F:1:{i}") } else { format!("F:2:{i}") }.into_bytes())
            .collect();
        let (arena, spans) =
            arena_of(&names.iter().map(|n| body_named(n, 0x41)).collect::<Vec<_>>());
        let mut builder = QnRefBuilder::with_guess(lane1);
        let mut lcp = QnBatchLcp::default();
        let mut refs: Vec<QnRef> = spans
            .iter()
            .map(|&(o, l)| {
                builder.build(&table, &mut lcp, arena_body(&arena, &(QnKey::default(), o, l)), o, l)
            })
            .collect();
        assert!(!builder.reclassify(&table, &arena, &mut refs, &lcp));
        assert_eq!(refs.iter().filter(|r| r.0.class == 0).count(), 3);
        assert_eq!(builder.guess(), lane1);
    }

    /// `build` reads the name only within the body: a body cut short of its
    /// declared name yields the empty name, never a panic.
    #[test]
    fn build_reads_nothing_past_the_body() {
        let full = body_named(b"A00132:53:HFHJ:1:1101:1:1", 0x41);
        let short = &full[..34];
        let table = QnClassTable::new();
        let (r, ..) =
            QnRefBuilder::default().build(&table, &mut QnBatchLcp::default(), short, 0, 34);
        assert_eq!(r.tail_len, 0);
    }

    /// The chunk orders like a sort of owned keys broken by offset — two
    /// flowcells interleaved, mates and secondaries with exact-name ties, long
    /// tails and names past the inline key capacity.
    #[rstest]
    fn queryname_chunk_from_qn_refs_orders_like_the_owned_key_sort(
        #[values(false, true)] natural: bool,
    ) {
        let mut names = Vec::new();
        for i in 0..3000u32 {
            let fc = if i % 2 == 0 { "A00132:53:HFHJ" } else { "B00777:9:XXYY" };
            let mut name = format!("{fc}:{}:{}:{}:{}", 1 + i % 4, 1100 + i % 97, i % 31, i % 7);
            if i % 10 < 3 {
                name.push_str(":ACGTACGTACGTACGTACGT");
            }
            if i % 20 == 0 {
                name.push_str(&"N".repeat(30));
            }
            names.push(name.into_bytes());
        }
        let bodies: Vec<Vec<u8>> = names
            .iter()
            .enumerate()
            .flat_map(|(i, n)| {
                let mut v = vec![body_named(n, 0x41), body_named(n, 0x81)];
                if i % 5 == 0 {
                    v.push(body_named(n, 0x141));
                    v.push(body_named(n, 0x141));
                }
                v
            })
            .collect();
        let (arena, spans) = arena_of(&bodies);
        let table = QnClassTable::new();
        let mut refs = Vec::new();
        let mut builder = QnRefBuilder::default();
        for batch in spans.chunks(QN_CLASS_BATCH_RECORDS) {
            let mut lcp = QnBatchLcp::default();
            let mut out: Vec<QnRef> = batch
                .iter()
                .map(|&(o, l)| {
                    builder.build(
                        &table,
                        &mut lcp,
                        arena_body(&arena, &(QnKey::default(), o, l)),
                        o,
                        l,
                    )
                })
                .collect();
            builder.reclassify(&table, &arena, &mut out, &lcp);
            refs.extend(out);
        }
        let mut expected = refs.clone();
        if natural {
            expected.sort_by(|a, b| owned_order::<QnNatural>(&arena, a, b));
        } else {
            expected.sort_by(|a, b| owned_order::<QnLex>(&arena, a, b));
        }
        let stats = QnCmpStats::new();
        let chunk = if natural {
            queryname_chunk_from_qn_refs::<QnNatural>(
                &table,
                Arc::clone(&arena),
                refs,
                Some(&stats),
            )
        } else {
            queryname_chunk_from_qn_refs::<QnLex>(&table, Arc::clone(&arena), refs, Some(&stats))
        };
        let got: Vec<&[u8]> = (0..chunk.len()).map(|i| chunk.record_bytes(i)).collect();
        let want: Vec<&[u8]> = expected.iter().map(|r| arena_body(&arena, r)).collect();
        assert!(got == want, "the chunk is not in owned-key order");
        let s = stats.snapshot();
        assert!(s.window + s.cross_class > 0, "the sort counted its compares");
        assert!(s.arena_fallbacks > 0, "names with tails past the window fall back: {s:?}");
    }

    /// The final run's pieces are exactly sized and together hold the run's
    /// refs once: the original allocation (here twice the refs) is freed, not
    /// kept by the first piece. Each piece is its range in owned-key order.
    #[rstest]
    fn pieces_hold_the_refs_once(#[values(2, 3, 8)] pieces: usize) {
        let bodies: Vec<Vec<u8>> = (0..1_001u32)
            .map(|i| body_named(format!("r{}", (i * 7_919) % 1_001).as_bytes(), 0x41))
            .collect();
        let (arena, spans) = arena_of(&bodies);
        let table = QnClassTable::new();
        let mut refs: Vec<QnRef> = Vec::with_capacity(2 * spans.len());
        refs.extend(spans.iter().map(|&(o, l)| {
            QnRefBuilder::default().build(
                &table,
                &mut QnBatchLcp::default(),
                arena_body(&arena, &(QnKey::default(), o, l)),
                o,
                l,
            )
        }));
        let size = crate::ref_sort::piece_len(refs.len(), pieces);
        let mut expected = refs.clone();
        for range in expected.chunks_mut(size) {
            range.sort_by(|a, b| owned_order::<QnLex>(&arena, a, b));
        }
        let chunks = queryname_chunks_from_qn_refs::<QnLex>(&table, &arena, refs, None, pieces);
        assert_eq!(chunks.len(), 1_001usize.div_ceil(size));
        let mut at = 0;
        for chunk in &chunks {
            assert_eq!(chunk.records_capacity(), chunk.len(), "a piece holds no spare capacity");
            for i in 0..chunk.len() {
                assert_eq!(chunk.record_bytes(i), arena_body(&arena, &expected[at + i]));
            }
            at += chunk.len();
        }
        assert_eq!(at, 1_001);
    }

    /// Exact ties — one name, one flag class — come out in ascending arena
    /// offset (ingest) order whatever order the refs arrive in, so the
    /// unstable parallel sort still yields the stable order.
    #[rstest]
    fn exact_ties_break_by_ascending_offset(#[values(false, true)] natural: bool) {
        let n = 2_000usize;
        let bodies: Vec<Vec<u8>> = (0..n)
            .map(|i| {
                let mut b = body_named(b"dup/name", 0x141);
                b[4..8].copy_from_slice(&u32::try_from(i).unwrap().to_le_bytes()); // pos = index
                b
            })
            .collect();
        let (arena, spans) = arena_of(&bodies);
        let table = QnClassTable::new();
        // Scrambled: odds then evens.
        let refs: Vec<QnRef> = (0..n)
            .filter(|i| i % 2 == 1)
            .chain((0..n).filter(|i| i % 2 == 0))
            .map(|i| {
                let (o, l) = spans[i];
                QnRefBuilder::default().build(&table, &mut QnBatchLcp::default(), &bodies[i], o, l)
            })
            .collect();
        let chunk = if natural {
            queryname_chunk_from_qn_refs::<QnNatural>(&table, Arc::clone(&arena), refs, None)
        } else {
            queryname_chunk_from_qn_refs::<QnLex>(&table, Arc::clone(&arena), refs, None)
        };
        let order: Vec<u32> = (0..chunk.len())
            .map(|i| u32::from_le_bytes(chunk.record_bytes(i)[4..8].try_into().unwrap()))
            .collect();
        assert_eq!(order, (0..u32::try_from(n).unwrap()).collect::<Vec<_>>());
    }

    #[test]
    fn the_summary_line_reports_classes_batches_and_compares() {
        let table = QnClassTable::new();
        let _ = table.class_for(b"A:");
        let stats = QnCmpStats::new();
        stats.note_batch(4, 3, true);
        let line = stats.summary_line(&table);
        assert_eq!(
            line,
            "Queryname keys: 2 classes; 1 of 1 batches re-windowed; 3 of 4 refs in class 0"
        );
        stats.window();
        stats.fallback();
        assert!(stats.summary_line(&table).ends_with(
            "; 1 compares (0.00% cross-class): decided without the arena 0.00%, arena \
             fallbacks 100.00%"
        ));
    }
}
