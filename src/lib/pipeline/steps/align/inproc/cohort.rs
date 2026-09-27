//! Cohort cutter, SE/PE layout, and per-sub-batch id-base arithmetic.
//!
//! These are the parity-critical *pure* functions of the in-process aligner:
//! they reproduce, exactly, three upstream bwa-mem3 rules so that the merged
//! output is byte-identical to the subprocess `bwa-mem3` preset: the `-K` cohort
//! cut, the SE/PE layout, and the per-sub-batch read-id bases. Each is pinned
//! by a proptest against a literal Rust port of the upstream rule in the test
//! module below.
//!
//! Nothing here calls FFI or wires a pipeline step.
//!
//! These pure primitives are the foundation the in-process backend is built on:
//! the pipeline steps that construct and thread the work items through the pool
//! (`AlignPrepareStep`, `AlignSeedExtendStep`, `CohortPeStatStep`,
//! `AlignPairEmitStep`) consume the types and functions here.

use std::sync::Arc;

use crate::pipeline::core::item::HeapSize;
use crate::template::Template;

use super::engine::{AlignEngine, IdBases};
use super::gate::CohortLease;

/// Reproduces bwa-mem3's `-K` even-parity cohort cut
/// (`fast_reader_bseq.c:129-137`).
///
/// bwa-mem3's fast reader appends reads to the current cohort one at a time,
/// accumulating `size += l_seq` after each, and closes the cohort as soon as
/// `size >= chunk_size && (n & 1) == 0` — i.e. the accumulated bases have
/// reached `-K` **and** an even number of reads has been buffered. The parity
/// test matters: with `-p` smart pairing an even read count keeps read pairs
/// together at a cohort boundary in the common all-paired case, but a cut can
/// still fall *between* the two reads of a pair when the running read count is
/// odd at the pair's start (a mixed SE/PE input) — bwa then classifies the two
/// reads as two separate singles in adjacent cohorts.
///
/// The cut is therefore evaluated **per read**, not per template, which is why
/// the cutter's state is a running byte size and read count rather than a
/// per-template tally.
pub(crate) struct CohortCutter {
    /// The `-K` chunk size in bases; a cohort closes once its accumulated read
    /// bases reach this and the read count is even.
    chunk_size: u64,
    /// Accumulated `l_seq` of the reads buffered in the current (open) cohort.
    size: u64,
    /// Number of reads buffered in the current (open) cohort.
    n_reads_in_cohort: u64,
}

impl CohortCutter {
    /// Create a cutter for the given `-K` chunk size (bases). Callers pass the
    /// validated aligner chunk size (`>= 1`), matching the CLI.
    pub(crate) fn new(chunk_size: u64) -> Self {
        Self { chunk_size, size: 0, n_reads_in_cohort: 0 }
    }

    /// Feed one template's surviving reads (1 for single-end, 2 for a pair)
    /// through the per-read cut rule, in order. Returns `true` iff the cohort
    /// boundary falls **exactly after this template** — i.e. the template's
    /// last read triggered the even-parity cut, so the current cohort closes
    /// here and the next template opens a new one.
    ///
    /// A cut that fires *between* a pair's two reads (possible only when the
    /// running read count is odd at the pair's start, in a mixed SE/PE input)
    /// leaves the boundary inside the template: the second read correctly opens
    /// the next cohort in the cutter's state, but this method returns `false`
    /// because the boundary is not after the template. Splitting such a
    /// template into two singles across the two cohorts is the caller's concern
    /// (`AlignPrepareStep` splits it); the cutter's own running state stays exact
    /// regardless.
    ///
    /// Production drives the cutter one read at a time via [`Self::push_read`] (so
    /// it can observe a mid-pair cut); this whole-template convenience is used
    /// only by the parity proptest, hence `#[cfg(test)]`.
    #[cfg(test)]
    pub(crate) fn push_template(&mut self, l_seqs: impl Iterator<Item = u32>) -> bool {
        let mut cut_after_last = false;
        for l_seq in l_seqs {
            cut_after_last = self.push_read(l_seq);
        }
        cut_after_last
    }

    /// Feed **one** read's `l_seq` through the per-read cut rule and return
    /// `true` iff the cohort closes **after this read** — i.e. the accumulated
    /// bases have reached `-K` and the running read count is now even
    /// (`fast_reader_bseq.c:129-137`). On a cut the running accumulators reset,
    /// so the next `push_read` opens a fresh cohort.
    ///
    /// This per-read granularity is what lets the caller
    /// (`AlignPrepareStep`) observe a
    /// **mid-pair cut**: feeding a pair's two reads individually, a `true` after
    /// the *first* read means the cohort boundary falls between the pair's reads,
    /// and bwa then classifies them as two separate singles in adjacent cohorts.
    /// `Self::push_template` cannot express that (its single
    /// `bool` reports only a cut after the *last* read), so the prepare step
    /// drives the cutter one read at a time via this method.
    pub(crate) fn push_read(&mut self, l_seq: u32) -> bool {
        self.size += u64::from(l_seq);
        self.n_reads_in_cohort += 1;
        if self.size >= self.chunk_size && self.n_reads_in_cohort.is_multiple_of(2) {
            // Close the cohort after this read and open a fresh one.
            self.size = 0;
            self.n_reads_in_cohort = 0;
            true
        } else {
            false
        }
    }
}

/// Where one template of a sub-batch sits in bwa-mem3's `-p` SE/PE
/// classification. `idx` is the template's ordinal **within its group** — its
/// index into bwa's `pairs[]` or `singles[]` array — which is what the id
/// formulas and the record-emission origin indices key off.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum TemplateSlot {
    /// A paired template (two reads with equal names); the `idx`-th pair.
    Pair { idx: u32 },
    /// A single-end template (one read); the `idx`-th single.
    Single { idx: u32 },
}

/// The SE/PE classification of a sub-batch's templates, in input order.
///
/// [`Layout::AllPairs`] is the production shape — every template is a full pair,
/// so a template's index is exactly its pair index and no per-template vector is
/// needed. [`Layout::Mixed`] carries one [`TemplateSlot`] per template for the
/// (rare) SE / mixed inputs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum Layout {
    /// Every template contributes two reads: template `i` is pair `i`.
    AllPairs,
    /// One slot per template, in input order.
    Mixed(Vec<TemplateSlot>),
}

impl Layout {
    /// Classify a sub-batch's templates into bwa-mem3's `-p` SE/PE groups from
    /// each template's surviving-read count (2 = pair, otherwise single),
    /// reproducing `bseq_classify` (`bwa.cpp:258-274`).
    ///
    /// Upstream pairs two *consecutive equal-named reads*; because
    /// `GroupByQueryname` never emits two consecutive same-name templates, a
    /// template is a pair exactly when it contributes two reads and a single
    /// otherwise. A non-empty sub-batch whose templates are
    /// all pairs collapses to [`Layout::AllPairs`]; any single (or an empty
    /// sub-batch) yields [`Layout::Mixed`].
    pub(crate) fn classify(read_counts: &[u8]) -> Self {
        if !read_counts.is_empty() && read_counts.iter().all(|&c| c == 2) {
            return Self::AllPairs;
        }
        let mut slots = Vec::with_capacity(read_counts.len());
        let mut pair_idx: u32 = 0;
        let mut single_idx: u32 = 0;
        for &count in read_counts {
            if count == 2 {
                slots.push(TemplateSlot::Pair { idx: pair_idx });
                pair_idx += 1;
            } else {
                slots.push(TemplateSlot::Single { idx: single_idx });
                single_idx += 1;
            }
        }
        Self::Mixed(slots)
    }
}

/// Cohort-level facts a sub-batch carries so the barrier can compute its id
/// bases. Offsets are counts *within this cohort* accumulated
/// over the earlier sub-batches, in input order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct CohortPos {
    /// Global reads before this cohort — bwa's `n_processed` at cohort start.
    pub(crate) cohort_read_base: u64,
    /// Single-end reads in earlier sub-batches of this cohort.
    pub(crate) se_offset: u64,
    /// Paired templates in earlier sub-batches of this cohort.
    pub(crate) pe_offset: u64,
    /// Single-end reads in this sub-batch.
    pub(crate) n_se: u32,
    /// Paired templates in this sub-batch.
    pub(crate) n_pe: u32,
}

/// Identity of a sub-batch. `serial` is dense from 0 across the whole run and
/// becomes `BamTemplateBatch::batch_serial`; `(cohort, index_in_cohort)` is for
/// the cohort barrier.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct SubBatchId {
    pub(crate) serial: u64,
    pub(crate) cohort: u32,
    pub(crate) index_in_cohort: u32,
}

/// Attached to the **last** sub-batch of a cohort (only known at cut time): the
/// cohort's sub-batch count and total SE/PE tallies, which the barrier needs to
/// know a cohort is complete and to compute `pe_n_processed`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct CohortCloser {
    pub(crate) n_sub_batches: u32,
    pub(crate) cohort_n_se: u64,
    pub(crate) cohort_n_pe: u64,
}

/// Global-ordinal id bases for one sub-batch, reproducing bwa-mem3's `worker_sam`
/// ids over a `-p` cohort (`fastmap.cpp:908-958` + `bwamem.cpp:2795-2884`).
///
/// bwa processes a cohort's SE group with `n_processed = cohort_read_base` and
/// its PE group with `n_processed = cohort_read_base + cohort_n_se`
/// (`fastmap.cpp:924-944`); `worker_sam` then assigns single `i` the id
/// `n_processed + i` and pair `i` the id `(n_processed >> 1) + i`
/// (`bwamem.cpp:2795-2884`). For a sub-batch that is preceded within its cohort
/// by `se_offset` singles and `pe_offset` pairs, the *first* single/pair ids are
/// therefore:
///
/// - `first_single_id = cohort_read_base + se_offset`
/// - `first_pair_id = ((cohort_read_base + cohort_n_se) >> 1) + pe_offset`
///
/// and the engine adds the per-read index `i` on top. The `>> 1` is bwa's own
/// truncating shift on the pair `n_processed`, reproduced verbatim.
pub(crate) fn id_bases(pos: &CohortPos, cohort_n_se: u64) -> IdBases {
    IdBases {
        first_single_id: pos.cohort_read_base + pos.se_offset,
        first_pair_id: ((pos.cohort_read_base + cohort_n_se) >> 1) + pos.pe_offset,
    }
}

/// Heap footprint shared by the in-process work items: the summed
/// `Template::heap_size` of the moved-out unmapped templates. The C-side read
/// copies and alignment regions are not charged here: they live in the cohort's
/// resident state (owned by the [`CohortLease`], freed with its last clone),
/// not in any one sub-batch, and the lease's gate reservation is what bounds
/// how many cohorts hold them at once.
fn work_heap_size(unmapped: &[Template]) -> usize {
    unmapped.iter().map(Template::heap_size).sum::<usize>()
}

/// A sub-batch prepared for seeding: its identity, cohort position, SE/PE
/// layout, and the moved-out unmapped templates in input order, carrying the
/// cohort's [`CohortLease`]. Produced by `AlignPrepareStep`; consumed by
/// `AlignSeedExtendStep`.
pub(crate) struct AlignWork {
    pub(crate) id: SubBatchId,
    pub(crate) pos: CohortPos,
    /// `Some` on the cohort's last sub-batch.
    pub(crate) closer: Option<CohortCloser>,
    pub(crate) layout: Layout,
    pub(crate) unmapped: Vec<Template>,
    pub(crate) lease: CohortLease,
}

impl HeapSize for AlignWork {
    fn heap_size(&self) -> usize {
        work_heap_size(&self.unmapped)
    }
}

/// [`AlignWork`] after seeding + SE-extension: adds the sub-batch's share of its
/// cohort's resident state ([`AlignEngine::Ranges`]; the reads and alignment
/// regions themselves stay in the cohort). Produced by
/// `AlignSeedExtendStep` and consumed
/// by `CohortPeStatStep`.
///
/// Generic over the engine so the parallel seed/extend step and its tests share
/// one type: production drives it with `BwaMem3Engine`, the scheduling/ordering
/// tests with `FakeEngine` (`Ranges = Vec<TaggedRead>`). Only [`AlignWork`] stays
/// non-generic — it carries no ranges.
pub(crate) struct ExtendedWork<E: AlignEngine> {
    pub(crate) id: SubBatchId,
    pub(crate) pos: CohortPos,
    pub(crate) closer: Option<CohortCloser>,
    pub(crate) layout: Layout,
    pub(crate) unmapped: Vec<Template>,
    pub(crate) lease: CohortLease,
    pub(crate) ranges: E::Ranges,
}

impl<E: AlignEngine> HeapSize for ExtendedWork<E> {
    fn heap_size(&self) -> usize {
        work_heap_size(&self.unmapped)
    }
}

/// A sub-batch past the cohort barrier, ready to pair + emit: carries the
/// cohort `pestat` (`None` only when the sub-batch has no pairs) and the
/// computed id bases. Produced by
/// `CohortPeStatStep`; consumed by
/// `AlignPairEmitStep`.
///
/// Generic over the engine for the same reason as [`ExtendedWork`]: production
/// drives it with `BwaMem3Engine` (`PeStat = MemPeStat`), the barrier's
/// scheduling/ordering tests with `FakeEngine`.
pub(crate) struct PairWork<E: AlignEngine> {
    pub(crate) id: SubBatchId,
    pub(crate) layout: Layout,
    pub(crate) unmapped: Vec<Template>,
    pub(crate) lease: CohortLease,
    pub(crate) ranges: E::Ranges,
    /// The cohort insert-size model, shared read-only across the cohort's
    /// sub-batches; `None` when this sub-batch classifies no pairs.
    pub(crate) pestat: Option<Arc<E::PeStat>>,
    pub(crate) ids: IdBases,
}

impl<E: AlignEngine> HeapSize for PairWork<E> {
    fn heap_size(&self) -> usize {
        work_heap_size(&self.unmapped)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    // -------------------------------------------------------------------------
    // CohortCutter vs. bseq_read_fast
    // -------------------------------------------------------------------------

    /// Literal port of bwa-mem3's `bseq_read_fast` inner loop
    /// (`fast_reader_bseq.c:129-137`): walk the flat read-length stream,
    /// accumulate `size` per read, and record a cut (by read index) whenever
    /// `size >= chunk_size && (n & 1) == 0`, resetting the accumulators on each
    /// cut. The final partial cohort at EOF is *not* a cut here — bwa's reader
    /// only cuts on the size/parity condition; the tail is closed separately by
    /// the caller at end-of-input.
    fn port_bseq_read_fast(reads: &[u32], chunk_size: u64) -> Vec<usize> {
        let mut cuts = Vec::new();
        let mut size: u64 = 0;
        let mut n: u64 = 0;
        for (i, &l_seq) in reads.iter().enumerate() {
            size += u64::from(l_seq);
            n += 1;
            if size >= chunk_size && n & 1 == 0 {
                cuts.push(i);
                size = 0;
                n = 0;
            }
        }
        cuts
    }

    /// Drive [`CohortCutter`] one single-read template at a time over `reads`,
    /// collecting the read indices after which `push_template` reported a cut.
    /// With one read per template, the per-read rule and the per-template return
    /// coincide, so the result must equal [`port_bseq_read_fast`].
    fn run_cutter_over_reads(reads: &[u32], chunk_size: u64) -> Vec<usize> {
        let mut cutter = CohortCutter::new(chunk_size);
        let mut cuts = Vec::new();
        for (i, &l_seq) in reads.iter().enumerate() {
            if cutter.push_template(std::iter::once(l_seq)) {
                cuts.push(i);
            }
        }
        cuts
    }

    proptest! {
        #[test]
        fn cutter_matches_bseq_read_fast(
            reads in prop::collection::vec(1u32..300, 0..5000),
            chunk in 1u64..2000,
        ) {
            let want = port_bseq_read_fast(&reads, chunk);
            let got = run_cutter_over_reads(&reads, chunk);
            prop_assert_eq!(got, want);
        }
    }

    // -------------------------------------------------------------------------
    // Layout vs. bseq_classify
    // -------------------------------------------------------------------------

    /// Group ordinal a template lands in, mirroring [`TemplateSlot`] but without
    /// the production `AllPairs` collapse, so the port and the impl compare
    /// element-for-element.
    #[derive(Debug, PartialEq, Eq)]
    enum PortSlot {
        Pair(u32),
        Single(u32),
    }

    /// Literal port of bwa-mem3's `bseq_classify` (`bwa.cpp:258-274`): walk the
    /// flat read-name stream, pairing two consecutive reads iff their names are
    /// equal, else emitting a single, tracking each group's ordinal.
    fn port_bseq_classify(names: &[&[u8]]) -> Vec<PortSlot> {
        let mut out = Vec::new();
        let mut pair_idx: u32 = 0;
        let mut single_idx: u32 = 0;
        let mut i = 0;
        while i < names.len() {
            if i + 1 < names.len() && names[i] == names[i + 1] {
                out.push(PortSlot::Pair(pair_idx));
                pair_idx += 1;
                i += 2;
            } else {
                out.push(PortSlot::Single(single_idx));
                single_idx += 1;
                i += 1;
            }
        }
        out
    }

    /// Normalize a [`Layout`] to the port's per-template slot list so the two are
    /// directly comparable (`AllPairs` expands to consecutive `Pair` slots).
    fn layout_to_port_slots(layout: &Layout, n_templates: usize) -> Vec<PortSlot> {
        match layout {
            Layout::AllPairs => {
                (0..u32::try_from(n_templates).unwrap()).map(PortSlot::Pair).collect()
            }
            Layout::Mixed(slots) => slots
                .iter()
                .map(|slot| match *slot {
                    TemplateSlot::Pair { idx } => PortSlot::Pair(idx),
                    TemplateSlot::Single { idx } => PortSlot::Single(idx),
                })
                .collect(),
        }
    }

    proptest! {
        #[test]
        fn layout_matches_bseq_classify(
            // `true` = paired template (2 reads), `false` = single (1 read).
            template_is_pair in prop::collection::vec(any::<bool>(), 1..64),
        ) {
            // Each template gets a globally unique name (its index), upholding
            // the `GroupByQueryname` invariant that no two *consecutive*
            // templates share a name — so `bseq_classify` pairs exactly the
            // two-read templates.
            let names: Vec<Vec<u8>> = (0..template_is_pair.len())
                .map(|i| i.to_le_bytes().to_vec())
                .collect();

            // Build the flat read-name stream (each template's name repeated per
            // surviving read) and the per-template read counts.
            let mut flat_names: Vec<&[u8]> = Vec::new();
            let mut read_counts: Vec<u8> = Vec::new();
            for (i, &is_pair) in template_is_pair.iter().enumerate() {
                let count: u8 = if is_pair { 2 } else { 1 };
                read_counts.push(count);
                for _ in 0..count {
                    flat_names.push(names[i].as_slice());
                }
            }

            let want = port_bseq_classify(&flat_names);
            let layout = Layout::classify(&read_counts);
            let got = layout_to_port_slots(&layout, template_is_pair.len());
            prop_assert_eq!(got, want);

            // The `AllPairs` collapse must fire exactly when every template pairs.
            let all_pairs = template_is_pair.iter().all(|&p| p);
            prop_assert_eq!(matches!(layout, Layout::AllPairs), all_pairs);
        }
    }

    // -------------------------------------------------------------------------
    // id_bases vs. worker_sam / fastmap id assignment
    // -------------------------------------------------------------------------

    #[derive(Clone, Copy, Debug)]
    enum Kind {
        Single,
        Pair,
    }

    /// Count the [`Kind::Single`] items in an iterator as a `u64`.
    fn count_singles<'a>(kinds: impl Iterator<Item = &'a Kind>) -> u64 {
        u64::try_from(kinds.filter(|k| matches!(k, Kind::Single)).count())
            .expect("cohort single count fits u64")
    }

    /// First single/pair ids a sub-batch would be assigned, computed by a
    /// literal port of the upstream rule (`fastmap.cpp:924-944` +
    /// `bwamem.cpp:2795-2884`): the cohort's SE group runs at
    /// `n_processed = base`, its PE group at `n_processed = base + cohort_n_se`,
    /// and `worker_sam` numbers single `k` as `n_processed + k` and pair `k` as
    /// `(n_processed >> 1) + k`, where `k` is the group-global ordinal (= the
    /// count of same-kind items in earlier sub-batches).
    fn port_id_bases(base: u64, sub_batches: &[Vec<Kind>]) -> Vec<IdBases> {
        let cohort_n_se = count_singles(sub_batches.iter().flatten());
        let se_n_processed = base;
        let pe_n_processed = base + cohort_n_se;

        let mut se_seen: u64 = 0;
        let mut pe_seen: u64 = 0;
        let mut out = Vec::with_capacity(sub_batches.len());
        for sub in sub_batches {
            out.push(IdBases {
                first_single_id: se_n_processed + se_seen,
                first_pair_id: (pe_n_processed >> 1) + pe_seen,
            });
            for kind in sub {
                match kind {
                    Kind::Single => se_seen += 1,
                    Kind::Pair => pe_seen += 1,
                }
            }
        }
        out
    }

    /// Build a [`CohortPos`] per sub-batch (tracking the same SE/PE offsets the
    /// prepare step would) and evaluate [`id_bases`] for each.
    // `n_se`/`n_pe` are the domain terms (they are also the `CohortPos` field
    // names); the pair reads more clearly as a pair than renamed apart.
    #[allow(clippy::similar_names)]
    fn id_bases_over_cohort(base: u64, sub_batches: &[Vec<Kind>]) -> Vec<IdBases> {
        let cohort_n_se = count_singles(sub_batches.iter().flatten());

        let mut se_offset: u64 = 0;
        let mut pe_offset: u64 = 0;
        let mut out = Vec::with_capacity(sub_batches.len());
        for sub in sub_batches {
            let n_se = u32::try_from(sub.iter().filter(|k| matches!(k, Kind::Single)).count())
                .expect("sub-batch single count fits u32");
            let n_pe = u32::try_from(sub.iter().filter(|k| matches!(k, Kind::Pair)).count())
                .expect("sub-batch pair count fits u32");
            let pos = CohortPos { cohort_read_base: base, se_offset, pe_offset, n_se, n_pe };
            out.push(id_bases(&pos, cohort_n_se));
            se_offset += u64::from(n_se);
            pe_offset += u64::from(n_pe);
        }
        out
    }

    proptest! {
        #[test]
        fn id_bases_matches_worker_sam_port(
            // Even base: bwa's `n_processed` at a cohort start is a sum of
            // earlier (even-parity) cohort read counts. The `>> 1` is exercised
            // by odd `base + cohort_n_se` values regardless.
            base_half in 0u64..500_000,
            sub_batches in prop::collection::vec(
                prop::collection::vec(
                    prop_oneof![Just(Kind::Single), Just(Kind::Pair)],
                    0..12,
                ),
                1..25,
            ),
        ) {
            let base = base_half * 2;
            let want = port_id_bases(base, &sub_batches);
            let got = id_bases_over_cohort(base, &sub_batches);
            prop_assert_eq!(got, want);
        }
    }

    // -------------------------------------------------------------------------
    // Item HeapSize
    // -------------------------------------------------------------------------

    /// The item `HeapSize` is the moved templates' heap footprint alone: the
    /// C-side reads and regions belong to the cohort's resident state, not to any
    /// one sub-batch. Tested through the shared `work_heap_size` helper the
    /// `ExtendedWork`/`PairWork` impls delegate to.
    #[test]
    fn work_heap_size_sums_the_templates() {
        let templates = vec![Template::new(), Template::new()];
        let templates_bytes: usize = templates.iter().map(Template::heap_size).sum();

        assert_eq!(work_heap_size(&templates), templates_bytes);
        assert_eq!(work_heap_size(&[]), 0);
    }

    /// Compile + wiring check for the pure work item: `AlignWork` (which carries
    /// no regs) is constructible from the pure types and its `HeapSize` is the
    /// templates' footprint.
    #[test]
    fn align_work_heap_size_matches_templates() {
        let work = AlignWork {
            id: SubBatchId { serial: 0, cohort: 0, index_in_cohort: 0 },
            pos: CohortPos { cohort_read_base: 0, se_offset: 0, pe_offset: 0, n_se: 0, n_pe: 0 },
            closer: Some(CohortCloser { n_sub_batches: 1, cohort_n_se: 0, cohort_n_pe: 1 }),
            layout: Layout::AllPairs,
            unmapped: vec![Template::new()],
            lease: CohortLease::for_test(),
        };
        let expected: usize = work.unmapped.iter().map(Template::heap_size).sum();
        assert_eq!(work.heap_size(), expected);
    }
}
