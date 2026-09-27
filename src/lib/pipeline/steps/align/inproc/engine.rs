//! The alignment engine adapter: the *single* module that names [`bwa_mem3_rs`].
//!
//! The scheduling steps of the in-process backend
//! ([`AlignSeedExtendStep`](super::seed_extend::AlignSeedExtendStep),
//! [`CohortPeStatStep`](super::pestat::CohortPeStatStep),
//! [`AlignPairEmitStep`](super::pair_emit::AlignPairEmitStep)) drive bwa-mem3's
//! three-phase `seed_extend -> infer_cohort -> pair_emit` pipeline. They do so
//! through the [`AlignEngine`] trait rather than calling `bwa_mem3_rs` directly,
//! for two reasons:
//!
//! 1. **Testability without FFI.** `FakeEngine` (a `#[cfg(test)]` sibling) is a
//!    deterministic, allocation-only stand-in that lets the scheduling/ordering
//!    tests exercise the real pipeline steps at many thread counts
//!    and sub-batch sizes without loading a bwa-mem3 index or compiling any C++.
//! 2. **One seam.** Every `bwa_mem3_rs` name the rest of the `inproc` module
//!    tree needs is re-exported from here: the types the steps pass
//!    ([`RecordSink`], [`RecordOrigin`], [`IdBases`]), the index and options the
//!    backend loads ([`BwaIndex`], [`MemOpts`]), and the linked [`version`] and
//!    [`shm`] staging probe. The rest of the tree is written against this module,
//!    not the binding crate.
//!
//! The real implementation, [`BwaMem3Engine`], is a thin wrapper over
//! `bwa_mem3_rs`'s [`ResidentCohort`]: one per `-K` cohort keeps every sub-batch's
//! reads and alignment regions resident from seed/extend to pair/emit, so no
//! sub-batch copies or frees them. The engine holds the shared, read-only
//! `Arc<BwaIndex>` and `Arc<MemOpts>` and maps each trait method onto the
//! corresponding `ResidentCohort` call. That API is entirely safe, so this module
//! adds no `unsafe`.
//!
//! The engine trait and both impls are consumed by the in-process backend's
//! steps ([`super::seed_extend`], [`super::pestat`], [`super::pair_emit`]) and
//! this module's tests.

use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use bwa_mem3_rs::{
    MemPeStat, PeOrientation, ReadPair as BwaReadPair, ResidentCohort, ResidentRange,
    SingleRead as BwaSingleRead,
};

use super::gate::CohortLease;

// Re-exports for the rest of the `inproc` module tree: these are the only
// `bwa_mem3_rs` names it references, and they go through here so `engine.rs`
// stays the sole module that names the binding crate.
pub(crate) use bwa_mem3_rs::{BwaIndex, IdBases, MemOpts, RecordOrigin, RecordSink, shm, version};

// ---------------------------------------------------------------------------
// Borrowed input batch (fgumi side)
// ---------------------------------------------------------------------------

/// One read as borrowed name/seq/qual spans, in the exact bytes bwa-mem3 must
/// see: `name` is the bare queryname (no `/1`, `/2`, no UMI),
/// `seq` is ASCII `ACGTN`, and `qual` is Phred+33 (or `None` for `'*'`).
///
/// The spans borrow storage owned elsewhere — in production a per-worker
/// `ReadArena` — so building an [`EngineBatch`] and converting it to a
/// [`bwa_mem3_rs::ReadBatch`] copies no sequence bytes.
#[derive(Clone, Copy, Debug)]
pub(crate) struct EngineRead<'a> {
    pub(crate) name: &'a [u8],
    pub(crate) seq: &'a [u8],
    pub(crate) qual: Option<&'a [u8]>,
}

impl<'a> EngineRead<'a> {
    /// Construct a read from its borrowed spans.
    pub(crate) fn new(name: &'a [u8], seq: &'a [u8], qual: Option<&'a [u8]>) -> Self {
        Self { name, seq, qual }
    }
}

/// A paired-end template's two reads (R1, R2), in input order.
#[derive(Clone, Copy, Debug)]
pub(crate) struct EnginePair<'a> {
    pub(crate) r1: EngineRead<'a>,
    pub(crate) r2: EngineRead<'a>,
}

/// A sub-batch of reads to seed+extend, already split into bwa-mem3's `-p` SE/PE
/// groups (the split is computed by [`super::cohort::Layout`] in the prepare
/// step). This is the fgumi-side, backend-neutral shape; [`BwaMem3Engine`]
/// converts it to a [`bwa_mem3_rs::ReadBatch`] borrowing the same spans, and
/// `FakeEngine` reads only the group lengths.
#[derive(Clone, Copy, Debug)]
pub(crate) struct EngineBatch<'a> {
    pub(crate) pairs: &'a [EnginePair<'a>],
    pub(crate) singles: &'a [EngineRead<'a>],
}

// ---------------------------------------------------------------------------
// The engine trait
// ---------------------------------------------------------------------------

/// The three-phase bwa-mem3 alignment interface the scheduling steps bind to,
/// over a **resident cohort**: one `-K` cohort's decoded reads and alignment regions stay
/// resident in a [`Self::Cohort`] from seed/extend through pair/emit, and each
/// sub-batch owns only its [`Self::Ranges`] of it. Nothing is copied or freed
/// per sub-batch; the cohort is freed once, when its last work item drops.
///
/// The cohort lives in the sub-batches' shared
/// [`CohortLease`] (see
/// [`cohort_of`]), so it is shared across pool workers and must be
/// `Send + Sync`. `Ranges` and `Scratch` are owned by exactly one work item /
/// worker at a time and only need `Send`. `PeStat` is inferred once per cohort
/// on a serial step and shared read-only across the cohort's `pair_emit` items.
pub(crate) trait AlignEngine: Send + Sync + 'static {
    /// Per-worker scratch, reused across calls. `Send`, not
    /// `Sync`: exactly one lives on each pool worker. `'static` because it is
    /// owned by the `'static` [`AlignSeedExtendStep`](super::seed_extend) /
    /// pair-emit steps.
    type Scratch: Send + 'static;
    /// One cohort's resident alignment state, shared by all its sub-batches.
    type Cohort: Send + Sync + 'static;
    /// One sub-batch's reserved share of its cohort: produced by
    /// [`Self::seed_extend`], consumed by [`Self::pair_emit`]. `'static` because
    /// it rides inside the `'static`
    /// [`ExtendedWork`](super::cohort::ExtendedWork)/`PairWork` pipeline items.
    type Ranges: Send + 'static;
    /// The cohort insert-size model. Shared read-only across a cohort's
    /// `pair_emit` items, hence `Send + Sync`; `'static` because it is `Arc`'d
    /// into the `'static` `PairWork` items.
    type PeStat: Send + Sync + 'static;

    /// Allocate a fresh per-worker scratch.
    ///
    /// # Errors
    /// Propagates the engine's allocation failure.
    fn new_scratch(&self) -> anyhow::Result<Self::Scratch>;

    /// Allocate an empty resident cohort. Called once per cohort, through
    /// [`cohort_of`].
    ///
    /// # Errors
    /// Propagates the engine's allocation failure.
    fn new_cohort(&self) -> anyhow::Result<Self::Cohort>;

    /// Phase 1: reserve `batch`'s share of `cohort`, copy its reads in, and
    /// seed + SE-extend them — the fused `worker_bwt_aln` work, per-read
    /// independent. Safe to call concurrently for different
    /// sub-batches of one cohort.
    ///
    /// # Errors
    /// Propagates an invalid batch or an engine failure.
    fn seed_extend(
        &self,
        scratch: &mut Self::Scratch,
        cohort: &Self::Cohort,
        batch: EngineBatch<'_>,
    ) -> anyhow::Result<Self::Ranges>;

    /// Phase 2 (cohort barrier): the `mem_pestat` insert-size model over every
    /// paired read of `cohort`. Every sub-batch's
    /// [`Self::seed_extend`] must have completed. The model is a histogram, so
    /// it does not depend on the order the sub-batches were extended in.
    ///
    /// # Errors
    /// Propagates an engine failure, including a sub-batch not yet extended.
    fn infer_cohort(&self, cohort: &Self::Cohort) -> anyhow::Result<Self::PeStat>;

    /// Phase 3: pairing + mate rescue + primary marking + emission for one
    /// sub-batch — the `worker_sam` work. Consumes `ranges`; emits one packed
    /// BAM record body per output record through `sink`, tagged with its
    /// sub-batch-local [`RecordOrigin`], in input order: pairs (R1 side then R2
    /// side, primary then supplementary), then singles.
    ///
    /// `pestat` may be `None` only for a sub-batch with no pairs; passing `None`
    /// for a sub-batch that classifies pairs is an error.
    ///
    /// # Errors
    /// Propagates a missing-pestat misuse or an engine failure.
    fn pair_emit(
        &self,
        scratch: &mut Self::Scratch,
        cohort: &Self::Cohort,
        ranges: Self::Ranges,
        pestat: Option<&Self::PeStat>,
        ids: IdBases,
        sink: &mut dyn RecordSink,
    ) -> anyhow::Result<()>;

    /// The number of PE templates (pairs) `ranges` holds — used by the barrier
    /// to decide whether a sub-batch needs the cohort pestat.
    fn n_pairs(ranges: &Self::Ranges) -> usize;

    /// A human-readable summary of `pestat`, logged by the cohort barrier at
    /// `info` once per cohort it computes a pestat for — the same information
    /// bwa's `[M::mem_pestat]` stderr lines carry on the subprocess route. The
    /// real engine reports the FR
    /// orientation's insert-size mean/std/low/high, matching upstream's log
    /// line; `FakeEngine` has no real insert-size model to report.
    fn pestat_summary(pestat: &Self::PeStat) -> String;
}

/// The resident cohort of the cohort `lease` belongs to, created through
/// `engine` by whichever of its sub-batches asks first and freed with the
/// lease's last clone.
///
/// # Errors
/// Returns an error if the cohort could not be created (the failure is kept,
/// so every sub-batch of that cohort reports it) or if the lease holds another
/// backend's state.
pub(crate) fn cohort_of<'l, E: AlignEngine>(
    engine: &E,
    lease: &'l CohortLease,
) -> io::Result<&'l E::Cohort> {
    let slot: &Result<E::Cohort, String> =
        lease.resident(|| engine.new_cohort().map_err(|e| format!("{e:#}")))?;
    slot.as_ref().map_err(|e| {
        io::Error::other(format!(
            "align-and-merge (in-process): creating the cohort's resident aligner state: {e}"
        ))
    })
}

// ---------------------------------------------------------------------------
// Real engine
// ---------------------------------------------------------------------------

/// The real in-process bwa-mem3 engine: a thin adapter over the merged
/// `bwa_mem3_rs` three-phase API. The `BwaIndex` and `MemOpts` are loaded/built
/// once at wire time and shared read-only across every worker (both are
/// `Send + Sync`).
pub(crate) struct BwaMem3Engine {
    idx: Arc<BwaIndex>,
    opts: Arc<MemOpts>,
    /// Wall time (ns) spent *inside* the three `bwa_mem3_rs` FFI calls, summed
    /// across every pool worker (the engine is shared behind an `Arc`, so these
    /// are atomics). The pipeline stats table already reports each *step's* total
    /// busy time (`AlignSeedExtend`/`CohortPeStat`/`AlignPairEmit`), which
    /// includes the fgumi-side arena fill, read-view conversion, and sink/template
    /// assembly around each FFI call. Subtracting these FFI totals from the step
    /// totals splits the cost into "bwa-mem3 C++" vs "fgumi-side overhead" — the
    /// question the `pair_emit`-vs-`seed_extend` finding raised. Reported once at
    /// [`Drop`] (i.e. pipeline end), at `debug`, so production runs stay quiet.
    seed_extend_ns: AtomicU64,
    seed_extend_calls: AtomicU64,
    infer_cohort_ns: AtomicU64,
    infer_cohort_calls: AtomicU64,
    pair_emit_ns: AtomicU64,
    pair_emit_calls: AtomicU64,
}

/// Nanoseconds since `start`, saturating into a `u64` (a `u64` of ns is ~584
/// years, so the cap is only a formality that avoids a `u128`-to-`u64` panic
/// path on an absurd clock).
fn elapsed_ns(start: Instant) -> u64 {
    u64::try_from(start.elapsed().as_nanos()).unwrap_or(u64::MAX)
}

impl BwaMem3Engine {
    /// Wrap a loaded index and built options. Both are shared (`Arc`) so worker
    /// copies of the steps clone the handles, not the data.
    pub(crate) fn new(idx: Arc<BwaIndex>, opts: Arc<MemOpts>) -> Self {
        Self {
            idx,
            opts,
            seed_extend_ns: AtomicU64::new(0),
            seed_extend_calls: AtomicU64::new(0),
            infer_cohort_ns: AtomicU64::new(0),
            infer_cohort_calls: AtomicU64::new(0),
            pair_emit_ns: AtomicU64::new(0),
            pair_emit_calls: AtomicU64::new(0),
        }
    }
}

impl Drop for BwaMem3Engine {
    fn drop(&mut self) {
        let se_ns = self.seed_extend_ns.load(Ordering::Relaxed);
        let ic_ns = self.infer_cohort_ns.load(Ordering::Relaxed);
        let pe_ns = self.pair_emit_ns.load(Ordering::Relaxed);
        if se_ns == 0 && ic_ns == 0 && pe_ns == 0 {
            return; // never ran (e.g. wiring failed, or a fake-index test path)
        }
        // Precision loss is irrelevant for a human-readable seconds display.
        #[allow(clippy::cast_precision_loss)]
        let secs = |ns: u64| ns as f64 / 1e9;
        log::debug!(
            "in-process bwa-mem3 FFI CPU (summed across workers): \
             seed_extend {:.2}s ({} calls), infer_cohort {:.3}s ({} calls), \
             pair_emit {:.2}s ({} calls). Each step's pipeline-stats total minus \
             this FFI time is fgumi-side overhead (arena fill / read-view build / \
             sink+template assembly).",
            secs(se_ns),
            self.seed_extend_calls.load(Ordering::Relaxed),
            secs(ic_ns),
            self.infer_cohort_calls.load(Ordering::Relaxed),
            secs(pe_ns),
            self.pair_emit_calls.load(Ordering::Relaxed),
        );
    }
}

/// A [`BwaMem3Engine`] sub-batch's share of its [`ResidentCohort`]: its pairs'
/// range and its singles' range, each absent when the sub-batch has none.
pub(crate) struct BwaRanges {
    pairs: Option<ResidentRange>,
    singles: Option<ResidentRange>,
}

impl AlignEngine for BwaMem3Engine {
    type Scratch = bwa_mem3_rs::AlignScratch;
    type Cohort = ResidentCohort;
    type Ranges = BwaRanges;
    type PeStat = MemPeStat;

    fn new_scratch(&self) -> anyhow::Result<Self::Scratch> {
        Ok(bwa_mem3_rs::AlignScratch::new()?)
    }

    fn new_cohort(&self) -> anyhow::Result<Self::Cohort> {
        Ok(ResidentCohort::new(self.opts.meth())?)
    }

    fn seed_extend(
        &self,
        scratch: &mut Self::Scratch,
        cohort: &Self::Cohort,
        batch: EngineBatch<'_>,
    ) -> anyhow::Result<Self::Ranges> {
        let t = Instant::now();
        // The C side copies each read into its resident slot during
        // `write_*`; the views here copy no sequence bytes.
        let mut ranges = BwaRanges { pairs: None, singles: None };
        if !batch.pairs.is_empty() {
            let mut range = cohort.reserve_pairs(batch.pairs.len())?;
            // One batch write takes the cohort's range lock once, not per read.
            // These borrowed views are rebuilt per call (one small allocation per
            // sub-batch, not per read): they borrow this sub-batch's arena, so
            // they cannot live in the `'static` per-worker scratch, and
            // `write_pairs` needs the whole range in one slice. bwa-mem3-rs then
            // collects its own same-sized Vec of C structs inside the call.
            let pairs: Vec<BwaReadPair<'_>> = batch
                .pairs
                .iter()
                .map(|p| BwaReadPair {
                    name_r1: p.r1.name,
                    seq_r1: p.r1.seq,
                    qual_r1: p.r1.qual,
                    name_r2: p.r2.name,
                    seq_r2: p.r2.seq,
                    qual_r2: p.r2.qual,
                })
                .collect();
            cohort.write_pairs(&mut range, &pairs)?;
            cohort.seed_extend(&self.idx, &self.opts, scratch, &mut range)?;
            ranges.pairs = Some(range);
        }
        if !batch.singles.is_empty() {
            let mut range = cohort.reserve_singles(batch.singles.len())?;
            let reads: Vec<BwaSingleRead<'_>> = batch
                .singles
                .iter()
                .map(|r| BwaSingleRead { name: r.name, seq: r.seq, qual: r.qual })
                .collect();
            cohort.write_singles(&mut range, &reads)?;
            cohort.seed_extend(&self.idx, &self.opts, scratch, &mut range)?;
            ranges.singles = Some(range);
        }
        self.seed_extend_ns.fetch_add(elapsed_ns(t), Ordering::Relaxed);
        self.seed_extend_calls.fetch_add(1, Ordering::Relaxed);
        Ok(ranges)
    }

    fn infer_cohort(&self, cohort: &Self::Cohort) -> anyhow::Result<Self::PeStat> {
        let t = Instant::now();
        let pestat = cohort.infer_cohort(&self.idx, &self.opts)?;
        self.infer_cohort_ns.fetch_add(elapsed_ns(t), Ordering::Relaxed);
        self.infer_cohort_calls.fetch_add(1, Ordering::Relaxed);
        Ok(pestat)
    }

    fn pair_emit(
        &self,
        scratch: &mut Self::Scratch,
        cohort: &Self::Cohort,
        ranges: Self::Ranges,
        pestat: Option<&Self::PeStat>,
        ids: IdBases,
        sink: &mut dyn RecordSink,
    ) -> anyhow::Result<()> {
        let t = Instant::now();
        // Pairs then singles, each tagged with its sub-batch-local origin
        // (`origin_base` 0) — the order and origins the legacy per-sub-batch
        // `pair_emit` produced.
        let BwaRanges { pairs, singles } = ranges;
        for mut range in [pairs, singles].into_iter().flatten() {
            cohort.pair_emit(
                &self.idx, &self.opts, scratch, &mut range, pestat, ids, 0, &mut *sink,
            )?;
        }
        self.pair_emit_ns.fetch_add(elapsed_ns(t), Ordering::Relaxed);
        self.pair_emit_calls.fetch_add(1, Ordering::Relaxed);
        Ok(())
    }

    fn n_pairs(ranges: &Self::Ranges) -> usize {
        ranges.pairs.as_ref().map_or(0, ResidentRange::n_pairs)
    }

    fn pestat_summary(pestat: &Self::PeStat) -> String {
        // FR is the orientation `mem_pestat`/`worker_sam` actually uses for
        // insert-size-aware mate rescue in the common (non-mate-pair) case;
        // upstream's `[M::mem_pestat]` line reports the same orientation.
        let fr = pestat.orientation(PeOrientation::Fr);
        format!(
            "FR insert mean={mean:.1} std={std:.1} low={low} high={high}{failed}",
            mean = fr.avg,
            std = fr.std,
            low = fr.low,
            high = fr.high,
            failed = if fr.failed { " (failed)" } else { "" },
        )
    }
}

#[cfg(test)]
pub(crate) mod fake {
    //! A deterministic, FFI-free [`AlignEngine`] for the scheduling/ordering tests.
    //!
    //! `seed_extend` records one [`TaggedRead`] per input read, in bwa-mem3's
    //! emission order — for each pair `i` an R1 then an R2 read tagged
    //! [`RecordOrigin::Pair(i)`], then for each single `j` one read tagged
    //! [`RecordOrigin::Single(j)`] (bwa emits pairs, R1 then R2, then singles).
    //! `pair_emit` replays that recorded order, emitting one
    //! **minimal-valid packed-BAM record body** per read through the sink tagged
    //! with its origin, so a scheduling test can assert records come back grouped
    //! by origin in input order — and, crucially for
    //! [`AlignPairEmitStep`](super::super::pair_emit), so
    //! [`Template::from_records`](crate::template::Template::from_records) can
    //! actually parse and group them and expose a `read_name()`/`flags()` for the
    //! name-match guard. Every read of a pair `i` gets the same synthetic name
    //! ([`fake_record_name`]) so its two records collapse into one template;
    //! `mate` picks `FIRST_SEGMENT`/`LAST_SEGMENT` so the pair has exactly one
    //! primary R1 and one primary R2. `infer_cohort` returns a unit `PeStat`; no
    //! FFI is called anywhere. The fake's cohort only counts the pairs seeded
    //! into it; its "ranges" are the recorded reads themselves.
    //!
    //! The body is a fixed 4-base record, so it is not usable to assert per-pair
    //! record *counts* (a real engine emits supplementaries; this fake never
    //! does) — only grouping, ordering, and guard-firing.

    use std::sync::atomic::{AtomicUsize, Ordering};

    use fgumi_raw_bam::{SamBuilder, flags};

    use super::{AlignEngine, EngineBatch, IdBases, RecordOrigin, RecordSink};

    /// The fake's resident cohort: how many pairs its sub-batches have seeded.
    #[derive(Debug, Default)]
    pub(crate) struct FakeCohort {
        pub(crate) seeded_pairs: AtomicUsize,
    }

    /// The synthetic queryname [`FakeEngine::pair_emit`] stamps on every record of
    /// an origin, so a pair's two records share one QNAME and collapse into one
    /// [`Template`](crate::template::Template) whose `name()` a test can predict
    /// (to line up — or deliberately mismatch — the unmapped half's name for the
    /// name-match guard). Both records of `Pair(i)` get the same name; each
    /// `Single(j)` gets its own.
    pub(crate) fn fake_record_name(origin: RecordOrigin) -> Vec<u8> {
        match origin {
            RecordOrigin::Pair(i) => format!("fakepair{i}").into_bytes(),
            RecordOrigin::Single(j) => format!("fakesingle{j}").into_bytes(),
        }
    }

    /// One recorded input read: its [`RecordOrigin`] (which input pair/single it
    /// came from) and `mate` (0 for R1 or a single, 1 for R2) so the body is
    /// deterministic and the R1/R2 order within a pair is inspectable.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub(crate) struct TaggedRead {
        pub(crate) origin: RecordOrigin,
        pub(crate) mate: u8,
    }

    /// The deterministic, FFI-free engine.
    #[derive(Clone, Copy, Debug, Default)]
    pub(crate) struct FakeEngine;

    impl AlignEngine for FakeEngine {
        type Scratch = ();
        type Cohort = FakeCohort;
        type Ranges = Vec<TaggedRead>;
        type PeStat = ();

        fn new_scratch(&self) -> anyhow::Result<Self::Scratch> {
            Ok(())
        }

        fn new_cohort(&self) -> anyhow::Result<Self::Cohort> {
            Ok(FakeCohort::default())
        }

        fn seed_extend(
            &self,
            _scratch: &mut Self::Scratch,
            cohort: &Self::Cohort,
            batch: EngineBatch<'_>,
        ) -> anyhow::Result<Self::Ranges> {
            cohort.seeded_pairs.fetch_add(batch.pairs.len(), Ordering::Relaxed);
            // Emission order: every pair (R1 then R2) in pair
            // index order, then every single in single index order.
            let mut reads = Vec::with_capacity(batch.pairs.len() * 2 + batch.singles.len());
            for i in 0..batch.pairs.len() {
                reads.push(TaggedRead { origin: RecordOrigin::Pair(i), mate: 0 });
                reads.push(TaggedRead { origin: RecordOrigin::Pair(i), mate: 1 });
            }
            for j in 0..batch.singles.len() {
                reads.push(TaggedRead { origin: RecordOrigin::Single(j), mate: 0 });
            }
            Ok(reads)
        }

        fn infer_cohort(&self, _cohort: &Self::Cohort) -> anyhow::Result<Self::PeStat> {
            Ok(())
        }

        fn pair_emit(
            &self,
            _scratch: &mut Self::Scratch,
            _cohort: &Self::Cohort,
            regs: Self::Ranges,
            _pestat: Option<&Self::PeStat>,
            _ids: IdBases,
            sink: &mut dyn RecordSink,
        ) -> anyhow::Result<()> {
            // Replay the recorded reads in emission order, one minimal-valid
            // packed-BAM record body per read tagged with its origin. Each record
            // is a real 4-base BAM record (built via `SamBuilder`, whose output is
            // a body with no `block_size` prefix — exactly a `RecordSink` body), so
            // `RawRecord::from(body)` + `Template::from_records` parse and group
            // them: a pair's two records share one name and split R1/R2 by `mate`.
            for read in &regs {
                let name = fake_record_name(read.origin);
                let record_flags = match read.origin {
                    RecordOrigin::Pair(_) => {
                        flags::PAIRED
                            | if read.mate == 0 {
                                flags::FIRST_SEGMENT
                            } else {
                                flags::LAST_SEGMENT
                            }
                    }
                    RecordOrigin::Single(_) => 0,
                };
                let mut builder = SamBuilder::new();
                builder
                    .read_name(&name)
                    .flags(record_flags)
                    .sequence(b"ACGT")
                    .qualities(&[30u8, 30, 30, 30]);
                let record = builder.build();
                sink.emit(read.origin, record.as_ref());
            }
            Ok(())
        }

        fn n_pairs(regs: &Self::Ranges) -> usize {
            // Each pair contributed two reads; count the R1 side (mate == 0) to
            // recover the number of PE templates.
            regs.iter().filter(|r| matches!(r.origin, RecordOrigin::Pair(_)) && r.mate == 0).count()
        }

        fn pestat_summary(_pestat: &Self::PeStat) -> String {
            // `PeStat = ()`: no real insert-size model to report.
            "n/a (FakeEngine)".to_string()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::fake::FakeEngine;
    use super::{AlignEngine, EngineBatch, EnginePair, EngineRead, IdBases};
    use bwa_mem3_rs::{RecordOrigin, RecordVec};

    /// Build a two-pair, one-single [`EngineBatch`] borrowing local buffers.
    fn sample_batch<'a>(
        pairs: &'a [EnginePair<'a>],
        singles: &'a [EngineRead<'a>],
    ) -> EngineBatch<'a> {
        EngineBatch { pairs, singles }
    }

    // `FakeEngine`'s `Scratch` and `PeStat` are the unit type, so binding the
    // values these driving calls return trips `let_unit_value`. Keeping the
    // bindings mirrors how a real step drives the engine, which is the point of
    // the test, so the lint is allowed locally.
    #[allow(clippy::let_unit_value)]
    #[test]
    fn fake_engine_emits_reads_grouped_by_origin_in_input_order() {
        let engine = FakeEngine;
        let mut scratch = engine.new_scratch().expect("scratch");

        let pairs = [
            EnginePair {
                r1: EngineRead::new(b"pair0", b"ACGT", None),
                r2: EngineRead::new(b"pair0", b"TGCA", None),
            },
            EnginePair {
                r1: EngineRead::new(b"pair1", b"AACC", None),
                r2: EngineRead::new(b"pair1", b"GGTT", None),
            },
        ];
        let singles = [EngineRead::new(b"sng0", b"NNNN", None)];
        let batch = sample_batch(&pairs, &singles);

        let cohort = engine.new_cohort().expect("cohort");
        let regs = engine.seed_extend(&mut scratch, &cohort, batch).expect("seed_extend");

        // n_pairs counts PE templates (not reads); the ranges hold one entry per read.
        assert_eq!(FakeEngine::n_pairs(&regs), 2);
        assert_eq!(regs.len(), 2 * 2 + 1);
        assert_eq!(cohort.seeded_pairs.load(std::sync::atomic::Ordering::Relaxed), 2);

        let pestat = engine.infer_cohort(&cohort).expect("infer_cohort");

        let ids = IdBases { first_single_id: 100, first_pair_id: 7 };
        let mut sink = RecordVec::default();
        engine
            .pair_emit(&mut scratch, &cohort, regs, Some(&pestat), ids, &mut sink)
            .expect("pair_emit");

        // Every input read produced exactly one record, in emission order:
        // pair 0 (R1, R2), pair 1 (R1, R2), then single 0.
        let origins: Vec<RecordOrigin> = sink.records.iter().map(|(o, _)| *o).collect();
        assert_eq!(
            origins,
            vec![
                RecordOrigin::Pair(0),
                RecordOrigin::Pair(0),
                RecordOrigin::Pair(1),
                RecordOrigin::Pair(1),
                RecordOrigin::Single(0),
            ],
        );
        // Bodies are non-empty and deterministic.
        assert!(sink.records.iter().all(|(_, body)| !body.is_empty()));
    }

    #[allow(clippy::let_unit_value)] // `FakeEngine::Scratch` is `()`; see the note above.
    #[test]
    fn fake_engine_handles_pure_single_and_empty_batches() {
        let engine = FakeEngine;
        let mut scratch = engine.new_scratch().expect("scratch");

        // Pure-single batch: no pairs, so a `None` pestat is valid.
        let singles =
            [EngineRead::new(b"sng0", b"ACG", None), EngineRead::new(b"sng1", b"GTA", None)];
        let batch = EngineBatch { pairs: &[], singles: &singles };
        let cohort = engine.new_cohort().expect("cohort");
        let regs = engine.seed_extend(&mut scratch, &cohort, batch).expect("seed_extend");
        assert_eq!(FakeEngine::n_pairs(&regs), 0);

        let mut sink = RecordVec::default();
        engine
            .pair_emit(&mut scratch, &cohort, regs, None, IdBases::default(), &mut sink)
            .expect("emit");
        let origins: Vec<RecordOrigin> = sink.records.iter().map(|(o, _)| *o).collect();
        assert_eq!(origins, vec![RecordOrigin::Single(0), RecordOrigin::Single(1)]);

        // Empty batch: no reads, no records.
        let empty = EngineBatch { pairs: &[], singles: &[] };
        let regs = engine.seed_extend(&mut scratch, &cohort, empty).expect("seed_extend");
        assert_eq!(regs.len(), 0);
        let mut sink = RecordVec::default();
        engine
            .pair_emit(&mut scratch, &cohort, regs, None, IdBases::default(), &mut sink)
            .expect("emit");
        assert!(sink.records.is_empty());
    }

    /// Smoke test for the real [`super::BwaMem3Engine`] over the full three-phase
    /// flow. Requires a prebuilt bwa-mem3 index at `FGUMI_BWA_MEM3_TEST_REF` (the
    /// index prefix). Skips cleanly when it is unset, unless
    /// `FGUMI_BWA_MEM3_REQUIRE_TOOLS` is set — as in the e2e-parity CI job — in
    /// which case a missing index is a hard failure rather than a silent pass.
    #[test]
    fn bwa_mem3_engine_smoke_pair_emit() {
        use std::sync::Arc;

        use super::BwaMem3Engine;

        let Ok(prefix) = std::env::var("FGUMI_BWA_MEM3_TEST_REF") else {
            assert!(
                std::env::var_os("FGUMI_BWA_MEM3_REQUIRE_TOOLS").is_none(),
                "FGUMI_BWA_MEM3_REQUIRE_TOOLS is set but FGUMI_BWA_MEM3_TEST_REF (a bwa-mem3 \
                 index prefix) is not"
            );
            eprintln!("skipping: FGUMI_BWA_MEM3_TEST_REF not set");
            return;
        };

        let idx = Arc::new(bwa_mem3_rs::BwaIndex::load(&prefix).expect("load index"));
        let mut opts = bwa_mem3_rs::MemOpts::new().expect("opts");
        opts.set_pe(true);
        let engine = BwaMem3Engine::new(idx, Arc::new(opts));

        let mut scratch = engine.new_scratch().expect("scratch");

        // One PE pair of short reads. The engine emits at least one record
        // (unmapped if nothing aligns) for each mate, all tagged `Pair(0)`.
        let pairs = [EnginePair {
            r1: EngineRead::new(b"read1", b"ACGTACGTACGTACGTACGT", Some(b"IIIIIIIIIIIIIIIIIIII")),
            r2: EngineRead::new(b"read1", b"TACGTACGTACGTACGTACG", Some(b"IIIIIIIIIIIIIIIIIIII")),
        }];
        let batch = EngineBatch { pairs: &pairs, singles: &[] };

        let cohort = engine.new_cohort().expect("cohort");
        let regs = engine.seed_extend(&mut scratch, &cohort, batch).expect("seed_extend");
        assert_eq!(BwaMem3Engine::n_pairs(&regs), 1, "one PE template was seeded");
        let pestat = engine.infer_cohort(&cohort).expect("infer_cohort");

        let ids = IdBases { first_single_id: 0, first_pair_id: 0 };
        let mut sink = RecordVec::default();
        engine
            .pair_emit(&mut scratch, &cohort, regs, Some(&pestat), ids, &mut sink)
            .expect("pair_emit");
        assert!(
            sink.records.len() >= 2,
            "expected >= 1 record per mate, got {}",
            sink.records.len()
        );
        assert!(
            sink.records.iter().all(|(origin, _)| *origin == RecordOrigin::Pair(0)),
            "every record is tagged with the one pair's origin"
        );
    }
}
