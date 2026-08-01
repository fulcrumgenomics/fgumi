//! Sort-order verification engine for `fgumi compare bams --command sort`.
//!
//! `sort` is different from every other comparison preset: the output *order* is the
//! payload, but two conforming sort implementations may legitimately break ties
//! differently. Coordinate sort leaves records with the same `(tid, pos, reverse)` as an
//! unordered set (samtools does not tie-break on read name), and fgumi's
//! template-coordinate sort tie-breaks the name lane with a hash while samtools/fgbio use
//! lexical order (the documented `SORT-01` residue). A plain index-positional byte compare
//! would report `DIFFER` on two files that are both correctly sorted but disagree only on
//! tie order — a false positive in the comparison oracle.
//!
//! This engine never re-sorts either input (`sort` is the one BAM comparison that must
//! not: re-sorting to force alignment would mask a genuine reordering regression). Instead
//! it:
//!
//! 1. Detects the shared sort order from both inputs' `@HD` `SO`/`GO`/`SS` header tags
//!    (see [`detect_sort_order`]) and errors if the two inputs disagree or declare an
//!    order this engine cannot verify.
//! 2. Reads each file exactly **once**, in file order, layering an [`OrderChecked`] stream
//!    adapter over that single keyed-record read. The adapter watches every key flow past
//!    and folds monotonicity into a per-file [`OrderTracker`] (violation count + first
//!    violation) as a side effect — this catches a genuine mis-sort in either file, and is
//!    the exact fold [`fgumi_sort::verify_sort_order`] performs in its own standalone pass.
//! 3. Over that same single pass, groups each file into maximal runs of records that share
//!    the same *core* sort key (the key without the name/hash tie-break lane — see
//!    [`fgumi_sort::TemplateKey::core_cmp`] for template-coordinate; coordinate has no
//!    tie-break lane beyond `(tid, pos, reverse)` in the first place), and asserts the
//!    multiset of full records is equal between the two files' corresponding runs. This
//!    tolerates intra-run reordering while still catching a missing/extra record or any
//!    content difference.
//!
//! The run comparison is exact — reporting the specific reads that went unmatched — for as
//! long as the two files' orderings stay within [`EXACT_WINDOW_RECORDS`] of each other. Past
//! that it degrades to an order-insensitive digest ([`MultisetDigest`]) that still decides
//! the same question in constant memory, but can no longer name individual reads. See
//! [`RunCanceller`] for why that fallback exists and what it costs.
//!
//! Reading each file once (rather than an order pass per file *plus* a reopened comparison
//! pass — four full BAM traversals) roughly halves the BAM I/O + BGZF decode this
//! benchmarking command spends on WGS-scale inputs. Because the comparison — including its
//! desync "drain" path — pulls records *through* the [`OrderChecked`] adapter, each file's
//! order is verified **completely**, even past a desync point, so the per-file violation
//! counts are identical to the old two-pass computation for every input.
//!
//! Both checks run independently and are `AND`ed together in
//! [`SortVerifyOutcome::is_match`]: an order violation in either file, or any run
//! mismatch, is a `DIFFER`.

use std::cmp::Ordering;
use std::path::Path;
use std::sync::LazyLock;

use anyhow::{Result, anyhow, bail};

use fgumi_raw_bam::{RawRecord, RawRecordView};
use fgumi_sort::{
    LibraryLookup, QuerynameComparator, RawQuerynameKey, RawQuerynameLexKey, RawSortKey,
    SortContext, SortOrder, TemplateKey, cb_hasher, extract_coordinate_key_inline,
    extract_template_key_inline,
};
use noodles::sam::Header;

use ahash::{AHashMap, RandomState};

use super::{CheckedRecords, OpenedInput};
use crate::sam::SamTag;

use super::super::raw_compare::content_key_exact;
use super::header::{compare_headers, fold_header_diffs};

/// Outcome of a [`sort_verify_compare`] run.
#[derive(Debug, Clone)]
pub struct SortVerifyOutcome {
    /// The sort order detected from both inputs' `@HD` header (they must agree).
    pub sort_order: SortOrder,
    /// Total number of records read from `bam1`.
    pub bam1_count: u64,
    /// Total number of records read from `bam2`.
    pub bam2_count: u64,
    /// Number of sort-order violations found in `bam1` alone (accumulated by the per-file
    /// order tracker, the inline fold of [`fgumi_sort::verify_sort_order`]).
    pub bam1_violations: u64,
    /// The first violation in `bam1`, if any: `(record_number, read_name)`.
    pub bam1_first_violation: Option<(u64, String)>,
    /// Number of sort-order violations found in `bam2` alone.
    pub bam2_violations: u64,
    /// The first violation in `bam2`, if any: `(record_number, read_name)`.
    pub bam2_first_violation: Option<(u64, String)>,
    /// Number of equal-core-sort-key runs whose record multiset differed between the two
    /// files (a content difference, or a missing/extra record within the run).
    pub run_mismatches: u64,
    /// Number of runs whose comparison fell back from the exact canceller to the
    /// order-insensitive digest because the two files' orderings diverged by more than
    /// this engine's exact-cancellation window (`EXACT_WINDOW_RECORDS`, enforced by
    /// `RunCanceller`). Both are private to this module, so they are named here rather
    /// than linked — an intra-doc link from this public field would not resolve.
    ///
    /// Purely informational, and deliberately *not* part of [`Self::is_match`]: degrading is
    /// not an error, and a run compared this way is still compared. It is reported so a
    /// reader can tell that an `IDENTICAL` verdict rested on a digest rather than on an
    /// exact multiset, and so a `DIFFER` verdict that names no reads is explicable.
    pub runs_compared_by_digest: u64,
    /// `true` if the two files' run sequences desynchronized — one file's run boundaries
    /// don't line up with the other's (a whole run present on only one side), or one file
    /// ran out of records before the other. Once this happens, no further runs are
    /// compared (no resync), though both streams are still drained for accurate counts.
    pub presence_mismatch: bool,
    /// `true` if the two inputs' `@HD`/`@SQ`/`@RG` headers disagreed on a field
    /// [`compare_headers`](super::header::compare_headers) considers significant (`@PG`/`@CO`
    /// are normalized and never contribute here). Note that `@HD` `SO`/`GO`/`SS` agreement is
    /// already implied by `detect_sort_order` succeeding on both inputs with a matching
    /// [`SortOrder`] (checked before this field is ever populated) *only when both writers use
    /// the same bare-vs-prefixed `SS` convention* — this field catches the case where they
    /// don't (see `R2-HDR-01`), plus any `@SQ`/`@RG` divergence `detect_sort_order` never looks
    /// at.
    pub header_mismatch: bool,
    /// Human-readable diff strings, capped at the caller-supplied `max_diffs`.
    pub diff_details: Vec<String>,
}

impl SortVerifyOutcome {
    /// Returns `true` iff both files are correctly ordered under the shared sort key, every
    /// equal-core-key run's record multiset matches between the two files, and no significant
    /// header divergence was found.
    #[must_use]
    pub fn is_match(&self) -> bool {
        self.bam1_violations == 0
            && self.bam2_violations == 0
            && self.run_mismatches == 0
            && !self.presence_mismatch
            && !self.header_mismatch
            && self.bam1_count == self.bam2_count
    }
}

/// Detect the sort order declared by a BAM header's `@HD` `SO`/`GO`/`SS` tags.
///
/// Accepts both the SAM-spec-compliant prefixed `SS` form (e.g. `queryname:natural`,
/// `queryname:lexicographical`, `unsorted:template-coordinate`) and fgumi's former
/// bare-form writer output (e.g. `natural`, `lexicographic`, `template-coordinate` —
/// see `R2-HDR-01`), so this detector keeps working across the writer fix. When a prefix
/// is present it is validated against the declared `SO` (via [`ss_subsort`]): a value
/// whose sort-order prefix *disagrees* with `SO` (e.g. `coordinate:natural` under
/// `SO:queryname`) is rejected rather than silently matched on its suffix. Both
/// `lexicographic` (fgumi's pre-fix bare spelling) and `lexicographical` (the
/// samtools/SAM-spec spelling fgumi now emits) are accepted as aliases for
/// [`QuerynameComparator::Lexicographic`].
///
/// # Errors
///
/// Returns an error if the header has no `@HD` line, no `SO` tag, or an `SO`/`GO`/`SS`
/// combination that isn't one of the three sort orders `fgumi sort` produces (i.e. the
/// file is unsorted, or declares an order this engine doesn't recognize).
/// Split an `@HD` `SS` value into its sub-sort token, validating any `<sort-order>:`
/// prefix against `expected_so`.
///
/// Accepts a bare sub-sort (`natural` — fgumi's R2-HDR-01 writer output) via `Ok(ss)`, or
/// the SAM-spec prefixed form (`queryname:natural`) via `Ok(rest)` when the prefix matches
/// the file's declared `SO`. Returns `Err(prefix)` when the value carries a sort-order
/// prefix that *disagrees* with `SO` — e.g. `coordinate:natural` under `SO:queryname` —
/// which the previous `rsplit(':')` suffix-only match silently accepted and then verified
/// with the wrong comparator.
fn ss_subsort<'a>(ss: &'a str, expected_so: &str) -> std::result::Result<&'a str, &'a str> {
    match ss.split_once(':') {
        None => Ok(ss),
        Some((prefix, rest)) if prefix == expected_so => Ok(rest),
        Some((prefix, _)) => Err(prefix),
    }
}

/// Map a BAM header's `@HD` `SO`/`GO`/`SS` tags to the single `SortOrder` enum value they
/// denote, normalizing fgbio's `unsorted:`-prefixed `SS` and fgumi's bare `SS` to the same
/// variant. This is the fgbio↔fgumi special-case builder.
pub(crate) fn sort_order_from_header(header: &Header) -> Result<SortOrder> {
    let hd = header
        .header()
        .ok_or_else(|| anyhow!("BAM header has no @HD line; cannot verify sort order"))?;
    let unknown = hd.other_fields();
    let so = unknown.get(b"SO").map(|v| String::from_utf8_lossy(v).into_owned());
    let go = unknown.get(b"GO").map(|v| String::from_utf8_lossy(v).into_owned());
    let ss = unknown.get(b"SS").map(|v| String::from_utf8_lossy(v).into_owned());

    match so.as_deref() {
        // This engine implements no coordinate sub-sort, so plain `SO:coordinate` with no
        // `SS` tag is accepted, but ANY declared `SS` sub-sort must be rejected rather than
        // silently validated as plain coordinate: `run_compare`'s equal-core-sort-key
        // grouping only ever verifies `(tid, pos, reverse)`, which is weaker than whatever a
        // declared sub-sort would promise. fgumi itself never emits a coordinate `SS` (see
        // `header_ss_tag` in `fgumi-sort`), so this cannot reject fgumi's own output.
        Some("coordinate") => match ss.as_deref().map(|s| ss_subsort(s, "coordinate")) {
            None => Ok(SortOrder::Coordinate),
            Some(Ok(bad_ss)) => bail!(
                "unrecognized/unsupported @HD SS sub-sort '{bad_ss}' for SO:coordinate \
                 (no coordinate sub-sort is currently implemented; expected no SS tag)"
            ),
            Some(Err(prefix)) => bail!(
                "@HD SS sort-order prefix '{prefix}' disagrees with SO:coordinate \
                 (expected no SS tag, since no coordinate sub-sort is currently implemented)"
            ),
        },
        // Validate the *whole* `SS` value, not just its suffix: a bare sub-sort or a
        // `queryname:`-prefixed one is accepted, but a value whose sort-order prefix
        // disagrees with `SO:queryname` (e.g. `coordinate:natural`) is rejected rather
        // than silently verified under the suffix's comparator.
        Some("queryname") => match ss.as_deref().map(|s| ss_subsort(s, "queryname")) {
            None => Ok(SortOrder::Queryname(QuerynameComparator::Lexicographic)),
            Some(Ok("natural")) => Ok(SortOrder::Queryname(QuerynameComparator::Natural)),
            Some(Ok("lexicographic" | "lexicographical")) => {
                Ok(SortOrder::Queryname(QuerynameComparator::Lexicographic))
            }
            Some(Ok(bad_ss)) => bail!(
                "unrecognized @HD SS sub-sort '{bad_ss}' for SO:queryname \
                 (expected 'lexicographical', 'lexicographic', or 'natural')"
            ),
            Some(Err(prefix)) => bail!(
                "@HD SS sort-order prefix '{prefix}' disagrees with SO:queryname \
                 (expected a bare sub-sort or a 'queryname:'-prefixed one)"
            ),
        },
        Some("unsorted") if go.as_deref() == Some("query") => {
            match ss.as_deref().map(|s| ss_subsort(s, "unsorted")) {
                Some(Ok("template-coordinate")) => Ok(SortOrder::TemplateCoordinate),
                _ => bail!(
                    "@HD SO:unsorted GO:query without SS:template-coordinate is not a sort \
                     order this engine can verify"
                ),
            }
        }
        Some(bad_so) => bail!(
            "unsupported @HD SO:'{bad_so}' for sort verification (expected 'coordinate', \
             'queryname', or 'unsorted' with GO:query and SS:template-coordinate)"
        ),
        None => bail!("BAM header has no @HD SO tag; cannot verify sort order"),
    }
}

/// Back-compat alias for existing internal callers.
pub(crate) fn detect_sort_order(header: &Header) -> Result<SortOrder> {
    sort_order_from_header(header)
}

/// A record together with its extracted sort key, read one step ahead of the run
/// currently being assembled (the lookahead needed to detect a run boundary).
type Keyed<K> = (K, RawRecord);

/// The per-file monotonic sort-order accumulator, mirroring the fold in
/// [`fgumi_sort::verify_sort_order`].
///
/// An [`OrderChecked`] adapter feeds every record's key through [`OrderTracker::observe`]
/// exactly once, in file order, so the accumulated violation count and first violation are
/// **complete** — even when the two compared streams desynchronize and the tail of one file
/// is drained rather than compared. `is_violation` is borrowed (both files' trackers share
/// one closure); the extracted key is cloned into `prev_key` to keep the yielded
/// `(key, record)` intact for the run comparison downstream.
struct OrderTracker<'a, K, IsViolation>
where
    IsViolation: Fn(&K, &K) -> bool,
{
    is_violation: &'a IsViolation,
    prev_key: Option<K>,
    total: u64,
    violations: u64,
    first_violation: Option<(u64, String)>,
}

impl<'a, K, IsViolation> OrderTracker<'a, K, IsViolation>
where
    K: Clone,
    IsViolation: Fn(&K, &K) -> bool,
{
    fn new(is_violation: &'a IsViolation) -> Self {
        Self { is_violation, prev_key: None, total: 0, violations: 0, first_violation: None }
    }

    /// Observe one record in file order: bump the running total and, if its key violates
    /// the sort invariant relative to the previous key, bump the violation count and (once)
    /// capture the first violation's `(record_number, read_name)`. This is exactly the
    /// per-record body of [`fgumi_sort::verify_sort_order`], including the 1-based record
    /// numbering and the `from_utf8_lossy` read-name extraction.
    fn observe(&mut self, key: &K, record_bytes: &[u8]) {
        self.total += 1;
        if let Some(prev) = &self.prev_key
            && (self.is_violation)(key, prev)
        {
            self.violations += 1;
            if self.first_violation.is_none() {
                let name = String::from_utf8_lossy(RawRecordView::new(record_bytes).read_name())
                    .to_string();
                self.first_violation = Some((self.total, name));
            }
        }
        self.prev_key = Some(key.clone());
    }
}

/// A keyed-record stream that order-checks every record it yields, as a side effect.
///
/// Wraps the single keyed-record read the run comparison already consumes, threading each
/// record's key through an [`OrderTracker`] before handing the `(key, record)` downstream
/// unchanged. Because the comparison — including its desync [`OrderChecked::drain`] path —
/// pulls records *through* this adapter, each file's order is verified completely even past
/// a desync (the semantic invariant the old end-to-end `verify_sort_order` pass guaranteed),
/// while the file is read only once.
struct OrderChecked<'a, K, ExtractKey, IsViolation>
where
    ExtractKey: Fn(&[u8]) -> K,
    IsViolation: Fn(&K, &K) -> bool,
{
    reader: CheckedRecords,
    extract_key: &'a ExtractKey,
    tracker: OrderTracker<'a, K, IsViolation>,
}

impl<'a, K, ExtractKey, IsViolation> OrderChecked<'a, K, ExtractKey, IsViolation>
where
    K: Clone,
    ExtractKey: Fn(&[u8]) -> K,
    IsViolation: Fn(&K, &K) -> bool,
{
    fn new(
        reader: CheckedRecords,
        extract_key: &'a ExtractKey,
        is_violation: &'a IsViolation,
    ) -> Self {
        Self { reader, extract_key, tracker: OrderTracker::new(is_violation) }
    }

    /// Read the next record, extract and order-check its key, and yield `(key, record)`;
    /// `None` at EOF.
    fn next_keyed(&mut self) -> Result<Option<Keyed<K>>> {
        match self.reader.next() {
            Some(rec) => {
                let rec = rec?;
                let key = (self.extract_key)(rec.as_ref());
                self.tracker.observe(&key, rec.as_ref());
                Ok(Some((key, rec)))
            }
            None => Ok(None),
        }
    }

    /// Drain the pending lookahead plus every remaining record *through* the order tracker.
    /// Used once the two streams have desynchronized so both the per-file violation counts
    /// and the tracker's record total (the authoritative count in [`SortVerifyOutcome`])
    /// still reflect each file end to end.
    ///
    /// The lookahead record was already read — and so already order-checked and counted —
    /// when [`Self::next_keyed`] produced it, so it is simply dropped here; only the
    /// still-unread tail is pulled through the tracker.
    fn drain(&mut self, next: &mut Option<Keyed<K>>) -> Result<()> {
        *next = None;
        while self.next_keyed()?.is_some() {}
        Ok(())
    }
}

/// Pull the next maximal run of records sharing the same core sort key from `stream`,
/// using `next` as one-record lookahead state (carried across calls). Returns `None`
/// once the stream (and lookahead) are exhausted. Every record read here flows through
/// `stream`'s order tracker exactly once.
fn next_in_run<K, ExtractKey, IsViolation>(
    stream: &mut OrderChecked<'_, K, ExtractKey, IsViolation>,
    next: &mut Option<Keyed<K>>,
    key: &K,
    same_core: &impl Fn(&K, &K) -> bool,
) -> Result<Option<RawRecord>>
where
    K: Clone,
    ExtractKey: Fn(&[u8]) -> K,
    IsViolation: Fn(&K, &K) -> bool,
{
    let Some((k, rec)) = next.take() else { return Ok(None) };
    if !same_core(&k, key) {
        // Belongs to the next run: put it back untouched so the run boundary is preserved.
        *next = Some((k, rec));
        return Ok(None);
    }
    *next = stream.next_keyed()?;
    Ok(Some(rec))
}

/// Consume the remainder of the run identified by `key`, returning how many records it held.
/// Counting only — used on the desync paths where one file ended first and the other's
/// current run just needs its length reported.
fn count_rest_of_run<K, ExtractKey, IsViolation>(
    stream: &mut OrderChecked<'_, K, ExtractKey, IsViolation>,
    next: &mut Option<Keyed<K>>,
    key: &K,
    same_core: &impl Fn(&K, &K) -> bool,
) -> Result<u64>
where
    K: Clone,
    ExtractKey: Fn(&[u8]) -> K,
    IsViolation: Fn(&K, &K) -> bool,
{
    let mut count = 0;
    while next_in_run(stream, next, key, same_core)?.is_some() {
        count += 1;
    }
    Ok(count)
}

/// Outcome of comparing one pair of corresponding equal-core-key runs.
struct RunPairOutcome {
    residual: RunResidual,
    count1: u64,
    count2: u64,
}

/// Compare the two files' records for the run identified by `key`, consuming both sides in
/// lockstep. A pair that is byte-identical is cancelled here directly and never reaches
/// [`RunCanceller`] — the common case, and the reason no content key is built for it.
/// Everything else (pairs that differ, and records arriving out of step) goes to the
/// canceller, which matches them by canonical content key as they arrive. Never
/// materializes the run, so an unmapped tail spanning the whole file costs memory
/// proportional only to how far the two orderings diverge.
///
/// Cancellation is confined to one run, but note *why* that is not what keeps it sound: a
/// record's run membership is a function of its content ([`content_key_exact`] covers the
/// core-field bytes carrying `ref_id`/`pos`/`flags`, from which the coordinate key
/// `(tid, pos, reverse)` — and likewise the queryname and template-coordinate core keys — is
/// derived). Two records with equal content keys therefore always fall in the same run, so a
/// record that moved to a different sort key can never find its old self to cancel against,
/// scoping or no scoping. Run scoping is what preserves the engine's *per-run* structure —
/// run-boundary misalignment detection, per-run record counts, the `run {n}` diff labels, and
/// the deliberate no-resync behavior — not what prevents a false MATCH.
fn compare_run<K, ExtractKey, IsViolation>(
    stream1: &mut OrderChecked<'_, K, ExtractKey, IsViolation>,
    next1: &mut Option<Keyed<K>>,
    stream2: &mut OrderChecked<'_, K, ExtractKey, IsViolation>,
    next2: &mut Option<Keyed<K>>,
    key: &K,
    same_core: &impl Fn(&K, &K) -> bool,
    exact_window: usize,
) -> Result<RunPairOutcome>
where
    K: Clone,
    ExtractKey: Fn(&[u8]) -> K,
    IsViolation: Fn(&K, &K) -> bool,
{
    let mut canceller = RunCanceller::new(exact_window);
    let (mut count1, mut count2) = (0u64, 0u64);

    loop {
        // Lockstep: taking one record from each side per iteration keeps pending at ~one
        // record when the two files agree on order, instead of draining one side first.
        let from1 = next_in_run(stream1, next1, key, same_core)?;
        let from2 = next_in_run(stream2, next2, key, same_core)?;
        if from1.is_none() && from2.is_none() {
            break;
        }
        // Byte-identical pair: cancel it here without building either record's
        // canonical content key.
        //
        // Sound because the run comparison decides on the symmetric difference of
        // the two sides' record multisets, and removing one element from each side
        // *with the same key* leaves that symmetric difference unchanged — whatever
        // else is pending. Byte equality implies content-key equality, and strictly
        // so: `content_key_exact` excludes `bin`, width-normalizes integer tags, and
        // sorts the tag multiset, so any two byte-identical records necessarily
        // share a key. The converse does not hold, which is why this is only a fast
        // path: content-equal records that differ in tag order, integer tag width,
        // or `bin` fall through to the key path below and still cancel there.
        //
        // Worth the special case because it is the overwhelmingly common one — two
        // files that agree, record for record — and `content_key_exact` is the
        // dominant remaining cost on the critical path (~32%: a `Vec` per aux tag,
        // a sort, a concatenation, and an `AHashMap` hash + insert, per record).
        // Here a `memcmp` replaces all of it.
        if let (Some(a), Some(b)) = (&from1, &from2)
            && a.as_ref() == b.as_ref()
        {
            count1 += 1;
            count2 += 1;
            continue;
        }
        if let Some(rec) = from1 {
            canceller.observe(Side::Bam1, rec.as_ref());
            count1 += 1;
        }
        if let Some(rec) = from2 {
            canceller.observe(Side::Bam2, rec.as_ref());
            count2 += 1;
        }
    }

    Ok(RunPairOutcome { residual: canceller.finish(), count1, count2 })
}

/// Records a single equal-core-key run comparison ([`RunCanceller`]) will hold pending,
/// across both sides, before it stops naming reads and switches to the digest.
///
/// Not a tuning knob for run *length*: a run of any length costs O(1) pending as long as the
/// two files list its records in the same order, because each record cancels against its
/// counterpart on arrival. Pending grows only with the *displacement* between the two
/// orderings, so this is reached when the inputs disagree about intra-run order beyond this
/// window — two independently-produced unmapped tails being the case that motivated it
/// (#699: a 77.3M-record both-ends-unmapped tail is one run under template-coordinate, and
/// two sorters order it differently).
///
/// This window buys *diagnostics*, not correctness: within it a mismatched run names the
/// reads that went unmatched, and past it [`RunCanceller`] still decides the same
/// question — in constant memory — but can only report counts. It is therefore sized to
/// cover any plausible genuine difference rather than any plausible displacement.
const EXACT_WINDOW_RECORDS: usize = 1_000_000;

/// Which input a record came from. Names the two symmetric halves of [`RunCanceller`] so
/// one implementation serves both directions.
#[derive(Debug, Clone, Copy)]
enum Side {
    Bam1,
    Bam2,
}

/// Two independently-seeded hashers whose outputs form the low and high halves of the
/// 128-bit per-record hash [`MultisetDigest`] accumulates.
///
/// The seeds are fixed constants rather than [`RandomState::new`]'s per-process random ones
/// so that repeated invocations of one binary agree. The verdict is order-insensitive
/// either way, since both sides are hashed by the same hashers within a run; what per-process
/// seeds would add is a re-rolled residual collision risk, letting the oracle change its mind
/// about the same pair of files between runs. Fixed seeds make that risk a property of the
/// inputs rather than of the run. (The digits below are the fractional part of pi, the usual
/// nothing-up-my-sleeve source; nothing depends on the particular values beyond the two
/// hashers being distinct.)
///
/// This is determinism, not a stable-hash guarantee. `ahash` explicitly does not fix its
/// output across releases or build configurations — its AES-backed and software paths are
/// selected at compile time and produce different values — so the same key hashes differently
/// under a different `ahash` version, target architecture, or feature set. Digest values are
/// accordingly only ever compared with each other inside a single process, and must never be
/// persisted, reported as an identity, or compared across builds.
static DIGEST_HASHERS: LazyLock<(RandomState, RandomState)> = LazyLock::new(|| {
    (
        RandomState::with_seeds(
            0x243f_6a88_85a3_08d3,
            0x1319_8a2e_0370_7344,
            0xa409_3822_299f_31d0,
            0x082e_fa98_ec4e_6c89,
        ),
        RandomState::with_seeds(
            0x4528_21e6_38d0_1377,
            0xbe54_66cf_34e9_0c6c,
            0xc0ac_29b7_c97c_50dd,
            0x3f84_d5b5_b547_0917,
        ),
    )
});

/// Hash a canonical content key to 128 bits by concatenating two independently-seeded
/// 64-bit hashes of it.
fn hash_content_key(key: &[u8]) -> u128 {
    let (low, high) = &*DIGEST_HASHERS;
    (u128::from(high.hash_one(key)) << 64) | u128::from(low.hash_one(key))
}

/// An order-insensitive fingerprint of one side's uncancelled records, in constant space.
///
/// Every field is a commutative fold over the per-record 128-bit content-key hash, which is
/// exactly what makes the accumulated value independent of the order records arrive in — the
/// property the exact canceller pays O(displacement) memory to achieve.
///
/// All three fields are load-bearing, and none is redundant:
///
/// - `count` alone catches a missing or extra record.
/// - `xor` alone does not: it is self-inverse, so `{A, A, B}` and `{B, C, C}` share both a
///   count and an xor while being different multisets. (Pinned by
///   `digest_distinguishes_xor_cancelling_multisets`.)
/// - `sum` covers that case, since duplicates add rather than annihilate, and `xor` in turn
///   covers sums that happen to coincide modulo 2¹²⁸.
///
/// Equality is therefore a probabilistic — not an exact — multiset equality: two genuinely
/// different runs agreeing on all three fields would be reported as equal. For BAM records
/// hashed this way that is a ~2⁻¹²⁸ event per run and is not reachable by construction from
/// non-adversarial input, but it is not zero, and it is the price of comparing a 77M-record
/// run at all rather than aborting on it.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
struct MultisetDigest {
    sum: u128,
    xor: u128,
    count: u64,
}

impl MultisetDigest {
    /// Fold one record's canonical content key into the digest.
    fn add(&mut self, key: &[u8]) {
        self.add_n(key, 1);
    }

    /// Fold `multiplicity` copies of one canonical content key into the digest, without
    /// hashing it more than once. Used when draining a pending map, whose entries already
    /// carry their multiplicity.
    fn add_n(&mut self, key: &[u8], multiplicity: usize) {
        let hash = hash_content_key(key);
        self.sum = self.sum.wrapping_add(hash.wrapping_mul(multiplicity as u128));
        // XOR is self-inverse, so an even number of copies contributes nothing.
        if multiplicity % 2 == 1 {
            self.xor ^= hash;
        }
        self.count += multiplicity as u64;
    }
}

/// One side's unmatched records for a single canonical content key, plus a read name for
/// diagnostics. `count` is the multiplicity (a BAM may legitimately hold several records
/// with one content key), and `name` is the read name of the first such record — kept once
/// per distinct key so a residual can be reported by name rather than by opaque key bytes.
#[derive(Debug)]
struct PendingEntry {
    count: usize,
    name: Vec<u8>,
}

/// What a run comparison had left over, in whichever of the two representations
/// [`RunCanceller`] finished in.
#[derive(Debug)]
enum RunResidual {
    /// The run stayed within the canceller's exact window: the exact symmetric difference of
    /// the two sides' record multisets, by read name.
    Exact { only_in_bam1: Vec<Vec<u8>>, only_in_bam2: Vec<Vec<u8>> },
    /// The run's displacement exceeded the window and it was finished by digest. The verdict
    /// is still decided; the individual unmatched records are not recoverable, because
    /// deciding this run at all meant not holding them. `window` is the window that was
    /// actually in force ([`EXACT_WINDOW_RECORDS`] in production), carried so the diagnostic
    /// quotes the threshold that fired rather than the production constant.
    Digest { matched: bool, window: usize },
}

/// Read names listed per side in a residual diagnostic before truncating. A mismatched run
/// can have arbitrarily many unmatched records; naming a few identifies the divergence
/// without letting one run emit an unbounded message.
const RESIDUAL_NAMES_SHOWN: usize = 5;

impl RunResidual {
    /// `true` when the two sides held equal record multisets — exactly, or (in the
    /// [`Self::Digest`] case) to the confidence [`MultisetDigest`] documents.
    fn is_match(&self) -> bool {
        match self {
            Self::Exact { only_in_bam1, only_in_bam2 } => {
                only_in_bam1.is_empty() && only_in_bam2.is_empty()
            }
            Self::Digest { matched, .. } => *matched,
        }
    }

    /// `true` when this run was finished by digest rather than by exact cancellation.
    fn compared_by_digest(&self) -> bool {
        matches!(self, Self::Digest { .. })
    }

    /// Render the unmatched read names per side, e.g. `; only in bam1: lonely_read`. Empty
    /// when nothing is left over. A digest run names no reads and says so, rather than
    /// rendering an empty string that would read as "nothing was left over".
    fn describe(&self) -> String {
        let (only_in_bam1, only_in_bam2) = match self {
            Self::Exact { only_in_bam1, only_in_bam2 } => (only_in_bam1, only_in_bam2),
            Self::Digest { window, .. } => {
                return format!(
                    "; compared by order-insensitive digest (the two inputs' orderings \
                     diverged by more than {window} records within this run), so the \
                     unmatched records cannot be named"
                );
            }
        };
        let side = |label: &str, names: &[Vec<u8>]| -> String {
            if names.is_empty() {
                return String::new();
            }
            let shown: Vec<String> = names
                .iter()
                .take(RESIDUAL_NAMES_SHOWN)
                .map(|name| String::from_utf8_lossy(name).into_owned())
                .collect();
            let more = names.len().saturating_sub(shown.len());
            let suffix = if more > 0 { format!(" (+{more} more)") } else { String::new() };
            format!("; only in {label}: {}{suffix}", shown.join(", "))
        };
        format!("{}{}", side("bam1", only_in_bam1), side("bam2", only_in_bam2))
    }
}

/// The exact half of [`RunCanceller`]: unmatched records held per side, keyed by canonical
/// content key, so a residual can name the reads involved.
#[derive(Debug, Default)]
struct ExactPending {
    pending1: AHashMap<Vec<u8>, PendingEntry>,
    pending2: AHashMap<Vec<u8>, PendingEntry>,
    pending1_len: usize,
    pending2_len: usize,
}

impl ExactPending {
    /// Records held across both sides.
    fn len(&self) -> usize {
        self.pending1_len + self.pending2_len
    }

    /// Cancel `record` against a waiting counterpart on the opposite side if one is there,
    /// otherwise hold it on its own side.
    fn cancel_or_hold(&mut self, side: Side, key: Vec<u8>, record: &[u8]) {
        let (opposite, opposite_len, own, own_len) = match side {
            Side::Bam1 => (
                &mut self.pending2,
                &mut self.pending2_len,
                &mut self.pending1,
                &mut self.pending1_len,
            ),
            Side::Bam2 => (
                &mut self.pending1,
                &mut self.pending1_len,
                &mut self.pending2,
                &mut self.pending2_len,
            ),
        };
        if let Some(entry) = opposite.get_mut(&key) {
            entry.count -= 1;
            if entry.count == 0 {
                opposite.remove(&key);
            }
            *opposite_len -= 1;
            return;
        }
        let name = RawRecordView::new(record).read_name().to_vec();
        own.entry(key)
            .and_modify(|entry| entry.count += 1)
            .or_insert(PendingEntry { count: 1, name });
        *own_len += 1;
    }

    /// Fold each side's held records into that side's digest, releasing the maps.
    fn into_digests(self) -> (MultisetDigest, MultisetDigest) {
        let fold = |pending: AHashMap<Vec<u8>, PendingEntry>| {
            let mut digest = MultisetDigest::default();
            for (key, entry) in pending {
                digest.add_n(&key, entry.count);
            }
            digest
        };
        (fold(self.pending1), fold(self.pending2))
    }

    /// The unmatched read names per side, sorted for deterministic diagnostics.
    fn into_residual(self) -> RunResidual {
        let names = |pending: AHashMap<Vec<u8>, PendingEntry>| -> Vec<Vec<u8>> {
            let mut names: Vec<Vec<u8>> = pending
                .into_values()
                .flat_map(|entry| std::iter::repeat_n(entry.name, entry.count))
                .collect();
            // Hash-map iteration order is unspecified; sort so diagnostics are deterministic
            // across runs and platforms.
            names.sort_unstable();
            names
        };
        RunResidual::Exact {
            only_in_bam1: names(self.pending1),
            only_in_bam2: names(self.pending2),
        }
    }
}

/// Whether a run comparison is still tracking individual records or has fallen back to
/// per-side digests.
#[derive(Debug)]
enum CancellerState {
    Exact(ExactPending),
    Digest { digest1: MultisetDigest, digest2: MultisetDigest },
}

/// Streaming multiset comparison of one equal-core-key run, in memory bounded by
/// [`EXACT_WINDOW_RECORDS`] rather than by the run's length or by how far the two files'
/// orderings diverge.
///
/// Each record that reaches the canceller is reduced to its canonical content key
/// (byte-identical pairs are cancelled by [`compare_run`] before this point — see its docs)
/// ([`content_key_exact`], whose byte-equality is *exactly* `Exact` content-equality) and
/// cancelled against the opposite side's pending set if a counterpart is waiting; otherwise
/// it joins its own side's pending set. Two files listing a run's records in the same order
/// therefore never hold more than one record pending, however long the run.
///
/// # Degrading to the digest
///
/// Files that list a run's records in *different* orders hold records pending until their
/// counterparts arrive, so pending grows with the displacement between the two orderings.
/// A BAM's both-ends-unmapped tail is a single run — every such record packs to one constant
/// coordinate key, and to one core template-coordinate key — and two independently-produced
/// sorts of a whole-genome BAM order that tail arbitrarily differently: #699 reports 77.3M
/// such records in one run of a 1.33B-record file. Holding that displacement exactly is not
/// affordable, and re-sorting the inputs into a common order (the remedy this cap used to
/// suggest) is *circular* when the two inputs are the outputs of the two sorters under
/// comparison — there is no third order to normalize to.
///
/// So on the observation that pushes pending past [`EXACT_WINDOW_RECORDS`], both pending sets
/// are folded into per-side [`MultisetDigest`]s, the maps are dropped, and the rest of the
/// run accumulates directly into those digests. This is sound because the comparison decides
/// the *symmetric difference* of the two sides' multisets: writing `A` and `B` for the two
/// sides' records and `C` for the pairs already cancelled — each of which removed one element
/// from `A` and an equal element from `B` — the digests cover exactly `A \ C` and `B \ C`, and
/// `A == B` iff `A \ C == B \ C`. Work done before the fallback is carried across it rather
/// than discarded, so a partial-displacement run still cancels everything it can exactly.
///
/// What is lost is *naming*: past the window a mismatched run can be reported only by count.
/// What is gained is that it can be reported at all. See [`MultisetDigest`] for the residual
/// (~2⁻¹²⁸) chance that a digest calls two different runs equal.
#[derive(Debug)]
struct RunCanceller {
    state: CancellerState,
    exact_window: usize,
    peak_pending: usize,
}

impl RunCanceller {
    fn new(exact_window: usize) -> Self {
        Self {
            state: CancellerState::Exact(ExactPending::default()),
            exact_window,
            peak_pending: 0,
        }
    }

    /// Observe one record from `side`: cancel it against a waiting counterpart, hold it, or —
    /// once this run has degraded — fold it into that side's digest.
    fn observe(&mut self, side: Side, record: &[u8]) {
        let key = content_key_exact(record);
        match &mut self.state {
            CancellerState::Exact(exact) => {
                exact.cancel_or_hold(side, key, record);
                self.peak_pending = self.peak_pending.max(exact.len());
                if exact.len() > self.exact_window {
                    self.degrade_to_digest();
                }
            }
            CancellerState::Digest { digest1, digest2 } => match side {
                Side::Bam1 => digest1.add(&key),
                Side::Bam2 => digest2.add(&key),
            },
        }
    }

    /// Fold both pending sets into per-side digests and release the maps. A no-op if this
    /// run has already degraded.
    ///
    /// The already-degraded arm has to hand the accumulated digests back, not let them
    /// drop: the placeholder `mem::replace` installs is empty, and two empty digests
    /// compare equal, so dropping them would turn a run of unmatched records into a
    /// `matched: true` from [`finish`](Self::finish) — a silent false `IDENTICAL`.
    fn degrade_to_digest(&mut self) {
        let placeholder = CancellerState::Digest {
            digest1: MultisetDigest::default(),
            digest2: MultisetDigest::default(),
        };
        self.state = match std::mem::replace(&mut self.state, placeholder) {
            CancellerState::Exact(exact) => {
                let (digest1, digest2) = exact.into_digests();
                CancellerState::Digest { digest1, digest2 }
            }
            already_degraded @ CancellerState::Digest { .. } => already_degraded,
        };
    }

    /// Highest number of records held pending at any point, across both sides. Exists so
    /// tests can assert the bounded-memory contract directly (an in-order run must peak at
    /// ~1 record however long it is, and no run may exceed the window by more than the one
    /// record that trips it) rather than inferring it from process RSS.
    #[cfg(test)]
    fn peak_pending(&self) -> usize {
        self.peak_pending
    }

    /// Consume the canceller, reporting what never found a counterpart.
    fn finish(self) -> RunResidual {
        let window = self.exact_window;
        match self.state {
            CancellerState::Exact(exact) => exact.into_residual(),
            CancellerState::Digest { digest1, digest2 } => {
                RunResidual::Digest { matched: digest1 == digest2, window }
            }
        }
    }
}

/// Result of the run-grouped multiset comparison pass (before the per-file sort-order
/// violation counts are merged in by [`run_full_verify`]).
#[derive(Debug, Default)]
struct RunCompareOutcome {
    run_mismatches: u64,
    runs_compared_by_digest: u64,
    presence_mismatch: bool,
    diff_details: Vec<String>,
}

/// Walk `stream1`/`stream2` in file order, grouping each into maximal equal-core-key runs
/// and asserting the record multiset matches within each pair of corresponding runs (see
/// the module docs). Never resyncs: once a run boundary or stream-length mismatch is
/// found, no further runs are compared, though both streams are drained *through their
/// order trackers* for accurate counts and complete per-file order verification.
fn run_compare<K, ExtractKey, IsViolation>(
    stream1: &mut OrderChecked<'_, K, ExtractKey, IsViolation>,
    stream2: &mut OrderChecked<'_, K, ExtractKey, IsViolation>,
    same_core: &impl Fn(&K, &K) -> bool,
    max_diffs: usize,
    exact_window: usize,
) -> Result<RunCompareOutcome>
where
    K: Clone,
    ExtractKey: Fn(&[u8]) -> K,
    IsViolation: Fn(&K, &K) -> bool,
{
    let mut out = RunCompareOutcome::default();

    let mut next1 = stream1.next_keyed()?;
    let mut next2 = stream2.next_keyed()?;

    let mut run_index = 0u64;

    loop {
        // Peek each side's current run key without consuming the lookahead record.
        let key1 = next1.as_ref().map(|(k, _)| k.clone());
        let key2 = next2.as_ref().map(|(k, _)| k.clone());

        match (key1, key2) {
            (None, None) => break,
            (Some(key1), None) => {
                let extra = count_rest_of_run(stream1, &mut next1, &key1, same_core)?;
                out.presence_mismatch = true;
                if out.diff_details.len() < max_diffs {
                    out.diff_details.push(format!(
                        "run {run_index}: bam1 has {extra} more record(s) than bam2 \
                         (bam2 exhausted first — no resync)"
                    ));
                }
                stream1.drain(&mut next1)?;
                break;
            }
            (None, Some(key2)) => {
                let extra = count_rest_of_run(stream2, &mut next2, &key2, same_core)?;
                out.presence_mismatch = true;
                if out.diff_details.len() < max_diffs {
                    out.diff_details.push(format!(
                        "run {run_index}: bam2 has {extra} more record(s) than bam1 \
                         (bam1 exhausted first — no resync)"
                    ));
                }
                stream2.drain(&mut next2)?;
                break;
            }
            (Some(key1), Some(key2)) => {
                if !same_core(&key1, &key2) {
                    out.presence_mismatch = true;
                    if out.diff_details.len() < max_diffs {
                        out.diff_details.push(format!(
                            "run {run_index}: sort-key run boundary mismatch — bam1 and \
                             bam2 groups do not align (no resync)"
                        ));
                    }
                    stream1.drain(&mut next1)?;
                    stream2.drain(&mut next2)?;
                    break;
                }

                let pair = compare_run(
                    stream1,
                    &mut next1,
                    stream2,
                    &mut next2,
                    &key1,
                    same_core,
                    exact_window,
                )?;
                if pair.residual.compared_by_digest() {
                    out.runs_compared_by_digest += 1;
                }
                if !pair.residual.is_match() {
                    out.run_mismatches += 1;
                    if out.diff_details.len() < max_diffs {
                        out.diff_details.push(format!(
                            "run {run_index}: record multiset differs ({} record(s) in \
                             bam1 vs {} in bam2){}",
                            pair.count1,
                            pair.count2,
                            pair.residual.describe()
                        ));
                    }
                }
            }
        }

        run_index += 1;
    }

    Ok(out)
}

/// Grouped parameters for [`run_full_verify`].
///
/// Bundled into a struct (rather than nine positional arguments) purely to keep the
/// function signature readable — this is a structural fix for `clippy::too_many_arguments`,
/// not a suppression of it.
struct VerifyContext<K, ExtractKey, IsViolation, SameCore>
where
    ExtractKey: Fn(&[u8]) -> K,
    IsViolation: Fn(&K, &K) -> bool,
    SameCore: Fn(&K, &K) -> bool,
{
    /// Record reader over the first input, already past its header. Carried rather
    /// than the path so the input is opened exactly once — see [`OpenedInput`].
    reader1: CheckedRecords,
    /// Record reader over the second input, already past its header.
    reader2: CheckedRecords,
    extract_key: ExtractKey,
    is_violation: IsViolation,
    same_core: SameCore,
    max_diffs: usize,
    order: SortOrder,
    /// Records a single run comparison may hold pending before degrading to the digest —
    /// [`EXACT_WINDOW_RECORDS`] in production, threaded through so tests can reach the
    /// fallback with a handful of records instead of a million.
    exact_window: usize,
}

/// Verify both files in a single synchronized streaming pass for a single key type `K` and
/// assemble the combined [`SortVerifyOutcome`].
///
/// Each file is read exactly once through an [`OrderChecked`] adapter that folds the
/// per-file monotonic-order check ([`OrderTracker`]) into the same read the run-grouped
/// multiset comparison consumes. The per-file violation counts and first violations are
/// read out of the two adapters' trackers after the comparison, and are complete even when
/// the streams desynchronize because [`run_compare`] drains the tail *through* the trackers
/// (see [`OrderChecked::drain`]).
fn run_full_verify<K, ExtractKey, IsViolation, SameCore>(
    ctx: VerifyContext<K, ExtractKey, IsViolation, SameCore>,
) -> Result<SortVerifyOutcome>
where
    K: Clone,
    ExtractKey: Fn(&[u8]) -> K,
    IsViolation: Fn(&K, &K) -> bool,
    SameCore: Fn(&K, &K) -> bool,
{
    let VerifyContext {
        reader1,
        reader2,
        extract_key,
        is_violation,
        same_core,
        max_diffs,
        order,
        exact_window,
    } = ctx;

    let mut stream1 = OrderChecked::new(reader1, &extract_key, &is_violation);
    let mut stream2 = OrderChecked::new(reader2, &extract_key, &is_violation);

    let runs = run_compare(&mut stream1, &mut stream2, &same_core, max_diffs, exact_window)?;

    Ok(SortVerifyOutcome {
        sort_order: order,
        bam1_count: stream1.tracker.total,
        bam2_count: stream2.tracker.total,
        bam1_violations: stream1.tracker.violations,
        bam1_first_violation: stream1.tracker.first_violation,
        bam2_violations: stream2.tracker.violations,
        bam2_first_violation: stream2.tracker.first_violation,
        run_mismatches: runs.run_mismatches,
        runs_compared_by_digest: runs.runs_compared_by_digest,
        presence_mismatch: runs.presence_mismatch,
        header_mismatch: false,
        diff_details: runs.diff_details,
    })
}

/// Verify `bam1` (fgumi's output) and `bam2` (the baseline) are both correctly sorted
/// under their shared, header-declared sort order, and compare them as a multiset grouped
/// by maximal equal-core-sort-key run (see the module docs).
///
/// Never re-sorts either input. Also compares the two inputs' headers via
/// [`compare_headers`](super::header::compare_headers) (`@HD`/`@SQ`/`@RG`, normalizing away
/// `@PG`/`@CO`); a significant divergence is folded into
/// [`SortVerifyOutcome::header_mismatch`]/[`SortVerifyOutcome::is_match`] alongside the
/// run-comparison findings. `max_diffs` caps the number of entries collected in
/// [`SortVerifyOutcome::diff_details`].
///
/// # Errors
///
/// Returns an error if either file cannot be opened or read, if the two files' `@HD`
/// headers declare different sort orders, or if the declared sort order isn't one this
/// engine can verify (see `detect_sort_order`).
pub fn sort_verify_compare(
    bam1: &Path,
    bam2: &Path,
    max_diffs: usize,
) -> Result<SortVerifyOutcome> {
    // One open per input, which is what the record readers below then stream —
    // see `OpenedInput::open`.
    sort_verify_compare_opened(OpenedInput::open(bam1)?, OpenedInput::open(bam2)?, max_diffs)
}

/// [`sort_verify_compare`] over inputs the caller already opened.
///
/// The entry point for `CompareBams::execute`, which has to read both headers to
/// check them for compatibility before dispatching here, and would otherwise have
/// to reopen both paths for their records — see [`OpenedInput`]. Handing the
/// streams over instead is what lets `compare bams --command sort` take a FIFO or
/// a process substitution.
///
/// # Errors
///
/// As [`sort_verify_compare`], minus the opens the caller already made.
pub(crate) fn sort_verify_compare_opened(
    input1: OpenedInput,
    input2: OpenedInput,
    max_diffs: usize,
) -> Result<SortVerifyOutcome> {
    sort_verify_compare_opened_with_window(input1, input2, max_diffs, EXACT_WINDOW_RECORDS)
}

/// [`sort_verify_compare_opened`] with an explicit exact-comparison window.
///
/// Exists so the digest fallback is reachable in a test: it is otherwise a million records
/// of displacement away, which no test fixture can afford to build. Production callers use
/// [`sort_verify_compare_opened`], which supplies [`EXACT_WINDOW_RECORDS`]; the window is
/// deliberately *not* a CLI flag, since exceeding it is no longer an error the caller has to
/// do anything about.
///
/// # Errors
///
/// As [`sort_verify_compare_opened`].
pub(crate) fn sort_verify_compare_opened_with_window(
    input1: OpenedInput,
    input2: OpenedInput,
    max_diffs: usize,
    exact_window: usize,
) -> Result<SortVerifyOutcome> {
    let OpenedInput { reader: reader1, header: header1, path: bam1 } = input1;
    let OpenedInput { reader: reader2, header: header2, path: bam2 } = input2;
    let (bam1, bam2) = (bam1.as_path(), bam2.as_path());

    let order1 = detect_sort_order(&header1)
        .map_err(|e| anyhow!("detecting sort order from {}: {e}", bam1.display()))?;
    let order2 = detect_sort_order(&header2)
        .map_err(|e| anyhow!("detecting sort order from {}: {e}", bam2.display()))?;
    if order1 != order2 {
        bail!(
            "{} and {} declare different sort orders ({order1:?} vs {order2:?}); cannot \
             verify sort order across a mismatched pair",
            bam1.display(),
            bam2.display()
        );
    }
    let order = order1;
    let header_diffs = compare_headers(&header1, &header2);

    // IMPORTANT: The per-SortOrder extractor/comparator selection below must stay
    // consistent with the matching arms in `OrderCheck::new`. Future changes to how an
    // arm selects its key extractor or comparator must be applied to both sites.
    let mut outcome = match order {
        SortOrder::Coordinate => {
            let nref = header1.reference_sequences().len() as u32;
            run_full_verify(VerifyContext {
                reader1,
                reader2,
                extract_key: |bam: &[u8]| extract_coordinate_key_inline(bam, nref),
                is_violation: |key: &u64, prev: &u64| key < prev,
                same_core: |a: &u64, b: &u64| a == b,
                max_diffs,
                order,
                exact_window,
            })
        }
        SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
            let ctx = SortContext::from_header(&header1);
            run_full_verify(VerifyContext {
                reader1,
                reader2,
                extract_key: move |bam: &[u8]| RawQuerynameLexKey::extract(bam, &ctx),
                is_violation: |key: &RawQuerynameLexKey, prev: &RawQuerynameLexKey| key < prev,
                same_core: |a: &RawQuerynameLexKey, b: &RawQuerynameLexKey| a == b,
                max_diffs,
                order,
                exact_window,
            })
        }
        SortOrder::Queryname(QuerynameComparator::Natural) => {
            let ctx = SortContext::from_header(&header1);
            run_full_verify(VerifyContext {
                reader1,
                reader2,
                extract_key: move |bam: &[u8]| RawQuerynameKey::extract(bam, &ctx),
                is_violation: |key: &RawQuerynameKey, prev: &RawQuerynameKey| key < prev,
                same_core: |a: &RawQuerynameKey, b: &RawQuerynameKey| a == b,
                max_diffs,
                order,
                exact_window,
            })
        }
        SortOrder::TemplateCoordinate => {
            let lib_lookup = LibraryLookup::from_header(&header1);
            let hasher = cb_hasher();
            // Matches `fgumi sort`'s own `--verify` (crate::commands::sort::parse_cell_tag):
            // template-coordinate always hashes the CB tag into the sort key when present.
            let cell_tag = Some(SamTag::CB);
            run_full_verify(VerifyContext {
                reader1,
                reader2,
                extract_key: move |bam: &[u8]| {
                    extract_template_key_inline(bam, &lib_lookup, cell_tag, &hasher)
                },
                is_violation: |key: &TemplateKey, prev: &TemplateKey| {
                    key.core_cmp(prev) == Ordering::Less
                },
                same_core: |a: &TemplateKey, b: &TemplateKey| a.core_cmp(b) == Ordering::Equal,
                max_diffs,
                order,
                exact_window,
            })
        }
    }?;

    fold_header_diffs(
        header_diffs,
        &mut outcome.header_mismatch,
        &mut outcome.diff_details,
        max_diffs,
    );

    Ok(outcome)
}

/// A per-file sort-order check folded into a comparison's own record pass, instead
/// of costing a dedicated traversal of the file.
///
/// `content` mode used to answer this question by opening each input and reading it
/// end to end *before* `positional_compare` streamed both files again — four full
/// BAM traversals for a comparison that needs two. This type lets the comparison
/// pass feed the records it is already decoding through the identical check.
///
/// The violation *count* is why this accumulates rather than failing on the first
/// bad record: the diagnostic reports how many records violated the order, which is
/// only known once the whole file has been seen. `positional_compare` already
/// drains both inputs to completion for its record counts, so the totals match a
/// dedicated pass exactly, and evaluating bam1's result before bam2's preserves
/// which file a mis-sorted pair names.
/// Extracts a record's sort key and reports whether it regressed against the
/// previous record, owning the previous key internally.
///
/// Boxed rather than generic so a comparison loop can hold one checker per input
/// without the loop itself becoming generic over the four sort-key types.
type OrderStep = Box<dyn FnMut(&[u8]) -> bool + Send>;

pub(crate) struct OrderCheck {
    order: SortOrder,
    step: OrderStep,
    total: u64,
    violations: u64,
    first_violation: Option<(u64, String)>,
}

impl OrderCheck {
    /// Build a checker for `order`, taking whatever context that order's key
    /// extractor needs (reference count, read-group library lookup) from `header`.
    ///
    /// IMPORTANT: the per-`SortOrder` extractor/comparator selection here must stay
    /// consistent with the matching arms in `sort_verify_compare`.
    pub(crate) fn new(header: &Header, order: SortOrder) -> Self {
        let step: OrderStep = match order {
            SortOrder::Coordinate => {
                let nref = header.reference_sequences().len() as u32;
                let mut prev: Option<u64> = None;
                Box::new(move |bam: &[u8]| {
                    let key = extract_coordinate_key_inline(bam, nref);
                    let violated = prev.is_some_and(|p| key < p);
                    prev = Some(key);
                    violated
                })
            }
            SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
                let ctx = SortContext::from_header(header);
                let mut prev: Option<RawQuerynameLexKey> = None;
                Box::new(move |bam: &[u8]| {
                    let key = RawQuerynameLexKey::extract(bam, &ctx);
                    let violated = prev.as_ref().is_some_and(|p| key < *p);
                    prev = Some(key);
                    violated
                })
            }
            SortOrder::Queryname(QuerynameComparator::Natural) => {
                let ctx = SortContext::from_header(header);
                let mut prev: Option<RawQuerynameKey> = None;
                Box::new(move |bam: &[u8]| {
                    let key = RawQuerynameKey::extract(bam, &ctx);
                    let violated = prev.as_ref().is_some_and(|p| key < *p);
                    prev = Some(key);
                    violated
                })
            }
            SortOrder::TemplateCoordinate => {
                let lib_lookup = LibraryLookup::from_header(header);
                let hasher = cb_hasher();
                let cell_tag = Some(SamTag::CB);
                let mut prev: Option<TemplateKey> = None;
                Box::new(move |bam: &[u8]| {
                    let key = extract_template_key_inline(bam, &lib_lookup, cell_tag, &hasher);
                    let violated = prev.as_ref().is_some_and(|p| key.core_cmp(p) == Ordering::Less);
                    prev = Some(key);
                    violated
                })
            }
        };
        Self { order, step, total: 0, violations: 0, first_violation: None }
    }

    /// Observe one record, in file order.
    pub(crate) fn observe(&mut self, bam: &[u8]) {
        self.total += 1;
        if (self.step)(bam) {
            self.violations += 1;
            if self.first_violation.is_none() {
                let name = String::from_utf8_lossy(RawRecordView::new(bam).read_name()).to_string();
                self.first_violation = Some((self.total, name));
            }
        }
    }

    /// Fail if any record violated the declared order.
    ///
    /// # Errors
    ///
    /// Returns an error naming the violation count and the first violating record's
    /// 1-based position and read name.
    pub(crate) fn into_result(self, path: &Path) -> Result<()> {
        if self.violations > 0 {
            let (record_num, name) = self
                .first_violation
                .expect("a nonzero violation count must carry a first violation");
            let (order, violations) = (self.order, self.violations);
            bail!(
                "{}: {violations} record(s) violate the declared {order:?} sort order (first at \
                 record {record_num}, read name '{name}'); records must actually be in {order:?} \
                 order, not merely declare it, for order-dependent comparison",
                path.display()
            );
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use bstr::BString;
    use fgumi_raw_bam::{SamBuilder, flags};
    use noodles::sam::alignment::io::Write as AlignmentWrite;
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::Header as HeaderRecord;
    use noodles::sam::header::record::value::map::header::tag as hd_tag;
    use rstest::rstest;

    use super::*;

    /// Builds a minimal header with the given `@HD` `SO`/`GO`/`SS` tags (and nothing else).
    fn header_with_hd(so: Option<&str>, go: Option<&str>, ss: Option<&str>) -> Header {
        let mut hd = Map::<HeaderRecord>::default();
        if let Some(so) = so {
            hd.other_fields_mut().insert(hd_tag::SORT_ORDER, BString::from(so));
        }
        if let Some(go) = go {
            hd.other_fields_mut().insert(hd_tag::GROUP_ORDER, BString::from(go));
        }
        if let Some(ss) = ss {
            hd.other_fields_mut().insert(hd_tag::SUBSORT_ORDER, BString::from(ss));
        }
        Header::builder().set_header(hd).build()
    }

    /// Every `@HD` `SS` spelling `detect_sort_order` must accept for `SO:queryname`:
    /// `queryname:lexicographical` (R2-HDR-01 corollary/Fix 1b — the SAM-spec/samtools
    /// spelling `detect_sort_order` must still accept, not `bail!`, once the writer emits
    /// it), `queryname:lexicographic` (back-compat: fgumi's own pre-fix bare spelling,
    /// R2-HDR-01), `queryname:natural` (fgbio/samtools/fgumi spelling), and no `SS` tag at
    /// all (defaults to lexicographic, matching `fgumi sort`'s CLI default).
    #[rstest]
    #[case::prefixed_lexicographical(
        Some("queryname:lexicographical"),
        SortOrder::Queryname(QuerynameComparator::Lexicographic)
    )]
    #[case::prefixed_lexicographic(
        Some("queryname:lexicographic"),
        SortOrder::Queryname(QuerynameComparator::Lexicographic)
    )]
    #[case::prefixed_natural(
        Some("queryname:natural"),
        SortOrder::Queryname(QuerynameComparator::Natural)
    )]
    #[case::bare_natural(Some("natural"), SortOrder::Queryname(QuerynameComparator::Natural))]
    #[case::bare_lexicographical(
        Some("lexicographical"),
        SortOrder::Queryname(QuerynameComparator::Lexicographic)
    )]
    #[case::defaults_to_lexicographic_when_ss_absent(
        None,
        SortOrder::Queryname(QuerynameComparator::Lexicographic)
    )]
    fn detect_sort_order_accepts_queryname_variants(
        #[case] ss: Option<&str>,
        #[case] expected: SortOrder,
    ) {
        let header = header_with_hd(Some("queryname"), None, ss);
        let order = detect_sort_order(&header).expect("SO:queryname variant must be accepted");
        assert_eq!(order, expected);
    }

    /// An unrecognized `SS` suffix for `SO:queryname` must still error, not silently fall
    /// back to lexicographic.
    #[test]
    fn detect_sort_order_rejects_unknown_queryname_subsort() {
        let header = header_with_hd(Some("queryname"), None, Some("queryname:bogus"));
        let err = detect_sort_order(&header).expect_err("unrecognized SS suffix must error");
        assert!(err.to_string().contains("bogus"), "error should name the bad suffix: {err}");
    }

    /// An `SS` whose sort-order prefix disagrees with `SO:queryname` (e.g.
    /// `coordinate:natural`) must be rejected — the previous suffix-only match accepted
    /// it and then verified the file with the `natural` comparator, which could return a
    /// spurious MATCH. The prefix, not just the suffix, is validated.
    #[rstest]
    #[case::coordinate_prefixed_natural("coordinate:natural")]
    #[case::unsorted_prefixed_lexicographical("unsorted:lexicographical")]
    fn detect_sort_order_rejects_ss_prefix_disagreeing_with_so(#[case] ss: &str) {
        let header = header_with_hd(Some("queryname"), None, Some(ss));
        let err = detect_sort_order(&header)
            .expect_err("SS with a sort-order prefix disagreeing with SO must error");
        assert!(
            err.to_string().contains("disagrees with SO:queryname"),
            "error should call out the prefix disagreement: {err}"
        );
    }

    /// Plain `SO:coordinate` with no `SS` tag must still be accepted as `SortOrder::Coordinate`
    /// (this engine implements no coordinate sub-sort, so the common case is bare `SO` alone).
    #[test]
    fn detect_sort_order_accepts_bare_coordinate() {
        let header = header_with_hd(Some("coordinate"), None, None);
        let order = detect_sort_order(&header).expect("bare SO:coordinate must be accepted");
        assert_eq!(order, SortOrder::Coordinate);
    }

    /// `SO:coordinate` with ANY `SS` sub-sort must be rejected, not silently validated as plain
    /// coordinate: this engine implements no coordinate sub-sort, so a declared
    /// `coordinate:<something>` sub-order is stronger than what `run_compare`'s
    /// equal-core-sort-key grouping actually verifies. Silently accepting it would let a file
    /// that only satisfies plain coordinate order falsely validate against a sub-sort it never
    /// promised.
    #[rstest]
    #[case::bare_unrecognized_subsort("bogus")]
    #[case::prefixed_unrecognized_subsort("coordinate:bogus")]
    fn detect_sort_order_rejects_coordinate_subsort(#[case] ss: &str) {
        let header = header_with_hd(Some("coordinate"), None, Some(ss));
        let err = detect_sort_order(&header)
            .expect_err("SO:coordinate with any SS sub-sort must be rejected");
        assert!(
            err.to_string().contains("coordinate"),
            "error should mention SO:coordinate: {err}"
        );
    }

    /// An `SS` whose sort-order prefix disagrees with `SO:coordinate` (e.g.
    /// `queryname:natural`) must be rejected with a message calling out the disagreement,
    /// mirroring the `SO:queryname` prefix-disagreement check above.
    #[test]
    fn detect_sort_order_rejects_coordinate_ss_prefix_disagreement() {
        let header = header_with_hd(Some("coordinate"), None, Some("queryname:natural"));
        let err = detect_sort_order(&header)
            .expect_err("SS with a sort-order prefix disagreeing with SO must error");
        assert!(
            err.to_string().contains("disagrees with SO:coordinate"),
            "error should call out the prefix disagreement: {err}"
        );
    }

    /// A record with a distinct read name, so each one is its own canonical content key.
    fn named(name: &[u8]) -> RawRecord {
        SamBuilder::new().read_name(name).flags(flags::FIRST_SEGMENT).build()
    }

    /// The memory contract for the whole engine: comparing an equal-core-key run must cost
    /// memory proportional to how far the two files' *orderings* diverge, not to the run's
    /// length. Every `tid = -1` record packs to one constant coordinate key
    /// (`extract_coordinate_key_inline`), so a BAM's unmapped tail is a single run whose
    /// length is the file's — buffering it whole is what made a 1M-read unmapped tail cost
    /// ~1 GB. Two files listing the same records in the same order cancel on arrival, so
    /// pending never exceeds one record per side regardless of run length.
    #[test]
    fn identical_run_cancels_on_arrival_and_keeps_pending_bounded() {
        let mut canceller = RunCanceller::new(EXACT_WINDOW_RECORDS);
        for i in 0..10_000u32 {
            let rec = named(format!("read{i}").as_bytes());
            canceller.observe(Side::Bam1, rec.as_ref());
            canceller.observe(Side::Bam2, rec.as_ref());
        }
        assert!(
            canceller.peak_pending() <= 1,
            "in-order identical runs must cancel on arrival, peaked at {}",
            canceller.peak_pending()
        );
        let residual = canceller.finish();
        assert!(residual.is_match(), "identical runs must leave no residual: {residual:?}");
        assert!(
            !residual.compared_by_digest(),
            "an in-order run never approaches the window, so it must stay exact: {residual:?}"
        );
    }

    /// Writes `records` to a temp BAM declaring `SO:coordinate`. All records are unmapped
    /// with no reference (`tid = -1`), so they pack to one constant coordinate key and form
    /// a single equal-core-key run — the shape this engine must compare without buffering.
    fn unmapped_tail_bam(records: &[RawRecord]) -> tempfile::NamedTempFile {
        let tmp = tempfile::NamedTempFile::new().expect("create temp BAM");
        let header = header_with_hd(Some("coordinate"), None, None);
        let mut writer = noodles::bam::io::Writer::new(
            std::fs::File::create(tmp.path()).expect("create test BAM"),
        );
        writer.write_header(&header).expect("write test header");
        for record in records {
            let record_buf = fgumi_raw_bam::raw_record_to_record_buf(record, &header)
                .expect("raw_record_to_record_buf should succeed in test");
            writer.write_alignment_record(&header, &record_buf).expect("write test record");
        }
        writer.try_finish().expect("finish test BAM");
        tmp
    }

    /// A mismatch inside the unmapped tail must name the read that differs. The whole tail
    /// is one run, so reporting only its record counts ("multiset differs (N vs M)") points
    /// at every unmapped read in the file at once and is useless for diagnosis; the
    /// cancellation residual knows exactly which read went unmatched.
    #[test]
    fn unmapped_tail_mismatch_names_the_offending_read() {
        let shared: Vec<RawRecord> =
            (0..64u32).map(|i| named(format!("shared{i:03}").as_bytes())).collect();
        let mut with_extra = shared.clone();
        with_extra.push(named(b"lonely_read"));

        let bam1 = unmapped_tail_bam(&with_extra);
        let bam2 = unmapped_tail_bam(&shared);

        let outcome = sort_verify_compare(bam1.path(), bam2.path(), 20).expect("compare runs");
        assert!(!outcome.is_match(), "an extra read in bam1 must not compare equal");
        assert!(
            outcome.diff_details.iter().any(|d| d.contains("lonely_read")),
            "diagnostics must name the unmatched read, got: {:?}",
            outcome.diff_details
        );
    }

    /// Cancellation tracks *multiplicity*, not mere presence. Read names are not unique in a
    /// BAM (mates, secondary/supplementary records share one), so a set-based canceller
    /// would let a second copy match the first and report a dropped duplicate as equal.
    #[rstest]
    #[case::bam1_has_the_extra_copy(2, 1, 1, 0)]
    #[case::bam2_has_the_extra_copy(1, 2, 0, 1)]
    #[case::equal_multiplicity_cancels(3, 3, 0, 0)]
    fn duplicate_records_cancel_by_multiplicity(
        #[case] copies_in_bam1: usize,
        #[case] copies_in_bam2: usize,
        #[case] expected_only_in_bam1: usize,
        #[case] expected_only_in_bam2: usize,
    ) {
        let mut canceller = RunCanceller::new(EXACT_WINDOW_RECORDS);
        let rec = named(b"duplicated");
        for _ in 0..copies_in_bam1 {
            canceller.observe(Side::Bam1, rec.as_ref());
        }
        for _ in 0..copies_in_bam2 {
            canceller.observe(Side::Bam2, rec.as_ref());
        }
        let RunResidual::Exact { only_in_bam1, only_in_bam2 } = canceller.finish() else {
            panic!("a handful of records is far inside the window; the run must stay exact");
        };
        assert_eq!(only_in_bam1.len(), expected_only_in_bam1, "bam1 residual");
        assert_eq!(only_in_bam2.len(), expected_only_in_bam2, "bam2 residual");
    }

    /// Records the two files list in opposite order still cancel completely — the run's
    /// multiset is what matters, not its internal order (that is the whole reason this
    /// engine compares runs as multisets). Memory is what degrades, not correctness.
    #[test]
    fn reversed_run_cancels_completely() {
        let records: Vec<RawRecord> =
            (0..256u32).map(|i| named(format!("read{i:03}").as_bytes())).collect();
        let mut canceller = RunCanceller::new(EXACT_WINDOW_RECORDS);
        for rec in &records {
            canceller.observe(Side::Bam1, rec.as_ref());
        }
        for rec in records.iter().rev() {
            canceller.observe(Side::Bam2, rec.as_ref());
        }
        assert!(canceller.finish().is_match(), "a reversed run is still the same multiset");
    }

    /// Feed `bam1`'s records then `bam2`'s (the reversed-order worst case: nothing cancels
    /// until the second side arrives, so pending grows to the run's length), through a
    /// canceller whose window is small enough to force the digest fallback.
    fn compare_via_digest(bam1: &[RawRecord], bam2: &[RawRecord], window: usize) -> RunCanceller {
        let mut canceller = RunCanceller::new(window);
        for rec in bam1 {
            canceller.observe(Side::Bam1, rec.as_ref());
        }
        for rec in bam2 {
            canceller.observe(Side::Bam2, rec.as_ref());
        }
        canceller
    }

    /// Distinct records numbered `0..n`.
    fn distinct_records(n: u32) -> Vec<RawRecord> {
        (0..n).map(|i| named(format!("read{i:03}").as_bytes())).collect()
    }

    /// #699: a run whose two orderings diverge past the window must still be *compared*,
    /// not aborted. The abort it replaced was unfixable in the case that motivated it —
    /// the two inputs are the outputs of the two sorters under comparison, so "re-sort both
    /// into a common order" destroys the very property being verified.
    #[test]
    fn a_run_past_the_window_is_compared_by_digest_rather_than_aborting() {
        let records = distinct_records(256);
        let reversed: Vec<RawRecord> = records.iter().rev().cloned().collect();
        let residual = compare_via_digest(&records, &reversed, 8).finish();

        assert!(residual.compared_by_digest(), "a fully-displaced run must degrade: {residual:?}");
        assert!(residual.is_match(), "a reversed run is still the same multiset: {residual:?}");
    }

    /// Degrading must not blind the comparison: the digest still decides the question a
    /// mismatched run is asked. Covers a record missing from either side, and the equal
    /// multiset that must still match. Every case here is settled by record count alone;
    /// the content the digest is keyed on is pinned separately, below.
    #[rstest]
    #[case::record_missing_from_bam2(255, 256, false)]
    #[case::record_missing_from_bam1(256, 255, false)]
    #[case::equal_multisets(256, 256, true)]
    fn the_digest_decides_a_mismatched_run(
        #[case] records_in_bam1: u32,
        #[case] records_in_bam2: u32,
        #[case] expected_match: bool,
    ) {
        let bam1 = distinct_records(records_in_bam1);
        let bam2: Vec<RawRecord> = distinct_records(records_in_bam2).into_iter().rev().collect();
        let residual = compare_via_digest(&bam1, &bam2, 8).finish();

        assert!(residual.compared_by_digest(), "the window must have been exceeded");
        assert_eq!(residual.is_match(), expected_match, "digest verdict: {residual:?}");
    }

    /// Equal counts, equal read names, one record's *content* changed.
    ///
    /// The cases above are all decided by `MultisetDigest::count`, so none of them says
    /// what the digest is keyed *on* — a digest that hashed the read name would pass every
    /// one. Here the two sides hold the same 64 names, so only a digest over the canonical
    /// content key can call it a mismatch, and a name-keyed one would report a false
    /// `IDENTICAL`.
    #[test]
    fn the_digest_catches_a_content_change_under_an_unchanged_name() {
        let bam1 = distinct_records(64);
        let mut bam2: Vec<RawRecord> = bam1.iter().rev().cloned().collect();
        // `named` builds a FIRST_SEGMENT record, so rebuilding `read000` as LAST_SEGMENT
        // changes its bytes and nothing else — the name it is found under is unchanged.
        let last = bam2.len() - 1;
        bam2[last] = SamBuilder::new().read_name(b"read000").flags(flags::LAST_SEGMENT).build();

        let residual = compare_via_digest(&bam1, &bam2, 8).finish();

        assert!(residual.compared_by_digest(), "the window must have been exceeded");
        assert!(
            !residual.is_match(),
            "equal counts and equal names, one record's content changed: {residual:?}"
        );
    }

    /// Degrading a run that has already degraded is documented as a no-op, and must be one.
    ///
    /// `degrade_to_digest` installs an *empty* placeholder before it matches, so an
    /// already-degraded arm that let the accumulated digests drop would leave two empty
    /// digests behind — equal to each other, and therefore a `matched: true` for a run that
    /// does not match. `observe` only ever calls it from `Exact`, so nothing reaches the
    /// second call but a direct one.
    #[test]
    fn degrading_an_already_degraded_run_keeps_its_digests() {
        let bam1 = distinct_records(64);
        let bam2: Vec<RawRecord> = distinct_records(63).into_iter().rev().collect();
        let mut canceller = compare_via_digest(&bam1, &bam2, 8);

        canceller.degrade_to_digest();

        let residual = canceller.finish();
        assert!(residual.compared_by_digest(), "the window must have been exceeded");
        assert!(
            !residual.is_match(),
            "the record missing from bam2 must survive a redundant degrade: {residual:?}"
        );
    }

    /// The digest must distinguish `{A, A, B}` from `{B, C, C}`. Both have three records and
    /// the *same* XOR, because XOR is self-inverse and each duplicated pair annihilates —
    /// so a digest built on count and XOR alone would report these different runs as equal.
    /// `MultisetDigest::sum` is what separates them, and this test is why it is there.
    #[test]
    fn digest_distinguishes_xor_cancelling_multisets() {
        let (a, b, c) = (named(b"aaa"), named(b"bbb"), named(b"ccc"));
        let bam1 = [a.clone(), a, b.clone()];
        let bam2 = [b, c.clone(), c];

        let residual = compare_via_digest(&bam1, &bam2, 1).finish();
        assert!(residual.compared_by_digest(), "window of 1 must force the digest");
        assert!(
            !residual.is_match(),
            "{{A,A,B}} and {{B,C,C}} share a count and an XOR but are different multisets"
        );
    }

    /// The fallback is a change of representation mid-run, not a restart: pairs cancelled
    /// exactly before the window is exceeded must stay cancelled after it. Here a prefix
    /// arrives in lockstep and cancels on arrival, then the two sides diverge far enough to
    /// degrade — and the run as a whole is still the same multiset.
    #[test]
    fn cancellations_made_before_degrading_survive_it() {
        let mut canceller = RunCanceller::new(8);
        for rec in &distinct_records(32) {
            canceller.observe(Side::Bam1, rec.as_ref());
            canceller.observe(Side::Bam2, rec.as_ref());
        }
        assert!(canceller.peak_pending() <= 1, "the lockstep prefix must cancel on arrival");

        let tail: Vec<RawRecord> =
            (0..64u32).map(|i| named(format!("tail{i:03}").as_bytes())).collect();
        for rec in &tail {
            canceller.observe(Side::Bam1, rec.as_ref());
        }
        for rec in tail.iter().rev() {
            canceller.observe(Side::Bam2, rec.as_ref());
        }

        let residual = canceller.finish();
        assert!(residual.compared_by_digest(), "the tail must have exceeded the window");
        assert!(
            residual.is_match(),
            "records cancelled before the fallback must not resurface after it: {residual:?}"
        );
    }

    /// The bounded-memory contract, asserted directly rather than inferred from RSS: however
    /// far the two orderings diverge, pending never exceeds the window by more than the one
    /// record whose arrival trips it. This is the whole point of #699 — a 77.3M-record
    /// unmapped tail must not be materialized.
    #[test]
    fn pending_never_exceeds_the_window() {
        const WINDOW: usize = 16;
        let records = distinct_records(1024);
        let reversed: Vec<RawRecord> = records.iter().rev().cloned().collect();
        let canceller = compare_via_digest(&records, &reversed, WINDOW);
        assert!(
            canceller.peak_pending() <= WINDOW + 1,
            "pending peaked at {} for a window of {WINDOW}",
            canceller.peak_pending()
        );
    }

    /// A mapped record at a given coordinate, so distinct positions form distinct runs.
    fn mapped(name: &[u8], pos: i32) -> RawRecord {
        SamBuilder::new().read_name(name).ref_id(0).pos(pos).flags(0).build()
    }

    /// Writes `records` to a temp BAM declaring `SO:coordinate` with a single reference
    /// sequence, so mapped records (`ref_id = 0`) can be written and distinct positions
    /// form distinct equal-core-key runs.
    fn mapped_bam(records: &[RawRecord]) -> tempfile::NamedTempFile {
        use noodles::sam::header::record::value::map::ReferenceSequence;

        let tmp = tempfile::NamedTempFile::new().expect("create temp BAM");
        let mut hd = Map::<HeaderRecord>::default();
        hd.other_fields_mut().insert(hd_tag::SORT_ORDER, BString::from("coordinate"));
        let header = Header::builder()
            .set_header(hd)
            .add_reference_sequence(
                BString::from("chr1"),
                Map::<ReferenceSequence>::new(
                    std::num::NonZeroUsize::new(100_000).expect("nonzero length"),
                ),
            )
            .build();

        let mut writer = noodles::bam::io::Writer::new(
            std::fs::File::create(tmp.path()).expect("create test BAM"),
        );
        writer.write_header(&header).expect("write test header");
        for record in records {
            let record_buf = fgumi_raw_bam::raw_record_to_record_buf(record, &header)
                .expect("raw_record_to_record_buf should succeed in test");
            writer.write_alignment_record(&header, &record_buf).expect("write test record");
        }
        writer.try_finish().expect("finish test BAM");
        tmp
    }

    /// A record that moved to a different coordinate must be reported, not absorbed — the
    /// regression a sort oracle exists to catch. It holds for a reason worth stating:
    /// [`content_key_exact`] covers the core-field bytes carrying `ref_id`/`pos`/`flags`,
    /// exactly the fields the coordinate run key `(tid, pos, reverse)` is derived from, so a
    /// moved record has a different content key from its old self and cannot cancel against
    /// it. See the note on [`compare_run`] for why that makes cross-run cancellation
    /// impossible independently of run scoping.
    #[test]
    fn a_record_that_moved_to_another_coordinate_is_reported() {
        let bam1 = mapped_bam(&[mapped(b"stationary", 100), mapped(b"wanderer", 200)]);
        let bam2 = mapped_bam(&[mapped(b"stationary", 100), mapped(b"wanderer", 300)]);

        let outcome = sort_verify_compare(bam1.path(), bam2.path(), 20).expect("compare runs");
        assert!(
            !outcome.is_match(),
            "a record moved to another coordinate must DIFFER: {:?}",
            outcome.diff_details
        );
    }

    /// One file having an entire *trailing run* the other lacks — unequal run counts, not a
    /// within-run difference — must be reported as a presence mismatch that names which side
    /// is longer and by how many records. This is the desync path where one stream is
    /// exhausted while the other still holds a run to count and drain
    /// ([`count_rest_of_run`]); the extra run has two records at one coordinate so the count
    /// is exercised past a single element.
    #[rstest]
    #[case::bam1_has_the_extra_trailing_run(true, "bam1 has 2 more")]
    #[case::bam2_has_the_extra_trailing_run(false, "bam2 has 2 more")]
    fn a_trailing_run_present_in_only_one_file_is_reported(
        #[case] extra_in_bam1: bool,
        #[case] expected_detail: &str,
    ) {
        let shared = [mapped(b"stationary", 100)];
        let longer =
            [mapped(b"stationary", 100), mapped(b"trailing_a", 200), mapped(b"trailing_b", 200)];

        let (bam1, bam2) = if extra_in_bam1 {
            (mapped_bam(&longer), mapped_bam(&shared))
        } else {
            (mapped_bam(&shared), mapped_bam(&longer))
        };

        let outcome = sort_verify_compare(bam1.path(), bam2.path(), 20).expect("compare runs");
        assert!(
            !outcome.is_match(),
            "an extra trailing run in one file must DIFFER: {:?}",
            outcome.diff_details
        );
        assert!(
            outcome.diff_details.iter().any(|d| d.contains(expected_detail)),
            "diagnostics must name the longer side and its extra record count, got: {:?}",
            outcome.diff_details
        );
    }

    /// Feed both sides through a canceller with the given window, in the interleaved order a
    /// run comparison would: one record from each side per step.
    fn verdict_with_window(left: &[u32], right: &[u32], window: usize) -> bool {
        let mut canceller = RunCanceller::new(window);
        for i in 0..left.len().max(right.len()) {
            if let Some(v) = left.get(i) {
                canceller.observe(Side::Bam1, named(format!("read{v}").as_bytes()).as_ref());
            }
            if let Some(v) = right.get(i) {
                canceller.observe(Side::Bam2, named(format!("read{v}").as_bytes()).as_ref());
            }
        }
        canceller.finish().is_match()
    }

    proptest::proptest! {
        /// The residual is empty exactly when the two sides hold equal record multisets — the
        /// property the whole run comparison rests on. Checked against a straightforward
        /// sorted content-key reference over randomized inputs.
        #[test]
        fn residual_is_empty_iff_content_multisets_are_equal(
            left in proptest::collection::vec(0u32..6, 0..12),
            right in proptest::collection::vec(0u32..6, 0..12),
        ) {
            let (mut sorted_left, mut sorted_right) = (left.clone(), right.clone());
            sorted_left.sort_unstable();
            sorted_right.sort_unstable();
            proptest::prop_assert_eq!(
                verdict_with_window(&left, &right, EXACT_WINDOW_RECORDS),
                sorted_left == sorted_right
            );
        }

        /// The correctness claim for the whole fallback: the digest reaches the *same*
        /// verdict the exact canceller would, on arbitrary input. A window of 0 degrades on
        /// the first held record, so every non-trivial case here runs entirely through the
        /// digest; `usize::MAX` never degrades. Any disagreement — a digest collision, a
        /// record lost across the fold, a multiplicity mishandled — fails this.
        #[test]
        fn digest_and_exact_verdicts_agree(
            left in proptest::collection::vec(0u32..6, 0..12),
            right in proptest::collection::vec(0u32..6, 0..12),
        ) {
            proptest::prop_assert_eq!(
                verdict_with_window(&left, &right, 0),
                verdict_with_window(&left, &right, usize::MAX)
            );
        }
    }

    /// [`sort_verify_compare`] against a deliberately tiny exact window, so the digest
    /// fallback is reachable from a fixture of a few dozen records rather than a million.
    fn compare_with_window(
        bam1: &Path,
        bam2: &Path,
        exact_window: usize,
    ) -> Result<SortVerifyOutcome> {
        sort_verify_compare_opened_with_window(
            OpenedInput::open(bam1)?,
            OpenedInput::open(bam2)?,
            20,
            exact_window,
        )
    }

    /// The #699 scenario end to end, through the real engine rather than the canceller alone:
    /// two BAMs holding the same unmapped tail in opposite orders — one equal-key run, fully
    /// displaced — must compare `IDENTICAL`, and must report that the verdict came from the
    /// digest. Before this change the same inputs aborted with an error whose suggested
    /// remedy could not be applied.
    #[test]
    fn a_reordered_unmapped_tail_compares_identical_via_the_digest() {
        let records = distinct_records(64);
        let reversed: Vec<RawRecord> = records.iter().rev().cloned().collect();
        let bam1 = unmapped_tail_bam(&records);
        let bam2 = unmapped_tail_bam(&reversed);

        let outcome = compare_with_window(bam1.path(), bam2.path(), 4).expect("compare runs");
        assert!(
            outcome.is_match(),
            "a reordered unmapped tail is the same multiset: {:?}",
            outcome.diff_details
        );
        assert_eq!(outcome.runs_compared_by_digest, 1, "the tail run must be reported as degraded");
        assert_eq!(outcome.bam1_count, 64, "both files must still be counted in full");
        assert_eq!(outcome.bam2_count, 64);
    }

    /// A genuine difference inside a degraded run must still be caught, and the diagnostic
    /// must say why it names no reads instead of reporting a bare count that reads as though
    /// the reads were checked and found anonymous.
    ///
    /// It must also quote the window that actually fired. The window is a field, not the
    /// constant, so a diagnostic built from [`EXACT_WINDOW_RECORDS`] would name a threshold
    /// this comparison never reached.
    #[test]
    fn a_difference_inside_a_degraded_run_is_reported_without_names() {
        const WINDOW: usize = 4;
        let records = distinct_records(64);
        let mut reversed: Vec<RawRecord> = records.iter().rev().cloned().collect();
        reversed.pop();
        let bam1 = unmapped_tail_bam(&records);
        let bam2 = unmapped_tail_bam(&reversed);

        let outcome = compare_with_window(bam1.path(), bam2.path(), WINDOW).expect("compare runs");
        assert!(!outcome.is_match(), "a dropped record must DIFFER");
        assert_eq!(outcome.run_mismatches, 1, "the tail run must be the mismatch");
        assert!(
            outcome.diff_details.iter().any(|d| d.contains("order-insensitive digest")),
            "the diagnostic must explain why no reads are named, got: {:?}",
            outcome.diff_details
        );
        assert!(
            outcome.diff_details.iter().any(|d| d.contains(&format!("more than {WINDOW} records"))),
            "the diagnostic must quote the window in force, not {EXACT_WINDOW_RECORDS}, got: {:?}",
            outcome.diff_details
        );
    }

    /// Runs that stay inside the window keep naming reads, even when another run in the same
    /// file degrades. The fallback is per-run, so a huge unmapped tail must not cost the
    /// mapped runs their diagnostics.
    #[test]
    fn a_degraded_run_does_not_cost_other_runs_their_diagnostics() {
        // One mapped run with a named difference, plus a displaced unmapped tail. Unmapped
        // records sort after mapped ones, so the tail is the file's last run either way.
        let tail = distinct_records(64);
        let mut bam1_records = vec![mapped(b"stationary", 100), mapped(b"only_in_bam1", 100)];
        bam1_records.extend(tail.iter().cloned());
        let mut bam2_records = vec![mapped(b"stationary", 100)];
        bam2_records.extend(tail.iter().rev().cloned());

        let bam1 = mapped_bam(&bam1_records);
        let bam2 = mapped_bam(&bam2_records);

        let outcome = compare_with_window(bam1.path(), bam2.path(), 4).expect("compare runs");
        assert!(!outcome.is_match(), "the extra mapped record must DIFFER");
        assert_eq!(outcome.runs_compared_by_digest, 1, "only the tail run should degrade");
        assert!(
            outcome.diff_details.iter().any(|d| d.contains("only_in_bam1")),
            "the exact run must still name its unmatched read, got: {:?}",
            outcome.diff_details
        );
    }

    /// A header carrying a single reference sequence, so `ref_id = 0` records are
    /// writable and `extract_coordinate_key_inline` sees a nonzero reference count.
    fn single_reference_header(so: &str) -> Header {
        use noodles::sam::header::record::value::map::ReferenceSequence;

        let mut hd = Map::<HeaderRecord>::default();
        hd.other_fields_mut().insert(hd_tag::SORT_ORDER, BString::from(so));
        Header::builder()
            .set_header(hd)
            .add_reference_sequence(
                BString::from("chr1"),
                Map::<ReferenceSequence>::new(
                    std::num::NonZeroUsize::new(100_000).expect("nonzero length"),
                ),
            )
            .build()
    }

    /// [`OrderCheck`] must apply the comparator its [`SortOrder`] names, not merely
    /// *a* comparator. The two queryname variants differ in nothing else, so picking
    /// the wrong one silently accepts a file that violates the order it declares —
    /// and content mode then pairs records positionally on that false assurance.
    ///
    /// Each queryname case is chosen so the *other* queryname comparator would
    /// disagree with it: `r1 < r10 < r2` holds lexicographically and not naturally,
    /// `r1 < r2 < r10` naturally and not lexicographically. A case table that only
    /// listed sequences both comparators accept would pass with the arms swapped.
    ///
    /// Every order carries a rejecting case for the same reason. Equal keys never
    /// violate under any comparator, so an accept-only entry — which is all
    /// `TemplateCoordinate` had — passes whichever extractor its arm is wired to,
    /// and that arm's extractor is the most involved of the four. An unpaired
    /// mapped read keys on `(tid, unclipped_pos, ...)`, so descending positions
    /// give it something to reject.
    #[rstest]
    #[case::coordinate_ascending(SortOrder::Coordinate, &[("a", 100), ("b", 200)], 0)]
    #[case::coordinate_backwards(SortOrder::Coordinate, &[("a", 200), ("b", 100)], 1)]
    #[case::coordinate_counts_every_violation(
        SortOrder::Coordinate,
        &[("a", 200), ("b", 100), ("c", 400), ("d", 300)],
        2
    )]
    #[case::lexicographic_accepts_lexicographic_order(
        SortOrder::Queryname(QuerynameComparator::Lexicographic),
        &[("r1", 100), ("r10", 100), ("r2", 100)],
        0
    )]
    #[case::lexicographic_rejects_natural_order(
        SortOrder::Queryname(QuerynameComparator::Lexicographic),
        &[("r2", 100), ("r10", 100)],
        1
    )]
    #[case::natural_accepts_natural_order(
        SortOrder::Queryname(QuerynameComparator::Natural),
        &[("r1", 100), ("r2", 100), ("r10", 100)],
        0
    )]
    #[case::natural_rejects_lexicographic_order(
        SortOrder::Queryname(QuerynameComparator::Natural),
        &[("r10", 100), ("r2", 100)],
        1
    )]
    #[case::template_coordinate_accepts_equal_keys(
        SortOrder::TemplateCoordinate,
        &[("a", 100), ("b", 100)],
        0
    )]
    #[case::template_coordinate_rejects_descending_positions(
        SortOrder::TemplateCoordinate,
        &[("a", 200), ("b", 100)],
        1
    )]
    fn order_check_applies_the_comparator_its_sort_order_names(
        #[case] order: SortOrder,
        #[case] records: &[(&str, i32)],
        #[case] expected_violations: u64,
    ) {
        let header = single_reference_header("coordinate");
        let mut check = OrderCheck::new(&header, order);
        for (name, pos) in records {
            check.observe(mapped(name.as_bytes(), *pos).as_ref());
        }

        let result = check.into_result(Path::new("input.bam"));
        if expected_violations == 0 {
            result.unwrap_or_else(|e| panic!("{order:?} must accept {records:?}: {e}"));
        } else {
            let err = result.expect_err(&format!("{order:?} must reject {records:?}"));
            let msg = err.to_string();
            assert!(
                msg.contains(&format!("{expected_violations} record(s) violate")),
                "must report every violation, not stop at the first: {msg}"
            );
            assert!(msg.contains("input.bam"), "must name the offending input: {msg}");
        }
    }

    /// The first violation names the record that broke the order, not merely the
    /// count — the report is what tells a user *where* to look.
    #[test]
    fn order_check_names_the_first_violating_record() {
        let header = single_reference_header("coordinate");
        let mut check = OrderCheck::new(&header, SortOrder::Coordinate);
        for (name, pos) in [("first", 100), ("second", 300), ("backwards", 200), ("last", 150)] {
            check.observe(mapped(name.as_bytes(), pos).as_ref());
        }

        let msg = check
            .into_result(Path::new("input.bam"))
            .expect_err("a backwards record must be rejected")
            .to_string();
        assert!(msg.contains("record 3"), "must give the 1-based record position: {msg}");
        assert!(msg.contains("backwards"), "must name the first violating read: {msg}");
    }
}
