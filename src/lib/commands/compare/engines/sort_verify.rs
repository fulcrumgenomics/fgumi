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
use std::fs::File;
use std::path::Path;

use anyhow::{Context, Result, anyhow, bail};

use fgumi_bam_io::create_raw_bam_reader;
use fgumi_raw_bam::{RawRecord, RawRecordView};
use fgumi_sort::{
    LibraryLookup, QuerynameComparator, RawBamRecordReader, RawQuerynameKey, RawQuerynameLexKey,
    RawSortKey, SortContext, SortOrder, TemplateKey, cb_hasher, extract_coordinate_key_inline,
    extract_template_key_inline, verify_sort_order,
};
use noodles::sam::Header;

use ahash::AHashMap;

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

/// Open a fresh raw-byte record reader over `path`, positioned just past the header, via
/// the shared [`super::open_raw_bam_reader`] helper.
fn open_raw_reader(path: &Path) -> Result<RawBamRecordReader<File>> {
    super::open_raw_bam_reader(path)
        .with_context(|| format!("opening raw BAM reader for {}", path.display()))
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
    reader: RawBamRecordReader<File>,
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
        reader: RawBamRecordReader<File>,
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
/// lockstep and cancelling each record against its counterpart as it arrives (see
/// [`RunCanceller`]). Never materializes the run, so an unmapped tail spanning the whole
/// file costs memory proportional only to how far the two orderings diverge.
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
) -> Result<RunPairOutcome>
where
    K: Clone,
    ExtractKey: Fn(&[u8]) -> K,
    IsViolation: Fn(&K, &K) -> bool,
{
    let mut canceller = RunCanceller::new(MAX_PENDING_RUN_RECORDS);
    let (mut count1, mut count2) = (0u64, 0u64);

    loop {
        // Lockstep: taking one record from each side per iteration keeps pending at ~one
        // record when the two files agree on order, instead of draining one side first.
        let from1 = next_in_run(stream1, next1, key, same_core)?;
        let from2 = next_in_run(stream2, next2, key, same_core)?;
        if from1.is_none() && from2.is_none() {
            break;
        }
        if let Some(rec) = from1 {
            canceller.observe_bam1(rec.as_ref())?;
            count1 += 1;
        }
        if let Some(rec) = from2 {
            canceller.observe_bam2(rec.as_ref())?;
            count2 += 1;
        }
    }

    Ok(RunPairOutcome { residual: canceller.finish(), count1, count2 })
}

/// Cap on records held pending inside a single equal-core-key run comparison
/// ([`RunCanceller`]), counted across both sides.
///
/// Not a tuning knob for run *length*: a run of any length costs O(1) pending as long as the
/// two files list its records in the same order, because each record cancels against its
/// counterpart on arrival. Pending grows only with the *displacement* between the two
/// orderings, so this fires when the inputs disagree about intra-run order beyond this
/// window — e.g. two independently-produced unmapped tails that were never re-sorted into a
/// common order. Re-sorting both inputs (`samtools sort`/`fgbio SortBam`) makes the
/// displacement zero.
const MAX_PENDING_RUN_RECORDS: usize = 5_000_000;

/// One side's unmatched records for a single canonical content key, plus a read name for
/// diagnostics. `count` is the multiplicity (a BAM may legitimately hold several records
/// with one content key), and `name` is the read name of the first such record — kept once
/// per distinct key so a residual can be reported by name rather than by opaque key bytes.
#[derive(Debug)]
struct PendingEntry {
    count: usize,
    name: Vec<u8>,
}

/// Records left unmatched on each side once a run is fully consumed: the exact symmetric
/// difference of the two sides' record multisets, by read name.
#[derive(Debug, Default)]
struct RunResidual {
    only_in_bam1: Vec<Vec<u8>>,
    only_in_bam2: Vec<Vec<u8>>,
}

/// Read names listed per side in a residual diagnostic before truncating. A mismatched run
/// can have arbitrarily many unmatched records; naming a few identifies the divergence
/// without letting one run emit an unbounded message.
const RESIDUAL_NAMES_SHOWN: usize = 5;

impl RunResidual {
    /// `true` when the two sides' record multisets were equal (nothing left over).
    fn is_empty(&self) -> bool {
        self.only_in_bam1.is_empty() && self.only_in_bam2.is_empty()
    }

    /// Render the unmatched read names per side, e.g.
    /// `; only in bam1: lonely_read`. Empty when nothing is left over.
    fn describe(&self) -> String {
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
        format!("{}{}", side("bam1", &self.only_in_bam1), side("bam2", &self.only_in_bam2))
    }
}

/// Streaming multiset comparison of one equal-core-key run, in memory proportional to how
/// far the two files' orderings diverge rather than to the run's length.
///
/// Each arriving record is reduced to its canonical content key ([`content_key_exact`],
/// whose byte-equality is *exactly* `Exact` content-equality) and cancelled against the
/// opposite side's pending set if a counterpart is waiting; otherwise it joins its own
/// side's pending set. Two files listing a run's records in the same order therefore never
/// hold more than one record pending, however long the run — which is what keeps a BAM's
/// unmapped tail (a single run spanning the whole file, since every `tid = -1` record packs
/// to one constant coordinate key) from being materialized in full.
///
/// Whatever remains at [`finish`](Self::finish) is the exact symmetric difference, so a
/// mismatch is reported as the specific reads involved instead of only a record count.
struct RunCanceller {
    pending1: AHashMap<Vec<u8>, PendingEntry>,
    pending2: AHashMap<Vec<u8>, PendingEntry>,
    pending1_len: usize,
    pending2_len: usize,
    peak_pending: usize,
    max_pending: usize,
}

impl RunCanceller {
    fn new(max_pending: usize) -> Self {
        Self {
            pending1: AHashMap::new(),
            pending2: AHashMap::new(),
            pending1_len: 0,
            pending2_len: 0,
            peak_pending: 0,
            max_pending,
        }
    }

    /// Observe a record from bam1: cancel it against a waiting bam2 record, else hold it.
    fn observe_bam1(&mut self, record: &[u8]) -> Result<()> {
        Self::observe(
            record,
            &mut self.pending2,
            &mut self.pending2_len,
            &mut self.pending1,
            &mut self.pending1_len,
        );
        self.after_observe()
    }

    /// Observe a record from bam2: cancel it against a waiting bam1 record, else hold it.
    fn observe_bam2(&mut self, record: &[u8]) -> Result<()> {
        Self::observe(
            record,
            &mut self.pending1,
            &mut self.pending1_len,
            &mut self.pending2,
            &mut self.pending2_len,
        );
        self.after_observe()
    }

    /// Cancel `record` against `opposite` if a counterpart waits there, otherwise add it to
    /// `own`. Side-agnostic so both directions share one implementation.
    fn observe(
        record: &[u8],
        opposite: &mut AHashMap<Vec<u8>, PendingEntry>,
        opposite_len: &mut usize,
        own: &mut AHashMap<Vec<u8>, PendingEntry>,
        own_len: &mut usize,
    ) {
        let key = content_key_exact(record);
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

    /// Refresh the peak watermark and enforce the pending cap.
    fn after_observe(&mut self) -> Result<()> {
        let pending = self.pending1_len + self.pending2_len;
        self.peak_pending = self.peak_pending.max(pending);
        if pending > self.max_pending {
            bail!(
                "sort-verify: more than {} records held pending while comparing one \
                 equal-sort-key run — the two inputs disagree about the order of records \
                 within that run by more than this window (an unmapped tail is one such \
                 run). Re-sort both inputs into a common order and retry.",
                self.max_pending
            );
        }
        Ok(())
    }

    /// Highest number of records held pending at any point, across both sides. Exists so
    /// tests can assert the bounded-memory contract directly (an in-order run must peak at
    /// ~1 record however long it is) rather than inferring it from process RSS.
    #[cfg(test)]
    fn peak_pending(&self) -> usize {
        self.peak_pending
    }

    /// Consume the canceller, reporting the records that never found a counterpart.
    fn finish(self) -> RunResidual {
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
        RunResidual { only_in_bam1: names(self.pending1), only_in_bam2: names(self.pending2) }
    }
}

/// Result of the run-grouped multiset comparison pass (before the per-file sort-order
/// violation counts are merged in by [`run_full_verify`]).
#[derive(Debug, Default)]
struct RunCompareOutcome {
    run_mismatches: u64,
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

                let pair = compare_run(stream1, &mut next1, stream2, &mut next2, &key1, same_core)?;
                if !pair.residual.is_empty() {
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
struct VerifyContext<'a, K, ExtractKey, IsViolation, SameCore>
where
    ExtractKey: Fn(&[u8]) -> K,
    IsViolation: Fn(&K, &K) -> bool,
    SameCore: Fn(&K, &K) -> bool,
{
    bam1: &'a Path,
    bam2: &'a Path,
    extract_key: ExtractKey,
    is_violation: IsViolation,
    same_core: SameCore,
    max_diffs: usize,
    order: SortOrder,
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
    ctx: VerifyContext<'_, K, ExtractKey, IsViolation, SameCore>,
) -> Result<SortVerifyOutcome>
where
    K: Clone,
    ExtractKey: Fn(&[u8]) -> K,
    IsViolation: Fn(&K, &K) -> bool,
    SameCore: Fn(&K, &K) -> bool,
{
    let VerifyContext { bam1, bam2, extract_key, is_violation, same_core, max_diffs, order } = ctx;

    let mut stream1 = OrderChecked::new(open_raw_reader(bam1)?, &extract_key, &is_violation);
    let mut stream2 = OrderChecked::new(open_raw_reader(bam2)?, &extract_key, &is_violation);

    let runs = run_compare(&mut stream1, &mut stream2, &same_core, max_diffs)?;

    Ok(SortVerifyOutcome {
        sort_order: order,
        bam1_count: stream1.tracker.total,
        bam2_count: stream2.tracker.total,
        bam1_violations: stream1.tracker.violations,
        bam1_first_violation: stream1.tracker.first_violation,
        bam2_violations: stream2.tracker.violations,
        bam2_first_violation: stream2.tracker.first_violation,
        run_mismatches: runs.run_mismatches,
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
    let (_, header1) = create_raw_bam_reader(bam1, 1)?;
    let (_, header2) = create_raw_bam_reader(bam2, 1)?;

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
    // consistent with the matching arms in verify_records_in_order. Future changes
    // to how an arm selects its key extractor or comparator must be applied to both sites.
    let mut outcome = match order {
        SortOrder::Coordinate => {
            let nref = header1.reference_sequences().len() as u32;
            run_full_verify(VerifyContext {
                bam1,
                bam2,
                extract_key: |bam: &[u8]| extract_coordinate_key_inline(bam, nref),
                is_violation: |key: &u64, prev: &u64| key < prev,
                same_core: |a: &u64, b: &u64| a == b,
                max_diffs,
                order,
            })
        }
        SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
            let ctx = SortContext::from_header(&header1);
            run_full_verify(VerifyContext {
                bam1,
                bam2,
                extract_key: move |bam: &[u8]| RawQuerynameLexKey::extract(bam, &ctx),
                is_violation: |key: &RawQuerynameLexKey, prev: &RawQuerynameLexKey| key < prev,
                same_core: |a: &RawQuerynameLexKey, b: &RawQuerynameLexKey| a == b,
                max_diffs,
                order,
            })
        }
        SortOrder::Queryname(QuerynameComparator::Natural) => {
            let ctx = SortContext::from_header(&header1);
            run_full_verify(VerifyContext {
                bam1,
                bam2,
                extract_key: move |bam: &[u8]| RawQuerynameKey::extract(bam, &ctx),
                is_violation: |key: &RawQuerynameKey, prev: &RawQuerynameKey| key < prev,
                same_core: |a: &RawQuerynameKey, b: &RawQuerynameKey| a == b,
                max_diffs,
                order,
            })
        }
        SortOrder::TemplateCoordinate => {
            let lib_lookup = LibraryLookup::from_header(&header1);
            let hasher = cb_hasher();
            // Matches `fgumi sort`'s own `--verify` (crate::commands::sort::parse_cell_tag):
            // template-coordinate always hashes the CB tag into the sort key when present.
            let cell_tag = Some(SamTag::CB);
            run_full_verify(VerifyContext {
                bam1,
                bam2,
                extract_key: move |bam: &[u8]| {
                    extract_template_key_inline(bam, &lib_lookup, cell_tag, &hasher)
                },
                is_violation: |key: &TemplateKey, prev: &TemplateKey| {
                    key.core_cmp(prev) == Ordering::Less
                },
                same_core: |a: &TemplateKey, b: &TemplateKey| a.core_cmp(b) == Ordering::Equal,
                max_diffs,
                order,
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

/// Verify every record in `path` is non-decreasing under `order`'s comparator. Err names the
/// first violating pair. Used as a universal precondition for order-dependent comparisons
/// (`CompareBams::execute`'s content-mode gate: an `@HD`-declared order must not be trusted
/// blindly — see `CompareMode::Content`'s doc comment).
///
/// This is a separate, fail-fast entry point from `sort_verify_compare`'s own per-file
/// violation *counting* (`run_full_verify`'s [`OrderChecked`]/[`OrderTracker`] fold, part of
/// the fuller `--command sort` report that also needs to keep going to compare run multisets
/// even after finding violations). Both call sites route through the same underlying
/// comparator machinery — the per-`SortOrder` key extractors
/// (`extract_coordinate_key_inline`, `RawQuerynameKey`/`RawQuerynameLexKey::extract`,
/// `extract_template_key_inline`) and their key comparisons (`<`, `TemplateKey::core_cmp`) —
/// so no comparator logic is reimplemented here, only the per-order dispatch. This entry
/// point calls [`verify_sort_order`] directly; the compare engine inlines that same
/// per-record fold ([`OrderTracker::observe`]) so it can verify order within its single
/// streaming pass.
///
/// # Errors
///
/// Returns an error if `path` cannot be opened/read, or if any record violates `order`
/// (naming the first violating record's 1-based position and read name).
pub(crate) fn verify_records_in_order(path: &Path, order: SortOrder) -> Result<()> {
    let (_, header) = create_raw_bam_reader(path, 1)
        .with_context(|| format!("opening {} to verify sort order", path.display()))?;

    // IMPORTANT: The per-SortOrder extractor/comparator selection below must stay
    // consistent with the matching arms in sort_verify_compare. Future changes
    // to how an arm selects its key extractor or comparator must be applied to both sites.
    let (_, violations, first_violation) = match order {
        SortOrder::Coordinate => {
            let nref = header.reference_sequences().len() as u32;
            let reader = open_raw_reader(path)?;
            verify_sort_order(
                reader,
                |bam: &[u8]| extract_coordinate_key_inline(bam, nref),
                |key: &u64, prev: &u64| key < prev,
            )?
        }
        SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
            let ctx = SortContext::from_header(&header);
            let reader = open_raw_reader(path)?;
            verify_sort_order(
                reader,
                move |bam: &[u8]| RawQuerynameLexKey::extract(bam, &ctx),
                |key: &RawQuerynameLexKey, prev: &RawQuerynameLexKey| key < prev,
            )?
        }
        SortOrder::Queryname(QuerynameComparator::Natural) => {
            let ctx = SortContext::from_header(&header);
            let reader = open_raw_reader(path)?;
            verify_sort_order(
                reader,
                move |bam: &[u8]| RawQuerynameKey::extract(bam, &ctx),
                |key: &RawQuerynameKey, prev: &RawQuerynameKey| key < prev,
            )?
        }
        SortOrder::TemplateCoordinate => {
            let lib_lookup = LibraryLookup::from_header(&header);
            let hasher = cb_hasher();
            // Matches `fgumi sort`'s own `--verify` (crate::commands::sort::parse_cell_tag):
            // template-coordinate always hashes the CB tag into the sort key when present.
            let cell_tag = Some(SamTag::CB);
            let reader = open_raw_reader(path)?;
            verify_sort_order(
                reader,
                move |bam: &[u8]| extract_template_key_inline(bam, &lib_lookup, cell_tag, &hasher),
                |key: &TemplateKey, prev: &TemplateKey| key.core_cmp(prev) == Ordering::Less,
            )?
        }
    };

    if violations > 0 {
        let (record_num, name) = first_violation
            .expect("verify_sort_order must report a first_violation when violations > 0");
        bail!(
            "{}: {violations} record(s) violate the declared {order:?} sort order (first at \
             record {record_num}, read name '{name}'); records must actually be in {order:?} \
             order, not merely declare it, for order-dependent comparison",
            path.display()
        );
    }

    Ok(())
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
        let mut canceller = RunCanceller::new(MAX_PENDING_RUN_RECORDS);
        for i in 0..10_000u32 {
            let rec = named(format!("read{i}").as_bytes());
            canceller.observe_bam1(rec.as_ref()).expect("bam1 observe within cap");
            canceller.observe_bam2(rec.as_ref()).expect("bam2 observe within cap");
        }
        assert!(
            canceller.peak_pending() <= 1,
            "in-order identical runs must cancel on arrival, peaked at {}",
            canceller.peak_pending()
        );
        let residual = canceller.finish();
        assert!(residual.is_empty(), "identical runs must leave no residual: {residual:?}");
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
        let mut canceller = RunCanceller::new(MAX_PENDING_RUN_RECORDS);
        let rec = named(b"duplicated");
        for _ in 0..copies_in_bam1 {
            canceller.observe_bam1(rec.as_ref()).expect("within cap");
        }
        for _ in 0..copies_in_bam2 {
            canceller.observe_bam2(rec.as_ref()).expect("within cap");
        }
        let residual = canceller.finish();
        assert_eq!(residual.only_in_bam1.len(), expected_only_in_bam1, "bam1 residual");
        assert_eq!(residual.only_in_bam2.len(), expected_only_in_bam2, "bam2 residual");
    }

    /// Records the two files list in opposite order still cancel completely — the run's
    /// multiset is what matters, not its internal order (that is the whole reason this
    /// engine compares runs as multisets). Memory is what degrades, not correctness.
    #[test]
    fn reversed_run_cancels_completely() {
        let records: Vec<RawRecord> =
            (0..256u32).map(|i| named(format!("read{i:03}").as_bytes())).collect();
        let mut canceller = RunCanceller::new(MAX_PENDING_RUN_RECORDS);
        for rec in &records {
            canceller.observe_bam1(rec.as_ref()).expect("within cap");
        }
        for rec in records.iter().rev() {
            canceller.observe_bam2(rec.as_ref()).expect("within cap");
        }
        assert!(canceller.finish().is_empty(), "a reversed run is still the same multiset");
    }

    /// The cap must fail with an actionable message rather than growing without bound. This
    /// is the reversed-order worst case: nothing cancels until the second side arrives, so
    /// pending grows to the run's length.
    #[test]
    fn exceeding_the_pending_cap_reports_an_actionable_error() {
        let mut canceller = RunCanceller::new(8);
        let err = (0..64u32)
            .find_map(|i| {
                let rec = named(format!("read{i:03}").as_bytes());
                canceller.observe_bam1(rec.as_ref()).err()
            })
            .expect("feeding only one side past the cap must error");
        let msg = err.to_string();
        assert!(msg.contains("held pending"), "error should explain the cap: {msg}");
        assert!(msg.contains("Re-sort both inputs"), "error should say how to fix it: {msg}");
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

    // The residual is empty exactly when the two sides hold equal record multisets — the
    // property the whole run comparison rests on. Checked against a straightforward sorted
    // content-key reference over randomized inputs.
    proptest::proptest! {
        #[test]
        fn residual_is_empty_iff_content_multisets_are_equal(
            left in proptest::collection::vec(0u32..6, 0..12),
            right in proptest::collection::vec(0u32..6, 0..12),
        ) {
            let mut canceller = RunCanceller::new(MAX_PENDING_RUN_RECORDS);
            for i in &left {
                canceller.observe_bam1(named(format!("read{i}").as_bytes()).as_ref()).unwrap();
            }
            for i in &right {
                canceller.observe_bam2(named(format!("read{i}").as_bytes()).as_ref()).unwrap();
            }
            let (mut sorted_left, mut sorted_right) = (left.clone(), right.clone());
            sorted_left.sort_unstable();
            sorted_right.sort_unstable();
            proptest::prop_assert_eq!(
                canceller.finish().is_empty(),
                sorted_left == sorted_right
            );
        }
    }
}
