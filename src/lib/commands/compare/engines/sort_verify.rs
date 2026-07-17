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
//! 2. Verifies each file independently is non-decreasing under that order's key, via
//!    [`fgumi_sort::verify_sort_order`] — this catches a genuine mis-sort in either file.
//! 3. Walks both files in file order, grouping each into maximal runs of records that
//!    share the same *core* sort key (the key without the name/hash tie-break lane — see
//!    [`fgumi_sort::TemplateKey::core_cmp`] for template-coordinate; coordinate has no
//!    tie-break lane beyond `(tid, pos, reverse)` in the first place), and asserts the
//!    multiset of full records is equal between the two files' corresponding runs. This
//!    tolerates intra-run reordering while still catching a missing/extra record or any
//!    content difference.
//!
//! Both checks run independently and are `AND`ed together in
//! [`SortVerifyOutcome::is_match`]: an order violation in either file, or any run
//! mismatch, is a `DIFFER`.

use std::cmp::Ordering;
use std::fs::File;
use std::path::Path;

use anyhow::{Context, Result, anyhow, bail};

use fgumi_bam_io::create_raw_bam_reader;
use fgumi_raw_bam::RawRecord;
use fgumi_sort::{
    LibraryLookup, QuerynameComparator, RawBamRecordReader, RawQuerynameKey, RawQuerynameLexKey,
    RawSortKey, SortContext, SortOrder, TemplateKey, VerifySummary, cb_hasher,
    extract_coordinate_key_inline, extract_template_key_inline, verify_sort_order,
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
    /// Number of sort-order violations found in `bam1` alone
    /// (via [`fgumi_sort::verify_sort_order`]).
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

/// Read one record and its key from `reader`, or `None` at EOF.
fn read_keyed<K>(
    reader: &mut RawBamRecordReader<File>,
    extract_key: &impl Fn(&[u8]) -> K,
) -> Result<Option<Keyed<K>>> {
    match reader.next() {
        Some(rec) => {
            let rec = rec?;
            let key = extract_key(rec.as_ref());
            Ok(Some((key, rec)))
        }
        None => Ok(None),
    }
}

/// Pull the next maximal run of records sharing the same core sort key from `reader`,
/// using `next` as one-record lookahead state (carried across calls). Returns `None`
/// once the stream (and lookahead) are exhausted.
fn pull_run<K>(
    reader: &mut RawBamRecordReader<File>,
    next: &mut Option<Keyed<K>>,
    extract_key: &impl Fn(&[u8]) -> K,
    same_core: &impl Fn(&K, &K) -> bool,
) -> Result<Option<(K, Vec<RawRecord>)>> {
    let Some((key, rec)) = next.take() else { return Ok(None) };
    let mut run = vec![rec];
    while let Some((k, r)) = read_keyed(reader, extract_key)? {
        if same_core(&k, &key) {
            run.push(r);
        } else {
            *next = Some((k, r));
            break;
        }
    }
    Ok(Some((key, run)))
}

/// Drain the remainder of `reader` (plus any pending lookahead), returning the number of
/// records consumed. Used once the two streams have desynchronized, so the final record
/// counts in [`SortVerifyOutcome`] still reflect the true totals.
fn drain_remaining<K>(
    reader: &mut RawBamRecordReader<File>,
    next: &mut Option<Keyed<K>>,
) -> Result<u64> {
    let mut count = u64::from(next.take().is_some());
    while reader.next().transpose()?.is_some() {
        count += 1;
    }
    Ok(count)
}

/// Returns `true` iff `a` and `b` contain the same multiset of records under
/// [`ContentPredicate::Exact`](super::content::ContentPredicate::Exact) (order-independent
/// tag comparison, all core fields exact).
///
/// Linear-time counting: each record is reduced to a canonical content key
/// ([`content_key_exact`]) whose byte-equality is *exactly* `Exact` content-equality (a
/// true equivalence relation — reflexive/symmetric/transitive, no epsilon tolerance), then
/// the two runs' key multisets are compared with a hash-count map. This replaces the
/// previous O(n²) greedy first-fit match, so a large tied-sort-key run — e.g. many records
/// piled on one genomic coordinate in high-depth data — no longer costs quadratic time.
///
/// Memory is O(run size): the run is already fully buffered by [`pull_run`], so counting
/// its keys adds no new unbounded buffering (and introduces no spill — see the design note
/// on why an external spill was judged unnecessary for this dev-only oracle path). The
/// canonical-key equivalence is locked to `content_diffs` by a proptest in
/// [`super::content`].
fn multiset_equal(a: &[RawRecord], b: &[RawRecord]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    // Net count per canonical content key: +1 per `a` record, -1 per `b` record. With equal
    // lengths, the two multisets are equal iff every net count ends at zero; a `b` key that
    // never appeared in `a` is an immediate mismatch.
    let mut counts: AHashMap<Vec<u8>, i64> = AHashMap::with_capacity(a.len());
    for rec in a {
        *counts.entry(content_key_exact(rec.as_ref())).or_insert(0) += 1;
    }
    for rec in b {
        match counts.get_mut(&content_key_exact(rec.as_ref())) {
            Some(count) => *count -= 1,
            None => return false,
        }
    }
    counts.values().all(|&count| count == 0)
}

/// Result of the run-grouped multiset comparison pass (before the per-file sort-order
/// violation counts are merged in by [`run_full_verify`]).
#[derive(Debug, Default)]
struct RunCompareOutcome {
    bam1_count: u64,
    bam2_count: u64,
    run_mismatches: u64,
    presence_mismatch: bool,
    diff_details: Vec<String>,
}

/// Walk `reader1`/`reader2` in file order, grouping each into maximal equal-core-key runs
/// and asserting the record multiset matches within each pair of corresponding runs (see
/// the module docs). Never resyncs: once a run boundary or stream-length mismatch is
/// found, no further runs are compared, though both readers are drained for accurate
/// counts.
fn run_compare<K>(
    mut reader1: RawBamRecordReader<File>,
    mut reader2: RawBamRecordReader<File>,
    extract_key: &impl Fn(&[u8]) -> K,
    same_core: &impl Fn(&K, &K) -> bool,
    max_diffs: usize,
) -> Result<RunCompareOutcome> {
    let mut out = RunCompareOutcome::default();

    let mut next1 = read_keyed(&mut reader1, extract_key)?;
    let mut next2 = read_keyed(&mut reader2, extract_key)?;

    let mut run_index = 0u64;

    loop {
        let run1 = pull_run(&mut reader1, &mut next1, extract_key, same_core)?;
        let run2 = pull_run(&mut reader2, &mut next2, extract_key, same_core)?;

        match (run1, run2) {
            (None, None) => break,
            (Some((_, recs1)), None) => {
                out.bam1_count += recs1.len() as u64;
                out.presence_mismatch = true;
                if out.diff_details.len() < max_diffs {
                    out.diff_details.push(format!(
                        "run {run_index}: bam1 has {} more record(s) than bam2 \
                         (bam2 exhausted first — no resync)",
                        recs1.len()
                    ));
                }
                out.bam1_count += drain_remaining(&mut reader1, &mut next1)?;
                break;
            }
            (None, Some((_, recs2))) => {
                out.bam2_count += recs2.len() as u64;
                out.presence_mismatch = true;
                if out.diff_details.len() < max_diffs {
                    out.diff_details.push(format!(
                        "run {run_index}: bam2 has {} more record(s) than bam1 \
                         (bam1 exhausted first — no resync)",
                        recs2.len()
                    ));
                }
                out.bam2_count += drain_remaining(&mut reader2, &mut next2)?;
                break;
            }
            (Some((k1, recs1)), Some((k2, recs2))) => {
                out.bam1_count += recs1.len() as u64;
                out.bam2_count += recs2.len() as u64;

                if !same_core(&k1, &k2) {
                    out.presence_mismatch = true;
                    if out.diff_details.len() < max_diffs {
                        out.diff_details.push(format!(
                            "run {run_index}: sort-key run boundary mismatch — bam1 and \
                             bam2 groups do not align (no resync)"
                        ));
                    }
                    out.bam1_count += drain_remaining(&mut reader1, &mut next1)?;
                    out.bam2_count += drain_remaining(&mut reader2, &mut next2)?;
                    break;
                }

                if !multiset_equal(&recs1, &recs2) {
                    out.run_mismatches += 1;
                    if out.diff_details.len() < max_diffs {
                        out.diff_details.push(format!(
                            "run {run_index}: record multiset differs ({} record(s) in \
                             bam1 vs {} in bam2)",
                            recs1.len(),
                            recs2.len()
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

/// Run both verification passes (per-file monotonic order, then run-grouped multiset
/// equality) for a single key type `K` and assemble the combined [`SortVerifyOutcome`].
fn run_full_verify<K, ExtractKey, IsViolation, SameCore>(
    ctx: VerifyContext<'_, K, ExtractKey, IsViolation, SameCore>,
) -> Result<SortVerifyOutcome>
where
    ExtractKey: Fn(&[u8]) -> K,
    IsViolation: Fn(&K, &K) -> bool,
    SameCore: Fn(&K, &K) -> bool,
{
    let VerifyContext { bam1, bam2, extract_key, is_violation, same_core, max_diffs, order } = ctx;

    let verify_one = |path: &Path| -> Result<VerifySummary> {
        let reader = open_raw_reader(path)?;
        verify_sort_order(reader, &extract_key, &is_violation)
    };

    let (_, bam1_violations, bam1_first_violation) = verify_one(bam1)?;
    let (_, bam2_violations, bam2_first_violation) = verify_one(bam2)?;

    let reader1 = open_raw_reader(bam1)?;
    let reader2 = open_raw_reader(bam2)?;
    let runs = run_compare(reader1, reader2, &extract_key, &same_core, max_diffs)?;

    Ok(SortVerifyOutcome {
        sort_order: order,
        bam1_count: runs.bam1_count,
        bam2_count: runs.bam2_count,
        bam1_violations,
        bam1_first_violation,
        bam2_violations,
        bam2_first_violation,
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
/// violation *counting* (`run_full_verify`'s `verify_one`, part of the fuller `--command
/// sort` report that also needs to keep going to compare run multisets even after finding
/// violations). Both call sites route through the same underlying comparator machinery —
/// the per-`SortOrder` key extractors (`extract_coordinate_key_inline`, `RawQuerynameKey`/
/// `RawQuerynameLexKey::extract`, `extract_template_key_inline`), their key comparisons
/// (`<`, `TemplateKey::core_cmp`), and [`verify_sort_order`] itself — so no comparator logic
/// is reimplemented here, only the per-order dispatch.
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
}
