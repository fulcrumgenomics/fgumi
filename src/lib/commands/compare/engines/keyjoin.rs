//! Key-join engine for `fgumi compare bams --command group`.
//!
//! `group` is the one BAM comparison that legitimately needs record order to differ
//! between the two tools under comparison (fgumi and fgbio number molecules
//! independently), so it is the sole exception to "`compare` never re-sorts a BAM".
//! This module holds the canonicalization half of that exception: given an
//! (arbitrarily ordered) BAM, produce a queryname-sorted copy via `fgumi_sort`'s
//! external (bounded-memory, spill-to-disk) sorter, so that a later merge-join can
//! walk both inputs in lockstep by read name. The merge-join itself (pairing records
//! on `RecordKey` and checking the fgumi-MI/fgbio-MI bijection) is added separately;
//! this module only provides the canonicalization primitive it depends on.
//!
//! # Why lexicographic, not natural, ordering
//!
//! The eventual merge-join only needs records with the *same name* to become
//! contiguous within each canonicalized stream — it never depends on the specific
//! string-ordering convention (natural vs. lexicographic), as long as **both** inputs
//! are canonicalized with the same comparator (true here: one call site, one
//! comparator, applied uniformly). `QuerynameComparator::Lexicographic` is
//! substantially faster than the natural-order comparator and carries no
//! samtools-compatibility requirement here (unlike `fgumi sort`'s user-facing
//! `--order queryname`), so it is used unconditionally.
//!
//! # Disk-backed scratch
//!
//! Per this repo's spill-validation guidance, the canonicalization sort's temp chunk
//! files must land on a disk-backed path — never a bare system temp directory, which
//! is tmpfs (RAM-backed) on many modern Linux distributions (Amazon Linux 2023,
//! RHEL9+, systemd-managed Ubuntu). [`resolve_sort_tmp_dirs`] resolves the CLI/env
//! precedence the same way `fgumi sort` does, but falls back to `/var/tmp` rather
//! than deferring to [`RawExternalSorter`]'s own system-temp default.

use std::cmp::Ordering;
use std::fs::File;
use std::path::{Path, PathBuf};

use ahash::AHashMap;
use anyhow::{Context, Result, ensure};
use fgumi_bam_io::create_raw_bam_reader;
use fgumi_raw_bam::RawRecord;
use fgumi_sort::{
    QuerynameComparator, RawBamRecordReader, RawExternalSorter, SortOrder, SortStats,
};
use noodles::sam::Header;

use super::super::bams::{MiKey, get_mi_tag_raw};
use super::super::record_key::{self, RecordKey};
use super::content::{ContentPredicate, content_diffs};
use super::header::{compare_headers, fold_header_diffs};
use super::push_diff;

/// Fallback scratch directory used when neither `--sort-tmp-dir` nor `FGUMI_TMP_DIRS`
/// resolve to anything.
///
/// On Unix this is `/var/tmp`: FHS-defined persistent temp storage that, unlike `/tmp` on
/// modern distros (Amazon Linux 2023, RHEL9+, systemd-managed Ubuntu), is not tmpfs by
/// default. Defaulting here — rather than letting [`RawExternalSorter`] fall through to its
/// own bare system-temp default — avoids silently spilling a large BAM sort into RAM (see
/// this repo's `CLAUDE.md` "Benchmark & Spill Validation" guidance).
///
/// On non-Unix platforms (e.g. Windows, where `/var/tmp` does not exist and `--command
/// group` would otherwise fail before comparing) this falls back to the platform temp
/// directory via [`std::env::temp_dir`], which is not RAM-backed there.
#[cfg(unix)]
fn default_disk_backed_tmp_dir() -> PathBuf {
    PathBuf::from("/var/tmp")
}

/// See [`default_disk_backed_tmp_dir`]'s Unix variant; this is the non-Unix fallback.
#[cfg(not(unix))]
fn default_disk_backed_tmp_dir() -> PathBuf {
    std::env::temp_dir()
}

/// Configuration for the key-join engine's internal canonicalization sort.
///
/// Threaded from `CompareBams`'s `--sort-memory`/`--sort-tmp-dir`/`--threads` flags.
pub struct KeyJoinConfig {
    /// Threads for the internal canonicalization sort (`RawExternalSorter::threads`).
    pub threads: usize,
    /// Memory budget (bytes) for the internal canonicalization sort
    /// (`RawExternalSorter::memory_limit`); exceeding it spills to `sort_tmp_dirs`.
    pub sort_memory: usize,
    /// Spill directories for the internal canonicalization sort.
    ///
    /// Should never be empty in production use — see [`resolve_sort_tmp_dirs`], which
    /// callers should use to build this field — but an empty vector is accepted here
    /// and simply defers to [`RawExternalSorter`]'s own default, which is convenient
    /// for unit tests that don't care about tmpfs safety.
    pub sort_tmp_dirs: Vec<PathBuf>,
    /// Maximum number of entries collected in [`KeyJoinOutcome::diff_details`].
    pub max_diffs: usize,
}

/// Resolve the final list of temp directories for the canonicalization sort.
///
/// Precedence matches `fgumi sort`/[`crate::commands::sort::resolve_tmp_dirs`]: CLI
/// flags (if non-empty) win, then the `FGUMI_TMP_DIRS` environment variable. Unlike
/// `fgumi sort`'s own resolution, an *empty* result here is replaced with
/// [`default_disk_backed_tmp_dir`] rather than left empty: `group`'s key-join
/// canonicalizes and spills two full-size BAMs, so silently deferring to whatever
/// [`RawExternalSorter`] picks on its own (a bare system temp directory, i.e. `/tmp`,
/// which is tmpfs/RAM on many modern hosts) risks an OOM that looks like a
/// disk-space problem rather than a memory one.
///
/// Called from `CompareBams::keyjoin_config` (`bams.rs`) to build the
/// [`KeyJoinConfig::sort_tmp_dirs`] threaded into [`keyjoin_compare`] for `--command
/// group`/`--mode grouping --ignore-order`.
pub(crate) fn resolve_sort_tmp_dirs(cli: &[PathBuf], env_value: Option<&str>) -> Vec<PathBuf> {
    let resolved = crate::commands::sort::resolve_tmp_dirs(cli, env_value);
    if resolved.is_empty() { vec![default_disk_backed_tmp_dir()] } else { resolved }
}

/// Canonicalize `input` to queryname order at `output`, via `fgumi_sort`'s external
/// (bounded-memory, spill-to-disk) merge-sort.
///
/// Uses [`QuerynameComparator::Lexicographic`] (see this module's doc comment for why
/// the comparator choice is safe for the merge-join) and fast (level 1) output
/// compression, since `output` is an internal scratch copy that the merge-join reads
/// once and then discards — unlike the sort's own spill chunks, which already default
/// to fast zstd compression internally.
///
/// # Errors
///
/// Returns an error if `input` cannot be read or the sort fails.
pub fn canonicalize_to_queryname(
    input: &Path,
    output: &Path,
    cfg: &KeyJoinConfig,
) -> Result<SortStats> {
    let mut sorter =
        RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::Lexicographic))
            .memory_limit(cfg.sort_memory)
            .threads(cfg.threads.max(1))
            .output_compression(1);
    if !cfg.sort_tmp_dirs.is_empty() {
        sorter = sorter.temp_dirs(cfg.sort_tmp_dirs.clone());
    }
    sorter
        .sort(input, output)
        .with_context(|| format!("canonicalizing {} to queryname order", input.display()))
}

/// Incrementally tracks whether the fgumi-MI <-> fgbio-MI mapping observed across the
/// key-join's matched pairs is a consistent bijection. This is the sole grouping-mode
/// engine: it folds together both the full-`MiKey` bijection check and the `base()`-level
/// duplex strand-pairing check (the two grouping invariants), fed one matched pair at a
/// time instead of from two fully-materialized `RecordKey -> MiKey` maps. Memory here is
/// O(number of distinct MI values), not O(number of records) — the property this phase
/// needs to stay off the O(records) path.
///
/// Only matched pairs where *both* sides carry an MI tag should be fed via [`Self::observe`];
/// records missing an MI on either side are excluded (the merge-join accounts for those
/// separately via `KeyJoinOutcome::missing_mi_bam1`/`missing_mi_bam2`).
#[derive(Default)]
struct MiBijectionTracker {
    /// For each MI seen in bam1, the first bam2 MI paired with it, and whether a later
    /// pairing disagreed (a genuine grouping split at full-`MiKey` granularity — this
    /// keeps `/A`/`/B` distinct, so it catches a strand *mislabel*).
    mi1_to_mi2: AHashMap<MiKey, (MiKey, bool)>,
    /// The mirror of `mi1_to_mi2`, in the bam2 -> bam1 direction.
    mi2_to_mi1: AHashMap<MiKey, (MiKey, bool)>,
    /// As above, but keyed by `MiKey::base()` (strand-suffix stripped) — catches a duplex
    /// strand-*pairing* split that the full-`MiKey` maps are structurally blind to (see
    /// [`MiKey::base`](super::super::bams::MiKey::base)'s doc-comment for the `X/A`+`X/B`
    /// vs. `X/A`+`Y/A` example this exists to catch).
    base1_to_base2: AHashMap<i64, (i64, bool)>,
    /// The mirror of `base1_to_base2`, in the bam2 -> bam1 direction.
    base2_to_base1: AHashMap<i64, (i64, bool)>,
}

impl MiBijectionTracker {
    /// Record one matched pair's `(bam1 MI, bam2 MI)` observation.
    fn observe(&mut self, mi1: MiKey, mi2: MiKey) {
        Self::observe_one(&mut self.mi1_to_mi2, mi1, mi2);
        Self::observe_one(&mut self.mi2_to_mi1, mi2, mi1);
        Self::observe_one(&mut self.base1_to_base2, mi1.base(), mi2.base());
        Self::observe_one(&mut self.base2_to_base1, mi2.base(), mi1.base());
    }

    /// Record a single `from -> to` observation in `map`, flagging a mismatch the first
    /// time `from` is seen paired with a `to` different from the one it was first seen
    /// paired with.
    fn observe_one<K: Eq + std::hash::Hash, V: Copy + PartialEq>(
        map: &mut AHashMap<K, (V, bool)>,
        from: K,
        to: V,
    ) {
        map.entry(from)
            .and_modify(|(first, has_mismatch)| {
                if *first != to {
                    *has_mismatch = true;
                }
            })
            .or_insert((to, false));
    }

    /// Total number of MI groups (either direction, either granularity) whose reads
    /// mapped to more than one counterpart group. Non-zero means the grouping is not a
    /// consistent bijection.
    fn mismatches(&self) -> u64 {
        fn count<K, V>(m: &AHashMap<K, (V, bool)>) -> u64 {
            m.values().filter(|(_, mismatch)| *mismatch).count() as u64
        }
        count(&self.mi1_to_mi2)
            + count(&self.mi2_to_mi1)
            + count(&self.base1_to_base2)
            + count(&self.base2_to_base1)
    }
}

/// Outcome of a [`keyjoin_compare`] run.
#[derive(Debug, Default)]
pub struct KeyJoinOutcome {
    /// Total number of records read from the canonicalized `bam1`.
    pub bam1_count: u64,
    /// Total number of records read from the canonicalized `bam2`.
    pub bam2_count: u64,
    /// Number of key-matched (equal `RecordKey`) pairs.
    pub matched: u64,
    /// Number of records whose key was present only in `bam1`.
    pub only_in_bam1: u64,
    /// Number of records whose key was present only in `bam2`.
    pub only_in_bam2: u64,
    /// Number of matched pairs whose content differed under
    /// [`ContentPredicate::ExactMinusMi`].
    pub content_diffs: u64,
    /// Number of matched pairs where `bam1`'s record was missing an MI tag.
    pub missing_mi_bam1: u64,
    /// Number of matched pairs where `bam2`'s record was missing an MI tag.
    pub missing_mi_bam2: u64,
    /// Number of MI groups (either direction, either granularity) that are not a
    /// consistent bijection — see [`MiBijectionTracker`].
    pub mi_bijection_mismatches: u64,
    /// `true` if the two *original* inputs' (not the internal queryname-canonicalized
    /// copies') `@HD`/`@SQ`/`@RG` headers disagreed on a field
    /// [`compare_headers`](super::header::compare_headers) considers significant (`@PG`/`@CO`
    /// are normalized and never contribute here).
    pub header_mismatch: bool,
    /// Human-readable diff strings (header, presence, content, and bijection mismatches),
    /// capped at the caller-supplied `max_diffs`.
    pub diff_details: Vec<String>,
}

impl KeyJoinOutcome {
    /// Returns `true` iff the two streams matched: no presence-only records on either
    /// side, no content diffs among matched pairs, no MI missing on a matched pair, the MI
    /// mapping observed across all matched pairs is a consistent bijection, and no
    /// significant header divergence was found.
    #[must_use]
    pub fn is_match(&self) -> bool {
        self.only_in_bam1 == 0
            && self.only_in_bam2 == 0
            && self.content_diffs == 0
            && self.missing_mi_bam1 == 0
            && self.missing_mi_bam2 == 0
            && self.mi_bijection_mismatches == 0
            && !self.header_mismatch
    }
}

/// Open a canonicalized BAM for sequential raw-record pulling, via the shared
/// [`super::open_raw_bam_reader`] helper (header via [`create_raw_bam_reader`], records via
/// [`RawBamRecordReader`] after `skip_header()`).
fn open_sorted_cursor(path: &Path) -> Result<RawBamRecordReader<File>> {
    super::open_raw_bam_reader(path)
}

/// F1 soundness guard: hard-fail if `key` is less than the last key pulled from this
/// side (`prev`), i.e. if the canonicalized stream is not non-decreasing under
/// [`RecordKey`]'s `Ord`. See [`keyjoin_compare`]'s doc comment for why this must be a
/// hard error (a real merge-join desync risk) rather than a `debug_assert!` (a no-op in
/// release builds).
///
/// # Errors
///
/// Returns an error if `key < *prev`.
fn ensure_non_decreasing(prev: Option<&RecordKey>, key: &RecordKey, side: &str) -> Result<()> {
    ensure!(
        prev.is_none_or(|p| p <= key),
        "{side} canonicalized stream violates RecordKey ordering: {key:?} came after {prev:?}"
    );
    Ok(())
}

/// Pull the complete run of consecutive records sharing `run_key` from `reader`, starting
/// with the already-read `first` (whose key is `run_key`). Returns the run and the first
/// record *past* it — the look-ahead whose key is strictly greater, or `None` at EOF.
///
/// Every additional record pulled is checked for `RecordKey` monotonicity against `run_key`
/// (the F1 soundness guard — see [`keyjoin_compare`]); a stream that regresses below the run
/// key is a hard error, not a silent desync.
///
/// # Errors
///
/// Returns an error if a record cannot be read, or if the stream violates `RecordKey`
/// ordering.
fn collect_equal_key_run(
    reader: &mut RawBamRecordReader<File>,
    first: RawRecord,
    run_key: &RecordKey,
    side: &str,
) -> Result<(Vec<RawRecord>, Option<RawRecord>)> {
    let mut run = vec![first];
    loop {
        match reader.next_record()? {
            None => return Ok((run, None)),
            Some(record) => {
                let key = record_key::record_key(&record);
                ensure_non_decreasing(Some(run_key), &key, side)?;
                if &key == run_key {
                    run.push(record);
                } else {
                    return Ok((run, Some(record)));
                }
            }
        }
    }
}

/// Fold one key-matched pair's `(bam1 MI, bam2 MI)` into `outcome`/`bijection`, exactly as
/// the merge-join's equal-key arm did per pair: observe the bijection when both records carry
/// an MI, and bump the missing-MI counters otherwise.
fn observe_pair_mi(
    outcome: &mut KeyJoinOutcome,
    bijection: &mut MiBijectionTracker,
    a: &RawRecord,
    b: &RawRecord,
) {
    match (get_mi_tag_raw(a), get_mi_tag_raw(b)) {
        (Some(mi1), Some(mi2)) => bijection.observe(mi1, mi2),
        (None, None) => {
            outcome.missing_mi_bam1 += 1;
            outcome.missing_mi_bam2 += 1;
        }
        (None, Some(_)) => outcome.missing_mi_bam1 += 1,
        (Some(_), None) => outcome.missing_mi_bam2 += 1,
    }
}

/// Match one equal-`RecordKey` run's content as a multiset under
/// [`ContentPredicate::ExactMinusMi`], folding the result into `outcome`/`bijection`.
///
/// `RecordKey` equality does not imply content equality, and the two canonicalized streams
/// may order an equal-key run differently (there is no intra-run tie-break — see
/// [`keyjoin_compare`]'s duplicate-key note). Pairing such a run positionally would invent
/// content diffs on a mere reorder; instead this matches by content in three passes:
///
/// 1. **Content multiset match** — greedily pair each `run1` record with a still-unpaired
///    `run2` record whose content is identical under `ExactMinusMi`. Content equality (minus
///    MI) is an equivalence relation, so greedy matching within each equivalence class is
///    optimal. These are clean matches (no content diff); each feeds the MI bijection.
/// 2. **Leftover positional pairing** — records with no content-equal partner are paired
///    positionally (up to the smaller leftover count) and reported as content diffs, exactly
///    as the unique-key path reports a single differing pair (one matched pair, one content
///    diff, one MI observation).
/// 3. **Multiplicity excess** — any records still unpaired on the longer side have no
///    counterpart at all within the run; each is reported once as a content diff.
///
/// For a unique key (`run1.len() == run2.len() == 1` — the overwhelmingly common case) this
/// reduces exactly to the previous single-pair behaviour: one match, plus one content diff
/// iff the pair differs.
#[allow(clippy::too_many_arguments)]
fn compare_equal_key_run(
    run1: &[RawRecord],
    run2: &[RawRecord],
    run_key: &RecordKey,
    header1: &Header,
    header2: &Header,
    outcome: &mut KeyJoinOutcome,
    bijection: &mut MiBijectionTracker,
    cfg: &KeyJoinConfig,
) {
    let content_matches = |a: &RawRecord, b: &RawRecord| {
        content_diffs(a.as_ref(), b.as_ref(), ContentPredicate::ExactMinusMi, header1, header2)
            .is_none()
    };

    // Pass 1: content-multiset match (order-independent), so an intra-run reorder is not
    // misreported as a content diff. `paired2[j]` marks a `run2` record already claimed.
    let mut paired2 = vec![false; run2.len()];
    let mut unmatched1: Vec<usize> = Vec::new();
    for (i, rec1) in run1.iter().enumerate() {
        match (0..run2.len()).find(|&j| !paired2[j] && content_matches(rec1, &run2[j])) {
            Some(j) => {
                paired2[j] = true;
                outcome.matched += 1;
                observe_pair_mi(outcome, bijection, rec1, &run2[j]);
            }
            None => unmatched1.push(i),
        }
    }
    let unmatched2: Vec<usize> = (0..run2.len()).filter(|&j| !paired2[j]).collect();

    // Pass 2: pair leftover records positionally — key-matched, but differing in content
    // (the same one-matched-pair-plus-one-content-diff accounting the unique-key path uses).
    let leftover_pairs = unmatched1.len().min(unmatched2.len());
    for p in 0..leftover_pairs {
        let rec1 = &run1[unmatched1[p]];
        let rec2 = &run2[unmatched2[p]];
        outcome.matched += 1;
        outcome.content_diffs += 1;
        if let Some(diffs) = content_diffs(
            rec1.as_ref(),
            rec2.as_ref(),
            ContentPredicate::ExactMinusMi,
            header1,
            header2,
        ) {
            push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
                format!("record {run_key:?}: {}", diffs.join("; "))
            });
        }
        observe_pair_mi(outcome, bijection, rec1, rec2);
    }

    // Pass 3: multiplicity excess — records with no counterpart at all in the other run.
    for _ in &unmatched1[leftover_pairs..] {
        outcome.content_diffs += 1;
        push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
            format!(
                "record {run_key:?}: extra record in bam1 with no counterpart in bam2's \
                 equal-key run"
            )
        });
    }
    for _ in &unmatched2[leftover_pairs..] {
        outcome.content_diffs += 1;
        push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
            format!(
                "record {run_key:?}: extra record in bam2 with no counterpart in bam1's \
                 equal-key run"
            )
        });
    }
}

/// Canonicalize both inputs to queryname order, then merge-join on `RecordKey`, comparing
/// matched pairs under [`ContentPredicate::ExactMinusMi`] and tracking the fgumi-MI/fgbio-MI
/// bijection via [`MiBijectionTracker`]. Also compares the two *original* inputs' headers via
/// [`compare_headers`](super::header::compare_headers) (`@HD`/`@SQ`/`@RG`, normalizing away
/// `@PG`/`@CO`) — a significant divergence there is folded into
/// [`KeyJoinOutcome::header_mismatch`]/[`KeyJoinOutcome::is_match`] alongside the merge-join's
/// own findings.
///
/// # Merge-join semantics (no resync)
///
/// The two canonicalized streams are walked in lockstep, always advancing whichever
/// cursor currently holds the **lower** [`RecordKey`] (per its `Ord`). When the keys are
/// equal, the pair is a match: its content is compared and, if both sides carry an MI
/// tag, the pair is fed to the bijection tracker. When the keys differ, the lower-keyed
/// record is recorded as present-only-on-its-side and *only that side's* cursor
/// advances — the two cursors are never buffered past a mismatch or re-aligned by
/// search-ahead. This mirrors the positional engine's "no resync" discipline (see
/// `super::positional`); a key-join is not license to loosen that principle, since a
/// spurious pairing across an actual discrepancy would be worse than reporting the
/// discrepancy directly.
///
/// # Soundness guard (F1)
///
/// The merge-join assumes each canonicalized stream is non-decreasing under
/// [`RecordKey`]'s `Ord`. `group` output R1/R2 satisfies this, but a record with *both*
/// `FIRST_OF_PAIR` and `LAST_OF_PAIR` set is a real counter-example:
/// [`record_key::record_key`] maps that flag combination to its own
/// [`Segment::FirstAndLast`](super::super::record_key::Segment::FirstAndLast) variant,
/// which is deliberately ordered *before* `First`/`Last` (so it sorts *early* among a
/// name's records), while `fgumi_sort`'s `queryname_flag_order` bit-packing sorts that
/// same flag combination *last* — so the canonicalized physical order and `RecordKey::Ord`
/// disagree, and a stream containing such a record is not actually non-decreasing under
/// the key the join advances on. (This intentional ordering disagreement is documented on
/// the `Segment` enum itself.) Rather than a
/// `debug_assert!` (a no-op in release builds — exactly the silent-desync failure mode
/// this effort exists to close), this is a hard `anyhow::ensure!` on every record pulled
/// from either side: a violation aborts the compare instead of silently producing a
/// wrong verdict.
///
/// # Errors
///
/// Returns an error if either input cannot be canonicalized or read, or if either
/// canonicalized stream is found to violate `RecordKey` monotonicity (see above).
pub fn keyjoin_compare(bam1: &Path, bam2: &Path, cfg: &KeyJoinConfig) -> Result<KeyJoinOutcome> {
    // Compare the *original* inputs' headers, not the internal queryname-canonicalized
    // copies': `canonicalize_to_queryname` always rewrites `@HD` to a fixed queryname sort
    // order (see `create_output_header`), so comparing the canonicalized copies' `@HD` would
    // trivially always agree regardless of what the two original inputs actually declared.
    let (_, orig_header1) = create_raw_bam_reader(bam1, 1)?;
    let (_, orig_header2) = create_raw_bam_reader(bam2, 1)?;
    let header_diffs = compare_headers(&orig_header1, &orig_header2);

    // F2: never fall back to a bare `tempfile::tempdir()` (system temp, i.e. `/tmp`,
    // which is tmpfs/RAM on many modern hosts) — always resolve to a disk-backed
    // directory first, via the same `resolve_sort_tmp_dirs` this module already uses for
    // the canonicalization sort's own spill chunks (env_value is `None` here: `cfg`'s
    // `sort_tmp_dirs` has already folded in any `FGUMI_TMP_DIRS` resolution upstream).
    let scratch_dirs = resolve_sort_tmp_dirs(&cfg.sort_tmp_dirs, None);
    let scratch_base = scratch_dirs.first().expect("non-empty by construction above");
    let scratch = tempfile::Builder::new()
        .prefix("fgumi-compare-keyjoin-")
        .tempdir_in(scratch_base)
        .with_context(|| {
            format!("creating key-join scratch directory under {}", scratch_base.display())
        })?;

    let canon1 = scratch.path().join("bam1.queryname.bam");
    let canon2 = scratch.path().join("bam2.queryname.bam");
    canonicalize_to_queryname(bam1, &canon1, cfg)?;
    canonicalize_to_queryname(bam2, &canon2, cfg)?;

    let (_, header1) = create_raw_bam_reader(&canon1, 1)?;
    let (_, header2) = create_raw_bam_reader(&canon2, 1)?;
    let mut reader1 = open_sorted_cursor(&canon1)?;
    let mut reader2 = open_sorted_cursor(&canon2)?;

    let mut outcome = KeyJoinOutcome::default();
    fold_header_diffs(
        header_diffs,
        &mut outcome.header_mismatch,
        &mut outcome.diff_details,
        cfg.max_diffs,
    );
    let mut bijection = MiBijectionTracker::default();
    // F1 soundness guard state: the last key pulled from each side, so every
    // subsequently pulled key can be checked for non-decreasing order (see this
    // function's doc comment).
    let mut prev1: Option<RecordKey> = None;
    let mut prev2: Option<RecordKey> = None;

    let mut next1 = reader1.next_record()?;
    let mut next2 = reader2.next_record()?;

    loop {
        let (a, b) = match (next1.take(), next2.take()) {
            (None, None) => break,
            (Some(a), None) => {
                outcome.bam1_count += 1;
                outcome.only_in_bam1 += 1;
                let key = record_key::record_key(&a);
                ensure_non_decreasing(prev1.as_ref(), &key, "bam1")?;
                push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
                    format!("record present only in bam1: {key:?}")
                });
                prev1 = Some(key);
                next1 = reader1.next_record()?;
                continue;
            }
            (None, Some(b)) => {
                outcome.bam2_count += 1;
                outcome.only_in_bam2 += 1;
                let key = record_key::record_key(&b);
                ensure_non_decreasing(prev2.as_ref(), &key, "bam2")?;
                push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
                    format!("record present only in bam2: {key:?}")
                });
                prev2 = Some(key);
                next2 = reader2.next_record()?;
                continue;
            }
            (Some(a), Some(b)) => (a, b),
        };

        let ka = record_key::record_key(&a);
        let kb = record_key::record_key(&b);

        ensure_non_decreasing(prev1.as_ref(), &ka, "bam1")?;
        ensure_non_decreasing(prev2.as_ref(), &kb, "bam2")?;

        match ka.cmp(&kb) {
            Ordering::Less => {
                outcome.bam1_count += 1;
                outcome.only_in_bam1 += 1;
                push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
                    format!("record present only in bam1: {ka:?}")
                });
                prev1 = Some(ka);
                next1 = reader1.next_record()?;
                next2 = Some(b);
            }
            Ordering::Greater => {
                outcome.bam2_count += 1;
                outcome.only_in_bam2 += 1;
                push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
                    format!("record present only in bam2: {kb:?}")
                });
                prev2 = Some(kb);
                next2 = reader2.next_record()?;
                next1 = Some(a);
            }
            Ordering::Equal => {
                // `RecordKey` is collision-resistant, not unique: a run of several records
                // can share one key (e.g. same-locus secondary/supplementary alignments, or
                // malformed duplicate primaries). The queryname canonicalization has no
                // tie-break within such a run, so the two streams may order it differently.
                // Consume the *whole* equal-key run from both sides and match its content as
                // a multiset, so an intra-run reorder is not misreported as a content diff
                // (and a real difference is not masked by an arbitrary positional pairing).
                let run_key = ka; // == kb
                let (run1, look1) = collect_equal_key_run(&mut reader1, a, &run_key, "bam1")?;
                let (run2, look2) = collect_equal_key_run(&mut reader2, b, &run_key, "bam2")?;
                outcome.bam1_count += run1.len() as u64;
                outcome.bam2_count += run2.len() as u64;
                compare_equal_key_run(
                    &run1,
                    &run2,
                    &run_key,
                    &header1,
                    &header2,
                    &mut outcome,
                    &mut bijection,
                    cfg,
                );
                prev1 = Some(run_key.clone());
                prev2 = Some(run_key);
                next1 = look1;
                next2 = look2;
            }
        }
    }

    outcome.mi_bijection_mismatches = bijection.mismatches();
    if outcome.mi_bijection_mismatches > 0 {
        push_diff(&mut outcome.diff_details, cfg.max_diffs, || {
            format!(
                "MI grouping is not a consistent bijection ({} conflicting mapping(s) detected)",
                outcome.mi_bijection_mismatches
            )
        });
    }
    Ok(outcome)
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::BString;
    use fgumi_raw_bam::{RawRecord, SamBuilder, flags, raw_record_to_record_buf};
    use noodles::bam;
    use noodles::sam::Header;
    use noodles::sam::alignment::io::Write as AlignmentWrite;
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReferenceSequence;
    use rstest::rstest;
    use std::fs;
    use std::num::NonZeroUsize;
    use tempfile::tempdir;

    /// Builds a minimal single-`@SQ` header, just enough to write and read back a BAM.
    fn minimal_header() -> Header {
        let reference_sequence =
            Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).expect("non-zero length"));
        Header::builder().add_reference_sequence(BString::from("chr1"), reference_sequence).build()
    }

    /// Builds a content-distinct record (name, flags, and a unique SEQ so each fingerprint is
    /// unique) for the canonicalization test.
    fn canon_record(name: &[u8], flags: u16, seq: &[u8]) -> RawRecord {
        SamBuilder::new()
            .read_name(name)
            .flags(flags)
            .sequence(seq)
            .qualities(&vec![30u8; seq.len()])
            .build()
    }

    /// Writes `records` to a BAM file at `path` under `header`.
    fn write_bam(path: &Path, header: &Header, records: &[RawRecord]) {
        let mut writer = bam::io::Writer::new(fs::File::create(path).expect("create test BAM"));
        writer.write_header(header).expect("write test header");
        for record in records {
            let buf = raw_record_to_record_buf(record, header).expect("decode test record");
            writer.write_alignment_record(header, &buf).expect("write test record");
        }
        writer.try_finish().expect("finish test BAM");
    }

    /// Zeroes a record's non-semantic BAM `bin` field (bytes 10-11) so it does not enter the
    /// fingerprint: `bin` is an index optimization the sort pipeline writes as 0, not part of
    /// the SAM data model (the same field `raw_core_fields_equal` excludes from comparison).
    fn strip_bin(mut bytes: Vec<u8>) -> Vec<u8> {
        if bytes.len() >= 12 {
            bytes[10..12].fill(0);
        }
        bytes
    }

    /// The samtools queryname flag-ordering tie-break (`bam_sort.c`): among records sharing a
    /// name, the order is READ1, READ2, then PRIMARY < SUPPLEMENTARY < SECONDARY. Written from
    /// the samtools formula, independently of `fgumi_sort`'s own `queryname_flag_order`.
    fn samtools_queryname_flag_order(flags: u16) -> u16 {
        ((flags & 0xc0) << 8) | ((flags & 0x100) << 3) | ((flags & 0x800) >> 3)
    }

    /// Reads back every record from a BAM file, in on-disk order, as
    /// `(bin-stripped raw bytes, read name, flags)`.
    fn read_all_records(path: &Path) -> Vec<(Vec<u8>, Vec<u8>, u16)> {
        let file = fs::File::open(path).expect("open test BAM");
        let mut reader = fgumi_sort::RawBamRecordReader::new(file).expect("open raw reader");
        reader.skip_header().expect("skip header");
        let mut records = Vec::new();
        while let Some(record) = reader.next_record().expect("read record") {
            records.push((
                strip_bin(record.as_ref().to_vec()),
                record.read_name().to_vec(),
                record.flags(),
            ));
        }
        records
    }

    /// Canonicalization must (a) preserve the exact *records* — every field, compared by an
    /// independent full-byte fingerprint, not just the read name — and (b) leave the output in
    /// non-decreasing order under the FULL queryname sort key (name plus the
    /// segment/secondary/supplementary flag tie lanes), whether or not the memory budget forces
    /// the sort to spill to disk.
    #[rstest]
    #[case::ample_memory(1024 * 1024, false)]
    #[case::tiny_memory_forces_spill(1, true)]
    fn canonicalize_sorts_by_queryname_and_preserves_records(
        #[case] sort_memory: usize,
        #[case] expect_spill: bool,
    ) {
        let tmp = tempdir().expect("tempdir");
        let header = minimal_header();
        // Scrambled input; `dup` repeats across four flag lanes (R1 primary, R1 supplementary,
        // R1 secondary, R2 primary) to exercise the tie-break, and every record has a distinct
        // SEQ so the fingerprint catches a payload landing on the wrong name.
        let records = vec![
            canon_record(b"readZ", flags::PAIRED | flags::FIRST_SEGMENT, b"AAAAAAAA"),
            canon_record(b"dup", flags::PAIRED | flags::LAST_SEGMENT, b"CCCCCCCC"),
            canon_record(
                b"dup",
                flags::PAIRED | flags::FIRST_SEGMENT | flags::SUPPLEMENTARY,
                b"GGGGGGGG",
            ),
            canon_record(b"aaa", 0, b"TTTTTTTT"),
            canon_record(b"dup", flags::PAIRED | flags::FIRST_SEGMENT, b"ACGTACGT"),
            canon_record(
                b"dup",
                flags::PAIRED | flags::FIRST_SEGMENT | flags::SECONDARY,
                b"TGCATGCA",
            ),
        ];
        let input = tmp.path().join("in.bam");
        write_bam(&input, &header, &records);

        let output = tmp.path().join("out.bam");
        let cfg = KeyJoinConfig { threads: 1, sort_memory, sort_tmp_dirs: vec![], max_diffs: 10 };
        let stats =
            canonicalize_to_queryname(&input, &output, &cfg).expect("canonicalize should succeed");

        assert_eq!(stats.total_records, records.len() as u64);
        assert_eq!(stats.output_records, records.len() as u64);
        assert_eq!(
            stats.chunks_written > 0,
            expect_spill,
            "chunks_written={} for sort_memory={sort_memory}",
            stats.chunks_written
        );

        let got = read_all_records(&output);

        // (b) Non-decreasing under the full sort key (name + flag tie-break).
        for pair in got.windows(2) {
            let key =
                |r: &(Vec<u8>, Vec<u8>, u16)| (r.1.clone(), samtools_queryname_flag_order(r.2));
            assert!(
                key(&pair[0]) <= key(&pair[1]),
                "canonicalized output must be non-decreasing by (name, flag-order): \
                 {:?}/{:#06x} then {:?}/{:#06x}",
                String::from_utf8_lossy(&pair[0].1),
                pair[0].2,
                String::from_utf8_lossy(&pair[1].1),
                pair[1].2,
            );
        }

        // (a) Output is an exact permutation of the input records, compared by full bytes.
        let mut want: Vec<Vec<u8>> =
            records.iter().map(|r| strip_bin(r.as_ref().to_vec())).collect();
        want.sort();
        let mut got_bytes: Vec<Vec<u8>> = got.into_iter().map(|(bytes, _, _)| bytes).collect();
        got_bytes.sort();
        assert_eq!(
            got_bytes, want,
            "canonicalization must preserve the exact record multiset (all fields), not just names"
        );
    }

    #[test]
    fn resolve_sort_tmp_dirs_cli_wins_over_env() {
        let cli = vec![PathBuf::from("/cli/dir")];
        let got = resolve_sort_tmp_dirs(&cli, Some("/env/a:/env/b"));
        assert_eq!(got, cli);
    }

    #[test]
    fn resolve_sort_tmp_dirs_falls_back_to_env_when_cli_empty() {
        #[cfg(unix)]
        let env = "/tmp/x:/tmp/y";
        #[cfg(windows)]
        let env = "C:/tmp/x;C:/tmp/y";

        let got = resolve_sort_tmp_dirs(&[], Some(env));
        assert_eq!(got.len(), 2);
        assert!(got[0].to_string_lossy().ends_with('x'));
        assert!(got[1].to_string_lossy().ends_with('y'));
    }

    /// The F2 fix: when CLI and env both resolve to nothing, this must NOT return
    /// an empty vector (which would defer to `RawExternalSorter`'s bare system-temp
    /// default, i.e. `/tmp` — tmpfs on many modern hosts). It must return an
    /// explicit, disk-backed default instead.
    #[test]
    fn resolve_sort_tmp_dirs_never_defaults_to_empty() {
        assert_eq!(resolve_sort_tmp_dirs(&[], None), vec![default_disk_backed_tmp_dir()]);
        assert_eq!(resolve_sort_tmp_dirs(&[], Some("")), vec![default_disk_backed_tmp_dir()]);
    }

    /// The platform-aware default must resolve to a directory that actually *exists* and
    /// into which a scratch temp dir can be created — not merely a plausible-looking path.
    /// On Windows the old hard-coded `/var/tmp` would not exist and `--command group` would
    /// fail here before ever comparing.
    #[test]
    fn default_disk_backed_tmp_dir_is_a_usable_scratch_root() {
        let base = default_disk_backed_tmp_dir();
        assert!(base.is_dir(), "default scratch root {base:?} must exist as a directory");
        let scratch = tempfile::Builder::new()
            .prefix("fgumi-compare-keyjoin-default-test-")
            .tempdir_in(&base)
            .expect("must be able to create a scratch directory under the default temp dir");
        assert!(scratch.path().is_dir());
    }

    // ---- MiBijectionTracker ---------------------------------------------

    /// `MiBijectionTracker` folds two checks fed one matched pair at a time: the full-`MiKey`
    /// bijection (`mi1_to_mi2`/`mi2_to_mi1`) and the `base()`-level duplex strand-pairing
    /// bijection (`base1_to_base2`/`base2_to_base1`). A pure MI relabel is not a mismatch, but
    /// one bam1 MI mapping to two distinct bam2 MIs is — and a duplex strand-*pairing* split
    /// (caught only via the `base()`-level maps) must be flagged just as reliably as a
    /// full-`MiKey` split, while a mere `/A`/`/B` relabel must not be flagged at all.
    ///
    /// The `split` case's expected count is `2`, not `1`: for a plain `MiKey::Int`,
    /// `base()` is the identity, so a real split is flagged independently by *both* the
    /// full-`MiKey` maps and the `base()`-level maps — the two checks are additive, not
    /// deduplicated. For this exact input the full-key check flags one `mi1_to_mi2` entry and
    /// the base-level check flags one `base1_to_base2` entry (both direction-1, since `Int`
    /// MIs make the base-level check redundant-but-not-wrong here), giving `2`, not `1`.
    ///
    /// `strand_pairing_split` (BAM1 pairs `0/A`+`0/B` into one molecule; BAM2 splits those
    /// same reads into `0/A`+`1/A`) is invisible to the full-`MiKey` maps (each direction is
    /// still a consistent bijection: `0/A<->0/A`, `0/B<->1/A`), so its expected count of `1`
    /// comes entirely from the `base()`-level maps. `ab_strand_swap` (swapping the `/A`/`/B`
    /// labels of a molecule) is a naming difference, not a pairing difference, and must not be
    /// flagged at all. `reverse_split` mirrors `split` in the opposite direction: a single
    /// bam2 MI mapping back to two distinct bam1 MIs, expected count `2` by the same full-key +
    /// base-level double-counting reasoning as `split`.
    #[rstest]
    #[case::consistent_relabel(
        &[(MiKey::Int(0), MiKey::Int(9)), (MiKey::Int(0), MiKey::Int(9))],
        0,
    )]
    #[case::split(&[(MiKey::Int(0), MiKey::Int(9)), (MiKey::Int(0), MiKey::Int(8))], 2)]
    #[case::strand_pairing_split(
        &[
            (MiKey::Strand { base: 0, strand: b'A' }, MiKey::Strand { base: 0, strand: b'A' }),
            (MiKey::Strand { base: 0, strand: b'B' }, MiKey::Strand { base: 1, strand: b'A' }),
        ],
        1,
    )]
    #[case::ab_strand_swap(
        &[
            (MiKey::Strand { base: 0, strand: b'A' }, MiKey::Strand { base: 0, strand: b'B' }),
            (MiKey::Strand { base: 0, strand: b'B' }, MiKey::Strand { base: 0, strand: b'A' }),
        ],
        0,
    )]
    #[case::reverse_split(&[(MiKey::Int(1), MiKey::Int(0)), (MiKey::Int(2), MiKey::Int(0))], 2)]
    fn mi_bijection_tracker_flags_splits(
        #[case] pairs: &[(MiKey, MiKey)],
        #[case] expected_mismatches: u64,
    ) {
        let mut tracker = MiBijectionTracker::default();
        for &(mi1, mi2) in pairs {
            tracker.observe(mi1, mi2);
        }
        assert_eq!(tracker.mismatches(), expected_mismatches);
    }
}
