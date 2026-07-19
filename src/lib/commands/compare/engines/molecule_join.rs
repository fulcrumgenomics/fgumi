//! Per-molecule comparison for the streaming grouping engine.
//!
//! [`compare_molecule`] compares two [`MoleculeRun`]s — the same molecule (matched by
//! canonical id) as read from each of the two BAMs under comparison — via three purely
//! local checks: record membership, per-record content (excluding MI), and the duplex
//! `/A`/`/B` strand partition. It has no notion of streaming or cross-molecule state;
//! that is the job of [`molecule_join_compare`], the streaming two-sided hash-join driver
//! that pairs up [`MoleculeRun`]s by canonical id (across the two whole files, independent
//! of file order or MI numbering) and feeds each matched pair here.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::path::Path;

use ahash::AHashMap;
use anyhow::{Context, Result, bail};
use fgumi_bam_io::create_raw_bam_reader;
use fgumi_raw_bam::RawRecord;
use fgumi_sort::RawBamRecordReader;
use noodles::sam::Header;

use super::super::bams::{MiKey, get_mi_tag_raw};
use super::super::molecule::{MoleculeRun, molecule_runs};
use super::super::record_key::{RecordKey, record_key};
use super::content::{ContentPredicate, content_diffs};
use super::header::require_compatible_headers;
use super::push_diff;

/// Build a `RecordKey -> [&RawRecord]` multiset index for a molecule's members.
///
/// `RecordKey` is collision-resistant but not collision-free (see its module docs): if two
/// members of one molecule genuinely share a key, both are kept here (in file order), rather
/// than one silently winning over the other. This makes multiplicity — how many physical
/// records share a given key — a first-class, comparable quantity: a key present a different
/// number of times on the two sides of a matched molecule is a genuine divergence (e.g. a
/// duplicated or dropped record), not something `compare_molecule` should be blind to.
fn index_by_key(members: &[RawRecord]) -> BTreeMap<RecordKey, Vec<&RawRecord>> {
    let mut index: BTreeMap<RecordKey, Vec<&RawRecord>> = BTreeMap::new();
    for r in members {
        index.entry(record_key(r)).or_default().push(r);
    }
    index
}

/// Returns `true` iff the two molecules' duplex `/A`/`/B` strand partitions are equivalent.
///
/// Groups each molecule's members into an `/A` and a `/B` multiset — a `RecordKey ->
/// occurrence count` map, not a plain set — so two records that happen to share a
/// `RecordKey` within the same strand are not silently collapsed into one (records with no
/// strand suffix — a plain [`MiKey::Int`] or no MI tag at all — contribute to neither
/// multiset and so trivially pass this check). Returns `true` iff `{a_A, a_B}` equals
/// `{b_A, b_B}` as an **unordered pair** — i.e. either the multisets line up directly, or
/// they line up after a single global A<->B relabel — so a molecule-wide strand-label swap
/// between the two files is accepted, while an actual strand-*pairing* split (a read moving
/// from one strand's set to the other, or between molecules) is not.
///
/// Known limitation (deliberately not addressed here): the per-strand multisets are keyed by
/// `RecordKey`, not content. If two *distinct* records collide on one `RecordKey` (e.g. a
/// chimeric read's several SUPPLEMENTARY alignments) and swap strands between the two files,
/// the counts still balance and this check reports equivalence — a false MATCH for that
/// contrived case. Making the strand multisets content-aware would close it; the reachable
/// half of the collision problem (a content false-DIFF in `compare_molecule`'s membership/
/// content check) is already handled there by content-multiset matching.
fn strand_partitions_equivalent(a: &[RawRecord], b: &[RawRecord]) -> bool {
    fn strand_multisets(
        members: &[RawRecord],
    ) -> (BTreeMap<RecordKey, usize>, BTreeMap<RecordKey, usize>) {
        let mut a_set: BTreeMap<RecordKey, usize> = BTreeMap::new();
        let mut b_set: BTreeMap<RecordKey, usize> = BTreeMap::new();
        for member in members {
            if let Some(MiKey::Strand { strand, .. }) = get_mi_tag_raw(member) {
                match strand {
                    b'A' => *a_set.entry(record_key(member)).or_insert(0) += 1,
                    b'B' => *b_set.entry(record_key(member)).or_insert(0) += 1,
                    _ => unreachable!("get_mi_tag_raw only returns b'A'/b'B' strand bytes"),
                }
            }
        }
        (a_set, b_set)
    }

    let (a_a, a_b) = strand_multisets(a);
    let (b_a, b_b) = strand_multisets(b);
    (a_a == b_a && a_b == b_b) || (a_a == b_b && a_b == b_a)
}

/// Compare two molecules matched by canonical id, pushing bounded diagnostics into
/// `diff_details` (capped at `max_diffs` via [`push_diff`]) and returning `true` iff the
/// molecules diverge on any of the three checks below.
///
/// The returned `bool` is deliberately **uncapped** — it reflects every divergence found,
/// even once `diff_details` has hit `max_diffs` and stopped accepting new lines — so the
/// caller's match/no-match verdict (`MoleculeJoinOutcome::matched`, via [`fold_matched_pair`])
/// stays correct regardless of `max_diffs` (including `max_diffs == 0`; see
/// [`MoleculeJoinOutcome::is_match`]'s doc comment for why that independence matters). Diff
/// strings themselves are still built lazily by [`push_diff`], so a fully-diverged molecule
/// under a small `max_diffs` never allocates more than `max_diffs` diagnostic strings.
///
/// Three local checks, each independent of the others (a molecule can fail more than one):
///
/// 1. **Membership** — the two molecules must contain the same [`RecordKey`] *multiset*: a
///    key present a different number of times on the two sides (including zero, i.e.
///    present on only one side) is a diff naming the molecule and the multiplicities.
/// 2. **Content** — every key present (the same number of times) on both sides has its
///    instances compared under [`ContentPredicate::ExactMinusMi`] (MI numbering may
///    legitimately differ across tools/runs; see [`super::super::bams::MiKey`]'s doc comment).
///    Multiplicity 1 (the overwhelmingly common case) is a direct pairwise comparison that
///    reports precise field diffs; a genuine `RecordKey` collision (multiplicity > 1) is
///    compared as a *content multiset* — instances are matched by content, not file order, so
///    colliding records that merely reordered between the files do not read as a diff.
/// 3. **Strand partition** — the duplex `/A`/`/B` split must be equivalent modulo a single
///    global A<->B relabel; see [`strand_partitions_equivalent`].
pub(crate) fn compare_molecule(
    a: &MoleculeRun,
    b: &MoleculeRun,
    ha: &Header,
    hb: &Header,
    max_diffs: usize,
    diff_details: &mut Vec<String>,
) -> bool {
    let mut differed = false;
    let canon = String::from_utf8_lossy(&a.canon);

    let keys_a = index_by_key(&a.members);
    let keys_b = index_by_key(&b.members);

    // 1 & 2. Multiset membership, then content of every key present (equally) on both sides.
    let mut all_keys: BTreeSet<RecordKey> = keys_a.keys().cloned().collect();
    all_keys.extend(keys_b.keys().cloned());
    for k in &all_keys {
        let ra: &[&RawRecord] = keys_a.get(k).map_or(&[], Vec::as_slice);
        let rb: &[&RawRecord] = keys_b.get(k).map_or(&[], Vec::as_slice);
        if ra.is_empty() {
            differed = true;
            push_diff(diff_details, max_diffs, || {
                format!("molecule {canon}: record {k:?} only in bam2 (x{})", rb.len())
            });
        } else if rb.is_empty() {
            differed = true;
            push_diff(diff_details, max_diffs, || {
                format!("molecule {canon}: record {k:?} only in bam1 (x{})", ra.len())
            });
        } else if ra.len() != rb.len() {
            differed = true;
            push_diff(diff_details, max_diffs, || {
                format!(
                    "molecule {canon}: record {k:?} multiplicity differs: {} in bam1 vs {} in bam2",
                    ra.len(),
                    rb.len()
                )
            });
        } else if ra.len() == 1 {
            // Common case: exactly one instance of this key on each side. Compare directly and
            // report the precise field-level diffs.
            if let Some(cd) = content_diffs(
                ra[0].as_ref(),
                rb[0].as_ref(),
                ContentPredicate::ExactMinusMi,
                ha,
                hb,
            ) {
                differed = true;
                push_diff(diff_details, max_diffs, || {
                    format!("molecule {canon}: record {k:?}: {}", cd.join("; "))
                });
            }
        } else {
            // Genuine `RecordKey` collision: multiple records share this key (e.g. a chimeric
            // read's several SUPPLEMENTARY alignments). `RecordKey` is not content, so pairing
            // the instances by *file order* would false-DIFF two colliding records that merely
            // reordered between the files. Compare the bucket as a CONTENT multiset instead:
            // greedily match each bam1 instance to an as-yet-unmatched, `ExactMinusMi`-equal
            // bam2 instance; any bam1 instance left unmatched is a real content divergence (and,
            // by equal multiplicity, so is some bam2 instance). N is the collision multiplicity
            // (tiny), so the O(N^2) greedy match is fine.
            let mut rb_matched = vec![false; rb.len()];
            for ia in ra {
                let hit = rb.iter().enumerate().position(|(j, ib)| {
                    !rb_matched[j]
                        && content_diffs(
                            ia.as_ref(),
                            ib.as_ref(),
                            ContentPredicate::ExactMinusMi,
                            ha,
                            hb,
                        )
                        .is_none()
                });
                if let Some(j) = hit {
                    rb_matched[j] = true;
                } else {
                    differed = true;
                    push_diff(diff_details, max_diffs, || {
                        format!(
                            "molecule {canon}: record {k:?}: a bam1 instance has no \
                             content-equivalent counterpart among the {} bam2 instance(s) \
                             sharing this key",
                            rb.len()
                        )
                    });
                }
            }
        }
    }

    // 3. Duplex strand partition (accepted modulo one global A<->B relabel).
    if !strand_partitions_equivalent(&a.members, &b.members) {
        differed = true;
        push_diff(diff_details, max_diffs, || {
            format!("molecule {canon}: duplex strand partition differs")
        });
    }

    differed
}

/// Outcome of a [`molecule_join_compare`] run.
#[derive(Debug, Default)]
pub struct MoleculeJoinOutcome {
    /// Total number of molecule runs read from `bam1`.
    pub bam1_molecules: u64,
    /// Total number of molecule runs read from `bam2`.
    pub bam2_molecules: u64,
    /// Number of molecules matched by canonical id across both files and found equivalent
    /// by `compare_molecule`.
    pub matched: u64,
    /// Human-readable diff strings, capped at the caller-supplied `max_diffs`: per-molecule
    /// `compare_molecule` divergences, plus one entry per residual molecule present in
    /// only one of the two files.
    pub diff_details: Vec<String>,
    /// `true` if the two inputs' headers disagreed (reserved; always `false` today — a
    /// header disagreement is instead a hard `Err` returned directly from
    /// [`molecule_join_compare`] via its own `require_compatible_headers` precondition
    /// check, so this outcome is never actually constructed for a header-incompatible pair;
    /// this field stays `false` by construction, not because it's checked and found clean).
    pub header_mismatch: bool,
}

impl MoleculeJoinOutcome {
    /// Returns `true` iff every molecule on both sides matched cleanly: `matched` (which
    /// increments in `fold_matched_pair` only when `compare_molecule` returns zero
    /// diffs for a canonical-id-matched pair) equals *both* `bam1_molecules` and
    /// `bam2_molecules`, and no header divergence was found.
    ///
    /// This is deliberately keyed off the `matched` counter rather than
    /// `diff_details.is_empty()`. `diff_details` is populated exclusively via `push_diff`,
    /// which caps at the caller-supplied `max_diffs` (`details.len() < max_diffs`) — with
    /// `max_diffs == 0` (a plausible CI invocation: `--max-diffs` is a plain `usize` with no
    /// range validator) `diff_details` is *always* empty regardless of how many real diffs
    /// were found. A verdict keyed off `diff_details.is_empty()` would therefore collapse to
    /// `bam1_molecules == bam2_molecules`, silently reporting EQUIVALENT for a genuine
    /// content/membership/strand divergence between matched molecules, or for a
    /// present-in-only-one-file residual on each side that happens to balance the counts.
    /// `matched` is never subject to the cap — it is incremented directly in
    /// `fold_matched_pair`, independent of `diff_details`/`max_diffs` entirely — so this
    /// check is sound for every `max_diffs` value, including `0`. A residual (only-in-one-
    /// file) molecule is counted in `bam1_molecules`/`bam2_molecules` but never in `matched`,
    /// so it correctly fails this equality too.
    #[must_use]
    pub fn is_match(&self) -> bool {
        self.matched == self.bam1_molecules
            && self.matched == self.bam2_molecules
            && !self.header_mismatch
    }
}

/// Open a fresh raw-byte record reader over `path`, positioned just past the header, via
/// the shared [`super::open_raw_bam_reader`] helper.
fn open_raw(path: &Path) -> Result<RawBamRecordReader<File>> {
    super::open_raw_bam_reader(path)
        .with_context(|| format!("opening raw BAM reader for {}", path.display()))
}

/// Fold one matched pair (already known bam1-run/bam2-run in that order) into `matched`/
/// `diff_details` via [`compare_molecule`]. Shared by both sides of [`molecule_join_compare`]'s
/// hash-join loop so the "compare, then either count as matched or record diffs" step isn't
/// duplicated between the two symmetric branches.
///
/// `compare_molecule` pushes its own (capped) diagnostic lines directly into `diff_details`
/// and returns an uncapped `differed` bool; `matched` increments only when that bool is
/// `false`, so the verdict stays sound even when `max_diffs` caps `diff_details` to nothing
/// (e.g. `max_diffs == 0`) — see [`compare_molecule`]'s doc comment.
fn fold_matched_pair(
    bam1_run: &MoleculeRun,
    bam2_run: &MoleculeRun,
    h1: &Header,
    h2: &Header,
    max_diffs: usize,
    matched: &mut u64,
    diff_details: &mut Vec<String>,
) {
    let differed = compare_molecule(bam1_run, bam2_run, h1, h2, max_diffs, diff_details);
    if !differed {
        *matched += 1;
    }
}

/// Two-sided streaming hash-join over `bam1`'s and `bam2`'s per-molecule runs (see
/// [`molecule_runs`]), comparing each pair matched by canonical id via [`compare_molecule`].
///
/// Unlike the retired key-join engine, this never re-sorts either input: both files are
/// already template-coordinate sorted with same-molecule reads consecutive (grouped output), so
/// [`molecule_runs`] can cut each into per-molecule runs in a single streaming pass. The two
/// streams are pulled in lockstep — one run from `bam1`, then one from `bam2`, each iteration —
/// each newly pulled run probed against the *other* side's `pending` map: a hit compares and
/// evicts the pair, a miss is buffered in the puller's own `pending` map awaiting its
/// partner. Because both inputs share the same underlying coordinate order, a molecule's
/// counterpart typically arrives at a nearby offset in the other file, so `pending` stays
/// small in practice. It is not windowed or silently dropped — a legitimately large transient
/// backlog is fully buffered — but the combined backlog is bounded by a deliberately high
/// last-resort cap ([`MAX_PENDING_MOLECULES`]): if two inputs are grouped in very different
/// orders, a molecule's partner never arrives at a nearby offset and the join would otherwise
/// buffer O(file-size) molecules and OOM, so exceeding the cap is turned into a clear error
/// rather than an eventual crash. Realistic co-sorted grouped BAMs never approach it. Once both
/// streams are exhausted, any molecule left in either `pending` map was present in only one file
/// and is reported as a diff naming its canonical id.
///
/// The two per-side steps below are written out separately, rather than driven from a shared
/// array of `(iterator, own map, other map, ...)` tuples: `pending1`/`pending2` each need to
/// be borrowed mutably as *both* "own" (for their side) and "other" (for the opposite side)
/// within the same iteration, which a single array literal can't express without borrowing
/// one of the two maps mutably twice at once.
///
/// Header compatibility: the CLI entry point (`CompareBams::execute`) already gates on
/// `require_compatible_headers` before dispatching to this driver, but `molecule_join_compare`
/// is a `pub` function and must be sound for any caller, not just the CLI — so it re-runs the
/// same check itself, immediately after opening both headers below. This is redundant (but
/// cheap — one more header parse) for the CLI path; for any other caller it is the only thing
/// standing between an incompatible `@SQ`/`@RG`/sort-order pair and a molecule-join that
/// silently pairs records against the wrong reference dictionary.
///
/// # Errors
///
/// Returns an error if either file cannot be opened/read, if the two inputs' headers are
/// incompatible (see [`require_compatible_headers`]), or if [`molecule_runs`] yields an
/// `Err` from either stream. On an `Err` from either side, the join stops immediately without
/// polling either iterator again: `molecule_runs` does not fuse after an `Err` (a stale
/// `pending` record inside it would otherwise silently leak into the next run's boundary), so
/// resuming past a first error would corrupt subsequent run boundaries.
///
/// Grouping mode compares *grouped* output (same-molecule reads consecutive under a shared
/// MI), so any non-empty input with a record that lacks a (parseable) MI is misuse of this
/// comparison, not a legitimate "no diffs" input. That precondition is enforced upstream, per
/// record, by [`molecule_runs`]: it yields an `Err` on the first MI-less record, which this
/// driver propagates via `?` (see above). Enforcing it per record — rather than with a
/// whole-input "never saw an MI" guard here — is what closes the *partial*-MI-less soundness
/// hole: a side with some records tagged and some not (or an entirely-MI-less side whose
/// single spanning run happened to canonical-id/content-match a real molecule on the grouped
/// side) would otherwise report a false MATCH despite never having verified grouping, because
/// the join compares content *excluding* MI. Two entirely *empty* BAMs (zero records) yield no
/// runs and no MI-less record, so they are never rejected — there is nothing to compare, a
/// vacuous MATCH.
/// Last-resort bound on the combined molecule backlog held in `pending1` + `pending2` (see
/// [`molecule_join_compare`]'s doc comment). This is **not** a tuning knob: for correctly
/// co-sorted grouped inputs the backlog stays small (a molecule's partner arrives at a nearby
/// offset), so the cap is set deliberately far above any legitimate transient backlog. Its only
/// job is to convert a pathological O(file-size) backlog — two inputs grouped in very different
/// orders — from an eventual OOM into a clear, actionable error. Realistic grouped BAMs, even
/// whole-genome ones, never approach 100M buffered molecules.
pub(crate) const MAX_PENDING_MOLECULES: usize = 100_000_000;

pub fn molecule_join_compare(
    bam1: &Path,
    bam2: &Path,
    max_diffs: usize,
) -> Result<MoleculeJoinOutcome> {
    molecule_join_compare_capped(bam1, bam2, max_diffs, MAX_PENDING_MOLECULES)
}

/// [`molecule_join_compare`] with an explicit backlog cap, so tests can exercise the
/// over-cap error path without materializing 100M molecules. Production always calls through
/// the public wrapper with [`MAX_PENDING_MOLECULES`].
pub(crate) fn molecule_join_compare_capped(
    bam1: &Path,
    bam2: &Path,
    max_diffs: usize,
    max_pending: usize,
) -> Result<MoleculeJoinOutcome> {
    let (_, h1) = create_raw_bam_reader(bam1, 1)?;
    let (_, h2) = create_raw_bam_reader(bam2, 1)?;
    // Make this public API sound on its own, independent of the CLI's own (redundant but
    // cheap) call — see this function's doc comment.
    require_compatible_headers(&h1, &h2)?;

    let mut it1 = molecule_runs(open_raw(bam1)?);
    let mut it2 = molecule_runs(open_raw(bam2)?);
    let mut pending1: AHashMap<Vec<u8>, MoleculeRun> = AHashMap::new();
    let mut pending2: AHashMap<Vec<u8>, MoleculeRun> = AHashMap::new();
    let mut done1 = false;
    let mut done2 = false;

    let mut bam1_molecules = 0u64;
    let mut bam2_molecules = 0u64;
    let mut matched = 0u64;
    let mut diff_details: Vec<String> = Vec::new();

    // Alternate pulling one run from each side; probe the opposite pending map on each pull.
    // Because both inputs are template-coordinate ordered, corresponding molecules arrive at
    // ~the same offset, so `pending` stays small; the combined backlog is bounded only by the
    // deliberately high `max_pending` last-resort cap checked at the end of each iteration.
    // On the first `Err` from either side, `?` returns immediately without polling either
    // iterator again (see this function's doc comment on why resuming would corrupt run
    // boundaries).
    while !done1 || !done2 {
        if !done1 {
            match it1.next() {
                None => done1 = true,
                Some(run) => {
                    let run = run?;
                    bam1_molecules += 1;
                    if let Some(partner) = pending2.remove(&run.canon) {
                        fold_matched_pair(
                            &run,
                            &partner,
                            &h1,
                            &h2,
                            max_diffs,
                            &mut matched,
                            &mut diff_details,
                        );
                    } else {
                        pending1.insert(run.canon.clone(), run);
                    }
                }
            }
        }
        if !done2 {
            match it2.next() {
                None => done2 = true,
                Some(run) => {
                    let run = run?;
                    bam2_molecules += 1;
                    if let Some(partner) = pending1.remove(&run.canon) {
                        fold_matched_pair(
                            &partner,
                            &run,
                            &h1,
                            &h2,
                            max_diffs,
                            &mut matched,
                            &mut diff_details,
                        );
                    } else {
                        pending2.insert(run.canon.clone(), run);
                    }
                }
            }
        }

        // Last-resort backlog bound (see `MAX_PENDING_MOLECULES`). Checked once per iteration
        // after both pulls: for co-sorted grouped input the combined backlog stays small, so
        // this never fires; it exists to turn a pathological O(file-size) backlog (two inputs
        // grouped in very different orders) into a clear error instead of an eventual OOM.
        if pending1.len() + pending2.len() > max_pending {
            bail!(
                "grouping comparison backlog exceeded {max_pending} buffered molecules: the two \
                 inputs appear not to be co-sorted (a molecule's partner is not arriving at a \
                 nearby offset in the other file), so the streaming molecule-join would buffer \
                 O(file-size) molecules. Ensure both inputs are template-coordinate sorted \
                 grouped output produced with the same ordering."
            );
        }
    }

    // The "every non-empty input must be MI-tagged grouped data" precondition is enforced
    // upstream, per record, by `molecule_runs` (its `Err` is propagated by the `?` in the pull
    // loop above), so there is no whole-input MI-less guard here — see this function's doc
    // comment on why per-record rejection, not a post-join guard, is what closes the
    // partial-MI-less false-MATCH hole.

    // Residual molecules present in only one file. `pending1`/`pending2` are `AHashMap`s, so
    // their own iteration order is not deterministic across runs; sort the canonical ids
    // before reporting them so both *which* ids survive the `max_diffs` cap and the order
    // they're reported in are stable run-to-run, not an artifact of hash iteration order.
    let mut residual1: Vec<Vec<u8>> = pending1.into_keys().collect();
    residual1.sort();
    for canon in residual1 {
        push_diff(&mut diff_details, max_diffs, || {
            format!("molecule {} only in bam1", String::from_utf8_lossy(&canon))
        });
    }
    let mut residual2: Vec<Vec<u8>> = pending2.into_keys().collect();
    residual2.sort();
    for canon in residual2 {
        push_diff(&mut diff_details, max_diffs, || {
            format!("molecule {} only in bam2", String::from_utf8_lossy(&canon))
        });
    }

    Ok(MoleculeJoinOutcome {
        bam1_molecules,
        bam2_molecules,
        matched,
        diff_details,
        header_mismatch: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::SamTag;
    use fgumi_raw_bam::{SamBuilder, flags};
    use rstest::rstest;

    /// A minimal header sufficient for `content_diffs`' RNAME/RNEXT resolution on the
    /// unmapped records these tests build.
    fn minimal_header() -> Header {
        Header::default()
    }

    /// Builds a `MoleculeRun` from `(name, mi)` pairs. `mi` follows the same string-tag
    /// convention as `molecule.rs`'s own test helper: a plain integer (`"5"`) or a duplex
    /// strand-suffixed id (`"5/A"`/`"5/B"`). `canon` is the lexicographically smallest name,
    /// matching `MoleculeRun::from_members`'s own convention.
    fn molecule_run_from(specs: &[(&str, &str)]) -> MoleculeRun {
        let members: Vec<RawRecord> = specs
            .iter()
            .map(|(name, mi)| {
                SamBuilder::new()
                    .read_name(name.as_bytes())
                    .flags(flags::FIRST_SEGMENT)
                    .add_string_tag(SamTag::MI, mi.as_bytes())
                    .build()
            })
            .collect();
        let canon = members.iter().map(|r| r.read_name().to_vec()).min().unwrap_or_default();
        MoleculeRun { canon, members }
    }

    #[rstest]
    // same reads, reordered within the molecule, MI renumbered → equivalent (empty)
    #[case::reorder_and_renumber(&[("r1","5"),("r2","5")], &[("r2","9"),("r1","9")], true)]
    // one file grouped an extra read into the molecule → membership diff
    #[case::membership_differs(&[("r1","5"),("r2","5")], &[("r1","9")], false)]
    // duplex: A/B strand partition disagrees → strand diff (modulo global A<->B is still equal)
    #[case::strand_relabel_ok(&[("r1","5/A"),("r2","5/B")], &[("r1","9/B"),("r2","9/A")], true)]
    #[case::strand_split(&[("r1","5/A"),("r2","5/B")], &[("r1","9/A"),("r2","9/A")], false)]
    // multiset regression: bam1 has TWO "dup"-named records (same RecordKey, since a
    // primary's identity is (name, segment) only), bam2 has ONE → a 2-vs-1 multiplicity
    // diff, not a silent BTreeMap-collapse MATCH.
    #[case::duplicate_key_multiplicity_differs(
        &[("r1", "5"), ("dup", "5"), ("dup", "5")],
        &[("r1", "9"), ("dup", "9")],
        false
    )]
    fn compare_molecule_cases(
        #[case] a: &[(&str, &str)],
        #[case] b: &[(&str, &str)],
        #[case] equal: bool,
    ) {
        let (ha, hb) = (minimal_header(), minimal_header());
        let ra = molecule_run_from(a);
        let rb = molecule_run_from(b);
        let mut diff_details = Vec::new();
        let differed = compare_molecule(&ra, &rb, &ha, &hb, usize::MAX, &mut diff_details);
        assert_eq!(!differed, equal, "diffs: {diff_details:?}");
    }

    /// `compare_molecule`'s returned `differed` bool must stay accurate even when
    /// `max_diffs == 0` suppresses every diagnostic line — the caller's match/no-match
    /// verdict must never depend on whether any diff string was actually kept.
    #[test]
    fn compare_molecule_differed_is_uncapped_by_max_diffs_zero() {
        let (ha, hb) = (minimal_header(), minimal_header());
        let ra = molecule_run_from(&[("r1", "5"), ("r2", "5")]);
        let rb = molecule_run_from(&[("r1", "9")]); // r2 missing on bam2's side
        let mut diff_details = Vec::new();
        let differed = compare_molecule(&ra, &rb, &ha, &hb, 0, &mut diff_details);
        assert!(differed, "a real membership diff must be reported even under max_diffs == 0");
        assert!(diff_details.is_empty(), "max_diffs == 0 must still suppress diagnostic lines");
    }

    /// Build one MI-tagged, single-record molecule (base MI `mi`, read name `name`).
    fn mi_rec(name: &[u8], mi: &str) -> RawRecord {
        SamBuilder::new()
            .read_name(name)
            .flags(flags::FIRST_SEGMENT)
            .add_string_tag(SamTag::MI, mi.as_bytes())
            .build()
    }

    /// Build an MI-tagged record with a specific `seq` (and matching all-30 quals). Two of
    /// these sharing `name` (and the `FIRST_SEGMENT` flag) but differing in `seq` collide on
    /// `RecordKey` — identity is (name, segment) only — while differing in content.
    fn mi_seq_rec(name: &[u8], mi: &str, seq: &[u8]) -> RawRecord {
        SamBuilder::new()
            .read_name(name)
            .flags(flags::FIRST_SEGMENT)
            .sequence(seq)
            .qualities(&vec![30u8; seq.len()])
            .add_string_tag(SamTag::MI, mi.as_bytes())
            .build()
    }

    /// Item 9a: two records that collide on `RecordKey` (same name+segment, different content)
    /// but merely REORDER between the two files must MATCH — the colliding bucket is compared
    /// as a content multiset, not zipped in file order. Before the fix this file-order zip
    /// paired `X` with `Y` and false-DIFFered.
    #[test]
    fn colliding_records_reordered_within_bucket_still_match() {
        let (ha, hb) = (minimal_header(), minimal_header());
        let a = MoleculeRun {
            canon: b"dup".to_vec(),
            members: vec![
                mi_seq_rec(b"dup", "5", b"AAAAAAAA"),
                mi_seq_rec(b"dup", "5", b"CCCCCCCC"),
            ],
        };
        let b = MoleculeRun {
            canon: b"dup".to_vec(),
            members: vec![
                mi_seq_rec(b"dup", "5", b"CCCCCCCC"),
                mi_seq_rec(b"dup", "5", b"AAAAAAAA"),
            ],
        };
        let mut diffs = Vec::new();
        let differed = compare_molecule(&a, &b, &ha, &hb, usize::MAX, &mut diffs);
        assert!(!differed, "reordered colliding records must match, got diffs: {diffs:?}");
    }

    /// The genuine-divergence counterpart: the same colliding key, but bam2's content multiset
    /// differs ({A, A} vs bam1's {A, C}). The `A` matches, the `C` finds no counterpart → DIFF.
    /// Confirms the content-multiset match still catches a real difference (no false MATCH).
    #[test]
    fn colliding_records_content_multiset_differs() {
        let (ha, hb) = (minimal_header(), minimal_header());
        let a = MoleculeRun {
            canon: b"dup".to_vec(),
            members: vec![
                mi_seq_rec(b"dup", "5", b"AAAAAAAA"),
                mi_seq_rec(b"dup", "5", b"CCCCCCCC"),
            ],
        };
        let b = MoleculeRun {
            canon: b"dup".to_vec(),
            members: vec![
                mi_seq_rec(b"dup", "5", b"AAAAAAAA"),
                mi_seq_rec(b"dup", "5", b"AAAAAAAA"),
            ],
        };
        let mut diffs = Vec::new();
        let differed = compare_molecule(&a, &b, &ha, &hb, usize::MAX, &mut diffs);
        assert!(differed, "a differing colliding-key content multiset must DIFFER");
        assert!(
            diffs.iter().any(|d| d.contains("no content-equivalent counterpart")),
            "expected the collision-bucket diagnostic, got: {diffs:?}"
        );
    }

    /// Write `records` to a temp BAM (unmapped, so a default header suffices) and return the
    /// handle, which must be kept alive for the file to persist. Modeled on `molecule.rs`'s
    /// `try_collect_runs` test writer.
    fn write_temp_bam(records: &[RawRecord]) -> tempfile::NamedTempFile {
        use noodles::sam::alignment::io::Write as _;
        let tmp = tempfile::NamedTempFile::new().expect("create temp BAM");
        let header = Header::default();
        let mut writer =
            noodles::bam::io::Writer::new(std::fs::File::create(tmp.path()).expect("create BAM"));
        writer.write_header(&header).expect("write header");
        for record in records {
            let buf = fgumi_raw_bam::raw_record_to_record_buf(record, &header)
                .expect("raw_record_to_record_buf");
            writer.write_alignment_record(&header, &buf).expect("write record");
        }
        writer.try_finish().expect("finish BAM");
        tmp
    }

    /// The last-resort backlog cap turns a pathological *non-co-sorted* pair — disjoint
    /// canonical ids on every molecule, so no partner ever arrives and both `pending` maps grow
    /// unbounded — into a clear error instead of an eventual OOM. A tiny cap makes the trip
    /// cheap to exercise; production uses [`MAX_PENDING_MOLECULES`].
    #[test]
    fn pending_backlog_over_cap_errors() {
        let bam1 = write_temp_bam(&[
            mi_rec(b"left1", "1"),
            mi_rec(b"left2", "2"),
            mi_rec(b"left3", "3"),
            mi_rec(b"left4", "4"),
        ]);
        let bam2 = write_temp_bam(&[
            mi_rec(b"right1", "1"),
            mi_rec(b"right2", "2"),
            mi_rec(b"right3", "3"),
            mi_rec(b"right4", "4"),
        ]);
        let err = molecule_join_compare_capped(bam1.path(), bam2.path(), 10, 2)
            .expect_err("disjoint-canon pair must exceed the tiny backlog cap");
        assert!(err.to_string().contains("backlog exceeded"), "got: {err}");
    }

    /// A co-sorted pair (identical canonical ids, each matched and evicted immediately) keeps
    /// the combined backlog at zero, so even a cap of 1 never trips: the cap targets ordering
    /// divergence, not input size.
    #[test]
    fn co_sorted_pair_does_not_trip_cap() {
        let recs = [mi_rec(b"mol1", "1"), mi_rec(b"mol2", "2"), mi_rec(b"mol3", "3")];
        let bam1 = write_temp_bam(&recs);
        let bam2 = write_temp_bam(&recs);
        let outcome = molecule_join_compare_capped(bam1.path(), bam2.path(), 10, 1)
            .expect("co-sorted identical pair must not trip the cap");
        assert!(outcome.is_match(), "identical grouped BAMs should MATCH: {outcome:?}");
    }
}
