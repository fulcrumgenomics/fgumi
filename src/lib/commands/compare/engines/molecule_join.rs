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
use super::push_diff;

/// Build a `RecordKey -> &RawRecord` index for a molecule's members.
///
/// `RecordKey` is collision-resistant but not collision-free (see its module docs); if two
/// members of one molecule genuinely share a key, the later one wins in this index. That
/// residual is accepted here, matching the key-join engine's own use of `RecordKey`.
fn index_by_key(members: &[RawRecord]) -> BTreeMap<RecordKey, &RawRecord> {
    members.iter().map(|r| (record_key(r), r)).collect()
}

/// Returns `true` iff the two molecules' duplex `/A`/`/B` strand partitions are equivalent.
///
/// Groups each molecule's members into an `/A` set and a `/B` set of [`RecordKey`]s (records
/// with no strand suffix — a plain [`MiKey::Int`] or no MI tag at all — contribute to neither
/// set and so trivially pass this check). Returns `true` iff `{a_A, a_B}` equals `{b_A, b_B}`
/// as an **unordered pair** — i.e. either the sets line up directly, or they line up after a
/// single global A<->B relabel — so a molecule-wide strand-label swap between the two files is
/// accepted, while an actual strand-*pairing* split (a read moving from one strand's set to
/// the other, or between molecules) is not.
fn strand_partitions_equivalent(a: &[RawRecord], b: &[RawRecord]) -> bool {
    fn strand_sets(members: &[RawRecord]) -> (BTreeSet<RecordKey>, BTreeSet<RecordKey>) {
        let mut a_set = BTreeSet::new();
        let mut b_set = BTreeSet::new();
        for member in members {
            if let Some(MiKey::Strand { strand, .. }) = get_mi_tag_raw(member) {
                match strand {
                    b'A' => {
                        a_set.insert(record_key(member));
                    }
                    b'B' => {
                        b_set.insert(record_key(member));
                    }
                    _ => unreachable!("get_mi_tag_raw only returns b'A'/b'B' strand bytes"),
                }
            }
        }
        (a_set, b_set)
    }

    let (a_a, a_b) = strand_sets(a);
    let (b_a, b_b) = strand_sets(b);
    (a_a == b_a && a_b == b_b) || (a_a == b_b && a_b == b_a)
}

/// Compare two molecules matched by canonical id. Returns a diff line per divergence;
/// an empty result means the molecules are equivalent.
///
/// Three local checks, each independent of the others (a molecule can fail more than one):
///
/// 1. **Membership** — the two molecules must contain the same set of [`RecordKey`]
///    identities. A key present on only one side is a diff naming the molecule.
/// 2. **Content** — every record present on both sides is compared under
///    [`ContentPredicate::ExactMinusMi`] (MI numbering may legitimately differ across
///    tools/runs; see [`super::super::bams::MiKey`]'s doc comment).
/// 3. **Strand partition** — the duplex `/A`/`/B` split must be equivalent modulo a single
///    global A<->B relabel; see [`strand_partitions_equivalent`].
pub(crate) fn compare_molecule(
    a: &MoleculeRun,
    b: &MoleculeRun,
    ha: &Header,
    hb: &Header,
) -> Vec<String> {
    let mut diffs = Vec::new();
    let canon = String::from_utf8_lossy(&a.canon);

    // 1. Membership: same set of record identities.
    let keys_a = index_by_key(&a.members);
    let keys_b = index_by_key(&b.members);
    for k in keys_a.keys() {
        if !keys_b.contains_key(k) {
            diffs.push(format!("molecule {canon}: record {k:?} only in bam1"));
        }
    }
    for k in keys_b.keys() {
        if !keys_a.contains_key(k) {
            diffs.push(format!("molecule {canon}: record {k:?} only in bam2"));
        }
    }

    // 2. Content of shared reads (ExactMinusMi — MI numbering may legitimately differ).
    for (k, ra) in &keys_a {
        if let Some(rb) = keys_b.get(k) {
            if let Some(cd) =
                content_diffs(ra.as_ref(), rb.as_ref(), ContentPredicate::ExactMinusMi, ha, hb)
            {
                diffs.push(format!("molecule {canon}: record {k:?}: {}", cd.join("; ")));
            }
        }
    }

    // 3. Duplex strand partition (accepted modulo one global A<->B relabel).
    if !strand_partitions_equivalent(&a.members, &b.members) {
        diffs.push(format!("molecule {canon}: duplex strand partition differs"));
    }

    diffs
}

/// Outcome of a [`molecule_join_compare`] run.
#[derive(Debug, Default)]
pub struct MoleculeJoinOutcome {
    /// Total number of molecule runs read from `bam1`.
    pub bam1_molecules: u64,
    /// Total number of molecule runs read from `bam2`.
    pub bam2_molecules: u64,
    /// Number of molecules matched by canonical id across both files and found equivalent
    /// by [`compare_molecule`].
    pub matched: u64,
    /// Human-readable diff strings, capped at the caller-supplied `max_diffs`: per-molecule
    /// [`compare_molecule`] divergences, plus one entry per residual molecule present in
    /// only one of the two files.
    pub diff_details: Vec<String>,
    /// `true` if the two inputs' headers disagreed (reserved for a future header check on
    /// this path; always `false` today — see the module-level note on the CLI entry
    /// already enforcing header compatibility before this driver runs).
    pub header_mismatch: bool,
}

impl MoleculeJoinOutcome {
    /// Returns `true` iff every matched molecule compared equivalent, no molecule was
    /// present in only one file, no header divergence was found, and both files contained
    /// the same number of molecules.
    #[must_use]
    pub fn is_match(&self) -> bool {
        self.diff_details.is_empty()
            && !self.header_mismatch
            && self.bam1_molecules == self.bam2_molecules
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
fn fold_matched_pair(
    bam1_run: &MoleculeRun,
    bam2_run: &MoleculeRun,
    h1: &Header,
    h2: &Header,
    max_diffs: usize,
    matched: &mut u64,
    diff_details: &mut Vec<String>,
) {
    let diffs = compare_molecule(bam1_run, bam2_run, h1, h2);
    if diffs.is_empty() {
        *matched += 1;
    } else {
        for d in diffs {
            push_diff(diff_details, max_diffs, || d.clone());
        }
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
/// small in practice — but it is never capped: an unbounded backlog (e.g. one file grouping
/// very differently from the other) is trusted and fully buffered rather than silently
/// dropped or windowed. Once both streams are exhausted, any molecule left in either
/// `pending` map was present in only one file and is reported as a diff naming its canonical
/// id.
///
/// The two per-side steps below are written out separately, rather than driven from a shared
/// array of `(iterator, own map, other map, ...)` tuples: `pending1`/`pending2` each need to
/// be borrowed mutably as *both* "own" (for their side) and "other" (for the opposite side)
/// within the same iteration, which a single array literal can't express without borrowing
/// one of the two maps mutably twice at once.
///
/// Header compatibility is not re-checked here: `require_compatible_headers` already gates
/// the CLI entry point (`CompareBams::execute`) before this driver ever runs, so by the time
/// `molecule_join_compare` is called the two headers are known compatible.
///
/// # Errors
///
/// Returns an error if either file cannot be opened/read, or if [`molecule_runs`] yields an
/// `Err` from either stream. On an `Err` from either side, the join stops immediately without
/// polling either iterator again: `molecule_runs` does not fuse after an `Err` (a stale
/// `pending` record inside it would otherwise silently leak into the next run's boundary), so
/// resuming past a first error would corrupt subsequent run boundaries.
///
/// Also returns an error if **neither** input ever saw an MI tag on any record. Grouping mode
/// compares *grouped* output (same-molecule reads consecutive under a shared MI), so a pair of
/// BAMs with zero MI tags between them is misuse of this comparison, not a legitimate "no
/// diffs" input — without this guard, two content-identical fully-MI-less BAMs would fold into
/// one giant MI-less "molecule" on each side and spuriously MATCH. A *partial* MI-less pair
/// (one side grouped, the other not) is unaffected by this guard: it is still caught downstream
/// as an ordinary molecule-count/membership asymmetry.
pub fn molecule_join_compare(
    bam1: &Path,
    bam2: &Path,
    max_diffs: usize,
) -> Result<MoleculeJoinOutcome> {
    let (_, h1) = create_raw_bam_reader(bam1, 1)?;
    let (_, h2) = create_raw_bam_reader(bam2, 1)?;
    // Header precondition already enforced at the CLI entry; header content is equal here.

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
    // A side "saw MI" iff any member of any run it produced carries an MI tag; see the
    // fully-MI-less guard below.
    let mut saw_mi1 = false;
    let mut saw_mi2 = false;

    // Alternate pulling one run from each side; probe the opposite pending map on each pull.
    // Because both inputs are template-coordinate ordered, corresponding molecules arrive at
    // ~the same offset, so `pending` stays small — but it is NOT capped (trust the data).
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
                    if run.members.iter().any(|m| get_mi_tag_raw(m).is_some()) {
                        saw_mi1 = true;
                    }
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
                    if run.members.iter().any(|m| get_mi_tag_raw(m).is_some()) {
                        saw_mi2 = true;
                    }
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
    }

    // Fully-MI-less guard (see this function's doc comment): if neither side ever saw an MI
    // tag, this is misuse of grouping comparison, not a legitimate "no diffs" verdict.
    if !saw_mi1 && !saw_mi2 {
        bail!(
            "grouping comparison requires MI-tagged (grouped) input; neither BAM contains any \
             MI tags"
        );
    }

    // Residual molecules present in only one file.
    for (canon, _) in pending1.drain() {
        push_diff(&mut diff_details, max_diffs, || {
            format!("molecule {} only in bam1", String::from_utf8_lossy(&canon))
        });
    }
    for (canon, _) in pending2.drain() {
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
    fn compare_molecule_cases(
        #[case] a: &[(&str, &str)],
        #[case] b: &[(&str, &str)],
        #[case] equal: bool,
    ) {
        let (ha, hb) = (minimal_header(), minimal_header());
        let ra = molecule_run_from(a);
        let rb = molecule_run_from(b);
        let diffs = compare_molecule(&ra, &rb, &ha, &hb);
        assert_eq!(diffs.is_empty(), equal, "diffs: {diffs:?}");
    }
}
