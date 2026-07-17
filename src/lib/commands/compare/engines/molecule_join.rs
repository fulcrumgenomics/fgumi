//! Per-molecule comparison for the streaming grouping engine.
//!
//! [`compare_molecule`] compares two [`MoleculeRun`]s — the same molecule (matched by
//! canonical id) as read from each of the two BAMs under comparison — via three purely
//! local checks: record membership, per-record content (excluding MI), and the duplex
//! `/A`/`/B` strand partition. It has no notion of streaming or cross-molecule state;
//! that is the job of the join driver (Task 6) that pairs up `MoleculeRun`s by canonical
//! id and feeds each matched pair here.

use std::collections::{BTreeMap, BTreeSet};

use fgumi_raw_bam::RawRecord;
use noodles::sam::Header;

use super::super::bams::{MiKey, get_mi_tag_raw};
use super::super::molecule::MoleculeRun;
use super::super::record_key::{RecordKey, record_key};
use super::content::{ContentPredicate, content_diffs};

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
// `#[allow(dead_code)]`: not yet consumed outside tests — wired up by the join driver in
// Task 6.
#[allow(dead_code)]
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
