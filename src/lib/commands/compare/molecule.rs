//! Per-molecule runs for the streaming grouping comparison.
//!
//! `fgumi group`/`fgbio group` output is template-coordinate sorted with same-MI reads
//! consecutive (a molecule's reads sit together). [`molecule_runs`] cuts a raw BAM stream
//! into those per-molecule runs so the streaming molecule-join (Task 6) can compare two
//! grouped BAMs without materializing either input in memory.

use anyhow::Result;
use fgumi_raw_bam::RawRecord;
use fgumi_sort::RawBamRecordReader;
use std::io::Read;

use super::bams::get_mi_tag_raw;

/// One molecule's consecutive records, plus its MI-invariant canonical id.
pub(crate) struct MoleculeRun {
    /// Lexicographically smallest read name across `members` — a stable, MI-independent,
    /// collision-free molecule identity (a QNAME belongs to exactly one molecule).
    pub canon: Vec<u8>,
    /// The molecule's records (both duplex strands when present).
    pub members: Vec<RawRecord>,
}

impl MoleculeRun {
    /// Build a run from its member records, computing `canon` as the lexicographically
    /// smallest read name.
    fn from_members(members: Vec<RawRecord>) -> Self {
        let canon = members.iter().map(|r| r.read_name().to_vec()).min().unwrap_or_default();
        MoleculeRun { canon, members }
    }
}

/// Stream `reader` into per-molecule runs, cutting a run when the base MI (see
/// [`super::bams::MiKey::base`]) changes. Assumes same-molecule reads are consecutive
/// (guaranteed by grouped output; see the design doc). Records with no MI tag
/// (`base == None`) form their own runs, cut on any change to/from `None`.
pub(crate) fn molecule_runs<R: Read>(
    mut reader: RawBamRecordReader<R>,
) -> impl Iterator<Item = Result<MoleculeRun>> {
    let mut pending: Vec<RawRecord> = Vec::new();
    let mut pending_base: Option<i64> = None;
    std::iter::from_fn(move || {
        loop {
            match reader.next_record() {
                Ok(Some(rec)) => {
                    let base = get_mi_tag_raw(&rec).map(|mi| mi.base());
                    if pending.is_empty() {
                        pending_base = base;
                        pending.push(rec);
                    } else if base == pending_base {
                        pending.push(rec);
                    } else {
                        let done = std::mem::take(&mut pending);
                        pending_base = base;
                        pending.push(rec);
                        return Some(Ok(MoleculeRun::from_members(done)));
                    }
                }
                Ok(None) => {
                    if pending.is_empty() {
                        return None;
                    }
                    return Some(Ok(MoleculeRun::from_members(std::mem::take(&mut pending))));
                }
                Err(e) => return Some(Err(e.into())),
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::SamTag;
    use fgumi_raw_bam::{SamBuilder, flags};
    use noodles::sam::Header;
    use noodles::sam::alignment::io::Write as AlignmentWrite;

    fn rec(name: &[u8], mi: &str) -> RawRecord {
        // MI is a string-typed aux tag; `add_string_tag` is the chained setter
        // (crates/fgumi-raw-bam/src/builder.rs:431). Duplex strands use `"7/A"`/`"7/B"`.
        SamBuilder::new()
            .read_name(name)
            .flags(flags::FIRST_SEGMENT)
            .add_string_tag(SamTag::MI, mi.as_bytes())
            .build()
    }

    /// Writes `records` to a temp BAM (unmapped, so a default/empty header suffices),
    /// re-opens it via the raw reader used by the sort/compare engines, and drains
    /// [`molecule_runs`] into a `Vec`. Modeled on `read_all_records` in
    /// `tests/integration/test_compare_bams.rs:864`.
    fn collect_runs(records: Vec<RawRecord>) -> Vec<MoleculeRun> {
        let tmp = tempfile::NamedTempFile::new().expect("create temp BAM");
        let header = Header::default();
        {
            let mut writer = noodles::bam::io::Writer::new(
                std::fs::File::create(tmp.path()).expect("create test BAM"),
            );
            writer.write_header(&header).expect("write test header");
            for record in &records {
                let record_buf = fgumi_raw_bam::raw_record_to_record_buf(record, &header)
                    .expect("raw_record_to_record_buf should succeed in test");
                writer.write_alignment_record(&header, &record_buf).expect("write test record");
            }
            writer.try_finish().expect("finish test BAM");
        }

        let file = std::fs::File::open(tmp.path()).expect("open temp BAM");
        let mut reader = RawBamRecordReader::new(file).expect("open raw reader");
        reader.skip_header().expect("skip header");

        molecule_runs(reader).collect::<Result<Vec<_>>>().expect("molecule_runs should succeed")
    }

    /// A run is the maximal block of consecutive records sharing a base MI; canon is the
    /// smallest read name within it.
    #[test]
    fn cuts_runs_on_base_mi_change_and_canonicalizes_by_min_name() {
        let recs = vec![
            rec(b"r_charlie", "1"),
            rec(b"r_alpha", "1"), // molecule 1
            rec(b"r_bravo", "2"), // molecule 2
        ];
        let runs = collect_runs(recs);
        assert_eq!(runs.len(), 2);
        assert_eq!(runs[0].canon, b"r_alpha"); // min of {charlie, alpha}
        assert_eq!(runs[0].members.len(), 2);
        assert_eq!(runs[1].canon, b"r_bravo");
    }

    /// Duplex strands X/A and X/B share a base MI → one molecule run.
    #[test]
    fn duplex_strands_of_one_molecule_form_a_single_run() {
        let recs = vec![rec(b"r_a", "7/A"), rec(b"r_b", "7/B"), rec(b"r_c", "8/A")];
        let runs = collect_runs(recs);
        assert_eq!(runs.len(), 2, "7/A and 7/B are one molecule; 8/A is another");
        assert_eq!(runs[0].members.len(), 2);
    }
}
