//! Per-molecule runs for the streaming grouping comparison.
//!
//! `fgumi group`/`fgbio group` output is template-coordinate sorted with same-MI reads
//! consecutive (a molecule's reads sit together). [`molecule_runs`] cuts a raw BAM stream
//! into those per-molecule runs so the streaming molecule-join (Task 6) can compare two
//! grouped BAMs without materializing either input in memory.

use anyhow::{Result, anyhow};
use fgumi_raw_bam::RawRecord;
use fgumi_sort::RawBamRecordReader;
use std::io::Read;

use super::bams::get_mi_tag_raw;

/// One molecule's consecutive records, plus its MI-invariant canonical id.
#[derive(Debug)]
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
/// (guaranteed by grouped output; see the design doc) and *enforces* two grouping invariants
/// in O(1) memory, yielding an `Err` on the first record that violates either:
///
/// 1. **Every record is MI-tagged.** `fgumi group`/`fgbio GroupReadsByUmi` tag every emitted
///    read with an MI, so a record with a missing or unparseable MI base (`base == None`)
///    means the input is not grouped. Such a record is rejected immediately rather than folded
///    into an MI-less run: because the molecule-join compares content *excluding* MI, an
///    MI-less run could otherwise still canonical-id/content-match its tagged counterpart,
///    silently passing a molecule that lost its MI (a false MATCH). Rejecting per-record is
///    strictly stronger than any whole-input "never saw an MI" guard, which a *partial*-MI-less
///    input (some records tagged, some not) slips through. An empty input yields no runs and no
///    error.
/// 2. **Base MI is strictly increasing across runs.** Both tools assign the MI as a
///    monotonically increasing counter (fgbio: `counter.getAndIncrement()`; fgumi:
///    `next_mi_base`) in template-coordinate scan order, so a valid grouped file's base MI is
///    *strictly increasing* across runs. A new run whose base is not greater than a base
///    already seen is either a reappearance (non-consecutive) or a non-monotonic id — i.e. not
///    grouped — and is rejected rather than silently starting a second, disconnected run (which
///    would corrupt the molecule-join's canonical-id matching downstream). Tracking only the
///    largest base seen makes this O(1), not O(molecules); the engine never materializes either
///    input's records.
/// 3. **A single run is bounded.** A run accumulates its members in `pending` before being
///    yielded, so a degenerate input where one MI spans (nearly) the whole file would grow
///    `pending` to O(file size) *before* the run is ever handed to the join — which the
///    downstream molecule-backlog cap ([`super::engines::molecule_join::MAX_PENDING_MOLECULES`])
///    cannot prevent, since that bounds the number of *unmatched molecules*, not the size of one.
///    A deliberately high per-run record ceiling ([`MAX_MOLECULE_RUN_RECORDS`]) turns that
///    pathological O(file) accumulation into a clear error instead of an eventual OOM; a
///    well-formed grouped input has bounded UMI families and never approaches it.
pub(crate) fn molecule_runs<R: Read>(
    reader: RawBamRecordReader<R>,
) -> impl Iterator<Item = Result<MoleculeRun>> {
    molecule_runs_capped(reader, MAX_MOLECULE_RUN_RECORDS)
}

/// Deliberately high last-resort ceiling on the number of records buffered for a *single*
/// molecule run before it is yielded (see invariant 3 in [`molecule_runs`]'s doc comment).
/// Not a tuning knob: legitimate UMI families are small (a handful to low thousands of reads),
/// so this only fires on degenerate ungrouped input (e.g. one MI spanning the file), converting
/// an O(file-size) buffer from an OOM into a clear error.
pub(crate) const MAX_MOLECULE_RUN_RECORDS: usize = 100_000_000;

/// [`molecule_runs`] with an explicit per-run record ceiling, so tests can exercise the
/// over-cap error path without materializing 100M records. Production always calls through the
/// public wrapper with [`MAX_MOLECULE_RUN_RECORDS`].
pub(crate) fn molecule_runs_capped<R: Read>(
    mut reader: RawBamRecordReader<R>,
    max_run_records: usize,
) -> impl Iterator<Item = Result<MoleculeRun>> {
    let mut pending: Vec<RawRecord> = Vec::new();
    let mut pending_base: Option<i64> = None;
    // The largest base MI of any run started so far — O(1), not an O(molecules) set (see
    // invariant 2 in the doc comment).
    let mut last_base: Option<i64> = None;
    std::iter::from_fn(move || {
        loop {
            match reader.next_record() {
                Ok(Some(rec)) => {
                    // Invariant 1: every record in a non-empty grouped input carries a
                    // parseable MI. `get_mi_tag_raw` returns `None` for both a missing tag and
                    // an unparseable one; reject the first such record (stricter than a
                    // whole-input guard, which a partial-MI-less side would slip through).
                    let base = match get_mi_tag_raw(&rec).map(|mi| mi.base()) {
                        Some(b) => b,
                        None => {
                            return Some(Err(anyhow!(
                                "input is not grouped: encountered a record with no (or an \
                                 unparseable) MI tag; grouping comparison requires every record \
                                 in a non-empty input to be MI-tagged (grouped output tags every \
                                 read)"
                            )));
                        }
                    };
                    if !pending.is_empty() && Some(base) == pending_base {
                        // Invariant 3: bound a single run's accumulation. `pending` already
                        // holds this molecule's members; refuse to grow it past the ceiling
                        // rather than buffering O(file-size) records for a degenerate one-MI
                        // input.
                        if pending.len() >= max_run_records {
                            return Some(Err(anyhow!(
                                "input is not grouped: a single molecule (MI base {base}) exceeds \
                                 {max_run_records} buffered records before yielding; a well-formed \
                                 grouped input has bounded UMI families, so this indicates \
                                 ungrouped input (e.g. one MI spanning the file) that would \
                                 otherwise buffer O(file-size) records and OOM"
                            )));
                        }
                        pending.push(rec);
                        continue;
                    }
                    // Starting a new run — either the very first record of the stream, or a
                    // cut from the previous base. Invariant 2: enforce strictly-increasing MI
                    // base.
                    if let Some(m) = last_base
                        && base <= m
                    {
                        return Some(Err(anyhow!(
                            "input is not grouped: MI base {base} is not greater than a \
                             previously seen base {m}; grouping comparison requires \
                             same-MI-consecutive (grouped) input with monotonically \
                             increasing molecule ids"
                        )));
                    }
                    last_base = Some(base);
                    if pending.is_empty() {
                        pending_base = Some(base);
                        pending.push(rec);
                    } else {
                        let done = std::mem::take(&mut pending);
                        pending_base = Some(base);
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
    use rstest::rstest;

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
    fn try_collect_runs(records: Vec<RawRecord>) -> Result<Vec<MoleculeRun>> {
        try_collect_runs_capped(records, MAX_MOLECULE_RUN_RECORDS)
    }

    fn try_collect_runs_capped(
        records: Vec<RawRecord>,
        max_run_records: usize,
    ) -> Result<Vec<MoleculeRun>> {
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

        molecule_runs_capped(reader, max_run_records).collect::<Result<Vec<_>>>()
    }

    fn collect_runs(records: Vec<RawRecord>) -> Vec<MoleculeRun> {
        try_collect_runs(records).expect("molecule_runs should succeed")
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

    /// A non-consecutive base MI (base 1's run is cut by base 2, then base 1 reappears)
    /// violates the same-MI-consecutive contract `molecule_runs` assumes; it must surface
    /// as an `Err`, not silently mis-group the reappearing base into a second run.
    #[test]
    fn base_mi_reappearing_after_run_closes_is_rejected() {
        let recs = vec![rec(b"r1", "1"), rec(b"r2", "1"), rec(b"r3", "2"), rec(b"r4", "1")];
        let err = try_collect_runs(recs)
            .expect_err("scattered (non-consecutive) MI base must be rejected");
        let msg = err.to_string();
        assert!(msg.contains("not grouped"), "got: {msg}");
        assert!(msg.contains('1'), "error should name the reappearing base: {msg}");
    }

    /// The valid counterpart: each base MI's records stay consecutive, so no error.
    #[test]
    fn base_mi_consecutive_runs_are_accepted() {
        let recs = vec![rec(b"r1", "1"), rec(b"r2", "1"), rec(b"r3", "2"), rec(b"r4", "2")];
        let runs = collect_runs(recs);
        assert_eq!(runs.len(), 2);
    }

    /// Contiguous but *non-monotonic* base MIs (`1`, then `3`, then `2`) are not valid grouped
    /// output — both grouping tools assign the MI as a strictly increasing counter — so the
    /// O(1) monotonicity check must reject it even though no base literally reappears (a plain
    /// "seen this base before" set would have missed the `3`→`2` regression).
    #[test]
    fn non_monotonic_base_mi_is_rejected() {
        let recs = vec![rec(b"r1", "1"), rec(b"r2", "3"), rec(b"r3", "2")];
        let err = try_collect_runs(recs).expect_err("non-monotonic MI base must be rejected");
        let msg = err.to_string();
        assert!(msg.contains("not grouped"), "got: {msg}");
        assert!(msg.contains('2') && msg.contains('3'), "error should name the bases: {msg}");
    }

    /// Invariant 1: a record with no MI tag in a non-empty input is rejected on the spot,
    /// rather than folded into an MI-less run (which — since the molecule-join excludes MI from
    /// content comparison — could otherwise let a molecule that lost its MI still MATCH its
    /// tagged counterpart). The reject fires on the *first* such record, so it applies to a
    /// *mixed* input (some records tagged, some not), not only a wholly-MI-less one — the case
    /// a whole-input "never saw MI" guard would slip through. The message names both "not
    /// grouped" and "MI-tagged" so the engine- and CLI-level tests can key off either.
    fn mi_less() -> RawRecord {
        SamBuilder::new().read_name(b"r_none").flags(flags::FIRST_SEGMENT).build()
    }

    #[rstest]
    #[case::leading_mi_less(vec![mi_less(), rec(b"r2", "1")])]
    #[case::trailing_mi_less(vec![rec(b"r1", "1"), mi_less()])]
    #[case::mixed_mid_stream(vec![rec(b"r1", "1"), rec(b"r2", "1"), mi_less(), rec(b"r4", "2")])]
    fn mi_less_record_in_non_empty_input_is_rejected(#[case] recs: Vec<RawRecord>) {
        let err = try_collect_runs(recs).expect_err("a record with no MI tag must be rejected");
        let msg = err.to_string();
        assert!(msg.contains("not grouped"), "got: {msg}");
        assert!(msg.contains("MI-tagged"), "got: {msg}");
    }

    /// Invariant 3: a single molecule that would buffer past the per-run ceiling is rejected
    /// rather than growing `pending` unboundedly — the degenerate "one MI spans the file" case
    /// that the downstream molecule-backlog cap cannot catch. A tiny cap makes the trip cheap;
    /// production uses [`MAX_MOLECULE_RUN_RECORDS`].
    #[test]
    fn single_molecule_over_run_cap_is_rejected() {
        // Three records sharing one MI base; a cap of 2 bails when the third would push.
        let recs = vec![rec(b"r_a", "1"), rec(b"r_b", "1"), rec(b"r_c", "1")];
        let err = try_collect_runs_capped(recs, 2)
            .expect_err("a single molecule exceeding the per-run cap must be rejected");
        let msg = err.to_string();
        assert!(msg.contains("not grouped"), "got: {msg}");
        assert!(msg.contains("single molecule"), "got: {msg}");
    }

    /// The boundary counterpart: a run of exactly the ceiling is accepted (the cap bounds the
    /// buffer, it does not reject a legitimately cap-sized family).
    #[test]
    fn single_molecule_at_run_cap_is_accepted() {
        let recs = vec![rec(b"r_a", "1"), rec(b"r_b", "1")];
        let runs = try_collect_runs_capped(recs, 2).expect("a run at the cap must be accepted");
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].members.len(), 2);
    }
}
