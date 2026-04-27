//! Shared fixture builders for MI-determinism regression tests.
//!
//! Provides the synthetic paired-end BAM fixture used by both
//! `test_group_determinism` and `test_runall_mi_determinism` to ensure the two
//! test suites exercise the same input data and that any future changes to the
//! fixture automatically apply to all callers.

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;

use crate::helpers::bam_generator::{create_minimal_header, to_record_buf};

/// Build a paired-end template for `--strategy paired` with explicit orientation.
///
/// `is_rf = false` produces FR orientation: R1 forward at `r1_pos`, R2 reverse
/// at `r2_pos`. `is_rf = true` flips the strand bits for both reads.
pub fn build_paired_template(
    name: &str,
    umi: &str,
    sequence: &str,
    quality: u8,
    r1_pos: i32,
    r2_pos: i32,
    is_rf: bool,
) -> (RawRecord, RawRecord) {
    let seq = sequence.as_bytes();
    let cigar_op = u32::try_from(seq.len()).expect("seq len fits u32") << 4;
    let cigar_str = format!("{}M", seq.len());
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_possible_wrap,
        reason = "small synthetic positions fit in i32"
    )]
    let tlen = (r2_pos - r1_pos + seq.len() as i32).abs();

    let r1_flags = flags::PAIRED
        | flags::FIRST_SEGMENT
        | if is_rf { flags::REVERSE } else { 0 }
        | if is_rf { 0 } else { flags::MATE_REVERSE };
    let r2_flags = flags::PAIRED
        | flags::LAST_SEGMENT
        | if is_rf { 0 } else { flags::REVERSE }
        | if is_rf { flags::MATE_REVERSE } else { 0 };

    let r1 = {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .ref_id(0)
            .pos(r1_pos - 1)
            .mapq(60)
            .flags(r1_flags)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(0)
            .mate_pos(r2_pos - 1)
            .template_length(if is_rf { -tlen } else { tlen })
            .sequence(seq)
            .qualities(&vec![quality; seq.len()])
            .add_string_tag(SamTag::RX, umi.as_bytes())
            .add_string_tag(SamTag::MC, cigar_str.as_bytes());
        b.build()
    };
    let r2 = {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .ref_id(0)
            .pos(r2_pos - 1)
            .mapq(60)
            .flags(r2_flags)
            .cigar_ops(&[cigar_op])
            .mate_ref_id(0)
            .mate_pos(r1_pos - 1)
            .template_length(if is_rf { tlen } else { -tlen })
            .sequence(seq)
            .qualities(&vec![quality; seq.len()])
            .add_string_tag(SamTag::RX, umi.as_bytes())
            .add_string_tag(SamTag::MC, cigar_str.as_bytes());
        b.build()
    };
    (r1, r2)
}

/// Read all `(QNAME, flags, MI:Z)` tuples from a BAM in input order.
///
/// Returning the tuple lets callers compare both per-output-position (which
/// catches sort-order differences induced by re-numbered MIs) and
/// per-(QNAME, flags) (which catches MI re-numbering on the same record).
pub fn read_qname_mi(path: &Path) -> Vec<(String, u16, String)> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).expect("open output BAM"));
    let _header = reader.read_header().expect("read header");
    let mi = SamTag::MI.to_noodles_tag();
    reader
        .records()
        .map(|r| {
            let r = r.expect("read record");
            let qname = r.name().expect("missing name").to_string();
            let flags = u16::from(r.flags());
            let mi_val = r.data().get(&mi).map(|v| format!("{v:?}")).unwrap_or_default();
            (qname, flags, mi_val)
        })
        .collect()
}

/// Assert that per-record `MI:Z` tags are identical across all provided runs.
///
/// `all_runs` is a slice of per-run record vectors, where each inner `Vec` was
/// produced by [`read_qname_mi`]. `total_runs` is the total number of runs
/// (used in assertion messages). `label` is interpolated into failure messages
/// to identify the command variant being tested (e.g. `"group --threads 4"` or
/// `"runall group --threads 4"`).
///
/// Performs two checks per run:
/// 1. **Keyed**: per-(QNAME, flags) MI tag equality — catches MI re-numbering
///    on the same logical record even if the output is re-sorted by MI.
/// 2. **Positional**: per-position equality — catches reordering induced by
///    different MI numbering across runs.
pub fn assert_mi_deterministic_across_runs(
    all_runs: &[Vec<(String, u16, String)>],
    total_runs: usize,
    label: &str,
) {
    let baseline = &all_runs[0];
    let baseline_by_key: HashMap<(String, u16), &str> =
        baseline.iter().map(|(q, f, m)| ((q.clone(), *f), m.as_str())).collect();

    for (run, run_records) in all_runs.iter().enumerate().skip(1) {
        assert_eq!(
            baseline.len(),
            run_records.len(),
            "{label}: produced different record counts across runs ({total_runs} total)",
        );

        // Per-(QNAME, flags) comparison: catches MI tag re-numbering on the
        // same logical record, even if the output is re-sorted by MI.
        let mut keyed_mismatches = 0usize;
        let mut sample_diffs: Vec<String> = Vec::new();
        for (qname, flags, mi_run) in run_records {
            let mi_base = baseline_by_key.get(&(qname.clone(), *flags));
            let Some(mi_base) = mi_base else {
                continue;
            };
            if *mi_base != mi_run.as_str() {
                keyed_mismatches += 1;
                if sample_diffs.len() < 5 {
                    sample_diffs
                        .push(format!("{qname} flags={flags}: run0={mi_base} run{run}={mi_run}",));
                }
            }
        }

        // Positional comparison: catches re-sorting that occurs because
        // templates are sorted by MI on the way out.
        let pos_mismatches =
            baseline.iter().zip(run_records.iter()).filter(|(x, y)| x != y).count();

        assert_eq!(
            keyed_mismatches,
            0,
            "{label}: per-(QNAME,flags) MI tags differ between run 0 and run {run} \
             ({keyed_mismatches}/{} records). Positional mismatches: {pos_mismatches}. \
             Examples: {sample_diffs:?}",
            baseline.len(),
        );
    }
}

/// Write a BAM containing the given templates against a minimal `chr1` header.
pub fn write_bam(path: &Path, templates: &[(RawRecord, RawRecord)]) {
    let header = create_minimal_header("chr1", 100_000);
    let mut writer =
        bam::io::Writer::new(fs::File::create(path).expect("Failed to create BAM file"));
    writer.write_header(&header).expect("Failed to write header");
    for (r1, r2) in templates {
        writer.write_alignment_record(&header, &to_record_buf(r1)).expect("Failed to write R1");
        writer.write_alignment_record(&header, &to_record_buf(r2)).expect("Failed to write R2");
    }
    writer.try_finish().expect("Failed to finish BAM");
}

/// Build a BAM with paired-end UMI families spread across multiple
/// orientations, multiple UMIs per orientation, and multiple positions.
///
/// The combination of distinct orientations and distinct paired UMIs is what
/// triggers the orientation-subgroup iteration-order bug: the assigner sees
/// more than one orientation subgroup, and within each subgroup more than one
/// UMI family, so `MoleculeId` numbering depends on subgroup iteration order.
///
/// We materialize a large number of distinct molecules per subgroup so that
/// even small variances in `AHashMap::values()` iteration order become
/// observable as different `MoleculeId` numbering across runs.
///
/// Produces: 40 position groups × 256 paired UMIs × 2 orientations
/// = 20480 templates = 40960 reads.
pub fn build_mixed_orientation_bam(path: &Path) {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut umis: Vec<String> = Vec::new();
    for &a in &bases {
        for &b in &bases {
            for &c in &bases {
                for &d in &bases {
                    let prefix = format!("{}{}{}{}", a as char, b as char, c as char, d as char,);
                    let suffix = format!("{}{}{}{}", d as char, c as char, b as char, a as char,);
                    umis.push(format!("{prefix}{prefix}-{suffix}{suffix}"));
                }
            }
        }
    }

    let mut templates = Vec::new();
    let mut counter = 0;
    let mut emit = |name: String, umi: &str, pos: i32, is_rf: bool| {
        templates.push(build_paired_template(
            &name,
            umi,
            "ACGTACGTACGTACGT",
            30,
            pos,
            pos + 100,
            is_rf,
        ));
    };

    // Many position groups, each carrying templates from BOTH orientations.
    // Position-grouping is per-coordinate, so distinct positions form
    // independent subgroup AHashMaps in the assigner — each one a separate
    // chance for iteration order to vary.
    for group_idx in 0..40 {
        let base_pos = 1_000 + group_idx * 1_000;
        for (u_idx, umi) in umis.iter().enumerate() {
            for is_rf in [false, true] {
                counter += 1;
                let name = format!("read_g{group_idx}_u{u_idx}_{is_rf}_{counter}");
                emit(name, umi, base_pos, is_rf);
            }
        }
    }
    write_bam(path, &templates);
}

/// Build a BAM with proper duplex pairs: every molecule contributes BOTH a
/// canonical-UMI strand and a swapped-UMI strand, so `--strategy paired`
/// assigns one of each to `PairedA(N)` and `PairedB(N)`.
///
/// The canonical strand uses UMI `"AAAA-CCCC"` in FR orientation; the
/// swapped strand uses UMI `"CCCC-AAAA"` in RF orientation, where the two
/// halves come from the same source molecule. The paired assigner detects
/// the reverse-pair relationship and stamps `PairedA(N)` on the canonical
/// strand and `PairedB(N)` on the reversed one. Downstream duplex consensus
/// then has both strands available and emits a duplex consensus record.
///
/// Produces: 40 position groups × 256 paired-UMI families × 2 strands × 2
/// supporting reads-per-strand = 40,960 templates / 81,920 reads.
pub fn build_duplex_pair_bam(path: &Path) {
    let bases = [b'A', b'C', b'G', b'T'];
    // Build 256 canonical UMI prefixes; suffix is reversed prefix.
    let mut canonical_halves: Vec<(String, String)> = Vec::new();
    for &a in &bases {
        for &b in &bases {
            for &c in &bases {
                for &d in &bases {
                    let prefix = format!("{}{}{}{}", a as char, b as char, c as char, d as char);
                    let suffix = format!("{}{}{}{}", d as char, c as char, b as char, a as char);
                    canonical_halves
                        .push((format!("{prefix}{prefix}"), format!("{suffix}{suffix}")));
                }
            }
        }
    }

    let mut templates = Vec::new();
    let mut counter: u64 = 0;

    // For each duplex pair: the top strand's R1 is at the low coord with the
    // canonical UMI; the bottom strand's R1 is at the HIGH coord with the
    // swapped UMI. The position swap flips `is_r1_genomically_earlier`, which
    // in turn flips the lower/higher prefix order in `umi_for_read_impl`. The
    // resulting prefixed UMI on the bottom strand is the lexical reverse of
    // the prefixed UMI on the top strand, which is exactly what
    // `PairedUmiAssigner` requires to detect a duplex pair and stamp
    // `PairedA(N)` / `PairedB(N)`.
    for group_idx in 0..40 {
        let low_pos = 1_000 + group_idx * 1_000;
        let high_pos = low_pos + 100;
        for (u_idx, (lhs, rhs)) in canonical_halves.iter().enumerate() {
            // Two reads per strand so the duplex caller's min-reads check
            // (default 1) clears comfortably and a few reads can be lost to
            // any per-read masking without dropping the whole pair.
            for rep in 0..2 {
                let canonical = format!("{lhs}-{rhs}");
                let swapped = format!("{rhs}-{lhs}");
                counter += 1;
                templates.push(build_paired_template(
                    &format!("read_g{group_idx}_u{u_idx}_top_r{rep}_{counter}"),
                    &canonical,
                    "ACGTACGTACGTACGT",
                    30,
                    low_pos,
                    high_pos,
                    false,
                ));
                counter += 1;
                templates.push(build_paired_template(
                    &format!("read_g{group_idx}_u{u_idx}_bot_r{rep}_{counter}"),
                    &swapped,
                    "ACGTACGTACGTACGT",
                    30,
                    high_pos,
                    low_pos,
                    true,
                ));
            }
        }
    }
    write_bam(path, &templates);
}
