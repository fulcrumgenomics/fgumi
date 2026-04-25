//! Determinism regression tests for the `group` and `dedup` commands.
//!
//! These tests guard against a class of non-determinism caused by iterating
//! `AHashMap`-keyed orientation subgroups (`(bool, bool) -> Vec<usize>`)
//! without sorting before assigning `MoleculeId`s. With Rust's per-process
//! random hasher, two runs of the same command on the same input would walk
//! orientation subgroups in different orders and produce identical molecule
//! groupings but mismatched `MI:Z` numbering across runs — which then
//! propagates into downstream consensus QNAMEs.
//!
//! The failing case requires multiple distinct orientations so that
//! `subgroups` has more than one entry; otherwise iteration order is trivially
//! stable.
//!
//! Both tests build a small paired-end BAM containing several UMI families
//! distributed across both FR (R1 forward, R2 reverse) and RF (R1 reverse,
//! R2 forward) orientations, run the command twice on the same input, and
//! assert that the per-record `MI:Z` tags are identical across runs.

use fgumi_lib::sam::SamTag;
use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
use noodles::bam;
use noodles::sam::alignment::io::Write as AlignmentWrite;
use std::fs;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::bam_generator::{create_minimal_header, to_record_buf};

/// Build a paired-end template for `--strategy paired` with explicit orientation.
///
/// `is_rf = false` produces FR orientation: R1 forward at `r1_pos`, R2 reverse
/// at `r2_pos`. `is_rf = true` flips the strand bits for both reads.
fn build_paired_template(
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

/// Write a BAM containing the given templates against a minimal `chr1` header.
fn write_bam(path: &Path, templates: &[(RawRecord, RawRecord)]) {
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
/// triggers the bug: the assigner sees more than one orientation subgroup,
/// and within each subgroup more than one UMI family, so `MoleculeId`
/// numbering depends on subgroup iteration order.
///
/// We materialize a large number of distinct molecules per subgroup so that
/// even small variances in `AHashMap::values()` iteration order become
/// observable as different `MoleculeId` numbering across runs.
fn build_mixed_orientation_bam(path: &Path) {
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

/// Read all `(QNAME, flags, MI:Z)` tuples from a BAM in input order.
///
/// Returning the tuple lets the test compare both per-output-position (which
/// catches sort-order differences induced by re-numbered MIs) and
/// per-(QNAME, flags) (which catches MI re-numbering on the same record).
fn read_qname_mi(path: &Path) -> Vec<(String, u16, String)> {
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

/// Run a fgumi subcommand `n_runs` times on the same input and assert that
/// the per-record MI tags are byte-identical across every pair of runs.
///
/// Each run is a fresh process invocation so `AHashMap`'s per-process random
/// hasher seed varies between runs; this is exactly the property that exposes
/// the orientation-subgroup iteration bug. Running more than two times and
/// comparing pairwise raises the probability of catching the variance even
/// when only a handful of subgroup keys are present.
fn assert_mi_deterministic(subcommand: &str, threads: &str, n_runs: usize) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    build_mixed_orientation_bam(&input_bam);

    let mut all_runs: Vec<Vec<(String, u16, String)>> = Vec::with_capacity(n_runs);
    for run in 0..n_runs {
        let out = temp_dir.path().join(format!("out_{run}.bam"));
        let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .args([
                subcommand,
                "--input",
                input_bam.to_str().unwrap(),
                "--output",
                out.to_str().unwrap(),
                "--strategy",
                "paired",
                "--threads",
                threads,
                "--compression-level",
                "1",
            ])
            .status()
            .expect("Failed to run command");
        assert!(status.success(), "{subcommand} run {run} failed");
        all_runs.push(read_qname_mi(&out));
    }

    let baseline = &all_runs[0];
    let baseline_by_key: std::collections::HashMap<(String, u16), &str> =
        baseline.iter().map(|(q, f, m)| ((q.clone(), *f), m.as_str())).collect();

    for (run, run_records) in all_runs.iter().enumerate().skip(1) {
        assert_eq!(
            baseline.len(),
            run_records.len(),
            "{subcommand} produced different record counts across runs",
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
            "{subcommand} --threads {threads}: per-(QNAME,flags) MI tags differ between run 0 and run {run} ({keyed_mismatches}/{} records). Positional mismatches: {pos_mismatches}. Examples: {sample_diffs:?}",
            baseline.len(),
        );
    }
}

/// Repeated runs of `fgumi group --strategy paired` must assign identical
/// `MI:Z` tags to every record.
#[test]
fn test_group_paired_mi_deterministic_threads_4() {
    assert_mi_deterministic("group", "4", 6);
}

/// Same as above for the single-threaded path.
#[test]
fn test_group_paired_mi_deterministic_threads_1() {
    assert_mi_deterministic("group", "1", 6);
}

/// Repeated runs of `fgumi dedup --strategy paired` must assign identical
/// `MI:Z` tags to every accepted record.
#[test]
fn test_dedup_paired_mi_deterministic_threads_4() {
    assert_mi_deterministic("dedup", "4", 6);
}
