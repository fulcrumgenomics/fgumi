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
use noodles::sam::alignment::record_buf::data::field::Value;
use rstest::rstest;
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
///
/// The MI tag is extracted via the typed `Value::String` accessor rather
/// than the `Debug` representation so the assertion compares the actual MI
/// string rather than coupling to `noodles`' `Value` Debug surface.
fn read_qname_mi(path: &Path) -> Vec<(String, u16, String)> {
    let mut reader = bam::io::Reader::new(fs::File::open(path).expect("open output BAM"));
    let header = reader.read_header().expect("read header");
    let mi = SamTag::MI.to_noodles_tag();
    reader
        .record_bufs(&header)
        .map(|r| {
            let r = r.expect("read record");
            let qname = r.name().expect("missing name").to_string();
            let flags = u16::from(r.flags());
            let mi_val = match r.data().get(&mi) {
                Some(Value::String(s)) => s.to_string(),
                Some(other) => panic!("expected MI tag to be a string, got {other:?}"),
                None => panic!("missing MI tag on output record {qname}"),
            };
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
///
/// `threads = None` omits the `--threads` flag entirely, which resolves to a
/// one-worker chain (`num_threads() == 1`); `Some("1")` requests one worker
/// explicitly. Both route through the same typed-step pipeline and must be
/// deterministic (and, per `group_no_threads_matches_threads_1`, identical).
fn assert_mi_deterministic(
    subcommand: &str,
    threads: Option<&str>,
    extra_args: &[&str],
    n_runs: usize,
) {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    build_mixed_orientation_bam(&input_bam);

    let mut all_runs: Vec<Vec<(String, u16, String)>> = Vec::with_capacity(n_runs);
    for run in 0..n_runs {
        let out = temp_dir.path().join(format!("out_{run}.bam"));
        let mut args: Vec<&str> = vec![
            subcommand,
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--strategy",
            "paired",
            "--compression-level",
            "1",
        ];
        if let Some(t) = threads {
            args.extend(["--threads", t]);
        }
        args.extend_from_slice(extra_args);
        let status = Command::new(env!("CARGO_BIN_EXE_fgumi"))
            .args(&args)
            .status()
            .expect("Failed to run command");
        assert!(status.success(), "{subcommand} run {run} failed");
        all_runs.push(read_qname_mi(&out));
    }

    let threads_label = threads.unwrap_or("<none>");
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
            pos_mismatches, 0,
            "{subcommand} --threads {threads_label}: output order differs between run 0 and run {run} ({pos_mismatches} positional mismatches)",
        );

        assert_eq!(
            keyed_mismatches,
            0,
            "{subcommand} --threads {threads_label}: per-(QNAME,flags) MI tags differ between run 0 and run {run} ({keyed_mismatches}/{} records). Positional mismatches: {pos_mismatches}. Examples: {sample_diffs:?}",
            baseline.len(),
        );
    }
}

/// Repeated runs of `fgumi group` and `fgumi dedup` with `--strategy paired`
/// must assign identical `MI:Z` tags to every record across runs and across
/// every supported threading mode.
///
/// The three threads cases exercise the pipeline at different worker counts:
///
/// * `None` — `--threads` omitted, which resolves to a one-worker chain.
/// * `Some("1")` — the unified pipeline with a single worker.
/// * `Some("4")` — the unified pipeline with four workers, which was the
///   originally failing case before the serial-ordered `MI Assign` hook.
#[rstest]
#[case::group_no_threads("group", None)]
#[case::group_threads_1("group", Some("1"))]
#[case::group_threads_4("group", Some("4"))]
#[case::dedup_threads_1("dedup", Some("1"))]
#[case::dedup_threads_4("dedup", Some("4"))]
fn test_paired_mi_deterministic(#[case] subcommand: &str, #[case] threads: Option<&str>) {
    assert_mi_deterministic(subcommand, threads, &[], 6);
}

/// B1 (audit): `fgumi group` with `--threads` unset must produce identical
/// grouping to `--threads 1` — same record order and same per-record `MI:Z`
/// tags. Historically the unset case took a bespoke `execute_single_threaded`
/// orchestration while `--threads 1` routed through the typed-step chain; this
/// pins that the two paths agree, which is the invariant that lets the
/// single-threaded path be folded onto the chain (deleting the second
/// orchestration). Uses the paired-orientation fixture — the case most sensitive
/// to MI-numbering drift.
///
/// Beyond the unset-vs-`--threads 1` self-consistency check, this test also
/// asserts the unset output against an INDEPENDENT oracle derived from the
/// fixture construction (`assert_paired_fixture_oracle`), per the repo
/// path-instruction that correctness-critical integration tests validate
/// against an expectation not taken from the tool's own output — so that a
/// grouping that is identical across the two paths but both wrong cannot pass.
#[test]
fn group_no_threads_matches_threads_1() {
    let temp_dir = TempDir::new().unwrap();
    let input_bam = temp_dir.path().join("input.bam");
    build_mixed_orientation_bam(&input_bam);

    let run = |out: &Path, threads: Option<&str>| {
        let mut args: Vec<&str> = vec![
            "group",
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--strategy",
            "paired",
            "--compression-level",
            "1",
        ];
        if let Some(t) = threads {
            args.extend(["--threads", t]);
        }
        let status =
            Command::new(env!("CARGO_BIN_EXE_fgumi")).args(&args).status().expect("run group");
        assert!(status.success(), "group (threads={threads:?}) failed");
    };

    let out_none = temp_dir.path().join("out_none.bam");
    let out_t1 = temp_dir.path().join("out_t1.bam");
    run(&out_none, None);
    run(&out_t1, Some("1"));

    let mi_none = read_qname_mi(&out_none);

    // Independent oracle (satisfies the path-instruction that correctness-critical
    // integration tests check against an expectation NOT derived from the tool's
    // own output, so "identical-but-both-wrong" grouping cannot pass): assert the
    // `--threads` unset output — the path PR 547 folds onto the chain — against
    // the MI-group structure derived purely from the fixture construction.
    assert_paired_fixture_oracle(&mi_none);

    // Self-consistency: the folded unset path must match `--threads 1` exactly.
    assert_eq!(
        mi_none,
        read_qname_mi(&out_t1),
        "group --threads unset must match --threads 1 (record order + MI tags)"
    );
}

/// Independent oracle for `build_mixed_orientation_bam` under `--strategy paired`.
///
/// Derived solely from the fixture construction (NOT the tool's output): the
/// fixture emits, at each of 40 positions, for each of 256 distinct paired UMIs
/// (4^4 prefixes; `{prefix}{prefix}-{suffix}{suffix}` with suffix = reversed
/// prefix), one template in each of the two orientations (FR = `is_rf` false,
/// RF = `is_rf` true).
///
/// Under the `paired` strategy the two orientations at the *same* position are
/// distinct orientation subgroups, NOT the two strands of one duplex molecule:
/// a duplex pair needs opposite-strand reads at *mirrored* coordinates, whereas
/// here both orientations share identical coordinates. So each (position, UMI,
/// orientation) is its own molecule and the orientations do NOT merge:
///
///   families = 40 positions × 256 UMIs × 2 orientations = 20480
///
/// (If they had merged it would be 40 × 256 = 10240; the empirical run confirms
/// 20480, i.e. the 40 × 512 branch.) Each template is 2 reads and no two
/// templates share a family, so every family holds exactly 2 records; fgumi
/// numbers molecule IDs contiguously with no gaps, so the distinct numeric MI
/// values form an exact gap-free range of length `families`.
fn assert_paired_fixture_oracle(records: &[(String, u16, String)]) {
    const POSITIONS: u64 = 40;
    const UMIS: u64 = 256; // 4^4 distinct 4-mer prefixes
    const ORIENTATIONS: u64 = 2; // FR + RF stay distinct under `paired`
    const EXPECTED_FAMILIES: u64 = POSITIONS * UMIS * ORIENTATIONS; // 20480
    const RECORDS_PER_FAMILY: usize = 2; // one 2-read template per family

    // Family id per record: the numeric part of `MI:Z`, stripping any
    // `/A`,`/B` strand suffix.
    let mut family_sizes: std::collections::HashMap<u64, usize> = std::collections::HashMap::new();
    for (qname, _flags, mi) in records {
        let numeric = mi.split('/').next().expect("non-empty MI");
        let id: u64 = numeric
            .parse()
            .unwrap_or_else(|_| panic!("MI {mi:?} on {qname} is not <n> or <n>/strand"));
        *family_sizes.entry(id).or_default() += 1;
    }

    // Family count matches the fixture-derived value (orientations do not merge).
    assert_eq!(
        family_sizes.len() as u64,
        EXPECTED_FAMILIES,
        "expected 40×256×2 = {EXPECTED_FAMILIES} MI families, got {}",
        family_sizes.len(),
    );

    // Total records = families × 2.
    assert_eq!(
        records.len() as u64,
        EXPECTED_FAMILIES * RECORDS_PER_FAMILY as u64,
        "expected {EXPECTED_FAMILIES} families × {RECORDS_PER_FAMILY} records each",
    );

    // Every family is uniformly sized (exactly one 2-read template).
    for (id, size) in &family_sizes {
        assert_eq!(
            *size, RECORDS_PER_FAMILY,
            "family {id} has {size} records, expected {RECORDS_PER_FAMILY}",
        );
    }

    // MI ids are numbered contiguously with no gaps: the distinct ids equal the
    // exact range [min, min + families).
    let min = family_sizes.keys().copied().min().expect("at least one family");
    let expected_ids: std::collections::HashSet<u64> = (min..min + EXPECTED_FAMILIES).collect();
    let actual_ids: std::collections::HashSet<u64> = family_sizes.keys().copied().collect();
    assert_eq!(
        actual_ids, expected_ids,
        "MI ids must be a contiguous, gap-free range of {EXPECTED_FAMILIES} values starting at {min}",
    );
}

/// `--parallel-group-min-templates` must route the chain `GroupProcess` step
/// through its per-group assigner-selection branch: when the threshold is set
/// there is no batch-wide assigner, so every group decides parallel-vs-sequential
/// by its own template count and builds its own assigner. At `--threads 1` the
/// `num_threads == 1` gate in `create_umi_assigner` takes the sequential
/// fallback, so this exercises that branch without paying per-group rayon-pool
/// construction (which the auto thresholds of 128–3072 exist to avoid; forcing
/// the threshold to 1 at high thread counts would spin up a pool per tiny group).
/// The assignment must stay deterministic across runs — the same lockstep
/// guarantee the non-threshold path provides.
#[test]
fn test_parallel_group_min_templates_chain_is_deterministic() {
    assert_mi_deterministic("group", Some("1"), &["--parallel-group-min-templates", "1"], 4);
}
