//! Parity gate for the `correct` command's legacy-path retirement (C4):
//! `CorrectUmis::execute` no longer has a no-`--threads` serial engine
//! (`execute_single_thread_mode`) or a `#[allow(dead_code)]`
//! `execute_threads_mode` — it *always* routes through the declarative chain
//! builder. This test proves that cutover lost nothing user-observable.
//!
//! Two independent things are checked here:
//!
//! 1. **The cutover actually happened** (`correct_no_threads_routes_through_chain`).
//!    A no-`--threads` run now emits the chain-only `"Using pipeline with N
//!    threads"` banner that `ChainBuilder::add_correct` logs and the retired
//!    serial engine never did. This is the genuine RED/GREEN discriminator:
//!    before the removal a no-`--threads` run took the serial path and printed no
//!    such line; after it, the chain does.
//!
//! 2. **Output parity with the pre-removal serial path**
//!    (`cutover_matches_baseline`, `cutover_min_corrected_failure_writes_metrics`).
//!    The current build's `correct` output — records (byte-identical, modulo the
//!    `@PG` line) and the `--metrics` TSV — must match the frozen serial baseline
//!    binary run with no `--threads` (its legacy single-threaded engine). The
//!    baseline path comes from `FGUMI_BASELINE_BIN`; when it is unset (or names a
//!    missing file) the case degrades to a self-consistency oracle rather than
//!    skipping — the exact fallback discipline of `test_sort_cutover_parity.rs`.
//!
//! The `--min-corrected`-failure case pins the fgbio-parity behavior the cutover
//! must preserve: a kept ratio below `--min-corrected` fails the run, but the
//! `--metrics` TSV is still written (fgbio CorrectUmis.scala:152 — "a ratio
//! below --min-corrected will cause a failure but all files will still be
//! written"). The retired legacy path wrote metrics before its ratio bail; the
//! chain matches it (metrics on a success-only hook, the gate on a later
//! success-only hook), and this test proves both the current build and the
//! baseline's legacy path error *and* leave identical metrics behind.
//!
//! **Not a RED/GREEN gate for the parity half.** Because the removed serial
//! engine and the chain were already output-equivalent, the baseline byte-parity
//! check passes on both sides of the change — like the sort/retag cutovers it
//! guards equivalence, it does not observe a regression the cutover introduces.
//! The chain-banner check in (1) is the part that flips RED→GREEN.

use std::ffi::OsStr;
use std::path::Path;
use std::process::Command;

use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::bam_generator::create_umi_family;
use crate::helpers::read_bam_output;
use crate::test_correct_command::create_umi_bam;

/// The fixed UMI set every case corrects against.
const WHITELIST: &[&str] = &["AAAAAAAA", "CCCCCCCC", "GGGGGGGG"];

/// Resolves the saved pre-removal serial baseline binary to compare against.
///
/// The path comes solely from `FGUMI_BASELINE_BIN`; when unset (or naming a
/// missing file) this returns `None` — "no baseline oracle available", never a
/// silent pass. Callers layer the baseline byte-parity check on top of the
/// always-available self-consistency oracle; a missing baseline drops only the
/// byte-parity half, it never skips the case. No hardcoded fallback: a baseline
/// binary is host-specific and must never be a path committed into the repo.
fn baseline_bin() -> Option<std::path::PathBuf> {
    let path = std::path::PathBuf::from(std::env::var_os("FGUMI_BASELINE_BIN")?);
    if path.is_file() {
        return Some(path);
    }
    eprintln!(
        "FGUMI_BASELINE_BIN={} does not name an existing file; baseline oracle unavailable",
        path.display()
    );
    None
}

/// Removes every `@PG` line from a SAM header text blob.
///
/// Both the current build and the baseline binary stamp a single `@PG` line
/// whose `VN` (git-describe version) and `CL` (command line, naming argv[0] — a
/// different binary path for each side — and the per-run output path) necessarily
/// differ between two independent invocations. Stripping the whole line is
/// correct here because neither side emits more than one `@PG` for these inputs
/// (the test BAMs carry no pre-existing `@PG`).
fn strip_pg_lines(text: &str) -> String {
    if text.is_empty() {
        return String::new();
    }
    let mut lines: Vec<&str> = text.split('\n').collect();
    let had_trailing_newline = lines.last() == Some(&"");
    if had_trailing_newline {
        lines.pop();
    }
    lines.retain(|line| !line.starts_with("@PG"));
    let mut out = lines.join("\n");
    if had_trailing_newline {
        out.push('\n');
    }
    out
}

/// Reads a BAM's raw BGZF stream, decompresses it, strips the `@PG` line(s) from
/// the embedded SAM header text, and returns the resulting bytes (new header +
/// unmodified `n_ref`/reference-list/record bytes).
///
/// Comparing the *decompressed* BAM binary — rather than `RecordBuf`-parsed
/// records or raw file bytes — keeps everything except the `@PG` line an exact,
/// uninterpreted byte comparison: BGZF block boundaries differ between the two
/// writers even for identical logical content (so raw file bytes never match),
/// while parsing into `RecordBuf` and re-encoding could mask a real tag-order or
/// binary-layout regression by normalizing it away. correct rewrites UMI/OX tag
/// bytes in place, so that masking risk is exactly what must be avoided.
fn decompressed_records_without_pg(path: &Path) -> Vec<u8> {
    let file = std::fs::File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    let mut reader = noodles::bgzf::io::Reader::new(std::io::BufReader::new(file));
    let mut raw = Vec::new();
    std::io::Read::read_to_end(&mut reader, &mut raw)
        .unwrap_or_else(|e| panic!("decompress BGZF stream for {}: {e}", path.display()));

    assert!(
        raw.len() >= 8 && &raw[0..4] == b"BAM\x01",
        "{} does not decompress to a BAM binary stream (missing magic)",
        path.display()
    );
    let l_text = i32::from_le_bytes(raw[4..8].try_into().expect("4 bytes"));
    let l_text = usize::try_from(l_text).expect("l_text is non-negative");
    let text_start = 8;
    let text_end = text_start + l_text;
    assert!(raw.len() >= text_end, "{} header text runs past end of stream", path.display());

    let text = String::from_utf8_lossy(&raw[text_start..text_end]);
    let stripped_text = strip_pg_lines(&text);

    let mut out = Vec::with_capacity(raw.len());
    out.extend_from_slice(b"BAM\x01");
    let new_l_text = i32::try_from(stripped_text.len()).expect("stripped header text fits i32");
    out.extend_from_slice(&new_l_text.to_le_bytes());
    out.extend_from_slice(stripped_text.as_bytes());
    out.extend_from_slice(&raw[text_end..]); // n_ref, reference list, and all records, untouched
    out
}

/// Writes a varied UMI-tagged corpus to `dir/in.bam` and returns its path.
///
/// Each family carries a distinct query-name prefix (so template grouping keeps
/// input order) and one of four RX shapes, exercising the exact-match,
/// single-mismatch-correction, and no-match-reject branches the chain must
/// reproduce (it does not cover the missing-UMI or wrong-length buckets — every
/// record here carries a well-formed 8-base RX):
/// - `AAAAAAAA` — exact match (kept, no OX rewrite)
/// - `AAAAAAAT` — 1 mismatch from `AAAAAAAA` (corrected, OX stashed)
/// - `CCCCCCCC` — exact match to a different whitelist UMI (kept)
/// - `TTTTTTTT` — matches nothing within `--max-mismatches 2` (rejected)
fn write_varied_input(dir: &Path) -> std::path::PathBuf {
    let input = dir.join("in.bam");
    let families = vec![
        create_umi_family("AAAAAAAA", 3, "exact_a", "ACGTACGT", 30),
        create_umi_family("AAAAAAAT", 2, "corr_a", "ACGTACGT", 30),
        create_umi_family("CCCCCCCC", 3, "exact_c", "ACGTACGT", 30),
        create_umi_family("TTTTTTTT", 2, "reject_t", "ACGTACGT", 30),
        create_umi_family("GGGGGGGG", 1, "exact_g", "ACGTACGT", 30),
    ];
    create_umi_bam(&input, families);
    input
}

/// Runs `<bin> correct` (no `--threads`, so the current build takes the
/// post-cutover chain path and the baseline takes its legacy serial path) with
/// `RUST_LOG=info`, and returns the process output for stderr/exit assertions.
fn run_correct(
    bin: &Path,
    input: &Path,
    output: &Path,
    metrics: Option<&Path>,
    min_corrected: Option<&str>,
) -> std::process::Output {
    let mut cmd = Command::new(bin);
    cmd.env("RUST_LOG", "info").args([
        OsStr::new("correct"),
        OsStr::new("-i"),
        input.as_os_str(),
        OsStr::new("-o"),
        output.as_os_str(),
        OsStr::new("--max-mismatches"),
        OsStr::new("2"),
        OsStr::new("--min-distance"),
        OsStr::new("1"),
        OsStr::new("--compression-level"),
        OsStr::new("1"),
    ]);
    for umi in WHITELIST {
        cmd.args([OsStr::new("--umis"), OsStr::new(umi)]);
    }
    if let Some(m) = metrics {
        cmd.args([OsStr::new("-M"), m.as_os_str()]);
    }
    if let Some(min) = min_corrected {
        cmd.args([OsStr::new("--min-corrected"), OsStr::new(min)]);
    }
    cmd.output().unwrap_or_else(|e| panic!("failed to spawn `{}` correct: {e}", bin.display()))
}

/// A no-`--threads` run now routes through the declarative chain, which logs the
/// `"Using pipeline with N threads"` banner from `ChainBuilder::add_correct`.
/// The retired serial engine logged no such line, so this is the RED
/// (pre-removal) → GREEN (post-removal) discriminator for the cutover.
#[test]
fn correct_no_threads_routes_through_chain() {
    let dir = TempDir::new().expect("temp dir");
    let input = write_varied_input(dir.path());
    let output = dir.path().join("out.bam");

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let out = run_correct(current_bin, &input, &output, None, None);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "no-threads correct must succeed; stderr:\n{stderr}");
    assert!(
        stderr.contains("Using pipeline with 1 threads"),
        "a no-`--threads` correct must route through the chain (which logs the pipeline banner); \
         the serial engine is retired. stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("Starting correct"),
        "the chain must still emit the `Starting correct` banner; stderr:\n{stderr}"
    );
}

/// Output parity of the post-cutover chain against the pre-removal serial
/// baseline binary (run with no `--threads`, i.e. its legacy engine), with and
/// without `--metrics` — plus the always-available self-consistency oracle when
/// no baseline is set.
#[rstest]
#[case::no_metrics(false)]
#[case::with_metrics(true)]
fn cutover_matches_baseline(#[case] with_metrics: bool) {
    let dir = TempDir::new().expect("temp dir");
    let input = write_varied_input(dir.path());

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let current_tsv = dir.path().join("current.tsv");
    let current_metrics = with_metrics.then(|| current_tsv.clone());
    let current = run_correct(current_bin, &input, &current_out, current_metrics.as_deref(), None);
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(current.status.success(), "current correct must succeed; stderr:\n{current_stderr}");

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let baseline_tsv = dir.path().join("baseline.tsv");
        let baseline_metrics = with_metrics.then(|| baseline_tsv.clone());
        let base = run_correct(&baseline, &input, &baseline_out, baseline_metrics.as_deref(), None);
        assert!(
            base.status.success(),
            "baseline correct failed; stderr:\n{}",
            String::from_utf8_lossy(&base.stderr)
        );

        assert_eq!(
            decompressed_records_without_pg(&current_out),
            decompressed_records_without_pg(&baseline_out),
            "chain correct output diverges from the pre-removal serial baseline binary ({}) \
             after stripping @PG — a real cutover parity bug, not something to relax",
            baseline.display(),
        );
        if with_metrics {
            assert_eq!(
                std::fs::read_to_string(&current_tsv).expect("current tsv"),
                std::fs::read_to_string(&baseline_tsv).expect("baseline tsv"),
                "chain --metrics TSV diverges from the serial baseline binary"
            );
        }
    } else {
        eprintln!(
            "SKIP baseline half of cutover_matches_baseline[with_metrics={with_metrics}]: \
             FGUMI_BASELINE_BIN is unset or does not name an existing file — running \
             self-consistency oracle instead"
        );
        assert_self_consistent(&current_out, current_metrics.as_deref());
    }
}

/// fgbio parity across the cutover: a kept ratio below `--min-corrected` fails
/// the run *and* still leaves the `--metrics` TSV behind (CorrectUmis.scala:152).
/// The retired legacy path wrote metrics before its ratio bail; the current
/// chain must match, and — when a baseline is available — the baseline's legacy
/// path must error the same way and produce a byte-identical metrics TSV.
#[test]
fn cutover_min_corrected_failure_writes_metrics() {
    let dir = TempDir::new().expect("temp dir");
    // Every read is `TTTTTTTT`, which matches nothing in WHITELIST within
    // `--max-mismatches 2`: kept/total = 0, below the 0.5 floor.
    let input = dir.path().join("all_reject.bam");
    create_umi_bam(&input, vec![create_umi_family("TTTTTTTT", 6, "reject", "ACGTACGT", 30)]);

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let current_tsv = dir.path().join("current.tsv");
    let current = run_correct(current_bin, &input, &current_out, Some(&current_tsv), Some("0.5"));
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(
        !current.status.success(),
        "correct must fail when the kept ratio is below --min-corrected; stderr:\n{current_stderr}"
    );
    assert!(
        current_stderr.contains("Final ratio of reads kept"),
        "the --min-corrected failure message must survive the cutover; stderr:\n{current_stderr}"
    );
    assert!(
        current_tsv.exists(),
        "the --metrics TSV must be written even on a --min-corrected failure (fgbio parity); \
         stderr:\n{current_stderr}"
    );

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let baseline_tsv = dir.path().join("baseline.tsv");
        let base = run_correct(&baseline, &input, &baseline_out, Some(&baseline_tsv), Some("0.5"));
        assert!(
            !base.status.success(),
            "baseline correct must also fail the --min-corrected gate; stderr:\n{}",
            String::from_utf8_lossy(&base.stderr)
        );
        assert!(
            baseline_tsv.exists(),
            "the baseline's legacy path must also leave a --metrics TSV behind on the failure"
        );
        assert_eq!(
            std::fs::read_to_string(&current_tsv).expect("current tsv"),
            std::fs::read_to_string(&baseline_tsv).expect("baseline tsv"),
            "chain --metrics TSV on a --min-corrected failure diverges from the serial baseline",
        );
    } else {
        eprintln!(
            "SKIP baseline half of cutover_min_corrected_failure_writes_metrics: \
             FGUMI_BASELINE_BIN is unset or does not name an existing file"
        );
        // Parse the TSV and assert the actual counts, not just a line count: all
        // 6 reads are TTTTTTTT (no match), so every whitelist UMI records zero
        // matches and the unmatched all-N bucket takes all 6. A numerically wrong
        // but structurally valid file fails here.
        let parsed = fgumi_lib::metrics::UmiCorrectionMetrics::read_metrics(&current_tsv)
            .expect("parse metrics TSV");
        let by_umi: std::collections::HashMap<&str, &fgumi_lib::metrics::UmiCorrectionMetrics> =
            parsed.iter().map(|m| (m.umi.as_str(), m)).collect();
        for umi in WHITELIST {
            let m = by_umi.get(umi).unwrap_or_else(|| {
                panic!("metrics TSV must include a row for whitelist UMI {umi}")
            });
            assert_eq!(m.total_matches, 0, "no read matches {umi}, so total_matches must be 0");
        }
        let unmatched =
            by_umi.get("NNNNNNNN").expect("metrics TSV must include the unmatched all-N row");
        assert_eq!(
            unmatched.total_matches, 6,
            "all 6 unmatched reads must credit the all-N bucket"
        );
        assert_eq!(
            parsed.iter().map(|m| m.total_matches).sum::<u64>(),
            6,
            "total matches across all rows must equal the 6 input reads"
        );
    }
}

/// Behavior change the cutover intentionally makes (a bug fix, not a
/// regression): a no-`--threads` run over EMPTY input under a positive
/// `--min-corrected` now ERRORS. The retired legacy serial engine silently
/// passed this case (its `0 / 0 = NaN` kept-ratio compared false against the
/// floor and exited 0); the chain — now the only path — treats "no reads
/// processed" as a failure, which is correct since an empty run can never meet a
/// positive minimum. When a baseline is available this is a genuine RED/GREEN
/// discriminator: the frozen pre-removal binary (legacy engine) SUCCEEDS on the
/// same input while the current build fails.
#[test]
fn cutover_empty_input_min_corrected_now_errors() {
    let dir = TempDir::new().expect("temp dir");
    // Header-only BAM: no records at all.
    let input = dir.path().join("empty.bam");
    create_umi_bam(&input, vec![]);

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let current = run_correct(current_bin, &input, &current_out, None, Some("0.5"));
    let current_stderr = String::from_utf8_lossy(&current.stderr);
    assert!(
        !current.status.success(),
        "the cutover makes a no-`--threads` empty-input --min-corrected run error \
         (was a legacy silent pass); stderr:\n{current_stderr}"
    );
    assert!(
        current_stderr.contains("No reads were processed"),
        "the empty-input failure must name the empty-run guard; stderr:\n{current_stderr}"
    );

    if let Some(baseline) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let base = run_correct(&baseline, &input, &baseline_out, None, Some("0.5"));
        assert!(
            base.status.success(),
            "the frozen pre-removal baseline (legacy serial engine) silently passed empty input \
             under --min-corrected; this asymmetry is exactly the bug the cutover fixes. \
             stderr:\n{}",
            String::from_utf8_lossy(&base.stderr)
        );
    }
}

/// Always-available oracle used when no baseline binary is set: the chain's
/// output must keep only records whose UMI is (or corrects to) a whitelist UMI,
/// rewrite the one 1-mismatch family's RX to its exact expected whitelist entry
/// (stashing the pre-correction UMI in `OX`), leave exact matches untouched, and
/// (with `--metrics`) carry a well-formed TSV whose per-UMI `total_matches` match
/// the fixture. Asserting the exact corrected values — not just "some whitelist
/// UMI" and a presence check — is what makes this catch a real correction
/// regression (e.g. a near-match corrected to the wrong entry, or a miscount).
fn assert_self_consistent(output: &Path, metrics: Option<&Path>) {
    use std::collections::HashMap;

    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;

    use fgumi_lib::sam::SamTag;

    let (_, records) = read_bam_output(output);
    // Kept families: exact_a(3) + corr_a(2) + exact_c(3) + exact_g(1) = 9; the
    // reject_t(2) family matches nothing and is dropped (no --rejects here).
    assert_eq!(records.len(), 9, "only correctable/exact records are kept");

    let tag_of =
        |name: &str| -> Tag { Tag::from(name.parse::<SamTag>().expect("valid two-char tag")) };
    let rx_string = |r: &noodles::sam::alignment::RecordBuf, name: &str| -> Option<String> {
        match r.data().get(&tag_of(name))? {
            Value::String(s) => Some(s.to_string()),
            _ => None,
        }
    };

    // Tally the post-correction RX distribution and, for every corrected record
    // (OX present), the exact original -> corrected mapping.
    let mut rx_counts: HashMap<String, usize> = HashMap::new();
    let mut correction_pairs: Vec<(String, String)> = Vec::new();
    for r in &records {
        let rx = rx_string(r, "RX").expect("kept record must carry RX");
        assert!(
            WHITELIST.contains(&rx.as_str()),
            "kept record RX {rx} must be a whitelist UMI after correction"
        );
        if let Some(ox) = rx_string(r, "OX") {
            correction_pairs.push((ox, rx.clone()));
        }
        *rx_counts.entry(rx).or_default() += 1;
    }

    // The 1-mismatch corr_a family (original RX AAAAAAAT) must correct to exactly
    // AAAAAAAA — not merely to "a whitelist UMI". Both of its records carry OX.
    assert_eq!(
        correction_pairs,
        vec![
            ("AAAAAAAT".to_string(), "AAAAAAAA".to_string()),
            ("AAAAAAAT".to_string(), "AAAAAAAA".to_string()),
        ],
        "the two corrected records must map original AAAAAAAT -> corrected AAAAAAAA"
    );

    // Exact per-UMI kept counts: exact_a(3) + corr_a(2, corrected to AAAAAAAA) = 5;
    // exact_c = 3; exact_g = 1. A miscount or a mis-corrected near-match fails here.
    let expected_rx: HashMap<String, usize> =
        [("AAAAAAAA".to_string(), 5), ("CCCCCCCC".to_string(), 3), ("GGGGGGGG".to_string(), 1)]
            .into_iter()
            .collect();
    assert_eq!(rx_counts, expected_rx, "post-correction RX distribution must match the fixture");

    if let Some(path) = metrics {
        let parsed = fgumi_lib::metrics::UmiCorrectionMetrics::read_metrics(path)
            .expect("parse metrics TSV");
        let by_umi: HashMap<&str, &fgumi_lib::metrics::UmiCorrectionMetrics> =
            parsed.iter().map(|m| (m.umi.as_str(), m)).collect();
        // Each whitelist UMI has a row; total_matches equals the kept-count for
        // that UMI (AAAAAAAA credits both the exact and the corrected family).
        for (umi, expected) in [("AAAAAAAA", 5u64), ("CCCCCCCC", 3), ("GGGGGGGG", 1)] {
            let m = by_umi.get(umi).unwrap_or_else(|| {
                panic!("metrics TSV must include a row for whitelist UMI {umi}")
            });
            assert_eq!(m.total_matches, expected, "total_matches for {umi}");
        }
        // The single 1-mismatch correction is recorded against AAAAAAAA.
        assert_eq!(
            by_umi["AAAAAAAA"].one_mismatch_matches, 2,
            "AAAAAAAA must record the two 1-mismatch corrections"
        );
    }
}
