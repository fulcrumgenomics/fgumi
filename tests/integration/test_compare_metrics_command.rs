//! Integration tests for the `fgumi compare metrics` command's CLI surface.
//!
//! Unit tests for the row key-join engine itself live alongside the
//! implementation in `src/lib/commands/compare/metrics.rs`; these tests instead
//! exercise the command end to end: argument parsing, the `RESULT: metrics
//! files ...` report line, `--quiet`/exit-code behavior, and the multi-column
//! `--key-columns` flag via the real CLI surface.

use clap::Parser;
use fgumi_lib::commands::command::Command as FgumiCommand;
use fgumi_lib::commands::compare::{CompareMetrics, CompareMismatch};
use fgumi_lib::variant_review::ConsensusVariantReviewInfo;
use proptest::prelude::*;
use std::path::Path;
use std::process::Command;
use tempfile::TempDir;

use crate::helpers::write_tsv;

/// Runs `CompareMetrics::execute()` in-process; returns `true` on match, `false`
/// on `CompareMismatch` (DIFFER). Panics on any other error.
fn run_compare_in_process(file1: &Path, file2: &Path, extra_args: &[&str]) -> bool {
    let mut args: Vec<&str> = vec!["metrics", file1.to_str().unwrap(), file2.to_str().unwrap()];
    args.extend(extra_args);
    let cmd = CompareMetrics::try_parse_from(args).expect("failed to parse compare metrics args");
    match cmd.execute("fgumi compare metrics") {
        Ok(()) => true,
        Err(e) if e.is::<CompareMismatch>() => false,
        Err(e) => panic!("compare metrics hit unexpected error: {e:#}"),
    }
}

/// Runs `fgumi compare metrics` via subprocess and returns `(exit_code, stdout, stderr)`.
///
/// Returns the *exact* exit code (not just success/failure) so callers can distinguish a
/// clean DIFFER (`Some(1)`, set explicitly for `CompareMismatch` in `main`) from a usage
/// error (`Some(2)`), a panic/abort (`Some(101)`/`Some(134)`), or a signal (`None`) — a
/// nonzero exit that is *not* 1 must never be mistaken for a difference result. `code()`
/// is `None` only when the process was terminated by a signal without exiting. `stderr` is
/// returned so a caller can additionally rule out a generic runtime error (which `main`
/// also exits 1 for, but via the `anyhow` `Error:` report on stderr) masquerading as a
/// clean DIFFER when stdout is suppressed by `--quiet`.
fn run_compare_subprocess(
    file1: &Path,
    file2: &Path,
    extra_args: &[&str],
) -> (Option<i32>, String, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_fgumi"))
        .args(["compare", "metrics"])
        .arg(file1)
        .arg(file2)
        .args(extra_args)
        .output()
        .expect("failed to run fgumi");
    (
        output.status.code(),
        String::from_utf8_lossy(&output.stdout).to_string(),
        String::from_utf8_lossy(&output.stderr).to_string(),
    )
}

#[test]
fn identical_files_match_in_process() {
    let dir = TempDir::new().expect("tempdir");
    let content = "umi\tcount\nAAA\t3\nCCC\t5\n";
    let f1 = write_tsv(dir.path(), "a.txt", content);
    let f2 = write_tsv(dir.path(), "b.txt", content);
    assert!(run_compare_in_process(&f1, &f2, &[]));
}

#[test]
fn reordered_rows_still_match_in_process() {
    let dir = TempDir::new().expect("tempdir");
    let f1 = write_tsv(dir.path(), "a.txt", "family_size\tcount\n1\t10\n2\t5\n3\t1\n");
    let f2 = write_tsv(dir.path(), "b.txt", "family_size\tcount\n3\t1\n1\t10\n2\t5\n");
    assert!(run_compare_in_process(&f1, &f2, &[]), "row order must not affect the key-join");
}

/// The composite key must use *both* named columns. The two rows share `ab_size` and
/// differ only in `ba_size`, so keying on `ab_size` alone would collide them into one
/// duplicate key (which the engine rejects) — only a genuine two-column key gives the two
/// rows distinct identities and lets the row-order-independent join succeed.
#[test]
fn multi_column_key_joins_on_both_columns() {
    let dir = TempDir::new().expect("tempdir");
    let f1 = write_tsv(dir.path(), "a.txt", "ab_size\tba_size\tcount\n1\t0\t100\n1\t5\t20\n");
    let f2 = write_tsv(dir.path(), "b.txt", "ab_size\tba_size\tcount\n1\t5\t20\n1\t0\t100\n");
    assert!(
        run_compare_in_process(&f1, &f2, &["--key-columns", "ab_size,ba_size"]),
        "rows sharing ab_size but differing in ba_size must join by the full (ab_size, ba_size) key"
    );
}

/// The negative direction of the above: with the same two-column key, a value change on
/// exactly one composite-keyed row must DIFFER and name that row's full `(ab_size, ba_size)`
/// identity — proving the join pairs rows by both columns, not just the first.
#[test]
fn multi_column_key_value_change_on_one_row_differs() {
    let dir = TempDir::new().expect("tempdir");
    let f1 = write_tsv(dir.path(), "a.txt", "ab_size\tba_size\tcount\n1\t0\t100\n1\t5\t20\n");
    let f2 = write_tsv(dir.path(), "b.txt", "ab_size\tba_size\tcount\n1\t0\t100\n1\t5\t21\n");
    let (code, stdout, _stderr) =
        run_compare_subprocess(&f1, &f2, &["--key-columns", "ab_size,ba_size"]);
    assert_eq!(code, Some(1), "a value change on one composite-keyed row must DIFFER: {stdout}");
    assert!(stdout.contains("RESULT: metrics files DIFFER"), "stdout was: {stdout}");
    assert!(
        stdout.contains("(1, 5)"),
        "the diff must name the offending composite key (ab_size, ba_size) = (1, 5): {stdout}"
    );
}

/// Exercises the CLI's exit code and the `RESULT: metrics files DIFFER` report
/// line, which downstream tooling (per the `compare bams` precedent) greps for.
#[test]
fn subprocess_reports_result_differ_and_exit_code_1() {
    let dir = TempDir::new().expect("tempdir");
    let f1 = write_tsv(dir.path(), "a.txt", "key\tvalue\nA\t1\n");
    let f2 = write_tsv(dir.path(), "b.txt", "key\tvalue\nA\t2\n");
    let (code, stdout, _stderr) = run_compare_subprocess(&f1, &f2, &[]);
    assert_eq!(
        code,
        Some(1),
        "a clean DIFFER must exit exactly 1, not crash/usage-error: {stdout}"
    );
    assert!(stdout.contains("RESULT: metrics files DIFFER"), "stdout was: {stdout}");
    assert!(!stdout.contains("only in"), "no presence-only rows in this fixture: {stdout}");
}

#[test]
fn subprocess_reports_result_identical_with_verbose() {
    let dir = TempDir::new().expect("tempdir");
    let content = "key\tvalue\nA\t1\n";
    let f1 = write_tsv(dir.path(), "a.txt", content);
    let f2 = write_tsv(dir.path(), "b.txt", content);
    let (code, stdout, _stderr) = run_compare_subprocess(&f1, &f2, &["--verbose"]);
    assert_eq!(code, Some(0), "identical metrics files must exit exactly 0: {stdout}");
    assert!(stdout.contains("RESULT: metrics files are IDENTICAL"), "stdout was: {stdout}");
}

#[test]
fn subprocess_quiet_suppresses_output_on_differ() {
    let dir = TempDir::new().expect("tempdir");
    let f1 = write_tsv(dir.path(), "a.txt", "key\tvalue\nA\t1\n");
    let f2 = write_tsv(dir.path(), "b.txt", "key\tvalue\nA\t2\n");
    let (code, stdout, stderr) = run_compare_subprocess(&f1, &f2, &["--quiet"]);
    assert_eq!(
        code,
        Some(1),
        "a quiet DIFFER must exit exactly 1, not crash/usage-error: {stderr}"
    );
    assert!(stdout.is_empty(), "quiet mode must suppress stdout: {stdout}");
    // The `--quiet` contract is "the exit code is the only *result* signal". stderr is not
    // asserted fully empty: `main` emits a process-wide startup banner (the fgumi version) at
    // `info` before dispatch, under `RUST_LOG` control — the same orthogonal-logging model as
    // fgbio's global `--log-level`. What quiet DOES guarantee is that nothing on stderr
    // communicates the comparison *result* or a runtime error:
    //   - no `anyhow` `Error:` report (a generic runtime error also exits 1, and would
    //     otherwise masquerade as a clean quiet DIFFER),
    //   - no result-determination text (the command's own "Metrics files differ"/"identical"
    //     `info!`s and the "Comparing metrics" timer line are gated behind `!quiet`).
    assert!(!stderr.contains("Error:"), "quiet DIFFER must not emit a runtime error: {stderr}");
    for banned in ["differ", "Differ", "DIFFER", "identical", "IDENTICAL", "Comparing metrics"] {
        assert!(
            !stderr.contains(banned),
            "quiet mode must not leak the result token {banned:?} to stderr: {stderr}"
        );
    }
}

#[test]
fn localized_differ_names_the_offending_key() {
    let dir = TempDir::new().expect("tempdir");
    let f1 = write_tsv(dir.path(), "a.txt", "family_size\tcount\n1\t10\n2\t5\n3\t1\n");
    let f2 = write_tsv(dir.path(), "b.txt", "family_size\tcount\n1\t10\n2\t5\n");
    let (code, stdout, _stderr) = run_compare_subprocess(&f1, &f2, &[]);
    assert_eq!(code, Some(1), "a clean DIFFER must exit exactly 1: {stdout}");
    assert!(stdout.contains("3 only in file1"), "stdout was: {stdout}");
}

// ---------------------------------------------------------------------------
// CMP3-02 (BS2): `compare metrics` is the review-`.txt` comparator.
//
// The `review` command writes its `<output>.txt` as a headered TSV (serde/csv
// serialize of `ConsensusVariantReviewInfo`), so the generic row key-join
// comparator handles it directly — one `<variant-site> × <consensus-read>` row
// per line, joined on `chrom,pos,consensus_read`. These two tests pin that the
// review file's own failure modes are visible to the oracle: row order between
// two `review` runs must be tolerated, and a per-read consensus base/quality
// divergence (the class the whole `review` command's parity findings live in)
// must DIFFER and be localized to the offending variant/read key.
// ---------------------------------------------------------------------------

/// A minimal `review`-`.txt` fixture: the real column header emitted by
/// `ConsensusVariantReviewInfo` (derived from the struct itself via
/// `tsv_header()`, so it can't drift if a field is added/renamed), followed by
/// the given data rows (each a full tab-joined line without its trailing
/// newline).
fn review_txt(rows: &[&str]) -> String {
    let mut out = ConsensusVariantReviewInfo::tsv_header();
    for row in rows {
        out.push('\n');
        out.push_str(row);
    }
    out.push('\n');
    out
}

/// A baseline `ConsensusVariantReviewInfo` for the review-`.txt` fixtures.
///
/// Each fixture starts from this and mutates only the fields it cares about via
/// [`review_row`], so every value stays bound to a *named* struct field. A
/// reorder of two same-typed fields (e.g. two `String`s) is then a compiler- and
/// header-visible change — the row moves in lockstep with the struct-derived
/// header — rather than a silent positional misassignment in a hand-typed
/// tab-joined literal.
fn baseline_review_info() -> ConsensusVariantReviewInfo {
    ConsensusVariantReviewInfo {
        chrom: "chr1".to_string(),
        pos: 100,
        ref_allele: "A".to_string(),
        genotype: "0/1".to_string(),
        filters: "PASS".to_string(),
        consensus_a: 3,
        consensus_c: 0,
        consensus_g: 1,
        consensus_t: 0,
        consensus_n: 0,
        consensus_read: "read".to_string(),
        consensus_insert: "chr1:100-200 | F1R2".to_string(),
        consensus_call: 'A',
        consensus_qual: 40,
        a: 5,
        c: 0,
        g: 2,
        t: 0,
        n: 0,
    }
}

/// Build one review-`.txt` row from the baseline, applying `overrides`, and
/// serialize it through `ConsensusVariantReviewInfo::to_tsv_row()` — the same
/// serde/csv mechanism as the struct-derived header — so the row can never drift
/// from the header's column order.
fn review_row(overrides: impl FnOnce(&mut ConsensusVariantReviewInfo)) -> String {
    let mut info = baseline_review_info();
    overrides(&mut info);
    info.to_tsv_row()
}

/// Build `n` review rows with pairwise-distinct join keys
/// (`chrom,pos,consensus_read`), so any reordering of the set is unambiguously a
/// pure permutation of identical content.
fn distinct_review_rows(n: usize) -> Vec<String> {
    (0..n)
        .map(|i| {
            review_row(|info| {
                info.chrom = format!("chr{i}");
                info.pos = 100 + i32::try_from(i).expect("row index fits i32");
                info.consensus_read = format!("read{i}");
                info.consensus_insert = format!("chr{i}:100-200 | F1R2");
            })
        })
        .collect()
}

/// `readX` — the baseline row, distinguished only by its read name.
fn review_row_x() -> String {
    review_row(|info| info.consensus_read = "readX".to_string())
}

/// The `readY`-specific fields shared by both the reference and the diverged
/// `readY` row *except* its consensus base call and quality. Both rows derive
/// from this, so the `consensus_call`/`consensus_qual` transition is guaranteed
/// to be their *only* difference — a later edit to readY's key/insert/raw counts
/// can't silently turn the divergence test into a multi-column diff.
fn apply_read_y_context(info: &mut ConsensusVariantReviewInfo) {
    info.consensus_read = "readY".to_string();
    info.consensus_insert = "chr1:100-200 | F2R1".to_string();
    info.a = 0;
    info.g = 0;
    info.t = 4;
}

/// `readY` — same variant site as `readX`, a distinct read on the F2R1 insert
/// calling `G@38`.
fn review_row_y() -> String {
    review_row(|info| {
        apply_read_y_context(info);
        info.consensus_call = 'G';
        info.consensus_qual = 38;
    })
}

/// `readZ` — a different variant site (`chr2:250`, ref `C`).
fn review_row_z() -> String {
    review_row(|info| {
        info.chrom = "chr2".to_string();
        info.pos = 250;
        info.ref_allele = "C".to_string();
        info.genotype = "1/1".to_string();
        info.consensus_a = 0;
        info.consensus_c = 2;
        info.consensus_g = 0;
        info.consensus_t = 1;
        info.consensus_read = "readZ".to_string();
        info.consensus_insert = "chr2:250-350 | F1R2".to_string();
        info.consensus_call = 'C';
        info.consensus_qual = 42;
        info.a = 0;
        info.c = 3;
        info.g = 0;
        info.t = 1;
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// Row order between two `review` runs is tolerated for ANY permutation, not
    /// just one hand-picked shuffle: two runs may legitimately emit identical
    /// rows in different order, and `compare metrics`' key-join on
    /// `chrom,pos,consensus_read` exists to tolerate exactly that. Generalizes
    /// the former single fixed shuffle to a property over a generated N-row set
    /// and an arbitrary permutation. (Coding guidelines: use proptest for
    /// property-based testing.)
    #[test]
    fn review_txt_row_order_is_tolerated(
        perm in (1usize..8).prop_flat_map(|n| Just((0..n).collect::<Vec<usize>>()).prop_shuffle())
    ) {
        let rows = distinct_review_rows(perm.len());
        let canonical: Vec<&str> = rows.iter().map(String::as_str).collect();
        let shuffled: Vec<&str> = perm.iter().map(|&i| rows[i].as_str()).collect();

        let dir = TempDir::new().expect("tempdir");
        let f1 = write_tsv(dir.path(), "a.txt", &review_txt(&canonical));
        let f2 = write_tsv(dir.path(), "b.txt", &review_txt(&shuffled));
        prop_assert!(
            run_compare_in_process(&f1, &f2, &["--key-columns", "chrom,pos,consensus_read"]),
            "reordered review rows with identical content must MATCH (perm={perm:?})"
        );
    }
}

#[test]
fn review_txt_consensus_call_and_qual_divergence_differs_and_is_localized() {
    let dir = TempDir::new().expect("tempdir");
    let (row_x, row_y, row_z) = (review_row_x(), review_row_y(), review_row_z());
    let f1 = write_tsv(
        dir.path(),
        "a.txt",
        &review_txt(&[row_x.as_str(), row_y.as_str(), row_z.as_str()]),
    );
    // readY's consensus base G→T and quality 38→20 — the per-read consensus
    // divergence class the `review` command's parity findings live in. Shares
    // `apply_read_y_context` with `review_row_y` so call+qual are the ONLY diff.
    let row_y_diff = review_row(|info| {
        apply_read_y_context(info);
        info.consensus_call = 'T';
        info.consensus_qual = 20;
    });
    let f2 = write_tsv(
        dir.path(),
        "b.txt",
        &review_txt(&[row_x.as_str(), row_y_diff.as_str(), row_z.as_str()]),
    );
    let (code, stdout, _stderr) =
        run_compare_subprocess(&f1, &f2, &["--key-columns", "chrom,pos,consensus_read"]);
    assert_eq!(code, Some(1), "a consensus base/qual divergence must exit non-zero (DIFFER)");
    assert!(stdout.contains("RESULT: metrics files DIFFER"), "stdout was: {stdout}");
    // Localized to the offending variant-site × read key (chr1:100, readY)...
    assert!(stdout.contains("chr1"), "diff must name the offending contig: {stdout}");
    assert!(stdout.contains("100"), "diff must name the offending position: {stdout}");
    assert!(stdout.contains("readY"), "diff must name the offending read: {stdout}");
    // ...naming both changed columns and their exact value transitions.
    assert!(stdout.contains("consensus_call"), "diff must name the call column: {stdout}");
    assert!(stdout.contains(r#""G" != "T""#), "diff must show the base transition: {stdout}");
    assert!(stdout.contains("consensus_qual"), "diff must name the qual column: {stdout}");
    assert!(stdout.contains(r#""38" != "20""#), "diff must show the qual transition: {stdout}");
    // The unchanged readX/readZ rows must NOT be reported (no cascade).
    assert!(!stdout.contains("readX"), "unchanged rows must not appear: {stdout}");
    assert!(!stdout.contains("readZ"), "unchanged rows must not appear: {stdout}");
}
