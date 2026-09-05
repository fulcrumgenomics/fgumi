//! Cutover parity test for `fgumi copy-umi`.
//!
//! `CopyUmi::execute` used to gate on `--threads`: unset ran a serial
//! read->copy-umi->write loop (`run_single_threaded`), set routed through the
//! declarative chain builder. The cutover deletes the serial path so `execute()`
//! always runs through the chain (`execute_chain`) regardless of `--threads`.
//! `copy-umi` itself is unchanged by the cutover -- only which code path
//! `execute()` reaches is -- so the pre-cutover binary is a valid oracle for the
//! post-cutover binary's output.
//!
//! This pins the current (chain-only) binary's output against the frozen
//! pre-cutover baseline binary named by `FGUMI_BASELINE_BIN`: BAM output must be
//! byte-identical modulo the `@PG` header line, and the `--metrics` TSV (which
//! embeds no timestamp) must be byte-identical outright. When `FGUMI_BASELINE_BIN`
//! is unset (the default in CI), the baseline half is skipped but the case still
//! runs a non-vacuous self-consistency oracle: it asserts the destination `RX`
//! tag is present and equals the *normalized* (not raw) source UMI, using fixed
//! expected values pinned independently in `umi::read_name`'s own unit tests
//! (`rAAAA+CCCC` -> `TTTT-CCCC`, etc.) -- a passthrough that copied the read-name
//! field verbatim, or dropped the `r`-prefix reverse-complement / `+`->`-`
//! translation, would fail this even with no baseline binary available.

use std::ffi::OsStr;
use std::path::Path;
use std::process::Command;

use fgumi_raw_bam::RawRecord;
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::{
    read_copy_umi_metrics, read_name_and_rx, record_named, record_named_with_rx, write_input,
};

/// Resolves the frozen pre-cutover baseline binary to compare against.
///
/// Comes solely from `FGUMI_BASELINE_BIN`; unset (or naming a missing file)
/// returns `None` -- "no baseline oracle available", never a passing result.
/// Callers layer the baseline byte-parity check on top of the always-available
/// self-consistency oracle; a missing baseline only drops the byte-parity half.
/// No hardcoded fallback: a baseline binary is host-specific and must never be a
/// path committed into the repo.
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
/// Both the cutover build and the baseline binary stamp a single `@PG` line
/// whose `VN` (git-describe version) and `CL` (command line, naming a different
/// binary path and output path per side) necessarily differ between two
/// independent invocations. Stripping the whole line is simplest and correct
/// here because neither side emits more than one `@PG` record for this input.
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

/// Reads a BAM file's raw BGZF stream, decompresses it, strips the `@PG`
/// line(s) from the embedded SAM header text, and returns the resulting bytes
/// (new header + unmodified `n_ref`/reference-list/record bytes).
///
/// Operates on the decompressed BAM binary format directly (not
/// `RecordBuf`-parsed records) so everything except the `@PG` line is an exact,
/// uninterpreted byte comparison -- BGZF block boundaries and compression
/// levels differ between the two writers even when the logical content is
/// identical, so raw file bytes never match, and a re-encode through noodles
/// risks masking a real tag-order or binary-layout regression by normalizing it
/// away. Mirrors `test_sort_cutover_parity.rs`'s helper of the same name.
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

/// Runs `<bin> copy-umi -i <input> -o <output> [extra...]`. Returns `Ok(())` on
/// success or `Err(stderr)` on failure -- never panics on a non-zero exit, so
/// callers can assert on either outcome.
fn run_copy_umi(bin: &Path, input: &Path, output: &Path, extra: &[&str]) -> Result<(), String> {
    let out = Command::new(bin)
        .args([OsStr::new("copy-umi"), OsStr::new("-i"), input.as_os_str()])
        .args([OsStr::new("-o"), output.as_os_str()])
        .args(extra.iter().map(OsStr::new))
        .output()
        .unwrap_or_else(|e| panic!("failed to spawn `{}`: {e}", bin.display()));
    if out.status.success() {
        Ok(())
    } else {
        Err(String::from_utf8_lossy(&out.stderr).into_owned())
    }
}

// ============================================================================
// Records + fixed expected (name, RX) values.
//
// The raw last-field UMI text (after the `:` delimiter) for each record, and
// the RX each must normalize to under `--reverse-complement-r-umis` true
// (fgbio default) vs. false, are pinned independently of this cutover: they
// are the exact values `umi::read_name::normalize_read_name_umi`'s own unit
// tests assert (`rAAAA+CCCC` -> `TTTT-CCCC` under revcomp, `AAAA-CCCC` under
// strip-only, etc.) A passthrough implementation that copied the raw field
// verbatim, or dropped the `+`->`-` translation or the `r`-prefix handling,
// diverges from these on every record but the first.
// ============================================================================

/// `(name_prefix, raw_umi_field)` pairs; the read name is `"{prefix}:{raw}"`.
const RECORD_SPECS: &[(&str, &str)] = &[
    ("read0", "ACGT"),
    ("read1", "rAAAA"),
    ("read2", "ACGT+CAGA"),
    ("read3", "rAAAA+CCCC"),
    ("read4", "rAAAA+rCCCC"),
];

/// Expected RX per `RECORD_SPECS` entry, in order, under
/// `--reverse-complement-r-umis` true (the default) or false (strip-only).
fn expected_rx(reverse_complement: bool) -> Vec<&'static str> {
    if reverse_complement {
        vec!["ACGT", "TTTT", "ACGT-CAGA", "TTTT-CCCC", "TTTT-GGGG"]
    } else {
        vec!["ACGT", "AAAA", "ACGT-CAGA", "AAAA-CCCC", "AAAA-CCCC"]
    }
}

fn parity_records() -> Vec<RawRecord> {
    RECORD_SPECS.iter().map(|(prefix, raw)| record_named(&format!("{prefix}:{raw}"))).collect()
}

// ============================================================================
// Main matrix: every non-default, output-affecting flag combination, byte-
// pinned against the baseline binary where available, and self-consistency-
// checked against fixed expected (name, RX) values regardless.
// ============================================================================

#[rstest]
#[case::default(&[][..], false, true)]
#[case::remove_umi(&["--remove-umi"][..], true, true)]
#[case::strip_r(&["--reverse-complement-r-umis", "false"][..], false, false)]
#[case::remove_and_strip(&["--remove-umi", "--reverse-complement-r-umis", "false"][..], true, false)]
fn cutover_matches_baseline_and_self_consistency(
    #[case] flags: &[&str],
    #[case] removes_umi: bool,
    #[case] reverse_complement: bool,
) {
    let dir = TempDir::new().expect("create temp dir");
    let input = write_input(dir.path(), &parity_records());

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let current_metrics = dir.path().join("current.tsv");
    let mut current_flags = flags.to_vec();
    let metrics_arg = current_metrics.display().to_string();
    current_flags.extend_from_slice(&["-M", &metrics_arg]);
    run_copy_umi(current_bin, &input, &current_out, &current_flags)
        .unwrap_or_else(|e| panic!("current binary copy-umi failed (flags={flags:?}): {e}"));

    // Self-consistency oracle: always runs, proving THIS run actually copied and
    // normalized the UMI (not a passthrough), regardless of baseline availability.
    let name_rx = read_name_and_rx(dir.path(), &current_out);
    let expected_rx_vals = expected_rx(reverse_complement);
    assert_eq!(name_rx.len(), RECORD_SPECS.len(), "record count must round-trip (flags={flags:?})");
    for (i, ((name, rx), (prefix, _raw))) in name_rx.iter().zip(RECORD_SPECS.iter()).enumerate() {
        let expected_name = if removes_umi {
            (*prefix).to_string()
        } else {
            format!("{prefix}:{}", RECORD_SPECS[i].1)
        };
        assert_eq!(name, &expected_name, "record {i} name mismatch (flags={flags:?})");
        assert_eq!(
            rx.as_deref(),
            Some(expected_rx_vals[i]),
            "record {i} RX mismatch (flags={flags:?}): expected {}, a passthrough or a dropped \
             normalization step would diverge here",
            expected_rx_vals[i]
        );
    }

    let metrics = read_copy_umi_metrics(&current_metrics);
    assert_eq!(metrics.total_records, 5, "flags={flags:?}");
    assert_eq!(metrics.rx_written, 5, "flags={flags:?}");
    assert_eq!(metrics.rx_overwritten, 0, "flags={flags:?}");
    let expected_trimmed: u64 = if removes_umi { 5 } else { 0 };
    assert_eq!(metrics.names_trimmed, expected_trimmed, "flags={flags:?}");

    match baseline_bin() {
        Some(baseline_bin_path) => {
            let baseline_out = dir.path().join("baseline.bam");
            let baseline_metrics = dir.path().join("baseline.tsv");
            let mut baseline_flags = flags.to_vec();
            let baseline_metrics_arg = baseline_metrics.display().to_string();
            baseline_flags.extend_from_slice(&["-M", &baseline_metrics_arg]);
            run_copy_umi(&baseline_bin_path, &input, &baseline_out, &baseline_flags)
                .unwrap_or_else(|e| {
                    panic!("baseline binary copy-umi failed (flags={flags:?}): {e}")
                });

            assert_eq!(
                decompressed_records_without_pg(&current_out),
                decompressed_records_without_pg(&baseline_out),
                "cutover output for flags={flags:?} diverges from the pre-cutover baseline binary \
                 ({}) after stripping the @PG header line -- this is a real cutover parity bug, \
                 not something to relax the assertion for",
                baseline_bin_path.display(),
            );
            assert_eq!(
                std::fs::read(&current_metrics).expect("read current metrics"),
                std::fs::read(&baseline_metrics).expect("read baseline metrics"),
                "the --metrics TSV must be byte-identical between the cutover and the baseline \
                 binary (flags={flags:?})"
            );
        }
        None => {
            eprintln!(
                "SKIP baseline half of cutover_matches_baseline_and_self_consistency[{flags:?}]: \
                 FGUMI_BASELINE_BIN is unset or does not name an existing file"
            );
        }
    }
}

// ============================================================================
// RX overwrite: a pre-existing RX tag must be overwritten (default policy) and
// counted, or reject the run under `--fail-if-tag-present` -- both byte-pinned
// against the baseline where available.
// ============================================================================

fn overwrite_records() -> Vec<RawRecord> {
    vec![record_named("plain:ACGT"), record_named_with_rx("stale:rAAAA", "GATTACA")]
}

#[test]
fn cutover_overwrites_existing_rx_and_counts_it() {
    let dir = TempDir::new().expect("create temp dir");
    let input = write_input(dir.path(), &overwrite_records());

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let current_metrics = dir.path().join("current.tsv");
    run_copy_umi(
        current_bin,
        &input,
        &current_out,
        &["-M", current_metrics.to_str().expect("valid utf8 path")],
    )
    .expect("current binary copy-umi (overwrite) should succeed");

    let name_rx = read_name_and_rx(dir.path(), &current_out);
    assert_eq!(
        name_rx,
        vec![
            ("plain:ACGT".to_string(), Some("ACGT".to_string())),
            ("stale:rAAAA".to_string(), Some("TTTT".to_string())),
        ],
        "the pre-existing RX on record 1 must be overwritten with the freshly-copied UMI"
    );

    let metrics = read_copy_umi_metrics(&current_metrics);
    assert_eq!(metrics.total_records, 2);
    assert_eq!(metrics.rx_overwritten, 1);

    if let Some(baseline_bin_path) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let baseline_metrics = dir.path().join("baseline.tsv");
        run_copy_umi(
            &baseline_bin_path,
            &input,
            &baseline_out,
            &["-M", baseline_metrics.to_str().expect("valid utf8 path")],
        )
        .expect("baseline binary copy-umi (overwrite) should succeed");
        assert_eq!(
            decompressed_records_without_pg(&current_out),
            decompressed_records_without_pg(&baseline_out),
            "overwrite-scenario output diverges from the baseline binary after stripping @PG"
        );
        assert_eq!(
            std::fs::read(&current_metrics).expect("read current metrics"),
            std::fs::read(&baseline_metrics).expect("read baseline metrics"),
            "the --metrics TSV must be byte-identical between the cutover and the baseline binary"
        );
    } else {
        eprintln!(
            "SKIP baseline half of cutover_overwrites_existing_rx_and_counts_it: \
             FGUMI_BASELINE_BIN is unset or does not name an existing file"
        );
    }
}

/// `--fail-if-tag-present` must reject a record already carrying an `RX` tag,
/// on both the cutover binary and (where available) the baseline binary, with
/// the same error substring.
#[test]
fn cutover_fail_if_tag_present_rejects_existing_rx() {
    let dir = TempDir::new().expect("create temp dir");
    let input = write_input(dir.path(), &overwrite_records());

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    let err = run_copy_umi(current_bin, &input, &current_out, &["--fail-if-tag-present"])
        .expect_err("--fail-if-tag-present must reject a pre-existing RX tag");
    assert!(
        err.contains("already has an RX tag"),
        "expected the fail-if-tag-present error to name the reason; got: {err}"
    );

    if let Some(baseline_bin_path) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        let baseline_err =
            run_copy_umi(&baseline_bin_path, &input, &baseline_out, &["--fail-if-tag-present"])
                .expect_err("baseline binary must also reject a pre-existing RX tag");
        assert!(
            baseline_err.contains("already has an RX tag"),
            "baseline binary's --fail-if-tag-present error differs; got: {baseline_err}"
        );
    } else {
        eprintln!(
            "SKIP baseline half of cutover_fail_if_tag_present_rejects_existing_rx: \
             FGUMI_BASELINE_BIN is unset or does not name an existing file"
        );
    }
}

// ============================================================================
// Custom `--field-delimiter`.
// ============================================================================

#[test]
fn cutover_custom_field_delimiter_matches_baseline_and_self_consistency() {
    let dir = TempDir::new().expect("create temp dir");
    let records: Vec<RawRecord> = [("read0", "ACGT"), ("read1", "rAAAA+CCCC")]
        .iter()
        .map(|(prefix, raw)| record_named(&format!("{prefix}_{raw}")))
        .collect();
    let input = write_input(dir.path(), &records);

    let current_bin = Path::new(env!("CARGO_BIN_EXE_fgumi"));
    let current_out = dir.path().join("current.bam");
    run_copy_umi(current_bin, &input, &current_out, &["--field-delimiter", "_"])
        .expect("current binary copy-umi (custom delimiter) should succeed");

    let name_rx = read_name_and_rx(dir.path(), &current_out);
    assert_eq!(
        name_rx,
        vec![
            ("read0_ACGT".to_string(), Some("ACGT".to_string())),
            ("read1_rAAAA+CCCC".to_string(), Some("TTTT-CCCC".to_string())),
        ],
        "a `_`-delimited last field must be located and normalized identically to `:`"
    );

    if let Some(baseline_bin_path) = baseline_bin() {
        let baseline_out = dir.path().join("baseline.bam");
        run_copy_umi(&baseline_bin_path, &input, &baseline_out, &["--field-delimiter", "_"])
            .expect("baseline binary copy-umi (custom delimiter) should succeed");
        assert_eq!(
            decompressed_records_without_pg(&current_out),
            decompressed_records_without_pg(&baseline_out),
            "custom-field-delimiter output diverges from the baseline binary after stripping @PG"
        );
    } else {
        eprintln!(
            "SKIP baseline half of cutover_custom_field_delimiter_matches_baseline_and_self_consistency: \
             FGUMI_BASELINE_BIN is unset or does not name an existing file"
        );
    }
}
