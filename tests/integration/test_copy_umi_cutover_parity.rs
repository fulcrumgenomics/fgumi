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

use fgumi_raw_bam::{RawRecord, SamBuilder};
use rstest::rstest;
use tempfile::TempDir;

use crate::helpers::{
    read_bam_output, read_copy_umi_metrics, read_name_and_rx, record_named, record_named_with_rx,
    write_input,
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

/// Asserts `copy-umi` preserved every alignment-bearing field of every record
/// against the `input` fixture: it rewrites only the read name and the `RX` tag,
/// so SEQ, QUAL, FLAG, RNAME, POS, MAPQ, CIGAR, the mate reference id, the mate
/// position, the template length (TLEN), and every non-`RX` auxiliary tag —
/// plus the header's reference sequences — must round-trip unchanged. This
/// upgrades the always-on (baseline-free) oracle from a name/`RX`-only check to
/// full-record identity, so a shared chain regression that corrupts a record
/// body is caught in CI where `FGUMI_BASELINE_BIN` is unset. Only the read name,
/// the `RX` tag, and the `@PG` header line are excluded — those are the parts
/// `copy-umi` is expected to change (the name/`RX` transformation itself is
/// asserted separately by each caller).
fn assert_alignment_fields_preserved(input: &Path, output: &Path) {
    // The one aux tag `copy-umi` writes; every other tag must survive verbatim.
    let rx_tag = fgumi_lib::sam::SamTag::RX.to_noodles_tag();
    let (mut in_header, in_recs) = read_bam_output(input);
    let (mut out_header, out_recs) = read_bam_output(output);
    assert_eq!(out_recs.len(), in_recs.len(), "record count must round-trip");
    // Every header record except `@PG` must round-trip: `copy-umi` stamps its
    // own `@PG` line, so clear the programs on both sides and compare the whole
    // header. A reference-sequences-only check would let a dropped or altered
    // `@HD` (or `@RG`/`@CO`) pass in the baseline-free CI path; comparing the
    // full header modulo `@PG` closes that gap.
    in_header.programs_mut().as_mut().clear();
    out_header.programs_mut().as_mut().clear();
    assert_eq!(
        out_header, in_header,
        "copy-umi must preserve every non-@PG header record (@HD, @SQ, @RG, @CO)"
    );
    for (i, (out, inp)) in out_recs.iter().zip(&in_recs).enumerate() {
        assert_eq!(out.sequence().as_ref(), inp.sequence().as_ref(), "record {i} SEQ must survive");
        assert_eq!(
            out.quality_scores().as_ref(),
            inp.quality_scores().as_ref(),
            "record {i} QUAL must survive"
        );
        assert_eq!(out.flags(), inp.flags(), "record {i} FLAG must survive");
        assert_eq!(
            out.reference_sequence_id(),
            inp.reference_sequence_id(),
            "record {i} RNAME must survive"
        );
        assert_eq!(out.alignment_start(), inp.alignment_start(), "record {i} POS must survive");
        assert_eq!(out.mapping_quality(), inp.mapping_quality(), "record {i} MAPQ must survive");
        let out_cigar: Vec<_> =
            out.cigar().as_ref().iter().map(|op| (op.kind(), op.len())).collect();
        let in_cigar: Vec<_> =
            inp.cigar().as_ref().iter().map(|op| (op.kind(), op.len())).collect();
        assert_eq!(out_cigar, in_cigar, "record {i} CIGAR must survive");
        assert_eq!(
            out.mate_reference_sequence_id(),
            inp.mate_reference_sequence_id(),
            "record {i} mate RNAME (mate reference id) must survive"
        );
        assert_eq!(
            out.mate_alignment_start(),
            inp.mate_alignment_start(),
            "record {i} mate POS must survive"
        );
        assert_eq!(
            out.template_length(),
            inp.template_length(),
            "record {i} TLEN (template length) must survive"
        );
        // Every auxiliary tag except `RX` (the one `copy-umi` writes) must
        // survive verbatim, in order.
        let out_aux: Vec<_> = out.data().iter().filter(|(tag, _)| *tag != rx_tag).collect();
        let in_aux: Vec<_> = inp.data().iter().filter(|(tag, _)| *tag != rx_tag).collect();
        assert_eq!(out_aux, in_aux, "record {i} non-RX auxiliary tags must survive");
    }
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

/// The 28-byte BGZF end-of-file marker (an empty gzip block). A fully written
/// BGZF stream ends with it; a run that aborts mid-stream (e.g. a
/// `--fail-if-tag-present` rejection) never finalizes the writer, so the partial
/// output it leaves lacks it.
const BGZF_EOF_MARKER: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// True if `path` ends with the BGZF EOF marker, i.e. the writer finalized the
/// stream. A partial BAM from an aborted run returns false.
fn has_bgzf_eof_marker(path: &Path) -> bool {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    bytes.ends_with(&BGZF_EOF_MARKER)
}

/// Reads whatever *complete* records a possibly-truncated BAM contains, tolerating
/// a missing EOF marker or a partial trailing block. A header-only partial BAM
/// (the shape a failed `copy-umi` run leaves) yields an empty vec rather than a
/// panic, so callers can assert the aborted run published fewer records than it read.
fn recover_records_lenient(path: &Path) -> Vec<noodles::sam::alignment::RecordBuf> {
    let mut reader = noodles::bam::io::Reader::new(std::io::BufReader::new(
        std::fs::File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display())),
    ));
    // A truncated or absent header is a legitimate shape for an aborted run's
    // partial BAM; report zero recovered records rather than panicking, so
    // `failed_output_state` stays comparable for zero-byte / header-truncated output.
    let Ok(header) = reader.read_header() else {
        return Vec::new();
    };
    let mut records = Vec::new();
    for result in reader.record_bufs(&header) {
        match result {
            Ok(record) => records.push(record),
            // A truncated tail is the expected shape of a failed run's partial
            // BAM; stop at the first unreadable record instead of failing.
            Err(_) => break,
        }
    }
    records
}

/// The observable post-failure state of a `copy-umi` output path: whether the
/// file exists, whether its BGZF stream was finalized (EOF marker present), and
/// how many complete records survived. Comparing this between the cutover and the
/// frozen baseline binary is how the failure test pins output-state parity, not
/// just the stderr message.
#[derive(Debug, PartialEq, Eq)]
struct FailedOutputState {
    exists: bool,
    finalized: bool,
    recovered_records: usize,
}

fn failed_output_state(path: &Path) -> FailedOutputState {
    if !path.exists() {
        return FailedOutputState { exists: false, finalized: false, recovered_records: 0 };
    }
    FailedOutputState {
        exists: true,
        finalized: has_bgzf_eof_marker(path),
        recovered_records: recover_records_lenient(path).len(),
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

/// Builds a parity fixture record: the same minimal alignment body as
/// [`record_named`], but with the mate reference id, mate position, template
/// length (TLEN), and a non-`RX` `NM:i` auxiliary tag populated with concrete,
/// per-record, non-default values.
///
/// `record_named` leaves all of these at their BAM defaults (mate ref/pos = -1,
/// TLEN = 0, no aux tags), so the always-on / baseline-free oracle
/// ([`assert_alignment_fields_preserved`]) could not observe a `copy-umi` chain
/// regression that dropped or corrupted them — a dropped mate ref/pos, TLEN, or
/// aux tag would compare "default == default" on both sides and pass silently.
/// Populating them here makes such a regression actually diverge input vs.
/// output. `copy-umi` rewrites only the read name and the `RX` tag, so every
/// field set here must round-trip unchanged.
fn record_with_mate_and_aux(name: &str, mate_pos: i32, tlen: i32, nm: i32) -> RawRecord {
    let mut b = SamBuilder::new();
    b.read_name(name.as_bytes())
        .sequence(b"ACGT")
        .qualities(&[30; 4])
        .flags(0)
        .ref_id(0)
        .pos(99)
        .mapq(60)
        .cigar_ops(&[4 << 4])
        .mate_ref_id(0)
        .mate_pos(mate_pos)
        .template_length(tlen)
        .add_int_tag(fgumi_lib::sam::SamTag::NM, nm);
    b.build()
}

fn parity_records() -> Vec<RawRecord> {
    RECORD_SPECS
        .iter()
        .enumerate()
        .map(|(i, (prefix, raw))| {
            let i = i32::try_from(i).expect("fixture index fits i32");
            // Distinct, non-default values per record so a dropped/corrupted
            // field diverges rather than aliasing another record's value.
            record_with_mate_and_aux(&format!("{prefix}:{raw}"), 200 + i, 300 + i, 1 + i)
        })
        .collect()
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

    // Full-record oracle (always on, baseline or not): copy-umi must leave every
    // alignment-bearing field and the header refs untouched — only name + RX change.
    assert_alignment_fields_preserved(&input, &current_out);

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

    // Full-record oracle (always on): only name + RX change on the overwrite path.
    assert_alignment_fields_preserved(&input, &current_out);

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

    // Output-state half of the contract (stderr is the other half): the rejection
    // aborts the run mid-stream, so the chain builder documents it leaves a partial
    // BAM. Pin that — the writer must not have finalized the BGZF stream, and must
    // not have published a complete copy of the input (which would mean the output
    // path silently ignored --fail-if-tag-present even though stderr reported it).
    let current_state = failed_output_state(&current_out);
    assert!(
        current_state.exists,
        "a failed copy-umi run must leave its (partial) output path behind, not remove it"
    );
    assert!(
        !current_state.finalized,
        "a failed copy-umi run must not finalize the BGZF stream (a partial BAM is \
         expected); {} ends with the BGZF EOF marker",
        current_out.display()
    );
    assert!(
        current_state.recovered_records < overwrite_records().len(),
        "a --fail-if-tag-present rejection must not publish a complete output; recovered \
         {} of {} records",
        current_state.recovered_records,
        overwrite_records().len()
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
        // Parity extends to the partial BAM the aborted run leaves behind, not just
        // the error message: the cutover must leave the same post-failure output
        // state (existence, finalization, surviving record count) as the frozen
        // pre-cutover binary.
        assert_eq!(
            failed_output_state(&baseline_out),
            current_state,
            "cutover and baseline must leave the same post-failure output state"
        );
    } else {
        eprintln!(
            "SKIP baseline half of cutover_fail_if_tag_present_rejects_existing_rx: \
             FGUMI_BASELINE_BIN is unset or does not name an existing file"
        );
    }
}

/// `failed_output_state` must stay total for every shape an aborted run can leave
/// — including a zero-byte or header-truncated output — rather than panicking in
/// `recover_records_lenient`. `failed_output_state` is the value both sides of the
/// baseline-parity `assert_eq!` need, so a panic here would take the whole parity
/// comparison down, not just the record count.
#[rstest]
#[case::empty(&[])]
#[case::truncated_header(b"BAM\x01\x00\x00")]
#[case::garbage(b"not a bam at all")]
fn failed_output_state_is_total_on_partial_output(#[case] bytes: &[u8]) {
    let dir = TempDir::new().expect("create temp dir");
    let path = dir.path().join("partial.bam");
    std::fs::write(&path, bytes).expect("write partial output");

    // No panic, and the header-less shape reports the expected comparable state:
    // the file exists, the stream was never finalized, and no complete record
    // survived.
    let state = failed_output_state(&path);
    assert_eq!(
        state,
        FailedOutputState { exists: true, finalized: false, recovered_records: 0 },
        "a zero-byte / truncated-header partial output must yield a comparable \
         header-less state, not a panic"
    );
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

    // Full-record oracle (always on): only name + RX change under a custom delimiter.
    assert_alignment_fields_preserved(&input, &current_out);

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
