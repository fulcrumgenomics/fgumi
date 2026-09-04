//! Shared helpers for command-cutover parity tests.
//!
//! A "cutover" test proves that routing a command through the declarative chain
//! builder produces output byte-identical to the pre-cutover owned-engine
//! baseline binary. Every such test needs the same three primitives: resolve the
//! baseline binary from the environment, and compare two BAMs' decompressed
//! record bytes with the necessarily-differing `@PG` line stripped. They live
//! here so each per-command cutover test (`test_sort_cutover_parity.rs`,
//! `test_retag_cutover_parity.rs`, and later ones) shares one copy.

use std::path::{Path, PathBuf};

/// Resolves the saved owned-engine baseline binary to compare a cutover against.
///
/// The baseline path comes solely from the `FGUMI_BASELINE_BIN` env var; when it
/// is unset (or names a missing file) this returns `None` — "no baseline oracle
/// available", never a passing result. Callers layer the baseline byte-parity
/// check on top of an always-available oracle (e.g. `fgumi sort --verify`, or a
/// self-consistency check); a missing baseline just drops the byte-parity half,
/// it never skips the case outright. No hardcoded fallback: a baseline binary is
/// host-specific and must never be a path committed into the repo.
#[must_use]
pub fn baseline_bin() -> Option<PathBuf> {
    let path = PathBuf::from(std::env::var_os("FGUMI_BASELINE_BIN")?);
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
/// Both a cutover build and the baseline binary stamp a single `@PG` line whose
/// `VN` (git-describe version) and `CL` (command line, naming argv[0] — a
/// different binary path for each side — and the per-run output path) necessarily
/// differ between two independent invocations. Stripping the whole line (not just
/// normalizing those two sub-fields) is simplest and correct because neither side
/// emits more than one `@PG` record for these inputs (the test BAMs carry no
/// pre-existing `@PG`).
#[must_use]
pub fn strip_pg_lines(text: &str) -> String {
    // Empty in, empty out: `"".split('\n')` yields `[""]`, which the
    // trailing-newline logic below would otherwise turn into `"\n"`.
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
    // Restore the trailing newline iff the source had one — never synthesize one
    // the source lacked.
    if had_trailing_newline {
        out.push('\n');
    }
    out
}

/// Reads a BAM file's raw BGZF stream, decompresses it, strips the `@PG` line(s)
/// from the embedded SAM header text, and returns the resulting bytes (new header
/// + unmodified `n_ref`/reference-list/record bytes).
///
/// This compares the two writers' *decompressed* output rather than
/// `RecordBuf`-parsed records or raw file bytes: BGZF block boundaries and
/// compression levels differ between the two writers even when the logical
/// content is identical, so raw file bytes never match; parsing into `RecordBuf`
/// and re-encoding risks masking a real tag-order or binary-layout regression by
/// normalizing it away. Operating directly on the decompressed BAM binary format
/// keeps everything except the `@PG` line an exact, uninterpreted byte comparison.
#[must_use]
pub fn decompressed_records_without_pg(path: &Path) -> Vec<u8> {
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
