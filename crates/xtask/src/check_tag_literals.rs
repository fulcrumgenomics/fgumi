//! Scans the workspace for bare 2-byte byte-literal SAM aux tag identifiers that
//! should use `SamTag::XX` constants instead.
//!
//! Exits non-zero if any findings remain outside the allowlist. Wired as
//! `cargo xtask check-tag-literals` and run as part of CI.
//!
//! # Allowlist strategy
//!
//! Three complementary allowlists prevent false positives:
//!
//! * **File allowlist** — files that are entirely exempt from the scan:
//!   - `crates/fgumi-tag/src/` — the tag definition crate itself (defines the constants).
//!   - `crates/xtask/src/check_tag_literals.rs` — this file.
//!   - `benches/` — microbenchmarks that deliberately exercise the raw byte API to
//!     measure its overhead without the `SamTag` abstraction layer.
//!
//! * **Payload allowlist** — specific 2-byte sequences that are known-good across the
//!   rest of the codebase because they are not SAM aux-tag identifiers:
//!   - SAM header sub-fields (`SO`, `GO`, `SS`) accepted by the noodles header API.
//!   - Read-name fixtures (`q1`, `q2`, `r0`–`r4`, `rd`).
//!   - 2-bp DNA sequence fragments (`AA`, `AC`, `AT`, …) used in sequence tests.
//!   - Opaque test-fixture tag names (`XA`–`XY`, `x0`–`x9`, `ZZ`, `Zz`, …) used to
//!     exercise SAM/BAM parsers without carrying semantic meaning.
//!
//! * **Path-scoped allowlist** — payloads that are known-good only within specific
//!   files (e.g. intentional compat-proof tests, format-specific constants):
//!   - `b"cd"` / `b"RX"` in the raw-bam builder/tags compat tests.
//!   - `b"BC"` in the BGZF writer (extra-subfield ID, not a SAM aux tag).

use anyhow::{Context, Result};
use std::fs;
use std::path::{Path, PathBuf};

/// A single finding: a bare 2-byte byte literal that should use a `SamTag` constant.
#[derive(Debug, Clone)]
pub struct Finding {
    /// Path relative to the workspace root.
    pub path: PathBuf,
    /// 1-based line number.
    pub line: usize,
    /// Trimmed source line (for display).
    pub snippet: String,
    /// The 2-byte payload of the literal.
    pub payload: [u8; 2],
}

/// Workspace-relative paths that are entirely exempt from scanning.
const FILE_ALLOWLIST: &[&str] = &[
    // Tag-definition crate — contains the constants themselves.
    "crates/fgumi-tag/src/tag.rs",
    "crates/fgumi-tag/src/lib.rs",
    // This scanner file.
    "crates/xtask/src/check_tag_literals.rs",
];

/// Path prefixes (relative to workspace root) that are entirely exempt.
const DIR_ALLOWLIST: &[&str] = &[
    // Benchmarks intentionally use the raw byte API to measure its overhead.
    "benches",
];

/// 2-byte payloads that are known-good across the workspace.
///
/// Each entry is annotated with the reason it belongs here.
const PAYLOAD_ALLOWLIST: &[&[u8; 2]] = &[
    // ---- SAM header sub-fields (noodles `other_fields` map key convention) ----
    b"SO", b"GO", b"SS", // ---- Read-name fixtures ----
    b"q1", b"q2", b"r0", b"r1", b"r2", b"r3", b"r4", b"rd",
    // ---- 2-bp DNA sequence fragments used in sequence/UMI tests ----
    b"AA", b"AC", b"AG", b"AT", b"CA", b"CC", b"CG", b"GA", b"GC", b"GT", b"NN", b"TA", b"TC",
    b"TT",
    // ---- Opaque test-fixture tags (no SamTag constant) ----
    // X-prefixed tags used in SAM/BAM parser exercisers.
    b"XA", b"XC", b"XF", b"XH", b"XI", b"XN", b"XU", b"XV", b"XX", b"XY", b"X0", b"X1", b"X2",
    b"X3", b"X4", b"X5", b"Xc", b"Xi", b"Xs", b"Xz", b"x1", b"x2", b"x3", b"x4", b"x5", b"x6",
    b"x7", b"x8", b"x9", b"YF", b"YI", b"ZZ", b"Zz", b"BB", b"CC", b"AF", b"AF", b"As", b"aa",
    b"bb", b"fa", b"fd", b"hi", b"id", b"pa", b"s1", b"sb", b"sc", b"sq", b"uI", b"a0", b"a1",
    b"a9", b"u1", b"xy",
];

/// Path-scoped allowlist: payloads that are only exempt within specific files.
///
/// Keys are workspace-relative file paths; values are the payloads exempt in that file.
const PATH_SCOPED_ALLOWLIST: &[(&str, &[&[u8; 2]])] = &[
    // compat-proof tests in raw-bam: explicitly verify b"cd"/b"RX" still compile alongside
    // SamTag::CD_BASES / SamTag::RX, proving backward-compatible byte API still works.
    ("crates/fgumi-raw-bam/src/builder.rs", &[b"cd", b"RX"]),
    ("crates/fgumi-raw-bam/src/tags.rs", &[b"RX"]),
    // BGZF extra-subfield ID field — not a SAM aux tag.
    ("crates/fgumi-bgzf/src/writer.rs", &[b"BC"]),
];

/// Scan the entire workspace for bare 2-byte SAM-tag byte literals outside the allowlist.
///
/// `root` must be the workspace root (the directory containing the top-level `Cargo.toml`).
pub fn scan_workspace(root: &Path) -> Result<Vec<Finding>> {
    let mut findings = Vec::new();
    for path in walk_rust_sources(root)? {
        let rel = path.strip_prefix(root).unwrap_or(&path).to_path_buf();

        // File-level allowlist check.
        if FILE_ALLOWLIST.iter().any(|p| rel == Path::new(p)) {
            continue;
        }

        // Directory-prefix allowlist check.
        if DIR_ALLOWLIST.iter().any(|prefix| {
            rel.components().next().and_then(|c| c.as_os_str().to_str()) == Some(*prefix)
        }) {
            continue;
        }

        // Collect any path-scoped payload exemptions for this file.
        let scoped: &[&[u8; 2]] = PATH_SCOPED_ALLOWLIST
            .iter()
            .find(|(p, _)| rel == Path::new(p))
            .map_or(&[], |(_, payloads)| *payloads);

        let body =
            fs::read_to_string(&path).with_context(|| format!("reading {}", path.display()))?;
        scan_body(&body, &rel, scoped, &mut findings);
    }
    Ok(findings)
}

/// Recursively collect all `.rs` files under `dir`, skipping generated/vendor directories.
fn walk_rust_sources(root: &Path) -> Result<Vec<PathBuf>> {
    fn recurse(dir: &Path, out: &mut Vec<PathBuf>) -> Result<()> {
        for entry in fs::read_dir(dir).with_context(|| format!("reading {}", dir.display()))? {
            let entry = entry?;
            let path = entry.path();
            if path.is_dir() {
                let name = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
                if matches!(name, "target" | ".git" | "node_modules" | ".cargo") {
                    continue;
                }
                recurse(&path, out)?;
            } else if path.extension().and_then(|s| s.to_str()) == Some("rs") {
                out.push(path);
            }
        }
        Ok(())
    }
    let mut out = Vec::new();
    recurse(root, &mut out)?;
    Ok(out)
}

/// Scanner state for the token-aware source walker.
#[derive(Clone, Copy, PartialEq, Eq)]
enum ScanState {
    Code,
    LineComment,
    BlockComment { depth: u32 },
    Str,                      // inside "..." or b"..."
    RawStr { hashes: usize }, // inside r"..." / r#"..."# / rb"..."
}

/// Scan the text body of a single source file for bare 2-byte SAM-tag byte literals.
///
/// Uses a character-level state machine to skip line comments, block comments
/// (including nested Rust `/* /* */ */`), regular string literals, raw string literals,
/// and char literals — preventing false positives from patterns like:
/// - `/* b"MI" */`
/// - `// b"RX"` (including inline)
/// - `"example b\"PG\" usage"` (inside a regular string)
/// - `r#"b"XX""#` (inside a raw string)
///
/// The byte-string literal `b"AB"` itself is NOT skipped; only its 2-byte *content*
/// is examined for tag-name patterns.
fn scan_body(body: &str, rel: &Path, scoped_allowlist: &[&[u8; 2]], out: &mut Vec<Finding>) {
    let bytes = body.as_bytes();
    let len = bytes.len();
    let mut state = ScanState::Code;
    let mut i = 0usize;
    let mut line = 1usize;

    while i < len {
        let b = bytes[i];
        if b == b'\n' {
            if let ScanState::LineComment = state {
                state = ScanState::Code;
            }
            line += 1;
            i += 1;
            continue;
        }
        match state {
            ScanState::Code => {
                i = advance_code(bytes, len, i, b, line, rel, scoped_allowlist, out, &mut state);
            }
            ScanState::LineComment => {
                i += 1;
            }
            ScanState::BlockComment { ref mut depth } => {
                if b == b'/' && i + 1 < len && bytes[i + 1] == b'*' {
                    *depth += 1;
                    i += 2;
                } else if b == b'*' && i + 1 < len && bytes[i + 1] == b'/' {
                    *depth -= 1;
                    if *depth == 0 {
                        state = ScanState::Code;
                    }
                    i += 2;
                } else {
                    i += 1;
                }
            }
            ScanState::Str => {
                if b == b'\\' {
                    i += 2; // skip escape sequence
                } else if b == b'"' {
                    state = ScanState::Code;
                    i += 1;
                } else {
                    i += 1;
                }
            }
            ScanState::RawStr { hashes } => {
                if b == b'"' {
                    let end_hashes = bytes[i + 1..].iter().take_while(|&&c| c == b'#').count();
                    if end_hashes >= hashes {
                        state = ScanState::Code;
                        i += 1 + hashes;
                        continue;
                    }
                }
                i += 1;
            }
        }
    }
}

/// Advance one step while in `Code` state, returning the new cursor position.
///
/// Handles comment starters, raw/byte/regular string openers, char literals, and
/// byte-string tag detection. Updates `state` in place when transitioning out of Code.
#[allow(clippy::too_many_arguments)]
fn advance_code(
    bytes: &[u8],
    len: usize,
    i: usize,
    b: u8,
    line: usize,
    rel: &Path,
    scoped_allowlist: &[&[u8; 2]],
    out: &mut Vec<Finding>,
    state: &mut ScanState,
) -> usize {
    // Line comment: //
    if b == b'/' && i + 1 < len && bytes[i + 1] == b'/' {
        *state = ScanState::LineComment;
        return i + 2;
    }
    // Block comment: /*
    if b == b'/' && i + 1 < len && bytes[i + 1] == b'*' {
        *state = ScanState::BlockComment { depth: 1 };
        return i + 2;
    }
    // Raw / raw-byte string literal: r"..." r#"..."# br"..." rb"..."
    let is_raw_prefix = (b == b'r') || (b == b'b' && i + 1 < len && bytes[i + 1] == b'r');
    if is_raw_prefix {
        let after_r = if b == b'r' { i + 1 } else { i + 2 };
        let mut j = after_r;
        let mut hashes = 0usize;
        while j < len && bytes[j] == b'#' {
            hashes += 1;
            j += 1;
        }
        if j < len && bytes[j] == b'"' {
            *state = ScanState::RawStr { hashes };
            return j + 1;
        }
    }
    // Byte string literal: b"..."
    if b == b'b' && i + 1 < len && bytes[i + 1] == b'"' {
        maybe_push_finding(bytes, len, i, line, rel, scoped_allowlist, out);
        *state = ScanState::Str;
        return i + 2; // skip `b"`
    }
    // Regular string literal: "..."
    if b == b'"' {
        *state = ScanState::Str;
        return i + 1;
    }
    // Char literal: '...' — can't hold b"XX", skip it inline
    if b == b'\'' {
        let mut j = i + 1;
        while j < len && bytes[j] != b'\'' && bytes[j] != b'\n' {
            if bytes[j] == b'\\' {
                j += 1; // skip escaped byte
            }
            j += 1;
        }
        return if j < len && bytes[j] == b'\'' { j + 1 } else { j };
    }
    i + 1
}

/// If the byte-string starting at `bytes[i]` (i.e. `b"..`) holds a 2-byte SAM tag that is
/// not in either allowlist, push a finding to `out`.
fn maybe_push_finding(
    bytes: &[u8],
    len: usize,
    i: usize,
    line: usize,
    rel: &Path,
    scoped_allowlist: &[&[u8; 2]],
    out: &mut Vec<Finding>,
) {
    let after = &bytes[i + 2..];
    if after.len() >= 3
        && after[2] == b'"'
        && is_first_tag_byte(after[0])
        && is_second_tag_byte(after[1])
    {
        let payload = [after[0], after[1]];
        let in_global = PAYLOAD_ALLOWLIST.iter().any(|t| **t == payload);
        let in_scoped = scoped_allowlist.iter().any(|t| **t == payload);
        if !in_global && !in_scoped {
            let line_start = bytes[..i].iter().rposition(|&c| c == b'\n').map_or(0, |p| p + 1);
            let line_end = bytes[i..].iter().position(|&c| c == b'\n').map_or(len, |p| i + p);
            let snippet = std::str::from_utf8(&bytes[line_start..line_end])
                .unwrap_or("")
                .trim_start()
                .to_string();
            out.push(Finding { path: rel.to_path_buf(), line, snippet, payload });
        }
    }
}

/// Returns `true` if `b` is valid as the first byte of a SAM aux-tag name.
#[inline]
fn is_first_tag_byte(b: u8) -> bool {
    b.is_ascii_alphabetic()
}

/// Returns `true` if `b` is valid as the second byte of a SAM aux-tag name.
#[inline]
fn is_second_tag_byte(b: u8) -> bool {
    b.is_ascii_alphanumeric()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Returns the workspace root (two levels above this crate's manifest directory).
    fn workspace_root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .canonicalize()
            .expect("workspace root should exist")
    }

    /// Convenience: scan a source snippet with no scoped allowlist.
    fn scan(src: &str) -> Vec<[u8; 2]> {
        let mut out = Vec::new();
        scan_body(src, Path::new("fake/src.rs"), &[], &mut out);
        out.into_iter().map(|f| f.payload).collect()
    }

    #[test]
    fn no_bare_sam_tag_byte_literals_in_workspace() {
        let root = workspace_root();
        let findings = scan_workspace(&root).expect("scanner should run without I/O error");
        if !findings.is_empty() {
            for f in &findings {
                eprintln!(
                    "{}:{}  {} (payload=b\"{}{}\")",
                    f.path.display(),
                    f.line,
                    f.snippet,
                    f.payload[0] as char,
                    f.payload[1] as char,
                );
            }
        }
        assert!(
            findings.is_empty(),
            "found {} bare SAM-tag byte literal(s) outside the allowlist (see stderr above)",
            findings.len(),
        );
    }

    // ------------------------------------------------------------------
    // Token-aware scanner unit tests
    // ------------------------------------------------------------------

    #[test]
    fn scanner_detects_bare_tag_in_code() {
        assert_eq!(scan("let x = b\"MI\";"), vec![[b'M', b'I']]);
    }

    #[test]
    fn scanner_skips_full_line_comments() {
        assert!(scan("// let x = b\"RG\";").is_empty(), "full-line comment should be skipped");
    }

    #[test]
    fn scanner_skips_inline_comments() {
        let src = "fn f() { let x = 1; /* b\"MI\" */ let y = 2; } // b\"RG\"";
        assert!(scan(src).is_empty(), "inline comments should be skipped");
    }

    #[test]
    fn scanner_skips_block_comments() {
        assert!(scan("/* b\"PG\" */\nlet z = 0;").is_empty(), "block comment should be skipped");
    }

    #[test]
    fn scanner_handles_nested_block_comments() {
        let src = "/* outer /* b\"PG\" */ still comment */ let z = 0;";
        assert!(scan(src).is_empty(), "nested block comment should be skipped");
    }

    #[test]
    fn scanner_skips_regular_string_contents() {
        // b"MI" appears inside a regular string literal — not a real byte-string.
        let src = r#"let s = "example b\"MI\" usage";"#;
        assert!(scan(src).is_empty(), "regular string contents should be skipped");
    }

    #[test]
    fn scanner_skips_raw_string_contents() {
        // b"MI" inside a raw string should not be flagged.
        let src = "let s = r#\"b\\\"MI\\\"\"#;";
        assert!(scan(src).is_empty(), "raw string contents should be skipped");
    }

    #[test]
    fn scanner_finds_tag_after_comment_block_ends() {
        let src = "/* comment */ let x = b\"MI\";";
        assert_eq!(scan(src), vec![[b'M', b'I']]);
    }

    #[test]
    fn scanner_respects_scoped_allowlist() {
        let src = "let x = b\"RX\";";
        let mut out = Vec::new();
        scan_body(src, Path::new("crates/fgumi-raw-bam/src/builder.rs"), &[b"RX"], &mut out);
        assert!(out.is_empty(), "scoped allowlist should suppress b\"RX\" in builder.rs");
    }

    #[test]
    fn scanner_flags_when_not_in_scoped_allowlist() {
        // b"RX" is NOT in the global allowlist; only in the scoped allowlist for builder.rs.
        // Other files should still flag it.
        let src = "let x = b\"RX\";";
        let mut out = Vec::new();
        scan_body(src, Path::new("src/lib/some_other_file.rs"), &[], &mut out);
        assert!(!out.is_empty(), "b\"RX\" should be flagged outside the scoped file");
    }
}
