//! Aligner-header synthesis for the in-process bwa-mem3 backend.
//!
//! The subprocess backend reads the aligner's emitted SAM/BAM header off its
//! stdout and merges it into the partial output header. The in-process backend
//! has no such stream, so it **synthesizes** the equivalent header from the
//! loaded index's contigs plus one `@PG ID:bwa-mem3` line, then feeds it through
//! the *same* [`validate_sq_consistency`](super::super::validate_sq_consistency)
//! and [`merge_aligner_header`](super::super::merge_aligner_header) path the
//! subprocess reader uses, so both backends share one `@SQ` check.
//!
//! What survives that merge: the synthesized `@SQ` table is validated against
//! the dict-derived partial header and then discarded (the partial's is
//! authoritative); the `@RG`/`@PG`/`@CO` lines of bwa-mem3's optional
//! per-index header sidecar (`<prefix>.hdr`, else `<baseprefix>.dict` — see
//! [`load_index_header_sidecar`]), which the CLI copies into its own header,
//! are merged exactly as the subprocess reader merges them; and the single
//! `@PG ID:bwa-mem3` line is appended. The `@PG CL` is deliberately **honest**
//! about being an in-process invocation rather than forging the subprocess
//! command line. That `CL` is the one header difference from the subprocess
//! preset; `tests/align_inproc_parity.rs` also normalizes the `@PG VN` git
//! suffix and each record's optional-tag order, and compares decoded BAM
//! records rather than raw bytes.

use std::ffi::OsString;
use std::io::{self, Read};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use noodles::sam::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::program::tag as pg_tag;
use noodles::sam::header::record::value::map::{Program, ReferenceSequence};

use super::engine;

/// The `@PG ID` of the synthesized bwa-mem3 program line, mirroring the CLI's
/// `@PG\tID:bwa-mem3` (`main.cpp:243`).
const BWA_MEM3_PG_ID: &str = "bwa-mem3";

/// The candidate sidecar paths bwa-mem3 tries for index `prefix`, in order:
/// `<prefix>.hdr`, then `<baseprefix>.dict`, where `baseprefix` drops a trailing
/// `.gz` and then the final dotted suffix inside the basename (`foo.fa` →
/// `foo`, `foo.fasta.gz` → `foo`, `GRCh38.p14` → `GRCh38`; no dot → unchanged).
/// An exact port of the path logic in bwa-mem3's `bwa_load_hdr_from_index`
/// (`bwa.cpp`).
fn index_header_sidecar_paths(prefix: &Path) -> [PathBuf; 2] {
    let mut hdr = OsString::from(prefix.as_os_str());
    hdr.push(".hdr");

    // A non-UTF-8 prefix falls back to a lossy copy for the `.dict` candidate.
    let text = prefix.to_string_lossy();
    let mut base: &str = text.strip_suffix(".gz").unwrap_or(&text);
    // Cut at the last '.' that comes after the last '/'.
    if let Some(i) = base.rfind(['.', '/'])
        && base.as_bytes()[i] == b'.'
    {
        base = &base[..i];
    }
    [PathBuf::from(hdr), PathBuf::from(format!("{base}.dict"))]
}

/// Load the `@RG`/`@PG`/`@CO` lines of bwa-mem3's optional per-index header
/// sidecar for index `prefix`, in file order, or `None` when neither sidecar
/// file exists (or it is empty).
///
/// The bwa-mem3 CLI reads `<prefix>.hdr`, else `<baseprefix>.dict`
/// ([`index_header_sidecar_paths`]), strips `\r`, and copies every line into
/// its output header ahead of its own `@PG`. The subprocess backend then merges
/// that header's `@RG`/`@PG`/`@CO` into the output
/// ([`merge_aligner_header`](super::super::merge_aligner_header)); these are the
/// only sidecar lines that survive the merge, so they are the only ones kept.
///
/// # Errors
/// Reading an existing sidecar fails, or its kept lines do not parse as SAM
/// header records.
pub(crate) fn load_index_header_sidecar(prefix: &Path) -> io::Result<Option<Header>> {
    // Like `fopen` in the CLI: a path that cannot be opened is simply absent.
    let Some(mut file) =
        index_header_sidecar_paths(prefix).into_iter().find_map(|p| std::fs::File::open(p).ok())
    else {
        return Ok(None);
    };
    let mut raw = Vec::new();
    file.read_to_end(&mut raw)?;
    raw.retain(|&b| b != b'\r');

    let mut kept = String::new();
    for line in raw.split(|&b| b == b'\n') {
        if line.starts_with(b"@RG\t") || line.starts_with(b"@PG\t") || line.starts_with(b"@CO\t") {
            kept.push_str(&String::from_utf8_lossy(line));
            kept.push('\n');
        }
    }
    if kept.is_empty() {
        return Ok(None);
    }
    kept.parse::<Header>().map(Some).map_err(|e| {
        io::Error::other(format!(
            "in-process bwa-mem3: cannot parse the @RG/@PG/@CO lines of the index header \
             sidecar for '{prefix}': {e}",
            prefix = prefix.display(),
        ))
    })
}

/// Synthesize the aligner-emitted header from the index's `contigs`, the
/// index header `sidecar`'s `@RG`/`@PG`/`@CO` lines (see
/// [`load_index_header_sidecar`]), and one `@PG ID:bwa-mem3` line.
///
/// Takes the contigs as `(name, length)` pairs — an iterator, not the whole
/// [`BwaIndex`](engine::BwaIndex) — so it is unit-testable with a hand-built
/// stand-in and never needs a real index. In production the caller passes
/// `idx.contigs()` (which yields `(&str, i64)`).
///
/// The result carries:
/// - one `@SQ` per contig, in index order, name + length (the only fields the
///   subprocess aligner would emit and the merge preserves after validation);
/// - the sidecar's `@RG`, `@PG` and `@CO` lines, in file order, with its `@PG`s
///   ahead of bwa-mem3's own, as the CLI writes them;
/// - one `@PG ID:bwa-mem3` line with `PN:bwa-mem3`, `VN` = the linked bwa-mem3
///   version ([`version`](engine::version)), and an honest in-process `CL`.
///
/// The `@HD` is intentionally omitted:
/// [`merge_aligner_header`](super::super::merge_aligner_header) takes the partial
/// header's `@HD`, never the aligner's, so synthesizing one would be
/// dead weight.
pub(crate) fn synthesize_aligner_header<'a, C>(
    contigs: C,
    sidecar: Option<&Header>,
    chunk_size: u64,
    reference: &Path,
) -> Header
where
    C: IntoIterator<Item = (&'a str, i64)>,
{
    let mut builder = Header::builder();

    for (name, length) in contigs {
        // A real index contig length is always >= 1; clamp defensively so the
        // `NonZeroUsize` conversion is total. A degenerate length would fail
        // `validate_sq_consistency` against the dict anyway (that mismatch is the
        // point), so clamping cannot mask a real inconsistency.
        let len = NonZeroUsize::new(usize::try_from(length).unwrap_or(1).max(1))
            .expect("clamped length is >= 1");
        builder = builder
            .add_reference_sequence(bstr::BString::from(name), Map::<ReferenceSequence>::new(len));
    }

    if let Some(sidecar) = sidecar {
        for (id, rg) in sidecar.read_groups() {
            builder = builder.add_read_group(id.clone(), rg.clone());
        }
        for (id, pg) in sidecar.programs().as_ref() {
            builder = builder.add_program(id.clone(), pg.clone());
        }
        for comment in sidecar.comments() {
            builder = builder.add_comment(comment.clone());
        }
    }

    let version = engine::version();
    // Honest in-process command line: describes the invocation shape
    // (`mem -p -K <K>` against the reference) and that it ran in-process via
    // bwa-mem3-rs, rather than forging the subprocess `bwa-mem3 mem ... /dev/stdin`.
    let command_line = format!(
        "bwa-mem3 mem -p -K {chunk_size} {reference} (in-process via bwa-mem3-rs {version})",
        reference = reference.display(),
    );
    let program = Map::<Program>::builder()
        .insert(pg_tag::NAME, BWA_MEM3_PG_ID)
        .insert(pg_tag::VERSION, version)
        .insert(pg_tag::COMMAND_LINE, command_line)
        .build()
        .expect("synthesized @PG map is valid");
    builder = builder.add_program(bstr::BString::from(BWA_MEM3_PG_ID), program);

    builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::steps::align::{merge_aligner_header, validate_sq_consistency};
    use noodles::sam::header::record::value::map::Header as HeaderMap;

    /// Build a dict-style partial header with an `@HD` and the given `@SQ`s.
    fn partial_header(refs: &[(&str, usize)]) -> Header {
        let mut b = Header::builder().set_header(Map::<HeaderMap>::default());
        for (name, length) in refs {
            let len = NonZeroUsize::new(*length).expect("nonzero ref length");
            b = b.add_reference_sequence(
                bstr::BString::from(*name),
                Map::<ReferenceSequence>::new(len),
            );
        }
        b.build()
    }

    /// The `PN`/`VN`/`CL` fields of the `@PG` with the given ID, if present.
    fn pg_fields(header: &Header, id: &str) -> Option<(String, String, String)> {
        let programs = header.programs();
        let pg = programs.as_ref().get(&bstr::BString::from(id))?;
        let other = pg.other_fields();
        Some((
            other.get(&pg_tag::NAME)?.to_string(),
            other.get(&pg_tag::VERSION)?.to_string(),
            other.get(&pg_tag::COMMAND_LINE)?.to_string(),
        ))
    }

    // (a) synthesize builds @SQ from the contigs + one @PG ID:bwa-mem3.
    #[test]
    fn synthesize_builds_sq_from_contigs_and_one_bwa_mem3_pg() {
        let contigs = vec![("chr1", 1000i64), ("chr2", 2000i64)];
        let synth =
            synthesize_aligner_header(contigs, None, 150_000_000, Path::new("/refs/genome.fa"));

        // @SQ built from the contigs, in order, name + length.
        let sq: Vec<(String, usize)> = synth
            .reference_sequences()
            .iter()
            .map(|(name, map)| (name.to_string(), map.length().get()))
            .collect();
        assert_eq!(sq, vec![("chr1".to_string(), 1000), ("chr2".to_string(), 2000)]);

        // Exactly one @PG, ID bwa-mem3, PN bwa-mem3, VN = the linked version, and
        // an honest in-process CL.
        assert_eq!(synth.programs().as_ref().len(), 1, "exactly one synthesized @PG");
        let (pn, vn, cl) = pg_fields(&synth, "bwa-mem3").expect("@PG ID:bwa-mem3 present");
        assert_eq!(pn, "bwa-mem3");
        assert_eq!(vn, engine::version(), "VN is the linked bwa-mem3 version");
        assert_eq!(
            cl,
            format!(
                "bwa-mem3 mem -p -K 150000000 /refs/genome.fa (in-process via bwa-mem3-rs {})",
                engine::version()
            ),
            "CL is the honest in-process invocation"
        );
    }

    // (b) merge keeps the partial @HD/@SQ and appends the synthesized @PG.
    #[test]
    fn merge_keeps_partial_hd_sq_and_appends_the_pg() {
        let partial = partial_header(&[("chr1", 1000), ("chr2", 2000)]);
        let synth = synthesize_aligner_header(
            vec![("chr1", 1000i64), ("chr2", 2000i64)],
            None,
            150_000_000,
            Path::new("/refs/genome.fa"),
        );
        validate_sq_consistency(&partial, &synth).expect("matching @SQ validates");
        let merged = merge_aligner_header(&partial, &synth);

        // Partial's @HD survives.
        assert!(merged.header().is_some(), "merged keeps the partial @HD");
        // Partial's (dict) @SQ survives unchanged.
        let sq: Vec<(String, usize)> = merged
            .reference_sequences()
            .iter()
            .map(|(name, map)| (name.to_string(), map.length().get()))
            .collect();
        assert_eq!(sq, vec![("chr1".to_string(), 1000), ("chr2".to_string(), 2000)]);
        // The synthesized @PG is appended.
        let (pn, _, _) = pg_fields(&merged, "bwa-mem3").expect("@PG bwa-mem3 appended");
        assert_eq!(pn, "bwa-mem3");
    }

    // (c) an @SQ mismatch (synth built from different contigs than the partial's
    // @SQ) fails validation with the existing message — the check that runs at
    // wire time before any input is read.
    #[test]
    fn sq_mismatch_fails_validation_at_wire_time() {
        let partial = partial_header(&[("chr1", 1000), ("chr2", 2000)]);
        // Synthesize from an index built against a DIFFERENT FASTA (wrong length).
        let synth = synthesize_aligner_header(
            vec![("chr1", 1000i64), ("chr2", 9999i64)],
            None,
            150_000_000,
            Path::new("/refs/other.fa"),
        );
        let err = validate_sq_consistency(&partial, &synth)
            .expect_err("a length mismatch must reject at wire time");
        let msg = err.to_string();
        assert!(msg.contains("length mismatch"), "message: {msg}");
        assert!(msg.contains("chr2"), "message names the mismatched contig: {msg}");
    }

    // (d) the sidecar paths bwa-mem3 tries: `<prefix>.hdr`, then
    // `<baseprefix>.dict` with a trailing `.gz` and the last dotted suffix of the
    // basename dropped.
    #[rstest::rstest]
    #[case::fasta("/refs/genome.fa", "/refs/genome.fa.hdr", "/refs/genome.dict")]
    #[case::gzipped("/refs/genome.fasta.gz", "/refs/genome.fasta.gz.hdr", "/refs/genome.dict")]
    #[case::patch_suffix("/refs/GRCh38.p14", "/refs/GRCh38.p14.hdr", "/refs/GRCh38.dict")]
    #[case::no_dot_in_basename("/refs.d/genome", "/refs.d/genome.hdr", "/refs.d/genome.dict")]
    #[case::relative("genome.fa", "genome.fa.hdr", "genome.dict")]
    fn sidecar_paths_match_bwa_mem3(
        #[case] prefix: &str,
        #[case] expected_hdr: &str,
        #[case] expected_dict: &str,
    ) {
        let [hdr, dict] = index_header_sidecar_paths(Path::new(prefix));
        assert_eq!(hdr, Path::new(expected_hdr));
        assert_eq!(dict, Path::new(expected_dict));
    }

    /// Write `contents` to `dir/name`.
    fn write_file(dir: &Path, name: &str, contents: &str) {
        std::fs::write(dir.join(name), contents).expect("write sidecar");
    }

    // (e) `.hdr` wins over `.dict`; `.dict` is the fallback; neither => `None`.
    #[test]
    fn sidecar_prefers_hdr_then_falls_back_to_dict() {
        let dir = tempfile::tempdir().expect("tempdir");
        let prefix = dir.path().join("genome.fa");

        assert!(load_index_header_sidecar(&prefix).expect("no sidecar is Ok").is_none());

        write_file(dir.path(), "genome.dict", "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n@CO\tfrom dict\n");
        let dict = load_index_header_sidecar(&prefix).expect("dict loads").expect("dict present");
        assert_eq!(dict.comments(), [bstr::BString::from("from dict")]);

        write_file(dir.path(), "genome.fa.hdr", "@CO\tfrom hdr\r\n");
        let hdr = load_index_header_sidecar(&prefix).expect("hdr loads").expect("hdr present");
        assert_eq!(hdr.comments(), [bstr::BString::from("from hdr")], "CR stripped, .hdr wins");
    }

    // (f) a sidecar with only @HD/@SQ contributes nothing that survives the merge.
    #[test]
    fn sidecar_with_only_hd_and_sq_is_none() {
        let dir = tempfile::tempdir().expect("tempdir");
        write_file(dir.path(), "genome.dict", "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n");
        let loaded = load_index_header_sidecar(&dir.path().join("genome.fa")).expect("dict loads");
        assert!(loaded.is_none());
    }

    // (g) the sidecar's @RG/@PG/@CO reach the merged header, its @PG ahead of
    // bwa-mem3's, as the subprocess reader would merge them from the CLI header.
    #[test]
    fn sidecar_rg_pg_co_are_merged_before_the_bwa_mem3_pg() {
        let dir = tempfile::tempdir().expect("tempdir");
        write_file(
            dir.path(),
            "genome.dict",
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:s1\n\
             @PG\tID:picard\tPN:picard\n@CO\tbuilt by picard\n",
        );
        let prefix = dir.path().join("genome.fa");
        let sidecar = load_index_header_sidecar(&prefix).expect("dict loads");
        let synth = synthesize_aligner_header(
            vec![("chr1", 1000i64)],
            sidecar.as_ref(),
            150_000_000,
            &prefix,
        );
        let merged = merge_aligner_header(&partial_header(&[("chr1", 1000)]), &synth);

        assert!(
            merged.read_groups().contains_key(&bstr::BString::from("rg1")),
            "sidecar @RG merged"
        );
        let pg_ids: Vec<String> =
            merged.programs().as_ref().keys().map(ToString::to_string).collect();
        assert_eq!(pg_ids, vec!["picard".to_string(), "bwa-mem3".to_string()]);
        assert_eq!(merged.comments(), [bstr::BString::from("built by picard")]);
    }
}
