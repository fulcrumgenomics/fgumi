//! Input-format detection shared by every fgumi reader.
//!
//! Detection is by content, never by file extension: one classifier so a given
//! file gets the same verdict — and so the same diagnosis — whichever command
//! opens it. Three independent sniffs (2, 4 and 18 bytes wide) previously
//! disagreed, so a plain-gzipped SAM was read as BAM by one command and as SAM
//! text by another, each failing with a different and misleading error.

/// gzip magic, which every BGZF block starts with.
///
/// SAM text cannot begin with these bytes: a SAM file starts either with `@`
/// (a header line) or with the first character of a read name, and the SAM
/// specification restricts read names to printable ASCII.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Bytes needed to tell BGZF from plain gzip.
///
/// A whole BGZF block header, because that is what the shared predicate reads:
/// BGZF is gzip with a mandatory `BC` extra subfield behind the `XLEN` at 10-11.
/// Anything shorter can only answer "gzip or not", which is how a plain-gzipped
/// file gets mistaken for BAM.
pub const FORMAT_PREFIX_LEN: usize = fgumi_bgzf::BGZF_HEADER_SIZE;

/// What an input's leading bytes say it is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
    /// BGZF — a BAM, or any other bgzip-compressed stream.
    Bgzf,
    /// gzip, but not BGZF that fgumi can block-decode.
    ///
    /// Plain gzip — one deflate stream, so not seekable or block-addressable —
    /// and also the rare gzip member that carries a `BC` subfield but does not
    /// match the fixed BGZF layout every fgumi decoder requires (see
    /// [`classify_input`]).
    ///
    /// This is a supported alignment input, not an error: the normalization
    /// boundary decompresses the member and re-frames what comes out, so the rest
    /// of the pipeline still sees only BGZF. What it costs is block-parallel
    /// decode — one deflate stream cannot be split across threads — so `bgzip` is
    /// a performance recommendation here rather than a requirement.
    Gzip,
    /// Uncompressed text; for alignment inputs, SAM.
    Text,
    /// No bytes at all.
    Empty,
}

/// Classify an input from its leading bytes.
///
/// `prefix` should hold up to [`FORMAT_PREFIX_LEN`] bytes (fewer only at end of
/// input). One classifier for every fgumi input sniff, so a given file gets the
/// same verdict — and so the same diagnosis — whichever command opens it.
///
/// Detection is by content, never by file extension.
///
/// # The BGZF verdict comes from the decoder's own predicate
///
/// The `Bgzf` arm delegates to [`fgumi_bgzf::is_bgzf_header`], which is the same
/// check [`fgumi_bgzf::reader`] applies to every block it frames. That shared
/// definition is the point: the verdict's job is to predict what the decoder will
/// accept, so answering it with a *second*, hand-rolled copy of the layout is how
/// the two came to disagree in the first place. See [`fgumi_bgzf::header`] for
/// why the layout is matched at fixed offsets rather than walked, and why the
/// predicate is the intersection of what fgumi's decoders accept rather than what
/// RFC 1952 permits.
///
/// Answering `Bgzf` for a member that misses this layout would hand it to a
/// decoder that rejects it — `invalid BGZF header` from noodles, or a named field
/// rejection from `fgumi_bgzf` — which is the misleading diagnosis this
/// classifier exists to prevent. Answering `Gzip` names it as
/// gzip-that-is-not-BGZF, which routes it to a spec-complete gzip decoder rather
/// than to a BGZF block reader that would reject it.
///
/// `Gzip` is a supported verdict on every input path, not an error. `extract`
/// reads FASTQ with `MultiGzDecoder`, and alignment inputs go through
/// `normalize_to_bgzf`, which decompresses the member and re-frames the result.
/// In both cases the only cost is block-parallel decode — one deflate stream
/// cannot be split across threads — so `bgzip` buys throughput, not readability.
#[must_use]
pub fn classify_input(prefix: &[u8]) -> InputFormat {
    if prefix.is_empty() {
        return InputFormat::Empty;
    }
    if !prefix.starts_with(&GZIP_MAGIC) {
        return InputFormat::Text;
    }

    // gzip. Whether it is BGZF is the block reader's question, so it answers it.
    if fgumi_bgzf::is_bgzf_header(prefix) { InputFormat::Bgzf } else { InputFormat::Gzip }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    /// A valid 18-byte BGZF header: gzip magic, deflate, FEXTRA, then the
    /// mandatory `BC` subfield that distinguishes BGZF from plain gzip.
    const BGZF: [u8; 18] = [
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00, 0x1b, 0x00,
    ];

    /// gzip magic and deflate, but no FEXTRA and no `BC` subfield.
    const PLAIN_GZIP: [u8; 18] =
        [0x1f, 0x8b, 0x08, 0x00, 0, 0, 0, 0, 0, 0xff, 0, 0, 0, 0, 0, 0, 0, 0];

    #[rstest]
    #[case::bgzf(&BGZF, InputFormat::Bgzf)]
    #[case::plain_gzip(&PLAIN_GZIP, InputFormat::Gzip)]
    #[case::sam_header(b"@HD\tVN:1.6\tSO:queryname\n", InputFormat::Text)]
    #[case::empty(b"", InputFormat::Empty)]
    // Too short to hold the BGZF extra field: gzip magic is all we can confirm,
    // and calling it BGZF would hand a truncated file to the block decoder.
    #[case::gzip_magic_but_truncated(&[0x1f, 0x8b, 0x08, 0x04], InputFormat::Gzip)]
    // A short *text* prefix is still unambiguously text.
    #[case::short_text(b"@H", InputFormat::Text)]
    // The per-field near misses — a foreign subfield, `BC` with a trailing
    // subfield, FEXTRA plus FNAME, a wrong SLEN, a non-deflate method — are the
    // shared predicate's contract and are pinned in `fgumi_bgzf::header`'s own
    // table. Duplicating them here would be a second copy of exactly the thing
    // this delegation removed. What belongs here is the dispatch: that a gzip
    // member which is not a BGZF block reads as `Gzip` rather than `Bgzf`.
    #[case::gzip_with_non_bc_extra(
        &[0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'X', b'Y', 0x02, 0x00, 0x1b, 0x00],
        InputFormat::Gzip
    )]
    fn classifies(#[case] prefix: &[u8], #[case] expected: InputFormat) {
        assert_eq!(classify_input(prefix), expected);
    }

    /// The `Bgzf` verdict must agree with the decoder that consumes the bytes.
    ///
    /// `classify_input` reads a fixed-offset window rather than walking the
    /// `FEXTRA` subfields, which is only sound while the decoder is at least as
    /// strict. This pins that: a real BGZF block, framed by the same writer
    /// fgumi uses, must classify as `Bgzf` and decode.
    #[test]
    fn bgzf_verdict_agrees_with_the_decoder() {
        use std::io::Write;

        let mut writer = noodles::bgzf::io::Writer::new(Vec::new());
        writer.write_all(b"payload").expect("write BGZF payload");
        let block = writer.finish().expect("finish BGZF stream");

        assert_eq!(
            classify_input(&block[..FORMAT_PREFIX_LEN]),
            InputFormat::Bgzf,
            "a block written by the BGZF writer must classify as BGZF"
        );

        let mut decoded = Vec::new();
        std::io::copy(&mut noodles::bgzf::io::Reader::new(&block[..]), &mut decoded)
            .expect("a block classified as BGZF must decode as BGZF");
        assert_eq!(decoded, b"payload");
    }
}
