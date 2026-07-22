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
/// BGZF is gzip with a mandatory `BC` extra subfield, which lives at bytes
/// 12-13 behind the `XLEN` at 10-11. Anything shorter can only answer "gzip or
/// not", which is how a plain-gzipped file gets mistaken for BAM.
pub const FORMAT_PREFIX_LEN: usize = 18;

/// The deflate compression method, the only one BGZF uses.
const DEFLATE: u8 = 0x08;

/// `FEXTRA`, and nothing else: BGZF sets exactly this one gzip flag.
const FEXTRA_ONLY: u8 = 0x04;

/// `XLEN`: the extra field is exactly the 6-byte `BC` subfield, no more.
const BGZF_XLEN: [u8; 2] = [0x06, 0x00];

/// The `BC` subfield identifier, which is what makes a gzip member BGZF.
const BC_SUBFIELD_ID: [u8; 2] = [b'B', b'C'];

/// `SLEN` for `BC`: two bytes, holding the total block size minus one.
const BC_SUBFIELD_LEN: [u8; 2] = [0x02, 0x00];

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
    /// [`classify_input`]). Either way fgumi's readers cannot consume it even
    /// though it decompresses, and the remedy is the same: recompress with
    /// `bgzip`.
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
/// # The BGZF check is at fixed offsets, deliberately
///
/// RFC 1952 permits several extra subfields per gzip member in any order, so
/// `BC` need not be first and a spec-literal check would have to walk the
/// `FEXTRA` area looking for it. This does not, and must not: the verdict's job
/// is to predict what the *decoder* will accept, and every decoder downstream of
/// it is keyed to these same fixed offsets.
///
/// * `fgumi_bgzf`'s `read_raw_block` (behind `read_raw_blocks`) — the sort
///   worker pool, the BAM writers, `fastq` — checks `BC` at bytes 12-13 and
///   takes `BSIZE` from bytes 16-17 to *frame* the block. The layout is
///   load-bearing for where one block ends and the next begins, not merely for
///   validation, and the deflate payload is taken as `[18 .. len - 8]`.
/// * `noodles_bgzf`'s `is_valid_header` additionally requires `FLG` to be
///   `FEXTRA` alone and `XLEN` to be exactly 6.
/// * htslib's `check_header` (`bgzf.c`) requires the same `XLEN == 6`, `BC` at
///   12-13, and `SLEN == 2`; it differs from noodles only in tolerating other
///   `FLG` bits alongside `FEXTRA`.
///
/// bgzip, samtools, htsjdk and noodles all *write* exactly this layout, so no
/// BGZF a user can produce with them is missed.
///
/// Answering `Bgzf` for a member that misses this layout would therefore hand it
/// to a decoder that rejects it — `invalid BGZF header` from noodles, or
/// `Invalid BGZF subfield ID` from `fgumi_bgzf` — which is the misleading
/// diagnosis this classifier exists to prevent. Answering `Gzip` names it as
/// gzip-that-is-not-BGZF and points at `bgzip`, which is the actual remedy — the
/// member has to be recompressed before fgumi can read it either way. For FASTQ,
/// where [`InputFormat::Gzip`] is a supported input rather than an error, the
/// verdict costs only block-parallel decode: `extract` reads it with
/// `MultiGzDecoder`, and BGZF is valid gzip.
#[must_use]
pub fn classify_input(prefix: &[u8]) -> InputFormat {
    if prefix.is_empty() {
        return InputFormat::Empty;
    }
    if !prefix.starts_with(&GZIP_MAGIC) {
        return InputFormat::Text;
    }

    // gzip. BGZF is gzip whose sole extra subfield is `BC`, at a fixed offset —
    // see this function's doc comment for why the offsets are not walked.
    if prefix.len() >= FORMAT_PREFIX_LEN
        && prefix[2] == DEFLATE
        && prefix[3] == FEXTRA_ONLY
        && prefix[10..12] == BGZF_XLEN
        && prefix[12..14] == BC_SUBFIELD_ID
        && prefix[14..16] == BC_SUBFIELD_LEN
    {
        InputFormat::Bgzf
    } else {
        InputFormat::Gzip
    }
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
    // FEXTRA set but the subfield is not `BC` — some other gzip extension.
    #[case::gzip_with_non_bc_extra(
        &[0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'X', b'Y', 0x02, 0x00, 0x1b, 0x00],
        InputFormat::Gzip
    )]
    // `BC` first but XLEN 10, so a second subfield trails it. RFC 1952 allows
    // this; `noodles_bgzf` and htslib both reject it (XLEN must be exactly 6),
    // so calling it BGZF would promise a decode that cannot happen — see
    // `classify_input`'s doc comment.
    #[case::bc_first_with_trailing_subfield(
        &[0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x0a, 0x00, b'B', b'C', 0x02, 0x00, 0x1b, 0x00],
        InputFormat::Gzip
    )]
    // A foreign subfield ahead of `BC`: same reasoning, and here the fixed
    // offsets cannot see `BC` at all within the 18-byte window.
    #[case::foreign_subfield_before_bc(
        &[0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x0a, 0x00, b'X', b'Y', 0x02, 0x00, 0x00, 0x00],
        InputFormat::Gzip
    )]
    // FEXTRA plus FNAME. The extra field precedes FNAME in a gzip member, so
    // `XLEN` and `BC` are still where they belong and htslib would take this
    // (`header[3] & 4`); `noodles_bgzf` requires FLG to be FEXTRA alone and
    // rejects it, and this verdict tracks our decoder rather than the spec.
    #[case::fextra_with_fname(
        &[0x1f, 0x8b, 0x08, 0x0c, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00, 0x1b, 0x00],
        InputFormat::Gzip
    )]
    // `BC` at the right offset but SLEN != 2, so the block size it carries is
    // not the two bytes every decoder reads.
    #[case::bc_with_wrong_slen(
        &[0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x04, 0x00, 0x1b, 0x00],
        InputFormat::Gzip
    )]
    // Compression method other than deflate: gzip's framing, nothing fgumi reads.
    #[case::non_deflate_method(
        &[0x1f, 0x8b, 0x00, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00, 0x1b, 0x00],
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
