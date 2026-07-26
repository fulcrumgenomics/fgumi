//! The one BGZF block-header layout every fgumi reader agrees on.
//!
//! A BGZF block is a gzip member whose extra field carries a `BC` subfield
//! holding the block's total size. Reading one means answering two questions —
//! "is this a BGZF block?" and "where does it end?" — and both were previously
//! answered by hand at four separate sites, which is how they came to disagree:
//! the input sniff accepted headers the block reader would reject, and the block
//! reader accepted headers `noodles_bgzf` would reject. This module is the single
//! answer to both.
//!
//! # Why the layout is fixed rather than walked
//!
//! RFC 1952 lets a gzip member carry several extra subfields in any order, so a
//! spec-literal reader would walk the extra field looking for `BC`. fgumi does
//! not, because the layout below is not merely *validated* here — it is what
//! [`block_size`] uses to find the end of the block, and what every consumer
//! assumes when it takes the deflate payload as `[18 .. len - 8]`.
//!
//! The predicate is therefore the **intersection** of what fgumi's decoders
//! accept, not the union and not the spec: a `true` from [`is_bgzf_header`] has
//! to be a promise that every downstream path can keep. `noodles_bgzf`'s
//! `is_valid_header` is the strictest of them and is not ours to change, so it
//! sets the bar — `FLG` must be `FEXTRA` alone and `XLEN` exactly 6. htslib's
//! `check_header` agrees on `XLEN`, `BC` and `SLEN`, differing only in tolerating
//! other `FLG` bits; a file that exercises that difference is readable by
//! samtools and not by fgumi, which this predicate reports honestly rather than
//! promising a decode that would fail later.
//!
//! Nothing is lost in practice: bgzip, samtools, htsjdk and noodles all *write*
//! exactly this layout, and so does fgumi — [`crate::BGZF_EOF`] is literally
//! these bytes.

use std::fmt;

/// Size of a BGZF block header.
pub const BGZF_HEADER_SIZE: usize = 18;

/// gzip magic (ID1, ID2), which every BGZF block starts with.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

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

/// Offset of the `BSIZE` value within the header.
const BSIZE_OFFSET: usize = 16;

/// Why a candidate header is not a BGZF block header fgumi can decode.
///
/// Carried rather than rendered so the streaming block reader can name the exact
/// field that failed — "invalid BGZF header" over a truncated file and over a
/// plain-gzip file send a reader looking in different places.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HeaderRejection {
    /// Fewer than [`BGZF_HEADER_SIZE`] bytes were available.
    TooShort {
        /// How many bytes there were.
        len: usize,
    },
    /// Not a gzip member at all.
    Magic {
        /// The first two bytes.
        got: [u8; 2],
    },
    /// gzip, but not deflate-compressed.
    CompressionMethod {
        /// The `CM` byte.
        got: u8,
    },
    /// `FLG` is not `FEXTRA` alone; see the module docs for why this is exact.
    Flags {
        /// The `FLG` byte.
        got: u8,
    },
    /// `XLEN` is not 6, so the extra field is not just the `BC` subfield.
    ExtraFieldLength {
        /// The `XLEN` value.
        got: u16,
    },
    /// The extra subfield is not `BC`.
    SubfieldId {
        /// The two subfield identifier bytes.
        got: [u8; 2],
    },
    /// `BC`'s `SLEN` is not 2, so it does not hold a `BSIZE`.
    SubfieldLength {
        /// The `SLEN` value.
        got: u16,
    },
}

impl fmt::Display for HeaderRejection {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::TooShort { len } => {
                write!(f, "truncated BGZF header: expected {BGZF_HEADER_SIZE} bytes, got {len}")
            }
            Self::Magic { got } => write!(
                f,
                "Invalid BGZF magic: expected 0x1f 0x8b, got 0x{:02x} 0x{:02x}",
                got[0], got[1]
            ),
            Self::CompressionMethod { got } => {
                write!(f, "Invalid compression method: expected 0x08, got 0x{got:02x}")
            }
            Self::Flags { got } => {
                write!(f, "Invalid BGZF flags: expected FEXTRA (0x04) alone, got 0x{got:02x}")
            }
            Self::ExtraFieldLength { got } => {
                write!(f, "Invalid BGZF extra-field length: expected 6, got {got}")
            }
            Self::SubfieldId { got } => write!(
                f,
                "Invalid BGZF subfield ID: expected 'BC', got '{}{}'",
                got[0] as char, got[1] as char
            ),
            Self::SubfieldLength { got } => {
                write!(f, "Invalid BGZF subfield length: expected 2, got {got}")
            }
        }
    }
}

impl std::error::Error for HeaderRejection {}

/// Check `prefix` against the BGZF block-header layout, naming the field that
/// fails.
///
/// `prefix` may be longer than a header; only the first [`BGZF_HEADER_SIZE`]
/// bytes are read.
///
/// # Errors
///
/// Returns the first [`HeaderRejection`] that applies, in field order, so the
/// message names the earliest thing that is wrong rather than an arbitrary one.
pub fn validate(prefix: &[u8]) -> Result<(), HeaderRejection> {
    if prefix.len() < BGZF_HEADER_SIZE {
        return Err(HeaderRejection::TooShort { len: prefix.len() });
    }
    if prefix[..2] != GZIP_MAGIC {
        return Err(HeaderRejection::Magic { got: [prefix[0], prefix[1]] });
    }
    if prefix[2] != DEFLATE {
        return Err(HeaderRejection::CompressionMethod { got: prefix[2] });
    }
    if prefix[3] != FEXTRA_ONLY {
        return Err(HeaderRejection::Flags { got: prefix[3] });
    }
    if prefix[10..12] != BGZF_XLEN {
        return Err(HeaderRejection::ExtraFieldLength {
            got: u16::from_le_bytes([prefix[10], prefix[11]]),
        });
    }
    if prefix[12..14] != BC_SUBFIELD_ID {
        return Err(HeaderRejection::SubfieldId { got: [prefix[12], prefix[13]] });
    }
    if prefix[14..16] != BC_SUBFIELD_LEN {
        return Err(HeaderRejection::SubfieldLength {
            got: u16::from_le_bytes([prefix[14], prefix[15]]),
        });
    }
    Ok(())
}

/// Whether `prefix` starts with a BGZF block header fgumi can decode.
///
/// The boolean form of [`validate`], for callers that only need the verdict —
/// the input classifier, which reports a format rather than a parse failure.
#[must_use]
pub fn is_bgzf_header(prefix: &[u8]) -> bool {
    validate(prefix).is_ok()
}

/// Total size of the block starting at the front of `block`, from its `BSIZE`.
///
/// `BSIZE` is stored as total size minus one. Returns `None` when `block` is too
/// short to hold a header; the value is *not* validated against
/// [`validate`], because the framing loops that call this per block run over a
/// stream whose header was already checked, and re-checking every block's fields
/// to compute an offset would cost a branch per block for no new information.
#[must_use]
pub fn block_size(block: &[u8]) -> Option<usize> {
    if block.len() < BGZF_HEADER_SIZE {
        return None;
    }
    let bsize = u16::from_le_bytes([block[BSIZE_OFFSET], block[BSIZE_OFFSET + 1]]);
    Some(usize::from(bsize) + 1)
}

/// The smallest `BSIZE` that can describe a real block: a header and a footer,
/// with an empty deflate payload between them.
pub const MIN_BLOCK_SIZE: usize = BGZF_HEADER_SIZE + crate::reader::BGZF_FOOTER_SIZE;

/// [`block_size`], rejecting a `BSIZE` too small to describe a block at all.
///
/// This is the form callers that *use* the size arithmetically want.
/// [`block_size`] reports `BSIZE` as stored, and a malformed or corrupt block can
/// store a value below [`MIN_BLOCK_SIZE`] — at which point `offset + size - 4`
/// underflows and the following index is wildly out of bounds. Every caller has to
/// enforce that floor, and each one enforcing its own copy is how the check came
/// to be missing from one of them: two framing loops in the FASTQ pipeline took
/// the raw value and panicked on input bytes alone.
///
/// Returns `None` when `block` is too short to hold a header *or* when the stored
/// `BSIZE` is below the floor; callers that need to tell those apart can check the
/// slice length first, which their loop guard has usually already done.
#[must_use]
pub fn block_size_checked(block: &[u8]) -> Option<usize> {
    block_size(block).filter(|size| *size >= MIN_BLOCK_SIZE)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    /// A valid 18-byte BGZF header: gzip magic, deflate, FEXTRA, then the
    /// mandatory `BC` subfield carrying `BSIZE`.
    const BGZF: [u8; 18] = [
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00, 0x1b, 0x00,
    ];

    /// Returns `BGZF` with one byte replaced, for building near-miss headers.
    fn bgzf_with(index: usize, value: u8) -> [u8; 18] {
        let mut header = BGZF;
        header[index] = value;
        header
    }

    #[rstest]
    #[case::valid(&BGZF, None)]
    #[case::not_gzip(b"@HD\tVN:1.6\tSO:queryname\n", Some(HeaderRejection::Magic { got: [b'@', b'H'] }))]
    #[case::truncated(&BGZF[..17], Some(HeaderRejection::TooShort { len: 17 }))]
    #[case::not_deflate(&bgzf_with(2, 0x00), Some(HeaderRejection::CompressionMethod { got: 0 }))]
    // No FEXTRA at all: plain gzip.
    #[case::plain_gzip(&bgzf_with(3, 0x00), Some(HeaderRejection::Flags { got: 0 }))]
    // FEXTRA plus FNAME. The extra field precedes FNAME in a gzip member, so
    // `XLEN` and `BC` are still where they belong and htslib would take this;
    // `noodles_bgzf` requires FEXTRA alone, so fgumi cannot promise a decode.
    #[case::fextra_with_fname(&bgzf_with(3, 0x0c), Some(HeaderRejection::Flags { got: 0x0c }))]
    // `BC` first but a second subfield trailing it. RFC 1952 permits this; no
    // decoder in the pipeline reads it.
    #[case::trailing_subfield(&bgzf_with(10, 0x0a), Some(HeaderRejection::ExtraFieldLength { got: 10 }))]
    #[case::foreign_subfield(&bgzf_with(12, b'X'), Some(HeaderRejection::SubfieldId { got: [b'X', b'C'] }))]
    #[case::wrong_subfield_len(&bgzf_with(14, 0x04), Some(HeaderRejection::SubfieldLength { got: 4 }))]
    fn validates(#[case] prefix: &[u8], #[case] expected: Option<HeaderRejection>) {
        assert_eq!(validate(prefix).err(), expected);
        assert_eq!(is_bgzf_header(prefix), expected.is_none());
    }

    /// Every rejection renders a message naming both the field that failed and the
    /// value it saw.
    ///
    /// One case per variant rather than a spot check: the whole point of carrying
    /// the rejection is that the message names the *specific* field, so a variant
    /// whose `Display` arm was never rendered is one whose message could be wrong
    /// -- or, in `SubfieldId`'s case, could render bytes that are not printable
    /// characters -- without any test noticing.
    #[rstest]
    #[case::too_short(HeaderRejection::TooShort { len: 17 }, "truncated", "17")]
    #[case::magic(HeaderRejection::Magic { got: [b'@', b'H'] }, "magic", "0x40")]
    #[case::compression_method(
        HeaderRejection::CompressionMethod { got: 0x00 },
        "compression method",
        "0x00"
    )]
    #[case::flags(HeaderRejection::Flags { got: 0x0c }, "FEXTRA", "0x0c")]
    #[case::extra_field_length(
        HeaderRejection::ExtraFieldLength { got: 10 },
        "extra-field length",
        "10"
    )]
    #[case::subfield_id(HeaderRejection::SubfieldId { got: [b'X', b'C'] }, "subfield ID", "XC")]
    #[case::subfield_length(HeaderRejection::SubfieldLength { got: 4 }, "subfield length", "4")]
    fn rejections_render(
        #[case] rejection: HeaderRejection,
        #[case] names_field: &str,
        #[case] names_value: &str,
    ) {
        let rendered = rejection.to_string();
        assert!(
            rendered.contains(names_field),
            "message must name the failing field {names_field:?}, got: {rendered}"
        );
        assert!(
            rendered.contains(names_value),
            "message must name the observed value {names_value:?}, got: {rendered}"
        );
    }

    #[rstest]
    // BSIZE is total size minus one: 0x1b == 27 means a 28-byte block.
    #[case::eof_block(&BGZF, Some(28))]
    #[case::too_short(&BGZF[..10], None)]
    fn reads_block_size(#[case] block: &[u8], #[case] expected: Option<usize>) {
        assert_eq!(block_size(block), expected);
    }

    /// The EOF marker the crate writes must satisfy the predicate it reads with.
    #[test]
    fn eof_marker_is_a_valid_header() {
        assert_eq!(validate(&crate::BGZF_EOF), Ok(()));
        assert_eq!(block_size(&crate::BGZF_EOF), Some(crate::BGZF_EOF.len()));
    }
}
