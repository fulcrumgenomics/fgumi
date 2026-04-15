//! [`SamTag`] newtype for two-byte SAM/BAM tag identifiers.
//!
//! SAM/BAM aux tags match `[A-Za-z][A-Za-z0-9]` (e.g. `RX`, `MI`, `CB`).
//! [`SamTag`] wraps `[u8; 2]` and enforces that syntax at construction time.
//! Use [`SamTag::new`] (const, panics on invalid bytes) for compile-time
//! constants, or [`std::str::FromStr`] (fallible, used by clap for CLI argument
//! parsing) for runtime input.
//!
//! # Named Constants
//!
//! Standard SAM-spec and fgumi-internal tags are available as associated constants:
//!
//! ```
//! use fgumi_sam::SamTag;
//!
//! assert_eq!(SamTag::RX.to_string(), "RX");
//! assert_eq!(SamTag::MI.to_string(), "MI");
//! assert_eq!(SamTag::CB.to_string(), "CB");
//! ```
//!
//! # CLI Integration
//!
//! [`SamTag`] implements [`std::str::FromStr`], so clap will parse and validate
//! tag arguments automatically:
//!
//! ```
//! use fgumi_sam::SamTag;
//!
//! let tag: SamTag = "XT".parse().expect("valid tag");
//! assert_eq!(*tag, [b'X', b'T']);
//!
//! let err = "ABC".parse::<SamTag>();
//! assert!(err.is_err());
//! ```

use std::fmt;
use std::ops::Deref;
use std::str::FromStr;

use noodles::sam::alignment::record::data::field::Tag;

/// A validated two-byte SAM/BAM tag identifier (e.g. `RX`, `MI`, `CB`).
///
/// Wraps `[u8; 2]` and guarantees the SAM aux-tag pattern `[A-Za-z][A-Za-z0-9]`. Use
/// [`SamTag::new`] for compile-time constants, or parse from a `&str` via
/// [`FromStr`] / [`TryFrom<&str>`] for runtime input.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SamTag([u8; 2]);

impl SamTag {
    /// Returns true if `[a, b]` matches the SAM aux-tag pattern `[A-Za-z][A-Za-z0-9]`.
    #[must_use]
    pub const fn is_valid_tag_bytes(a: u8, b: u8) -> bool {
        a.is_ascii_alphabetic() && (b.is_ascii_alphabetic() || b.is_ascii_digit())
    }

    /// Construct a `SamTag` from two raw bytes.
    ///
    /// This is a `const fn` intended for use in named constant definitions.
    /// Both bytes must match the SAM aux-tag syntax `[A-Za-z][A-Za-z0-9]`;
    /// the assertion is checked at compile time in `const` contexts and panics
    /// at runtime otherwise.
    ///
    /// # Panics
    ///
    /// Panics if the bytes do not match `[A-Za-z][A-Za-z0-9]`.
    #[must_use]
    pub const fn new(a: u8, b: u8) -> Self {
        assert!(Self::is_valid_tag_bytes(a, b), "SAM tag must match [A-Za-z][A-Za-z0-9]");
        Self([a, b])
    }

    /// Convert to a noodles [`Tag`] in `const` context.
    #[must_use]
    pub const fn to_noodles_tag(self) -> Tag {
        Tag::new(self.0[0], self.0[1])
    }

    // ── SAM spec standard tags ──────────────────────────────────────────────

    /// Raw UMI sequence (SAM spec).
    pub const RX: SamTag = SamTag::new(b'R', b'X');
    /// UMI base quality scores (SAM spec).
    pub const QX: SamTag = SamTag::new(b'Q', b'X');
    /// Molecule identifier assigned by grouping (SAM spec).
    pub const MI: SamTag = SamTag::new(b'M', b'I');
    /// Cell barcode sequence (SAM spec).
    pub const CB: SamTag = SamTag::new(b'C', b'B');
    /// Cell barcode quality scores (SAM spec).
    pub const CY: SamTag = SamTag::new(b'C', b'Y');
    /// Sample/library barcode bases (SAM spec).
    pub const BC: SamTag = SamTag::new(b'B', b'C');
    /// Sample/library barcode quality scores (SAM spec).
    pub const QT: SamTag = SamTag::new(b'Q', b'T');
    /// Original UMI bases, before error correction (SAM spec).
    pub const OX: SamTag = SamTag::new(b'O', b'X');
    /// Original UMI base quality scores, before error correction (SAM spec).
    pub const BZ: SamTag = SamTag::new(b'B', b'Z');
    /// Read group (SAM spec).
    pub const RG: SamTag = SamTag::new(b'R', b'G');
    /// Number of mismatches (edit distance) to the reference (SAM spec).
    pub const NM: SamTag = SamTag::new(b'N', b'M');
    /// String encoding of the reference-relative alignment (SAM spec).
    pub const MD: SamTag = SamTag::new(b'M', b'D');
    /// Mate's mapping quality (SAM spec).
    pub const MQ: SamTag = SamTag::new(b'M', b'Q');
    /// Mate's CIGAR string (SAM spec).
    pub const MC: SamTag = SamTag::new(b'M', b'C');
    /// Mate's alignment score, lowercase per SAM spec convention (SAM spec).
    pub const MS: SamTag = SamTag::new(b'm', b's');
    /// Alignment score (SAM spec / aligner-defined).
    pub const AS: SamTag = SamTag::new(b'A', b'S');
    /// Suboptimal alignment score (aligner-defined, widely used).
    pub const XS: SamTag = SamTag::new(b'X', b'S');
    /// Program record identifier (SAM spec).
    pub const PG: SamTag = SamTag::new(b'P', b'G');
    /// Adapter clipping position (e.g. from Picard `MarkIlluminaAdapters`).
    pub const XT: SamTag = SamTag::new(b'X', b'T');

    // ── fgumi-internal tags ─────────────────────────────────────────────────

    /// Template-coordinate sort key for secondary/supplementary reads.
    ///
    /// Stores the primary alignments' template-coordinate sort key so
    /// secondary/supplementary reads can sort adjacent to their primaries.
    /// Lowercase per the SAM specification convention for non-standard
    /// program tags. Written by `fgumi zipper`; validated by `fgumi dedup`.
    pub const TC: SamTag = SamTag::new(b't', b'c');
}

impl Deref for SamTag {
    type Target = [u8; 2];

    fn deref(&self) -> &[u8; 2] {
        &self.0
    }
}

impl AsRef<[u8; 2]> for SamTag {
    fn as_ref(&self) -> &[u8; 2] {
        &self.0
    }
}

impl From<SamTag> for [u8; 2] {
    fn from(tag: SamTag) -> [u8; 2] {
        tag.0
    }
}

impl From<SamTag> for Tag {
    fn from(tag: SamTag) -> Tag {
        Tag::from(tag.0)
    }
}

impl PartialEq<[u8; 2]> for SamTag {
    fn eq(&self, other: &[u8; 2]) -> bool {
        self.0 == *other
    }
}

impl PartialEq<SamTag> for [u8; 2] {
    fn eq(&self, other: &SamTag) -> bool {
        *self == other.0
    }
}

impl fmt::Display for SamTag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.0[0] as char, self.0[1] as char)
    }
}

impl FromStr for SamTag {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> anyhow::Result<Self> {
        let bytes = s.as_bytes();
        anyhow::ensure!(
            bytes.len() == 2 && Self::is_valid_tag_bytes(bytes[0], bytes[1]),
            "SAM tag must match [A-Za-z][A-Za-z0-9], got: {s:?}"
        );
        Ok(SamTag([bytes[0], bytes[1]]))
    }
}

impl TryFrom<&str> for SamTag {
    type Error = anyhow::Error;

    fn try_from(s: &str) -> anyhow::Result<Self> {
        s.parse()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constants_are_correct() {
        assert_eq!(*SamTag::RX, [b'R', b'X']);
        assert_eq!(*SamTag::MI, [b'M', b'I']);
        assert_eq!(*SamTag::CB, [b'C', b'B']);
        assert_eq!(*SamTag::QX, [b'Q', b'X']);
        assert_eq!(*SamTag::CY, [b'C', b'Y']);
        assert_eq!(*SamTag::TC, [b't', b'c']);
        assert_eq!(*SamTag::MS, [b'm', b's']);
        assert_eq!(*SamTag::NM, [b'N', b'M']);
    }

    #[test]
    fn test_display() {
        assert_eq!(SamTag::RX.to_string(), "RX");
        assert_eq!(SamTag::MI.to_string(), "MI");
        assert_eq!(SamTag::CB.to_string(), "CB");
        assert_eq!(SamTag::TC.to_string(), "tc");
    }

    #[test]
    fn test_parse_valid() {
        assert_eq!("RX".parse::<SamTag>().unwrap(), SamTag::RX);
        assert_eq!("XT".parse::<SamTag>().unwrap(), SamTag::XT);
    }

    #[test]
    fn test_parse_invalid_length() {
        assert!("R".parse::<SamTag>().is_err());
        assert!("RXX".parse::<SamTag>().is_err());
        assert!("".parse::<SamTag>().is_err());
    }

    #[test]
    fn test_parse_invalid_non_ascii() {
        assert!("R\x01".parse::<SamTag>().is_err());
        assert!(" X".parse::<SamTag>().is_err());
    }

    #[test]
    fn test_parse_rejects_non_aux_tag_syntax() {
        // Punctuation in second byte
        assert!("A!".parse::<SamTag>().is_err());
        // Digit in first byte
        assert!("1A".parse::<SamTag>().is_err());
        // Punctuation in first byte
        assert!("!A".parse::<SamTag>().is_err());
        // Digit in second byte is allowed
        assert!("A0".parse::<SamTag>().is_ok());
    }

    #[test]
    fn test_into_noodles_tag() {
        let t: Tag = SamTag::CB.into();
        assert_eq!(t, Tag::from([b'C', b'B']));
    }

    #[test]
    fn test_to_noodles_tag_const() {
        const T: Tag = SamTag::TC.to_noodles_tag();
        assert_eq!(T, Tag::from([b't', b'c']));
    }

    #[test]
    fn test_partial_eq_with_bytes() {
        assert_eq!(SamTag::RX, *b"RX");
        assert_eq!(*b"MI", SamTag::MI);
    }

    #[test]
    fn test_deref_to_byte_array() {
        let tag = SamTag::CB;
        let bytes: &[u8; 2] = &tag;
        assert_eq!(bytes, b"CB");
    }
}
