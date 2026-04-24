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
//! use fgumi_tag::SamTag;
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
//! use fgumi_tag::SamTag;
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
    /// Bisulfite strand (written by `bwameth` and related tools): `Z:f` for the
    /// forward/top strand, `Z:r` for the reverse/bottom strand. Used by `fgumi
    /// zipper --restore-unconverted-bases` to decide which bases to restore.
    pub const YD: SamTag = SamTag::new(b'Y', b'D');
    /// Base modification calls, e.g. `MM:Z:C+m,5,12,0;` (SAM spec, "Base modifications").
    pub const MM: SamTag = SamTag::new(b'M', b'M');
    /// Per-base modification probabilities, paired with `MM` (SAM spec).
    pub const ML: SamTag = SamTag::new(b'M', b'L');
    /// Phred-likelihood that the segment is incorrect (SAM spec).
    pub const UQ: SamTag = SamTag::new(b'U', b'Q');
    /// Original base quality scores, before recalibration (SAM spec).
    pub const OQ: SamTag = SamTag::new(b'O', b'Q');
    /// Per-base offset to base alignment quality (BAQ), `B:C` array (SAM spec).
    ///
    /// Note: the Rust identifier uses a descriptive suffix to avoid a name
    /// collision with [`SamTag::BD`], which is the fgumi-internal per-read BA
    /// single-strand depth tag (`bD`).
    pub const BAQ_DELTA: SamTag = SamTag::new(b'B', b'D');

    // ── fgumi-internal tags ─────────────────────────────────────────────────

    /// Template-coordinate sort key for secondary/supplementary reads.
    ///
    /// Stores the primary alignments' template-coordinate sort key so
    /// secondary/supplementary reads can sort adjacent to their primaries.
    /// Lowercase per the SAM specification convention for non-standard
    /// program tags. Written by `fgumi zipper`; validated by `fgumi dedup`.
    pub const TC: SamTag = SamTag::new(b't', b'c');

    // ── fgbio-style consensus tags (per-read) ──────────────────────────────
    //
    // First letter encodes the strand: `c` = combined consensus, `a` = AB
    // single-strand, `b` = BA single-strand. Second letter (uppercase = per-read,
    // lowercase = per-base array) encodes the metric.

    /// Per-read max raw-read depth contributing to the combined consensus (`cD`).
    pub const CD: SamTag = SamTag::new(b'c', b'D');
    /// Per-read min raw-read depth contributing to the combined consensus (`cM`).
    pub const CM: SamTag = SamTag::new(b'c', b'M');
    /// Per-read raw-read disagreement rate vs. the combined consensus (`cE`).
    pub const CE: SamTag = SamTag::new(b'c', b'E');
    /// Per-read max raw-read depth for the AB single-strand consensus (`aD`).
    pub const AD: SamTag = SamTag::new(b'a', b'D');
    /// Per-read min raw-read depth for the AB single-strand consensus (`aM`).
    pub const AM: SamTag = SamTag::new(b'a', b'M');
    /// Per-read raw-read disagreement rate vs. the AB single-strand consensus (`aE`).
    pub const AE: SamTag = SamTag::new(b'a', b'E');
    /// Per-read max raw-read depth for the BA single-strand consensus (`bD`).
    pub const BD: SamTag = SamTag::new(b'b', b'D');
    /// Per-read min raw-read depth for the BA single-strand consensus (`bM`).
    pub const BM: SamTag = SamTag::new(b'b', b'M');
    /// Per-read raw-read disagreement rate vs. the BA single-strand consensus (`bE`).
    pub const BE: SamTag = SamTag::new(b'b', b'E');

    // ── fgbio-style consensus tags (per-base arrays) ───────────────────────
    //
    // Lowercase second byte distinguishes per-base i16 arrays from the per-read
    // scalars above. The `_BASES` suffix on the Rust identifier disambiguates
    // from the matching per-read constant (case-insensitive identifier collision).

    /// Per-base raw-read depth array for the combined consensus (`cd`, `B,s`).
    pub const CD_BASES: SamTag = SamTag::new(b'c', b'd');
    /// Per-base raw-read disagreement count array for the combined consensus (`ce`, `B,s`).
    pub const CE_BASES: SamTag = SamTag::new(b'c', b'e');
    /// Per-base raw-read depth array for the AB single-strand consensus (`ad`, `B,s`).
    pub const AD_BASES: SamTag = SamTag::new(b'a', b'd');
    /// Per-base raw-read disagreement count array for the AB single-strand consensus (`ae`, `B,s`).
    pub const AE_BASES: SamTag = SamTag::new(b'a', b'e');
    /// Per-base raw-read depth array for the BA single-strand consensus (`bd`, `B,s`).
    pub const BD_BASES: SamTag = SamTag::new(b'b', b'd');
    /// Per-base raw-read disagreement count array for the BA single-strand consensus (`be`, `B,s`).
    pub const BE_BASES: SamTag = SamTag::new(b'b', b'e');

    /// AB single-strand consensus base sequence (`ac`, `Z`).
    pub const AC: SamTag = SamTag::new(b'a', b'c');
    /// BA single-strand consensus base sequence (`bc`, `Z`).
    pub const BC_BASES: SamTag = SamTag::new(b'b', b'c');
    /// AB single-strand consensus Phred+33 quality string (`aq`, `Z`).
    pub const AQ: SamTag = SamTag::new(b'a', b'q');
    /// BA single-strand consensus Phred+33 quality string (`bq`, `Z`).
    pub const BQ: SamTag = SamTag::new(b'b', b'q');

    // ── fgbio-style EM-Seq / methylation consensus tags ────────────────────

    /// Per-base unconverted-cytosine count for the combined consensus (`cu`, `B,s`).
    pub const CU: SamTag = SamTag::new(b'c', b'u');
    /// Per-base converted-cytosine count for the combined consensus (`ct`, `B,s`).
    pub const CT: SamTag = SamTag::new(b'c', b't');
    /// Per-base unconverted-cytosine count for the AB single-strand consensus (`au`, `B,s`).
    pub const AU: SamTag = SamTag::new(b'a', b'u');
    /// Per-base converted-cytosine count for the AB single-strand consensus (`at`, `B,s`).
    pub const AT: SamTag = SamTag::new(b'a', b't');
    /// Per-base unconverted-cytosine count for the BA single-strand consensus (`bu`, `B,s`).
    pub const BU: SamTag = SamTag::new(b'b', b'u');
    /// Per-base converted-cytosine count for the BA single-strand consensus (`bt`, `B,s`).
    pub const BT: SamTag = SamTag::new(b'b', b't');

    /// Per-base methylation-modification string for the AB single-strand consensus (`am`, `Z`).
    ///
    /// Stores an MM-format annotation (e.g. `C+m,5,12,0;`) for the top strand.
    /// Uses the `_BASES` suffix to disambiguate from per-read [`SamTag::AM`]
    /// (bytes `aM`, min depth scalar), which shares the Rust identifier stem.
    pub const AM_BASES: SamTag = SamTag::new(b'a', b'm');

    /// Per-base methylation-modification string for the BA single-strand consensus (`bm`, `Z`).
    ///
    /// Stores an MM-format annotation (e.g. `G-m,...`) for the bottom strand.
    /// Uses the `_BASES` suffix to disambiguate from per-read [`SamTag::BM`]
    /// (bytes `bM`, min depth scalar), which shares the Rust identifier stem.
    pub const BM_BASES: SamTag = SamTag::new(b'b', b'm');

    // ── Slice helpers ──────────────────────────────────────────────────────

    /// All per-read consensus tags (combined, AB strand, BA strand).
    pub const PER_READ_CONSENSUS_TAGS: &'static [SamTag] =
        &[Self::CD, Self::CM, Self::CE, Self::AD, Self::AM, Self::AE, Self::BD, Self::BM, Self::BE];

    /// All per-base consensus tags, including methylation conversion counts.
    pub const PER_BASE_CONSENSUS_TAGS: &'static [SamTag] = &[
        Self::CD_BASES,
        Self::CE_BASES,
        Self::AD_BASES,
        Self::AE_BASES,
        Self::BD_BASES,
        Self::BE_BASES,
        Self::AC,
        Self::BC_BASES,
        Self::AQ,
        Self::BQ,
        Self::CU,
        Self::CT,
        Self::AU,
        Self::AT,
        Self::BU,
        Self::BT,
    ];

    /// Per-base tags that must be **reversed** (not reverse-complemented) when
    /// the read is on the negative strand. These are scalar/quality arrays
    /// keyed by read position.
    pub const PER_BASE_TAGS_TO_REVERSE: &'static [SamTag] = &[
        Self::CD_BASES,
        Self::CE_BASES,
        Self::AD_BASES,
        Self::AE_BASES,
        Self::BD_BASES,
        Self::BE_BASES,
        Self::AQ,
        Self::BQ,
        Self::CU,
        Self::CT,
        Self::AU,
        Self::AT,
        Self::BU,
        Self::BT,
    ];

    /// Per-base tags that must be **reverse-complemented** when the read is on
    /// the negative strand. These are base-sequence strings.
    pub const PER_BASE_TAGS_TO_REVCOMP: &'static [SamTag] = &[Self::AC, Self::BC_BASES];
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

    #[test]
    fn test_sam_spec_extension_constants() {
        assert_eq!(*SamTag::MM, [b'M', b'M']);
        assert_eq!(*SamTag::ML, [b'M', b'L']);
        assert_eq!(*SamTag::UQ, [b'U', b'Q']);
        assert_eq!(*SamTag::OQ, [b'O', b'Q']);
    }

    #[test]
    fn test_per_read_consensus_constants() {
        assert_eq!(*SamTag::CD, [b'c', b'D']);
        assert_eq!(*SamTag::CM, [b'c', b'M']);
        assert_eq!(*SamTag::CE, [b'c', b'E']);
        assert_eq!(*SamTag::AD, [b'a', b'D']);
        assert_eq!(*SamTag::AM, [b'a', b'M']);
        assert_eq!(*SamTag::AE, [b'a', b'E']);
        assert_eq!(*SamTag::BD, [b'b', b'D']);
        assert_eq!(*SamTag::BM, [b'b', b'M']);
        assert_eq!(*SamTag::BE, [b'b', b'E']);
    }

    #[test]
    fn test_per_base_consensus_constants() {
        // depth/error count arrays
        assert_eq!(*SamTag::CD_BASES, [b'c', b'd']);
        assert_eq!(*SamTag::CE_BASES, [b'c', b'e']);
        assert_eq!(*SamTag::AD_BASES, [b'a', b'd']);
        assert_eq!(*SamTag::AE_BASES, [b'a', b'e']);
        assert_eq!(*SamTag::BD_BASES, [b'b', b'd']);
        assert_eq!(*SamTag::BE_BASES, [b'b', b'e']);
        // single-strand consensus bases / qualities
        assert_eq!(*SamTag::AC, [b'a', b'c']);
        assert_eq!(*SamTag::BC_BASES, [b'b', b'c']);
        assert_eq!(*SamTag::AQ, [b'a', b'q']);
        assert_eq!(*SamTag::BQ, [b'b', b'q']);
        // EM-Seq methylation conversion counts
        assert_eq!(*SamTag::CU, [b'c', b'u']);
        assert_eq!(*SamTag::CT, [b'c', b't']);
        assert_eq!(*SamTag::AU, [b'a', b'u']);
        assert_eq!(*SamTag::AT, [b'a', b't']);
        assert_eq!(*SamTag::BU, [b'b', b'u']);
        assert_eq!(*SamTag::BT, [b'b', b't']);
    }

    #[test]
    fn test_per_base_constants_are_distinct_from_per_read() {
        // Lowercase second byte = per-base array, uppercase = per-read scalar.
        assert_ne!(SamTag::CD, SamTag::CD_BASES);
        assert_ne!(SamTag::CE, SamTag::CE_BASES);
        assert_ne!(SamTag::AD, SamTag::AD_BASES);
        assert_ne!(SamTag::AE, SamTag::AE_BASES);
        assert_ne!(SamTag::BD, SamTag::BD_BASES);
        assert_ne!(SamTag::BE, SamTag::BE_BASES);
        // BC_BASES collides with the SAM-spec sample-barcode BC tag.
        assert_ne!(SamTag::BC, SamTag::BC_BASES);
    }

    #[test]
    fn test_per_base_methylation_annotation_constants() {
        assert_eq!(*SamTag::AM_BASES, [b'a', b'm']);
        assert_eq!(*SamTag::BM_BASES, [b'b', b'm']);
        // These must NOT collide with the per-read AM/BM (which have uppercase M).
        assert_ne!(SamTag::AM, SamTag::AM_BASES);
        assert_ne!(SamTag::BM, SamTag::BM_BASES);
    }

    #[test]
    fn test_per_base_consensus_tag_set() {
        let tags = SamTag::PER_BASE_CONSENSUS_TAGS;
        // exact set & length — guards against accidental drift
        let expected = [
            SamTag::CD_BASES,
            SamTag::CE_BASES,
            SamTag::AD_BASES,
            SamTag::AE_BASES,
            SamTag::BD_BASES,
            SamTag::BE_BASES,
            SamTag::AC,
            SamTag::BC_BASES,
            SamTag::AQ,
            SamTag::BQ,
            SamTag::CU,
            SamTag::CT,
            SamTag::AU,
            SamTag::AT,
            SamTag::BU,
            SamTag::BT,
        ];
        assert_eq!(tags.len(), expected.len());
        for t in &expected {
            assert!(tags.contains(t), "missing {t} from PER_BASE_CONSENSUS_TAGS");
        }
    }

    #[test]
    fn test_per_read_consensus_tag_set() {
        let tags = SamTag::PER_READ_CONSENSUS_TAGS;
        let expected = [
            SamTag::CD,
            SamTag::CM,
            SamTag::CE,
            SamTag::AD,
            SamTag::AM,
            SamTag::AE,
            SamTag::BD,
            SamTag::BM,
            SamTag::BE,
        ];
        assert_eq!(tags.len(), expected.len());
        for t in &expected {
            assert!(tags.contains(t), "missing {t} from PER_READ_CONSENSUS_TAGS");
        }
    }

    #[test]
    fn test_per_base_tags_to_reverse_set() {
        let tags = SamTag::PER_BASE_TAGS_TO_REVERSE;
        for t in tags {
            assert!(SamTag::PER_BASE_CONSENSUS_TAGS.contains(t));
        }
        // Bases tags are reverse-COMPLEMENTED, not just reversed; they must NOT appear here.
        assert!(!tags.contains(&SamTag::AC));
        assert!(!tags.contains(&SamTag::BC_BASES));
    }

    #[test]
    fn test_per_base_tags_to_revcomp_set() {
        assert_eq!(SamTag::PER_BASE_TAGS_TO_REVCOMP, &[SamTag::AC, SamTag::BC_BASES]);
    }

    #[test]
    fn test_per_base_partitions_are_exhaustive() {
        // REVERSE and REVCOMP must together cover every tag in PER_BASE_CONSENSUS_TAGS.
        assert_eq!(
            SamTag::PER_BASE_TAGS_TO_REVERSE.len() + SamTag::PER_BASE_TAGS_TO_REVCOMP.len(),
            SamTag::PER_BASE_CONSENSUS_TAGS.len(),
            "REVERSE + REVCOMP must partition PER_BASE_CONSENSUS_TAGS exactly",
        );
        for t in SamTag::PER_BASE_TAGS_TO_REVCOMP {
            assert!(
                SamTag::PER_BASE_CONSENSUS_TAGS.contains(t),
                "REVCOMP tag {t} missing from PER_BASE_CONSENSUS_TAGS",
            );
        }
        for t in SamTag::PER_BASE_TAGS_TO_REVERSE {
            assert!(
                SamTag::PER_BASE_CONSENSUS_TAGS.contains(t),
                "REVERSE tag {t} missing from PER_BASE_CONSENSUS_TAGS",
            );
        }
    }
}
