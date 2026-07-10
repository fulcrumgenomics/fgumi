//! A local read-structure implementation that supports a `+` (any-length) segment in **any**
//! position, matching fgbio 4.1.0's `ReadStructure` (PR #1157).
//!
//! # Why not the `read-structure` crate?
//! fgumi depends on `read-structure = "0.2.0"`, which rejects a non-terminal `+` at parse time
//! and, for a fully-fixed structure, silently ignores trailing bases on an over-long read. fgbio
//! 4.1.0 does neither: it accepts a `+` anywhere (resolving segments after the `+` by walking
//! back from the read end) and *requires* a fully-fixed structure to match the read length
//! exactly. Bumping the crate to 0.3.0 would fix the first gap but change the trailing-`+`
//! semantics from zero-or-more to one-or-more (a new parity break) and force a breaking API
//! migration; the crate pin is therefore held at 0.2.0 and this module reimplements the parsing,
//! offset, and length-validation logic faithfully to fgbio. The crate's [`SegmentType`] enum is
//! reused as-is.
//!
//! The model mirrors `fgbio` `com.fulcrumgenomics.util.ReadStructure`:
//! - At most one segment may be indefinite-length (`+`); it may sit at any index.
//! - Segments before/at the `+` (or all segments, when there is no `+`) have forward offsets
//!   measured from the read start. Segments strictly after the `+` have offsets measured from the
//!   read end (stored as positive distances from the read end).
//! - The `+` segment absorbs `read_len - fixed_length_sum` bases, i.e. it is **zero-or-more**.

use std::fmt;
use std::str::FromStr;

pub use read_structure::SegmentType;

/// The `+` character marking an any-length segment.
const ANY_LENGTH_CHAR: char = '+';

/// A single segment of a [`ReadStructure`]: a kind and an optional fixed length (`None` == `+`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadSegment {
    /// The kind of the segment (template, molecular/sample/cellular barcode, or skip).
    pub kind: SegmentType,
    /// The fixed length of the segment, or `None` for the any-length (`+`) segment.
    length: Option<usize>,
}

impl ReadSegment {
    /// Returns the fixed length of the segment, or `None` for the any-length (`+`) segment.
    #[must_use]
    pub fn length(&self) -> Option<usize> {
        self.length
    }

    /// Returns true if the segment has a fixed length (i.e. is not the `+` segment).
    #[must_use]
    pub fn has_length(&self) -> bool {
        self.length.is_some()
    }
}

impl fmt::Display for ReadSegment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.length {
            Some(l) => write!(f, "{l}")?,
            None => write!(f, "{ANY_LENGTH_CHAR}")?,
        }
        write!(f, "{}", self.kind.value())
    }
}

/// The outcome of validating a read length against a [`ReadStructure`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LengthCheck {
    /// The read length is acceptable for this structure.
    Ok,
    /// The read is shorter than the fixed bases the structure requires (`need` bases minimum).
    TooFew {
        /// The minimum number of bases the structure requires.
        need: usize,
    },
    /// The structure is fully fixed and the read is longer than its fixed length (`need` bases).
    ///
    /// fgbio treats silent truncation of trailing bases as almost always a bug, so an over-long
    /// read against a fully-fixed structure is an error rather than a truncation.
    OverLong {
        /// The exact number of bases a fully-fixed structure requires.
        need: usize,
    },
}

/// Errors that can occur while parsing a [`ReadStructure`] from a string.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadStructureError {
    /// The read structure had no segments.
    Empty,
    /// The read structure had more than one `+` (any-length) segment.
    MultiplePlus(String),
    /// A segment was given an explicit length of zero.
    ZeroLength(String),
    /// A length token was expected but not found.
    MissingLength(String),
    /// A segment operator (kind) was expected but not found.
    MissingOperator(String),
    /// A segment operator (kind) was not one of the recognized types.
    UnknownType(String),
    /// A segment length or the accumulated fixed length overflowed `usize`.
    LengthOverflow(String),
}

impl fmt::Display for ReadStructureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "Read structure contained no segments"),
            Self::MultiplePlus(rs) => {
                write!(f, "Read structure contains more than one any-length (+) segment: {rs}")
            }
            Self::ZeroLength(rs) => {
                write!(f, "Read structure contains a zero-length segment: {rs}")
            }
            Self::MissingLength(rs) => {
                write!(f, "Read structure is missing a length before an operator: {rs}")
            }
            Self::MissingOperator(rs) => {
                write!(f, "Read structure is missing a segment operator: {rs}")
            }
            Self::UnknownType(rs) => {
                write!(f, "Read structure contains an unknown segment type: {rs}")
            }
            Self::LengthOverflow(rs) => {
                write!(f, "Read structure segment length overflowed: {rs}")
            }
        }
    }
}

impl std::error::Error for ReadStructureError {}

/// A segment's start position within a read, anchored either to the read start or the read end.
///
/// Segments before/at the `+` (or all segments, when there is no `+`) are anchored to the read
/// start; segments strictly after the `+` are anchored to the read end (their start depends on the
/// read length). Modeling the anchor explicitly avoids signed-integer offset arithmetic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SegOffset {
    /// Start is this many bases from the read start.
    FromStart(usize),
    /// Start is this many bases before the read end (start = `read_len - dist`).
    FromEnd(usize),
}

/// A read structure: an ordered list of [`ReadSegment`]s with at most one any-length (`+`) segment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadStructure {
    /// The segments, in order.
    segments: Vec<ReadSegment>,
    /// The index of the any-length (`+`) segment, if any.
    plus_index: Option<usize>,
    /// Per-segment start anchor (see [`SegOffset`]).
    offsets: Vec<SegOffset>,
    /// The sum of the fixed lengths of the segments strictly after the `+` segment.
    post_plus_len: usize,
    /// The sum of all fixed segment lengths (the `+` contributes zero).
    fixed_length_sum: usize,
}

impl ReadStructure {
    /// Builds a [`ReadStructure`] from its segments, computing offsets and validating the `+` rule.
    fn new(segments: Vec<ReadSegment>, rendered: &str) -> Result<Self, ReadStructureError> {
        if segments.is_empty() {
            return Err(ReadStructureError::Empty);
        }

        let plus_indices: Vec<usize> =
            segments.iter().enumerate().filter(|(_, s)| !s.has_length()).map(|(i, _)| i).collect();
        if plus_indices.len() > 1 {
            return Err(ReadStructureError::MultiplePlus(rendered.to_string()));
        }
        let plus_index = plus_indices.first().copied();

        // Accumulate the total fixed length with checked arithmetic so a pathological structure
        // (many large segments, or a single over-large length token) is rejected rather than
        // wrapping in release builds and later underflowing in `span_of`.
        let overflow = || ReadStructureError::LengthOverflow(rendered.to_string());
        let fixed_length_sum: usize = segments
            .iter()
            .filter_map(ReadSegment::length)
            .try_fold(0usize, usize::checked_add)
            .ok_or_else(overflow)?;

        // Bases occupied by fixed segments strictly after the `+`.
        let post_plus_len: usize = match plus_index {
            Some(p) => segments[p + 1..]
                .iter()
                .filter_map(ReadSegment::length)
                .try_fold(0usize, usize::checked_add)
                .ok_or_else(overflow)?,
            None => 0,
        };

        // Offsets. Forward pass (anchored to the read start) up to and including the `+` (or all
        // segments if no `+`); backward pass (anchored to the read end) for segments after the `+`.
        // Both `off` and `dist_from_end` are partial sums of the fixed lengths, so they are bounded
        // by `fixed_length_sum` (checked above) and cannot overflow.
        let n = segments.len();
        let mut offsets = vec![SegOffset::FromStart(0); n];
        let forward_end = plus_index.map_or(n, |p| p + 1);
        let mut off: usize = 0;
        for (i, seg) in segments.iter().enumerate().take(forward_end) {
            offsets[i] = SegOffset::FromStart(off);
            off += seg.length.unwrap_or(0);
        }
        if let Some(p) = plus_index {
            let mut dist_from_end: usize = 0;
            for i in (p + 1..n).rev() {
                // Post-`+` segments always have a fixed length: only one segment may be indefinite.
                dist_from_end += segments[i].length.unwrap_or(0);
                offsets[i] = SegOffset::FromEnd(dist_from_end);
            }
        }

        Ok(Self { segments, plus_index, offsets, post_plus_len, fixed_length_sum })
    }

    /// Returns the segments of the read structure.
    #[must_use]
    pub fn segments(&self) -> &[ReadSegment] {
        &self.segments
    }

    /// Returns an iterator over the segments of the read structure.
    pub fn iter(&self) -> impl Iterator<Item = &ReadSegment> {
        self.segments.iter()
    }

    /// Returns the number of segments in the read structure.
    #[must_use]
    pub fn number_of_segments(&self) -> usize {
        self.segments.len()
    }

    /// Returns the segments of the given kind, in order.
    pub fn segments_by_type(&self, kind: SegmentType) -> impl Iterator<Item = &ReadSegment> {
        self.segments.iter().filter(move |s| s.kind == kind)
    }

    /// Returns true if the structure has a fixed (non-variable) length (no `+` segment).
    #[must_use]
    pub fn has_fixed_length(&self) -> bool {
        self.plus_index.is_none()
    }

    /// Returns the `[start, end)` span within a read of length `read_len` for the segment at `index`.
    ///
    /// Mirrors fgbio `ReadStructure.spanOf`: the `+` segment runs from its forward start offset to
    /// just before the fixed post-`+` region; a post-`+` segment's start is measured from the read
    /// end.
    ///
    /// # Panics
    /// Panics if `index` is out of range, or if `read_len` is too short to contain the segment
    /// (callers must validate the length via [`ReadStructure::check_read_length`] first).
    #[must_use]
    pub fn span_of(&self, index: usize, read_len: usize) -> (usize, usize) {
        let start = match self.offsets[index] {
            SegOffset::FromStart(s) => s,
            SegOffset::FromEnd(dist) => read_len - dist,
        };
        if self.plus_index == Some(index) {
            // The `+` segment absorbs everything up to the fixed post-`+` region (zero-or-more).
            (start, read_len - self.post_plus_len)
        } else {
            (start, start + self.segments[index].length.expect("fixed segment has a length"))
        }
    }

    /// Validates a read of `read_len` bases against this structure (fgbio `validateReadLength`).
    ///
    /// A structure with a `+` accepts any read with at least the fixed-segment length in bases
    /// (the `+` is zero-or-more). A fully-fixed structure requires the read length to match its
    /// fixed length exactly; an over-long read is [`LengthCheck::OverLong`] rather than being
    /// silently truncated.
    #[must_use]
    pub fn check_read_length(&self, read_len: usize) -> LengthCheck {
        if read_len < self.fixed_length_sum {
            LengthCheck::TooFew { need: self.fixed_length_sum }
        } else if self.has_fixed_length() && read_len > self.fixed_length_sum {
            LengthCheck::OverLong { need: self.fixed_length_sum }
        } else {
            LengthCheck::Ok
        }
    }
}

impl fmt::Display for ReadStructure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for seg in &self.segments {
            write!(f, "{seg}")?;
        }
        Ok(())
    }
}

impl FromStr for ReadStructure {
    type Err = ReadStructureError;

    fn from_str(rs: &str) -> Result<Self, Self::Err> {
        let chars: Vec<char> = rs.to_uppercase().chars().filter(|c| !c.is_whitespace()).collect();
        let rendered: String = chars.iter().collect();
        let mut segments: Vec<ReadSegment> = Vec::new();
        let mut i = 0;
        while i < chars.len() {
            // Length token: '+' (any-length) or one-or-more digits.
            let length = if chars[i] == ANY_LENGTH_CHAR {
                i += 1;
                None
            } else if chars[i].is_ascii_digit() {
                let mut len: usize = 0;
                while i < chars.len() && chars[i].is_ascii_digit() {
                    // Reject an over-large length token instead of wrapping in release builds.
                    len = len
                        .checked_mul(10)
                        .and_then(|v| v.checked_add(chars[i] as usize - '0' as usize))
                        .ok_or_else(|| ReadStructureError::LengthOverflow(rendered.clone()))?;
                    i += 1;
                }
                Some(len)
            } else {
                return Err(ReadStructureError::MissingLength(rendered));
            };

            // Operator token: the segment kind.
            if i >= chars.len() {
                return Err(ReadStructureError::MissingOperator(rendered));
            }
            let kind = SegmentType::try_from(chars[i])
                .map_err(|_| ReadStructureError::UnknownType(rendered.clone()))?;
            if length == Some(0) {
                return Err(ReadStructureError::ZeroLength(rendered));
            }
            i += 1;
            segments.push(ReadSegment { kind, length });
        }

        Self::new(segments, &rendered)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    fn rs(s: &str) -> ReadStructure {
        ReadStructure::from_str(s).expect("valid read structure")
    }

    #[test]
    fn parses_terminal_plus_and_fully_fixed() {
        assert_eq!(rs("5M+T").to_string(), "5M+T");
        assert_eq!(rs("+T").to_string(), "+T");
        assert_eq!(rs("10T").to_string(), "10T");
        assert_eq!(rs("8B8B75T").to_string(), "8B8B75T");
    }

    #[test]
    fn parses_non_terminal_plus() {
        // These are exactly the structures the 0.2.0 crate rejects.
        assert_eq!(rs("8B+M10T").to_string(), "8B+M10T");
        assert_eq!(rs("+M70T").to_string(), "+M70T");
        assert_eq!(rs("2M1S2M+T").to_string(), "2M1S2M+T");
    }

    #[rstest]
    #[case::plus_then_plus_operator("++M")]
    #[case::plus_operator_after_segment("5M++T")]
    #[case::two_any_length_segments("+M+T")]
    #[case::explicit_zero_length("0T")]
    #[case::unknown_segment_type("9R")]
    #[case::operator_without_length("T")]
    #[case::operator_without_following_length("23TT")]
    #[case::trailing_length_without_operator("23T2")]
    fn rejects_malformed_structures(#[case] bad: &str) {
        assert!(ReadStructure::from_str(bad).is_err(), "expected error for {bad}");
    }

    #[test]
    fn bare_plus_without_operator_is_missing_operator() {
        // A bare `+` with no operator is missing its operator.
        assert!(matches!(
            ReadStructure::from_str("8B+"),
            Err(ReadStructureError::MissingOperator(_))
        ));
    }

    #[test]
    fn rejects_length_token_that_overflows_usize() {
        // 20 nines exceed usize::MAX on 64-bit targets, so parsing must reject rather than wrap.
        let overflowing = "9".repeat(20) + "T";
        assert!(matches!(
            ReadStructure::from_str(&overflowing),
            Err(ReadStructureError::LengthOverflow(_))
        ));
    }

    #[test]
    fn multiple_plus_is_a_specific_error() {
        assert!(matches!(
            ReadStructure::from_str("+M+T"),
            Err(ReadStructureError::MultiplePlus(_))
        ));
    }

    #[test]
    fn extracts_non_terminal_plus_by_walking_back_from_end() {
        // 8B + variable M + 10T on a 30bp read: M = bases[8..20], T = bases[20..30].
        let structure = rs("8B+M10T");
        let read_len = 30;
        assert_eq!(structure.span_of(0, read_len), (0, 8)); // 8B
        assert_eq!(structure.span_of(1, read_len), (8, 20)); // +M (absorbs 30 - 8 - 10 = 12)
        assert_eq!(structure.span_of(2, read_len), (20, 30)); // 10T
    }

    #[test]
    fn terminal_plus_absorbs_the_remainder() {
        let structure = rs("4M+T");
        assert_eq!(structure.span_of(0, 10), (0, 4));
        assert_eq!(structure.span_of(1, 10), (4, 10));
        // Zero-or-more: a read that only satisfies the fixed prefix yields an empty `+` segment.
        assert_eq!(structure.span_of(1, 4), (4, 4));
    }

    #[test]
    fn length_check_fully_fixed_requires_exact_length() {
        let structure = rs("8M2T"); // fixed length 10
        assert_eq!(structure.check_read_length(10), LengthCheck::Ok);
        assert_eq!(structure.check_read_length(12), LengthCheck::OverLong { need: 10 });
        assert_eq!(structure.check_read_length(8), LengthCheck::TooFew { need: 10 });
    }

    #[test]
    fn length_check_variable_is_zero_or_more() {
        let structure = rs("4M+T"); // min length 4
        assert_eq!(structure.check_read_length(4), LengthCheck::Ok); // empty +T is allowed
        assert_eq!(structure.check_read_length(100), LengthCheck::Ok);
        assert_eq!(structure.check_read_length(3), LengthCheck::TooFew { need: 4 });
    }

    #[test]
    fn segments_by_type_and_counts() {
        let structure = rs("8B+M10T");
        assert_eq!(structure.number_of_segments(), 3);
        assert_eq!(structure.segments_by_type(SegmentType::Template).count(), 1);
        assert_eq!(structure.segments_by_type(SegmentType::MolecularBarcode).count(), 1);
        assert_eq!(structure.segments_by_type(SegmentType::SampleBarcode).count(), 1);
    }

    proptest::proptest! {
        /// Any well-formed read structure round-trips through `from_str`/`Display`, and for a
        /// read length the structure accepts, the segment spans tile `0..read_len` contiguously
        /// with `start <= end <= read_len`. This exercises the forward/backward offset arithmetic
        /// (including a non-terminal `+`) far beyond the enumerated cases above.
        #[test]
        fn parse_round_trips_and_spans_tile_the_read(
            // 1..=6 segments, each a `(kind index, fixed length >= 1)`.
            segs in proptest::collection::vec((0usize..5, 1usize..=20), 1..=6),
            // Optionally promote one segment to the any-length `+` segment.
            plus_pick in proptest::option::of(0usize..6),
            // Extra bases the `+` absorbs when one is present.
            extra in 0usize..=50,
        ) {
            const KINDS: [char; 5] = ['T', 'B', 'M', 'C', 'S'];
            let plus_index = plus_pick.filter(|&i| i < segs.len());

            // Render the canonical structure string (uppercase digits + kind chars).
            let mut rendered = String::new();
            for (i, (kind, len)) in segs.iter().enumerate() {
                if Some(i) == plus_index {
                    rendered.push(ANY_LENGTH_CHAR);
                } else {
                    rendered.push_str(&len.to_string());
                }
                rendered.push(KINDS[*kind]);
            }

            let structure =
                ReadStructure::from_str(&rendered).expect("generated read structure must parse");

            // Round-trip: `Display` is the canonical rendering.
            proptest::prop_assert_eq!(structure.to_string(), rendered.clone());

            // A read length the structure accepts: a fixed structure needs its exact fixed sum;
            // a `+` structure accepts the fixed sum plus any number of absorbed bases.
            let fixed_sum: usize = segs
                .iter()
                .enumerate()
                .filter(|(i, _)| Some(*i) != plus_index)
                .map(|(_, (_, len))| *len)
                .sum();
            let read_len = if plus_index.is_some() { fixed_sum + extra } else { fixed_sum };
            proptest::prop_assert_eq!(structure.check_read_length(read_len), LengthCheck::Ok);

            // Spans tile `0..read_len` contiguously with no gaps or overlaps.
            let mut cursor = 0usize;
            for i in 0..structure.number_of_segments() {
                let (start, end) = structure.span_of(i, read_len);
                proptest::prop_assert_eq!(start, cursor, "segment {} start is not contiguous", i);
                proptest::prop_assert!(start <= end, "segment {} has start > end", i);
                proptest::prop_assert!(end <= read_len, "segment {} end exceeds read_len", i);
                cursor = end;
            }
            proptest::prop_assert_eq!(cursor, read_len, "segments must cover the whole read");
        }
    }
}
