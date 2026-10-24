//! The `rw` (realign window) tag value.
//!
//! The `rw` tag records the genomic window spanned by the raw reads that were
//! combined into a consensus read, so that a later re-alignment of the consensus
//! can prefer the alignment that falls inside that window. The value has the form
//! `<contig>:<start>-<end>:<strand>` with 1-based inclusive coordinates and a
//! strand of `+` or `-`.

use std::io::Write;

use anyhow::{Result, anyhow, bail};

use crate::xa::parse_digits;

/// A genomic window (1-based, inclusive) on a given strand.
///
/// The contig name is borrowed from the parsed tag value, so parsing allocates
/// nothing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RealignWindow<'a> {
    /// Reference contig name (may itself contain `:`).
    pub contig: &'a str,
    /// 1-based start position.
    pub start: u32,
    /// 1-based inclusive end position.
    pub end: u32,
    /// True if the window is on the reverse strand.
    pub reverse: bool,
}

impl<'a> RealignWindow<'a> {
    /// Parses an `rw` tag value of the form `<contig>:<start>-<end>:<strand>`.
    ///
    /// Fields are split from the right, so the contig name may contain `:`.
    ///
    /// # Errors
    /// Returns an error if the value is not valid UTF-8, is missing a field,
    /// has a coordinate that is not a positive integer (they are 1-based), has
    /// `start > end`, or has a strand other than `+` or `-`.
    pub fn parse(value: &'a [u8]) -> Result<Self> {
        let text = std::str::from_utf8(value).map_err(|e| anyhow!("rw value is not UTF-8: {e}"))?;
        let (rest, strand) =
            text.rsplit_once(':').ok_or_else(|| anyhow!("invalid rw value '{text}'"))?;
        let (contig, range) =
            rest.rsplit_once(':').ok_or_else(|| anyhow!("invalid rw value '{text}'"))?;
        let (start, end) =
            range.split_once('-').ok_or_else(|| anyhow!("invalid rw range in '{text}'"))?;
        let reverse = match strand {
            "+" => false,
            "-" => true,
            _ => bail!("invalid rw strand '{strand}' in '{text}'"),
        };
        if contig.is_empty() {
            bail!("empty contig in rw value '{text}'");
        }
        let start = parse_position(start).ok_or_else(|| anyhow!("invalid rw start in '{text}'"))?;
        let end = parse_position(end).ok_or_else(|| anyhow!("invalid rw end in '{text}'"))?;
        if start > end {
            bail!("rw start {start} is greater than end {end} in '{text}'");
        }
        Ok(Self { contig, start, end, reverse })
    }

    /// Appends the `rw` tag value for this window to `out` (inverse of [`Self::parse`]).
    pub fn write(&self, out: &mut Vec<u8>) {
        write_value(self.contig.as_bytes(), self.start, self.end, self.reverse, out);
    }

    /// Returns the union of two windows on the same contig, or `None` if the
    /// contigs differ. The result spans the minimum start to the maximum end and
    /// keeps `self`'s strand.
    #[must_use]
    pub fn union(&self, other: &Self) -> Option<Self> {
        if self.contig != other.contig {
            return None;
        }
        Some(Self {
            contig: self.contig,
            start: self.start.min(other.start),
            end: self.end.max(other.end),
            reverse: self.reverse,
        })
    }

    /// Returns true iff the interval `[start, end]` is on the same contig and
    /// strand as this window and intersects it (closed intervals).
    #[must_use]
    pub fn overlaps(&self, contig: &str, start: u32, end: u32, reverse: bool) -> bool {
        self.contig == contig && self.reverse == reverse && start <= self.end && end >= self.start
    }
}

/// Appends the `rw` tag value `<contig>:<start>-<end>:<strand>` to `out`: the one
/// writer of the format, for a [`RealignWindow`] or for a window held as parts.
pub fn write_value(contig: &[u8], start: u32, end: u32, reverse: bool, out: &mut Vec<u8>) {
    out.extend_from_slice(contig);
    let strand = if reverse { '-' } else { '+' };
    // Writing to a `Vec` cannot fail.
    write!(out, ":{start}-{end}:{strand}").expect("writing to a Vec cannot fail");
}

/// Parses a 1-based coordinate: ASCII digits only (no sign), and not 0.
fn parse_position(field: &str) -> Option<u32> {
    parse_digits(field).filter(|&p| p > 0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    fn win(contig: &str, start: u32, end: u32, reverse: bool) -> RealignWindow<'_> {
        RealignWindow { contig, start, end, reverse }
    }

    #[rstest]
    #[case::forward(b"chr1:100-250:+", win("chr1", 100, 250, false))]
    #[case::reverse(b"chr1:100-250:-", win("chr1", 100, 250, true))]
    #[case::contig_with_colons(
        b"HLA-A*01:01:01:01:100-200:-",
        win("HLA-A*01:01:01:01", 100, 200, true)
    )]
    fn roundtrip(#[case] value: &[u8], #[case] expected: RealignWindow<'_>) {
        let parsed = RealignWindow::parse(value).unwrap();
        assert_eq!(parsed, expected);
        let mut out = Vec::new();
        parsed.write(&mut out);
        assert_eq!(out, value);
    }

    #[rstest]
    #[case::missing_strand(b"chr1:100-250", "invalid rw value")]
    #[case::start_gt_end(b"chr1:300-250:+", "rw start 300 is greater than end 250")]
    #[case::bad_strand(b"chr1:100-250:x", "invalid rw strand 'x'")]
    #[case::non_numeric(b"chr1:a-250:+", "invalid rw start")]
    #[case::zero_start(b"chr1:0-250:+", "invalid rw start")]
    #[case::signed_end(b"chr1:1-+250:+", "invalid rw end")]
    #[case::empty_contig(b":1-250:+", "empty contig")]
    fn parse_rejects(#[case] value: &[u8], #[case] message: &str) {
        let err = RealignWindow::parse(value).expect_err("malformed rw").to_string();
        assert!(err.contains(message), "{err}");
    }

    #[rstest]
    #[case::overlapping(
        win("chr1", 100, 200, false),
        win("chr1", 150, 300, false),
        Some(win("chr1", 100, 300, false))
    )]
    #[case::disjoint(
        win("chr1", 100, 200, false),
        win("chr1", 400, 500, false),
        Some(win("chr1", 100, 500, false))
    )]
    #[case::keeps_self_strand(
        win("chr1", 100, 200, true),
        win("chr1", 50, 120, false),
        Some(win("chr1", 50, 200, true))
    )]
    #[case::other_contig(win("chr1", 100, 200, false), win("chr2", 100, 200, false), None)]
    fn union_cases(
        #[case] a: RealignWindow<'_>,
        #[case] b: RealignWindow<'_>,
        #[case] expected: Option<RealignWindow<'_>>,
    ) {
        assert_eq!(a.union(&b), expected);
    }

    #[rstest]
    #[case::inside("chr1", 120, 180, false, true)]
    #[case::partial_left("chr1", 50, 100, false, true)]
    #[case::adjacent_left_no_overlap("chr1", 50, 99, false, false)]
    #[case::touches_right_edge("chr1", 200, 300, false, true)]
    #[case::adjacent_no_overlap("chr1", 201, 300, false, false)]
    #[case::other_contig("chr2", 120, 180, false, false)]
    #[case::other_strand("chr1", 120, 180, true, false)]
    fn overlaps_cases(
        #[case] contig: &str,
        #[case] start: u32,
        #[case] end: u32,
        #[case] reverse: bool,
        #[case] expected: bool,
    ) {
        assert_eq!(win("chr1", 100, 200, false).overlaps(contig, start, end, reverse), expected);
    }
}
