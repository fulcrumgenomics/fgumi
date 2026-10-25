//! Parser and writer for bwa's `XA` (alternative hits) tag.
//!
//! The value is a list of hits, each terminated by `;`, of the form
//! `chr,[+-]pos,CIGAR,NM;` where `pos` is 1-based. `bwa-mem3 mem -u` appends
//! each hit's alignment score and mapping quality (`chr,[+-]pos,CIGAR,NM,score,mapq;`)
//! and still names the tag `XA`; both forms are accepted and written back as read.

use std::fmt::Write as _;

use anyhow::{Result, anyhow, bail};

/// A single alternative alignment hit from an `XA` tag.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct XaHit {
    /// Reference contig name.
    pub contig: String,
    /// 1-based leftmost reference position.
    pub start: u32,
    /// True if the hit is on the reverse strand.
    pub reverse: bool,
    /// CIGAR string of the hit.
    pub cigar: String,
    /// Edit distance of the hit.
    pub nm: u32,
    /// The hit's alignment score and mapping quality, present when the aligner
    /// wrote the six-field form (`bwa-mem3 mem -u`).
    pub score_mapq: Option<(i32, u32)>,
}

impl XaHit {
    /// The hit's CIGAR as BAM-packed ops, or `None` if it is not a valid CIGAR
    /// for a mapped alignment (see [`fgumi_raw_bam::parse_cigar_string_ops`]:
    /// unknown or zero-length operations, a run length past BAM's 28-bit field,
    /// misplaced clips, or no reference-consuming operation).
    #[must_use]
    pub fn cigar_ops(&self) -> Option<Vec<u32>> {
        fgumi_raw_bam::parse_cigar_string_ops(&self.cigar)
    }

    /// The 1-based inclusive reference end of the hit, given its parsed CIGAR
    /// `ops` ([`Self::cigar_ops`]). Saturates rather than overflowing.
    #[must_use]
    pub fn end(&self, ops: &[u32]) -> u32 {
        let span = u32::try_from(fgumi_raw_bam::saturating_reference_length(ops)).unwrap_or(0);
        self.start.saturating_add(span).saturating_sub(1).max(self.start)
    }
}

/// Parses an unsigned decimal field made only of ASCII digits (no sign, no spaces).
pub(crate) fn parse_digits<T: std::str::FromStr>(field: &str) -> Option<T> {
    if field.is_empty() || !field.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    field.parse().ok()
}

/// Parses an `XA` tag value into its hits.
///
/// # Errors
/// Returns an error if the value is not UTF-8 or any hit is malformed: a field
/// count other than 4 or 6, a missing strand sign, a position that is not a
/// positive integer (`XA` positions are 1-based), or a non-numeric `NM`, score
/// or mapping quality. CIGARs are not checked here; see [`XaHit::cigar_ops`].
pub fn parse_xa(value: &[u8]) -> Result<Vec<XaHit>> {
    let text = std::str::from_utf8(value).map_err(|e| anyhow!("XA value is not UTF-8: {e}"))?;
    let mut hits = Vec::new();
    for entry in text.split(';').filter(|e| !e.is_empty()) {
        let fields: Vec<&str> = entry.split(',').collect();
        let (contig, pos, cigar, nm, extra) = match fields[..] {
            [contig, pos, cigar, nm] => (contig, pos, cigar, nm, None),
            [contig, pos, cigar, nm, score, mapq] => (contig, pos, cigar, nm, Some((score, mapq))),
            _ => bail!("invalid XA hit '{entry}': expected 4 or 6 comma-separated fields"),
        };
        let (reverse, digits) = match pos.as_bytes().first() {
            Some(b'+') => (false, &pos[1..]),
            Some(b'-') => (true, &pos[1..]),
            _ => bail!("invalid XA strand in '{entry}'"),
        };
        let start: u32 =
            parse_digits(digits).ok_or_else(|| anyhow!("invalid XA position in '{entry}'"))?;
        if start == 0 {
            bail!("invalid XA position 0 in '{entry}' (XA positions are 1-based)");
        }
        let nm: u32 = parse_digits(nm).ok_or_else(|| anyhow!("invalid XA NM in '{entry}'"))?;
        let score_mapq = match extra {
            None => None,
            Some((score, mapq)) => {
                let score: i32 =
                    score.parse().map_err(|_| anyhow!("invalid XA score in '{entry}'"))?;
                let mapq: u32 =
                    parse_digits(mapq).ok_or_else(|| anyhow!("invalid XA mapq in '{entry}'"))?;
                Some((score, mapq))
            }
        };
        hits.push(XaHit {
            contig: contig.to_string(),
            start,
            reverse,
            cigar: cigar.to_string(),
            nm,
            score_mapq,
        });
    }
    Ok(hits)
}

/// Appends hits to `out` in `XA` format (inverse of [`parse_xa`]); every hit,
/// including the last, is terminated by `;`. A hit with a score and mapping
/// quality is written in the six-field form.
pub fn write_xa(hits: &[XaHit], out: &mut Vec<u8>) {
    for h in hits {
        let mut line = format!(
            "{},{}{},{},{}",
            h.contig,
            if h.reverse { '-' } else { '+' },
            h.start,
            h.cigar,
            h.nm
        );
        if let Some((score, mapq)) = h.score_mapq {
            // Writing to a `String` cannot fail.
            let _ = write!(line, ",{score},{mapq}");
        }
        line.push(';');
        out.extend_from_slice(line.as_bytes());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    fn hit(contig: &str, start: u32, reverse: bool, cigar: &str, nm: u32) -> XaHit {
        XaHit {
            contig: contig.to_string(),
            start,
            reverse,
            cigar: cigar.to_string(),
            nm,
            score_mapq: None,
        }
    }

    fn hit6(
        contig: &str,
        start: u32,
        reverse: bool,
        cigar: &str,
        nm: u32,
        score_mapq: (i32, u32),
    ) -> XaHit {
        XaHit { score_mapq: Some(score_mapq), ..hit(contig, start, reverse, cigar, nm) }
    }

    #[rstest]
    #[case::one_hit(b"chr2,+1000,150M,1;", vec![hit("chr2", 1000, false, "150M", 1)])]
    #[case::two_hits(
        b"chr2,+1000,150M,1;chr7,+55,20S130M,2;",
        vec![hit("chr2", 1000, false, "150M", 1), hit("chr7", 55, false, "20S130M", 2)]
    )]
    #[case::reverse_hit(b"chr3,-500,50M2D100M,3;", vec![hit("chr3", 500, true, "50M2D100M", 3)])]
    #[case::six_field(
        b"chr1,+2001,300M,0,300,0;",
        vec![hit6("chr1", 2001, false, "300M", 0, (300, 0))]
    )]
    #[case::six_field_negative_score(
        b"chr2,-7,10M,1,-3,60;",
        vec![hit6("chr2", 7, true, "10M", 1, (-3, 60))]
    )]
    fn xa_roundtrip(#[case] value: &[u8], #[case] expected: Vec<XaHit>) {
        let parsed = parse_xa(value).unwrap();
        assert_eq!(parsed, expected);
        let mut out = Vec::new();
        write_xa(&parsed, &mut out);
        assert_eq!(out, value);
    }

    #[rstest]
    #[case::match_only("150M", Some(1149))]
    #[case::deletion("50M2D100M", Some(1151))]
    #[case::soft_clip("20S130M", Some(1129))]
    #[case::saturates(&"268435455M".repeat(17), Some(u32::try_from(i32::MAX).unwrap() + 999))]
    #[case::oversized_run_length("268435456M", None)]
    #[case::trailing_length("10M5", None)]
    #[case::unknown_operator("10Q5M", None)]
    #[case::zero_length("0M", None)]
    #[case::internal_soft_clip("10M5S10M", None)]
    fn xa_end(#[case] cigar: &str, #[case] expected: Option<u32>) {
        let h = hit("chr1", 1000, false, cigar, 0);
        assert_eq!(h.cigar_ops().map(|ops| h.end(&ops)), expected);
    }

    #[rstest]
    #[case::missing_nm(b"chr2,+1000,150M;", "expected 4 or 6 comma-separated fields")]
    #[case::five_fields(b"chr2,+1000,150M,1,150;", "expected 4 or 6 comma-separated fields")]
    #[case::bad_strand(b"chr2,1000,150M,1;", "invalid XA strand")]
    #[case::zero_position(b"chr2,+0,150M,1;", "invalid XA position 0")]
    #[case::signed_position(b"chr2,++5,150M,1;", "invalid XA position in")]
    #[case::bad_nm(b"chr2,+5,150M,x;", "invalid XA NM")]
    #[case::bad_score(b"chr2,+5,150M,1,x,0;", "invalid XA score")]
    #[case::bad_mapq(b"chr2,+5,150M,1,150,-1;", "invalid XA mapq")]
    fn xa_rejects(#[case] value: &[u8], #[case] message: &str) {
        let err = parse_xa(value).expect_err("malformed XA").to_string();
        assert!(err.contains(message), "{err}");
    }
}
