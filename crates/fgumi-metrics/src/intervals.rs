//! Genomic interval parsing and overlap checking.
//!
//! Supports both BED (0-based half-open) and Picard interval list (1-based closed)
//! file formats with automatic format detection.

use anyhow::Result;
use std::path::Path;

/// Genomic interval for filtering, stored as 0-based half-open coordinates.
#[derive(Clone, Debug)]
pub struct Interval {
    /// Reference sequence name (e.g. "chr1").
    pub ref_name: String,
    /// 0-based start position (inclusive).
    pub start: i32,
    /// 0-based end position (exclusive).
    pub end: i32,
}

/// Parses an intervals file in BED or Picard interval list format.
///
/// Auto-detects the format: if any line starts with `@`, the file is treated as a
/// Picard interval list (1-based closed coordinates with a SAM header); otherwise
/// it is treated as BED (0-based half-open coordinates).
///
/// Intervals are stored internally using BED conventions (0-based half-open).
///
/// # Errors
///
/// Returns an error if the file cannot be read or lines cannot be parsed.
///
/// # Panics
///
/// This function will not panic in practice. The only `.unwrap()` call is on
/// `fields.next()` after splitting a non-empty, non-comment, non-header line,
/// which is guaranteed to yield at least one element.
pub fn parse_intervals(path: &Path) -> Result<Vec<Interval>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut intervals = Vec::new();
    let mut is_interval_list = false;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // Skip SAM header lines (interval list format)
        if line.starts_with('@') {
            is_interval_list = true;
            continue;
        }

        let mut fields = line.splitn(4, '\t');
        let ref_name = fields.next().unwrap(); // guaranteed non-empty (line isn't empty)
        let start_str = fields.next();
        let end_str = fields.next();

        let (Some(start_str), Some(end_str)) = (start_str, end_str) else {
            let fmt = if is_interval_list { "interval list" } else { "BED" };
            anyhow::bail!("Invalid {fmt} line (needs at least 3 fields): {line}");
        };

        if is_interval_list {
            // Picard interval list: chr start end strand name (1-based, closed)
            let start: i32 = start_str
                .parse::<i32>()
                .map_err(|_| anyhow::anyhow!("Invalid start position: {start_str}"))?
                - 1; // Convert 1-based to 0-based
            let end: i32 =
                end_str.parse().map_err(|_| anyhow::anyhow!("Invalid end position: {end_str}"))?;
            // end stays the same: 1-based closed end == 0-based half-open end
            intervals.push(Interval { ref_name: ref_name.to_string(), start, end });
        } else {
            // BED format: chr start end [name] [score] [strand] (0-based, half-open)
            let start: i32 = start_str
                .parse()
                .map_err(|_| anyhow::anyhow!("Invalid start position: {start_str}"))?;
            let end: i32 =
                end_str.parse().map_err(|_| anyhow::anyhow!("Invalid end position: {end_str}"))?;
            intervals.push(Interval { ref_name: ref_name.to_string(), start, end });
        }
    }

    Ok(intervals)
}

/// Checks if a genomic region overlaps any provided interval.
///
/// Determines whether a genomic region overlaps with any of the specified intervals.
/// If no intervals are provided, all regions are considered to overlap (no filtering).
///
/// Positions are 1-based inclusive (SAM convention). Intervals are stored as
/// 0-based half-open (BED convention). The overlap test converts the 1-based
/// region `[start, end]` to 0-based half-open `[start-1, end)`, then checks
/// `(start-1) < iv.end && iv.start < end`, which simplifies (for integers) to
/// `start <= iv.end && iv.start < end`.
///
/// # Returns
///
/// `true` if the region overlaps any interval or if no intervals are provided,
/// `false` if the region is unmapped or does not overlap any interval.
#[must_use]
pub fn overlaps_intervals(
    ref_name: Option<&str>,
    position: Option<i32>,
    end_position: Option<i32>,
    intervals: &[Interval],
) -> bool {
    if intervals.is_empty() {
        return true; // No filtering if no intervals provided
    }

    if let (Some(ref_name), Some(start), Some(end)) = (ref_name, position, end_position) {
        intervals.iter().any(|interval| {
            interval.ref_name == ref_name && start <= interval.end && interval.start < end
        })
    } else {
        false // Unmapped reads or reads without proper coordinates don't overlap any interval
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_parse_bed_intervals() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.bed");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "chr1\t100\t200").unwrap();
        writeln!(f, "chr2\t300\t400").unwrap();

        let intervals = parse_intervals(&path).unwrap();
        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].ref_name, "chr1");
        assert_eq!(intervals[0].start, 100);
        assert_eq!(intervals[0].end, 200);
        assert_eq!(intervals[1].ref_name, "chr2");
        assert_eq!(intervals[1].start, 300);
        assert_eq!(intervals[1].end, 400);
    }

    #[test]
    fn test_parse_picard_interval_list() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.interval_list");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "@HD\tVN:1.0").unwrap();
        writeln!(f, "@SQ\tSN:chr1\tLN:1000").unwrap();
        writeln!(f, "chr1\t101\t200\t+\ttarget1").unwrap();

        let intervals = parse_intervals(&path).unwrap();
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].ref_name, "chr1");
        // 1-based 101 -> 0-based 100
        assert_eq!(intervals[0].start, 100);
        // 1-based closed 200 == 0-based half-open 200
        assert_eq!(intervals[0].end, 200);
    }

    #[test]
    fn test_overlaps_intervals_no_intervals() {
        assert!(overlaps_intervals(Some("chr1"), Some(100), Some(200), &[]));
    }

    #[test]
    fn test_overlaps_intervals_hit() {
        let intervals = vec![Interval { ref_name: "chr1".to_string(), start: 50, end: 150 }];
        assert!(overlaps_intervals(Some("chr1"), Some(100), Some(200), &intervals));
    }

    #[test]
    fn test_overlaps_intervals_miss() {
        let intervals = vec![Interval { ref_name: "chr1".to_string(), start: 300, end: 400 }];
        assert!(!overlaps_intervals(Some("chr1"), Some(100), Some(200), &intervals));
    }

    #[test]
    fn test_overlaps_intervals_wrong_chrom() {
        let intervals = vec![Interval { ref_name: "chr2".to_string(), start: 50, end: 150 }];
        assert!(!overlaps_intervals(Some("chr1"), Some(100), Some(200), &intervals));
    }

    #[test]
    fn test_overlaps_intervals_unmapped() {
        let intervals = vec![Interval { ref_name: "chr1".to_string(), start: 50, end: 150 }];
        assert!(!overlaps_intervals(None, None, None, &intervals));
    }
}
