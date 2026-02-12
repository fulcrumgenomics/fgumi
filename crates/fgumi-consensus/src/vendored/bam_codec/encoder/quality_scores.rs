use noodles::sam::alignment::record::QualityScores;
use std::{io, iter};

use super::num::write_u8;

const MAX_SCORE: u8 = 93;
const MISSING: u8 = 255;

pub(super) fn write_quality_scores<S>(
    dst: &mut Vec<u8>,
    base_count: usize,
    quality_scores: S,
) -> io::Result<()>
where
    S: QualityScores,
{

    if quality_scores.len() == base_count {
        for result in quality_scores.iter() {
            let n = result?;

            if is_valid_score(n) {
                write_u8(dst, n);
            } else {
                return Err(io::Error::from(io::ErrorKind::InvalidInput));
            }
        }
    } else if quality_scores.is_empty() {
        dst.extend(iter::repeat_n(MISSING, base_count));
    } else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "sequence-quality scores length mismatch: expected {}, got {}",
                base_count,
                quality_scores.len()
            ),
        ));
    }

    Ok(())
}

/// Fast path: encode quality scores from a slice (no dynamic dispatch).
///
/// This bypasses the trait-based iterator and uses bulk operations:
/// 1. Validates all scores in one pass (LLVM can auto-vectorize this)
/// 2. Uses `extend_from_slice` for a single memcpy
#[inline]
pub(super) fn write_quality_scores_from_slice(
    dst: &mut Vec<u8>,
    base_count: usize,
    scores: &[u8],
) -> io::Result<()> {

    if scores.len() == base_count {
        // Validate entire slice first - this can be auto-vectorized by LLVM
        if scores.iter().all(|&s| s <= MAX_SCORE) {
            dst.extend_from_slice(scores); // Single memcpy!
            Ok(())
        } else {
            Err(io::Error::new(io::ErrorKind::InvalidInput, "quality score out of range [0, 93]"))
        }
    } else if scores.is_empty() {
        // Fill with MISSING (0xFF)
        dst.resize(dst.len() + base_count, MISSING);
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "sequence-quality scores length mismatch: expected {}, got {}",
                base_count,
                scores.len()
            ),
        ))
    }
}

fn is_valid_score(score: u8) -> bool {
    // § 4.2.3 "SEQ and QUAL encoding" (2023-05-24): "Base qualities are stored as bytes in the
    // range [0, 93]..."
    score <= MAX_SCORE
}

#[cfg(test)]
mod tests {
    use noodles::sam::alignment::record_buf::QualityScores as QualityScoresBuf;

    use super::*;

    #[test]
    fn test_write_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            buf: &mut Vec<u8>,
            base_count: usize,
            quality_scores: &QualityScoresBuf,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_quality_scores(buf, base_count, quality_scores)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 0, &QualityScoresBuf::default(), &[])?;
        t(&mut buf, 4, &QualityScoresBuf::default(), &[0xff, 0xff, 0xff, 0xff])?;

        let quality_scores = QualityScoresBuf::from(vec![45, 35, 43, 50]);
        t(&mut buf, 4, &quality_scores, &[0x2d, 0x23, 0x2b, 0x32])?;

        let quality_scores = QualityScoresBuf::from(vec![45, 35, 43, 50]);
        buf.clear();
        assert!(matches!(
            write_quality_scores(&mut buf, 3, &quality_scores),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        // Out-of-range score (> 93) should be rejected
        let quality_scores = QualityScoresBuf::from(vec![45, 94, 43, 50]);
        buf.clear();
        assert!(matches!(
            write_quality_scores(&mut buf, 4, &quality_scores),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_quality_scores_from_slice() {
        let mut buf = Vec::new();

        // Valid scores
        write_quality_scores_from_slice(&mut buf, 4, &[45, 35, 43, 50]).unwrap();
        assert_eq!(buf, [45, 35, 43, 50]);

        // Empty scores (fill with 0xFF)
        buf.clear();
        write_quality_scores_from_slice(&mut buf, 3, &[]).unwrap();
        assert_eq!(buf, [0xff, 0xff, 0xff]);

        // Zero base count with empty scores
        buf.clear();
        write_quality_scores_from_slice(&mut buf, 0, &[]).unwrap();
        assert!(buf.is_empty());

        // Invalid score (> 93)
        buf.clear();
        assert!(write_quality_scores_from_slice(&mut buf, 2, &[45, 94]).is_err());

        // Length mismatch
        buf.clear();
        assert!(write_quality_scores_from_slice(&mut buf, 3, &[45, 35]).is_err());
    }
}
