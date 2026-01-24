//! Metrics for the `correct` command.
//!
//! This module provides metrics for tracking UMI correction operations.

use anyhow::{Context, Result};
use fgoxide::io::DelimFile;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::path::PathBuf;

use super::Metric;

/// Serializes an f64 as an integer string if it's a whole number.
///
/// This matches fgbio's output format where 0.0 is output as "0".
#[allow(clippy::trivially_copy_pass_by_ref)] // serde requires &T signature
fn serialize_f64_as_int<S>(value: &f64, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let value = *value;
    // Check for NaN and Infinity first
    if value.is_nan() {
        return serializer.serialize_str("NaN");
    }
    if value.is_infinite() {
        return serializer.serialize_str(if value > 0.0 { "Infinity" } else { "-Infinity" });
    }

    // If the value is a whole number, serialize without decimal point
    if value.fract() == 0.0 && value.abs() < (i64::MAX as f64) {
        serializer.serialize_str(&format!("{}", value as i64))
    } else {
        serializer.serialize_str(&format!("{value}"))
    }
}

/// Deserializes an f64 from a string, handling special values like NaN and Infinity.
fn deserialize_f64_from_str<'de, D>(deserializer: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    let s: String = String::deserialize(deserializer)?;
    match s.as_str() {
        "NaN" => Ok(f64::NAN),
        "Infinity" => Ok(f64::INFINITY),
        "-Infinity" => Ok(f64::NEG_INFINITY),
        _ => s.parse().map_err(serde::de::Error::custom),
    }
}

/// Metrics tracking how well observed UMIs match expected UMI sequences.
///
/// These metrics are generated per-UMI and track the distribution of match types
/// (perfect matches, single mismatches, etc.) for each expected UMI.
///
/// # Fields
///
/// * `umi` - The expected/corrected UMI sequence (or all Ns for unmatched)
/// * `total_matches` - Total UMI sequences matched/corrected to this UMI
/// * `perfect_matches` - Number of reads with zero mismatches
/// * `one_mismatch_matches` - Number of reads with exactly one mismatch
/// * `two_mismatch_matches` - Number of reads with exactly two mismatches
/// * `other_matches` - Number of reads with three or more mismatches
/// * `fraction_of_matches` - Proportion of all reads matching this UMI
/// * `representation` - Ratio of this UMI's count to the mean count across all UMIs
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UmiCorrectionMetrics {
    /// The corrected UMI sequence (or all Ns for unmatched).
    pub umi: String,

    /// The number of UMI sequences that matched/were corrected to this UMI.
    pub total_matches: u64,

    /// The number of UMI sequences that were perfect matches to this UMI.
    pub perfect_matches: u64,

    /// The number of UMI sequences that matched with a single mismatch.
    pub one_mismatch_matches: u64,

    /// The number of UMI sequences that matched with two mismatches.
    pub two_mismatch_matches: u64,

    /// The number of UMI sequences that matched with three or more mismatches.
    pub other_matches: u64,

    /// The fraction of all UMIs that matched or were corrected to this UMI.
    #[serde(
        serialize_with = "serialize_f64_as_int",
        deserialize_with = "deserialize_f64_from_str"
    )]
    pub fraction_of_matches: f64,

    /// The `total_matches` for this UMI divided by the mean `total_matches` for all UMIs.
    #[serde(
        serialize_with = "serialize_f64_as_int",
        deserialize_with = "deserialize_f64_from_str"
    )]
    pub representation: f64,
}

impl UmiCorrectionMetrics {
    /// Creates a new metrics instance for the given UMI.
    ///
    /// All count fields are initialized to zero, and calculated fields to 0.0.
    ///
    /// # Arguments
    ///
    /// * `umi` - The UMI sequence these metrics are for
    #[must_use]
    pub fn new(umi: String) -> Self {
        Self {
            umi,
            total_matches: 0,
            perfect_matches: 0,
            one_mismatch_matches: 0,
            two_mismatch_matches: 0,
            other_matches: 0,
            fraction_of_matches: 0.0,
            representation: 0.0,
        }
    }

    /// Writes metrics to a TSV file.
    ///
    /// The output format is tab-separated with a header row followed by
    /// one row per UMI.
    ///
    /// # Arguments
    ///
    /// * `metrics` - Slice of metrics to write
    /// * `path` - Path to the output file
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be created or written to.
    pub fn write_metrics(metrics: &[Self], path: &PathBuf) -> Result<()> {
        DelimFile::default()
            .write_tsv(path, metrics)
            .with_context(|| format!("Failed to write UMI correction metrics: {}", path.display()))
    }

    /// Reads metrics from a TSV file.
    ///
    /// Parses a metrics file created by `write_metrics`.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the metrics file to read
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The file cannot be read
    /// - The file format is invalid
    /// - Any values cannot be parsed
    pub fn read_metrics(path: &PathBuf) -> Result<Vec<Self>> {
        DelimFile::default()
            .read_tsv(path)
            .with_context(|| format!("Failed to read UMI correction metrics: {}", path.display()))
    }
}

impl Default for UmiCorrectionMetrics {
    fn default() -> Self {
        Self::new(String::new())
    }
}

impl Metric for UmiCorrectionMetrics {
    fn metric_name() -> &'static str {
        "UMI correction"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_umi_correction_metrics_new() {
        let metrics = UmiCorrectionMetrics::new("ACGT".to_string());
        assert_eq!(metrics.umi, "ACGT");
        assert_eq!(metrics.total_matches, 0);
        assert_eq!(metrics.perfect_matches, 0);
        assert_eq!(metrics.one_mismatch_matches, 0);
        assert_eq!(metrics.two_mismatch_matches, 0);
        assert_eq!(metrics.other_matches, 0);
        assert!(metrics.fraction_of_matches.abs() < f64::EPSILON);
        assert!(metrics.representation.abs() < f64::EPSILON);
    }

    #[test]
    fn test_umi_correction_metrics_default() {
        let metrics = UmiCorrectionMetrics::default();
        assert!(metrics.umi.is_empty());
        assert_eq!(metrics.total_matches, 0);
    }

    #[test]
    fn test_write_and_read_metrics() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let path = temp_file.path().to_path_buf();

        let mut m1 = UmiCorrectionMetrics::new("AAAA".to_string());
        m1.total_matches = 100;
        m1.perfect_matches = 80;
        m1.one_mismatch_matches = 15;
        m1.two_mismatch_matches = 5;
        m1.fraction_of_matches = 0.5;
        m1.representation = 1.2;

        let mut m2 = UmiCorrectionMetrics::new("TTTT".to_string());
        m2.total_matches = 100;
        m2.perfect_matches = 90;
        m2.one_mismatch_matches = 10;
        m2.fraction_of_matches = 0.5;
        m2.representation = 0.8;

        let metrics = vec![m1, m2];
        UmiCorrectionMetrics::write_metrics(&metrics, &path)?;

        let read_metrics = UmiCorrectionMetrics::read_metrics(&path)?;
        assert_eq!(read_metrics.len(), 2);
        assert_eq!(read_metrics[0].umi, "AAAA");
        assert_eq!(read_metrics[0].total_matches, 100);
        assert_eq!(read_metrics[1].umi, "TTTT");
        assert_eq!(read_metrics[1].perfect_matches, 90);

        Ok(())
    }

    #[test]
    fn test_metric_trait_impl() {
        assert_eq!(UmiCorrectionMetrics::metric_name(), "UMI correction");
    }

    #[test]
    fn test_nan_infinity_serialization() -> Result<()> {
        let temp_file = NamedTempFile::new()?;
        let path = temp_file.path().to_path_buf();

        // Create metrics with NaN and Infinity values (edge cases)
        let mut m1 = UmiCorrectionMetrics::new("AAAA".to_string());
        m1.total_matches = 0;
        m1.fraction_of_matches = f64::NAN; // 0/0 case
        m1.representation = f64::NAN;

        let mut m2 = UmiCorrectionMetrics::new("TTTT".to_string());
        m2.total_matches = 100;
        m2.fraction_of_matches = 1.0;
        m2.representation = f64::INFINITY; // n/0 case

        let metrics = vec![m1, m2];
        UmiCorrectionMetrics::write_metrics(&metrics, &path)?;

        // Read back and verify NaN/Infinity are preserved
        let read_metrics = UmiCorrectionMetrics::read_metrics(&path)?;
        assert_eq!(read_metrics.len(), 2);

        // NaN should remain NaN (NaN != NaN, so use is_nan())
        assert!(read_metrics[0].fraction_of_matches.is_nan());
        assert!(read_metrics[0].representation.is_nan());

        // Infinity should remain Infinity
        assert!(read_metrics[1].representation.is_infinite());
        assert!(read_metrics[1].representation > 0.0); // Positive infinity

        Ok(())
    }
}
