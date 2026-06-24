//! Input validation utilities
//!
//! This module provides common validation functions for command-line parameters,
//! file paths, and SAM tags with consistent error messages.
//!
//! All validation functions now use structured error types from [`crate::errors`] to provide
//! rich contextual information when validation fails.
//!
//! `validate_file_exists` and `parse_memory_size` are re-exported from `fgumi-cli-common`
//! (single source of truth); umbrella-only validators are defined here.

pub use fgumi_cli_common::{parse_memory_size, validate_file_exists};

use crate::errors::{FgumiError, Result};
use std::fmt::Display;
use std::path::Path;

/// Validate that multiple files exist
///
/// # Arguments
/// * `files` - Slice of (path, description) tuples
///
/// # Errors
/// Returns an error for the first file that doesn't exist
///
/// # Example
/// ```no_run
/// use fgumi_lib::validation::validate_files_exist;
/// use std::path::PathBuf;
///
/// let files = vec![
///     (PathBuf::from("input.bam"), "Input BAM"),
///     (PathBuf::from("ref.fa"), "Reference"),
/// ];
/// validate_files_exist(&files).unwrap();
/// ```
pub fn validate_files_exist<P: AsRef<Path>>(files: &[(P, &str)]) -> Result<()> {
    for (path, desc) in files {
        validate_file_exists(path, desc)?;
    }
    Ok(())
}

/// Validate that max >= min for optional max values
///
/// # Arguments
/// * `min_val` - Minimum value
/// * `max_val` - Optional maximum value
/// * `min_name` - Name of minimum parameter for error messages
/// * `max_name` - Name of maximum parameter for error messages
///
/// # Errors
/// Returns an error if max < min
///
/// # Example
/// ```
/// use fgumi_lib::validation::validate_min_max;
///
/// // Valid: max >= min
/// validate_min_max(1, Some(10), "min-reads", "max-reads").unwrap();
///
/// // Valid: max is None
/// validate_min_max(1, None, "min-reads", "max-reads").unwrap();
///
/// // Invalid: max < min
/// let result = validate_min_max(10, Some(5), "min-reads", "max-reads");
/// assert!(result.is_err());
/// ```
#[allow(clippy::needless_pass_by_value)]
pub fn validate_min_max<T: Ord + Display>(
    min_val: T,
    max_val: Option<T>,
    min_name: &str,
    max_name: &str,
) -> Result<()> {
    if let Some(max) = max_val {
        if max < min_val {
            return Err(FgumiError::InvalidParameter {
                parameter: max_name.to_string(),
                reason: format!("{max_name} ({max}) must be >= {min_name} ({min_val})"),
            });
        }
    }
    Ok(())
}

/// Validate that an error rate is in the valid range [0.0, 1.0]
///
/// # Arguments
/// * `rate` - Error rate to validate
/// * `_name` - Parameter name, accepted for call-site symmetry with
///   `validate_min_max`/`validate_positive` but currently unused: the returned
///   `FgumiError::InvalidFrequency` does not interpolate it.
///
/// # Errors
/// Returns an error if the rate is not in [0.0, 1.0]
///
/// # Example
/// ```
/// use fgumi_lib::validation::validate_error_rate;
///
/// validate_error_rate(0.01, "error-rate-pre-umi").unwrap();
/// validate_error_rate(1.0, "error-rate-post-umi").unwrap();
///
/// let result = validate_error_rate(1.5, "error-rate");
/// assert!(result.is_err());
/// ```
pub fn validate_error_rate(rate: f64, _name: &str) -> Result<()> {
    if !(0.0..=1.0).contains(&rate) {
        return Err(FgumiError::InvalidFrequency { value: rate, min: 0.0, max: 1.0 });
    }
    Ok(())
}

/// Validate that a quality score is in the valid Phred range [0, 93]
///
/// # Arguments
/// * `quality` - Quality score to validate
/// * `_name` - Parameter name, accepted for call-site symmetry with
///   `validate_min_max`/`validate_positive` but currently unused: the returned
///   `FgumiError::InvalidQuality` does not interpolate it.
///
/// # Errors
/// Returns an error if the quality is not in [0, 93]
///
/// # Example
/// ```
/// use fgumi_lib::validation::validate_quality_score;
///
/// validate_quality_score(30, "min-base-quality").unwrap();
///
/// let result = validate_quality_score(100, "min-base-quality");
/// assert!(result.is_err());
/// ```
pub fn validate_quality_score(quality: u8, _name: &str) -> Result<()> {
    if quality > 93 {
        return Err(FgumiError::InvalidQuality { value: quality, max: 93 });
    }
    Ok(())
}

/// Validate that a value is positive (> 0)
///
/// # Arguments
/// * `value` - Value to validate
/// * `name` - Name of the parameter for error messages
///
/// # Errors
/// Returns an error if the value is not positive
///
/// # Example
/// ```
/// use fgumi_lib::validation::validate_positive;
///
/// validate_positive(10, "min-reads").unwrap();
///
/// let result = validate_positive(0, "min-reads");
/// assert!(result.is_err());
/// ```
#[allow(clippy::needless_pass_by_value)]
pub fn validate_positive<T: Ord + Display + Default>(value: T, name: &str) -> Result<()> {
    if value <= T::default() {
        return Err(FgumiError::InvalidParameter {
            parameter: name.to_string(),
            reason: format!("Must be positive (> 0), got: {value}"),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use std::path::PathBuf;
    use tempfile::NamedTempFile;

    #[test]
    fn test_validate_file_exists_valid() {
        let temp_file = NamedTempFile::new().expect("creating temp file/dir should succeed");
        validate_file_exists(temp_file.path(), "Test file")
            .expect("file validation should succeed");
    }

    #[test]
    fn test_validate_file_exists_invalid() {
        let result = validate_file_exists("/nonexistent/file.bam", "Input file");
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("Input file"));
        assert!(err_msg.contains("does not exist"));
    }

    #[test]
    fn test_validate_files_exist_all_valid() {
        let temp1 = NamedTempFile::new().expect("creating temp file/dir should succeed");
        let temp2 = NamedTempFile::new().expect("creating temp file/dir should succeed");

        let files =
            vec![(temp1.path().to_path_buf(), "File 1"), (temp2.path().to_path_buf(), "File 2")];

        validate_files_exist(&files).expect("file validation should succeed");
    }

    #[test]
    fn test_validate_files_exist_one_invalid() {
        let temp1 = NamedTempFile::new().expect("creating temp file/dir should succeed");

        let files = vec![
            (temp1.path().to_path_buf(), "File 1"),
            (PathBuf::from("/nonexistent.bam"), "File 2"),
        ];

        let result = validate_files_exist(&files);
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("File 2"));
    }

    #[test]
    fn test_validate_min_max_valid() -> Result<()> {
        // max > min
        validate_min_max(1, Some(10), "min-reads", "max-reads")?;

        // max == min
        validate_min_max(5, Some(5), "min-reads", "max-reads")?;

        // max is None
        validate_min_max(1, None, "min-reads", "max-reads")?;

        Ok(())
    }

    #[test]
    fn test_validate_min_max_invalid() {
        let result = validate_min_max(10, Some(5), "min-reads", "max-reads");
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("max-reads"));
        assert!(err_msg.contains("min-reads"));
        assert!(err_msg.contains(">="));
    }

    #[rstest]
    #[case(0.0, true, "minimum valid rate")]
    #[case(0.01, true, "typical low rate")]
    #[case(0.5, true, "middle rate")]
    #[case(1.0, true, "maximum valid rate")]
    #[case(-0.1, false, "negative rate")]
    #[case(1.5, false, "above maximum")]
    #[case(2.0, false, "far above maximum")]
    fn test_validate_error_rate(
        #[case] rate: f64,
        #[case] should_succeed: bool,
        #[case] description: &str,
    ) {
        let result = validate_error_rate(rate, "error-rate");
        if should_succeed {
            assert!(result.is_ok(), "Failed for: {description}");
        } else {
            assert!(result.is_err(), "Should have failed for: {description}");
            let err_msg = result.unwrap_err().to_string();
            assert!(
                err_msg.contains("Invalid frequency threshold"),
                "Missing expected error for: {description}"
            );
            assert!(err_msg.contains("between 0 and 1"), "Missing range info for: {description}");
        }
    }

    #[rstest]
    #[case(0, true, "minimum valid quality")]
    #[case(30, true, "typical quality")]
    #[case(93, true, "maximum valid quality")]
    #[case(60, true, "high quality")]
    #[case(94, false, "just above maximum")]
    #[case(100, false, "far above maximum")]
    fn test_validate_quality_score(
        #[case] score: u8,
        #[case] should_succeed: bool,
        #[case] description: &str,
    ) {
        let result = validate_quality_score(score, "min-base-quality");
        if should_succeed {
            assert!(result.is_ok(), "Failed for: {description}");
        } else {
            assert!(result.is_err(), "Should have failed for: {description}");
            let err_msg = result.unwrap_err().to_string();
            assert!(
                err_msg.contains("Invalid quality threshold"),
                "Missing expected error for: {description}"
            );
            assert!(err_msg.contains("between 0 and 93"), "Missing range info for: {description}");
        }
    }

    #[test]
    fn test_validate_positive_valid() -> Result<()> {
        validate_positive(1, "min-reads")?;
        validate_positive(100, "min-reads")?;
        validate_positive(1_usize, "threshold")?;
        Ok(())
    }

    #[test]
    fn test_validate_positive_zero() {
        let result = validate_positive(0, "min-reads");
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("Invalid parameter 'min-reads'"));
        assert!(err_msg.contains("Must be positive"));
        assert!(err_msg.contains("got: 0"));
    }

    #[test]
    fn test_validate_positive_negative() {
        let result = validate_positive(-5, "threshold");
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("Invalid parameter 'threshold'"));
        assert!(err_msg.contains("Must be positive"));
        assert!(err_msg.contains("got: -5"));
    }

    // ========== Tests for memory size parsing ==========

    #[test]
    fn test_parse_memory_size_plain_numbers() {
        assert_eq!(
            parse_memory_size("768").expect("parse '768' should succeed"),
            768 * 1024 * 1024
        );
        assert_eq!(parse_memory_size("1").expect("parse '1' should succeed"), 1024 * 1024);
        assert_eq!(
            parse_memory_size("4096").expect("parse '4096' should succeed"),
            4096 * 1024 * 1024
        );
    }

    #[test]
    fn test_parse_memory_size_human_readable() {
        assert_eq!(
            parse_memory_size("2GB").expect("parse '2GB' should succeed"),
            2 * 1000 * 1000 * 1000
        );
        assert_eq!(
            parse_memory_size("2G").expect("parse '2G' should succeed"),
            2 * 1000 * 1000 * 1000
        );
        assert_eq!(
            parse_memory_size("1024MB").expect("parse '1024MB' should succeed"),
            1024 * 1000 * 1000
        );
        assert_eq!(
            parse_memory_size("1024M").expect("parse '1024M' should succeed"),
            1024 * 1000 * 1000
        );
        assert_eq!(
            parse_memory_size("1GiB").expect("parse '1GiB' should succeed"),
            1024 * 1024 * 1024
        );
        assert_eq!(
            parse_memory_size("512MiB").expect("parse '512MiB' should succeed"),
            512 * 1024 * 1024
        );
    }

    #[test]
    fn test_parse_memory_size_invalid() {
        assert!(parse_memory_size("invalid").is_err());
        assert!(parse_memory_size("").is_err());
        assert!(parse_memory_size("GB2").is_err());
    }

    #[test]
    fn test_parse_memory_size_zero() {
        assert!(parse_memory_size("0").is_err());
        assert!(parse_memory_size("0MB").is_err());
        assert!(parse_memory_size("0GB").is_err());
    }

    #[test]
    fn test_parse_memory_size_edge_cases() {
        assert!(parse_memory_size("").is_err());
        assert!(parse_memory_size("   ").is_err());
        assert!(parse_memory_size("-100").is_err());
        assert!(parse_memory_size("-1GB").is_err());
        assert!(parse_memory_size("1.5").is_err());
        assert!(parse_memory_size("1.5GB").is_ok());
        assert!(parse_memory_size("2.5GB").is_ok());
        assert!(parse_memory_size("1e3").is_err());
        assert!(parse_memory_size("1E6").is_err());
        assert!(parse_memory_size("9999999").is_err());
    }

    #[test]
    fn test_parse_memory_size_whitespace_handling() {
        assert_eq!(
            parse_memory_size("  768  ").expect("parse trimmed '768' should succeed"),
            768 * 1024 * 1024
        );
        assert_eq!(
            parse_memory_size("\t1GB\n").expect("parse trimmed '1GB' should succeed"),
            1000 * 1000 * 1000
        );
    }

    #[test]
    fn test_parse_memory_size_overflow() {
        let very_large = format!("{}", u64::MAX / 1024);
        assert!(parse_memory_size(&very_large).is_err());
    }
}
