//! Input validation utilities
//!
//! This module provides common validation functions for command-line parameters,
//! file paths, and SAM tags with consistent error messages.
//!
//! All validation functions now use structured error types from [`crate::errors`] to provide
//! rich contextual information when validation fails.

use crate::errors::{FgumiError, Result};
use noodles::sam::alignment::record::data::field::Tag;
use std::fmt::Display;
use std::path::Path;

/// Validate that a file exists
///
/// # Arguments
/// * `path` - Path to validate
/// * `description` - Human-readable description of the file (e.g., "Input file", "Reference")
///
/// # Errors
/// Returns an error if the file does not exist
///
/// # Example
/// ```
/// use fgumi_lib::validation::validate_file_exists;
/// use std::path::Path;
///
/// let result = validate_file_exists("/nonexistent/file.bam", "Input file");
/// assert!(result.is_err());
/// ```
pub fn validate_file_exists<P: AsRef<Path>>(path: P, description: &str) -> Result<()> {
    let path_ref = path.as_ref();
    if !path_ref.exists() {
        return Err(FgumiError::InvalidFileFormat {
            file_type: description.to_string(),
            path: path_ref.display().to_string(),
            reason: "File does not exist".to_string(),
        });
    }
    Ok(())
}

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

/// Validate that a SAM tag is exactly 2 characters
///
/// # Arguments
/// * `tag` - Tag string to validate
/// * `name` - Name of the parameter for error messages
///
/// # Returns
/// A 2-byte array representing the tag
///
/// # Errors
/// Returns an error if the tag is not exactly 2 characters
///
/// # Example
/// ```
/// use fgumi_lib::validation::validate_tag;
///
/// let tag = validate_tag("MI", "UMI tag").unwrap();
/// assert_eq!(tag, [b'M', b'I']);
///
/// let result = validate_tag("ABC", "UMI tag");
/// assert!(result.is_err());
/// ```
pub fn validate_tag(tag: &str, name: &str) -> Result<[u8; 2]> {
    if tag.len() != 2 {
        return Err(FgumiError::InvalidParameter {
            parameter: name.to_string(),
            reason: format!("Tag must be exactly 2 characters, got: '{tag}'"),
        });
    }
    let bytes = tag.as_bytes();
    Ok([bytes[0], bytes[1]])
}

/// Convert a validated string tag to noodles Tag type
///
/// This combines validation and conversion for convenience.
///
/// # Arguments
/// * `tag` - Tag string to validate and convert
/// * `name` - Name of the parameter for error messages
///
/// # Returns
/// A noodles `Tag` object
///
/// # Errors
/// Returns an error if the tag is not exactly 2 characters
///
/// # Example
/// ```
/// use fgumi_lib::validation::string_to_tag;
///
/// let tag = string_to_tag("MI", "UMI tag").unwrap();
/// ```
pub fn string_to_tag(tag: &str, name: &str) -> Result<Tag> {
    let tag_array = validate_tag(tag, name)?;
    Ok(Tag::from(tag_array))
}

/// Convert an optional string tag to an optional noodles Tag type
///
/// This is a convenience function for optional CLI arguments like cell barcode tags.
///
/// # Arguments
/// * `tag` - Optional tag string to validate and convert
/// * `name` - Name of the parameter for error messages
///
/// # Returns
/// `Some(Tag)` if the input is `Some` and valid, `None` if the input is `None`
///
/// # Errors
/// Returns an error if the tag is `Some` but not exactly 2 characters
///
/// # Example
/// ```
/// use fgumi_lib::validation::optional_string_to_tag;
///
/// // Some tag - validates and converts
/// let tag = optional_string_to_tag(Some("CB"), "cell tag").unwrap();
/// assert!(tag.is_some());
///
/// // None - returns None
/// let tag = optional_string_to_tag(None, "cell tag").unwrap();
/// assert!(tag.is_none());
///
/// // Invalid tag - returns error
/// let result = optional_string_to_tag(Some("ABC"), "cell tag");
/// assert!(result.is_err());
/// ```
pub fn optional_string_to_tag(tag: Option<&str>, name: &str) -> Result<Option<Tag>> {
    tag.map(|t| string_to_tag(t, name)).transpose()
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
/// * `name` - Name of the parameter for error messages
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
/// * `name` - Name of the parameter for error messages
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
        let temp_file = NamedTempFile::new().unwrap();
        validate_file_exists(temp_file.path(), "Test file").unwrap();
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
        let temp1 = NamedTempFile::new().unwrap();
        let temp2 = NamedTempFile::new().unwrap();

        let files =
            vec![(temp1.path().to_path_buf(), "File 1"), (temp2.path().to_path_buf(), "File 2")];

        validate_files_exist(&files).unwrap();
    }

    #[test]
    fn test_validate_files_exist_one_invalid() {
        let temp1 = NamedTempFile::new().unwrap();

        let files = vec![
            (temp1.path().to_path_buf(), "File 1"),
            (PathBuf::from("/nonexistent.bam"), "File 2"),
        ];

        let result = validate_files_exist(&files);
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("File 2"));
    }

    #[rstest]
    #[case("MI", true, Some([b'M', b'I']), "valid MI tag")]
    #[case("RX", true, Some([b'R', b'X']), "valid RX tag")]
    #[case("AB", true, Some([b'A', b'B']), "valid AB tag")]
    #[case("M", false, None, "too short")]
    #[case("ABC", false, None, "too long")]
    #[case("", false, None, "empty string")]
    fn test_validate_tag(
        #[case] input: &str,
        #[case] should_succeed: bool,
        #[case] expected: Option<[u8; 2]>,
        #[case] description: &str,
    ) {
        let result = validate_tag(input, "test tag");
        if should_succeed {
            assert!(result.is_ok(), "Failed for: {description}");
            assert_eq!(result.unwrap(), expected.unwrap(), "Failed for: {description}");
        } else {
            assert!(result.is_err(), "Should have failed for: {description}");
            let err_msg = result.unwrap_err().to_string();
            assert!(
                err_msg.contains("must be exactly 2 characters"),
                "Missing expected error message for: {description}"
            );
        }
    }

    #[test]
    fn test_string_to_tag_valid() -> Result<()> {
        let tag = string_to_tag("MI", "UMI tag")?;
        assert_eq!(tag, Tag::from([b'M', b'I']));
        Ok(())
    }

    #[test]
    fn test_string_to_tag_invalid_length() {
        let result = string_to_tag("ABC", "UMI tag");
        assert!(result.is_err());
    }

    #[test]
    fn test_optional_string_to_tag_some_valid() -> Result<()> {
        let tag = optional_string_to_tag(Some("CB"), "cell tag")?;
        assert!(tag.is_some());
        assert_eq!(tag.unwrap(), Tag::from([b'C', b'B']));
        Ok(())
    }

    #[test]
    fn test_optional_string_to_tag_none() -> Result<()> {
        let tag = optional_string_to_tag(None, "cell tag")?;
        assert!(tag.is_none());
        Ok(())
    }

    #[test]
    fn test_optional_string_to_tag_some_invalid() {
        let result = optional_string_to_tag(Some("ABC"), "cell tag");
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("must be exactly 2 characters"));
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
}
