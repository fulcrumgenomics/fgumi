//! Input validation utilities
//!
//! This module provides common validation functions for command-line parameters,
//! file paths, and SAM tags with consistent error messages.
//!
//! All validation functions now use structured error types from [`crate::errors`] to provide
//! rich contextual information when validation fails.

use crate::errors::{FgumiError, Result};
use bytesize::ByteSize;
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
    if let Some(max) = max_val
        && max < min_val
    {
        return Err(FgumiError::InvalidParameter {
            parameter: max_name.to_string(),
            reason: format!("{max_name} ({max}) must be >= {min_name} ({min_val})"),
        });
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

/// Does `s` contain a genuine scientific-notation exponent?
///
/// True only for a digit or `.` followed by `e`/`E` and then an optional sign
/// and at least one digit -- `1e5`, `1.5e10`, `2E+3`. A bare `e`/`E` test would
/// also match the exabyte units `EB` and `EiB`, which are the only units
/// `bytesize` accepts that contain an `E`, so `1EB` would be rejected as
/// scientific notation.
fn looks_like_scientific_notation(s: &str) -> bool {
    let bytes = s.as_bytes();
    bytes.iter().enumerate().skip(1).any(|(i, &c)| {
        if c != b'e' && c != b'E' {
            return false;
        }
        // Mantissa: the character before must end a number.
        let prev = bytes[i - 1];
        if !prev.is_ascii_digit() && prev != b'.' {
            return false;
        }
        // Exponent: an optional sign then at least one digit.
        let mut rest = bytes[i + 1..].iter();
        match rest.next() {
            Some(&b'+' | &b'-') => rest.next().is_some_and(u8::is_ascii_digit),
            Some(c) => c.is_ascii_digit(),
            None => false,
        }
    })
}

/// Parses a memory size string into bytes.
///
/// Accepts both plain numbers (interpreted as MiB) and human-readable formats like:
/// - "2GB", "2G" -> 2 gigabytes (decimal: 2,000,000,000)
/// - "1.5GB" -> 1.5 gigabytes
/// - "1024MB", "1024M" -> 1024 megabytes (decimal)
/// - "512MiB" -> 512 mebibytes (binary: 536,870,912)
/// - "768" -> 768 MiB (plain numbers are interpreted as mebibytes, i.e. `n * 1024 * 1024`)
///
/// # Examples
///
/// ```
/// # use fgumi_lib::validation::parse_memory_size;
/// assert_eq!(parse_memory_size("768").unwrap(), 768 * 1024 * 1024);
/// assert_eq!(parse_memory_size("2GB").unwrap(), 2 * 1000 * 1000 * 1000);
/// assert_eq!(parse_memory_size("1024MiB").unwrap(), 1024 * 1024 * 1024);
/// ```
///
/// # Errors
///
/// Returns [`FgumiError::InvalidMemorySize`] if the string cannot be parsed as a valid size.
pub fn parse_memory_size(size_str: &str) -> Result<u64> {
    let trimmed = size_str.trim();
    if trimmed.is_empty() {
        return Err(FgumiError::InvalidMemorySize {
            reason: "Memory size cannot be empty".to_string(),
        });
    }

    // Handle negative values early
    if trimmed.starts_with('-') {
        return Err(FgumiError::InvalidMemorySize {
            reason: format!("Memory size cannot be negative: '{trimmed}'"),
        });
    }

    // First try parsing as a plain integer in MiB (backward compatibility)
    // Only accept simple integers, not floats or scientific notation
    if let Ok(mb_value) = trimmed.parse::<u64>() {
        if mb_value == 0 {
            return Err(FgumiError::InvalidMemorySize {
                reason: "Memory size cannot be zero".to_string(),
            });
        }
        if mb_value > 1_000_000 {
            // Sanity guard: >1TB as a plain number likely means the user forgot a unit suffix.
            return Err(FgumiError::InvalidMemorySize {
                reason: format!(
                    "Plain number memory size too large: {} MiB. Use human-readable format like '{}GB' instead.",
                    mb_value,
                    mb_value / 1000
                ),
            });
        }

        return mb_value.checked_mul(1024 * 1024).ok_or_else(|| FgumiError::InvalidMemorySize {
            reason: format!("Memory size calculation overflow for {mb_value} MiB"),
        });
    }

    // Reject scientific notation (e.g. "1e3") but allow decimals in human-readable
    // sizes (e.g. "1.5GB"). Testing for a bare 'e'/'E' would also reject the
    // exabyte units "EB"/"EiB" -- the only supported units containing 'E' -- so
    // match a genuine mantissa-and-exponent instead.
    if looks_like_scientific_notation(trimmed) {
        return Err(FgumiError::InvalidMemorySize {
            reason: format!(
                "Scientific notation not supported: '{trimmed}'. Use integer values or human-readable formats like '2GB'."
            ),
        });
    }

    // Reject bare decimal numbers without a unit suffix (e.g. "1.5") since plain numbers are MiB
    if trimmed.contains('.') && trimmed.chars().all(|c| c.is_ascii_digit() || c == '.') {
        return Err(FgumiError::InvalidMemorySize {
            reason: format!(
                "Plain decimal numbers not supported: '{trimmed}'. Use an integer for MiB (e.g. '768') or a human-readable format (e.g. '1.5GB')."
            ),
        });
    }

    // Fall back to parsing as a human-readable size (like "2GB", "1024MiB")
    match trimmed.parse::<ByteSize>() {
        Ok(size) => {
            if size.0 == 0 {
                return Err(FgumiError::InvalidMemorySize {
                    reason: format!("Memory size cannot be zero: '{trimmed}'"),
                });
            }
            Ok(size.0)
        }
        Err(_) => Err(FgumiError::InvalidMemorySize {
            reason: format!(
                "Invalid memory size '{trimmed}'. Valid formats:\n\
                 - Plain numbers (interpreted as MiB): '768', '4096'\n\
                 - Human-readable (decimal): '2GB', '1024MB'\n\
                 - Human-readable (binary): '1GiB', '512MiB'"
            ),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `EB`/`EiB` are the only units `bytesize` accepts that contain an `E`, so
    /// a bare `contains('e')` check rejected them as scientific notation. The
    /// replacement must still reject real exponents.
    #[rstest]
    #[case::exabyte("1EB", true)]
    #[case::exbibyte("1EiB", true)]
    #[case::exabyte_lower("1eb", true)]
    #[case::gigabyte_unaffected("2GB", true)]
    #[case::decimal_with_unit("1.5GB", true)]
    #[case::plain_integer_mib("768", true)]
    #[case::simple_exponent("1e5", false)]
    #[case::decimal_exponent("1.5e10", false)]
    #[case::signed_exponent("2E+3", false)]
    #[case::negative_exponent("2e-3", false)]
    fn test_parse_memory_size_exabytes_are_not_scientific_notation(
        #[case] input: &str,
        #[case] expect_ok: bool,
    ) {
        let result = parse_memory_size(input);
        assert_eq!(result.is_ok(), expect_ok, "{input}: expected ok={expect_ok}, got {result:?}",);
        if !expect_ok {
            let msg = result.unwrap_err().to_string();
            assert!(
                msg.contains("Scientific notation"),
                "{input} should be rejected as scientific notation, got: {msg}",
            );
        }
    }

    /// The exponent detector on its own, including the shapes that are *not*
    /// exponents and must be left for `bytesize` to parse.
    #[rstest]
    #[case::plain_exponent("1e5", true)]
    #[case::uppercase_exponent("1E5", true)]
    #[case::after_decimal_point("1.e5", true)]
    #[case::signed("3e+2", true)]
    #[case::exabyte_unit("1EB", false)]
    #[case::exbibyte_unit("1EiB", false)]
    #[case::leading_e("e5", false)]
    #[case::trailing_e("5e", false)]
    #[case::sign_without_digit("5e+", false)]
    #[case::letter_before_e("Xe5", false)]
    #[case::no_e_at_all("1024MB", false)]
    fn test_looks_like_scientific_notation(#[case] input: &str, #[case] expected: bool) {
        assert_eq!(looks_like_scientific_notation(input), expected, "input: {input}");
    }
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
