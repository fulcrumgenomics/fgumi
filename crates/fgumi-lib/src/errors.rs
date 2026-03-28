//! Custom error types for fgumi operations.

use thiserror::Error;

/// Result type alias for fgumi operations
pub type Result<T> = std::result::Result<T, FgumiError>;

/// Error type for fgumi operations
#[derive(Error, Debug)]
pub enum FgumiError {
    /// Invalid parameter value provided
    #[error("Invalid parameter '{parameter}': {reason}")]
    InvalidParameter {
        /// The parameter name
        parameter: String,
        /// Explanation of why it's invalid
        reason: String,
    },

    /// Invalid frequency threshold
    #[error("Invalid frequency threshold: {value} (must be between {min} and {max})")]
    InvalidFrequency {
        /// The invalid frequency value
        value: f64,
        /// Minimum valid value
        min: f64,
        /// Maximum valid value
        max: f64,
    },

    /// Invalid quality threshold
    #[error("Invalid quality threshold: {value} (must be between 0 and {max})")]
    InvalidQuality {
        /// The invalid quality value
        value: u8,
        /// Maximum valid value (usually 93 for SAM/BAM)
        max: u8,
    },

    /// File format error
    #[error("Invalid {file_type} file '{path}': {reason}")]
    InvalidFileFormat {
        /// Type of file (e.g., "BAM", "FASTQ")
        file_type: String,
        /// Path to the file
        path: String,
        /// Explanation of the problem
        reason: String,
    },

    /// Required reference sequence not found
    #[error("Reference sequence '{ref_name}' not found in header")]
    ReferenceNotFound {
        /// The reference sequence name
        ref_name: String,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_invalid_parameter() {
        let error = FgumiError::InvalidParameter {
            parameter: "min-reads".to_string(),
            reason: "must be >= 1".to_string(),
        };
        let msg = format!("{error}");
        assert!(msg.contains("Invalid parameter 'min-reads'"));
        assert!(msg.contains("must be >= 1"));
    }

    #[test]
    fn test_invalid_frequency() {
        let error = FgumiError::InvalidFrequency { value: 1.5, min: 0.0, max: 1.0 };
        let msg = format!("{error}");
        assert!(msg.contains("1.5"));
        assert!(msg.contains("between 0 and 1"));
    }

    #[test]
    fn test_invalid_file_format() {
        let error = FgumiError::InvalidFileFormat {
            file_type: "BAM".to_string(),
            path: "/path/to/file.bam".to_string(),
            reason: "truncated file".to_string(),
        };
        let msg = format!("{error}");
        assert!(msg.contains("Invalid BAM file"));
        assert!(msg.contains("truncated file"));
    }

    #[test]
    fn test_reference_not_found() {
        let error = FgumiError::ReferenceNotFound { ref_name: "chr1".to_string() };
        let msg = format!("{error}");
        assert!(msg.contains("Reference sequence 'chr1' not found"));
    }
}
