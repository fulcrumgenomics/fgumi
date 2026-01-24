//! Compare two TSV metrics files with optional float precision rounding.
//!
//! This is useful for comparing metrics files produced by fgbio and fgumi,
//! which may have slightly different floating-point representations.

use anyhow::Result;
use clap::Parser;
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::validation::validate_file_exists;
use log::info;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use crate::commands::command::Command;

/// Compare two TSV metrics files with optional float precision rounding.
///
/// Compares TSV files line by line, field by field. Floats can be rounded to
/// a specified precision before comparison to account for floating-point
/// representation differences between tools.
#[derive(Debug, Parser)]
#[command(
    name = "metrics",
    about = "Compare two TSV metrics files",
    long_about = r#"
Compare two TSV metrics files with optional float precision rounding.

This tool compares TSV files line by line and field by field:
- Integers are compared exactly
- Floats are optionally rounded to a specified precision before comparison
- Float comparison uses both relative and absolute tolerance
- Strings are compared exactly

The tool is designed to compare metrics files from fgbio and fgumi, which
may produce slightly different floating-point representations.

Example usage:
  fgumi compare metrics file1.txt file2.txt
  fgumi compare metrics file1.txt file2.txt --precision 6
  fgumi compare metrics file1.txt file2.txt --rel-tol 1e-6 --abs-tol 1e-9
  fgumi compare metrics file1.txt file2.txt --quiet
"#
)]
pub struct CompareMetrics {
    /// First TSV file
    #[arg(index = 1)]
    pub file1: PathBuf,

    /// Second TSV file
    #[arg(index = 2)]
    pub file2: PathBuf,

    /// Round floats to this many decimal places (-1 to disable rounding)
    #[arg(long = "precision", default_value = "6")]
    pub precision: i32,

    /// Relative tolerance for float comparison
    #[arg(long = "rel-tol", default_value = "1e-9")]
    pub rel_tol: f64,

    /// Absolute tolerance for float comparison
    #[arg(long = "abs-tol", default_value = "1e-9")]
    pub abs_tol: f64,

    /// Maximum number of differences to display
    #[arg(short = 'm', long = "max-diffs", default_value = "20")]
    pub max_diffs: usize,

    /// Quiet mode - only exit code indicates result (0=equal, 1=different)
    #[arg(short = 'q', long = "quiet")]
    pub quiet: bool,

    /// Verbose mode - print success message when files match
    #[arg(short = 'v', long = "verbose")]
    pub verbose: bool,
}

/// Parsed value from a TSV field.
#[derive(Debug, Clone)]
enum ParsedValue {
    Int(i64),
    Float(f64),
    String(String),
}

impl std::fmt::Display for ParsedValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParsedValue::Int(i) => write!(f, "{i}"),
            ParsedValue::Float(fl) => write!(f, "{fl}"),
            ParsedValue::String(s) => write!(f, "{s}"),
        }
    }
}

/// Parse a string value as int, float (optionally rounded), or keep as string.
fn parse_value(value: &str, precision: Option<u32>) -> ParsedValue {
    // Try to parse as integer first
    if let Ok(i) = value.parse::<i64>() {
        return ParsedValue::Int(i);
    }

    // Try to parse as float
    if let Ok(f) = value.parse::<f64>() {
        let f = if let Some(p) = precision {
            let multiplier = 10_f64.powi(p as i32);
            (f * multiplier).round() / multiplier
        } else {
            f
        };
        return ParsedValue::Float(f);
    }

    // Keep as string
    ParsedValue::String(value.to_string())
}

/// Compare two values, using tolerance for floats.
fn values_equal(val1: &ParsedValue, val2: &ParsedValue, rel_tol: f64, abs_tol: f64) -> bool {
    match (val1, val2) {
        (ParsedValue::Int(i1), ParsedValue::Int(i2)) => i1 == i2,
        (ParsedValue::Float(f1), ParsedValue::Float(f2)) => {
            // Treat NaN == NaN (unlike IEEE 754 semantics where NaN != NaN)
            if f1.is_nan() && f2.is_nan() {
                return true;
            }
            // Handle infinities: +inf == +inf, -inf == -inf, but +inf != -inf
            if f1.is_infinite() && f2.is_infinite() {
                return f1.signum() == f2.signum();
            }
            // If only one is infinite, they are not equal
            if f1.is_infinite() || f2.is_infinite() {
                return false;
            }
            let diff = (f1 - f2).abs();
            let max_val = f1.abs().max(f2.abs());
            diff <= (rel_tol * max_val).max(abs_tol)
        }
        // Try to compare as floats if one is int and one is float
        (ParsedValue::Int(i), ParsedValue::Float(f))
        | (ParsedValue::Float(f), ParsedValue::Int(i)) => {
            let f1 = *i as f64;
            let f2 = *f;
            let diff = (f1 - f2).abs();
            let max_val = f1.abs().max(f2.abs());
            diff <= (rel_tol * max_val).max(abs_tol)
        }
        (ParsedValue::String(s1), ParsedValue::String(s2)) => s1 == s2,
        _ => false,
    }
}

impl Command for CompareMetrics {
    fn execute(&self, _command_line: &str) -> Result<()> {
        validate_file_exists(&self.file1, "First file")?;
        validate_file_exists(&self.file2, "Second file")?;

        let timer = OperationTimer::new("Comparing metrics");

        let precision = if self.precision >= 0 { Some(self.precision as u32) } else { None };

        // Open files as streaming iterators
        let file1 = File::open(&self.file1)?;
        let file2 = File::open(&self.file2)?;
        let reader1 = BufReader::new(file1);
        let reader2 = BufReader::new(file2);

        let mut lines1 = reader1.lines();
        let mut lines2 = reader2.lines();

        let mut differences: Vec<String> = Vec::new();
        let mut line_num = 0usize;
        let progress = ProgressTracker::new("Compared lines").with_interval(100_000);

        // Stream through both files simultaneously
        loop {
            let line1 = lines1.next();
            let line2 = lines2.next();

            match (line1, line2) {
                (Some(Ok(l1)), Some(Ok(l2))) => {
                    line_num += 1;
                    progress.log_if_needed(1);

                    let fields1: Vec<&str> =
                        l1.trim_end_matches(['\n', '\r']).split('\t').collect();
                    let fields2: Vec<&str> =
                        l2.trim_end_matches(['\n', '\r']).split('\t').collect();

                    if fields1.len() != fields2.len() {
                        differences.push(format!(
                            "Line {line_num}: field count mismatch ({} vs {})",
                            fields1.len(),
                            fields2.len()
                        ));
                        continue;
                    }

                    for (col_num, (f1, f2)) in fields1.iter().zip(fields2.iter()).enumerate() {
                        let col_num = col_num + 1; // 1-based
                        let val1 = parse_value(f1, precision);
                        let val2 = parse_value(f2, precision);

                        if !values_equal(&val1, &val2, self.rel_tol, self.abs_tol) {
                            differences.push(format!(
                                "Line {line_num}, col {col_num}: {f1:?} != {f2:?} (parsed: {val1} != {val2})"
                            ));
                        }
                    }
                }
                (None, None) => {
                    // Both files ended at the same time - success
                    break;
                }
                (Some(_), None) => {
                    // File1 has more lines
                    // Count remaining lines in file1
                    let extra_lines = 1 + lines1.count();
                    if !self.quiet {
                        println!(
                            "Line count mismatch: {} has {} extra lines after line {}",
                            self.file1.display(),
                            extra_lines,
                            line_num
                        );
                    }
                    std::process::exit(1);
                }
                (None, Some(_)) => {
                    // File2 has more lines
                    let extra_lines = 1 + lines2.count();
                    if !self.quiet {
                        println!(
                            "Line count mismatch: {} has {} extra lines after line {}",
                            self.file2.display(),
                            extra_lines,
                            line_num
                        );
                    }
                    std::process::exit(1);
                }
                (Some(Err(e)), _) | (_, Some(Err(e))) => {
                    return Err(anyhow::anyhow!(
                        "Error reading file at line {}: {}",
                        line_num + 1,
                        e
                    ));
                }
            }
        }

        progress.log_final();
        let is_equal = differences.is_empty();

        if !self.quiet {
            if !differences.is_empty() {
                println!("Found {} difference(s) between files:", differences.len());
                for diff in differences.iter().take(self.max_diffs) {
                    println!("  {diff}");
                }
                if differences.len() > self.max_diffs {
                    println!("  ... and {} more", differences.len() - self.max_diffs);
                }
            } else if self.verbose {
                println!("Files match ({line_num} lines compared)");
            }
        }

        if is_equal {
            info!("Metrics files are identical");
            timer.log_completion(line_num as u64);
            Ok(())
        } else {
            info!("Metrics files differ");
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nan_equality() {
        // NaN values should be considered equal
        let nan1 = ParsedValue::Float(f64::NAN);
        let nan2 = ParsedValue::Float(f64::NAN);
        assert!(values_equal(&nan1, &nan2, 1e-9, 1e-9));
    }

    #[test]
    fn test_positive_infinity_equality() {
        // +Infinity values should be considered equal
        let inf1 = ParsedValue::Float(f64::INFINITY);
        let inf2 = ParsedValue::Float(f64::INFINITY);
        assert!(values_equal(&inf1, &inf2, 1e-9, 1e-9));
    }

    #[test]
    fn test_negative_infinity_equality() {
        // -Infinity values should be considered equal
        let neg_inf1 = ParsedValue::Float(f64::NEG_INFINITY);
        let neg_inf2 = ParsedValue::Float(f64::NEG_INFINITY);
        assert!(values_equal(&neg_inf1, &neg_inf2, 1e-9, 1e-9));
    }

    #[test]
    fn test_mixed_infinity_not_equal() {
        // +Infinity and -Infinity should NOT be equal
        let pos_inf = ParsedValue::Float(f64::INFINITY);
        let neg_inf = ParsedValue::Float(f64::NEG_INFINITY);
        assert!(!values_equal(&pos_inf, &neg_inf, 1e-9, 1e-9));
    }

    #[test]
    fn test_nan_not_equal_to_number() {
        // NaN and a regular number should NOT be equal
        let nan = ParsedValue::Float(f64::NAN);
        let num = ParsedValue::Float(1.0);
        assert!(!values_equal(&nan, &num, 1e-9, 1e-9));
    }

    #[test]
    fn test_parse_nan_string() {
        // "NaN" string should parse to NaN float
        let parsed = parse_value("NaN", Some(6));
        match parsed {
            ParsedValue::Float(f) => assert!(f.is_nan()),
            _ => panic!("Expected Float variant for NaN"),
        }
    }

    #[test]
    fn test_parse_infinity_string() {
        // "Infinity" and "inf" should parse to infinity
        let parsed = parse_value("Infinity", Some(6));
        match parsed {
            ParsedValue::Float(f) => assert!(f.is_infinite() && f.is_sign_positive()),
            _ => panic!("Expected Float variant for Infinity"),
        }

        let parsed = parse_value("inf", Some(6));
        match parsed {
            ParsedValue::Float(f) => assert!(f.is_infinite() && f.is_sign_positive()),
            _ => panic!("Expected Float variant for inf"),
        }
    }
}
