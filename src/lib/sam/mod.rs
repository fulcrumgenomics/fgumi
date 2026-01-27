//! SAM/BAM file utilities and header manipulation.
//!
//! This module provides utilities for working with SAM/BAM files, including:
//! - Checking and validating SAM header sort orders
//! - Reversing and reverse-complementing per-base tag values
//! - Template-coordinate sorting validation
//! - Test utilities for building SAM/BAM records
//! - Alignment tag regeneration (NM, MD, UQ)
//! - Record-level utilities (position mapping, FR pair detection, CIGAR parsing)
//!
//! # Sort Orders
//!
//! The module supports several important sort orders:
//! - **queryname** - Reads sorted by query name (required for grouping by UMI)
//! - **template-coordinate** - Special sort order from fgbio where reads are grouped
//!   by query name but ordered by genomic position within each template
//!
//! # Tag Manipulation
//!
//! Functions are provided to reverse or reverse-complement per-base tag values when
//! reads are mapped to the negative strand, ensuring tag values match the orientation
//! of the read sequence.
//!
//! # Record Utilities
//!
//! The [`record_utils`] submodule provides utilities for working with individual records:
//! - [`record_utils::read_pos_at_ref_pos`] - Map reference position to read position
//! - [`record_utils::is_fr_pair_from_tags`] - Check if read is part of FR pair using tags
//! - [`record_utils::mate_unclipped_start`] / [`record_utils::mate_unclipped_end`] - Get mate boundaries from MC tag
//! - [`record_utils::num_bases_extending_past_mate`] - Calculate overlap with mate
//! - [`record_utils::parse_cigar_string`] - Parse CIGAR string to operations

pub mod alignment_tags;
pub mod builder;
pub mod record_utils;

// Re-export commonly used items from submodules for convenience
pub use alignment_tags::regenerate_alignment_tags;
pub use builder::{
    ConsensusTagsBuilder, FragBuilder, MAPPED_PG_ID, PairBuilder, REFERENCE_LENGTH, RecordBuilder,
    SamBuilder, Strand, create_default_test_fasta, create_ref_dict, create_test_fasta,
    degrading_qualities, parse_cigar, repeat_n, uniform_qualities,
};
pub use record_utils::{
    PairOrientation, alignment_end, cigar_reference_length, get_pair_orientation, is_fr_pair,
    is_fr_pair_from_tags, leading_clipping, leading_soft_clipping, mate_unclipped_end,
    mate_unclipped_start, num_bases_extending_past_mate, parse_cigar_string, read_pos_at_ref_pos,
    reference_length, trailing_clipping, trailing_soft_clipping, unclipped_end,
    unclipped_five_prime_position, unclipped_start,
};

use bstr::ByteSlice;
use log::warn;
use noodles::sam::Header;

use crate::dna::complement_base;
use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
use noodles::sam::header::record::value::map::header::sort_order::{QUERY_NAME, UNSORTED};
use std::path::Path;

/// Checks if a BAM file has a specified sort order according to its header.
///
/// Examines the SAM header for the SO (sort order) tag and compares it to the
/// specified sort order. This is useful for validating that input files are
/// properly sorted before processing.
///
/// # Arguments
///
/// * `header` - SAM header to check
/// * `sort_order` - The expected sort order (e.g., `QUERY_NAME`, `COORDINATE`)
///
/// # Returns
///
/// `true` if the header's SO tag matches the specified sort order, `false` otherwise
///
/// # Examples
///
/// ```rust,ignore
/// use noodles::sam::header::record::value::map::header::sort_order::QUERY_NAME;
///
/// if is_sorted(&header, QUERY_NAME) {
///     // Process queryname-sorted data
/// }
/// ```
#[must_use]
pub fn is_sorted(header: &Header, sort_order: &[u8]) -> bool {
    if let Some(hdr_map) = header.header() {
        hdr_map
            .other_fields()
            .get(b"SO")
            .is_some_and(|so| <_ as AsRef<[u8]>>::as_ref(so) == sort_order)
    } else {
        false
    }
}

/// Checks if a BAM file is template-coordinate sorted according to its header.
///
/// Template-coordinate sorting is a special sort order used by fgbio where reads are:
/// 1. Grouped by query name (all reads with same name together)
/// 2. Ordered by the genomic position of the template's lower coordinate
///
/// This sort order is indicated in the SAM header by:
/// - SO:unsorted (not coordinate sorted)
/// - GO:query (grouped by query name)
/// - SS:template-coordinate (optional, but must match if present)
///
/// This matches the behavior of fgbio's `GroupReadsByUmi` output format.
///
/// # Arguments
///
/// * `header` - SAM header to check
///
/// # Returns
///
/// `true` if the header indicates template-coordinate sorting, `false` otherwise
///
/// # Examples
///
/// ```rust,ignore
/// if is_template_coordinate_sorted(&header) {
///     // Process template-coordinate sorted data efficiently
/// }
/// ```
#[must_use]
pub fn is_template_coordinate_sorted(header: &Header) -> bool {
    if let Some(hdr_map) = header.header() {
        let other_fields = hdr_map.other_fields();

        // Check SO tag - must be "unsorted"
        let is_unsorted =
            other_fields.get(b"SO").is_some_and(|so| <_ as AsRef<[u8]>>::as_ref(so) == UNSORTED);

        // Check GO tag - must be "query"
        let is_query_grouped =
            other_fields.get(b"GO").is_some_and(|go| <_ as AsRef<[u8]>>::as_ref(go) == b"query");

        // Check SS tag - if present, must be "template-coordinate" (or "SO:template-coordinate"), but it's optional
        // The SS tag may be prefixed with the sort order (e.g., "unsorted:template-coordinate")
        // per the SAM spec, so we need to extract the part after the colon
        let ss_matches = other_fields.get(b"SS").is_none_or(|ss| {
            let ss_bytes = <_ as AsRef<[u8]>>::as_ref(ss);
            // Find the last colon and take everything after it, or use the whole value if no colon
            if let Some(colon_pos) = ss_bytes.iter().position(|&b| b == b':') {
                &ss_bytes[colon_pos + 1..] == b"template-coordinate"
            } else {
                ss_bytes == b"template-coordinate"
            }
        }); // If SS is missing, that's acceptable

        is_unsorted && is_query_grouped && ss_matches
    } else {
        false
    }
}

/// Checks if a BAM file is queryname sorted and logs a warning if not.
///
/// This function validates that the input file has the correct sort order for
/// processing (typically queryname). If the sort order is incorrect, it logs
/// a warning but does not fail, allowing processing to continue with potentially
/// incorrect results.
///
/// Use this for non-critical validation where you want to warn users but not
/// prevent them from proceeding.
///
/// # Arguments
///
/// * `header` - SAM header to check
/// * `path` - Path to the BAM file (used in warning messages)
/// * `name` - Descriptive name of the file (e.g., "unmapped", "input")
///
/// # Examples
///
/// ```rust,ignore
/// check_sort(&header, input_path, "input");
/// // Will log: "input file 'foo.bam' does not appear to be queryname sorted..."
/// ```
pub fn check_sort(header: &Header, path: &Path, name: &str) {
    if !is_sorted(header, QUERY_NAME) {
        warn!(
            "{name} file {} does not appear to be queryname sorted per the SAM header.",
            path.display()
        );
        warn!("Continuing, but your output may be incorrect.");
    }
}

/// Reverses a `BufValue` (array or string).
///
/// This function reverses per-base tag values to match read orientation when reads
/// are mapped to the negative strand. Arrays and strings are reversed element-wise.
///
/// Supported types:
/// - Arrays of any numeric type (i8, u8, i16, u16, i32, u32, f32)
/// - Strings (character order reversed)
///
/// Other value types (integers, floats, characters) are returned unchanged.
///
/// # Arguments
///
/// * `value` - The value to reverse
///
/// # Returns
///
/// A new `BufValue` with reversed contents, or a clone if not reversible
///
/// # Examples
///
/// ```rust,ignore
/// // Reverse per-base quality scores for negative strand read
/// let quals = BufValue::Array(UInt8(vec![30, 25, 20, 15]));
/// let reversed = reverse_buf_value(&quals);
/// // Result: Array(UInt8(vec![15, 20, 25, 30]))
/// ```
#[must_use]
pub fn reverse_buf_value(value: &BufValue) -> BufValue {
    use noodles::sam::alignment::record_buf::data::field::value::Array::{
        Float, Int8, Int16, Int32, UInt8, UInt16, UInt32,
    };
    match value {
        BufValue::Array(arr) => {
            use noodles::sam::alignment::record_buf::data::field::value::Array;
            let new_arr: Array = match arr {
                Int8(vec) => {
                    let mut values = vec.clone();
                    values.reverse();
                    Int8(values)
                }
                UInt8(vec) => {
                    let mut values = vec.clone();
                    values.reverse();
                    UInt8(values)
                }
                Int16(vec) => {
                    let mut values = vec.clone();
                    values.reverse();
                    Int16(values)
                }
                UInt16(vec) => {
                    let mut values = vec.clone();
                    values.reverse();
                    UInt16(values)
                }
                Int32(vec) => {
                    let mut values = vec.clone();
                    values.reverse();
                    Int32(values)
                }
                UInt32(vec) => {
                    let mut values = vec.clone();
                    values.reverse();
                    UInt32(values)
                }
                Float(vec) => {
                    let mut values = vec.clone();
                    values.reverse();
                    Float(values)
                }
            };
            BufValue::Array(new_arr)
        }
        BufValue::String(s) => {
            let mut bytes = s.as_bytes().to_vec();
            bytes.reverse();
            BufValue::from(String::from_utf8_lossy(&bytes).to_string())
        }
        _ => value.clone(),
    }
}

/// Reverse complements a DNA sequence `BufValue`.
///
/// This function performs a full reverse complement operation on DNA sequences stored
/// in per-base tags. It reverses the order of bases and complements each base:
/// - A <-> T
/// - C <-> G
/// - Preserves case (a <-> t, c <-> g)
/// - Other characters are left unchanged
///
/// This is essential for per-base sequence tags (like consensus bases) when reads
/// are mapped to the negative strand, ensuring the tag sequences match the read orientation.
///
/// # Arguments
///
/// * `value` - The value to reverse complement (must be a String)
///
/// # Returns
///
/// A new `BufValue` with reverse complemented sequence, or a clone if not a string
///
/// # Examples
///
/// ```rust,ignore
/// // Reverse complement per-base consensus sequences for negative strand
/// let bases = BufValue::String("ACGT".to_string());
/// let revcomp = revcomp_buf_value(&bases);
/// // Result: String("ACGT") [reverse of TGCA]
/// ```
#[must_use]
pub fn revcomp_buf_value(value: &BufValue) -> BufValue {
    match value {
        BufValue::String(s) => {
            let revcomp: Vec<u8> = s.as_bytes().iter().rev().map(|&b| complement_base(b)).collect();
            BufValue::from(String::from_utf8_lossy(&revcomp).to_string())
        }
        _ => value.clone(),
    }
}

/// Converts an integer value to the smallest signed integer `BufValue` that fits.
///
/// This is used to ensure integer tags are written with the same type encoding as fgbio.
/// fgbio uses signed integer types for tags like MQ, AS, XS, and ms, while noodles
/// may choose unsigned types when writing positive values.
///
/// The function chooses the smallest signed type that fits:
/// - Int8 for values in [-128, 127]
/// - Int16 for values in [-32768, 32767]
/// - Int32 for larger values
///
/// # Arguments
///
/// * `value` - The integer value to convert
///
/// # Returns
///
/// A `BufValue` using the smallest signed type that fits the value
///
/// # Examples
///
/// ```rust,ignore
/// let mq_value = to_smallest_signed_int(60);
/// // Result: BufValue::Int8(60)
///
/// let as_value = to_smallest_signed_int(1000);
/// // Result: BufValue::Int16(1000)
/// ```
#[must_use]
pub fn to_smallest_signed_int(value: i32) -> BufValue {
    if let Ok(v) = i8::try_from(value) {
        BufValue::Int8(v)
    } else if let Ok(v) = i16::try_from(value) {
        BufValue::Int16(v)
    } else {
        BufValue::Int32(value)
    }
}

/// Converts an integer `BufValue` from any integer type to the smallest signed type.
///
/// This is used to normalize integer tags that may have been read with varying types
/// (`Int32`, `UInt8`, etc.) to signed integer format that matches fgbio's encoding.
///
/// # Arguments
///
/// * `value` - The `BufValue` to convert
///
/// # Returns
///
/// `Some(BufValue)` with the smallest signed type if the value is an integer type,
/// `None` if not an integer type
#[must_use]
pub fn buf_value_to_smallest_signed_int(value: &BufValue) -> Option<BufValue> {
    let int_value = match value {
        BufValue::Int8(i) => i32::from(*i),
        BufValue::Int16(i) => i32::from(*i),
        BufValue::Int32(i) => *i,
        BufValue::UInt8(i) => i32::from(*i),
        BufValue::UInt16(i) => i32::from(*i),
        BufValue::UInt32(i) => i32::try_from(*i).ok()?,
        _ => return None,
    };
    Some(to_smallest_signed_int(int_value))
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::record_buf::data::field::value::Array;

    #[test]
    fn test_reverse_buf_value_int8_array() {
        let value = BufValue::Array(Array::Int8(vec![1, 2, 3, 4, 5]));
        let reversed = reverse_buf_value(&value);

        if let BufValue::Array(Array::Int8(vals)) = reversed {
            assert_eq!(vals, vec![5, 4, 3, 2, 1]);
        } else {
            panic!("Expected Int8 array");
        }
    }

    #[test]
    fn test_reverse_buf_value_uint8_array() {
        let value = BufValue::Array(Array::UInt8(vec![10, 20, 30]));
        let reversed = reverse_buf_value(&value);

        if let BufValue::Array(Array::UInt8(vals)) = reversed {
            assert_eq!(vals, vec![30, 20, 10]);
        } else {
            panic!("Expected UInt8 array");
        }
    }

    #[test]
    fn test_reverse_buf_value_int16_array() {
        let value = BufValue::Array(Array::Int16(vec![100, 200, 300]));
        let reversed = reverse_buf_value(&value);

        if let BufValue::Array(Array::Int16(vals)) = reversed {
            assert_eq!(vals, vec![300, 200, 100]);
        } else {
            panic!("Expected Int16 array");
        }
    }

    #[test]
    fn test_reverse_buf_value_uint16_array() {
        let value = BufValue::Array(Array::UInt16(vec![1000, 2000]));
        let reversed = reverse_buf_value(&value);

        if let BufValue::Array(Array::UInt16(vals)) = reversed {
            assert_eq!(vals, vec![2000, 1000]);
        } else {
            panic!("Expected UInt16 array");
        }
    }

    #[test]
    fn test_reverse_buf_value_int32_array() {
        let value = BufValue::Array(Array::Int32(vec![10000, 20000, 30000]));
        let reversed = reverse_buf_value(&value);

        if let BufValue::Array(Array::Int32(vals)) = reversed {
            assert_eq!(vals, vec![30000, 20000, 10000]);
        } else {
            panic!("Expected Int32 array");
        }
    }

    #[test]
    fn test_reverse_buf_value_uint32_array() {
        let value = BufValue::Array(Array::UInt32(vec![100_000, 200_000]));
        let reversed = reverse_buf_value(&value);

        if let BufValue::Array(Array::UInt32(vals)) = reversed {
            assert_eq!(vals, vec![200_000, 100_000]);
        } else {
            panic!("Expected UInt32 array");
        }
    }

    #[test]
    fn test_reverse_buf_value_float_array() {
        let value = BufValue::Array(Array::Float(vec![1.1, 2.2, 3.3]));
        let reversed = reverse_buf_value(&value);

        if let BufValue::Array(Array::Float(vals)) = reversed {
            assert!((vals[0] - 3.3).abs() < 0.001);
            assert!((vals[1] - 2.2).abs() < 0.001);
            assert!((vals[2] - 1.1).abs() < 0.001);
        } else {
            panic!("Expected Float array");
        }
    }

    #[test]
    fn test_reverse_buf_value_string() {
        let value = BufValue::from("abcde".to_string());
        let reversed = reverse_buf_value(&value);

        if let BufValue::String(s) = reversed {
            assert_eq!(s.to_string(), "edcba");
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_revcomp_buf_value_simple() {
        let value = BufValue::from("ACGT".to_string());
        let revcomp = revcomp_buf_value(&value);

        if let BufValue::String(s) = revcomp {
            assert_eq!(s.to_string(), "ACGT"); // reverse of "TGCA" -> "ACGT"
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_revcomp_buf_value_lowercase() {
        let value = BufValue::from("acgt".to_string());
        let revcomp = revcomp_buf_value(&value);

        // Lowercase is normalized to uppercase
        if let BufValue::String(s) = revcomp {
            assert_eq!(s.to_string(), "ACGT");
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_revcomp_buf_value_mixed_case() {
        let value = BufValue::from("AcGt".to_string());
        let revcomp = revcomp_buf_value(&value);

        // Mixed case is normalized to uppercase
        if let BufValue::String(s) = revcomp {
            assert_eq!(s.to_string(), "ACGT");
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_revcomp_buf_value_with_n() {
        let value = BufValue::from("ACGTN".to_string());
        let revcomp = revcomp_buf_value(&value);

        if let BufValue::String(s) = revcomp {
            assert_eq!(s.to_string(), "NACGT");
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_revcomp_buf_value_complex() {
        let value = BufValue::from("AAAGG".to_string());
        let revcomp = revcomp_buf_value(&value);

        if let BufValue::String(s) = revcomp {
            assert_eq!(s.to_string(), "CCTTT");
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_revcomp_buf_value_array() {
        let value = BufValue::Array(Array::UInt8(vec![1, 2, 3]));
        let revcomp = revcomp_buf_value(&value);

        // Should return unchanged for non-string types
        if let BufValue::Array(Array::UInt8(vals)) = revcomp {
            assert_eq!(vals, vec![1, 2, 3]);
        } else {
            panic!("Expected UInt8 array");
        }
    }

    #[test]
    fn test_reverse_and_revcomp_combined() {
        // Test that reverse and revcomp work correctly together
        let original = BufValue::from("ACGT".to_string());

        // Reverse: "ACGT" -> "TGCA"
        let reversed = reverse_buf_value(&original);
        if let BufValue::String(s) = reversed {
            assert_eq!(s.to_string(), "TGCA");
        }

        // Revcomp: "ACGT" -> "ACGT"
        let revcomped = revcomp_buf_value(&original);
        if let BufValue::String(s) = revcomped {
            assert_eq!(s.to_string(), "ACGT");
        }
    }

    #[test]
    fn test_empty_string_operations() {
        let value = BufValue::from(String::new());

        let reversed = reverse_buf_value(&value);
        if let BufValue::String(s) = reversed {
            assert_eq!(s.to_string(), "");
        }

        let revcomped = revcomp_buf_value(&value);
        if let BufValue::String(s) = revcomped {
            assert_eq!(s.to_string(), "");
        }
    }

    #[test]
    fn test_empty_array_operations() {
        let value = BufValue::Array(Array::UInt8(vec![]));
        let reversed = reverse_buf_value(&value);

        if let BufValue::Array(Array::UInt8(vals)) = reversed {
            assert!(vals.is_empty());
        } else {
            panic!("Expected empty UInt8 array");
        }
    }

    // =========================================================================
    // Tests for is_sorted()
    // =========================================================================

    fn create_header_with_so(sort_order: &str) -> Header {
        let header_str = format!("@HD\tVN:1.6\tSO:{sort_order}\n");
        header_str.parse().unwrap()
    }

    fn create_header_without_so() -> Header {
        let header_str = "@HD\tVN:1.6\n";
        header_str.parse().unwrap()
    }

    fn create_empty_header() -> Header {
        Header::default()
    }

    #[test]
    fn test_is_sorted_queryname_matches() {
        use noodles::sam::header::record::value::map::header::sort_order::QUERY_NAME;
        let header = create_header_with_so("queryname");
        assert!(is_sorted(&header, QUERY_NAME));
    }

    #[test]
    fn test_is_sorted_coordinate_matches() {
        use noodles::sam::header::record::value::map::header::sort_order::COORDINATE;
        let header = create_header_with_so("coordinate");
        assert!(is_sorted(&header, COORDINATE));
    }

    #[test]
    fn test_is_sorted_unsorted_matches() {
        let header = create_header_with_so("unsorted");
        assert!(is_sorted(&header, UNSORTED));
    }

    #[test]
    fn test_is_sorted_mismatch() {
        use noodles::sam::header::record::value::map::header::sort_order::COORDINATE;
        let header = create_header_with_so("queryname");
        assert!(!is_sorted(&header, COORDINATE));
    }

    #[test]
    fn test_is_sorted_no_so_tag() {
        use noodles::sam::header::record::value::map::header::sort_order::COORDINATE;
        let header = create_header_without_so();
        assert!(!is_sorted(&header, COORDINATE));
    }

    #[test]
    fn test_is_sorted_empty_header() {
        use noodles::sam::header::record::value::map::header::sort_order::COORDINATE;
        let header = create_empty_header();
        assert!(!is_sorted(&header, COORDINATE));
    }

    // =========================================================================
    // Tests for is_template_coordinate_sorted()
    // =========================================================================

    fn create_template_coord_header(so: &str, go: Option<&str>, ss: Option<&str>) -> Header {
        use std::fmt::Write;
        let mut header_str = format!("@HD\tVN:1.6\tSO:{so}");
        if let Some(go_val) = go {
            write!(header_str, "\tGO:{go_val}").unwrap();
        }
        if let Some(ss_val) = ss {
            write!(header_str, "\tSS:{ss_val}").unwrap();
        }
        header_str.push('\n');
        header_str.parse().unwrap()
    }

    #[test]
    fn test_is_template_coordinate_sorted_valid_minimal() {
        // SO:unsorted + GO:query (no SS)
        let header = create_template_coord_header("unsorted", Some("query"), None);
        assert!(is_template_coordinate_sorted(&header));
    }

    #[test]
    fn test_is_template_coordinate_sorted_valid_with_ss() {
        // SO:unsorted + GO:query + SS:template-coordinate
        let header =
            create_template_coord_header("unsorted", Some("query"), Some("template-coordinate"));
        assert!(is_template_coordinate_sorted(&header));
    }

    #[test]
    fn test_is_template_coordinate_sorted_valid_with_prefixed_ss() {
        // SO:unsorted + GO:query + SS:unsorted:template-coordinate
        let header = create_template_coord_header(
            "unsorted",
            Some("query"),
            Some("unsorted:template-coordinate"),
        );
        assert!(is_template_coordinate_sorted(&header));
    }

    #[test]
    fn test_is_template_coordinate_sorted_invalid_so_coordinate() {
        // SO:coordinate is wrong
        let header = create_template_coord_header("coordinate", Some("query"), None);
        assert!(!is_template_coordinate_sorted(&header));
    }

    #[test]
    fn test_is_template_coordinate_sorted_invalid_so_queryname() {
        // SO:queryname is wrong
        let header = create_template_coord_header("queryname", Some("query"), None);
        assert!(!is_template_coordinate_sorted(&header));
    }

    #[test]
    fn test_is_template_coordinate_sorted_invalid_go_none() {
        // GO:none is wrong
        let header = create_template_coord_header("unsorted", Some("none"), None);
        assert!(!is_template_coordinate_sorted(&header));
    }

    #[test]
    fn test_is_template_coordinate_sorted_invalid_go_reference() {
        // GO:reference is wrong
        let header = create_template_coord_header("unsorted", Some("reference"), None);
        assert!(!is_template_coordinate_sorted(&header));
    }

    #[test]
    fn test_is_template_coordinate_sorted_missing_go() {
        // Missing GO tag
        let header = create_template_coord_header("unsorted", None, None);
        assert!(!is_template_coordinate_sorted(&header));
    }

    #[test]
    fn test_is_template_coordinate_sorted_invalid_ss() {
        // Wrong SS value
        let header = create_template_coord_header("unsorted", Some("query"), Some("coordinate"));
        assert!(!is_template_coordinate_sorted(&header));
    }

    #[test]
    fn test_is_template_coordinate_sorted_empty_header() {
        let header = create_empty_header();
        assert!(!is_template_coordinate_sorted(&header));
    }

    // =========================================================================
    // Tests for reverse_buf_value with scalar types
    // =========================================================================

    #[test]
    fn test_reverse_buf_value_integer_unchanged() {
        let value = BufValue::from(42_i32);
        let reversed = reverse_buf_value(&value);
        assert_eq!(reversed, BufValue::from(42_i32));
    }

    #[test]
    fn test_reverse_buf_value_character_unchanged() {
        let value = BufValue::Character(b'X');
        let reversed = reverse_buf_value(&value);
        assert_eq!(reversed, BufValue::Character(b'X'));
    }

    // =========================================================================
    // Tests for revcomp_buf_value edge cases
    // =========================================================================

    #[test]
    fn test_revcomp_buf_value_ambiguous_bases_unchanged() {
        // IUPAC ambiguity codes are not complemented
        let value = BufValue::from("RYSWKM".to_string());
        let revcomp = revcomp_buf_value(&value);
        if let BufValue::String(s) = revcomp {
            // Reversed but not complemented (our function only complements ACGTN)
            assert_eq!(s.to_string(), "MKWSYR");
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_revcomp_buf_value_integer_unchanged() {
        let value = BufValue::from(42_i32);
        let revcomp = revcomp_buf_value(&value);
        assert_eq!(revcomp, BufValue::from(42_i32));
    }

    // =========================================================================
    // Tests for to_smallest_signed_int()
    // =========================================================================

    #[test]
    fn test_to_smallest_signed_int_fits_in_i8() {
        // Values that fit in i8 (-128 to 127) should be Int8
        assert_eq!(to_smallest_signed_int(0), BufValue::Int8(0));
        assert_eq!(to_smallest_signed_int(60), BufValue::Int8(60));
        assert_eq!(to_smallest_signed_int(127), BufValue::Int8(127));
        assert_eq!(to_smallest_signed_int(-128), BufValue::Int8(-128));
        assert_eq!(to_smallest_signed_int(-1), BufValue::Int8(-1));
    }

    #[test]
    fn test_to_smallest_signed_int_fits_in_i16() {
        // Values that fit in i16 but not i8 should be Int16
        assert_eq!(to_smallest_signed_int(128), BufValue::Int16(128));
        assert_eq!(to_smallest_signed_int(1000), BufValue::Int16(1000));
        assert_eq!(to_smallest_signed_int(32767), BufValue::Int16(32767));
        assert_eq!(to_smallest_signed_int(-129), BufValue::Int16(-129));
        assert_eq!(to_smallest_signed_int(-32768), BufValue::Int16(-32768));
    }

    #[test]
    fn test_to_smallest_signed_int_requires_i32() {
        // Values that require i32 should be Int32
        assert_eq!(to_smallest_signed_int(32768), BufValue::Int32(32768));
        assert_eq!(to_smallest_signed_int(100_000), BufValue::Int32(100_000));
        assert_eq!(to_smallest_signed_int(-32769), BufValue::Int32(-32769));
    }

    // =========================================================================
    // Tests for buf_value_to_smallest_signed_int()
    // =========================================================================

    #[test]
    fn test_buf_value_to_smallest_signed_int_from_uint8() {
        // UInt8 values should be converted to the smallest signed type
        let value = BufValue::UInt8(60);
        assert_eq!(buf_value_to_smallest_signed_int(&value), Some(BufValue::Int8(60)));

        let value = BufValue::UInt8(200);
        assert_eq!(buf_value_to_smallest_signed_int(&value), Some(BufValue::Int16(200)));
    }

    #[test]
    fn test_buf_value_to_smallest_signed_int_from_int8() {
        // Int8 values should stay as Int8
        let value = BufValue::Int8(60);
        assert_eq!(buf_value_to_smallest_signed_int(&value), Some(BufValue::Int8(60)));

        let value = BufValue::Int8(-50);
        assert_eq!(buf_value_to_smallest_signed_int(&value), Some(BufValue::Int8(-50)));
    }

    #[test]
    fn test_buf_value_to_smallest_signed_int_from_int16() {
        // Int16 values should be converted appropriately
        let value = BufValue::Int16(60);
        assert_eq!(buf_value_to_smallest_signed_int(&value), Some(BufValue::Int8(60)));

        let value = BufValue::Int16(1000);
        assert_eq!(buf_value_to_smallest_signed_int(&value), Some(BufValue::Int16(1000)));
    }

    #[test]
    fn test_buf_value_to_smallest_signed_int_from_int32() {
        // Int32 values should be converted appropriately
        let value = BufValue::Int32(60);
        assert_eq!(buf_value_to_smallest_signed_int(&value), Some(BufValue::Int8(60)));

        let value = BufValue::Int32(1000);
        assert_eq!(buf_value_to_smallest_signed_int(&value), Some(BufValue::Int16(1000)));

        let value = BufValue::Int32(100_000);
        assert_eq!(buf_value_to_smallest_signed_int(&value), Some(BufValue::Int32(100_000)));
    }

    #[test]
    fn test_buf_value_to_smallest_signed_int_non_integer_returns_none() {
        // Non-integer types should return None
        let value = BufValue::from("hello".to_string());
        assert_eq!(buf_value_to_smallest_signed_int(&value), None);

        let value = BufValue::Character(b'X');
        assert_eq!(buf_value_to_smallest_signed_int(&value), None);
    }
}
