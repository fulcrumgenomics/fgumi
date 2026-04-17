//! Raw byte-level BAM record comparison helpers.
//!
//! Provides functions to compare BAM records at the byte level, avoiding the overhead
//! of full record decoding. Supports full-record comparison, core-field-only comparison,
//! tag comparison (both order-dependent and order-independent), and a structured
//! comparison that reports which parts differ.

use fgumi_raw_bam::fields::{aux_data_offset_from_record, tag_value_size};

/// Result of a structured raw BAM record comparison.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RawCompareResult {
    /// Whether the core fields (everything before aux data) are identical.
    pub core_match: bool,
    /// Whether the aux tags are byte-identical (same order and content).
    pub tags_match: bool,
    /// Whether the aux tags are semantically identical regardless of order.
    pub tag_order_match: bool,
}

/// Returns `true` if two raw BAM records are byte-identical.
#[must_use]
pub fn raw_records_byte_equal(r1: &[u8], r2: &[u8]) -> bool {
    r1 == r2
}

/// Offset of the BAM `bin` field in the fixed header.
const BIN_OFFSET: usize = 10;
/// Length of the BAM `bin` field (u16).
const BIN_LEN: usize = 2;

/// Returns `true` if the core fields (everything up to but not including aux data) are identical.
///
/// The BAM `bin` field (bytes 10-11) is excluded from comparison because it is a
/// BAM-specific index optimization that is not part of the SAM data model.  Different
/// tools may compute different bin values for the same alignment, so it should not
/// cause records to be reported as differing.
///
/// Returns `false` if either record is too short to determine the aux data offset.
#[must_use]
pub fn raw_core_fields_equal(r1: &[u8], r2: &[u8]) -> bool {
    let (Some(off1), Some(off2)) =
        (aux_data_offset_from_record(r1), aux_data_offset_from_record(r2))
    else {
        return false;
    };
    let end1 = off1.min(r1.len());
    let end2 = off2.min(r2.len());
    if end1 != end2 {
        return false;
    }
    // Compare bytes before the bin field, then after, skipping bytes 10-11.
    r1[..BIN_OFFSET] == r2[..BIN_OFFSET]
        && r1[BIN_OFFSET + BIN_LEN..end1] == r2[BIN_OFFSET + BIN_LEN..end2]
}

/// Returns `true` if the aux data regions are byte-identical (order-sensitive).
#[must_use]
pub fn raw_tags_byte_equal(r1: &[u8], r2: &[u8]) -> bool {
    let aux1 = fgumi_raw_bam::fields::aux_data_slice(r1);
    let aux2 = fgumi_raw_bam::fields::aux_data_slice(r2);
    aux1 == aux2
}

/// Collects individual tag entries from an aux data byte slice.
///
/// Each entry is `(tag_name, value_start, value_end)` where the offsets are relative
/// to the start of the aux data slice. `value_start` is the offset of the type byte,
/// and `value_end` is one past the last byte of the tag value.
///
/// Returns `None` if the aux data is malformed (truncated tag or unknown type).
fn collect_tag_entries(aux_data: &[u8]) -> Option<Vec<([u8; 2], usize, usize)>> {
    let mut entries = Vec::new();
    let mut i = 0;
    while i < aux_data.len() {
        if i + 3 > aux_data.len() {
            return None;
        }
        let tag_name = [aux_data[i], aux_data[i + 1]];
        let val_type = aux_data[i + 2];
        let val_start = i + 2; // include type byte in the "value" region for comparison
        let data_start = i + 3;
        let size = tag_value_size(val_type, aux_data.get(data_start..)?)?;
        let val_end = data_start + size;
        entries.push((tag_name, val_start, val_end));
        i = val_end;
    }
    Some(entries)
}

/// Returns `true` if the given BAM aux type byte is an integer type.
///
/// BAM integer types: `c` (i8), `C` (u8), `s` (i16), `S` (u16), `i` (i32), `I` (u32).
fn is_int_type(t: u8) -> bool {
    matches!(t, b'c' | b'C' | b's' | b'S' | b'i' | b'I')
}

/// Decodes a BAM integer tag value to `i64` for semantic comparison.
///
/// `type_byte` is the BAM type character and `data` is the value bytes (excluding the type byte).
/// Returns `None` if the data is too short for the given type.
fn decode_int_tag(type_byte: u8, data: &[u8]) -> Option<i64> {
    match type_byte {
        b'c' => data.first().map(|&b| i64::from(b as i8)),
        b'C' => data.first().map(|&b| i64::from(b)),
        b's' => data.get(..2).map(|b| i64::from(i16::from_le_bytes([b[0], b[1]]))),
        b'S' => data.get(..2).map(|b| i64::from(u16::from_le_bytes([b[0], b[1]]))),
        b'i' => data.get(..4).map(|b| i64::from(i32::from_le_bytes([b[0], b[1], b[2], b[3]]))),
        b'I' => data.get(..4).map(|b| i64::from(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))),
        _ => None,
    }
}

/// Returns `true` if the aux data regions contain the same tags with the same values,
/// regardless of the order in which tags appear.
///
/// Integer tag values are compared semantically: if both tags have integer types
/// (`c`/`C`/`s`/`S`/`i`/`I`) but different encodings (e.g., `C` for u8 vs `s` for i16),
/// they are decoded to `i64` and compared by value.
///
/// Note: Assumes tag names are unique within each record per the BAM specification.
/// If duplicate tag names exist, comparison results are undefined.
///
/// Returns `false` if either record has malformed aux data.
#[must_use]
pub fn raw_tags_equal_order_independent(r1: &[u8], r2: &[u8]) -> bool {
    let aux1 = fgumi_raw_bam::fields::aux_data_slice(r1);
    let aux2 = fgumi_raw_bam::fields::aux_data_slice(r2);

    let (Some(entries1), Some(entries2)) = (collect_tag_entries(aux1), collect_tag_entries(aux2))
    else {
        return false;
    };

    if entries1.len() != entries2.len() {
        return false;
    }

    // For each tag in r1, find the matching tag in r2 and compare values.
    for &(tag1, start1, end1) in &entries1 {
        let Some(&(_, start2, end2)) = entries2.iter().find(|(t, _, _)| *t == tag1) else {
            return false;
        };
        // Fast path: byte-identical values (includes type byte)
        if aux1[start1..end1] == aux2[start2..end2] {
            continue;
        }
        // If both are integer types, compare semantically by decoding to i64
        let type1 = aux1[start1];
        let type2 = aux2[start2];
        if is_int_type(type1) && is_int_type(type2) {
            let data1 = &aux1[start1 + 1..end1];
            let data2 = &aux2[start2 + 1..end2];
            match (decode_int_tag(type1, data1), decode_int_tag(type2, data2)) {
                (Some(v1), Some(v2)) if v1 == v2 => continue,
                _ => return false,
            }
        } else {
            return false;
        }
    }

    true
}

/// Performs a structured comparison of two raw BAM records, reporting which parts match.
///
/// Returns a [`RawCompareResult`] with `core_match`, `tags_match`, and `tag_order_match`.
/// If core fields cannot be parsed (record too short), `core_match` is `false`.
/// If aux data is malformed, `tag_order_match` is `false`.
#[must_use]
pub fn raw_compare_structured(r1: &[u8], r2: &[u8]) -> RawCompareResult {
    let core_match = raw_core_fields_equal(r1, r2);
    let tags_match = raw_tags_byte_equal(r1, r2);
    let tag_order_match = if tags_match { true } else { raw_tags_equal_order_independent(r1, r2) };
    RawCompareResult { core_match, tags_match, tag_order_match }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::SamTag;
    use fgumi_raw_bam::testutil::make_bam_bytes;

    /// Helper to build aux bytes for a Z-type string tag.
    fn make_z_tag(tag: [u8; 2], value: &[u8]) -> Vec<u8> {
        let mut aux = vec![tag[0], tag[1], b'Z'];
        aux.extend_from_slice(value);
        aux.push(0); // null terminator
        aux
    }

    /// Helper to build aux bytes for an i-type (i32) integer tag.
    fn make_i_tag(tag: [u8; 2], value: i32) -> Vec<u8> {
        let mut aux = vec![tag[0], tag[1], b'i'];
        aux.extend_from_slice(&value.to_le_bytes());
        aux
    }

    /// Helper to build aux bytes for a C-type (u8) integer tag.
    fn make_c_tag(tag: [u8; 2], value: u8) -> Vec<u8> {
        vec![tag[0], tag[1], b'C', value]
    }

    /// Helper to build aux bytes for an s-type (i16) integer tag.
    fn make_s_tag(tag: [u8; 2], value: i16) -> Vec<u8> {
        let mut aux = vec![tag[0], tag[1], b's'];
        aux.extend_from_slice(&value.to_le_bytes());
        aux
    }

    /// Helper to build aux bytes for an S-type (u16) integer tag.
    fn make_upper_s_tag(tag: [u8; 2], value: u16) -> Vec<u8> {
        let mut aux = vec![tag[0], tag[1], b'S'];
        aux.extend_from_slice(&value.to_le_bytes());
        aux
    }

    /// Helper to build aux bytes for an I-type (u32) integer tag.
    fn make_upper_i_tag(tag: [u8; 2], value: u32) -> Vec<u8> {
        let mut aux = vec![tag[0], tag[1], b'I'];
        aux.extend_from_slice(&value.to_le_bytes());
        aux
    }

    /// Helper to build aux bytes for an f-type (f32) float tag.
    fn make_f_tag(tag: [u8; 2], value: f32) -> Vec<u8> {
        let mut aux = vec![tag[0], tag[1], b'f'];
        aux.extend_from_slice(&value.to_le_bytes());
        aux
    }

    fn base_record(aux: &[u8]) -> Vec<u8> {
        make_bam_bytes(1, 100, 0x41, b"rea", &[], 4, 1, 200, aux)
    }

    // ========================================================================
    // raw_records_byte_equal tests
    // ========================================================================

    #[test]
    fn test_byte_equal_identical() {
        let bytes = base_record(&[]);
        assert!(raw_records_byte_equal(&bytes, &bytes));
    }

    #[test]
    fn test_byte_equal_different() {
        let r1 = base_record(&[]);
        let r2 = make_bam_bytes(2, 100, 0x41, b"rea", &[], 4, 1, 200, &[]);
        assert!(!raw_records_byte_equal(&r1, &r2));
    }

    #[test]
    fn test_byte_equal_empty() {
        let r1: &[u8] = &[];
        let r2: &[u8] = &[];
        assert!(raw_records_byte_equal(r1, r2));
    }

    // ========================================================================
    // raw_core_fields_equal tests
    // ========================================================================

    #[test]
    fn test_core_equal_identical_records() {
        let rec = base_record(&[]);
        assert!(raw_core_fields_equal(&rec, &rec));
    }

    #[test]
    fn test_core_equal_same_core_different_tags() {
        let aux1 = make_z_tag(*SamTag::RG, b"sample1");
        let aux2 = make_z_tag(*SamTag::RG, b"sample2");
        let r1 = base_record(&aux1);
        let r2 = base_record(&aux2);
        assert!(raw_core_fields_equal(&r1, &r2));
    }

    #[test]
    fn test_core_not_equal_different_pos() {
        let r1 = make_bam_bytes(1, 100, 0, b"rea", &[], 4, 1, 200, &[]);
        let r2 = make_bam_bytes(1, 200, 0, b"rea", &[], 4, 1, 200, &[]);
        assert!(!raw_core_fields_equal(&r1, &r2));
    }

    #[test]
    fn test_core_not_equal_different_flags() {
        let r1 = make_bam_bytes(1, 100, 0, b"rea", &[], 4, 1, 200, &[]);
        let r2 = make_bam_bytes(1, 100, 0x10, b"rea", &[], 4, 1, 200, &[]);
        assert!(!raw_core_fields_equal(&r1, &r2));
    }

    #[test]
    fn test_core_equal_short_record_returns_false() {
        let short = vec![0u8; 10];
        let rec = base_record(&[]);
        assert!(!raw_core_fields_equal(&short, &rec));
    }

    #[test]
    fn test_core_equal_different_bin_values() {
        // The BAM `bin` field (bytes 10-11) is not semantically meaningful and may
        // differ between tools. Core comparison should ignore it.
        let mut r1 = base_record(&[]);
        let mut r2 = base_record(&[]);
        // Set different bin values
        r1[10..12].copy_from_slice(&100u16.to_le_bytes());
        r2[10..12].copy_from_slice(&200u16.to_le_bytes());
        assert!(raw_core_fields_equal(&r1, &r2));
    }

    // ========================================================================
    // raw_tags_byte_equal tests
    // ========================================================================

    #[test]
    fn test_tags_byte_equal_no_tags() {
        let r1 = base_record(&[]);
        let r2 = base_record(&[]);
        assert!(raw_tags_byte_equal(&r1, &r2));
    }

    #[test]
    fn test_tags_byte_equal_same_tags() {
        let aux = make_z_tag(*SamTag::RG, b"sample1");
        let r1 = base_record(&aux);
        let r2 = base_record(&aux);
        assert!(raw_tags_byte_equal(&r1, &r2));
    }

    #[test]
    fn test_tags_byte_not_equal_different_values() {
        let r1 = base_record(&make_z_tag(*SamTag::RG, b"sample1"));
        let r2 = base_record(&make_z_tag(*SamTag::RG, b"sample2"));
        assert!(!raw_tags_byte_equal(&r1, &r2));
    }

    #[test]
    fn test_tags_byte_not_equal_different_order() {
        let mut aux1 = make_z_tag(*SamTag::RG, b"s1");
        aux1.extend_from_slice(&make_i_tag(*b"NM", 5));
        let mut aux2 = make_i_tag(*b"NM", 5);
        aux2.extend_from_slice(&make_z_tag(*SamTag::RG, b"s1"));
        let r1 = base_record(&aux1);
        let r2 = base_record(&aux2);
        assert!(!raw_tags_byte_equal(&r1, &r2));
    }

    // ========================================================================
    // raw_tags_equal_order_independent tests
    // ========================================================================

    #[test]
    fn test_order_independent_same_order() {
        let mut aux = make_z_tag(*SamTag::RG, b"s1");
        aux.extend_from_slice(&make_i_tag(*b"NM", 5));
        let r1 = base_record(&aux);
        let r2 = base_record(&aux);
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_different_order() {
        let mut aux1 = make_z_tag(*SamTag::RG, b"s1");
        aux1.extend_from_slice(&make_i_tag(*b"NM", 5));
        let mut aux2 = make_i_tag(*b"NM", 5);
        aux2.extend_from_slice(&make_z_tag(*SamTag::RG, b"s1"));
        let r1 = base_record(&aux1);
        let r2 = base_record(&aux2);
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_different_values() {
        let r1 = base_record(&make_i_tag(*b"NM", 5));
        let r2 = base_record(&make_i_tag(*b"NM", 10));
        assert!(!raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_missing_tag() {
        let mut aux1 = make_z_tag(*SamTag::RG, b"s1");
        aux1.extend_from_slice(&make_i_tag(*b"NM", 5));
        let aux2 = make_z_tag(*SamTag::RG, b"s1");
        let r1 = base_record(&aux1);
        let r2 = base_record(&aux2);
        assert!(!raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_extra_tag() {
        let aux1 = make_i_tag(*b"NM", 5);
        let mut aux2 = make_i_tag(*b"NM", 5);
        aux2.extend_from_slice(&make_z_tag(*SamTag::RG, b"s1"));
        let r1 = base_record(&aux1);
        let r2 = base_record(&aux2);
        assert!(!raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_no_tags() {
        let r1 = base_record(&[]);
        let r2 = base_record(&[]);
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_multiple_types() {
        // Mix of Z, i, C, and f tags in different orders
        let mut aux1 = make_z_tag(*SamTag::RG, b"grp");
        aux1.extend_from_slice(&make_i_tag(*b"NM", 3));
        aux1.extend_from_slice(&make_c_tag(*b"MQ", 42));
        aux1.extend_from_slice(&make_f_tag(*b"GC", 0.75));

        let mut aux2 = make_f_tag(*b"GC", 0.75);
        aux2.extend_from_slice(&make_c_tag(*b"MQ", 42));
        aux2.extend_from_slice(&make_z_tag(*SamTag::RG, b"grp"));
        aux2.extend_from_slice(&make_i_tag(*b"NM", 3));

        let r1 = base_record(&aux1);
        let r2 = base_record(&aux2);
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_short_record() {
        let short = vec![0u8; 10];
        let rec = base_record(&make_i_tag(*b"NM", 5));
        // Both have empty aux slices due to short record, but the short record
        // returns empty aux, so it won't match a record with tags.
        assert!(!raw_tags_equal_order_independent(&short, &rec));
    }

    // ========================================================================
    // raw_compare_structured tests
    // ========================================================================

    #[test]
    fn test_structured_identical() {
        let aux = make_i_tag(*b"NM", 5);
        let rec = base_record(&aux);
        let result = raw_compare_structured(&rec, &rec);
        assert_eq!(
            result,
            RawCompareResult { core_match: true, tags_match: true, tag_order_match: true }
        );
    }

    #[test]
    fn test_structured_different_core() {
        let r1 = make_bam_bytes(1, 100, 0, b"rea", &[], 4, 1, 200, &[]);
        let r2 = make_bam_bytes(1, 999, 0, b"rea", &[], 4, 1, 200, &[]);
        let result = raw_compare_structured(&r1, &r2);
        assert!(!result.core_match);
        assert!(result.tags_match); // both have no tags
    }

    #[test]
    fn test_structured_same_core_different_tag_order() {
        let mut aux1 = make_z_tag(*SamTag::RG, b"s1");
        aux1.extend_from_slice(&make_i_tag(*b"NM", 5));
        let mut aux2 = make_i_tag(*b"NM", 5);
        aux2.extend_from_slice(&make_z_tag(*SamTag::RG, b"s1"));
        let r1 = base_record(&aux1);
        let r2 = base_record(&aux2);
        let result = raw_compare_structured(&r1, &r2);
        assert!(result.core_match);
        assert!(!result.tags_match);
        assert!(result.tag_order_match);
    }

    #[test]
    fn test_structured_different_tags() {
        let r1 = base_record(&make_i_tag(*b"NM", 5));
        let r2 = base_record(&make_i_tag(*b"NM", 10));
        let result = raw_compare_structured(&r1, &r2);
        assert!(result.core_match);
        assert!(!result.tags_match);
        assert!(!result.tag_order_match);
    }

    #[test]
    fn test_structured_short_records() {
        let short = vec![0u8; 10];
        let rec = base_record(&[]);
        let result = raw_compare_structured(&short, &rec);
        assert!(!result.core_match);
    }

    #[test]
    fn test_structured_different_bin_different_tag_order() {
        // Records with different bin values and different tag order should report:
        // core_match=true (bin is not semantically meaningful),
        // tags_match=false (byte order differs), tag_order_match=true (values match)
        let mut aux1 = make_z_tag(*SamTag::RG, b"s1");
        aux1.extend_from_slice(&make_i_tag(*b"NM", 5));
        let mut aux2 = make_i_tag(*b"NM", 5);
        aux2.extend_from_slice(&make_z_tag(*SamTag::RG, b"s1"));
        let mut r1 = base_record(&aux1);
        let mut r2 = base_record(&aux2);
        // Set different bin values
        r1[10..12].copy_from_slice(&100u16.to_le_bytes());
        r2[10..12].copy_from_slice(&200u16.to_le_bytes());
        let result = raw_compare_structured(&r1, &r2);
        assert!(result.core_match);
        assert!(!result.tags_match);
        assert!(result.tag_order_match);
    }

    // ========================================================================
    // collect_tag_entries tests
    // ========================================================================

    #[test]
    fn test_collect_entries_empty() {
        let entries =
            collect_tag_entries(&[]).expect("collect_tag_entries should succeed for empty input");
        assert!(entries.is_empty());
    }

    #[test]
    fn test_collect_entries_single_int() {
        let aux = make_i_tag(*b"NM", 42);
        let entries = collect_tag_entries(&aux)
            .expect("collect_tag_entries should succeed for single int tag");
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].0, *b"NM");
    }

    #[test]
    fn test_collect_entries_multiple() {
        let mut aux = make_z_tag(*SamTag::RG, b"grp");
        aux.extend_from_slice(&make_i_tag(*b"NM", 5));
        let entries = collect_tag_entries(&aux)
            .expect("collect_tag_entries should succeed for multiple tags");
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].0, *SamTag::RG);
        assert_eq!(entries[1].0, *b"NM");
    }

    #[test]
    fn test_collect_entries_truncated_returns_none() {
        // Just a tag name with no type byte
        let aux = [b'N', b'M'];
        assert!(collect_tag_entries(&aux).is_none());
    }

    #[test]
    fn test_collect_entries_b_array_tag() {
        let aux = fgumi_raw_bam::testutil::make_b_int_array_tag(*b"XA", &[1, 2, 3]);
        let entries =
            collect_tag_entries(&aux).expect("collect_tag_entries should succeed for B-array tag");
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].0, *b"XA");
        // The entry should span the entire aux data
        assert_eq!(entries[0].2, aux.len());
    }

    // ========================================================================
    // Semantic integer comparison tests
    // ========================================================================

    #[test]
    fn test_order_independent_same_value_different_int_types_u8_vs_i16() {
        // cD tag: u8 (C) value 158 vs i16 (s) value 158
        // These are semantically equal but have different byte encodings.
        let r1 = base_record(&make_c_tag(*b"cD", 158));
        let r2 = base_record(&make_s_tag(*b"cD", 158));
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_same_value_different_int_types_u8_vs_i32() {
        // NM tag: u8 (C) value 42 vs i32 (i) value 42
        let r1 = base_record(&make_c_tag(*b"NM", 42));
        let r2 = base_record(&make_i_tag(*b"NM", 42));
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_same_value_different_int_types_i16_vs_i32() {
        // MQ tag: i16 (s) value 300 vs i32 (i) value 300
        let r1 = base_record(&make_s_tag(*b"MQ", 300));
        let r2 = base_record(&make_i_tag(*b"MQ", 300));
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_same_value_different_int_types_u32_vs_i32() {
        // NM tag: u32 (I) value 1000 vs i32 (i) value 1000
        let r1 = base_record(&make_upper_i_tag(*b"NM", 1000));
        let r2 = base_record(&make_i_tag(*b"NM", 1000));
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_same_value_different_int_types_u32_vs_u8() {
        // NM tag: u32 (I) value 42 vs u8 (C) value 42
        let r1 = base_record(&make_upper_i_tag(*b"NM", 42));
        let r2 = base_record(&make_c_tag(*b"NM", 42));
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_different_values_different_int_types() {
        // Same tag name, different int types, different values => should NOT match
        let r1 = base_record(&make_c_tag(*b"NM", 5));
        let r2 = base_record(&make_i_tag(*b"NM", 10));
        assert!(!raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_int_vs_non_int_same_tag_not_equal() {
        // Same tag name, one is int (C), other is string (Z) => should NOT match
        let r1 = base_record(&make_c_tag(*b"XY", 65)); // 65 = 'A'
        let r2 = base_record(&make_z_tag(*b"XY", b"A"));
        assert!(!raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_structured_semantic_int_tags_match() {
        // Same value with different int encodings: structured comparison should
        // report tags_match=false (bytes differ) but tag_order_match=true (semantically equal)
        let r1 = base_record(&make_c_tag(*b"cD", 158));
        let r2 = base_record(&make_s_tag(*b"cD", 158));
        let result = raw_compare_structured(&r1, &r2);
        assert!(result.core_match);
        assert!(!result.tags_match); // raw bytes differ
        assert!(result.tag_order_match); // semantically equal
    }

    #[test]
    fn test_order_independent_semantic_int_with_reordered_tags() {
        // Multiple tags, different order, AND different integer encodings
        let mut aux1 = make_z_tag(*SamTag::RG, b"grp");
        aux1.extend_from_slice(&make_c_tag(*b"cD", 200));
        aux1.extend_from_slice(&make_i_tag(*b"NM", 5));

        let mut aux2 = make_i_tag(*b"NM", 5);
        aux2.extend_from_slice(&make_upper_s_tag(*b"cD", 200));
        aux2.extend_from_slice(&make_z_tag(*SamTag::RG, b"grp"));

        let r1 = base_record(&aux1);
        let r2 = base_record(&aux2);
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_negative_value_across_types() {
        // -100 as i8 (c) vs -100 as i16 (s) => should match
        let r1 = base_record(&[b'X', b'N', b'c', (-100i8) as u8]);
        let r2 = base_record(&make_s_tag(*b"XN", -100));
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_same_bytes_different_sign_interpretation() {
        // i8(-1) = 0xFF vs u8(255) = 0xFF => semantically different, should NOT match
        let r1 = base_record(&[b'X', b'V', b'c', 0xFF]); // -1 as i8
        let r2 = base_record(&make_c_tag(*b"XV", 255)); // 255 as u8
        assert!(!raw_tags_equal_order_independent(&r1, &r2));
    }
}
