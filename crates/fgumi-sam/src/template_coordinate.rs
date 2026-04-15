//! TC (Template-Coordinate) sort-key tag for secondary/supplementary reads.
//!
//! `fgumi zipper` writes a `tc:B:i,<6 int32s>` tag on secondary/supplementary
//! records. The tag carries the primary alignments' template-coordinate sort
//! key so downstream tools (e.g. `fgumi sort --order template-coordinate`,
//! `fgumi dedup`) can keep those records adjacent to their primaries.
//!
//! The tag name is lowercase per the SAM specification convention for
//! non-standard program tags, and is deliberately distinct from bwa-mem's
//! `pa:f:<float>` (primary-alignment score fraction).

use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;

use crate::SamTag;

/// The TC (Template-Coordinate) tag for secondary/supplementary reads.
///
/// Stores the template-coordinate sort key as a `B:i` array (6 int32s).
/// Added by `fgumi zipper` to enable correct downstream handling of
/// secondary/supplementary reads.
pub const TC_TAG: Tag = SamTag::TC.to_noodles_tag();

/// Template-coordinate sort-key info for secondary/supplementary reads.
///
/// Stores the template-coordinate sort key extracted from the primary
/// alignments, enabling secondary/supplementary reads to sort adjacent
/// to their primaries.
///
/// # Binary Format
///
/// Stored as a `B:i` (int32 array) BAM tag with 6 elements:
/// `[tid1, pos1, neg1, tid2, pos2, neg2]`
/// where `neg1`/`neg2` are 0 for forward, 1 for reverse strand.
///
/// This format is faster to parse than a string representation and ensures
/// supplementary reads get the exact same sort key as their primary reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TemplateCoordinateInfo {
    /// Reference ID of the earlier mate (lower position).
    pub tid1: i32,
    /// Unclipped 5' position of the earlier mate.
    pub pos1: i32,
    /// True if earlier mate is on reverse strand.
    pub neg1: bool,
    /// Reference ID of the later mate.
    pub tid2: i32,
    /// Unclipped 5' position of the later mate.
    pub pos2: i32,
    /// True if later mate is on reverse strand.
    pub neg2: bool,
}

impl TemplateCoordinateInfo {
    /// Creates a new `TemplateCoordinateInfo`.
    #[must_use]
    pub const fn new(tid1: i32, pos1: i32, neg1: bool, tid2: i32, pos2: i32, neg2: bool) -> Self {
        Self { tid1, pos1, neg1, tid2, pos2, neg2 }
    }

    /// Serializes to a BAM tag value (B:i array with 6 elements).
    #[must_use]
    pub fn to_tag_value(&self) -> noodles::sam::alignment::record_buf::data::field::Value {
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        let values: Vec<i32> = vec![
            self.tid1,
            self.pos1,
            i32::from(self.neg1),
            self.tid2,
            self.pos2,
            i32::from(self.neg2),
        ];
        Value::Array(Array::Int32(values))
    }

    /// Deserializes from a BAM tag value (B:i array with 6 int32 elements).
    ///
    /// Optimized with a fast path for Int32 arrays (the expected format) that
    /// avoids heap allocation by directly indexing the array.
    #[must_use]
    pub fn from_tag_value(
        value: &noodles::sam::alignment::record_buf::data::field::Value,
    ) -> Option<Self> {
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        match value {
            Value::Array(arr) => {
                // Fast path: Int32 array (expected format from zipper) - no allocation
                if let Array::Int32(v) = arr {
                    if v.len() == 6 {
                        return Some(Self {
                            tid1: v[0],
                            pos1: v[1],
                            neg1: v[2] != 0,
                            tid2: v[3],
                            pos2: v[4],
                            neg2: v[5] != 0,
                        });
                    }
                    return None;
                }

                // Slow path: other array types (rare) - requires allocation
                let values: Vec<i32> = match arr {
                    Array::Int8(v) => v.iter().map(|&x| i32::from(x)).collect(),
                    Array::UInt8(v) => v.iter().map(|&x| i32::from(x)).collect(),
                    Array::Int16(v) => v.iter().map(|&x| i32::from(x)).collect(),
                    Array::UInt16(v) => v.iter().map(|&x| i32::from(x)).collect(),
                    Array::Int32(_) => unreachable!(), // Handled above
                    Array::UInt32(v) => {
                        // Use try_from to avoid wrapping for values > i32::MAX
                        let result: Result<Vec<i32>, _> =
                            v.iter().map(|&x| i32::try_from(x)).collect();
                        match result {
                            Ok(vals) => vals,
                            Err(_) => return None,
                        }
                    }
                    Array::Float(_) => return None,
                };

                if values.len() != 6 {
                    return None;
                }

                Some(Self {
                    tid1: values[0],
                    pos1: values[1],
                    neg1: values[2] != 0,
                    tid2: values[3],
                    pos2: values[4],
                    neg2: values[5] != 0,
                })
            }
            _ => None,
        }
    }

    /// Extracts from a BAM record's TC tag.
    #[must_use]
    pub fn from_record(record: &RecordBuf) -> Option<Self> {
        let value = record.data().get(&TC_TAG)?;
        Self::from_tag_value(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_template_coordinate_info_roundtrip() {
        let info = TemplateCoordinateInfo::new(5, 1000, false, 3, 2000, true);
        let value = info.to_tag_value();
        let parsed = TemplateCoordinateInfo::from_tag_value(&value);

        assert!(parsed.is_some());
        let parsed = parsed.expect("parsing should succeed");
        assert_eq!(parsed.tid1, 5);
        assert_eq!(parsed.pos1, 1000);
        assert!(!parsed.neg1);
        assert_eq!(parsed.tid2, 3);
        assert_eq!(parsed.pos2, 2000);
        assert!(parsed.neg2);
    }

    #[test]
    fn test_template_coordinate_info_from_record() {
        use crate::builder::RecordBuilder;
        use noodles::sam::alignment::record::Flags;

        // Create a supplementary record with tc tag
        let info = TemplateCoordinateInfo::new(0, 100, false, 0, 200, true);
        let record = RecordBuilder::new()
            .name("test")
            .sequence("ACGT")
            .flags(Flags::SUPPLEMENTARY)
            .tag("tc", info.to_tag_value())
            .build();

        let result = TemplateCoordinateInfo::from_record(&record);
        assert!(result.is_some());
        let result = result.expect("result should be Some");
        assert_eq!(result.tid1, 0);
        assert_eq!(result.pos1, 100);
        assert!(!result.neg1);
        assert_eq!(result.tid2, 0);
        assert_eq!(result.pos2, 200);
        assert!(result.neg2);
    }

    #[test]
    fn test_template_coordinate_info_missing() {
        use crate::builder::RecordBuilder;
        use noodles::sam::alignment::record::Flags;

        let record =
            RecordBuilder::new().name("test").sequence("ACGT").flags(Flags::SUPPLEMENTARY).build();

        let result = TemplateCoordinateInfo::from_record(&record);
        assert!(result.is_none());
    }

    // ========================================================================
    // from_tag_value tests (fast path coverage)
    // ========================================================================

    #[test]
    fn test_from_tag_value_int32_fast_path() {
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        // Create Int32 array (the expected format from zipper)
        let values: Vec<i32> = vec![5, 1000, 0, 3, 2000, 1];
        let value = Value::Array(Array::Int32(values));

        let result = TemplateCoordinateInfo::from_tag_value(&value);
        assert!(result.is_some());
        let info = result.expect("result should be Some");
        assert_eq!(info.tid1, 5);
        assert_eq!(info.pos1, 1000);
        assert!(!info.neg1);
        assert_eq!(info.tid2, 3);
        assert_eq!(info.pos2, 2000);
        assert!(info.neg2);
    }

    #[test]
    fn test_from_tag_value_int32_wrong_length() {
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        // Int32 array with wrong number of elements
        let values: Vec<i32> = vec![5, 1000, 0]; // Only 3 elements
        let value = Value::Array(Array::Int32(values));

        let result = TemplateCoordinateInfo::from_tag_value(&value);
        assert!(result.is_none());
    }

    #[test]
    fn test_from_tag_value_int16_fallback() {
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        // Int16 array (rare, but should work via fallback path)
        let values: Vec<i16> = vec![5, 1000, 0, 3, 2000, 1];
        let value = Value::Array(Array::Int16(values));

        let result = TemplateCoordinateInfo::from_tag_value(&value);
        assert!(result.is_some());
        let info = result.expect("result should be Some");
        assert_eq!(info.tid1, 5);
        assert_eq!(info.pos1, 1000);
    }

    #[test]
    fn test_from_tag_value_non_array_returns_none() {
        use noodles::sam::alignment::record_buf::data::field::Value;

        // String value instead of array
        let value = Value::String("not_an_array".into());

        let result = TemplateCoordinateInfo::from_tag_value(&value);
        assert!(result.is_none());
    }

    #[test]
    fn test_from_tag_value_float_array_returns_none() {
        use noodles::sam::alignment::record_buf::data::field::Value;
        use noodles::sam::alignment::record_buf::data::field::value::Array;

        // Float array is not supported
        let values: Vec<f32> = vec![5.0, 1000.0, 0.0, 3.0, 2000.0, 1.0];
        let value = Value::Array(Array::Float(values));

        let result = TemplateCoordinateInfo::from_tag_value(&value);
        assert!(result.is_none());
    }

    // ========================================================================
    // TemplateCoordinateInfo::new edge cases
    // ========================================================================

    #[test]
    fn test_template_coordinate_info_new_stores_values_unchanged() {
        // Verify new() stores values exactly as provided (no normalization)
        let info = TemplateCoordinateInfo::new(5, 1000, true, 3, 500, false);

        // Values should be stored exactly as provided
        assert_eq!(info.tid1, 5);
        assert_eq!(info.pos1, 1000);
        assert!(info.neg1);
        assert_eq!(info.tid2, 3);
        assert_eq!(info.pos2, 500);
        assert!(!info.neg2);
    }

    #[test]
    fn test_template_coordinate_info_new_with_negative_positions() {
        // Edge case: negative positions (can happen with soft clips before position 0)
        let info = TemplateCoordinateInfo::new(0, -5, false, 0, -10, true);

        assert_eq!(info.tid1, 0);
        assert_eq!(info.pos1, -5);
        assert!(!info.neg1);
        assert_eq!(info.tid2, 0);
        assert_eq!(info.pos2, -10);
        assert!(info.neg2);
    }

    #[test]
    fn test_template_coordinate_info_new_with_max_values() {
        // Edge case: maximum i32 values
        let info = TemplateCoordinateInfo::new(i32::MAX, i32::MAX, true, i32::MAX, i32::MAX, true);

        assert_eq!(info.tid1, i32::MAX);
        assert_eq!(info.pos1, i32::MAX);
        assert!(info.neg1);
        assert_eq!(info.tid2, i32::MAX);
        assert_eq!(info.pos2, i32::MAX);
        assert!(info.neg2);
    }

    #[test]
    fn test_template_coordinate_info_new_with_zero_values() {
        // Edge case: all zeros (unmapped or start of reference)
        let info = TemplateCoordinateInfo::new(0, 0, false, 0, 0, false);

        assert_eq!(info.tid1, 0);
        assert_eq!(info.pos1, 0);
        assert!(!info.neg1);
        assert_eq!(info.tid2, 0);
        assert_eq!(info.pos2, 0);
        assert!(!info.neg2);
    }

    #[test]
    fn test_template_coordinate_info_roundtrip_with_negative_positions() {
        // Verify negative positions survive roundtrip
        let info = TemplateCoordinateInfo::new(0, -10, false, 0, -5, true);
        let value = info.to_tag_value();
        let parsed =
            TemplateCoordinateInfo::from_tag_value(&value).expect("parsing should succeed");

        assert_eq!(parsed.tid1, 0);
        assert_eq!(parsed.pos1, -10);
        assert!(!parsed.neg1);
        assert_eq!(parsed.tid2, 0);
        assert_eq!(parsed.pos2, -5);
        assert!(parsed.neg2);
    }

    #[test]
    fn test_template_coordinate_info_position_order_is_caller_responsibility() {
        // Verify that new() does NOT enforce pos1 < pos2 ordering
        // (ordering is done by the caller in zipper.rs)
        let info = TemplateCoordinateInfo::new(0, 2000, false, 0, 1000, false);

        // Values are stored as provided, even if "out of order"
        assert_eq!(info.pos1, 2000); // pos1 > pos2 is allowed
        assert_eq!(info.pos2, 1000);
    }
}
