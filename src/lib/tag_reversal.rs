//! Per-base tag reversal for negative-strand reads.
//!
//! When consensus reads are mapped to negative strand, their per-base tags need to be
//! reversed to match the orientation of the sequence in the BAM file.

use anyhow::Result;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value;

use crate::consensus_tags::per_base;
use crate::dna::reverse_complement;

/// Reverses per-base tags for a negative-strand read
///
/// This function reverses array-based per-base tags (cd, ce, ad, bd, ae, be, aq, bq)
/// to match the reversed sequence orientation in negative-strand reads.
///
/// # Arguments
/// * `record` - The record to reverse tags for (modified in place)
///
/// # Returns
/// True if tags were reversed, false if not needed
pub fn reverse_per_base_tags(record: &mut RecordBuf) -> Result<bool> {
    // Check if read is mapped to negative strand
    let flags = record.flags();
    if !flags.is_reverse_complemented() {
        return Ok(false);
    }

    // Get list of tags to reverse
    let tags_to_reverse = per_base::tags_to_reverse();

    for tag_str in tags_to_reverse {
        let tag = per_base::tag(tag_str);

        if let Some(value) = record.data().get(&tag) {
            if let Some(reversed_value) = reverse_array_value(value)? {
                record.data_mut().insert(tag, reversed_value);
            }
        }
    }

    // Handle tags that need reverse complement (ac, bc)
    let tags_to_revcomp = per_base::tags_to_reverse_complement();

    for tag_str in tags_to_revcomp {
        let tag = per_base::tag(tag_str);

        if let Some(Value::String(bases)) = record.data().get(&tag) {
            // Convert BString to Vec<u8>
            let bases_vec: Vec<u8> = bases.iter().copied().collect();
            let revcomp = reverse_complement(&bases_vec);
            record
                .data_mut()
                .insert(tag, Value::from(String::from_utf8_lossy(&revcomp).to_string()));
        }
    }

    Ok(true)
}

/// Reverses an array value
fn reverse_array_value(value: &Value) -> Result<Option<Value>> {
    match value {
        Value::Array(arr) => {
            use noodles::sam::alignment::record_buf::data::field::value::Array;

            let reversed = match arr {
                Array::Int8(values) => {
                    let mut reversed: Vec<i8> = values.clone();
                    reversed.reverse();
                    Value::from(reversed)
                }
                Array::UInt8(values) => {
                    let mut reversed: Vec<u8> = values.clone();
                    reversed.reverse();
                    Value::from(reversed)
                }
                Array::Int16(values) => {
                    let mut reversed: Vec<i16> = values.clone();
                    reversed.reverse();
                    Value::from(reversed)
                }
                Array::UInt16(values) => {
                    let mut reversed: Vec<u16> = values.clone();
                    reversed.reverse();
                    Value::from(reversed)
                }
                Array::Int32(values) => {
                    let mut reversed: Vec<i32> = values.clone();
                    reversed.reverse();
                    Value::from(reversed)
                }
                Array::UInt32(values) => {
                    let mut reversed: Vec<u32> = values.clone();
                    reversed.reverse();
                    Value::from(reversed)
                }
                Array::Float(values) => {
                    let mut reversed: Vec<f32> = values.clone();
                    reversed.reverse();
                    Value::from(reversed)
                }
            };

            Ok(Some(reversed))
        }
        _ => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::builder::RecordBuilder;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::value::Array;

    #[test]
    fn test_reverse_per_base_tags_positive_strand() {
        let mut record = RecordBuilder::new().sequence("ACGT").build();

        let tag = Tag::from([b'c', b'd']);
        record.data_mut().insert(tag, Value::from(vec![1u16, 2, 3]));

        let reversed = reverse_per_base_tags(&mut record).unwrap();
        assert!(!reversed); // Should not reverse for positive strand

        if let Some(Value::Array(Array::UInt16(arr))) = record.data().get(&tag) {
            assert_eq!(arr.clone(), vec![1, 2, 3]);
        }
    }

    #[test]
    fn test_reverse_per_base_tags_negative_strand() {
        let mut record = RecordBuilder::new().sequence("ACGT").reverse_complement(true).build();

        let tag = Tag::from([b'c', b'd']);
        record.data_mut().insert(tag, Value::from(vec![1u16, 2, 3]));

        let reversed = reverse_per_base_tags(&mut record).unwrap();
        assert!(reversed); // Should reverse for negative strand

        if let Some(Value::Array(Array::UInt16(arr))) = record.data().get(&tag) {
            assert_eq!(arr.clone(), vec![3, 2, 1]);
        }
    }

    #[test]
    fn test_reverse_complement_tag() {
        let mut record = RecordBuilder::new()
            .sequence("ACGT")
            .reverse_complement(true)
            .tag("ac", "ACGT")
            .build();

        let reversed = reverse_per_base_tags(&mut record).unwrap();
        assert!(reversed);

        let tag = Tag::from([b'a', b'c']);
        if let Some(Value::String(bases)) = record.data().get(&tag) {
            let bases_vec: Vec<u8> = bases.iter().copied().collect();
            assert_eq!(bases_vec, b"ACGT");
        }
    }

    #[test]
    fn test_reverse_array_value_all_types() {
        // Test all array types with same pattern
        let test_cases: Vec<(Value, Value)> = vec![
            (Value::from(vec![1i8, 2, 3]), Value::from(vec![3i8, 2, 1])),
            (Value::from(vec![1u8, 2, 3]), Value::from(vec![3u8, 2, 1])),
            (Value::from(vec![1i16, 2, 3]), Value::from(vec![3i16, 2, 1])),
            (Value::from(vec![1u16, 2, 3]), Value::from(vec![3u16, 2, 1])),
            (Value::from(vec![1i32, 2, 3]), Value::from(vec![3i32, 2, 1])),
            (Value::from(vec![1u32, 2, 3]), Value::from(vec![3u32, 2, 1])),
            (Value::from(vec![1.0f32, 2.0, 3.0]), Value::from(vec![3.0f32, 2.0, 1.0])),
        ];

        for (input, expected) in test_cases {
            let reversed = reverse_array_value(&input).unwrap().unwrap();
            assert_eq!(reversed, expected);
        }
    }

    #[test]
    fn test_reverse_array_value_non_array() {
        let value = Value::from(42i32);
        let reversed = reverse_array_value(&value).unwrap();
        assert!(reversed.is_none());
    }
}
