//! Per-base tag reversal for negative-strand reads.
//!
//! When consensus reads are mapped to negative strand, their per-base tags need to be
//! reversed to match the orientation of the sequence in the BAM file.

use anyhow::{Result, bail};
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value;

use crate::consensus_tags::per_base;
use crate::dna::reverse_complement;
use crate::sort::bam_fields;

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
///
/// # Errors
///
/// Currently infallible; returns `Result` for API consistency with
/// [`reverse_per_base_tags_raw`].
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
            match value {
                Value::String(s) => {
                    // Z-type string tags (e.g. aq, bq quality strings) — reverse bytes
                    let mut reversed: Vec<u8> = s.iter().copied().collect();
                    reversed.reverse();
                    record
                        .data_mut()
                        .insert(tag, Value::from(String::from_utf8_lossy(&reversed).to_string()));
                }
                _ => {
                    if let Some(reversed_value) = reverse_array_value(value) {
                        record.data_mut().insert(tag, reversed_value);
                    }
                }
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

/// Reverses per-base tags for a negative-strand read using raw BAM bytes.
///
/// Operates directly on the binary record without parsing into `RecordBuf`.
/// Array tags (cd, ce, ad, ae, bd, be, aq, bq) are reversed in-place.
/// String tags (ac, bc) are reverse-complemented in-place.
///
/// # Returns
/// `Ok(true)` if tags were reversed, `Ok(false)` if not on reverse strand.
///
/// # Errors
///
/// Returns an error if the record is too short to be a valid BAM record.
pub fn reverse_per_base_tags_raw(record: &mut [u8]) -> Result<bool> {
    if record.len() < bam_fields::MIN_BAM_RECORD_LEN {
        bail!(
            "BAM record too short ({} bytes, minimum {})",
            record.len(),
            bam_fields::MIN_BAM_RECORD_LEN
        );
    }
    let flg = bam_fields::flags(record);
    if (flg & bam_fields::flags::REVERSE) == 0 {
        return Ok(false);
    }

    let aux_off = bam_fields::aux_data_offset_from_record(record).unwrap_or(record.len());
    if aux_off >= record.len() {
        return Ok(true);
    }

    // Tags to reverse: cd, ce, ad, ae, bd, be, aq, bq
    // These may be B-type arrays or Z-type strings (aq, bq are Phred+33 strings)
    for tag_str in per_base::tags_to_reverse() {
        let tag_bytes: [u8; 2] = [tag_str.as_bytes()[0], tag_str.as_bytes()[1]];
        // Check tag type to determine reversal method
        let tag_type = bam_fields::find_tag_type(&record[aux_off..], &tag_bytes);
        match tag_type {
            Some(b'B') => {
                bam_fields::reverse_array_tag_in_place(record, aux_off, &tag_bytes);
            }
            Some(b'Z') => {
                bam_fields::reverse_string_tag_in_place(record, aux_off, &tag_bytes);
            }
            _ => {} // Tag not found or unsupported type — skip
        }
    }

    // Tags to reverse-complement: ac, bc
    for tag_str in per_base::tags_to_reverse_complement() {
        let tag_bytes: [u8; 2] = [tag_str.as_bytes()[0], tag_str.as_bytes()[1]];
        bam_fields::reverse_complement_string_tag_in_place(record, aux_off, &tag_bytes);
    }

    Ok(true)
}

/// Reverses an array value
fn reverse_array_value(value: &Value) -> Option<Value> {
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

            Some(reversed)
        }
        _ => None,
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
            let reversed = reverse_array_value(&input).unwrap();
            assert_eq!(reversed, expected);
        }
    }

    #[test]
    fn test_reverse_array_value_non_array() {
        let value = Value::from(42i32);
        let reversed = reverse_array_value(&value);
        assert!(reversed.is_none());
    }

    #[test]
    fn test_reverse_string_tag_aq() {
        // aq is a Z-type quality string — should be reversed on negative strand
        let mut record = RecordBuilder::new()
            .sequence("ACGT")
            .reverse_complement(true)
            .tag("aq", "IIHG")
            .build();

        let reversed = reverse_per_base_tags(&mut record).unwrap();
        assert!(reversed);

        let tag = Tag::from([b'a', b'q']);
        if let Some(Value::String(s)) = record.data().get(&tag) {
            let bytes: Vec<u8> = s.iter().copied().collect();
            assert_eq!(bytes, b"GHII");
        } else {
            panic!("Expected aq tag to be present as String");
        }
    }

    #[test]
    fn test_reverse_per_base_tags_raw_string_tag_aq() {
        // Verify aq Z-type tag is reversed via the raw path too
        use crate::sort::bam_fields;
        use crate::vendored::bam_codec::encoder::encode_record_buf;
        use noodles::sam::Header;

        let record_buf = RecordBuilder::new()
            .sequence("ACGT")
            .reverse_complement(true)
            .tag("aq", "IIHG")
            .build();

        let header = Header::default();
        let mut raw = Vec::new();
        encode_record_buf(&mut raw, &header, &record_buf).unwrap();

        let result = reverse_per_base_tags_raw(&mut raw).unwrap();
        assert!(result);

        let aux = bam_fields::aux_data_slice(&raw);
        let s = bam_fields::find_string_tag(aux, b"aq").expect("aq tag should exist");
        assert_eq!(s, b"GHII");
    }

    #[test]
    fn test_reverse_per_base_tags_raw_short_record() {
        // A record shorter than MIN_BAM_RECORD_LEN (32 bytes) should return an error
        let mut short_record = vec![0u8; 16];
        let result = reverse_per_base_tags_raw(&mut short_record);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("too short"));
    }

    #[test]
    fn test_reverse_per_base_tags_raw_positive_strand() {
        // Build a valid raw record on positive strand — should return Ok(false)
        use crate::vendored::bam_codec::encoder::encode_record_buf;
        use noodles::sam::Header;

        let record_buf = RecordBuilder::new().sequence("ACGT").build();
        let header = Header::default();
        let mut raw = Vec::new();
        encode_record_buf(&mut raw, &header, &record_buf).unwrap();

        let result = reverse_per_base_tags_raw(&mut raw);
        assert!(result.is_ok());
        assert!(!result.unwrap()); // positive strand => false
    }

    #[test]
    fn test_reverse_per_base_tags_raw_negative_strand() {
        // Build a valid raw record on negative strand with a cd tag — should reverse it
        use crate::sort::bam_fields;
        use crate::vendored::bam_codec::encoder::encode_record_buf;
        use noodles::sam::Header;

        let mut record_buf = RecordBuilder::new().sequence("ACGT").reverse_complement(true).build();
        let tag = Tag::from([b'c', b'd']);
        record_buf.data_mut().insert(tag, Value::from(vec![1u16, 2, 3, 4]));

        let header = Header::default();
        let mut raw = Vec::new();
        encode_record_buf(&mut raw, &header, &record_buf).unwrap();

        let result = reverse_per_base_tags_raw(&mut raw);
        assert!(result.is_ok());
        assert!(result.unwrap()); // negative strand => true

        // Verify the cd tag was reversed: find it in the raw record's aux data
        let aux = bam_fields::aux_data_slice(&raw);
        let arr = bam_fields::find_array_tag(aux, b"cd").expect("cd tag should exist");
        let values = bam_fields::array_tag_to_vec_u16(&arr);
        assert_eq!(values, vec![4, 3, 2, 1]);
    }
}
