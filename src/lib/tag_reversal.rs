//! Per-base tag reversal for negative-strand reads.
//!
//! When consensus reads are mapped to negative strand, their per-base tags need to be
//! reversed to match the orientation of the sequence in the BAM file.

use anyhow::{Result, bail};

use crate::consensus_tags::per_base;
use crate::sort::bam_fields;
use fgumi_raw_bam::RawRecordView;

/// Reverses per-base tags for a negative-strand read using raw BAM bytes.
///
/// Operates directly on the binary record without parsing into `RecordBuf`.
/// Tags in [`per_base::tags_to_reverse`] are reversed in place (B-arrays element-wise,
/// preserving each element's internal byte order; Z-strings such as `aq`/`bq` have their
/// bytes reversed). Tags in
/// [`per_base::tags_to_reverse_complement`] (`ac`, `bc`) are reverse-complemented
/// in place. See those helpers for the authoritative tag list.
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
    let flg = RawRecordView::new(record).flags();
    if (flg & bam_fields::flags::REVERSE) == 0 {
        return Ok(false);
    }

    let aux_off = bam_fields::aux_data_offset_from_record(record).unwrap_or(record.len());
    if aux_off >= record.len() {
        return Ok(true);
    }

    // Tags to reverse: cd, ce, ad, ae, bd, be, aq, bq
    // These may be B-type arrays or Z-type strings (aq, bq are Phred+33 strings)
    for tag in per_base::tags_to_reverse() {
        // Check tag type to determine reversal method
        let tag_type = bam_fields::find_tag_type(&record[aux_off..], tag);
        match tag_type {
            Some(b'B') => {
                bam_fields::reverse_array_tag_in_place(record, aux_off, tag);
            }
            Some(b'Z') => {
                bam_fields::reverse_string_tag_in_place(record, aux_off, tag);
            }
            _ => {} // Tag not found or unsupported type — skip
        }
    }

    // Tags to reverse-complement: ac, bc
    for tag in per_base::tags_to_reverse_complement() {
        bam_fields::reverse_complement_string_tag_in_place(record, aux_off, tag);
    }

    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::{
        SamBuilder as RawSamBuilder, SamTag, flags as raw_flags, raw_record_to_record_buf,
    };
    use noodles::sam::Header;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    use rstest::rstest;

    fn to_record_buf(raw: &fgumi_raw_bam::RawRecord) -> noodles::sam::alignment::RecordBuf {
        raw_record_to_record_buf(raw, &Header::default())
            .expect("raw_record_to_record_buf failed in test")
    }

    #[rstest]
    #[case(Value::from(vec![1i8, 2, 3]), Value::from(vec![3i8, 2, 1]))]
    #[case(Value::from(vec![1u8, 2, 3]), Value::from(vec![3u8, 2, 1]))]
    #[case(Value::from(vec![1i16, 2, 3]), Value::from(vec![3i16, 2, 1]))]
    #[case(Value::from(vec![1u16, 2, 3]), Value::from(vec![3u16, 2, 1]))]
    #[case(Value::from(vec![1i32, 2, 3]), Value::from(vec![3i32, 2, 1]))]
    #[case(Value::from(vec![1u32, 2, 3]), Value::from(vec![3u32, 2, 1]))]
    #[case(Value::from(vec![1.0f32, 2.0, 3.0]), Value::from(vec![3.0f32, 2.0, 1.0]))]
    fn test_reverse_buf_value_all_array_types(#[case] input: Value, #[case] expected: Value) {
        let reversed = fgumi_sam::reverse_buf_value(&input);
        assert_eq!(reversed, expected);
    }

    #[test]
    fn test_reverse_buf_value_non_array() {
        let value = Value::from(42i32);
        let reversed = fgumi_sam::reverse_buf_value(&value);
        assert_eq!(reversed, value); // Non-array values are returned unchanged
    }

    #[test]
    fn test_reverse_per_base_tags_raw_string_tag_aq() {
        // Verify aq Z-type tag is reversed via the raw path too
        use crate::sort::bam_fields;

        let mut b = RawSamBuilder::new();
        b.sequence(b"ACGT").qualities(&[30; 4]).flags(raw_flags::REVERSE);
        b.add_string_tag(SamTag::AQ, b"IIHG");
        let record_buf = to_record_buf(&b.build());

        let header = Header::default();
        let raw_rec = fgumi_raw_bam::encode_record_buf_to_raw(&record_buf, &header)
            .expect("encoding record should succeed");
        let mut raw: Vec<u8> = raw_rec.as_ref().to_vec();

        let result =
            reverse_per_base_tags_raw(&mut raw).expect("reverse_per_base_tags_raw should succeed");
        assert!(result);

        let aux = bam_fields::aux_data_slice(&raw);
        let s = bam_fields::find_string_tag(aux, SamTag::AQ).expect("aq tag should exist");
        assert_eq!(s, b"GHII");
    }

    #[rstest]
    #[case(SamTag::AC)]
    #[case(SamTag::BC_BASES)]
    fn test_reverse_per_base_tags_raw_revcomp_string_tags(#[case] tag: SamTag) {
        // Verify ac/bc Z-type tags are reverse-complemented on the raw path.
        use crate::sort::bam_fields;

        let mut b = RawSamBuilder::new();
        b.sequence(b"ACGT").qualities(&[30; 4]).flags(raw_flags::REVERSE);
        b.add_string_tag(tag, b"ACGA");
        let record_buf = to_record_buf(&b.build());

        let header = Header::default();
        let raw_rec = fgumi_raw_bam::encode_record_buf_to_raw(&record_buf, &header)
            .expect("encoding record should succeed");
        let mut raw: Vec<u8> = raw_rec.as_ref().to_vec();

        let result =
            reverse_per_base_tags_raw(&mut raw).expect("reverse_per_base_tags_raw should succeed");
        assert!(result);

        let aux = bam_fields::aux_data_slice(&raw);
        let s = bam_fields::find_string_tag(aux, tag).expect("tag should exist");
        assert_eq!(s, b"TCGT");
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
        let mut b = RawSamBuilder::new();
        b.sequence(b"ACGT").qualities(&[30; 4]).flags(0);
        let record_buf = to_record_buf(&b.build());
        let header = Header::default();
        let raw_rec = fgumi_raw_bam::encode_record_buf_to_raw(&record_buf, &header)
            .expect("encoding record should succeed");
        let mut raw: Vec<u8> = raw_rec.as_ref().to_vec();

        let result = reverse_per_base_tags_raw(&mut raw);
        assert!(result.is_ok());
        assert!(!result.expect("result should be Ok")); // positive strand => false
    }

    #[test]
    fn test_reverse_per_base_tags_raw_negative_strand() {
        // Build a valid raw record on negative strand with a cd tag — should reverse it
        use crate::sort::bam_fields;

        let mut b = RawSamBuilder::new();
        b.sequence(b"ACGT").qualities(&[30; 4]).flags(raw_flags::REVERSE);
        let mut record_buf = to_record_buf(&b.build());
        let tag = Tag::from([b'c', b'd']);
        record_buf.data_mut().insert(tag, Value::from(vec![1u16, 2, 3, 4]));

        let header = Header::default();
        let raw_rec = fgumi_raw_bam::encode_record_buf_to_raw(&record_buf, &header)
            .expect("encoding record should succeed");
        let mut raw: Vec<u8> = raw_rec.as_ref().to_vec();

        let result = reverse_per_base_tags_raw(&mut raw);
        assert!(result.is_ok());
        assert!(result.expect("result should be Ok")); // negative strand => true

        // Verify the cd tag was reversed: find it in the raw record's aux data
        let aux = bam_fields::aux_data_slice(&raw);
        let arr = bam_fields::find_array_tag(aux, SamTag::CD_BASES).expect("cd tag should exist");
        let values = bam_fields::array_tag_to_vec_u16(&arr);
        assert_eq!(values, vec![4, 3, 2, 1]);
    }
}
