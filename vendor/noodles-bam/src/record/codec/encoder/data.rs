//! BAM record data field writer.

pub mod field;

use std::io;

use noodles_sam::alignment::record::{Data, data::field::Tag};
use noodles_sam::alignment::record_buf::data::field::Value as RecordBufValue;

use self::field::write_field;
use super::num::{write_f32_le, write_i16_le, write_u32_le};

pub(super) fn write_data<D>(dst: &mut Vec<u8>, data: D) -> io::Result<()>
where
    D: Data,
{
    for result in data.iter() {
        let (tag, value) = result?;

        if &tag == Tag::CIGAR.as_ref() {
            continue;
        }

        write_field(dst, tag, &value)?;
    }

    Ok(())
}

/// Fast path for writing data fields from a RecordBuf.
///
/// This function provides optimized encoding for common tags:
/// - Standard SAM tags: NH, NM, AS, RG
/// - UMI tags: RX, MI, QX, OX
/// - Consensus per-read tags: cD, cM, cE, aD, bD, aM, bM, aE, bE
/// - Consensus per-base arrays: cd, ce, ad, ae, bd, be (bulk memcpy for Int16)
/// - Consensus strings: ac, bc, aq, bq
pub(super) fn write_data_fast(
    dst: &mut Vec<u8>,
    data: &noodles_sam::alignment::record_buf::Data,
) -> io::Result<()> {
    for (tag, value) in data.iter() {
        // Skip CIGAR tag (handled separately for oversized CIGARs)
        if tag == Tag::CIGAR {
            continue;
        }

        // Try fast path first
        if !try_write_common_tag(dst, &tag, value) {
            // Fall back to generic path
            write_field_from_value(dst, tag, value)?;
        }
    }

    Ok(())
}

/// Try to write a common tag using the fast path.
/// Returns true if handled, false if caller should use generic path.
#[inline]
fn try_write_common_tag(dst: &mut Vec<u8>, tag: &Tag, value: &RecordBufValue) -> bool {
    let tag_bytes = tag.as_ref();

    match (tag_bytes, value) {
        // =====================================================================
        // Standard SAM tags
        // =====================================================================

        // NH:C (alignment hit count) - UInt8
        ([b'N', b'H'], RecordBufValue::UInt8(n)) => {
            dst.extend([b'N', b'H', b'C', *n]);
            true
        }

        // NM:i (edit distance) - Int32
        ([b'N', b'M'], RecordBufValue::Int32(n)) => {
            dst.extend([b'N', b'M', b'i']);
            dst.extend(n.to_le_bytes());
            true
        }

        // AS:i (alignment score) - Int32
        ([b'A', b'S'], RecordBufValue::Int32(n)) => {
            dst.extend([b'A', b'S', b'i']);
            dst.extend(n.to_le_bytes());
            true
        }

        // =====================================================================
        // UMI and string tags (bulk copy)
        // =====================================================================

        // RX, MI, QX, OX, RG - String tags with bulk copy
        ([b'R', b'X'], RecordBufValue::String(s))
        | ([b'M', b'I'], RecordBufValue::String(s))
        | ([b'Q', b'X'], RecordBufValue::String(s))
        | ([b'O', b'X'], RecordBufValue::String(s))
        | ([b'B', b'Z'], RecordBufValue::String(s))
        | ([b'R', b'G'], RecordBufValue::String(s)) => {
            dst.extend([tag_bytes[0], tag_bytes[1], b'Z']);
            dst.extend_from_slice(s.as_ref());
            dst.push(0); // NUL terminator
            true
        }

        // =====================================================================
        // Consensus per-read tags
        // =====================================================================

        // cD, cM, aD, bD, aM, bM - Int16 (depth/min depth metrics)
        ([b'c' | b'a' | b'b', b'D' | b'M'], RecordBufValue::Int16(n)) => {
            dst.extend([tag_bytes[0], tag_bytes[1], b's']);
            write_i16_le(dst, *n);
            true
        }

        // cE, aE, bE - Float (error rate)
        ([b'c' | b'a' | b'b', b'E'], RecordBufValue::Float(n)) => {
            dst.extend([tag_bytes[0], tag_bytes[1], b'f']);
            write_f32_le(dst, *n);
            true
        }

        // =====================================================================
        // Consensus per-base arrays (BIGGEST OPPORTUNITY)
        // =====================================================================

        // cd, ce, ad, ae, bd, be - Int16 arrays (per-base depths/errors)
        ([b'c' | b'a' | b'b', b'd' | b'e'], RecordBufValue::Array(arr)) => {
            if let noodles_sam::alignment::record_buf::data::field::value::Array::Int16(values) =
                arr
            {
                // Write header: tag + 'B' + 's' + count
                dst.extend([tag_bytes[0], tag_bytes[1], b'B', b's']);
                write_u32_le(dst, values.len() as u32);

                // Bulk copy Int16 array (works on little-endian systems)
                // SAFETY: i16 has no padding and same byte representation on little-endian
                #[cfg(target_endian = "little")]
                {
                    let bytes: &[u8] = unsafe {
                        std::slice::from_raw_parts(
                            values.as_ptr() as *const u8,
                            values.len() * std::mem::size_of::<i16>(),
                        )
                    };
                    dst.extend_from_slice(bytes);
                }

                // Fallback for big-endian (rare)
                #[cfg(target_endian = "big")]
                {
                    for &v in values.iter() {
                        write_i16_le(dst, v);
                    }
                }

                true
            } else {
                false
            }
        }

        // =====================================================================
        // Consensus string tags (ac, bc, aq, bq)
        // =====================================================================
        ([b'a' | b'b', b'c' | b'q'], RecordBufValue::String(s)) => {
            dst.extend([tag_bytes[0], tag_bytes[1], b'Z']);
            dst.extend_from_slice(s.as_ref());
            dst.push(0);
            true
        }

        // Not a common tag - use generic path
        _ => false,
    }
}

/// Write a field from a RecordBufValue (used in fast path fallback).
#[inline]
fn write_field_from_value(
    dst: &mut Vec<u8>,
    tag: Tag,
    value: &RecordBufValue,
) -> io::Result<()> {
    use noodles_sam::alignment::record::data::field::Value;

    // Convert RecordBufValue to Value reference for write_field
    // This is a bit awkward but necessary for compatibility
    let value_ref: Value<'_> = value.into();
    write_field(dst, tag, &value_ref)
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record_buf::{Data as DataBuf, data::field::Value};

    use super::*;

    #[test]
    fn test_write_data() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, data: &DataBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_data(buf, data)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let data = DataBuf::default();
        t(&mut buf, &data, &[])?;

        let data = [(Tag::ALIGNMENT_HIT_COUNT, Value::from(1))]
            .into_iter()
            .collect();
        t(&mut buf, &data, &[b'N', b'H', b'C', 0x01])?;

        let data = [
            (Tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
            (Tag::READ_GROUP, Value::from("rg0")),
        ]
        .into_iter()
        .collect();
        t(
            &mut buf,
            &data,
            &[
                b'N', b'H', b'C', 0x01, // NH:C:1
                b'R', b'G', b'Z', b'r', b'g', b'0', 0x00, // RG:Z:rg0
            ],
        )?;

        Ok(())
    }
}
