use crate::sequence::pack_sequence_into;
use crate::tags::{append_string_tag, append_int_tag, append_float_tag, append_i16_array_tag, append_phred33_string_tag};

// ============================================================================
// Unmapped BAM Record Builder
// ============================================================================

/// Reusable builder for constructing unmapped BAM records as raw bytes.
///
/// Produces BAM record bytes directly, bypassing noodles `RecordBuf`
/// serialization. Designed for the consensus output hot path where
/// records are always unmapped (`ref_id=-1`, `pos=-1`, no CIGAR).
///
/// The internal buffer is reused across calls to [`Self::build_record`] +
/// [`Self::clear`] to avoid repeated allocation.
///
/// # Usage
///
/// ```rust,ignore
/// let mut builder = UnmappedBamRecordBuilder::new();
///
/// builder.build_record(b"cons:1:ACG-TCG", my_flags, &bases, &quals);
/// builder.append_string_tag(b"RG", b"sample1");
/// builder.append_int_tag(b"cD", max_depth);
/// builder.append_float_tag(b"cE", error_rate);
/// builder.append_i16_array_tag(b"cd", &depth_array);
/// builder.write_with_block_size(&mut output);
///
/// builder.clear();
/// // ... build next record ...
/// ```
pub struct UnmappedBamRecordBuilder {
    buf: Vec<u8>,
    sealed: bool,
}

/// Unmapped BAM bin (SAM spec §4.2.1: `reg2bin(-1, 0)` = 4680).
const UNMAPPED_BIN: u16 = 4680;

impl UnmappedBamRecordBuilder {
    /// Create a new builder with default capacity (512 bytes).
    #[must_use]
    pub fn new() -> Self {
        Self::with_capacity(512)
    }

    /// Create a new builder with the given initial capacity.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self { buf: Vec::with_capacity(capacity), sealed: false }
    }

    /// Clear the builder for reuse without deallocating.
    pub fn clear(&mut self) {
        self.buf.clear();
        self.sealed = false;
    }

    /// Write the 32-byte fixed header, read name, packed sequence, and
    /// quality scores for an unmapped record.
    ///
    /// After calling this, append tags via the `append_*_tag` methods,
    /// then call [`Self::write_with_block_size`] or [`Self::as_bytes`].
    ///
    /// Hardcoded unmapped fields: `ref_id=-1`, `pos=-1`, `mapq=0`,
    /// `bin=4680`, `n_cigar_op=0`, `next_ref_id=-1`, `next_pos=-1`,
    /// `tlen=0`.
    ///
    /// # Arguments
    ///
    /// * `name`  — read name **without** null terminator
    /// * `flag`  — BAM flags (`u16`)
    /// * `bases` — ASCII sequence (e.g. `b"ACGT"`)
    /// * `quals` — raw Phred quality scores (0-93, **not** Phred+33)
    ///
    /// # Panics
    ///
    /// Panics if `quals` is non-empty and `bases.len() != quals.len()`.
    pub fn build_record(&mut self, name: &[u8], flag: u16, bases: &[u8], quals: &[u8]) {
        assert!(
            bases.len() == quals.len() || quals.is_empty(),
            "bases.len() ({}) != quals.len() ({})",
            bases.len(),
            quals.len(),
        );

        self.buf.clear();
        self.sealed = false;

        assert!(name.len() < 255, "read name too long ({} bytes, max 254)", name.len());
        let l_read_name = u8::try_from(name.len() + 1).expect("read name length overflow"); // +1 for NUL
        let l_seq = u32::try_from(bases.len()).expect("sequence length overflow");
        let packed_seq_len = bases.len().div_ceil(2);

        // Pre-reserve for header + name + seq + qual + typical tag overhead
        self.buf.reserve(32 + l_read_name as usize + packed_seq_len + bases.len() + 100);

        // === Fixed 32-byte header ===
        self.buf.extend_from_slice(&(-1i32).to_le_bytes()); // ref_id = -1
        self.buf.extend_from_slice(&(-1i32).to_le_bytes()); // pos = -1
        self.buf.push(l_read_name); // l_read_name
        self.buf.push(0); // mapq = 0
        self.buf.extend_from_slice(&UNMAPPED_BIN.to_le_bytes()); // bin = 4680
        self.buf.extend_from_slice(&0u16.to_le_bytes()); // n_cigar_op = 0
        self.buf.extend_from_slice(&flag.to_le_bytes()); // flags
        self.buf.extend_from_slice(&l_seq.to_le_bytes()); // l_seq
        self.buf.extend_from_slice(&(-1i32).to_le_bytes()); // next_ref_id = -1
        self.buf.extend_from_slice(&(-1i32).to_le_bytes()); // next_pos = -1
        self.buf.extend_from_slice(&0i32.to_le_bytes()); // tlen = 0

        // === Read name + NUL ===
        self.buf.extend_from_slice(name);
        self.buf.push(0);

        // === Packed sequence ===
        pack_sequence_into(&mut self.buf, bases);

        // === Quality scores ===
        if quals.is_empty() && !bases.is_empty() {
            // Missing quality: 0xFF per SAM spec
            self.buf.resize(self.buf.len() + bases.len(), 0xFF);
        } else {
            self.buf.extend_from_slice(quals);
        }

        self.sealed = true;
    }

    /// Append a string (`Z`-type) tag.
    #[inline]
    pub fn append_string_tag(&mut self, tag: &[u8; 2], value: &[u8]) {
        debug_assert!(self.sealed, "must call build_record before appending tags");
        append_string_tag(&mut self.buf, tag, value);
    }

    /// Append an integer tag using the smallest signed type that fits.
    #[inline]
    pub fn append_int_tag(&mut self, tag: &[u8; 2], value: i32) {
        debug_assert!(self.sealed, "must call build_record before appending tags");
        append_int_tag(&mut self.buf, tag, value);
    }

    /// Append a float (`f`-type) tag.
    #[inline]
    pub fn append_float_tag(&mut self, tag: &[u8; 2], value: f32) {
        debug_assert!(self.sealed, "must call build_record before appending tags");
        append_float_tag(&mut self.buf, tag, value);
    }

    /// Append an `i16` array (`B:s`-type) tag.
    #[inline]
    pub fn append_i16_array_tag(&mut self, tag: &[u8; 2], values: &[i16]) {
        debug_assert!(self.sealed, "must call build_record before appending tags");
        append_i16_array_tag(&mut self.buf, tag, values);
    }

    /// Append a Phred+33 encoded quality string (`Z`-type) tag.
    #[inline]
    pub fn append_phred33_string_tag(&mut self, tag: &[u8; 2], quals: &[u8]) {
        debug_assert!(self.sealed, "must call build_record before appending tags");
        append_phred33_string_tag(&mut self.buf, tag, quals);
    }

    /// Get the completed record bytes (**without** the 4-byte `block_size` prefix).
    #[inline]
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        debug_assert!(self.sealed, "must call build_record first");
        &self.buf
    }

    /// Write the record with a 4-byte little-endian `block_size` prefix to `output`.
    ///
    /// This matches the output format expected by the pipeline serialize
    /// functions and the BGZF compression step.
    ///
    /// # Panics
    ///
    /// Panics if the internal buffer length exceeds `u32::MAX`.
    #[inline]
    pub fn write_with_block_size(&self, output: &mut Vec<u8>) {
        debug_assert!(self.sealed, "must call build_record first");
        let block_size = u32::try_from(self.buf.len()).expect("record too large for BAM block_size");
        output.extend_from_slice(&block_size.to_le_bytes());
        output.extend_from_slice(&self.buf);
    }
}

impl Default for UnmappedBamRecordBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
#[allow(clippy::identity_op)]
mod tests {
    use super::*;
    use crate::fields::*;
    use crate::tags::*;
    use crate::sequence::*;
    use crate::testutil::*;

    #[test]
    fn test_builder_basic_unmapped_record() {
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"read1", flags::UNMAPPED, b"ACGT", &[30, 25, 35, 40]);

        let bam = builder.as_bytes();

        // Verify fixed header fields via read primitives
        assert_eq!(ref_id(bam), -1);
        assert_eq!(pos(bam), -1);
        assert_eq!(l_read_name(bam), 6); // "read1" + NUL
        assert_eq!(mapq(bam), 0);
        assert_eq!(n_cigar_op(bam), 0);
        assert_eq!(flags(bam), flags::UNMAPPED);
        assert_eq!(l_seq(bam), 4);
        assert_eq!(mate_ref_id(bam), -1);
        assert_eq!(mate_pos(bam), -1);
        assert_eq!(read_name(bam), b"read1");

        // Verify packed sequence via read primitives
        let so = seq_offset(bam);
        assert_eq!(get_base(bam, so, 0), 1); // A
        assert_eq!(get_base(bam, so, 1), 2); // C
        assert_eq!(get_base(bam, so, 2), 4); // G
        assert_eq!(get_base(bam, so, 3), 8); // T

        // Verify quality scores via read primitives
        let qo = qual_offset(bam);
        assert_eq!(get_qual(bam, qo, 0), 30);
        assert_eq!(get_qual(bam, qo, 1), 25);
        assert_eq!(get_qual(bam, qo, 2), 35);
        assert_eq!(get_qual(bam, qo, 3), 40);

        // Verify no aux data
        assert!(aux_data_slice(bam).is_empty());
    }

    #[test]
    fn test_builder_with_tags() {
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"cons:1:ACG", flags::UNMAPPED, b"ACGT", &[30, 30, 30, 30]);
        builder.append_string_tag(b"RG", b"sample1");
        builder.append_string_tag(b"MI", b"42");
        builder.append_int_tag(b"cD", 10);
        builder.append_int_tag(b"cM", 3);
        builder.append_float_tag(b"cE", 0.05);
        builder.append_i16_array_tag(b"cd", &[10, 8, 12, 9]);

        let bam = builder.as_bytes();
        let aux = aux_data_slice(bam);

        // Verify string tags
        assert_eq!(find_string_tag(aux, b"RG"), Some(b"sample1" as &[u8]));
        assert_eq!(find_string_tag(aux, b"MI"), Some(b"42" as &[u8]));

        // Verify integer tags
        assert_eq!(find_int_tag(aux, b"cD"), Some(10));
        assert_eq!(find_int_tag(aux, b"cM"), Some(3));

        // Verify float tag
        let ce = find_float_tag(aux, b"cE").unwrap();
        assert!((ce - 0.05).abs() < 1e-7);
    }

    #[test]
    fn test_builder_reuse() {
        let mut builder = UnmappedBamRecordBuilder::new();

        // Build first record
        builder.build_record(b"r1", flags::UNMAPPED, b"ACGT", &[30, 30, 30, 30]);
        builder.append_string_tag(b"MI", b"1");
        let len1 = builder.as_bytes().len();
        assert_eq!(l_seq(builder.as_bytes()), 4);
        assert_eq!(read_name(builder.as_bytes()), b"r1");

        // Reuse for different-length record
        builder.clear();
        builder.build_record(b"r2", flags::UNMAPPED, b"AC", &[30, 30]);
        builder.append_string_tag(b"MI", b"2");

        let bam = builder.as_bytes();
        assert_ne!(bam.len(), len1);
        assert_eq!(l_seq(bam), 2);
        assert_eq!(read_name(bam), b"r2");
        assert_eq!(find_string_tag(aux_data_slice(bam), b"MI"), Some(b"2" as &[u8]));
    }

    #[test]
    fn test_builder_write_with_block_size() {
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"r1", flags::UNMAPPED, b"ACGT", &[30, 30, 30, 30]);

        let mut output = Vec::new();
        builder.write_with_block_size(&mut output);

        // First 4 bytes = block_size
        let block_size = u32::from_le_bytes([output[0], output[1], output[2], output[3]]);
        assert_eq!(block_size as usize, builder.as_bytes().len());
        assert_eq!(&output[4..], builder.as_bytes());
    }

    #[test]
    fn test_builder_paired_consensus_flags() {
        let mut builder = UnmappedBamRecordBuilder::new();

        // R1
        let r1_flags =
            flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | flags::FIRST_SEGMENT;
        builder.build_record(b"cons:1:ACG", r1_flags, b"ACGT", &[30, 30, 30, 30]);
        assert_eq!(flags(builder.as_bytes()), r1_flags);

        // R2
        builder.clear();
        let r2_flags = flags::PAIRED | flags::UNMAPPED | flags::MATE_UNMAPPED | flags::LAST_SEGMENT;
        builder.build_record(b"cons:1:ACG", r2_flags, b"TGCA", &[25, 35, 30, 40]);
        assert_eq!(flags(builder.as_bytes()), r2_flags);
    }

    #[test]
    fn test_builder_empty_sequence() {
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"empty", flags::UNMAPPED, b"", &[]);
        builder.append_string_tag(b"RG", b"rg0");

        let bam = builder.as_bytes();
        assert_eq!(l_seq(bam), 0);
        assert_eq!(find_string_tag(aux_data_slice(bam), b"RG"), Some(b"rg0" as &[u8]));
    }

    #[test]
    fn test_builder_missing_quality() {
        let mut builder = UnmappedBamRecordBuilder::new();
        // Empty quals with non-empty bases -> fills with 0xFF
        builder.build_record(b"r1", flags::UNMAPPED, b"ACGT", &[]);

        let bam = builder.as_bytes();
        let qo = qual_offset(bam);
        for i in 0..4 {
            assert_eq!(get_qual(bam, qo, i), 0xFF);
        }
    }

    #[test]
    fn test_builder_header_matches_vendored_default() {
        // Cross-validate header layout against known bytes from
        // vendored encoder test_encode_with_default_fields.
        // Default RecordBuf: name="*", flags=UNMAPPED, mapq=255, seq=empty
        // Our builder: custom name, flags, mapq=0 (hardcoded for consensus)
        //
        // Verify shared unmapped constants: ref_id=-1, pos=-1, bin=4680,
        // n_cigar_op=0, next_ref_id=-1, next_pos=-1, tlen=0
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"*", flags::UNMAPPED, b"", &[]);

        let bam = builder.as_bytes();
        // ref_id = -1
        assert_eq!(&bam[0..4], &0xFF_FF_FF_FFu32.to_le_bytes());
        // pos = -1
        assert_eq!(&bam[4..8], &0xFF_FF_FF_FFu32.to_le_bytes());
        // l_read_name = 2 ("*" + NUL)
        assert_eq!(bam[8], 0x02);
        // mapq = 0 (consensus-specific, vendored default is 255)
        assert_eq!(bam[9], 0x00);
        // bin = 4680 = 0x1248
        assert_eq!(&bam[10..12], &[0x48, 0x12]);
        // n_cigar_op = 0
        assert_eq!(&bam[12..14], &[0x00, 0x00]);
        // flags = UNMAPPED = 4
        assert_eq!(&bam[14..16], &[0x04, 0x00]);
        // l_seq = 0
        assert_eq!(&bam[16..20], &[0x00, 0x00, 0x00, 0x00]);
        // next_ref_id = -1
        assert_eq!(&bam[20..24], &0xFF_FF_FF_FFu32.to_le_bytes());
        // next_pos = -1
        assert_eq!(&bam[24..28], &0xFF_FF_FF_FFu32.to_le_bytes());
        // tlen = 0
        assert_eq!(&bam[28..32], &[0x00, 0x00, 0x00, 0x00]);
        // name = "*\0"
        assert_eq!(&bam[32..34], &[b'*', 0x00]);
    }

    #[test]
    fn test_builder_phred33_string_tag() {
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"test", flags::UNMAPPED, b"ACGT", &[30, 30, 30, 30]);
        builder.append_phred33_string_tag(b"aq", &[0, 10, 30, 40]);
        let bytes = builder.as_bytes();
        let parsed = ParsedBamRecord::from_bytes(bytes);
        let aq = parsed.get_string_tag(b"aq").unwrap();
        assert_eq!(aq, b"!+?I");
    }

    #[test]
    fn test_parsed_bam_record_roundtrip() {
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(
            b"myread",
            flags::UNMAPPED | flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_UNMAPPED,
            b"ACGTNN",
            &[30, 25, 20, 15, 10, 5],
        );
        builder.append_string_tag(b"RG", b"sample1");
        builder.append_int_tag(b"cD", 5);
        builder.append_float_tag(b"cE", 0.01);
        builder.append_i16_array_tag(b"cd", &[3, 4, 5, 3, 2, 1]);

        let mut output = Vec::new();
        builder.write_with_block_size(&mut output);

        let records = ParsedBamRecord::parse_all(&output);
        assert_eq!(records.len(), 1);
        let rec = &records[0];
        assert_eq!(rec.name, b"myread");
        assert_eq!(
            rec.flag,
            flags::UNMAPPED | flags::PAIRED | flags::FIRST_SEGMENT | flags::MATE_UNMAPPED
        );
        assert_eq!(rec.bases, b"ACGTNN");
        assert_eq!(rec.quals, vec![30, 25, 20, 15, 10, 5]);
        assert_eq!(rec.get_string_tag(b"RG").unwrap(), b"sample1");
        assert_eq!(rec.get_int_tag(b"cD").unwrap(), 5);
        let ce = rec.get_float_tag(b"cE").unwrap();
        assert!((ce - 0.01).abs() < 0.001);
        assert_eq!(rec.get_i16_array_tag(b"cd").unwrap(), vec![3, 4, 5, 3, 2, 1]);
    }

    // ========================================================================
    // UnmappedBamRecordBuilder additional tests
    // ========================================================================

    #[test]
    fn test_builder_default() {
        let builder = UnmappedBamRecordBuilder::default();
        // Default builder has empty buffer and is not sealed
        assert!(builder.buf.is_empty());
        assert!(!builder.sealed);
    }

    #[test]
    fn test_builder_with_capacity() {
        let builder = UnmappedBamRecordBuilder::with_capacity(1024);
        assert!(builder.buf.capacity() >= 1024);
    }

    #[test]
    fn test_builder_odd_length_sequence() {
        // Odd-length sequence to exercise the packing edge case
        let mut builder = UnmappedBamRecordBuilder::new();
        builder.build_record(b"read1", flags::UNMAPPED, b"ACG", &[30, 25, 20]);

        let bam = builder.as_bytes();
        assert_eq!(l_seq(bam), 3);
        let so = seq_offset(bam);
        assert_eq!(get_base(bam, so, 0), 1); // A
        assert_eq!(get_base(bam, so, 1), 2); // C
        assert_eq!(get_base(bam, so, 2), 4); // G
    }

    #[test]
    fn test_builder_long_sequence() {
        // Longer sequence to exercise multi-byte packing
        let mut builder = UnmappedBamRecordBuilder::new();
        let bases = b"ACGTACGTACGTACGT"; // 16 bases
        let quals = vec![30u8; 16];
        builder.build_record(b"r1", flags::UNMAPPED, bases, &quals);

        let bam = builder.as_bytes();
        assert_eq!(l_seq(bam), 16);
        let so = seq_offset(bam);
        for (i, &base) in bases.iter().enumerate() {
            let expected = match base {
                b'A' => 1,
                b'C' => 2,
                b'G' => 4,
                b'T' => 8,
                _ => 0,
            };
            assert_eq!(get_base(bam, so, i), expected, "base mismatch at position {i}");
        }
    }

    #[test]
    fn test_builder_multiple_records_via_write_with_block_size() {
        // Build two records into the same output buffer
        let mut builder = UnmappedBamRecordBuilder::new();
        let mut output = Vec::new();

        builder.build_record(b"r1", flags::UNMAPPED, b"ACGT", &[30, 30, 30, 30]);
        builder.append_string_tag(b"MI", b"1");
        builder.write_with_block_size(&mut output);

        builder.clear();
        builder.build_record(b"r2", flags::UNMAPPED, b"TGCA", &[25, 25, 25, 25]);
        builder.append_string_tag(b"MI", b"2");
        builder.write_with_block_size(&mut output);

        let records = ParsedBamRecord::parse_all(&output);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].name, b"r1");
        assert_eq!(records[0].bases, b"ACGT");
        assert_eq!(records[0].get_string_tag(b"MI").unwrap(), b"1");
        assert_eq!(records[1].name, b"r2");
        assert_eq!(records[1].bases, b"TGCA");
        assert_eq!(records[1].get_string_tag(b"MI").unwrap(), b"2");
    }
}
