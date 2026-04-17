use crate::raw_bam_record::RawRecord;
use crate::sequence::pack_sequence_into;
use crate::tags::{
    append_f32_array_tag, append_float_tag, append_i16_array_tag, append_i32_array_tag,
    append_int_tag, append_phred33_string_tag, append_string_tag, append_u8_array_tag,
    append_u16_array_tag,
};

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
/// let mut builder = UnmappedSamBuilder::new();
///
/// builder.build_record(b"cons:1:ACG-TCG", my_flags, &bases, &quals);
/// builder.append_string_tag(b"RG", b"sample1");
/// builder.append_int_tag(b"cD", max_depth);
/// builder.append_float_tag(b"cE", error_rate);
/// builder.append_i16_array_tag(b"cd", &depth_array);
///
/// // Return the record as a RawRecord (moves the buffer out).
/// let record: RawRecord = builder.build();
///
/// // -- or write directly without taking ownership --
/// builder.write_with_block_size(&mut output);
///
/// builder.clear();
/// // ... build next record ...
/// ```
pub struct UnmappedSamBuilder {
    buf: Vec<u8>,
    sealed: bool,
}

/// Unmapped BAM bin (SAM spec §4.2.1: `reg2bin(-1, 0)` = 4680).
const UNMAPPED_BIN: u16 = 4680;

impl UnmappedSamBuilder {
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
        #[expect(
            clippy::cast_possible_truncation,
            reason = "assert above guarantees name.len() < 255"
        )]
        let l_read_name = (name.len() + 1) as u8; // +1 for NUL
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

    /// Append an integer tag using the smallest type that fits the value.
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

    /// Append a `u8` array (`B:C`-type) tag.
    #[inline]
    pub fn append_u8_array_tag(&mut self, tag: &[u8; 2], values: &[u8]) {
        debug_assert!(self.sealed, "must call build_record before appending tags");
        append_u8_array_tag(&mut self.buf, tag, values);
    }

    /// Append a Phred+33 encoded quality string (`Z`-type) tag.
    #[inline]
    pub fn append_phred33_string_tag(&mut self, tag: &[u8; 2], quals: &[u8]) {
        debug_assert!(self.sealed, "must call build_record before appending tags");
        append_phred33_string_tag(&mut self.buf, tag, quals);
    }

    /// Consume the completed record bytes and return them as a [`RawRecord`].
    ///
    /// This moves the internal buffer out of the builder. After calling
    /// `build()`, the builder's buffer is empty; call [`Self::build_record`]
    /// again before appending tags or writing.
    ///
    /// Callers that need the raw bytes can use [`RawRecord::into_inner`].
    #[inline]
    #[must_use]
    pub fn build(&mut self) -> RawRecord {
        debug_assert!(self.sealed, "must call build_record first");
        self.sealed = false;
        RawRecord::from(std::mem::take(&mut self.buf))
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
        let block_size =
            u32::try_from(self.buf.len()).expect("record too large for BAM block_size");
        output.extend_from_slice(&block_size.to_le_bytes());
        output.extend_from_slice(&self.buf);
    }
}

impl Default for UnmappedSamBuilder {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Mapped BAM Record Builder
// ============================================================================

/// Builder for constructing mapped (or unmapped) BAM records as raw bytes.
///
/// Unlike [`UnmappedSamBuilder`], this builder supports any combination
/// of `ref_id`, `pos`, CIGAR operations, `mapq`, flags, mate fields, sequence,
/// quality scores, and auxiliary tags. It uses a fluent setter API; call
/// [`Self::build`] to assemble the final [`RawRecord`].
///
/// The internal state is reused across calls via [`Self::clear`] to avoid
/// repeated allocations.
///
/// # Default field values
///
/// | Field           | Default  |
/// |-----------------|----------|
/// | `ref_id`        | `-1`     |
/// | `pos`           | `-1`     |
/// | `mapq`          | `0`      |
/// | `flags`         | `0`      |
/// | `mate_ref_id`   | `-1`     |
/// | `mate_pos`      | `-1`     |
/// | `template_length` | `0`  |
/// | `bin`           | `0`      |
/// | `read_name`     | `b"*"`   |
/// | `cigar_ops`     | empty    |
/// | `sequence`      | empty    |
/// | `qualities`     | empty    |
///
/// # Usage
///
/// ```rust,ignore
/// let mut b = SamBuilder::new();
/// b.ref_id(0)
///  .pos(999)
///  .mapq(60)
///  .flags(flags::PAIRED | flags::FIRST_SEGMENT)
///  .read_name(b"read1")
///  .cigar_ops(&[/* 100M */ (100u32 << 4) | 0])
///  .sequence(b"ACGT")
///  .qualities(&[30, 25, 35, 40]);
/// b.add_string_tag(b"RG", b"sample1");
/// let record: RawRecord = b.build();
/// ```
pub struct SamBuilder {
    ref_id: i32,
    pos: i32,
    mapq: u8,
    flags: u16,
    mate_ref_id: i32,
    mate_pos: i32,
    tlen: i32,
    bin: u16,
    read_name: Vec<u8>,
    cigar_ops: Vec<u32>,
    sequence: Vec<u8>,
    qualities: Vec<u8>,
    aux: Vec<u8>,
}

impl SamBuilder {
    /// Create a new builder with default field values.
    #[must_use]
    pub fn new() -> Self {
        Self {
            ref_id: -1,
            pos: -1,
            mapq: 0,
            flags: 0,
            mate_ref_id: -1,
            mate_pos: -1,
            tlen: 0,
            bin: 0,
            read_name: b"*".to_vec(),
            cigar_ops: Vec::new(),
            sequence: Vec::new(),
            qualities: Vec::new(),
            aux: Vec::new(),
        }
    }

    /// Reset all fields to their defaults, retaining allocated capacity.
    pub fn clear(&mut self) {
        self.ref_id = -1;
        self.pos = -1;
        self.mapq = 0;
        self.flags = 0;
        self.mate_ref_id = -1;
        self.mate_pos = -1;
        self.tlen = 0;
        self.bin = 0;
        self.read_name.clear();
        self.read_name.extend_from_slice(b"*");
        self.cigar_ops.clear();
        self.sequence.clear();
        self.qualities.clear();
        self.aux.clear();
    }

    // -------------------------------------------------------------------------
    // Fluent setters
    // -------------------------------------------------------------------------

    /// Set the reference sequence ID (`ref_id`). Use `-1` for unmapped.
    pub fn ref_id(&mut self, v: i32) -> &mut Self {
        self.ref_id = v;
        self
    }

    /// Set the 0-based leftmost alignment position. Use `-1` for unmapped.
    pub fn pos(&mut self, v: i32) -> &mut Self {
        self.pos = v;
        self
    }

    /// Set the mapping quality.
    pub fn mapq(&mut self, v: u8) -> &mut Self {
        self.mapq = v;
        self
    }

    /// Set the bitwise flags.
    pub fn flags(&mut self, v: u16) -> &mut Self {
        self.flags = v;
        self
    }

    /// Set the mate reference sequence ID. Use `-1` for unmapped mate.
    pub fn mate_ref_id(&mut self, v: i32) -> &mut Self {
        self.mate_ref_id = v;
        self
    }

    /// Set the mate 0-based alignment position. Use `-1` for unmapped mate.
    pub fn mate_pos(&mut self, v: i32) -> &mut Self {
        self.mate_pos = v;
        self
    }

    /// Set the template length (`tlen`).
    pub fn template_length(&mut self, v: i32) -> &mut Self {
        self.tlen = v;
        self
    }

    /// Set the BAM bin (SAM spec §4.2.1). Defaults to `0`; callers may compute
    /// with `reg2bin` or leave as-is for unmapped records.
    pub fn bin(&mut self, v: u16) -> &mut Self {
        self.bin = v;
        self
    }

    /// Set the read name (without null terminator). Maximum 254 bytes.
    ///
    /// # Panics
    ///
    /// Panics in [`Self::build`] if `name.len() >= 255`.
    pub fn read_name(&mut self, name: &[u8]) -> &mut Self {
        self.read_name.clear();
        self.read_name.extend_from_slice(name);
        self
    }

    /// Set the CIGAR operations as BAM-encoded `u32` words
    /// (`(len << 4) | op_code`).
    pub fn cigar_ops(&mut self, ops: &[u32]) -> &mut Self {
        self.cigar_ops.clear();
        self.cigar_ops.extend_from_slice(ops);
        self
    }

    /// Set the ASCII sequence (e.g. `b"ACGT"`).
    pub fn sequence(&mut self, seq: &[u8]) -> &mut Self {
        self.sequence.clear();
        self.sequence.extend_from_slice(seq);
        self
    }

    /// Set the raw Phred quality scores (0–93, not Phred+33 ASCII).
    pub fn qualities(&mut self, quals: &[u8]) -> &mut Self {
        self.qualities.clear();
        self.qualities.extend_from_slice(quals);
        self
    }

    // -------------------------------------------------------------------------
    // Aux tag helpers
    // -------------------------------------------------------------------------

    /// Append a string (`Z`-type) aux tag.
    pub fn add_string_tag(&mut self, tag: &[u8; 2], value: &[u8]) -> &mut Self {
        append_string_tag(&mut self.aux, tag, value);
        self
    }

    /// Append an integer aux tag (smallest type that fits).
    pub fn add_int_tag(&mut self, tag: &[u8; 2], value: i32) -> &mut Self {
        append_int_tag(&mut self.aux, tag, value);
        self
    }

    /// Append a float (`f`-type) aux tag.
    pub fn add_float_tag(&mut self, tag: &[u8; 2], value: f32) -> &mut Self {
        append_float_tag(&mut self.aux, tag, value);
        self
    }

    /// Append a `u8` array (`B:C`-type) aux tag.
    pub fn add_array_u8(&mut self, tag: &[u8; 2], values: &[u8]) -> &mut Self {
        append_u8_array_tag(&mut self.aux, tag, values);
        self
    }

    /// Append a `u16` array (`B:S`-type) aux tag.
    pub fn add_array_u16(&mut self, tag: &[u8; 2], values: &[u16]) -> &mut Self {
        append_u16_array_tag(&mut self.aux, tag, values);
        self
    }

    /// Append an `i16` array (`B:s`-type) aux tag.
    pub fn add_array_i16(&mut self, tag: &[u8; 2], values: &[i16]) -> &mut Self {
        append_i16_array_tag(&mut self.aux, tag, values);
        self
    }

    /// Append an `i32` array (`B:i`-type) aux tag.
    pub fn add_array_i32(&mut self, tag: &[u8; 2], values: &[i32]) -> &mut Self {
        append_i32_array_tag(&mut self.aux, tag, values);
        self
    }

    /// Append an `f32` array (`B:f`-type) aux tag.
    pub fn add_array_f32(&mut self, tag: &[u8; 2], values: &[f32]) -> &mut Self {
        append_f32_array_tag(&mut self.aux, tag, values);
        self
    }

    // -------------------------------------------------------------------------
    // Finalize
    // -------------------------------------------------------------------------

    /// Assemble all fields into a [`RawRecord`] and return it.
    ///
    /// The builder's state is **not** cleared after calling `build()`. Call
    /// [`Self::clear`] explicitly before reusing for a different record.
    ///
    /// # Panics
    ///
    /// - If `read_name.len() >= 255` (BAM spec limit is 254 characters + NUL).
    /// - If `sequence` and `qualities` are both non-empty and their lengths differ.
    /// - If `sequence.len()` overflows `u32`.
    /// - If `cigar_ops.len()` overflows `u16`.
    #[must_use]
    pub fn build(&mut self) -> RawRecord {
        assert!(
            self.read_name.len() < 255,
            "read name too long ({} bytes, max 254)",
            self.read_name.len(),
        );
        assert!(
            self.qualities.is_empty() || self.sequence.len() == self.qualities.len(),
            "sequence length ({}) != qualities length ({})",
            self.sequence.len(),
            self.qualities.len(),
        );

        #[expect(
            clippy::cast_possible_truncation,
            reason = "assert above guarantees read_name.len() < 255"
        )]
        let l_read_name = (self.read_name.len() + 1) as u8; // +1 for NUL

        let n_cigar_op =
            u16::try_from(self.cigar_ops.len()).expect("too many CIGAR operations (max 65535)");

        let l_seq = u32::try_from(self.sequence.len()).expect("sequence length overflow");

        let packed_seq_len = self.sequence.len().div_ceil(2);
        let record_len = 32
            + l_read_name as usize
            + self.cigar_ops.len() * 4
            + packed_seq_len
            + self.sequence.len() // qual bytes (same length as seq)
            + self.aux.len();

        let mut buf = Vec::with_capacity(record_len);

        // === Fixed 32-byte header ===
        buf.extend_from_slice(&self.ref_id.to_le_bytes());
        buf.extend_from_slice(&self.pos.to_le_bytes());
        buf.push(l_read_name);
        buf.push(self.mapq);
        buf.extend_from_slice(&self.bin.to_le_bytes());
        buf.extend_from_slice(&n_cigar_op.to_le_bytes());
        buf.extend_from_slice(&self.flags.to_le_bytes());
        buf.extend_from_slice(&l_seq.to_le_bytes());
        buf.extend_from_slice(&self.mate_ref_id.to_le_bytes());
        buf.extend_from_slice(&self.mate_pos.to_le_bytes());
        buf.extend_from_slice(&self.tlen.to_le_bytes());

        // === Read name + NUL ===
        buf.extend_from_slice(&self.read_name);
        buf.push(0);

        // === CIGAR ops (4 bytes each) ===
        for &op in &self.cigar_ops {
            buf.extend_from_slice(&op.to_le_bytes());
        }

        // === Packed sequence (4 bits per base) ===
        if !self.sequence.is_empty() {
            pack_sequence_into(&mut buf, &self.sequence);

            // === Quality scores ===
            if self.qualities.is_empty() {
                // Absent quality: 0xFF per SAM spec
                buf.resize(buf.len() + self.sequence.len(), 0xFF);
            } else {
                buf.extend_from_slice(&self.qualities);
            }
        }

        // === Aux tags ===
        buf.extend_from_slice(&self.aux);

        RawRecord::from(buf)
    }
}

impl Default for SamBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
#[allow(clippy::identity_op)]
mod tests {
    use super::*;
    use crate::fields::*;
    use crate::sequence::*;
    use crate::tags::*;
    use crate::testutil::*;

    #[test]
    fn test_builder_basic_unmapped_record() {
        let mut builder = UnmappedSamBuilder::new();
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
        let mut builder = UnmappedSamBuilder::new();
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
        let ce = find_float_tag(aux, b"cE").expect("cE float tag should be present");
        assert!((ce - 0.05).abs() < 1e-7);
    }

    #[test]
    fn test_builder_reuse() {
        let mut builder = UnmappedSamBuilder::new();

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
        let mut builder = UnmappedSamBuilder::new();
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
        let mut builder = UnmappedSamBuilder::new();

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
        let mut builder = UnmappedSamBuilder::new();
        builder.build_record(b"empty", flags::UNMAPPED, b"", &[]);
        builder.append_string_tag(b"RG", b"rg0");

        let bam = builder.as_bytes();
        assert_eq!(l_seq(bam), 0);
        assert_eq!(find_string_tag(aux_data_slice(bam), b"RG"), Some(b"rg0" as &[u8]));
    }

    #[test]
    fn test_builder_missing_quality() {
        let mut builder = UnmappedSamBuilder::new();
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
        let mut builder = UnmappedSamBuilder::new();
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
        let mut builder = UnmappedSamBuilder::new();
        builder.build_record(b"test", flags::UNMAPPED, b"ACGT", &[30, 30, 30, 30]);
        builder.append_phred33_string_tag(b"aq", &[0, 10, 30, 40]);
        let bytes = builder.as_bytes();
        let parsed = ParsedBamRecord::from_bytes(bytes);
        let aq =
            parsed.get_string_tag(b"aq").expect("aq string tag should be present in parsed record");
        assert_eq!(aq, b"!+?I");
    }

    #[test]
    fn test_parsed_bam_record_roundtrip() {
        let mut builder = UnmappedSamBuilder::new();
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
        assert_eq!(rec.get_string_tag(b"RG").expect("RG tag should be present"), b"sample1");
        assert_eq!(rec.get_int_tag(b"cD").expect("cD tag should be present"), 5);
        let ce = rec.get_float_tag(b"cE").expect("cE float tag should be present");
        assert!((ce - 0.01).abs() < 0.001);
        assert_eq!(
            rec.get_i16_array_tag(b"cd").expect("cd i16 array tag should be present"),
            vec![3, 4, 5, 3, 2, 1]
        );
    }

    // ========================================================================
    // UnmappedSamBuilder additional tests
    // ========================================================================

    #[test]
    fn test_builder_build_returns_raw_record() {
        let mut builder = UnmappedSamBuilder::new();
        builder.build_record(b"r1", flags::UNMAPPED, b"ACGT", &[30, 30, 30, 30]);
        builder.append_string_tag(b"MI", b"7");

        let record = builder.build();

        // Should be a valid RawRecord with the expected fields.
        assert_eq!(l_seq(record.as_ref()), 4);
        assert_eq!(read_name(record.as_ref()), b"r1");
        assert_eq!(find_string_tag(aux_data_slice(record.as_ref()), b"MI"), Some(b"7" as &[u8]));

        // After build(), the builder's internal buffer should be empty (moved out).
        assert!(builder.buf.is_empty());
        assert!(!builder.sealed);

        // The builder can be reused.
        builder.build_record(b"r2", flags::UNMAPPED, b"AC", &[20, 25]);
        assert_eq!(l_seq(builder.as_bytes()), 2);
        assert_eq!(read_name(builder.as_bytes()), b"r2");
    }

    #[test]
    fn test_builder_build_into_inner() {
        let mut builder = UnmappedSamBuilder::new();
        builder.build_record(b"r3", flags::UNMAPPED, b"ACGT", &[10, 20, 30, 40]);

        let bytes_before = builder.as_bytes().to_vec();
        let record = builder.build();

        // into_inner() should yield the same bytes.
        assert_eq!(record.into_inner(), bytes_before);
    }

    #[test]
    fn test_builder_default() {
        let builder = UnmappedSamBuilder::default();
        // Default builder has empty buffer and is not sealed
        assert!(builder.buf.is_empty());
        assert!(!builder.sealed);
    }

    #[test]
    fn test_builder_with_capacity() {
        let builder = UnmappedSamBuilder::with_capacity(1024);
        assert!(builder.buf.capacity() >= 1024);
    }

    #[test]
    fn test_builder_odd_length_sequence() {
        // Odd-length sequence to exercise the packing edge case
        let mut builder = UnmappedSamBuilder::new();
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
        let mut builder = UnmappedSamBuilder::new();
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
        let mut builder = UnmappedSamBuilder::new();
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
        assert_eq!(
            records[0].get_string_tag(b"MI").expect("MI tag should be present on record 0"),
            b"1"
        );
        assert_eq!(records[1].name, b"r2");
        assert_eq!(records[1].bases, b"TGCA");
        assert_eq!(
            records[1].get_string_tag(b"MI").expect("MI tag should be present on record 1"),
            b"2"
        );
    }

    // ========================================================================
    // SamBuilder tests
    // ========================================================================

    /// Encode a BAM CIGAR operation: `(len << 4) | op_code`.
    /// Op codes: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
    const fn cigar_op(len: u32, op: u32) -> u32 {
        (len << 4) | op
    }

    #[test]
    fn test_mapped_builder_default_values() {
        // Default build should produce a valid minimal record.
        let mut b = SamBuilder::new();
        let rec = b.build();
        let bam = rec.as_ref();

        assert_eq!(ref_id(bam), -1);
        assert_eq!(pos(bam), -1);
        assert_eq!(mapq(bam), 0);
        assert_eq!(flags(bam), 0);
        assert_eq!(mate_ref_id(bam), -1);
        assert_eq!(mate_pos(bam), -1);
        assert_eq!(template_length(bam), 0);
        assert_eq!(l_seq(bam), 0);
        assert_eq!(read_name(bam), b"*");
        assert_eq!(n_cigar_op(bam), 0);
        assert!(aux_data_slice(bam).is_empty());
    }

    #[test]
    fn test_mapped_builder_all_header_fields() {
        // Set every header field and verify round-trip.
        let mut b = SamBuilder::new();
        b.ref_id(3)
            .pos(999)
            .mapq(60)
            .flags(flags::PAIRED | flags::REVERSE)
            .mate_ref_id(3)
            .mate_pos(1499)
            .template_length(500)
            .bin(4681)
            .read_name(b"myread")
            .cigar_ops(&[cigar_op(100, 0)]) // 100M
            .sequence(b"ACGTACGT")
            .qualities(&[30, 25, 35, 40, 30, 25, 35, 40]);
        let rec = b.build();
        let bam = rec.as_ref();

        assert_eq!(ref_id(bam), 3);
        assert_eq!(pos(bam), 999);
        assert_eq!(mapq(bam), 60);
        assert_eq!(flags(bam), flags::PAIRED | flags::REVERSE);
        assert_eq!(mate_ref_id(bam), 3);
        assert_eq!(mate_pos(bam), 1499);
        assert_eq!(template_length(bam), 500);
        // bin
        assert_eq!(u16::from_le_bytes([bam[10], bam[11]]), 4681u16);
        assert_eq!(read_name(bam), b"myread");
        assert_eq!(n_cigar_op(bam), 1);
        assert_eq!(l_seq(bam), 8);

        // Verify packed sequence
        let so = seq_offset(bam);
        let expected = [1u8, 2, 4, 8, 1, 2, 4, 8]; // A C G T A C G T
        for (i, &exp) in expected.iter().enumerate() {
            assert_eq!(get_base(bam, so, i), exp, "base mismatch at position {i}");
        }

        // Verify quality scores
        let qo = qual_offset(bam);
        let expected_quals = [30u8, 25, 35, 40, 30, 25, 35, 40];
        for (i, &exp) in expected_quals.iter().enumerate() {
            assert_eq!(get_qual(bam, qo, i), exp, "qual mismatch at position {i}");
        }

        assert!(aux_data_slice(bam).is_empty());
    }

    #[test]
    fn test_mapped_builder_with_all_tag_types() {
        // Verify every add_*_tag method appends bytes that can be read back.
        let mut b = SamBuilder::new();
        b.read_name(b"tagged").sequence(b"ACGT").qualities(&[30, 30, 30, 30]);
        b.add_string_tag(b"RG", b"lane1");
        b.add_int_tag(b"NM", 2);
        b.add_float_tag(b"AS", 98.5);
        b.add_array_u8(b"BC", &[10, 20, 30]);
        b.add_array_u16(b"QT", &[100, 200, 300]);
        b.add_array_i16(b"cd", &[-1, 0, 1]);
        b.add_array_i32(b"id", &[-100_000, 0, 100_000]);
        b.add_array_f32(b"fd", &[1.0_f32, 2.5, -3.0]);
        let rec = b.build();
        let bam = rec.as_ref();
        let aux = aux_data_slice(bam);

        assert_eq!(find_string_tag(aux, b"RG"), Some(b"lane1" as &[u8]));
        assert_eq!(find_int_tag(aux, b"NM"), Some(2));
        let as_val = find_float_tag(aux, b"AS").expect("AS float tag present");
        assert!((as_val - 98.5).abs() < 1e-4);

        // Verify array tags exist (type byte check)
        assert_eq!(find_tag_type(aux, b"BC"), Some(b'B'));
        assert_eq!(find_tag_type(aux, b"QT"), Some(b'B'));
        assert_eq!(find_tag_type(aux, b"cd"), Some(b'B'));
        assert_eq!(find_tag_type(aux, b"id"), Some(b'B'));
        assert_eq!(find_tag_type(aux, b"fd"), Some(b'B'));
    }

    #[test]
    fn test_mapped_builder_cigar_multple_ops() {
        // Build a record with a multi-op CIGAR: 10S 80M 10S
        let cigar = [cigar_op(10, 4), cigar_op(80, 0), cigar_op(10, 4)];
        let mut b = SamBuilder::new();
        b.ref_id(0).pos(100).mapq(30).read_name(b"clip").cigar_ops(&cigar);
        let rec = b.build();
        let bam = rec.as_ref();

        assert_eq!(n_cigar_op(bam), 3);
        // Verify raw CIGAR words stored in order at the cigar offset
        let cigar_off = 32 + l_read_name(bam) as usize;
        for (i, &op) in cigar.iter().enumerate() {
            let stored = u32::from_le_bytes([
                bam[cigar_off + i * 4],
                bam[cigar_off + i * 4 + 1],
                bam[cigar_off + i * 4 + 2],
                bam[cigar_off + i * 4 + 3],
            ]);
            assert_eq!(stored, op, "CIGAR op mismatch at index {i}");
        }
    }

    #[test]
    fn test_mapped_builder_empty_cigar_unmapped_like() {
        // A record with ref_id=-1, pos=-1, no cigar is valid (unmapped-like).
        let mut b = SamBuilder::new();
        b.flags(flags::UNMAPPED).read_name(b"unmap").sequence(b"ACGT").qualities(&[30; 4]);
        let rec = b.build();
        let bam = rec.as_ref();

        assert_eq!(ref_id(bam), -1);
        assert_eq!(pos(bam), -1);
        assert_eq!(n_cigar_op(bam), 0);
        assert_eq!(l_seq(bam), 4);
    }

    #[test]
    fn test_mapped_builder_absent_qualities_filled_with_ff() {
        // When sequence is set but qualities are omitted, emit 0xFF per SAM spec.
        let mut b = SamBuilder::new();
        b.read_name(b"r1").sequence(b"ACGT"); // no qualities()
        let rec = b.build();
        let bam = rec.as_ref();

        assert_eq!(l_seq(bam), 4);
        let qo = qual_offset(bam);
        for i in 0..4 {
            assert_eq!(get_qual(bam, qo, i), 0xFF, "expected 0xFF sentinel at pos {i}");
        }
    }

    #[test]
    fn test_mapped_builder_empty_seq_and_qual() {
        // Both empty: l_seq == 0, no seq/qual bytes emitted.
        let mut b = SamBuilder::new();
        b.read_name(b"noSeq");
        let rec = b.build();
        let bam = rec.as_ref();
        assert_eq!(l_seq(bam), 0);
        // The record length should be exactly: 32 (header) + 6 (name+NUL) + 0 (CIGAR) + 0 (seq/qual)
        assert_eq!(bam.len(), 32 + 6);
    }

    #[test]
    fn test_mapped_builder_clear_and_reuse() {
        let mut b = SamBuilder::new();

        // First record
        b.ref_id(1).pos(100).read_name(b"r1").sequence(b"ACGT").qualities(&[30; 4]);
        b.add_string_tag(b"RG", b"sample1");
        let rec1 = b.build();
        let bam1 = rec1.as_ref();
        assert_eq!(ref_id(bam1), 1);
        assert_eq!(pos(bam1), 100);
        assert_eq!(read_name(bam1), b"r1");
        assert_eq!(find_string_tag(aux_data_slice(bam1), b"RG"), Some(b"sample1" as &[u8]));

        // Clear and build a second, different record
        b.clear();
        b.ref_id(2).pos(500).read_name(b"r2").sequence(b"TTTT").qualities(&[25; 4]);
        b.add_int_tag(b"NM", 0);
        let rec2 = b.build();
        let bam2 = rec2.as_ref();

        assert_eq!(ref_id(bam2), 2);
        assert_eq!(pos(bam2), 500);
        assert_eq!(read_name(bam2), b"r2");
        // Old tag should be gone
        assert!(find_string_tag(aux_data_slice(bam2), b"RG").is_none());
        assert_eq!(find_int_tag(aux_data_slice(bam2), b"NM"), Some(0));
    }

    #[test]
    fn test_mapped_builder_default_impl() {
        let b = SamBuilder::default();
        assert_eq!(b.ref_id, -1);
        assert_eq!(b.pos, -1);
        assert_eq!(b.read_name, b"*");
    }

    #[cfg(feature = "noodles")]
    #[test]
    fn test_mapped_builder_noodles_roundtrip() {
        use crate::noodles_compat::raw_record_to_record_buf;
        use noodles::sam::Header;
        use noodles::sam::alignment::record::Flags;
        use noodles::sam::alignment::record::Sequence as _;

        // Build a mapped record and decode it via noodles to cross-validate.
        let header = Header::default();
        let mut b = SamBuilder::new();
        b.flags(flags::PAIRED | flags::FIRST_SEGMENT)
            .read_name(b"noodles_read")
            .sequence(b"ACGTACGT")
            .qualities(&[30, 25, 35, 40, 30, 25, 35, 40]);
        b.add_string_tag(b"RG", b"sg1");
        b.add_int_tag(b"NM", 1);
        let rec = b.build();

        let record_buf = raw_record_to_record_buf(&rec, &header)
            .expect("noodles should decode the record without error");

        assert_eq!(
            record_buf.flags(),
            Flags::SEGMENTED | Flags::FIRST_SEGMENT,
            "flags round-trip mismatch"
        );
        // Sequence
        let seq_bytes: Vec<u8> = record_buf.sequence().iter().collect();
        assert_eq!(seq_bytes, b"ACGTACGT");
        // Quality scores
        let qual_bytes: Vec<u8> = record_buf.quality_scores().as_ref().to_vec();
        assert_eq!(qual_bytes, [30, 25, 35, 40, 30, 25, 35, 40]);
    }
}
