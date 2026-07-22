use crate::fields::{
    flags,
    flags::{REVERSE, UNMAPPED},
    l_read_name, n_cigar_op, pos,
};

/// Converts a raw BAM CIGAR u32 op word to a noodles `Kind`.
///
/// Extracts the 4-bit op type from the raw u32 and maps it to the corresponding
/// `Kind` variant. Returns `Kind::Pad` for unknown op types.
#[cfg(feature = "noodles")]
#[inline]
#[must_use]
pub fn cigar_op_kind(raw_op: u32) -> noodles::sam::alignment::record::cigar::op::Kind {
    use noodles::sam::alignment::record::cigar::op::Kind;
    let op = raw_op & 0xF;
    match op {
        0 => Kind::Match,
        1 => Kind::Insertion,
        2 => Kind::Deletion,
        3 => Kind::Skip,
        4 => Kind::SoftClip,
        5 => Kind::HardClip,
        7 => Kind::SequenceMatch,
        8 => Kind::SequenceMismatch,
        // 6 => Pad; 9..=15 are reserved/invalid per the BAM spec and fall back
        // to Pad as the documented default so callers don't need to branch.
        _ => Kind::Pad,
    }
}

/// BAM CIGAR op-type -> (`consume_query`, `consume_ref`) packed bitmask.
///
/// Encodes both properties for every op type (0..8) in a single u32. Bit
/// `(op << 1)` is `consume_query`; bit `(op << 1 | 1)` is `consume_ref`. Bits
/// for op values >= 9 are zero. Matches htslib's `BAM_CIGAR_TYPE` in
/// `sam.h`, giving branchless `consumes_*` queries.
///
/// | op | mnemonic | consume_query | consume_ref |
/// |----|----------|---------------|-------------|
/// | 0  | M        | yes           | yes         |
/// | 1  | I        | yes           | no          |
/// | 2  | D        | no            | yes         |
/// | 3  | N        | no            | yes         |
/// | 4  | S        | yes           | no          |
/// | 5  | H        | no            | no          |
/// | 6  | P        | no            | no          |
/// | 7  | =        | yes           | yes         |
/// | 8  | X        | yes           | yes         |
const BAM_CIGAR_TYPE: u32 = 0x3C1A7;

/// Returns true if the CIGAR op type consumes the reference (M, D, N, =, X).
///
/// Values outside the 0..=15 spec range are masked with `& 0xF` — keeps the
/// shift amount bounded in case a caller passes a raw CIGAR word or a noodles
/// conversion bug leaks a larger value.
#[inline]
#[must_use]
pub const fn consumes_ref(op_type: u32) -> bool {
    (BAM_CIGAR_TYPE >> ((op_type & 0xF) << 1)) & 2 != 0
}

/// Returns true if the CIGAR op type consumes the query including soft clips (M, I, S, =, X).
///
/// See [`consumes_ref`] for the `& 0xF` rationale.
#[inline]
#[must_use]
pub const fn consumes_query(op_type: u32) -> bool {
    (BAM_CIGAR_TYPE >> ((op_type & 0xF) << 1)) & 1 != 0
}

/// Returns true if the CIGAR op type consumes the read/query excluding soft clips (M, I, =, X).
#[inline]
#[must_use]
pub const fn consumes_read(op_type: u32) -> bool {
    matches!(op_type, 0 | 1 | 7 | 8)
}

/// Extract CIGAR operations from BAM record.
#[inline]
#[must_use]
pub fn get_cigar_ops(bam: &[u8]) -> Vec<u32> {
    let l_read_name = l_read_name(bam) as usize;
    let n_cigar_op = n_cigar_op(bam) as usize;

    if n_cigar_op == 0 {
        return Vec::new();
    }

    let cigar_start = 32 + l_read_name;
    let cigar_end = cigar_start + n_cigar_op * 4;

    if cigar_end > bam.len() {
        return Vec::new();
    }

    // Read CIGAR ops bytewise to avoid alignment issues — the CIGAR data offset
    // (32 + l_read_name) is not guaranteed to be 4-byte aligned.
    let cigar_bytes = &bam[cigar_start..cigar_end];
    cigar_bytes.chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect()
}

/// CIGAR operation type to SAM character mapping.
const CIGAR_OP_CHARS: [u8; 9] = [b'M', b'I', b'D', b'N', b'S', b'H', b'P', b'=', b'X'];

/// Formats CIGAR operations from a raw BAM record as a SAM-style CIGAR string.
///
/// Returns an empty string if the record has no CIGAR operations, or `"*"` is not used
/// (callers should check `n_cigar_op == 0` separately if they need `"*"`).
///
/// # Panics
///
/// Panics if a CIGAR operation has an invalid type (>= 9).
#[must_use]
pub fn cigar_to_string_from_raw(bam: &[u8]) -> String {
    let ops = get_cigar_ops(bam);
    if ops.is_empty() {
        return String::new();
    }
    let mut result = String::with_capacity(ops.len() * 4);
    for &op in &ops {
        let op_len = op >> 4;
        let op_type = (op & 0xF) as usize;
        assert!(op_type < CIGAR_OP_CHARS.len(), "invalid CIGAR op type: {op_type}");
        result.push_str(&op_len.to_string());
        result.push(CIGAR_OP_CHARS[op_type] as char);
    }
    result
}

/// Calculate reference-consuming length from CIGAR operations.
///
/// This is the sum of M/D/N/=/X operations, which represents how many
/// reference bases the alignment spans. Used for BAM index generation.
#[inline]
#[must_use]
pub fn reference_length_from_cigar(cigar_ops: &[u32]) -> i32 {
    let mut ref_len = 0i32;

    for &op in cigar_ops {
        let op_len = (op >> 4).cast_signed();
        let op_type = op & 0xF;

        if consumes_ref(op_type) {
            ref_len += op_len;
        }
    }

    ref_len
}

/// Compute reference length directly from raw CIGAR bytes (zero allocation).
///
/// Returns `0` for a record with no CIGAR **and** for a malformed record whose
/// declared read-name/CIGAR region extends past the buffer. Callers that must
/// distinguish those two cases (e.g. to fail closed on corrupt input) should use
/// [`reference_length_from_raw_bam_checked`] instead.
#[inline]
#[must_use]
pub fn reference_length_from_raw_bam(bam: &[u8]) -> i32 {
    reference_length_from_raw_bam_checked(bam).unwrap_or(0)
}

/// Like [`reference_length_from_raw_bam`], but validates the record's read-name
/// and CIGAR bounds against the buffer and distinguishes a legitimately
/// zero-length span from a truncated/malformed record.
///
/// Returns:
/// - `Some(len)` for a well-formed record — `Some(0)` when it genuinely has no
///   CIGAR (`n_cigar_op == 0`), otherwise the reference-consuming span.
/// - `None` when the record is malformed: shorter than the 32-byte fixed core, or
///   its declared read-name/CIGAR region extends past the buffer (validated with
///   checked arithmetic so a corrupt `l_read_name`/`n_cigar_op` cannot overflow the
///   offset math and wrap back into bounds). Callers must reject such records
///   rather than treat them as a zero-length (point) alignment at their start.
#[inline]
#[must_use]
pub fn reference_length_from_raw_bam_checked(bam: &[u8]) -> Option<i32> {
    use crate::fields::MIN_BAM_RECORD_LEN;

    // The 32-byte fixed core must be present to read `l_read_name`/`n_cigar_op`.
    if bam.len() < MIN_BAM_RECORD_LEN {
        return None;
    }
    let n_cigar_op = n_cigar_op(bam) as usize;
    let l_read_name = l_read_name(bam) as usize;
    let cigar_start = MIN_BAM_RECORD_LEN.checked_add(l_read_name)?;
    let cigar_end = cigar_start.checked_add(n_cigar_op.checked_mul(4)?)?;
    if cigar_end > bam.len() {
        return None;
    }
    // Validate the read-name/CIGAR bounds *before* short-circuiting on a zero-op
    // CIGAR, so a core-only record declaring an oversized `l_read_name` fails closed
    // rather than being accepted as a valid zero-length (point) alignment.
    if n_cigar_op == 0 {
        return Some(0);
    }

    let mut ref_len = 0i32;
    for i in 0..n_cigar_op {
        let op = cigar_op_at(bam, cigar_start + i * 4);
        let op_type = op & 0xF;
        if consumes_ref(op_type) {
            // Fail closed on a malformed-but-in-bounds CIGAR whose ref-consuming
            // lengths overflow `i32`: return `None` rather than a wrapped (negative)
            // span that downstream span math would trust as non-negative.
            ref_len = ref_len.checked_add((op >> 4).cast_signed())?;
        }
    }
    Some(ref_len)
}

/// Compute the query-consuming length of CIGAR operations (the "read length").
///
/// This is the sum of M/I/S/=/X operations.
#[inline]
#[must_use]
pub fn query_length_from_cigar(cigar_ops: &[u32]) -> usize {
    let mut len = 0usize;
    for &op in cigar_ops {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        if consumes_query(op_type) {
            len += op_len;
        }
    }
    len
}

/// Compute unclipped 5' position directly from raw BAM bytes (zero allocation).
///
/// This is the primary entry point for raw-byte callers, replacing the pattern:
/// `let cigar = get_cigar_ops(bam); unclipped_5prime_1based(pos, reverse, unmapped, &cigar)`
#[inline]
#[must_use]
pub fn unclipped_5prime_from_raw_bam(bam: &[u8]) -> i32 {
    let flg = flags(bam);
    let unmapped = (flg & UNMAPPED) != 0;
    if unmapped {
        return 0;
    }

    let n_cigar_op = n_cigar_op(bam) as usize;
    if n_cigar_op == 0 {
        return i32::MAX;
    }

    let l_read_name = l_read_name(bam) as usize;
    let cigar_start = 32 + l_read_name;
    let cigar_end = cigar_start + n_cigar_op * 4;
    if cigar_end > bam.len() {
        return i32::MAX;
    }

    let reverse = (flg & REVERSE) != 0;
    let pos_1based = pos(bam) + 1;

    if reverse {
        unclipped_end_from_raw_cigar(pos_1based, bam, cigar_start, n_cigar_op)
    } else {
        unclipped_start_from_raw_cigar(pos_1based, bam, cigar_start, n_cigar_op)
    }
}

/// Compute unclipped 5' position directly from raw BAM bytes (zero allocation).
///
/// Like [`unclipped_5prime_from_raw_bam`] but uses pre-extracted `pos` (0-based) and
/// `is_reverse` flag, avoiding redundant field reads when the caller already has them.
/// Returns the position in the same basis as the input `pos`.
#[inline]
#[must_use]
pub fn unclipped_5prime_raw(bam: &[u8], pos: i32, is_reverse: bool) -> i32 {
    let n_cigar_op = n_cigar_op(bam) as usize;
    if n_cigar_op == 0 {
        return pos;
    }

    let l_read_name = l_read_name(bam) as usize;
    let cigar_start = 32 + l_read_name;
    let cigar_end = cigar_start + n_cigar_op * 4;
    if cigar_end > bam.len() {
        return pos;
    }

    if is_reverse {
        unclipped_end_from_raw_cigar(pos, bam, cigar_start, n_cigar_op)
    } else {
        unclipped_start_from_raw_cigar(pos, bam, cigar_start, n_cigar_op)
    }
}

/// Calculate mate's unclipped 5' coordinate from MC tag.
///
/// For forward strand mate: `unclipped_start`
/// For reverse strand mate: `unclipped_end`
#[inline]
#[must_use]
pub fn mate_unclipped_5prime(mate_pos: i32, mate_reverse: bool, mc_cigar: &str) -> i32 {
    if mate_reverse {
        unclipped_other_end(mate_pos, mc_cigar)
    } else {
        unclipped_other_start(mate_pos, mc_cigar)
    }
}

/// Calculate mate's unclipped 5' coordinate using 1-based positions (matching noodles).
///
/// Adds 1 to the 0-based BAM mate position before applying clip adjustment.
#[inline]
#[must_use]
pub fn mate_unclipped_5prime_1based(
    mate_pos_0based: i32,
    mate_reverse: bool,
    mc_cigar: &str,
) -> i32 {
    let mate_pos_1based = mate_pos_0based.saturating_add(1);
    if mate_reverse {
        unclipped_other_end(mate_pos_1based, mc_cigar)
    } else {
        unclipped_other_start(mate_pos_1based, mc_cigar)
    }
}

/// Calculate mate's unclipped start from MC tag CIGAR string.
///
/// For forward strand mates, this is the 5' position.
///
/// The result is deliberately *not* clamped to the first coordinate: a read
/// whose leading clips extend past the start of the contig has an unclipped
/// start below it, and that is true of real data, not just malformed MC values.
/// The record's own side ([`unclipped_start_from_raw_cigar`]) does not clamp
/// either, and the template-coordinate key compares the two against each other
/// (see `external.rs`), so clamping only here would make the two records of one
/// template canonicalize differently — and would diverge from fgbio and
/// `samtools sort --template-coordinate`.
#[inline]
#[must_use]
pub(crate) fn unclipped_other_start(mate_pos: i32, mc_cigar: &str) -> i32 {
    mate_pos.saturating_sub(parse_leading_clips(mc_cigar.as_bytes()))
}

/// Calculate mate's unclipped end from MC tag CIGAR string.
///
/// For reverse strand mates, this is the 5' position.
///
/// The `-1` that converts a length to an inclusive end is applied *before* the
/// saturating additions: doing it last would turn a sum that legitimately lands
/// on [`i32::MAX`] into `i32::MAX - 1`, because the addition saturates first and
/// the subtraction then walks the saturated value back.
#[inline]
#[must_use]
pub(crate) fn unclipped_other_end(mate_pos: i32, mc_cigar: &str) -> i32 {
    let (ref_len, trailing_clips) = parse_ref_len_and_trailing_clips(mc_cigar.as_bytes());
    mate_pos.saturating_sub(1).saturating_add(ref_len).saturating_add(trailing_clips)
}

/// Compute 1-based alignment end position from raw BAM bytes.
///
/// Returns `pos + ref_len` (1-based inclusive end, matching noodles convention).
/// Returns `None` if the record is unmapped or has no CIGAR.
#[inline]
#[must_use]
pub fn alignment_end_from_raw(bam: &[u8]) -> Option<usize> {
    let p = pos(bam);
    if p < 0 {
        return None;
    }
    let ref_len = reference_length_from_raw_bam(bam);
    if ref_len == 0 {
        return None;
    }
    // pos is 0-based, convert to 1-based and add ref_len - 1 for inclusive end
    Some((p + ref_len).cast_unsigned() as usize)
}

/// Compute 1-based alignment start position from raw BAM bytes.
///
/// Returns `None` if `pos < 0` (unmapped).
#[inline]
#[must_use]
pub fn alignment_start_from_raw(bam: &[u8]) -> Option<usize> {
    let p = pos(bam);
    if p < 0 { None } else { Some((p + 1).cast_unsigned() as usize) }
}

/// Virtual CIGAR clipping on raw u32 CIGAR ops.
///
/// Computes clipped CIGAR ops and the number of reference bases consumed by clipping,
/// without modifying any record. This is the raw-byte equivalent of
/// `SamRecordClipper::clip_start_of_read` / `clip_end_of_read` with Hard clipping mode.
///
/// Like the Clipper, this first accounts for existing H+S clips at the relevant end.
/// If `clip_amount` <= existing clips, only soft clips are upgraded to hard clips
/// (no alignment change). Otherwise, `clip_amount - existing_clips` bases are clipped
/// from the alignment.
///
/// Returns `(new_cigar_ops, ref_bases_consumed)`.
/// `ref_bases_consumed` is used to adjust `alignment_start` for start-clipping.
///
/// # Panics
///
/// Panics if a CIGAR operation length exceeds `u32::MAX`.
#[must_use]
pub fn clip_cigar_ops_raw(
    cigar_ops: &[u32],
    clip_amount: usize,
    from_start: bool,
) -> (Vec<u32>, usize) {
    if clip_amount == 0 || cigar_ops.is_empty() {
        return (cigar_ops.to_vec(), 0);
    }

    // Helper to encode a CIGAR op as raw u32
    let encode_op = |op_type: u32, len: usize| -> u32 {
        (u32::try_from(len).expect("CIGAR op length overflows u32") << 4) | op_type
    };

    // Count existing H+S clips at the relevant end (matching clip_*_of_read)
    let existing_clip: usize = if from_start {
        cigar_ops
            .iter()
            .take_while(|&&op| matches!(op & 0xF, 4 | 5))
            .map(|&op| (op >> 4) as usize)
            .sum()
    } else {
        cigar_ops
            .iter()
            .rev()
            .take_while(|&&op| matches!(op & 0xF, 4 | 5))
            .map(|&op| (op >> 4) as usize)
            .sum()
    };

    if clip_amount <= existing_clip {
        // Just upgrade soft clips to hard clips (no alignment change)
        upgrade_clipping_raw(cigar_ops, clip_amount, from_start, encode_op)
    } else {
        // Clip into alignment: clip (clip_amount - existing_clip) additional bases
        let alignment_clip = clip_amount - existing_clip;
        if from_start {
            clip_cigar_start_raw(cigar_ops, alignment_clip, encode_op)
        } else {
            clip_cigar_end_raw(cigar_ops, alignment_clip, encode_op)
        }
    }
}

/// Returns the query position (1-based) at a given reference position, from raw CIGAR ops.
///
/// This is the raw-byte equivalent of `CodecConsensusCaller::read_pos_at_ref_pos`.
///
/// # Arguments
/// * `cigar_ops` - Raw u32 CIGAR operations
/// * `alignment_start` - 1-based alignment start position
/// * `ref_pos` - 1-based reference position to query
/// * `return_last_base_if_deleted` - If true, returns last query position before deletion
///
/// Returns `None` if the position falls outside the alignment or in a deletion
/// (when `return_last_base_if_deleted` is false).
#[must_use]
pub fn read_pos_at_ref_pos_raw(
    cigar_ops: &[u32],
    alignment_start: usize,
    ref_pos: usize,
    return_last_base_if_deleted: bool,
) -> Option<usize> {
    if ref_pos < alignment_start {
        return None;
    }

    let mut ref_offset = 0usize;
    let mut query_offset = 0usize;

    for &op in cigar_ops {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;

        let op_ref_start = alignment_start + ref_offset;

        if consumes_ref(op_type) {
            let op_ref_end = op_ref_start + op_len - 1;

            if ref_pos >= op_ref_start && ref_pos <= op_ref_end {
                if consumes_query(op_type) {
                    // M, =, X: we have a base at this position
                    let offset_in_op = ref_pos - op_ref_start;
                    return Some(query_offset + offset_in_op + 1); // 1-based
                }
                // D, N: position falls in a deletion
                if return_last_base_if_deleted {
                    return Some(if query_offset > 0 { query_offset } else { 1 });
                }
                return None;
            }
        }

        if consumes_ref(op_type) {
            ref_offset += op_len;
        }
        if consumes_query(op_type) {
            query_offset += op_len;
        }
    }

    None
}

/// Read a single CIGAR op directly from raw BAM bytes at the given byte offset.
#[inline]
fn cigar_op_at(bam: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes([bam[offset], bam[offset + 1], bam[offset + 2], bam[offset + 3]])
}

/// Compute unclipped start by iterating raw CIGAR bytes directly (zero allocation).
///
/// Equivalent to `unclipped_start_from_cigar(pos, &get_cigar_ops(bam))` but avoids
/// the `Vec<u32>` allocation.
#[inline]
fn unclipped_start_from_raw_cigar(pos: i32, bam: &[u8], cigar_start: usize, n_ops: usize) -> i32 {
    let mut clipped = 0i32;
    for i in 0..n_ops {
        let op = cigar_op_at(bam, cigar_start + i * 4);
        let op_type = op & 0xF;
        match op_type {
            4 | 5 => clipped += (op >> 4).cast_signed(), // S or H
            _ => break,
        }
    }
    pos - clipped
}

/// Compute unclipped end by iterating raw CIGAR bytes directly (zero allocation).
///
/// Equivalent to `unclipped_end_from_cigar(pos, &get_cigar_ops(bam))` but avoids
/// the `Vec<u32>` allocation.
#[inline]
fn unclipped_end_from_raw_cigar(pos: i32, bam: &[u8], cigar_start: usize, n_ops: usize) -> i32 {
    let mut ref_len = 0i32;
    let mut trailing_clips = 0i32;
    let mut saw_ref_op = false;

    for i in 0..n_ops {
        let op = cigar_op_at(bam, cigar_start + i * 4);
        let op_len = (op >> 4).cast_signed();
        let op_type = op & 0xF;

        if consumes_ref(op_type) {
            ref_len += op_len;
            trailing_clips = 0;
            saw_ref_op = true;
        } else if matches!(op_type, 4 | 5) && saw_ref_op {
            trailing_clips += op_len;
        }
    }

    pos + ref_len + trailing_clips - 1
}

/// Accumulator for one CIGAR run-length, matching `str::parse::<i32>()`'s
/// behavior on the digit runs these parsers can produce.
///
/// The MC tag is *text*, so its run lengths are not bounded by BAM's 28-bit
/// CIGAR op field and its bytes are not guaranteed to be ASCII. Accumulating
/// here — rather than slicing the `&str` by byte offset and calling `parse` —
/// matters for three reasons:
///
/// 1. **No panic.** Slicing a `&str` at an offset landing inside a multi-byte
///    UTF-8 character panics. `std::str::from_utf8` admits such values (they are
///    valid UTF-8, merely not ASCII), so the previous form could panic in
///    release on a malformed MC reaching `fgumi sort`'s template key or the
///    decode step's `GroupKey`.
/// 2. **No re-scan.** The caller is already iterating those digits.
/// 3. **Overflow parity.** `str::parse::<i32>()` *fails* on a run wider than
///    `i32`, which `unwrap_or(0)` turns into `0` — it does not saturate. The
///    `u64` accumulator plus `i32::try_from(..).unwrap_or(0)` reproduces that
///    exactly, including runs too wide even for `u64`, and the empty run
///    (`"".parse()` is likewise an error, hence `0`).
///
/// `parse_mi_bytes` in `tags.rs` uses the same idiom.
struct DigitRun(u64);

impl DigitRun {
    const fn new() -> Self {
        Self(0)
    }

    /// Feeds one byte. Returns `true` when it was a digit and was consumed.
    fn push(&mut self, byte: u8) -> bool {
        if byte.is_ascii_digit() {
            self.0 = self.0.saturating_mul(10).saturating_add(u64::from(byte - b'0'));
            true
        } else {
            false
        }
    }

    /// Returns the accumulated run and resets, mirroring
    /// `slice.parse::<i32>().unwrap_or(0)`.
    fn take(&mut self) -> i32 {
        let value = std::mem::replace(&mut self.0, 0);
        i32::try_from(value).unwrap_or(0)
    }
}

/// Parse leading clips (S or H) from CIGAR bytes.
///
/// Returns the total length of leading soft/hard clips, saturating rather than
/// wrapping if a malformed MC tag declares runs that overflow `i32`.
#[inline]
fn parse_leading_clips(cigar: &[u8]) -> i32 {
    let mut clipped = 0i32;
    let mut run = DigitRun::new();

    for &c in cigar {
        if run.push(c) {
            continue;
        }
        let num = run.take();

        if c == b'S' || c == b'H' {
            clipped = clipped.saturating_add(num);
        } else {
            // Non-clip operation, stop processing
            break;
        }
    }

    clipped
}

/// Parse reference length and trailing clips from CIGAR bytes.
///
/// Returns `(reference_length, trailing_clips)`.
/// Reference length is the sum of M/D/N/=/X operations.
/// Trailing clips are S/H operations after the last reference-consuming op.
#[inline]
fn parse_ref_len_and_trailing_clips(cigar: &[u8]) -> (i32, i32) {
    let mut ref_len = 0i32;
    let mut trailing_clips = 0i32;
    let mut run = DigitRun::new();
    let mut saw_ref_op = false;

    for &c in cigar {
        if run.push(c) {
            continue;
        }
        let num = run.take();

        match c {
            b'M' | b'D' | b'N' | b'=' | b'X' => {
                ref_len = ref_len.saturating_add(num);
                trailing_clips = 0; // Reset trailing clips
                saw_ref_op = true;
            }
            b'S' | b'H' if saw_ref_op => {
                trailing_clips = trailing_clips.saturating_add(num);
            }
            _ => {}
        }
    }

    (ref_len, trailing_clips)
}

/// Upgrade existing soft clips to hard clips without changing alignment.
///
/// Matches `SamRecordClipper::upgrade_clipping` in Hard mode.
/// Used when `clip_amount` <= existing H+S clips.
fn upgrade_clipping_raw(
    cigar_ops: &[u32],
    clip_amount: usize,
    from_start: bool,
    encode_op: impl Fn(u32, usize) -> u32,
) -> (Vec<u32>, usize) {
    if from_start {
        let mut existing_hard = 0usize;
        let mut existing_soft = 0usize;
        let mut skip_count = 0usize;

        for &op in cigar_ops {
            if op & 0xF == 5 {
                existing_hard += (op >> 4) as usize;
                skip_count += 1;
            } else {
                break;
            }
        }
        for &op in &cigar_ops[skip_count..] {
            if op & 0xF == 4 {
                existing_soft += (op >> 4) as usize;
                skip_count += 1;
            } else {
                break;
            }
        }

        let length_to_upgrade = existing_soft.min(clip_amount.saturating_sub(existing_hard));
        let new_hard = existing_hard + length_to_upgrade;
        let remaining_soft = existing_soft - length_to_upgrade;

        let mut result = Vec::new();
        result.push(encode_op(5, new_hard));
        if remaining_soft > 0 {
            result.push(encode_op(4, remaining_soft));
        }
        result.extend_from_slice(&cigar_ops[skip_count..]);

        (result, 0)
    } else {
        let mut existing_hard = 0usize;
        let mut existing_soft = 0usize;
        let mut skip_count = 0usize;

        for &op in cigar_ops.iter().rev() {
            if op & 0xF == 5 {
                existing_hard += (op >> 4) as usize;
                skip_count += 1;
            } else {
                break;
            }
        }
        let end_idx = cigar_ops.len() - skip_count;
        for &op in cigar_ops[..end_idx].iter().rev() {
            if op & 0xF == 4 {
                existing_soft += (op >> 4) as usize;
                skip_count += 1;
            } else {
                break;
            }
        }

        let length_to_upgrade = existing_soft.min(clip_amount.saturating_sub(existing_hard));
        let new_hard = existing_hard + length_to_upgrade;
        let remaining_soft = existing_soft - length_to_upgrade;

        let end_content = cigar_ops.len() - skip_count;
        let mut result = cigar_ops[..end_content].to_vec();
        if remaining_soft > 0 {
            result.push(encode_op(4, remaining_soft));
        }
        result.push(encode_op(5, new_hard));

        (result, 0)
    }
}

/// Clip from the start of alignment (hard clip mode).
fn clip_cigar_start_raw(
    cigar_ops: &[u32],
    clip_amount: usize,
    encode_op: impl Fn(u32, usize) -> u32,
) -> (Vec<u32>, usize) {
    // Extract existing hard and soft clips from the start
    let mut existing_hard_clip = 0usize;
    let mut existing_soft_clip = 0usize;
    let mut skip_count = 0usize;

    for &op in cigar_ops {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        if op_type == 5 {
            // H
            existing_hard_clip += op_len;
            skip_count += 1;
        } else {
            break;
        }
    }
    for &op in &cigar_ops[skip_count..] {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        if op_type == 4 {
            // S
            existing_soft_clip += op_len;
            skip_count += 1;
        } else {
            break;
        }
    }

    let post_clip_ops = &cigar_ops[skip_count..];

    let mut read_bases_clipped = 0usize;
    let mut ref_bases_clipped = 0usize;
    let mut new_ops: Vec<u32> = Vec::new();
    let mut idx = 0;

    // Clip operations from the start
    while idx < post_clip_ops.len() {
        let op = post_clip_ops[idx];
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;

        // Check if we've clipped enough but need to skip trailing deletions
        if read_bases_clipped == clip_amount && new_ops.is_empty() && op_type == 2 {
            // Deletion: skip it (ref-consuming only)
            ref_bases_clipped += op_len;
            idx += 1;
            continue;
        }

        if read_bases_clipped >= clip_amount {
            break;
        }

        // Note: S (soft clip) is NOT read-consuming here, matching the Clipper.
        // After stripping leading H+S, remaining S ops (e.g. trailing) should not
        // count toward read bases clipped.
        let is_read = consumes_read(op_type); // M, I, =, X
        let is_ref = consumes_ref(op_type); // M, D, N, =, X

        if is_read && op_len > (clip_amount - read_bases_clipped) {
            if op_type == 1 {
                // Insertion: consume entire at clip boundary
                read_bases_clipped += op_len;
            } else {
                // Split the operation
                let remaining_clip = clip_amount - read_bases_clipped;
                let remaining_length = op_len - remaining_clip;
                read_bases_clipped += remaining_clip;
                if is_ref {
                    ref_bases_clipped += remaining_clip;
                }
                new_ops.push(encode_op(op_type, remaining_length));
            }
        } else {
            if is_read {
                read_bases_clipped += op_len;
            }
            if is_ref {
                ref_bases_clipped += op_len;
            }
        }

        idx += 1;
    }

    // Add remaining operations
    new_ops.extend_from_slice(&post_clip_ops[idx..]);

    // Hard clip mode: convert all existing soft clips to hard clips
    let added_hard_clip = existing_soft_clip + read_bases_clipped;
    let total_hard_clip = existing_hard_clip + added_hard_clip;
    let mut result = Vec::with_capacity(1 + new_ops.len());
    result.push(encode_op(5, total_hard_clip)); // H
    result.extend(new_ops);

    (result, ref_bases_clipped)
}

/// Clip from the end of alignment (hard clip mode).
fn clip_cigar_end_raw(
    cigar_ops: &[u32],
    clip_amount: usize,
    encode_op: impl Fn(u32, usize) -> u32,
) -> (Vec<u32>, usize) {
    // Extract existing hard and soft clips from the end
    let mut existing_hard_clip = 0usize;
    let mut existing_soft_clip = 0usize;
    let mut skip_count = 0usize;

    for &op in cigar_ops.iter().rev() {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        if op_type == 5 {
            // H
            existing_hard_clip += op_len;
            skip_count += 1;
        } else {
            break;
        }
    }
    let end_idx = cigar_ops.len() - skip_count;
    for &op in cigar_ops[..end_idx].iter().rev() {
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;
        if op_type == 4 {
            // S
            existing_soft_clip += op_len;
            skip_count += 1;
        } else {
            break;
        }
    }

    let post_clip_end = cigar_ops.len() - skip_count;
    let post_clip_ops = &cigar_ops[..post_clip_end];

    let mut read_bases_clipped = 0usize;
    let mut new_ops: Vec<u32> = Vec::new();
    let mut idx = post_clip_ops.len();

    // Clip operations from the end (working backwards)
    while idx > 0 {
        let op = post_clip_ops[idx - 1];
        let op_type = op & 0xF;
        let op_len = (op >> 4) as usize;

        // Check if we've clipped enough but need to skip adjacent deletions
        if read_bases_clipped == clip_amount && new_ops.is_empty() && op_type == 2 {
            // Deletion: skip it
            idx -= 1;
            continue;
        }

        if read_bases_clipped >= clip_amount {
            break;
        }

        // Note: S (soft clip) is NOT read-consuming here, matching the Clipper.
        let is_read = consumes_read(op_type); // M, I, =, X

        if is_read && op_len > (clip_amount - read_bases_clipped) {
            if op_type == 1 {
                // Insertion: consume entire at clip boundary
                read_bases_clipped += op_len;
            } else {
                // Split the operation
                let remaining_clip = clip_amount - read_bases_clipped;
                let remaining_length = op_len - remaining_clip;
                read_bases_clipped += remaining_clip;
                new_ops.push(encode_op(op_type, remaining_length));
            }
        } else if is_read {
            read_bases_clipped += op_len;
        }

        idx -= 1;
    }

    // Add remaining operations (new_ops is in reverse, remaining are forward)
    let remaining: Vec<u32> = post_clip_ops[..idx].to_vec();
    let mut result = remaining;
    // new_ops collected in reverse order, need to reverse
    new_ops.reverse();
    result.extend(new_ops);

    // Hard clip mode: convert all existing soft clips to hard clips
    let added_hard_clip = existing_soft_clip + read_bases_clipped;
    let total_hard_clip = existing_hard_clip + added_hard_clip;
    result.push(encode_op(5, total_hard_clip)); // H

    (result, 0) // ref_bases_consumed is 0 for end clipping (no position adjustment)
}

// ============================================================================
// CigarKind and CigarOp — typed, zero-cost CIGAR iteration
// ============================================================================

/// BAM CIGAR operation kind. Matches the wire-format op-type low 4 bits.
///
/// Variant ordering is intentionally identical to the BAM spec so that
/// `as u8` gives the correct op-type byte (useful in low-level code).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum CigarKind {
    /// Alignment match (M, op-type 0). Consumes both query and reference.
    Match = 0,
    /// Insertion to reference (I, op-type 1). Consumes query only.
    Insertion = 1,
    /// Deletion from reference (D, op-type 2). Consumes reference only.
    Deletion = 2,
    /// Skipped region from reference (N, op-type 3). Consumes reference only.
    Skip = 3,
    /// Soft clip (S, op-type 4). Consumes query only.
    SoftClip = 4,
    /// Hard clip (H, op-type 5). Consumes neither.
    HardClip = 5,
    /// Padding (P, op-type 6). Consumes neither.
    Pad = 6,
    /// Sequence match (=, op-type 7). Consumes both query and reference.
    SequenceMatch = 7,
    /// Sequence mismatch (X, op-type 8). Consumes both query and reference.
    SequenceMismatch = 8,
}

impl CigarKind {
    /// Single-char SAM representation (M/I/D/N/S/H/P/=/X).
    #[inline]
    #[must_use]
    pub const fn as_char(self) -> char {
        match self {
            Self::Match => 'M',
            Self::Insertion => 'I',
            Self::Deletion => 'D',
            Self::Skip => 'N',
            Self::SoftClip => 'S',
            Self::HardClip => 'H',
            Self::Pad => 'P',
            Self::SequenceMatch => '=',
            Self::SequenceMismatch => 'X',
        }
    }

    /// Parse from the wire-format 4-bit op-type value.
    ///
    /// Returns `None` for values >= 9 (reserved/invalid in the BAM spec).
    #[inline]
    #[must_use]
    pub const fn from_op_type(v: u32) -> Option<Self> {
        match v {
            0 => Some(Self::Match),
            1 => Some(Self::Insertion),
            2 => Some(Self::Deletion),
            3 => Some(Self::Skip),
            4 => Some(Self::SoftClip),
            5 => Some(Self::HardClip),
            6 => Some(Self::Pad),
            7 => Some(Self::SequenceMatch),
            8 => Some(Self::SequenceMismatch),
            _ => None,
        }
    }
}

/// A single typed BAM CIGAR operation — `Copy` newtype over the raw u32 packing.
///
/// Provides the same ergonomic shape as noodles `record::cigar::op::Op` but
/// backed by a single `u32` (no allocation, no trait dispatch, inline methods).
///
/// # Wire format
///
/// BAM encodes each CIGAR op as a u32 where:
/// - low 4 bits (bits 0–3) encode the op type (0 = M, 1 = I, …, 8 = X)
/// - high 28 bits (bits 4–31) encode the op length
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct CigarOp(u32);

impl CigarOp {
    /// Construct from the raw 4-bit-op-type + 28-bit-length packing.
    #[inline]
    #[must_use]
    pub const fn from_raw(raw: u32) -> Self {
        Self(raw)
    }

    /// Construct from a (kind, len) pair.
    ///
    /// Rejects zero-length ops: the BAM spec treats them as invalid (see
    /// [`CigarOp::len`]), and the typed constructor should not mint what the
    /// type claims to forbid. Use [`CigarOp::from_raw`] on a raw `u32` word if
    /// you must preserve malformed wire data.
    ///
    /// # Panics
    ///
    /// - `len == 0`
    /// - `len > 0x0FFF_FFFF` (28-bit limit)
    #[inline]
    #[must_use]
    pub fn new(kind: CigarKind, len: u32) -> Self {
        assert!(len > 0, "CIGAR op length must be non-zero");
        assert!(len <= 0x0FFF_FFFF, "CIGAR op length exceeds 28-bit limit");
        Self((len << 4) | (kind as u32))
    }

    /// The raw u32 packing, for byte-level code paths.
    #[inline]
    #[must_use]
    pub const fn raw(self) -> u32 {
        self.0
    }

    /// Op kind (M, I, D, N, S, H, P, =, X).
    ///
    /// Decodes the low 4 bits of the raw u32 and returns the corresponding
    /// [`CigarKind`] variant. Returns [`CigarKind::Pad`] for reserved/invalid
    /// op-type values (>= 9) — same fallback used by [`cigar_op_kind`].
    #[inline]
    #[must_use]
    pub fn kind(self) -> CigarKind {
        CigarKind::from_op_type(self.0 & 0xF).unwrap_or(CigarKind::Pad)
    }

    /// Op length (number of bases).
    ///
    /// Note: zero-length CIGAR ops are invalid per the BAM spec, so there is
    /// intentionally no companion `is_empty` method.
    #[allow(clippy::len_without_is_empty)]
    #[inline]
    #[must_use]
    pub const fn len(self) -> u32 {
        self.0 >> 4
    }

    /// True if this op type consumes query bases (M, I, S, =, X).
    ///
    /// Uses the `BAM_CIGAR_TYPE` bitmask for a branchless test.
    #[inline]
    #[must_use]
    pub fn consumes_query(self) -> bool {
        consumes_query(self.0 & 0xF)
    }

    /// True if this op type consumes reference bases (M, D, N, =, X).
    ///
    /// Uses the `BAM_CIGAR_TYPE` bitmask for a branchless test.
    #[inline]
    #[must_use]
    pub fn consumes_ref(self) -> bool {
        consumes_ref(self.0 & 0xF)
    }
}

#[cfg(feature = "noodles")]
impl From<CigarKind> for noodles::sam::alignment::record::cigar::op::Kind {
    fn from(k: CigarKind) -> Self {
        use noodles::sam::alignment::record::cigar::op::Kind;
        match k {
            CigarKind::Match => Kind::Match,
            CigarKind::Insertion => Kind::Insertion,
            CigarKind::Deletion => Kind::Deletion,
            CigarKind::Skip => Kind::Skip,
            CigarKind::SoftClip => Kind::SoftClip,
            CigarKind::HardClip => Kind::HardClip,
            CigarKind::Pad => Kind::Pad,
            CigarKind::SequenceMatch => Kind::SequenceMatch,
            CigarKind::SequenceMismatch => Kind::SequenceMismatch,
        }
    }
}

#[cfg(feature = "noodles")]
impl From<noodles::sam::alignment::record::cigar::op::Kind> for CigarKind {
    fn from(k: noodles::sam::alignment::record::cigar::op::Kind) -> Self {
        use noodles::sam::alignment::record::cigar::op::Kind;
        match k {
            Kind::Match => CigarKind::Match,
            Kind::Insertion => CigarKind::Insertion,
            Kind::Deletion => CigarKind::Deletion,
            Kind::Skip => CigarKind::Skip,
            Kind::SoftClip => CigarKind::SoftClip,
            Kind::HardClip => CigarKind::HardClip,
            Kind::Pad => CigarKind::Pad,
            Kind::SequenceMatch => CigarKind::SequenceMatch,
            Kind::SequenceMismatch => CigarKind::SequenceMismatch,
        }
    }
}

use crate::fields::RawRecordView;

impl<'a> RawRecordView<'a> {
    /// Raw CIGAR bytes (4 bytes per op), zero-allocation.
    #[inline]
    #[must_use]
    pub fn cigar_raw_bytes(&self) -> &'a [u8] {
        let bam = self.as_bytes();
        let lrn = l_read_name(bam) as usize;
        let n = n_cigar_op(bam) as usize;
        let start = 32 + lrn;
        let end = start + n * 4;
        if end <= bam.len() { &bam[start..end] } else { &[] }
    }

    /// Zero-allocation iterator yielding decoded CIGAR ops as raw `u32`.
    #[inline]
    pub fn cigar_ops_iter(&self) -> impl Iterator<Item = u32> + 'a {
        self.cigar_raw_bytes().chunks_exact(4).map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
    }

    /// Typed iteration over CIGAR ops. Ergonomically equivalent to
    /// noodles `cigar().iter()` but zero-cost: each [`CigarOp`] is a `Copy`
    /// newtype over the raw u32.
    ///
    /// The returned iterator has lifetime `'a` (tied to the underlying byte
    /// slice, not to `&self`), so it can outlive a temporary view.
    #[inline]
    pub fn cigar_ops_typed(&self) -> impl Iterator<Item = CigarOp> + 'a {
        self.cigar_ops_iter().map(CigarOp::from_raw)
    }

    /// Convenience: collect CIGAR ops into a `Vec<u32>`.
    #[inline]
    #[must_use]
    pub fn cigar_ops_vec(&self) -> Vec<u32> {
        get_cigar_ops(self.as_bytes())
    }

    /// Format CIGAR as a SAM-style string.
    #[inline]
    #[must_use]
    pub fn cigar_to_string(&self) -> String {
        cigar_to_string_from_raw(self.as_bytes())
    }

    /// Sum of M/D/N/=/X op lengths (reference-consuming).
    #[inline]
    #[must_use]
    pub fn reference_length(&self) -> i32 {
        reference_length_from_raw_bam(self.as_bytes())
    }

    /// Sum of M/I/S/=/X op lengths (query-consuming).
    #[inline]
    #[must_use]
    pub fn query_length(&self) -> usize {
        self.cigar_ops_iter()
            .filter(|&op| consumes_query(op & 0xF))
            .map(|op| (op >> 4) as usize)
            .sum()
    }

    /// Compute reference-consuming and query-consuming lengths in one pass.
    ///
    /// Equivalent to `(self.reference_length(), self.query_length())` but
    /// walks the CIGAR only once. Use this when both values are needed
    /// (e.g. record validation, alignment-end computation).
    ///
    /// Returns `(ref_length, query_length)`.
    #[inline]
    #[must_use]
    pub fn cigar_lengths(&self) -> (i32, usize) {
        let mut ref_len = 0i32;
        let mut q_len = 0usize;
        for op in self.cigar_ops_iter() {
            let ty = op & 0xF;
            let len_u32 = op >> 4;
            if consumes_ref(ty) {
                ref_len += len_u32.cast_signed();
            }
            if consumes_query(ty) {
                q_len += len_u32 as usize;
            }
        }
        (ref_len, q_len)
    }

    /// Compute 1-based alignment end position.
    #[inline]
    #[must_use]
    pub fn alignment_end_1based(&self) -> Option<usize> {
        alignment_end_from_raw(self.as_bytes())
    }

    /// Compute 1-based alignment start position.
    #[inline]
    #[must_use]
    pub fn alignment_start_1based(&self) -> Option<usize> {
        alignment_start_from_raw(self.as_bytes())
    }

    /// Computes the unclipped 5' coordinate (1-based) for grouping keys.
    #[inline]
    #[must_use]
    pub fn unclipped_5prime_1based(&self) -> i32 {
        unclipped_5prime_from_raw_bam(self.as_bytes())
    }

    /// Virtual clipping — see [`clip_cigar_ops_raw`].
    #[inline]
    #[must_use]
    pub fn clip_cigar_ops(&self, clip_amount: usize, from_start: bool) -> (Vec<u32>, usize) {
        clip_cigar_ops_raw(&self.cigar_ops_vec(), clip_amount, from_start)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::*;

    // ========================================================================
    // BAM_CIGAR_TYPE bitmask tests
    // ========================================================================

    #[test]
    fn test_cigar_type_bitmask_matches_spec() {
        // Spec: (consumes_query, consumes_ref) for each op type 0..8
        let expected = [
            (true, true),   // 0: M
            (true, false),  // 1: I
            (false, true),  // 2: D
            (false, true),  // 3: N
            (true, false),  // 4: S
            (false, false), // 5: H
            (false, false), // 6: P
            (true, true),   // 7: =
            (true, true),   // 8: X
        ];
        for (op, &(eq, er)) in expected.iter().enumerate() {
            let op = u32::try_from(op).unwrap();
            assert_eq!(consumes_query(op), eq, "consumes_query({op})");
            assert_eq!(consumes_ref(op), er, "consumes_ref({op})");
        }
        // Op values >= 9 should return false for both
        for op in 9u32..=15 {
            assert!(!consumes_query(op), "consumes_query({op}) should be false");
            assert!(!consumes_ref(op), "consumes_ref({op}) should be false");
        }
    }

    #[test]
    fn test_consumes_fns_mask_values_beyond_4_bits() {
        // Without the `& 0xF` guard, op_type >= 16 overshifts (debug panic /
        // release wrap). With the guard, the input is normalised to 0..=15 and
        // the result matches the low-nibble equivalent. For example op=16 aliases
        // to op=0 (M) and op=22 aliases to op=6 (P).
        for (bad, good) in [(16u32, 0u32), (22u32, 6u32), (257u32, 1u32), (u32::MAX, 15u32)] {
            assert_eq!(
                consumes_query(bad),
                consumes_query(good),
                "consumes_query({bad}) should equal consumes_query({good})",
            );
            assert_eq!(
                consumes_ref(bad),
                consumes_ref(good),
                "consumes_ref({bad}) should equal consumes_ref({good})",
            );
        }
    }

    // ========================================================================
    // unclipped_start_from_cigar / unclipped_end_from_cigar tests
    // ========================================================================

    // ========================================================================
    // query_length_from_cigar tests
    // ========================================================================

    #[test]
    fn test_query_length_from_cigar_simple_match() {
        // 50M: query length = 50
        let cigar = &[encode_op(0, 50)];
        assert_eq!(query_length_from_cigar(cigar), 50);
    }

    #[test]
    fn test_query_length_from_cigar_with_insertion() {
        // 10M5I10M: M and I consume query
        // query_len = 10 + 5 + 10 = 25
        let cigar = &[encode_op(0, 10), encode_op(1, 5), encode_op(0, 10)];
        assert_eq!(query_length_from_cigar(cigar), 25);
    }

    #[test]
    fn test_query_length_from_cigar_with_soft_clip() {
        // 5S10M3S: S and M consume query
        // query_len = 5 + 10 + 3 = 18
        let cigar = &[encode_op(4, 5), encode_op(0, 10), encode_op(4, 3)];
        assert_eq!(query_length_from_cigar(cigar), 18);
    }

    #[test]
    fn test_query_length_from_cigar_with_deletion() {
        // 10M3D5M: D does not consume query
        // query_len = 10 + 5 = 15
        let cigar = &[encode_op(0, 10), encode_op(2, 3), encode_op(0, 5)];
        assert_eq!(query_length_from_cigar(cigar), 15);
    }

    #[test]
    fn test_query_length_from_cigar_with_eq_and_x() {
        // 10=3X2I: =, X, and I consume query
        // query_len = 10 + 3 + 2 = 15
        let cigar = &[encode_op(7, 10), encode_op(8, 3), encode_op(1, 2)];
        assert_eq!(query_length_from_cigar(cigar), 15);
    }

    #[test]
    fn test_query_length_from_cigar_empty() {
        assert_eq!(query_length_from_cigar(&[]), 0);
    }

    #[test]
    fn test_query_length_from_cigar_hard_clip_only() {
        // 5H: hard clips do not consume query
        let cigar = &[encode_op(5, 5)];
        assert_eq!(query_length_from_cigar(cigar), 0);
    }

    // ========================================================================
    // read_pos_at_ref_pos_raw tests
    // ========================================================================

    #[test]
    fn test_read_pos_at_ref_pos_raw_simple_match() {
        // 10M starting at ref pos 100 (1-based)
        let cigar = &[encode_op(0, 10)];
        // Position 100: first query base
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 100, false), Some(1));
        // Position 105: 6th query base
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 105, false), Some(6));
        // Position 109: last query base
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 109, false), Some(10));
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_before_alignment() {
        let cigar = &[encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 99, false), None);
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_past_alignment() {
        let cigar = &[encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 110, false), None);
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_in_deletion() {
        // 5M3D5M: deletion at ref 105-107
        let cigar = &[encode_op(0, 5), encode_op(2, 3), encode_op(0, 5)];
        // Position 106 is in deletion
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 106, false), None);
        // With return_last_base_if_deleted=true
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 106, true), Some(5));
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_with_insertion() {
        // 5M3I5M: insertion between ref 104 and 105
        let cigar = &[encode_op(0, 5), encode_op(1, 3), encode_op(0, 5)];
        // Position 104: 5th query base (before insertion)
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 104, false), Some(5));
        // Position 105: 9th query base (after 5M + 3I)
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 105, false), Some(9));
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_with_soft_clip() {
        // 3S10M: soft clip consumes 3 query bases, then 10M
        let cigar = &[encode_op(4, 3), encode_op(0, 10)];
        // Position 100: 4th query base (after 3S)
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 100, false), Some(4));
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_deletion_at_start() {
        // 2D5M: deletion at ref 100-101
        let cigar = &[encode_op(2, 2), encode_op(0, 5)];
        // Position 100 is in the deletion with no prior query bases
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 100, true), Some(1));
    }

    // ========================================================================
    // mate_unclipped_5prime_1based tests
    // ========================================================================

    #[test]
    fn test_mate_unclipped_5prime_1based_forward() {
        // MC=5S10M: forward 5' = pos+1 - 5 = 96
        assert_eq!(mate_unclipped_5prime_1based(100, false, "5S10M"), 96);
    }

    #[test]
    fn test_mate_unclipped_5prime_1based_reverse() {
        // MC=10M5S: reverse 5' = pos+1 + 10 + 5 - 1 = 115
        assert_eq!(mate_unclipped_5prime_1based(100, true, "10M5S"), 115);
    }

    // ========================================================================
    // parse_leading_clips tests
    // ========================================================================

    #[test]
    fn test_parse_leading_clips_soft() {
        // "5S10M" -> leading clips = 5
        assert_eq!(parse_leading_clips(b"5S10M"), 5);
    }

    #[test]
    fn test_parse_leading_clips_hard_and_soft() {
        // "3H5S10M" -> leading clips = 3 + 5 = 8
        assert_eq!(parse_leading_clips(b"3H5S10M"), 8);
    }

    #[test]
    fn test_parse_leading_clips_no_clips() {
        // "10M" -> leading clips = 0
        assert_eq!(parse_leading_clips(b"10M"), 0);
    }

    #[test]
    fn test_parse_leading_clips_all_clips() {
        // "5S3H" — all clips, no ref ops
        assert_eq!(parse_leading_clips(b"5S3H"), 8);
    }

    // ========================================================================
    // parse_ref_len_and_trailing_clips tests
    // ========================================================================

    #[test]
    fn test_parse_ref_len_and_trailing_clips_basic() {
        // "10M5S" -> ref_len=10, trailing_clips=5
        assert_eq!(parse_ref_len_and_trailing_clips(b"10M5S"), (10, 5));
    }

    #[test]
    fn test_parse_ref_len_and_trailing_clips_complex() {
        // "5S10M2I3D5M3S2H"
        // ref consuming: 10M + 3D + 5M = 18
        // trailing clips: 3S + 2H = 5
        assert_eq!(parse_ref_len_and_trailing_clips(b"5S10M2I3D5M3S2H"), (18, 5));
    }

    /// Malformed MC values must not panic and must not produce a negative
    /// coordinate.
    ///
    /// MC is a *text* tag, so its run lengths are not bounded by BAM's 28-bit
    /// CIGAR op field, and its bytes are not guaranteed to be ASCII. Both
    /// parsers are reachable in release builds from `fgumi sort`'s template
    /// key (`external.rs`) and from the decode step's `GroupKey` construction.
    ///
    /// Expectations are literal, not computed from the previous implementation
    /// -- for these inputs the old code panicked or wrapped, so it cannot serve
    /// as an oracle.
    #[rstest]
    // Non-ASCII but valid UTF-8: `str::from_utf8` admits it, so it reaches the
    // parser. Byte 4 falls inside 'e-acute', which panicked when the old code
    // sliced the &str there.
    #[case::non_ascii_utf8("10M\u{e9}", (10, 0))]
    #[case::non_ascii_leading("\u{e9}10M", (10, 0))]
    // Unknown ASCII operators are ignored, as before.
    #[case::unknown_operator("10Mfoo5S", (10, 5))]
    #[case::placeholder_star("*", (0, 0))]
    // A single run wider than i32: `str::parse::<i32>` returns Err, and
    // `unwrap_or(0)` makes it 0. Saturation would wrongly give i32::MAX.
    #[case::run_overflows_i32("99999999999M", (0, 0))]
    // Totals that overflow across operators must clamp, not wrap to negative.
    #[case::sum_overflows_i32("2000000000M2000000000M", (i32::MAX, 0))]
    #[case::trailing_sum_overflows("1M2000000000S2000000000S", (1, i32::MAX))]
    // A non-ref, non-clip operator between the last ref op and the trailing
    // clips must not reset the run (matches htsjdk's getUnclippedEnd).
    #[case::insertion_before_trailing_clip("10M5I3S", (10, 3))]
    #[case::padding_before_trailing_clip("10M5P3S", (10, 3))]
    // An unknown operator AFTER a trailing clip must not clear the run --
    // only a reference-consuming operator resets it.
    #[case::garbage_after_trailing_clip("10M3Sfoo", (10, 3))]
    // A run too wide for u64 itself: the accumulator must saturate, not wrap
    // back into i32 range. 18446744073709551617 == u64::MAX + 2.
    #[case::run_overflows_u64("18446744073709551617M", (0, 0))]
    #[case::leading_zeros("0007M", (7, 0))]
    #[case::empty("", (0, 0))]
    #[case::digits_only("10", (0, 0))]
    fn test_parse_ref_len_and_trailing_clips_malformed(
        #[case] cigar: &str,
        #[case] expected: (i32, i32),
    ) {
        assert_eq!(parse_ref_len_and_trailing_clips(cigar.as_bytes()), expected);
    }

    /// Same contract for the leading-clip parser.
    #[rstest]
    #[case::non_ascii_utf8("\u{e9}10M", 0)]
    #[case::run_overflows_i32("99999999999S", 0)]
    #[case::run_overflows_u64("18446744073709551617S", 0)]
    #[case::sum_overflows_i32("2000000000S2000000000S", i32::MAX)]
    #[case::leading_zeros("0007S", 7)]
    #[case::empty("", 0)]
    #[case::digits_only("10", 0)]
    #[case::stops_at_first_non_clip("5S10M5S", 5)]
    fn test_parse_leading_clips_malformed(#[case] cigar: &str, #[case] expected: i32) {
        assert_eq!(parse_leading_clips(cigar.as_bytes()), expected);
    }

    /// One token of a generated MC value: a run of digits, a CIGAR operator, a
    /// whole non-ASCII character, or an arbitrary byte.
    ///
    /// Generating *tokens* rather than independent bytes is what makes this
    /// property worth running, for two reasons.
    ///
    /// A per-byte digit generator has to emit ~20 consecutive digits by chance
    /// to reach the `u64` accumulator's ceiling and ~10 to reach `i32`'s, which
    /// effectively never happens -- it would exercise little beyond the
    /// "unknown byte" arm. Drawing a run length up to 24 instead puts both
    /// overflow ceilings, and the saturating totals that span several
    /// operators, in reach on nearly every case.
    ///
    /// The multi-byte-character arm matters for the same reason in the other
    /// direction: a *lone* high byte is not valid UTF-8, so an `any::<u8>()`
    /// arm alone can essentially never build a value that survives
    /// `str::from_utf8` while still being non-ASCII -- exactly the shape that
    /// panicked the old `&str`-slicing parsers. Emitting a whole `char`'s
    /// encoding keeps that case reachable. The `any::<u8>()` arm still covers
    /// arbitrary and invalid-UTF-8 bytes for the byte-slice parsers.
    fn cigar_ish_token() -> impl proptest::strategy::Strategy<Value = Vec<u8>> {
        use proptest::strategy::Strategy as _;
        proptest::prop_oneof![
            6 => proptest::collection::vec(b'0'..=b'9', 1..24),
            6 => proptest::sample::select(vec![b'M', b'I', b'D', b'N', b'S', b'H', b'P', b'=', b'X'])
                .prop_map(|op| vec![op]),
            3 => proptest::prelude::any::<char>()
                .prop_map(|ch| ch.to_string().into_bytes()),
            2 => proptest::prelude::any::<u8>().prop_map(|byte| vec![byte]),
        ]
    }

    /// An MC-tag-shaped byte string built from [`cigar_ish_token`].
    fn cigar_ish_bytes() -> impl proptest::strategy::Strategy<Value = Vec<u8>> {
        use proptest::strategy::Strategy as _;
        proptest::collection::vec(cigar_ish_token(), 0..12).prop_map(|tokens| tokens.concat())
    }

    proptest::proptest! {
        /// The parsers must tolerate the full input space, not just the cases
        /// the tables above name.
        ///
        /// MC is a text tag copied verbatim from another record, so the bytes
        /// reaching these parsers are untrusted: the contract is that *any*
        /// byte string returns non-negative lengths without panicking or
        /// wrapping. The tables pin specific expectations; this pins the
        /// contract itself.
        #[test]
        fn test_parse_leading_clips_never_panics_or_goes_negative(
            cigar in cigar_ish_bytes(),
        ) {
            proptest::prop_assert!(parse_leading_clips(&cigar) >= 0);
        }

        #[test]
        fn test_parse_ref_len_and_trailing_clips_never_panics_or_goes_negative(
            cigar in cigar_ish_bytes(),
        ) {
            let (ref_len, trailing_clips) = parse_ref_len_and_trailing_clips(&cigar);
            proptest::prop_assert!(ref_len >= 0);
            proptest::prop_assert!(trailing_clips >= 0);
        }

        /// The public entry point must not panic either, on any MC that
        /// survives `str::from_utf8` -- the filter these values pass through.
        #[test]
        fn test_mate_unclipped_5prime_never_panics(
            cigar in cigar_ish_bytes(),
            mate_pos in proptest::prelude::any::<i32>(),
            mate_reverse in proptest::prelude::any::<bool>(),
        ) {
            if let Ok(mc) = std::str::from_utf8(&cigar) {
                // Only asserting it returns: the panic and the debug-build
                // overflow are what this is guarding against.
                let _ = mate_unclipped_5prime(mate_pos, mate_reverse, mc);
                let _ = mate_unclipped_5prime_1based(mate_pos, mate_reverse, mc);
            }
        }
    }

    /// The public entry points must not panic on a non-ASCII MC either.
    #[rstest]
    #[case::forward(false, 96)]
    #[case::reverse(true, 110)]
    fn test_mate_unclipped_5prime_1based_non_ascii_mc(
        #[case] mate_reverse: bool,
        #[case] expected: i32,
    ) {
        // "5S10M\u{e9}": leading clips 5, ref_len 10, no trailing clips.
        assert_eq!(mate_unclipped_5prime_1based(100, mate_reverse, "5S10M\u{e9}"), expected);
    }

    #[test]
    fn test_parse_ref_len_and_trailing_clips_no_ref_ops() {
        // "5S3H" — no ref-consuming ops
        assert_eq!(parse_ref_len_and_trailing_clips(b"5S3H"), (0, 0));
    }

    #[test]
    fn test_parse_ref_len_and_trailing_clips_eq_and_x() {
        // "10=3X" → ref_len=13, trailing_clips=0
        assert_eq!(parse_ref_len_and_trailing_clips(b"10=3X"), (13, 0));
    }

    // ========================================================================
    // reference_length_from_cigar tests
    // ========================================================================

    #[test]
    fn test_reference_length_simple_match() {
        // 50M
        let cigar = &[(50 << 4)];
        assert_eq!(reference_length_from_cigar(cigar), 50);
    }

    #[test]
    fn test_reference_length_with_deletions() {
        // 10M3D5M2N8M
        // M (0), D (2), N (3) all consume reference
        // ref_len = 10 + 3 + 5 + 2 + 8 = 28
        let cigar = &[
            (10 << 4),    // 10M
            (3 << 4) | 2, // 3D
            (5 << 4),     // 5M
            (2 << 4) | 3, // 2N
            (8 << 4),     // 8M
        ];
        assert_eq!(reference_length_from_cigar(cigar), 28);
    }

    #[test]
    fn test_reference_length_with_insertions() {
        // 10M5I10M: insertions don't consume reference
        // ref_len = 10 + 10 = 20
        let cigar = &[
            (10 << 4),    // 10M
            (5 << 4) | 1, // 5I
            (10 << 4),    // 10M
        ];
        assert_eq!(reference_length_from_cigar(cigar), 20);
    }

    // ========================================================================
    // reference_length_from_raw_bam tests
    // ========================================================================

    #[test]
    fn test_reference_length_from_raw_bam_no_cigar() {
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(reference_length_from_raw_bam(&rec), 0);
    }

    #[test]
    fn test_reference_length_from_raw_bam_simple() {
        // 50M
        let cigar = &[(50 << 4)];
        let rec = make_bam_bytes(0, 100, 0, b"rea", cigar, 50, -1, -1, &[]);
        assert_eq!(reference_length_from_raw_bam(&rec), 50);
    }

    #[test]
    fn test_reference_length_from_raw_bam_with_insertions() {
        // 10M5I10M: insertions don't consume reference, ref_len = 20
        let cigar = &[(10 << 4), (5 << 4) | 1, (10 << 4)];
        let rec = make_bam_bytes(0, 100, 0, b"rea", cigar, 25, -1, -1, &[]);
        assert_eq!(reference_length_from_raw_bam(&rec), 20);
    }

    #[test]
    fn test_reference_length_from_raw_bam_truncated_cigar() {
        let mut rec = make_bam_bytes(0, 100, 0, b"rea", &[(10 << 4)], 10, -1, -1, &[]);
        rec[12..14].copy_from_slice(&100u16.to_le_bytes());
        assert_eq!(reference_length_from_raw_bam(&rec), 0);
    }

    #[test]
    fn test_reference_length_from_raw_bam_checked_no_cigar_is_valid_zero() {
        // A record that genuinely has no CIGAR is well-formed with a zero span.
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(reference_length_from_raw_bam_checked(&rec), Some(0));
    }

    #[test]
    fn test_reference_length_from_raw_bam_checked_simple() {
        let cigar = &[(50 << 4)]; // 50M
        let rec = make_bam_bytes(0, 100, 0, b"rea", cigar, 50, -1, -1, &[]);
        assert_eq!(reference_length_from_raw_bam_checked(&rec), Some(50));
    }

    #[test]
    fn test_reference_length_from_raw_bam_checked_truncated_cigar_is_none() {
        // n_cigar_op declares 100 ops but the buffer holds only one: malformed,
        // must report `None` rather than the `0` the infallible variant returns.
        let mut rec = make_bam_bytes(0, 100, 0, b"rea", &[(10 << 4)], 10, -1, -1, &[]);
        rec[12..14].copy_from_slice(&100u16.to_le_bytes());
        assert_eq!(reference_length_from_raw_bam_checked(&rec), None);
        // The infallible variant still collapses the malformed record to `0`.
        assert_eq!(reference_length_from_raw_bam(&rec), 0);
    }

    #[test]
    fn test_reference_length_from_raw_bam_checked_too_short_is_none() {
        // Shorter than the 32-byte fixed core: cannot be a valid record.
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[(10 << 4)], 10, -1, -1, &[]);
        assert_eq!(reference_length_from_raw_bam_checked(&rec[..20]), None);
    }

    #[test]
    fn test_reference_length_from_raw_bam_checked_zero_cigar_oversized_read_name_is_none() {
        // Core-only record with no CIGAR but l_read_name corrupted to 255: the
        // declared read name overruns the buffer, so it must fail closed rather than
        // be accepted as a zero-length point alignment.
        let mut rec = make_bam_bytes(0, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        rec[8] = 255; // l_read_name (byte 8) -> 255, far past the 36-byte buffer
        assert_eq!(reference_length_from_raw_bam_checked(&rec), None);
        assert_eq!(reference_length_from_raw_bam(&rec), 0);
    }

    #[test]
    fn test_reference_length_from_raw_bam_checked_overflow_is_none() {
        // Nine ref-consuming M ops of length 2^28-1 each sum to ~2.4B, past
        // i32::MAX (~2.1B). The CIGAR is in-bounds but malformed: the accumulator
        // must fail closed (None) rather than wrap negative (or panic in debug).
        let big_m = 0x0FFF_FFFF_u32 << 4; // length = 2^28-1, op = M (0, consumes ref)
        let cigar = vec![big_m; 9];
        let rec = make_bam_bytes(0, 100, 0, b"rea", &cigar, 0, -1, -1, &[]);
        assert_eq!(reference_length_from_raw_bam_checked(&rec), None);
        // The infallible wrapper collapses the overflow to 0 like any other malform.
        assert_eq!(reference_length_from_raw_bam(&rec), 0);
    }

    #[test]
    fn test_reference_length_from_raw_bam_mixed_ops() {
        // 5S10M2D3N5M3I8X2H: ref consuming = 10M+2D+3N+5M+8X = 28
        let cigar = &[
            (5 << 4) | 4, // 5S
            (10 << 4),    // 10M
            (2 << 4) | 2, // 2D
            (3 << 4) | 3, // 3N
            (5 << 4),     // 5M
            (3 << 4) | 1, // 3I
            (8 << 4) | 8, // 8X (= mismatch, consumes ref)
            (2 << 4) | 5, // 2H
        ];
        let rec = make_bam_bytes(0, 100, 0, b"rea", cigar, 31, -1, -1, &[]);
        assert_eq!(reference_length_from_raw_bam(&rec), 28);
    }

    // ========================================================================
    // unclipped_5prime_from_raw_bam tests
    // ========================================================================

    #[test]
    fn test_unclipped_5prime_from_raw_bam_unmapped() {
        let rec = make_bam_bytes(0, 100, flags::UNMAPPED, b"rea", &[(10 << 4)], 10, -1, -1, &[]);
        assert_eq!(unclipped_5prime_from_raw_bam(&rec), 0);
    }

    #[test]
    fn test_unclipped_5prime_from_raw_bam_no_cigar() {
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(unclipped_5prime_from_raw_bam(&rec), i32::MAX);
    }

    #[test]
    fn test_unclipped_5prime_from_raw_bam_forward() {
        // 5S10M: forward, pos_0based=100 → 1-based=101, unclipped = 101 - 5 = 96
        let cigar = &[(5 << 4) | 4, (10 << 4)]; // 5S10M
        let rec = make_bam_bytes(0, 100, 0, b"rea", cigar, 15, -1, -1, &[]);
        assert_eq!(unclipped_5prime_from_raw_bam(&rec), 96);
    }

    #[test]
    fn test_unclipped_5prime_from_raw_bam_reverse() {
        // 10M5S: reverse, pos_0based=100 → 1-based=101, end = 101 + 10 + 5 - 1 = 115
        let cigar = &[(10 << 4), (5 << 4) | 4]; // 10M5S
        let rec = make_bam_bytes(0, 100, flags::REVERSE, b"rea", cigar, 15, -1, -1, &[]);
        assert_eq!(unclipped_5prime_from_raw_bam(&rec), 115);
    }

    #[test]
    fn test_unclipped_5prime_from_raw_bam_truncated_cigar() {
        // Record with n_cigar_op > 0 but record too short to hold CIGAR data
        let mut rec = make_bam_bytes(0, 100, 0, b"rea", &[(10 << 4)], 10, -1, -1, &[]);
        rec[12..14].copy_from_slice(&100u16.to_le_bytes()); // claim 100 CIGAR ops
        assert_eq!(unclipped_5prime_from_raw_bam(&rec), i32::MAX);
    }

    #[test]
    fn test_unclipped_5prime_from_raw_bam_hard_and_soft_clip() {
        // 3H5S10M4S2H: forward → start = 101 - 3 - 5 = 93
        let cigar = &[
            (3 << 4) | 5, // 3H
            (5 << 4) | 4, // 5S
            (10 << 4),    // 10M
            (4 << 4) | 4, // 4S
            (2 << 4) | 5, // 2H
        ];
        let rec = make_bam_bytes(0, 100, 0, b"rea", cigar, 15, -1, -1, &[]);
        // Forward: pos_1based=101, leading clips=3+5=8, unclipped_start = 101 - 8 = 93
        assert_eq!(unclipped_5prime_from_raw_bam(&rec), 93);
    }

    #[test]
    fn test_unclipped_5prime_from_raw_bam_reverse_hard_and_soft_clip() {
        // 3H5S10M4S2H: reverse → end = 101 + 10 + 4 + 2 - 1 = 116
        let cigar = &[
            (3 << 4) | 5, // 3H
            (5 << 4) | 4, // 5S
            (10 << 4),    // 10M
            (4 << 4) | 4, // 4S
            (2 << 4) | 5, // 2H
        ];
        let rec = make_bam_bytes(0, 100, flags::REVERSE, b"rea", cigar, 15, -1, -1, &[]);
        // Reverse: pos_1based=101, ref_len=10, trailing_clips=4+2=6, end = 101 + 10 + 6 - 1 = 116
        assert_eq!(unclipped_5prime_from_raw_bam(&rec), 116);
    }

    // ========================================================================
    // get_cigar_ops tests
    // ========================================================================

    #[test]
    fn test_get_cigar_ops_no_cigar() {
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        assert!(get_cigar_ops(&rec).is_empty());
    }

    #[test]
    fn test_get_cigar_ops_truncated_record() {
        // Record claims n_cigar_op=5 but is too short to hold that many
        let mut rec = make_bam_bytes(0, 100, 0, b"rea", &[(10 << 4)], 10, -1, -1, &[]);
        // Override n_cigar_op to a large value
        rec[12..14].copy_from_slice(&100u16.to_le_bytes());
        assert!(get_cigar_ops(&rec).is_empty());
    }

    #[test]
    fn test_get_cigar_ops_multiple_ops() {
        let cigar = &[(5 << 4) | 4, (10 << 4), (3 << 4) | 4]; // 5S10M3S
        let rec = make_bam_bytes(0, 100, 0, b"rea", cigar, 18, -1, -1, &[]);
        let ops = get_cigar_ops(&rec);
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0], (5 << 4) | 4);
        assert_eq!(ops[1], (10 << 4));
        assert_eq!(ops[2], (3 << 4) | 4);
    }

    // ========================================================================
    // alignment_start_from_raw / alignment_end_from_raw tests
    // ========================================================================

    #[test]
    fn test_alignment_start_from_raw_mapped() {
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[encode_op(0, 10)], 10, -1, -1, &[]);
        assert_eq!(alignment_start_from_raw(&rec), Some(101));
    }

    #[test]
    fn test_alignment_start_from_raw_unmapped() {
        let rec = make_bam_bytes(-1, -1, flags::UNMAPPED, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(alignment_start_from_raw(&rec), None);
    }

    #[test]
    fn test_alignment_end_from_raw_mapped() {
        // 10M at pos 100 (0-based): end = 100 + 10 = 110 (exclusive)
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[encode_op(0, 10)], 10, -1, -1, &[]);
        assert_eq!(alignment_end_from_raw(&rec), Some(110));
    }

    #[test]
    fn test_alignment_end_from_raw_unmapped() {
        let rec = make_bam_bytes(-1, -1, flags::UNMAPPED, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(alignment_end_from_raw(&rec), None);
    }

    #[test]
    fn test_alignment_end_from_raw_no_cigar() {
        // Mapped but no CIGAR => ref_len=0 => None
        let rec = make_bam_bytes(0, 100, 0, b"rea", &[], 0, -1, -1, &[]);
        assert_eq!(alignment_end_from_raw(&rec), None);
    }

    #[test]
    fn test_alignment_end_from_raw_with_deletion() {
        // 5M2D5M at pos=0 → ref_len = 5+2+5 = 12
        let cigar = &[encode_op(0, 5), encode_op(2, 2), encode_op(0, 5)];
        let rec = make_bam_bytes(0, 0, 0, b"rd", cigar, 10, -1, -1, &[]);
        assert_eq!(alignment_end_from_raw(&rec), Some(12));
    }

    // ========================================================================
    // clip_cigar_ops_raw tests
    // ========================================================================

    #[test]
    fn test_clip_cigar_ops_raw_zero_clip() {
        let cigar = &[encode_op(0, 10)]; // 10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 0, true);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result, cigar);
    }

    #[test]
    fn test_clip_cigar_ops_raw_empty_cigar() {
        let (result, ref_consumed) = clip_cigar_ops_raw(&[], 5, true);
        assert_eq!(ref_consumed, 0);
        assert!(result.is_empty());
    }

    #[test]
    fn test_clip_cigar_ops_raw_upgrade_path_from_start() {
        // 5S10M: clip_amount=3, existing_clip=5, 3<=5 => upgrade path
        let cigar = &[encode_op(4, 5), encode_op(0, 10)]; // 5S10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 3, true);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result[0], encode_op(5, 3)); // 3H
        assert_eq!(result[1], encode_op(4, 2)); // 2S
        assert_eq!(result[2], encode_op(0, 10)); // 10M
    }

    #[test]
    fn test_clip_cigar_ops_raw_upgrade_path_from_end() {
        // 10M5S: clip_amount=3, existing_clip=5, 3<=5 => upgrade path
        let cigar = &[encode_op(0, 10), encode_op(4, 5)]; // 10M5S
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 3, false);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result[0], encode_op(0, 10)); // 10M
        assert_eq!(result[1], encode_op(4, 2)); // 2S
        assert_eq!(result[2], encode_op(5, 3)); // 3H
    }

    #[test]
    fn test_clip_cigar_ops_raw_alignment_clip_from_start() {
        // 10M: clip_amount=3, no existing clips => clips 3 bases from alignment
        let cigar = &[encode_op(0, 10)]; // 10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 3, true);
        assert_eq!(ref_consumed, 3);
        assert_eq!(result.len(), 2); // 3H, 7M
        assert_eq!(result[0], encode_op(5, 3)); // 3H
        assert_eq!(result[1], encode_op(0, 7)); // 7M
    }

    #[test]
    fn test_clip_cigar_ops_raw_alignment_clip_from_end() {
        // 10M: clip_amount=3, no existing clips => clips 3 bases from alignment end
        let cigar = &[encode_op(0, 10)]; // 10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 3, false);
        assert_eq!(ref_consumed, 0); // end-clipping doesn't consume ref for position adjustment
        assert_eq!(result.len(), 2); // 7M, 3H
        assert_eq!(result[0], encode_op(0, 7)); // 7M
        assert_eq!(result[1], encode_op(5, 3)); // 3H
    }

    #[test]
    fn test_clip_cigar_ops_raw_clip_past_existing_from_start() {
        // 2S10M: clip_amount=5, existing_clip=2, alignment_clip=3
        let cigar = &[encode_op(4, 2), encode_op(0, 10)]; // 2S10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 5, true);
        assert_eq!(ref_consumed, 3);
        assert_eq!(result.len(), 2); // 5H, 7M
        assert_eq!(result[0], encode_op(5, 5)); // 5H (2S + 3 from alignment)
        assert_eq!(result[1], encode_op(0, 7)); // 7M
    }

    #[test]
    fn test_clip_cigar_ops_raw_clip_past_existing_from_end() {
        // 10M2S: clip_amount=5, existing_clip=2, alignment_clip=3
        let cigar = &[encode_op(0, 10), encode_op(4, 2)]; // 10M2S
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 5, false);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 2); // 7M, 5H
        assert_eq!(result[0], encode_op(0, 7)); // 7M
        assert_eq!(result[1], encode_op(5, 5)); // 5H (2S + 3 from alignment)
    }

    #[test]
    fn test_clip_cigar_ops_raw_with_insertion_at_boundary_start() {
        // 10M3I5M: clip 10 from start should clip through 10M and consume entire 3I
        let cigar = &[encode_op(0, 10), encode_op(1, 3), encode_op(0, 5)]; // 10M3I5M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 10, true);
        // Clips 10M (10 read bases, 10 ref bases), then hits I:
        // After 10M, read_bases_clipped=10 == clip_amount=10, done
        assert_eq!(ref_consumed, 10);
        assert_eq!(result.len(), 3); // 10H, 3I, 5M
        assert_eq!(result[0], encode_op(5, 10)); // 10H
        assert_eq!(result[1], encode_op(1, 3)); // 3I preserved
        assert_eq!(result[2], encode_op(0, 5)); // 5M preserved
    }

    #[test]
    fn test_clip_cigar_ops_raw_with_deletion_at_boundary_start() {
        // 5M2D10M: clip 5 from start should clip 5M and skip the 2D
        let cigar = &[encode_op(0, 5), encode_op(2, 2), encode_op(0, 10)]; // 5M2D10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 5, true);
        // Clips 5M (5 read + 5 ref), then hits deletion at boundary (read_bases_clipped==5),
        // skip 2D (ref_consumed += 2)
        assert_eq!(ref_consumed, 7); // 5 from M + 2 from D
        assert_eq!(result.len(), 2); // 5H, 10M
        assert_eq!(result[0], encode_op(5, 5)); // 5H
        assert_eq!(result[1], encode_op(0, 10)); // 10M
    }

    #[test]
    fn test_clip_cigar_ops_raw_with_deletion_at_boundary_end() {
        // 10M2D5M: clip 5 from end should clip 5M and skip the 2D
        let cigar = &[encode_op(0, 10), encode_op(2, 2), encode_op(0, 5)]; // 10M2D5M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 5, false);
        // Clips 5M from end (5 read bases), then skip 2D adjacent
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 2); // 10M, 5H
        assert_eq!(result[0], encode_op(0, 10)); // 10M
        assert_eq!(result[1], encode_op(5, 5)); // 5H
    }

    #[test]
    fn test_clip_cigar_ops_raw_split_match_from_start() {
        // 10M: clip 4, splits the M
        let cigar = &[encode_op(0, 10)]; // 10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 4, true);
        assert_eq!(ref_consumed, 4);
        assert_eq!(result.len(), 2); // 4H, 6M
        assert_eq!(result[0], encode_op(5, 4));
        assert_eq!(result[1], encode_op(0, 6));
    }

    #[test]
    fn test_clip_cigar_ops_raw_split_match_from_end() {
        // 10M: clip 4, splits the M from end
        let cigar = &[encode_op(0, 10)]; // 10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 4, false);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 2); // 6M, 4H
        assert_eq!(result[0], encode_op(0, 6));
        assert_eq!(result[1], encode_op(5, 4));
    }

    #[test]
    fn test_clip_cigar_ops_raw_insertion_consumed_at_boundary_start() {
        // 5M3I5M: clip 6 from start
        // First 5M: 5 read bases clipped, 5 ref bases consumed
        // Then 3I: need 1 more, but insertion is consumed entirely
        // read_bases_clipped = 5 + 3 = 8 > 6 (entire I consumed at boundary)
        let cigar = &[encode_op(0, 5), encode_op(1, 3), encode_op(0, 5)]; // 5M3I5M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 6, true);
        // 5M consumed (5 read, 5 ref), then 3I consumed entirely (3 read)
        // total read bases = 8, ref_consumed = 5
        assert_eq!(ref_consumed, 5);
        assert_eq!(result.len(), 2); // 8H, 5M
        assert_eq!(result[0], encode_op(5, 8)); // 8H
        assert_eq!(result[1], encode_op(0, 5)); // 5M
    }

    #[test]
    fn test_clip_cigar_ops_raw_with_eq_and_x_ops() {
        // 5=3X: clip 4 from start splits the = op, consuming 4 ref bases
        let cigar = &[encode_op(7, 5), encode_op(8, 3)]; // 5= 3X
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 4, true);
        assert_eq!(ref_consumed, 4);
        assert_eq!(result.len(), 3); // 4H, 1=, 3X
        assert_eq!(result[0], encode_op(5, 4));
        assert_eq!(result[1], encode_op(7, 1));
        assert_eq!(result[2], encode_op(8, 3));
    }

    #[test]
    fn test_clip_cigar_ops_raw_clip_entire_alignment_from_start() {
        // 10M: clip all 10 bases from start
        let cigar = &[encode_op(0, 10)]; // 10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 10, true);
        assert_eq!(ref_consumed, 10);
        assert_eq!(result.len(), 1); // 10H only
        assert_eq!(result[0], encode_op(5, 10));
    }

    #[test]
    fn test_clip_cigar_ops_raw_clip_entire_alignment_from_end() {
        // 10M: clip all 10 bases from end
        let cigar = &[encode_op(0, 10)]; // 10M
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 10, false);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 1); // 10H only
        assert_eq!(result[0], encode_op(5, 10));
    }

    #[test]
    fn test_clip_cigar_ops_raw_complex_cigar_start() {
        // 3S10M2I5M4S: clip_amount=8
        // existing_clip = 3 (3S), alignment_clip = 8 - 3 = 5
        // After stripping S: 10M2I5M4S
        // Clip 5 from 10M: split into 5 clipped + 5M remaining, ref=5
        let cigar = &[
            encode_op(4, 3),  // 3S
            encode_op(0, 10), // 10M
            encode_op(1, 2),  // 2I
            encode_op(0, 5),  // 5M
            encode_op(4, 4),  // 4S
        ];
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 8, true);
        assert_eq!(ref_consumed, 5);
        // Result: 8H (3S+5 from alignment), 5M, 2I, 5M, 4S
        assert_eq!(result[0], encode_op(5, 8)); // 8H
        assert_eq!(result[1], encode_op(0, 5)); // 5M (remaining from 10M)
        assert_eq!(result[2], encode_op(1, 2)); // 2I
        assert_eq!(result[3], encode_op(0, 5)); // 5M
        assert_eq!(result[4], encode_op(4, 4)); // 4S (trailing)
    }

    #[test]
    fn test_clip_cigar_ops_raw_complex_cigar_end() {
        // 3S10M2I5M4S: clip_amount=8
        // existing_clip at end = 4 (4S), alignment_clip = 8 - 4 = 4
        // After stripping trailing S: 3S10M2I5M
        // Clip 4 from end of 5M: split into 1M remaining + 4 clipped
        let cigar = &[
            encode_op(4, 3),  // 3S
            encode_op(0, 10), // 10M
            encode_op(1, 2),  // 2I
            encode_op(0, 5),  // 5M
            encode_op(4, 4),  // 4S
        ];
        let (result, ref_consumed) = clip_cigar_ops_raw(cigar, 8, false);
        assert_eq!(ref_consumed, 0);
        // Result: 3S, 10M, 2I, 1M, 8H (4S+4 from alignment)
        assert_eq!(result[0], encode_op(4, 3)); // 3S (leading)
        assert_eq!(result[1], encode_op(0, 10)); // 10M
        assert_eq!(result[2], encode_op(1, 2)); // 2I
        assert_eq!(result[3], encode_op(0, 1)); // 1M (remaining from 5M)
        assert_eq!(result[4], encode_op(5, 8)); // 8H (4S+4)
    }

    // ========================================================================
    // upgrade_clipping_raw tests
    // ========================================================================

    #[test]
    fn test_upgrade_clipping_raw_from_start_soft_only() {
        // 5S10M: clip_amount=3, should upgrade 3S to H
        // existing: H=0, S=5, clip_amount=3 => upgrade 3 => new_hard=3, remaining_soft=2
        let cigar = &[encode_op(4, 5), encode_op(0, 10)]; // 5S10M
        let (result, ref_consumed) = upgrade_clipping_raw(cigar, 3, true, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 3); // 3H, 2S, 10M
        assert_eq!(result[0], encode_op(5, 3)); // 3H
        assert_eq!(result[1], encode_op(4, 2)); // 2S
        assert_eq!(result[2], encode_op(0, 10)); // 10M
    }

    #[test]
    fn test_upgrade_clipping_raw_from_start_hard_and_soft() {
        // 2H5S10M: clip_amount=4, existing_hard=2, existing_soft=5
        // length_to_upgrade = min(5, 4-2) = 2
        // new_hard = 2+2 = 4, remaining_soft = 5-2 = 3
        let cigar = &[encode_op(5, 2), encode_op(4, 5), encode_op(0, 10)]; // 2H5S10M
        let (result, ref_consumed) = upgrade_clipping_raw(cigar, 4, true, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 3); // 4H, 3S, 10M
        assert_eq!(result[0], encode_op(5, 4)); // 4H
        assert_eq!(result[1], encode_op(4, 3)); // 3S
        assert_eq!(result[2], encode_op(0, 10)); // 10M
    }

    #[test]
    fn test_upgrade_clipping_raw_from_start_all_soft_to_hard() {
        // 5S10M: clip_amount=5, upgrade all 5S to H
        let cigar = &[encode_op(4, 5), encode_op(0, 10)]; // 5S10M
        let (result, ref_consumed) = upgrade_clipping_raw(cigar, 5, true, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 2); // 5H, 10M (no remaining S)
        assert_eq!(result[0], encode_op(5, 5)); // 5H
        assert_eq!(result[1], encode_op(0, 10)); // 10M
    }

    #[test]
    fn test_upgrade_clipping_raw_from_start_clip_amount_equals_hard() {
        // 3H5S10M: clip_amount=3, existing_hard=3
        // length_to_upgrade = min(5, 3-3) = 0
        // new_hard=3, remaining_soft=5
        let cigar = &[encode_op(5, 3), encode_op(4, 5), encode_op(0, 10)]; // 3H5S10M
        let (result, ref_consumed) = upgrade_clipping_raw(cigar, 3, true, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 3); // 3H, 5S, 10M
        assert_eq!(result[0], encode_op(5, 3));
        assert_eq!(result[1], encode_op(4, 5));
        assert_eq!(result[2], encode_op(0, 10));
    }

    #[test]
    fn test_upgrade_clipping_raw_from_end_soft_only() {
        // 10M5S: clip_amount=3, upgrade 3S at end to H
        let cigar = &[encode_op(0, 10), encode_op(4, 5)]; // 10M5S
        let (result, ref_consumed) = upgrade_clipping_raw(cigar, 3, false, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 3); // 10M, 2S, 3H
        assert_eq!(result[0], encode_op(0, 10)); // 10M
        assert_eq!(result[1], encode_op(4, 2)); // 2S
        assert_eq!(result[2], encode_op(5, 3)); // 3H
    }

    #[test]
    fn test_upgrade_clipping_raw_from_end_hard_and_soft() {
        // 10M5S2H: clip_amount=5, existing_hard=2, existing_soft=5
        // length_to_upgrade = min(5, 5-2) = 3
        // new_hard = 2+3 = 5, remaining_soft = 5-3 = 2
        let cigar = &[encode_op(0, 10), encode_op(4, 5), encode_op(5, 2)]; // 10M5S2H
        let (result, ref_consumed) = upgrade_clipping_raw(cigar, 5, false, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 3); // 10M, 2S, 5H
        assert_eq!(result[0], encode_op(0, 10));
        assert_eq!(result[1], encode_op(4, 2));
        assert_eq!(result[2], encode_op(5, 5));
    }

    #[test]
    fn test_upgrade_clipping_raw_from_end_all_soft_to_hard() {
        // 10M5S: clip_amount=5, upgrade all soft
        let cigar = &[encode_op(0, 10), encode_op(4, 5)]; // 10M5S
        let (result, ref_consumed) = upgrade_clipping_raw(cigar, 5, false, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 2); // 10M, 5H (no remaining S)
        assert_eq!(result[0], encode_op(0, 10));
        assert_eq!(result[1], encode_op(5, 5));
    }

    // ========================================================================
    // clip_cigar_start_raw / clip_cigar_end_raw tests
    // ========================================================================

    #[test]
    fn test_clip_cigar_start_raw_with_existing_hard_clip() {
        // 3H10M: clip 5 from start
        // existing_hard=3, existing_soft=0, post_clip_ops=[10M]
        // alignment_clip = 5 (called from clip_cigar_ops_raw: 5 - 3 = 2)
        // But directly calling clip_cigar_start_raw with clip_amount=2:
        // clips 2 from 10M => 8M remains, ref_consumed=2
        // total_hard = 3 + 0 + 2 = 5
        let cigar = &[encode_op(5, 3), encode_op(0, 10)]; // 3H10M
        let (result, ref_consumed) = clip_cigar_start_raw(cigar, 2, encode_op);
        assert_eq!(ref_consumed, 2);
        assert_eq!(result.len(), 2); // 5H, 8M
        assert_eq!(result[0], encode_op(5, 5)); // 3 existing + 2 new
        assert_eq!(result[1], encode_op(0, 8));
    }

    #[test]
    fn test_clip_cigar_end_raw_with_existing_hard_clip() {
        // 10M3H: clip 2 from end
        // existing_hard=3, existing_soft=0, post_clip_ops=[10M]
        // clips 2 from end of 10M => 8M remains
        // total_hard = 3 + 0 + 2 = 5
        let cigar = &[encode_op(0, 10), encode_op(5, 3)]; // 10M3H
        let (result, ref_consumed) = clip_cigar_end_raw(cigar, 2, encode_op);
        assert_eq!(ref_consumed, 0);
        assert_eq!(result.len(), 2); // 8M, 5H
        assert_eq!(result[0], encode_op(0, 8));
        assert_eq!(result[1], encode_op(5, 5)); // 3 existing + 2 new
    }

    #[test]
    fn test_clip_cigar_start_raw_insertion_at_exact_boundary() {
        // 3I10M: clip 1 from start
        // The insertion has 3 read bases. clip_amount=1, but I consumes entire (3 bases)
        // total read bases clipped = 3, ref_consumed = 0
        let cigar = &[encode_op(1, 3), encode_op(0, 10)]; // 3I10M
        let (result, ref_consumed) = clip_cigar_start_raw(cigar, 1, encode_op);
        assert_eq!(ref_consumed, 0);
        // 3I consumed entirely => 3H, 10M
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], encode_op(5, 3));
        assert_eq!(result[1], encode_op(0, 10));
    }

    #[test]
    fn test_clip_cigar_end_raw_insertion_at_exact_boundary() {
        // 10M3I: clip 1 from end
        // The insertion has 3 read bases. clip_amount=1, but I consumes entire (3 bases)
        let cigar = &[encode_op(0, 10), encode_op(1, 3)]; // 10M3I
        let (result, ref_consumed) = clip_cigar_end_raw(cigar, 1, encode_op);
        assert_eq!(ref_consumed, 0);
        // 3I consumed entirely => 10M, 3H
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], encode_op(0, 10));
        assert_eq!(result[1], encode_op(5, 3));
    }

    // ========================================================================
    // read_pos_at_ref_pos_raw tests
    // ========================================================================

    #[test]
    fn test_read_pos_at_ref_pos_raw_simple() {
        // 10M starting at position 100: ref pos 102 => query pos 3
        let cigar = &[encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 102, false), Some(3));
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_after_alignment() {
        // ref_pos after alignment => None
        let cigar = &[encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 110, false), None);
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_in_deletion_return_last() {
        // 5M3D5M: deletion at ref 105-107
        // ref_pos=106 is in deletion, return_last_base_if_deleted=true => returns 5
        let cigar = &[encode_op(0, 5), encode_op(2, 3), encode_op(0, 5)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 106, true), Some(5));
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_at_exact_start() {
        // 10M at position 100: query at ref 100 = 1
        let cigar = &[encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 100, false), Some(1));
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_at_exact_end() {
        // 10M at position 100: query at ref 109 = 10
        let cigar = &[encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 109, false), Some(10));
    }

    // ========================================================================
    // query_length_from_cigar tests
    // ========================================================================

    #[test]
    fn test_query_length_from_cigar_simple() {
        let cigar = &[encode_op(0, 10)]; // 10M
        assert_eq!(query_length_from_cigar(cigar), 10);
    }

    #[test]
    fn test_query_length_from_cigar_with_clips_and_insertion() {
        // 3S10M2I5M4S: query consuming = 3+10+2+5+4 = 24
        let cigar = &[
            encode_op(4, 3),  // 3S
            encode_op(0, 10), // 10M
            encode_op(1, 2),  // 2I
            encode_op(0, 5),  // 5M
            encode_op(4, 4),  // 4S
        ];
        assert_eq!(query_length_from_cigar(cigar), 24);
    }

    #[test]
    fn test_query_length_from_cigar_deletion_not_counted() {
        // 10M3D5M: deletion doesn't consume query
        let cigar = &[encode_op(0, 10), encode_op(2, 3), encode_op(0, 5)];
        assert_eq!(query_length_from_cigar(cigar), 15);
    }

    #[test]
    fn test_query_length_from_cigar_hard_clip_not_counted() {
        // 3H10M2H: hard clips don't consume query
        let cigar = &[encode_op(5, 3), encode_op(0, 10), encode_op(5, 2)];
        assert_eq!(query_length_from_cigar(cigar), 10);
    }

    #[test]
    fn test_query_length_from_cigar_eq_and_x() {
        // 5=3X: both consume query, total 8
        let cigar = &[encode_op(7, 5), encode_op(8, 3)];
        assert_eq!(query_length_from_cigar(cigar), 8);
    }

    // ========================================================================
    // unclipped_other_start / unclipped_other_end tests
    // ========================================================================

    #[test]
    fn test_unclipped_other_start_no_clips() {
        assert_eq!(unclipped_other_start(100, "10M"), 100);
    }

    #[test]
    fn test_unclipped_other_end_no_trailing_clips() {
        // 10M: end = 100 + 10 + 0 - 1 = 109
        assert_eq!(unclipped_other_end(100, "10M"), 109);
    }

    #[test]
    fn test_unclipped_other_start_complex() {
        // 3H5S10M: leading = 3+5=8, start = 100 - 8 = 92
        assert_eq!(unclipped_other_start(100, "3H5S10M"), 92);
    }

    #[test]
    fn test_unclipped_other_end_complex() {
        // 10M5S3H: ref=10, trailing=5+3=8, end = 100 + 10 + 8 - 1 = 117
        assert_eq!(unclipped_other_end(100, "10M5S3H"), 117);
    }

    /// The inclusive-end `-1` must be applied before the saturating additions.
    /// Applying it last saturates the sum to `i32::MAX` and then walks it back,
    /// reporting `i32::MAX - 1` for a coordinate whose exact value is `i32::MAX`.
    #[rstest::rstest]
    // 1M at MAX: MAX + 1 - 1 == MAX exactly, so no saturation is warranted.
    #[case::exactly_max(i32::MAX, "1M", i32::MAX)]
    // Genuinely past the top: saturation is correct here.
    #[case::past_max(i32::MAX, "10M", i32::MAX)]
    // A malformed MC whose spans saturate must still clamp at the top, not wrap.
    #[case::saturated_spans(100, "2000000000M2000000000M", i32::MAX)]
    fn test_unclipped_other_end_saturates_at_the_coordinate_ceiling(
        #[case] mate_pos: i32,
        #[case] mc_cigar: &str,
        #[case] expected: i32,
    ) {
        assert_eq!(unclipped_other_end(mate_pos, mc_cigar), expected);
    }

    /// An unclipped start below the first coordinate is real data, not a bug:
    /// a read whose leading clips run off the front of the contig has one.
    ///
    /// The record's own side (`unclipped_start_from_raw_cigar`) does not clamp,
    /// and the template-coordinate key compares the two against each other, so
    /// clamping only the mate side would make the two records of one template
    /// canonicalize differently -- and would diverge from fgbio and
    /// `samtools sort --template-coordinate`. This pins the choice against a
    /// well-meaning `.max(0)`.
    #[rstest::rstest]
    // 1-based pos 3 with 10 leading clips: 3 - 10 = -7, before the contig start.
    #[case::clips_past_contig_start(3, "10S5M", -7)]
    #[case::clips_exactly_to_zero(10, "10S5M", 0)]
    // Even a malformed, saturated clip stays a (very negative) number rather
    // than wrapping: `saturating_sub` bottoms out at `i32::MIN`.
    #[case::saturated_clip_does_not_wrap(100, "2000000000S2000000000S", 100 - i32::MAX)]
    fn test_unclipped_other_start_is_not_clamped_to_the_contig_start(
        #[case] mate_pos: i32,
        #[case] mc_cigar: &str,
        #[case] expected: i32,
    ) {
        assert_eq!(unclipped_other_start(mate_pos, mc_cigar), expected);
    }

    // ========================================================================
    // unclipped_5prime_raw tests
    // ========================================================================

    use rstest::rstest;

    #[rstest]
    // Forward strand: 5' = start. 5S10M at pos 100 → unclipped start = 100 - 5 = 95
    #[case(false, &[encode_op(4, 5), encode_op(0, 10)], 100, 95)]
    // Reverse strand: 5' = end. 10M5S at pos 100 → ref end = 100+10-1=109, unclipped end = 109+5=114
    #[case(true, &[encode_op(0, 10), encode_op(4, 5)], 100, 114)]
    // Forward strand: no clipping. 10M at pos 50 → 50
    #[case(false, &[encode_op(0, 10)], 50, 50)]
    // Reverse strand: no clipping. 10M at pos 50 → 50+10-1=59
    #[case(true, &[encode_op(0, 10)], 50, 59)]
    fn test_unclipped_5prime_raw_basic(
        #[case] is_reverse: bool,
        #[case] cigar: &[u32],
        #[case] pos: i32,
        #[case] expected: i32,
    ) {
        let flags = if is_reverse { flags::REVERSE } else { 0 };
        let rec = make_bam_bytes(0, pos, flags, b"read", cigar, 15, -1, -1, &[]);
        assert_eq!(unclipped_5prime_raw(&rec, pos, is_reverse), expected);
    }

    #[test]
    fn test_unclipped_5prime_raw_no_cigar() {
        // n_cigar_op == 0 → returns pos unchanged
        let rec = make_bam_bytes(0, 100, 0, b"read", &[], 15, -1, -1, &[]);
        assert_eq!(unclipped_5prime_raw(&rec, 100, false), 100);
    }

    #[test]
    fn test_unclipped_5prime_raw_truncated_record() {
        // Record too short for CIGAR → falls back to returning pos
        let flags = 0u16;
        let mut rec = make_bam_bytes(0, 100, flags, b"read", &[encode_op(0, 10)], 10, -1, -1, &[]);
        // Truncate so cigar_end > bam.len()
        rec.truncate(36);
        assert_eq!(unclipped_5prime_raw(&rec, 100, false), 100);
    }

    // ========================================================================
    // cigar_to_string_from_raw tests
    // ========================================================================

    #[test]
    fn test_cigar_to_string_simple() {
        let rec = make_bam_bytes(0, 0, 0, b"r1", &[encode_op(0, 100)], 100, -1, -1, &[]);
        assert_eq!(cigar_to_string_from_raw(&rec), "100M");
    }

    #[test]
    fn test_cigar_to_string_complex() {
        let ops = [encode_op(4, 5), encode_op(0, 50), encode_op(1, 3), encode_op(0, 42)];
        let rec = make_bam_bytes(0, 0, 0, b"r1", &ops, 100, -1, -1, &[]);
        assert_eq!(cigar_to_string_from_raw(&rec), "5S50M3I42M");
    }

    #[test]
    fn test_cigar_to_string_no_cigar() {
        let rec = make_bam_bytes(0, 0, 0, b"r1", &[], 4, -1, -1, &[]);
        assert_eq!(cigar_to_string_from_raw(&rec), "");
    }

    #[cfg(feature = "noodles")]
    #[test]
    fn test_cigar_op_kind_all_ops() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        assert_eq!(cigar_op_kind(0), Kind::Match);
        assert_eq!(cigar_op_kind(1), Kind::Insertion);
        assert_eq!(cigar_op_kind(2), Kind::Deletion);
        assert_eq!(cigar_op_kind(3), Kind::Skip);
        assert_eq!(cigar_op_kind(4), Kind::SoftClip);
        assert_eq!(cigar_op_kind(5), Kind::HardClip);
        assert_eq!(cigar_op_kind(7), Kind::SequenceMatch);
        assert_eq!(cigar_op_kind(8), Kind::SequenceMismatch);
        // Unknown ops (6 = Pad, 9+ = Pad fallback)
        assert_eq!(cigar_op_kind(6), Kind::Pad);
        assert_eq!(cigar_op_kind(9), Kind::Pad);
        assert_eq!(cigar_op_kind(15), Kind::Pad);
    }

    #[cfg(feature = "noodles")]
    #[test]
    fn test_cigar_op_kind_with_length_bits() {
        use noodles::sam::alignment::record::cigar::op::Kind;
        // Raw BAM encodes op_type in low 4 bits, length in upper 28 bits
        // Verify cigar_op_kind masks correctly: 50M = (50 << 4) | 0
        assert_eq!(cigar_op_kind(50 << 4), Kind::Match);
        // 10I = (10 << 4) | 1
        assert_eq!(cigar_op_kind(10 << 4 | 1), Kind::Insertion);
        // 5D = (5 << 4) | 2
        assert_eq!(cigar_op_kind(5 << 4 | 2), Kind::Deletion);
    }

    #[test]
    fn test_view_cigar_methods() {
        use crate::fields::RawRecordView;
        use crate::testutil::*;
        let rec = make_bam_bytes(
            0,
            100,
            0,
            b"r",
            &[encode_op(0, 50), encode_op(2, 5), encode_op(0, 50)],
            100,
            -1,
            -1,
            &[],
        );
        let v = RawRecordView::new(&rec);
        assert_eq!(v.cigar_ops_vec(), vec![encode_op(0, 50), encode_op(2, 5), encode_op(0, 50)]);
        assert_eq!(v.cigar_to_string(), "50M5D50M");
        assert_eq!(v.reference_length(), 105);
        assert_eq!(v.query_length(), 100);
        assert_eq!(v.cigar_raw_bytes().len(), 12);
        assert_eq!(v.cigar_ops_iter().count(), 3);
    }

    // ========================================================================
    // cigar_lengths tests
    // ========================================================================

    #[test]
    fn test_cigar_lengths_matches_separate_calls() {
        use crate::fields::RawRecordView;

        // Record with mixed CIGAR: 50M 5D 10I 35M = ref 90, query 95
        let rec = make_bam_bytes(
            0,
            100,
            0,
            b"r",
            &[
                encode_op(0, 50), // 50M
                encode_op(2, 5),  // 5D (ref only)
                encode_op(1, 10), // 10I (query only)
                encode_op(0, 35), // 35M
            ],
            95,
            -1,
            -1,
            &[],
        );
        let v = RawRecordView::new(&rec);
        let (ref_len, q_len) = v.cigar_lengths();
        assert_eq!(ref_len, 90);
        assert_eq!(q_len, 95);
        // Agree with the single-purpose methods
        assert_eq!(ref_len, v.reference_length());
        assert_eq!(q_len, v.query_length());
    }

    #[test]
    fn test_cigar_lengths_empty_cigar() {
        use crate::fields::RawRecordView;
        let rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        let v = RawRecordView::new(&rec);
        assert_eq!(v.cigar_lengths(), (0, 0));
    }

    // ========================================================================
    // CigarKind tests
    // ========================================================================

    #[test]
    fn test_cigar_kind_as_char_all_variants() {
        assert_eq!(CigarKind::Match.as_char(), 'M');
        assert_eq!(CigarKind::Insertion.as_char(), 'I');
        assert_eq!(CigarKind::Deletion.as_char(), 'D');
        assert_eq!(CigarKind::Skip.as_char(), 'N');
        assert_eq!(CigarKind::SoftClip.as_char(), 'S');
        assert_eq!(CigarKind::HardClip.as_char(), 'H');
        assert_eq!(CigarKind::Pad.as_char(), 'P');
        assert_eq!(CigarKind::SequenceMatch.as_char(), '=');
        assert_eq!(CigarKind::SequenceMismatch.as_char(), 'X');
    }

    #[test]
    fn test_cigar_kind_from_op_type_all_valid() {
        assert_eq!(CigarKind::from_op_type(0), Some(CigarKind::Match));
        assert_eq!(CigarKind::from_op_type(1), Some(CigarKind::Insertion));
        assert_eq!(CigarKind::from_op_type(2), Some(CigarKind::Deletion));
        assert_eq!(CigarKind::from_op_type(3), Some(CigarKind::Skip));
        assert_eq!(CigarKind::from_op_type(4), Some(CigarKind::SoftClip));
        assert_eq!(CigarKind::from_op_type(5), Some(CigarKind::HardClip));
        assert_eq!(CigarKind::from_op_type(6), Some(CigarKind::Pad));
        assert_eq!(CigarKind::from_op_type(7), Some(CigarKind::SequenceMatch));
        assert_eq!(CigarKind::from_op_type(8), Some(CigarKind::SequenceMismatch));
    }

    #[test]
    fn test_cigar_kind_from_op_type_reserved_returns_none() {
        for v in 9u32..=15 {
            assert_eq!(CigarKind::from_op_type(v), None, "op_type {v} should be None");
        }
    }

    // ========================================================================
    // CigarOp tests
    // ========================================================================

    #[test]
    fn test_cigar_op_round_trip_from_raw() {
        // All 9 op types with a representative length each.
        let cases: &[(u32, u32)] = &[
            (0, 50), // 50M
            (1, 5),  // 5I
            (2, 3),  // 3D
            (3, 10), // 10N
            (4, 4),  // 4S
            (5, 2),  // 2H
            (6, 1),  // 1P
            (7, 20), // 20=
            (8, 7),  // 7X
        ];
        for &(op_type, op_len) in cases {
            let raw = (op_len << 4) | op_type;
            let cigar_op = CigarOp::from_raw(raw);
            assert_eq!(cigar_op.raw(), raw, "raw round-trip for op_type {op_type}");
            assert_eq!(cigar_op.len(), op_len, "len for op_type {op_type}");
            let expected_kind = CigarKind::from_op_type(op_type).unwrap();
            assert_eq!(cigar_op.kind(), expected_kind, "kind for op_type {op_type}");
        }
    }

    #[test]
    #[should_panic(expected = "CIGAR op length must be non-zero")]
    fn test_cigar_op_new_rejects_zero_length() {
        let _ = CigarOp::new(CigarKind::Match, 0);
    }

    #[test]
    fn test_cigar_op_from_raw_allows_zero_length_for_malformed_data() {
        // from_raw intentionally preserves malformed wire data, so zero-length
        // must still be constructable via the raw path even though CigarOp::new
        // rejects it.
        let op = CigarOp::from_raw(CigarKind::Match as u32); // len = 0
        assert_eq!(op.len(), 0);
    }

    #[test]
    fn test_cigar_op_new_round_trip() {
        // Construct via CigarOp::new and verify kind/len match.
        let cases: &[(CigarKind, u32)] = &[
            (CigarKind::Match, 100),
            (CigarKind::Insertion, 3),
            (CigarKind::Deletion, 8),
            (CigarKind::Skip, 1000),
            (CigarKind::SoftClip, 5),
            (CigarKind::HardClip, 2),
            (CigarKind::Pad, 1),
            (CigarKind::SequenceMatch, 42),
            (CigarKind::SequenceMismatch, 7),
        ];
        for &(kind, len) in cases {
            let op = CigarOp::new(kind, len);
            assert_eq!(op.kind(), kind, "kind mismatch for {kind:?}");
            assert_eq!(op.len(), len, "len mismatch for {kind:?}");
        }
    }

    #[test]
    fn test_cigar_op_consumes_query_matches_free_fn() {
        // For all 9 op types, CigarOp::consumes_query should match consumes_query(op_type).
        for op_type in 0u32..9 {
            let raw = (10 << 4) | op_type; // length=10, arbitrary
            let cigar_op = CigarOp::from_raw(raw);
            assert_eq!(
                cigar_op.consumes_query(),
                consumes_query(op_type),
                "consumes_query mismatch for op_type {op_type}"
            );
        }
    }

    #[test]
    fn test_cigar_op_consumes_ref_matches_free_fn() {
        // For all 9 op types, CigarOp::consumes_ref should match consumes_ref(op_type).
        for op_type in 0u32..9 {
            let raw = (10 << 4) | op_type;
            let cigar_op = CigarOp::from_raw(raw);
            assert_eq!(
                cigar_op.consumes_ref(),
                consumes_ref(op_type),
                "consumes_ref mismatch for op_type {op_type}"
            );
        }
    }

    // ========================================================================
    // cigar_ops_typed on RawRecordView
    // ========================================================================

    #[test]
    fn test_cigar_ops_typed_known_cigar() {
        // 5S10M3I2D4= 2X 1H => 8 ops
        let raw_cigar = &[
            encode_op(4, 5),  // 5S
            encode_op(0, 10), // 10M
            encode_op(1, 3),  // 3I
            encode_op(2, 2),  // 2D
            encode_op(7, 4),  // 4=
            encode_op(8, 2),  // 2X
            encode_op(5, 1),  // 1H
        ];
        let seq_len = 5 + 10 + 3 + 4 + 2; // query-consuming ops
        let rec = make_bam_bytes(0, 100, 0, b"rd", raw_cigar, seq_len, -1, -1, &[]);
        let view = RawRecordView::new(&rec);

        let typed: Vec<CigarOp> = view.cigar_ops_typed().collect();
        assert_eq!(typed.len(), 7);
        assert_eq!(typed[0].kind(), CigarKind::SoftClip);
        assert_eq!(typed[0].len(), 5);
        assert_eq!(typed[1].kind(), CigarKind::Match);
        assert_eq!(typed[1].len(), 10);
        assert_eq!(typed[2].kind(), CigarKind::Insertion);
        assert_eq!(typed[2].len(), 3);
        assert_eq!(typed[3].kind(), CigarKind::Deletion);
        assert_eq!(typed[3].len(), 2);
        assert_eq!(typed[4].kind(), CigarKind::SequenceMatch);
        assert_eq!(typed[4].len(), 4);
        assert_eq!(typed[5].kind(), CigarKind::SequenceMismatch);
        assert_eq!(typed[5].len(), 2);
        assert_eq!(typed[6].kind(), CigarKind::HardClip);
        assert_eq!(typed[6].len(), 1);
    }

    #[test]
    fn test_cigar_ops_typed_empty_cigar() {
        let rec = make_bam_bytes(0, 0, 0, b"r", &[], 0, -1, -1, &[]);
        let view = RawRecordView::new(&rec);
        let typed: Vec<CigarOp> = view.cigar_ops_typed().collect();
        assert!(typed.is_empty());
    }

    #[test]
    fn test_cigar_ops_typed_matches_cigar_ops_iter() {
        // cigar_ops_typed().map(|op| op.raw()) should equal cigar_ops_iter()
        let raw_cigar = &[
            encode_op(0, 30), // 30M
            encode_op(1, 2),  // 2I
            encode_op(0, 20), // 20M
            encode_op(4, 5),  // 5S
        ];
        let rec = make_bam_bytes(0, 50, 0, b"rd", raw_cigar, 57, -1, -1, &[]);
        let view = RawRecordView::new(&rec);

        let from_raw: Vec<u32> = view.cigar_ops_iter().collect();
        let from_typed: Vec<u32> = view.cigar_ops_typed().map(CigarOp::raw).collect();
        assert_eq!(from_raw, from_typed);
    }

    #[cfg(feature = "noodles")]
    mod noodles_roundtrip {
        use super::*;

        #[test]
        fn test_cigar_kind_noodles_roundtrip_all_variants() {
            use noodles::sam::alignment::record::cigar::op::Kind;

            let kinds = [
                CigarKind::Match,
                CigarKind::Insertion,
                CigarKind::Deletion,
                CigarKind::Skip,
                CigarKind::SoftClip,
                CigarKind::HardClip,
                CigarKind::Pad,
                CigarKind::SequenceMatch,
                CigarKind::SequenceMismatch,
            ];

            for original in kinds {
                // CigarKind -> noodles Kind -> CigarKind
                let noodles_kind: Kind = original.into();
                let round_tripped: CigarKind = noodles_kind.into();
                assert_eq!(round_tripped, original, "round-trip failed for {original:?}");
            }
        }

        #[test]
        fn test_noodles_kind_to_cigar_kind_all_variants() {
            use noodles::sam::alignment::record::cigar::op::Kind;

            let noodles_kinds = [
                Kind::Match,
                Kind::Insertion,
                Kind::Deletion,
                Kind::Skip,
                Kind::SoftClip,
                Kind::HardClip,
                Kind::Pad,
                Kind::SequenceMatch,
                Kind::SequenceMismatch,
            ];
            let expected = [
                CigarKind::Match,
                CigarKind::Insertion,
                CigarKind::Deletion,
                CigarKind::Skip,
                CigarKind::SoftClip,
                CigarKind::HardClip,
                CigarKind::Pad,
                CigarKind::SequenceMatch,
                CigarKind::SequenceMismatch,
            ];

            for (nk, ck) in noodles_kinds.iter().zip(expected.iter()) {
                let converted: CigarKind = (*nk).into();
                assert_eq!(converted, *ck, "noodles->CigarKind failed for {nk:?}");
            }
        }
    }
}
