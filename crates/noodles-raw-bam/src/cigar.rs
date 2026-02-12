use crate::fields::{
    flags,
    flags::{REVERSE, UNMAPPED},
    l_read_name, n_cigar_op, pos,
};

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

        // M (0), D (2), N (3), = (7), X (8) consume reference bases
        if matches!(op_type, 0 | 2 | 3 | 7 | 8) {
            ref_len += op_len;
        }
    }

    ref_len
}

/// Compute reference length directly from raw CIGAR bytes (zero allocation).
#[inline]
#[must_use]
pub fn reference_length_from_raw_bam(bam: &[u8]) -> i32 {
    let n_cigar_op = n_cigar_op(bam) as usize;
    if n_cigar_op == 0 {
        return 0;
    }
    let l_read_name = l_read_name(bam) as usize;
    let cigar_start = 32 + l_read_name;
    let cigar_end = cigar_start + n_cigar_op * 4;
    if cigar_end > bam.len() {
        return 0;
    }

    let mut ref_len = 0i32;
    for i in 0..n_cigar_op {
        let op = cigar_op_at(bam, cigar_start + i * 4);
        let op_type = op & 0xF;
        if matches!(op_type, 0 | 2 | 3 | 7 | 8) {
            ref_len += (op >> 4).cast_signed();
        }
    }
    ref_len
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
        if matches!(op_type, 0 | 1 | 4 | 7 | 8) {
            len += op_len;
        }
    }
    len
}

/// Calculate unclipped start position from BAM record.
///
/// For forward strand reads, this is the 5' position.
/// Returns: `pos - leading_clips` (0-based like the input position).
#[inline]
#[must_use]
pub fn unclipped_start_from_cigar(pos: i32, cigar_ops: &[u32]) -> i32 {
    let mut clipped = 0i32;

    for &op in cigar_ops {
        let op_len = (op >> 4).cast_signed();
        let op_type = op & 0xF;

        match op_type {
            4 | 5 => clipped += op_len, // S (4) or H (5)
            _ => break,                 // Non-clip operation
        }
    }

    pos - clipped
}

/// Calculate unclipped end position from BAM record.
///
/// For reverse strand reads, this is the 5' position.
/// Returns: `pos + ref_length + trailing_clips - 1` (0-based).
#[inline]
#[must_use]
pub fn unclipped_end_from_cigar(pos: i32, cigar_ops: &[u32]) -> i32 {
    let mut ref_len = 0i32;
    let mut trailing_clips = 0i32;
    let mut saw_ref_op = false;

    for &op in cigar_ops {
        let op_len = (op >> 4).cast_signed();
        let op_type = op & 0xF;

        match op_type {
            0 | 2 | 3 | 7 | 8 => {
                // M (0), D (2), N (3), = (7), X (8)
                ref_len += op_len;
                trailing_clips = 0;
                saw_ref_op = true;
            }
            4 | 5 if saw_ref_op => {
                // S (4) or H (5) after ref-consuming ops
                trailing_clips += op_len;
            }
            _ => {}
        }
    }

    pos + ref_len + trailing_clips - 1
}

/// Calculate unclipped 5' coordinate for this read.
///
/// For forward strand: `unclipped_start` (leftmost including clips)
/// For reverse strand: `unclipped_end` (rightmost including clips)
#[inline]
#[must_use]
pub fn unclipped_5prime(pos: i32, reverse: bool, cigar_ops: &[u32]) -> i32 {
    if reverse {
        unclipped_end_from_cigar(pos, cigar_ops)
    } else {
        unclipped_start_from_cigar(pos, cigar_ops)
    }
}

/// Calculate unclipped 5' coordinate using 1-based positions (matching noodles).
///
/// Adds 1 to the 0-based BAM position before applying clip adjustment,
/// producing values identical to `record_utils::unclipped_five_prime_position()`.
///
/// Returns 0 for unmapped reads (matching `get_unclipped_position_for_groupkey`).
/// Returns `i32::MAX` for mapped reads with no CIGAR operations.
#[inline]
#[must_use]
pub fn unclipped_5prime_1based(
    pos_0based: i32,
    reverse: bool,
    unmapped: bool,
    cigar_ops: &[u32],
) -> i32 {
    if unmapped {
        return 0;
    }
    if cigar_ops.is_empty() {
        return i32::MAX;
    }
    // Convert to 1-based, then apply clip adjustment
    let pos_1based = pos_0based + 1;
    if reverse {
        unclipped_end_from_cigar(pos_1based, cigar_ops)
    } else {
        unclipped_start_from_cigar(pos_1based, cigar_ops)
    }
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
    let mate_pos_1based = mate_pos_0based + 1;
    if mate_reverse {
        unclipped_other_end(mate_pos_1based, mc_cigar)
    } else {
        unclipped_other_start(mate_pos_1based, mc_cigar)
    }
}

/// Calculate mate's unclipped start from MC tag CIGAR string.
///
/// For forward strand mates, this is the 5' position.
#[inline]
#[must_use]
pub fn unclipped_other_start(mate_pos: i32, mc_cigar: &str) -> i32 {
    mate_pos - parse_leading_clips(mc_cigar)
}

/// Calculate mate's unclipped end from MC tag CIGAR string.
///
/// For reverse strand mates, this is the 5' position.
#[inline]
#[must_use]
pub fn unclipped_other_end(mate_pos: i32, mc_cigar: &str) -> i32 {
    let (ref_len, trailing_clips) = parse_ref_len_and_trailing_clips(mc_cigar);
    mate_pos + ref_len + trailing_clips - 1
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

        let consumes_ref = matches!(op_type, 0 | 2 | 3 | 7 | 8); // M, D, N, =, X
        let consumes_query = matches!(op_type, 0 | 1 | 4 | 7 | 8); // M, I, S, =, X

        let op_ref_start = alignment_start + ref_offset;

        if consumes_ref {
            let op_ref_end = op_ref_start + op_len - 1;

            if ref_pos >= op_ref_start && ref_pos <= op_ref_end {
                if consumes_query {
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

        if consumes_ref {
            ref_offset += op_len;
        }
        if consumes_query {
            query_offset += op_len;
        }
    }

    None
}

/// Zero-allocation unclipped 5' position from raw BAM bytes (0-based, for sorting).
///
/// Unlike [`unclipped_5prime_from_raw_bam`] which adds 1 for 1-based group keys,
/// this uses the raw 0-based BAM position directly — matching what the sort
/// comparator needs.
#[inline]
pub(crate) fn unclipped_5prime_sort(bam: &[u8], pos_0based: i32, reverse: bool) -> i32 {
    let n_cigar_op = n_cigar_op(bam) as usize;
    if n_cigar_op == 0 {
        return pos_0based;
    }
    let l_read_name = l_read_name(bam) as usize;
    let cigar_start = 32 + l_read_name;
    if cigar_start + n_cigar_op * 4 > bam.len() {
        return pos_0based;
    }
    if reverse {
        unclipped_end_from_raw_cigar(pos_0based, bam, cigar_start, n_cigar_op)
    } else {
        unclipped_start_from_raw_cigar(pos_0based, bam, cigar_start, n_cigar_op)
    }
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

        match op_type {
            0 | 2 | 3 | 7 | 8 => {
                ref_len += op_len;
                trailing_clips = 0;
                saw_ref_op = true;
            }
            4 | 5 if saw_ref_op => {
                trailing_clips += op_len;
            }
            _ => {}
        }
    }

    pos + ref_len + trailing_clips - 1
}

/// Parse leading clips (S or H) from a CIGAR string.
///
/// Returns the total length of leading soft/hard clips.
#[inline]
fn parse_leading_clips(cigar: &str) -> i32 {
    let mut clipped = 0i32;
    let mut num_start = 0;
    let bytes = cigar.as_bytes();

    for (i, &c) in bytes.iter().enumerate() {
        if c.is_ascii_digit() {
            continue;
        }
        // Parse the number before this operation
        let num: i32 = cigar[num_start..i].parse().unwrap_or(0);

        if c == b'S' || c == b'H' {
            clipped += num;
            num_start = i + 1;
        } else {
            // Non-clip operation, stop processing
            break;
        }
    }

    clipped
}

/// Parse reference length and trailing clips from a CIGAR string.
///
/// Returns `(reference_length, trailing_clips)`.
/// Reference length is the sum of M/D/N/=/X operations.
/// Trailing clips are S/H operations after the last reference-consuming op.
#[inline]
fn parse_ref_len_and_trailing_clips(cigar: &str) -> (i32, i32) {
    let mut ref_len = 0i32;
    let mut trailing_clips = 0i32;
    let mut num_start = 0;
    let mut saw_ref_op = false;
    let bytes = cigar.as_bytes();

    for (i, &c) in bytes.iter().enumerate() {
        if c.is_ascii_digit() {
            continue;
        }

        let num: i32 = cigar[num_start..i].parse().unwrap_or(0);
        num_start = i + 1;

        match c {
            b'M' | b'D' | b'N' | b'=' | b'X' => {
                ref_len += num;
                trailing_clips = 0; // Reset trailing clips
                saw_ref_op = true;
            }
            b'S' | b'H' if saw_ref_op => {
                trailing_clips += num;
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
        let consumes_read = matches!(op_type, 0 | 1 | 7 | 8); // M, I, =, X
        let consumes_ref = matches!(op_type, 0 | 2 | 3 | 7 | 8); // M, D, N, =, X

        if consumes_read && op_len > (clip_amount - read_bases_clipped) {
            if op_type == 1 {
                // Insertion: consume entire at clip boundary
                read_bases_clipped += op_len;
            } else {
                // Split the operation
                let remaining_clip = clip_amount - read_bases_clipped;
                let remaining_length = op_len - remaining_clip;
                read_bases_clipped += remaining_clip;
                if consumes_ref {
                    ref_bases_clipped += remaining_clip;
                }
                new_ops.push(encode_op(op_type, remaining_length));
            }
        } else {
            if consumes_read {
                read_bases_clipped += op_len;
            }
            if consumes_ref {
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
        let consumes_read = matches!(op_type, 0 | 1 | 7 | 8); // M, I, =, X

        if consumes_read && op_len > (clip_amount - read_bases_clipped) {
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
        } else if consumes_read {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::*;

    // ========================================================================
    // unclipped_start_from_cigar / unclipped_end_from_cigar tests
    // ========================================================================

    #[test]
    fn test_unclipped_start_from_cigar_no_clips() {
        // 10M = (10 << 4)
        let cigar = &[(10 << 4)];
        assert_eq!(unclipped_start_from_cigar(100, cigar), 100);
    }

    #[test]
    fn test_unclipped_start_from_cigar_soft_clip() {
        // 5S10M: soft clip (op type 4), then match
        let cigar = &[(5 << 4) | 4, (10 << 4)];
        assert_eq!(unclipped_start_from_cigar(100, cigar), 95);
    }

    #[test]
    fn test_unclipped_start_from_cigar_hard_clip() {
        // 3H10M: hard clip (op type 5), then match
        let cigar = &[(3 << 4) | 5, (10 << 4)];
        assert_eq!(unclipped_start_from_cigar(100, cigar), 97);
    }

    #[test]
    fn test_unclipped_end_from_cigar_no_clips() {
        // 10M: end = pos + 10 - 1 = 109
        let cigar = &[(10 << 4)];
        assert_eq!(unclipped_end_from_cigar(100, cigar), 109);
    }

    #[test]
    fn test_unclipped_end_from_cigar_trailing_clips() {
        // 10M5S3H: end = 100 + 10 + 5 + 3 - 1 = 117
        let cigar = &[(10 << 4), (5 << 4) | 4, (3 << 4) | 5];
        assert_eq!(unclipped_end_from_cigar(100, cigar), 117);
    }

    // ========================================================================
    // unclipped_5prime tests
    // ========================================================================

    #[test]
    fn test_unclipped_5prime_forward_vs_reverse() {
        // 5S10M3S: forward uses start, reverse uses end
        let cigar = &[(5 << 4) | 4, (10 << 4), (3 << 4) | 4];
        // Forward: unclipped_start = 100 - 5 = 95
        assert_eq!(unclipped_5prime(100, false, cigar), 95);
        // Reverse: unclipped_end = 100 + 10 + 3 - 1 = 112
        assert_eq!(unclipped_5prime(100, true, cigar), 112);
    }

    // ========================================================================
    // unclipped_5prime_1based tests
    // ========================================================================

    #[test]
    fn test_unclipped_5prime_1based_unmapped() {
        assert_eq!(unclipped_5prime_1based(100, false, true, &[(10 << 4)]), 0);
    }

    #[test]
    fn test_unclipped_5prime_1based_no_cigar() {
        assert_eq!(unclipped_5prime_1based(100, false, false, &[]), i32::MAX);
    }

    #[test]
    fn test_unclipped_5prime_1based_forward() {
        // 5S10M: forward 5' = pos+1 - 5 = 96
        let cigar = &[(5 << 4) | 4, (10 << 4)];
        assert_eq!(unclipped_5prime_1based(100, false, false, cigar), 96);
    }

    #[test]
    fn test_unclipped_5prime_1based_reverse() {
        // 10M5S: reverse 5' = pos+1 + 10 + 5 - 1 = 115
        let cigar = &[(10 << 4), (5 << 4) | 4];
        assert_eq!(unclipped_5prime_1based(100, true, false, cigar), 115);
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
        assert_eq!(parse_leading_clips("5S10M"), 5);
    }

    #[test]
    fn test_parse_leading_clips_hard_and_soft() {
        // "3H5S10M" -> leading clips = 3 + 5 = 8
        assert_eq!(parse_leading_clips("3H5S10M"), 8);
    }

    #[test]
    fn test_parse_leading_clips_no_clips() {
        // "10M" -> leading clips = 0
        assert_eq!(parse_leading_clips("10M"), 0);
    }

    #[test]
    fn test_parse_leading_clips_all_clips() {
        // "5S3H" — all clips, no ref ops
        assert_eq!(parse_leading_clips("5S3H"), 8);
    }

    // ========================================================================
    // parse_ref_len_and_trailing_clips tests
    // ========================================================================

    #[test]
    fn test_parse_ref_len_and_trailing_clips_basic() {
        // "10M5S" -> ref_len=10, trailing_clips=5
        assert_eq!(parse_ref_len_and_trailing_clips("10M5S"), (10, 5));
    }

    #[test]
    fn test_parse_ref_len_and_trailing_clips_complex() {
        // "5S10M2I3D5M3S2H"
        // ref consuming: 10M + 3D + 5M = 18
        // trailing clips: 3S + 2H = 5
        assert_eq!(parse_ref_len_and_trailing_clips("5S10M2I3D5M3S2H"), (18, 5));
    }

    #[test]
    fn test_parse_ref_len_and_trailing_clips_no_ref_ops() {
        // "5S3H" — no ref-consuming ops
        assert_eq!(parse_ref_len_and_trailing_clips("5S3H"), (0, 0));
    }

    #[test]
    fn test_parse_ref_len_and_trailing_clips_eq_and_x() {
        // "10=3X" → ref_len=13, trailing_clips=0
        assert_eq!(parse_ref_len_and_trailing_clips("10=3X"), (13, 0));
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
    fn test_read_pos_at_ref_pos_raw_before_alignment() {
        // ref_pos < alignment_start => None
        let cigar = &[encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 99, false), None);
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_after_alignment() {
        // ref_pos after alignment => None
        let cigar = &[encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 110, false), None);
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_in_deletion() {
        // 5M3D5M: deletion at ref 105-107
        // ref_pos=106 is in deletion => None (without return_last_base_if_deleted)
        let cigar = &[encode_op(0, 5), encode_op(2, 3), encode_op(0, 5)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 106, false), None);
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_in_deletion_return_last() {
        // 5M3D5M: deletion at ref 105-107
        // ref_pos=106 is in deletion, return_last_base_if_deleted=true => returns 5
        let cigar = &[encode_op(0, 5), encode_op(2, 3), encode_op(0, 5)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 106, true), Some(5));
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_with_insertion() {
        // 5M3I5M: insertion doesn't consume reference
        // ref 105 is first base of second 5M, query pos = 5 (from first M) + 3 (from I) + 1 = 9
        let cigar = &[encode_op(0, 5), encode_op(1, 3), encode_op(0, 5)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 105, false), Some(9));
    }

    #[test]
    fn test_read_pos_at_ref_pos_raw_with_soft_clip() {
        // 3S10M: soft clip consumes query but not reference
        // ref 100 = first ref base in 10M, query pos = 3 (from S) + 1 = 4
        let cigar = &[encode_op(4, 3), encode_op(0, 10)];
        assert_eq!(read_pos_at_ref_pos_raw(cigar, 100, 100, false), Some(4));
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

    #[test]
    fn test_query_length_from_cigar_empty() {
        assert_eq!(query_length_from_cigar(&[]), 0);
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
}
