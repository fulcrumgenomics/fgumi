//! Raw byte-level BAM record comparison helpers.
//!
//! Provides functions to compare BAM records at the byte level, avoiding the overhead
//! of full record decoding. Supports full-record comparison, core-field-only comparison,
//! tag comparison (both order-dependent and order-independent), and a structured
//! comparison that reports which parts differ.

use fgumi_raw_bam::fields::{aux_data_offset_from_record, tag_value_size};

use crate::sam::SamTag;

/// Returns `true` for aux tags that fgumi emits but fgbio never persists, and which are
/// therefore ignored (both presence and value, on either side) in every fgumi-vs-fgbio
/// content comparison — an accepted divergence, see the compare-hardening design spec
/// §"Accepted divergences".
///
/// Currently just `tc`, the template-coordinate sort key that `fgumi zipper` writes on
/// secondary/supplementary reads (and that `fgumi dedup` later consumes); fgbio computes
/// the same key transiently but never writes it to the BAM. The carve-out is deliberately
/// narrow — only `tc` — so a genuine difference in any other tag still reports a DIFFER.
#[must_use]
pub(crate) fn is_ignored_fgumi_only_tag(tag: [u8; 2]) -> bool {
    tag == *SamTag::TC
}

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
        // `tag_value_size` returns the declared size without confirming the buffer
        // is long enough; reject entries that would slice past the end of the aux
        // data so downstream comparisons cannot panic on truncated records.
        if val_end > aux_data.len() {
            return None;
        }
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
/// Duplicate tag names are handled as a *multiset*: each `r1` entry consumes a distinct
/// unused matching `r2` entry, so `[NM=1, NM=2]` never compares equal to `[NM=2]` even
/// though both share the same tag *name* — matching the SAM spec's uniqueness expectation
/// while still failing loudly on a malformed record that violates it.
///
/// Returns `false` if either record has malformed aux data.
#[must_use]
pub fn raw_tags_equal_order_independent(r1: &[u8], r2: &[u8]) -> bool {
    raw_tags_equal_multiset(r1, r2, is_ignored_fgumi_only_tag)
}

/// [`raw_tags_equal_order_independent`], but additionally excludes the `MI` tag (its value
/// *and* its presence) from the comparison — the multiset equivalent of
/// [`super::engines::content::ContentPredicate::ExactMinusMi`]'s tag rule. `group`
/// legitimately renumbers MI across tools, and MI equivalence is verified separately by the
/// molecule-join engine's matching of molecules on an MI-invariant canonical id (plus its
/// record-membership/strand-partition checks), so this predicate must not also react to an MI
/// difference. Duplicate non-MI tags are still matched as a multiset (a duplicate `NM` does
/// not collapse into one), so a genuine duplicate-tag divergence still reports a difference.
#[must_use]
pub fn raw_tags_equal_order_independent_excluding_mi(r1: &[u8], r2: &[u8]) -> bool {
    raw_tags_equal_multiset(r1, r2, |tag| is_ignored_fgumi_only_tag(tag) || tag == *SamTag::MI)
}

/// Shared multiset tag comparison for the order-independent predicates: drops every tag for
/// which `exclude(tag)` is `true`, then requires the two remaining tag *multisets* to be
/// equal (each `r1` entry pairs with a distinct unused `r2` entry, preserving duplicates).
///
/// Integer tags are compared semantically (decoded to `i64`); every other type is compared
/// by raw bytes. Returns `false` if either record has malformed aux data.
fn raw_tags_equal_multiset(r1: &[u8], r2: &[u8], exclude: impl Fn([u8; 2]) -> bool) -> bool {
    let aux1 = fgumi_raw_bam::fields::aux_data_slice(r1);
    let aux2 = fgumi_raw_bam::fields::aux_data_slice(r2);

    let (Some(entries1), Some(entries2)) = (collect_tag_entries(aux1), collect_tag_entries(aux2))
    else {
        return false;
    };

    // Drop excluded tags (always the fgumi-only `tc`; optionally `MI`) so a record carrying
    // one still compares equal to a record that lacks it.
    let entries1: Vec<_> = entries1.into_iter().filter(|(t, _, _)| !exclude(*t)).collect();
    let entries2: Vec<_> = entries2.into_iter().filter(|(t, _, _)| !exclude(*t)).collect();

    if entries1.len() != entries2.len() {
        return false;
    }

    // For each tag in r1, consume an unused r2 entry that matches it by *both* name and
    // value. Searching for a value match (not just the first same-name entry) is required
    // for duplicate tag names: `[NM=1, NM=2]` must still match a reordered `[NM=2, NM=1]`,
    // which a first-same-name-wins match would wrongly reject (it would pair `NM=1` with
    // `NM=2`, mismatch, and bail). `matched` prevents two r1 entries sharing one r2 entry.
    let mut matched = vec![false; entries2.len()];
    for &(tag1, start1, end1) in &entries1 {
        let candidate = entries2.iter().enumerate().find(|&(i, &(tag2, start2, end2))| {
            !matched[i] && tag2 == tag1 && tag_values_equal(aux1, start1, end1, aux2, start2, end2)
        });
        match candidate {
            Some((idx2, _)) => matched[idx2] = true,
            None => return false,
        }
    }

    true
}

/// Returns `true` if the tag value at `aux1[start1..end1]` equals the one at
/// `aux2[start2..end2]` under the same rule the multiset match uses: byte-identical
/// (type byte included), or — when *both* are integer types — semantically equal after
/// decoding to `i64` (so different on-disk widths of the same value compare equal).
fn tag_values_equal(
    aux1: &[u8],
    start1: usize,
    end1: usize,
    aux2: &[u8],
    start2: usize,
    end2: usize,
) -> bool {
    // Fast path: byte-identical values (includes the type byte).
    if aux1[start1..end1] == aux2[start2..end2] {
        return true;
    }
    // Otherwise, only two integer tags can still be equal — by decoded value.
    let (type1, type2) = (aux1[start1], aux2[start2]);
    if is_int_type(type1) && is_int_type(type2) {
        let v1 = decode_int_tag(type1, &aux1[start1 + 1..end1]);
        let v2 = decode_int_tag(type2, &aux2[start2 + 1..end2]);
        matches!((v1, v2), (Some(v1), Some(v2)) if v1 == v2)
    } else {
        false
    }
}

/// A canonical byte encoding of a record's content under
/// [`ContentPredicate::Exact`](super::engines::content::ContentPredicate::Exact) equality:
/// two records produce byte-equal keys **iff**
/// [`content_diffs(.., Exact, ..)`](super::engines::content::content_diffs) reports them
/// content-equal. This lets an equal-sort-key run be multiset-compared by cancelling each
/// record against its counterpart as it arrives (hashing the keys) instead of buffering and
/// pairwise-matching the run — see the run canceller in
/// [`sort_verify`](super::engines::sort_verify), whose per-key byte-equality cancellation is
/// exactly this `Exact` relation, bounded by order divergence rather than run length.
///
/// The encoding mirrors `Exact`'s equality decision
/// (`raw_core_fields_equal && raw_tags_equal_order_independent`) and its three tolerances
/// exactly:
/// - the BAM `bin` field (bytes 10-11) is excluded, matching [`raw_core_fields_equal`];
/// - aux tags are emitted sorted by their canonical `(name, value)` form, so tag *order*
///   is not significant — as a **multiset** (duplicate tag names are preserved, never
///   collapsed, matching [`raw_tags_equal_order_independent`]);
/// - integer tags are width-normalized to their decoded `i64` value (so `c`/`C`/`s`/`S`/
///   `i`/`I` encodings of the same value collide); the fgumi-only `tc` tag is dropped.
///
/// A record whose core offset or aux data cannot be parsed is *malformed*: under `Exact`
/// it is only ever equal to a byte-identical record (the parse-failure path in
/// `raw_tags_equal_order_independent` never reports two distinct malformed records equal,
/// and `content_diffs`' one match is its byte-equal fast path). Such a record is encoded
/// by its full raw bytes under a distinct leading discriminant, so it collides only with a
/// byte-identical twin and never with a well-formed record's structured key.
///
/// The `content_key_exact_matches_content_diffs_*` proptests in this module's tests lock
/// this equivalence in.
#[must_use]
pub(crate) fn content_key_exact(record: &[u8]) -> Vec<u8> {
    content_key_exact_structured(record).unwrap_or_else(|| {
        // Malformed record: encode by full raw bytes under the raw-fallback discriminant so
        // it collides only with a byte-identical twin.
        let mut key = Vec::with_capacity(record.len() + 1);
        key.push(1u8);
        key.extend_from_slice(record);
        key
    })
}

/// Builds the structured (well-formed) content key for [`content_key_exact`], or `None`
/// when the record's core offset or aux data cannot be parsed (caller then uses the raw
/// fallback).
fn content_key_exact_structured(record: &[u8]) -> Option<Vec<u8>> {
    let aux_off = aux_data_offset_from_record(record)?;
    let end = aux_off.min(record.len());
    // Guard the bin-field slice (bytes 10-11); a record too short for it is degenerate and
    // routes to the raw fallback (where it matches only byte-identical records, which is
    // exactly what `raw_core_fields_equal` — returning `false` here — would decide anyway).
    if record.len() < BIN_OFFSET + BIN_LEN || end < BIN_OFFSET + BIN_LEN {
        return None;
    }

    let aux = fgumi_raw_bam::fields::aux_data_slice(record);
    let entries = collect_tag_entries(aux)?;

    // One canonical (name, normalized value) tuple per non-`tc` tag; sorting the tuples
    // makes the comparison order-independent while preserving duplicates (a multiset).
    let mut tags: Vec<Vec<u8>> = Vec::with_capacity(entries.len());
    for (tag, val_start, val_end) in entries {
        if is_ignored_fgumi_only_tag(tag) {
            continue;
        }
        let type_byte = aux[val_start];
        let mut tuple = Vec::with_capacity(3 + (val_end - val_start));
        tuple.extend_from_slice(&tag);
        if is_int_type(type_byte) {
            // Width-normalize integer tags to their decoded i64 (matching the semantic
            // integer comparison in `raw_tags_equal_multiset`).
            let value = decode_int_tag(type_byte, &aux[val_start + 1..val_end])?;
            tuple.push(1u8);
            tuple.extend_from_slice(&value.to_le_bytes());
        } else {
            // Non-integer tags compare by raw bytes including the type byte (as
            // `raw_tags_equal_multiset` does): `aux[val_start..val_end]`.
            tuple.push(0u8);
            tuple.extend_from_slice(&aux[val_start..val_end]);
        }
        tags.push(tuple);
    }
    tags.sort_unstable();

    // Core (bin excluded) then the tag multiset, every variable-length region explicitly
    // length-prefixed so distinct (core, tags) pairs can never serialize to equal bytes.
    let core_len = BIN_OFFSET + (end - (BIN_OFFSET + BIN_LEN));
    let mut key = Vec::with_capacity(1 + 4 + core_len + 4 + tags.len() * 8);
    key.push(0u8);
    key.extend_from_slice(&(core_len as u32).to_le_bytes());
    key.extend_from_slice(&record[..BIN_OFFSET]);
    key.extend_from_slice(&record[BIN_OFFSET + BIN_LEN..end]);
    key.extend_from_slice(&(tags.len() as u32).to_le_bytes());
    for tuple in &tags {
        key.extend_from_slice(&(tuple.len() as u32).to_le_bytes());
        key.extend_from_slice(tuple);
    }
    Some(key)
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
        aux1.extend_from_slice(&make_i_tag(*SamTag::NM, 5));
        let mut aux2 = make_i_tag(*SamTag::NM, 5);
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
        aux.extend_from_slice(&make_i_tag(*SamTag::NM, 5));
        let r1 = base_record(&aux);
        let r2 = base_record(&aux);
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn tc_carveout_ignores_fgumi_only_tc_tag() {
        // r1 (fgumi) carries the fgumi-only `tc` tag that r2 (fgbio) lacks — they must
        // still compare equal, since fgbio never persists `tc`.
        let mut aux1 = make_i_tag(*SamTag::NM, 5);
        aux1.extend_from_slice(&make_i_tag(*SamTag::TC, 99));
        let r1 = base_record(&aux1);
        let r2 = base_record(&make_i_tag(*SamTag::NM, 5));
        assert!(raw_tags_equal_order_independent(&r1, &r2));
        assert!(is_ignored_fgumi_only_tag(*SamTag::TC));
        assert!(!is_ignored_fgumi_only_tag(*SamTag::NM));
    }

    #[test]
    fn tc_carveout_is_narrow_other_tag_diff_still_differs() {
        // The carve-out drops only `tc`: a real difference in any other tag still DIFFERs,
        // even when an ignored `tc` is also present.
        let mut aux1 = make_i_tag(*SamTag::NM, 5);
        aux1.extend_from_slice(&make_i_tag(*SamTag::TC, 99));
        let r1 = base_record(&aux1);
        let r2 = base_record(&make_i_tag(*SamTag::NM, 6));
        assert!(!raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_different_order() {
        let mut aux1 = make_z_tag(*SamTag::RG, b"s1");
        aux1.extend_from_slice(&make_i_tag(*SamTag::NM, 5));
        let mut aux2 = make_i_tag(*SamTag::NM, 5);
        aux2.extend_from_slice(&make_z_tag(*SamTag::RG, b"s1"));
        let r1 = base_record(&aux1);
        let r2 = base_record(&aux2);
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_different_values() {
        let r1 = base_record(&make_i_tag(*SamTag::NM, 5));
        let r2 = base_record(&make_i_tag(*SamTag::NM, 10));
        assert!(!raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_duplicate_tag_is_multiset_not_collapsed() {
        // Regression: a name-keyed map would collapse `[NM=1, NM=2]` to a single `NM`, so it
        // would compare equal to `[NM=2]` (a false MATCH). The multiset match keeps both
        // `NM` entries distinct: the two records have different tag counts and must DIFFER.
        let mut aux_dup = make_i_tag(*SamTag::NM, 1);
        aux_dup.extend_from_slice(&make_i_tag(*SamTag::NM, 2));
        let r_dup = base_record(&aux_dup);
        let r_single = base_record(&make_i_tag(*SamTag::NM, 2));
        assert!(!raw_tags_equal_order_independent(&r_dup, &r_single));
        // Two records that both carry the same duplicate multiset are still equal.
        let r_dup2 = base_record(&aux_dup);
        assert!(raw_tags_equal_order_independent(&r_dup, &r_dup2));
    }

    #[test]
    fn test_order_independent_reordered_duplicate_values_match() {
        // Regression: matching the *first* unused same-name entry (before comparing values)
        // would pair `NM=1` with `NM=2`, mismatch, and wrongly report DIFFER. `[NM=1, NM=2]`
        // must match a reordered `[NM=2, NM=1]` — the match searches for a value-equal
        // same-name entry, not merely the first same-name one.
        let ab = base_record(&[make_i_tag(*SamTag::NM, 1), make_i_tag(*SamTag::NM, 2)].concat());
        let ba = base_record(&[make_i_tag(*SamTag::NM, 2), make_i_tag(*SamTag::NM, 1)].concat());
        assert!(raw_tags_equal_order_independent(&ab, &ba));
        // But `[NM=1, NM=2]` must still NOT match `[NM=1, NM=1]` (a genuine multiset diff).
        let aa = base_record(&[make_i_tag(*SamTag::NM, 1), make_i_tag(*SamTag::NM, 1)].concat());
        assert!(!raw_tags_equal_order_independent(&ab, &aa));
    }

    // ========================================================================
    // raw_tags_equal_order_independent_excluding_mi tests
    // ========================================================================

    #[test]
    fn excluding_mi_ignores_mi_value_and_presence() {
        // MI value differs -> still equal (MI is excluded entirely).
        let mut aux1 = make_i_tag(*SamTag::NM, 5);
        aux1.extend_from_slice(&make_i_tag(*SamTag::MI, 1));
        let mut aux2 = make_i_tag(*SamTag::NM, 5);
        aux2.extend_from_slice(&make_i_tag(*SamTag::MI, 2));
        let r1 = base_record(&aux1);
        let r2 = base_record(&aux2);
        assert!(raw_tags_equal_order_independent_excluding_mi(&r1, &r2));
        // MI present on one side only -> still equal.
        let r3 = base_record(&make_i_tag(*SamTag::NM, 5));
        assert!(raw_tags_equal_order_independent_excluding_mi(&r1, &r3));
    }

    #[test]
    fn excluding_mi_still_flags_non_mi_tag_difference() {
        // A real non-MI difference must still DIFFER even when MI is present and equal.
        let mut aux1 = make_i_tag(*SamTag::NM, 5);
        aux1.extend_from_slice(&make_i_tag(*SamTag::MI, 1));
        let mut aux2 = make_i_tag(*SamTag::NM, 6);
        aux2.extend_from_slice(&make_i_tag(*SamTag::MI, 1));
        let r1 = base_record(&aux1);
        let r2 = base_record(&aux2);
        assert!(!raw_tags_equal_order_independent_excluding_mi(&r1, &r2));
    }

    #[test]
    fn excluding_mi_duplicate_non_mi_tag_is_multiset_not_collapsed() {
        // The MI-excluding predicate must keep the same multiset discipline for non-MI tags:
        // `[NM=1, NM=2]` must not collapse to a single `NM` and false-MATCH `[NM=2]`.
        let mut aux_dup = make_i_tag(*SamTag::NM, 1);
        aux_dup.extend_from_slice(&make_i_tag(*SamTag::NM, 2));
        let r_dup = base_record(&aux_dup);
        let r_single = base_record(&make_i_tag(*SamTag::NM, 2));
        assert!(!raw_tags_equal_order_independent_excluding_mi(&r_dup, &r_single));
    }

    #[test]
    fn test_order_independent_missing_tag() {
        let mut aux1 = make_z_tag(*SamTag::RG, b"s1");
        aux1.extend_from_slice(&make_i_tag(*SamTag::NM, 5));
        let aux2 = make_z_tag(*SamTag::RG, b"s1");
        let r1 = base_record(&aux1);
        let r2 = base_record(&aux2);
        assert!(!raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_extra_tag() {
        let aux1 = make_i_tag(*SamTag::NM, 5);
        let mut aux2 = make_i_tag(*SamTag::NM, 5);
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
        aux1.extend_from_slice(&make_i_tag(*SamTag::NM, 3));
        aux1.extend_from_slice(&make_c_tag(*SamTag::MQ, 42));
        aux1.extend_from_slice(&make_f_tag(*b"GC", 0.75));

        let mut aux2 = make_f_tag(*b"GC", 0.75);
        aux2.extend_from_slice(&make_c_tag(*SamTag::MQ, 42));
        aux2.extend_from_slice(&make_z_tag(*SamTag::RG, b"grp"));
        aux2.extend_from_slice(&make_i_tag(*SamTag::NM, 3));

        let r1 = base_record(&aux1);
        let r2 = base_record(&aux2);
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_short_record() {
        let short = vec![0u8; 10];
        let rec = base_record(&make_i_tag(*SamTag::NM, 5));
        // Both have empty aux slices due to short record, but the short record
        // returns empty aux, so it won't match a record with tags.
        assert!(!raw_tags_equal_order_independent(&short, &rec));
    }

    // ========================================================================
    // raw_compare_structured tests
    // ========================================================================

    #[test]
    fn test_structured_identical() {
        let aux = make_i_tag(*SamTag::NM, 5);
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
        aux1.extend_from_slice(&make_i_tag(*SamTag::NM, 5));
        let mut aux2 = make_i_tag(*SamTag::NM, 5);
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
        let r1 = base_record(&make_i_tag(*SamTag::NM, 5));
        let r2 = base_record(&make_i_tag(*SamTag::NM, 10));
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
        aux1.extend_from_slice(&make_i_tag(*SamTag::NM, 5));
        let mut aux2 = make_i_tag(*SamTag::NM, 5);
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
        let aux = make_i_tag(*SamTag::NM, 42);
        let entries = collect_tag_entries(&aux)
            .expect("collect_tag_entries should succeed for single int tag");
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].0, *SamTag::NM);
    }

    #[test]
    fn test_collect_entries_multiple() {
        let mut aux = make_z_tag(*SamTag::RG, b"grp");
        aux.extend_from_slice(&make_i_tag(*SamTag::NM, 5));
        let entries = collect_tag_entries(&aux)
            .expect("collect_tag_entries should succeed for multiple tags");
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].0, *SamTag::RG);
        assert_eq!(entries[1].0, *SamTag::NM);
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

    #[test]
    fn test_collect_entries_truncated_fixed_size_returns_none() {
        // `i` type (i32) declares 4 bytes of value but only 2 are present in the slice.
        // Without bounds enforcement on val_end, this produces an entry with
        // val_end > aux_data.len() and slicing by a downstream caller panics.
        let aux = [b'N', b'M', b'i', 0x00, 0x01];
        assert!(collect_tag_entries(&aux).is_none());
    }

    #[test]
    fn test_collect_entries_truncated_b_array_returns_none() {
        // B-array declares 5 i32 elements but the buffer only has room for 2 after the header.
        // `tag_value_size` returns 5 + 5*4 = 25 which overshoots the slice.
        let mut aux = vec![b'X', b'A', b'B', b'i']; // B array of i32
        aux.extend_from_slice(&5u32.to_le_bytes()); // claim 5 elements
        aux.extend_from_slice(&1i32.to_le_bytes()); // only 2 actually follow
        aux.extend_from_slice(&2i32.to_le_bytes());
        assert!(collect_tag_entries(&aux).is_none());
    }

    #[test]
    fn test_order_independent_truncated_fixed_size_returns_false() {
        // Regression: `raw_tags_equal_order_independent` used to panic when aux data was
        // truncated for a fixed-size type because `collect_tag_entries` produced entries
        // whose `val_end` exceeded `aux_data.len()`. It must return `false` instead.
        let aux_good = make_i_tag(*SamTag::NM, 42);
        let aux_bad = [b'N', b'M', b'i', 0x00, 0x01]; // truncated i32 (only 2 bytes)
        let r1 = base_record(&aux_good);
        let r2 = base_record(&aux_bad);
        assert!(!raw_tags_equal_order_independent(&r1, &r2));
        // Also when both sides are truncated: malformed aux => false.
        let r3 = base_record(&aux_bad);
        assert!(!raw_tags_equal_order_independent(&r2, &r3));
    }

    // ========================================================================
    // Semantic integer comparison tests
    // ========================================================================

    #[test]
    fn test_order_independent_same_value_different_int_types_u8_vs_i16() {
        // CD tag: u8 (C) value 158 vs i16 (s) value 158
        // These are semantically equal but have different byte encodings.
        let r1 = base_record(&make_c_tag(*SamTag::CD, 158));
        let r2 = base_record(&make_s_tag(*SamTag::CD, 158));
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_same_value_different_int_types_u8_vs_i32() {
        // NM tag: u8 (C) value 42 vs i32 (i) value 42
        let r1 = base_record(&make_c_tag(*SamTag::NM, 42));
        let r2 = base_record(&make_i_tag(*SamTag::NM, 42));
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_same_value_different_int_types_i16_vs_i32() {
        // MQ tag: i16 (s) value 300 vs i32 (i) value 300
        let r1 = base_record(&make_s_tag(*SamTag::MQ, 300));
        let r2 = base_record(&make_i_tag(*SamTag::MQ, 300));
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_same_value_different_int_types_u32_vs_i32() {
        // NM tag: u32 (I) value 1000 vs i32 (i) value 1000
        let r1 = base_record(&make_upper_i_tag(*SamTag::NM, 1000));
        let r2 = base_record(&make_i_tag(*SamTag::NM, 1000));
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_same_value_different_int_types_u32_vs_u8() {
        // NM tag: u32 (I) value 42 vs u8 (C) value 42
        let r1 = base_record(&make_upper_i_tag(*SamTag::NM, 42));
        let r2 = base_record(&make_c_tag(*SamTag::NM, 42));
        assert!(raw_tags_equal_order_independent(&r1, &r2));
    }

    #[test]
    fn test_order_independent_different_values_different_int_types() {
        // Same tag name, different int types, different values => should NOT match
        let r1 = base_record(&make_c_tag(*SamTag::NM, 5));
        let r2 = base_record(&make_i_tag(*SamTag::NM, 10));
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
        let r1 = base_record(&make_c_tag(*SamTag::CD, 158));
        let r2 = base_record(&make_s_tag(*SamTag::CD, 158));
        let result = raw_compare_structured(&r1, &r2);
        assert!(result.core_match);
        assert!(!result.tags_match); // raw bytes differ
        assert!(result.tag_order_match); // semantically equal
    }

    #[test]
    fn test_order_independent_semantic_int_with_reordered_tags() {
        // Multiple tags, different order, AND different integer encodings
        let mut aux1 = make_z_tag(*SamTag::RG, b"grp");
        aux1.extend_from_slice(&make_c_tag(*SamTag::CD, 200));
        aux1.extend_from_slice(&make_i_tag(*SamTag::NM, 5));

        let mut aux2 = make_i_tag(*SamTag::NM, 5);
        aux2.extend_from_slice(&make_upper_s_tag(*SamTag::CD, 200));
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

    // ========================================================================
    // content_key_exact: canonical content key mirrors content_diffs(Exact)
    // ========================================================================

    use crate::commands::compare::engines::content::{ContentPredicate, content_diffs};
    use noodles::sam::Header;
    use proptest::prelude::*;

    /// `content_diffs(a, b, Exact, ..).is_none()` — the exact equality relation
    /// `content_key_exact` must mirror. Headers are irrelevant to the Exact *decision*
    /// (they only render diff strings), so a default header is used on both sides.
    fn exact_equal(a: &[u8], b: &[u8]) -> bool {
        content_diffs(a, b, ContentPredicate::Exact, &Header::default(), &Header::default())
            .is_none()
    }

    #[test]
    fn content_key_equal_iff_exact_equal_for_reordered_tags() {
        // Same two tags in different on-record order: Exact-equal, so keys must match.
        let r1 =
            base_record(&[make_i_tag(*SamTag::NM, 1), make_z_tag(*SamTag::RG, b"grp")].concat());
        let r2 =
            base_record(&[make_z_tag(*SamTag::RG, b"grp"), make_i_tag(*SamTag::NM, 1)].concat());
        assert!(exact_equal(&r1, &r2));
        assert_eq!(content_key_exact(&r1), content_key_exact(&r2));
    }

    #[test]
    fn content_key_normalizes_integer_tag_width() {
        // Same NM value, different on-disk integer width (i32 vs i16): Exact-equal, so the
        // canonical keys must collide even though the records are not byte-equal.
        let wide = base_record(&make_i_tag(*SamTag::NM, 5));
        let narrow = base_record(&make_s_tag(*SamTag::NM, 5));
        assert_ne!(wide, narrow, "records must differ on disk (different widths)");
        assert!(exact_equal(&wide, &narrow));
        assert_eq!(content_key_exact(&wide), content_key_exact(&narrow));
    }

    #[test]
    fn content_key_drops_fgumi_only_tc_tag() {
        // The fgumi-only `tc` tag is ignored by Exact; a record carrying it keys the same
        // as one without it.
        let with_tc =
            base_record(&[make_i_tag(*SamTag::NM, 1), make_i_tag(*SamTag::TC, 9)].concat());
        let without_tc = base_record(&make_i_tag(*SamTag::NM, 1));
        assert!(exact_equal(&with_tc, &without_tc));
        assert_eq!(content_key_exact(&with_tc), content_key_exact(&without_tc));
    }

    #[test]
    fn content_key_excludes_bin_field() {
        // The BAM `bin` field (bytes 10-11) is excluded from Exact core comparison; two
        // records differing only there must key identically.
        let mut a = base_record(&make_i_tag(*SamTag::NM, 1));
        let mut b = a.clone();
        b[10..12].copy_from_slice(&0xABCDu16.to_le_bytes());
        assert_ne!(a, b, "records must differ in the bin field");
        assert!(exact_equal(&a, &b));
        assert_eq!(content_key_exact(&a), content_key_exact(&b));
        // Sanity: a real core difference (pos) still separates the keys.
        a[4..8].copy_from_slice(&999i32.to_le_bytes());
        assert!(!exact_equal(&a, &b));
        assert_ne!(content_key_exact(&a), content_key_exact(&b));
    }

    #[test]
    fn content_key_duplicate_tags_are_a_multiset() {
        // [NM=1, NM=2] must not collapse: it is Exact-equal to a reordered [NM=2, NM=1] but
        // NOT to a single [NM=2].
        let ab = base_record(&[make_i_tag(*SamTag::NM, 1), make_i_tag(*SamTag::NM, 2)].concat());
        let ba = base_record(&[make_i_tag(*SamTag::NM, 2), make_i_tag(*SamTag::NM, 1)].concat());
        let single = base_record(&make_i_tag(*SamTag::NM, 2));
        assert_eq!(content_key_exact(&ab), content_key_exact(&ba));
        assert!(exact_equal(&ab, &ba));
        assert_ne!(content_key_exact(&ab), content_key_exact(&single));
        assert!(!exact_equal(&ab, &single));
    }

    #[test]
    fn content_key_malformed_matches_only_byte_identical() {
        // A truncated (malformed) aux region: content_diffs can only call two such records
        // equal via its byte-equal fast path, and content_key_exact's raw fallback matches
        // the same — byte-identical twins collide, a differing malformed record does not.
        let malformed = base_record(&[b'N', b'M', b'i', 0x01]); // 'i' claims 4 bytes, only 1 present
        let same = malformed.clone();
        let mut different = malformed.clone();
        *different.last_mut().unwrap() ^= 0xFF;
        assert!(exact_equal(&malformed, &same));
        assert_eq!(content_key_exact(&malformed), content_key_exact(&same));
        assert!(!exact_equal(&malformed, &different));
        assert_ne!(content_key_exact(&malformed), content_key_exact(&different));
    }

    // ---- proptest: content_key_exact byte-equality <=> content_diffs(Exact) equality ----

    /// One generated aux tag. Encoded to raw BAM aux bytes by [`encode_tag`]; integer tags
    /// can be re-encoded at a wider type by [`encode_tag_wide`] to exercise
    /// width-normalization.
    #[derive(Clone, Debug)]
    enum TagGen {
        /// Integer tag (small value so it fits `s`/`i16`), name from a tiny pool.
        Int { name: [u8; 2], value: i16 },
        /// Z-type string tag.
        Str { name: [u8; 2], value: Vec<u8> },
        /// f-type float tag (compared by raw bytes under Exact).
        Float { name: [u8; 2], bits: u32 },
    }

    fn tag_gen_strategy() -> impl Strategy<Value = TagGen> {
        // A tiny tag-name pool, deliberately including duplicate-prone names, `MI` (compared
        // under Exact) and the ignored fgumi-only `tc`. Built at runtime because `*SamTag`
        // deref coercion is not permitted in a `const`.
        let names: [[u8; 2]; 6] =
            [*SamTag::NM, *SamTag::MD, *SamTag::RG, *SamTag::MI, *SamTag::TC, *b"XA"];
        let name = (0..names.len()).prop_map(move |i| names[i]);
        prop_oneof![
            (name.clone(), any::<i16>()).prop_map(|(name, value)| TagGen::Int { name, value }),
            (
                name.clone(),
                proptest::collection::vec(any::<u8>().prop_filter("no NUL", |b| *b != 0), 0..4)
            )
                .prop_map(|(name, value)| TagGen::Str { name, value }),
            (name, any::<u32>()).prop_map(|(name, bits)| TagGen::Float { name, bits }),
        ]
    }

    /// Encode a tag to raw aux bytes at its "natural" width (integers as `s`/i16).
    fn encode_tag(t: &TagGen) -> Vec<u8> {
        match t {
            TagGen::Int { name, value } => make_s_tag(*name, *value),
            TagGen::Str { name, value } => make_z_tag(*name, value),
            TagGen::Float { name, bits } => make_f_tag(*name, f32::from_bits(*bits)),
        }
    }

    /// Encode a tag, widening integers to `i`/i32 — a different on-disk width for the same
    /// value, so a re-encoded record is content-equal but not byte-equal.
    fn encode_tag_wide(t: &TagGen) -> Vec<u8> {
        match t {
            TagGen::Int { name, value } => make_i_tag(*name, i32::from(*value)),
            other => encode_tag(other),
        }
    }

    fn aux_from(tags: &[TagGen], wide: bool) -> Vec<u8> {
        tags.iter().flat_map(|t| if wide { encode_tag_wide(t) } else { encode_tag(t) }).collect()
    }

    proptest! {
        /// The core equivalence: for two independently generated records (same core), the
        /// canonical keys are byte-equal **iff** `content_diffs(.., Exact, ..)` calls them
        /// content-equal. This is the property the sort-verify run canceller's key-based
        /// cancellation relies on.
        #[test]
        fn content_key_exact_matches_content_diffs_random(
            tags_a in proptest::collection::vec(tag_gen_strategy(), 0..6),
            tags_b in proptest::collection::vec(tag_gen_strategy(), 0..6),
            wide_a in any::<bool>(),
            wide_b in any::<bool>(),
        ) {
            let ra = base_record(&aux_from(&tags_a, wide_a));
            let rb = base_record(&aux_from(&tags_b, wide_b));
            prop_assert_eq!(
                content_key_exact(&ra) == content_key_exact(&rb),
                exact_equal(&ra, &rb),
                "key-equality must match content_diffs(Exact) for {:?} vs {:?}",
                tags_a,
                tags_b
            );
        }

        /// A record and a *reordered, integer-widened* copy of the same tags are always
        /// content-equal, so their keys must always collide — the reorder- and
        /// width-independence guarantees, fuzzed.
        #[test]
        fn content_key_exact_stable_under_reorder_and_width(
            tags in proptest::collection::vec(tag_gen_strategy(), 0..6),
            rotate in 0usize..6,
        ) {
            let mut reordered = tags.clone();
            if !reordered.is_empty() {
                let shift = rotate % reordered.len();
                reordered.rotate_right(shift);
            }
            let base = base_record(&aux_from(&tags, false));
            let variant = base_record(&aux_from(&reordered, true));
            prop_assert!(exact_equal(&base, &variant));
            prop_assert_eq!(content_key_exact(&base), content_key_exact(&variant));
        }
    }
}
