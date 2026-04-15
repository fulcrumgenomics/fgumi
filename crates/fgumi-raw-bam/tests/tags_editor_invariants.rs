#![cfg(feature = "test-utils")]

use fgumi_raw_bam::RawRecord;
use fgumi_raw_bam::testutil::{encode_op, make_bam_bytes};
use proptest::prelude::*;

/// Build a minimal valid [`RawRecord`] with the given aux bytes.
///
/// Uses a 4-base Match CIGAR op and a short name whose length + NUL is
/// word-aligned so the aux section starts at a clean offset.
fn base_record(aux: &[u8]) -> RawRecord {
    RawRecord::from(make_bam_bytes(0, 0, 0, b"r", &[encode_op(0, 4)], 4, -1, -1, aux))
}

proptest! {
    // =========================================================================
    // INVARIANT 1a: append_string(tag, v); remove(tag) → bytes unchanged
    // when `tag` was absent before the append.
    // =========================================================================

    /// A record whose aux section is built exclusively from well-typed integer
    /// tags is left byte-for-byte identical after appending and then removing
    /// a "Xz" tag.
    ///
    /// We generate valid aux data (zero or more i32 tags with distinct
    /// two-uppercase-letter names) rather than random bytes, so the scanner
    /// can always navigate past the existing tags to reach the newly appended
    /// one.  The invariant requires the aux section to be parseable.
    #[test]
    fn append_then_remove_string_is_noop(
        // 0–4 existing i32 tags under uppercase two-letter names (not "Xz")
        existing in proptest::collection::vec(
            ("[A-WY-Z]{2}", any::<i32>()).prop_map(|(t, v)| (t.into_bytes(), v)),
            0..4,
        ),
        val in "[A-Za-z]{0,40}",
    ) {
        let tag = b"Xz";

        // Build valid aux bytes from the generated tags (deduplicating names).
        let mut aux: Vec<u8> = Vec::new();
        let mut seen: Vec<[u8; 2]> = Vec::new();
        for (name, val_i) in &existing {
            let k = [name[0], name[1]];
            // Skip duplicates and any accidental collision with our test tag.
            if seen.contains(&k) || &k == tag {
                continue;
            }
            seen.push(k);
            // Write as i32 (type 'i', 4-byte LE value).
            aux.push(k[0]);
            aux.push(k[1]);
            aux.push(b'i');
            aux.extend_from_slice(&val_i.to_le_bytes());
        }

        let mut rec = base_record(&aux);
        let before = rec.as_ref().to_vec();

        {
            let mut ed = rec.tags_editor();
            ed.append_string(tag, val.as_bytes());
            ed.remove(tag);
        }

        prop_assert_eq!(rec.as_ref(), before.as_slice());
    }

    // =========================================================================
    // INVARIANT 1b: append_int(tag, v); remove(tag) → bytes unchanged
    // when `tag` was absent before the append.
    // =========================================================================

    /// Same noop invariant as 1a, exercised over the full i32 domain.
    #[test]
    fn append_then_remove_int_is_noop(v in any::<i32>()) {
        let tag = b"Xz";
        let mut rec = base_record(&[]);
        let before = rec.as_ref().to_vec();

        {
            let mut ed = rec.tags_editor();
            ed.append_int(tag, v);
            ed.remove(tag);
        }

        prop_assert_eq!(rec.as_ref(), before.as_slice());
    }

    // =========================================================================
    // INVARIANT 2a: update_int round-trip across i32 domain.
    // =========================================================================

    /// Writing an integer tag and reading it back must yield the same value
    /// (widened to i64) on every i32 input.
    #[test]
    fn update_int_roundtrip(v in any::<i32>()) {
        let tag = b"NM";
        let mut rec = base_record(&[]);
        {
            let mut ed = rec.tags_editor();
            ed.update_int(tag, v);
        }
        let found = rec.tags().find_int(tag);
        prop_assert_eq!(found, Some(i64::from(v)));
    }

    // =========================================================================
    // INVARIANT 2b: update_string round-trip.
    // =========================================================================

    /// Writing a string tag and reading it back must yield the identical bytes.
    #[test]
    fn update_string_roundtrip(s in "[A-Za-z0-9 _-]{0,200}") {
        let tag = b"Zz";
        let mut rec = base_record(&[]);
        {
            let mut ed = rec.tags_editor();
            ed.update_string(tag, s.as_bytes());
        }
        let found = rec.tags().find_string(tag);
        prop_assert_eq!(found, Some(s.as_bytes()));
    }

    // =========================================================================
    // INVARIANT 2c: update_float round-trip across finite f32 domain.
    // =========================================================================

    /// A float is stored as 4 LE bytes and read back verbatim; the bit pattern
    /// must be identical (NaN is excluded because NaN != NaN).
    #[test]
    fn update_float_roundtrip(
        v in any::<f32>().prop_filter("finite only", |x| x.is_finite()),
    ) {
        let tag = b"As";
        let mut rec = base_record(&[]);
        {
            let mut ed = rec.tags_editor();
            ed.update_float(tag, v);
        }
        let found = rec.tags().find_float(tag);
        prop_assert!(found.is_some(), "tag not found after update_float");
        prop_assert_eq!(
            found.unwrap().to_bits(),
            v.to_bits(),
            "float round-trip bit mismatch: wrote {:?}, read back {:?}",
            v,
            found
        );
    }

    // =========================================================================
    // INVARIANT 3: aux_offset stays valid across N mixed ops.
    //
    // After every mutation, tags().iter() must terminate without panic, and
    // every tag entry must have value_bytes within the aux section bounds.
    // =========================================================================

    /// Any sequence of up to 20 update_int / update_string / remove operations
    /// must leave the record in a state where aux-tag iteration terminates and
    /// every reported entry is within bounds.
    #[test]
    fn mixed_ops_preserve_iterability(
        ops in proptest::collection::vec(
            prop_oneof![
                // update_int: kind='I', tag, int_val, str_val unused
                ("[A-Z]{2}", any::<i32>()).prop_map(|(t, v)| {
                    let tb = t.into_bytes();
                    (b'I', [tb[0], tb[1]], v, Vec::<u8>::new())
                }),
                // update_string: kind='S', tag, str_val, int_val unused
                ("[A-Z]{2}", "[a-z]{0,20}").prop_map(|(t, s)| {
                    let tb = t.into_bytes();
                    (b'S', [tb[0], tb[1]], 0i32, s.into_bytes())
                }),
                // remove: kind='R', tag, both values unused
                "[A-Z]{2}".prop_map(|t| {
                    let tb = t.into_bytes();
                    (b'R', [tb[0], tb[1]], 0i32, Vec::<u8>::new())
                }),
            ],
            1..20,
        ),
    ) {
        let mut rec = base_record(&[]);
        let aux_len = rec.as_ref().len();

        for (kind, tag, int_val, str_val) in &ops {
            {
                let mut ed = rec.tags_editor();
                match kind {
                    b'I' => ed.update_int(tag, *int_val),
                    b'S' => ed.update_string(tag, str_val),
                    b'R' => ed.remove(tag),
                    _ => unreachable!(),
                }
            }

            // After the mutation, iteration must terminate and stay in-bounds.
            let tags = rec.tags();
            let total_aux = tags.as_bytes().len();
            let mut count = 0usize;
            for entry in &tags {
                prop_assert!(
                    entry.value_bytes.len() <= total_aux,
                    "entry value_bytes.len() {} exceeds aux section length {}",
                    entry.value_bytes.len(),
                    total_aux
                );
                count += 1;
                prop_assert!(count < 1000, "iter produced >1000 entries — likely infinite loop");
            }

            // The record must not be shorter than it was at the start (header + fixed fields).
            prop_assert!(
                rec.as_ref().len() >= aux_len,
                "record shrank below initial size: {} < {}",
                rec.as_ref().len(),
                aux_len
            );
        }
    }
}
