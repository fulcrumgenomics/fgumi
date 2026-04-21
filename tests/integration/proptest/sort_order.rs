//! Property: the template-coordinate sort key induces a total order.
//!
//! `TemplateKey` is the fixed-width sort key used by fgumi's external sorter
//! when writing template-coordinate BAMs. A sort key must, at minimum, be a
//! total order: for every pair `(a, b)`, exactly one of `a < b`, `a == b`,
//! or `a > b` must hold, and `cmp` must be consistent with `eq`.

use std::cmp::Ordering;

use fgumi_lib::sort::TemplateKey;
use proptest::prelude::*;

/// Build a `TemplateKey` from a tuple of the fields that the production
/// `extract_template_key_inline` would pack out of a BAM record.
#[allow(clippy::too_many_arguments, clippy::fn_params_excessive_bools)]
fn make_key(
    tid1: i32,
    pos1: i32,
    neg1: bool,
    tid2: i32,
    pos2: i32,
    neg2: bool,
    cb_hash: u64,
    library: u32,
    mi_value: u64,
    mi_is_a: bool,
    name_hash: u64,
    is_upper: bool,
) -> TemplateKey {
    TemplateKey::new(
        tid1,
        pos1,
        neg1,
        tid2,
        pos2,
        neg2,
        cb_hash,
        library,
        (mi_value, mi_is_a),
        name_hash,
        is_upper,
    )
}

/// Strategy covering the full relevant domain of `TemplateKey` fields,
/// including edge cases: unmapped sentinels (`tid == i32::MAX`), both strands,
/// full `u64` ranges for hashes, and both MI strand suffixes.
fn template_key_strategy() -> impl Strategy<Value = TemplateKey> {
    (
        prop_oneof![Just(i32::MAX), 0i32..8],
        -1i32..10_000,
        any::<bool>(),
        prop_oneof![Just(i32::MAX), 0i32..8],
        -1i32..10_000,
        any::<bool>(),
        any::<u64>(),
        0u32..16,
        0u64..(1 << 40),
        any::<bool>(),
        any::<u64>(),
        any::<bool>(),
    )
        .prop_map(
            |(
                tid1,
                pos1,
                neg1,
                tid2,
                pos2,
                neg2,
                cb_hash,
                library,
                mi_value,
                mi_is_a,
                name_hash,
                is_upper,
            )| {
                make_key(
                    tid1, pos1, neg1, tid2, pos2, neg2, cb_hash, library, mi_value, mi_is_a,
                    name_hash, is_upper,
                )
            },
        )
}

proptest! {
    /// Trichotomy: exactly one of `<`, `==`, `>` holds.
    #[test]
    fn sort_key_cmp_is_trichotomous(
        a in template_key_strategy(),
        b in template_key_strategy(),
    ) {
        let ord = a.cmp(&b);
        let less = ord == Ordering::Less;
        let equal = ord == Ordering::Equal;
        let greater = ord == Ordering::Greater;
        // Exactly one is true.
        prop_assert_eq!(u8::from(less) + u8::from(equal) + u8::from(greater), 1);

        // `cmp == Equal` iff the structs are `==`.
        prop_assert_eq!(equal, a == b);
    }

    /// Antisymmetry: `cmp(a, b)` is the reverse of `cmp(b, a)`.
    #[test]
    fn sort_key_cmp_is_antisymmetric(
        a in template_key_strategy(),
        b in template_key_strategy(),
    ) {
        prop_assert_eq!(a.cmp(&b), b.cmp(&a).reverse());
    }

    /// Transitivity: if `a <= b` and `b <= c` then `a <= c`.
    #[test]
    fn sort_key_cmp_is_transitive(
        a in template_key_strategy(),
        b in template_key_strategy(),
        c in template_key_strategy(),
    ) {
        if a.cmp(&b) != Ordering::Greater && b.cmp(&c) != Ordering::Greater {
            prop_assert!(a.cmp(&c) != Ordering::Greater);
        }
    }

    /// Reflexivity: every key is equal to itself.
    #[test]
    fn sort_key_cmp_is_reflexive(a in template_key_strategy()) {
        prop_assert_eq!(a.cmp(&a), Ordering::Equal);
    }
}
