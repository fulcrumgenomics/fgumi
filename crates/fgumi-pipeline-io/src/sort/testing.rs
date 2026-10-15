//! The sort's test oracle: the stable in-memory reference sort, and the
//! tie-heavy fixture the parity tests feed it.
//!
//! One home for both, so every crate's sort tests check the same contract:
//! `sort/tests.rs` (this crate), the root crate's integration tests
//! (`fgumi merge`, the chain sort) through the `test-utils` feature.
#![allow(
    clippy::missing_panics_doc,
    clippy::must_use_candidate,
    clippy::cast_possible_truncation,
    clippy::should_implement_trait
)]

use rand::RngExt as _;
use std::cmp::Ordering;

use fgumi_raw_bam::RawRecord;
use fgumi_raw_bam::testutil::make_bam_bytes;
use fgumi_sort::{QuerynameComparator, RawSortKey, SortOrder};
use noodles::sam::Header;

/// Every record's bytes, in output order.
pub type RecordBytes = Vec<Vec<u8>>;

/// Number of reference sequences the synthetic header declares.
pub const N_TEST_REFS: usize = 4;

// ── Reference: the stable in-memory sort ────────────────────────────────────

/// One record's key in a sort order — the definition of that order.
pub enum OrderKey {
    Coordinate(fgumi_sort::RawCoordinateKey),
    QuerynameLex(fgumi_sort::RawQuerynameLexKey),
    QuerynameNatural(fgumi_sort::RawQuerynameKey),
    Template(fgumi_sort::TemplateKey),
}

impl OrderKey {
    /// Compare two keys of one order.
    ///
    /// # Panics
    ///
    /// If the keys belong to different orders.
    #[must_use]
    pub fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Self::Coordinate(a), Self::Coordinate(b)) => a.cmp(b),
            (Self::QuerynameLex(a), Self::QuerynameLex(b)) => a.cmp(b),
            (Self::QuerynameNatural(a), Self::QuerynameNatural(b)) => a.cmp(b),
            (Self::Template(a), Self::Template(b)) => a.cmp(b),
            _ => unreachable!("keys of one order"),
        }
    }
}

/// Every record's key in `order`, in record order. Template-coordinate uses
/// the full key: a narrowed lane set orders exactly like the full key.
pub fn order_keys<R: AsRef<[u8]>>(
    order: SortOrder,
    header: &Header,
    records: &[R],
    cell_tag: Option<fgumi_raw_bam::SamTag>,
) -> Vec<OrderKey> {
    let ctx = fgumi_sort::SortContext::from_header(header);
    let lib = fgumi_sort::LibraryLookup::from_header(header);
    let hasher = fgumi_sort::cb_hasher();
    records
        .iter()
        .map(|r| {
            let r = r.as_ref();
            match order {
                SortOrder::Coordinate => {
                    OrderKey::Coordinate(fgumi_sort::RawCoordinateKey::extract(r, &ctx))
                }
                SortOrder::Queryname(QuerynameComparator::Lexicographic) => {
                    OrderKey::QuerynameLex(fgumi_sort::RawQuerynameLexKey::extract(r, &ctx))
                }
                SortOrder::Queryname(QuerynameComparator::Natural) => {
                    OrderKey::QuerynameNatural(fgumi_sort::RawQuerynameKey::extract(r, &ctx))
                }
                SortOrder::TemplateCoordinate => OrderKey::Template(
                    fgumi_sort::extract_template_key_inline(r, &lib, cell_tag, &hasher),
                ),
            }
        })
        .collect()
}

/// In-memory STABLE sort of `records` by the order's key — the definition of
/// the output contract: equal keys keep input order.
pub fn reference_sort(
    order: SortOrder,
    header: &Header,
    records: &[RawRecord],
    cell_tag: Option<fgumi_raw_bam::SamTag>,
) -> RecordBytes {
    let keys = order_keys(order, header, records, cell_tag);
    let mut idx: Vec<usize> = (0..records.len()).collect();
    // `sort_by` is stable.
    idx.sort_by(|&a, &b| keys[a].cmp(&keys[b]));
    idx.into_iter().map(|i| records[i].as_ref().to_vec()).collect()
}

// ── Tie-heavy fixture ───────────────────────────────────────────────────────

/// The input index a [`tie_heavy_records`] record carries in its last aux
/// field, `XI:i:<index>`.
pub fn xi_of(record: &[u8]) -> u32 {
    let n = record.len();
    assert!(n >= 7 && &record[n - 7..n - 4] == b"XIi", "record must end with an XI:i tag");
    u32::from_le_bytes(record[n - 4..].try_into().expect("4 bytes"))
}

/// Adjacent runs of `out` whose records share one key in `order`, as index
/// ranges into `out`.
pub fn equal_key_groups(
    order: SortOrder,
    header: &Header,
    out: &[Vec<u8>],
) -> Vec<std::ops::Range<usize>> {
    let keys = order_keys(order, header, out, None);
    let mut groups = Vec::new();
    let mut start = 0;
    for i in 1..=keys.len() {
        if i == keys.len() || keys[i].cmp(&keys[start]) != Ordering::Equal {
            groups.push(start..i);
            start = i;
        }
    }
    groups
}

/// Within every group of equal keys in `out`, the records keep input order:
/// their `XI` indices strictly ascend.
pub fn assert_input_order_within_equal_keys(order: SortOrder, header: &Header, out: &[Vec<u8>]) {
    for group in equal_key_groups(order, header, out) {
        let xis: Vec<u32> = out[group.clone()].iter().map(|r| xi_of(r)).collect();
        assert!(
            xis.windows(2).all(|w| w[0] < w[1]),
            "{order:?}: equal keys out of input order at output {group:?}: XI {xis:?}"
        );
    }
}

/// The size of the largest group of equal keys in `out`.
pub fn max_equal_key_group(order: SortOrder, header: &Header, out: &[Vec<u8>]) -> usize {
    equal_key_groups(order, header, out).iter().map(ExactSizeIterator::len).max().unwrap_or(0)
}

/// Whether some group of equal keys in `out` has members from two different
/// runs, given the runs' boundaries as input record counts.
pub fn equal_key_group_spans_a_boundary(
    order: SortOrder,
    header: &Header,
    out: &[Vec<u8>],
    spill_boundaries: &[u64],
) -> bool {
    let run_of = |record: &Vec<u8>| {
        let xi = u64::from(xi_of(record));
        spill_boundaries.iter().filter(|&&b| b <= xi).count()
    };
    equal_key_groups(order, header, out).into_iter().any(|group| {
        let first = run_of(&out[group.start]);
        out[group].iter().any(|r| run_of(r) != first)
    })
}

/// The seeded generator every fixture draws from, so fixtures are reproducible.
#[must_use]
pub fn fixture_rng(seed: u64) -> rand::rngs::StdRng {
    rand::SeedableRng::seed_from_u64(seed)
}

/// `true` with probability `percent / 100`.
pub fn chance(rng: &mut impl rand::Rng, percent: u64) -> bool {
    rand::RngExt::random_range(rng, 0..100) < percent
}

/// The header every synthetic fixture uses: `N_TEST_REFS` references of 1 Mb.
pub fn synthetic_header() -> Header {
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReferenceSequence;
    use std::num::NonZeroUsize;
    let len = NonZeroUsize::new(1_000_000).expect("nonzero");
    let mut builder = Header::builder();
    for i in 0..N_TEST_REFS {
        builder =
            builder.add_reference_sequence(format!("chr{i}"), Map::<ReferenceSequence>::new(len));
    }
    builder.build()
}

/// About `n` records, unsorted in every order and rich in equal keys, each
/// ending with a unique `XI:i:<input index>` so equal-key records are never
/// byte-identical and their input order is observable in the output.
///
/// Templates are read pairs with Illumina-like names. Ties come from:
/// - secondaries (30 %) and supplementaries (20 %) — one or two each, same
///   name and flag class; secondaries copy the first read's alignment, so they
///   tie with it in coordinate and template-coordinate order, and with each
///   other in queryname order;
/// - pairs forced onto one shared coordinate (5 %);
/// - whole templates duplicated under a new name at the same coordinates (3 %);
/// - unmapped pairs (2 %), which all share the unmapped coordinate key.
///
/// The records are shuffled before the `XI` tags are assigned, so equal-key
/// records are spread over the input, and therefore over the runs.
pub fn tie_heavy_records(n: usize, seed: u64) -> (Header, Vec<RawRecord>) {
    /// `(tid, pos, flags, name, mate_tid, mate_pos)` of one record, before its tag.
    type Spec = (i32, i32, u16, Vec<u8>, i32, i32);
    const PAIRED: u16 = 0x1;
    const UNMAPPED: u16 = 0x4;
    const MATE_UNMAPPED: u16 = 0x8;
    const REVERSE: u16 = 0x10;
    const MATE_REVERSE: u16 = 0x20;
    const READ1: u16 = 0x40;
    const READ2: u16 = 0x80;
    const SECONDARY: u16 = 0x100;
    const SUPPLEMENTARY: u16 = 0x800;

    let mut rng = fixture_rng(seed);
    let mut specs: Vec<Spec> = Vec::with_capacity(n + n / 4);
    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    while specs.len() < n {
        let name = format!(
            "A00132:53:HFHJ:{}:{}:{}:{}",
            1 + rng.random_range(0..4),
            1101 + rng.random_range(0..1578),
            rng.random_range(0..40_000),
            rng.random_range(0..40_000)
        )
        .into_bytes();
        if chance(&mut rng, 2) {
            for read in [READ1, READ2] {
                let flags = PAIRED | UNMAPPED | MATE_UNMAPPED | read;
                specs.push((-1, -1, flags, name.clone(), -1, -1));
            }
            continue;
        }
        let shared = chance(&mut rng, 5);
        let (t1, p1, r1) = if shared {
            (0, 5_000, false)
        } else {
            (
                rng.random_range(0..N_TEST_REFS as u64) as i32,
                rng.random_range(0..1_000_000) as i32,
                chance(&mut rng, 50),
            )
        };
        let (t2, p2, r2) = (t1, p1 + 100 + rng.random_range(0..400) as i32, !r1);
        let strand = |reverse: bool, mate_reverse: bool| {
            (if reverse { REVERSE } else { 0 }) | if mate_reverse { MATE_REVERSE } else { 0 }
        };
        let first: Spec = (t1, p1, PAIRED | READ1 | strand(r1, r2), name.clone(), t2, p2);
        let second: Spec = (t2, p2, PAIRED | READ2 | strand(r2, r1), name.clone(), t1, p1);
        let copies = if chance(&mut rng, 3) { 2 } else { 1 };
        for copy in 0..copies {
            let name = if copy == 0 { name.clone() } else { [&name[..], b"d"].concat() };
            let (mut a, mut b) = (first.clone(), second.clone());
            a.3.clone_from(&name);
            b.3.clone_from(&name);
            if chance(&mut rng, 30) {
                for _ in 0..=rng.random_range(0..2) {
                    specs.push((a.0, a.1, a.2 | SECONDARY, name.clone(), a.4, a.5));
                }
            }
            if chance(&mut rng, 20) {
                let pos = rng.random_range(0..1_000_000) as i32;
                for _ in 0..=rng.random_range(0..2) {
                    specs.push((b.0, pos, b.2 | SUPPLEMENTARY, name.clone(), b.4, b.5));
                }
            }
            specs.push(a);
            specs.push(b);
        }
    }
    rand::seq::SliceRandom::shuffle(specs.as_mut_slice(), &mut rng);
    let records = specs
        .into_iter()
        .enumerate()
        .map(|(i, (tid, pos, flags, name, mate_tid, mate_pos))| {
            let mut aux = b"XIi".to_vec();
            aux.extend_from_slice(&u32::try_from(i).expect("index fits u32").to_le_bytes());
            RawRecord::from(make_bam_bytes(
                tid,
                pos,
                flags,
                &name,
                &[],
                40,
                mate_tid,
                mate_pos,
                &aux,
            ))
        })
        .collect();
    (synthetic_header(), records)
}

/// [`reference_sort`] in template-coordinate order: the stable sort by the
/// full template key (every lane, `cell_tag` included), whatever lanes a
/// sort's `--key-types` would narrow to. A merge of independently sorted
/// inputs must order on it.
pub fn reference_sort_full_template(
    header: &Header,
    records: &[RawRecord],
    cell_tag: Option<fgumi_raw_bam::SamTag>,
) -> RecordBytes {
    reference_sort(SortOrder::TemplateCoordinate, header, records, cell_tag)
}
