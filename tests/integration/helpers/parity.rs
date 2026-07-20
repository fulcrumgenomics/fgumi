//! BAM-equivalence assertion for runall parity tests.
//!
//! Two BAM files are "parity-equivalent" when they contain the same
//! sequence of alignment records (`RecordBuf` byte-equal). Header
//! differences are tolerated for fields that legitimately diverge
//! between invocations (`@PG` carries the literal command-line string;
//! `@HD` may carry slight SO/GO/SS variation depending on which
//! pipeline produced the file).
//!
//! For the parity tests' purposes, the RECORDS are what defines the
//! contract: same input + same parameters must yield the same record
//! stream regardless of CLI entry point.

#![allow(dead_code)]

use std::fs;
use std::path::Path;

use noodles::bam;
use noodles::sam::alignment::record_buf::RecordBuf;

/// Read all records from a BAM file into a `Vec<RecordBuf>`, in file
/// order. Used by `assert_bams_record_equivalent` and by tests that
/// need to inspect record streams directly.
pub fn read_bam_records(path: &Path) -> Vec<RecordBuf> {
    let mut reader = bam::io::Reader::new(
        fs::File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display())),
    );
    let header =
        reader.read_header().unwrap_or_else(|e| panic!("read header from {}: {e}", path.display()));
    reader
        .record_bufs(&header)
        .collect::<std::io::Result<Vec<_>>>()
        .unwrap_or_else(|e| panic!("read records from {}: {e}", path.display()))
}

/// Assert that two BAMs contain the same record stream (count + each
/// record byte-equal as `RecordBuf`).
///
/// Header content is intentionally NOT compared — `@PG` lines record
/// the literal command-line string, which differs between runall and
/// the equivalent standalone/staged chain. The parity contract is
/// about the record stream, not the metadata.
///
/// # Panics
///
/// Panics with a diagnostic on the first mismatch (record count or
/// first divergent record index).
pub fn assert_bams_record_equivalent(a: &Path, b: &Path) {
    let recs_a = read_bam_records(a);
    let recs_b = read_bam_records(b);
    assert_eq!(
        recs_a.len(),
        recs_b.len(),
        "record count differs: {} has {}, {} has {}",
        a.display(),
        recs_a.len(),
        b.display(),
        recs_b.len(),
    );
    for (i, (ra, rb)) in recs_a.iter().zip(recs_b.iter()).enumerate() {
        assert_eq!(ra, rb, "record {i} differs between {} and {}", a.display(), b.display());
    }
}

/// Like [`assert_bams_record_equivalent`], but additionally asserts the
/// record stream is NON-empty.
///
/// `assert_bams_record_equivalent` is vacuously satisfied when both BAMs
/// emit zero records (`0 == 0`, the per-record loop never runs), so a
/// parity test whose chain is *expected* to retain records could stay
/// green even if the fused pipeline silently emitted a header-only BAM
/// (S9a-001). Use this variant for every parity chain that must carry
/// records; reserve the bare `assert_bams_record_equivalent` (paired with
/// an explicit `assert_eq!(read_bam_records(..).len(), 0, ..)`) only for
/// chains that legitimately produce zero records by design.
///
/// # Panics
///
/// Panics if the BAMs are not record-equivalent, or if `a` (and therefore
/// `b`, since the counts were just asserted equal) contains zero records.
pub fn assert_bams_record_equivalent_nonempty(a: &Path, b: &Path) {
    assert_bams_record_equivalent(a, b);
    let n = read_bam_records(a).len();
    assert!(
        n > 0,
        "expected a non-empty record stream but {} has 0 records — \
         a header-only BAM would satisfy the equivalence vacuously",
        a.display(),
    );
}

/// Assert the two BAM headers agree on the correctness-relevant fields,
/// ignoring only `@PG` (which records the literal command line and so
/// legitimately diverges between a fused runall and the equivalent staged
/// chain).
///
/// Compares (S9a-005):
///   * `@HD` sort order / grouping / sub-sort (`SO`/`GO`/`SS`),
///   * the `@SQ` reference-sequence dictionary (names + lengths, in order),
///   * the `@RG` read-group records (id + every field, e.g. `LB`/`SM`, in order).
///
/// A runall regression that emitted `SO:unsorted` after a sort stage,
/// dropped an `@SQ`, or failed to collapse read groups would be invisible
/// to `assert_bams_record_equivalent` (which ignores the whole header);
/// this closes that gap without re-asserting the `@PG` command-line text.
///
/// # Panics
///
/// Panics on the first divergent field.
/// One `@RG` record reduced to `(id, sorted (tag, value) fields)` — the
/// canonical form compared by [`assert_bam_headers_equivalent_ignoring_pg`].
type ReadGroupRecord = (Vec<u8>, Vec<(Vec<u8>, bstr::BString)>);

/// Project a header's `@RG` lines to comparable [`ReadGroupRecord`]s (id plus
/// every field, e.g. `LB`/`SM`), preserving file order.
fn read_group_records(h: &noodles::sam::Header) -> Vec<ReadGroupRecord> {
    h.read_groups()
        .iter()
        .map(|(id, map)| {
            let mut fields: Vec<(Vec<u8>, bstr::BString)> = map
                .other_fields()
                .iter()
                .map(|(tag, value)| (tag.as_ref().to_vec(), value.clone()))
                .collect();
            fields.sort();
            (id.to_vec(), fields)
        })
        .collect()
}

pub fn assert_bam_headers_equivalent_ignoring_pg(a: &Path, b: &Path) {
    fn read_header(path: &Path) -> noodles::sam::Header {
        let mut reader = bam::io::Reader::new(
            fs::File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display())),
        );
        reader.read_header().unwrap_or_else(|e| panic!("read header from {}: {e}", path.display()))
    }

    use noodles::sam::header::record::value::map::header::tag;
    use noodles::sam::header::record::value::map::header::tag::Standard;
    use noodles::sam::header::record::value::map::tag::Other;

    let ha = read_header(a);
    let hb = read_header(b);

    // @HD — SO/GO/SS sort-order fields (read from `other_fields`, the only
    // place noodles exposes them). VN is intentionally NOT compared (it can
    // differ harmlessly across writers).
    let hd_field = |h: &noodles::sam::Header, t: Other<Standard>| -> Option<bstr::BString> {
        h.header().and_then(|m| m.other_fields().get(&t).cloned())
    };
    for (t, name) in [(tag::SORT_ORDER, "SO"), (tag::GROUP_ORDER, "GO"), (tag::SUBSORT_ORDER, "SS")]
    {
        assert_eq!(
            hd_field(&ha, t),
            hd_field(&hb, t),
            "@HD {name} differs between {} and {}",
            a.display(),
            b.display(),
        );
    }

    // @SQ — reference dictionary (names + lengths, in file order).
    let sq = |h: &noodles::sam::Header| -> Vec<(Vec<u8>, usize)> {
        h.reference_sequences()
            .iter()
            .map(|(name, map)| (name.to_vec(), map.length().get()))
            .collect()
    };
    assert_eq!(
        sq(&ha),
        sq(&hb),
        "@SQ dictionary differs between {} and {}",
        a.display(),
        b.display(),
    );

    // @RG — full read-group records (id plus every field, e.g. LB/SM/PL), in
    // file order. Comparing only ids would bless a fused path that altered an
    // @RG field that later stages read — `LibraryIndex::from_header` consumes
    // @RG LB — so a changed library would diverge downstream yet pass here.
    assert_eq!(
        read_group_records(&ha),
        read_group_records(&hb),
        "@RG read-group records differ between {} and {}",
        a.display(),
        b.display(),
    );
}
