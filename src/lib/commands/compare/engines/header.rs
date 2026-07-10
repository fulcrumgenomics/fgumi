//! `@HD`/`@SQ`/`@RG` header-equality comparison for `fgumi compare bams`.
//!
//! The comparison engines (`content`/`positional`, `keyjoin`, `sort_verify`) already parse
//! both inputs' `noodles::sam::Header`, but historically only ever used it to resolve
//! RNAME/RNEXT reference names for diff rendering — the header itself was never compared,
//! so two BAMs with genuinely different reference dictionaries, read groups, or sort-order
//! metadata could still report `IDENTICAL` as long as their records matched (R2 gap, see the
//! compare-hardening design spec's "Compare the header" cross-cutting fix).
//!
//! [`compare_headers`] closes that gap for the fields where a divergence is a real content
//! difference, while normalizing away the fields that legitimately differ between fgumi and
//! fgbio (or between two independent tool invocations) even on functionally equivalent
//! output:
//!
//! - **`@HD`** — the `SO` (sort order), `GO` (group order), and `SS` (sub-sort) tags are
//!   compared byte-for-byte. This is what surfaces `R2-HDR-01` (fgumi previously wrote the
//!   `SS` sub-sort tag without its `<sort-order>:` prefix, e.g. `SS:natural` instead of
//!   `SS:queryname:natural`): before that writer fix, comparing a fgumi-sorted BAM against
//!   an equivalently-sorted fgbio/samtools BAM would have spuriously flagged `SS` as
//!   different. The `VN` (format version) tag is deliberately not compared here — it is a
//!   SAM-spec version marker, not a content field.
//! - **`@SQ`** — the reference sequence dictionary is compared as an ordered list of every
//!   field: a name, length, *order*, **or any other field** (`M5`, `UR`, `AS`, `SP`, …)
//!   difference is a genuine divergence (two BAMs aligned to different references, or a
//!   mis-sorted/edited dictionary, are not functionally equivalent even if every record
//!   happens to match by chance). In particular a differing `M5` checksum under the same
//!   name and length means the two dictionaries describe *different* reference bases — so,
//!   exactly as `@RG` compares every non-`ID` field, `@SQ` compares every non-`SN`/`LN`
//!   field, not just name and length.
//! - **`@RG`** — read groups are compared by `ID`, with every other field (`SM`, `LB`, `PL`,
//!   etc.) checked for that `ID`. Order between distinct read group `ID`s does not matter
//!   (the SAM spec does not require `@RG` order to be significant, unlike `@SQ`), but a
//!   missing/extra `ID` or a field difference under a shared `ID` is a genuine divergence.
//! - **`@PG` and `@CO` are NOT compared.** Program records (name, version, command line) and
//!   free-text comments are tool-invocation metadata: fgumi and fgbio each stamp their own
//!   `@PG` chain, and even two runs of the *same* tool can differ in `@PG`'s command-line
//!   text (e.g. absolute vs. relative input paths) or `@CO` timestamps without the BAM
//!   content itself differing at all. Diffing them would produce a `DIFFER` on every
//!   cross-tool (and many same-tool) comparison regardless of actual content equivalence, so
//!   they are normalized away entirely rather than compared.

use std::collections::BTreeMap;

use noodles::sam::Header;
use noodles::sam::header::record::value::map::header::tag as header_tag;

/// Compares two BAM headers, returning human-readable diff strings for every genuine
/// divergence (see the module docs for exactly which fields are compared vs. normalized).
///
/// Returns `None` if the two headers agree on every compared field, or `Some(diffs)` with
/// one line per disagreement (order: `@HD`, then `@SQ`, then `@RG`) otherwise. Never
/// returns `Some(vec![])` — an empty diff list is represented as `None`.
#[must_use]
pub fn compare_headers(h1: &Header, h2: &Header) -> Option<Vec<String>> {
    let mut diffs = compare_hd(h1, h2);
    diffs.extend(compare_sq(h1, h2));
    diffs.extend(compare_rg(h1, h2));
    if diffs.is_empty() { None } else { Some(diffs) }
}

/// Fold an optional [`compare_headers`] result into an engine's `header_mismatch` flag and
/// `diff_details`, capping entries at `max_diffs` (via [`super::push_diff`]).
///
/// Shared by every engine (`positional`, `keyjoin`, `sort_verify`) that calls
/// [`compare_headers`] and needs to merge its result into its own outcome/diff-details
/// state identically.
pub(crate) fn fold_header_diffs(
    diffs: Option<Vec<String>>,
    mismatch: &mut bool,
    details: &mut Vec<String>,
    max_diffs: usize,
) {
    if let Some(diffs) = diffs {
        *mismatch = true;
        for diff in diffs {
            super::push_diff(details, max_diffs, || format!("header: {diff}"));
        }
    }
}

/// Fetches a single `@HD` tag's raw bytes, or `None` if the header has no `@HD` line at all
/// or the tag is absent from it.
fn hd_tag_value(header: &Header, tag: [u8; 2]) -> Option<&[u8]> {
    header.header().and_then(|hd| hd.other_fields().get(&tag)).map(AsRef::as_ref)
}

/// Renders an `Option<&[u8]>` `@HD`/`@RG` field value for a diff message.
fn format_opt_bytes(v: Option<&[u8]>) -> String {
    match v {
        Some(b) => format!("{:?}", String::from_utf8_lossy(b)),
        None => "<absent>".to_string(),
    }
}

/// Compares the `@HD` `SO`/`GO`/`SS` tags (see the module docs for why exactly these three,
/// and not `VN`).
fn compare_hd(h1: &Header, h2: &Header) -> Vec<String> {
    let fields: [(&str, [u8; 2]); 3] = [
        ("SO", *header_tag::SORT_ORDER.as_ref()),
        ("GO", *header_tag::GROUP_ORDER.as_ref()),
        ("SS", *header_tag::SUBSORT_ORDER.as_ref()),
    ];
    fields
        .iter()
        .filter_map(|(name, tag)| {
            let v1 = hd_tag_value(h1, *tag);
            let v2 = hd_tag_value(h2, *tag);
            (v1 != v2).then(|| {
                format!("@HD {name}: {} vs {}", format_opt_bytes(v1), format_opt_bytes(v2))
            })
        })
        .collect()
}

/// One reference sequence's compared identity: its name (`SN`), `LN` length, and every
/// other field keyed by 2-byte tag. noodles types only `SN`/`LN` on `@SQ`, so `M5`, `UR`,
/// `AS`, `SP`, `DS`, … all land in `other_fields()` and are captured here.
type SqEntry = (String, usize, BTreeMap<[u8; 2], Vec<u8>>);

/// Renders one [`SqEntry`] as `SN:<name> LN:<len> {M5:…, UR:…}` for a diff message
/// (the `{…}` field block is omitted when the sequence has no non-`SN`/`LN` fields).
fn format_sq_entry((name, length, other): &SqEntry) -> String {
    let base = format!("SN:{name} LN:{length}");
    if other.is_empty() {
        return base;
    }
    let fields: Vec<String> = other
        .iter()
        .map(|(tag, value)| {
            format!(
                "{}{}:{}",
                char::from(tag[0]),
                char::from(tag[1]),
                String::from_utf8_lossy(value)
            )
        })
        .collect();
    format!("{base} {{{}}}", fields.join(", "))
}

/// Compares the `@SQ` reference sequence dictionary as an ordered list of every field —
/// order is significant here (unlike `@RG`), since two dictionaries that agree on every field
/// but disagree on order describe different `tid` numbering. Every non-`SN`/`LN` field
/// (`M5`, `UR`, `AS`, `SP`, …) is compared too, not just name and length: two sequences with
/// the same name and length but a different `M5` checksum describe *different* reference
/// bases and must not match — exactly as `@RG` compares every non-`ID` field.
fn compare_sq(h1: &Header, h2: &Header) -> Vec<String> {
    let seqs = |h: &Header| -> Vec<SqEntry> {
        h.reference_sequences()
            .iter()
            .map(|(name, seq)| {
                let other: BTreeMap<[u8; 2], Vec<u8>> = seq
                    .other_fields()
                    .iter()
                    .map(|(tag, value)| (*tag.as_ref(), value.to_vec()))
                    .collect();
                (name.to_string(), seq.length().get(), other)
            })
            .collect()
    };
    let seqs1 = seqs(h1);
    let seqs2 = seqs(h2);
    if seqs1 == seqs2 {
        Vec::new()
    } else {
        let render = |seqs: &[SqEntry]| -> String {
            format!("[{}]", seqs.iter().map(format_sq_entry).collect::<Vec<_>>().join(", "))
        };
        vec![format!(
            "@SQ reference sequences (name, length, fields, order): {} vs {}",
            render(&seqs1),
            render(&seqs2)
        )]
    }
}

/// One read group's non-`ID` fields, keyed by 2-byte tag, sorted for stable diff rendering.
type RgFields = BTreeMap<[u8; 2], Vec<u8>>;

/// Renders an `Option<&RgFields>` `@RG` field-set for a diff message.
fn format_opt_rg_fields(fields: Option<&RgFields>) -> String {
    match fields {
        None => "<absent>".to_string(),
        Some(f) => {
            let rendered: Vec<String> = f
                .iter()
                .map(|(tag, value)| {
                    format!(
                        "{}{}:{}",
                        char::from(tag[0]),
                        char::from(tag[1]),
                        String::from_utf8_lossy(value)
                    )
                })
                .collect();
            format!("{{{}}}", rendered.join(", "))
        }
    }
}

/// Compares `@RG` read groups by `ID`: every other field (`SM`, `LB`, `PL`, etc.) is checked
/// under each shared `ID`, and an `ID` present in only one header is reported directly.
/// `ID` order between distinct read groups is not compared (only within-`ID` field content).
fn compare_rg(h1: &Header, h2: &Header) -> Vec<String> {
    let rg_map = |h: &Header| -> BTreeMap<String, RgFields> {
        h.read_groups()
            .iter()
            .map(|(id, rg)| {
                let fields: RgFields = rg
                    .other_fields()
                    .iter()
                    .map(|(tag, value)| (*tag.as_ref(), value.to_vec()))
                    .collect();
                (id.to_string(), fields)
            })
            .collect()
    };

    let rgs1 = rg_map(h1);
    let rgs2 = rg_map(h2);

    let mut ids: Vec<&String> = rgs1.keys().chain(rgs2.keys()).collect();
    ids.sort();
    ids.dedup();

    ids.into_iter()
        .filter_map(|id| {
            let f1 = rgs1.get(id);
            let f2 = rgs2.get(id);
            if f1 == f2 {
                return None;
            }
            Some(format!("@RG {id}: {} vs {}", format_opt_rg_fields(f1), format_opt_rg_fields(f2)))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::BString;
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::header::tag as hd_tag;
    use noodles::sam::header::record::value::map::{
        Header as HeaderRecord, ReadGroup, ReferenceSequence,
    };
    use rstest::rstest;
    use std::num::NonZeroUsize;

    /// Builds a minimal header with a single `@SQ` (`chr1`, length `1000`) and no `@HD`,
    /// `@RG`, `@PG`, or `@CO` records.
    fn base_header() -> Header {
        let seq = Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).expect("non-zero"));
        Header::builder().add_reference_sequence(BString::from("chr1"), seq).build()
    }

    /// Builds a header with an `@HD` line carrying the given `SO`/`GO`/`SS` tags (any of
    /// which may be omitted).
    fn header_with_hd(so: Option<&str>, go: Option<&str>, ss: Option<&str>) -> Header {
        let mut hd = Map::<HeaderRecord>::default();
        if let Some(so) = so {
            hd.other_fields_mut().insert(hd_tag::SORT_ORDER, BString::from(so));
        }
        if let Some(go) = go {
            hd.other_fields_mut().insert(hd_tag::GROUP_ORDER, BString::from(go));
        }
        if let Some(ss) = ss {
            hd.other_fields_mut().insert(hd_tag::SUBSORT_ORDER, BString::from(ss));
        }
        let seq = Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).expect("non-zero"));
        Header::builder().set_header(hd).add_reference_sequence(BString::from("chr1"), seq).build()
    }

    /// Builds a header with a single `@SQ` reference sequence.
    fn header_with_sq(name: &str, length: usize) -> Header {
        let seq = Map::<ReferenceSequence>::new(NonZeroUsize::new(length).expect("non-zero"));
        Header::builder().add_reference_sequence(BString::from(name), seq).build()
    }

    /// Builds a header with a single `@SQ` reference sequence carrying an `M5` (MD5
    /// checksum) field in addition to `SN`/`LN`.
    fn header_with_sq_m5(name: &str, length: usize, m5: &str) -> Header {
        use noodles::sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;
        let mut seq = Map::<ReferenceSequence>::new(NonZeroUsize::new(length).expect("non-zero"));
        seq.other_fields_mut().insert(MD5_CHECKSUM, BString::from(m5));
        Header::builder().add_reference_sequence(BString::from(name), seq).build()
    }

    /// Builds a header with two `@SQ` reference sequences in the given order.
    fn header_with_sq_pair(first: (&str, usize), second: (&str, usize)) -> Header {
        let seq1 = Map::<ReferenceSequence>::new(NonZeroUsize::new(first.1).expect("non-zero"));
        let seq2 = Map::<ReferenceSequence>::new(NonZeroUsize::new(second.1).expect("non-zero"));
        Header::builder()
            .add_reference_sequence(BString::from(first.0), seq1)
            .add_reference_sequence(BString::from(second.0), seq2)
            .build()
    }

    /// Builds a header with one `@RG` read group carrying an `SM` tag.
    fn header_with_rg(id: &str, sample: &str) -> Header {
        let mut rg = Map::<ReadGroup>::default();
        rg.other_fields_mut().insert(
            noodles::sam::header::record::value::map::read_group::tag::SAMPLE,
            BString::from(sample),
        );
        Header::builder().add_read_group(BString::from(id), rg).build()
    }

    /// Builds a header with one `@PG` program record carrying the given version.
    fn header_with_pg(id: &str, version: &str) -> Header {
        use noodles::sam::header::record::value::map::Program;
        let mut pg = Map::<Program>::default();
        pg.other_fields_mut().insert(
            noodles::sam::header::record::value::map::program::tag::VERSION,
            BString::from(version),
        );
        Header::builder().add_program(BString::from(id), pg).build()
    }

    #[test]
    fn identical_headers_have_no_diff() {
        let h = base_header();
        assert!(compare_headers(&h, &h).is_none());
    }

    #[test]
    fn identical_but_distinct_header_objects_have_no_diff() {
        assert!(compare_headers(&base_header(), &base_header()).is_none());
    }

    #[rstest]
    #[case::name_differs(("chr1", 1000), ("chr2", 1000))]
    #[case::length_differs(("chr1", 1000), ("chr1", 2000))]
    fn differing_sq_name_or_length_is_a_diff(#[case] a: (&str, usize), #[case] b: (&str, usize)) {
        let h1 = header_with_sq(a.0, a.1);
        let h2 = header_with_sq(b.0, b.1);
        let diffs = compare_headers(&h1, &h2).expect("@SQ divergence must be reported");
        assert!(diffs.iter().any(|d| d.starts_with("@SQ")), "diffs: {diffs:?}");
    }

    #[test]
    fn differing_sq_order_is_a_diff() {
        let h1 = header_with_sq_pair(("chr1", 1000), ("chr2", 2000));
        let h2 = header_with_sq_pair(("chr2", 2000), ("chr1", 1000));
        let diffs = compare_headers(&h1, &h2).expect("@SQ order difference must be reported");
        assert!(diffs.iter().any(|d| d.starts_with("@SQ")), "diffs: {diffs:?}");
    }

    /// A differing `M5` checksum under an identical `SN`/`LN` is a genuine `@SQ` divergence:
    /// the two dictionaries describe different reference bases even though name and length
    /// agree. Regression guard for comparing every `@SQ` field, not just name and length.
    #[test]
    fn differing_sq_m5_is_a_diff() {
        let h1 = header_with_sq_m5("chr1", 1000, "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
        let h2 = header_with_sq_m5("chr1", 1000, "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb");
        let diffs = compare_headers(&h1, &h2).expect("@SQ M5 divergence must be reported");
        assert!(diffs.iter().any(|d| d.starts_with("@SQ")), "diffs: {diffs:?}");
    }

    /// The converse: an identical `M5` (and name/length) is *not* a diff — comparing the
    /// extra fields must not spuriously flag two genuinely-equal dictionaries.
    #[test]
    fn matching_sq_m5_is_not_a_diff() {
        let m5 = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
        let h1 = header_with_sq_m5("chr1", 1000, m5);
        let h2 = header_with_sq_m5("chr1", 1000, m5);
        assert!(compare_headers(&h1, &h2).is_none());
    }

    /// An `M5` present on only one side (same name/length) is still a divergence — a field
    /// added or dropped changes the dictionary's identity.
    #[test]
    fn sq_m5_present_only_on_one_side_is_a_diff() {
        let h1 = header_with_sq_m5("chr1", 1000, "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
        let h2 = header_with_sq("chr1", 1000);
        let diffs = compare_headers(&h1, &h2).expect("one-sided @SQ M5 must be reported");
        assert!(diffs.iter().any(|d| d.starts_with("@SQ")), "diffs: {diffs:?}");
    }

    #[rstest]
    #[case::sort_order(Some("coordinate"), None, None, Some("queryname"), None, None)]
    #[case::group_order(Some("unsorted"), Some("query"), None, Some("unsorted"), None, None)]
    #[case::subsort(
        Some("queryname"),
        None,
        Some("queryname:natural"),
        Some("queryname"),
        None,
        Some("queryname:lexicographical")
    )]
    fn differing_hd_sort_fields_is_a_diff(
        #[case] so1: Option<&str>,
        #[case] go1: Option<&str>,
        #[case] ss1: Option<&str>,
        #[case] so2: Option<&str>,
        #[case] go2: Option<&str>,
        #[case] ss2: Option<&str>,
    ) {
        let h1 = header_with_hd(so1, go1, ss1);
        let h2 = header_with_hd(so2, go2, ss2);
        let diffs = compare_headers(&h1, &h2).expect("@HD divergence must be reported");
        assert!(diffs.iter().any(|d| d.starts_with("@HD")), "diffs: {diffs:?}");
    }

    /// The R2-HDR-01 regression check: once the writer emits the prefixed `SS` form on both
    /// sides (fix applied), two headers that both declare `queryname:natural` must compare
    /// equal — this is the case that would have spuriously reported `DIFFER` before the writer fix
    /// (fgumi's bare `SS:natural` vs. fgbio/samtools' `SS:queryname:natural`).
    #[test]
    fn matching_prefixed_ss_tag_is_not_a_diff() {
        let h1 = header_with_hd(Some("queryname"), None, Some("queryname:natural"));
        let h2 = header_with_hd(Some("queryname"), None, Some("queryname:natural"));
        assert!(compare_headers(&h1, &h2).is_none());
    }

    #[test]
    fn differing_rg_sample_is_a_diff() {
        let h1 = header_with_rg("rg1", "sampleA");
        let h2 = header_with_rg("rg1", "sampleB");
        let diffs = compare_headers(&h1, &h2).expect("@RG SM divergence must be reported");
        assert!(diffs.iter().any(|d| d.starts_with("@RG")), "diffs: {diffs:?}");
    }

    #[test]
    fn rg_present_only_on_one_side_is_a_diff() {
        let h1 = header_with_rg("rg1", "sampleA");
        let h2 = base_header();
        let diffs = compare_headers(&h1, &h2).expect("a one-sided @RG must be reported");
        assert!(diffs.iter().any(|d| d.starts_with("@RG")), "diffs: {diffs:?}");
    }

    #[test]
    fn differing_pg_version_is_not_a_diff() {
        let h1 = header_with_pg("fgumi", "1.0.0");
        let h2 = header_with_pg("fgbio", "2.3.4");
        assert!(compare_headers(&h1, &h2).is_none(), "@PG must be normalized, not compared");
    }

    #[test]
    fn differing_comments_are_not_a_diff() {
        let h1 = Header::builder().add_comment(BString::from("built by fgumi")).build();
        let h2 = Header::builder().add_comment(BString::from("built by fgbio")).build();
        assert!(compare_headers(&h1, &h2).is_none(), "@CO must be normalized, not compared");
    }
}
