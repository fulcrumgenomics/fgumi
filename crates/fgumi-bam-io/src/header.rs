//! Utilities for adding @PG (program) records to SAM headers.
//!
//! This module provides functions for managing @PG records in SAM/BAM headers,
//! including automatic PP (previous program) chaining and ID collision handling.

use anyhow::Result;
use bstr::BString;
use noodles::sam::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::Program;
use noodles::sam::header::record::value::map::program::tag;
use std::collections::{HashMap, HashSet};

/// Get the ID of the last program in the @PG chain (for PP chaining).
///
/// Finds the program that is not referenced by any other program's PP tag,
/// i.e., the "leaf" of the chain. When the header holds several chains, the
/// leaf that appears last in header order is chosen, which is the newest one
/// only when header order is chronological. When every program is referenced
/// (a PP cycle), the last program in header order is returned.
///
/// # Arguments
///
/// * `header` - The SAM header to search
///
/// # Returns
///
/// The ID of the last program in the chain, or `None` if there are no programs.
#[must_use]
pub fn get_last_program_id(header: &Header) -> Option<String> {
    let programs = header.programs();
    let program_map = programs.as_ref();

    // Collect all program IDs that are referenced as PP by other programs
    let mut referenced: HashSet<&[u8]> = HashSet::new();
    for (_id, pg) in program_map {
        if let Some(pp) = pg.other_fields().get(&tag::PREVIOUS_PROGRAM_ID) {
            referenced.insert(pp.as_ref());
        }
    }

    program_map
        .keys()
        .rev()
        .find(|id| !referenced.contains(id.as_slice()))
        .or_else(|| program_map.keys().next_back())
        .map(|id| String::from_utf8_lossy(id).to_string())
}

/// Create a unique program ID by appending .1, .2, etc. if needed.
///
/// # Arguments
///
/// * `header` - The SAM header to check for existing IDs
/// * `base_id` - The base program ID to use (e.g., "fgumi")
///
/// # Returns
///
/// A program ID not already in the header, either the base ID or with a numeric suffix.
#[must_use]
pub fn make_unique_program_id(header: &Header, base_id: &str) -> String {
    let program_map = header.programs();
    let program_map = program_map.as_ref();

    // Check if base ID is available
    if !program_map.contains_key(base_id.as_bytes()) {
        return base_id.to_string();
    }

    suffixed_id(base_id.as_bytes(), |candidate| program_map.contains_key(candidate)).to_string()
}

/// Returns the first `{id}.{n}` (n = 1, 2, ...) for which `is_taken` returns false.
///
/// This is the scheme fgumi uses whenever a header record needs a fresh ID: its own
/// @PG ([`make_unique_program_id`]) and conflicting @RG/@PG records renamed by `zipper`
/// and `merge`. `{id}.{n}` splits back uniquely at its last dot, so distinct `id`s never
/// yield the same candidate. `is_taken` must return false for some `n`, as it does for
/// any finite set of taken IDs.
#[must_use]
pub fn suffixed_id(id: &[u8], is_taken: impl Fn(&[u8]) -> bool) -> BString {
    let mut suffix = 1_usize;
    loop {
        let mut candidate = BString::from(id);
        candidate.extend_from_slice(format!(".{suffix}").as_bytes());
        if !is_taken(&candidate) {
            return candidate;
        }
        suffix += 1;
    }
}

/// Returns `record` with the ID in its `tag` field (e.g. a @PG `PP` or an @RG `PG`)
/// rewritten through `renames`, if that field names a renamed ID.
#[must_use]
pub fn with_renamed_reference<T, S>(
    record: &Map<T>,
    tag: noodles::sam::header::record::value::map::tag::Other<T::StandardTag>,
    renames: &HashMap<BString, BString, S>,
) -> Map<T>
where
    T: noodles::sam::header::record::value::map::Inner,
    Map<T>: Clone,
    S: std::hash::BuildHasher,
{
    let mut record = record.clone();
    let fields = record.other_fields_mut();
    if let Some(renamed) = fields.get(&tag).and_then(|id| renames.get(id)).cloned() {
        fields.insert(tag, renamed);
    }
    record
}

/// Build a @PG record with all standard fields.
///
/// # Arguments
///
/// * `version` - Program version string
/// * `command_line` - Full command line invocation
/// * `previous_program` - Optional ID of previous program for PP chaining
///
/// # Returns
///
/// A `Map<Program>` ready to add to a header.
/// # Errors
///
/// Returns an error if the program record cannot be built.
pub fn build_program_record(
    version: &str,
    command_line: &str,
    previous_program: Option<&str>,
) -> Result<Map<Program>> {
    let mut builder = Map::<Program>::builder()
        .insert(tag::NAME, "fgumi")
        .insert(tag::VERSION, version)
        .insert(tag::COMMAND_LINE, header_safe_value(command_line));

    if let Some(pp) = previous_program {
        builder = builder.insert(tag::PREVIOUS_PROGRAM_ID, pp);
    }

    Ok(builder.build()?)
}

/// Returns `value` with every character a SAM header field cannot hold replaced, so a
/// record built from it can always be written.
///
/// The SAM spec (§1.3) limits header field values to printable ASCII (`[ -~]+`), and
/// noodles refuses to write anything else, so one tab or non-ASCII byte in a command line
/// would otherwise fail the whole output. Tabs, newlines and carriage returns become a
/// space, as `samtools` does for tabs; any other character becomes its `\u{..}` escape,
/// which keeps a non-ASCII path legible.
#[must_use]
pub fn header_safe_value(value: &str) -> String {
    let mut safe = String::with_capacity(value.len());
    for c in value.chars() {
        match c {
            ' '..='~' => safe.push(c),
            '\t' | '\n' | '\r' => safe.push(' '),
            _ => safe.extend(c.escape_unicode()),
        }
    }
    safe
}

/// Add a @PG record to an existing header with automatic PP chaining.
///
/// This function:
/// 1. Finds the last program in the existing @PG chain (see [`get_last_program_id`])
/// 2. Creates a unique ID (appending .1, .2 if "fgumi" exists)
/// 3. Adds exactly one new @PG with PP pointing to the previous program, even when the
///    header holds several program chains
///
/// This intentionally differs from samtools (and noodles' `Programs::add`), which add
/// one @PG per chain leaf: repeated steps would then grow the header exponentially.
/// Existing @PG records are never modified, and a PP cycle or a PP naming a missing
/// program is tolerated rather than rejected.
///
/// # Arguments
///
/// * `header` - The header to modify
/// * `version` - Program version string
/// * `command_line` - Full command line invocation
///
/// # Returns
///
/// The modified header with the new @PG record.
/// # Errors
///
/// Returns an error if the program record cannot be built.
pub fn add_pg_record(mut header: Header, version: &str, command_line: &str) -> Result<Header> {
    let previous_program = get_last_program_id(&header);
    let unique_id = make_unique_program_id(&header, "fgumi");
    let pg_record = build_program_record(version, command_line, previous_program.as_deref())?;

    header.programs_mut().as_mut().insert(BString::from(unique_id), pg_record);

    Ok(header)
}

/// Ensure the header carries an `@HD` line, synthesizing `@HD VN:1.6 SO:unsorted`
/// when it is absent.
///
/// The SAM spec makes `@HD` optional, but fgbio (and essentially every tool)
/// synthesizes `@HD VN:1.6 SO:unsorted` when reading input that lacks one, and
/// downstream readers can choke on its absence. Commands that read an input
/// header and pass it through (e.g. `correct`, `filter`) otherwise propagate a
/// missing `@HD` — producing a BAM whose header starts at `@PG`, which diverges
/// from fgbio. This normalizes that case.
///
/// An existing `@HD` (and its sort order) is left untouched.
///
/// # Arguments
///
/// * `header` - The header to normalize
///
/// # Returns
///
/// The header, guaranteed to carry an `@HD` line.
///
/// # Errors
///
/// Returns an error if the synthesized `@HD` map cannot be built.
pub fn ensure_hd_record(mut header: Header) -> Result<Header> {
    use noodles::sam::header::record::value::map::Header as HeaderMap;
    use noodles::sam::header::record::value::map::header::{Version, tag as header_tag};

    if header.header().is_none() {
        // Set VN:1.6 explicitly (rather than relying on the map's default) so the
        // synthesized `@HD VN:1.6 SO:unsorted` matches fgbio regardless of the
        // noodles default version.
        let hd = Map::<HeaderMap>::builder()
            .set_version(Version::new(1, 6))
            .insert(header_tag::SORT_ORDER, BString::from("unsorted"))
            .build()?;
        *header.header_mut() = Some(hd);
    }

    Ok(header)
}

/// Add a @PG record to a header builder (for commands creating new headers).
///
/// Use this when building a header from scratch (no PP chaining needed).
///
/// # Arguments
///
/// * `builder` - The header builder to modify
/// * `version` - Program version string
/// * `command_line` - Full command line invocation
///
/// # Returns
///
/// The modified header builder.
/// # Errors
///
/// Returns an error if the program record cannot be built.
pub fn add_pg_to_builder(
    builder: noodles::sam::header::Builder,
    version: &str,
    command_line: &str,
) -> Result<noodles::sam::header::Builder> {
    let pg_record = build_program_record(version, command_line, None)?;
    Ok(builder.add_program("fgumi", pg_record))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_last_program_id_empty() {
        let header = Header::default();
        assert_eq!(get_last_program_id(&header), None);
    }

    #[test]
    fn test_get_last_program_id_single() {
        let mut header = Header::default();
        let pg = Map::<Program>::default();
        header
            .programs_mut()
            .add(BString::from("bwa"), pg)
            .expect("adding program to header should succeed");
        assert_eq!(get_last_program_id(&header), Some("bwa".to_string()));
    }

    #[test]
    fn test_get_last_program_id_chained() {
        let mut header = Header::default();

        // Add first program
        let pg1 = Map::<Program>::default();
        header
            .programs_mut()
            .add(BString::from("bwa"), pg1)
            .expect("adding program to header should succeed");

        // Add second program that references the first
        let pg2 = Map::<Program>::builder()
            .insert(tag::PREVIOUS_PROGRAM_ID, "bwa")
            .build()
            .expect("build should succeed");
        header
            .programs_mut()
            .add(BString::from("samtools"), pg2)
            .expect("adding program to header should succeed");

        // The last program should be samtools (not referenced by anyone)
        assert_eq!(get_last_program_id(&header), Some("samtools".to_string()));
    }

    #[test]
    fn test_get_last_program_id_multiple_chains_picks_newest_leaf() {
        assert_eq!(get_last_program_id(&header_with_two_chains()), Some("bwa-mem3".to_string()));

        let header = Header::builder()
            .add_program("samtools", Map::<Program>::default())
            .add_program("bwa-mem3", Map::<Program>::default())
            .add_program("fgumi", program_with_pp("samtools"))
            .build();
        assert_eq!(get_last_program_id(&header), Some("fgumi".to_string()));
    }

    #[test]
    fn test_get_last_program_id_cycle_picks_last_program() {
        assert_eq!(get_last_program_id(&header_with_pp_cycle()), Some("b".to_string()));
    }

    #[rstest::rstest]
    #[case::printable_ascii_unchanged(
        "fgumi merge -o out.bam a.bam",
        "fgumi merge -o out.bam a.bam"
    )]
    #[case::whitespace_controls_become_spaces("a\tb\nc\rd", "a b c d")]
    #[case::other_controls_escaped("a\u{1}b\u{7f}", "a\\u{1}b\\u{7f}")]
    #[case::non_ascii_escaped("caf\u{e9}.bam", "caf\\u{e9}.bam")]
    fn test_header_safe_value(#[case] value: &str, #[case] expected: &str) {
        assert_eq!(header_safe_value(value), expected);
    }

    /// A command line with characters a header cannot hold still yields a writable @PG.
    #[test]
    fn test_build_program_record_with_unprintable_command_line_is_writable() {
        let pg = build_program_record("1.0", "fgumi merge -o caf\u{e9}.bam\ta.bam", None)
            .expect("program record");
        let header = Header::builder().add_program("fgumi", pg).build();
        let mut writer = noodles::sam::io::Writer::new(Vec::new());
        writer.write_header(&header).expect("header with a sanitized CL must be writable");
        let text = String::from_utf8(writer.get_ref().clone()).expect("UTF-8 header");
        assert!(text.contains("CL:fgumi merge -o caf\\u{e9}.bam a.bam"), "{text}");
    }

    #[rstest::rstest]
    #[case::first_suffix_free("A", &[], "A.1")]
    #[case::skips_taken_suffixes("A", &["A.1", "A.2"], "A.3")]
    #[case::ignores_unrelated_ids("A", &["B.1", "A.1.1"], "A.1")]
    #[case::dotted_id("A.1", &["A.1.1"], "A.1.2")]
    fn test_suffixed_id(#[case] id: &str, #[case] taken: &[&str], #[case] expected: &str) {
        let taken: HashSet<&[u8]> = taken.iter().map(|t| t.as_bytes()).collect();
        assert_eq!(suffixed_id(id.as_bytes(), |candidate| taken.contains(candidate)), expected);
    }

    #[rstest::rstest]
    #[case::renamed(Some("bwa"), Some("bwa.1"))]
    #[case::not_renamed(Some("samtools"), Some("samtools"))]
    #[case::no_reference(None, None)]
    fn test_with_renamed_reference(#[case] previous: Option<&str>, #[case] expected: Option<&str>) {
        let mut builder = Map::<Program>::builder().insert(tag::NAME, "fgumi");
        if let Some(pp) = previous {
            builder = builder.insert(tag::PREVIOUS_PROGRAM_ID, pp);
        }
        let pg = builder.build().expect("valid @PG");
        let renames = HashMap::from([(BString::from("bwa"), BString::from("bwa.1"))]);

        let renamed = with_renamed_reference(&pg, tag::PREVIOUS_PROGRAM_ID, &renames);

        let pp = renamed.other_fields().get(&tag::PREVIOUS_PROGRAM_ID).map(ToString::to_string);
        assert_eq!(pp.as_deref(), expected);
        assert_eq!(renamed.other_fields().get(&tag::NAME), pg.other_fields().get(&tag::NAME));
    }

    #[test]
    fn test_make_unique_program_id_no_collision() {
        let header = Header::default();
        assert_eq!(make_unique_program_id(&header, "fgumi"), "fgumi");
    }

    #[test]
    fn test_make_unique_program_id_with_collision() {
        let mut header = Header::default();
        let pg = Map::<Program>::default();
        header
            .programs_mut()
            .add(BString::from("fgumi"), pg)
            .expect("adding program to header should succeed");

        assert_eq!(make_unique_program_id(&header, "fgumi"), "fgumi.1");
    }

    #[test]
    fn test_make_unique_program_id_multiple_collisions() {
        let mut header = Header::default();

        let pg1 = Map::<Program>::default();
        header
            .programs_mut()
            .add(BString::from("fgumi"), pg1)
            .expect("adding program to header should succeed");

        let pg2 = Map::<Program>::default();
        header
            .programs_mut()
            .add(BString::from("fgumi.1"), pg2)
            .expect("adding program to header should succeed");

        assert_eq!(make_unique_program_id(&header, "fgumi"), "fgumi.2");
    }

    #[test]
    fn test_add_pg_record_empty_header() {
        let header = Header::default();
        let result =
            add_pg_record(header, "1.0.0", "fgumi test").expect("add_pg_record should succeed");
        let programs = result.programs();
        assert_eq!(programs.as_ref().len(), 1);
        assert!(programs.as_ref().contains_key(b"fgumi".as_slice()));

        // Verify the program has expected fields
        let pg =
            programs.as_ref().get(b"fgumi".as_slice()).expect("expected key should be present");
        assert_eq!(
            pg.other_fields().get(&tag::NAME).map(std::convert::AsRef::as_ref),
            Some(b"fgumi".as_slice())
        );
        assert_eq!(
            pg.other_fields().get(&tag::VERSION).map(std::convert::AsRef::as_ref),
            Some(b"1.0.0".as_slice())
        );
        assert_eq!(
            pg.other_fields().get(&tag::COMMAND_LINE).map(std::convert::AsRef::as_ref),
            Some(b"fgumi test".as_slice())
        );
        assert!(pg.other_fields().get(&tag::PREVIOUS_PROGRAM_ID).is_none());
    }

    #[test]
    fn test_add_pg_record_with_existing_fgumi() {
        let mut header = Header::default();
        let pg = Map::<Program>::default();
        header
            .programs_mut()
            .add(BString::from("fgumi"), pg)
            .expect("adding program to header should succeed");

        let result =
            add_pg_record(header, "1.0.0", "fgumi test2").expect("add_pg_record should succeed");
        let programs = result.programs();
        assert_eq!(programs.as_ref().len(), 2);
        assert!(programs.as_ref().contains_key(b"fgumi.1".as_slice()));

        // Verify PP chaining
        let pg =
            programs.as_ref().get(b"fgumi.1".as_slice()).expect("expected key should be present");
        assert_eq!(
            pg.other_fields().get(&tag::PREVIOUS_PROGRAM_ID).map(std::convert::AsRef::as_ref),
            Some(b"fgumi".as_slice())
        );
    }

    #[test]
    fn test_add_pg_record_chains_to_non_fgumi() {
        let mut header = Header::default();

        // Add a BWA program first
        let bwa_pg = Map::<Program>::builder()
            .insert(tag::NAME, "bwa")
            .insert(tag::VERSION, "0.7.17")
            .build()
            .expect("building program map should succeed");
        header
            .programs_mut()
            .add(BString::from("bwa"), bwa_pg)
            .expect("adding program to header should succeed");

        let result = add_pg_record(header, "1.0.0", "fgumi group -i in.bam")
            .expect("add_pg_record should succeed");
        let programs = result.programs();

        // fgumi should chain to bwa
        let pg =
            programs.as_ref().get(b"fgumi".as_slice()).expect("expected key should be present");
        assert_eq!(
            pg.other_fields().get(&tag::PREVIOUS_PROGRAM_ID).map(std::convert::AsRef::as_ref),
            Some(b"bwa".as_slice())
        );
    }

    fn program_with_pp(pp: &str) -> Map<Program> {
        Map::<Program>::builder()
            .insert(tag::PREVIOUS_PROGRAM_ID, pp)
            .build()
            .expect("building program map should succeed")
    }

    fn previous_program(header: &Header, id: &str) -> Option<String> {
        header
            .programs()
            .as_ref()
            .get(id.as_bytes())
            .and_then(|pg| pg.other_fields().get(&tag::PREVIOUS_PROGRAM_ID))
            .map(ToString::to_string)
    }

    /// Two root programs, as in `zipper`'s merge of an unmapped and a mapped BAM header.
    fn header_with_two_chains() -> Header {
        Header::builder()
            .add_program("samtools", Map::<Program>::default())
            .add_program("bwa-mem3", Map::<Program>::default())
            .build()
    }

    fn header_with_pp_cycle() -> Header {
        Header::builder()
            .add_program("a", program_with_pp("b"))
            .add_program("b", program_with_pp("a"))
            .build()
    }

    #[test]
    fn test_add_pg_record_single_chain_chains_to_leaf() {
        let header = Header::builder()
            .add_program("bwa", Map::<Program>::default())
            .add_program("samtools", program_with_pp("bwa"))
            .build();

        let result =
            add_pg_record(header, "1.0.0", "fgumi sort").expect("add_pg_record should succeed");

        assert_eq!(result.programs().as_ref().len(), 3);
        assert_eq!(previous_program(&result, "fgumi").as_deref(), Some("samtools"));
    }

    #[test]
    fn test_add_pg_record_adds_one_record_with_multiple_chains() {
        let result = add_pg_record(header_with_two_chains(), "1.0.0", "fgumi zipper")
            .expect("add_pg_record should succeed");

        let ids: Vec<String> = result.programs().as_ref().keys().map(ToString::to_string).collect();
        assert_eq!(ids, ["samtools", "bwa-mem3", "fgumi"]);
        assert_eq!(previous_program(&result, "fgumi").as_deref(), Some("bwa-mem3"));
    }

    #[rstest::rstest]
    #[case::pp_cycle(header_with_pp_cycle(), "b")]
    #[case::dangling_pp(Header::builder().add_program("bwa", program_with_pp("missing")).build(), "bwa")]
    fn test_add_pg_record_tolerates_malformed_chain(
        #[case] header: Header,
        #[case] expected_pp: &str,
    ) {
        let before = header.programs().clone();

        let result =
            add_pg_record(header, "1.0.0", "fgumi sort").expect("add_pg_record should succeed");

        let after = result.programs().as_ref();
        assert_eq!(after.len(), before.as_ref().len() + 1);
        assert!(before.as_ref().iter().all(|(id, pg)| after.get(id) == Some(pg)));
        assert_eq!(previous_program(&result, "fgumi").as_deref(), Some(expected_pp));
    }

    #[test]
    fn test_add_pg_record_never_overwrites_existing_id() {
        let mut builder = Header::builder().add_program("fgumi", Map::<Program>::default());
        for i in 1..=1000 {
            builder = builder.add_program(format!("fgumi.{i}"), Map::<Program>::default());
        }
        let header = builder.build();
        let before = header.programs().clone();

        let result =
            add_pg_record(header, "1.0.0", "fgumi sort").expect("add_pg_record should succeed");

        let after = result.programs().as_ref();
        assert_eq!(after.len(), 1002);
        assert!(before.as_ref().iter().all(|(id, pg)| after.get(id) == Some(pg)));
        assert_eq!(previous_program(&result, "fgumi.1001").as_deref(), Some("fgumi.1000"));
    }

    #[test]
    fn test_add_pg_record_repeated_adds_grow_by_one() {
        let mut header = header_with_two_chains();
        let mut expected_pp = "bwa-mem3".to_string();

        for i in 0..6 {
            let expected_id = if i == 0 { "fgumi".to_string() } else { format!("fgumi.{i}") };
            let before = header.programs().as_ref().len();

            header =
                add_pg_record(header, "1.0.0", "fgumi sort").expect("add_pg_record should succeed");

            let programs = header.programs();
            assert_eq!(programs.as_ref().len(), before + 1, "each add must insert one @PG");
            let (last_id, _) = programs.as_ref().last().expect("header has programs");
            assert_eq!(last_id.to_string(), expected_id);
            assert_eq!(previous_program(&header, &expected_id), Some(expected_pp));
            expected_pp = expected_id;
        }

        assert_eq!(header.programs().as_ref().len(), 8);
        assert_eq!(header.programs().leaves().expect("no cycles").count(), 2);
    }

    #[test]
    fn test_add_pg_to_builder() {
        let builder = Header::builder();
        let builder = add_pg_to_builder(builder, "1.0.0", "fgumi extract")
            .expect("add_pg_to_builder should succeed");
        let header = builder.build();

        let programs = header.programs();
        assert_eq!(programs.as_ref().len(), 1);

        let pg =
            programs.as_ref().get(b"fgumi".as_slice()).expect("expected key should be present");
        assert_eq!(
            pg.other_fields().get(&tag::NAME).map(std::convert::AsRef::as_ref),
            Some(b"fgumi".as_slice())
        );
        assert!(pg.other_fields().get(&tag::PREVIOUS_PROGRAM_ID).is_none());
    }

    #[test]
    fn test_add_pg_record_empty_command_line() {
        let header = Header::default();
        let result = add_pg_record(header, "1.0.0", "").expect("add_pg_record should succeed");
        let programs = result.programs();
        assert_eq!(programs.as_ref().len(), 1);
        assert!(programs.as_ref().contains_key(b"fgumi".as_slice()));
    }

    #[test]
    fn test_add_pg_record_write_to_bam() {
        use crate::writer::create_bam_writer;
        use tempfile::TempDir;

        let dir = TempDir::new().expect("creating temp file/dir should succeed");
        let output_path = dir.path().join("test.bam");

        let header = Header::default();
        let result =
            add_pg_record(header, "1.0.0", "fgumi test").expect("add_pg_record should succeed");

        // Try to write the header to a BAM file
        let _writer = create_bam_writer(&output_path, &result, 1, 6)
            .expect("creating BAM writer should succeed");
    }

    #[test]
    fn test_add_pg_record_chains_to_empty_program() {
        use crate::writer::create_bam_writer;
        use tempfile::TempDir;

        // Simulate what SamBuilder does - adds an empty/default program
        let pg_map = Map::<Program>::default();
        let header = Header::builder().add_program("SamBuilder", pg_map).build();

        // Now add our fgumi @PG record
        let result =
            add_pg_record(header, "1.0.0", "fgumi test").expect("add_pg_record should succeed");
        let programs = result.programs();
        assert_eq!(programs.as_ref().len(), 2);

        // fgumi should chain to SamBuilder
        let pg =
            programs.as_ref().get(b"fgumi".as_slice()).expect("expected key should be present");
        assert_eq!(
            pg.other_fields().get(&tag::PREVIOUS_PROGRAM_ID).map(std::convert::AsRef::as_ref),
            Some(b"SamBuilder".as_slice())
        );

        // Try to write to BAM
        let dir = TempDir::new().expect("creating temp file/dir should succeed");
        let output_path = dir.path().join("test.bam");
        let _writer = create_bam_writer(&output_path, &result, 1, 6)
            .expect("creating BAM writer should succeed");
    }

    #[test]
    fn test_ensure_hd_record_synthesizes_when_absent() {
        use noodles::sam::header::record::value::map::header::{Version, tag as header_tag};

        // A header built without an explicit @HD (only @SQ etc.) reports None.
        let header = Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<noodles::sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZero::new(1000).expect("non-zero reference length"),
                ),
            )
            .build();
        assert!(header.header().is_none(), "precondition: no @HD");

        let header = ensure_hd_record(header).expect("ensure_hd_record should succeed");
        let hd = header.header().expect("@HD must now be present");
        assert_eq!(hd.version(), Version::new(1, 6), "@HD VN must be 1.6");
        assert_eq!(
            hd.other_fields().get(&header_tag::SORT_ORDER).map(|v| v.to_vec()),
            Some(b"unsorted".to_vec()),
            "@HD SO must be 'unsorted' to match fgbio"
        );
    }

    #[test]
    fn test_ensure_hd_record_preserves_existing() {
        use noodles::sam::header::record::value::map::Header as HeaderMap;
        use noodles::sam::header::record::value::map::header::{Version, tag as header_tag};

        // Existing @HD with a NON-default version and SO:coordinate must be left
        // entirely untouched (not just the SO tag).
        let existing = Map::<HeaderMap>::builder()
            .set_version(Version::new(1, 5))
            .insert(header_tag::SORT_ORDER, BString::from("coordinate"))
            .build()
            .expect("building @HD map should succeed");
        let header = Header::builder().set_header(existing.clone()).build();

        let header = ensure_hd_record(header).expect("ensure_hd_record should succeed");
        let hd = header.header().expect("@HD present");
        assert_eq!(*hd, existing, "an existing @HD (version + all fields) must be left untouched");
    }
}
