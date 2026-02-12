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
use std::collections::HashSet;

/// Get the ID of the last program in the @PG chain (for PP chaining).
///
/// Finds the program that is not referenced by any other program's PP tag,
/// i.e., the "leaf" of the chain.
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

    if program_map.is_empty() {
        return None;
    }

    // Collect all program IDs that are referenced as PP by other programs
    let mut referenced: HashSet<&[u8]> = HashSet::new();
    for (_id, pg) in program_map {
        if let Some(pp) = pg.other_fields().get(&tag::PREVIOUS_PROGRAM_ID) {
            referenced.insert(pp.as_ref());
        }
    }

    // Find a program that is NOT referenced (the leaf/end of chain)
    for (id, _pg) in program_map {
        if !referenced.contains(id.as_slice()) {
            return Some(String::from_utf8_lossy(id).to_string());
        }
    }

    // Fallback: return any program ID (shouldn't happen with valid headers)
    program_map.keys().next().map(|id| String::from_utf8_lossy(id).to_string())
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
/// A unique program ID, either the base ID or with a numeric suffix.
#[must_use]
pub fn make_unique_program_id(header: &Header, base_id: &str) -> String {
    let programs = header.programs();
    let program_map = programs.as_ref();

    // Check if base ID is available
    if !program_map.contains_key(base_id.as_bytes()) {
        return base_id.to_string();
    }

    // Append numeric suffix until unique
    for i in 1..=1000 {
        let candidate = format!("{base_id}.{i}");
        if !program_map.contains_key(candidate.as_bytes()) {
            return candidate;
        }
    }

    // Extremely unlikely fallback
    format!("{base_id}.{}", std::process::id())
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
        .insert(tag::COMMAND_LINE, command_line);

    if let Some(pp) = previous_program {
        builder = builder.insert(tag::PREVIOUS_PROGRAM_ID, pp);
    }

    Ok(builder.build()?)
}

/// Add a @PG record to an existing header with automatic PP chaining.
///
/// This function:
/// 1. Finds the last program in the existing @PG chain
/// 2. Creates a unique ID (appending .1, .2 if "fgumi" exists)
/// 3. Adds the new @PG with PP pointing to the previous program
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
/// Returns an error if the program record cannot be added to the header.
pub fn add_pg_record(mut header: Header, version: &str, command_line: &str) -> Result<Header> {
    let previous_program = get_last_program_id(&header);
    let unique_id = make_unique_program_id(&header, "fgumi");
    let pg_record = build_program_record(version, command_line, previous_program.as_deref())?;

    header.programs_mut().add(BString::from(unique_id), pg_record)?;

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
        header.programs_mut().add(BString::from("bwa"), pg).unwrap();
        assert_eq!(get_last_program_id(&header), Some("bwa".to_string()));
    }

    #[test]
    fn test_get_last_program_id_chained() {
        let mut header = Header::default();

        // Add first program
        let pg1 = Map::<Program>::default();
        header.programs_mut().add(BString::from("bwa"), pg1).unwrap();

        // Add second program that references the first
        let pg2 =
            Map::<Program>::builder().insert(tag::PREVIOUS_PROGRAM_ID, "bwa").build().unwrap();
        header.programs_mut().add(BString::from("samtools"), pg2).unwrap();

        // The last program should be samtools (not referenced by anyone)
        assert_eq!(get_last_program_id(&header), Some("samtools".to_string()));
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
        header.programs_mut().add(BString::from("fgumi"), pg).unwrap();

        assert_eq!(make_unique_program_id(&header, "fgumi"), "fgumi.1");
    }

    #[test]
    fn test_make_unique_program_id_multiple_collisions() {
        let mut header = Header::default();

        let pg1 = Map::<Program>::default();
        header.programs_mut().add(BString::from("fgumi"), pg1).unwrap();

        let pg2 = Map::<Program>::default();
        header.programs_mut().add(BString::from("fgumi.1"), pg2).unwrap();

        assert_eq!(make_unique_program_id(&header, "fgumi"), "fgumi.2");
    }

    #[test]
    fn test_add_pg_record_empty_header() {
        let header = Header::default();
        let result = add_pg_record(header, "1.0.0", "fgumi test").unwrap();
        let programs = result.programs();
        assert_eq!(programs.as_ref().len(), 1);
        assert!(programs.as_ref().contains_key(b"fgumi".as_slice()));

        // Verify the program has expected fields
        let pg = programs.as_ref().get(b"fgumi".as_slice()).unwrap();
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
        header.programs_mut().add(BString::from("fgumi"), pg).unwrap();

        let result = add_pg_record(header, "1.0.0", "fgumi test2").unwrap();
        let programs = result.programs();
        assert_eq!(programs.as_ref().len(), 2);
        assert!(programs.as_ref().contains_key(b"fgumi.1".as_slice()));

        // Verify PP chaining
        let pg = programs.as_ref().get(b"fgumi.1".as_slice()).unwrap();
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
            .unwrap();
        header.programs_mut().add(BString::from("bwa"), bwa_pg).unwrap();

        let result = add_pg_record(header, "1.0.0", "fgumi group -i in.bam").unwrap();
        let programs = result.programs();

        // fgumi should chain to bwa
        let pg = programs.as_ref().get(b"fgumi".as_slice()).unwrap();
        assert_eq!(
            pg.other_fields().get(&tag::PREVIOUS_PROGRAM_ID).map(std::convert::AsRef::as_ref),
            Some(b"bwa".as_slice())
        );
    }

    #[test]
    fn test_add_pg_to_builder() {
        let builder = Header::builder();
        let builder = add_pg_to_builder(builder, "1.0.0", "fgumi extract").unwrap();
        let header = builder.build();

        let programs = header.programs();
        assert_eq!(programs.as_ref().len(), 1);

        let pg = programs.as_ref().get(b"fgumi".as_slice()).unwrap();
        assert_eq!(
            pg.other_fields().get(&tag::NAME).map(std::convert::AsRef::as_ref),
            Some(b"fgumi".as_slice())
        );
        assert!(pg.other_fields().get(&tag::PREVIOUS_PROGRAM_ID).is_none());
    }

    #[test]
    fn test_add_pg_record_empty_command_line() {
        let header = Header::default();
        let result = add_pg_record(header, "1.0.0", "").unwrap();
        let programs = result.programs();
        assert_eq!(programs.as_ref().len(), 1);
        assert!(programs.as_ref().contains_key(b"fgumi".as_slice()));
    }

    #[test]
    fn test_add_pg_record_write_to_bam() {
        use crate::bam_io::create_bam_writer;
        use tempfile::TempDir;

        let dir = TempDir::new().unwrap();
        let output_path = dir.path().join("test.bam");

        let header = Header::default();
        let result = add_pg_record(header, "1.0.0", "fgumi test").unwrap();

        // Try to write the header to a BAM file
        let _writer = create_bam_writer(&output_path, &result, 1, 6).unwrap();
    }

    #[test]
    fn test_add_pg_record_chains_to_empty_program() {
        use crate::bam_io::create_bam_writer;
        use tempfile::TempDir;

        // Simulate what SamBuilder does - adds an empty/default program
        let pg_map = Map::<Program>::default();
        let header = Header::builder().add_program("SamBuilder", pg_map).build();

        // Now add our fgumi @PG record
        let result = add_pg_record(header, "1.0.0", "fgumi test").unwrap();
        let programs = result.programs();
        assert_eq!(programs.as_ref().len(), 2);

        // fgumi should chain to SamBuilder
        let pg = programs.as_ref().get(b"fgumi".as_slice()).unwrap();
        assert_eq!(
            pg.other_fields().get(&tag::PREVIOUS_PROGRAM_ID).map(std::convert::AsRef::as_ref),
            Some(b"SamBuilder".as_slice())
        );

        // Try to write to BAM
        let dir = TempDir::new().unwrap();
        let output_path = dir.path().join("test.bam");
        let _writer = create_bam_writer(&output_path, &result, 1, 6).unwrap();
    }
}
