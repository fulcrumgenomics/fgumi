//! Merge pre-sorted BAM files into a single sorted BAM.
//!
//! Performs a k-way merge of BAM files that are already sorted in the same
//! order, producing a single merged output that preserves the sort order.
//!
//! Similar to `samtools merge`, but supports template-coordinate order and
//! uses the same high-performance merge infrastructure as `fgumi sort`.

use std::collections::HashSet;
use std::path::PathBuf;

use crate::logging::OperationTimer;
use crate::validation::validate_file_exists;
use anyhow::{Result, bail};
use clap::Parser;
#[cfg(test)]
use fgumi_bam_io::create_bam_reader;
use fgumi_bam_io::create_raw_bam_reader;
use fgumi_sort::RawExternalSorter;
use log::info;
use noodles::sam::Header;

use crate::commands::command::Command;
use crate::commands::sort::SortOrderArg;

/// Merge pre-sorted BAM files.
///
/// Performs a k-way merge of multiple BAM files that are already sorted in
/// the same order, similar to `samtools merge`. Input files must all be
/// sorted in the specified order.
#[derive(Debug, Parser)]
#[command(
    name = "merge",
    about = "\x1b[38;5;72m[ALIGNMENT]\x1b[0m      \x1b[36mMerge pre-sorted BAM files into a single sorted BAM\x1b[0m",
    long_about = r#"
Merge pre-sorted BAM files into a single sorted BAM.

Performs a k-way merge of multiple BAM files that are already sorted in the
same order, producing a single merged output that preserves the sort order.
Similar to `samtools merge`, but supports template-coordinate order.

Input files must all be sorted in the specified sort order.

EXAMPLES:

  # Merge coordinate-sorted BAMs
  fgumi merge -o merged.bam sorted1.bam sorted2.bam sorted3.bam

  # Merge template-coordinate sorted BAMs
  fgumi merge -o merged.bam --order template-coordinate tc1.bam tc2.bam

  # Merge from a file listing input BAMs (one per line)
  fgumi merge -o merged.bam -b input_list.txt --order queryname

  # Merge with multiple threads
  fgumi merge -o merged.bam -@ 4 sorted1.bam sorted2.bam

"#
)]
pub struct Merge {
    /// Output BAM file.
    #[arg(short = 'o', long = "output")]
    pub output: PathBuf,

    /// Input BAM files to merge (positional).
    #[arg(required_unless_present = "input_list")]
    pub inputs: Vec<PathBuf>,

    /// File containing a list of input BAM paths, one per line.
    ///
    /// Can be combined with positional inputs.
    #[arg(short = 'b', long = "input-list")]
    pub input_list: Option<PathBuf>,

    /// Sort order of the input files.
    #[arg(long = "order", default_value = "template-coordinate", value_parser = SortOrderArg::parse)]
    pub order: SortOrderArg,

    /// Number of threads for parallel operations.
    ///
    /// Used for multi-threaded BGZF compression.
    #[arg(short = '@', short_alias = 't', long = "threads", default_value = "1")]
    pub threads: usize,

    /// Compression level for output BAM (1-12).
    ///
    /// Level 1 is fastest with larger files.
    /// Level 6 (default) balances speed and file size.
    /// Level 12 produces smallest files but is slowest.
    #[arg(long = "compression-level", default_value_t = 6)]
    pub compression_level: u32,
}

impl Command for Merge {
    fn execute(&self, _command_line: &str) -> Result<()> {
        let mut input_paths: Vec<PathBuf> = self.inputs.clone();

        if let Some(ref list_path) = self.input_list {
            validate_file_exists(list_path, "Input list")?;
            let contents = std::fs::read_to_string(list_path)?;
            for line in contents.lines() {
                let line = line.trim();
                if !line.is_empty() && !line.starts_with('#') {
                    input_paths.push(PathBuf::from(line));
                }
            }
        }

        if input_paths.is_empty() {
            bail!("No input files specified");
        }

        for path in &input_paths {
            validate_file_exists(path, "Input BAM")?;
        }

        // Check output doesn't alias any input
        if let Ok(output_canon) = std::fs::canonicalize(&self.output) {
            for path in &input_paths {
                if let Ok(input_canon) = std::fs::canonicalize(path) {
                    if output_canon == input_canon {
                        bail!(
                            "Output file '{}' is the same as input file '{}'",
                            self.output.display(),
                            path.display()
                        );
                    }
                }
            }
        }

        let cell_tag = crate::commands::sort::parse_cell_tag(self.order)?;

        let timer = OperationTimer::new("Merging BAMs");

        info!("Starting Merge");
        info!("Inputs: {} files", input_paths.len());
        for path in &input_paths {
            info!("  {}", path.display());
        }
        info!("Output: {}", self.output.display());
        info!("Sort order: {:?}", self.order);
        if let Some(ct) = cell_tag {
            let ct_bytes = *ct;
            info!("Cell tag: {}{}", ct_bytes[0] as char, ct_bytes[1] as char);
        }
        info!("Threads: {}", self.threads);

        // Read and merge headers from all inputs (also validates each input's
        // declared sort order against --order; MERGE3-01 fast check).
        let header = merge_headers(&input_paths, self.order)?;

        let mut sorter = RawExternalSorter::new(self.order.into())
            .threads(self.threads)
            .output_compression(self.compression_level);

        if let Some(ct) = cell_tag {
            sorter = sorter.cell_tag(ct);
        }

        let records_merged = sorter.merge_bams(&input_paths, &header, &self.output)?;

        info!("=== Summary ===");
        info!("Records merged: {records_merged}");
        info!("Output: {}", self.output.display());

        timer.log_completion(records_merged);
        Ok(())
    }
}

/// A BAM header's *declared* sort order, as far as merge validation cares.
///
/// `classify_declared_order` returns `None` when the header asserts no usable
/// order (`SO` absent, or `SO:unsorted` without the template-coordinate
/// sub-sort). Those inputs pass the fast header check and are instead verified
/// record-by-record during the merge (the streaming monotonicity check).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DeclaredOrder {
    Coordinate,
    Queryname,
    TemplateCoordinate,
}

impl DeclaredOrder {
    fn as_str(self) -> &'static str {
        match self {
            DeclaredOrder::Coordinate => "coordinate",
            DeclaredOrder::Queryname => "queryname",
            DeclaredOrder::TemplateCoordinate => "template-coordinate",
        }
    }
}

/// Classifies a header's declared sort order for merge-input validation.
fn classify_declared_order(header: &Header) -> Option<DeclaredOrder> {
    use noodles::sam::header::record::value::map::header::sort_order::{COORDINATE, QUERY_NAME};
    if fgumi_sam::is_template_coordinate_sorted(header) {
        Some(DeclaredOrder::TemplateCoordinate)
    } else if fgumi_sam::is_sorted(header, COORDINATE) {
        Some(DeclaredOrder::Coordinate)
    } else if fgumi_sam::is_sorted(header, QUERY_NAME) {
        Some(DeclaredOrder::Queryname)
    } else {
        None
    }
}

/// The declared order an input must carry to be compatible with merge `order`.
fn expected_declared_order(order: SortOrderArg) -> DeclaredOrder {
    match order {
        SortOrderArg::Coordinate => DeclaredOrder::Coordinate,
        SortOrderArg::Queryname | SortOrderArg::QuerynameNatural => DeclaredOrder::Queryname,
        SortOrderArg::TemplateCoordinate => DeclaredOrder::TemplateCoordinate,
    }
}

/// The `fgumi sort --order` value string for an order (used in error hints).
fn order_flag_value(order: SortOrderArg) -> &'static str {
    match order {
        SortOrderArg::Coordinate => "coordinate",
        SortOrderArg::Queryname => "queryname",
        SortOrderArg::QuerynameNatural => "queryname::natural",
        SortOrderArg::TemplateCoordinate => "template-coordinate",
    }
}

/// MERGE3-01 (fast header check): reject an input whose header *declares* a sort
/// order that conflicts with the requested merge `order`.
///
/// The k-way merge only yields globally-sorted output if every input is already
/// sorted in `order`; a coordinate-sorted BAM fed to the default
/// `--order template-coordinate` merge (or any declared mismatch) silently
/// corrupts the output with a success exit. This catches the common footgun
/// before any records are read. Inputs that declare no usable order pass here
/// and are verified record-by-record during the merge instead.
///
/// # Errors
///
/// Returns an error if the header declares an order that conflicts with `order`.
fn check_input_declared_order(header: &Header, order: SortOrderArg, source: &str) -> Result<()> {
    let expected = expected_declared_order(order);
    if let Some(declared) = classify_declared_order(header) {
        if declared != expected {
            // `{source}` names the input for identification only; the remediation
            // command uses a `<input>` placeholder (matching the streaming-verify
            // error) so a path with shell metacharacters can't turn the
            // copy-pasteable hint into something unexpected.
            bail!(
                "Input '{source}' is sorted by {declared} but merge was asked for --order \
                 {requested}. Every input must already be sorted in the merge order, or the \
                 k-way merge silently corrupts the output.\n\nEither merge in the inputs' \
                 existing order:\n  fgumi merge --order {declared} ...\nor sort the inputs to \
                 {requested} first:\n  fgumi sort -i <input> -o sorted.bam --order {requested}",
                declared = declared.as_str(),
                requested = order_flag_value(order),
            );
        }
    }
    Ok(())
}

/// Merge headers from multiple BAM files.
///
/// Uses the first input's reference sequences and header line as the base.
/// Combines read groups and program records from all inputs, with earlier
/// inputs taking precedence for duplicate IDs (matching `samtools merge`).
/// Validates that all inputs share the same reference sequences (names and order)
/// and that each input's declared sort order is compatible with `order`
/// (MERGE3-01 fast check).
fn merge_headers(input_paths: &[PathBuf], order: SortOrderArg) -> Result<Header> {
    if input_paths.is_empty() {
        bail!("No input files to merge headers from");
    }

    // Read all headers once (raw reader; records aren't consumed here).
    let headers: Vec<Header> = input_paths
        .iter()
        .map(|path| {
            let (_, header) = create_raw_bam_reader(path, 1)?;
            Ok(header)
        })
        .collect::<Result<Vec<_>>>()?;

    // MERGE3-01: reject any input whose declared order conflicts with the merge
    // order (validated for every input, including the single-input case, since
    // the output header is stamped with the requested order regardless).
    for (path, header) in input_paths.iter().zip(headers.iter()) {
        check_input_declared_order(header, order, &path.display().to_string())?;
    }

    let first_header = &headers[0];

    if headers.len() == 1 {
        return Ok(first_header.clone());
    }

    // Verify reference sequences match across all inputs
    let first_refs = first_header.reference_sequences();
    for (i, (path, header)) in input_paths[1..].iter().zip(headers[1..].iter()).enumerate() {
        let other_refs = header.reference_sequences();
        if first_refs.len() != other_refs.len() {
            bail!(
                "Reference sequence count mismatch: {} has {} references, {} has {}",
                input_paths[0].display(),
                first_refs.len(),
                path.display(),
                other_refs.len()
            );
        }
        for ((name1, _), (name2, _)) in first_refs.iter().zip(other_refs.iter()) {
            if name1 != name2 {
                bail!(
                    "Reference sequence mismatch at input {}: '{}' has '{}', '{}' has '{}'",
                    i + 2,
                    input_paths[0].display(),
                    String::from_utf8_lossy(name1.as_ref()),
                    path.display(),
                    String::from_utf8_lossy(name2.as_ref()),
                );
            }
        }
    }

    let mut builder = Header::builder();

    // Reference sequences from first input
    for (name, seq) in first_header.reference_sequences() {
        builder = builder.add_reference_sequence(name.clone(), seq.clone());
    }

    // Header line from first input
    if let Some(hdr) = first_header.header() {
        builder = builder.set_header(hdr.clone());
    }

    // Collect read groups and programs from all inputs (first input wins ties)
    let mut rg_ids = HashSet::new();
    let mut pg_ids = HashSet::new();

    for header in &headers {
        for (id, rg) in header.read_groups() {
            if rg_ids.insert(id.clone()) {
                builder = builder.add_read_group(id.clone(), rg.clone());
            }
        }

        for (id, pg) in header.programs().as_ref() {
            if pg_ids.insert(id.clone()) {
                builder = builder.add_program(id.clone(), pg.clone());
            }
        }
    }

    // Comments from first input only (matching samtools behavior)
    for comment in first_header.comments() {
        builder = builder.add_comment(comment.clone());
    }

    Ok(builder.build())
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::BString;
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReadGroup;
    use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
    use std::num::NonZeroUsize;

    /// Create a BAM with the given read group IDs and return its path.
    fn write_bam_with_read_groups(dir: &std::path::Path, name: &str, rg_ids: &[&str]) -> PathBuf {
        use noodles::sam::header::record::value::map::{Program, ReferenceSequence};

        let mut header_builder = Header::builder();

        // Add one reference sequence so we can write mapped records
        let map = Map::<ReferenceSequence>::new(
            NonZeroUsize::new(200_000_000).expect("non-zero reference length"),
        );
        header_builder = header_builder.add_reference_sequence(BString::from("chr1"), map);

        for rg_id in rg_ids {
            let rg = Map::<ReadGroup>::builder()
                .insert(rg_tag::LIBRARY, format!("Lib_{rg_id}"))
                .build()
                .expect("valid RG");
            header_builder = header_builder.add_read_group(BString::from(*rg_id), rg);
        }

        // Add a program record
        let pg = Map::<Program>::default();
        header_builder = header_builder.add_program(BString::from("test"), pg);

        let header = header_builder.build();

        // Write a minimal BAM with this header (no records needed for header tests)
        let path = dir.join(format!("{name}.bam"));
        let file = std::fs::File::create(&path).expect("failed to create BAM file");
        let mut writer = noodles::bam::io::Writer::new(file);
        writer.write_header(&header).expect("failed to write BAM header");
        // write EOF
        drop(writer);
        path
    }

    #[test]
    fn test_merge_headers_single_input() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let bam = write_bam_with_read_groups(dir.path(), "single", &["RG1", "RG2"]);

        let header = merge_headers(std::slice::from_ref(&bam), SortOrderArg::Coordinate)
            .expect("merge_headers should succeed");

        // Should be the same header (same RGs)
        let rg_ids: Vec<String> =
            header.read_groups().iter().map(|(id, _)| id.to_string()).collect();
        assert_eq!(rg_ids.len(), 2);
        assert!(rg_ids.contains(&"RG1".to_string()));
        assert!(rg_ids.contains(&"RG2".to_string()));
    }

    #[test]
    fn test_merge_headers_combines_read_groups() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let bam_a = write_bam_with_read_groups(dir.path(), "a", &["RG1"]);
        let bam_b = write_bam_with_read_groups(dir.path(), "b", &["RG2"]);

        let header = merge_headers(&[bam_a, bam_b], SortOrderArg::Coordinate)
            .expect("merge_headers should succeed");

        let rg_ids: Vec<String> =
            header.read_groups().iter().map(|(id, _)| id.to_string()).collect();
        assert_eq!(rg_ids.len(), 2, "expected 2 read groups, got {rg_ids:?}");
        assert!(rg_ids.contains(&"RG1".to_string()));
        assert!(rg_ids.contains(&"RG2".to_string()));
    }

    #[test]
    fn test_merge_headers_deduplicates_read_groups() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        // Both BAMs have RG1, but with different library names
        let bam_a = write_bam_with_read_groups(dir.path(), "a", &["RG1"]);
        let bam_b = write_bam_with_read_groups(dir.path(), "b", &["RG1", "RG2"]);

        let header = merge_headers(&[bam_a, bam_b], SortOrderArg::Coordinate)
            .expect("merge_headers should succeed");

        let rg_ids: Vec<String> =
            header.read_groups().iter().map(|(id, _)| id.to_string()).collect();
        // RG1 from first input wins, RG2 from second input is added => 2 total
        assert_eq!(rg_ids.len(), 2, "expected 2 unique read groups, got {rg_ids:?}");
        assert!(rg_ids.contains(&"RG1".to_string()));
        assert!(rg_ids.contains(&"RG2".to_string()));

        // Verify RG1's library comes from the first input
        let (_, rg1) = header
            .read_groups()
            .iter()
            .find(|(id, _)| <BString as AsRef<[u8]>>::as_ref(id) == b"RG1")
            .expect("RG1 read group not found in merged header");
        let lib = rg1.other_fields().get(&rg_tag::LIBRARY).map(|v| v.to_string());
        assert_eq!(lib, Some("Lib_RG1".to_string()));
    }

    /// Create a BAM with the given reference sequence names and return its path.
    fn write_bam_with_refs(dir: &std::path::Path, name: &str, ref_names: &[&str]) -> PathBuf {
        use noodles::sam::header::record::value::map::ReferenceSequence;

        let mut header_builder = Header::builder();

        for ref_name in ref_names {
            let map = Map::<ReferenceSequence>::new(
                NonZeroUsize::new(200_000_000).expect("non-zero reference length"),
            );
            header_builder = header_builder.add_reference_sequence(BString::from(*ref_name), map);
        }

        let header = header_builder.build();

        let path = dir.join(format!("{name}.bam"));
        let file = std::fs::File::create(&path).expect("failed to create BAM file");
        let mut writer = noodles::bam::io::Writer::new(file);
        writer.write_header(&header).expect("failed to write BAM header");
        drop(writer);
        path
    }

    #[test]
    fn test_merge_headers_rejects_different_ref_count() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let bam_a = write_bam_with_refs(dir.path(), "a", &["chr1", "chr2"]);
        let bam_b = write_bam_with_refs(dir.path(), "b", &["chr1"]);

        let result = merge_headers(&[bam_a, bam_b], SortOrderArg::Coordinate);
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("Reference sequence count mismatch"), "unexpected error: {msg}");
    }

    #[test]
    fn test_merge_headers_rejects_different_ref_names() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let bam_a = write_bam_with_refs(dir.path(), "a", &["chr1", "chr2"]);
        let bam_b = write_bam_with_refs(dir.path(), "b", &["chr1", "chrX"]);

        let result = merge_headers(&[bam_a, bam_b], SortOrderArg::Coordinate);
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("Reference sequence mismatch"), "unexpected error: {msg}");
    }

    #[test]
    fn test_merge_output_aliasing_input() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let bam = write_bam_with_read_groups(dir.path(), "input", &["RG1"]);

        let merge = Merge {
            output: bam.clone(),
            inputs: vec![bam.clone()],
            input_list: None,
            order: SortOrderArg::Coordinate,
            threads: 1,
            compression_level: 6,
        };

        let result = merge.execute("test");
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("is the same as input file"), "unexpected error: {msg}");
    }

    #[test]
    fn test_merge_empty_bams_produces_output() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let bam_a = write_bam_with_read_groups(dir.path(), "a", &["RG1"]);
        let bam_b = write_bam_with_read_groups(dir.path(), "b", &["RG2"]);
        let output = dir.path().join("merged.bam");

        let merge = Merge {
            output: output.clone(),
            inputs: vec![bam_a, bam_b],
            input_list: None,
            order: SortOrderArg::Coordinate,
            threads: 1,
            compression_level: 6,
        };

        let result = merge.execute("test");
        assert!(result.is_ok(), "merge failed: {:?}", result.unwrap_err());
        assert!(output.exists(), "output BAM was not created");

        // Verify the output is a valid BAM by reading its header
        let (_, header) = create_bam_reader(&output, 1).expect("failed to read merged BAM output");
        let rg_ids: Vec<String> =
            header.read_groups().iter().map(|(id, _)| id.to_string()).collect();
        assert_eq!(rg_ids.len(), 2);
    }

    /// MERGE3-01 fast header check: an input whose header *declares* an order
    /// conflicting with `--order` is rejected; matching or undeclared inputs pass
    /// (undeclared ones are verified record-by-record during the merge instead).
    #[rstest::rstest]
    // Declared coordinate.
    #[case::coord_ok("@HD\tVN:1.6\tSO:coordinate\n", SortOrderArg::Coordinate, true)]
    #[case::coord_into_tc("@HD\tVN:1.6\tSO:coordinate\n", SortOrderArg::TemplateCoordinate, false)]
    #[case::coord_into_qname("@HD\tVN:1.6\tSO:coordinate\n", SortOrderArg::Queryname, false)]
    // Declared queryname (both lex and natural merges accept a queryname header).
    #[case::qname_ok("@HD\tVN:1.6\tSO:queryname\n", SortOrderArg::Queryname, true)]
    #[case::qname_natural_ok("@HD\tVN:1.6\tSO:queryname\n", SortOrderArg::QuerynameNatural, true)]
    #[case::qname_into_coord("@HD\tVN:1.6\tSO:queryname\n", SortOrderArg::Coordinate, false)]
    #[case::qname_into_tc("@HD\tVN:1.6\tSO:queryname\n", SortOrderArg::TemplateCoordinate, false)]
    // Declared template-coordinate.
    #[case::tc_ok(
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\tSS:template-coordinate\n",
        SortOrderArg::TemplateCoordinate,
        true
    )]
    #[case::tc_into_coord(
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\tSS:template-coordinate\n",
        SortOrderArg::Coordinate,
        false
    )]
    #[case::tc_into_qname(
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\tSS:template-coordinate\n",
        SortOrderArg::Queryname,
        false
    )]
    // Undeclared: bare, plain unsorted, or query-grouped-without-SS → pass here.
    #[case::bare_any("@HD\tVN:1.6\n", SortOrderArg::Coordinate, true)]
    #[case::unsorted_any("@HD\tVN:1.6\tSO:unsorted\n", SortOrderArg::TemplateCoordinate, true)]
    #[case::query_grouped_no_ss(
        "@HD\tVN:1.6\tSO:unsorted\tGO:query\n",
        SortOrderArg::TemplateCoordinate,
        true
    )]
    fn test_check_input_declared_order(
        #[case] header_str: &str,
        #[case] order: SortOrderArg,
        #[case] expect_ok: bool,
    ) {
        let header: Header = header_str.parse().expect("parse header");
        let result = check_input_declared_order(&header, order, "in.bam");
        assert_eq!(result.is_ok(), expect_ok, "header {header_str:?} into --order {order:?}");
        if let Err(e) = result {
            let msg = e.to_string();
            assert!(msg.contains("in.bam"), "message missing source: {msg}");
            assert!(msg.contains(order_flag_value(order)), "message missing order hint: {msg}");
        }
    }

    /// MERGE3-01 hardening: the copy-pasteable remediation command must use a
    /// fixed `<input>` placeholder, never the interpolated source path, so a
    /// filename with shell metacharacters can't turn the hint into an unintended
    /// command. The source is still named in the diagnostic text for identification.
    #[test]
    fn test_declared_order_error_does_not_interpolate_source_into_shell_command() {
        // A coordinate header into a queryname merge triggers the declared-order
        // conflict; the source path carries shell metacharacters.
        let header: Header = "@HD\tVN:1.6\tSO:coordinate\n".parse().expect("parse header");
        let malicious = "a'; rm -rf ~; echo '.bam";
        let err = check_input_declared_order(&header, SortOrderArg::Queryname, malicious)
            .expect_err("coordinate header into queryname merge must be rejected");
        let msg = err.to_string();

        // The remediation command uses the fixed placeholder...
        assert!(
            msg.contains("fgumi sort -i <input> -o sorted.bam"),
            "remediation must use the <input> placeholder: {msg}"
        );
        // ...and the source path is never spliced into any `-i` argument.
        assert!(
            !msg.contains(&format!("-i '{malicious}'")),
            "source path must not be interpolated into the sort command: {msg}"
        );
        assert!(
            !msg.contains(&format!("-i {malicious}")),
            "source path must not be interpolated into the sort command: {msg}"
        );
        // The source is still named for identification (in the descriptive text).
        assert!(msg.contains(malicious), "message should still name the source: {msg}");
    }
}
