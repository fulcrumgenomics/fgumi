//! Merge pre-sorted BAM files into a single sorted BAM.
//!
//! Performs a k-way merge of BAM files that are already sorted in the same
//! order, producing a single merged output that preserves the sort order.
//!
//! Similar to `samtools merge`, but supports template-coordinate order and
//! uses the same high-performance merge infrastructure as `fgumi sort`.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use crate::logging::OperationTimer;
use crate::sam::SamTag;
use crate::validation::validate_file_exists;
use anyhow::{Result, bail};
use bstr::{BStr, BString};
use clap::Parser;
#[cfg(test)]
use fgumi_bam_io::create_bam_reader;
use fgumi_bam_io::create_raw_bam_reader;
use fgumi_bam_io::header::{suffixed_id, with_renamed_reference};
use fgumi_sort::RawExternalSorter;
use indexmap::IndexMap;
use log::{info, warn};
use noodles::sam::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::program::tag as program_tag;
use noodles::sam::header::record::value::map::read_group::tag as read_group_tag;
use noodles::sam::header::record::value::map::{Program, ReadGroup};

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

@RG and @PG records from all inputs are combined. A record that appears with
identical content in several inputs is written once; records count as identical
only if the @PG their `PP` (@PG) or `PG` (@RG) names is identical too. When
inputs reuse an ID for different records (e.g. read group `A` with different
libraries), each distinct later record is written under a fresh ID (`A.1`,
`A.2`, ...), and the `RG`/`PG` tags on that input's reads, the `PP` of its @PG
records and the `PG` of its @RG records are rewritten to match, so every read
keeps its own read group and program. Each rename is logged. Unlike `samtools
merge`, identical records are combined without `-c`/`-p`, and fresh IDs are
deterministic. The output gets one @PG for this merge, chained to the last
input's program chain; `samtools merge` adds one per input chain.

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
    fn execute(&self, command_line: &str) -> Result<()> {
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
                if let Ok(input_canon) = std::fs::canonicalize(path)
                    && output_canon == input_canon
                {
                    bail!(
                        "Output file '{}' is the same as input file '{}'",
                        self.output.display(),
                        path.display()
                    );
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
        let MergedHeader { header, input_renames } = merge_headers(&input_paths, self.order)?;
        log_renames(&input_paths, &input_renames);
        // Record this merge with one @PG, as every fgumi command does, chained to the last
        // program chain end in merged-header order, which is the last input's when inputs
        // carry separate chains. `samtools merge` instead adds one @PG per chain end; fgumi
        // adds one per command (see `fgumi_bam_io::header::add_pg_record`).
        let header = crate::commands::common::add_pg_record(header, command_line)?;
        let mut rewriter = ReadTagRewriter::new(&input_paths, input_renames);

        let mut sorter = RawExternalSorter::new(self.order.into())
            .threads(self.threads)
            .output_compression(self.compression_level);

        if let Some(ct) = cell_tag {
            sorter = sorter.cell_tag(ct);
        }

        let records_merged =
            sorter.merge_bams_rewriting(&input_paths, &header, &self.output, |input, record| {
                rewriter.rewrite(input, record);
            })?;
        rewriter.log_undeclared_summary();

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
    if let Some(declared) = classify_declared_order(header)
        && declared != expected
    {
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
    Ok(())
}

/// The merged output header, and how each input's IDs map into it.
#[derive(Debug)]
struct MergedHeader {
    header: Header,
    /// One entry per input, in input order.
    input_renames: Vec<InputIdRenames>,
}

/// The @RG and @PG IDs one input's records are written under in the merged header,
/// for the IDs that changed, plus the IDs that input's header declares.
#[derive(Debug, Default)]
struct InputIdRenames {
    read_groups: HashMap<BString, BString>,
    programs: HashMap<BString, BString>,
    declared_read_groups: HashSet<BString>,
    declared_programs: HashSet<BString>,
}

impl InputIdRenames {
    /// Whether no ID of this input was renamed, so its reads pass through unchanged.
    fn is_empty(&self) -> bool {
        self.read_groups.is_empty() && self.programs.is_empty()
    }
}

/// Merge headers from multiple BAM files.
///
/// Uses the first input's reference sequences and header line as the base.
/// Combines read groups and program records from all inputs, in input order: a
/// record with the same ID and content as one already written is written once,
/// and a record whose ID is already used for different content is written under a
/// fresh ID (see [`merge_records`]). Returns, per input, the IDs that were renamed,
/// so its reads can be rewritten to follow them. Validates that all inputs share the
/// same reference sequences (names and order) and that each input's declared sort
/// order is compatible with `order` (MERGE3-01 fast check).
fn merge_headers(input_paths: &[PathBuf], order: SortOrderArg) -> Result<MergedHeader> {
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
        return Ok(MergedHeader {
            header: first_header.clone(),
            input_renames: vec![InputIdRenames::default()],
        });
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

    // Fresh IDs avoid every ID declared by any input, so a later input's own ID can
    // never collide with one generated for an earlier input.
    let mut reserved_read_group_ids: HashSet<BString> =
        headers.iter().flat_map(|h| h.read_groups().keys().cloned()).collect();
    let mut reserved_program_ids: HashSet<BString> =
        headers.iter().flat_map(|h| h.programs().as_ref().keys().cloned()).collect();

    let mut read_groups: IndexMap<BString, Map<ReadGroup>> = IndexMap::new();
    let mut programs: IndexMap<BString, Map<Program>> = IndexMap::new();
    let mut minted_read_group_ids = MintedIds::new();
    let mut minted_program_ids = MintedIds::new();
    let mut input_renames = Vec::with_capacity(headers.len());
    for header in &headers {
        let input_programs: Vec<(&BString, &Map<Program>)> =
            header.programs().as_ref().iter().collect();
        let program_renames = merge_records(
            &input_programs,
            &mut programs,
            &mut minted_program_ids,
            &mut reserved_program_ids,
            |pg, renames| with_renamed_reference(pg, program_tag::PREVIOUS_PROGRAM_ID, renames),
        );
        let input_read_groups: Vec<(&BString, &Map<ReadGroup>)> =
            header.read_groups().iter().collect();
        let read_group_renames = merge_records(
            &input_read_groups,
            &mut read_groups,
            &mut minted_read_group_ids,
            &mut reserved_read_group_ids,
            |rg, _| with_renamed_reference(rg, read_group_tag::PROGRAM, &program_renames),
        );
        input_renames.push(InputIdRenames {
            read_groups: read_group_renames,
            programs: program_renames,
            declared_read_groups: header.read_groups().keys().cloned().collect(),
            declared_programs: header.programs().as_ref().keys().cloned().collect(),
        });
    }
    for (id, rg) in read_groups {
        builder = builder.add_read_group(id, rg);
    }
    for (id, pg) in programs {
        builder = builder.add_program(id, pg);
    }

    // Comments from first input only (matching samtools behavior)
    for comment in first_header.comments() {
        builder = builder.add_comment(comment.clone());
    }

    Ok(MergedHeader { header: builder.build(), input_renames })
}

/// For each original @RG/@PG ID, the fresh IDs earlier inputs' records were written
/// under, so a later input holding the same record can reuse one instead of minting
/// another.
type MintedIds = HashMap<BString, Vec<BString>>;

/// Adds one input's @RG or @PG records to the merged `output` records, and returns the
/// IDs it renamed.
///
/// `translate` rewrites a record's cross-references (`PP` for @PG, `PG` for @RG) through
/// this input's renames. A record whose translated content equals the output record with
/// its ID is combined with it. Otherwise, if its ID is already used in `output`, it is
/// combined with an identical record an earlier input wrote under a fresh ID minted from
/// the same original ID (`minted`), or written under a new fresh `{id}.{n}` not in
/// `reserved`. A record whose ID is not yet used keeps it.
///
/// Renaming a @PG changes the translation of a record whose `PP` names it, which can
/// change that record's decision in turn, so every decision is recomputed from the
/// previous pass's renames until they stop changing. Records are written only after
/// that, so nothing written by an earlier input is ever replaced. A `PP` cycle can keep
/// the decisions from settling; then every record whose ID is used is written under a
/// fresh ID, which is always consistent. Records are appended in input order.
fn merge_records<T>(
    records: &[(&BString, &Map<T>)],
    output: &mut IndexMap<BString, Map<T>>,
    minted: &mut MintedIds,
    reserved: &mut HashSet<BString>,
    translate: impl Fn(&Map<T>, &HashMap<BString, BString>) -> Map<T>,
) -> HashMap<BString, BString>
where
    T: noodles::sam::header::record::value::map::Inner,
    Map<T>: PartialEq,
{
    // Fresh IDs are only reserved once used, so this is stable across passes, and
    // distinct original IDs never share one (see `suffixed_id`).
    let fresh = |id: &BString| suffixed_id(id, |candidate| reserved.contains(candidate));

    let mut renames: HashMap<BString, BString> = HashMap::new();
    let mut settled = false;
    for _ in 0..=records.len() {
        let mut next: HashMap<BString, BString> = HashMap::new();
        for (id, record) in records {
            let Some(existing) = output.get(*id) else { continue };
            let translated = translate(record, &renames);
            if translated == *existing {
                continue;
            }
            let reused = minted
                .get(*id)
                .and_then(|ids| ids.iter().find(|m| output.get(*m) == Some(&translated)));
            next.insert((*id).clone(), reused.cloned().unwrap_or_else(|| fresh(id)));
        }
        if next == renames {
            settled = true;
            break;
        }
        renames = next;
    }
    if !settled {
        renames = records
            .iter()
            .filter(|(id, _)| output.contains_key(*id))
            .map(|(id, _)| ((*id).clone(), fresh(id)))
            .collect();
    }

    for (id, record) in records {
        let translated = translate(record, &renames);
        let target = renames.get(*id).unwrap_or(id);
        if let Some(existing) = output.get(target) {
            // Settled decisions combine only with an identical record.
            debug_assert!(*existing == translated, "combined records must be identical");
            continue;
        }
        if renames.contains_key(*id) {
            reserved.insert(target.clone());
            minted.entry((*id).clone()).or_default().push(target.clone());
        }
        output.insert(target.clone(), translated);
    }
    renames
}

/// Logs, per input, each @RG/@PG ID written under a different ID in the merged header.
fn log_renames(input_paths: &[PathBuf], input_renames: &[InputIdRenames]) {
    for (path, renames) in input_paths.iter().zip(input_renames) {
        for (record_type, ids) in [("@RG", &renames.read_groups), ("@PG", &renames.programs)] {
            let mut ids: Vec<_> = ids.iter().collect();
            ids.sort();
            for (from, to) in ids {
                info!(
                    "{}: {record_type} {from} differs from an earlier input's {from}; \
                     writing it, and this input's reads' tags, as {to}",
                    path.display()
                );
            }
        }
    }
}

/// Rewrites the `RG` and `PG` tags of each input's reads to the IDs that input's read
/// groups and programs are written under in the merged header.
///
/// When no input has a renamed ID, every read passes through untouched, at no
/// per-record cost, so an undeclared `PG` tag that happens to equal the ID of the
/// merge's own @PG (`fgumi`, ...) goes unreported and resolves to it. Otherwise every read's tags are checked: a tag naming a renamed ID
/// is rewritten, and a tag naming an ID its input's header does not declare is left
/// unchanged and reported (once per input and tag, then as a count). When that ID is a
/// fresh ID minted for a renamed record, the read now joins that record, so the
/// collision is reported separately.
struct ReadTagRewriter<'a> {
    input_paths: &'a [PathBuf],
    input_renames: Vec<InputIdRenames>,
    any_renames: bool,
    /// The fresh @RG IDs minted for any input's renamed read groups.
    fresh_read_group_ids: HashSet<BString>,
    /// The fresh @PG IDs minted for any input's renamed programs.
    fresh_program_ids: HashSet<BString>,
    /// Reads seen with an undeclared tag value, per `(input, tag)`.
    undeclared: HashMap<(usize, SamTag), u64>,
    /// Reads seen with an undeclared tag value that is a fresh ID, per `(input, tag)`.
    fresh_id_collisions: HashMap<(usize, SamTag), u64>,
}

impl<'a> ReadTagRewriter<'a> {
    fn new(input_paths: &'a [PathBuf], input_renames: Vec<InputIdRenames>) -> Self {
        let any_renames = input_renames.iter().any(|renames| !renames.is_empty());
        let fresh_read_group_ids = input_renames
            .iter()
            .flat_map(|renames| renames.read_groups.values().cloned())
            .collect();
        let fresh_program_ids =
            input_renames.iter().flat_map(|renames| renames.programs.values().cloned()).collect();
        Self {
            input_paths,
            input_renames,
            any_renames,
            fresh_read_group_ids,
            fresh_program_ids,
            undeclared: HashMap::new(),
            fresh_id_collisions: HashMap::new(),
        }
    }

    /// Rewrites the `RG` and `PG` tags of `record`, a raw BAM record from input `input`.
    fn rewrite(&mut self, input: usize, record: &mut Vec<u8>) {
        if !self.any_renames {
            return;
        }
        self.rewrite_tag(input, record, SamTag::RG);
        self.rewrite_tag(input, record, SamTag::PG);
    }

    fn rewrite_tag(&mut self, input: usize, record: &mut Vec<u8>, tag: SamTag) {
        let renames = &self.input_renames[input];
        let (ids, declared, fresh, record_type) = if tag == SamTag::RG {
            (&renames.read_groups, &renames.declared_read_groups, &self.fresh_read_group_ids, "@RG")
        } else {
            (&renames.programs, &renames.declared_programs, &self.fresh_program_ids, "@PG")
        };
        let aux = fgumi_raw_bam::aux_data_slice(record);
        let Some(value) = fgumi_raw_bam::find_string_tag(aux, tag) else { return };
        let value = BStr::new(value);
        if let Some(renamed) = ids.get(value) {
            fgumi_raw_bam::update_string_tag(record, tag, renamed);
        } else if !declared.contains(value) {
            let [t0, t1] = *tag;
            if fresh.contains(value) {
                let count = self.fresh_id_collisions.entry((input, tag)).or_insert(0);
                *count += 1;
                if *count == 1 {
                    warn!(
                        "{}: a read's {}{}:Z:{value} tag names no {record_type} in that input's \
                         header, but {value} is the fresh ID a renamed {record_type} is written \
                         under; leaving it unchanged, so the read joins that {record_type}",
                        self.input_paths[input].display(),
                        t0 as char,
                        t1 as char,
                    );
                }
            }
            let count = self.undeclared.entry((input, tag)).or_insert(0);
            *count += 1;
            if *count == 1 {
                warn!(
                    "{}: a read's {}{}:Z:{value} tag names no {record_type} in that input's \
                     header; leaving it unchanged",
                    self.input_paths[input].display(),
                    t0 as char,
                    t1 as char,
                );
            }
        }
    }

    /// Logs how many reads, per input and tag, named an undeclared ID, and how many of
    /// those named a fresh ID, where more than the one already warned about did.
    fn log_undeclared_summary(&self) {
        for (counts, what) in [
            (&self.undeclared, "naming no record in that input's header"),
            (&self.fresh_id_collisions, "naming a renamed record's fresh ID"),
        ] {
            let mut counts: Vec<_> = counts.iter().filter(|(_, count)| **count > 1).collect();
            counts.sort_by_key(|((input, tag), _)| (*input, **tag));
            for ((input, tag), count) in counts {
                let [t0, t1] = **tag;
                warn!(
                    "{}: {count} read(s) had an {}{} tag {what}",
                    self.input_paths[*input].display(),
                    t0 as char,
                    t1 as char,
                );
            }
        }
    }
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
            .expect("merge_headers should succeed")
            .header;

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
            .expect("merge_headers should succeed")
            .header;

        let rg_ids: Vec<String> =
            header.read_groups().iter().map(|(id, _)| id.to_string()).collect();
        assert_eq!(rg_ids.len(), 2, "expected 2 read groups, got {rg_ids:?}");
        assert!(rg_ids.contains(&"RG1".to_string()));
        assert!(rg_ids.contains(&"RG2".to_string()));
    }

    #[test]
    fn test_merge_headers_deduplicates_read_groups() {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        // Both BAMs have an identical RG1 (`write_bam_with_read_groups` derives the
        // library from the ID), so it is written once
        let bam_a = write_bam_with_read_groups(dir.path(), "a", &["RG1"]);
        let bam_b = write_bam_with_read_groups(dir.path(), "b", &["RG1", "RG2"]);

        let header = merge_headers(&[bam_a, bam_b], SortOrderArg::Coordinate)
            .expect("merge_headers should succeed")
            .header;

        let rg_ids: Vec<String> =
            header.read_groups().iter().map(|(id, _)| id.to_string()).collect();
        // The shared RG1 is combined, RG2 from the second input is added => 2 total
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

    /// Writes a header-only BAM whose header is one `chr1` reference plus `lines`, the
    /// tab-separated @RG/@PG lines, and returns its path.
    fn write_bam_with_header_lines(dir: &std::path::Path, name: &str, lines: &[&str]) -> PathBuf {
        let mut text = String::from("@SQ\tSN:chr1\tLN:1000\n");
        for line in lines {
            text.push_str(line);
            text.push('\n');
        }
        let header: Header = text.parse().expect("valid SAM header text");
        let path = dir.join(format!("{name}.bam"));
        let file = std::fs::File::create(&path).expect("failed to create BAM file");
        let mut writer = noodles::bam::io::Writer::new(file);
        writer.write_header(&header).expect("failed to write BAM header");
        drop(writer);
        path
    }

    /// The @RG and @PG lines of `header`, in header order, as SAM text.
    fn read_group_and_program_lines(header: &Header) -> Vec<String> {
        let mut writer = noodles::sam::io::Writer::new(Vec::new());
        writer.write_header(header).expect("failed to write SAM header");
        String::from_utf8(writer.get_ref().clone())
            .expect("SAM header is UTF-8")
            .lines()
            .filter(|line| line.starts_with("@RG") || line.starts_with("@PG"))
            .map(str::to_string)
            .collect()
    }

    /// One input's expected renamed `(from, to)` @RG IDs, then @PG IDs.
    type RenamedIds = (Vec<(&'static str, &'static str)>, Vec<(&'static str, &'static str)>);

    /// Sorted `(from, to)` ID pairs.
    type IdPairs = Vec<(String, String)>;

    /// The renamed @RG and @PG IDs of one input, as sorted `(from, to)` pairs.
    fn renamed_ids(renames: &InputIdRenames) -> (IdPairs, IdPairs) {
        let sorted = |ids: &HashMap<BString, BString>| {
            let mut pairs: Vec<(String, String)> =
                ids.iter().map(|(from, to)| (from.to_string(), to.to_string())).collect();
            pairs.sort();
            pairs
        };
        (sorted(&renames.read_groups), sorted(&renames.programs))
    }

    /// Inputs that reuse an @RG/@PG ID: identical records are written once, and a record
    /// whose ID is already used for different content gets a fresh `{id}.{n}` ID, with
    /// its input's `PP` (@PG) and `PG` (@RG) references rewritten to follow it.
    /// `expected_renames` lists, per input, the renamed `(from, to)` @RG and @PG IDs.
    #[rstest::rstest]
    #[case::identical_records_combined(
        vec![
            vec!["@RG\tID:A\tLB:x", "@PG\tID:bwa\tCL:bwa ref.fa"],
            vec!["@RG\tID:A\tLB:x", "@PG\tID:bwa\tCL:bwa ref.fa"],
        ],
        vec!["@RG\tID:A\tLB:x", "@PG\tID:bwa\tCL:bwa ref.fa"],
        vec![(vec![], vec![]), (vec![], vec![])],
    )]
    #[case::conflicting_read_group_renamed(
        vec![vec!["@RG\tID:A\tLB:x"], vec!["@RG\tID:A\tLB:y", "@RG\tID:B\tLB:z"]],
        vec!["@RG\tID:A\tLB:x", "@RG\tID:A.1\tLB:y", "@RG\tID:B\tLB:z"],
        vec![(vec![], vec![]), (vec![("A", "A.1")], vec![])],
    )]
    #[case::fresh_id_skips_ids_of_the_same_input(
        vec![vec!["@RG\tID:A\tLB:x"], vec!["@RG\tID:A\tLB:y", "@RG\tID:A.1\tLB:z"]],
        vec!["@RG\tID:A\tLB:x", "@RG\tID:A.2\tLB:y", "@RG\tID:A.1\tLB:z"],
        vec![(vec![], vec![]), (vec![("A", "A.2")], vec![])],
    )]
    #[case::fresh_id_skips_ids_of_later_inputs(
        vec![vec!["@RG\tID:A\tLB:x"], vec!["@RG\tID:A\tLB:y"], vec!["@RG\tID:A.1\tLB:z"]],
        vec!["@RG\tID:A\tLB:x", "@RG\tID:A.2\tLB:y", "@RG\tID:A.1\tLB:z"],
        vec![(vec![], vec![]), (vec![("A", "A.2")], vec![]), (vec![], vec![])],
    )]
    #[case::program_rename_rewrites_pp(
        vec![
            vec!["@PG\tID:bwa\tCL:bwa ref1.fa"],
            vec!["@PG\tID:bwa\tCL:bwa ref2.fa", "@PG\tID:fgumi\tPP:bwa\tCL:fgumi zipper"],
        ],
        vec![
            "@PG\tID:bwa\tCL:bwa ref1.fa",
            "@PG\tID:bwa.1\tCL:bwa ref2.fa",
            "@PG\tID:fgumi\tPP:bwa.1\tCL:fgumi zipper",
        ],
        vec![(vec![], vec![]), (vec![], vec![("bwa", "bwa.1")])],
    )]
    // `samtools` is listed before the `bwa` it follows, so a single pass in header order
    // would compare it before `bwa` is renamed and wrongly combine it.
    #[case::program_rename_cascades_through_pp(
        vec![
            vec!["@PG\tID:bwa\tCL:bwa ref1.fa", "@PG\tID:samtools\tPP:bwa\tCL:samtools view"],
            vec!["@PG\tID:samtools\tPP:bwa\tCL:samtools view", "@PG\tID:bwa\tCL:bwa ref2.fa"],
        ],
        vec![
            "@PG\tID:bwa\tCL:bwa ref1.fa",
            "@PG\tID:samtools\tPP:bwa\tCL:samtools view",
            "@PG\tID:samtools.1\tPP:bwa.1\tCL:samtools view",
            "@PG\tID:bwa.1\tCL:bwa ref2.fa",
        ],
        vec![(vec![], vec![]), (vec![], vec![("bwa", "bwa.1"), ("samtools", "samtools.1")])],
    )]
    #[case::read_group_follows_renamed_program(
        vec![
            vec!["@RG\tID:A\tLB:x\tPG:bwa", "@PG\tID:bwa\tCL:bwa ref1.fa"],
            vec!["@RG\tID:A\tLB:x\tPG:bwa", "@PG\tID:bwa\tCL:bwa ref2.fa"],
        ],
        vec![
            "@RG\tID:A\tLB:x\tPG:bwa",
            "@RG\tID:A.1\tLB:x\tPG:bwa.1",
            "@PG\tID:bwa\tCL:bwa ref1.fa",
            "@PG\tID:bwa.1\tCL:bwa ref2.fa",
        ],
        vec![(vec![], vec![]), (vec![("A", "A.1")], vec![("bwa", "bwa.1")])],
    )]
    // Later inputs that share a definition conflicting with the first input's reuse one
    // fresh ID rather than each minting their own.
    #[case::later_inputs_share_a_fresh_id(
        vec![vec!["@RG\tID:A\tLB:x"], vec!["@RG\tID:A\tLB:y"], vec!["@RG\tID:A\tLB:y"]],
        vec!["@RG\tID:A\tLB:x", "@RG\tID:A.1\tLB:y"],
        vec![(vec![], vec![]), (vec![("A", "A.1")], vec![]), (vec![("A", "A.1")], vec![])],
    )]
    #[case::third_definition_gets_its_own_fresh_id(
        vec![vec!["@RG\tID:A\tLB:x"], vec!["@RG\tID:A\tLB:y"], vec!["@RG\tID:A\tLB:z"]],
        vec!["@RG\tID:A\tLB:x", "@RG\tID:A.1\tLB:y", "@RG\tID:A.2\tLB:z"],
        vec![(vec![], vec![]), (vec![("A", "A.1")], vec![]), (vec![("A", "A.2")], vec![])],
    )]
    // The third input repeats the second's renamed chain, `samtools` listed before the
    // `bwa` it follows: `samtools` matches `samtools.1` only once `bwa` resolves to
    // `bwa.1`, so the reuse must be decided after the cascade, not before it.
    #[case::later_inputs_share_a_renamed_pp_chain(
        vec![
            vec!["@PG\tID:bwa\tCL:bwa ref1.fa", "@PG\tID:samtools\tPP:bwa\tCL:samtools view"],
            vec!["@PG\tID:bwa\tCL:bwa ref2.fa", "@PG\tID:samtools\tPP:bwa\tCL:samtools view"],
            vec!["@PG\tID:samtools\tPP:bwa\tCL:samtools view", "@PG\tID:bwa\tCL:bwa ref2.fa"],
        ],
        vec![
            "@PG\tID:bwa\tCL:bwa ref1.fa",
            "@PG\tID:samtools\tPP:bwa\tCL:samtools view",
            "@PG\tID:bwa.1\tCL:bwa ref2.fa",
            "@PG\tID:samtools.1\tPP:bwa.1\tCL:samtools view",
        ],
        vec![
            (vec![], vec![]),
            (vec![], vec![("bwa", "bwa.1"), ("samtools", "samtools.1")]),
            (vec![], vec![("bwa", "bwa.1"), ("samtools", "samtools.1")]),
        ],
    )]
    fn test_merge_headers_renames_conflicting_ids(
        #[case] inputs: Vec<Vec<&str>>,
        #[case] expected_lines: Vec<&str>,
        #[case] expected_renames: Vec<RenamedIds>,
    ) {
        let dir = tempfile::tempdir().expect("failed to create temp dir");
        let paths: Vec<PathBuf> = inputs
            .iter()
            .enumerate()
            .map(|(i, lines)| write_bam_with_header_lines(dir.path(), &format!("in{i}"), lines))
            .collect();

        let merged =
            merge_headers(&paths, SortOrderArg::Coordinate).expect("merge_headers should succeed");

        assert_eq!(read_group_and_program_lines(&merged.header), expected_lines);
        let owned = |pairs: &[(&str, &str)]| -> Vec<(String, String)> {
            pairs.iter().map(|(from, to)| ((*from).to_string(), (*to).to_string())).collect()
        };
        let actual: Vec<_> = merged.input_renames.iter().map(renamed_ids).collect();
        let expected: Vec<_> =
            expected_renames.iter().map(|(rg, pg)| (owned(rg), owned(pg))).collect();
        assert_eq!(actual, expected);
    }

    /// When the rename decisions never settle (as a `PP` cycle could make them), every
    /// record whose ID is taken is written under a fresh ID, with the content translated
    /// through those renames, rather than combined with a record it may not match.
    #[test]
    fn test_merge_records_falls_back_to_fresh_ids_when_decisions_do_not_settle() {
        let record = |library: &str| {
            Map::<ReadGroup>::builder().insert(rg_tag::LIBRARY, library).build().expect("valid RG")
        };
        let id = BString::from("A");
        let mut output: IndexMap<BString, Map<ReadGroup>> =
            IndexMap::from([(id.clone(), record("x"))]);
        let mut minted = MintedIds::new();
        let mut reserved: HashSet<BString> = HashSet::from([id.clone()]);
        let input = record("in");
        // Differs from the output record with no renames and matches it with any, so
        // every pass flips the previous pass's decision.
        let translate = |_: &Map<ReadGroup>, renames: &HashMap<BString, BString>| {
            if renames.is_empty() { record("y") } else { record("x") }
        };

        let renames =
            merge_records(&[(&id, &input)], &mut output, &mut minted, &mut reserved, translate);

        assert_eq!(renames, HashMap::from([(id.clone(), BString::from("A.1"))]));
        let ids: Vec<&BString> = output.keys().collect();
        assert_eq!(ids, [&id, &BString::from("A.1")]);
        assert_eq!(output[&BString::from("A.1")], record("x"));
        assert!(reserved.contains(&BString::from("A.1")));
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

    /// A read whose tag names an ID its own input does not declare is left unchanged,
    /// and counted as a collision only when that ID is a fresh ID minted for a renamed
    /// record, since the read then silently joins that record.
    #[rstest::rstest]
    #[case::rg_names_another_inputs_fresh_id(0, SamTag::RG, b"A.1", true)]
    #[case::pg_names_another_inputs_fresh_id(0, SamTag::PG, b"P.1", true)]
    #[case::rg_names_its_own_inputs_fresh_id(1, SamTag::RG, b"A.1", true)]
    #[case::pg_names_its_own_inputs_fresh_id(1, SamTag::PG, b"P.1", true)]
    #[case::rg_names_an_unrelated_id(0, SamTag::RG, b"B", false)]
    #[case::pg_names_an_unrelated_id(0, SamTag::PG, b"Q", false)]
    fn test_read_tag_rewriter_reports_undeclared_tags_naming_fresh_ids(
        #[case] input: usize,
        #[case] tag: SamTag,
        #[case] value: &[u8],
        #[case] collides: bool,
    ) {
        let input_paths = [PathBuf::from("a.bam"), PathBuf::from("b.bam")];
        let input_renames = vec![
            InputIdRenames {
                declared_read_groups: HashSet::from([BString::from("A")]),
                declared_programs: HashSet::from([BString::from("P")]),
                ..InputIdRenames::default()
            },
            InputIdRenames {
                read_groups: HashMap::from([(BString::from("A"), BString::from("A.1"))]),
                programs: HashMap::from([(BString::from("P"), BString::from("P.1"))]),
                declared_read_groups: HashSet::from([BString::from("A")]),
                declared_programs: HashSet::from([BString::from("P")]),
            },
        ];
        let mut rewriter = ReadTagRewriter::new(&input_paths, input_renames);
        let mut record =
            fgumi_raw_bam::SamBuilder::new().add_string_tag(tag, value).build().into_inner();

        rewriter.rewrite(input, &mut record);

        let aux = fgumi_raw_bam::aux_data_slice(&record);
        assert_eq!(fgumi_raw_bam::find_string_tag(aux, tag), Some(value), "tag left unchanged");
        assert_eq!(rewriter.undeclared.get(&(input, tag)), Some(&1));
        assert_eq!(rewriter.fresh_id_collisions.get(&(input, tag)).copied(), collides.then_some(1));
    }
}
