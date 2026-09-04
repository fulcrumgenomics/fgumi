//! Downsample BAM files by UMI family using a single-pass streaming algorithm.
//!
//! This tool reads a BAM file that has been processed by fgumi group (or fgbio `GroupReadsByUmi`),
//! uniformly samples UMI families based on the MI tag, and outputs kept reads directly to a BAM file.
//! By default the two duplex strands of a molecule (`<base>/A` and `<base>/B`) are sampled as one
//! family; `--per-strand` restores the legacy behavior of sampling each raw MI tag independently.
//!
//! Requires input BAM to be in template-coordinate order (from group).

use crate::logging::OperationTimer;
use crate::sam::{SamTag, is_template_coordinate_sorted};
use anyhow::{Result, bail};
use clap::Parser;
use fgumi_bam_io::ProgressTracker;
use fgumi_bam_io::{RawBamWriter, create_raw_bam_reader_with_opts, create_raw_bam_writer};
use fgumi_raw_bam::{
    RawBamReader, RawRecord, aux_data_slice, find_int_tag, find_string_tag, read_name,
};
use log::info;
use rand::SeedableRng;
use rand::rngs::StdRng;
use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::commands::command::Command;
use crate::commands::common::{BamIoOptions, CompressionOptions, reject_output_collisions};
use crate::umi::extract_mi_base;

/// Downsample a BAM file by UMI family using streaming.
///
/// Reads consecutive records sharing the same MI tag, randomly samples families
/// based on the keep fraction, and writes kept families directly to output BAM.
#[derive(Debug, Parser)]
#[command(
    name = "downsample",
    about = "\x1b[38;5;166m[UTILITIES]\x1b[0m      \x1b[36mDownsample BAM by UMI family using streaming\x1b[0m",
    long_about = r#"
Downsample a BAM file by UMI family using a single-pass streaming algorithm.

This tool reads a BAM file that has been processed by fgumi group (or fgbio GroupReadsByUmi)
containing MI tags, uniformly samples UMI families, and outputs kept reads directly to a BAM file.

Requires input BAM to be in template-coordinate order:
  - SO:unsorted (or not set)
  - GO:query
  - SS:unsorted:template-coordinate or SS:template-coordinate

The tool processes families in streaming fashion by grouping consecutive reads with the same
MI tag value. For each family, a random decision is made based on the fraction parameter to
either keep or reject all reads in that family.

NOTE: this command is NOT a port of Picard `DownsampleSam`. It is a UMI-family sampler by
design, and differs from `DownsampleSam` in several deliberate ways:
  - Sampling unit: whole UMI families (MI tag), not individual read templates. Every record
    must carry an MI tag and the input must be group-produced, template-coordinate-sorted BAM.
    `DownsampleSam` downsamples any BAM in any order by read name.
  - Decision function: one sequential random draw per family (order-dependent), rather than
    `DownsampleSam`'s stateless, name+seed hash (order-independent, cross-tool reproducible).
  - Determinism: without --seed, a run is non-deterministic (seeded from entropy), whereas
    `DownsampleSam` defaults to a fixed seed of 1. Pass --seed for reproducible output.
  - --fraction must be in (0.0, 1.0]; 0.0 is rejected (unlike `DownsampleSam`, which allows
    PROBABILITY=0 for an empty pass). NaN and infinities are also rejected.
Because the sampling unit, decision function, and RNG all differ, output is intentionally not
bit-identical to `DownsampleSam`; only a statistically-equivalent fraction of families is kept.

DUPLEX DATA: the `paired` grouping strategy tags the two strands of a molecule `<base>/A` and
`<base>/B`. By default downsample reduces each MI to its molecule base (the same last-`/`
truncation `group` and `duplex` use) and samples whole molecules — both strands kept or dropped
together — so duplex families survive at the intended fraction. Simplex MIs (no `/A`,`/B` suffix)
are a no-op under this rule and unaffected. Pass `--per-strand` to instead sample each raw MI tag
independently (legacy behavior); on duplex data this collapses duplex families, because a molecule
then survives as a duplex only with probability fraction^2 — use it only when you deliberately want
strand-level sampling.

Example usage:
  fgumi downsample -i grouped.bam -o downsampled.bam -f 0.1 --seed 42
  fgumi downsample -i grouped.bam -o kept.bam -f 0.5 --rejects rejected.bam
  fgumi downsample -i grouped.bam -o kept.bam -f 0.1 --histogram-kept kept_hist.txt
  fgumi downsample -i duplex.grouped.bam -o kept.bam -f 0.1 --per-strand   # legacy per-strand
"#
)]
pub struct Downsample {
    /// Input/output BAM options
    #[command(flatten)]
    pub io: BamIoOptions,

    /// Fraction of UMI families to keep (0.0 exclusive to 1.0 inclusive)
    #[arg(short = 'f', long = "fraction")]
    pub fraction: f64,

    /// Optional output BAM file for rejected reads
    #[arg(long = "rejects")]
    pub rejects: Option<PathBuf>,

    /// Random seed for reproducibility
    #[arg(long = "seed")]
    pub seed: Option<u64>,

    /// Validate that MI tags appear in consecutive groups (error if seen non-consecutively)
    #[arg(long = "validate-mi-order", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub validate_mi_order: bool,

    /// Sample by raw MI tag rather than by molecule, sampling each duplex strand independently.
    ///
    /// By default downsample groups by molecule: each MI is reduced to its molecule base
    /// before grouping (the same last-`/` truncation `group`/`duplex` use), so the two
    /// strands the `paired` (duplex) strategy tags `<base>/A` and `<base>/B` form ONE family
    /// sharing a single keep/reject draw and are kept or dropped together. Simplex MIs (no
    /// `/A`,`/B` suffix) have no `/` to strip and are unaffected either way. With this flag
    /// each raw MI value is a separate family with its own draw — so a duplex molecule
    /// survives as a duplex only if BOTH strand draws hit (probability `fraction`^2),
    /// collapsing duplex families at low fractions. Use only when you deliberately want
    /// strand-level sampling. Off by default.
    #[arg(long = "per-strand", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub per_strand: bool,

    /// Output file for kept family size histogram
    #[arg(long = "histogram-kept")]
    pub histogram_kept: Option<PathBuf>,

    /// Output file for rejected family size histogram
    #[arg(long = "histogram-rejected")]
    pub histogram_rejected: Option<PathBuf>,

    /// Compression options for output BAM.
    #[command(flatten)]
    pub compression: CompressionOptions,
}

impl Downsample {
    /// Validates that `--fraction` is a finite value in `(0.0, 1.0]`.
    ///
    /// The non-finite check must come first: `NaN` compares `false` against both
    /// `<= 0.0` and `> 1.0`, so a bare range test admits it. A `NaN` fraction then
    /// makes every `rng < fraction` sampling decision false, silently producing an
    /// empty output BAM with exit code 0.
    ///
    /// # Errors
    ///
    /// Returns an error if `fraction` is not finite or falls outside `(0.0, 1.0]`.
    fn validate_fraction(fraction: f64) -> Result<()> {
        if !fraction.is_finite() {
            bail!("--fraction must be a finite number, got {fraction}");
        }
        if fraction <= 0.0 || fraction > 1.0 {
            bail!("--fraction must be between 0.0 (exclusive) and 1.0 (inclusive), got {fraction}");
        }
        Ok(())
    }
}

impl Command for Downsample {
    fn execute(&self, command_line: &str) -> Result<()> {
        // Validate inputs. Use the shared validator so stdin (`-` / `/dev/stdin`)
        // is exempt from the file-existence check, matching every other streaming
        // command; the BAM reader already handles stdin.
        self.io.validate()?;
        let mut outputs: Vec<(&Path, &str)> = vec![(self.io.output.as_path(), "--output")];
        if let Some(path) = &self.rejects {
            outputs.push((path.as_path(), "--rejects"));
        }
        if let Some(path) = &self.histogram_kept {
            outputs.push((path.as_path(), "--histogram-kept"));
        }
        if let Some(path) = &self.histogram_rejected {
            outputs.push((path.as_path(), "--histogram-rejected"));
        }
        reject_output_collisions(&outputs)?;

        // Validate fraction
        Self::validate_fraction(self.fraction)?;

        let timer = OperationTimer::new("Downsampling reads");

        info!("Starting Downsample");
        info!("Input: {}", self.io.input.display());
        info!("Output: {}", self.io.output.display());
        info!("Target fraction: {}", self.fraction);
        if let Some(seed) = self.seed {
            info!("Random seed: {seed}");
        }
        if self.validate_mi_order {
            info!("MI order validation: enabled");
        }
        if self.per_strand {
            info!(
                "Sampling unit: strand (legacy per-MI sampling; duplex strands sampled independently)"
            );
        } else {
            info!("Sampling unit: molecule (duplex strands /A,/B grouped and sampled together)");
        }
        // downsample is single-threaded and reads through fgumi-bgzf's decoder
        // (`create_raw_bam_reader_with_opts` with threads=1), so it honors
        // --check-crc/--no-check-crc (#800).
        self.io.log_effective_check_crc();

        // Initialize RNG
        let mut rng = match self.seed {
            Some(seed) => StdRng::seed_from_u64(seed),
            None => rand::make_rng(),
        };

        let (reader, header) =
            create_raw_bam_reader_with_opts(&self.io.input, 1, self.io.pipeline_reader_opts())?;

        // Validate header - input must be template-coordinate sorted (output from group)
        if !is_template_coordinate_sorted(&header) {
            bail!(
                "Input BAM must be template-coordinate sorted (output from group).\n\n\
                Expected header fields: SO:unsorted, GO:query, SS:template-coordinate\n\n\
                The input to this tool should be the output of fgumi group or fgbio GroupReadsByUmi."
            );
        }
        info!("Header validation passed (template-coordinate order confirmed)");

        // Add @PG record with PP chaining to input's last program
        let header = crate::commands::common::add_pg_record(header, command_line)?;

        // Create output BAM writer (single-threaded, downsample doesn't have threads parameter)
        let mut writer =
            create_raw_bam_writer(&self.io.output, &header, 1, self.compression.compression_level)?;

        // Create optional rejects writer
        let mut rejects_writer: Option<RawBamWriter> = self
            .rejects
            .as_ref()
            .map(|path| create_raw_bam_writer(path, &header, 1, self.compression.compression_level))
            .transpose()?;

        // Statistics
        let mut total_families: u64 = 0;
        let mut kept_families: u64 = 0;
        let mut kept_reads: u64 = 0;
        let mut rejected_reads: u64 = 0;
        let mut record_count: usize = 0;
        let progress = ProgressTracker::new("Processed records").with_interval(1_000_000);

        // Histograms: family_size -> count
        let mut hist_kept: BTreeMap<usize, u64> = BTreeMap::new();
        let mut hist_rejected: BTreeMap<usize, u64> = BTreeMap::new();

        // For MI order validation
        let mut seen_mis: HashSet<String> = HashSet::new();

        info!("Processing reads...");

        let mut family_iter = FamilyIterator::new(raw_record_iter(reader), self.per_strand);

        while let Some(family_result) = family_iter.next_family()? {
            let (mi, family) = family_result;
            total_families += 1;
            let family_size = family.len();

            // Validate MI order if requested
            if self.validate_mi_order {
                if seen_mis.contains(&mi) {
                    bail!(
                        "MI tag '{mi}' seen non-consecutively. Input BAM may not be properly grouped by MI."
                    );
                }
                seen_mis.insert(mi);
            }

            // Random sampling decision
            let keep = rand::RngExt::random::<f64>(&mut rng) < self.fraction;

            // Count all records processed (kept + rejected)
            record_count += family_size;

            if keep {
                kept_families += 1;
                kept_reads += family_size as u64;
                *hist_kept.entry(family_size).or_insert(0) += 1;

                for record in &family {
                    writer.write_raw_record(record.as_ref())?;
                }
            } else {
                rejected_reads += family_size as u64;
                *hist_rejected.entry(family_size).or_insert(0) += 1;

                if let Some(ref mut rw) = rejects_writer {
                    for record in &family {
                        rw.write_raw_record(record.as_ref())?;
                    }
                }
            }
            progress.log_if_needed(family_size as u64);
        }

        progress.log_final();

        // Write histograms
        if let Some(ref path) = self.histogram_kept {
            write_histogram(&hist_kept, path)?;
            info!("Wrote kept histogram to: {}", path.display());
        }

        if let Some(ref path) = self.histogram_rejected {
            write_histogram(&hist_rejected, path)?;
            info!("Wrote rejected histogram to: {}", path.display());
        }

        // Finalize writers before summary so any I/O failure surfaces here rather
        // than silently on drop (flushes buffered records + writes BGZF EOF).
        writer.finish()?;
        if let Some(rw) = rejects_writer {
            rw.finish()?;
        }

        // Summary
        info!("=== Summary ===");
        info!("Total reads processed: {}", kept_reads + rejected_reads);
        info!("Input families: {total_families}");
        if total_families > 0 {
            let kept_pct = 100.0 * kept_families as f64 / total_families as f64;
            info!("Kept families: {kept_families} ({kept_pct:.2}%)");
        } else {
            info!("Kept families: 0");
        }
        info!("Kept reads: {kept_reads}");
        info!("Rejected reads: {rejected_reads}");
        info!("Output BAM: {}", self.io.output.display());
        if let Some(ref rejects) = self.rejects {
            info!("Rejects BAM: {}", rejects.display());
        }

        timer.log_completion(record_count as u64);
        Ok(())
    }
}

/// Write family size histogram in fgbio-compatible TSV format.
fn write_histogram(histogram: &BTreeMap<usize, u64>, path: &PathBuf) -> Result<()> {
    let mut file = File::create(path)?;
    writeln!(file, "family_size\tcount")?;
    for (size, count) in histogram {
        writeln!(file, "{size}\t{count}")?;
    }
    Ok(())
}

/// Wrap a [`RawBamReader`] as a streaming iterator of [`RawRecord`]s.
///
/// Allocates a fresh `RawRecord` per iteration so downstream code can buffer
/// records by family without stepping on each other's storage.
fn raw_record_iter<R: std::io::Read>(
    mut reader: RawBamReader<R>,
) -> impl Iterator<Item = Result<RawRecord>> {
    let mut exhausted = false;
    std::iter::from_fn(move || {
        if exhausted {
            return None;
        }
        let mut rec = RawRecord::new();
        match reader.read_record(&mut rec) {
            Ok(0) => {
                exhausted = true;
                None
            }
            Ok(_) => Some(Ok(rec)),
            Err(e) => {
                exhausted = true;
                Some(Err(anyhow::Error::from(e)))
            }
        }
    })
}

/// Iterator that groups consecutive records by MI tag.
///
/// This provides streaming iteration over UMI families without loading the entire BAM into memory.
struct FamilyIterator<I>
where
    I: Iterator<Item = Result<RawRecord>>,
{
    records: std::iter::Peekable<I>,
    /// When true (`--per-strand`), group by the raw MI tag (legacy). When false
    /// (default), group by the molecule base MI so a duplex molecule's `<base>/A`
    /// and `<base>/B` strands fall into one family. See `Downsample::per_strand`.
    per_strand: bool,
}

impl<I> FamilyIterator<I>
where
    I: Iterator<Item = Result<RawRecord>>,
{
    fn new(records: I, per_strand: bool) -> Self {
        Self { records: records.peekable(), per_strand }
    }

    /// Get the next family of records sharing the same family key.
    ///
    /// The key is the molecule base MI (any trailing `/A`/`/B` duplex-strand suffix
    /// stripped, so both strands of a molecule fall into one family), or — when
    /// `per_strand` is set — the raw MI tag. Returns `Ok(Some((key, records)))` for
    /// each family, or `Ok(None)` when exhausted.
    fn next_family(&mut self) -> Result<Option<(String, Vec<RawRecord>)>> {
        // Peek at the first record to get the family key
        let key = match self.records.peek() {
            Some(Ok(record)) => family_key(record, self.per_strand)?,
            Some(Err(_)) => {
                // Consume the error — next() is Some because peek() was Some(Err(_))
                return Err(self.records.next().expect("peek() returned Some").unwrap_err());
            }
            None => return Ok(None),
        };

        // Collect all records with the same family key
        let mut family = Vec::new();

        while let Some(peek_result) = self.records.peek() {
            match peek_result {
                Ok(record) => {
                    // Compare the key against the family key without allocating a
                    // `String` per record — the family key was already materialized
                    // once above.
                    if !family_key_equals(record, key.as_bytes(), self.per_strand)? {
                        break;
                    }
                    // Consume the record — next() is Some because peek() was Some
                    family.push(self.records.next().expect("peek() returned Some")?);
                }
                Err(_) => {
                    // Consume and return the error — next() is Some because peek() was Some(Err(_))
                    return Err(self.records.next().expect("peek() returned Some").unwrap_err());
                }
            }
        }

        Ok(Some((key, family)))
    }
}

/// The family key for a record: its molecule base by default, or the raw MI tag
/// when `per_strand` is set. The molecule base is computed by the canonical
/// [`extract_mi_base`], which truncates at the last `/` exactly as `group`,
/// `duplex`, and `duplex_metrics` do when collapsing strands to their source
/// molecule — so grouping here matches the rest of the pipeline rather than a
/// second, divergent strip rule.
fn family_key(record: &RawRecord, per_strand: bool) -> Result<String> {
    let mut mi = get_mi_tag(record)?;
    if !per_strand {
        // `mi` is already valid UTF-8, and `extract_mi_base` truncates at an ASCII
        // '/' boundary, so truncating in place stays on a char boundary and reuses
        // the existing allocation rather than round-tripping through a new String.
        let base_len = extract_mi_base(&mi).len();
        mi.truncate(base_len);
    }
    Ok(mi)
}

/// Whether a record's family key equals `target`, without allocating.
///
/// Mirrors [`family_key`]: compares its molecule base by default, or the raw MI
/// tag when `per_strand` is set. Delegates to [`mi_tag_equals`] in the raw case;
/// in the molecule case it reduces the record's MI to its base via the canonical
/// [`extract_mi_base`] (the Z-typed path) and compares to `target`, which is
/// itself already a base key.
fn family_key_equals(record: &RawRecord, target: &[u8], per_strand: bool) -> Result<bool> {
    if per_strand {
        return mi_tag_equals(record, target);
    }

    let aux = aux_data_slice(record.as_ref());

    if let Some(bytes) = find_string_tag(aux, SamTag::MI) {
        let mi = std::str::from_utf8(bytes)
            .map_err(|e| anyhow::anyhow!("MI tag is not valid UTF-8: {e}"))?;
        return Ok(extract_mi_base(mi).as_bytes() == target);
    }

    if let Some(v) = find_int_tag(aux, SamTag::MI) {
        // Integer MIs never carry a strand suffix, so the base equals the rendering.
        return Ok(v.to_string().as_bytes() == target);
    }

    let name = String::from_utf8_lossy(read_name(record.as_ref())).into_owned();
    let display_name = if name.is_empty() { "<unknown>".to_string() } else { name };
    bail!("Read '{display_name}' is missing required MI tag")
}

/// Extract the MI tag value from a raw BAM record.
///
/// MI is Z-typed per SAM spec, but fgbio historically also writes it as an integer.
/// The raw path supports both: we look for a Z-type first, then fall back to any
/// integer encoding (c/C/s/S/i/I).
fn get_mi_tag(record: &RawRecord) -> Result<String> {
    let aux = aux_data_slice(record.as_ref());

    if let Some(bytes) = find_string_tag(aux, SamTag::MI) {
        return std::str::from_utf8(bytes)
            .map(str::to_string)
            .map_err(|e| anyhow::anyhow!("MI tag is not valid UTF-8: {e}"));
    }

    if let Some(v) = find_int_tag(aux, SamTag::MI) {
        return Ok(v.to_string());
    }

    let name = String::from_utf8_lossy(read_name(record.as_ref())).into_owned();
    let display_name = if name.is_empty() { "<unknown>".to_string() } else { name };
    bail!("Read '{display_name}' is missing required MI tag")
}

/// Whether a record's MI tag equals `target`, without allocating.
///
/// Mirrors [`get_mi_tag`]'s type handling (Z-typed first, then an integer
/// encoding), but compares against `target` in place rather than building a
/// `String` for every record — the family-key `String` is materialized only
/// once per family. The Z path still validates UTF-8 so an invalid MI is an
/// error exactly as in [`get_mi_tag`]; the legacy integer path is rare and falls
/// back to rendering.
fn mi_tag_equals(record: &RawRecord, target: &[u8]) -> Result<bool> {
    let aux = aux_data_slice(record.as_ref());

    if let Some(bytes) = find_string_tag(aux, SamTag::MI) {
        std::str::from_utf8(bytes)
            .map_err(|e| anyhow::anyhow!("MI tag is not valid UTF-8: {e}"))?;
        return Ok(bytes == target);
    }

    if let Some(v) = find_int_tag(aux, SamTag::MI) {
        return Ok(v.to_string().as_bytes() == target);
    }

    let name = String::from_utf8_lossy(read_name(record.as_ref())).into_owned();
    let display_name = if name.is_empty() { "<unknown>".to_string() } else { name };
    bail!("Read '{display_name}' is missing required MI tag")
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::SamBuilder as RawSamBuilder;
    use rstest::rstest;

    #[rstest]
    #[case::typical(0.5, true)]
    #[case::lower_bound_exclusive(0.0, false)]
    #[case::upper_bound_inclusive(1.0, true)]
    #[case::negative(-0.1, false)]
    #[case::above_one(1.5, false)]
    #[case::nan(f64::NAN, false)]
    #[case::infinity(f64::INFINITY, false)]
    #[case::neg_infinity(f64::NEG_INFINITY, false)]
    fn test_validate_fraction(#[case] fraction: f64, #[case] expected_ok: bool) {
        assert_eq!(
            Downsample::validate_fraction(fraction).is_ok(),
            expected_ok,
            "unexpected validation result for fraction {fraction}"
        );
    }

    /// Create a test record with a string MI tag.
    fn create_test_record(name: &str, mi: &str) -> RawRecord {
        let mut b = RawSamBuilder::new();
        b.read_name(name.as_bytes());
        b.add_string_tag(SamTag::MI, mi.as_bytes());
        b.build()
    }

    /// Create a test record with an integer MI tag.
    fn create_test_record_int_mi(name: &str, mi: i32) -> RawRecord {
        let mut b = RawSamBuilder::new();
        b.read_name(name.as_bytes());
        b.add_int_tag(SamTag::MI, mi);
        b.build()
    }

    /// Create a test record with a raw-bytes (Z-typed) MI tag, allowing invalid
    /// UTF-8 that `&str` cannot express.
    fn create_test_record_bytes_mi(name: &str, mi: &[u8]) -> RawRecord {
        let mut b = RawSamBuilder::new();
        b.read_name(name.as_bytes());
        b.add_string_tag(SamTag::MI, mi);
        b.build()
    }

    /// Create a test record without an MI tag.
    fn create_test_record_no_mi(name: &str) -> RawRecord {
        let mut b = RawSamBuilder::new();
        b.read_name(name.as_bytes());
        b.build()
    }

    /// `mi_tag_equals` must agree with a `get_mi_tag` + compare for both Z and
    /// integer MI encodings, and error on a missing MI exactly as `get_mi_tag`
    /// does — it is the allocation-free replacement for that comparison.
    #[test]
    fn test_mi_tag_equals_matches_get_mi_tag() {
        // Z-typed MI: match and mismatch.
        let z = create_test_record("read1", "12345");
        assert!(mi_tag_equals(&z, b"12345").unwrap());
        assert!(!mi_tag_equals(&z, b"12346").unwrap());
        assert!(!mi_tag_equals(&z, b"1234").unwrap()); // prefix is not a match
        assert!(mi_tag_equals(&z, get_mi_tag(&z).unwrap().as_bytes()).unwrap());

        // Integer-typed MI renders the same way get_mi_tag does.
        let i = create_test_record_int_mi("read2", 421);
        assert!(mi_tag_equals(&i, b"421").unwrap());
        assert!(!mi_tag_equals(&i, b"420").unwrap());
        assert!(mi_tag_equals(&i, get_mi_tag(&i).unwrap().as_bytes()).unwrap());

        // A Z-typed MI with invalid UTF-8 is an error, and mi_tag_equals must
        // surface the same UTF-8 error as get_mi_tag rather than reporting a
        // (mis)match against the raw bytes.
        let bad_utf8 = create_test_record_bytes_mi("read3", &[0xff, 0xfe]);
        let equals_err =
            mi_tag_equals(&bad_utf8, &[0xff, 0xfe]).expect_err("invalid UTF-8 MI must be an error");
        let get_err = get_mi_tag(&bad_utf8).expect_err("invalid UTF-8 MI must be an error");
        assert!(equals_err.to_string().contains("not valid UTF-8"));
        assert!(get_err.to_string().contains("not valid UTF-8"));
        assert_eq!(equals_err.to_string(), get_err.to_string());

        // Missing MI is an error, matching get_mi_tag.
        let none = create_test_record_no_mi("read4");
        assert!(mi_tag_equals(&none, b"anything").is_err());
        assert!(get_mi_tag(&none).is_err());
    }

    /// Canonical `BamIoOptions` used by `Downsample` test constructors.
    fn test_bam_io_options() -> BamIoOptions {
        BamIoOptions {
            input: PathBuf::from("input.bam"),
            output: PathBuf::from("output.bam"),
            async_reader: false,
            check_crc: false,
            no_check_crc: false,
        }
    }

    #[test]
    fn test_get_mi_tag_string() {
        let record = create_test_record("read1", "12345");
        let mi = get_mi_tag(&record).expect("get_mi_tag should succeed for string MI");
        assert_eq!(mi, "12345");
    }

    #[test]
    fn test_get_mi_tag_integer() {
        let record = create_test_record_int_mi("read1", 42);
        let mi = get_mi_tag(&record).expect("get_mi_tag should succeed for integer MI");
        assert_eq!(mi, "42");
    }

    #[test]
    fn test_get_mi_tag_missing() {
        let record = create_test_record_no_mi("read1");
        let result = get_mi_tag(&record);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("missing required MI tag"));
    }

    #[test]
    fn test_family_iterator_single_family() {
        let records = vec![
            Ok(create_test_record("r1", "100")),
            Ok(create_test_record("r2", "100")),
            Ok(create_test_record("r3", "100")),
        ];

        let mut iter = FamilyIterator::new(records.into_iter(), false);

        let family1 =
            iter.next_family().expect("next_family should succeed").expect("expected a family");
        assert_eq!(family1.0, "100");
        assert_eq!(family1.1.len(), 3);

        let family2 = iter.next_family().expect("next_family should succeed");
        assert!(family2.is_none());
    }

    #[test]
    fn test_family_iterator_multiple_families() {
        let records = vec![
            Ok(create_test_record("r1", "100")),
            Ok(create_test_record("r2", "100")),
            Ok(create_test_record("r3", "200")),
            Ok(create_test_record("r4", "200")),
            Ok(create_test_record("r5", "200")),
            Ok(create_test_record("r6", "300")),
        ];

        let mut iter = FamilyIterator::new(records.into_iter(), false);

        let family1 =
            iter.next_family().expect("next_family should succeed").expect("expected family 1");
        assert_eq!(family1.0, "100");
        assert_eq!(family1.1.len(), 2);

        let family2 =
            iter.next_family().expect("next_family should succeed").expect("expected family 2");
        assert_eq!(family2.0, "200");
        assert_eq!(family2.1.len(), 3);

        let family3 =
            iter.next_family().expect("next_family should succeed").expect("expected family 3");
        assert_eq!(family3.0, "300");
        assert_eq!(family3.1.len(), 1);

        let family4 = iter.next_family().expect("next_family should succeed");
        assert!(family4.is_none());
    }

    #[test]
    fn test_family_iterator_empty() {
        let records: Vec<Result<RawRecord>> = vec![];
        let mut iter = FamilyIterator::new(records.into_iter(), false);

        let family = iter.next_family().expect("next_family should succeed");
        assert!(family.is_none());
    }

    #[test]
    fn test_write_histogram() {
        use tempfile::NamedTempFile;

        let mut hist = BTreeMap::new();
        hist.insert(1, 10);
        hist.insert(2, 20);
        hist.insert(5, 5);

        let temp_file = NamedTempFile::new().expect("failed to create temp file");
        write_histogram(&hist, &temp_file.path().to_path_buf())
            .expect("write_histogram should succeed");

        let contents =
            std::fs::read_to_string(temp_file.path()).expect("failed to read histogram file");
        assert!(contents.contains("family_size\tcount"));
        assert!(contents.contains("1\t10"));
        assert!(contents.contains("2\t20"));
        assert!(contents.contains("5\t5"));
    }

    #[test]
    #[allow(clippy::float_cmp)] // Testing exact value assignment, not computation
    fn test_downsample_parameters() {
        let cmd = Downsample {
            io: test_bam_io_options(),
            fraction: 0.1,
            rejects: Some(PathBuf::from("rejects.bam")),
            seed: Some(42),
            validate_mi_order: true,
            per_strand: true,
            histogram_kept: Some(PathBuf::from("kept.txt")),
            histogram_rejected: Some(PathBuf::from("rejected.txt")),
            compression: CompressionOptions { compression_level: 1 },
        };

        assert_eq!(cmd.fraction, 0.1);
        assert_eq!(cmd.seed, Some(42));
        assert!(cmd.validate_mi_order);
        assert!(cmd.per_strand);
        assert!(cmd.rejects.is_some());
        assert!(cmd.histogram_kept.is_some());
        assert!(cmd.histogram_rejected.is_some());
    }

    // Note: Header validation tests are in crate::sam::tests for is_template_coordinate_sorted

    #[test]
    fn test_deterministic_sampling_with_seed() {
        // Test that the same seed produces the same results
        use rand::RngExt;

        let seed = 12345u64;

        // First run
        let mut rng1 = StdRng::seed_from_u64(seed);
        let results1: Vec<bool> = (0..100).map(|_| rng1.random::<f64>() < 0.5).collect();

        // Second run with same seed
        let mut rng2 = StdRng::seed_from_u64(seed);
        let results2: Vec<bool> = (0..100).map(|_| rng2.random::<f64>() < 0.5).collect();

        assert_eq!(results1, results2);
    }

    #[test]
    fn test_histogram_sorted_by_family_size() {
        let mut hist = BTreeMap::new();
        hist.insert(5, 10);
        hist.insert(1, 20);
        hist.insert(3, 15);

        // BTreeMap maintains sorted order
        let sizes: Vec<usize> = hist.keys().copied().collect();
        assert_eq!(sizes, vec![1, 3, 5]);
    }

    #[test]
    fn test_family_key_respects_per_strand() {
        let a = create_test_record("r", "7/A");
        // Default (molecule) strips the strand suffix; --per-strand keeps the raw tag.
        assert_eq!(family_key(&a, false).unwrap(), "7");
        assert_eq!(family_key(&a, true).unwrap(), "7/A");

        // An integer MI has no `/`, so both modes agree.
        let i = create_test_record_int_mi("r", 7);
        assert_eq!(family_key(&i, false).unwrap(), "7");
        assert_eq!(family_key(&i, true).unwrap(), "7");
    }

    /// The default molecule key uses the canonical `extract_mi_base` rule (strip at
    /// the last `/`), NOT a `/A`,`/B`-only strip: any final `/`-suffix is removed, and
    /// a leading `/` (empty base) is preserved. Pinned deliberately so a future change
    /// cannot silently narrow the rule to `/A`,`/B` and diverge downsample's grouping
    /// from `group`/`duplex`, which collapse strands with this same canonical rule.
    #[test]
    fn test_family_key_default_strips_any_suffix_canonically() {
        // A non-/A,/B suffix is still stripped (matches group/duplex, not narrowed).
        assert_eq!(family_key(&create_test_record("r", "7/C"), false).unwrap(), "7");
        // A leading '/' is an empty base and is preserved verbatim (extract_mi_base).
        assert_eq!(family_key(&create_test_record("r", "/A"), false).unwrap(), "/A");
    }

    #[test]
    fn test_family_key_equals_per_strand() {
        let a = create_test_record("r", "7/A");
        let b = create_test_record("r", "7/B");

        // Default (molecule): both strands match the base key "7".
        assert!(family_key_equals(&a, b"7", false).unwrap());
        assert!(family_key_equals(&b, b"7", false).unwrap());
        assert!(!family_key_equals(&a, b"8", false).unwrap());

        // --per-strand: /A and /B are different families.
        assert!(family_key_equals(&a, b"7/A", true).unwrap());
        assert!(!family_key_equals(&b, b"7/A", true).unwrap());

        // A missing MI is an error in both modes (mirrors get_mi_tag).
        let none = create_test_record_no_mi("r");
        assert!(family_key_equals(&none, b"7", false).is_err());
        assert!(family_key_equals(&none, b"7", true).is_err());
    }

    /// `family_key_equals` must agree with `family_key` + compare, in BOTH modes and
    /// across Z-typed and integer MIs — they are the two sides of the same grouping
    /// decision and drift between them would split or merge families incorrectly.
    #[test]
    fn test_family_key_equals_matches_family_key() {
        let records = [
            create_test_record("r", "7/A"),
            create_test_record("r", "7/B"),
            create_test_record("r", "7"), // simplex MI, no suffix
            create_test_record_int_mi("r", 7),
        ];
        for per_strand in [false, true] {
            for rec in &records {
                let key = family_key(rec, per_strand).unwrap();
                assert!(
                    family_key_equals(rec, key.as_bytes(), per_strand).unwrap(),
                    "family_key_equals disagreed with family_key (per_strand={per_strand})"
                );
            }
        }
    }

    #[test]
    fn test_family_iterator_duplex_strands_grouped_by_default_and_split_by_per_strand() {
        // One molecule's two strands (7/A, 7/B), then a second molecule (8/A).
        // Built fresh per iterator: Vec<Result<_, anyhow::Error>> is not Clone.
        let records = || {
            vec![
                Ok(create_test_record("r1", "7/A")),
                Ok(create_test_record("r2", "7/A")),
                Ok(create_test_record("r3", "7/B")),
                Ok(create_test_record("r4", "8/A")),
            ]
        };

        // Default (per_strand = false): 7/A + 7/B collapse into ONE family keyed "7".
        let mut iter = FamilyIterator::new(records().into_iter(), false);
        let fam1 = iter.next_family().unwrap().expect("family 1");
        assert_eq!(fam1.0, "7");
        assert_eq!(fam1.1.len(), 3); // both strands of molecule 7
        let fam2 = iter.next_family().unwrap().expect("family 2");
        assert_eq!(fam2.0, "8");
        assert_eq!(fam2.1.len(), 1);
        assert!(iter.next_family().unwrap().is_none());

        // --per-strand (true): 7/A and 7/B are separate families.
        let mut iter = FamilyIterator::new(records().into_iter(), true);
        let f1 = iter.next_family().unwrap().expect("family 7/A");
        assert_eq!((f1.0.as_str(), f1.1.len()), ("7/A", 2));
        let f2 = iter.next_family().unwrap().expect("family 7/B");
        assert_eq!((f2.0.as_str(), f2.1.len()), ("7/B", 1));
        let f3 = iter.next_family().unwrap().expect("family 8/A");
        assert_eq!((f3.0.as_str(), f3.1.len()), ("8/A", 1));
        assert!(iter.next_family().unwrap().is_none());
    }
}
