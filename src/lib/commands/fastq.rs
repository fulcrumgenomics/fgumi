//! Convert BAM to FASTQ format.
//!
//! This tool reads a BAM file and outputs interleaved FASTQ to stdout for piping to aligners.
//! Input should be queryname-sorted or template-coordinate sorted.

use crate::commands::common::parse_bool;
use crate::logging::OperationTimer;
use crate::sam::SamTag;
use crate::validation::validate_file_exists;
use anyhow::Result;
use clap::Parser;
use fgumi_bam_io::{BgzfWriterEnum, create_bgzf_writer, create_raw_bam_reader};
use fgumi_raw_bam::{
    RawRecord, extract_sequence_into, quality_scores_slice, read_name as raw_read_name,
};
use log::info;
use std::fs::File;
use std::io::{BufWriter, StdoutLock, Write, stdout};
use std::path::{Path, PathBuf};

use crate::commands::command::Command;

/// Lookup table for Phred to Phred+33 ASCII conversion (clamped to 126)
static QUAL_TO_ASCII: [u8; 256] = {
    let mut table = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        let val = (i as u8).saturating_add(33);
        table[i] = if val > 126 { 126 } else { val };
        i += 1;
    }
    table
};

/// Lookup table for base complement (A<->T, C<->G, others->N)
static COMPLEMENT: [u8; 256] = {
    let mut table = [b'N'; 256];
    table[b'A' as usize] = b'T';
    table[b'a' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b't' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'c' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'g' as usize] = b'C';
    table[b'N' as usize] = b'N';
    table[b'n' as usize] = b'N';
    table
};

/// Parse a comma-separated list of two-character SAM tags (e.g. `"RX,OX"`).
///
/// Whitespace around each tag is trimmed. Returns an error if any token is not
/// a valid SAM aux-tag identifier, which also rejects empty tokens (`""`,
/// `"RX,"`).
fn parse_umi_tags(spec: &str) -> Result<Vec<SamTag>> {
    spec.split(',').map(|token| token.trim().parse::<SamTag>()).collect()
}

/// Configuration for embedding a UMI into the FASTQ read header, mirroring
/// `samtools fastq -U`.
///
/// The UMI is read from the first present tag in `tags` (default `RX` then
/// `OX`) and appended to the read name as `<name><name_delim><umi>`, with each
/// `-` in the stored UMI (fgumi's duplex join) replaced by `umi_sep` (default
/// `+`, the delimiter Illumina DRAGEN expects between duplex UMIs).
pub(crate) struct UmiHeader {
    /// Tags to search, in priority order; the first present wins.
    tags: Vec<SamTag>,
    /// Bytes inserted between the read name and the UMI.
    name_delim: Vec<u8>,
    /// Bytes that replace each `-` (duplex join) in the stored UMI.
    umi_sep: Vec<u8>,
}

impl UmiHeader {
    /// Build from raw CLI strings, validating the tag list.
    pub(crate) fn from_args(umi_tag: &str, name_delim: &str, umi_sep: &str) -> Result<Self> {
        Ok(Self {
            tags: parse_umi_tags(umi_tag)?,
            name_delim: name_delim.as_bytes().to_vec(),
            umi_sep: umi_sep.as_bytes().to_vec(),
        })
    }

    /// Look up the UMI value from the record's tags (first present, in priority order).
    fn lookup<'a>(&self, record: &'a RawRecord) -> Option<&'a [u8]> {
        let tags = record.tags();
        self.tags.iter().find_map(|tag| tags.find_string(*tag))
    }

    /// Write `<name_delim><umi>` to `writer`, replacing each `-` in `umi` with
    /// `umi_sep`. Maximal runs between separators are written in one call.
    fn write_annotation<W: Write>(&self, writer: &mut W, umi: &[u8]) -> std::io::Result<()> {
        writer.write_all(&self.name_delim)?;
        let mut start = 0;
        for (i, &b) in umi.iter().enumerate() {
            if b == b'-' {
                writer.write_all(&umi[start..i])?;
                writer.write_all(&self.umi_sep)?;
                start = i + 1;
            }
        }
        writer.write_all(&umi[start..])?;
        Ok(())
    }
}

/// Which FASTQ stream a record belongs to, based on its segment flags.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Segment {
    /// First segment of a pair (`FIRST_SEGMENT` set, `LAST_SEGMENT` clear) → R1.
    Read1,
    /// Last segment of a pair (`LAST_SEGMENT` set, `FIRST_SEGMENT` clear) → R2.
    Read2,
    /// Neither or both segment bits set (single-end or ambiguous) → "other".
    Other,
}

/// Classify a record into R1 / R2 / other from its flags, matching how
/// `samtools fastq` routes reads to `-1` / `-2` / `-0`.
fn classify_segment(flags: u16) -> Segment {
    use fgumi_raw_bam::flags as flag_bits;
    let is_first = (flags & flag_bits::FIRST_SEGMENT) != 0;
    let is_last = (flags & flag_bits::LAST_SEGMENT) != 0;
    match (is_first, is_last) {
        (true, false) => Segment::Read1,
        (false, true) => Segment::Read2,
        _ => Segment::Other,
    }
}

/// Returns `true` if `path` should be written compressed (ends in `.gz`/`.bgz`).
fn path_is_gzip(path: &Path) -> bool {
    matches!(path.extension().and_then(|e| e.to_str()), Some("gz" | "bgz" | "bgzf"))
}

/// Default BGZF/DEFLATE compression level for `.gz` FASTQ output (matches the
/// samtools default).
const FASTQ_GZIP_LEVEL: u32 = 6;

/// A FASTQ output stream. Compressed file outputs use BGZF (block gzip) via the
/// same multi-threaded backend fgumi uses for BAM — gzip-compatible for
/// downstream tools yet block-indexable. Plain files and stdout are written
/// uncompressed.
enum FastqSink {
    /// Uncompressed file.
    Plain(BufWriter<File>),
    /// BGZF-compressed file (chosen for `.gz`/`.bgz` paths).
    Bgzf(BgzfWriterEnum),
    /// Uncompressed stdout (for piping to an aligner).
    Stdout(BufWriter<StdoutLock<'static>>),
}

impl FastqSink {
    /// Output buffer capacity, sized for efficient pipe/file throughput.
    const BUF_CAPACITY: usize = 64 * 1024 * 1024;

    /// Open a file sink, writing BGZF when the path ends in `.gz`/`.bgz`. BGZF
    /// compression is parallelized across `threads` worker threads.
    fn open_file(path: &Path, threads: usize) -> Result<Self> {
        let buf = BufWriter::with_capacity(Self::BUF_CAPACITY, File::create(path)?);
        if path_is_gzip(path) {
            let boxed: Box<dyn Write + Send> = Box::new(buf);
            Ok(FastqSink::Bgzf(create_bgzf_writer(boxed, threads, FASTQ_GZIP_LEVEL)))
        } else {
            Ok(FastqSink::Plain(buf))
        }
    }

    /// A sink writing uncompressed FASTQ to stdout.
    fn stdout() -> Self {
        FastqSink::Stdout(BufWriter::with_capacity(Self::BUF_CAPACITY, stdout().lock()))
    }

    /// Finish the stream: for BGZF this flushes all workers and writes the EOF
    /// block; all variants flush their underlying writer.
    fn finish(self) -> std::io::Result<()> {
        match self {
            FastqSink::Plain(mut w) => w.flush(),
            FastqSink::Stdout(mut w) => w.flush(),
            FastqSink::Bgzf(w) => w.finish(),
        }
    }
}

impl Write for FastqSink {
    #[inline]
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            FastqSink::Plain(w) => w.write(buf),
            FastqSink::Stdout(w) => w.write(buf),
            FastqSink::Bgzf(w) => w.write(buf),
        }
    }

    #[inline]
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            FastqSink::Plain(w) => w.flush(),
            FastqSink::Stdout(w) => w.flush(),
            FastqSink::Bgzf(w) => w.flush(),
        }
    }
}

/// Convert BAM to FASTQ format.
#[derive(Debug, Parser)]
#[command(
    name = "fastq",
    about = "\x1b[38;5;72m[ALIGNMENT]\x1b[0m      \x1b[36mConvert BAM to FASTQ format\x1b[0m",
    long_about = r#"
Convert a BAM file to FASTQ.

By default, reads BAM records and writes interleaved FASTQ to stdout for piping
to aligners. Can also write paired split files (--out1/--out2, plus --out0),
gzip (BGZF) output for any .gz/.bgz path, and embed the UMI in the read name
(--annotate-read-names). Input should be queryname-sorted or template-coordinate
sorted.

EXAMPLES:

  # Pipe to bwa mem for alignment
  fgumi fastq -i unmapped.bam | bwa mem -t 16 -p -K 150000000 -Y ref.fa -

  # With multi-threaded BAM decompression
  fgumi fastq -i unmapped.bam -@ 4 | bwa mem -t 16 -p ref.fa -

  # Exclude secondary and supplementary alignments (default)
  fgumi fastq -i aligned.bam -F 0x900 | bwa mem ...

  # Embed the UMI (RX tag) in the read name for DRAGEN, like `samtools fastq -U`.
  # A duplex UMI stored as `AAA-CCC` is emitted as `readname:AAA+CCC`.
  fgumi fastq -i extracted.bam --annotate-read-names -n > umi_in_name.fq

  # Custom delimiters for other platforms (underscore separator, no duplex join)
  fgumi fastq -i extracted.bam --annotate-read-names --umi-name-delim _ --umi-sep '' ...
"#
)]
pub struct Fastq {
    /// Input BAM file.
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Output FASTQ file. If omitted, the FASTQ stream is written to stdout
    /// (the default, intended for piping straight to an aligner).
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Don't append /1 and /2 to read names.
    #[arg(short = 'n', long = "no-read-suffix", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub no_suffix: bool,

    /// Exclude reads with any of these flags present [0x900 = secondary|supplementary].
    #[arg(short = 'F', long = "exclude-flags", default_value_t = 0x900, value_parser = parse_flags)]
    pub exclude_flags: u16,

    /// Only include reads with all of these flags present.
    #[arg(short = 'f', long = "require-flags", default_value_t = 0, value_parser = parse_flags)]
    pub require_flags: u16,

    /// Number of threads for BAM decompression and, for `.gz`/`.bgz` outputs,
    /// BGZF compression (per output). Use several threads for compressed output.
    #[arg(short = '@', short_alias = 't', long = "threads", default_value = "1")]
    pub threads: usize,

    /// BWA -K parameter value (bases per batch). Sizes output buffer to match bwa's
    /// batch size for optimal pipe throughput. Default matches common bwa mem usage.
    #[arg(short = 'K', long = "bwa-chunk-size", default_value = "150000000")]
    pub bwa_chunk_size: u64,

    /// Append the UMI (from a BAM tag) to the read name, matching
    /// `fgumi extract --annotate-read-names` and `samtools fastq -U`. The UMI is
    /// inserted before any /1 or /2 suffix.
    #[arg(short = 'a', visible_short_alias = 'U', long = "annotate-read-names", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = parse_bool)]
    pub annotate_read_names: bool,

    /// Comma-separated tag(s) to source the UMI from; the first present wins.
    #[arg(long = "umi-tag", default_value = "RX,OX")]
    pub umi_tag: String,

    /// Delimiter inserted between the read name and the UMI.
    #[arg(long = "umi-name-delim", default_value = ":")]
    pub umi_name_delim: String,

    /// Delimiter joining multiple (duplex) UMIs, replacing the stored '-'
    /// [default '+' matches Illumina DRAGEN].
    #[arg(long = "umi-sep", default_value = "+")]
    pub umi_sep: String,

    /// Write read 1 (R1) to this file instead of the interleaved stream; requires
    /// --out2. A `.gz`/`.bgz` path is written as BGZF.
    #[arg(short = '1', long = "out1", requires = "out2", conflicts_with = "output")]
    pub out1: Option<PathBuf>,

    /// Write read 2 (R2) to this file; requires --out1.
    /// A `.gz`/`.bgz` path is written as BGZF.
    #[arg(short = '2', long = "out2", requires = "out1", conflicts_with = "output")]
    pub out2: Option<PathBuf>,

    /// Write reads that are neither cleanly R1 nor R2 (single-end / ambiguous)
    /// here; requires --out1. If omitted, such reads go to stdout, matching
    /// `samtools fastq` without -0. A `.gz`/`.bgz` path is written as BGZF.
    #[arg(short = '0', long = "out0", requires = "out1")]
    pub out0: Option<PathBuf>,
}

/// Parse flag values supporting both decimal and hex (0x) notation.
fn parse_flags(s: &str) -> Result<u16, String> {
    if s.starts_with("0x") || s.starts_with("0X") {
        u16::from_str_radix(&s[2..], 16).map_err(|e| e.to_string())
    } else {
        s.parse().map_err(|e: std::num::ParseIntError| e.to_string())
    }
}

impl Fastq {
    /// Build the optional UMI-in-read-name config, validating the tag list up front.
    fn build_umi_header(&self) -> Result<Option<UmiHeader>> {
        if self.annotate_read_names {
            info!("UMI in read name: from tag(s) {}", self.umi_tag);
            Ok(Some(UmiHeader::from_args(&self.umi_tag, &self.umi_name_delim, &self.umi_sep)?))
        } else {
            Ok(None)
        }
    }

    /// Log the shared conversion configuration.
    fn log_config(&self) {
        info!("Input: {}", self.input.display());
        info!("Threads: {}", self.threads);
        info!("Exclude flags: 0x{:X}", self.exclude_flags);
        info!("Require flags: 0x{:X}", self.require_flags);
        info!("Read name suffix: {}", if self.no_suffix { "disabled" } else { "enabled" });
    }

    /// Returns `true` if the record passes the include/exclude flag filters.
    #[inline]
    fn passes_filters(&self, flags: u16) -> bool {
        (flags & self.exclude_flags) == 0 && (flags & self.require_flags) == self.require_flags
    }

    /// Interleaved conversion: every record is written to a single sink (stdout
    /// for piping, or one `--output` file).
    fn run_interleaved(&self, mut sink: FastqSink) -> Result<()> {
        let timer = OperationTimer::new("Converting BAM to FASTQ");
        self.log_config();
        info!("BWA chunk size: {} bases", self.bwa_chunk_size);
        let umi_header = self.build_umi_header()?;

        let (mut reader, _header) = create_raw_bam_reader(&self.input, self.threads)?;

        let mut total_records: u64 = 0;
        let mut written_records: u64 = 0;

        // Batch tracking (mirrors bwa's logic)
        let mut bases_this_batch: u64 = 0;
        let mut records_this_batch: usize = 0;

        // Reusable buffers to avoid per-record allocations
        let mut seq_buf: Vec<u8> = Vec::with_capacity(512);
        let mut qual_buf: Vec<u8> = Vec::with_capacity(512);
        let mut record = RawRecord::new();

        loop {
            let n = reader.read_record(&mut record)?;
            if n == 0 {
                break; // EOF
            }
            total_records += 1;

            let flags = record.flags();
            if !self.passes_filters(flags) {
                continue;
            }

            // Get sequence length for batch tracking
            let seq_len = record.l_seq() as usize;

            write_fastq_record(
                &mut sink,
                &record,
                flags,
                self.no_suffix,
                umi_header.as_ref(),
                &mut seq_buf,
                &mut qual_buf,
            )?;
            written_records += 1;

            // Track batch progress
            bases_this_batch += seq_len as u64;
            records_this_batch += 1;

            // Flush at batch boundary (like bwa: bases >= K AND record count is even)
            if bases_this_batch >= self.bwa_chunk_size && records_this_batch.is_multiple_of(2) {
                sink.flush()?;
                bases_this_batch = 0;
                records_this_batch = 0;
            }
        }

        sink.finish()?;

        info!("Read {total_records} records, wrote {written_records} FASTQ records");
        timer.log_completion(written_records);
        Ok(())
    }

    /// Paired conversion: route R1/R2/other reads to the `--out1`/`--out2`/`--out0`
    /// sinks by segment flag, mirroring `samtools fastq -1/-2/-0`. Assumes input
    /// is name-grouped (the standard queryname / template-coordinate order);
    /// a singleton mate follows its own segment flag, as `samtools fastq` does
    /// without `-s`.
    fn run_paired(&self) -> Result<()> {
        let timer = OperationTimer::new("Converting BAM to paired FASTQ");
        self.log_config();
        let umi_header = self.build_umi_header()?;

        // clap `requires` guarantees out2 is present whenever out1 is.
        let out1 = self.out1.as_ref().expect("out1 present in paired mode");
        let out2 = self.out2.as_ref().expect("out2 present in paired mode");
        info!("R1 -> {}", out1.display());
        info!("R2 -> {}", out2.display());

        // Each compressed (BGZF) output gets its own pool of `--threads`
        // compression workers. Surface the effective total so a large `--threads`
        // with several gzipped outputs doesn't silently oversubscribe the CPU
        // (e.g. `--threads 16` with three `.gz` outputs = 48 compression workers,
        // plus the reader's own decompression threads).
        let threads = self.threads;
        let gzip_output_count = [Some(out1), Some(out2), self.out0.as_ref()]
            .into_iter()
            .flatten()
            .filter(|p| path_is_gzip(p))
            .count();
        if threads > 1 && gzip_output_count > 1 {
            info!(
                "Compressing {} BGZF outputs with {} workers each ({} compression workers total)",
                gzip_output_count,
                threads,
                threads * gzip_output_count
            );
        }
        let mut sink1 = FastqSink::open_file(out1, threads)?;
        let mut sink2 = FastqSink::open_file(out2, threads)?;
        // Reads that are neither cleanly R1 nor R2 go to --out0 if given, else to
        // stdout (matching `samtools fastq` without -0).
        let mut sink_other = match &self.out0 {
            Some(path) => {
                info!("other -> {}", path.display());
                FastqSink::open_file(path, threads)?
            }
            None => FastqSink::stdout(),
        };

        let (mut reader, _header) = create_raw_bam_reader(&self.input, self.threads)?;

        let mut total_records: u64 = 0;
        let mut written_records: u64 = 0;
        let mut seq_buf: Vec<u8> = Vec::with_capacity(512);
        let mut qual_buf: Vec<u8> = Vec::with_capacity(512);
        let mut record = RawRecord::new();

        // In split mode the R1/R2 file already identifies the mate, so the
        // `/1` `/2` suffix is always omitted — matching `samtools fastq -1/-2`,
        // which keeps R1 and R2 read names identical for downstream pairing.
        let no_suffix = true;

        loop {
            let n = reader.read_record(&mut record)?;
            if n == 0 {
                break; // EOF
            }
            total_records += 1;

            let flags = record.flags();
            if !self.passes_filters(flags) {
                continue;
            }

            let sink = match classify_segment(flags) {
                Segment::Read1 => &mut sink1,
                Segment::Read2 => &mut sink2,
                Segment::Other => &mut sink_other,
            };
            write_fastq_record(
                sink,
                &record,
                flags,
                no_suffix,
                umi_header.as_ref(),
                &mut seq_buf,
                &mut qual_buf,
            )?;
            written_records += 1;
        }

        sink1.finish()?;
        sink2.finish()?;
        sink_other.finish()?;

        info!("Read {total_records} records, wrote {written_records} FASTQ records");
        timer.log_completion(written_records);
        Ok(())
    }
}

impl Command for Fastq {
    fn execute(&self, _command_line: &str) -> Result<()> {
        validate_file_exists(&self.input, "Input BAM")?;

        // Refuse to clobber the input BAM or write two streams to the same file.
        // `File::create(path)` truncates before the input is opened, so an output
        // path pointing at the input silently destroys it, and two outputs sharing
        // a path silently corrupt each other (two writer threads appending to one
        // truncated file). Compare by canonical path where the file exists (which
        // also catches `./x` vs `x`, symlinks, and absolute vs relative). For
        // outputs that don't exist yet, canonicalize the (existing) parent
        // directory and re-join the file name so aliases that only differ through
        // a symlinked intermediate directory (or a leading `./`) still collide in
        // the `seen` dedup — otherwise two sink threads could truncate the same
        // not-yet-existing file.
        let output_paths: Vec<&PathBuf> =
            [&self.output, &self.out1, &self.out2, &self.out0].into_iter().flatten().collect();
        // Nested `if let` (not a let-chain) to stay compatible with the 1.87.0 MSRV.
        let canonical = |p: &Path| -> Result<PathBuf> {
            if p.exists() {
                return Ok(std::fs::canonicalize(p)?);
            }
            if let Some(parent) = p.parent().filter(|d| !d.as_os_str().is_empty()) {
                if let Ok(canon_parent) = std::fs::canonicalize(parent) {
                    if let Some(file_name) = p.file_name() {
                        return Ok(canon_parent.join(file_name));
                    }
                }
            }
            Ok(p.to_path_buf())
        };
        let input_key = canonical(&self.input)?;
        let mut seen: std::collections::HashSet<PathBuf> = std::collections::HashSet::new();
        for output in &output_paths {
            // The literal `/dev/null` is exempt — routing several streams to it
            // is intentional, as with `samtools fastq` (e.g. `--out0 /dev/null`).
            // Only the exact path is exempt, not symlinks to it or `/dev/stdout`.
            if output.as_os_str() == "/dev/null" {
                continue;
            }
            let key = canonical(output)?;
            if key == input_key {
                anyhow::bail!(
                    "output path {} must differ from --input {} (would truncate the input BAM)",
                    output.display(),
                    self.input.display()
                );
            }
            if !seen.insert(key) {
                anyhow::bail!(
                    "output path {} is used more than once; each FASTQ output must be a distinct file",
                    output.display()
                );
            }
        }

        // Paired split output (`-1`/`-2`, optionally `-0`). clap guarantees
        // `--out1`/`--out2` come together and conflict with `--output`.
        if self.out1.is_some() {
            return self.run_paired();
        }

        let sink = match &self.output {
            Some(path) => FastqSink::open_file(path, self.threads)?,
            None => FastqSink::stdout(),
        };
        self.run_interleaved(sink)
    }
}

/// Write a single FASTQ record to the writer.
///
/// `pub(crate)` so the runall AAM chain (`crate::align_and_merge`)
/// can reuse the same BAM→FASTQ logic the standalone `fgumi fastq`
/// command uses, instead of duplicating the encoding.
///
/// When `umi_header` is `Some`, the UMI from the configured tag is appended to
/// the read name (before any `/1`/`/2` suffix), mirroring `samtools fastq -U`.
#[inline]
pub(crate) fn write_fastq_record<W: Write>(
    writer: &mut W,
    record: &RawRecord,
    flags: u16,
    no_suffix: bool,
    umi_header: Option<&UmiHeader>,
    seq_buf: &mut Vec<u8>,
    qual_buf: &mut Vec<u8>,
) -> Result<()> {
    use fgumi_raw_bam::flags as flag_bits;

    // Get read name (without null terminator)
    let name = raw_read_name(record);

    // Determine read suffix (/1 or /2)
    let is_first = (flags & flag_bits::FIRST_SEGMENT) != 0;
    let is_last = (flags & flag_bits::LAST_SEGMENT) != 0;
    let suffix: &[u8] = if no_suffix {
        b""
    } else if is_first && !is_last {
        b"/1"
    } else if is_last && !is_first {
        b"/2"
    } else {
        b"" // Single-end or both flags set
    };

    // Decode sequence from 4-bit BAM encoding to ASCII bases (reuses seq_buf).
    extract_sequence_into(record, seq_buf);

    // Copy quality bytes and transform to Phred+33 ASCII
    qual_buf.clear();
    qual_buf.extend(quality_scores_slice(record).iter().map(|&s| QUAL_TO_ASCII[s as usize]));

    if (flags & flag_bits::REVERSE) != 0 {
        // Reverse complement sequence in place using lookup table
        seq_buf.reverse();
        for base in seq_buf.iter_mut() {
            *base = COMPLEMENT[*base as usize];
        }
        // Reverse quality in place
        qual_buf.reverse();
    }

    // Write all parts
    writer.write_all(b"@")?;
    writer.write_all(name)?;
    // Append the UMI to the read name (before the /1 /2 suffix), if requested.
    // Nested `if let` (not a let-chain) to stay compatible with the 1.87.0 MSRV.
    if let Some(uh) = umi_header {
        if let Some(umi) = uh.lookup(record) {
            uh.write_annotation(writer, umi)?;
        }
    }
    writer.write_all(suffix)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq_buf)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(qual_buf)?;
    writer.write_all(b"\n")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Write reverse complement of sequence bytes to a buffer (test helper).
    fn write_reverse_complement_bytes<W: Write>(writer: &mut W, bases: &[u8]) -> Result<()> {
        for &base in bases.iter().rev() {
            let comp = match base {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                b'N' | b'n' => b'N',
                _ => b'N',
            };
            writer.write_all(&[comp])?;
        }
        Ok(())
    }

    /// Write quality scores as Phred+33 ASCII (test helper).
    fn write_quality_bytes<W: Write>(writer: &mut W, quals: &[u8]) -> Result<()> {
        for &score in quals {
            let ascii = score.saturating_add(33).min(126);
            writer.write_all(&[ascii])?;
        }
        Ok(())
    }

    /// Write reversed quality scores as Phred+33 ASCII (test helper).
    fn write_reversed_quality_bytes<W: Write>(writer: &mut W, quals: &[u8]) -> Result<()> {
        for &score in quals.iter().rev() {
            let ascii = score.saturating_add(33).min(126);
            writer.write_all(&[ascii])?;
        }
        Ok(())
    }

    #[test]
    fn test_parse_flags_decimal() {
        assert_eq!(parse_flags("2304").expect("parse decimal '2304' should succeed"), 2304);
        assert_eq!(parse_flags("0").expect("parse decimal '0' should succeed"), 0);
        assert_eq!(parse_flags("65535").expect("parse decimal '65535' should succeed"), 65535);
    }

    #[test]
    fn test_parse_flags_hex() {
        assert_eq!(parse_flags("0x900").expect("parse hex '0x900' should succeed"), 0x900);
        assert_eq!(parse_flags("0X900").expect("parse hex '0X900' should succeed"), 0x900);
        assert_eq!(parse_flags("0xff").expect("parse hex '0xff' should succeed"), 0xff);
        assert_eq!(parse_flags("0xFFFF").expect("parse hex '0xFFFF' should succeed"), 0xFFFF);
    }

    #[test]
    fn test_parse_flags_invalid() {
        assert!(parse_flags("invalid").is_err());
        assert!(parse_flags("0xGGGG").is_err());
        assert!(parse_flags("-1").is_err());
    }

    #[test]
    fn test_qual_to_ascii_lookup_table() {
        // Quality 0 -> ASCII 33 ('!')
        assert_eq!(QUAL_TO_ASCII[0], 33);
        // Quality 30 -> ASCII 63 ('?')
        assert_eq!(QUAL_TO_ASCII[30], 63);
        // Quality 40 -> ASCII 73 ('I')
        assert_eq!(QUAL_TO_ASCII[40], 73);
        // Quality 93 -> ASCII 126 ('~') - max valid
        assert_eq!(QUAL_TO_ASCII[93], 126);
        // Quality > 93 should be clamped to 126
        assert_eq!(QUAL_TO_ASCII[94], 126);
        assert_eq!(QUAL_TO_ASCII[255], 126);
    }

    #[test]
    fn test_complement_lookup_table() {
        // Standard bases
        assert_eq!(COMPLEMENT[b'A' as usize], b'T');
        assert_eq!(COMPLEMENT[b'T' as usize], b'A');
        assert_eq!(COMPLEMENT[b'C' as usize], b'G');
        assert_eq!(COMPLEMENT[b'G' as usize], b'C');
        assert_eq!(COMPLEMENT[b'N' as usize], b'N');

        // Lowercase
        assert_eq!(COMPLEMENT[b'a' as usize], b'T');
        assert_eq!(COMPLEMENT[b't' as usize], b'A');
        assert_eq!(COMPLEMENT[b'c' as usize], b'G');
        assert_eq!(COMPLEMENT[b'g' as usize], b'C');
        assert_eq!(COMPLEMENT[b'n' as usize], b'N');

        // Unknown bases map to N
        assert_eq!(COMPLEMENT[b'X' as usize], b'N');
        assert_eq!(COMPLEMENT[0], b'N');
    }

    #[test]
    fn test_write_reverse_complement_bytes() {
        let mut output = Vec::new();
        // ACGT reversed = TGCA, complemented = ACGT
        write_reverse_complement_bytes(&mut output, b"ACGT")
            .expect("write_reverse_complement_bytes should succeed");
        assert_eq!(output, b"ACGT");

        output.clear();
        write_reverse_complement_bytes(&mut output, b"AAAA")
            .expect("write_reverse_complement_bytes should succeed");
        assert_eq!(output, b"TTTT");

        output.clear();
        // ATCG reversed = GCTA, complemented = CGAT
        write_reverse_complement_bytes(&mut output, b"ATCG")
            .expect("write_reverse_complement_bytes should succeed");
        assert_eq!(output, b"CGAT");

        output.clear();
        // Test with N
        write_reverse_complement_bytes(&mut output, b"ANCG")
            .expect("write_reverse_complement_bytes should succeed");
        assert_eq!(output, b"CGNT");
    }

    #[test]
    fn test_write_quality_bytes() {
        let mut output = Vec::new();
        // Quality 0 -> ASCII 33 ('!')
        // Quality 30 -> ASCII 63 ('?')
        write_quality_bytes(&mut output, &[0, 30, 40]).expect("write_quality_bytes should succeed");
        assert_eq!(output, vec![33, 63, 73]);
    }

    #[test]
    fn test_write_reversed_quality_bytes() {
        let mut output = Vec::new();
        write_reversed_quality_bytes(&mut output, &[0, 30, 40])
            .expect("write_reversed_quality_bytes should succeed");
        // Reversed: [40, 30, 0] -> [73, 63, 33]
        assert_eq!(output, vec![73, 63, 33]);
    }

    #[test]
    fn test_quality_encoding_edge_cases() {
        let mut output = Vec::new();
        // Test max valid quality (93)
        write_quality_bytes(&mut output, &[93]).expect("write_quality_bytes should succeed");
        assert_eq!(output, vec![126]); // '~'

        output.clear();
        // Test overflow clamping (94+ should clamp to 126)
        write_quality_bytes(&mut output, &[94, 100, 255])
            .expect("write_quality_bytes should succeed");
        assert_eq!(output, vec![126, 126, 126]);
    }

    // ── UMI-in-header (mirrors `samtools fastq -U`) ─────────────────────────

    use fgumi_raw_bam::flags as fb;
    use rstest::rstest;
    use std::path::Path;

    #[rstest]
    #[case::default_list("RX,OX", vec![SamTag::RX, SamTag::OX])]
    #[case::single("RX", vec![SamTag::RX])]
    #[case::trims_whitespace(" RX , OX ", vec![SamTag::RX, SamTag::OX])]
    fn parse_umi_tags_accepts_valid(#[case] spec: &str, #[case] expected: Vec<SamTag>) {
        assert_eq!(parse_umi_tags(spec).expect("valid"), expected);
    }

    #[rstest]
    #[case::one_char("R")]
    #[case::empty("")]
    #[case::trailing_comma("RX,")]
    fn parse_umi_tags_rejects_malformed(#[case] spec: &str) {
        assert!(parse_umi_tags(spec).is_err(), "spec {spec:?} must be rejected");
    }

    #[rstest]
    // name_delim, umi_sep, stored UMI -> annotation
    #[case::single_umi(":", "+", b"AAAA".as_slice(), b":AAAA".as_slice())]
    #[case::duplex_hyphen_to_plus(":", "+", b"AAAA-GGGG", b":AAAA+GGGG")]
    #[case::multi_segment(":", "+", b"A-B-C", b":A+B+C")]
    // Underscore name separator (umi-tools style) with an empty UMI separator
    // (platforms that want the duplex UMIs concatenated with no delimiter).
    #[case::custom_delims_empty_sep("_", "", b"AAAA-GGGG", b"_AAAAGGGG")]
    fn umi_annotation_writes_expected(
        #[case] name_delim: &str,
        #[case] umi_sep: &str,
        #[case] umi: &[u8],
        #[case] expected: &[u8],
    ) {
        let uh = UmiHeader::from_args("RX", name_delim, umi_sep).expect("valid config");
        let mut out = Vec::new();
        uh.write_annotation(&mut out, umi).expect("write");
        assert_eq!(out, expected);
    }

    #[test]
    fn umi_header_rejects_bad_tag() {
        assert!(UmiHeader::from_args("nope", ":", "+").is_err());
    }

    // ── Paired split output (mirrors `samtools fastq -1/-2/-0`) ─────────────

    #[rstest]
    #[case::read1(fb::PAIRED | fb::FIRST_SEGMENT, Segment::Read1)]
    #[case::read2(fb::PAIRED | fb::LAST_SEGMENT, Segment::Read2)]
    // Both segment bits set, or neither, is "other" (single-end / ambiguous).
    #[case::both_bits(fb::PAIRED | fb::FIRST_SEGMENT | fb::LAST_SEGMENT, Segment::Other)]
    #[case::no_bits(0, Segment::Other)]
    fn classify_segment_by_flags(#[case] flags: u16, #[case] expected: Segment) {
        assert_eq!(classify_segment(flags), expected);
    }

    #[rstest]
    #[case::fq_gz("out_R1.fq.gz", true)]
    #[case::gz("out.gz", true)]
    #[case::bgz("out.bgz", true)]
    #[case::plain_fq("out_R1.fq", false)]
    #[case::fastq("out.fastq", false)]
    fn gz_extension_is_detected(#[case] path: &str, #[case] expected: bool) {
        assert_eq!(path_is_gzip(Path::new(path)), expected);
    }
}
