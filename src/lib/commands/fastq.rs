//! Convert BAM to FASTQ format.
//!
//! This tool reads a BAM file and outputs interleaved FASTQ to stdout for piping to aligners.
//! Input should be queryname-sorted or template-coordinate sorted.

use crate::logging::OperationTimer;
use crate::sam::SamTag;
use crate::validation::validate_input_exists;
use anyhow::Result;
use clap::Parser;
use fgumi_bam_io::{create_raw_bam_reader, is_stdin_path, is_stdout_path};
use fgumi_raw_bam::{
    RawRecord, aux_data_slice, extract_sequence_into, find_string_tag, quality_scores_slice,
    read_name as raw_read_name,
};
use log::{info, warn};
use std::io::{self, BufWriter, Write, stdout};
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
    // Values are shared from `fgumi_dna::COMPLEMENT` (IUPAC-aware, case-preserving),
    // but this table keeps the FASTQ-validity policy: any byte outside the valid
    // IUPAC alphabet folds to 'N', since an emitted FASTQ base must be a valid
    // nucleotide code (unlike the consensus/tag paths, which pass unknowns through).
    let mut table = [b'N'; 256];
    const VALID: &[u8] = b"ACGTURYSWKMBVDHNacgturyswkmbvdhn";
    let mut i = 0;
    while i < VALID.len() {
        let base = VALID[i];
        table[base as usize] = fgumi_dna::COMPLEMENT[base as usize];
        i += 1;
    }
    table
};

/// Convert BAM to FASTQ format.
#[derive(Debug, Parser)]
#[command(
    name = "fastq",
    about = "\x1b[38;5;72m[ALIGNMENT]\x1b[0m      \x1b[36mConvert BAM to FASTQ format\x1b[0m",
    long_about = r#"
Convert a BAM file to interleaved FASTQ format.

Reads BAM records and outputs FASTQ to stdout for piping to aligners.
Input should be queryname-sorted or template-coordinate sorted.

EXAMPLES:

  # Pipe to bwa mem for alignment
  fgumi fastq -i unmapped.bam | bwa mem -t 16 -p -K 150000000 -Y ref.fa -

  # With multi-threaded BAM decompression
  fgumi fastq -i unmapped.bam -@ 4 | bwa mem -t 16 -p ref.fa -

  # Exclude secondary and supplementary alignments (default)
  fgumi fastq -i aligned.bam -F 0x900 | bwa mem ...

NOTES:

  Read-name suffixes (/1, /2): appended from the FLAG unless --no-read-suffix,
  matching `samtools fastq -N`. A QNAME that already carries a mate suffix
  (a re-imported BAM; the SAM spec forbids it) is not stripped, so the output
  would double it (`name/1/1`) -- same as samtools -N. Pass --no-read-suffix
  to leave names untouched.

  Missing base qualities: a read with no stored quality (all 0xFF, the SAM
  no-quality sentinel) is emitted with a fixed Q33 ('B') per base, matching
  `samtools fastq` -- not Q93, which would falsely claim near-perfect quality.

  Interleaved pairing: output is interleaved for `bwa mem -p`. If a template
  contributes a lone mate (its pair is absent or dropped by --exclude-flags /
  --require-flags), the R1/R2 stream desyncs; fgumi warns when the paired R1/R2
  counts differ.
"#
)]
pub struct Fastq {
    /// Input BAM file.
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Output FASTQ file, or `-` / `/dev/stdout` for stdout. If omitted, the
    /// FASTQ stream is written to stdout (the default, intended for piping
    /// straight to an aligner). A `.gz`/`.bgz` path is BGZF-compressed; stdout
    /// is always plain text.
    #[arg(short = 'o', long = "output")]
    pub output: Option<PathBuf>,

    /// Don't append /1 and /2 to read names.
    #[arg(short = 'n', long = "no-read-suffix", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub no_suffix: bool,

    /// Exclude reads with any of these flags present [0x900 = secondary|supplementary].
    #[arg(short = 'F', long = "exclude-flags", default_value_t = 0x900, value_parser = parse_flags)]
    pub exclude_flags: u16,

    /// Only include reads with all of these flags present.
    #[arg(short = 'f', long = "require-flags", default_value_t = 0, value_parser = parse_flags)]
    pub require_flags: u16,

    /// Number of threads for BAM decompression.
    #[arg(short = '@', short_alias = 't', long = "threads", default_value = "1")]
    pub threads: usize,

    /// BWA -K parameter value (bases per batch). Sizes output buffer to match bwa's
    /// batch size for optimal pipe throughput. Default matches common bwa mem usage.
    #[arg(short = 'K', long = "bwa-chunk-size", default_value = "150000000")]
    pub bwa_chunk_size: u64,

    /// Append the record's UMI to the read name, before any /1 or /2 suffix.
    ///
    /// With the default delimiters this matches `samtools fastq -U`
    /// (`readname:AAAA+CCCC`), the layout DRAGEN expects.
    #[arg(short = 'a', short_alias = 'U', long = "annotate-read-names", value_name = "true|false", default_value = "false", num_args = 0..=1, default_missing_value = "true", action = clap::ArgAction::Set, value_parser = clap::builder::BoolishValueParser::new(), hide_possible_values = true)]
    pub annotate_read_names: bool,

    /// Tags to read the UMI from, in priority order; the first present wins.
    #[arg(long = "umi-tag", default_value = "RX,OX", value_delimiter = ',')]
    pub umi_tag: Vec<String>,

    /// Delimiter between the read name and the UMI.
    #[arg(long = "umi-name-delim", default_value = ":")]
    pub umi_name_delim: String,

    /// Separator between the two halves of a duplex UMI in the read name.
    ///
    /// fgumi stores duplex UMIs as `AAAA-CCCC`; `samtools fastq -U` and DRAGEN
    /// expect `AAAA+CCCC`, so the stored `-` is rewritten to this value.
    #[arg(long = "umi-sep", default_value = "+")]
    pub umi_sep: String,
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
    /// Run the BAM-to-FASTQ conversion loop against the given writer.
    ///
    /// Extracted so `execute` can dispatch between stdout (the default piping
    /// path) and a file (`--output`) without duplicating the conversion loop.
    fn run_with_writer<W: Write>(&self, writer: &mut W) -> Result<()> {
        use fgumi_raw_bam::flags as raw_flag_bits;
        let timer = OperationTimer::new("Converting BAM to FASTQ");

        info!("Input: {}", self.input.display());
        info!("Threads: {}", self.threads);
        info!("Exclude flags: 0x{:X}", self.exclude_flags);
        info!("Require flags: 0x{:X}", self.require_flags);
        info!("Read name suffix: {}", if self.no_suffix { "disabled" } else { "enabled" });
        info!("BWA chunk size: {} bases", self.bwa_chunk_size);

        let (mut reader, _header) = create_raw_bam_reader(&self.input, self.threads)?;

        let mut total_records: u64 = 0;
        let mut written_records: u64 = 0;
        // Track paired R1/R2 balance to detect singletons that would desync the
        // interleaved stream under `bwa mem -p` (FASTQ3-02).
        let mut paired_r1_written: u64 = 0;
        let mut paired_r2_written: u64 = 0;

        // Batch tracking (mirrors bwa's logic)
        let mut bases_this_batch: u64 = 0;
        let mut records_this_batch: usize = 0;

        // Reusable buffers to avoid per-record allocations
        let mut buffers = FastqRecordBuffers::with_capacity(512);
        let mut record = RawRecord::new();

        // Resolve the UMI annotation config once, outside the per-record loop.
        let umi_annotation = if self.annotate_read_names {
            Some(UmiNameAnnotation::new(&self.umi_tag, &self.umi_name_delim, &self.umi_sep)?)
        } else {
            None
        };

        loop {
            let n = reader.read_record(&mut record)?;
            if n == 0 {
                break; // EOF
            }
            total_records += 1;

            let flags = record.flags();

            // Filter by flags
            if (flags & self.exclude_flags) != 0 {
                continue;
            }
            if (flags & self.require_flags) != self.require_flags {
                continue;
            }

            // Get sequence length for batch tracking
            let seq_len = record.l_seq() as usize;

            // Write FASTQ record
            write_fastq_record(
                &mut *writer,
                &record,
                flags,
                self.no_suffix,
                &mut buffers,
                umi_annotation.as_ref(),
            )?;
            written_records += 1;
            if (flags & raw_flag_bits::PAIRED) != 0 {
                if (flags & raw_flag_bits::FIRST_SEGMENT) != 0 {
                    paired_r1_written += 1;
                }
                if (flags & raw_flag_bits::LAST_SEGMENT) != 0 {
                    paired_r2_written += 1;
                }
            }

            // Track batch progress
            bases_this_batch += seq_len as u64;
            records_this_batch += 1;

            // Flush at batch boundary (like bwa: bases >= K AND record count is even)
            if bases_this_batch >= self.bwa_chunk_size && records_this_batch.is_multiple_of(2) {
                writer.flush()?;
                bases_this_batch = 0;
                records_this_batch = 0;
            }
        }

        // Final flush
        writer.flush()?;

        // FASTQ3-02: this is interleaved output. If the paired R1/R2 counts differ,
        // at least one template contributed a lone mate (its pair was absent or
        // dropped by `--exclude-flags`/`--require-flags`), which desyncs R1/R2 pairing
        // for every subsequent read under `bwa mem -p`. Warn rather than fail —
        // single-end input and deliberate mate filtering are legitimate uses.
        //
        // NOTE: this is a best-effort NET-count heuristic, not per-QNAME pairing. Two
        // singletons that cancel — one template drops its R2, another drops its R1 —
        // leave R1==R2 and are NOT detected here (a robust check would need per-read-name
        // pairing, impractical for a streaming converter). So the absence of a warning
        // does not guarantee a perfectly-synced interleaved stream.
        if paired_r1_written != paired_r2_written {
            // Net imbalance, a lower bound on the true singleton count: canceling
            // singletons (one template drops its R2, another its R1) leave R1==R2
            // and are not counted here (see NOTE above).
            let net_imbalance = paired_r1_written.abs_diff(paired_r2_written);
            warn!(
                "Interleaved FASTQ paired-read counts differ (R1={paired_r1_written}, \
                 R2={paired_r2_written}): net R1/R2 imbalance of {net_imbalance} (at least \
                 {net_imbalance} singleton mate(s)) will desync R1/R2 pairing under \
                 `bwa mem -p`. Ensure the input is queryname/template-coordinate sorted with \
                 both mates present, or reconsider `--exclude-flags` (0x{:X}) / `--require-flags` \
                 (0x{:X}).",
                self.exclude_flags, self.require_flags
            );
        }

        info!("Read {total_records} records, wrote {written_records} FASTQ records");
        timer.log_completion(written_records);
        Ok(())
    }
}

/// Refuse an `--output` that would clobber the `--input` BAM.
///
/// `File::create(path)` truncates before the input is opened in
/// `run_with_writer`, so an `--output` naming the same file silently destroys
/// the input data. Lexical equality catches the obvious case; canonicalising
/// both sides also catches `./in.bam` vs `in.bam`, symlinks, and absolute vs
/// relative paths pointing at the same file.
///
/// # Errors
///
/// Returns an error if `output` and `input` name the same file, or if either
/// path exists but cannot be canonicalised.
fn reject_output_clobbering_input(input: &Path, output: Option<&PathBuf>) -> Result<()> {
    let Some(output) = output else { return Ok(()) };

    // Stdin has no filesystem entity to canonicalise -- `canonicalize("-")`
    // fails with `NotFound` -- and cannot name the output file, so there is
    // nothing to clobber. Without this the guard turns `--input -` into a
    // confusing IO error whenever `--output` happens to already exist.
    if is_stdin_path(input) {
        return Ok(());
    }

    // Likewise stdout: `-` names a stream, not a file, so it cannot clobber the
    // input — and a stray `-` file left in the working directory must not make
    // this guard canonicalise something that is not the destination.
    if is_stdout_path(output) {
        return Ok(());
    }

    let same_path = output == input
        || (output.exists() && std::fs::canonicalize(output)? == std::fs::canonicalize(input)?);
    if same_path {
        anyhow::bail!(
            "--output {} must differ from --input {} (would truncate the input BAM)",
            output.display(),
            input.display()
        );
    }
    Ok(())
}

impl Command for Fastq {
    fn execute(&self, _command_line: &str) -> Result<()> {
        validate_input_exists(&self.input, "Input BAM")?;
        reject_output_clobbering_input(&self.input, self.output.as_ref())?;

        match &self.output {
            // `-` and `/dev/stdout` name the stream that omitting `--output`
            // already writes. Spelling it out is what makes this command follow
            // the same `-` convention as the rest of the CLI; before, `-o -`
            // created a regular file called `-`.
            Some(path) if is_stdout_path(path) => self.write_to_stdout(),
            // A `.gz`/`.bgz` output path must actually be compressed. Writing plain
            // text under a `.gz` name produces a file every downstream tool rejects.
            // BGZF is gzip-compatible, so `zcat`/`gzip -d` read it, and it stays
            // block-indexable.
            Some(path) if is_gzip_output_path(path) => {
                let file = std::fs::File::create(path)?;
                let mut writer = BgzfFastqWriter::new(file, BGZF_OUTPUT_COMPRESSION_LEVEL);
                self.run_with_writer(&mut writer)?;
                // Append the BGZF EOF marker so the stream is not seen as truncated.
                writer.finish()?;
                Ok(())
            }
            Some(path) => {
                let file = std::fs::File::create(path)?;
                let mut writer = BufWriter::with_capacity(FASTQ_BUF_CAPACITY, file);
                self.run_with_writer(&mut writer)?;
                // `BufWriter`'s `Drop` flushes but discards the error, so a failed
                // final write — a full disk here — would leave a truncated FASTQ
                // behind an exit code of zero. Flush explicitly to report it, as
                // the `.gz` arm's `finish` already does.
                writer.flush()?;
                Ok(())
            }
            None => self.write_to_stdout(),
        }
    }
}

/// Output buffer size — 64MB, for efficient pipe throughput.
const FASTQ_BUF_CAPACITY: usize = 64 * 1024 * 1024;

impl Fastq {
    /// Writes uncompressed FASTQ to stdout.
    ///
    /// Uncompressed regardless of how stdout was asked for: this stream is meant
    /// to be piped into an aligner, which wants plain FASTQ.
    ///
    /// # Errors
    ///
    /// Returns an error if a record cannot be written, or if the final flush
    /// fails — `EPIPE` when the aligner on the other end exits early. Dropping
    /// the [`BufWriter`] would swallow that flush error and report success on a
    /// truncated stream.
    fn write_to_stdout(&self) -> Result<()> {
        let mut writer = BufWriter::with_capacity(FASTQ_BUF_CAPACITY, stdout().lock());
        self.run_with_writer(&mut writer)?;
        writer.flush()?;
        Ok(())
    }
}

/// ASCII quality emitted for a base whose quality is entirely absent.
///
/// Per the SAM spec, a record with no quality stores `0xFF` for every base.
/// fgumi previously mapped `0xFF` through [`QUAL_TO_ASCII`] to `~` (Q93), which
/// falsely claims near-perfect quality for a read that has none. Instead emit
/// Q33 (`B`) for an entirely-absent quality string, matching `samtools fastq`
/// (which `fgumi fastq` is modeled on).
const MISSING_QUALITY_ASCII: u8 = b'B';

/// Encode a record's raw quality bytes as Phred+33 ASCII into `out`.
///
/// When the quality is entirely absent (all `0xFF`), fills with
/// [`MISSING_QUALITY_ASCII`] rather than the misleading Q93 that
/// `QUAL_TO_ASCII[0xFF]` would produce. A partially-`0xFF` string is left to the
/// per-byte mapping (a genuinely malformed record, not the "no quality" sentinel).
fn encode_quality_into(quals: &[u8], out: &mut Vec<u8>) {
    out.clear();
    if !quals.is_empty() && quals.iter().all(|&q| q == 0xFF) {
        out.resize(quals.len(), MISSING_QUALITY_ASCII);
    } else {
        out.extend(quals.iter().map(|&s| QUAL_TO_ASCII[s as usize]));
    }
}

/// BGZF compression level used for `.gz`/`.bgz` FASTQ output.
///
/// 6 is the zlib/bgzip default: a middle trade-off between output size and CPU,
/// appropriate for an intermediate file handed straight to an aligner.
const BGZF_OUTPUT_COMPRESSION_LEVEL: u32 = 6;

/// Returns `true` if `path` names a gzip-family output that must be compressed.
fn is_gzip_output_path(path: &std::path::Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("gz") || e.eq_ignore_ascii_case("bgz"))
}

/// `io::Write` adapter that block-compresses FASTQ text as BGZF.
///
/// BGZF is gzip-compatible, so any consumer that reads `.gz` reads this, while the
/// output stays block-indexable. The BGZF EOF marker is appended on drop-free
/// finalization (`finish`), which `Write::flush` alone must not do — flush can be
/// called mid-stream.
struct BgzfFastqWriter<W: Write> {
    inner: W,
    compressor: fgumi_bgzf::InlineBgzfCompressor,
}

impl<W: Write> BgzfFastqWriter<W> {
    fn new(inner: W, compression_level: u32) -> Self {
        Self { inner, compressor: fgumi_bgzf::InlineBgzfCompressor::new(compression_level) }
    }

    /// Flushes remaining buffered bytes and appends the BGZF EOF marker so the
    /// stream is complete and tools do not report a truncated file.
    fn finish(&mut self) -> io::Result<()> {
        self.compressor.flush()?;
        self.compressor.write_blocks_to(&mut self.inner)?;
        self.inner.write_all(&fgumi_bgzf::BGZF_EOF)?;
        self.inner.flush()
    }
}

impl<W: Write> Write for BgzfFastqWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.compressor.write_all(buf)?;
        // Drain any blocks the compressor completed so memory stays bounded.
        self.compressor.write_blocks_to(&mut self.inner)?;
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.compressor.write_blocks_to(&mut self.inner)?;
        self.inner.flush()
    }
}

/// How to append a record's UMI to its read name.
///
/// Built once per run from the CLI options so the per-record path does no
/// string parsing or allocation beyond the name buffer itself.
#[derive(Debug, Clone)]
pub struct UmiNameAnnotation {
    /// Tags to consult, in priority order; the first present wins.
    tags: Vec<SamTag>,
    /// Delimiter between the read name and the UMI.
    name_delim: Vec<u8>,
    /// Separator written between duplex UMI halves (replaces the stored `-`).
    umi_sep: Vec<u8>,
}

impl UmiNameAnnotation {
    /// Builds the annotation config from raw CLI strings.
    ///
    /// # Errors
    ///
    /// Returns an error if a tag name is not exactly two characters.
    pub fn new(tags: &[String], name_delim: &str, umi_sep: &str) -> Result<Self> {
        let tags = tags
            .iter()
            .map(|t| {
                let bytes = t.as_bytes();
                if bytes.len() != 2 || !SamTag::is_valid_tag_bytes(bytes[0], bytes[1]) {
                    anyhow::bail!(
                        "--umi-tag values must be a valid two-character SAM tag, got '{t}'"
                    );
                }
                Ok(SamTag::new(bytes[0], bytes[1]))
            })
            .collect::<Result<Vec<_>>>()?;
        if tags.is_empty() {
            anyhow::bail!("--umi-tag requires at least one tag");
        }
        Ok(Self {
            tags,
            name_delim: name_delim.as_bytes().to_vec(),
            umi_sep: umi_sep.as_bytes().to_vec(),
        })
    }

    /// Appends `<delim><umi>` to `out` if the record carries one of the configured
    /// tags, rewriting the stored duplex `-` separator to the configured one.
    ///
    /// Records without any of the tags are left unannotated rather than failing —
    /// a BAM can legitimately mix UMI-bearing and UMI-free reads.
    fn append_to(&self, record: &RawRecord, out: &mut Vec<u8>) {
        let aux = aux_data_slice(record);
        let Some(umi) = self.tags.iter().find_map(|&tag| find_string_tag(aux, tag)) else {
            return;
        };
        out.extend_from_slice(&self.name_delim);
        for &byte in umi {
            if byte == b'-' {
                out.extend_from_slice(&self.umi_sep);
            } else {
                out.push(byte);
            }
        }
    }
}

/// Reusable per-record scratch buffers, kept across the conversion loop so the
/// hot path does not allocate.
#[derive(Debug, Default)]
pub(crate) struct FastqRecordBuffers {
    /// Decoded sequence bases.
    seq: Vec<u8>,
    /// Phred+33 encoded qualities.
    qual: Vec<u8>,
    /// Read name, used only when a UMI is appended.
    name: Vec<u8>,
}

impl FastqRecordBuffers {
    pub(crate) fn with_capacity(capacity: usize) -> Self {
        Self {
            seq: Vec::with_capacity(capacity),
            qual: Vec::with_capacity(capacity),
            name: Vec::with_capacity(256),
        }
    }
}

/// Write a single FASTQ record to the writer.
#[inline]
pub(crate) fn write_fastq_record<W: Write>(
    writer: &mut W,
    record: &RawRecord,
    flags: u16,
    no_suffix: bool,
    buffers: &mut FastqRecordBuffers,
    umi_annotation: Option<&UmiNameAnnotation>,
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

    // Decode sequence from 4-bit BAM encoding to ASCII bases (reuses the scratch buffer).
    extract_sequence_into(record, &mut buffers.seq);

    // Copy quality bytes and transform to Phred+33 ASCII (absent quality → default)
    encode_quality_into(quality_scores_slice(record), &mut buffers.qual);

    if (flags & flag_bits::REVERSE) != 0 {
        // Reverse complement sequence in place using lookup table
        buffers.seq.reverse();
        for base in buffers.seq.iter_mut() {
            *base = COMPLEMENT[*base as usize];
        }
        // Reverse quality in place
        buffers.qual.reverse();
    }

    // Write all parts. The UMI goes between the name and the /1 /2 suffix, matching
    // `samtools fastq -U`.
    writer.write_all(b"@")?;
    if let Some(annotation) = umi_annotation {
        buffers.name.clear();
        buffers.name.extend_from_slice(name);
        annotation.append_to(record, &mut buffers.name);
        writer.write_all(&buffers.name)?;
    } else {
        writer.write_all(name)?;
    }
    writer.write_all(suffix)?;
    writer.write_all(b"\n")?;
    writer.write_all(&buffers.seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(&buffers.qual)?;
    writer.write_all(b"\n")?;

    Ok(())
}

// ==== ported from feat-runall for the chain builder (R2) ====
/// Which FASTQ stream a record belongs to, based on its segment flags.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Segment {
    /// First segment of a pair (`FIRST_SEGMENT` set, `LAST_SEGMENT` clear) → R1.
    Read1,
    /// Last segment of a pair (`LAST_SEGMENT` set, `FIRST_SEGMENT` clear) → R2.
    Read2,
    /// Neither or both segment bits set (single-end or ambiguous) → "other".
    Other,
}

/// Classify a record into R1 / R2 / other from its flags, matching how
/// `samtools fastq` routes reads to `-1` / `-2` / `-0`.
pub(crate) fn classify_segment(flags: u16) -> Segment {
    use fgumi_raw_bam::flags as flag_bits;
    let is_first = (flags & flag_bits::FIRST_SEGMENT) != 0;
    let is_last = (flags & flag_bits::LAST_SEGMENT) != 0;
    match (is_first, is_last) {
        (true, false) => Segment::Read1,
        (false, true) => Segment::Read2,
        _ => Segment::Other,
    }
}

/// Returns `true` if `path` should be written compressed (ends in
/// `.gz`/`.bgz`/`.bgzf`).
///
/// KNOWN DIVERGENCE — reconcile in the fastq WIRING PR: this is case-SENSITIVE
/// and also matches `.bgzf`, whereas the live standalone equivalent
/// [`is_gzip_output_path`] is case-INSENSITIVE and matches only `gz`/`bgz`. The
/// two decide the same thing ("is this output gzip-compressed?") and should be a
/// single shared function once the fastq stage is wired onto the chain; they
/// cannot diverge in observable behavior today because this path is dormant.
pub(crate) fn path_is_gzip(path: &Path) -> bool {
    matches!(path.extension().and_then(|e| e.to_str()), Some("gz" | "bgz" | "bgzf"))
}

// ==== ported from feat-runall for the chain builder (R2) ====
/// Per-stage options for [`crate::pipeline::chains::Stage::Fastq`] — the
/// BAM→FASTQ encode. Carries the flag filters, read-name suffix behavior, and
/// the optional UMI-in-read-name config the encode step needs. The output
/// destination(s) and compression are carried by the chain's `SinkSpec`, not
/// here.
#[derive(Clone)]
pub struct FastqOptions {
    /// Exclude reads with any of these flag bits set.
    pub exclude_flags: u16,
    /// Only include reads with all of these flag bits set.
    pub require_flags: u16,
    /// When `true`, omit the `/1` `/2` read-name suffix.
    pub no_suffix: bool,
    /// When `Some`, append the UMI (from the configured tag) to the read name.
    pub umi_header: Option<UmiNameAnnotation>,
}

// ==== impls ported from feat-runall for the chain builder (R2) ====
impl FastqOptions {
    /// Returns `true` if a record with `flags` passes the include/exclude filters.
    #[must_use]
    pub(crate) fn passes_filters(&self, flags: u16) -> bool {
        (flags & self.exclude_flags) == 0 && (flags & self.require_flags) == self.require_flags
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    /// Build a single-end record carrying an optional UMI tag.
    fn record_with_umi(name: &[u8], umi: Option<(&str, &[u8])>) -> RawRecord {
        let mut b = fgumi_raw_bam::SamBuilder::new();
        b.read_name(name);
        if let Some((tag, value)) = umi {
            let t = tag.as_bytes();
            b.add_string_tag(SamTag::new(t[0], t[1]), value);
        }
        b.build()
    }

    #[rstest]
    #[case::gz_lower("out.fq.gz", true)]
    #[case::gz_upper("OUT.FQ.GZ", true)]
    #[case::bgz("out.fq.bgz", true)]
    #[case::plain_fq("out.fq", false)]
    #[case::plain_fastq("out.fastq", false)]
    #[case::no_extension("out", false)]
    #[case::gz_in_stem("out.gz.fq", false)]
    fn test_is_gzip_output_path(#[case] path: &str, #[case] expected: bool) {
        assert_eq!(is_gzip_output_path(std::path::Path::new(path)), expected);
    }

    #[rstest]
    // Default delimiters reproduce `samtools fastq -U`: name:UMI, duplex `-` -> `+`.
    #[case::simplex_umi(Some(("RX", &b"ACGT"[..])), ":", "+", "read1:ACGT")]
    #[case::duplex_umi(Some(("RX", &b"ACGT-TTTT"[..])), ":", "+", "read1:ACGT+TTTT")]
    #[case::custom_delims(Some(("RX", &b"ACGT-TTTT"[..])), "_", "|", "read1_ACGT|TTTT")]
    // A record without any of the configured tags is left unannotated rather than failing.
    #[case::no_umi_tag(None, ":", "+", "read1")]
    // OX is consulted when RX is absent (tag priority order).
    #[case::fallback_tag(Some(("OX", &b"GGGG"[..])), ":", "+", "read1:GGGG")]
    fn test_umi_name_annotation(
        #[case] umi: Option<(&str, &[u8])>,
        #[case] name_delim: &str,
        #[case] umi_sep: &str,
        #[case] expected: &str,
    ) {
        let record = record_with_umi(b"read1", umi);
        let annotation =
            UmiNameAnnotation::new(&["RX".to_string(), "OX".to_string()], name_delim, umi_sep)
                .expect("valid annotation config");

        let mut out = b"read1".to_vec();
        annotation.append_to(&record, &mut out);
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[rstest]
    #[case::too_short("R")]
    #[case::too_long("RXX")]
    #[case::invalid_first_char("1X")]
    fn test_umi_name_annotation_rejects_bad_tag(#[case] tag: &str) {
        assert!(UmiNameAnnotation::new(&[tag.to_string()], ":", "+").is_err());
    }

    #[test]
    fn test_umi_name_annotation_rejects_empty_tag_list() {
        assert!(UmiNameAnnotation::new(&[], ":", "+").is_err());
    }

    /// The BGZF writer must emit a real, decompressible BGZF stream terminated by
    /// the EOF marker — a plain-text file under a `.gz` name is exactly the bug
    /// this path exists to prevent.
    #[test]
    fn test_bgzf_fastq_writer_round_trips() {
        let payload = b"@read1\nACGT\n+\nIIII\n@read2\nTTTT\n+\nJJJJ\n";
        let mut out: Vec<u8> = Vec::new();
        {
            let mut writer = BgzfFastqWriter::new(&mut out, BGZF_OUTPUT_COMPRESSION_LEVEL);
            writer.write_all(payload).expect("write");
            writer.finish().expect("finish");
        }

        // Real BGZF/gzip magic, not plain text.
        assert_eq!(&out[..2], &[0x1f, 0x8b], "output must carry gzip magic");
        assert_ne!(&out[..1], b"@", "output must not be plain FASTQ text");

        // Terminated by the BGZF EOF marker so readers do not see a truncated file.
        assert!(out.ends_with(&fgumi_bgzf::BGZF_EOF), "BGZF stream must end with the EOF marker");

        // And it decompresses back to exactly what went in.
        let mut cursor = std::io::Cursor::new(&out);
        let blocks = fgumi_bgzf::read_raw_blocks(&mut cursor, 64).expect("read blocks");
        let mut decompressor = libdeflater::Decompressor::new();
        let mut decoded = Vec::new();
        for block in &blocks {
            decoded.extend_from_slice(
                &fgumi_bgzf::decompress_block(block, &mut decompressor).expect("decompress block"),
            );
        }
        assert_eq!(decoded, payload);
    }

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

        // Lowercase (case preserved)
        assert_eq!(COMPLEMENT[b'a' as usize], b't');
        assert_eq!(COMPLEMENT[b't' as usize], b'a');
        assert_eq!(COMPLEMENT[b'c' as usize], b'g');
        assert_eq!(COMPLEMENT[b'g' as usize], b'c');
        assert_eq!(COMPLEMENT[b'n' as usize], b'n'); // case preserved (shared table)

        // Unknown bases map to N
        assert_eq!(COMPLEMENT[b'X' as usize], b'N');
        assert_eq!(COMPLEMENT[0], b'N');
    }

    #[test]
    fn test_complement_lookup_table_iupac_and_case() {
        // IUPAC-aware (values shared from fgumi_dna::COMPLEMENT), case-preserving,
        // but keeps the FASTQ-validity policy of folding invalid bytes to N.
        assert_eq!(COMPLEMENT[b'R' as usize], b'Y');
        assert_eq!(COMPLEMENT[b'Y' as usize], b'R');
        assert_eq!(COMPLEMENT[b'K' as usize], b'M');
        assert_eq!(COMPLEMENT[b'S' as usize], b'S');
        assert_eq!(COMPLEMENT[b'a' as usize], b't'); // case preserved
        assert_eq!(COMPLEMENT[b'r' as usize], b'y');
        assert_eq!(COMPLEMENT[b'.' as usize], b'N'); // invalid -> N (FASTQ validity)
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

    /// FASTQ3-03: a record whose quality is entirely absent (all `0xFF`, the SAM
    /// no-quality sentinel) must not be emitted as `~` (Q93, near-perfect); it
    /// gets the `MISSING_QUALITY_ASCII` default (`B`, Q33, matching samtools). A
    /// present quality maps per-byte; a partially-`0xFF` (malformed) string is NOT
    /// the sentinel and is left to the map; empty stays empty.
    #[rstest]
    #[case::absent_all_0xff(&[0xFF, 0xFF, 0xFF, 0xFF], b"BBBB")]
    #[case::present_phred33(&[0, 30, 40, 93], &[33, 63, 73, 126])]
    #[case::partial_0xff_is_not_sentinel(&[30, 0xFF], &[63, 126])]
    #[case::empty(&[], b"")]
    fn encode_quality_maps_present_and_defaults_absent(
        #[case] quals: &[u8],
        #[case] expected: &[u8],
    ) {
        let mut out = Vec::new();
        encode_quality_into(quals, &mut out);
        assert_eq!(out.as_slice(), expected);
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

    /// Resolve a case-table path spec against `dir`. `-` and `/dev/stdin` are
    /// returned verbatim; anything else names a file inside `dir`.
    fn resolve_path(dir: &std::path::Path, spec: &str) -> PathBuf {
        if spec == "-" || spec == "/dev/stdin" { PathBuf::from(spec) } else { dir.join(spec) }
    }

    /// The clobber guard must reject an `--output` that names the `--input`
    /// file, and must exempt stdin: `canonicalize("-")` fails with `NotFound`,
    /// so canonicalising a stdin input would turn `--input -` into a confusing
    /// IO error the moment `--output` happens to already exist.
    #[rstest]
    #[case::stdin_dash_with_existing_output("-", "out.fq", true, true)]
    #[case::stdin_dev_stdin_with_existing_output("/dev/stdin", "out.fq", true, true)]
    #[case::stdin_dash_with_new_output("-", "out.fq", false, true)]
    #[case::distinct_paths_output_exists("in.bam", "out.fq", true, true)]
    #[case::distinct_paths_output_missing("in.bam", "out.fq", false, true)]
    #[case::output_is_input("in.bam", "in.bam", true, false)]
    fn test_reject_output_clobbering_input(
        #[case] input: &str,
        #[case] output: &str,
        #[case] output_exists: bool,
        #[case] accepted: bool,
    ) {
        let dir = tempfile::TempDir::new().expect("create temp dir");
        std::fs::write(dir.path().join("in.bam"), b"not a real bam").expect("write input");

        let input_path = resolve_path(dir.path(), input);
        let output_path = resolve_path(dir.path(), output);
        if output_exists && !output_path.exists() {
            std::fs::write(&output_path, b"prior run").expect("write output");
        }

        let result = reject_output_clobbering_input(&input_path, Some(&output_path));
        assert_eq!(
            result.is_ok(),
            accepted,
            "unexpected verdict for --input {input} --output {output}: {result:?}"
        );
        if !accepted {
            let err = result.expect_err("case expects rejection").to_string();
            assert!(err.contains("must differ"), "unexpected error message: {err}");
        }
    }

    /// No `--output` means FASTQ goes to stdout, so there is nothing to clobber.
    #[test]
    fn test_reject_output_clobbering_input_allows_no_output() {
        reject_output_clobbering_input(&PathBuf::from("in.bam"), None)
            .expect("stdout output has nothing to clobber");
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

    /// `classify_segment` routes by the FIRST/LAST segment bits exactly as
    /// `samtools fastq` does: first-only → R1, last-only → R2, neither/both → other.
    #[rstest]
    #[case::read1(fgumi_raw_bam::flags::FIRST_SEGMENT, Segment::Read1)]
    #[case::read2(fgumi_raw_bam::flags::LAST_SEGMENT, Segment::Read2)]
    #[case::both_bits(
        fgumi_raw_bam::flags::FIRST_SEGMENT | fgumi_raw_bam::flags::LAST_SEGMENT,
        Segment::Other
    )]
    #[case::neither(0, Segment::Other)]
    fn classify_segment_routes_by_first_last_bits(#[case] flags: u16, #[case] expected: Segment) {
        assert_eq!(classify_segment(flags), expected);
    }

    /// `passes_filters` requires ALL `require_flags` set and NONE of
    /// `exclude_flags` set (an inverted bit-mask here would pass silently).
    #[rstest]
    // exclude FIRST_SEGMENT: an R1 read is excluded, others pass.
    #[case::excluded(
        fgumi_raw_bam::flags::FIRST_SEGMENT,
        0,
        fgumi_raw_bam::flags::FIRST_SEGMENT,
        false
    )]
    #[case::not_excluded(fgumi_raw_bam::flags::FIRST_SEGMENT, 0, 0, true)]
    // require PAIRED: only paired reads pass.
    #[case::required_present(0, fgumi_raw_bam::flags::PAIRED, fgumi_raw_bam::flags::PAIRED, true)]
    #[case::required_absent(0, fgumi_raw_bam::flags::PAIRED, 0, false)]
    // both: must satisfy require AND avoid exclude.
    #[case::require_met_exclude_hit(
        fgumi_raw_bam::flags::FIRST_SEGMENT,
        fgumi_raw_bam::flags::PAIRED,
        fgumi_raw_bam::flags::PAIRED | fgumi_raw_bam::flags::FIRST_SEGMENT,
        false
    )]
    #[case::require_met_exclude_clear(
        fgumi_raw_bam::flags::FIRST_SEGMENT,
        fgumi_raw_bam::flags::PAIRED,
        fgumi_raw_bam::flags::PAIRED,
        true
    )]
    fn passes_filters_applies_exclude_and_require(
        #[case] exclude_flags: u16,
        #[case] require_flags: u16,
        #[case] flags: u16,
        #[case] expected: bool,
    ) {
        let opts =
            FastqOptions { exclude_flags, require_flags, no_suffix: false, umi_header: None };
        assert_eq!(opts.passes_filters(flags), expected);
    }

    /// `path_is_gzip` matches `gz`/`bgz`/`bgzf`, case-sensitively. See the fn's
    /// KNOWN DIVERGENCE note vs the live `is_gzip_output_path` (the uppercase
    /// case below pins that divergence, to be reconciled at fastq wiring).
    #[rstest]
    #[case::gz("out.fq.gz", true)]
    #[case::bgz("out.fq.bgz", true)]
    #[case::bgzf("out.fq.bgzf", true)]
    #[case::plain_fq("out.fq", false)]
    #[case::no_extension("out", false)]
    #[case::gz_in_stem("out.gz.fq", false)]
    #[case::uppercase_not_matched("OUT.FQ.GZ", false)]
    fn path_is_gzip_matches_gz_bgz_bgzf(#[case] path: &str, #[case] expected: bool) {
        assert_eq!(path_is_gzip(std::path::Path::new(path)), expected);
    }
}
