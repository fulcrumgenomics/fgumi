//! Convert BAM to FASTQ format.
//!
//! This tool reads a BAM file and outputs interleaved FASTQ to stdout for piping to aligners.
//! Input should be queryname-sorted or template-coordinate sorted.

use crate::commands::common::parse_bool;
use crate::logging::OperationTimer;
use crate::validation::validate_file_exists;
use anyhow::Result;
use clap::Parser;
use fgumi_bam_io::create_raw_bam_reader;
use fgumi_raw_bam::{
    RawRecord, extract_sequence_into, quality_scores_slice, read_name as raw_read_name,
};
use log::{info, warn};
use std::io::{BufWriter, Write, stdout};
use std::path::PathBuf;

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

    /// Number of threads for BAM decompression.
    #[arg(short = '@', short_alias = 't', long = "threads", default_value = "1")]
    pub threads: usize,

    /// BWA -K parameter value (bases per batch). Sizes output buffer to match bwa's
    /// batch size for optimal pipe throughput. Default matches common bwa mem usage.
    #[arg(short = 'K', long = "bwa-chunk-size", default_value = "150000000")]
    pub bwa_chunk_size: u64,
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
    fn run_with_writer<W: Write>(&self, mut writer: W) -> Result<()> {
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
                &mut writer,
                &record,
                flags,
                self.no_suffix,
                &mut seq_buf,
                &mut qual_buf,
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

impl Command for Fastq {
    fn execute(&self, _command_line: &str) -> Result<()> {
        validate_file_exists(&self.input, "Input BAM")?;

        // Refuse to clobber the input BAM. `File::create(path)` truncates
        // before the input is opened in `run_with_writer`, so an `--output`
        // pointing at the same file silently destroys the input data.
        // Lexical equality catches the obvious case; canonicalising both
        // sides also catches `./in.bam` vs `in.bam`, symlinks, and absolute
        // vs relative paths pointing at the same file.
        if let Some(output) = &self.output {
            let same_path = output == &self.input
                || (output.exists()
                    && std::fs::canonicalize(output)? == std::fs::canonicalize(&self.input)?);
            if same_path {
                anyhow::bail!(
                    "--output {} must differ from --input {} (would truncate the input BAM)",
                    output.display(),
                    self.input.display()
                );
            }
        }

        // Use 64MB buffer for efficient pipe throughput.
        const BUF_CAPACITY: usize = 64 * 1024 * 1024;
        match &self.output {
            Some(path) => {
                let file = std::fs::File::create(path)?;
                self.run_with_writer(BufWriter::with_capacity(BUF_CAPACITY, file))
            }
            None => self.run_with_writer(BufWriter::with_capacity(BUF_CAPACITY, stdout().lock())),
        }
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

/// Write a single FASTQ record to the writer.
#[inline]
fn write_fastq_record<W: Write>(
    writer: &mut W,
    record: &RawRecord,
    flags: u16,
    no_suffix: bool,
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

    // Copy quality bytes and transform to Phred+33 ASCII (absent quality → default)
    encode_quality_into(quality_scores_slice(record), qual_buf);

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
    use rstest::rstest;

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
}
