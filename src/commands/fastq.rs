//! Convert BAM to FASTQ format.
//!
//! This tool reads a BAM file and outputs interleaved FASTQ to stdout for piping to aligners.
//! Input should be queryname-sorted or template-coordinate sorted.

use anyhow::{Context, Result};
use clap::Parser;
use fgumi_lib::bam_io::create_bam_reader;
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::validation::validate_file_exists;
use log::info;
use noodles::bam;
use noodles::sam::alignment::record::Flags;
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
"#
)]
pub struct Fastq {
    /// Input BAM file.
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Don't append /1 and /2 to read names.
    #[arg(short = 'n', long = "no-read-suffix", default_value = "false")]
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

impl Command for Fastq {
    fn execute(&self, _command_line: &str) -> Result<()> {
        validate_file_exists(&self.input, "Input BAM")?;

        let timer = OperationTimer::new("Converting BAM to FASTQ");

        info!("Input: {}", self.input.display());
        info!("Threads: {}", self.threads);
        info!("Exclude flags: 0x{:X}", self.exclude_flags);
        info!("Require flags: 0x{:X}", self.require_flags);
        info!("Read name suffix: {}", if self.no_suffix { "disabled" } else { "enabled" });
        info!("BWA chunk size: {} bases", self.bwa_chunk_size);

        let (mut reader, _header) = create_bam_reader(&self.input, self.threads)?;

        // Use 64MB buffer for efficient pipe throughput
        let mut writer = BufWriter::with_capacity(64 * 1024 * 1024, stdout().lock());

        let mut total_records: u64 = 0;
        let mut written_records: u64 = 0;

        // Batch tracking (mirrors bwa's logic)
        let mut bases_this_batch: u64 = 0;
        let mut records_this_batch: usize = 0;

        // Reusable buffers to avoid per-record allocations
        let mut seq_buf: Vec<u8> = Vec::with_capacity(512);
        let mut qual_buf: Vec<u8> = Vec::with_capacity(512);

        for result in reader.records() {
            let record = result.context("Failed to read BAM record")?;
            total_records += 1;

            let flags = record.flags();

            // Filter by flags
            if (flags.bits() & self.exclude_flags) != 0 {
                continue;
            }
            if (flags.bits() & self.require_flags) != self.require_flags {
                continue;
            }

            // Get sequence length for batch tracking
            let seq_len = record.sequence().len();

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

        info!("Read {total_records} records, wrote {written_records} FASTQ records");
        timer.log_completion(written_records);
        Ok(())
    }
}

/// Write a single FASTQ record to the writer.
#[inline]
fn write_fastq_record<W: Write>(
    writer: &mut W,
    record: &bam::Record,
    flags: Flags,
    no_suffix: bool,
    seq_buf: &mut Vec<u8>,
    qual_buf: &mut Vec<u8>,
) -> Result<()> {
    // Get read name - use concrete method for direct access
    let name = record.name().context("Record missing name")?;

    // Determine read suffix (/1 or /2)
    let suffix: &[u8] = if no_suffix {
        b""
    } else if flags.is_first_segment() && !flags.is_last_segment() {
        b"/1"
    } else if flags.is_last_segment() && !flags.is_first_segment() {
        b"/2"
    } else {
        b"" // Single-end or both flags set
    };

    // Get sequence and quality - use concrete methods for direct slice access
    let seq = record.sequence();
    let qual = record.quality_scores();

    // Collect sequence bytes into reusable buffer
    seq_buf.clear();
    seq_buf.extend(seq.iter());

    // Collect and transform quality bytes using direct slice access (no Result unwrapping)
    qual_buf.clear();
    qual_buf.extend(qual.as_ref().iter().map(|&s| QUAL_TO_ASCII[s as usize]));

    if flags.is_reverse_complemented() {
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
    writer.write_all(name.as_ref())?;
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
        assert_eq!(parse_flags("2304").unwrap(), 2304);
        assert_eq!(parse_flags("0").unwrap(), 0);
        assert_eq!(parse_flags("65535").unwrap(), 65535);
    }

    #[test]
    fn test_parse_flags_hex() {
        assert_eq!(parse_flags("0x900").unwrap(), 0x900);
        assert_eq!(parse_flags("0X900").unwrap(), 0x900);
        assert_eq!(parse_flags("0xff").unwrap(), 0xff);
        assert_eq!(parse_flags("0xFFFF").unwrap(), 0xFFFF);
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
        write_reverse_complement_bytes(&mut output, b"ACGT").unwrap();
        assert_eq!(output, b"ACGT");

        output.clear();
        write_reverse_complement_bytes(&mut output, b"AAAA").unwrap();
        assert_eq!(output, b"TTTT");

        output.clear();
        // ATCG reversed = GCTA, complemented = CGAT
        write_reverse_complement_bytes(&mut output, b"ATCG").unwrap();
        assert_eq!(output, b"CGAT");

        output.clear();
        // Test with N
        write_reverse_complement_bytes(&mut output, b"ANCG").unwrap();
        assert_eq!(output, b"CGNT");
    }

    #[test]
    fn test_write_quality_bytes() {
        let mut output = Vec::new();
        // Quality 0 -> ASCII 33 ('!')
        // Quality 30 -> ASCII 63 ('?')
        write_quality_bytes(&mut output, &[0, 30, 40]).unwrap();
        assert_eq!(output, vec![33, 63, 73]);
    }

    #[test]
    fn test_write_reversed_quality_bytes() {
        let mut output = Vec::new();
        write_reversed_quality_bytes(&mut output, &[0, 30, 40]).unwrap();
        // Reversed: [40, 30, 0] -> [73, 63, 33]
        assert_eq!(output, vec![73, 63, 33]);
    }

    #[test]
    fn test_quality_encoding_edge_cases() {
        let mut output = Vec::new();
        // Test max valid quality (93)
        write_quality_bytes(&mut output, &[93]).unwrap();
        assert_eq!(output, vec![126]); // '~'

        output.clear();
        // Test overflow clamping (94+ should clamp to 126)
        write_quality_bytes(&mut output, &[94, 100, 255]).unwrap();
        assert_eq!(output, vec![126, 126, 126]);
    }
}
