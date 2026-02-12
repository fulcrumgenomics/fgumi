//! Convert BAM to FASTQ format.
//!
//! This tool reads a BAM file and outputs interleaved FASTQ to stdout for piping to aligners.
//! Input should be queryname-sorted or template-coordinate sorted.

use anyhow::{Context, Result};
use clap::Parser;
use fgumi_lib::bam_io::create_raw_bam_reader;
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::sort::bam_fields;
use fgumi_lib::validation::validate_file_exists;
use fgumi_lib::vendored::RawRecord;
use log::info;
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

        let (mut raw_reader, _header) = create_raw_bam_reader(&self.input, self.threads)?;

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
        let mut raw_record = RawRecord::new();

        loop {
            let n = raw_reader.read_record(&mut raw_record).with_context(|| {
                format!("Failed to read BAM record after {total_records} records")
            })?;
            if n == 0 {
                break;
            }

            let bam = raw_record.as_ref();
            total_records += 1;

            let flags = bam_fields::flags(bam);

            // Filter by flags
            if (flags & self.exclude_flags) != 0 {
                continue;
            }
            if (flags & self.require_flags) != self.require_flags {
                continue;
            }

            // Get sequence length for batch tracking
            let seq_len = bam_fields::l_seq(bam) as usize;

            // Write FASTQ record
            write_fastq_record_raw(
                &mut writer,
                bam,
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

/// Write a single FASTQ record from raw BAM bytes to the writer.
#[inline]
fn write_fastq_record_raw<W: Write>(
    writer: &mut W,
    bam: &[u8],
    flags: u16,
    no_suffix: bool,
    seq_buf: &mut Vec<u8>,
    qual_buf: &mut Vec<u8>,
) -> Result<()> {
    // Get read name
    let name = bam_fields::read_name(bam);

    // Determine read suffix (/1 or /2)
    let suffix: &[u8] = if no_suffix {
        b""
    } else if (flags & bam_fields::flags::FIRST_SEGMENT) != 0
        && (flags & bam_fields::flags::LAST_SEGMENT) == 0
    {
        b"/1"
    } else if (flags & bam_fields::flags::LAST_SEGMENT) != 0
        && (flags & bam_fields::flags::FIRST_SEGMENT) == 0
    {
        b"/2"
    } else {
        b"" // Single-end or both flags set
    };

    // Decode packed 4-bit sequence into reusable buffer
    let seq_len = bam_fields::l_seq(bam) as usize;
    let seq_off = bam_fields::seq_offset(bam);
    seq_buf.clear();
    seq_buf.extend((0..seq_len).map(|i| {
        let code = bam_fields::get_base(bam, seq_off, i);
        bam_fields::BAM_BASE_TO_ASCII[code as usize]
    }));

    // Transform quality scores using direct slice access
    let qual_scores = bam_fields::quality_scores_slice(bam);
    qual_buf.clear();
    qual_buf.extend(qual_scores.iter().map(|&s| QUAL_TO_ASCII[s as usize]));

    if (flags & bam_fields::flags::REVERSE) != 0 {
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

    /// Create a minimal raw BAM record with sequence and quality for FASTQ testing.
    #[allow(clippy::cast_possible_truncation)]
    fn create_fastq_test_record(name: &str, seq: &[u8], quals: &[u8], flags: u16) -> Vec<u8> {
        let l_read_name = (name.len() + 1) as u8;
        let l_seq = seq.len() as u32;
        let mut record = vec![0u8; 32];

        // refID = -1 (unmapped)
        record[0..4].copy_from_slice(&(-1i32).to_le_bytes());
        // pos = -1
        record[4..8].copy_from_slice(&(-1i32).to_le_bytes());
        // l_read_name
        record[8] = l_read_name;
        // flags
        record[14..16].copy_from_slice(&flags.to_le_bytes());
        // l_seq
        record[16..20].copy_from_slice(&l_seq.to_le_bytes());
        // next_refID = -1
        record[20..24].copy_from_slice(&(-1i32).to_le_bytes());
        // next_pos = -1
        record[24..28].copy_from_slice(&(-1i32).to_le_bytes());

        // Read name + null terminator
        record.extend_from_slice(name.as_bytes());
        record.push(0);

        // No CIGAR (n_cigar_op = 0)

        // Pack sequence (4-bit encoding, 2 bases per byte)
        bam_fields::pack_sequence_into(&mut record, seq);

        // Quality scores (raw Phred, not +33)
        record.extend_from_slice(quals);

        record
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

    #[test]
    fn test_write_fastq_record_raw_forward() {
        let record = create_fastq_test_record("read1", b"ACGT", &[30, 25, 20, 15], 0x41); // paired, first segment
        let mut output = Vec::new();
        let mut seq_buf = Vec::new();
        let mut qual_buf = Vec::new();

        write_fastq_record_raw(
            &mut output,
            &record,
            bam_fields::flags(&record),
            false,
            &mut seq_buf,
            &mut qual_buf,
        )
        .unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();
        assert_eq!(lines.len(), 4);
        assert_eq!(lines[0], "@read1/1");
        assert_eq!(lines[1], "ACGT");
        assert_eq!(lines[2], "+");
        // quals [30,25,20,15] → Phred+33 → [63,58,53,48] → "?:50"
        assert_eq!(lines[3], "?:50");
    }

    #[test]
    fn test_write_fastq_record_raw_reverse_complement() {
        // ACGT reverse complemented = ACGT (palindrome)
        let record = create_fastq_test_record("read1", b"AACG", &[30, 25, 20, 15], 0x10); // reverse
        let mut output = Vec::new();
        let mut seq_buf = Vec::new();
        let mut qual_buf = Vec::new();

        write_fastq_record_raw(
            &mut output,
            &record,
            bam_fields::flags(&record),
            true, // no suffix
            &mut seq_buf,
            &mut qual_buf,
        )
        .unwrap();

        let result = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = result.lines().collect();
        assert_eq!(lines[0], "@read1");
        // AACG reversed = GCAA, complemented = CGTT
        assert_eq!(lines[1], "CGTT");
    }

    #[test]
    fn test_write_fastq_record_raw_no_suffix() {
        let record = create_fastq_test_record("read1", b"ACGT", &[30, 25, 20, 15], 0x41);
        let mut output = Vec::new();
        let mut seq_buf = Vec::new();
        let mut qual_buf = Vec::new();

        write_fastq_record_raw(
            &mut output,
            &record,
            bam_fields::flags(&record),
            true, // no suffix
            &mut seq_buf,
            &mut qual_buf,
        )
        .unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.starts_with("@read1\n"));
    }

    #[test]
    fn test_write_fastq_record_raw_read2_suffix() {
        let record = create_fastq_test_record("read1", b"ACGT", &[30, 25, 20, 15], 0x81); // paired, last segment
        let mut output = Vec::new();
        let mut seq_buf = Vec::new();
        let mut qual_buf = Vec::new();

        write_fastq_record_raw(
            &mut output,
            &record,
            bam_fields::flags(&record),
            false,
            &mut seq_buf,
            &mut qual_buf,
        )
        .unwrap();

        let result = String::from_utf8(output).unwrap();
        assert!(result.starts_with("@read1/2\n"));
    }
}
