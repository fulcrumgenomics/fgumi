//! FASTQ parsing stages — gzip variant ([`FastqParse`]) and BGZF
//! variant ([`FastqBlockParse`]).
//!
//! Both are parallel pool stages that turn decompressed chunk bytes
//! into parsed records. They differ in boundary handling:
//! `FastqParse` relies on pre-computed offsets attached to the
//! gzip-path chunk, while `FastqBlockParse` scans the block for
//! record boundaries with SIMD and emits prefix/suffix fragments for
//! cross-block record stitching.
//!
//! ## Input
//!
//! [`PerStreamChunk`] — decompressed bytes tagged by `stream_idx` and
//! `batch_num`. `FastqParse` requires `offsets` to be populated (by
//! the gzip-path `FastqFileRead`); `FastqBlockParse` requires
//! `offsets == None` and discovers boundaries itself.
//!
//! ## Output
//!
//! - `FastqParse` → [`FastqParsedStream`]: fully parsed records with
//!   per-stream ordering metadata.
//! - `FastqBlockParse` → [`BlockParsed`]: middle records plus prefix
//!   and suffix byte fragments that the downstream merge stitches
//!   across block boundaries.
//!
//! ## Ordering guarantees
//!
//! Each stage preserves the per-stream `(stream_idx, batch_num)` /
//! `(stream_idx, block_idx)` identifiers verbatim; within a chunk,
//! records are emitted in physical byte order. Cross-stream / cross-
//! chunk ordering is imposed by the downstream
//! [`crate::runall::engine::stages::fastq_pair::FastqPair`] or
//! [`crate::runall::engine::stages::fastq_pair::FastqBlockMerge`] stage.
//!
//! ## Memory model
//!
//! 1:1 with input: one output per input chunk/block. Record vectors
//! are sized from the discovered/provided offset count.
//!
//! ## Determinism
//!
//! Byte-identical: parsing is a pure function of the input bytes.
//! SIMD boundary detection produces the same offsets across runs.

use anyhow::Result;
use fgumi_simd_fastq::find_record_offsets;

use crate::runall::engine::fastq_types::{
    BlockParsed, FastqParsedStream, OwnedFastqRecord, PerStreamChunk,
};
use crate::runall::engine::stage::{Parallelism, Stage};

/// Parses gzip-path chunks. Input has `offsets` pre-computed by `FastqFileRead`.
#[derive(Clone, Default)]
pub struct FastqParse;

impl Stage for FastqParse {
    type Input = PerStreamChunk;
    type Output = FastqParsedStream;

    #[tracing::instrument(name = "fastq_parse", skip_all)]
    fn process(&mut self, input: Self::Input, out: &mut dyn FnMut(Self::Output)) -> Result<()> {
        let offsets = input
            .offsets
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("FastqParse requires gzip-path chunk with offsets"))?;
        // Sentinel convention: offsets has len = records + 1 with leading 0
        // and trailing data.len(). Iterating windows(2) yields (start, end)
        // pairs for every record.
        let mut records = Vec::with_capacity(offsets.len().saturating_sub(1));
        for w in offsets.windows(2) {
            let (start, end) = (w[0], w[1]);
            if start == end {
                continue;
            }
            records.push(parse_one_record_owned(&input.data[start..end])?);
        }
        out(FastqParsedStream {
            batch_num: input.batch_num,
            stream_idx: input.stream_idx,
            records,
            is_last: input.is_last,
        });
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }
    fn output_memory_estimate(&self, out: &Self::Output) -> usize {
        out.records.iter().map(|r| r.name.len() + r.sequence.len() + r.quality.len()).sum()
    }
    fn name(&self) -> &'static str {
        "FastqParse"
    }
}

/// BGZF path: finds boundaries inside a decompressed block and parses.
///
/// Emits prefix/suffix fragments for cross-block record handling.
#[derive(Clone, Default)]
pub struct FastqBlockParse;

impl Stage for FastqBlockParse {
    type Input = PerStreamChunk;
    type Output = BlockParsed;

    #[tracing::instrument(name = "fastq_block_parse", skip_all)]
    fn process(&mut self, input: Self::Input, out: &mut dyn FnMut(Self::Output)) -> Result<()> {
        // EOF sentinel: empty data + is_last=true passes through.
        if input.is_last && input.data.is_empty() {
            out(BlockParsed {
                block_idx: input.batch_num,
                stream_idx: input.stream_idx,
                records: Vec::new(),
                prefix_bytes: Vec::new(),
                suffix_bytes: Vec::new(),
                is_last: true,
            });
            return Ok(());
        }

        let data = &input.data;
        let prefix_end = detect_prefix_end(data);
        let suffix_start = if prefix_end >= data.len() {
            data.len()
        } else {
            detect_suffix_start(&data[prefix_end..]) + prefix_end
        };

        let prefix_bytes = data[..prefix_end].to_vec();
        let suffix_bytes = data[suffix_start..].to_vec();
        let middle = &data[prefix_end..suffix_start];

        // Parse complete records from the middle using SIMD-accelerated offset detection.
        let offsets = find_record_offsets(middle);
        let num_records = offsets.len().saturating_sub(1);
        let mut records = Vec::with_capacity(num_records);
        for i in 0..num_records {
            let start = offsets[i];
            let end = offsets[i + 1];
            if start >= end || start >= middle.len() {
                continue;
            }
            let end = end.min(middle.len());
            let record = parse_one_record_owned(&middle[start..end])?;
            records.push(record);
        }

        out(BlockParsed {
            block_idx: input.batch_num,
            stream_idx: input.stream_idx,
            records,
            prefix_bytes,
            suffix_bytes,
            is_last: input.is_last,
        });
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }
    fn output_memory_estimate(&self, out: &Self::Output) -> usize {
        out.records.iter().map(|r| r.name.len() + r.sequence.len() + r.quality.len()).sum::<usize>()
            + out.prefix_bytes.len()
            + out.suffix_bytes.len()
    }
    fn name(&self) -> &'static str {
        "FastqBlockParse"
    }
}

/// Detect where the first complete FASTQ record starts in `data`.
///
/// Returns the byte offset where the first complete record begins (bytes
/// `data[..prefix_end]` form an incomplete record fragment carried over from
/// the previous block). Returns `data.len()` if no complete record boundary
/// can be detected — the whole buffer is treated as a prefix fragment.
///
/// Mirrors `src/lib/unified_pipeline/fastq.rs::detect_prefix_end`.
fn detect_prefix_end(data: &[u8]) -> usize {
    if data.is_empty() {
        return 0;
    }

    // Collect the first 8 newline positions.
    let mut newlines = [0usize; 8];
    let mut count = 0;
    for (i, &b) in data.iter().enumerate() {
        if b == b'\n' {
            newlines[count] = i;
            count += 1;
            if count == 8 {
                break;
            }
        }
    }

    if count < 4 {
        return data.len();
    }

    for start in 0..count.saturating_sub(3) {
        let line0_start = if start == 0 { 0 } else { newlines[start - 1] + 1 };
        let line0_end = newlines[start];
        let line1_end = newlines[start + 1];
        let line2_end = newlines[start + 2];
        let line3_end = newlines[start + 3];

        if data[line0_start] != b'@' {
            continue;
        }
        let line2_start = line1_end + 1;
        if line2_start >= data.len() || data[line2_start] != b'+' {
            continue;
        }
        let seq_len = line1_end - (line0_end + 1);
        let qual_len = line3_end - (line2_end + 1);
        if seq_len != qual_len {
            continue;
        }
        return line0_start;
    }

    data.len()
}

/// Detect where the last complete FASTQ record ends in `data`.
///
/// Scans backward for newlines and tests sliding windows against FASTQ
/// invariants. Returns `data.len()` if all data forms complete records;
/// returns `0` if no complete trailing record can be located.
///
/// Mirrors `src/lib/unified_pipeline/fastq.rs::detect_suffix_start`.
fn detect_suffix_start(data: &[u8]) -> usize {
    if data.is_empty() {
        return 0;
    }

    let mut newlines = [0usize; 8];
    let mut count = 0;
    let mut i = data.len();
    while i > 0 {
        i -= 1;
        if data[i] == b'\n' {
            newlines[count] = i;
            count += 1;
            if count == 8 {
                break;
            }
        }
    }

    if count < 4 {
        return 0;
    }

    newlines[..count].reverse();

    let window_start = count.saturating_sub(4);
    for start in (0..=window_start).rev() {
        if start + 3 >= count {
            continue;
        }
        let line0_start = if start == 0 { 0 } else { newlines[start - 1] + 1 };
        let line0_end = newlines[start];
        let line1_end = newlines[start + 1];
        let line2_end = newlines[start + 2];
        let line3_end = newlines[start + 3];

        if data[line0_start] != b'@' {
            continue;
        }
        let line2_start = line1_end + 1;
        if line2_start >= data.len() || data[line2_start] != b'+' {
            continue;
        }
        let seq_len = line1_end - (line0_end + 1);
        let qual_len = line3_end - (line2_end + 1);
        if seq_len != qual_len {
            continue;
        }
        return line3_end + 1;
    }

    0
}

/// Parse a single FASTQ record from a byte slice. The slice must contain
/// exactly 4 lines: `@name\n`, `sequence\n`, `+\n`, `quality\n`. A trailing
/// newline is optional on the quality line.
fn parse_one_record_owned(bytes: &[u8]) -> Result<OwnedFastqRecord> {
    let mut lines = bytes.split(|&b| b == b'\n');
    let name_line = lines.next().ok_or_else(|| anyhow::anyhow!("missing name line"))?;
    let seq_line = lines.next().ok_or_else(|| anyhow::anyhow!("missing sequence line"))?;
    let _plus = lines.next().ok_or_else(|| anyhow::anyhow!("missing + line"))?;
    let qual_line = lines.next().ok_or_else(|| anyhow::anyhow!("missing quality line"))?;
    let name = name_line
        .strip_prefix(b"@")
        .ok_or_else(|| anyhow::anyhow!("name line missing @ prefix: {name_line:?}"))?;
    Ok(OwnedFastqRecord {
        name: name.to_vec(),
        sequence: seq_line.to_vec(),
        quality: qual_line.to_vec(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk_chunk_gzip(data: Vec<u8>, offsets: Vec<usize>, is_last: bool) -> PerStreamChunk {
        PerStreamChunk { stream_idx: 0, batch_num: 0, data, offsets: Some(offsets), is_last }
    }

    fn mk_chunk_bgzf(data: Vec<u8>, is_last: bool) -> PerStreamChunk {
        PerStreamChunk { stream_idx: 0, batch_num: 0, data, offsets: None, is_last }
    }

    #[test]
    fn test_parse_one_record_owned() {
        let bytes: &[u8] = b"@r1\nACGT\n+\nIIII\n";
        let rec = parse_one_record_owned(bytes).unwrap();
        assert_eq!(rec.name, b"r1");
        assert_eq!(rec.sequence, b"ACGT");
        assert_eq!(rec.quality, b"IIII");
    }

    #[test]
    fn test_fastq_parse_gzip_path_two_records() {
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nCGTA\n+\nJJJJ\n".to_vec();
        // Sentinel convention: offsets has len = records + 1, with leading 0
        // and trailing data.len().
        let second_start = data.windows(3).position(|w| w == b"@r2").unwrap();
        let data_len = data.len();
        let chunk = mk_chunk_gzip(data, vec![0, second_start, data_len], true);
        let mut stage = FastqParse;
        let mut captured = None;
        stage
            .process(chunk, &mut |v| {
                captured = Some(v);
            })
            .unwrap();
        let out = captured.expect("stage must emit");
        assert_eq!(out.records.len(), 2);
        assert_eq!(&out.records[0].name, b"r1");
        assert_eq!(&out.records[1].name, b"r2");
    }

    #[test]
    fn test_fastq_block_parse_complete_records() {
        let data = b"@r1\nACGT\n+\nIIII\n@r2\nCGTA\n+\nJJJJ\n".to_vec();
        let chunk = mk_chunk_bgzf(data, true);
        let mut stage = FastqBlockParse;
        let mut captured = None;
        stage
            .process(chunk, &mut |v| {
                captured = Some(v);
            })
            .unwrap();
        let out = captured.expect("stage must emit");
        assert_eq!(out.records.len(), 2);
        assert_eq!(out.prefix_bytes, Vec::<u8>::new());
        assert_eq!(out.suffix_bytes, Vec::<u8>::new());
    }

    #[test]
    fn test_fastq_block_parse_with_prefix_fragment() {
        // Prefix fragment (tail of previous block's quality line) + 1 complete record.
        // The fragment "GT\n+\nIIII\n" looks like seq/plus/qual of a previous record;
        // our detector walks the first 8 newlines looking for a valid 4-line window
        // and locates the real record start.
        let mut data = b"GT\n+\nIIII\n".to_vec();
        data.extend_from_slice(b"@r1\nACGT\n+\nIIII\n");
        let prefix_len = b"GT\n+\nIIII\n".len();
        let chunk = mk_chunk_bgzf(data, true);
        let mut stage = FastqBlockParse;
        let mut captured = None;
        stage
            .process(chunk, &mut |v| {
                captured = Some(v);
            })
            .unwrap();
        let out = captured.expect("stage must emit");
        assert_eq!(out.records.len(), 1);
        assert_eq!(&out.records[0].name, b"r1");
        assert_eq!(out.prefix_bytes.len(), prefix_len);
        assert_eq!(out.suffix_bytes, Vec::<u8>::new());
    }

    #[test]
    fn test_fastq_block_parse_with_suffix_fragment() {
        // One complete record, then a partial record missing the + and quality lines.
        let mut data = b"@r1\nACGT\n+\nIIII\n".to_vec();
        data.extend_from_slice(b"@r2\nCGTA\n");
        let chunk = mk_chunk_bgzf(data, false);
        let mut stage = FastqBlockParse;
        let mut captured = None;
        stage
            .process(chunk, &mut |v| {
                captured = Some(v);
            })
            .unwrap();
        let out = captured.expect("stage must emit");
        assert_eq!(out.records.len(), 1);
        assert_eq!(&out.records[0].name, b"r1");
        // The partial "@r2\nCGTA\n" becomes the suffix fragment.
        assert_eq!(out.suffix_bytes, b"@r2\nCGTA\n".to_vec());
        assert_eq!(out.prefix_bytes, Vec::<u8>::new());
    }
}
