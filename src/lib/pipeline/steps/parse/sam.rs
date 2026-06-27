//! `ParseSamChunk` mid-step. `Parallel + ByItemOrdinal`. Consumes a
//! `SamChunk` from the SAM source, parses each line via noodles' SAM
//! reader, re-encodes the record into its BAM body bytes via
//! `bam::io::Writer`, computes a `GroupKey`, and emits a
//! `DecodedRecordBatch` — the same output type `DecodeRecords`
//! produces for the BAM chain.
//!
//! Each chunk is independent (the upstream Read step splits at line
//! boundaries, so records never cross chunks). Workers process
//! different chunks concurrently — this is the parallel-parse phase
//! that lets the SAM ingest path keep up with the multi-MB/s aligner
//! output stream.

use std::io;
use std::sync::Arc;

use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::io::Write as _;

use crate::pipeline::core::Unpushed;
use crate::pipeline::core::held::HeldSlot;
use crate::pipeline::core::outputs::OrderedBytesSingle;
use crate::pipeline::core::queues::QueueSpec;
use crate::pipeline::core::reorder::BranchOrdering;
use crate::pipeline::core::step::{Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::{DecodedRecordBatch, SamChunk};
use fgumi_bam_io::{DecodedRecord, GroupKeyConfig, compute_group_key_from_raw, name_hash_key};

/// Parse every line in `chunk` and produce a `DecodedRecordBatch` carrying
/// the chunk's batch serial. SAM lines are encoded into BAM record body
/// bytes via `bam::io::Writer` (which writes `u32 LE block_size` +
/// record body to its inner buffer); the 4-byte length prefix is
/// stripped so the resulting `DecodedRecord` matches the shape produced
/// by the BAM-side `DecodeRecords`.
///
/// # Errors
///
/// Returns I/O errors from SAM line parsing or BAM encoding.
pub fn parse_sam_chunk_into_decoded(
    chunk: SamChunk,
    header: &sam::Header,
    key_config: &GroupKeyConfig,
) -> io::Result<DecodedRecordBatch> {
    let SamChunk { batch_serial, bytes, line_offsets } = chunk;
    let n_records = line_offsets.len().saturating_sub(1);
    let mut decoded: Vec<DecodedRecord> = Vec::with_capacity(n_records);

    let library_index = &key_config.library_index;
    let cell_tag = key_config.cell_tag;
    let name_hash_only = key_config.name_hash_only;

    let mut sam_record = sam::Record::default();
    let mut encoder = bam::io::Writer::from(Vec::<u8>::with_capacity(4096));

    // Walk lines explicitly so we can skip empty ones (e.g. a trailing
    // `\n\n` in the input shows up as a length-1 empty offset slot;
    // feeding `\n` to noodles' `sam::io::Reader` errors with
    // "unexpected EOL"). Each non-empty line is read with a fresh
    // single-line reader so we don't have to coordinate state across
    // skipped lines.
    for i in 0..n_records {
        let start = line_offsets[i] as usize;
        let end = line_offsets[i + 1] as usize;
        // Length 1 means just the `\n`; length 0 would be an
        // ill-formed offset (defensive — shouldn't happen with the
        // current splitter, but skipping is the right behavior either way).
        if end.saturating_sub(start) <= 1 {
            continue;
        }
        let line = &bytes[start..end];
        let mut sam_reader = sam::io::Reader::new(line);
        let n = sam_reader.read_record(&mut sam_record)?;
        if n == 0 {
            continue;
        }
        encoder.get_mut().clear();
        encoder.write_alignment_record(header, &sam_record)?;
        // `DecodedRecord` takes ownership of the body bytes, so swap the
        // encoder's inner buffer out (replacing it with a pre-sized fresh one
        // for the next record) and strip the 4-byte little-endian `block_size`
        // prefix that `bam::io::Writer` prepends in place. This avoids the extra
        // full-body memcpy that `get_ref()[4..].to_vec()` performed per record
        // (this is the hot SAM parse loop). The replacement is sized to the
        // outgoing buffer's capacity so the encoder keeps room for the next
        // record instead of re-growing a zero-capacity `Vec` from scratch (which
        // `std::mem::take` would leave behind).
        let next_capacity = encoder.get_mut().capacity().max(4096);
        let mut body_bytes =
            std::mem::replace(encoder.get_mut(), Vec::with_capacity(next_capacity));
        body_bytes.drain(..4);
        // Match the BAM decode path (`DecodeFromRecords`): under
        // `name_hash_only` the key is the name hash alone, so SAM- and
        // BAM-ingested records group identically.
        let key = if name_hash_only {
            name_hash_key(&body_bytes)
        } else {
            compute_group_key_from_raw(&body_bytes, library_index, cell_tag)
        };
        decoded.push(DecodedRecord::from_raw_bytes(body_bytes, key));
    }

    Ok(DecodedRecordBatch::new(batch_serial, decoded))
}

/// `Parallel + ByItemOrdinal` SAM parser. `SamChunk → DecodedRecordBatch`.
/// Sibling of `parse::decode::DecodeRecords` for the SAM ingest path —
/// converging at `DecodedRecordBatch` so every downstream step
/// (`GroupBam`/`GroupByPosition`/`process_ordered`/…) is reused
/// unchanged.
pub struct ParseSamChunk {
    /// Reference header. Shared across all worker clones — encoding
    /// needs it for reference-sequence ID resolution. `Arc` so cloning
    /// the step is cheap.
    header: Arc<sam::Header>,
    /// Group-key config (library index + cell tag). Shared via the
    /// `Arc<LibraryIndex>` inside `GroupKeyConfig`.
    key_config: GroupKeyConfig,
    held: HeldSlot<Unpushed<DecodedRecordBatch>>,
    output_byte_limit: u64,
}

impl ParseSamChunk {
    #[must_use]
    pub fn new(
        header: Arc<sam::Header>,
        key_config: GroupKeyConfig,
        output_byte_limit: u64,
    ) -> Self {
        Self { header, key_config, held: HeldSlot::new(), output_byte_limit }
    }
}

impl Clone for ParseSamChunk {
    fn clone(&self) -> Self {
        Self {
            header: Arc::clone(&self.header),
            key_config: self.key_config.clone(),
            held: HeldSlot::new(),
            output_byte_limit: self.output_byte_limit,
        }
    }
}

impl Step for ParseSamChunk {
    type Input = SamChunk;
    type Outputs = OrderedBytesSingle<DecodedRecordBatch>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "ParseSamChunk",
            kind: StepKind::Parallel,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: self.output_byte_limit }],
            branch_ordering: vec![BranchOrdering::ByItemOrdinal],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        if let Some(unpushed) = self.held.take() {
            match ctx.outputs.retry(unpushed) {
                Ok(()) => {}
                Err(again) => {
                    self.held.put(again);
                    return Ok(StepOutcome::Contention);
                }
            }
        }

        let Some(chunk) = ctx.input.pop() else {
            // No input this call. If upstream is drained, every item has been
            // processed (held output was flushed by the Contention preamble
            // above) and this step will never push again — report Finished.
            // For a Parallel step only the last clone to finish closes the
            // shared output (gated by the StepDrainCounter in the driver).
            if ctx.input.is_drained() {
                return Ok(StepOutcome::Finished);
            }
            return Ok(StepOutcome::NoProgress);
        };

        let batch = parse_sam_chunk_into_decoded(chunk, &self.header, &self.key_config)?;

        match ctx.outputs.push(batch) {
            Ok(()) => Ok(StepOutcome::Progress),
            Err(unpushed) => {
                self.held.put(unpushed);
                Ok(StepOutcome::Progress)
            }
        }
    }

    fn new_worker_copy(&self) -> Self {
        self.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::steps::types::SamChunk;
    use fgumi_bam_io::GroupKeyConfig;
    use noodles::sam;

    const SAM_TEXT: &str = "\
@HD\tVN:1.6\tSO:unsorted\n\
@SQ\tSN:chr1\tLN:1000\n\
read1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n\
read2\t16\tchr1\t20\t60\t5M\t*\t0\t0\tTGCAT\t!!!!!\n";

    fn parse_header(text: &str) -> sam::Header {
        let mut reader = sam::io::Reader::new(text.as_bytes());
        reader.read_header().expect("header parses")
    }

    fn key_config_for(header: &sam::Header) -> GroupKeyConfig {
        let library_index = fgumi_bam_io::LibraryIndex::from_header(header);
        GroupKeyConfig::new_raw_no_cell(library_index)
    }

    /// Build a `SamChunk` from the two record lines in `SAM_TEXT` (skipping
    /// the header lines). This is what `ReadSamChunks` would emit downstream
    /// of header consumption.
    fn sam_chunk_from_records(serial: u64) -> SamChunk {
        let records = "\
read1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n\
read2\t16\tchr1\t20\t60\t5M\t*\t0\t0\tTGCAT\t!!!!!\n";
        let bytes = records.as_bytes().to_vec();
        // Two lines: offsets at [0, end_of_line_0, end_of_line_1].
        let mut line_offsets = vec![0u32];
        for (i, &b) in bytes.iter().enumerate() {
            if b == b'\n' {
                line_offsets.push(u32::try_from(i + 1).unwrap());
            }
        }
        SamChunk { batch_serial: serial, bytes, line_offsets }
    }

    #[test]
    fn parse_sam_chunk_emits_one_decoded_record_per_input_line() {
        let header = parse_header(SAM_TEXT);
        let key_config = key_config_for(&header);
        let chunk = sam_chunk_from_records(7);

        let batch = parse_sam_chunk_into_decoded(chunk, &header, &key_config).expect("parses");
        assert_eq!(batch.batch_serial, 7, "serial propagated from input chunk");
        assert_eq!(batch.records.len(), 2);
    }

    #[test]
    fn parse_sam_chunk_emits_records_with_correct_read_names() {
        // Round-trip test: SAM text → BAM body bytes → decode names. Proves
        // the encode step preserves identity (per SAM spec, BAM record body
        // layout: refID/pos/l_read_name/.../read_name@offset 32).
        let header = parse_header(SAM_TEXT);
        let key_config = key_config_for(&header);
        let chunk = sam_chunk_from_records(0);

        let batch = parse_sam_chunk_into_decoded(chunk, &header, &key_config).expect("parses");

        let mut names: Vec<Vec<u8>> = Vec::new();
        for rec in &batch.records {
            let body = rec.raw_bytes();
            // BAM record body layout: bytes 0..32 are fixed fields;
            // byte 8 is l_read_name (including NUL); read_name starts at
            // byte 32 and is `l_read_name - 1` bytes (excluding NUL).
            let l_read_name = body[8] as usize;
            names.push(body[32..32 + l_read_name - 1].to_vec());
        }
        assert_eq!(names, vec![b"read1".to_vec(), b"read2".to_vec()]);
    }

    #[test]
    fn parse_sam_chunk_skips_empty_lines() {
        // A SAM file ending in `record\n\n` produces offsets that include an
        // empty line slot. The parser must skip such empty offsets — feeding
        // an empty line to noodles' SAM reader raises "unexpected EOL".
        let header = parse_header(SAM_TEXT);
        let key_config = key_config_for(&header);
        // Construct a chunk with a real record followed by an empty line.
        let records = b"read1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n\n";
        let bytes = records.to_vec();
        let mut line_offsets = vec![0u32];
        for (i, &b) in bytes.iter().enumerate() {
            if b == b'\n' {
                line_offsets.push(u32::try_from(i + 1).unwrap());
            }
        }
        // Three offsets => two "lines", second is empty (the `\n` after `\n`).
        assert_eq!(line_offsets.len(), 3);
        let chunk = SamChunk { batch_serial: 0, bytes, line_offsets };

        let batch = parse_sam_chunk_into_decoded(chunk, &header, &key_config).expect("parses");
        assert_eq!(batch.records.len(), 1, "empty line should be skipped, not parsed");
    }

    #[test]
    fn parse_sam_chunk_with_empty_record_count_returns_empty_batch() {
        let header = parse_header(SAM_TEXT);
        let key_config = key_config_for(&header);
        let chunk = SamChunk { batch_serial: 11, bytes: Vec::new(), line_offsets: Vec::new() };

        let batch = parse_sam_chunk_into_decoded(chunk, &header, &key_config).expect("parses");
        assert_eq!(batch.batch_serial, 11);
        assert!(batch.records.is_empty());
    }

    /// Full-field parity: BAM-encode a record via noodles + `DecodeRecords`
    /// produces the same `DecodedRecord.raw_bytes()` as SAM-encode of the
    /// equivalent text via `ParseSamChunk`. This is the strongest correctness
    /// assertion for the SAM parse path — every BAM record field (flags,
    /// refID, pos, MAPQ, CIGAR, sequence, quality, mate info, tags) must
    /// round-trip identically.
    #[test]
    fn parse_sam_chunk_record_bytes_match_bam_decode_path() {
        use fgumi_bam_io::DecodedRecord as _DecodedRecord;
        use noodles::bam;
        use noodles::sam::alignment::io::Write as _;

        // Header with one ref + a read group.
        let header_text = "\
@HD\tVN:1.6\tSO:coordinate\n\
@SQ\tSN:chr1\tLN:1000\n\
@RG\tID:rg1\tSM:sample\n";
        let header = parse_header(header_text);
        let key_config = key_config_for(&header);

        // Record with the full menagerie of fields: mapped, paired with mate,
        // MAPQ 60, multi-op CIGAR (M/I/D/S/H), full sequence + qualities,
        // and several optional tag types (Z string, i int, B array, A char).
        // `RX` is a UMI tag the group-key extractor reads (relevant via
        // GroupKeyConfig).
        let sam_line = "rec1\t99\tchr1\t100\t60\t3M1I2M1D2S1H\t=\t200\t300\tACGTGCAC\tIIIIIIII\tRX:Z:ACGT\tNM:i:2\tBQ:Z:GGGGGGGG\tBI:B:S,10,20,30\n";

        // ── BAM path: parse to noodles record, encode body bytes via the
        // same encoder DecodeRecords would consume after BgzfDecompress +
        // FindBamBoundaries.
        let bam_body = {
            let mut reader = noodles::sam::io::Reader::new(sam_line.as_bytes());
            let mut rec = noodles::sam::Record::default();
            reader.read_record(&mut rec).expect("sam parse");
            let mut buf = bam::io::Writer::from(Vec::<u8>::new());
            buf.write_alignment_record(&header, &rec).expect("encode");
            // Strip the 4-byte block_size prefix.
            buf.get_ref()[4..].to_vec()
        };
        let bam_decoded = _DecodedRecord::from_raw_bytes(
            bam_body.clone(),
            fgumi_bam_io::compute_group_key_from_raw(
                &bam_body,
                &key_config.library_index,
                key_config.cell_tag,
            ),
        );

        // ── SAM path: same line through ParseSamChunk.
        let bytes = sam_line.as_bytes().to_vec();
        let mut line_offsets = vec![0u32];
        for (i, &b) in bytes.iter().enumerate() {
            if b == b'\n' {
                line_offsets.push(u32::try_from(i + 1).unwrap());
            }
        }
        let chunk = SamChunk { batch_serial: 0, bytes, line_offsets };
        let batch = parse_sam_chunk_into_decoded(chunk, &header, &key_config).expect("parses");
        assert_eq!(batch.records.len(), 1);
        let sam_decoded = &batch.records[0];

        // Byte-identical raw record bytes (the SAM and BAM encoder paths
        // both go through `bam::io::Writer::write_alignment_record`).
        assert_eq!(
            sam_decoded.raw_bytes(),
            bam_decoded.raw_bytes(),
            "SAM and BAM paths must produce byte-identical DecodedRecord.raw_bytes()",
        );
    }

    /// Under `name_hash_only`, the SAM parse path must key records by the read
    /// name hash alone — identical to the BAM `DecodeFromRecords` path — so
    /// SAM- and BAM-ingested reads group together. Regression for the SAM path
    /// previously ignoring `name_hash_only` and always computing the full
    /// position/library key.
    #[test]
    fn parse_sam_chunk_honors_name_hash_only() {
        let header = parse_header(SAM_TEXT);
        let library_index = fgumi_bam_io::LibraryIndex::from_header(&header);
        let key_config = GroupKeyConfig::name_hash_only(library_index);

        let chunk = sam_chunk_from_records(0);
        let batch = parse_sam_chunk_into_decoded(chunk, &header, &key_config).expect("parses");
        assert_eq!(batch.records.len(), 2);

        for rec in &batch.records {
            let raw = rec.raw_bytes();
            assert_eq!(
                rec.key,
                fgumi_bam_io::name_hash_key(raw),
                "name_hash_only SAM key must be the name-hash key, not the full group key",
            );
            // The fixture reads are mapped (chr1:10 / chr1:20), so the full
            // position/library key differs from the name-hash key — confirming
            // the branch actually changed paths.
            assert_ne!(
                rec.key,
                fgumi_bam_io::compute_group_key_from_raw(
                    raw,
                    &key_config.library_index,
                    key_config.cell_tag,
                ),
                "name_hash_only key must differ from the full position/library key for a mapped read",
            );
        }
    }

    #[test]
    fn parse_sam_chunk_step_profile_is_parallel_byordinal() {
        let header = std::sync::Arc::new(parse_header(SAM_TEXT));
        let key_config = key_config_for(&header);
        let step = ParseSamChunk::new(header, key_config, 1024 * 1024);
        let profile = step.profile();
        assert_eq!(profile.name, "ParseSamChunk");
        assert_eq!(profile.kind, crate::pipeline::core::step::StepKind::Parallel);
        assert!(!profile.sticky);
        assert_eq!(
            profile.branch_ordering,
            vec![crate::pipeline::core::reorder::BranchOrdering::ByItemOrdinal]
        );
        assert!(matches!(
            profile.output_queues[0],
            crate::pipeline::core::queues::QueueSpec::ByteBounded { .. }
        ));
    }
}
