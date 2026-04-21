//! BAM file source that emits `SerializedBatch` items (not per-record).
//!
//! Reads a BGZF-compressed BAM file, batches records into `SerializedBatch`
//! chunks of at most [`MAX_RECORDS_PER_BATCH`], and emits them. Each
//! `SerializedBatch.primary.data` is a sequence of length-prefixed BAM records
//! in the same format emitted by `ExtractStage` and consumed by `CorrectStage`
//! and `BgzfCompress`.
//!
//! This source is the `--start-from correct` entry point for `pipeline`:
//! it bridges a BAM file to the Correct → Compress → Write pipeline without
//! needing a separate batching stage.

use std::path::PathBuf;

use anyhow::{Context, Result};
use fgumi_raw_bam::RawRecord;
use noodles::sam::Header;

use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::output_types::{RawBytes, SerializedBatch};
use crate::runall::engine::source::{OutputQueue, Source};
use crate::runall::engine::stage::SequencedItem;

/// Maximum records per emitted `SerializedBatch`. Matches the cap used by
/// `FastqBlockMerge` / `FastqPair` so downstream stages see consistent
/// batch sizes regardless of entry point.
const MAX_RECORDS_PER_BATCH: usize = 256;

/// Capacity hint for the per-batch byte buffer. Typical short-read BAM records
/// run ~150-300 bytes after decoding, so a 64 KiB buffer covers `MAX_RECORDS_PER_BATCH`
/// without reallocating in the common case.
const BATCH_BUFFER_CAPACITY: usize = 64 * 1024;

/// BAM file source emitting length-prefixed record batches.
pub struct BamBatchedSource {
    input_path: PathBuf,
    threads: usize,
}

impl BamBatchedSource {
    /// Construct a new source pointing at a BAM file.
    #[must_use]
    pub fn new(input_path: PathBuf, threads: usize) -> Self {
        Self { input_path, threads }
    }

    /// Read and return the BAM header. Opens the file to parse the header
    /// but does not consume records.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened or the header cannot
    /// be parsed.
    pub fn read_header(&self) -> Result<Header> {
        let (_reader, header) =
            crate::bam_io::create_raw_bam_reader(&self.input_path, self.threads)
                .with_context(|| format!("failed to open BAM {}", self.input_path.display()))?;
        Ok(header)
    }
}

impl Source for BamBatchedSource {
    type Output = SerializedBatch;

    fn run(
        self: Box<Self>,
        output: Box<dyn OutputQueue<Self::Output>>,
        cancel: CancelToken,
    ) -> Result<()> {
        let this = *self;
        let (mut reader, _header) =
            crate::bam_io::create_raw_bam_reader(&this.input_path, this.threads)
                .with_context(|| format!("failed to open BAM {}", this.input_path.display()))?;

        let mut record = RawRecord::new();
        let mut ordinal: u64 = 0;
        let mut batch_data: Vec<u8> = Vec::with_capacity(BATCH_BUFFER_CAPACITY);
        let mut batch_count: u64 = 0;
        let mut total_records: u64 = 0;

        loop {
            if cancel.is_cancelled() {
                break;
            }
            let n = reader.read_record(&mut record)?;
            if n == 0 {
                break;
            }

            let record_bytes = record.as_ref();
            let block_size = u32::try_from(record_bytes.len())
                .context("BamBatchedSource: record too large for u32 block_size")?;
            batch_data.extend_from_slice(&block_size.to_le_bytes());
            batch_data.extend_from_slice(record_bytes);
            batch_count += 1;

            if batch_count >= MAX_RECORDS_PER_BATCH as u64 {
                let next_buffer = Vec::with_capacity(BATCH_BUFFER_CAPACITY);
                let batch = SerializedBatch {
                    primary: RawBytes {
                        data: std::mem::replace(&mut batch_data, next_buffer),
                        record_count: batch_count,
                    },
                    secondary: None,
                    ordinal,
                };
                let mem = batch.primary.data.len();
                let item = SequencedItem::new(ordinal, batch, mem);
                if output.push_until_cancelled(item, &cancel).is_err() {
                    output.close();
                    return Ok(());
                }
                crate::progress::records_read(batch_count);
                total_records += batch_count;
                ordinal += 1;
                batch_count = 0;
            }
        }

        // Emit the trailing partial batch (if any).
        if batch_count > 0 {
            let batch = SerializedBatch {
                primary: RawBytes { data: batch_data, record_count: batch_count },
                secondary: None,
                ordinal,
            };
            let mem = batch.primary.data.len();
            let item = SequencedItem::new(ordinal, batch, mem);
            if output.push_until_cancelled(item, &cancel).is_ok() {
                crate::progress::records_read(batch_count);
                total_records += batch_count;
            }
        }

        if !cancel.is_cancelled() {
            crate::progress::reads_complete(total_records);
        }

        output.close();
        Ok(())
    }

    fn name(&self) -> &'static str {
        "BamBatchedSource"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::memory::MemoryTracker;
    use crate::runall::engine::queue::StageQueue;
    use std::sync::Arc;

    fn write_empty_bam(path: &std::path::Path) {
        use noodles::sam::header::record::value::Map;
        let header = noodles::sam::Header::builder()
            .set_header(Map::<noodles::sam::header::record::value::map::Header>::new(
                noodles::sam::header::record::value::map::header::Version::new(1, 6),
            ))
            .build();
        let writer = crate::bam_io::create_raw_bam_writer(path, &header, 1, 6).unwrap();
        writer.finish().unwrap();
    }

    #[test]
    fn test_bam_batched_source_empty_file() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_empty_bam(tmp.path());

        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let queue: Arc<StageQueue<SerializedBatch>> =
            Arc::new(StageQueue::new("bam_batched_out", 16, 100_000, tracker));

        let source = Box::new(BamBatchedSource::new(tmp.path().to_path_buf(), 1));
        source.run(Box::new(queue.clone()), CancelToken::new()).unwrap();

        assert!(queue.is_closed());
        assert!(queue.pop().is_none());
    }

    #[test]
    fn test_bam_batched_source_name() {
        let source = BamBatchedSource::new(PathBuf::from("/tmp/nonexistent.bam"), 1);
        assert_eq!(Source::name(&source), "BamBatchedSource");
    }
}
