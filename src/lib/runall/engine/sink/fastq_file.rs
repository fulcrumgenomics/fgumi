//! Sequential FASTQ-file sink with reorder buffer.
//!
//! Consumes [`SerializedFastqBatch`] items in arbitrary ordinal order and
//! appends them to a single FASTQ file in strict ordinal order. FASTQ is
//! plain text with no header, so the sink simply appends each batch's
//! `data` buffer to the output file as the reorder buffer produces
//! in-order batches.
//!
//! Mirrors [`crate::runall::engine::sink::bam_file_write::BamFileWrite`] but
//! simpler: no compression, no secondary output, no EOF marker.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};

use crate::runall::engine::backoff::Backoff;
use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::output_types::SerializedFastqBatch;
use crate::runall::engine::reorder::ReorderBuffer;
use crate::runall::engine::sink::{InputQueue, Sink};

/// Sequential sink that reorders [`SerializedFastqBatch`] items by their
/// `ordinal` and appends their bytes to a plain FASTQ file.
pub struct FastqFileSink {
    path: PathBuf,
    /// Currently unused. Reserved for backpressure/queue tuning if the
    /// engine later allows sinks to publish their own capacity hints.
    _queue_capacity: usize,
}

impl FastqFileSink {
    /// Construct a new sink.
    #[must_use]
    pub fn new(path: PathBuf, queue_capacity: usize) -> Self {
        Self { path, _queue_capacity: queue_capacity }
    }
}

impl Sink for FastqFileSink {
    type Input = SerializedFastqBatch;

    fn run(
        self: Box<Self>,
        input: Box<dyn InputQueue<Self::Input>>,
        cancel: CancelToken,
    ) -> Result<()> {
        let this = *self;
        let mut writer = BufWriter::with_capacity(
            64 * 1024,
            File::create(&this.path)
                .with_context(|| format!("open FASTQ output {}", this.path.display()))?,
        );

        let mut reorder: ReorderBuffer<SerializedFastqBatch> = ReorderBuffer::new();
        let mut backoff = Backoff::new();

        loop {
            if cancel.is_cancelled() {
                return Ok(());
            }

            if let Some(item) = input.pop() {
                let batch = item.item;
                let ordinal = batch.ordinal;
                reorder.push(ordinal, batch);
                while let Some(ready) = reorder.pop_ready() {
                    writer.write_all(&ready.data).with_context(|| {
                        format!("append FASTQ batch to {}", this.path.display())
                    })?;
                }
                backoff.reset();
            } else if input.is_drained() {
                break;
            } else {
                backoff.snooze();
            }
        }

        // Drain any remaining in-order items (shouldn't be any if queue is
        // drained and ordinals are dense, but pop_ready is safe to call on
        // an empty buffer).
        while let Some(ready) = reorder.pop_ready() {
            writer
                .write_all(&ready.data)
                .with_context(|| format!("drain: append FASTQ batch to {}", this.path.display()))?;
        }

        if !reorder.is_empty() {
            anyhow::bail!(
                "FastqFileSink: {} batches still pending in reorder buffer \
                 (next_expected={}); upstream stage produced non-dense ordinals",
                reorder.len(),
                reorder.next_expected(),
            );
        }

        writer.flush().with_context(|| format!("flush FASTQ output {}", this.path.display()))?;

        Ok(())
    }

    fn name(&self) -> &'static str {
        "FastqFileSink"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::memory::MemoryTracker;
    use crate::runall::engine::queue::StageQueue;
    use crate::runall::engine::stage::SequencedItem;
    use std::sync::Arc;

    fn new_queue() -> Arc<StageQueue<SerializedFastqBatch>> {
        let tracker = Arc::new(MemoryTracker::new(1_000_000_000));
        Arc::new(StageQueue::new("test_fastq_sink_in", 64, 100_000_000, tracker))
    }

    fn make_batch(data: &[u8], ordinal: u64) -> SerializedFastqBatch {
        SerializedFastqBatch { data: data.to_vec(), record_count: 1, ordinal }
    }

    #[test]
    fn test_fastq_sink_out_of_order_writes_in_order() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let q = new_queue();
        q.push(SequencedItem::new(2, make_batch(b"@B2\n", 2), 4)).unwrap();
        q.push(SequencedItem::new(0, make_batch(b"@B0\n", 0), 4)).unwrap();
        q.push(SequencedItem::new(1, make_batch(b"@B1\n", 1), 4)).unwrap();
        q.close();

        let sink = Box::new(FastqFileSink::new(tmp.path().to_path_buf(), 16));
        sink.run(Box::new(q), CancelToken::new()).expect("sink must succeed");

        let bytes = std::fs::read(tmp.path()).unwrap();
        assert_eq!(bytes, b"@B0\n@B1\n@B2\n");
    }

    #[test]
    fn test_fastq_sink_completion_validation_detects_gap() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let q = new_queue();
        // Push 0, 2 — skip 1.
        q.push(SequencedItem::new(0, make_batch(b"@B0\n", 0), 4)).unwrap();
        q.push(SequencedItem::new(2, make_batch(b"@B2\n", 2), 4)).unwrap();
        q.close();

        let sink = Box::new(FastqFileSink::new(tmp.path().to_path_buf(), 16));
        let err = sink
            .run(Box::new(q), CancelToken::new())
            .expect_err("sink must fail on non-dense ordinals");
        let msg = err.to_string();
        assert!(
            msg.contains("pending") && msg.contains("next_expected=1"),
            "unexpected error: {msg}"
        );
    }

    #[test]
    fn test_fastq_sink_empty_input_produces_empty_file() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let q = new_queue();
        q.close();

        let sink = Box::new(FastqFileSink::new(tmp.path().to_path_buf(), 16));
        sink.run(Box::new(q), CancelToken::new()).expect("sink must succeed");

        let bytes = std::fs::read(tmp.path()).unwrap();
        assert!(bytes.is_empty());
    }

    #[test]
    fn test_fastq_sink_cancellation_skips_completion_validation() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let q = new_queue();
        let cancel = CancelToken::new();
        cancel.cancel();

        let sink = Box::new(FastqFileSink::new(tmp.path().to_path_buf(), 16));
        sink.run(Box::new(q), cancel).expect("cancelled sink must return Ok");
    }
}
