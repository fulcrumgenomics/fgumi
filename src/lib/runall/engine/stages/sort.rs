//! Sort stage: template-coordinate sort with external merge.
//!
//! [`SpecialStage`] barrier: consumes **all** input before producing
//! any output. Wraps `RawExternalSorter` from `src/lib/sort/raw.rs`:
//! the accumulation phase may spill chunks to disk under the supplied
//! memory limit; the merge phase performs a k-way merge and pushes
//! records in sort order, re-batched into [`SerializedBatch`] chunks.
//!
//! Use via `Pipeline::builder().special_stage(SortStage::new(...))`.
//!
//! ## Input
//!
//! [`SerializedBatch`] — length-prefixed BAM record bytes in any
//! order.
//!
//! ## Output
//!
//! [`SerializedBatch`] — up to 256 records per batch
//! (`MAX_RECORDS_PER_BATCH`) in template-coordinate sort order, with
//! a fresh monotonic `ordinal` counter (0..N) stamped at merge time.
//! Each record is re-framed with a 4-byte little-endian `block_size`
//! prefix matching the format emitted by `ExtractStage`,
//! `CorrectStage`, and `BamBatchedSource`.
//!
//! ## Ordering guarantees
//!
//! Strict template-coordinate order (see
//! [`crate::sort::SortOrder::TemplateCoordinate`]). Input ordering is
//! discarded; output ordering is canonical.
//!
//! ## Memory model
//!
//! In-memory accumulation up to `memory_limit` bytes, then spill to
//! disk under `temp_dir`. Phase 2 merge streams records one at a
//! time through an internal `BatchingQueueSink`, batching up to 256
//! records per output chunk. A bounded `sync_channel(1024)` between
//! the drainer thread and the sorter caps in-flight per-record
//! buffers.
//!
//! ## Determinism
//!
//! Sort order is deterministic (stable under equal keys). The
//! accumulate-then-merge split across threads does not affect output:
//! final ordering is imposed by k-way merge on the sorted chunks.

use std::path::PathBuf;

use anyhow::Context;
use noodles::sam::Header;

use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::grouping_types::iter_length_prefixed;
use crate::runall::engine::output_types::{RawBytes, SerializedBatch};
use crate::runall::engine::sink::InputQueue;
use crate::runall::engine::source::OutputQueue;
use crate::runall::engine::special_stage::SpecialStage;
use crate::runall::engine::stage::SequencedItem;
use crate::sort::sink::SortedRecordSink;
use crate::sort::{RawExternalSorter, SortOrder, TemplateKey, merge_sorted_chunks_to_sink};

/// Target records per emitted `SerializedBatch`. Matches the cap used by
/// `BamBatchedSource`, `FastqBlockMerge`, and `FastqPair` so downstream stages
/// see consistent batch sizes regardless of entry point.
const MAX_RECORDS_PER_BATCH: u64 = 256;

/// Template-coordinate sort stage (barrier).
pub struct SortStage {
    /// Memory limit for in-memory accumulation before spilling.
    pub memory_limit: usize,
    /// Temp directory for spill files (None uses the system default).
    pub temp_dir: Option<PathBuf>,
    /// Threads for parallel sorting.
    pub threads: usize,
    /// BAM header (used for library lookups and sort-key extraction).
    pub header: Header,
}

impl SortStage {
    /// Construct a new `SortStage`.
    #[must_use]
    pub fn new(
        memory_limit: usize,
        temp_dir: Option<PathBuf>,
        threads: usize,
        header: Header,
    ) -> Self {
        Self { memory_limit, temp_dir, threads, header }
    }
}

/// `SortedRecordSink` that accumulates sorted records into
/// [`SerializedBatch`] chunks of up to [`MAX_RECORDS_PER_BATCH`] records and
/// pushes them to an `OutputQueue<SerializedBatch>` with monotonic ordinals.
///
/// Each emitted `record_bytes` is the BAM record body (no `block_size`
/// prefix); this sink prepends the 4-byte little-endian `block_size` before
/// appending to the per-batch buffer so downstream consumers see the same
/// length-prefixed format emitted by `ExtractStage`, `CorrectStage`, and
/// `BamBatchedSource`.
struct BatchingQueueSink<'a> {
    output: &'a dyn OutputQueue<SerializedBatch>,
    cancel: &'a CancelToken,
    data: Vec<u8>,
    record_count: u64,
    ordinal: u64,
}

impl<'a> BatchingQueueSink<'a> {
    fn new(output: &'a dyn OutputQueue<SerializedBatch>, cancel: &'a CancelToken) -> Self {
        Self { output, cancel, data: Vec::with_capacity(256 * 1024), record_count: 0, ordinal: 0 }
    }

    // `flush` keeps the `Result` return type to match the
    // `SortedRecordSink::finish` contract (which calls `flush` directly) and
    // to leave room for future flush-time errors; today the only failure
    // path was a "cancelled during merge" bail that was shadowing the real
    // upstream error, so it now exits cleanly.
    #[allow(clippy::unnecessary_wraps)]
    fn flush(&mut self) -> anyhow::Result<()> {
        if self.record_count == 0 {
            return Ok(());
        }
        let batch = SerializedBatch {
            primary: RawBytes {
                data: std::mem::take(&mut self.data),
                record_count: self.record_count,
            },
            secondary: None,
            ordinal: self.ordinal,
        };
        let mem = batch.primary.data.len();
        let item = SequencedItem::new(self.ordinal, batch, mem);
        // Cancellation here means another stage already errored; exit
        // cleanly so the real error wins in the driver's first-error-in-
        // join-order selection rather than shadowing it with our own
        // generic "sort cancelled" message.
        if self.output.push_until_cancelled(item, self.cancel).is_err() {
            return Ok(());
        }
        self.ordinal += 1;
        self.record_count = 0;
        self.data = Vec::with_capacity(256 * 1024);
        Ok(())
    }
}

impl SortedRecordSink for BatchingQueueSink<'_> {
    fn emit(&mut self, _key: &TemplateKey, record_bytes: Vec<u8>) -> anyhow::Result<()> {
        let block_size = u32::try_from(record_bytes.len())
            .map_err(|_| anyhow::anyhow!("SortStage: record exceeds u32::MAX bytes"))?;
        self.data.extend_from_slice(&block_size.to_le_bytes());
        self.data.extend_from_slice(&record_bytes);
        self.record_count += 1;
        if self.record_count >= MAX_RECORDS_PER_BATCH {
            self.flush()?;
        }
        Ok(())
    }

    fn finish(&mut self) -> anyhow::Result<()> {
        self.flush()
    }
}

impl SpecialStage for SortStage {
    type Input = SerializedBatch;
    type Output = SerializedBatch;

    #[tracing::instrument(name = "sort", skip_all)]
    fn run(
        self: Box<Self>,
        input: Box<dyn InputQueue<SerializedBatch>>,
        output: Box<dyn OutputQueue<SerializedBatch>>,
        cancel: CancelToken,
    ) -> anyhow::Result<()> {
        // Phase 1: accumulate.
        let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<u8>>(1024);

        let cancel_for_drainer = cancel.clone();
        let drainer = std::thread::Builder::new().name("sort_drainer".into()).spawn(
            move || -> anyhow::Result<u64> {
                let mut backoff = crate::runall::engine::backoff::Backoff::new();
                let mut count: u64 = 0;
                loop {
                    if cancel_for_drainer.is_cancelled() {
                        break;
                    }
                    if let Some(item) = input.pop() {
                        for record in iter_length_prefixed(&item.item.primary.data) {
                            let record = record.context("SortStage: framed-BAM parse error")?;
                            if tx.send(record.to_vec()).is_err() {
                                return Ok(count);
                            }
                            count += 1;
                        }
                        backoff.reset();
                    } else if input.is_drained() {
                        break;
                    } else {
                        backoff.snooze();
                    }
                }
                drop(tx);
                Ok(count)
            },
        )?;

        let sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
            .memory_limit(self.memory_limit)
            .threads(self.threads);
        let sorter =
            if let Some(path) = self.temp_dir.clone() { sorter.temp_dir(path) } else { sorter };

        let sort_result = sorter.sort_accumulate_from_iter(rx.into_iter(), &self.header);
        if sort_result.is_err() {
            cancel.cancel();
        }
        let drainer_result =
            drainer.join().map_err(|_| anyhow::anyhow!("SortStage: drainer thread panicked"))?;
        let sorted_chunks = match (sort_result, drainer_result) {
            (Ok(chunks), Ok(_count)) => chunks,
            (Err(sort_err), _) => {
                output.close();
                return Err(sort_err);
            }
            (Ok(_), Err(drain_err)) => {
                output.close();
                return Err(drain_err);
            }
        };

        if cancel.is_cancelled() {
            output.close();
            anyhow::bail!("sort cancelled");
        }

        // Phase 2: merge — stream sorted records into re-batched output.
        let mut sink = BatchingQueueSink::new(&*output, &cancel);
        merge_sorted_chunks_to_sink(sorted_chunks, &mut sink)?;
        sink.finish()?;

        output.close();
        Ok(())
    }

    fn name(&self) -> &'static str {
        "SortStage"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::memory::MemoryTracker;

    fn empty_header() -> Header {
        use noodles::sam::header::record::value::Map;
        Header::builder()
            .set_header(Map::<noodles::sam::header::record::value::map::Header>::new(
                noodles::sam::header::record::value::map::header::Version::new(1, 6),
            ))
            .build()
    }

    #[test]
    fn test_sort_stage_special_empty_input() {
        use std::sync::Arc;

        use crate::runall::engine::queue::StageQueue;

        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let input: Arc<StageQueue<SerializedBatch>> =
            Arc::new(StageQueue::new("sort_spec_in", 16, 100_000, tracker.clone()));
        let output: Arc<StageQueue<SerializedBatch>> =
            Arc::new(StageQueue::new("sort_spec_out", 16, 100_000, tracker));

        input.close();

        let stage: Box<SortStage> =
            Box::new(SortStage::new(64 * 1024 * 1024, None, 1, empty_header()));
        <SortStage as SpecialStage>::run(
            stage,
            Box::new(input),
            Box::new(output.clone()),
            CancelToken::new(),
        )
        .unwrap();

        assert!(OutputQueue::is_closed(&output));
        assert_eq!(output.pop().map(|_| ()), None);
    }
}
