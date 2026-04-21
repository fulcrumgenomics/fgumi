//! BAM file Sink: consumes `ConsensusOutput` items (pre-serialized BAM records
//! with 4-byte `block_size` prefixes), reorders them by sequence number, and
//! writes them through a `RawBamWriter` which handles BGZF compression.

use std::path::PathBuf;
use std::sync::{Arc, OnceLock};

use anyhow::{Context, Result};
use fgumi_consensus::caller::ConsensusOutput;
use noodles::sam::Header;

use crate::bam_io::{RawBamWriter, create_raw_bam_writer};
use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::reorder::ReorderBuffer;
use crate::runall::engine::sink::{InputQueue, Sink};

/// Header handle shared between `AlignAndMerge` and [`BamSink`]: the sink
/// waits for the aligner-aware output header to arrive before opening the
/// writer, matching what standalone `fgumi zipper` writes.
///
/// When `None`, the sink uses its plan-time `header`. When `Some`, the sink
/// blocks until `AlignAndMerge::run` populates the `OnceLock` with the header
/// built from `(unmapped_header, mapped_header, dict_path)`.
pub type DeferredHeader = Arc<OnceLock<Header>>;

/// BAM file Sink.
///
/// Consumes `ConsensusOutput` items in arbitrary order, reorders them by
/// sequence number via a [`ReorderBuffer`], and writes each item's BAM record
/// bytes via [`RawBamWriter`] (which applies BGZF compression).
pub struct BamSink {
    output_path: PathBuf,
    header: Header,
    threads: usize,
    compression_level: u32,
    /// Optional runtime-populated header override. When present, the sink
    /// waits for this `OnceLock` to be set before opening the writer. Used
    /// by plans that include `AlignAndMerge`, whose output header requires
    /// the aligner's mapped header (only available at runtime).
    deferred_header: Option<DeferredHeader>,
}

impl BamSink {
    /// Construct a new `BamSink` writing to `output_path`.
    #[must_use]
    pub fn new(
        output_path: PathBuf,
        header: Header,
        threads: usize,
        compression_level: u32,
    ) -> Self {
        Self { output_path, header, threads, compression_level, deferred_header: None }
    }

    /// Attach a deferred header that will be populated at runtime (by
    /// [`AlignAndMerge`]). The sink will wait for it before opening the
    /// output writer, falling back to the plan-time header only if the
    /// deferred slot stays empty.
    #[must_use]
    pub fn with_deferred_header(mut self, deferred_header: DeferredHeader) -> Self {
        self.deferred_header = Some(deferred_header);
        self
    }
}

impl Sink for BamSink {
    type Input = ConsensusOutput;

    fn run(
        self: Box<Self>,
        input: Box<dyn InputQueue<Self::Input>>,
        cancel: CancelToken,
    ) -> Result<()> {
        let this = *self;
        // Resolve the effective header. If a deferred slot is attached, wait
        // for AlignAndMerge to populate it; fall back to the plan-time header
        // if the producer never sets it (e.g. cancelled before first batch).
        let effective_header = if let Some(ref slot) = this.deferred_header {
            wait_for_deferred_header(slot, &cancel, &this.header)
        } else {
            this.header.clone()
        };
        let mut writer: RawBamWriter = create_raw_bam_writer(
            &this.output_path,
            &effective_header,
            this.threads,
            this.compression_level,
        )?;

        let mut reorder: ReorderBuffer<ConsensusOutput> = ReorderBuffer::new();
        let mut backoff = crate::runall::engine::backoff::Backoff::new();

        loop {
            if cancel.is_cancelled() {
                return Ok(());
            }

            if let Some(item) = input.pop() {
                reorder.push(item.seq, item.item);
                while let Some(co) = reorder.pop_ready() {
                    let count = co.count as u64;
                    writer.write_raw_bytes(&co.data).context("failed to write consensus batch")?;
                    if count > 0 {
                        crate::progress::records_written(count);
                    }
                }
                backoff.reset();
            } else if input.is_drained() {
                break;
            } else {
                backoff.snooze();
            }
        }

        // Drain any remaining reorder-buffered items (may exist if sequence
        // numbers had gaps that were filled as the last items arrived).
        while let Some(co) = reorder.pop_ready() {
            let count = co.count as u64;
            writer.write_raw_bytes(&co.data).context("failed to write consensus batch")?;
            if count > 0 {
                crate::progress::records_written(count);
            }
        }

        writer.finish().context("failed to finalize BAM output")?;
        Ok(())
    }

    fn name(&self) -> &'static str {
        "BamSink"
    }
}

/// Block until `slot` has been populated by the producer (e.g. `AlignAndMerge`),
/// or until cancellation. Falls back to `fallback` if the wait is cancelled
/// before the slot is filled. A stuck producer is the pipeline watchdog's
/// problem, not this function's.
fn wait_for_deferred_header(
    slot: &OnceLock<Header>,
    cancel: &CancelToken,
    fallback: &Header,
) -> Header {
    let mut backoff = crate::runall::engine::backoff::Backoff::new();
    loop {
        if let Some(h) = slot.get() {
            return h.clone();
        }
        if cancel.is_cancelled() {
            // Open the writer with the plan-time header so an empty cancelled
            // output is still a valid BAM.
            return fallback.clone();
        }
        backoff.snooze();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::memory::MemoryTracker;
    use crate::runall::engine::queue::StageQueue;
    use crate::runall::engine::stage::SequencedItem;
    use std::sync::Arc;

    fn empty_header() -> Header {
        use noodles::sam::header::record::value::Map;
        Header::builder()
            .set_header(Map::<noodles::sam::header::record::value::map::Header>::new(
                noodles::sam::header::record::value::map::header::Version::new(1, 6),
            ))
            .build()
    }

    #[test]
    fn test_bamsink_reorders_and_writes_empty_header() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let header = empty_header();

        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let queue: Arc<StageQueue<ConsensusOutput>> =
            Arc::new(StageQueue::new("sink_in", 16, 100_000, tracker));

        // Push out of order.
        let co_a = ConsensusOutput { data: vec![], count: 0 };
        let co_b = ConsensusOutput { data: vec![], count: 0 };
        // `ConsensusOutput` does not implement `Debug`, so use `ok()` + `expect`.
        assert!(queue.push(SequencedItem::new(1, co_b, 1)).is_ok());
        assert!(queue.push(SequencedItem::new(0, co_a, 1)).is_ok());
        queue.close();

        let sink = Box::new(BamSink::new(tmp.path().to_path_buf(), header, 1, 6));
        sink.run(Box::new(queue), CancelToken::new()).unwrap();

        assert!(tmp.path().exists());
    }

    #[test]
    fn test_bamsink_name() {
        let sink = BamSink::new(PathBuf::from("/tmp/out.bam"), empty_header(), 1, 6);
        assert_eq!(Sink::name(&sink), "BamSink");
    }
}
