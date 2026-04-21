//! Position-batch stage: groups sorted BAM records by template-
//! coordinate position key.
//!
//! [`SpecialStage`] — single-threaded because boundary detection needs
//! in-order per-record scanning. The generic engine drives it
//! alongside pool-parallel stages.
//!
//! ## Input
//!
//! [`SerializedBatch`] — concat-byte BAM records in template-coordinate
//! sort order, as emitted by
//! [`crate::runall::engine::stages::sort::SortStage`].
//!
//! ## Output
//!
//! [`PositionGroupBatch`] — one batch per distinct
//! `(primary, secondary, cb_hash)` sort key encountered; `data` holds
//! all records for that key in arrival order, with a fresh monotonic
//! `ordinal` counter (0..N) stamped on emit. `position_key` is
//! `(primary, secondary)` (cell-barcode hash is consumed internally).
//!
//! ## Ordering guarantees
//!
//! Strictly in-order by upstream sort key: the loop tracks the current
//! key and flushes whenever it changes, guaranteeing one output batch
//! per key and that batches are emitted in sort-key order. Requires
//! that the upstream `SortStage` emits records in canonical sort
//! order (it does).
//!
//! ## Memory model
//!
//! Bounded by the size of the largest position bucket. The in-flight
//! buffer is preallocated to 256 KB and reset after each flush.
//!
//! ## Determinism
//!
//! Fully deterministic given a deterministic (sorted) input stream.
//! Cell-barcode hashing (`cb_hasher()`) is a keyed but deterministic
//! hash.

use anyhow::Result;
use noodles::sam::Header;

use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::grouping_types::{PositionGroupBatch, iter_length_prefixed};
use crate::runall::engine::output_types::SerializedBatch;
use crate::runall::engine::sink::InputQueue;
use crate::runall::engine::source::OutputQueue;
use crate::runall::engine::special_stage::SpecialStage;
use crate::runall::engine::stage::SequencedItem;
use crate::sam::SamTag;
use crate::sort::{LibraryLookup, cb_hasher, extract_template_key_inline};

/// Boundary-detecting position-batch stage.
pub struct PositionBatchStage {
    header: Header,
    cell_tag: Option<SamTag>,
}

impl PositionBatchStage {
    /// Construct a new `PositionBatchStage` with no cell-barcode hashing.
    #[must_use]
    pub fn new(header: Header) -> Self {
        Self { header, cell_tag: None }
    }

    /// Construct with cell-barcode hashing for single-cell workflows.
    #[must_use]
    pub fn with_cell_tag(header: Header, cell_tag: SamTag) -> Self {
        Self { header, cell_tag: Some(cell_tag) }
    }
}

impl SpecialStage for PositionBatchStage {
    type Input = SerializedBatch;
    type Output = PositionGroupBatch;

    #[tracing::instrument(name = "position_batch", skip_all)]
    fn run(
        self: Box<Self>,
        input: Box<dyn InputQueue<SerializedBatch>>,
        output: Box<dyn OutputQueue<PositionGroupBatch>>,
        cancel: CancelToken,
    ) -> Result<()> {
        let lib_lookup = LibraryLookup::from_header(&self.header);
        let hasher = cb_hasher();

        let mut current_key: Option<(u64, u64, u64)> = None;
        let mut current_data: Vec<u8> = Vec::with_capacity(256 * 1024);
        let mut current_count: u64 = 0;
        let mut out_ordinal: u64 = 0;

        let mut backoff = crate::runall::engine::backoff::Backoff::new();

        let flush = |current_key: &mut Option<(u64, u64, u64)>,
                     current_data: &mut Vec<u8>,
                     current_count: &mut u64,
                     out_ordinal: &mut u64|
         -> Result<()> {
            if *current_count == 0 {
                return Ok(());
            }
            let Some((primary, secondary, _cb_hash)) = *current_key else {
                return Ok(());
            };
            let data = std::mem::take(current_data);
            let batch = PositionGroupBatch {
                data,
                record_count: *current_count,
                ordinal: *out_ordinal,
                position_key: (primary, secondary),
            };
            let mem = batch.data.len();
            let item = SequencedItem::new(*out_ordinal, batch, mem);
            if output.push_until_cancelled(item, &cancel).is_err() {
                anyhow::bail!("PositionBatchStage: cancelled during push");
            }
            *out_ordinal += 1;
            *current_count = 0;
            *current_data = Vec::with_capacity(256 * 1024);
            *current_key = None;
            Ok(())
        };

        loop {
            if cancel.is_cancelled() {
                break;
            }
            if let Some(item) = input.pop() {
                // Iterate length-prefixed records within the input batch.
                // Upstream `SortStage` emits `SerializedBatch` items with
                // monotonic ordinals in sort-key order, so records within
                // and across batches arrive in sort-key order already.
                for rec_result in iter_length_prefixed(&item.item.primary.data) {
                    let rec = rec_result?;
                    let tkey =
                        extract_template_key_inline(rec, &lib_lookup, self.cell_tag, &hasher);
                    let pkey = (tkey.primary, tkey.secondary, tkey.cb_hash);

                    match current_key {
                        Some(ck) if ck == pkey => {
                            let bs = u32::try_from(rec.len()).map_err(|_| {
                                anyhow::anyhow!("PositionBatchStage: record exceeds u32::MAX bytes")
                            })?;
                            current_data.extend_from_slice(&bs.to_le_bytes());
                            current_data.extend_from_slice(rec);
                            current_count += 1;
                        }
                        _ => {
                            flush(
                                &mut current_key,
                                &mut current_data,
                                &mut current_count,
                                &mut out_ordinal,
                            )?;
                            current_key = Some(pkey);
                            let bs = u32::try_from(rec.len()).map_err(|_| {
                                anyhow::anyhow!("PositionBatchStage: record exceeds u32::MAX bytes")
                            })?;
                            current_data.extend_from_slice(&bs.to_le_bytes());
                            current_data.extend_from_slice(rec);
                            current_count += 1;
                        }
                    }
                }
                backoff.reset();
            } else if input.is_drained() {
                break;
            } else {
                backoff.snooze();
            }
        }

        // Flush the trailing batch.
        flush(&mut current_key, &mut current_data, &mut current_count, &mut out_ordinal)?;
        output.close();
        Ok(())
    }

    fn name(&self) -> &'static str {
        "PositionBatchStage"
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::runall::engine::memory::MemoryTracker;
    use crate::runall::engine::queue::StageQueue;

    fn empty_header() -> Header {
        use noodles::sam::header::record::value::Map;
        Header::builder()
            .set_header(Map::<noodles::sam::header::record::value::map::Header>::new(
                noodles::sam::header::record::value::map::header::Version::new(1, 6),
            ))
            .build()
    }

    #[test]
    fn test_position_batch_empty_input() {
        let tracker = Arc::new(MemoryTracker::new(1_000_000));
        let input: Arc<StageQueue<SerializedBatch>> =
            Arc::new(StageQueue::new("pb_in", 16, 100_000, tracker.clone()));
        let output: Arc<StageQueue<PositionGroupBatch>> =
            Arc::new(StageQueue::new("pb_out", 16, 100_000, tracker));

        input.close();

        let stage = Box::new(PositionBatchStage::new(empty_header()));
        <PositionBatchStage as SpecialStage>::run(
            stage,
            Box::new(input),
            Box::new(output.clone()),
            CancelToken::new(),
        )
        .unwrap();

        assert!(OutputQueue::is_closed(&output));
        assert!(output.pop().is_none());
    }
}
