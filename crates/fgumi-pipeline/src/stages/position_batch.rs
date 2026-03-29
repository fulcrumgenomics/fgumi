//! `PositionBatcher` — bridges the sort/merge phase into Zone 3 processing.
//!
//! Implements [`SortedRecordSink`] to consume sorted BAM records one at a time,
//! detects position-group boundaries from [`TemplateKey`] changes, and pushes
//! completed [`PositionGroupBatch`] items into a downstream [`WorkQueue`].

use anyhow::Result;
use fgumi_lib::sort::inline_buffer::TemplateKey;
use fgumi_lib::sort::sink::SortedRecordSink;

use super::group_assign::PositionGroupBatch;
use crate::pipeline::queue::WorkQueue;
use crate::stage::SequencedItem;

/// Buffers sorted BAM records by position group and forwards completed batches
/// to a Zone 3 work queue.
///
/// Two records belong to the same position group when their `(primary, secondary,
/// cb_hash)` triple is identical. When any of those fields change, the current
/// batch is flushed as a [`PositionGroupBatch`] and a new batch begins.
///
/// Call [`finish`](PositionBatcher::finish) after the last record to flush any
/// in-progress batch before dropping the batcher.
pub struct PositionBatcher {
    output_queue: WorkQueue<SequencedItem<PositionGroupBatch>>,
    current_batch: Vec<Vec<u8>>,
    /// `(primary, secondary, cb_hash)` of the batch currently being accumulated.
    current_pos_key: Option<(u64, u64, u64)>,
    /// Monotonically increasing sequence number assigned to each flushed batch.
    seq_counter: u64,
}

impl PositionBatcher {
    /// Create a new `PositionBatcher` that pushes completed batches to `output_queue`.
    #[must_use]
    pub fn new(output_queue: WorkQueue<SequencedItem<PositionGroupBatch>>) -> Self {
        Self { output_queue, current_batch: Vec::new(), current_pos_key: None, seq_counter: 0 }
    }

    /// Return the number of batches pushed to the output queue so far.
    ///
    /// Each batch corresponds to one position group; the counter increments on
    /// every successful flush.
    #[must_use]
    pub fn batches_pushed(&self) -> u64 {
        self.seq_counter
    }

    /// Flush the current batch to the output queue if it is non-empty.
    ///
    /// The batch is tagged with the next sequence number and the accumulated
    /// `current_batch` buffer is replaced with an empty `Vec`.
    ///
    /// # Errors
    ///
    /// Returns an error if the output queue channel is disconnected.
    fn flush_batch(&mut self) -> Result<()> {
        if self.current_batch.is_empty() {
            return Ok(());
        }

        let pos_key =
            self.current_pos_key.expect("current_pos_key must be set when batch is non-empty");
        // The batch's `position_key` carries only `(primary, secondary)` — `cb_hash` (pos_key.2)
        // is used solely for boundary detection here (different `cb_hash` values split batches)
        // and is not needed downstream since all records in a batch share the same `cb_hash`.
        let batch = PositionGroupBatch {
            records: std::mem::take(&mut self.current_batch),
            position_key: (pos_key.0, pos_key.1),
        };
        let size_bytes = batch.records.iter().map(Vec::len).sum::<usize>();
        self.output_queue
            .push_blocking(SequencedItem::new(self.seq_counter, batch, size_bytes), size_bytes)?;
        self.seq_counter += 1;
        Ok(())
    }
}

impl SortedRecordSink for PositionBatcher {
    /// Append `record_bytes` to the current position-group batch.
    ///
    /// If `key`'s `(primary, secondary, cb_hash)` differs from the previous record,
    /// the accumulated batch is flushed first and a new batch is started.
    ///
    /// # Errors
    ///
    /// Returns an error if flushing the previous batch fails (channel disconnected).
    fn emit(&mut self, key: &TemplateKey, record_bytes: Vec<u8>) -> Result<()> {
        let pos_key = (key.primary, key.secondary, key.cb_hash);

        if self.current_pos_key != Some(pos_key) {
            self.flush_batch()?;
            self.current_pos_key = Some(pos_key);
        }

        self.current_batch.push(record_bytes);
        Ok(())
    }

    /// Flush any in-progress batch and finalize the sink.
    ///
    /// Must be called once after all records have been emitted to ensure the
    /// last position group is delivered to the output queue.
    ///
    /// # Errors
    ///
    /// Returns an error if the final flush fails (channel disconnected).
    fn finish(&mut self) -> Result<()> {
        self.flush_batch()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_key(primary: u64, secondary: u64, cb_hash: u64) -> TemplateKey {
        TemplateKey { primary, secondary, cb_hash, tertiary: 0, name_hash_upper: 0 }
    }

    #[test]
    fn test_single_position_group() {
        let queue = WorkQueue::new(100, 1_000_000, "test");
        let mut batcher = PositionBatcher::new(queue.clone());

        let key = make_key(1, 2, 0);
        batcher.emit(&key, vec![1, 2, 3]).unwrap();
        batcher.emit(&key, vec![4, 5, 6]).unwrap();
        batcher.finish().unwrap();

        let item = queue.pop().unwrap();
        assert_eq!(item.seq, 0);
        assert_eq!(item.item.records.len(), 2);
        assert_eq!(item.item.records[0], vec![1, 2, 3]);
        assert_eq!(item.item.records[1], vec![4, 5, 6]);
        assert!(queue.try_pop().is_none());
        assert_eq!(batcher.batches_pushed(), 1);
    }

    #[test]
    fn test_multiple_position_groups() {
        let queue = WorkQueue::new(100, 1_000_000, "test");
        let mut batcher = PositionBatcher::new(queue.clone());

        let key1 = make_key(1, 2, 0);
        let key2 = make_key(3, 4, 0);

        batcher.emit(&key1, vec![1]).unwrap();
        batcher.emit(&key1, vec![2]).unwrap();
        batcher.emit(&key2, vec![3]).unwrap();
        batcher.finish().unwrap();

        let batch1 = queue.pop().unwrap();
        assert_eq!(batch1.seq, 0);
        assert_eq!(batch1.item.records.len(), 2);
        assert_eq!(batch1.item.position_key, (1, 2));

        let batch2 = queue.pop().unwrap();
        assert_eq!(batch2.seq, 1);
        assert_eq!(batch2.item.records.len(), 1);
        assert_eq!(batch2.item.position_key, (3, 4));

        assert!(queue.try_pop().is_none());
        assert_eq!(batcher.batches_pushed(), 2);
    }

    #[test]
    fn test_position_key_uses_cb_hash() {
        // Two records that differ only in cb_hash must land in separate batches.
        let queue = WorkQueue::new(100, 1_000_000, "test");
        let mut batcher = PositionBatcher::new(queue.clone());

        let key_a = make_key(1, 1, 0xaaa);
        let key_b = make_key(1, 1, 0xbbb);

        batcher.emit(&key_a, vec![10]).unwrap();
        batcher.emit(&key_b, vec![20]).unwrap();
        batcher.finish().unwrap();

        let batch_a = queue.pop().unwrap();
        assert_eq!(batch_a.item.records.len(), 1);
        assert_eq!(batch_a.item.records[0], vec![10]);

        let batch_b = queue.pop().unwrap();
        assert_eq!(batch_b.item.records.len(), 1);
        assert_eq!(batch_b.item.records[0], vec![20]);
    }

    #[test]
    fn test_tertiary_ignored_for_boundary() {
        // Records with the same (primary, secondary, cb_hash) but different tertiary
        // (i.e. different MI) must be batched together — MI grouping happens downstream.
        let queue = WorkQueue::new(100, 1_000_000, "test");
        let mut batcher = PositionBatcher::new(queue.clone());

        let key1 =
            TemplateKey { primary: 5, secondary: 6, cb_hash: 0, tertiary: 1, name_hash_upper: 0 };
        let key2 =
            TemplateKey { primary: 5, secondary: 6, cb_hash: 0, tertiary: 2, name_hash_upper: 100 };

        batcher.emit(&key1, vec![1]).unwrap();
        batcher.emit(&key2, vec![2]).unwrap();
        batcher.finish().unwrap();

        let batch = queue.pop().unwrap();
        assert_eq!(
            batch.item.records.len(),
            2,
            "different tertiary/name_hash should not split the batch"
        );
        assert!(queue.try_pop().is_none());
    }

    #[test]
    fn test_empty_no_batches_pushed() {
        let queue = WorkQueue::new(100, 1_000_000, "test");
        let mut batcher = PositionBatcher::new(queue.clone());
        batcher.finish().unwrap();
        assert!(queue.try_pop().is_none());
        assert_eq!(batcher.batches_pushed(), 0);
    }

    #[test]
    fn test_sequence_numbers_are_monotonic() {
        let queue = WorkQueue::new(100, 1_000_000, "test");
        let mut batcher = PositionBatcher::new(queue.clone());

        for i in 0u64..5 {
            let key = make_key(i, 0, 0);
            #[allow(clippy::cast_possible_truncation)]
            batcher.emit(&key, vec![i as u8]).unwrap();
        }
        batcher.finish().unwrap();

        for expected_seq in 0u64..5 {
            let item = queue.pop().unwrap();
            assert_eq!(item.seq, expected_seq);
        }
        assert!(queue.try_pop().is_none());
    }
}
