//! Batch flow control for aligner subprocess coordination.
//!
//! When bwa mem is run with `-K <bases>`, it processes reads in chunks of approximately
//! K bases. `ToFastq` and `ReadSam` coordinate through [`AlignerBatchState`] to prevent pipe
//! deadlocks: `ToFastq` sends one batch, waits for `ReadSam` to consume the SAM output,
//! then sends the next batch.

use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

/// Shared state for coordinating FASTQ send / SAM receive batch boundaries.
///
/// `ToFastq` writes FASTQ records to the aligner's stdin, tracking bases and record
/// counts. When a batch boundary is reached (enough bases with an even record count
/// to keep read pairs together), it marks the batch as sent and waits. `ReadSam` reads
/// SAM records from the aligner's stdout, and once all records for the batch have been
/// received, it calls [`advance_batch`](Self::advance_batch) to reset counters and
/// unblock the next batch.
pub struct AlignerBatchState {
    /// bwa's -K value in bases.
    k_chunk_bases: u64,
    /// Bases written to stdin in the current batch.
    bases_sent: AtomicU64,
    /// Records sent in the current batch.
    records_sent: AtomicU64,
    /// True when the current batch has been fully sent to stdin.
    batch_sent: AtomicBool,
    /// Records received from stdout for the current batch.
    records_received: AtomicU64,
    /// Batch sequence counter for debugging/logging.
    batch_number: AtomicU64,
}

impl AlignerBatchState {
    /// Create a new batch state with the given `-K` value.
    #[must_use]
    pub fn new(k_chunk_bases: u64) -> Self {
        Self {
            k_chunk_bases,
            bases_sent: AtomicU64::new(0),
            records_sent: AtomicU64::new(0),
            batch_sent: AtomicBool::new(false),
            records_received: AtomicU64::new(0),
            batch_number: AtomicU64::new(0),
        }
    }

    /// Called by `ToFastq`: can we send more FASTQ data?
    ///
    /// Returns `true` if the current batch has not yet been fully sent, meaning
    /// `ToFastq` can continue writing records to the aligner's stdin.
    #[must_use]
    pub fn can_send(&self) -> bool {
        !self.batch_sent.load(Ordering::Acquire)
    }

    /// Called by `ToFastq`: record that `bases` bases and `records` records were written
    /// to stdin.
    ///
    /// When the accumulated bases reach `k_chunk_bases` AND the record count is even
    /// (to keep read pairs together, matching bwa's behavior), marks the batch as sent.
    pub fn record_send(&self, bases: u64, records: u64) {
        let total_bases = self.bases_sent.fetch_add(bases, Ordering::AcqRel) + bases;
        let total_records = self.records_sent.fetch_add(records, Ordering::AcqRel) + records;
        // Batch boundary: enough bases AND even record count (read pairs stay together)
        if total_bases >= self.k_chunk_bases && total_records.is_multiple_of(2) {
            self.batch_sent.store(true, Ordering::Release);
        }
    }

    /// Called by `ReadSam`: record that `records` records were received from stdout.
    pub fn record_receive(&self, records: u64) {
        self.records_received.fetch_add(records, Ordering::AcqRel);
    }

    /// Called by `ReadSam`: check if the current batch is complete.
    ///
    /// A batch is complete when:
    /// 1. `batch_sent` is true (`ToFastq` finished sending)
    /// 2. `records_received >= records_sent` (all SAM output consumed)
    #[must_use]
    pub fn is_batch_complete(&self) -> bool {
        if !self.batch_sent.load(Ordering::Acquire) {
            return false;
        }
        let sent = self.records_sent.load(Ordering::Acquire);
        let received = self.records_received.load(Ordering::Acquire);
        received >= sent
    }

    /// Called by `ReadSam`: signal that the current batch's SAM output is fully consumed.
    /// Resets counters for the next batch.
    pub fn advance_batch(&self) {
        self.bases_sent.store(0, Ordering::Release);
        self.records_sent.store(0, Ordering::Release);
        self.records_received.store(0, Ordering::Release);
        self.batch_sent.store(false, Ordering::Release);
        self.batch_number.fetch_add(1, Ordering::Relaxed);
    }

    /// Current batch number (for logging/debugging).
    #[must_use]
    pub fn batch_number(&self) -> u64 {
        self.batch_number.load(Ordering::Relaxed)
    }

    /// The `-K` value in bases.
    #[must_use]
    pub fn k_chunk_bases(&self) -> u64 {
        self.k_chunk_bases
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_can_send_initially() {
        let state = AlignerBatchState::new(1000);
        assert!(state.can_send(), "should be able to send initially");
        assert_eq!(state.batch_number(), 0);
    }

    #[test]
    fn test_batch_boundary() {
        let state = AlignerBatchState::new(100);

        // Send 100 bases in 2 records (even) — should trigger batch boundary.
        state.record_send(100, 2);

        assert!(!state.can_send(), "batch_sent should be true after reaching k_chunk_bases");
        assert!(
            !state.is_batch_complete(),
            "batch should not be complete until records are received"
        );
    }

    #[test]
    fn test_odd_records_delay_boundary() {
        let state = AlignerBatchState::new(100);

        // Send 150 bases but only 1 record (odd) — should NOT trigger batch boundary.
        state.record_send(150, 1);
        assert!(state.can_send(), "odd record count should delay batch boundary");

        // Send 1 more record with 0 bases — now 2 records (even), boundary should trigger.
        state.record_send(0, 1);
        assert!(!state.can_send(), "even record count should trigger batch boundary");
    }

    #[test]
    fn test_batch_complete_detection() {
        let state = AlignerBatchState::new(100);

        // Not complete when nothing has been sent.
        assert!(!state.is_batch_complete());

        // Send a full batch: 100 bases, 4 records.
        state.record_send(100, 4);
        assert!(!state.is_batch_complete(), "not complete until records received");

        // Receive 3 of 4 records.
        state.record_receive(3);
        assert!(!state.is_batch_complete(), "not complete with fewer records than sent");

        // Receive the 4th record.
        state.record_receive(1);
        assert!(state.is_batch_complete(), "complete when records_received >= records_sent");
    }

    #[test]
    fn test_batch_complete_with_extra_receives() {
        let state = AlignerBatchState::new(50);

        state.record_send(50, 2);
        // Receive more records than sent (e.g. supplementary alignments).
        state.record_receive(4);
        assert!(state.is_batch_complete(), "complete when records_received > records_sent");
    }

    #[test]
    fn test_advance_batch_resets() {
        let state = AlignerBatchState::new(100);

        // Complete a batch.
        state.record_send(100, 2);
        state.record_receive(2);
        assert!(state.is_batch_complete());
        assert_eq!(state.batch_number(), 0);

        // Advance to next batch.
        state.advance_batch();
        assert!(state.can_send(), "should be able to send after advance");
        assert!(!state.is_batch_complete(), "new batch should not be complete");
        assert_eq!(state.batch_number(), 1);
    }

    #[test]
    fn test_batch_number_increments() {
        let state = AlignerBatchState::new(10);

        for i in 0..5 {
            assert_eq!(state.batch_number(), i);
            state.record_send(10, 2);
            state.record_receive(2);
            state.advance_batch();
        }
        assert_eq!(state.batch_number(), 5);
    }

    #[test]
    fn test_k_chunk_bases_accessor() {
        let state = AlignerBatchState::new(150_000_000);
        assert_eq!(state.k_chunk_bases(), 150_000_000);
    }

    #[test]
    fn test_incremental_sends_accumulate() {
        let state = AlignerBatchState::new(100);

        // Send in small increments — should not trigger boundary until total >= k.
        state.record_send(30, 2);
        assert!(state.can_send());
        state.record_send(30, 2);
        assert!(state.can_send());
        // Now at 80 bases, 6 records (even) — still under k.
        state.record_send(20, 2);
        assert!(state.can_send()); // 80 < 100, no trigger yet
        // Now send enough to exceed k: 80 + 30 = 110 bases, 8 records (even).
        state.record_send(30, 2);
        assert!(!state.can_send());
    }
}
