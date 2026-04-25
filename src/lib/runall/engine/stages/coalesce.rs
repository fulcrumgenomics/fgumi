//! `Coalesce` [`SpecialStage`]: reorder incoming batches by ordinal and
//! coalesce their bytes into larger chunks sized to match ~N full BGZF
//! blocks.
//!
//! Generic over [`CoalesceInput`] so the same reorder-and-concatenate
//! engine handles both byte-level batches ([`SerializedBatch`] — emitted
//! by Extract/Correct/AlignAndMerge/Sort) and structured batches
//! ([`MiGroupBatch`](crate::runall::engine::grouping_types::MiGroupBatch) —
//! emitted by `GroupAssignStage`). Each input type flattens to bytes via
//! `CoalesceInput::into_parts`.
//!
//! ## Input
//!
//! Any `T: CoalesceInput` — currently [`SerializedBatch`] or
//! [`MiGroupBatch`](crate::runall::engine::grouping_types::MiGroupBatch).
//! Input batches carry an upstream ordinal and either concat-byte or
//! structured payloads.
//!
//! ## Output
//!
//! [`SerializedBatch`] — concatenated primary (and optional secondary)
//! bytes with a fresh, monotonic 0..N `ordinal` counter stamped by this
//! stage. Record counts are summed across the contributing batches.
//!
//! ## Ordering guarantees
//!
//! Strictly in-order by input `ordinal`: a [`ReorderBuffer`] consumes
//! input ordinals from 0..N and emits concatenations in that exact
//! order. Structured inputs (e.g. MI groups) are flattened in the
//! order their container presents them — `GroupAssignStage` sorts MI
//! groups ascending by MI before emitting, so the final byte stream is
//! canonical.
//!
//! ## Memory model
//!
//! Bounded by the `ReorderBuffer` (holds at most "gap-size" batches
//! waiting for the next ordinal) plus one in-flight output chunk. Each
//! emitted chunk is a fresh [`Vec<u8>`] of capacity `target_chunk_size`
//! (default 4 × `BGZF_MAX_BLOCK_SIZE`); the previous chunk's storage
//! is handed to the output queue via `std::mem::take`.
//!
//! ## Determinism
//!
//! Fully deterministic given a deterministic upstream input: reorder
//! imposes a canonical ordinal sequence, and concatenation is
//! commutative-free under that order.
//!
//! ## Rationale for coalescing
//!
//! Upstream pool-parallel stages emit at the cadence of their internal
//! batching (e.g., 256 records ≈ 87 KB for 150 bp reads). Compressed
//! individually by `BgzfCompress`, each batch yields one full 64 KB
//! BGZF block plus one tiny trailing block — doubling compressed file
//! size. Coalescing to ~4 × `BGZF_MAX_BLOCK_SIZE` per chunk amortizes
//! that per-batch trailing-block overhead while keeping chunks small
//! enough that many pool workers can compress concurrently.

use std::marker::PhantomData;

use anyhow::Result;
use fgumi_bgzf::BGZF_MAX_BLOCK_SIZE;

use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::grouping_types::CoalesceInput;
use crate::runall::engine::output_types::{RawBytes, SerializedBatch};
use crate::runall::engine::reorder::ReorderBuffer;
use crate::runall::engine::sink::InputQueue;
use crate::runall::engine::source::OutputQueue;
use crate::runall::engine::special_stage::SpecialStage;
use crate::runall::engine::stage::SequencedItem;

/// Default target chunk size = 4 BGZF max blocks. Chosen to amortize the
/// per-batch trailing-block overhead while still keeping chunks small enough
/// that many pool workers can compress concurrently.
pub const DEFAULT_TARGET_CHUNK_SIZE: usize = 4 * BGZF_MAX_BLOCK_SIZE;

/// Builder for [`Coalesce`]. Construct via [`Coalesce::builder`].
#[derive(Debug, Clone)]
pub struct CoalesceBuilder<T: CoalesceInput> {
    target_chunk_size: usize,
    _marker: PhantomData<fn() -> T>,
}

impl<T: CoalesceInput> Default for CoalesceBuilder<T> {
    fn default() -> Self {
        Self { target_chunk_size: DEFAULT_TARGET_CHUNK_SIZE, _marker: PhantomData }
    }
}

impl<T: CoalesceInput> CoalesceBuilder<T> {
    /// Start a fresh builder with all tunables at their defaults.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Minimum primary-buffer size, in uncompressed bytes, before a chunk is
    /// emitted. The actual emitted chunk may be slightly larger because
    /// full upstream batches are not split — a batch that pushes the buffer
    /// past the threshold is the final batch in that chunk.
    ///
    /// Default: [`DEFAULT_TARGET_CHUNK_SIZE`] (4 × `BGZF_MAX_BLOCK_SIZE`).
    ///
    /// # Panics
    ///
    /// Panics if `size == 0`.
    #[must_use]
    pub fn target_chunk_size(mut self, size: usize) -> Self {
        assert!(size > 0, "Coalesce target_chunk_size must be > 0");
        self.target_chunk_size = size;
        self
    }

    /// Finalize the builder.
    #[must_use]
    pub fn build(self) -> Coalesce<T> {
        Coalesce { target_chunk_size: self.target_chunk_size, _marker: PhantomData }
    }
}

/// Reorder + coalesce stage. Generic over [`CoalesceInput`]. See module docs.
#[derive(Debug, Clone)]
pub struct Coalesce<T: CoalesceInput> {
    target_chunk_size: usize,
    _marker: PhantomData<fn() -> T>,
}

impl<T: CoalesceInput> Coalesce<T> {
    /// Start constructing a `Coalesce` with default tunables.
    #[must_use]
    pub fn builder() -> CoalesceBuilder<T> {
        CoalesceBuilder::new()
    }
}

impl<T: CoalesceInput> SpecialStage for Coalesce<T> {
    type Input = T;
    type Output = SerializedBatch;

    #[allow(
        clippy::too_many_lines,
        reason = "single cohesive loop — reorder + coalesce + emit; splitting hurts clarity"
    )]
    #[tracing::instrument(name = "coalesce", skip_all)]
    fn run(
        self: Box<Self>,
        input: Box<dyn InputQueue<T>>,
        output: Box<dyn OutputQueue<SerializedBatch>>,
        cancel: CancelToken,
    ) -> Result<()> {
        let target = self.target_chunk_size;

        // Ordinal-keyed reorder buffer.
        let mut reorder: ReorderBuffer<T> = ReorderBuffer::new();

        // Accumulating buffers for the currently building chunk.
        let mut chunk_primary: Vec<u8> = Vec::with_capacity(target);
        let mut chunk_primary_records: u64 = 0;
        let mut chunk_secondary: Vec<u8> = Vec::new();
        let mut chunk_secondary_records: u64 = 0;
        let mut any_secondary_seen = false;

        // Fresh monotonic ordinal for emitted chunks.
        let mut out_ordinal: u64 = 0;

        let mut backoff = crate::runall::engine::backoff::Backoff::new();

        let emit_chunk = |chunk_primary: &mut Vec<u8>,
                          chunk_primary_records: &mut u64,
                          chunk_secondary: &mut Vec<u8>,
                          chunk_secondary_records: &mut u64,
                          any_secondary_seen: bool,
                          out_ordinal: &mut u64|
         -> Result<bool> {
            if *chunk_primary_records == 0 && !any_secondary_seen {
                return Ok(false);
            }
            let primary_data = std::mem::take(chunk_primary);
            let primary_record_count = *chunk_primary_records;
            let secondary = if any_secondary_seen {
                Some(RawBytes {
                    data: std::mem::take(chunk_secondary),
                    record_count: *chunk_secondary_records,
                })
            } else {
                None
            };
            let batch = SerializedBatch {
                primary: RawBytes { data: primary_data, record_count: primary_record_count },
                secondary,
                ordinal: *out_ordinal,
            };
            let mem =
                batch.primary.data.len() + batch.secondary.as_ref().map_or(0, |s| s.data.len());
            let item = SequencedItem::new(*out_ordinal, batch, mem);
            if output.push_until_cancelled(item, &cancel).is_err() {
                anyhow::bail!("Coalesce: output queue dropped during push");
            }
            *out_ordinal += 1;
            *chunk_primary_records = 0;
            *chunk_secondary_records = 0;
            *chunk_primary = Vec::with_capacity(target);
            Ok(true)
        };

        loop {
            if cancel.is_cancelled() {
                break;
            }

            if let Some(item) = input.pop() {
                // Reorder by the typed batch's ordinal. `CoalesceInput::ordinal`
                // is contiguous 0..N from the upstream batcher (SortStage,
                // AlignAndMerge, PositionBatchStage, etc.).
                let ordinal = item.item.ordinal();
                reorder.push(ordinal, item.item);
                backoff.reset();

                while let Some(batch) = reorder.pop_ready() {
                    let (primary, secondary) = batch.into_parts();
                    chunk_primary.extend_from_slice(&primary.data);
                    chunk_primary_records += primary.record_count;
                    if let Some(sec) = secondary {
                        any_secondary_seen = true;
                        chunk_secondary.extend_from_slice(&sec.data);
                        chunk_secondary_records += sec.record_count;
                    }

                    if chunk_primary.len() >= target {
                        emit_chunk(
                            &mut chunk_primary,
                            &mut chunk_primary_records,
                            &mut chunk_secondary,
                            &mut chunk_secondary_records,
                            any_secondary_seen,
                            &mut out_ordinal,
                        )?;
                    }
                }
            } else if input.is_drained() {
                break;
            } else {
                backoff.snooze();
            }
        }

        // Drain any remaining in-order items and emit the final partial chunk.
        while let Some(batch) = reorder.pop_ready() {
            let (primary, secondary) = batch.into_parts();
            chunk_primary.extend_from_slice(&primary.data);
            chunk_primary_records += primary.record_count;
            if let Some(sec) = secondary {
                any_secondary_seen = true;
                chunk_secondary.extend_from_slice(&sec.data);
                chunk_secondary_records += sec.record_count;
            }
        }

        if !cancel.is_cancelled() {
            emit_chunk(
                &mut chunk_primary,
                &mut chunk_primary_records,
                &mut chunk_secondary,
                &mut chunk_secondary_records,
                any_secondary_seen,
                &mut out_ordinal,
            )?;
        }

        output.close();
        Ok(())
    }

    fn name(&self) -> &'static str {
        "Coalesce"
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::runall::engine::grouping_types::{MiGroup, MiGroupBatch};
    use crate::runall::engine::memory::MemoryTracker;
    use crate::runall::engine::queue::StageQueue;

    fn raw(data: &[u8], records: u64) -> RawBytes {
        RawBytes { data: data.to_vec(), record_count: records }
    }

    fn batch(ordinal: u64, primary: &[u8], records: u64) -> SerializedBatch {
        SerializedBatch { primary: raw(primary, records), secondary: None, ordinal }
    }

    fn seq<T>(ordinal: u64, item: T, mem: usize) -> SequencedItem<T> {
        SequencedItem::new(ordinal, item, mem)
    }

    fn make_serialized_queues()
    -> (Arc<StageQueue<SerializedBatch>>, Arc<StageQueue<SerializedBatch>>) {
        let tracker = Arc::new(MemoryTracker::new(10_000_000));
        let input: Arc<StageQueue<SerializedBatch>> =
            Arc::new(StageQueue::new("coalesce_in", 64, 1_000_000, tracker.clone()));
        let output: Arc<StageQueue<SerializedBatch>> =
            Arc::new(StageQueue::new("coalesce_out", 64, 1_000_000, tracker));
        (input, output)
    }

    fn make_migroup_queues() -> (Arc<StageQueue<MiGroupBatch>>, Arc<StageQueue<SerializedBatch>>) {
        let tracker = Arc::new(MemoryTracker::new(10_000_000));
        let input: Arc<StageQueue<MiGroupBatch>> =
            Arc::new(StageQueue::new("coalesce_mi_in", 64, 1_000_000, tracker.clone()));
        let output: Arc<StageQueue<SerializedBatch>> =
            Arc::new(StageQueue::new("coalesce_mi_out", 64, 1_000_000, tracker));
        (input, output)
    }

    fn run_coalesce_serialized(
        stage: Coalesce<SerializedBatch>,
        input: Arc<StageQueue<SerializedBatch>>,
        output: Arc<StageQueue<SerializedBatch>>,
    ) {
        <Coalesce<SerializedBatch> as SpecialStage>::run(
            Box::new(stage),
            Box::new(input),
            Box::new(output),
            CancelToken::new(),
        )
        .unwrap();
    }

    fn run_coalesce_migroup(
        stage: Coalesce<MiGroupBatch>,
        input: Arc<StageQueue<MiGroupBatch>>,
        output: Arc<StageQueue<SerializedBatch>>,
    ) {
        <Coalesce<MiGroupBatch> as SpecialStage>::run(
            Box::new(stage),
            Box::new(input),
            Box::new(output),
            CancelToken::new(),
        )
        .unwrap();
    }

    #[test]
    fn test_coalesce_empty_input() {
        let (input, output) = make_serialized_queues();
        input.close();

        let stage: Coalesce<SerializedBatch> = Coalesce::builder().build();
        run_coalesce_serialized(stage, input, output.clone());

        assert!(OutputQueue::is_closed(&output));
        assert!(output.pop().is_none());
    }

    #[test]
    fn test_coalesce_single_small_batch_emits_on_drain() {
        let (input, output) = make_serialized_queues();
        let b = batch(0, b"ABCDE", 2);
        input.push(seq(0, b, 5)).unwrap();
        input.close();

        let stage: Coalesce<SerializedBatch> = Coalesce::builder().target_chunk_size(1024).build();
        run_coalesce_serialized(stage, input, output.clone());

        let popped = output.pop().expect("one chunk expected");
        assert_eq!(popped.item.ordinal, 0);
        assert_eq!(popped.item.primary.data, b"ABCDE");
        assert_eq!(popped.item.primary.record_count, 2);
        assert!(output.pop().is_none());
    }

    #[test]
    fn test_coalesce_reorders_by_ordinal_before_concatenation() {
        let (input, output) = make_serialized_queues();
        input.push(seq(1, batch(1, b"B", 1), 1)).unwrap();
        input.push(seq(2, batch(2, b"C", 1), 1)).unwrap();
        input.push(seq(0, batch(0, b"A", 1), 1)).unwrap();
        input.close();

        let stage: Coalesce<SerializedBatch> = Coalesce::builder().target_chunk_size(1024).build();
        run_coalesce_serialized(stage, input, output.clone());

        let popped = output.pop().expect("one chunk expected");
        assert_eq!(popped.item.primary.data, b"ABC");
        assert_eq!(popped.item.primary.record_count, 3);
    }

    #[test]
    fn test_coalesce_emits_when_target_size_reached() {
        let (input, output) = make_serialized_queues();
        input.push(seq(0, batch(0, b"AAAA", 1), 4)).unwrap();
        input.push(seq(1, batch(1, b"BBBB", 1), 4)).unwrap();
        input.push(seq(2, batch(2, b"CCCC", 1), 4)).unwrap();
        input.close();

        let stage: Coalesce<SerializedBatch> = Coalesce::builder().target_chunk_size(6).build();
        run_coalesce_serialized(stage, input, output.clone());

        let first = output.pop().expect("first chunk");
        assert_eq!(first.item.ordinal, 0);
        assert_eq!(first.item.primary.data, b"AAAABBBB");
        assert_eq!(first.item.primary.record_count, 2);

        let second = output.pop().expect("second chunk");
        assert_eq!(second.item.ordinal, 1);
        assert_eq!(second.item.primary.data, b"CCCC");
        assert_eq!(second.item.primary.record_count, 1);

        assert!(output.pop().is_none());
    }

    #[test]
    fn test_coalesce_propagates_secondary_stream() {
        let (input, output) = make_serialized_queues();
        let b0 =
            SerializedBatch { primary: raw(b"PP1", 1), secondary: Some(raw(b"S1", 1)), ordinal: 0 };
        let b1 =
            SerializedBatch { primary: raw(b"PP2", 1), secondary: Some(raw(b"S2", 1)), ordinal: 1 };
        input.push(seq(0, b0, 5)).unwrap();
        input.push(seq(1, b1, 5)).unwrap();
        input.close();

        let stage: Coalesce<SerializedBatch> = Coalesce::builder().target_chunk_size(1024).build();
        run_coalesce_serialized(stage, input, output.clone());

        let popped = output.pop().expect("one chunk expected");
        assert_eq!(popped.item.primary.data, b"PP1PP2");
        assert_eq!(popped.item.primary.record_count, 2);
        let sec = popped.item.secondary.expect("secondary present");
        assert_eq!(sec.data, b"S1S2");
        assert_eq!(sec.record_count, 2);
    }

    #[test]
    fn test_coalesce_builder_defaults() {
        let builder: CoalesceBuilder<SerializedBatch> = Coalesce::builder();
        assert_eq!(builder.target_chunk_size, DEFAULT_TARGET_CHUNK_SIZE);
    }

    /// `Coalesce<MiGroupBatch>` flattens MI-grouped records into concat bytes.
    #[test]
    fn test_coalesce_migroup_flattens_and_coalesces() {
        let (input, output) = make_migroup_queues();

        let g0 = MiGroup { data: b"\x01A".to_vec(), record_count: 1, mi: 0, local_mi: None };
        let g1 = MiGroup { data: b"\x01B".to_vec(), record_count: 1, mi: 1, local_mi: None };
        let mg0 = MiGroupBatch {
            groups: vec![g0, g1],
            ordinal: 0,
            position_key: (0, 0),
            assign_tag: *b"MI",
        };

        let g2 = MiGroup { data: b"\x01C".to_vec(), record_count: 1, mi: 2, local_mi: None };
        let mg1 =
            MiGroupBatch { groups: vec![g2], ordinal: 1, position_key: (1, 0), assign_tag: *b"MI" };

        input.push(seq(1, mg1, 2)).unwrap();
        input.push(seq(0, mg0, 4)).unwrap();
        input.close();

        let stage: Coalesce<MiGroupBatch> = Coalesce::builder().target_chunk_size(1024).build();
        run_coalesce_migroup(stage, input, output.clone());

        let popped = output.pop().expect("one chunk");
        // Expected primary: g0.data ++ g1.data ++ g2.data in MI order, ordinal 0.
        assert_eq!(popped.item.primary.data, b"\x01A\x01B\x01C");
        assert_eq!(popped.item.primary.record_count, 3);
        assert!(popped.item.secondary.is_none());
        assert!(output.pop().is_none());
    }
}
