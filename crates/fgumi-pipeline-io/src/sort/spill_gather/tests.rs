//! Unit tests for `SpillGather`'s `StepCtx`-free core: chunk fan-out into raw
//! blocks (`frame_chunk_into_blocks`) and event staging (`stage_event`).

use super::*;
use crate::sort::protocol::MemoryChunkErased;
use fgumi_sort::{InMemoryChunk, RawCoordinateKey};
use proptest::prelude::*;

/// Build a coordinate chunk from raw payloads with distinct keys, so the framed
/// blocks exercise key serialization (`RawCoordinateKey` is embedded, so no key
/// prefix is written — but the chunk path is still the production one).
fn coord_chunk(payloads: Vec<Vec<u8>>) -> MemoryChunkErased {
    let recs = payloads
        .into_iter()
        .enumerate()
        .map(|(i, b)| (RawCoordinateKey { sort_key: i as u64 }, b))
        .collect();
    MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(recs))
}

/// Concatenate the framed-record bytes of all blocks (drops the per-block
/// boundaries) — the decompressed-stream-equivalent the readers see.
fn concat(blocks: &[Vec<u8>]) -> Vec<u8> {
    blocks.iter().flat_map(|b| b.iter().copied()).collect()
}

/// Drive `SpillGather` over a sequence of input events exactly as `try_run`
/// does — stage each event, then fully frame any active chunk (draining produced
/// blocks) before staging the next — and return every emitted event in order.
/// This preserves the dense-ordinal invariant (a chunk is fully framed, and its
/// block ordinals minted, before the next event is staged).
fn drive(step: &mut SpillGather, events: Vec<SortChunkEvent>) -> Vec<SpillBlockEvent> {
    let mut out = Vec::new();
    for event in events {
        step.stage_event(event);
        // Frame the active chunk to completion, collecting blocks as we go.
        while step.active.is_some() {
            step.produce_blocks().unwrap();
            while let Some(ev) = step.pending.pop_front() {
                out.push(ev);
            }
        }
        // Drain any passthrough (Residual / AllAnnounced) events.
        while let Some(ev) = step.pending.pop_front() {
            out.push(ev);
        }
    }
    out
}

#[test]
fn frame_empty_chunk_yields_no_blocks() {
    let chunk = coord_chunk(Vec::new());
    let blocks = frame_chunk_into_blocks(&chunk, BGZF_MAX_BLOCK_SIZE).unwrap();
    assert!(blocks.is_empty(), "empty chunk must yield zero blocks");
}

#[test]
fn frame_cuts_blocks_at_threshold_without_splitting_records() {
    // 10 records of 100 payload bytes each; embedded key adds only the 4-byte
    // length prefix → 104 framed bytes/record. A 250-byte block threshold cuts
    // *before* the record that would exceed 250, so each block holds 2 records
    // (2×104=208 ≤ 250; a 3rd would be 312 > 250) → 5 blocks.
    let chunk = coord_chunk((0u8..10).map(|i| vec![i; 100]).collect());
    let blocks = frame_chunk_into_blocks(&chunk, 250).unwrap();
    assert_eq!(blocks.len(), 5, "250-byte threshold must pack 2 records/block over 10 records");
    // Every block stays within the threshold and holds a whole number of
    // 104-byte records (no record split).
    for b in &blocks {
        assert!(b.len() <= 250, "block exceeded threshold: {} bytes", b.len());
        assert_eq!(b.len() % 104, 0, "block boundary split a record: {} bytes", b.len());
    }
    // Reassembling the blocks reproduces a single 10-record stream.
    assert_eq!(concat(&blocks).len(), 10 * 104);
}

proptest! {
    /// The "no record split, ≤ threshold, exact reconstruction" packing invariant
    /// of `frame_chunk_into_blocks` holds for arbitrary record sizes and
    /// thresholds — far more of the space than the hand-picked cases above.
    #[test]
    fn frame_chunk_into_blocks_packing_invariant(
        payloads in prop::collection::vec(prop::collection::vec(any::<u8>(), 0..40), 0..30),
        threshold in 1usize..512,
    ) {
        let chunk = coord_chunk(payloads.clone());
        let blocks = frame_chunk_into_blocks(&chunk, threshold).unwrap();

        // Content is threshold-independent: block boundaries move, bytes do not.
        // Compare against both the all-in-one framing and the one-record-per-block
        // framing (threshold 1), pinning exact reconstruction of the framed stream.
        let one_block = frame_chunk_into_blocks(&chunk, BGZF_MAX_BLOCK_SIZE).unwrap();
        let per_record = frame_chunk_into_blocks(&chunk, 1).unwrap();
        prop_assert_eq!(concat(&blocks), concat(&one_block));
        prop_assert_eq!(concat(&blocks), concat(&per_record));

        // Empty chunk → no blocks; any record (even a zero-length payload, which
        // still frames its length prefix) → at least one, and never an empty one.
        prop_assert_eq!(blocks.is_empty(), payloads.is_empty());
        prop_assert!(blocks.iter().all(|b| !b.is_empty()), "no empty blocks");

        // No record is split: a block may exceed `threshold` only when a single
        // record is itself larger than it (the "cut before exceeding" rule can't
        // shrink one record), so every block is bounded by
        // `max(threshold, largest single framed record)`.
        let max_record = per_record.iter().map(Vec::len).max().unwrap_or(0);
        let bound = threshold.max(max_record);
        for b in &blocks {
            prop_assert!(b.len() <= bound, "block {} exceeds bound {bound}", b.len());
        }
    }
}

#[test]
fn frame_one_block_when_under_threshold() {
    let chunk = coord_chunk((0u8..5).map(|i| vec![i; 100]).collect());
    let blocks = frame_chunk_into_blocks(&chunk, BGZF_MAX_BLOCK_SIZE).unwrap();
    assert_eq!(blocks.len(), 1, "5 small records fit one block");
    assert_eq!(blocks[0].len(), 5 * 104);
}

#[test]
fn stage_spill_mints_dense_ordinals_and_marks_last_block() {
    let mut step = SpillGather::new(1 << 20);
    step.block_size = 250; // force several blocks
    let chunk = coord_chunk((0u8..10).map(|i| vec![i; 100]).collect());
    let events = drive(
        &mut step,
        vec![SortChunkEvent::Spill { seq: 3, chunk, records_ingested_so_far: 10 }],
    );

    assert!(events.len() > 1, "expected multiple block events");
    // Ordinals are dense 0..n.
    for (i, ev) in events.iter().enumerate() {
        assert_eq!(ev.ordinal(), i as u64, "ordinal not dense at index {i}");
    }
    // Exactly the final block is flagged is_last_in_file; all carry file_id=3.
    for (i, ev) in events.iter().enumerate() {
        let SpillBlockEvent::Block { file_id, is_last_in_file, records_ingested_so_far, .. } = ev
        else {
            panic!("expected Block event");
        };
        assert_eq!(*file_id, 3);
        assert_eq!(*records_ingested_so_far, 10);
        assert_eq!(*is_last_in_file, i == events.len() - 1, "is_last wrong at {i}");
    }
    // The chunk was freed once fully framed.
    assert!(step.active.is_none(), "active chunk must be cleared after framing");
}

#[test]
fn produce_blocks_keeps_pending_bounded() {
    // With many small records and a tiny block size, a single produce_blocks call
    // must not materialize the whole chunk — it tops up at most
    // MAX_EVENTS_PER_LOCK blocks.
    let mut step = SpillGather::new(1 << 20);
    step.block_size = 120; // ~1 record/block over 104-byte records
    let chunk = coord_chunk((0u8..100).map(|i| vec![i; 100]).collect());
    step.stage_event(SortChunkEvent::Spill { seq: 0, chunk, records_ingested_so_far: 100 });
    step.produce_blocks().unwrap();
    assert_eq!(step.pending.len(), 8, "one produce call tops up to MAX_EVENTS_PER_LOCK blocks");
    assert!(step.active.is_some(), "chunk still mid-framing (not fully drained in one call)");
}

#[test]
fn stage_residual_and_announced_pass_through_with_ordinals() {
    let mut step = SpillGather::new(1 << 20);
    // A spill (2 small records → 1 block, ordinal 0), then residual, then
    // AllAnnounced — ordinals must stay dense across the variant boundary, with
    // the spill's block ordinal minted *before* the later events (the chunk is
    // fully framed before the next event is staged).
    let events = drive(
        &mut step,
        vec![
            SortChunkEvent::Spill {
                seq: 0,
                chunk: coord_chunk(vec![vec![1u8; 10], vec![2u8; 10]]),
                records_ingested_so_far: 2,
            },
            SortChunkEvent::Residual {
                chunk: coord_chunk(vec![vec![3u8; 10]]),
                records_ingested_so_far: 3,
            },
            SortChunkEvent::AllAnnounced { slot_count: 1, memory_chunk_count: 1, total_records: 3 },
        ],
    );

    let ords: Vec<u64> = events.iter().map(SpillBlockEvent::ordinal).collect();
    assert_eq!(ords, vec![0, 1, 2], "ordinals must be dense across variants");

    assert!(matches!(events[0], SpillBlockEvent::Block { is_last_in_file: true, .. }));
    assert!(matches!(events[1], SpillBlockEvent::Residual { .. }));
    assert!(matches!(
        events[2],
        SpillBlockEvent::AllAnnounced {
            slot_count: 1,
            memory_chunk_count: 1,
            total_records: 3,
            ..
        }
    ));
}
