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

/// Build a coordinate chunk with explicit (already-sorted) `sort_key`s and dummy
/// record bodies, so the run-former's min/max bounds are controllable. `keys`
/// must be non-decreasing (chunks are pre-sorted at seal).
fn coord_chunk_keyed(keys: &[u64]) -> MemoryChunkErased {
    let recs = keys.iter().map(|&k| (RawCoordinateKey { sort_key: k }, vec![0u8; 8])).collect();
    MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(recs))
}

/// The `(file_id, is_last_in_file)` of every spill block — a planned `Range` or
/// the run-final raw `Block` — in order. Panics on a passthrough event, so
/// callers separate spills from passthroughs.
fn block_files(events: &[SpillBlockEvent]) -> Vec<(u32, bool)> {
    events
        .iter()
        .map(|e| match e {
            SpillBlockEvent::Range { file_id, is_last_in_file, .. }
            | SpillBlockEvent::Block { file_id, is_last_in_file, .. } => {
                (*file_id, *is_last_in_file)
            }
            _ => panic!("expected only Range/Block events (a passthrough event was present)"),
        })
        .collect()
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
///
/// After the last event, it closes the final run by flushing any withheld block
/// with `is_last_in_file = true` — mirroring `try_run`'s drain path (in
/// production `AllAnnounced` closes it; a `drive` without one relies on this).
fn drive(step: &mut SpillGather, events: Vec<SortChunkEvent>) -> Vec<SpillBlockEvent> {
    let mut out = Vec::new();
    for event in events {
        step.stage_event(event);
        // Frame the active chunk to completion, collecting blocks as we go.
        while step.active.is_some() {
            step.plan_ranges().unwrap();
            while let Some(ev) = step.pending.pop_front() {
                out.push(ev);
            }
        }
        // Drain any passthrough (Residual / AllAnnounced) events.
        while let Some(ev) = step.pending.pop_front() {
            out.push(ev);
        }
    }
    // Close the final run (no-op if a Residual/AllAnnounced already flushed it).
    step.flush_held_block(true);
    while let Some(ev) = step.pending.pop_front() {
        out.push(ev);
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
    // Every block but the last is a planned range; the last is the raw
    // run-final block, and only it is flagged is_last_in_file; all carry
    // file_id=3.
    let last = events.len() - 1;
    for (i, ev) in events.iter().enumerate() {
        let (file_id, is_last_in_file, records_ingested_so_far) = match ev {
            SpillBlockEvent::Range {
                file_id, is_last_in_file, records_ingested_so_far, ..
            } => {
                assert!(i < last, "the run-final block must be a raw Block, not a Range");
                (file_id, is_last_in_file, records_ingested_so_far)
            }
            SpillBlockEvent::Block {
                file_id, is_last_in_file, records_ingested_so_far, ..
            } => {
                assert_eq!(i, last, "only the run-final block is framed by the gather");
                (file_id, is_last_in_file, records_ingested_so_far)
            }
            _ => panic!("expected spill blocks"),
        };
        assert_eq!(*file_id, 3);
        assert_eq!(*records_ingested_so_far, 10);
        assert_eq!(*is_last_in_file, i == last, "is_last wrong at {i}");
    }
    // The chunk was freed once fully framed.
    assert!(step.active.is_none(), "active chunk must be cleared after framing");
}

#[test]
fn plan_ranges_keeps_pending_bounded() {
    // With many small records and a tiny block size, a single plan_ranges call
    // must not materialize the whole chunk — it tops up at most
    // MAX_EVENTS_PER_LOCK blocks.
    let mut step = SpillGather::new(1 << 20);
    step.block_size = 120; // ~1 record/block over 104-byte records
    let chunk = coord_chunk((0u8..100).map(|i| vec![i; 100]).collect());
    step.stage_event(SortChunkEvent::Spill { seq: 0, chunk, records_ingested_so_far: 100 });
    step.plan_ranges().unwrap();
    assert_eq!(step.pending.len(), 8, "one plan call tops up to MAX_EVENTS_PER_LOCK ranges");
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

// ── Step wiring + ordinal minting ────────────────────────────────────────────

#[test]
fn profile_runs_off_pool_on_the_coordination_driver() {
    let step = SpillGather::new(4096);
    let profile = step.profile();
    assert_eq!(profile.name, "SpillGather");
    // Detached keeps the serial framing + ordinal minting off a pool worker, so
    // it cannot starve the parallel spill compressors downstream.
    assert_eq!(profile.kind, StepKind::Detached);
    assert!(!profile.sticky);
    assert_eq!(profile.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    match profile.output_queues.as_slice() {
        [QueueSpec::ByteBounded { limit_bytes }] => assert_eq!(*limit_bytes, 4096),
        other => panic!("expected one byte-bounded queue, got {other:?}"),
    }
    assert_eq!(step.detached_group(), DetachedGroup::Shared(crate::sort::SORT_COORD_GROUP));
}

#[test]
fn ordinals_are_dense_and_monotonic_across_event_kinds() {
    // Downstream `ByItemOrdinal` reordering depends on a gap-free sequence, and
    // the counter is shared by every emitted event, not per-variant.
    let mut step = SpillGather::new(1 << 20);

    step.stage_event(SortChunkEvent::Residual {
        chunk: coord_chunk(vec![vec![1u8; 8]]),
        records_ingested_so_far: 1,
    });
    step.stage_event(SortChunkEvent::AllAnnounced {
        slot_count: 1,
        memory_chunk_count: 1,
        total_records: 1,
    });
    step.stage_event(SortChunkEvent::Residual {
        chunk: coord_chunk(vec![vec![2u8; 8]]),
        records_ingested_so_far: 2,
    });

    let ordinals: Vec<u64> = step.pending.iter().map(SpillBlockEvent::ordinal).collect();
    assert_eq!(ordinals, vec![0, 1, 2], "ordinals must be dense across variants");
    assert_eq!(step.next_ordinal, 3);
}

#[test]
fn an_empty_spill_chunk_is_skipped_rather_than_opening_a_file() {
    // A zero-block file would never receive an `is_last_in_file` terminator, so
    // `SpillWrite` would be left with a dangling open file forever.
    let mut step = SpillGather::new(1 << 20);
    step.stage_event(SortChunkEvent::Spill {
        seq: 0,
        chunk: coord_chunk(Vec::new()),
        records_ingested_so_far: 0,
    });
    assert!(step.active.is_none(), "an empty chunk must not become the active spill");
    assert!(step.pending.is_empty(), "and must emit nothing");
    assert_eq!(step.next_ordinal, 0, "and must not consume an ordinal");
}

// ── Run extension (spill coalescing) ─────────────────────────────────────────

/// The `slot_count` an `AllAnnounced` event carries (its rewritten run count).
/// Panics unless the last event is an `AllAnnounced`.
fn announced_slot_count(events: &[SpillBlockEvent]) -> u32 {
    match events.last().expect("at least one event") {
        SpillBlockEvent::AllAnnounced { slot_count, .. } => *slot_count,
        _ => panic!("last event must be AllAnnounced"),
    }
}

/// Split the emitted stream into its `Block` events and the trailing
/// `AllAnnounced` (asserting there are no other event kinds).
fn blocks_and_announced(events: &[SpillBlockEvent]) -> (Vec<(u32, bool)>, u32) {
    let (last, blocks) = events.split_last().expect("at least an AllAnnounced");
    let slot_count = match last {
        SpillBlockEvent::AllAnnounced { slot_count, .. } => *slot_count,
        _ => panic!("last event must be AllAnnounced"),
    };
    (block_files(blocks), slot_count)
}

#[test]
fn contiguous_chunks_coalesce_into_one_run() {
    // Two chunks, the second starting strictly after the first ends → one run:
    // both blocks share file_id 0, and only the final block closes the file.
    let mut step = SpillGather::new(1 << 20);
    let events = drive(
        &mut step,
        vec![
            SortChunkEvent::Spill {
                seq: 0,
                chunk: coord_chunk_keyed(&[10, 20]),
                records_ingested_so_far: 2,
            },
            SortChunkEvent::Spill {
                seq: 1,
                chunk: coord_chunk_keyed(&[30, 40]),
                records_ingested_so_far: 4,
            },
            SortChunkEvent::AllAnnounced { slot_count: 2, memory_chunk_count: 0, total_records: 4 },
        ],
    );

    let (blocks, slot_count) = blocks_and_announced(&events);
    assert_eq!(blocks, vec![(0, false), (0, true)], "both blocks in one file, last closes it");
    assert_eq!(slot_count, 1, "slot_count is the RUN count, not the chunk count");
    assert_eq!(step.runs_written, 1);
    // Ordinals stay dense across the withheld-block flush.
    let ords: Vec<u64> = events.iter().map(SpillBlockEvent::ordinal).collect();
    assert_eq!(ords, vec![0, 1, 2]);
}

#[test]
fn exact_boundary_tie_extends_the_run_for_content_keys() {
    // Coordinate keys are content-based: an exact boundary tie (20 == 20) still
    // satisfies precedes-or-equal, so the chunks coalesce.
    let mut step = SpillGather::new(1 << 20);
    let events = drive(
        &mut step,
        vec![
            SortChunkEvent::Spill {
                seq: 0,
                chunk: coord_chunk_keyed(&[10, 20]),
                records_ingested_so_far: 2,
            },
            SortChunkEvent::Spill {
                seq: 1,
                chunk: coord_chunk_keyed(&[20, 30]),
                records_ingested_so_far: 4,
            },
            SortChunkEvent::AllAnnounced { slot_count: 2, memory_chunk_count: 0, total_records: 4 },
        ],
    );
    let (blocks, slot_count) = blocks_and_announced(&events);
    assert_eq!(blocks, vec![(0, false), (0, true)]);
    assert_eq!(slot_count, 1);
}

#[test]
fn non_contiguous_chunk_starts_a_new_run() {
    // The second chunk's min (20) falls inside the first run (…40), so it cannot
    // extend it: two runs, each a self-contained file.
    let mut step = SpillGather::new(1 << 20);
    let events = drive(
        &mut step,
        vec![
            SortChunkEvent::Spill {
                seq: 0,
                chunk: coord_chunk_keyed(&[10, 40]),
                records_ingested_so_far: 2,
            },
            SortChunkEvent::Spill {
                seq: 1,
                chunk: coord_chunk_keyed(&[20, 30]),
                records_ingested_so_far: 4,
            },
            SortChunkEvent::AllAnnounced { slot_count: 2, memory_chunk_count: 0, total_records: 4 },
        ],
    );
    let (blocks, slot_count) = blocks_and_announced(&events);
    assert_eq!(blocks, vec![(0, true), (1, true)], "each chunk closes its own file");
    assert_eq!(slot_count, 2);
    assert_eq!(step.runs_written, 2);
}

#[test]
fn many_presorted_chunks_collapse_to_a_single_run() {
    // The already-sorted-input case: every chunk is contiguous, so N chunks
    // become one run — the whole point of run extension.
    let mut step = SpillGather::new(1 << 20);
    let events = drive(
        &mut step,
        vec![
            SortChunkEvent::Spill {
                seq: 0,
                chunk: coord_chunk_keyed(&[0, 1]),
                records_ingested_so_far: 2,
            },
            SortChunkEvent::Spill {
                seq: 1,
                chunk: coord_chunk_keyed(&[2, 3]),
                records_ingested_so_far: 4,
            },
            SortChunkEvent::Spill {
                seq: 2,
                chunk: coord_chunk_keyed(&[4, 5]),
                records_ingested_so_far: 6,
            },
            SortChunkEvent::AllAnnounced { slot_count: 3, memory_chunk_count: 0, total_records: 6 },
        ],
    );
    let (blocks, slot_count) = blocks_and_announced(&events);
    assert_eq!(blocks, vec![(0, false), (0, false), (0, true)], "one file, only the last closes");
    assert_eq!(slot_count, 1);
}

#[test]
fn an_empty_chunk_does_not_close_the_open_run() {
    // An empty spill is a pure no-op — it neither closes the run nor breaks the
    // contiguity of the chunks on either side of it.
    let mut step = SpillGather::new(1 << 20);
    let events = drive(
        &mut step,
        vec![
            SortChunkEvent::Spill {
                seq: 0,
                chunk: coord_chunk_keyed(&[10, 20]),
                records_ingested_so_far: 2,
            },
            SortChunkEvent::Spill {
                seq: 1,
                chunk: coord_chunk_keyed(&[]),
                records_ingested_so_far: 2,
            },
            SortChunkEvent::Spill {
                seq: 2,
                chunk: coord_chunk_keyed(&[30, 40]),
                records_ingested_so_far: 4,
            },
            SortChunkEvent::AllAnnounced { slot_count: 3, memory_chunk_count: 0, total_records: 4 },
        ],
    );
    let (blocks, slot_count) = blocks_and_announced(&events);
    assert_eq!(blocks, vec![(0, false), (0, true)], "the empty chunk contributes nothing");
    assert_eq!(slot_count, 1, "still one run — the empty chunk did not split it");
}

#[test]
fn a_residual_closes_the_open_run_before_a_later_spill() {
    // A residual (in-memory tail) ends the spill run; a spill after it starts a
    // fresh run even if its keys would otherwise have been contiguous.
    let mut step = SpillGather::new(1 << 20);
    let events = drive(
        &mut step,
        vec![
            SortChunkEvent::Spill {
                seq: 0,
                chunk: coord_chunk_keyed(&[10, 20]),
                records_ingested_so_far: 2,
            },
            SortChunkEvent::Residual {
                chunk: coord_chunk_keyed(&[25]),
                records_ingested_so_far: 3,
            },
            SortChunkEvent::Spill {
                seq: 2,
                chunk: coord_chunk_keyed(&[30, 40]),
                records_ingested_so_far: 5,
            },
            SortChunkEvent::AllAnnounced { slot_count: 2, memory_chunk_count: 1, total_records: 5 },
        ],
    );
    // Stream: Block(file0,last), Residual, Block(file2,last), AllAnnounced.
    assert!(matches!(events[0], SpillBlockEvent::Block { file_id: 0, is_last_in_file: true, .. }));
    assert!(matches!(events[1], SpillBlockEvent::Residual { .. }));
    assert!(matches!(events[2], SpillBlockEvent::Block { file_id: 2, is_last_in_file: true, .. }));
    assert_eq!(announced_slot_count(&events), 2, "two runs, split by the residual");
    let ords: Vec<u64> = events.iter().map(SpillBlockEvent::ordinal).collect();
    assert_eq!(ords, vec![0, 1, 2, 3], "ordinals dense across the residual boundary");
}

#[test]
fn run_formation_summary_matches_the_owned_engine_wording() {
    // Three contiguous chunks → one run, so 2 of the 3 chunks extended it. The
    // line must read exactly like the retired engine's run-formation line.
    let mut step = SpillGather::new(1 << 20);
    drive(
        &mut step,
        vec![
            SortChunkEvent::Spill {
                seq: 0,
                chunk: coord_chunk_keyed(&[0, 1]),
                records_ingested_so_far: 2,
            },
            SortChunkEvent::Spill {
                seq: 1,
                chunk: coord_chunk_keyed(&[2, 3]),
                records_ingested_so_far: 4,
            },
            SortChunkEvent::Spill {
                seq: 2,
                chunk: coord_chunk_keyed(&[4, 5]),
                records_ingested_so_far: 6,
            },
            SortChunkEvent::AllAnnounced { slot_count: 3, memory_chunk_count: 0, total_records: 6 },
        ],
    );
    assert_eq!(step.chunks_spilled, 3);
    assert_eq!(step.runs_written, 1);
    assert_eq!(
        step.run_formation_summary().as_deref(),
        Some("Spill runs: 1 from 3 chunks (2 extended an existing run)"),
    );
}

#[test]
fn run_formation_summary_counts_each_run_when_nothing_coalesces() {
    // Reverse-ordered chunks never extend: 3 chunks → 3 runs, 0 extended.
    let mut step = SpillGather::new(1 << 20);
    drive(
        &mut step,
        vec![
            SortChunkEvent::Spill {
                seq: 0,
                chunk: coord_chunk_keyed(&[40, 50]),
                records_ingested_so_far: 2,
            },
            SortChunkEvent::Spill {
                seq: 1,
                chunk: coord_chunk_keyed(&[20, 30]),
                records_ingested_so_far: 4,
            },
            SortChunkEvent::Spill {
                seq: 2,
                chunk: coord_chunk_keyed(&[0, 10]),
                records_ingested_so_far: 6,
            },
            SortChunkEvent::AllAnnounced { slot_count: 3, memory_chunk_count: 0, total_records: 6 },
        ],
    );
    assert_eq!(step.chunks_spilled, 3);
    assert_eq!(step.runs_written, 3);
    assert_eq!(
        step.run_formation_summary().as_deref(),
        Some("Spill runs: 3 from 3 chunks (0 extended an existing run)"),
    );
}

#[test]
fn run_formation_summary_is_none_when_nothing_spilled() {
    // A residual-only sort (no spills) reports no run-formation line: the line
    // is logged only when something spilled.
    let mut step = SpillGather::new(1 << 20);
    drive(
        &mut step,
        vec![
            SortChunkEvent::Residual { chunk: coord_chunk_keyed(&[1]), records_ingested_so_far: 1 },
            SortChunkEvent::AllAnnounced { slot_count: 0, memory_chunk_count: 1, total_records: 1 },
        ],
    );
    assert_eq!(step.chunks_spilled, 0);
    assert_eq!(step.run_formation_summary(), None);
}

#[test]
fn the_run_final_block_is_withheld_until_the_run_is_closed() {
    // Frame a single spill chunk with no terminating event yet: its final block
    // is withheld (the run is still open), which is exactly the state the drain
    // guard fails closed on — `held_last_block.is_held()` at drain.
    let mut step = SpillGather::new(1 << 20);
    step.stage_event(SortChunkEvent::Spill {
        seq: 0,
        chunk: coord_chunk_keyed(&[10, 20]),
        records_ingested_so_far: 2,
    });
    step.plan_ranges().unwrap();
    assert!(step.active.is_none(), "the one-block chunk is fully framed");
    assert!(step.pending.is_empty(), "the sole (run-final) block is withheld, not emitted");
    assert!(
        step.held_last_block.is_held(),
        "a run-final block is withheld until close — a drain here would trip the guard"
    );

    // Closing the run (AllAnnounced) flushes the withheld block with is_last=true
    // and leaves nothing held, so the drain guard is clear.
    step.stage_event(SortChunkEvent::AllAnnounced {
        slot_count: 1,
        memory_chunk_count: 0,
        total_records: 2,
    });
    assert!(!step.held_last_block.is_held(), "closing the run flushes the withheld block");
    assert!(
        matches!(
            step.pending.front(),
            Some(SpillBlockEvent::Block { is_last_in_file: true, file_id: 0, .. })
        ),
        "the flushed block closes its file"
    );
}

#[test]
fn a_non_empty_spill_chunk_becomes_active_without_emitting_yet() {
    let mut step = SpillGather::new(1 << 20);
    step.stage_event(SortChunkEvent::Spill {
        seq: 7,
        chunk: coord_chunk(vec![vec![3u8; 32], vec![4u8; 32]]),
        records_ingested_so_far: 2,
    });
    let active = step.active.as_ref().expect("chunk must be parked for incremental framing");
    assert_eq!(active.file_id, 7, "file_id comes from the spill seq, not write order");
    assert_eq!(active.next_idx, 0);
    assert!(step.pending.is_empty(), "planning happens in plan_ranges, not stage_event");
}

/// The `TemplateCoordinate` variant is framed by a *different* function than the
/// other three: `frame_record_at` routes it to `c.frame_record_into`, while
/// `Coordinate` / `QuerynameLex` / `QuerynameNatural` share `frame_keyed_record_into`.
/// Every other test here uses `coord_chunk`, so that branch never ran. A layout
/// divergence produces spill files that `SortMerge` misreads — wrong output, not a
/// crash — so the fourth variant needs its own framing coverage.
#[test]
fn template_coordinate_chunks_frame_through_their_own_path() {
    use fgumi_sort::{TemplateKey24, TemplateMemChunk};

    let payloads = [vec![0xA1u8; 24], vec![0xB2u8; 40], vec![0xC3u8; 8]];
    let recs = payloads.iter().map(|b| (TemplateKey24::default(), b.clone())).collect::<Vec<_>>();
    let chunk = MemoryChunkErased::TemplateCoordinate(TemplateMemChunk::K24(
        InMemoryChunk::from_owned_records(recs),
    ));

    // Frame the whole chunk; every record must be emitted exactly once.
    let mut blocks = Vec::new();
    let mut next = 0usize;
    while next < chunk.len() {
        let mut out = Vec::new();
        let advanced =
            frame_one_block(&chunk, next, BGZF_MAX_BLOCK_SIZE, &mut out).expect("template framing");
        assert!(advanced > next, "framing must make progress on every call");
        next = advanced;
        blocks.push(out);
    }
    assert_eq!(next, payloads.len(), "every template record is framed");

    // Each payload appears in the framed stream, so the template layout carries
    // the record bodies through unchanged.
    let framed = concat(&blocks);
    for (i, p) in payloads.iter().enumerate() {
        assert!(
            framed.windows(p.len()).any(|w| w == p.as_slice()),
            "record {i}'s body must survive template framing"
        );
    }
}

// ── Range planning (the cut rule, run-final framing, arena lifetime) ─────────

/// A chunk of `kind` whose records have the given body lengths (bodies are
/// opaque to framing: only the template lanes write a key, which defaults).
fn chunk_of_kind(kind: SpillKeyKind, lens: &[usize]) -> MemoryChunkErased {
    use fgumi_sort::{CbKey32, QnKey, TemplateKey, TemplateKey24, TemplateMemChunk, TertKey32};
    fn records<K: Default>(lens: &[usize]) -> Vec<(K, Vec<u8>)> {
        lens.iter()
            .enumerate()
            .map(|(i, &n)| (K::default(), vec![u8::try_from(i % 251).unwrap(); n]))
            .collect()
    }
    match kind {
        SpillKeyKind::Coordinate => {
            MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(records::<
                RawCoordinateKey,
            >(lens)))
        }
        SpillKeyKind::QuerynameLex => MemoryChunkErased::QuerynameLex(
            InMemoryChunk::from_owned_records(records::<QnKey>(lens)),
        ),
        SpillKeyKind::QuerynameNatural => MemoryChunkErased::QuerynameNatural(
            InMemoryChunk::from_owned_records(records::<QnKey>(lens)),
        ),
        SpillKeyKind::TemplateK24 => MemoryChunkErased::TemplateCoordinate(TemplateMemChunk::K24(
            InMemoryChunk::from_owned_records(records::<TemplateKey24>(lens)),
        )),
        SpillKeyKind::TemplateCb32 => MemoryChunkErased::TemplateCoordinate(
            TemplateMemChunk::Cb32(InMemoryChunk::from_owned_records(records::<CbKey32>(lens))),
        ),
        SpillKeyKind::TemplateTert32 => MemoryChunkErased::TemplateCoordinate(
            TemplateMemChunk::Tert32(InMemoryChunk::from_owned_records(records::<TertKey32>(lens))),
        ),
        SpillKeyKind::TemplateK40 => MemoryChunkErased::TemplateCoordinate(TemplateMemChunk::K40(
            InMemoryChunk::from_owned_records(records::<TemplateKey>(lens)),
        )),
    }
}

/// The planner's cut rule reproduces `frame_one_block`'s blocks exactly — the
/// same `[start, end)` and the same framed size — for every spill key kind
/// (all four template lanes included), block
/// size and an oversized record, so a planned range frames to the block the
/// serial framer would have produced.
#[rstest::rstest]
fn plan_next_block_matches_frame_one_block(
    #[values(
        SpillKeyKind::Coordinate,
        SpillKeyKind::QuerynameLex,
        SpillKeyKind::QuerynameNatural,
        SpillKeyKind::TemplateK24,
        SpillKeyKind::TemplateCb32,
        SpillKeyKind::TemplateTert32,
        SpillKeyKind::TemplateK40
    )]
    kind: SpillKeyKind,
    #[values(1usize, 100, 64 * 1024)] block_size: usize,
) {
    let mut lens: Vec<usize> = (0..300).map(|i| 1 + (i * 37) % 300).collect();
    lens.insert(120, 70_000); // larger than every block size
    lens.insert(7, 0);
    let chunk = chunk_of_kind(kind, &lens);
    let mut start = 0;
    while start < chunk.len() {
        let (end, framed) = plan_next_block(&chunk, start, block_size);
        let mut out = Vec::new();
        let framed_end = frame_one_block(&chunk, start, block_size, &mut out).unwrap();
        assert_eq!((end, framed as usize), (framed_end, out.len()), "{kind:?} block at {start}");
        assert!(end > start, "every block holds at least one record");
        start = end;
    }
}

/// A run-final block framed differently from the planner's cut fails the
/// sort, naming both ranges, rather than dropping the skipped records.
#[test]
fn a_run_final_framing_disagreement_is_an_error() {
    assert!(super::check_run_final_block(5, (9, 400), (9, 400)).is_ok());
    for framed in [(8, 300), (9, 399), (10, 500)] {
        let err = super::check_run_final_block(5, (9, 400), framed).unwrap_err();
        assert_eq!(
            err.to_string(),
            format!(
                "SpillGather: the run-final block framed records [5, {}) ({} bytes) but the \
                 planner cut [5, 9) (400 bytes); refusing to drop records",
                framed.0, framed.1
            )
        );
    }
}

/// A run smaller than one block has only its final block: it is framed
/// serially and withheld, and no range is planned.
#[test]
fn single_block_run_is_framed_serially() {
    let mut step = SpillGather::new(1 << 20);
    step.stage_event(SortChunkEvent::Spill {
        seq: 0,
        chunk: coord_chunk(vec![vec![1u8; 30], vec![2u8; 40]]),
        records_ingested_so_far: 2,
    });
    step.plan_ranges().unwrap();
    assert!(step.pending.is_empty(), "no range is planned for a one-block run");
    assert!(step.active.is_none(), "the step released the chunk");
    let held = step.held_last_block.take().expect("the run-final block is withheld");
    assert_eq!(
        held.bytes,
        concat(
            &frame_chunk_into_blocks(
                &coord_chunk(vec![vec![1u8; 30], vec![2u8; 40]]),
                BGZF_MAX_BLOCK_SIZE
            )
            .unwrap()
        )
    );
}

/// With a capacity-1 arena pool, run k's arena comes back only once every
/// range planned from it has been compressed — and does come back then, so
/// run k+1 can never wait on a chunk the planner still holds.
#[test]
fn arena_returns_only_after_the_last_range_is_compressed() {
    use fgumi_sort::{ArenaPool, RecordRef, SpillCodec};

    let pool = ArenaPool::new(1, 1 << 20);
    let mut arena = pool.try_acquire().expect("the pool's one arena");
    let refs: Vec<RecordRef> = (0..60u64)
        .map(|i| {
            let off = arena.extend_from_slice(&[u8::try_from(i).unwrap(); 100]) as u64;
            RecordRef::new(i, off, 100)
        })
        .collect();
    let chunk = MemoryChunkErased::Coordinate(fgumi_sort::coordinate_chunk_from_refs(
        Arc::new(arena),
        refs,
        1,
    ));
    let mut step = SpillGather::new(1 << 20);
    step.block_size = 250; // two records per block → many ranges
    step.stage_event(SortChunkEvent::Spill { seq: 0, chunk, records_ingested_so_far: 60 });
    let mut ranges = Vec::new();
    while step.active.is_some() {
        step.plan_ranges().unwrap();
        ranges.extend(step.pending.drain(..));
    }
    assert!(ranges.len() > 1, "the chunk was planned into several ranges");
    assert!(step.held_last_block.is_held(), "the run-final block was framed and withheld");
    let mut compress = crate::sort::SpillBlockCompress::new(SpillCodec::Zstd, 1, 1 << 20);
    let last = ranges.pop().expect("ranges");
    for range in ranges {
        compress.compress_event(range).unwrap();
        assert!(pool.try_acquire().is_none(), "a range still pending holds the arena");
    }
    compress.compress_event(last).unwrap();
    assert!(pool.try_acquire().is_some(), "the last compressed range returns the arena");
}

/// Every spill file the planned path writes (`SpillGather` ranges framed and
/// compressed by `SpillBlockCompress`, written by `SpillWrite`) is
/// byte-identical to the file the serial framing writes (the same chunks cut
/// by `frame_chunk_into_blocks`, compressed, written), for both codecs and a
/// run extended across two chunks.
#[rstest::rstest]
fn spill_files_are_byte_identical_to_serial_framing(
    #[values(fgumi_sort::SpillCodec::Zstd, fgumi_sort::SpillCodec::Bgzf)]
    codec: fgumi_sort::SpillCodec,
) {
    use crate::sort::{SpillBlockCompress, SpillWrite};
    use fgumi_sort::TmpDirAllocator;

    const BLOCK: usize = 700;
    let run = |keys: std::ops::Range<u64>| {
        let recs = keys
            .map(|k| {
                (
                    RawCoordinateKey { sort_key: k },
                    vec![
                        u8::try_from(k % 251).unwrap();
                        37 + usize::try_from(k * 13 % 90).unwrap()
                    ],
                )
            })
            .collect();
        MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(recs))
    };
    // Run 0 is two contiguous chunks (one file); run 1 starts lower (a new file).
    let chunks = || vec![run(0..40), run(40..90), run(5..30)];
    let write = |events: Vec<SpillBlockEvent>| -> (tempfile::TempDir, Vec<Vec<u8>>) {
        let dir = tempfile::TempDir::new().unwrap();
        let alloc = TmpDirAllocator::new(vec![dir.path().to_path_buf()]).unwrap();
        let mut writer = SpillWrite::new(
            Arc::new(parking_lot::Mutex::new(alloc)),
            codec,
            1 << 20,
            Arc::new(Vec::new()),
        );
        let mut compress = SpillBlockCompress::new(codec, 1, 1 << 20);
        for event in events {
            writer.process_event(compress.compress_event(event).unwrap()).unwrap();
        }
        let mut names: Vec<_> =
            std::fs::read_dir(dir.path()).unwrap().map(|e| e.unwrap().path()).collect();
        names.sort();
        let files = names.iter().map(|p| std::fs::read(p).unwrap()).collect();
        (dir, files)
    };

    // Planned: through SpillGather.
    let mut step = SpillGather::new(1 << 20);
    step.block_size = BLOCK;
    let mut input: Vec<SortChunkEvent> = chunks()
        .into_iter()
        .enumerate()
        .map(|(seq, chunk)| SortChunkEvent::Spill {
            seq: u32::try_from(seq).unwrap(),
            chunk,
            records_ingested_so_far: 0,
        })
        .collect();
    input.push(SortChunkEvent::AllAnnounced {
        slot_count: 0,
        memory_chunk_count: 0,
        total_records: 0,
    });
    let planned = drive(&mut step, input);
    assert!(
        planned.iter().any(|e| matches!(e, SpillBlockEvent::Range { .. })),
        "ranges were planned"
    );
    assert!(
        planned.iter().any(
            |e| matches!(e, SpillBlockEvent::Range { framed_bytes, .. } if *framed_bytes as usize > BLOCK)
        ),
        "some range spans several frames"
    );
    let (_d1, planned_files) = write(planned);

    // Serial: the same runs framed whole, as raw blocks.
    let mut serial = Vec::new();
    let file_of_chunk = [0u32, 0, 2];
    let mut ordinal = 0;
    let all = chunks();
    for (c, chunk) in all.iter().enumerate() {
        let blocks = frame_chunk_into_blocks(chunk, BLOCK).unwrap();
        let n = blocks.len();
        for (b, bytes) in blocks.into_iter().enumerate() {
            let last_of_file =
                b + 1 == n && (c + 1 == all.len() || file_of_chunk[c + 1] != file_of_chunk[c]);
            serial.push(SpillBlockEvent::Block {
                ordinal,
                file_id: file_of_chunk[c],
                key_kind: SpillKeyKind::Coordinate,
                is_last_in_file: last_of_file,
                records_ingested_so_far: 0,
                bytes,
            });
            ordinal += 1;
        }
    }
    serial.push(SpillBlockEvent::AllAnnounced {
        ordinal,
        slot_count: 2,
        memory_chunk_count: 0,
        total_records: 0,
    });
    let (_d2, serial_files) = write(serial);

    assert_eq!(planned_files.len(), 2, "two runs, two files");
    assert!(
        planned_files == serial_files,
        "spill files differ from the serial framing ({codec:?})"
    );
}

/// Each planned range spans up to `FRAMES_PER_RANGE` frames, never the
/// run-final one, and carries the frame size `SpillBlockCompress` re-cuts it
/// at: the ranges' frames, then the withheld final block, are exactly the
/// frames the planner cuts one at a time.
#[test]
fn ranges_span_up_to_frames_per_range_frames() {
    const BLOCK: usize = 250;
    let chunk = coord_chunk((0..70u8).map(|i| vec![i; 100]).collect());
    let mut frames = vec![0usize];
    while *frames.last().unwrap() < chunk.len() {
        frames.push(plan_next_block(&chunk, *frames.last().unwrap(), BLOCK).0);
    }
    let n_frames = frames.len() - 1;
    let mut step = SpillGather::new(1 << 20);
    step.block_size = BLOCK;
    step.stage_event(SortChunkEvent::Spill { seq: 0, chunk, records_ingested_so_far: 70 });
    let mut ranges = Vec::new();
    while step.active.is_some() {
        step.plan_ranges().unwrap();
        ranges.extend(step.pending.drain(..));
    }
    let mut at = 0usize;
    let mut spans = Vec::new();
    for range in &ranges {
        let SpillBlockEvent::Range { start, end, frame_bytes, .. } = range else {
            panic!("only ranges are planned before the run closes");
        };
        assert_eq!(*frame_bytes as usize, BLOCK);
        assert_eq!(usize::try_from(*start).unwrap(), frames[at], "ranges are contiguous");
        let span = frames[at..]
            .iter()
            .position(|&f| f == usize::try_from(*end).unwrap())
            .expect("ends on a frame");
        spans.push(span);
        at += span;
    }
    assert_eq!(at, n_frames - 1, "every frame but the final one is in a range");
    assert!(step.held_last_block.is_held(), "the run-final frame is withheld");
    let full = n_frames.saturating_sub(1) / FRAMES_PER_RANGE;
    assert!(spans[..full].iter().all(|&s| s == FRAMES_PER_RANGE), "{spans:?}");
    assert!(spans.iter().all(|&s| (1..=FRAMES_PER_RANGE).contains(&s)), "{spans:?}");
}
