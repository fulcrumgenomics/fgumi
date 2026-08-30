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

/// The `(file_id, is_last_in_file)` of every `Block` event, in order. Panics if a
/// non-`Block` event is present, so callers separate spills from passthroughs.
fn block_files(events: &[SpillBlockEvent]) -> Vec<(u32, bool)> {
    events
        .iter()
        .map(|e| match e {
            SpillBlockEvent::Block { file_id, is_last_in_file, .. } => (*file_id, *is_last_in_file),
            _ => panic!("expected only Block events (a non-Block event was present)"),
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

    let ordinals: Vec<u64> = step
        .pending
        .iter()
        .map(|e| match e {
            SpillBlockEvent::Block { ordinal, .. }
            | SpillBlockEvent::Residual { ordinal, .. }
            | SpillBlockEvent::AllAnnounced { ordinal, .. } => *ordinal,
        })
        .collect();
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
    // line must read exactly like fgumi_sort::external::log_run_formation.
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
    // A residual-only sort (no spills) reports no run-formation line, matching
    // the owned engine (which only logs when it spilled).
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
    assert!(step.pending.is_empty(), "framing happens in produce_blocks, not stage_event");
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
