//! Unit tests for `SortAndSpill::finalize_into_pending` (S2-003).
//!
//! The end-to-end three-step chain tests live in `sort/tests.rs`; those exercise
//! the happy path (in-memory / multi-spill matches the legacy sorter) but never
//! assert the *values* inside the terminal `AllAnnounced` sentinel — i.e. the
//! `slot_count` / `memory_chunk_count` derivation and the empty-chunk-skip logic
//! the #449 overflow guard protects. This module covers `finalize_into_pending`
//! directly:
//!
//! - the announced `slot_count` equals the number of emitted `SpillReady` events;
//! - the announced `memory_chunk_count` equals the number of emitted
//!   (non-empty) `MemoryChunk` events;
//! - an all-in-memory stream announces zero slots;
//! - a stream with nothing buffered finalizes to `None` (no events at all).

use super::*;

use fgumi_raw_bam::RawRecord;
use fgumi_raw_bam::testutil::make_bam_bytes;
use fgumi_sort::{RawExternalSorter, SortOrder};
use noodles::sam::Header;

/// Build `n` coordinate-sortable records, each carrying a `payload_len`-base
/// (and `payload_len`-quality) body so the in-memory footprint per record is
/// large enough that a small `memory_limit` forces spills.
fn sized_records(n: usize, payload_len: usize) -> Vec<RawRecord> {
    (0..n)
        .map(|i| {
            let name = format!("rec_{i:06}");
            // tid 0, ascending positions; coordinate sort handles these directly.
            let pos = i32::try_from(i).expect("pos fits i32");
            let bytes = make_bam_bytes(0, pos, 0, name.as_bytes(), &[], payload_len, -1, -1, &[]);
            RawRecord::from(bytes)
        })
        .collect()
}

/// Drive a coordinate `SortStream` to completion via `finalize_into_pending` and
/// return `(spill_ready_events, non_empty_memory_chunk_events, announced)` where
/// `announced` is the trailing `AllAnnounced { slot_count, memory_chunk_count }`,
/// or `None` if `finalize_into_pending` produced no events.
fn finalize_counts(memory_limit: usize, records: &[RawRecord]) -> Option<(usize, usize, u32, u32)> {
    let header = Header::default();
    let sorter = RawExternalSorter::new(SortOrder::Coordinate)
        .memory_limit(memory_limit)
        .threads(2)
        .output_compression(1)
        .temp_compression(1);

    let mut stream = build_stream(sorter, &header).expect("build coordinate stream");
    stream.push_records(records.iter().map(RawRecord::as_ref)).expect("push records into stream");

    let pending = SortAndSpill::finalize_into_pending(stream, "test")
        .expect("finalize_into_pending succeeds")?;
    let (events, _temp_dirs) = pending;

    let mut spill_ready = 0usize;
    let mut memory_chunks = 0usize;
    let mut announced: Option<(u32, u32)> = None;
    for event in &events {
        match event {
            SortPhase1Event::SpillReady { .. } => spill_ready += 1,
            SortPhase1Event::MemoryChunk { chunk, .. } => {
                // Empty chunks must never be emitted (they are skipped during
                // finalize); assert that invariant directly.
                assert!(!chunk.is_empty(), "an empty MemoryChunk leaked into the event stream");
                memory_chunks += 1;
            }
            SortPhase1Event::AllAnnounced { slot_count, memory_chunk_count, .. } => {
                assert!(announced.is_none(), "more than one AllAnnounced emitted");
                announced = Some((*slot_count, *memory_chunk_count));
            }
        }
    }

    // AllAnnounced must be the LAST event of the stream.
    assert!(
        matches!(events.back(), Some(SortPhase1Event::AllAnnounced { .. })),
        "AllAnnounced must be the final event"
    );

    let (slot_count, memory_chunk_count) = announced.expect("AllAnnounced present");
    Some((spill_ready, memory_chunks, slot_count, memory_chunk_count))
}

#[test]
fn announced_counts_match_emitted_events_with_spills() {
    // A small memory limit against many sizeable records forces at least one
    // spill AND leaves an in-memory residual chunk, so both announced counts are
    // exercised against actual emitted events.
    let records = sized_records(4_000, 200);
    let (spill_ready, memory_chunks, slot_count, memory_chunk_count) =
        finalize_counts(64 * 1024, &records).expect("non-empty finalize");

    assert!(slot_count >= 1, "small memory limit should have forced at least one spill slot");
    assert_eq!(
        usize::try_from(slot_count).unwrap(),
        spill_ready,
        "announced slot_count must equal the number of SpillReady events"
    );
    assert_eq!(
        usize::try_from(memory_chunk_count).unwrap(),
        memory_chunks,
        "announced memory_chunk_count must equal the number of (non-empty) MemoryChunk events"
    );
}

#[test]
fn all_in_memory_announces_zero_slots() {
    // A generous memory limit keeps everything in memory: zero spill slots,
    // and the residual is announced purely as memory chunks.
    let records = sized_records(500, 50);
    let (spill_ready, memory_chunks, slot_count, memory_chunk_count) =
        finalize_counts(256 * 1024 * 1024, &records).expect("non-empty finalize");

    assert_eq!(slot_count, 0, "no spills expected under a large memory limit");
    assert_eq!(spill_ready, 0, "no SpillReady events expected");
    assert!(memory_chunk_count >= 1, "the in-memory residual must be announced as a chunk");
    assert_eq!(usize::try_from(memory_chunk_count).unwrap(), memory_chunks);
}

#[test]
fn empty_stream_finalizes_to_none() {
    // No records pushed -> nothing to announce -> finalize_into_pending returns
    // None (no events, not an AllAnnounced with zeroed counts).
    let header = Header::default();
    let sorter = RawExternalSorter::new(SortOrder::Coordinate)
        .memory_limit(256 * 1024 * 1024)
        .threads(2)
        .output_compression(1)
        .temp_compression(1);
    let stream = build_stream(sorter, &header).expect("build coordinate stream");
    let pending =
        SortAndSpill::finalize_into_pending(stream, "test").expect("finalize_into_pending");
    assert!(pending.is_none(), "an empty stream must finalize to None (no events emitted)");
}
