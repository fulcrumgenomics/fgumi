// End-to-end merge behavior is covered by the integration tests in sort/tests.rs;
// the unit tests here pin the memory-lane fail-closed guard.

use super::*;

/// A residual chunk stranded in a lane that does not match the sort order must
/// fail closed: `build_driver` and the single-chunk fast path consume only the
/// selected lane, so such a chunk would otherwise be dropped silently while
/// `total_len()` still counted it toward setup completion.
#[test]
fn mismatched_memory_lane_fails_closed() {
    let mut chunks = MemoryChunksByKind::default();
    let chunk =
        InMemoryChunk::from_owned_records(vec![(RawCoordinateKey { sort_key: 1 }, vec![9u8; 8])]);
    chunks.push(MemoryChunkErased::Coordinate(chunk)).expect("coordinate lane never mismatches");

    // The coordinate lane matches a Coordinate sort → accepted.
    chunks.ensure_single_lane(SortOrder::Coordinate).expect("matching lane is accepted");

    // The same chunk under a Queryname sort is a lane mismatch → fail closed.
    let err = chunks
        .ensure_single_lane(SortOrder::Queryname(QuerynameComparator::Natural))
        .expect_err("stray coordinate chunk under a queryname sort must error");
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
}

/// Template-coordinate spill slots present but no residual chunk to identify the
/// `--key-types` lane must fail closed: defaulting to K40 would mis-decode narrow
/// (K24/Cb32/Tert32) spill files. Unreachable for valid input (Phase-1 always
/// emits a variant-tagged residual), so this guards against a seal-logic
/// regression. With no slots, the empty-input case is still accepted.
#[test]
fn empty_template_lane_with_spill_slots_fails_closed() {
    let slot = Arc::new(SortMergeSlot::for_test(0, fgumi_sort::SpillCodec::Bgzf));
    // `Box<dyn MergeDriverDyn>` isn't `Debug`, so match rather than `expect_err`.
    match build_driver(SortOrder::TemplateCoordinate, vec![slot], MemoryChunksByKind::default(), 1)
    {
        Err(e) => assert_eq!(e.kind(), std::io::ErrorKind::InvalidData),
        Ok(_) => panic!("empty template lane with spill slots must fail closed"),
    }

    // No slots → empty input; any key width is safe (nothing to merge).
    assert!(
        build_driver(SortOrder::TemplateCoordinate, Vec::new(), MemoryChunksByKind::default(), 0)
            .is_ok(),
        "empty template lane with no slots is valid",
    );
}

/// The `--key-types` narrowed-lane variant is chosen once per sort and is global
/// to the run. A template chunk arriving with a different variant means phase 1
/// and the merge disagree about the key width; merging on would compare keys of
/// different layouts and emit silently mis-ordered output. It must fail closed,
/// like the sibling `ensure_single_lane` / `build_driver` violations — not panic.
#[test]
fn template_variant_change_mid_sort_fails_closed() {
    use fgumi_sort::{TemplateKey24, TemplateMemChunk, TertKey32};

    let mut chunks = MemoryChunksByKind::default();

    let k24 = InMemoryChunk::from_owned_records(vec![(TemplateKey24::default(), vec![1u8; 8])]);
    chunks
        .push(MemoryChunkErased::TemplateCoordinate(TemplateMemChunk::K24(k24)))
        .expect("the first template chunk establishes the variant");

    // A second chunk in a different lane width is the protocol violation.
    let tert = InMemoryChunk::from_owned_records(vec![(TertKey32::default(), vec![2u8; 8])]);
    let err = chunks
        .push(MemoryChunkErased::TemplateCoordinate(TemplateMemChunk::Tert32(tert)))
        .expect_err("a variant change must be rejected");

    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    let msg = err.to_string();
    assert!(msg.contains("variant changed mid-sort"), "unexpected message: {msg}");
    // Both variants are named so the failure is diagnosable from the log alone.
    assert!(msg.contains("K24"), "error names the accumulated variant: {msg}");
    assert!(msg.contains("Tert32"), "error names the offending variant: {msg}");
}

/// Repeated chunks of the SAME variant are the normal path and must keep working.
#[test]
fn repeated_template_chunks_of_one_variant_accumulate() {
    use fgumi_sort::{TemplateKey24, TemplateMemChunk};

    let mut chunks = MemoryChunksByKind::default();
    for i in 0..3u8 {
        let c = InMemoryChunk::from_owned_records(vec![(TemplateKey24::default(), vec![i; 8])]);
        chunks
            .push(MemoryChunkErased::TemplateCoordinate(TemplateMemChunk::K24(c)))
            .expect("same-variant chunks accumulate");
    }
    assert_eq!(chunks.total_len(), 3, "all three chunks are retained");
}

/// A dispatch that popped input reports `Progress` whatever the emit after it
/// reports: the pop freed a slot on the input edge, and the plan's reverse wake
/// to an upstream producer holding an item for it runs only on `Progress`. Here
/// the dispatch absorbs a spill slot's announcement and the final
/// `AllAnnounced`, completes the setup, and the first merge pass stalls at once
/// (the slot has no decompressed block yet), which the emit reports as
/// `Contention`; the dispatch must still be `Progress`.
#[test]
fn a_dispatch_that_absorbed_setup_input_is_progress() {
    use fgumi_pipeline_core::testing::StepProbe;
    let mut step = SortMerge::<RecordBatchOutput>::new(SortOrder::Coordinate, 1 << 20);
    let probe = StepProbe::new(&step);
    let slot = Arc::new(SortMergeSlot::for_test(0, fgumi_sort::SpillCodec::Bgzf));
    probe.push_input(SortPhase2Event::SpillReady {
        slot,
        path: std::path::PathBuf::from("spill-0"),
        records_ingested_so_far: 1,
    });
    probe.push_input(SortPhase2Event::AllAnnounced {
        slot_count: 1,
        memory_chunk_count: 0,
        total_records: 1,
    });
    let outcome = probe.try_run(&mut step).expect("try_run");
    assert!(probe.input_is_empty(), "both announcements were absorbed");
    assert_eq!(outcome, StepOutcome::Progress, "a dispatch that took input is not idle");
}

/// A one-record spill block at `pos` on reference 0 (the coordinate key is
/// embedded in the record, so the frame is `[u32 len][record]`).
fn one_record_block(pos: usize) -> Vec<u8> {
    let pos = i32::try_from(pos).expect("small test positions");
    let record = fgumi_raw_bam::testutil::make_bam_bytes(
        0,
        pos,
        0,
        format!("r{pos:04}").as_bytes(),
        &[],
        4,
        -1,
        -1,
        &[],
    );
    let mut block = Vec::new();
    fgumi_sort::frame_keyed_record_into(&mut block, &RawCoordinateKey::default(), &record)
        .expect("frame record");
    block
}

/// A stall whose partial batch is refused by a full output edge leaves the
/// merge blocked on its output, not on the spill supply: it lowers `starved`
/// (so the pool's refill walk stops putting the supply ahead of the steps that
/// drain that edge), withdraws its registration, and the dispatch that then
/// parks on the held batch is not a parking registration.
#[test]
fn a_partial_batch_held_on_a_full_output_lowers_starved() {
    use fgumi_pipeline_core::testing::StepProbe;
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    // A 256-byte edge: a few one-record batches (their buffers are sized up to
    // the 256-byte cap) fill it, and one record never fills a batch, so each
    // stall flushes a partial batch.
    let mut step =
        SortMerge::<RecordBatchOutput>::with_target_batch_count(SortOrder::Coordinate, 256, 1024)
            .with_merge_demand(Arc::clone(&demand));
    let probe = StepProbe::new(&step);
    let slots: Vec<Arc<SortMergeSlot>> = (0..2u32)
        .map(|file_id| {
            let slot = Arc::new(SortMergeSlot::for_test(file_id, fgumi_sort::SpillCodec::Bgzf));
            slot.push_decompressed_for_test(one_record_block(file_id as usize));
            slot
        })
        .collect();
    for slot in &slots {
        probe.push_input(SortPhase2Event::SpillReady {
            slot: Arc::clone(slot),
            path: std::path::PathBuf::from("spill"),
            records_ingested_so_far: 1000,
        });
    }
    probe.push_input(SortPhase2Event::AllAnnounced {
        slot_count: 2,
        memory_chunk_count: 0,
        total_records: 1000,
    });
    let starved = demand.starved_signal();
    // Each dispatch merges one record and runs dry on that record's slot,
    // flushing it as a partial batch; feeding the slot's next record first
    // keeps every dispatch a stall with one record merged. Nothing drains the
    // edge, so a partial batch is eventually refused and held.
    let mut next = 2;
    loop {
        assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Progress);
        if step.held.is_held() {
            break;
        }
        assert!(next < 1000, "the output edge never filled");
        slots[next % 2].push_decompressed_for_test(one_record_block(next));
        next += 1;
    }
    assert!(
        !starved.load(std::sync::atomic::Ordering::Relaxed),
        "a merge blocked on its output must not keep the refill walk raised"
    );
    assert_eq!(demand.awaited(), None, "the stall's registration was withdrawn");
    // The edge is still full: this dispatch parks on the held batch.
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Contention);
    assert_eq!(demand.snapshot().parking_registrations, 0, "it parks on output, not a slot");
}

/// A stall books the awaited slot's state: a slot with a read in progress is
/// `issued`, and one with nothing read, stashed or in flight is `starved`. A
/// slot whose read blocks wait in the stash is not a stall: the merge serves
/// itself (`a_stashed_slot_is_served_not_parked_on`).
#[rstest::rstest]
#[case::starved(0, (1, 0, 0))]
#[case::issued(4096, (0, 1, 0))]
fn a_stall_books_the_awaited_slots_state(#[case] issued: u64, #[case] want: (u64, u64, u64)) {
    use fgumi_pipeline_core::testing::StepProbe;
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let mut step = SortMerge::<RecordBatchOutput>::with_target_batch_count(
        SortOrder::Coordinate,
        1 << 20,
        1024,
    )
    .with_merge_demand(Arc::clone(&demand));
    let probe = StepProbe::new(&step);
    // Slot 0 holds a record; slot 1 has nothing decompressed, so priming
    // stalls on it.
    let ready = Arc::new(SortMergeSlot::for_test(0, fgumi_sort::SpillCodec::Bgzf));
    ready.push_decompressed_for_test(one_record_block(0));
    let awaited = Arc::new(SortMergeSlot::for_test(1, fgumi_sort::SpillCodec::Bgzf));
    awaited.bp_note_issued(issued);
    for slot in [&ready, &awaited] {
        probe.push_input(SortPhase2Event::SpillReady {
            slot: Arc::clone(slot),
            path: std::path::PathBuf::from("spill"),
            records_ingested_so_far: 2,
        });
    }
    probe.push_input(SortPhase2Event::AllAnnounced {
        slot_count: 2,
        memory_chunk_count: 0,
        total_records: 2,
    });
    let _ = probe.try_run(&mut step).unwrap();
    let s = demand.snapshot();
    assert_eq!(s.stall_episodes, 1, "{s:?}");
    assert_eq!((s.awaited_starved, s.awaited_issued, s.awaited_decompressing), want);
}

/// A stall on a slot whose read blocks wait in the stash is served by the
/// merge itself, not booked as a stall at all: the only stall episode is
/// the later one on an emptied slot (`starved`), and the self-served
/// registration is not a parking one.
#[test]
fn a_stashed_slot_is_served_not_parked_on() {
    use fgumi_pipeline_core::testing::StepProbe;
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    let mut step = SortMerge::<RecordBatchOutput>::with_target_batch_count(
        SortOrder::Coordinate,
        1 << 20,
        1024,
    )
    .with_merge_demand(Arc::clone(&demand));
    let probe = StepProbe::new(&step);
    let ready = Arc::new(SortMergeSlot::for_test(0, fgumi_sort::SpillCodec::Bgzf));
    ready.push_decompressed_for_test(one_record_block(0));
    let stashed = Arc::new(SortMergeSlot::for_test(1, fgumi_sort::SpillCodec::Bgzf));
    let frame = fgumi_sort::SpillBlockCompressor::new(fgumi_sort::SpillCodec::Bgzf, 1)
        .unwrap()
        .compress_block(&one_record_block(1))
        .unwrap();
    stashed.bp_stash_frames_for_test(vec![frame], false);
    for slot in [&ready, &stashed] {
        probe.push_input(SortPhase2Event::SpillReady {
            slot: Arc::clone(slot),
            path: std::path::PathBuf::from("spill"),
            records_ingested_so_far: 2,
        });
    }
    probe.push_input(SortPhase2Event::AllAnnounced {
        slot_count: 2,
        memory_chunk_count: 0,
        total_records: 2,
    });
    let _ = probe.try_run(&mut step).unwrap();
    let s = demand.snapshot();
    assert_eq!((s.self_served_episodes, s.self_served_blocks), (1, 1), "{s:?}");
    assert_eq!(s.stall_episodes, s.awaited_starved, "a served stash is not a stall: {s:?}");
    assert!(
        s.parking_registrations < s.registrations,
        "the self-served registration is not a parking one: {s:?}"
    );
}
