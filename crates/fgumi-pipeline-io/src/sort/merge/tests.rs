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
    match build_driver(
        SortOrder::TemplateCoordinate,
        vec![slot],
        MemoryChunksByKind::default(),
        1,
        None,
        None,
    ) {
        Err(e) => assert_eq!(e.kind(), std::io::ErrorKind::InvalidData),
        Ok(_) => panic!("empty template lane with spill slots must fail closed"),
    }

    // No slots → empty input; any key width is safe (nothing to merge).
    assert!(
        build_driver(
            SortOrder::TemplateCoordinate,
            Vec::new(),
            MemoryChunksByKind::default(),
            0,
            None,
            None
        )
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

/// A full batch refused by a full output edge leaves the merge blocked on its
/// output, not on the spill supply: `starved` is low (the record that closed
/// the batch ended the stall episode, so the pool's refill walk does not put
/// the supply ahead of the steps that drain that edge), the last stall's
/// registration was withdrawn, and the dispatch that then parks on the held
/// batch is not a parking registration.
#[test]
fn a_batch_held_on_a_full_output_leaves_starved_low() {
    use fgumi_pipeline_core::testing::StepProbe;
    let demand = Arc::new(fgumi_sort::MergeDemand::new());
    // A 256-byte edge and one-record batches: every produced record closes a
    // full batch (its buffer sized up to the 256-byte cap), so a few of them
    // fill the edge, and the `Batch` arm's refused push is what holds one.
    let mut step =
        SortMerge::<RecordBatchOutput>::with_target_batch_count(SortOrder::Coordinate, 256, 1)
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
    // Each dispatch merges one record, delivers it as a full batch and runs
    // dry on that record's slot; feeding the slot's next record first keeps
    // every dispatch one record then a stall. Nothing drains the edge, so a
    // full batch is eventually refused and held.
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

/// The merge of `fgumi merge`'s inputs checks, before it accepts `Done`, that
/// every input was read to its clean end: a merge whose driver retired a
/// source early (simulated here by a driver over only the first of the two
/// inputs) is an error naming the unread input, not a short clean exit. With
/// every input drained, the same merge completes.
#[rstest::rstest]
#[case::an_unread_input_fails(false)]
#[case::every_input_drained_completes(true)]
fn merge_inputs_done_requires_every_input_read(#[case] second_drained: bool) {
    use fgumi_pipeline_core::testing::StepProbe;
    let mode = MergeInputs {
        inputs: vec!["a.bam".to_string(), "b.bam".to_string()],
        order_name: "coordinate",
        template_ctx: None,
        rewrite: None,
    };
    let mut step = SortMerge::<RecordBatchOutput>::with_target_batch_count(
        SortOrder::Coordinate,
        1 << 20,
        1024,
    )
    .with_merge_inputs(mode.clone());
    let probe = StepProbe::new(&step);
    let slots: Vec<Arc<SortMergeSlot>> = (0..2u32)
        .map(|file_id| Arc::new(SortMergeSlot::for_test(file_id, fgumi_sort::SpillCodec::Bgzf)))
        .collect();
    for slot in &slots {
        probe.push_input(SortPhase2Event::SpillReady {
            slot: Arc::clone(slot),
            path: std::path::PathBuf::from("input"),
            records_ingested_so_far: 0,
        });
    }
    probe.push_input(SortPhase2Event::AllAnnounced {
        slot_count: 2,
        memory_chunk_count: 0,
        total_records: 0,
    });
    // Nothing is decompressed yet: the merge sets up and stalls.
    let _ = probe.try_run(&mut step).unwrap();
    assert!(matches!(step.state, SortMergeState::Merging { .. }), "the merge set up");
    assert_eq!(step.merge_input_slots.len(), 2, "both inputs are kept for the check");
    // Swap in a driver that has lost the second input, then finish the first.
    let short = build_merge_inputs_driver(
        SortOrder::Coordinate,
        vec![Arc::clone(&slots[0])],
        None,
        None,
        &mode,
    )
    .expect("driver");
    let SortMergeState::Merging { driver, .. } = &mut step.state else { unreachable!() };
    *driver = short;
    slots[0].push_decompressed_for_test(one_record_block(0));
    slots[0].set_queue_eof_for_test();
    if second_drained {
        slots[1].set_queue_eof_for_test();
    }
    let mut result = Ok(StepOutcome::Progress);
    for _ in 0..8 {
        result = probe.try_run(&mut step);
        if result.is_err() || matches!(step.state, SortMergeState::Done) {
            break;
        }
    }
    if second_drained {
        result.expect("every input drained: the merge completes");
        assert!(matches!(step.state, SortMergeState::Done));
    } else {
        let err = result.expect_err("an unread input must fail the merge");
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("b.bam"), "names the unread input: {err}");
    }
}

/// `ensure_merge_inputs_read` rejects a failed input and one with blocks left
/// or no EOF, and accepts inputs that each reached a clean, empty EOF.
#[test]
fn ensure_merge_inputs_read_needs_a_clean_end_on_every_input() {
    let names = vec!["a.bam".to_string()];
    let slot = || Arc::new(SortMergeSlot::for_test(0, fgumi_sort::SpillCodec::Bgzf));
    let not_at_eof = slot();
    let err = ensure_merge_inputs_read(&[not_at_eof], &names).expect_err("an input not at EOF");
    assert!(err.to_string().contains("lost records"), "{err}");
    assert!(err.to_string().contains("a.bam"), "{err}");
    let left_over = slot();
    left_over.push_decompressed_for_test(one_record_block(0));
    left_over.set_queue_eof_for_test();
    let err = ensure_merge_inputs_read(&[left_over], &names).expect_err("blocks left over");
    assert!(err.to_string().contains("lost records"), "{err}");
    assert!(err.to_string().contains("a.bam"), "{err}");
    let failed = slot();
    failed.mark_failed();
    let err = ensure_merge_inputs_read(&[failed], &names).expect_err("a failed input");
    assert!(err.to_string().contains("failed reading input a.bam"), "{err}");
    let clean = slot();
    clean.set_queue_eof_for_test();
    ensure_merge_inputs_read(&[clean], &names).expect("a clean end");
}

/// A single in-memory chunk of `n` coordinate records, one per position.
fn gather_test_chunk(n: usize) -> MemoryChunkErased {
    MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(
        (0..n)
            .map(|i| {
                let pos = i32::try_from(i).expect("small test positions");
                let record = fgumi_raw_bam::testutil::make_bam_bytes(
                    0,
                    pos,
                    0,
                    format!("r{i:05}").as_bytes(),
                    &[],
                    4,
                    -1,
                    -1,
                    &[],
                );
                (RawCoordinateKey { sort_key: i as u64 }, record)
            })
            .collect(),
    ))
}

/// Run `posted` gather jobs and deliver their results to `inbox` out of order:
/// the awaited (lowest) ordinal last, with a dispatch in between that must
/// park and emit nothing ahead of it.
fn deliver_out_of_order(
    probe: &fgumi_pipeline_core::testing::StepProbe<SortMerge<BlockOutput>>,
    step: &mut SortMerge<BlockOutput>,
    inbox: &fgumi_pipeline_core::OrderedInbox<crate::sort::merge_job::MergeResult<BlockOutput>>,
    blocks: &fgumi_pipeline_core::BranchInputHandle<DecompressedBlock>,
    posted: Vec<crate::sort::merge_job::MergeTicket<BlockOutput>>,
) {
    use crate::sort::merge_job::MERGE_SLICE;
    use fgumi_pipeline_core::{InputHandle, Job, JobStep};

    let mut done: Vec<_> = posted
        .into_iter()
        .map(|mut job| {
            let ordinal = job.ordinal();
            match job.run_slice(MERGE_SLICE).expect("gather") {
                JobStep::Done(out) => (ordinal, out),
                JobStep::Continue | JobStep::Spent => panic!("a gather is one slice"),
            }
        })
        .collect();
    done.reverse();
    let (awaited, last) = done.pop().expect("one result at least");
    if !done.is_empty() {
        for (ordinal, out) in done {
            let _ = inbox.deliver(ordinal, out);
        }
        assert_eq!(
            probe.try_run(step).expect("try_run"),
            StepOutcome::NoProgress,
            "the next ordinal is still missing"
        );
        assert!(blocks.pop().is_none(), "nothing leaves ahead of the awaited ordinal");
    }
    let _ = inbox.deliver(awaited, last);
}

/// Gather mode posts one `MergeJob::Gather` ticket per planned block on branch 1, with
/// ordinals `0..k` in order and never more than `window` (`fanout_width * 4`)
/// ahead of the results taken back; it parks (`NoProgress`) while the next
/// ordinal is missing, and results delivered out of order leave branch 0 in
/// ordinal order, byte-identical to gathering `plan_fast_path_blocks`' plan
/// serially.
#[test]
fn gather_mode_posts_jobs_within_window() {
    use fgumi_pipeline_core::testing::StepProbe;
    use fgumi_pipeline_core::{InputHandle, Job, OrderedInbox};

    use crate::sort::merge_job::{MergeResult, MergeTicket};

    const N: usize = 2000;
    const TARGET: usize = 64;
    const LIMIT: u64 = 1 << 20;
    const WIDTH: usize = 2;
    let window = (WIDTH * FAST_PATH_BLOCKS_PER_WORKER_WINDOW) as u64;

    let expected: Vec<DecompressedBlock> = {
        let chunk = gather_test_chunk(N);
        let plan = plan_fast_path_blocks::<BlockOutput>(&chunk, TARGET, LIMIT);
        plan.iter()
            .enumerate()
            .map(|(i, &b)| {
                gather_fast_block::<BlockOutput>(
                    &chunk,
                    b,
                    i as u64,
                    INITIAL_OUTPUT_BUFFER_BYTES,
                    TARGET,
                )
                .unwrap()
            })
            .collect()
    };
    let k = expected.len() as u64;
    assert!(k > 2 * window, "fixture must span several windows ({k} blocks)");

    let inbox = OrderedInbox::<MergeResult<BlockOutput>>::new(true);
    let mut step =
        SortMerge::<BlockOutput>::with_target_batch_count(SortOrder::Coordinate, LIMIT, TARGET)
            .with_inbox(Arc::clone(&inbox))
            .with_fanout_width(WIDTH)
            .with_fast_path_min_records(0);
    let mut probe = StepProbe::new(&step);
    let jobs = probe.take_output::<MergeTicket<BlockOutput>>(1);
    let blocks = probe.take_output::<DecompressedBlock>(0);
    probe.push_input(SortPhase2Event::MemoryChunk {
        chunk: Arc::new(gather_test_chunk(N)),
        records_ingested_so_far: N as u64,
    });
    probe.push_input(SortPhase2Event::AllAnnounced {
        slot_count: 0,
        memory_chunk_count: 1,
        total_records: N as u64,
    });
    probe.close_input();

    let mut posted: Vec<MergeTicket<BlockOutput>> = Vec::new();
    let mut next_ordinal = 0u64;
    let mut emitted: Vec<DecompressedBlock> = Vec::new();
    let mut finished = false;
    for _ in 0..10_000 {
        let outcome = probe.try_run(&mut step).expect("try_run");
        while let Some(job) = jobs.pop() {
            assert_eq!(job.ordinal(), next_ordinal, "jobs are posted in ordinal order");
            next_ordinal += 1;
            posted.push(job);
        }
        while let Some(block) = blocks.pop() {
            emitted.push(block);
        }
        let taken = emitted.len() as u64;
        assert!(
            next_ordinal <= taken + window,
            "posted {next_ordinal} with {taken} taken exceeds the window {window}"
        );
        if outcome == StepOutcome::Finished {
            finished = true;
            break;
        }
        if outcome == StepOutcome::NoProgress {
            assert!(!posted.is_empty(), "it parks only with jobs in flight");
            deliver_out_of_order(&probe, &mut step, &inbox, &blocks, std::mem::take(&mut posted));
        }
    }
    assert!(finished, "gather mode finished");
    assert_eq!(next_ordinal, k, "one job per planned block");
    assert!(jobs.is_empty() && blocks.is_empty());
    assert_eq!(emitted.len(), expected.len());
    for (i, (got, want)) in emitted.iter().zip(&expected).enumerate() {
        assert_eq!(got.batch_serial, i as u64, "block {i} ordinal");
        assert_eq!(got.bytes, want.bytes, "block {i} bytes");
    }
    assert_eq!(inbox.in_flight_bytes(), 0, "every result taken");
    assert_eq!(probe.counter(RECORDS), N as u64);
}

/// A gather-mode merge over `N` records at fan-out width `WIDTH`, its input
/// queued and closed, with the job branch's byte bound optionally shrunk.
/// Returns the merge, its probe and the job branch.
fn gather_mode_fixture(
    job_edge_limit: Option<u64>,
) -> (
    SortMerge<BlockOutput>,
    fgumi_pipeline_core::testing::StepProbe<SortMerge<BlockOutput>>,
    fgumi_pipeline_core::BranchInputHandle<crate::sort::merge_job::MergeTicket<BlockOutput>>,
) {
    use fgumi_pipeline_core::OrderedInbox;
    use fgumi_pipeline_core::testing::StepProbe;

    use crate::sort::merge_job::{MergeResult, MergeTicket};

    const N: usize = 2000;
    const TARGET: usize = 64;
    const LIMIT: u64 = 1 << 20;
    const WIDTH: usize = 2;

    let inbox = OrderedInbox::<MergeResult<BlockOutput>>::new(true);
    let mut step =
        SortMerge::<BlockOutput>::with_target_batch_count(SortOrder::Coordinate, LIMIT, TARGET)
            .with_inbox(inbox)
            .with_fanout_width(WIDTH)
            .with_fast_path_min_records(0);
    if let Some(limit) = job_edge_limit {
        step = step.with_job_edge_limit(limit);
    }
    let mut probe = StepProbe::new(&step);
    let jobs = probe.take_output::<MergeTicket<BlockOutput>>(1);
    probe.push_input(SortPhase2Event::MemoryChunk {
        chunk: Arc::new(gather_test_chunk(N)),
        records_ingested_so_far: N as u64,
    });
    probe.push_input(SortPhase2Event::AllAnnounced {
        slot_count: 0,
        memory_chunk_count: 1,
        total_records: N as u64,
    });
    probe.close_input();
    (step, probe, jobs)
}

/// Dispatch `step` until it parks (`NoProgress`), failing on any other
/// outcome or if it never parks.
fn run_until_parked(
    probe: &fgumi_pipeline_core::testing::StepProbe<SortMerge<BlockOutput>>,
    step: &mut SortMerge<BlockOutput>,
) {
    for _ in 0..1000 {
        match probe.try_run(step).expect("try_run") {
            StepOutcome::NoProgress => return,
            StepOutcome::Progress => {}
            other => panic!("unexpected outcome {other:?} before parking"),
        }
    }
    panic!("the merge never parked");
}

/// With nothing popping the job branch, gather mode posts exactly `window`
/// jobs (the branch, sized for twice the window, never refuses one) and then
/// parks waiting on ordinal 0.
#[test]
fn gather_mode_posts_full_window_unrefused() {
    use fgumi_pipeline_core::{InputHandle, Job};

    let (mut step, probe, jobs) = gather_mode_fixture(None);
    let window = step.gather_window() as u64;
    run_until_parked(&probe, &mut step);
    assert!(!step.held_job.is_held(), "no post within the window is refused");
    let ordinals: Vec<u64> = std::iter::from_fn(|| jobs.pop()).map(|j| j.ordinal()).collect();
    assert_eq!(ordinals, (0..window).collect::<Vec<_>>(), "exactly the window is queued");
    assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::NoProgress);
}

/// A post the job branch refuses is held and re-posted before any new job,
/// so jobs leave in ordinal order.
#[test]
fn gather_mode_reposts_refused_job_first() {
    use fgumi_pipeline_core::{InputHandle, Job};

    // Room for one job descriptor: the second post is refused.
    let (mut step, probe, jobs) = gather_mode_fixture(Some(1));
    run_until_parked(&probe, &mut step);
    assert!(step.held_job.is_held(), "the second post was refused and held");
    assert_eq!(jobs.pop().map(|j| j.ordinal()), Some(0));
    assert!(jobs.pop().is_none(), "the refused job is not on the branch");
    for expected in 1..4u64 {
        // The edge has room again: the held job goes first, and the next new
        // post is refused and held in its place.
        assert_eq!(probe.try_run(&mut step).expect("try_run"), StepOutcome::Progress);
        assert_eq!(jobs.pop().map(|j| j.ordinal()), Some(expected), "held job re-posted first");
        assert!(jobs.pop().is_none());
        assert!(step.held_job.is_held(), "ordinal {} is held", expected + 1);
    }
}

/// The job branch keeps a fixed limit of twice the gather window of
/// descriptors, outside the queue-memory budget; the output branch is a
/// budgeted byte-bounded edge.
#[rstest::rstest]
#[case::width_one(1)]
#[case::width_four(4)]
fn job_branch_is_a_fixed_byte_bound(#[case] width: usize) {
    let merge = SortMerge::<RecordBatchOutput>::new(SortOrder::Coordinate, 1 << 20)
        .with_fanout_width(width);
    let queues = merge.profile().output_queues;
    let per_job = std::mem::size_of::<crate::sort::merge_job::MergeTicket<RecordBatchOutput>>();
    let expected_limit = (width * FAST_PATH_BLOCKS_PER_WORKER_WINDOW * 2 * per_job) as u64;
    assert_eq!(merge.job_edge_limit(), expected_limit);
    assert_eq!(
        queues,
        vec![
            QueueSpec::ByteBounded { limit_bytes: 1 << 20 },
            QueueSpec::FixedByteBounded { limit_bytes: expected_limit },
        ]
    );
    assert!(queues[0].is_budgeted());
    assert!(!queues[1].is_budgeted());
}

// ── Range mode ──────────────────────────────────────────────────────────────

use crate::sort::merge_job::{MERGE_SLICE, MergeTicket};
use crate::sort::merge_mode::{MergeMode, SerialReason};
use crate::sort::supply_ledger::SpillSupply;
use fgumi_pipeline_core::builder::{Pipeline, PipelineConfig};
use fgumi_pipeline_core::{Ordered, PipelineError};
use fgumi_sort::range_merge_test_support::{serial_reference, write_indexed_run};
use fgumi_sort::{
    CbKey32, MAX_RANGES, RangeSizing, RawSortKey, SpillKeyKind, TemplateKey, TertKey32,
};

/// Output batch record cap of the range-mode tests.
const RANGE_TARGET: usize = 64;
/// Output batch byte cap (and output edge bound) of the range-mode tests:
/// small enough that some batches close on bytes, not records.
const RANGE_LIMIT: u64 = 8 * 1024;
/// Spill frame size of the range-mode runs: many frames, so many splitters.
const RANGE_FRAME_BYTES: usize = 2048;
/// Plan sizing of the range-mode tests: ranges of ~16 KiB, so a few hundred
/// KiB of input cuts into a few dozen ranges.
const RANGE_SIZING: RangeSizing =
    RangeSizing { target_range_bytes: 16 * 1024, max_ranges: MAX_RANGES };
/// The range-mode tests' batch caps: `(target_batch_count, output_byte_limit)`.
const RANGE_BATCH: (usize, u64) = (RANGE_TARGET, RANGE_LIMIT);
/// Ranges of ~1 KiB (about one spill frame each, a dozen records)...
const TINY_RANGE_SIZING: RangeSizing =
    RangeSizing { target_range_bytes: 1024, max_ranges: MAX_RANGES };
/// ...under batches of 2000 records with no byte cap: well over a window of
/// consecutive ranges fits in one batch, so the cutter's open batch pins the
/// window's floor while nothing is in flight.
const TINY_RANGE_BATCH: (usize, u64) = (2000, 1 << 30);

/// A merge key the range-mode tests run over.
trait Lane: RawSortKey + Default + Copy + Send + Sync + std::fmt::Debug + 'static {
    /// The spill kind whose runs this key merges.
    const KIND: SpillKeyKind;
    /// The sort order that merges this key.
    const ORDER: SortOrder;
    /// A record with key ordinal `ord` (equal ordinals give equal keys) and a
    /// body unique to `id`.
    fn record(ord: u64, id: u64) -> (Self, Vec<u8>);
    /// `records` (sorted) as a residual chunk.
    fn chunk(records: Vec<(Self, Vec<u8>)>) -> MemoryChunkErased;
}

/// A BAM record at `pos` on reference 0, named and sized by `id`.
fn range_body(pos: i32, id: u64) -> Vec<u8> {
    fgumi_raw_bam::testutil::make_bam_bytes(
        0,
        pos,
        0,
        format!("r{id:06}").as_bytes(),
        &[],
        8 + usize::try_from(id % 64).unwrap(),
        -1,
        -1,
        &[],
    )
}

impl Lane for RawCoordinateKey {
    const KIND: SpillKeyKind = SpillKeyKind::Coordinate;
    const ORDER: SortOrder = SortOrder::Coordinate;

    fn record(ord: u64, id: u64) -> (Self, Vec<u8>) {
        let body = range_body(i32::try_from(ord).unwrap(), id);
        (Self::extract_from_record(&body), body)
    }

    fn chunk(records: Vec<(Self, Vec<u8>)>) -> MemoryChunkErased {
        MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(records))
    }
}

impl Lane for TemplateKey24 {
    const KIND: SpillKeyKind = SpillKeyKind::TemplateK24;
    const ORDER: SortOrder = SortOrder::TemplateCoordinate;

    fn record(ord: u64, id: u64) -> (Self, Vec<u8>) {
        let key = Self { primary: ord >> 3, secondary: ord & 7, name_hash_upper: 0 };
        (key, range_body(0, id))
    }

    fn chunk(records: Vec<(Self, Vec<u8>)>) -> MemoryChunkErased {
        MemoryChunkErased::TemplateCoordinate(fgumi_sort::TemplateMemChunk::K24(
            InMemoryChunk::from_owned_records(records),
        ))
    }
}

/// Fill a 32-byte key's lanes from `ord` (`set` receives primary, secondary and
/// the `opt` lane).
fn key32_lanes(ord: u64, mut set: impl FnMut(u64, u64, u64)) {
    set(ord >> 3, ord & 7, ord % 3);
}

impl Lane for CbKey32 {
    const KIND: SpillKeyKind = SpillKeyKind::TemplateCb32;
    const ORDER: SortOrder = SortOrder::TemplateCoordinate;

    fn record(ord: u64, id: u64) -> (Self, Vec<u8>) {
        let mut key = Self::default();
        key32_lanes(ord, |p, s, o| (key.0.primary, key.0.secondary, key.0.opt) = (p, s, o));
        (key, range_body(0, id))
    }

    fn chunk(records: Vec<(Self, Vec<u8>)>) -> MemoryChunkErased {
        MemoryChunkErased::TemplateCoordinate(fgumi_sort::TemplateMemChunk::Cb32(
            InMemoryChunk::from_owned_records(records),
        ))
    }
}

impl Lane for TertKey32 {
    const KIND: SpillKeyKind = SpillKeyKind::TemplateTert32;
    const ORDER: SortOrder = SortOrder::TemplateCoordinate;

    fn record(ord: u64, id: u64) -> (Self, Vec<u8>) {
        let mut key = Self::default();
        key32_lanes(ord, |p, s, o| (key.0.primary, key.0.secondary, key.0.opt) = (p, s, o));
        (key, range_body(0, id))
    }

    fn chunk(records: Vec<(Self, Vec<u8>)>) -> MemoryChunkErased {
        MemoryChunkErased::TemplateCoordinate(fgumi_sort::TemplateMemChunk::Tert32(
            InMemoryChunk::from_owned_records(records),
        ))
    }
}

impl Lane for TemplateKey {
    const KIND: SpillKeyKind = SpillKeyKind::TemplateK40;
    const ORDER: SortOrder = SortOrder::TemplateCoordinate;

    fn record(ord: u64, id: u64) -> (Self, Vec<u8>) {
        let key = Self {
            primary: ord >> 3,
            secondary: ord & 7,
            cb_hash: 0,
            tertiary: ord % 3,
            name_hash_upper: 0,
        };
        (key, range_body(0, id))
    }

    fn chunk(records: Vec<(Self, Vec<u8>)>) -> MemoryChunkErased {
        MemoryChunkErased::TemplateCoordinate(fgumi_sort::TemplateMemChunk::K40(
            InMemoryChunk::from_owned_records(records),
        ))
    }
}

/// The merge's leaves: spill runs (by `file_id`) and the residual, each sorted.
/// Range mode merges spill runs only, so every range test generates an empty
/// residual, except the one that checks range mode refuses a residual chunk.
struct RangeInputs<K> {
    runs: Vec<Vec<(K, Vec<u8>)>>,
    residual: Vec<(K, Vec<u8>)>,
}

impl<K: Lane> RangeInputs<K> {
    /// `runs` runs of `per_run` records and a residual of `residual`; with
    /// `ties`, keys repeat within and across leaves (the leaf-index tie-break
    /// decides their order), otherwise every key is distinct.
    fn generate(runs: usize, per_run: usize, residual: usize, ties: bool) -> Self {
        let total = (runs * per_run + residual) as u64;
        let mut state = 0x5EED_u64;
        let mut next = move || {
            // splitmix64: the data is a pure function of the seed.
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        };
        let mut id = 0u64;
        let mut leaf = |n: usize| {
            let mut recs: Vec<(K, Vec<u8>)> = (0..n)
                .map(|_| {
                    id += 1;
                    let ord = if ties { next() % (total / 2) } else { next() % (1 << 28) };
                    K::record(ord, id)
                })
                .collect();
            recs.sort_by(|a, b| a.0.cmp(&b.0));
            recs
        };
        let runs = (0..runs).map(|_| leaf(per_run)).collect();
        let residual = leaf(residual);
        Self { runs, residual }
    }

    fn total(&self) -> u64 {
        (self.runs.iter().map(Vec::len).sum::<usize>() + self.residual.len()) as u64
    }

    /// Every run as an indexed zstd spill run under `dir`, `file_id` = its
    /// position.
    fn write_runs(&self, dir: &std::path::Path) -> Vec<Arc<SortMergeSlot>> {
        let slots = self.write_runs_without_kind(dir);
        for slot in &slots {
            slot.attach_key_kind(K::KIND);
        }
        slots
    }

    /// [`Self::write_runs`] without the slots' key kind attached.
    fn write_runs_without_kind(&self, dir: &std::path::Path) -> Vec<Arc<SortMergeSlot>> {
        self.runs
            .iter()
            .enumerate()
            .map(|(i, run)| {
                write_indexed_run(dir, u32::try_from(i).unwrap(), K::KIND, run, RANGE_FRAME_BYTES)
            })
            .collect()
    }

    /// The merge's phase-2 input over `slots` and the residual.
    fn events(&self, slots: &[Arc<SortMergeSlot>]) -> Vec<SortPhase2Event> {
        let total = self.total();
        let mut events: Vec<SortPhase2Event> = slots
            .iter()
            .map(|slot| SortPhase2Event::SpillReady {
                slot: Arc::clone(slot),
                path: std::path::PathBuf::from("spill"),
                records_ingested_so_far: total,
            })
            .collect();
        let chunks = u32::from(!self.residual.is_empty());
        if chunks > 0 {
            events.push(SortPhase2Event::MemoryChunk {
                chunk: Arc::new(K::chunk(self.residual.clone())),
                records_ingested_so_far: total,
            });
        }
        events.push(SortPhase2Event::AllAnnounced {
            slot_count: u32::try_from(slots.len()).unwrap(),
            memory_chunk_count: chunks,
            total_records: total,
        });
        events
    }

    /// The serial merge's record order (leaf-index tie-break).
    fn expected_bodies(&self) -> Vec<Vec<u8>> {
        let mut leaves = self.runs.clone();
        leaves.push(self.residual.clone());
        serial_reference(&leaves)
    }
}

/// Read `slot`'s whole run into its stash, as the read planner would, so the
/// serial merge can serve itself every frame without a planner or decompress
/// step.
fn stash_whole_run(slot: &SortMergeSlot) {
    let pool = fgumi_bam_io::pread::SliceBufferPool::new(0);
    let start = slot.body_start();
    let mut body = vec![0u8; usize::try_from(slot.len() - start).unwrap()];
    fgumi_bam_io::pread::read_at_exact(&**slot.source(), &mut body, start).unwrap();
    slot.bp_note_issued(body.len() as u64);
    slot.bp_ingest_slice(0, pool.lease(body), true).expect("ingest the run");
}

/// A spill supply whose merge mode is already published as `mode`.
fn supply_with(mode: MergeMode) -> SpillSupply {
    let supply = SpillSupply::new();
    supply.merge_mode.set(mode).expect("a fresh supply has no mode");
    supply
}

/// The range-mode decision the tests publish for a pool `threads` wide: the
/// fit of three small runs within 256 MiB, a window of two ranges per worker
/// and at least four (the target only sizes the plan, which the tests
/// override).
fn range_mode(threads: usize) -> MergeMode {
    let elig = crate::sort::RangeEligibility::for_test(threads, 256 << 20, 256 << 20);
    let inputs = crate::sort::FitInputs {
        k: 3,
        frames: 100,
        key_width: 8,
        raw_bytes: 1 << 20,
        comp_bytes: 1 << 19,
        max_frame_raw: 0,
    };
    let fit = crate::sort::range_fit(&elig, &inputs).expect("three small runs fit");
    assert_eq!(fit.window, (2 * threads).max(4), "{fit:?}");
    MergeMode::Range { budget: elig.budget(), inputs, fit, index_bytes: 0 }
}

/// Emitted items, each as its ordinal and record bodies (so record
/// boundaries compare too).
type ItemsRecords = Vec<(u64, Vec<Vec<u8>>)>;

/// An emitted item's ordinal and record bodies, so record boundaries compare
/// too.
trait ItemRecords {
    fn ordinal_and_records(&self) -> (u64, Vec<Vec<u8>>);
}

impl ItemRecords for DecompressedBlock {
    fn ordinal_and_records(&self) -> (u64, Vec<Vec<u8>>) {
        let mut records = Vec::new();
        let mut span = self.bytes.as_slice();
        while let Some((prefix, rest)) = span.split_first_chunk::<4>() {
            let (body, rest) = rest.split_at(u32::from_le_bytes(*prefix) as usize);
            records.push(body.to_vec());
            span = rest;
        }
        assert!(span.is_empty(), "a block ends on a record boundary");
        (self.batch_serial, records)
    }
}

impl ItemRecords for RecordBatch {
    fn ordinal_and_records(&self) -> (u64, Vec<Vec<u8>>) {
        (self.ordinal(), self.iter_record_bytes().map(<[u8]>::to_vec).collect())
    }
}

/// A serial sink that collects every item it receives.
struct CollectSink<T> {
    received: Arc<parking_lot::Mutex<Vec<T>>>,
}

impl<T: Send + HeapSize + 'static> Step for CollectSink<T> {
    type Input = T;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "CollectSink",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        match ctx.input.pop() {
            Some(item) => {
                self.received.lock().push(item);
                Ok(StepOutcome::Progress)
            }
            None if ctx.input.is_drained() => Ok(StepOutcome::Finished),
            None => Ok(StepOutcome::NoProgress),
        }
    }
}

/// Run `f` on its own thread and fail the test if it has not returned within
/// a minute (a wedged merge).
fn within_watchdog<T: Send + 'static>(f: impl FnOnce() -> T + Send + 'static) -> T {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(f());
    });
    rx.recv_timeout(std::time::Duration::from_secs(60)).expect("WEDGED: the merge did not finish")
}

/// What a range-mode test merge run is given.
struct MergeRun {
    events: Vec<SortPhase2Event>,
    supply: SpillSupply,
    order: SortOrder,
    threads: usize,
    sizing: Option<RangeSizing>,
    /// `(target_batch_count, output_byte_limit)`.
    batch: (usize, u64),
    sort_stats: bool,
    /// The range-mode claim policy (`None`: the environment's, i.e. `auto`).
    claim: Option<super::claim::ClaimPolicy>,
    /// The job (ticket) edge's byte bound (`None`: the default).
    job_edge: Option<u64>,
}

/// Run `SortMerge<O>` (built `with_pooled_gather(threads)` and joined to
/// `supply`) and its `MergeWork` through `Pipeline::run` on `threads` pool
/// threads, under a watchdog, and return every emitted item's ordinal and
/// records.
fn run_merge<O: MergeOutput>(run: MergeRun) -> Result<ItemsRecords, PipelineError>
where
    O::Item: ItemRecords,
{
    within_watchdog(move || {
        let received = Arc::new(parking_lot::Mutex::new(Vec::<O::Item>::new()));
        let (target, limit) = run.batch;
        let mut merge = SortMerge::<O>::with_target_batch_count(run.order, limit, target)
            .with_spill_supply(&run.supply)
            .with_sort_stats(run.sort_stats);
        if let Some(sizing) = run.sizing {
            merge = merge.with_range_sizing_for_test(sizing);
        }
        if let Some(policy) = run.claim {
            merge = merge.with_claim_policy(policy);
        }
        if let Some(limit) = run.job_edge {
            merge = merge.with_job_edge_limit(limit);
        }
        let (merge, work) = merge.with_pooled_gather(run.threads, None);
        let builder = Pipeline::builder();
        let fan = builder
            .chain(crate::sort::tests::Phase2EventSource::new(run.events, limit))
            .chain(merge)
            .into_multi();
        fan.b0.chain(CollectSink { received: Arc::clone(&received) }).into_sink_marker();
        fan.b1.chain(work).into_sink_marker();
        builder.build().expect("build the pipeline").run(
            PipelineConfig { threads: run.threads, ..Default::default() }.with_deadlock_timeout(10),
        )?;
        let items = std::mem::take(&mut *received.lock());
        Ok(items.iter().map(ItemRecords::ordinal_and_records).collect())
    })
}

/// Range mode and the serial merge (mode `Serial(ForcedByEnv)`, every run
/// stashed so it serves itself) over the same leaves; returns both outputs and
/// the oracle's record order.
fn range_and_serial<K: Lane, O: MergeOutput>(
    threads: usize,
    sizing: RangeSizing,
    batch: (usize, u64),
) -> (ItemsRecords, ItemsRecords, Vec<Vec<u8>>)
where
    O::Item: ItemRecords,
{
    range_and_serial_claiming::<K, O>(threads, sizing, batch, None)
}

/// [`range_and_serial`] with range mode under claim policy `claim`.
fn range_and_serial_claiming<K: Lane, O: MergeOutput>(
    threads: usize,
    sizing: RangeSizing,
    batch: (usize, u64),
    claim: Option<super::claim::ClaimPolicy>,
) -> (ItemsRecords, ItemsRecords, Vec<Vec<u8>>)
where
    O::Item: ItemRecords,
{
    let inputs = RangeInputs::<K>::generate(3, 1200, 0, true);
    let range_dir = tempfile::tempdir().unwrap();
    let serial_dir = tempfile::tempdir().unwrap();
    let range_slots = inputs.write_runs(range_dir.path());
    let serial_slots = inputs.write_runs(serial_dir.path());
    for slot in &serial_slots {
        stash_whole_run(slot);
    }
    let serial = run_merge::<O>(MergeRun {
        events: inputs.events(&serial_slots),
        supply: supply_with(MergeMode::Serial(SerialReason::ForcedByEnv)),
        order: K::ORDER,
        threads,
        sizing: None,
        batch,
        sort_stats: false,
        claim: None,
        job_edge: None,
    })
    .expect("serial merge");
    // Range mode cannot fall back to the serial merge here: nothing reads or
    // stashes these slots, so a serial merge over them would never finish.
    let range = run_merge::<O>(MergeRun {
        events: inputs.events(&range_slots),
        supply: supply_with(range_mode(threads)),
        order: K::ORDER,
        threads,
        sizing: Some(sizing),
        batch,
        sort_stats: false,
        claim,
        job_edge: None,
    })
    .expect("range merge");
    (range, serial, inputs.expected_bodies())
}

/// The lane a range-mode parity case runs.
#[derive(Clone, Copy, Debug)]
enum RangeLane {
    Coordinate,
    K24,
    K40,
}

/// The output framing a range-mode parity case emits.
#[derive(Clone, Copy, Debug)]
enum RangeOutput {
    Block,
    Batch,
}

/// Range mode emits exactly the serial merge's items: the same records, the
/// same batch boundaries and the same dense ordinals, for every lane, pool
/// width and output framing; and the serial merge emits the oracle's order.
#[rstest::rstest]
#[case::coordinate(RangeLane::Coordinate)]
#[case::k24(RangeLane::K24)]
#[case::k40(RangeLane::K40)]
fn range_mode_matches_serial(
    #[case] lane: RangeLane,
    #[values(2, 4, 8)] threads: usize,
    #[values(RangeOutput::Block, RangeOutput::Batch)] output: RangeOutput,
) {
    let (range, serial, expected) = match (lane, output) {
        (RangeLane::Coordinate, RangeOutput::Block) => {
            range_and_serial::<RawCoordinateKey, BlockOutput>(threads, RANGE_SIZING, RANGE_BATCH)
        }
        (RangeLane::Coordinate, RangeOutput::Batch) => range_and_serial::<
            RawCoordinateKey,
            RecordBatchOutput,
        >(threads, RANGE_SIZING, RANGE_BATCH),
        (RangeLane::K24, RangeOutput::Block) => {
            range_and_serial::<TemplateKey24, BlockOutput>(threads, RANGE_SIZING, RANGE_BATCH)
        }
        (RangeLane::K24, RangeOutput::Batch) => {
            range_and_serial::<TemplateKey24, RecordBatchOutput>(threads, RANGE_SIZING, RANGE_BATCH)
        }
        (RangeLane::K40, RangeOutput::Block) => {
            range_and_serial::<TemplateKey, BlockOutput>(threads, RANGE_SIZING, RANGE_BATCH)
        }
        (RangeLane::K40, RangeOutput::Batch) => {
            range_and_serial::<TemplateKey, RecordBatchOutput>(threads, RANGE_SIZING, RANGE_BATCH)
        }
    };
    assert!(serial.len() > 10, "the fixture spans many batches ({})", serial.len());
    let serial_bodies: Vec<Vec<u8>> = serial.iter().flat_map(|(_, r)| r.clone()).collect();
    assert_eq!(serial_bodies, expected, "the serial merge emits the oracle's order");
    assert_eq!(range.len(), serial.len(), "batch count");
    for (i, (got, want)) in range.iter().zip(&serial).enumerate() {
        assert_eq!(got.0, i as u64, "batch {i} ordinal is dense");
        assert_eq!(got, want, "batch {i} differs from the serial merge");
    }
}

/// Range mode emits exactly the serial merge's items under every claim
/// policy (the merge thread runs every job, the pool runs every job, or both
/// claim) at one, two and four phase-2 threads: the two extremes bracket
/// every interleaving of who runs what.
#[rstest::rstest]
#[case::coordinate_blocks(RangeLane::Coordinate, RangeOutput::Block)]
#[case::k24_batches(RangeLane::K24, RangeOutput::Batch)]
fn range_mode_matches_serial_under_each_claim_policy(
    #[case] lane: RangeLane,
    #[case] output: RangeOutput,
    #[values(
        super::claim::ClaimPolicy::Auto,
        super::claim::ClaimPolicy::Merge,
        super::claim::ClaimPolicy::Pool
    )]
    claim: super::claim::ClaimPolicy,
    #[values(1, 2, 4)] threads: usize,
) {
    let claim = Some(claim);
    let (range, serial, expected) = match (lane, output) {
        (RangeLane::Coordinate, RangeOutput::Block) => range_and_serial_claiming::<
            RawCoordinateKey,
            BlockOutput,
        >(
            threads, RANGE_SIZING, RANGE_BATCH, claim
        ),
        (RangeLane::K24, RangeOutput::Batch) => range_and_serial_claiming::<
            TemplateKey24,
            RecordBatchOutput,
        >(threads, RANGE_SIZING, RANGE_BATCH, claim),
        other => unreachable!("no case runs {other:?}"),
    };
    let serial_bodies: Vec<Vec<u8>> = serial.iter().flat_map(|(_, r)| r.clone()).collect();
    assert_eq!(serial_bodies, expected, "the serial merge emits the oracle's order");
    assert!(range == serial, "range mode under {claim:?} at {threads} threads differs");
}

/// The batches of ranges cut from merged buffers and of ranges merged
/// straight into builders (fused), in any mix, are exactly the serial merge's
/// batches: the open batch carries across a fused range as segments copied
/// into its lead builder, and out of one as the lead of the next batch plan.
fn check_fused_mix<O: MergeOutput>(ranges: &[Vec<u32>], fused: &[bool], target: usize)
where
    O::Item: ItemRecords,
{
    use std::collections::VecDeque;

    use fgumi_sort::RecordSink;

    use super::range_mode::{FusedSink, begin_fused};
    use crate::sort::merge_job::{RangeBufPool, assemble_batch};
    use crate::sort::range_cut::{BatchCutter, range_out_for_test};

    const LIMIT: u64 = 600_000;
    let body = |r: usize, i: usize, len: u32| -> Vec<u8> {
        (0..len as usize).map(|k| (r * 31 + i * 7 + k).to_le_bytes()[0]).collect()
    };
    let overhead = O::Builder::FRAME_OVERHEAD_PER_RECORD as u64;

    // The serial merge's batches: every record in order, closed by `batch_full`.
    let mut want: Vec<O::Item> = Vec::new();
    let (mut builder, mut count, mut bytes) = (None::<O::Builder>, 0usize, 0u64);
    for (r, lens) in ranges.iter().enumerate() {
        for (i, &len) in lens.iter().enumerate() {
            let ordinal = want.len() as u64;
            builder
                .get_or_insert_with(|| O::Builder::with_capacity(ordinal, 0, target))
                .push_record_bytes(&body(r, i, len))
                .unwrap();
            count += 1;
            bytes += u64::from(len) + overhead;
            if crate::sort::merge::batch_full(count, bytes, target, LIMIT) {
                want.push(builder.take().unwrap().build());
                (count, bytes) = (0, 0);
            }
        }
    }
    want.extend(builder.map(MergeBatchBuilder::build));

    let pool = RangeBufPool::new(4, usize::MAX);
    let mut cutter = BatchCutter::new(target, LIMIT, O::Builder::FRAME_OVERHEAD_PER_RECORD);
    let mut lead: Option<O::Builder> = None;
    let mut got: Vec<O::Item> = Vec::new();
    let assemble = |plans: VecDeque<_>, lead: &mut Option<O::Builder>, got: &mut Vec<_>| {
        for plan in plans {
            let plan: crate::sort::range_cut::BatchPlan = plan;
            let l =
                if plan.lead { Some(lead.take().expect("a lead batch has a lead")) } else { None };
            got.push(assemble_batch::<O>(&plan, target, l).unwrap());
        }
    };
    for (r, lens) in ranges.iter().enumerate() {
        if fused[r] {
            begin_fused(&mut cutter, &mut lead, 0, target, &pool).unwrap();
            let mut closed = Vec::new();
            let mut sink = FusedSink {
                cutter: &mut cutter,
                lead: &mut lead,
                closed: &mut closed,
                bytes_cap: 0,
                target_batch_count: target,
            };
            for (i, &len) in lens.iter().enumerate() {
                sink.push_record(&body(r, i, len)).unwrap();
            }
            got.extend(closed.into_iter().map(MergeBatchBuilder::build));
        } else {
            let out = range_out_for_test(lens, |i, len| body(r, i, len));
            let mut plans = VecDeque::new();
            cutter.cut(r, &out, &mut plans);
            assemble(plans, &mut lead, &mut got);
        }
    }
    assemble(cutter.finish().into_iter().collect(), &mut lead, &mut got);
    assert!(lead.is_none(), "the lead went into the last batch");
    let got: Vec<_> = got.iter().map(ItemRecords::ordinal_and_records).collect();
    let want: Vec<_> = want.iter().map(ItemRecords::ordinal_and_records).collect();
    assert_eq!(got, want);
}

proptest::proptest! {
    #![proptest_config(proptest::prelude::ProptestConfig::with_cases(64))]

    #[test]
    fn fused_and_cut_ranges_match_the_serial_batches(
        ranges in proptest::collection::vec(proptest::collection::vec(1u32..9000, 0..300), 0..8),
        fused in proptest::collection::vec(proptest::prelude::any::<bool>(), 8),
        target in 1usize..400,
        blocks in proptest::prelude::any::<bool>(),
    ) {
        if blocks {
            check_fused_mix::<BlockOutput>(&ranges, &fused, target);
        } else {
            check_fused_mix::<RecordBatchOutput>(&ranges, &fused, target);
        }
    }
}

/// Under the `merge` claim policy the merge thread runs every job in order,
/// so it fuses every range: no range goes through a merged buffer. It keeps
/// taking results while the output (the test's small output edge) refuses
/// blocks, so a range merged ahead still fuses, and with no assembly
/// outstanding every fused batch is streamed to the output during its range.
#[test]
fn merge_policy_fuses_every_range() {
    let logs = capture_logs();
    let inputs = RangeInputs::<RawCoordinateKey>::generate(3, 1200, 0, true);
    let dir = tempfile::tempdir().unwrap();
    let slots = inputs.write_runs(dir.path());
    let items = run_merge::<BlockOutput>(MergeRun {
        events: inputs.events(&slots),
        supply: supply_with(range_mode(2)),
        order: SortOrder::Coordinate,
        threads: 2,
        sizing: Some(RANGE_SIZING),
        batch: RANGE_BATCH,
        sort_stats: true,
        claim: Some(super::claim::ClaimPolicy::Merge),
        job_edge: None,
    })
    .expect("merge");
    let bodies: Vec<Vec<u8>> = items.iter().flat_map(|(_, r)| r.clone()).collect();
    assert_eq!(bodies, inputs.expected_bodies());
    let lines = logs.lines();
    let fused = lines.iter().find(|l| l.starts_with("Range fused: ")).expect("the fused line");
    let numbers: Vec<usize> =
        fused.split(|c: char| !c.is_ascii_digit()).filter_map(|n| n.parse().ok()).collect();
    let merge_line = lines.iter().find(|l| l.starts_with("Range merge: ")).expect("summary");
    let ranges: usize = merge_line["Range merge: ".len()..]
        .split_once(' ')
        .and_then(|(n, _)| n.parse().ok())
        .expect("a range count");
    assert!(ranges > 10, "{merge_line}");
    assert_eq!(numbers[..2], [ranges, ranges], "{fused}");
    // Every batch but the last partial one (posted by the cutter's finish).
    assert!(items.len() - numbers[2] <= 1, "{} batches, {fused}", items.len());
    assert_eq!(numbers[3], numbers[2], "every fused batch streamed: {fused}");
    // No compressor counts pops here, so the gate reads the refused-block
    // backlog.
    let depth = lines
        .iter()
        .find(|l| l.starts_with("Range output depth: "))
        .expect("the output depth line");
    assert!(depth.contains("not counted; claims refused by the gate "), "{depth}");
    assert!(depth.contains("(gate on: workers claim only with no refused block)"), "{depth}");
}

/// While the output refuses blocks, the merge thread keeps merging ahead in
/// order (taking results, so every range still fuses) until the window is
/// full: a block waiting for the output keeps its range live, so the merge
/// parks with ranges left rather than queueing the whole sort. Each time it
/// parks, the output drains; the output is the serial order and every range
/// fused.
#[test]
fn held_output_merges_ahead_within_the_window() {
    use fgumi_pipeline_core::InputHandle;
    use fgumi_pipeline_core::testing::StepProbe;

    let inputs = RangeInputs::<RawCoordinateKey>::generate(3, 1500, 0, false);
    let dir = tempfile::tempdir().unwrap();
    let slots = inputs.write_runs(dir.path());
    let supply = supply_with(range_mode(2));
    let (target, limit) = RANGE_BATCH;
    let (mut step, _work) =
        SortMerge::<BlockOutput>::with_target_batch_count(SortOrder::Coordinate, limit, target)
            .with_spill_supply(&supply)
            .with_range_sizing_for_test(RANGE_SIZING)
            .with_claim_policy(super::claim::ClaimPolicy::Merge)
            .with_pooled_gather(2, None);
    let mut probe = StepProbe::new(&step);
    let _jobs = probe.take_output::<MergeTicket<BlockOutput>>(1);
    let blocks = probe.take_output::<DecompressedBlock>(0);
    for event in inputs.events(&slots) {
        probe.push_input(event);
    }
    probe.close_input();

    let mut emitted = Vec::new();
    let (mut held_parks, mut finished) = (0usize, false);
    let (mut max_live, mut window, mut ranges) = (0usize, 0usize, 0usize);
    for _ in 0..1_000_000 {
        let outcome = probe.try_run(&mut step).expect("try_run");
        if let SortMergeState::RangeMode(rm) = &step.state {
            max_live = max_live.max(rm.stats.max_live_ranges);
            (_, window, _, ranges) = rm.window_state_for_test();
        }
        if outcome == StepOutcome::Finished {
            finished = true;
            break;
        }
        if outcome == StepOutcome::Progress {
            continue;
        }
        // Parked: only on a refused block with the window full (the merge
        // thread runs every job itself, so nothing else can stop it).
        let SortMergeState::RangeMode(rm) = &step.state else {
            panic!("parked outside range mode: {outcome:?}")
        };
        let (live, window, posted, ranges) = rm.window_state_for_test();
        assert_eq!(outcome, StepOutcome::Contention, "a park on a refused block");
        assert!(step.held.is_held(), "parked with nothing held");
        if posted < ranges {
            assert!(live >= window, "parked with {live} of {window} live ranges");
            held_parks += 1;
        }
        emitted.extend(std::iter::from_fn(|| blocks.pop()));
    }
    assert!(finished, "range mode finished");
    assert!(held_parks > 0, "the merge never filled the window behind a refused block");
    assert!(ranges > 10, "the plan spans several windows ({ranges} ranges)");
    assert!(max_live <= window, "peak live ranges {max_live} above the window {window}");
    emitted.extend(std::iter::from_fn(|| blocks.pop()));
    let bodies: Vec<Vec<u8>> = emitted.iter().flat_map(|b| b.ordinal_and_records().1).collect();
    assert_eq!(bodies, inputs.expected_bodies());
    let ordinals: Vec<u64> = emitted.iter().map(|b| b.batch_serial).collect();
    assert_eq!(ordinals, (0..emitted.len() as u64).collect::<Vec<_>>(), "dense batch order");
}

/// Range mode over sparse run indexes (one entry per `stride` frames, as the
/// phase-1 policy leaves them when the dense index would not fit) emits
/// exactly the serial merge's items: a range bound may fall inside an entry,
/// whose records below the bound are decoded and dropped.
#[rstest::rstest]
fn range_mode_matches_serial_on_sparse_indexes(#[values(2, 4)] stride: u32) {
    let inputs = RangeInputs::<RawCoordinateKey>::generate(3, 1200, 0, true);
    let (range_dir, serial_dir) = (tempfile::tempdir().unwrap(), tempfile::tempdir().unwrap());
    let dense = inputs.write_runs(range_dir.path());
    let range_slots: Vec<_> = dense
        .iter()
        .map(|d| {
            let path = range_dir.path().join(format!("run-{}.zsp", d.file_id));
            let slot = fgumi_sort::open_spill_slot(&path, d.file_id, fgumi_sort::FramesEnd::AtEof)
                .unwrap();
            let mut index = fgumi_sort::RunIndex::clone(d.run_index().unwrap());
            index.sparsify(stride).unwrap();
            slot.attach_run_index(Arc::new(index));
            slot.attach_key_kind(SpillKeyKind::Coordinate);
            slot
        })
        .collect();
    assert!(range_slots.iter().all(|s| s.run_index().unwrap().stride() == stride));
    let serial_slots = inputs.write_runs(serial_dir.path());
    for slot in &serial_slots {
        stash_whole_run(slot);
    }
    let run = |events, mode, sizing| {
        run_merge::<BlockOutput>(MergeRun {
            events,
            supply: supply_with(mode),
            order: SortOrder::Coordinate,
            threads: 4,
            sizing,
            batch: RANGE_BATCH,
            sort_stats: false,
            claim: None,
            job_edge: None,
        })
        .expect("merge")
    };
    let serial =
        run(inputs.events(&serial_slots), MergeMode::Serial(SerialReason::ForcedByEnv), None);
    let range = run(inputs.events(&range_slots), range_mode(4), Some(RANGE_SIZING));
    assert!(range == serial, "sparse range mode differs from the serial merge");
}

/// A pool 2 wide has a window of 4 live ranges (posted, merged but not yet
/// cut, or referenced by an unfinished batch). Under the pool claim policy
/// (the merge thread runs no job) jobs run one at a time, oldest first, only
/// when the merge parks, so the merge sits at its limit as long as possible.
///
/// - `batches_close`: 200 ranges, each closing batches; live ranges never
///   exceed the window and fill it.
/// - `ranges_smaller_than_a_batch`: ranges far smaller than one batch, so the
///   cutter's open batch pins the window's floor; the merge must post past
///   the window (one range at a time, with nothing in flight) or it parks
///   with no job left to wake it. The strict window rule panics here with
///   "the merge parks only with jobs in flight".
#[rstest::rstest]
#[case::batches_close(RangeSizing { target_range_bytes: 1, max_ranges: 200 }, (4, RANGE_LIMIT), Some(200), 4, 4)]
// Past the window's floor every live range lies inside the open batch (fewer
// than its 2000 records, each range at least one), plus one range in flight.
#[case::ranges_smaller_than_a_batch(TINY_RANGE_SIZING, TINY_RANGE_BATCH, None, 5, TINY_RANGE_BATCH.0 + 4 + 1)]
fn range_mode_window_bounds_live_ranges(
    #[case] sizing: RangeSizing,
    #[case] batch: (usize, u64),
    #[case] expected_ranges: Option<usize>,
    #[case] min_max_live: usize,
    #[case] max_max_live: usize,
) {
    use fgumi_pipeline_core::testing::StepProbe;
    use fgumi_pipeline_core::{InputHandle, Job, JobStep};

    let inputs = RangeInputs::<RawCoordinateKey>::generate(3, 1500, 0, false);
    let dir = tempfile::tempdir().unwrap();
    let slots = inputs.write_runs(dir.path());
    let supply = supply_with(range_mode(2));
    let (target, limit) = batch;
    let (mut step, _work) =
        SortMerge::<BlockOutput>::with_target_batch_count(SortOrder::Coordinate, limit, target)
            .with_spill_supply(&supply)
            .with_range_sizing_for_test(sizing)
            .with_claim_policy(super::claim::ClaimPolicy::Pool)
            .with_pooled_gather(2, None);
    let inbox = Arc::clone(step.inbox.as_ref().expect("pooled"));
    let mut probe = StepProbe::new(&step);
    let jobs = probe.take_output::<MergeTicket<BlockOutput>>(1);
    let blocks = probe.take_output::<DecompressedBlock>(0);
    for event in inputs.events(&slots) {
        probe.push_input(event);
    }
    probe.close_input();

    let mut pending = std::collections::VecDeque::new();
    let mut emitted = Vec::new();
    let (mut max_live, mut range_count) = (0usize, 0usize);
    let mut finished = false;
    for _ in 0..1_000_000 {
        let outcome = probe.try_run(&mut step).expect("try_run");
        if let SortMergeState::RangeMode(rm) = &step.state {
            max_live = max_live.max(rm.stats.max_live_ranges);
            range_count = rm.plan.range_count();
        }
        pending.extend(std::iter::from_fn(|| jobs.pop()));
        emitted.extend(std::iter::from_fn(|| blocks.pop()));
        if outcome == StepOutcome::Finished {
            finished = true;
            break;
        }
        if outcome == StepOutcome::NoProgress {
            let mut job: MergeTicket<BlockOutput> =
                pending.pop_front().expect("the merge parks only with jobs in flight");
            // A ticket claims the oldest unstarted job and runs it to the end.
            let JobStep::Done(result) = job.run_slice(MERGE_SLICE).expect("job") else {
                panic!("under the pool policy every ticket finds its job");
            };
            let _ = inbox.deliver(job.ordinal(), result);
        }
    }
    assert!(finished, "range mode finished");
    if let Some(expected) = expected_ranges {
        assert_eq!(range_count, expected, "ranges cut");
    }
    assert!(range_count > 4 * 4, "the plan spans several windows ({range_count} ranges)");
    assert!(
        (min_max_live..=max_max_live).contains(&max_live),
        "peak live ranges {max_live} outside {min_max_live}..={max_max_live}"
    );
    let bodies: Vec<Vec<u8>> = emitted.iter().flat_map(|b| b.ordinal_and_records().1).collect();
    assert_eq!(bodies, inputs.expected_bodies());
    assert_eq!(probe.counter(RECORDS), inputs.total());
}

/// Ranges far smaller than a batch through `Pipeline::run`: range mode
/// finishes (under the watchdog) and emits exactly the serial merge's items.
/// With the strict window rule this run wedges: the open batch pins the floor
/// and nothing is left in flight to wake the merge.
#[test]
fn range_mode_tiny_ranges_match_serial() {
    let (range, serial, expected) =
        range_and_serial::<RawCoordinateKey, BlockOutput>(4, TINY_RANGE_SIZING, TINY_RANGE_BATCH);
    let serial_bodies: Vec<Vec<u8>> = serial.iter().flat_map(|(_, r)| r.clone()).collect();
    assert_eq!(serial_bodies, expected, "the serial merge emits the oracle's order");
    assert!(serial.len() > 1, "several batches ({})", serial.len());
    assert_eq!(range, serial);
}

/// A range job that fails (a spill run truncated under its index) fails the
/// sort with the job's error (`job N: ...`), and the run ends promptly rather
/// than wedging, whichever side ran it: under the `merge` policy the merge's
/// own thread runs (and fuses) every range, so `SortMerge` reports it; under
/// `pool` only `MergeWork` workers run jobs; under `auto` either may.
/// Columns: the claim policy, the steps that may report the error.
#[rstest::rstest]
#[case::merge_thread(Some(super::claim::ClaimPolicy::Merge), vec!["SortMerge"])]
#[case::worker(Some(super::claim::ClaimPolicy::Pool), vec!["MergeWork"])]
#[case::either(None, vec!["MergeWork", "SortMerge"])]
fn range_job_error_fails_the_sort(
    #[case] claim: Option<super::claim::ClaimPolicy>,
    #[case] expected_steps: Vec<&str>,
) {
    let inputs = RangeInputs::<RawCoordinateKey>::generate(3, 1200, 0, true);
    let dir = tempfile::tempdir().unwrap();
    let slots = inputs.write_runs(dir.path());
    let path = dir.path().join("run-1.zsp");
    let len = std::fs::metadata(&path).unwrap().len();
    std::fs::OpenOptions::new().write(true).open(&path).unwrap().set_len(len / 2).unwrap();
    let result = run_merge::<BlockOutput>(MergeRun {
        events: inputs.events(&slots),
        supply: supply_with(range_mode(4)),
        order: SortOrder::Coordinate,
        threads: 4,
        sizing: Some(RANGE_SIZING),
        batch: RANGE_BATCH,
        sort_stats: false,
        claim,
        job_edge: None,
    });
    match result {
        Err(PipelineError::Io { step, source }) => {
            assert!(expected_steps.contains(&step), "{step} reported the error");
            let msg = source.to_string();
            assert!(msg.contains("job "), "{msg}");
            assert!(msg.contains("spill slot 1 truncated"), "{msg}");
        }
        other => panic!("expected the range job's error, got {other:?}"),
    }
}

/// A merge told to range-merge without the pooled inbox (nothing would run
/// its jobs) fails at the transition instead of merging serially over slots
/// the planner no longer feeds.
#[test]
fn told_range_without_inbox_errors() {
    use fgumi_pipeline_core::testing::StepProbe;

    let inputs = RangeInputs::<RawCoordinateKey>::generate(1, 100, 0, false);
    let dir = tempfile::tempdir().unwrap();
    let slots = inputs.write_runs(dir.path());
    let supply = supply_with(range_mode(1));
    let mut step = SortMerge::<RecordBatchOutput>::new(SortOrder::Coordinate, RANGE_LIMIT)
        .with_spill_supply(&supply);
    let probe = StepProbe::new(&step);
    for event in inputs.events(&slots) {
        probe.push_input(event);
    }
    let err = probe.try_run(&mut step).expect_err("range mode without an inbox must fail");
    assert!(err.to_string().contains("told to range-merge but cannot"), "{err}");
}

/// Every formatted `log` record, for the `--sort-stats` line checks.
static CAPTURED_LOGS: std::sync::OnceLock<std::sync::Mutex<Vec<String>>> =
    std::sync::OnceLock::new();

/// A `log::Log` that appends each message to [`CAPTURED_LOGS`].
struct CaptureLogger;

impl log::Log for CaptureLogger {
    fn enabled(&self, _: &log::Metadata<'_>) -> bool {
        true
    }

    fn log(&self, record: &log::Record<'_>) {
        if let Some(buf) = CAPTURED_LOGS.get() {
            buf.lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .push(format!("{}", record.args()));
        }
    }

    fn flush(&self) {}
}

/// A test's exclusive hold on [`CAPTURED_LOGS`]: the logger and its buffer
/// are process-global, and plain `cargo test` runs tests as threads of one
/// process, so two log-reading tests at once would read each other's lines.
/// Only these tests turn on `--sort-stats`, so while one holds the capture
/// the lines it reads are its own.
struct LogCapture {
    _exclusive: std::sync::MutexGuard<'static, ()>,
    buffer: &'static std::sync::Mutex<Vec<String>>,
}

impl LogCapture {
    /// Every line logged since the capture began.
    fn lines(&self) -> Vec<String> {
        self.buffer.lock().unwrap_or_else(std::sync::PoisonError::into_inner).clone()
    }
}

/// The capture buffer, created on first use.
fn captured_buffer() -> &'static std::sync::Mutex<Vec<String>> {
    CAPTURED_LOGS.get_or_init(|| std::sync::Mutex::new(Vec::new()))
}

/// Install [`CaptureLogger`] (once per process), wait for any other test's
/// capture to end, and start this one with the buffer cleared. Hold the
/// returned capture until the test has read its lines.
fn capture_logs() -> LogCapture {
    static INIT: std::sync::Once = std::sync::Once::new();
    static LOGGER: CaptureLogger = CaptureLogger;
    static EXCLUSIVE: std::sync::Mutex<()> = std::sync::Mutex::new(());
    // A test that panicked holding the capture poisons it; the next one
    // starts over with a cleared buffer, so the poison is harmless.
    let exclusive = EXCLUSIVE.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    INIT.call_once(|| {
        let _ = log::set_logger(&LOGGER);
        log::set_max_level(log::LevelFilter::Info);
    });
    let buffer = captured_buffer();
    buffer.lock().unwrap_or_else(std::sync::PoisonError::into_inner).clear();
    LogCapture { _exclusive: exclusive, buffer }
}

/// Under `--sort-stats` the merge logs its mode at the transition and, in
/// range mode, the range merge's summary at the end. Each case lists the
/// line prefixes and the substrings each must contain.
#[rstest::rstest]
#[case::range(
    range_mode(4),
    vec![
        ("Sort merge mode: range (", vec![" ranges, window 8, 4 workers, ~0 MiB per range, index stride 1, index 0 MiB, mem 80 MiB of range budget 256 MiB (merge memory 256 MiB, floor R 256 MiB), planned decode amplification 1.050, plan amplification ", ", rho 2.50, final run spilled)"]),
        ("Range merge: ", vec![" on the merge thread, which also assembled ", " batches), 3600 records, largest range ", "% of records, decode amplification ", " MiB read"]),
        ("Range merge thread: busy ", vec!["(own ranges ", ", own assemblies ", " batches, compression ", " blocks, other ", "; idle ", ": waiting on worker ", " s, output held ", " s, window full ", " s, draining ", " s, waiting on nothing (self-check) ", "; between busy dispatches "]),
        ("Range window occupancy (live ranges: share of time): ", vec!["%"]),
        ("Range output backlog: blocks refused by the output edge ", vec!["% of the time"]),
        ("Range output depth: not counted; claims refused by the gate ", vec!["(gate on: workers claim only with no refused block)"]),
        ("Range claims (policy auto): merge thread ", vec![" ranges, ", " assemblies; workers ", "; spent tickets ", " refused while compression was behind, gate on); table high-water "]),
        ("Range fused: ", vec![" of the merge thread's ", " output batches (", " streamed to the output during their range)"]),
        ("Range output: merge thread compressed ", vec![" blocks (", "%; merge-thread compression not configured)"]),
        ("Range output buffers: ", vec![" ranges reused a recycled buffer, ", " buffers recycled, ", " pooled at the end (at most 8 of "]),
    ]
)]
#[case::serial(
    MergeMode::Serial(SerialReason::ForcedByEnv),
    vec![("Sort merge mode: serial (forced by FGUMI_SORT_SERIAL_MERGE)", vec![])]
)]
fn range_mode_logs_stats(#[case] mode: MergeMode, #[case] expected: Vec<(&str, Vec<&str>)>) {
    let logs = capture_logs();
    let inputs = RangeInputs::<RawCoordinateKey>::generate(3, 1200, 0, true);
    let dir = tempfile::tempdir().unwrap();
    let slots = inputs.write_runs(dir.path());
    if matches!(mode, MergeMode::Serial(_)) {
        for slot in &slots {
            stash_whole_run(slot);
        }
    }
    run_merge::<BlockOutput>(MergeRun {
        events: inputs.events(&slots),
        supply: supply_with(mode),
        order: SortOrder::Coordinate,
        threads: 4,
        sizing: Some(RANGE_SIZING),
        batch: RANGE_BATCH,
        sort_stats: true,
        claim: None,
        job_edge: None,
    })
    .expect("merge");
    let lines = logs.lines();
    for (prefix, parts) in expected {
        let line = lines
            .iter()
            .find(|l| l.starts_with(prefix))
            .unwrap_or_else(|| panic!("no line starting {prefix:?} in {lines:#?}"));
        for part in parts {
            assert!(line.contains(part), "{line:?} lacks {part:?}");
        }
    }
    // The fixture's plan really cuts many ranges (it is the parity tests'
    // fixture and sizing), so range mode is exercised across range edges.
    if let Some(line) = lines.iter().find(|l| l.starts_with("Range merge: ")) {
        let ranges: usize = line["Range merge: ".len()..]
            .split_once(' ')
            .and_then(|(n, _)| n.parse().ok())
            .unwrap_or_else(|| panic!("no range count in {line:?}"));
        assert!(ranges > 10, "only {ranges} ranges: {line:?}");
        // The merge's own thread merges ranges beside the workers (it takes
        // the first range before posting any), and never all of them here.
        let own: usize = line
            .split_once(" ranges (")
            .and_then(|(_, rest)| rest.split_once(' '))
            .and_then(|(n, _)| n.parse().ok())
            .unwrap_or_else(|| panic!("no merge-thread range count in {line:?}"));
        assert!((1..=ranges).contains(&own), "{own} of {ranges} ranges on the merge thread");
        // Every kept record was decoded, so the realized amplification is at
        // least 1 by construction; the check that can fail is that it stays
        // within the plan's (each range decodes at most its entry spans, and a
        // cursor stops at the range's upper bound).
        let figure = |line: &str, label: &str| -> f64 {
            line.split_once(label)
                .and_then(|(_, rest)| rest.split_once(','))
                .and_then(|(x, _)| x.parse().ok())
                .unwrap_or_else(|| panic!("no {label:?} in {line:?}"))
        };
        let amplification = figure(line, "decode amplification ");
        let mode_line = lines
            .iter()
            .find(|l| l.starts_with("Sort merge mode: range ("))
            .expect("the range mode line");
        let planned = figure(mode_line, "plan amplification ");
        assert!(
            (1.0..=planned + 1e-3).contains(&amplification),
            "amplification {amplification} outside 1..={planned}: {line:?}"
        );
    }
}

/// Merge `inputs` (no residual) from `slots` in `mode`, stashing the runs first
/// for the serial merge.
fn merge_slots_only<K: Lane>(
    inputs: &RangeInputs<K>,
    slots: &[Arc<SortMergeSlot>],
    mode: MergeMode,
    sizing: Option<RangeSizing>,
) -> Result<ItemsRecords, PipelineError> {
    if matches!(mode, MergeMode::Serial(_)) {
        for slot in slots {
            stash_whole_run(slot);
        }
    }
    run_merge::<BlockOutput>(MergeRun {
        events: inputs.events(slots),
        supply: supply_with(mode),
        order: K::ORDER,
        threads: 4,
        sizing,
        batch: RANGE_BATCH,
        sort_stats: false,
        claim: None,
        job_edge: None,
    })
}

/// With every residual chunk spilled, the template lane comes from the slots'
/// key kind: both the serial merge and range mode emit the reference order.
fn template_lane_case<K: Lane>() {
    let inputs = RangeInputs::<K>::generate(2, 800, 0, true);
    let serial_dir = tempfile::tempdir().unwrap();
    let serial = merge_slots_only(
        &inputs,
        &inputs.write_runs(serial_dir.path()),
        MergeMode::Serial(SerialReason::ForcedByEnv),
        None,
    )
    .expect("serial merge");
    let serial_bodies: Vec<Vec<u8>> = serial.iter().flat_map(|(_, r)| r.clone()).collect();
    assert_eq!(serial_bodies, inputs.expected_bodies());

    let range_dir = tempfile::tempdir().unwrap();
    let range = merge_slots_only(
        &inputs,
        &inputs.write_runs(range_dir.path()),
        range_mode(4),
        Some(RANGE_SIZING),
    )
    .expect("range merge");
    assert_eq!(range, serial);
}

#[rstest::rstest]
#[case::k24(SpillKeyKind::TemplateK24)]
#[case::cb32(SpillKeyKind::TemplateCb32)]
#[case::tert32(SpillKeyKind::TemplateTert32)]
#[case::k40(SpillKeyKind::TemplateK40)]
fn template_lane_from_slots(#[case] kind: SpillKeyKind) {
    match kind {
        SpillKeyKind::TemplateK24 => template_lane_case::<TemplateKey24>(),
        SpillKeyKind::TemplateCb32 => template_lane_case::<CbKey32>(),
        SpillKeyKind::TemplateTert32 => template_lane_case::<TertKey32>(),
        SpillKeyKind::TemplateK40 => template_lane_case::<TemplateKey>(),
        other => unreachable!("not a template kind: {other:?}"),
    }
}

/// Slots with no key kind and no residual cannot name the lane: fail closed.
#[rstest::rstest]
#[case::serial(MergeMode::Serial(SerialReason::ForcedByEnv), None)]
#[case::range(range_mode(4), Some(RANGE_SIZING))]
fn template_lane_unknown_without_kind(
    #[case] mode: MergeMode,
    #[case] sizing: Option<RangeSizing>,
) {
    let inputs = RangeInputs::<TemplateKey24>::generate(2, 200, 0, true);
    let dir = tempfile::tempdir().unwrap();
    let slots = inputs.write_runs_without_kind(dir.path());
    let err = merge_slots_only(&inputs, &slots, mode, sizing).expect_err("lane is unknown");
    assert!(err.to_string().contains("no residual chunk to identify"), "{err}");
}

/// Range mode never sees a residual chunk (the planner fails closed to
/// `Serial(ResidualInMemory)` first); should one arrive anyway the merge
/// fails rather than drop it.
#[test]
fn range_mode_with_memory_chunks_errors() {
    let inputs = RangeInputs::<RawCoordinateKey>::generate(3, 200, 50, true);
    let dir = tempfile::tempdir().unwrap();
    let slots = inputs.write_runs(dir.path());
    let result = run_merge::<BlockOutput>(MergeRun {
        events: inputs.events(&slots),
        supply: supply_with(range_mode(4)),
        order: SortOrder::Coordinate,
        threads: 4,
        sizing: Some(RANGE_SIZING),
        batch: RANGE_BATCH,
        sort_stats: false,
        claim: None,
        job_edge: None,
    });
    match result {
        Err(PipelineError::Io { source, .. }) => {
            assert!(source.to_string().contains("range mode with"), "{source}");
        }
        other => panic!("expected the residual error, got {other:?}"),
    }
}

/// When the phase-1 index policy dropped every index after the front spilled
/// the final run, the merge runs serially over the spill runs alone (no
/// residual) and emits the reference order.
#[test]
fn serial_merge_over_runs_whose_indexes_were_dropped() {
    let inputs = RangeInputs::<RawCoordinateKey>::generate(3, 600, 0, true);
    let dir = tempfile::tempdir().unwrap();
    let indexed = inputs.write_runs(dir.path());
    let slots: Vec<Arc<SortMergeSlot>> = indexed
        .iter()
        .map(|s| {
            let path = dir.path().join(format!("run-{}.zsp", s.file_id));
            let slot = fgumi_sort::open_spill_slot(&path, s.file_id, fgumi_sort::FramesEnd::AtEof)
                .unwrap();
            slot.attach_key_kind(SpillKeyKind::Coordinate);
            slot
        })
        .collect();
    assert!(slots.iter().all(|s| s.run_index().is_none()));
    let budget = crate::sort::RangeBudget::new(1 << 20, 1 << 20);
    let mode = MergeMode::Serial(SerialReason::IndexesDropped { k: 3, budget });
    let out = merge_slots_only(&inputs, &slots, mode, None).expect("serial merge");
    let bodies: Vec<Vec<u8>> = out.iter().flat_map(|(_, r)| r.clone()).collect();
    assert_eq!(bodies, inputs.expected_bodies());
}

// ─────────────────────────────────────────────────────────────────────────────
// Compression on the merge thread (range mode, `CompressibleBlockOutput`).
// ─────────────────────────────────────────────────────────────────────────────

/// Inflate a BGZF byte stream.
fn inflate(bytes: &[u8]) -> Vec<u8> {
    let mut cursor = std::io::Cursor::new(bytes.to_vec());
    let mut decompressor = fgumi_bgzf::Decompressor::new();
    let mut out = Vec::new();
    loop {
        let blocks = fgumi_bgzf::reader::read_raw_blocks(&mut cursor, 64).expect("valid BGZF");
        if blocks.is_empty() {
            return out;
        }
        for block in &blocks {
            out.extend(fgumi_bgzf::reader::decompress_block(block, &mut decompressor).unwrap());
        }
    }
}

impl ItemRecords for crate::output_compress::OutputBlock {
    fn ordinal_and_records(&self) -> (u64, Vec<Vec<u8>>) {
        use crate::output_compress::OutputBlock;
        match self {
            OutputBlock::Raw(b) => b.ordinal_and_records(),
            OutputBlock::Compressed(b) => {
                let raw =
                    DecompressedBlock { batch_serial: b.batch_serial, bytes: inflate(&b.bytes) };
                raw.ordinal_and_records()
            }
        }
    }
}

/// A test stand-in for the sink's `BgzfCompress<OutputBlock>` (which lives
/// in the binary crate): compresses raw blocks with its own
/// [`OutputCompressor`] and passes the merge's compressed blocks through,
/// as that step does.
struct PassThroughCompress {
    compressor: crate::output_compress::OutputCompressor,
    held: HeldSlot<Unpushed<crate::types::BgzfBlock>>,
    /// The merge's queue-depth counter, popped as `BgzfCompress` pops it.
    depth: Arc<crate::output_compress::OutputDepth>,
}

impl Step for PassThroughCompress {
    type Input = crate::output_compress::OutputBlock;
    type Outputs = fgumi_pipeline_core::outputs::Single<crate::types::BgzfBlock>;

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "PassThroughCompress",
            kind: StepKind::Serial,
            sticky: false,
            output_queues: vec![QueueSpec::ByteBounded { limit_bytes: 1 << 30 }],
            branch_ordering: vec![BranchOrdering::None],
        }
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        use crate::output_compress::OutputBlock;
        if let Some(unpushed) = self.held.take()
            && let Err(again) = ctx.outputs.retry(unpushed)
        {
            self.held.put(again);
            return Ok(StepOutcome::Contention);
        }
        let Some(item) = ctx.input.pop() else {
            return Ok(if ctx.input.is_drained() {
                StepOutcome::Finished
            } else {
                StepOutcome::NoProgress
            });
        };
        self.depth.note_popped();
        let block = match item {
            OutputBlock::Raw(raw) => self.compressor.compress(raw)?,
            OutputBlock::Compressed(block) => block,
        };
        if let Err(unpushed) = ctx.outputs.push(block) {
            self.held.put(unpushed);
        }
        Ok(StepOutcome::Progress)
    }
}

/// An emitted BGZF block's serial, bytes and BAI manifest (physical lengths
/// and each record's `uoffset`, `len` and alignment context).
type CompressedOut = (u64, Vec<u8>, Option<(Vec<u32>, Vec<String>)>);

/// How a compressed range-mode test merge runs.
#[derive(Clone, Copy)]
struct CompressRun {
    order: SortOrder,
    threads: usize,
    /// `Some`: range mode with the merge thread's compression on or off;
    /// `None`: the serial merge.
    compress: Option<bool>,
    index_bam: bool,
    batch: (usize, u64),
}

/// What [`run_compressed`] returns.
struct CompressedRun {
    /// The BGZF blocks the sink received, with their manifests.
    blocks: Vec<CompressedOut>,
    /// The blocks the merge thread compressed, from the `Range output:` line
    /// (`None`: no such line, as in a serial merge).
    compressed: Option<u64>,
    /// Every line the run logged.
    lines: Vec<String>,
}

/// Run `inputs` through `SortMerge<CompressibleBlockOutput>` (with output
/// compression at level 1) and [`PassThroughCompress`] at the same
/// settings. The shared [`OutputDepth`](crate::output_compress::OutputDepth)
/// must count one pop per block the sink received and, in range mode, one
/// push per pop.
fn run_compressed<K: Lane>(inputs: &RangeInputs<K>, run: CompressRun) -> CompressedRun {
    use super::{CompressibleBlockOutput, OutputCompression};
    use crate::output_compress::{OutputCompressor, OutputDepth};

    let logs = capture_logs();
    let depth = Arc::new(OutputDepth::new());
    let run_depth = Arc::clone(&depth);
    let dir = tempfile::tempdir().unwrap();
    let slots = inputs.write_runs(dir.path());
    let mode = if run.compress.is_some() {
        range_mode(run.threads)
    } else {
        for slot in &slots {
            stash_whole_run(slot);
        }
        MergeMode::Serial(SerialReason::ForcedByEnv)
    };
    let supply = supply_with(mode);
    let events = inputs.events(&slots);
    let blocks = within_watchdog(move || {
        let received = Arc::new(parking_lot::Mutex::new(Vec::<crate::types::BgzfBlock>::new()));
        let (target, limit) = run.batch;
        let settings = OutputCompression { level: 1, index_bam: run.index_bam };
        let depth = run_depth;
        let mut merge =
            SortMerge::<CompressibleBlockOutput>::with_target_batch_count(K::ORDER, limit, target)
                .with_spill_supply(&supply)
                .with_sort_stats(true)
                .with_output_depth(Arc::clone(&depth));
        // The merge compresses only blocks of an output it was given the
        // compression of.
        if run.compress != Some(false) {
            merge = merge.with_output_compression(settings);
        }
        if run.compress.is_some() {
            merge = merge.with_range_sizing_for_test(RANGE_SIZING);
        }
        let (merge, work) = merge.with_pooled_gather(run.threads, None);
        let builder = Pipeline::builder();
        let fan = builder
            .chain(crate::sort::tests::Phase2EventSource::new(events, limit))
            .chain(merge)
            .into_multi();
        fan.b0
            .chain(PassThroughCompress {
                compressor: OutputCompressor::new(1, run.index_bam),
                held: HeldSlot::new(),
                depth,
            })
            .chain(CollectSink { received: Arc::clone(&received) })
            .into_sink_marker();
        fan.b1.chain(work).into_sink_marker();
        builder.build().expect("build the pipeline").run(
            PipelineConfig { threads: run.threads, ..Default::default() }.with_deadlock_timeout(10),
        )?;
        let blocks = std::mem::take(&mut *received.lock());
        Ok::<_, PipelineError>(blocks)
    })
    .expect("merge");
    let _ = run.order;
    let out: Vec<CompressedOut> = blocks
        .iter()
        .map(|b| {
            let manifest = b.index.as_ref().map(|m| {
                let records =
                    m.records.iter().map(|r| format!("{} {} {:?}", r.uoffset, r.len, r.ctx));
                (m.phys_comp_len.clone(), records.collect())
            });
            (b.batch_serial, b.bytes.to_vec(), manifest)
        })
        .collect();
    assert_eq!(depth.popped(), out.len() as u64, "one counted pop per block");
    // A range merge counts every block it pushes, so the raw counts match (a
    // saturating `depth() == 0` would also pass with pushes left uncounted);
    // the serial merge counts none.
    let pushed = if run.compress.is_some() { depth.popped() } else { 0 };
    assert_eq!(depth.pushed_for_test(), pushed, "one counted push per block popped");
    let lines = logs.lines();
    let compressed = lines.iter().find_map(|l| {
        let rest = l.strip_prefix("Range output: merge thread compressed ")?;
        let (n, _) = rest.split_once(' ')?;
        Some(n.parse().unwrap_or_else(|e| panic!("unparsable count in {l:?}: {e}")))
    });
    CompressedRun { blocks: out, compressed, lines }
}

/// The BGZF stream and the BAI manifests the sink builds from are
/// byte-identical whether the merge thread compressed blocks itself or left
/// every block to the sink's compressor (a merge given no output compression),
/// and equal the serial merge's, at one, two and four phase-2 threads, for a coordinate
/// merge (indexed) and a template-coordinate one. The output edge is the
/// tests' 8 KiB, so it refuses blocks and the merge thread compresses, at
/// one phase-2 thread too.
#[rstest::rstest]
#[case::coordinate_indexed(SortOrder::Coordinate, true)]
#[case::template_coordinate(SortOrder::TemplateCoordinate, false)]
fn merge_thread_compression_is_byte_identical(
    #[case] order: SortOrder,
    #[case] index_bam: bool,
    #[values(1, 2, 4)] threads: usize,
) {
    fn check<K: Lane>(order: SortOrder, index_bam: bool, threads: usize) {
        let inputs = RangeInputs::<K>::generate(3, 1200, 0, true);
        let run =
            |compress| CompressRun { order, threads, compress, index_bam, batch: RANGE_BATCH };
        let serial = run_compressed(&inputs, run(None));
        let off = run_compressed(&inputs, run(Some(false)));
        let on = run_compressed(&inputs, run(Some(true)));
        assert_eq!(serial.compressed, None, "the serial merge logs no range output line");
        let off_compressed = off.compressed.expect("the range output line (compression off)");
        let on_compressed = on.compressed.expect("the range output line (compression on)");
        assert_eq!(off_compressed, 0, "compression off: the merge thread never compresses");
        assert!(on_compressed > 0, "the merge thread compressed nothing at {threads} threads");
        let (serial, off, on) = (serial.blocks, off.blocks, on.blocks);
        assert!(serial.len() > 10, "the fixture spans many blocks ({})", serial.len());
        assert!(serial.iter().all(|b| b.2.is_some() == index_bam), "manifests iff indexing");
        assert!(off == serial, "range (compression off) differs from serial");
        assert!(on == serial, "range (merge-thread compression) differs from serial");
        let records: Vec<Vec<u8>> = on
            .iter()
            .flat_map(|(s, bytes, _)| {
                DecompressedBlock { batch_serial: *s, bytes: inflate(bytes) }
                    .ordinal_and_records()
                    .1
            })
            .collect();
        assert_eq!(records, inputs.expected_bodies());
    }
    match order {
        SortOrder::Coordinate => check::<RawCoordinateKey>(order, index_bam, threads),
        _ => check::<TemplateKey24>(order, index_bam, threads),
    }
}

/// With an output edge that never refuses (a byte limit far above the
/// output), the merge thread never compresses: compression only fills time
/// it would otherwise spend parked on a refused block.
#[test]
fn merge_thread_does_not_compress_an_accepted_output() {
    let inputs = RangeInputs::<RawCoordinateKey>::generate(3, 1200, 0, true);
    let run = CompressRun {
        order: SortOrder::Coordinate,
        threads: 2,
        compress: Some(true),
        index_bam: true,
        batch: (RANGE_TARGET, 1 << 30),
    };
    let CompressedRun { blocks, compressed, lines } = run_compressed(&inputs, run);
    assert!(blocks.len() > 10, "{} blocks", blocks.len());
    assert_eq!(compressed, Some(0), "the range output line, with nothing compressed");
    let depth = lines
        .iter()
        .find(|l| l.starts_with("Range output depth: compression queue >= 2 blocks "))
        .unwrap_or_else(|| panic!("no counted depth line in {lines:#?}"));
    assert!(depth.contains("(gate on: workers claim only below 2)"), "{depth}");
}

/// An intermediate (`RecordBatchOutput`) range merge never compresses, even
/// given output compression settings and an output that refuses batches:
/// its items are record batches for a decode step, not BGZF blocks.
#[test]
fn record_batch_output_is_never_compressed() {
    let logs = capture_logs();
    let inputs = RangeInputs::<RawCoordinateKey>::generate(3, 1200, 0, true);
    let dir = tempfile::tempdir().unwrap();
    let slots = inputs.write_runs(dir.path());
    let supply = supply_with(range_mode(2));
    let events = inputs.events(&slots);
    let items = within_watchdog(move || {
        let received = Arc::new(parking_lot::Mutex::new(Vec::<RecordBatch>::new()));
        let (target, limit) = RANGE_BATCH;
        let (merge, work) = SortMerge::<RecordBatchOutput>::with_target_batch_count(
            SortOrder::Coordinate,
            limit,
            target,
        )
        .with_spill_supply(&supply)
        .with_sort_stats(true)
        .with_range_sizing_for_test(RANGE_SIZING)
        .with_output_compression(super::OutputCompression { level: 1, index_bam: true })
        .with_pooled_gather(2, None);
        let builder = Pipeline::builder();
        let fan = builder
            .chain(crate::sort::tests::Phase2EventSource::new(events, limit))
            .chain(merge)
            .into_multi();
        fan.b0.chain(CollectSink { received: Arc::clone(&received) }).into_sink_marker();
        fan.b1.chain(work).into_sink_marker();
        builder
            .build()
            .expect("build the pipeline")
            .run(PipelineConfig { threads: 2, ..Default::default() }.with_deadlock_timeout(10))
            .expect("merge");
        let items = std::mem::take(&mut *received.lock());
        items.iter().map(ItemRecords::ordinal_and_records).collect::<Vec<_>>()
    });
    let bodies: Vec<Vec<u8>> = items.iter().flat_map(|(_, r)| r.clone()).collect();
    assert_eq!(bodies, inputs.expected_bodies());
    let lines = logs.lines();
    let line = lines
        .iter()
        .find(|l| l.starts_with("Range output: merge thread compressed "))
        .expect("the range output line");
    assert!(line.starts_with("Range output: merge thread compressed 0 of "), "{line}");
    let held = lines
        .iter()
        .find(|l| l.starts_with("Range merge thread: busy "))
        .expect("the merge thread line");
    assert!(held.contains("compression 0.0 s over 0 blocks"), "{held}");
}

/// With its output stalled and the window full, the merge thread compresses
/// every waiting block (the held one first) before it parks; it then parks
/// on `Contention` with nothing left to compress, and resumes when the
/// output drains. The output, part compressed, is the serial order.
#[test]
fn merge_thread_compresses_before_parking() {
    use fgumi_pipeline_core::InputHandle;
    use fgumi_pipeline_core::testing::StepProbe;

    use super::{CompressibleBlockOutput, OutputCompression};
    use crate::output_compress::OutputBlock;

    let inputs = RangeInputs::<RawCoordinateKey>::generate(3, 1500, 0, false);
    let dir = tempfile::tempdir().unwrap();
    let slots = inputs.write_runs(dir.path());
    let supply = supply_with(range_mode(2));
    let (target, limit) = RANGE_BATCH;
    let (mut step, _work) = SortMerge::<CompressibleBlockOutput>::with_target_batch_count(
        SortOrder::Coordinate,
        limit,
        target,
    )
    .with_spill_supply(&supply)
    .with_range_sizing_for_test(RANGE_SIZING)
    .with_claim_policy(super::claim::ClaimPolicy::Merge)
    .with_output_compression(OutputCompression { level: 1, index_bam: false })
    .with_pooled_gather(2, None);
    let mut probe = StepProbe::new(&step);
    let _jobs = probe.take_output::<MergeTicket<CompressibleBlockOutput>>(1);
    let blocks = probe.take_output::<OutputBlock>(0);
    for event in inputs.events(&slots) {
        probe.push_input(event);
    }
    probe.close_input();

    let mut emitted = Vec::new();
    let (mut parks, mut finished) = (0usize, false);
    for _ in 0..1_000_000 {
        let outcome = probe.try_run(&mut step).expect("try_run");
        if outcome == StepOutcome::Finished {
            finished = true;
            break;
        }
        if outcome == StepOutcome::Progress {
            continue;
        }
        assert_eq!(outcome, StepOutcome::Contention, "a park on a refused block");
        let held = step.held.take().expect("parked with nothing held");
        assert!(
            matches!(held.item(), OutputBlock::Compressed(_)),
            "parked with the held block raw"
        );
        step.held.put(held);
        let SortMergeState::RangeMode(rm) = &step.state else {
            panic!("parked outside range mode")
        };
        assert_eq!(
            rm.raw_waiting_for_test(|b| matches!(b, OutputBlock::Raw(_))),
            0,
            "parked with a waiting block left to compress"
        );
        parks += 1;
        emitted.extend(std::iter::from_fn(|| blocks.pop()));
    }
    assert!(finished, "range mode finished");
    assert!(parks > 0, "the output never stalled the merge");
    emitted.extend(std::iter::from_fn(|| blocks.pop()));
    let compressed = emitted.iter().filter(|b| matches!(b, OutputBlock::Compressed(_))).count();
    assert!(
        compressed > 0 && compressed < emitted.len(),
        "{compressed} of {} compressed",
        emitted.len()
    );
    let bodies: Vec<Vec<u8>> = emitted.iter().flat_map(|b| b.ordinal_and_records().1).collect();
    assert_eq!(bodies, inputs.expected_bodies());
    let ordinals: Vec<u64> = emitted.iter().map(OutputBlock::batch_serial).collect();
    assert_eq!(ordinals, (0..emitted.len() as u64).collect::<Vec<_>>(), "dense batch order");
}

// ─────────────────────────────────────────────────────────────────────────────
// Posting independent of the ticket edge; the idle reasons.
// ─────────────────────────────────────────────────────────────────────────────

/// At one phase-2 thread with a ticket edge that holds one ticket and is
/// never popped (no worker runs), the merge thread still posts every range
/// up to its window and runs every job itself: posting does not wait on the
/// ticket edge. (Posting only while no ticket was held stopped after the
/// first range here, with nothing left to wake the merge.) The ticket-edge
/// idle self-check reads zero.
#[test]
fn posting_runs_ahead_of_an_unpopped_ticket_edge() {
    use fgumi_pipeline_core::InputHandle;
    use fgumi_pipeline_core::testing::StepProbe;

    let logs = capture_logs();
    let inputs = RangeInputs::<RawCoordinateKey>::generate(3, 1500, 0, false);
    let dir = tempfile::tempdir().unwrap();
    let slots = inputs.write_runs(dir.path());
    let supply = supply_with(range_mode(1));
    let (target, limit) = RANGE_BATCH;
    let one_ticket = std::mem::size_of::<MergeTicket<BlockOutput>>() as u64;
    let (mut step, _work) =
        SortMerge::<BlockOutput>::with_target_batch_count(SortOrder::Coordinate, limit, target)
            .with_spill_supply(&supply)
            .with_sort_stats(true)
            .with_range_sizing_for_test(RANGE_SIZING)
            .with_job_edge_limit(one_ticket)
            .with_pooled_gather(1, None);
    let mut probe = StepProbe::new(&step);
    let jobs = probe.take_output::<MergeTicket<BlockOutput>>(1);
    let blocks = probe.take_output::<DecompressedBlock>(0);
    for event in inputs.events(&slots) {
        probe.push_input(event);
    }
    probe.close_input();

    let mut emitted = Vec::new();
    let (mut max_posted_ahead, mut ranges, mut finished) = (0usize, 0usize, false);
    for _ in 0..1_000_000 {
        let outcome = probe.try_run(&mut step).expect("try_run");
        if let SortMergeState::RangeMode(rm) = &step.state {
            let (live, _, posted, count) = rm.window_state_for_test();
            max_posted_ahead = max_posted_ahead.max(live);
            ranges = count;
            assert!(posted <= count);
        }
        emitted.extend(std::iter::from_fn(|| blocks.pop()));
        if outcome == StepOutcome::Finished {
            finished = true;
            break;
        }
    }
    assert!(finished, "range mode wedged behind the unpopped ticket edge");
    assert!(ranges > 4, "{ranges} ranges");
    // The window at one phase-2 thread is the floor, 4.
    assert_eq!(max_posted_ahead, 4, "posting stopped short of the window");
    let unpopped = std::iter::from_fn(|| jobs.pop()).count();
    assert!(unpopped <= 1, "the edge holds one ticket ({unpopped})");
    let bodies: Vec<Vec<u8>> = emitted.iter().flat_map(|b| b.ordinal_and_records().1).collect();
    assert_eq!(bodies, inputs.expected_bodies());
    let lines = logs.lines();
    let thread = lines
        .iter()
        .find(|l| l.starts_with("Range merge thread: busy "))
        .expect("the merge thread line");
    assert!(thread.contains("ticket edge (self-check) 0.0 s"), "{thread}");
}

/// At one phase-2 thread with a one-ticket edge, in a real pipeline (the
/// pool pops tickets only between compressions), range mode finishes under
/// the watchdog with the serial merge's output.
#[test]
fn range_mode_is_live_at_one_thread_with_a_tiny_ticket_edge() {
    let inputs = RangeInputs::<RawCoordinateKey>::generate(3, 1200, 0, true);
    let dir = tempfile::tempdir().unwrap();
    let slots = inputs.write_runs(dir.path());
    let items = run_merge::<BlockOutput>(MergeRun {
        events: inputs.events(&slots),
        supply: supply_with(range_mode(1)),
        order: SortOrder::Coordinate,
        threads: 1,
        sizing: Some(RANGE_SIZING),
        batch: RANGE_BATCH,
        sort_stats: false,
        claim: None,
        job_edge: Some(std::mem::size_of::<MergeTicket<BlockOutput>>() as u64),
    })
    .expect("merge");
    let bodies: Vec<Vec<u8>> = items.iter().flat_map(|(_, r)| r.clone()).collect();
    assert_eq!(bodies, inputs.expected_bodies());
    let ordinals: Vec<u64> = items.iter().map(|(o, _)| *o).collect();
    assert_eq!(ordinals, (0..items.len() as u64).collect::<Vec<_>>());
}

/// When workers merge ranges (the `pool` claim policy, where the merge
/// thread fuses none), later ranges build their output in earlier ranges'
/// recycled buffers, and the output is still the oracle's.
#[test]
fn worker_ranges_reuse_recycled_buffers() {
    let logs = capture_logs();
    let inputs = RangeInputs::<RawCoordinateKey>::generate(3, 1200, 0, true);
    let dir = tempfile::tempdir().unwrap();
    let slots = inputs.write_runs(dir.path());
    let items = run_merge::<BlockOutput>(MergeRun {
        events: inputs.events(&slots),
        supply: supply_with(range_mode(4)),
        order: SortOrder::Coordinate,
        threads: 4,
        sizing: Some(RANGE_SIZING),
        batch: RANGE_BATCH,
        sort_stats: true,
        claim: Some(super::claim::ClaimPolicy::Pool),
        job_edge: None,
    })
    .expect("merge");
    let bodies: Vec<Vec<u8>> = items.iter().flat_map(|(_, r)| r.clone()).collect();
    assert_eq!(bodies, inputs.expected_bodies());
    let lines = logs.lines();
    let buffers = lines
        .iter()
        .find(|l| l.starts_with("Range output buffers: "))
        .expect("the output buffers line");
    let reused: usize = buffers["Range output buffers: ".len()..]
        .split_once(' ')
        .and_then(|(n, _)| n.parse().ok())
        .unwrap_or_else(|| panic!("no reuse count in {buffers:?}"));
    assert!(reused > 0, "no range reused a buffer: {buffers:?}");
}
