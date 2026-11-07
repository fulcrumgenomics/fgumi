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
    posted: Vec<crate::sort::merge_job::MergeJob<BlockOutput>>,
) {
    use crate::sort::merge_job::MERGE_SLICE;
    use fgumi_pipeline_core::{InputHandle, Job, JobStep};

    let mut done: Vec<_> = posted
        .into_iter()
        .map(|mut job| {
            let ordinal = job.ordinal();
            match job.run_slice(MERGE_SLICE).expect("gather") {
                JobStep::Done(out) => (ordinal, out),
                JobStep::Continue => panic!("a gather is one slice"),
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

/// Gather mode posts one `MergeJob::Gather` per planned block on branch 1, with
/// ordinals `0..k` in order and never more than `window` (`fanout_width * 4`)
/// ahead of the results taken back; it parks (`NoProgress`) while the next
/// ordinal is missing, and results delivered out of order leave branch 0 in
/// ordinal order, byte-identical to gathering `plan_fast_path_blocks`' plan
/// serially.
#[test]
fn gather_mode_posts_jobs_within_window() {
    use fgumi_pipeline_core::testing::StepProbe;
    use fgumi_pipeline_core::{InputHandle, Job, OrderedInbox};

    use crate::sort::merge_job::{MergeJob, MergeResult};

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
    let jobs = probe.take_output::<MergeJob<BlockOutput>>(1);
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

    let mut posted: Vec<MergeJob<BlockOutput>> = Vec::new();
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
    fgumi_pipeline_core::BranchInputHandle<crate::sort::merge_job::MergeJob<BlockOutput>>,
) {
    use fgumi_pipeline_core::OrderedInbox;
    use fgumi_pipeline_core::testing::StepProbe;

    use crate::sort::merge_job::{MergeJob, MergeResult};

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
    let jobs = probe.take_output::<MergeJob<BlockOutput>>(1);
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
    let per_job = std::mem::size_of::<crate::sort::merge_job::MergeJob<RecordBatchOutput>>();
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
