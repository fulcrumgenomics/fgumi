//! `ConsolidateRuns` alone (framing, ordinals, per-call budgets) and the whole
//! consolidation loop through [`KChain`]: `SpillWrite` → planner →
//! `PreadSpillSlices` → `SortSpillDecompress` ⇢ slots; `ConsolidateRuns` →
//! `SpillBlockCompressConsolidation` → `WriteConsolidatedRun`. Oracles are
//! computed without the code under test (`merged_reference`, the gather's
//! `frame_one_block`, stable sorts).

use std::io::{self, Write as _};
use std::ops::Range;
use std::path::PathBuf;
use std::sync::Arc;

use fgumi_bam_io::pread::SliceBufferPool;
use fgumi_pipeline_core::step::StepOutcome;
use fgumi_pipeline_core::testing::StepProbe;
use fgumi_pipeline_core::{CapPhase, PhaseCap, PhasedCap, StepCap};
use fgumi_sort::{
    CONSOLIDATION_SLOT_IDS_START, MergeDemand, RawCoordinateKey, RawSortKey, SpillBlockCompressor,
    SpillCodec, SpillKeyKind, TemplateKey, open_spill_slot, spill_magic, spill_trailer,
};
use rstest::rstest;
use tempfile::TempDir;

use super::ConsolidateRuns;
use crate::pread::{PreadSlices, ReadRequest, ReadSlice, ReadTarget, SpillClass};
use crate::sort::protocol::{ConsolidationJob, SpillBlockEvent};
use crate::sort::spill_block_compress::SpillBlockCompress;
use crate::sort::spill_decompress::SortSpillDecompress;
use crate::sort::spill_fixtures::{
    FrameShape, KChain, KOpts, coordinate_body, embedded_frame, framed_template_records,
    keyed_frame, keyed_run, loaded_slot, merged_reference, read_frames, template_key,
    write_frames_run, write_run_file,
};
use crate::sort::spill_gather::frame_one_block;
use crate::sort::supply_ledger::SpillSupply;

/// One emitted `Block`: its job, ordinal, last-in-file flag and raw bytes.
struct Emitted {
    job: u32,
    ordinal: u64,
    last: bool,
    bytes: Vec<u8>,
}

/// Run `step` over `jobs` to `Finished` through a `StepProbe` and return the
/// emitted blocks in order. Checks on the way: `records_ingested_so_far == 0`;
/// exactly one `is_last_in_file` per job, on that job's final block.
fn drive_consolidate(mut step: ConsolidateRuns, jobs: Vec<ConsolidationJob>) -> Vec<Emitted> {
    let n_jobs = jobs.len();
    let mut probe = StepProbe::new(&step);
    let out = probe.take_output::<SpillBlockEvent>(0);
    for j in jobs {
        probe.push_input(j);
    }
    probe.close_input();
    let mut events = Vec::new();
    loop {
        let o = probe.try_run(&mut step).unwrap();
        while let Some(ev) = fgumi_pipeline_core::step::InputHandle::pop(&out) {
            events.push(ev);
        }
        if o == StepOutcome::Finished {
            break;
        }
    }
    let mut emitted = Vec::new();
    let mut lasts = Vec::new();
    for ev in events {
        let SpillBlockEvent::Block {
            ordinal,
            file_id,
            is_last_in_file,
            bytes,
            records_ingested_so_far,
            ..
        } = ev
        else {
            panic!("only blocks")
        };
        assert_eq!(records_ingested_so_far, 0);
        assert!(lasts.last() != Some(&file_id), "block after job {file_id}'s last block");
        if is_last_in_file {
            lasts.push(file_id);
        }
        emitted.push(Emitted { job: file_id, ordinal, last: is_last_in_file, bytes });
    }
    assert_eq!(lasts.len(), n_jobs, "one last block per job");
    emitted
}

/// Review S4 (r1): consolidation output frames exactly like the gather. Oracle:
/// the production `frame_one_block` over the same sorted records. Catches: a
/// different cut rule (split records, an off-by-one at the threshold, no
/// rollback of an oversized record).
#[rstest]
fn try_step_framed_frames_like_the_gather(
    #[values(FrameShape::Coordinate, FrameShape::TemplateK40)] shape: FrameShape,
    #[values(1usize, 100, 65_536)] block_size: usize,
) {
    let chunk = shape.sorted_chunk(500);
    let mut expected = Vec::new();
    let mut i = 0;
    while i < chunk.len() {
        let mut block = Vec::new();
        i = frame_one_block(&chunk, i, block_size, &mut block).unwrap();
        expected.push(block);
    }
    // Two interleaved runs (even / odd records), merged back by the driver.
    let runs: Vec<Vec<u8>> = (0..2)
        .map(|parity| {
            let mut raw = Vec::new();
            for j in (parity..chunk.len()).step_by(2) {
                chunk.frame_record_into(j, &mut raw).unwrap();
            }
            raw
        })
        .collect();
    let job = ConsolidationJob {
        job_id: 3,
        slots: vec![loaded_slot(0, 0, runs[0].clone()), loaded_slot(1, 1, runs[1].clone())],
        key_kind: shape.key_kind(),
        input_paths: vec![PathBuf::from("run-0"), PathBuf::from("run-1")],
    };
    let got = drive_consolidate(
        ConsolidateRuns::new(block_size, u64::MAX).with_merge_demand(Arc::new(MergeDemand::new())),
        vec![job],
    );
    // A job that ends on a full block emits an empty last block; the gather has
    // no such block, so drop it before comparing.
    let blocks: Vec<Vec<u8>> =
        got.into_iter().filter(|e| !(e.last && e.bytes.is_empty())).map(|e| e.bytes).collect();
    assert_eq!(blocks, expected);
}

/// Review I1 (r1): ordinals are dense across jobs, so the ordered SBC edge never
/// waits on an ordinal that will not come. Catches: `next_ordinal` reset per job.
#[test]
fn consecutive_jobs_get_dense_ordinals() {
    let job = |id: u32| ConsolidationJob {
        job_id: id,
        slots: vec![
            loaded_slot(0, 10 * id, framed_template_records(0..50)),
            loaded_slot(1, 10 * id + 1, framed_template_records(25..75)),
        ],
        key_kind: SpillKeyKind::TemplateK40,
        input_paths: vec![PathBuf::from("a"), PathBuf::from("b")],
    };
    let got = drive_consolidate(
        ConsolidateRuns::new(256, u64::MAX).with_merge_demand(Arc::new(MergeDemand::new())),
        vec![job(1), job(2)],
    );
    let ordinals: Vec<u64> = got.iter().map(|e| e.ordinal).collect();
    assert_eq!(
        ordinals,
        (0..got.len() as u64).collect::<Vec<_>>(),
        "dense from 0 across both jobs"
    );
    assert!(got.iter().any(|e| e.job == 1) && got.iter().any(|e| e.job == 2));
    assert!(got.windows(2).all(|w| w[0].job <= w[1].job), "job 1's blocks precede job 2's");
}

/// Overlapping key ranges, so merges interleave.
fn runs(n: u32) -> Vec<(u32, Range<u64>)> {
    (0..n).map(|id| (id, u64::from(id) * 3..u64::from(id) * 3 + 10)).collect()
}

fn feed_and_settle(k: &mut KChain, runs: &[(u32, Range<u64>)]) {
    for (id, keys) in runs {
        k.feed(keyed_run(SpillCodec::Zstd, *id, keys.clone()));
        k.settle().unwrap();
    }
    k.feed(SpillBlockEvent::AllAnnounced {
        ordinal: 0,
        slot_count: u32::try_from(runs.len()).unwrap(),
        memory_chunk_count: 0,
        total_records: 0,
    });
    k.close();
    k.settle().unwrap();
}

/// Review Focus 1 (r1 B1). Limit 2 → fan-in clamp(2/2, 2, 32) = 2, so every run
/// after the first triggers one pass and the merged run keeps `file_id` 0: the
/// cascading reuse. Passes in closed form: ⌈(N − L + 1)/(f − 1)⌉ = N − 1.
/// Catches: file_id-keyed supply maps (two slots for file 0 collide), retired
/// slots kept alive (descriptor leak), lost or duplicated records.
#[test]
fn consolidated_run_reuses_no_supply_state() {
    const N: u32 = 5;
    const L: u32 = 2;
    const F: u32 = 2;
    let mut k = KChain::new(KOpts { limit: L as usize, ..KOpts::default() });
    let input = runs(N);
    for (id, keys) in &input {
        k.feed(keyed_run(SpillCodec::Zstd, *id, keys.clone()));
        k.settle().unwrap();
        let live = k.live_slot_ids();
        for (job_id, at) in &k.retired_at {
            let (_, slot_ids, _, sources) = k.jobs.iter().find(|j| j.0 == *job_id).unwrap();
            for (sid, src) in slot_ids.iter().zip(sources) {
                assert_eq!(
                    Arc::strong_count(src),
                    1,
                    "job {job_id} slot {sid}: only the test's clone may remain"
                );
                assert!(!live.contains(sid), "retired slot {sid} still referenced: {live:?}");
            }
            assert!(
                k.requests[*at..].iter().all(|(s, _)| !slot_ids.contains(s)),
                "reads after retire of job {job_id}"
            );
        }
    }
    let passes = (N - L + 1).div_ceil(F - 1);
    assert!(passes >= 3);
    assert_eq!(u32::try_from(k.jobs.len()).unwrap(), passes);
    let all_slot_ids: Vec<u32> = k.jobs.iter().flat_map(|j| j.1.clone()).collect();
    let mut dedup = all_slot_ids.clone();
    dedup.sort_unstable();
    dedup.dedup();
    assert_eq!(dedup.len(), all_slot_ids.len(), "every opened slot got a fresh id");
    assert!(all_slot_ids.iter().all(|&s| s >= CONSOLIDATION_SLOT_IDS_START));
    k.feed(SpillBlockEvent::AllAnnounced {
        ordinal: 0,
        slot_count: N,
        memory_chunk_count: 0,
        total_records: 0,
    });
    k.close();
    k.settle().unwrap();
    assert_eq!(
        k.surviving_runs(),
        vec![merged_reference(&input)],
        "one surviving run = the stable merge of all inputs"
    );
}

/// R2 / review I5 (r1): consolidation reads only through the planner. When a job
/// is forwarded, a decoy run (same names, different records) is renamed over
/// each input path, so any path-based open reads the decoy and the record oracle
/// fails. (A plain unlink would also break `SpillWrite`'s own removal of the
/// absorbed inputs.) Also: each slot is read exactly once, `len − body_start`
/// bytes.
#[test]
fn consolidation_reads_go_through_the_planner() {
    let mut k = KChain::new(KOpts::default());
    let decoy_dir = TempDir::new().unwrap();
    let mut n = 0u32;
    k.on_job = Some(Box::new(move |job: &ConsolidationJob| {
        for p in &job.input_paths {
            let decoy = decoy_dir.path().join(format!("decoy-{n}"));
            n += 1;
            write_run_file(&decoy, SpillCodec::Zstd, 999, 1_000..1_010);
            std::fs::rename(&decoy, p).unwrap();
        }
    }));
    let input = runs(4);
    feed_and_settle(&mut k, &input);
    assert_eq!(k.surviving_runs(), vec![merged_reference(&input)]);
    for (_, slot_ids, _, sources) in &k.jobs {
        for (sid, src) in slot_ids.iter().zip(sources) {
            let read: u64 =
                k.requests.iter().filter(|(s, _)| s == sid).map(|(_, l)| u64::from(*l)).sum();
            let body_start = spill_magic(SpillCodec::Zstd).len() as u64;
            // The open descriptor's length (fstat): the original file, not the
            // decoy now at its path.
            assert_eq!(read, src.metadata().unwrap().len() - body_start, "slot {sid}");
        }
    }
}

/// Review r2 I6 (and Design issue 2): the phase switch happens at `AllAnnounced`
/// and never with a job in flight. Catches: an early switch (consolidation
/// supply escaping `--sort-threads`), and a switch while a job runs.
#[test]
fn phased_cap_switches_only_after_retire() {
    let phased = PhasedCap::new(
        Some(PhaseCap::new("sort-phase1", 2)),
        Some(PhaseCap::new("sort-phase2", 2)),
    )
    .unwrap();
    let mut k = KChain::new(KOpts { phased: Some(Arc::clone(&phased)), ..KOpts::default() });
    for (id, keys) in runs(3) {
        k.feed(keyed_run(SpillCodec::Zstd, id, keys));
        k.settle().unwrap();
        assert_eq!(phased.phase(), CapPhase::One, "after run {id}");
    }
    k.feed(SpillBlockEvent::AllAnnounced {
        ordinal: 0,
        slot_count: 3,
        memory_chunk_count: 0,
        total_records: 0,
    });
    k.close();
    k.settle().unwrap();
    assert_eq!(phased.phase(), CapPhase::Two);
}

/// R2: the three planner-fed steps admit through the active arm. Catches:
/// consolidation supply taking phase-2 permits (escaping `--sort-threads`), and
/// an admit that pops before refusing.
#[test]
fn supply_steps_take_the_phase1_cap_during_consolidation() {
    let (p1, p2) = (PhaseCap::new("sort-phase1", 1), PhaseCap::new("sort-phase2", 4));
    let phased = PhasedCap::new(Some(Arc::clone(&p1)), Some(Arc::clone(&p2))).unwrap();
    let cap = StepCap::Phased(Arc::clone(&phased));
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("r0");
    write_run_file(&path, SpillCodec::Zstd, 0, 0..10);
    let slot = open_spill_slot(&path, 0).unwrap();
    let body_len = u32::try_from(slot.len() - slot.body_start()).unwrap();
    // As the planner does when it issues the read the slice answers.
    slot.bp_note_issued(u64::from(body_len));
    let request = || ReadRequest {
        ordinal: 0,
        stream: 0,
        seq: 0,
        source: Arc::clone(slot.source()) as Arc<dyn fgumi_bam_io::pread::PositionalSource>,
        offset: slot.body_start(),
        len: body_len,
        last: true,
        target: ReadTarget::Spill { slot: Arc::clone(&slot), class: SpillClass::ConsolidationCold },
    };
    let pool = SliceBufferPool::new(2);
    let body =
        std::fs::read(&path).unwrap()[usize::try_from(slot.body_start()).unwrap()..].to_vec();
    let slice = || ReadSlice {
        ordinal: 0,
        stream: 0,
        seq: 0,
        offset: slot.body_start(),
        bytes: pool.lease(body.clone()),
        last: true,
        target: ReadTarget::Spill { slot: Arc::clone(&slot), class: SpillClass::ConsolidationCold },
    };
    let block = || SpillBlockEvent::Block {
        ordinal: 0,
        file_id: 1,
        key_kind: SpillKeyKind::TemplateK40,
        is_last_in_file: true,
        records_ingested_so_far: 0,
        bytes: framed_template_records(0..10),
    };

    let mut pread = PreadSlices::spill(SliceBufferPool::new(2), 1 << 20).with_step_cap(cap.clone());
    let mut dec = SortSpillDecompress::new(1 << 20, &SpillSupply::new()).with_step_cap(cap.clone());
    let mut sbc = SpillBlockCompress::new(SpillCodec::Zstd, 1, 1 << 20)
        .named("SpillBlockCompressConsolidation")
        .with_step_cap(cap);
    let (pp, dp, sp) = (StepProbe::new(&pread), StepProbe::new(&dec), StepProbe::new(&sbc));
    pp.push_input(request());
    dp.push_input(slice());
    sp.push_input(block());

    let held = p1.try_acquire().expect("the phase-1 cap's only permit");
    assert_eq!(pp.try_run(&mut pread).unwrap(), StepOutcome::Capped);
    assert_eq!(dp.try_run(&mut dec).unwrap(), StepOutcome::Capped);
    assert_eq!(sp.try_run(&mut sbc).unwrap(), StepOutcome::Capped);
    assert!(
        !pp.input_is_empty() && !dp.input_is_empty() && !sp.input_is_empty(),
        "refusals pop nothing"
    );
    assert_eq!(p2.refused(), 0, "phase 2 was never consulted");

    phased.switch_to(CapPhase::Two);
    assert_eq!(pp.try_run(&mut pread).unwrap(), StepOutcome::Progress);
    assert_eq!(dp.try_run(&mut dec).unwrap(), StepOutcome::Progress);
    assert_eq!(sp.try_run(&mut sbc).unwrap(), StepOutcome::Progress);
    drop(held);
}

// ── The retired `run_consolidate` tests, through the real tail ──────────────

/// Stable merge in run order; ties resolved by run. Catches: a tie broken by
/// slot id or by arrival order instead of run order.
#[rstest]
fn merge_is_the_stable_merge_of_its_inputs_in_run_order(
    #[values(SpillCodec::Zstd, SpillCodec::Bgzf)] codec: SpillCodec,
) {
    let mut k = KChain::new(KOpts { codec, block_size: 29, ..KOpts::default() });
    let runs: Vec<Vec<(TemplateKey, Vec<u8>)>> = (0..3u8)
        .map(|run| {
            (0..40u64)
                .map(|i| {
                    (
                        template_key(i / 2 + u64::from(run)),
                        vec![run, u8::try_from(i).unwrap(), 0xAB],
                    )
                })
                .collect()
        })
        .collect();
    let inputs: Vec<PathBuf> = runs
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let p = k.dir.path().join(format!("chunk_{i}"));
            write_frames_run(
                &p,
                codec,
                &r.iter().map(|(key, b)| keyed_frame(key, b)).collect::<Vec<_>>(),
                23,
            );
            p
        })
        .collect();
    let output = k.dir.path().join("merged");
    k.start_job(1, SpillKeyKind::TemplateK40, &inputs, &output);
    k.settle().unwrap();
    let mut expected: Vec<(TemplateKey, Vec<u8>)> = runs.concat();
    expected.sort_by(|a, b| a.0.cmp(&b.0)); // stable
    assert_eq!(
        read_frames(&output, codec, 40),
        expected.iter().map(|(key, b)| keyed_frame(key, b)).collect::<Vec<_>>()
    );
    assert_eq!(k.registry.take_completed().map(|c| c.0), Some(1));
}

/// Embedded keys (coordinate): re-extracted from the body. Catches: a
/// body-keyed slot parsed as if it had a key prefix, or compared on garbage.
#[test]
fn merge_of_embedded_keys_re_extracts_them_from_the_body() {
    let mut k = KChain::new(KOpts::default());
    let runs: Vec<Vec<Vec<u8>>> = (0..4u32)
        .map(|run| {
            (0..25u32)
                .map(|i| coordinate_body(0, i32::try_from(i * 2 + run % 2).unwrap(), run * 100 + i))
                .collect()
        })
        .collect();
    let inputs: Vec<PathBuf> = runs
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let p = k.dir.path().join(format!("chunk_{i}"));
            write_frames_run(
                &p,
                SpillCodec::Zstd,
                &r.iter().map(|b| embedded_frame(b)).collect::<Vec<_>>(),
                50,
            );
            p
        })
        .collect();
    let output = k.dir.path().join("merged");
    k.start_job(1, SpillKeyKind::Coordinate, &inputs, &output);
    k.settle().unwrap();
    let mut expected: Vec<Vec<u8>> = runs.concat();
    expected.sort_by_key(|b| RawCoordinateKey::extract_from_record(b)); // stable
    assert_eq!(
        read_frames(&output, SpillCodec::Zstd, 0),
        expected.iter().map(|b| embedded_frame(b)).collect::<Vec<_>>()
    );
}

/// Review Focus (items §5.8): empty inputs are skipped, and an all-empty job
/// writes an empty run (magic + trailer) and still completes. Catches: an
/// all-empty job that never emits its last block (`SpillWrite` waits forever).
#[test]
fn empty_inputs_are_skipped_and_an_all_empty_job_writes_an_empty_run() {
    let mut k = KChain::new(KOpts { codec: SpillCodec::Bgzf, block_size: 64, ..KOpts::default() });
    let path = |n: &str| k.dir.path().join(n);
    let records = vec![keyed_frame(&template_key(1), &[1]), keyed_frame(&template_key(2), &[2])];
    write_frames_run(&path("a"), SpillCodec::Bgzf, &[], 64);
    write_frames_run(&path("b"), SpillCodec::Bgzf, &records, 64);
    write_frames_run(&path("c"), SpillCodec::Bgzf, &[], 64);
    let (a, b, c) = (path("a"), path("b"), path("c"));
    let (merged, empty) = (path("merged"), path("empty"));
    k.start_job(1, SpillKeyKind::TemplateK40, &[a.clone(), b, c.clone()], &merged);
    k.settle().unwrap();
    assert_eq!(read_frames(&merged, SpillCodec::Bgzf, 40), records);
    k.start_job(2, SpillKeyKind::TemplateK40, &[a, c], &empty);
    k.settle().unwrap();
    assert!(read_frames(&empty, SpillCodec::Bgzf, 40).is_empty());
    let bytes = std::fs::read(&empty).unwrap();
    assert_eq!(bytes, [spill_magic(SpillCodec::Bgzf), spill_trailer(SpillCodec::Bgzf)].concat());
    let done: Vec<u32> = std::iter::from_fn(|| k.registry.take_completed().map(|c| c.0)).collect();
    assert_eq!(done, [1, 2]);
}

/// Review S3 (r1): a truncated input names its path. Catches: an error that
/// names a slot id or nothing.
#[test]
fn a_run_truncated_mid_record_is_an_unexpected_eof_naming_the_input() {
    let mut k = KChain::new(KOpts::default());
    let path = k.dir.path().join("chunk");
    let mut raw = keyed_frame(&template_key(1), &[1, 2, 3]);
    raw.extend(keyed_frame(&template_key(2), &[4, 5, 6, 7]));
    raw.truncate(raw.len() - 2);
    let mut f = std::fs::File::create(&path).unwrap();
    f.write_all(spill_magic(SpillCodec::Zstd)).unwrap();
    f.write_all(
        &SpillBlockCompressor::new(SpillCodec::Zstd, 1).unwrap().compress_block(&raw).unwrap(),
    )
    .unwrap();
    drop(f);
    let merged = k.dir.path().join("merged");
    k.start_job(1, SpillKeyKind::TemplateK40, std::slice::from_ref(&path), &merged);
    let err = k.settle().expect_err("the second record is truncated");
    assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    assert_eq!(err.to_string(), format!("consolidation input {} is truncated", path.display()));
}

/// The output must not already exist (`create_new`). Catches: truncating a live
/// input named as the output.
#[test]
fn the_output_path_must_not_already_exist() {
    let mut k = KChain::new(KOpts::default());
    let input = k.dir.path().join("chunk");
    write_frames_run(&input, SpillCodec::Zstd, &[keyed_frame(&template_key(1), &[1])], 64);
    let before = std::fs::read(&input).unwrap();
    k.start_job(1, SpillKeyKind::TemplateK40, std::slice::from_ref(&input), &input);
    let err = k.settle().expect_err("creating the output over a live input must fail");
    assert_eq!(err.kind(), io::ErrorKind::AlreadyExists);
    assert_eq!(std::fs::read(&input).unwrap(), before, "the input is untouched");
}

/// Per-call budgets bound one `try_run` (records and bytes independently).
/// Catches: a call that drains a whole job (a long `sort-io` driver stall).
#[test]
fn per_call_budgets_bound_one_try_run() {
    let raw: Vec<u8> =
        (0..10u64).flat_map(|i| keyed_frame(&template_key(i), &[7u8; 100])).collect();
    let job = || ConsolidationJob {
        job_id: 1,
        slots: vec![loaded_slot(0, 0, raw.clone())],
        key_kind: SpillKeyKind::TemplateK40,
        input_paths: vec![PathBuf::from("a")],
    };
    // 250 bytes of 100-byte bodies: the third record crosses the byte budget.
    let mut step = ConsolidateRuns::new(1 << 20, u64::MAX)
        .with_merge_demand(Arc::new(MergeDemand::new()))
        .with_per_call(usize::MAX, 250);
    let probe = StepProbe::new(&step);
    probe.push_input(job());
    assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Progress); // pops the job and merges
    assert_eq!(step.records_merged_for_test(), 3);
    // A record budget of 1 merges exactly one record per call.
    let mut step = ConsolidateRuns::new(1 << 20, u64::MAX)
        .with_merge_demand(Arc::new(MergeDemand::new()))
        .with_per_call(1, usize::MAX);
    let probe = StepProbe::new(&step);
    probe.push_input(job());
    for want in 1..=3 {
        assert_eq!(probe.try_run(&mut step).unwrap(), StepOutcome::Progress);
        assert_eq!(step.records_merged_for_test(), want);
    }
}

/// The retired `SpillWrite` test, now through the real tail: limit 4, 25
/// overlapping runs. Catches: live runs above the limit, a lost or duplicated
/// record, unsorted surviving runs, absorbed inputs left on disk.
#[test]
fn consolidation_keeps_live_runs_under_the_limit_and_preserves_every_record() {
    const LIMIT: usize = 4;
    const RUNS: u32 = 25;
    let mut k = KChain::new(KOpts { limit: LIMIT, ..KOpts::default() });
    let input: Vec<(u32, Range<u64>)> =
        (0..RUNS).map(|id| (id, u64::from(id) * 3..u64::from(id) * 3 + 10)).collect();
    for (id, keys) in &input {
        k.feed(keyed_run(SpillCodec::Zstd, *id, keys.clone()));
        k.settle().unwrap();
        assert!(k.live_runs() < LIMIT, "live runs {} after run {id}", k.live_runs());
    }
    assert!(!k.jobs.is_empty(), "25 runs under a limit of 4 must consolidate");
    k.feed(SpillBlockEvent::AllAnnounced {
        ordinal: 0,
        slot_count: RUNS,
        memory_chunk_count: 0,
        total_records: 0,
    });
    k.close();
    k.settle().unwrap();
    let survivors = k.surviving_runs();
    assert!(survivors.len() < LIMIT);
    for run in &survivors {
        assert!(
            run.windows(2).all(|w| TemplateKey::read_from(&mut &w[0].0[..]).unwrap()
                <= TemplateKey::read_from(&mut &w[1].0[..]).unwrap()),
            "a surviving run is unsorted"
        );
    }
    let mut got: Vec<_> = survivors.concat();
    let mut want = merged_reference(&input);
    got.sort_by(|a, b| a.1.cmp(&b.1));
    want.sort_by(|a, b| a.1.cmp(&b.1));
    assert_eq!(got, want, "records lost, duplicated or altered");
    assert_eq!(
        std::fs::read_dir(k.dir.path()).unwrap().count(),
        survivors.len(),
        "absorbed inputs deleted"
    );
}
