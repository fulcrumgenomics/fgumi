//! Unit tests of the merge driver: the slot-backed `MergeDriver::from_slots`
//! merge, the kind dispatch, body-keyed BAM slots and slot ingest.

use super::*;
use crate::template_key::{KeyTypesSpec, extract_template_key_inline};
use fgumi_bam_io::create_raw_bam_reader;
use fgumi_raw_bam::SamTag;
use noodles::sam::Header;
use std::io::Read;
use std::path::{Path, PathBuf};

mod from_slots_merge_tests {
    //! Direct unit tests for the slot-backed `MergeDriver::from_slots` merge
    //! driver: cross-block record parsing (incl. records spanning many blocks),
    //! non-blocking stall/resume on empty-non-EOF slots, the embedded and
    //! non-embedded sort-key arms, and memory-chunk mixing. Ported verbatim
    //! from the issue-#330 source branch (the slimmed `from_slots` API is
    //! identical at this layer). These exercise the parser/state-machine paths
    //! that the end-to-end `three_step_chain_*` tests do not isolate.
    use super::*;
    use std::io::{Read, Write};
    use std::sync::Arc as StdArc;

    #[derive(Default, PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Debug)]
    struct TestKey(u64);

    impl RawSortKey for TestKey {
        type MergeCtx = ();

        type MergeHint = ();

        const SERIALIZED_SIZE: Option<usize> = Some(8);
        const EMBEDDED_IN_RECORD: bool = false;
        fn extract(_bam: &[u8], _ctx: &crate::keys::SortContext) -> Self {
            unimplemented!(
                "TestKey is only used in MergeDriver::from_slots tests, not for ingestion"
            )
        }
        fn write_to<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
            w.write_all(&self.0.to_le_bytes())
        }
        fn read_from<R: Read>(r: &mut R) -> std::io::Result<Self> {
            let mut buf = [0u8; 8];
            r.read_exact(&mut buf)?;
            Ok(TestKey(u64::from_le_bytes(buf)))
        }
    }

    /// Serialize a sequence of `(TestKey, record_bytes)` pairs in the spill
    /// file's wire format: `[key(8)][len(4)][record(len)]` per entry.
    fn serialize_records(records: &[(TestKey, Vec<u8>)]) -> Vec<u8> {
        let mut out = Vec::new();
        for (key, rec) in records {
            key.write_to(&mut out).unwrap();
            #[allow(clippy::cast_possible_truncation)]
            let len = rec.len() as u32;
            out.write_all(&len.to_le_bytes()).unwrap();
            out.write_all(rec).unwrap();
        }
        out
    }

    /// Build a populated `SortMergeSlot` with `file_id` whose decompressed
    /// queue contains `records` partitioned across `n_blocks` "blocks"
    /// (so the slot parser exercises cross-block reads). Sets `queue_eof`
    /// so the slot reports drained once consumer empties it.
    fn populated_slot(
        file_id: u32,
        records: &[(TestKey, Vec<u8>)],
        n_blocks: usize,
    ) -> StdArc<SortMergeSlot> {
        assert!(n_blocks >= 1, "n_blocks must be >= 1");
        let bytes = serialize_records(records);
        let block_size = bytes.len().div_ceil(n_blocks);
        let slot = StdArc::new(SortMergeSlot::for_test(file_id, crate::codec::SpillCodec::Bgzf));
        {
            let mut dec = slot.decompressed.lock().unwrap();
            for chunk in bytes.chunks(block_size) {
                dec.push_back(chunk.to_vec());
            }
            slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        }
        slot
    }

    /// Drive a `MergeDriver` to exhaustion via the non-blocking `try_step`,
    /// returning the emitted record bytes in emission order.
    ///
    /// Tests validate ordering by inspecting the human-readable record bytes
    /// (e.g. `b"AAA-1"`, `b"BBB-2"`) — the keys themselves are consumed by
    /// the parser before the record is exposed, so we have no direct access to
    /// them at this layer. These slots are pre-populated and EOF-marked before
    /// the merge, so `try_step` never returns `Stalled` here.
    fn drain_merge_driver_bytes<K: RawSortKey + Default + Send + 'static, R: ChunkRef<K>>(
        mut driver: MergeDriver<K, R>,
    ) -> Vec<Vec<u8>> {
        let mut out = Vec::new();
        loop {
            match driver.try_step().expect("try_step") {
                MergeStep::Produced(bytes) => out.push(bytes.to_vec()),
                MergeStep::Done => break,
                MergeStep::Stalled => panic!("unexpected Stalled on pre-populated EOF slots"),
            }
        }
        out
    }

    #[test]
    fn from_slots_merges_three_pre_populated_slots() {
        // Three slots, each with sorted records by TestKey. The merged output
        // should be globally sorted.
        let file0: Vec<(TestKey, Vec<u8>)> = vec![
            (TestKey(1), b"AAA-1".to_vec()),
            (TestKey(4), b"AAA-4".to_vec()),
            (TestKey(7), b"AAA-7".to_vec()),
        ];
        let file1: Vec<(TestKey, Vec<u8>)> = vec![
            (TestKey(2), b"BBB-2".to_vec()),
            (TestKey(5), b"BBB-5".to_vec()),
            (TestKey(8), b"BBB-8".to_vec()),
        ];
        let file2: Vec<(TestKey, Vec<u8>)> = vec![
            (TestKey(3), b"CCC-3".to_vec()),
            (TestKey(6), b"CCC-6".to_vec()),
            (TestKey(9), b"CCC-9".to_vec()),
        ];

        let slots = vec![
            populated_slot(0, &file0, 1),
            populated_slot(1, &file1, 1),
            populated_slot(2, &file2, 1),
        ];

        let driver = MergeDriver::<TestKey>::from_slots(slots, Vec::new(), 9, ());
        let emitted_bytes = drain_merge_driver_bytes(driver);
        let expected: Vec<Vec<u8>> = vec![
            b"AAA-1".to_vec(),
            b"BBB-2".to_vec(),
            b"CCC-3".to_vec(),
            b"AAA-4".to_vec(),
            b"BBB-5".to_vec(),
            b"CCC-6".to_vec(),
            b"AAA-7".to_vec(),
            b"BBB-8".to_vec(),
            b"CCC-9".to_vec(),
        ];
        assert_eq!(emitted_bytes, expected, "merged sequence not globally sorted by key");
    }

    #[test]
    fn from_slots_handles_records_spanning_block_boundaries() {
        // Same data, but each slot's bytes are split across 4 "decompressed
        // blocks" so the slot parser's read_exact + advance_to_next_block
        // path is exercised. The key+len header itself may straddle a block
        // boundary depending on byte counts.
        let records: Vec<(TestKey, Vec<u8>)> = (0u64..20)
            .map(|i| (TestKey(i), format!("record-{i:02}-with-some-padding").into_bytes()))
            .collect();

        let slots = vec![populated_slot(0, &records, 4)];
        let driver = MergeDriver::<TestKey>::from_slots(slots, Vec::new(), 20, ());
        let emitted_bytes = drain_merge_driver_bytes(driver);
        let expected: Vec<Vec<u8>> = records.into_iter().map(|(_, b)| b).collect();
        assert_eq!(emitted_bytes, expected, "cross-block parse corrupted record sequence");
    }

    #[test]
    fn from_slots_empty_sources_drain_to_zero_records() {
        let empty_slot = StdArc::new(SortMergeSlot::for_test(0, crate::codec::SpillCodec::Bgzf));
        // Mark drained without inserting any blocks.
        empty_slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        // Lazy priming: `from_slots` always yields a driver; the empty
        // determination happens on the first `try_step`, which returns `Done`.
        let driver = MergeDriver::<TestKey>::from_slots(vec![empty_slot], Vec::new(), 0, ());
        let emitted = drain_merge_driver_bytes(driver);
        assert!(emitted.is_empty(), "all-empty merge must emit zero records");
    }

    /// The length prefix is untrusted: it comes straight off a decompressed
    /// spill block. The slow path used to hand it to `out.reserve(len)` with no
    /// check, so a corrupt slot could request up to 4 GiB — and `reserve`
    /// answers an allocation failure by aborting the process, not by returning
    /// the truncation error the next loop iteration would have produced.
    ///
    /// The length here is over the bound but still allocatable, so this test
    /// discriminates on the *error message* rather than on surviving: before
    /// the bound existed it reserved happily and then failed as a truncated
    /// body, which is the wrong diagnosis for a corrupt length.
    #[test]
    fn implausible_record_length_is_rejected_before_reserving() {
        let bogus_len = u32::try_from(crate::inline::SORT_SEGMENT_SIZE).unwrap() + 1;
        let mut bytes = Vec::new();
        TestKey(1).write_to(&mut bytes).unwrap();
        bytes.write_all(&bogus_len.to_le_bytes()).unwrap();
        bytes.write_all(b"only-a-few-body-bytes").unwrap();

        let slot = StdArc::new(SortMergeSlot::for_test(5, crate::codec::SpillCodec::Bgzf));
        slot.decompressed.lock().unwrap().push_back(bytes);
        slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);

        let mut driver = MergeDriver::<TestKey>::from_slots(vec![slot], Vec::new(), 1, ());
        let err = driver.try_step().expect_err("an implausible length must surface as Err");
        assert!(
            err.to_string().contains("implausible record length") && err.to_string().contains('5'),
            "the error must name the bad length and the slot; got: {err}",
        );
    }

    /// An embedded key is read out of the record body itself, so an empty body
    /// cannot carry one. `RawCoordinateKey::extract_from_record` reads `ref_id`
    /// off the front, which on an empty body indexes out of bounds — a panic in
    /// the middle of a merge rather than a corrupt-slot error. A zero length is
    /// reachable from a corrupt length prefix: the body loop sees
    /// `out.len() >= 0` and breaks immediately with `out` still empty.
    #[test]
    fn embedded_key_on_an_empty_body_is_an_error_not_a_panic() {
        let err = slot_parse_key::<crate::keys::RawCoordinateKey>(&(), &mut (), &[], &[])
            .expect_err("an empty body cannot carry an embedded key");
        assert!(
            err.to_string().contains("empty and cannot carry an embedded sort key"),
            "the error must name the empty body; got: {err}",
        );
    }

    /// Truncation at EOF must be an error, not a stall — at **every** stage of
    /// the framer, not just the last one.
    ///
    /// Every other test either completes its record or stalls on a non-EOF slot,
    /// so the three `BlockLoad::Eof` arms in `slot_collect_pending` were
    /// unexercised. They are siblings: swapping any one for `WouldBlock` hangs
    /// the pipeline forever on a truncated spill while the suite stays green, so
    /// covering only the body arm left two of the three live.
    ///
    /// A serialized record is `[key(8)][len(4)][body(len)]`, so keeping 4 and 10
    /// bytes lands mid-key and mid-length; `None` means "all but the last 5",
    /// which lands mid-body.
    #[rstest::rstest]
    #[case::mid_key(Some(4), "truncated record key in slot 3")]
    #[case::mid_length(Some(10), "truncated record length in slot 3")]
    #[case::mid_body(None, "truncated record body in slot 3")]
    fn truncation_at_eof_is_an_error_not_a_stall_at_every_stage(
        #[case] keep: Option<usize>,
        #[case] expected: &str,
    ) {
        let serialized = serialize_records(&[(TestKey(1), b"a-body-that-is-cut-short".to_vec())]);
        let keep = keep.unwrap_or(serialized.len() - 5);
        assert!(keep < serialized.len(), "the case must actually truncate");

        let slot = StdArc::new(SortMergeSlot::for_test(3, crate::codec::SpillCodec::Bgzf));
        slot.decompressed.lock().unwrap().push_back(serialized[..keep].to_vec());
        slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);

        let mut driver = MergeDriver::<TestKey>::from_slots(vec![slot], Vec::new(), 1, ());
        let err = driver.try_step().expect_err("truncation must surface as Err");
        assert!(err.to_string().contains(expected), "expected `{expected}`, got: {err}");
    }

    /// `from_slots` sorts slots by `file_id` so equal-key ties resolve the same
    /// way the legacy chunk-files path resolves them — the `LoserTree` breaks
    /// ties by leaf index, so leaf order must follow `file_id` and not the
    /// caller's `Vec` order. Nothing pinned that, so the sort could vanish in a
    /// later cleanup and the divergence would be invisible except on inputs
    /// that contain ties. Pass the slots reversed and require `file_id` to win.
    #[test]
    fn key_ties_break_by_file_id_not_by_caller_push_order() {
        let file0 = populated_slot(0, &[(TestKey(1), b"from-file-0".to_vec())], 1);
        let file1 = populated_slot(1, &[(TestKey(1), b"from-file-1".to_vec())], 1);
        let driver = MergeDriver::<TestKey>::from_slots(
            vec![file1, file0], // reversed on purpose
            Vec::new(),
            2,
            (),
        );
        assert_eq!(
            drain_merge_driver_bytes(driver),
            vec![b"from-file-0".to_vec(), b"from-file-1".to_vec()],
            "equal-key ties must resolve by file_id, not by caller push order",
        );
    }

    /// Predictions: published only on a source switch, naming the runner-up's
    /// file, never a memory source. The expected switch count and hits are
    /// read off the output (each record names its source), independently of
    /// the driver.
    #[test]
    fn predictions_skip_memory_sources_and_track_switches() {
        // Slot 0: keys 1,2 ; slot 1: keys 3,4 ; memory chunk: key 5.
        let s0 =
            populated_slot(0, &[(TestKey(1), b"src0".to_vec()), (TestKey(2), b"src0".to_vec())], 1);
        let s1 =
            populated_slot(1, &[(TestKey(3), b"src1".to_vec()), (TestKey(4), b"src1".to_vec())], 1);
        let mem =
            crate::inline::InMemoryChunk::from_owned_records(vec![(TestKey(5), b"mem".to_vec())]);
        let demand = StdArc::new(crate::MergeDemand::new());
        let mut d = MergeDriver::<TestKey>::from_slots(vec![s0, s1], vec![mem], 5, ())
            .with_demand(StdArc::clone(&demand));
        let mut seen: Vec<(Vec<u8>, Option<u32>)> = Vec::new();
        loop {
            match d.try_step().unwrap() {
                MergeStep::Produced(b) => seen.push((b.to_vec(), demand.predicted())),
                MergeStep::Done => break,
                MergeStep::Stalled => unreachable!("pre-populated EOF slots never stall"),
            }
        }
        let labels: Vec<&[u8]> = seen.iter().map(|(b, _)| b.as_slice()).collect();
        let expected: [&[u8]; 5] = [b"src0", b"src0", b"src1", b"src1", b"mem"];
        assert_eq!(labels, expected);
        // While draining slot 0 the runner-up is slot 1's file.
        assert_eq!(seen[0].1, Some(1));
        assert_eq!(seen[1].1, Some(1), "no switch, no new prediction");
        // While draining slot 1 the runner-up is the memory source: not a file.
        assert_eq!(seen[2].1, None, "a memory source is never predicted");
        assert_eq!(seen[4].1, None, "nothing remains behind the last source");
        let switches = labels.windows(2).filter(|w| w[0] != w[1]).count() as u64;
        let s = demand.snapshot();
        assert_eq!(s.predictions, switches, "each switch after the first scores one prediction");
        // s0→s1 was predicted (Some(1)), and so was the switch to memory
        // (None: a memory source).
        assert_eq!(s.prediction_hits, 2);
    }

    /// The frontier is the lowest undrained slot's file id: published when the
    /// demand is attached, advanced on each drain, and `None` once every slot
    /// is drained.
    #[test]
    fn frontier_advances_on_drain_and_clears_at_the_end() {
        let s0 = populated_slot(0, &[(TestKey(1), b"a".to_vec())], 1);
        let s1 = populated_slot(1, &[(TestKey(2), b"b".to_vec())], 1);
        let demand = StdArc::new(crate::MergeDemand::new());
        let mut d = MergeDriver::<TestKey>::from_slots(vec![s0, s1], vec![], 2, ())
            .with_demand(StdArc::clone(&demand));
        assert_eq!(demand.frontier(), Some(0), "published before the first step");
        assert!(matches!(d.try_step().unwrap(), MergeStep::Produced(_)));
        assert_eq!(demand.frontier(), Some(0));
        assert!(matches!(d.try_step().unwrap(), MergeStep::Produced(_)));
        assert_eq!(demand.frontier(), Some(1), "slot 0 drained on the refill");
        assert!(matches!(d.try_step().unwrap(), MergeStep::Done));
        assert_eq!(demand.frontier(), None);
    }

    /// A slot that is empty at EOF from the start is drained during priming:
    /// the frontier skips it before the first record.
    #[test]
    fn frontier_skips_a_slot_drained_while_priming() {
        let empty = StdArc::new(SortMergeSlot::for_test(0, crate::codec::SpillCodec::Bgzf));
        {
            let _g = empty.decompressed.lock().unwrap();
            empty.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        }
        let s1 = populated_slot(1, &[(TestKey(2), b"b".to_vec())], 1);
        let demand = StdArc::new(crate::MergeDemand::new());
        let mut d = MergeDriver::<TestKey>::from_slots(vec![empty, s1], vec![], 1, ())
            .with_demand(StdArc::clone(&demand));
        assert!(matches!(d.try_step().unwrap(), MergeStep::Produced(_)));
        assert_eq!(demand.frontier(), Some(1));
    }

    /// `stalled_slot` names the slot whose empty FIFO stalled the merge, in
    /// priming and in mid-merge refill, and is cleared by the next `Produced`.
    #[test]
    fn stalled_slot_names_the_awaited_slot() {
        let slot_a = populated_slot(0, &[(TestKey(1), b"A-1".to_vec())], 1);
        // Slot 1: empty, not EOF → priming stalls on it.
        let slot_b = StdArc::new(SortMergeSlot::for_test(1, crate::codec::SpillCodec::Bgzf));
        let mut driver = MergeDriver::<TestKey>::from_slots(
            vec![StdArc::clone(&slot_a), StdArc::clone(&slot_b)],
            Vec::new(),
            2,
            (),
        );
        assert!(matches!(driver.try_step().unwrap(), MergeStep::Stalled));
        assert_eq!(driver.stalled_slot().map(|s| s.file_id), Some(1), "priming stall");
        // Feed slot 1 one record, still open; the merge resumes.
        slot_b
            .decompressed
            .lock()
            .unwrap()
            .push_back(serialize_records(&[(TestKey(2), b"B-2".to_vec())]));
        assert!(matches!(driver.try_step().unwrap(), MergeStep::Produced(b) if b == b"A-1"));
        assert_eq!(driver.stalled_slot().map(|s| s.file_id), None, "cleared by Produced");
        assert!(matches!(driver.try_step().unwrap(), MergeStep::Produced(b) if b == b"B-2"));
        // Refilling slot 1's winner finds its FIFO empty and the slot open.
        assert!(matches!(driver.try_step().unwrap(), MergeStep::Stalled));
        assert_eq!(driver.stalled_slot().map(|s| s.file_id), Some(1), "mid-merge refill stall");
        {
            let _g = slot_b.decompressed.lock().unwrap();
            slot_b.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        }
        assert!(matches!(driver.try_step().unwrap(), MergeStep::Done));
    }

    /// The core BUG #3 regression: a slot whose decompressed queue is empty
    /// and `!queue_eof` must make `try_step` return `Stalled` (yield), NOT
    /// block the caller. The old blocking consumer (`block_ready.wait()`)
    /// would park here forever at a single worker. After the producer fills
    /// the slot and marks EOF, the merge resumes and drains to completion.
    #[test]
    fn from_slots_try_step_stalls_on_empty_non_eof_slot_then_resumes() {
        let slot = StdArc::new(SortMergeSlot::for_test(0, crate::codec::SpillCodec::Bgzf));
        // Empty queue, NOT eof — the producer is "still feeding".
        let mut driver =
            MergeDriver::<TestKey>::from_slots(vec![StdArc::clone(&slot)], Vec::new(), 0, ());

        // Priming cannot read the slot's first record → Stalled, not a hang.
        assert!(
            matches!(driver.try_step().expect("try_step"), MergeStep::Stalled),
            "empty non-eof slot must Stall, not block"
        );
        // Repeated calls keep stalling (idempotent yield).
        assert!(matches!(driver.try_step().expect("try_step"), MergeStep::Stalled));

        // Producer fills the slot and marks EOF.
        {
            let bytes = serialize_records(&[(TestKey(1), b"only-1".to_vec())]);
            let mut dec = slot.decompressed.lock().unwrap();
            dec.push_back(bytes);
            slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        }

        // The merge now resumes: emits the record, then is Done.
        match driver.try_step().expect("try_step") {
            MergeStep::Produced(bytes) => assert_eq!(bytes, b"only-1"),
            MergeStep::Stalled => panic!("should have resumed, not Stalled"),
            MergeStep::Done => panic!("should have produced the record, not Done"),
        }
        assert!(matches!(driver.try_step().expect("try_step"), MergeStep::Done));
    }

    /// `try_step` must also stall (not block) when a slot drains to empty
    /// **mid-merge** (after priming) while still `!queue_eof`, then resume
    /// once the producer pushes the next block and marks EOF.
    #[test]
    fn from_slots_try_step_stalls_mid_merge_then_resumes() {
        let slot = StdArc::new(SortMergeSlot::for_test(0, crate::codec::SpillCodec::Bgzf));
        // One record present, NOT eof (more "coming").
        {
            let bytes = serialize_records(&[(TestKey(1), b"rec-A".to_vec())]);
            slot.decompressed.lock().unwrap().push_back(bytes);
        }
        let mut driver =
            MergeDriver::<TestKey>::from_slots(vec![StdArc::clone(&slot)], Vec::new(), 0, ());

        // Prime + emit the first record.
        match driver.try_step().expect("try_step") {
            MergeStep::Produced(bytes) => assert_eq!(bytes, b"rec-A"),
            other => panic!("expected Produced(rec-A), got {other:?}"),
        }
        // Deferred refill finds the slot drained + not eof → Stalled.
        assert!(
            matches!(driver.try_step().expect("try_step"), MergeStep::Stalled),
            "drained mid-merge non-eof slot must Stall"
        );

        // Producer pushes the next record and marks EOF.
        {
            let bytes = serialize_records(&[(TestKey(2), b"rec-B".to_vec())]);
            let mut dec = slot.decompressed.lock().unwrap();
            dec.push_back(bytes);
            slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        }

        match driver.try_step().expect("try_step") {
            MergeStep::Produced(bytes) => assert_eq!(bytes, b"rec-B"),
            other => panic!("expected Produced(rec-B), got {other:?}"),
        }
        assert!(matches!(driver.try_step().expect("try_step"), MergeStep::Done));
    }

    /// Resumable-framer regression: feed a single **non-embedded** record one
    /// byte per block, asserting `try_step` `Stalled`s (retaining
    /// `parser.pending`) until the whole record is available, then reassembles
    /// the exact bytes. With 1-byte blocks the framer takes the slow path and
    /// resumes a `WouldBlock` across the key prefix, the length prefix, AND the
    /// body — the partial-record path no other test covers byte-for-byte.
    #[test]
    fn from_slots_try_step_resumes_record_fed_one_byte_at_a_time() {
        let rec_body = b"multi-stage-resumable-record-body".to_vec();
        let serialized = serialize_records(&[(TestKey(7), rec_body.clone())]);
        assert!(serialized.len() > 12, "must span the 8-byte key + 4-byte len header");

        let slot = StdArc::new(SortMergeSlot::for_test(0, crate::codec::SpillCodec::Bgzf));
        slot.decompressed.lock().unwrap().push_back(vec![serialized[0]]);

        let mut driver =
            MergeDriver::<TestKey>::from_slots(vec![StdArc::clone(&slot)], Vec::new(), 1, ());

        // Each step before the last byte consumes one byte into the pending
        // record and stalls (queue empty, not EOF).
        for (i, &b) in serialized.iter().enumerate().skip(1) {
            assert!(
                matches!(driver.try_step().expect("try_step"), MergeStep::Stalled),
                "expected Stalled with only {i} of {} bytes available",
                serialized.len(),
            );
            slot.decompressed.lock().unwrap().push_back(vec![b]);
        }

        // The step that consumes the final byte completes the record.
        match driver.try_step().expect("try_step") {
            MergeStep::Produced(bytes) => assert_eq!(bytes, rec_body.as_slice()),
            other => panic!("expected Produced(reassembled), got {other:?}"),
        }
        // No more bytes and not yet EOF → Stalled; EOF → Done.
        assert!(matches!(driver.try_step().expect("try_step"), MergeStep::Stalled));
        slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        assert!(matches!(driver.try_step().expect("try_step"), MergeStep::Done));
    }

    /// Same resumable-framer property for an **embedded**-key record split mid
    /// body across two blocks (the key lives in the body, so this exercises the
    /// embedded `extract_from_record` path resuming after a `WouldBlock`).
    #[test]
    fn from_slots_try_step_resumes_embedded_record_split_mid_body() {
        // Embedded record: [len:4][body], key = first 8 body bytes.
        let record = embedded_record(123, b"-embedded-tail-bytes");
        let serialized = serialize_embedded(std::slice::from_ref(&record));
        // Cut mid-body: 4-byte len + 10 body bytes in the first block.
        let split = 4 + 10;
        assert!(split < serialized.len() && split > 4 + 8, "split mid-body, past the key");
        let (head, tail) = serialized.split_at(split);

        let slot = StdArc::new(SortMergeSlot::for_test(0, crate::codec::SpillCodec::Bgzf));
        slot.decompressed.lock().unwrap().push_back(head.to_vec());

        let mut driver = MergeDriver::<TestEmbeddedKey>::from_slots(
            vec![StdArc::clone(&slot)],
            Vec::new(),
            1,
            (),
        );

        // Body incomplete in the first block → Stalled mid-record.
        assert!(matches!(driver.try_step().expect("try_step"), MergeStep::Stalled));

        // Deliver the rest + EOF; the record reassembles exactly.
        {
            let mut dec = slot.decompressed.lock().unwrap();
            dec.push_back(tail.to_vec());
            slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        }
        match driver.try_step().expect("try_step") {
            MergeStep::Produced(bytes) => assert_eq!(bytes, record.as_slice()),
            other => panic!("expected Produced(reassembled embedded), got {other:?}"),
        }
        assert!(matches!(driver.try_step().expect("try_step"), MergeStep::Done));
    }

    /// A merge of a slot with an arena-backed in-memory chunk: the driver takes
    /// each chunk record's key by index and emits `chunk.record_bytes(idx)` in
    /// place. Catches: a defect in that index, or in the `idx`/`len` pairing.
    #[test]
    fn from_slots_merges_slots_with_shared_memory_chunks() {
        let slot_records: Vec<(TestKey, Vec<u8>)> = vec![
            (TestKey(1), b"slot-1".to_vec()),
            (TestKey(3), b"slot-3".to_vec()),
            (TestKey(5), b"slot-5".to_vec()),
        ];
        let memory_chunk = crate::inline::InMemoryChunk::from_owned_records(vec![
            (TestKey(2), b"mem-2".to_vec()),
            (TestKey(4), b"mem-4".to_vec()),
            (TestKey(6), b"mem-6".to_vec()),
        ]);

        let slots = vec![populated_slot(0, &slot_records, 2)];
        let driver = MergeDriver::<TestKey>::from_slots(slots, vec![memory_chunk], 6, ());
        let expected: Vec<Vec<u8>> = vec![
            b"slot-1".to_vec(),
            b"mem-2".to_vec(),
            b"slot-3".to_vec(),
            b"mem-4".to_vec(),
            b"slot-5".to_vec(),
            b"mem-6".to_vec(),
        ];
        assert_eq!(
            drain_merge_driver_bytes(driver),
            expected,
            "slot + shared-memory merge not globally sorted",
        );
    }

    // (v4: `from_slots_bails_on_unfilled_gap_in_decompressed` removed.
    //  v3.1's `decompressed: Mutex<ReorderBuffer<Vec<u8>>>` is now a
    //  `Mutex<VecDeque<Vec<u8>>>` (FIFO). There is no concept of
    //  "ordinals" or "gaps" — the producer pushes in order under the
    //  reader lock, so a missing-block-in-the-middle state is no
    //  longer representable.)

    // ------------------------------------------------------------------------
    // Coverage for the EMBEDDED_IN_RECORD = true parser arm.
    //
    // All three production sort keys (RawCoordinateKey, RawQuerynameKey,
    // TemplateKey) set EMBEDDED_IN_RECORD = true, but TestKey above uses
    // the non-embedded format. This test pins the embedded-format slot path
    // with a synthetic embedded key (TestEmbeddedKey) whose value lives at
    // bytes 0..8 of the record itself.

    /// Embedded-key test type: serialized format is `[len(4)][record(len)]`
    /// (no separate key prefix). The key is the first 8 LE bytes of the record.
    #[derive(Default, PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Debug)]
    struct TestEmbeddedKey(u64);

    impl RawSortKey for TestEmbeddedKey {
        type MergeCtx = ();

        type MergeHint = ();

        const SERIALIZED_SIZE: Option<usize> = Some(0);
        const EMBEDDED_IN_RECORD: bool = true;
        fn extract(_bam: &[u8], _ctx: &crate::keys::SortContext) -> Self {
            unimplemented!("TestEmbeddedKey is only used in MergeDriver::from_slots tests")
        }
        fn extract_from_record(bam: &[u8]) -> Self {
            assert!(bam.len() >= 8, "embedded-key test record must have ≥8 bytes");
            let mut buf = [0u8; 8];
            buf.copy_from_slice(&bam[..8]);
            TestEmbeddedKey(u64::from_le_bytes(buf))
        }
        fn write_to<W: Write>(&self, _w: &mut W) -> std::io::Result<()> {
            // No-op: embedded keys are not written separately; they live
            // inside the record bytes.
            Ok(())
        }
        fn read_from<R: Read>(_r: &mut R) -> std::io::Result<Self> {
            unreachable!(
                "EMBEDDED_IN_RECORD = true keys take the slot_try_next_record \
                 embedded arm, which never calls read_from"
            )
        }
    }

    /// Build a record whose first 8 bytes encode the embedded sort key.
    fn embedded_record(key: u64, tail: &[u8]) -> Vec<u8> {
        let mut v = Vec::with_capacity(8 + tail.len());
        v.extend_from_slice(&key.to_le_bytes());
        v.extend_from_slice(tail);
        v
    }

    /// Serialize records in the embedded-key spill format: `[len(4)][record(len)]`.
    fn serialize_embedded(records: &[Vec<u8>]) -> Vec<u8> {
        let mut out = Vec::new();
        for rec in records {
            #[allow(clippy::cast_possible_truncation)]
            let len = rec.len() as u32;
            out.write_all(&len.to_le_bytes()).unwrap();
            out.write_all(rec).unwrap();
        }
        out
    }

    fn embedded_populated_slot(
        file_id: u32,
        records: &[Vec<u8>],
        n_blocks: usize,
    ) -> StdArc<SortMergeSlot> {
        let bytes = serialize_embedded(records);
        let block_size = bytes.len().div_ceil(n_blocks.max(1));
        let slot = StdArc::new(SortMergeSlot::for_test(file_id, crate::codec::SpillCodec::Bgzf));
        {
            let mut dec = slot.decompressed.lock().unwrap();
            for chunk in bytes.chunks(block_size) {
                dec.push_back(chunk.to_vec());
            }
            slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        }
        slot
    }

    #[test]
    fn from_slots_embedded_key_path_merges_correctly() {
        // Three slots, each with sorted records by embedded key. The merged
        // output must be globally sorted. Exercises `slot_try_next_record`'s
        // EMBEDDED_IN_RECORD arm + the `extract_from_record` key-from-bytes
        // path.
        let file0 =
            vec![embedded_record(1, b"-A"), embedded_record(4, b"-A"), embedded_record(7, b"-A")];
        let file1 =
            vec![embedded_record(2, b"-B"), embedded_record(5, b"-B"), embedded_record(8, b"-B")];
        let file2 =
            vec![embedded_record(3, b"-C"), embedded_record(6, b"-C"), embedded_record(9, b"-C")];
        let slots = vec![
            embedded_populated_slot(0, &file0, 1),
            embedded_populated_slot(1, &file1, 1),
            embedded_populated_slot(2, &file2, 1),
        ];

        let driver = MergeDriver::<TestEmbeddedKey>::from_slots(slots, Vec::new(), 9, ());
        let emitted = drain_merge_driver_bytes(driver);

        // Verify keys are emitted in sorted order. Tail bytes prove the
        // identity of each source.
        let expected = vec![
            embedded_record(1, b"-A"),
            embedded_record(2, b"-B"),
            embedded_record(3, b"-C"),
            embedded_record(4, b"-A"),
            embedded_record(5, b"-B"),
            embedded_record(6, b"-C"),
            embedded_record(7, b"-A"),
            embedded_record(8, b"-B"),
            embedded_record(9, b"-C"),
        ];
        assert_eq!(emitted, expected);
    }

    #[test]
    fn from_slots_skips_empty_slots_mixed_with_populated() {
        // Verify a slot whose decompressed queue is empty (and reader marked EOF)
        // is correctly elided from the LoserTree. Without correct handling
        // the priming step would either panic or produce a phantom source.
        let populated = vec![
            (TestKey(1), b"one".to_vec()),
            (TestKey(2), b"two".to_vec()),
            (TestKey(3), b"three".to_vec()),
        ];

        let empty_slot = StdArc::new(SortMergeSlot::for_test(99, crate::codec::SpillCodec::Bgzf));
        empty_slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        // Insert empty-slot first, then a populated slot. The driver should
        // skip the empty one during priming and only emit records from the
        // populated slot.
        let slots = vec![empty_slot, populated_slot(0, &populated, 1)];

        let driver = MergeDriver::<TestKey>::from_slots(slots, Vec::new(), 3, ());
        let emitted = drain_merge_driver_bytes(driver);
        assert_eq!(emitted, vec![b"one".to_vec(), b"two".to_vec(), b"three".to_vec()]);
    }

    /// `slot_try_load_block` checks `decomp_error` BEFORE `queue_eof`, so a slot
    /// whose decompression failed surfaces `Err` rather than the clean EOF that
    /// an empty queue plus `queue_eof` looks like. Nothing exercised the error
    /// arm: every existing test drives `queue_eof`, so swapping the two checks
    /// would silently convert a failed spill into a short read — records
    /// dropped, exit status zero — with the whole suite green.
    #[test]
    fn decomp_error_surfaces_as_an_error_not_a_clean_eof() {
        let slot = StdArc::new(SortMergeSlot::for_test(7, crate::codec::SpillCodec::Bgzf));
        // The shape a failed decompression leaves behind: queue drained, EOF
        // set, and the error flag raised. Ordering matters — both are set, so
        // only the check order decides the outcome.
        slot.decomp_error.store(true, std::sync::atomic::Ordering::Release);
        slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);

        let mut parser = SlotParserState::new();
        let msg = match slot_try_load_block(&slot, &mut parser) {
            Err(e) => e.to_string(),
            Ok(_) => panic!("a slot with decomp_error must not report a clean EOF"),
        };
        // The whole phrase, not a bare `7`: a lone digit matches any number that
        // happens to appear in the message, which is a weaker claim than "this
        // error named the slot that failed".
        assert!(
            msg.contains("spill decompression error on slot 7"),
            "the error must name the failure and the slot; got: {msg}",
        );
    }

    #[test]
    fn from_slots_records_merged_equals_emitted_total() {
        // `records_merged()` (surfaced by `SortMerge` at completion) must equal
        // the number of records actually emitted across all sources.
        let f0: Vec<(TestKey, Vec<u8>)> =
            (0u64..15).map(|i| (TestKey(i * 2), format!("a-{i:02}").into_bytes())).collect();
        let f1: Vec<(TestKey, Vec<u8>)> =
            (0u64..15).map(|i| (TestKey(i * 2 + 1), format!("b-{i:02}").into_bytes())).collect();
        let total = (f0.len() + f1.len()) as u64;

        let slots = vec![populated_slot(0, &f0, 3), populated_slot(1, &f1, 2)];
        let mut driver = MergeDriver::<TestKey>::from_slots(slots, Vec::new(), total, ());

        let mut emitted = 0u64;
        loop {
            match driver.try_step().expect("try_step") {
                MergeStep::Produced(_) => emitted += 1,
                MergeStep::Done => break,
                MergeStep::Stalled => panic!("unexpected Stalled on pre-populated EOF slots"),
            }
        }
        assert_eq!(emitted, total, "should emit every input record");
        assert_eq!(
            driver.records_merged(),
            total,
            "records_merged must equal the emitted record count"
        );
    }

    /// The two spill-frame layouts the slot parser distinguishes: a separate
    /// key prefix (`TestKey`, `[key(8)][len(4)][body]`) or a key embedded in
    /// the body (`TestEmbeddedKey`, `[len(4)][body]`, the layout of every
    /// production key), whose fast path parses the key from the borrowed body.
    trait Layout: RawSortKey + Default + Copy + Send + Sync + 'static {
        /// The key `k`.
        fn key(k: u64) -> Self;
        /// The record body for key `k` with `tail`.
        fn body(k: u64, tail: &[u8]) -> Vec<u8>;
        /// The spill frame of a record `body` with key `k`.
        fn frame(k: u64, body: &[u8]) -> Vec<u8>;
    }

    impl Layout for TestKey {
        fn key(k: u64) -> Self {
            TestKey(k)
        }
        fn body(_k: u64, tail: &[u8]) -> Vec<u8> {
            tail.to_vec()
        }
        fn frame(k: u64, body: &[u8]) -> Vec<u8> {
            serialize_records(&[(TestKey(k), body.to_vec())])
        }
    }

    impl Layout for TestEmbeddedKey {
        fn key(k: u64) -> Self {
            TestEmbeddedKey(k)
        }
        fn body(k: u64, tail: &[u8]) -> Vec<u8> {
            embedded_record(k, tail)
        }
        fn frame(_k: u64, body: &[u8]) -> Vec<u8> {
            serialize_embedded(&[body.to_vec()])
        }
    }

    /// A run of `(key, body)` records.
    type KeyedBodies = Vec<(u64, Vec<u8>)>;

    /// `(key, L::body(key, tail))` for each `(key, tail)`.
    fn bodies<L: Layout>(records: &[(u64, &[u8])]) -> KeyedBodies {
        records.iter().map(|&(k, tail)| (k, L::body(k, tail))).collect()
    }

    /// A populated, EOF-marked slot holding `records` framed in `L`'s layout
    /// as one block.
    fn layout_slot<L: Layout>(file_id: u32, records: &KeyedBodies) -> StdArc<SortMergeSlot> {
        let bytes: Vec<u8> = records.iter().flat_map(|(k, b)| L::frame(*k, b)).collect();
        let slot = StdArc::new(SortMergeSlot::for_test(file_id, crate::codec::SpillCodec::Bgzf));
        slot.decompressed.lock().unwrap().push_back(bytes);
        slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        slot
    }

    /// A record wholly inside its block is emitted by reference to that block,
    /// with no copy into the per-leaf scratch, in both key layouts.
    #[test]
    fn produced_bytes_are_borrowed_from_the_block() {
        fn check<L: Layout>() {
            let slot = layout_slot::<L>(0, &bodies::<L>(&[(1, b"AAAA"), (2, b"BBBB")]));
            let mut d = MergeDriver::<L>::from_slots(vec![slot], vec![], 2, Default::default());
            let MergeStep::Produced(b) = d.try_step().unwrap() else { panic!("a record") };
            assert_eq!(b, L::body(1, b"AAAA"));
            let p = b.as_ptr();
            assert!(d.block_contains_for_test(0, p), "the record must point into slot 0's block");
            assert_eq!(d.scratch_capacity_for_test(0), 0, "the scratch was never used");
            assert_eq!(d.record_counts(), (1, 0));
        }
        check::<TestKey>();
        check::<TestEmbeddedKey>();
    }

    /// Two runs fed one block at a time; on `Stalled`, the stalled slot gets
    /// its next block and the driver steps again. Blocks are cut at arbitrary
    /// byte boundaries, so records spanning blocks are completed across a
    /// `WouldBlock` while the other source's record stays borrowed. The output
    /// equals a stable merge computed here, and the `(borrowed, assembled)`
    /// split equals the one computed from the cuts: exactly the records whose
    /// frame crosses a block boundary are assembled.
    fn check_incremental_feed<L: Layout>(
        runs: &[KeyedBodies; 2],
        block_records: [usize; 2],
    ) -> Result<(), proptest::test_runner::TestCaseError> {
        let total: usize = runs.iter().map(Vec::len).sum();
        let mut oracle: Vec<&(u64, Vec<u8>)> = runs.iter().flatten().collect();
        oracle.sort_by_key(|(k, _)| *k);
        let oracle: Vec<Vec<u8>> = oracle.into_iter().map(|(_, b)| b.clone()).collect();
        let (fa, spans_a) = incremental_slot::<L>(0, &runs[0], block_records[0]);
        let (fb, spans_b) = incremental_slot::<L>(1, &runs[1], block_records[1]);
        let mut d = MergeDriver::<L>::from_slots(
            vec![StdArc::clone(&fa.slot), StdArc::clone(&fb.slot)],
            vec![],
            total as u64,
            Default::default(),
        );
        let mut out = Vec::new();
        loop {
            match d.try_step().unwrap() {
                MergeStep::Produced(bytes) => out.push(bytes.to_vec()),
                MergeStep::Stalled => {
                    let id = d.stalled_slot().expect("a slot stalled").file_id;
                    let fed = if id == 0 { fa.feed_next() } else { fb.feed_next() };
                    proptest::prop_assert!(fed, "fed past the end");
                }
                MergeStep::Done => break,
            }
        }
        proptest::prop_assert_eq!(out, oracle);
        let assembled = spans_a + spans_b;
        proptest::prop_assert_eq!(d.record_counts(), (total as u64 - assembled, assembled));
        Ok(())
    }

    /// One sorted run: 1..60 records, keys in 0..64 (so the runs tie), bodies
    /// built by `L` from a tail of `tag` then 0..300 arbitrary bytes.
    fn sorted_run<L: Layout>(tag: u8) -> impl proptest::strategy::Strategy<Value = KeyedBodies> {
        use proptest::prelude::*;
        prop::collection::vec((0u64..64, prop::collection::vec(any::<u8>(), 0..300)), 1..60)
            .prop_map(move |mut records| {
                records.sort_by_key(|(k, _)| *k);
                records
                    .into_iter()
                    .map(|(k, rest)| {
                        let mut tail = vec![tag];
                        tail.extend(rest);
                        (k, L::body(k, &tail))
                    })
                    .collect()
            })
    }

    proptest::proptest! {
        /// [`check_incremental_feed`] with a separate key prefix.
        #[test]
        fn incremental_feed_merges_identically(
            a in sorted_run::<TestKey>(b'A'), b in sorted_run::<TestKey>(b'B'),
            block_a in 1usize..9, block_b in 1usize..9,
        ) {
            check_incremental_feed::<TestKey>(&[a, b], [block_a, block_b])?;
        }

        /// [`check_incremental_feed`] with the key embedded in the body.
        #[test]
        fn incremental_feed_merges_identically_with_an_embedded_key(
            a in sorted_run::<TestEmbeddedKey>(b'A'), b in sorted_run::<TestEmbeddedKey>(b'B'),
            block_a in 1usize..9, block_b in 1usize..9,
        ) {
            check_incremental_feed::<TestEmbeddedKey>(&[a, b], [block_a, block_b])?;
        }
    }

    /// `[len][body]` frames: the slot format of an embedded key.
    fn embedded_frames(bodies: &[Vec<u8>]) -> Vec<u8> {
        let mut out = Vec::new();
        for body in bodies {
            out.write_all(&u32::try_from(body.len()).unwrap().to_le_bytes()).unwrap();
            out.write_all(body).unwrap();
        }
        out
    }

    /// A queryname record named `A00132:53:HFHJ:1:1101:{x}:{y}` with the R1
    /// (`0x41`) or R2 (`0x81`) flags.
    fn qn_record(x: u8, y: u8, r1: bool) -> Vec<u8> {
        let name = format!("A00132:53:HFHJ:1:1101:{x}:{y}");
        let flags = if r1 { 0x41 } else { 0x81 };
        fgumi_raw_bam::testutil::make_bam_bytes(0, 1, flags, name.as_bytes(), &[], 4, -1, -1, &[])
    }

    /// `records` sorted by the owned lexicographic key.
    fn qn_sorted(mut records: Vec<Vec<u8>>) -> Vec<Vec<u8>> {
        records.sort_by_key(|b| crate::keys::RawQuerynameLexKey::extract_from_record(b));
        records
    }

    /// One run of 1..60 queryname records sorted by the owned lexicographic
    /// key. Names are drawn from 21 values, so equal names occur within and
    /// across runs; all share a class prefix the merge's table holds.
    fn qn_sorted_run() -> impl proptest::strategy::Strategy<Value = Vec<Vec<u8>>> {
        use proptest::prelude::*;
        prop::collection::vec((0u8..7, 0u8..3, any::<bool>()), 1..60).prop_map(|names| {
            qn_sorted(names.into_iter().map(|(x, y, r1)| qn_record(x, y, r1)).collect())
        })
    }

    /// The sort's class table, holding the runs' common prefix.
    fn qn_table() -> StdArc<crate::queryname_ref::QnClassTable> {
        let table = crate::queryname_ref::QnClassTable::new();
        table.class_for(b"A00132:53:HFHJ:1:");
        table
    }

    /// The `QnMergeKey` path borrows a record wholly inside its block too.
    #[test]
    fn qn_merge_key_produced_bytes_are_borrowed_from_the_block() {
        use crate::queryname_ref::{QnLex, QnMergeKey};
        let a = qn_sorted(vec![qn_record(3, 1, true), qn_record(5, 0, false)]);
        let slot = StdArc::new(SortMergeSlot::for_test(0, crate::codec::SpillCodec::Bgzf));
        slot.decompressed.lock().unwrap().push_back(embedded_frames(&a));
        slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        let mut d = MergeDriver::<QnMergeKey<QnLex>, crate::queryname_ref::QnKey>::from_slots(
            vec![slot],
            vec![],
            2,
            qn_table(),
        );
        let MergeStep::Produced(b) = d.try_step().unwrap() else { panic!("a record") };
        let p = b.as_ptr();
        assert!(d.block_contains_for_test(0, p), "the record must point into slot 0's block");
        assert_eq!(d.scratch_capacity_for_test(0), 0, "the scratch was never used");
        assert_eq!(d.record_counts(), (1, 0));
    }

    proptest::proptest! {
        /// `incremental_feed_merges_identically` with the queryname merge key:
        /// two runs fed one cut block at a time, records spanning blocks across
        /// a `WouldBlock`, merged through the class window against the sort's
        /// table. The output equals a stable merge by the owned key.
        #[test]
        fn qn_merge_key_incremental_feed_merges_identically(
            a in qn_sorted_run(), b in qn_sorted_run(), block_a in 1usize..9, block_b in 1usize..9,
        ) {
            use crate::queryname_ref::{QnLex, QnMergeKey};
            let mut oracle: Vec<Vec<u8>> = a.iter().chain(&b).cloned().collect();
            oracle.sort_by_cached_key(|r| crate::keys::RawQuerynameLexKey::extract_from_record(r));
            let feeder = |id: u32, bodies: &[Vec<u8>], per: usize| {
                let bytes = embedded_frames(bodies);
                let target = (bytes.len() / bodies.len().max(1)).max(1) * per;
                Feeder {
                    slot: StdArc::new(SortMergeSlot::for_test(id, crate::codec::SpillCodec::Bgzf)),
                    blocks: std::cell::RefCell::new(
                        bytes.chunks(target.max(1)).map(<[u8]>::to_vec).collect(),
                    ),
                    seq: std::cell::Cell::new(0),
                }
            };
            let (fa, fb) = (feeder(0, &a, block_a), feeder(1, &b, block_b));
            let mut d = MergeDriver::<QnMergeKey<QnLex>, crate::queryname_ref::QnKey>::from_slots(
                vec![StdArc::clone(&fa.slot), StdArc::clone(&fb.slot)],
                vec![],
                (a.len() + b.len()) as u64,
                qn_table(),
            );
            let mut out = Vec::new();
            loop {
                match d.try_step().unwrap() {
                    MergeStep::Produced(bytes) => out.push(bytes.to_vec()),
                    MergeStep::Stalled => {
                        let id = d.stalled_slot().expect("a slot stalled").file_id;
                        let fed = if id == 0 { fa.feed_next() } else { fb.feed_next() };
                        proptest::prop_assert!(fed, "fed past the end");
                    }
                    MergeStep::Done => break,
                }
            }
            proptest::prop_assert_eq!(out, oracle);
        }
    }

    /// A slot rebuilds an embedded key against the merge's context and carries
    /// the source's hint: a queryname key from a slot takes the sort table's
    /// class, and the hint follows it.
    #[test]
    fn slot_keys_are_rebuilt_against_the_merge_context() {
        use crate::queryname_ref::{QnLex, QnMergeKey};
        let table = qn_table();
        let body = fgumi_raw_bam::testutil::make_bam_bytes(
            0,
            1,
            0x41,
            b"A00132:53:HFHJ:1:1101:5:7",
            &[],
            4,
            -1,
            -1,
            &[],
        );
        let mut hint = 0u8;
        let key = slot_parse_key::<QnMergeKey<QnLex>>(&table, &mut hint, &[], &body).unwrap();
        assert_eq!(key.key.class, 1, "the table's class");
        assert_eq!(hint, 1, "the hint follows the class");
    }

    /// The shared residual (`MemoryShared`) is borrowed from the chunk's arena
    /// (pointer equal to `record_bytes(i)` taken before the chunk moves), with
    /// no copy; the owned residual (`Memory`) is swapped into the scratch. Both
    /// outputs equal the oracle, in both key layouts.
    #[test]
    fn memory_sources_borrow_or_swap() {
        fn check<L: Layout>() {
            let chunk_recs = bodies::<L>(&[(1, b"c-1"), (3, b"c-3"), (5, b"c-5")]);
            let slot_recs = bodies::<L>(&[(2, b"s-2"), (4, b"s-4")]);
            let oracle: Vec<Vec<u8>> =
                [(1, b"c-1"), (2, b"s-2"), (3, b"c-3"), (4, b"s-4"), (5, b"c-5")]
                    .iter()
                    .map(|&(k, t)| L::body(k, t))
                    .collect();
            let keyed = |recs: &KeyedBodies| -> Vec<(L, Vec<u8>)> {
                recs.iter().map(|(k, b)| (L::key(*k), b.clone())).collect()
            };
            // Chunk records carry `c` three bytes from the end of the body.
            let from_chunk = |b: &[u8]| b[b.len() - 3] == b'c';
            // Shared (production): borrowed from the chunk arena.
            let chunk = InMemoryChunk::from_owned_records(keyed(&chunk_recs));
            let ptrs: Vec<usize> =
                (0..3).map(|i| chunk.record_bytes(i).as_ptr() as usize).collect();
            let slot = layout_slot::<L>(0, &slot_recs);
            let mut d =
                MergeDriver::<L>::from_slots(vec![slot], vec![chunk], 5, Default::default());
            let (mut out, mut chunk_ptrs) = (Vec::new(), Vec::new());
            loop {
                match d.try_step().unwrap() {
                    MergeStep::Produced(b) => {
                        if from_chunk(b) {
                            chunk_ptrs.push(b.as_ptr() as usize);
                        }
                        out.push(b.to_vec());
                    }
                    MergeStep::Stalled => panic!("populated slots never stall"),
                    MergeStep::Done => break,
                }
            }
            assert_eq!(out, oracle);
            assert_eq!(chunk_ptrs, ptrs, "chunk records must be borrowed from the arena");
            assert_eq!(d.record_counts(), (5, 0), "two block-borrowed + three chunk-borrowed");
        }
        check::<TestKey>();
        check::<TestEmbeddedKey>();
    }

    /// The driver-level twin of
    /// `embedded_key_on_an_empty_body_is_an_error_not_a_panic`: an
    /// embedded-key `RawCoordinateKey` slot whose single record has length
    /// prefix 0. The borrowed path returns the same corrupt-slot error, not a
    /// panic.
    #[test]
    fn empty_body_on_the_borrowed_path_is_the_same_error() {
        let slot = StdArc::new(SortMergeSlot::for_test(0, crate::codec::SpillCodec::Bgzf));
        slot.bp_commit_read(1, true);
        // [len = 0][no body]
        slot.bp_insert_drain_finalize(0, vec![0u32.to_le_bytes().to_vec()], 1);
        let mut d =
            MergeDriver::<crate::keys::RawCoordinateKey>::from_slots(vec![slot], vec![], 1, ());
        let err = match d.try_step() {
            Err(e) => e,
            Ok(step) => panic!("an empty embedded-key body must be an error, got {step:?}"),
        };
        assert!(err.to_string().contains("empty and cannot carry an embedded sort key"), "{err}");
    }

    /// A slot plus its pre-cut blocks, inserted one at a time on demand.
    struct Feeder {
        slot: StdArc<SortMergeSlot>,
        blocks: std::cell::RefCell<std::collections::VecDeque<Vec<u8>>>,
        seq: std::cell::Cell<u64>,
    }

    impl Feeder {
        /// Insert the next block; `false` when none is left.
        fn feed_next(&self) -> bool {
            let Some(block) = self.blocks.borrow_mut().pop_front() else { return false };
            let last = self.blocks.borrow().is_empty();
            self.slot.bp_commit_read(1, last);
            self.slot.bp_insert_drain_finalize(self.seq.get(), vec![block], 1);
            self.seq.set(self.seq.get() + 1);
            true
        }
    }

    /// An empty slot and a feeder holding `records` framed in `L`'s layout and
    /// cut into blocks of about `block_records` records each, at arbitrary byte
    /// boundaries so records span blocks; plus how many records' frames cross
    /// a block boundary (the parser's slow path, so the driver's `assembled`).
    fn incremental_slot<L: Layout>(
        id: u32,
        records: &KeyedBodies,
        block_records: usize,
    ) -> (Feeder, u64) {
        let frames: Vec<Vec<u8>> = records.iter().map(|(k, b)| L::frame(*k, b)).collect();
        let bytes = frames.concat();
        let target = ((bytes.len() / records.len().max(1)).max(1) * block_records).max(1);
        let (mut spanning, mut start) = (0u64, 0usize);
        for frame in &frames {
            let end = start + frame.len();
            if start / target != (end - 1) / target {
                spanning += 1;
            }
            start = end;
        }
        let blocks = bytes.chunks(target).map(<[u8]>::to_vec).collect();
        let slot = StdArc::new(SortMergeSlot::for_test(id, crate::codec::SpillCodec::Bgzf));
        let feeder =
            Feeder { slot, blocks: std::cell::RefCell::new(blocks), seq: std::cell::Cell::new(0) };
        (feeder, spanning)
    }

    /// A queryname residual record's merge key (rebuilt from its body through
    /// `QnKey`'s `ChunkRef`) equals the key the merge parses for the same body
    /// read from a slot, for both orders.
    #[test]
    fn memory_residual_key_equals_slot_key() {
        use crate::keys::{RawQuerynameKey, RawQuerynameLexKey};
        use crate::queryname_ref::{QnBatchLcp, QnClassTable, QnRefBuilder};
        let table = QnClassTable::new();
        let class = table.class_for(b"A00132:53:HFHJ:1:");
        for name in [b"A00132:53:HFHJ:1:1101:23:45".as_slice(), b"r7x", b""] {
            let body =
                fgumi_raw_bam::testutil::make_bam_bytes(0, 1, 0x81, name, &[], 6, -1, -1, &[]);
            let (key, ..) = QnRefBuilder::with_guess(class).build(
                &table,
                &mut QnBatchLcp::default(),
                &body,
                0,
                0,
            );
            let lex: RawQuerynameLexKey = key.merge_key(&body);
            let natural: RawQuerynameKey = key.merge_key(&body);
            assert_eq!(
                lex.cmp(&slot_parse_key::<RawQuerynameLexKey>(&(), &mut (), &[], &body).unwrap()),
                std::cmp::Ordering::Equal
            );
            assert_eq!(
                natural.cmp(&slot_parse_key::<RawQuerynameKey>(&(), &mut (), &[], &body).unwrap()),
                std::cmp::Ordering::Equal
            );
        }
    }

    /// A queryname residual held as `InMemoryChunk<QnKey>` merges with a slot
    /// into the owned-key order, and every record is borrowed: the residual's
    /// from its arena, the slot's from its single block.
    #[test]
    fn memory_shared_queryname_residual_is_borrowed() {
        use crate::keys::RawQuerynameLexKey;
        use crate::queryname_ref::{
            QnBatchLcp, QnClassTable, QnKey, QnLex, QnRefBuilder, queryname_chunk_from_qn_refs,
        };
        let body = |name: String| {
            fgumi_raw_bam::testutil::make_bam_bytes(
                0,
                1,
                0x41,
                name.as_bytes(),
                &[],
                6,
                -1,
                -1,
                &[],
            )
        };
        let residual: Vec<Vec<u8>> = (0..40).map(|i| body(format!("F:1:{}", i * 2))).collect();
        let spilled: Vec<Vec<u8>> = (0..30).map(|i| body(format!("F:1:{}", i * 2 + 1))).collect();
        // The slot: one block of `[len][body]` records, sorted.
        let mut sorted_spill = spilled.clone();
        sorted_spill.sort_by(|a, b| {
            RawQuerynameLexKey::extract_from_record(a)
                .cmp(&crate::keys::RawQuerynameLexKey::extract_from_record(b))
        });
        let mut block = Vec::new();
        for b in &sorted_spill {
            crate::spill_block::frame_embedded_record_into(&mut block, b).unwrap();
        }
        let slot = StdArc::new(SortMergeSlot::for_test(0, crate::codec::SpillCodec::Bgzf));
        slot.decompressed.lock().unwrap().push_back(block);
        slot.queue_eof.store(true, std::sync::atomic::Ordering::Release);
        // The residual: QnKey refs over an arena, sorted by the window comparator.
        let mut arena = crate::segmented_buf::SegmentedBuf::new();
        let spans: Vec<(u64, u32)> = residual
            .iter()
            .map(|b| (arena.extend_from_slice(b) as u64, u32::try_from(b.len()).unwrap()))
            .collect();
        let arena = StdArc::new(crate::arena_pool::PooledSegmentedBuf::unpooled(arena));
        let table = QnClassTable::new();
        let builder = QnRefBuilder::default();
        let refs = residual
            .iter()
            .zip(&spans)
            .map(|(b, &(o, l))| builder.build(&table, &mut QnBatchLcp::default(), b, o, l))
            .collect();
        let chunk: InMemoryChunk<QnKey> =
            queryname_chunk_from_qn_refs::<QnLex>(&table, StdArc::clone(&arena), refs, None);
        let n = (residual.len() + spilled.len()) as u64;
        let driver =
            MergeDriver::<RawQuerynameLexKey, QnKey>::from_slots(vec![slot], vec![chunk], n, ());
        let mut expected: Vec<Vec<u8>> = residual.iter().chain(&spilled).cloned().collect();
        expected.sort_by(|a, b| {
            RawQuerynameLexKey::extract_from_record(a)
                .cmp(&crate::keys::RawQuerynameLexKey::extract_from_record(b))
        });
        let mut d = driver;
        let mut out = Vec::new();
        while let MergeStep::Produced(b) = d.try_step().unwrap() {
            out.push(b.to_vec());
        }
        assert_eq!(out, expected);
        assert_eq!(d.record_counts(), (n, 0), "every record borrowed, none assembled");
    }
}

mod kind_dispatch_tests {
    //! `merge_driver_for_kind` and `try_step_framed`: the consolidation merge's
    //! routing and its published demand.
    use super::*;
    use crate::{MergeDemand, SpillKeyKind, merge_slots::CONSOLIDATION_SLOT_IDS_START};

    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    fn template_key(primary: u64) -> crate::TemplateKey {
        crate::TemplateKey::new(
            0,
            primary as i32,
            false,
            i32::MAX,
            i32::MAX,
            false,
            0,
            0,
            (0, false),
            0,
            false,
        )
    }

    /// A slot (`file_id`, supply id `slot_id`) already holding `records` framed
    /// as `TemplateKey` records in one decompressed block, at EOF: the
    /// `fgumi-sort`-side twin of pipeline-io's `spill_fixtures::loaded_slot`.
    fn loaded_slot_raw(file_id: u32, slot_id: u32, records: &[(u64, &[u8])]) -> Arc<SortMergeSlot> {
        let mut raw = Vec::new();
        for (k, body) in records {
            crate::frame_keyed_record_into(&mut raw, &template_key(*k), body).unwrap();
        }
        let slot = Arc::new(
            SortMergeSlot::new(
                file_id,
                Arc::new(tempfile::tempfile().unwrap()),
                0,
                0,
                crate::SpillCodec::Zstd,
            )
            .with_slot_id(slot_id),
        );
        slot.bp_commit_read(1, true);
        slot.bp_insert_drain_finalize(0, vec![raw], 1);
        slot
    }

    /// Every kind is routed to its own key type: the driver's kind is derived
    /// from the key type the arm instantiates (`spill_kind_of`), not echoed from
    /// the request. Catches: a copy-pasted arm (e.g. `TemplateCb32` built as
    /// `TemplateTert32`), which would merge with the wrong comparator and frame
    /// the wrong key width.
    #[rstest::rstest]
    fn merge_driver_for_kind_dispatches_every_kind(
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
    ) {
        let empty = |id| {
            let s = Arc::new(SortMergeSlot::new(
                id,
                Arc::new(tempfile::tempfile().unwrap()),
                0,
                0,
                crate::SpillCodec::Zstd,
            ));
            s.bp_commit_read(0, true);
            s.bp_drain_and_finalize();
            s
        };
        let ctx = if matches!(kind, SpillKeyKind::QuerynameLex | SpillKeyKind::QuerynameNatural) {
            MergeCtxErased::Queryname(crate::QnClassTable::new())
        } else {
            MergeCtxErased::None
        };
        let mut d = merge_driver_for_kind(
            kind,
            vec![empty(0), empty(1)],
            ctx,
            Arc::new(MergeDemand::new()),
        );
        assert_eq!(d.key_kind(), Some(kind));
        assert_eq!(d.try_step_framed(&mut Vec::new()).unwrap(), MergeStepFramed::Done);
    }

    /// A slot holding one queryname record named `name` (embedded framing), at
    /// EOF.
    fn queryname_slot(file_id: u32, name: &[u8]) -> Arc<SortMergeSlot> {
        let body = fgumi_raw_bam::testutil::make_bam_bytes(-1, -1, 0x4, name, &[], 0, -1, -1, &[]);
        let mut raw = Vec::new();
        crate::frame_embedded_record_into(&mut raw, &body).unwrap();
        let slot = Arc::new(SortMergeSlot::new(
            file_id,
            Arc::new(tempfile::tempfile().unwrap()),
            0,
            0,
            crate::SpillCodec::Zstd,
        ));
        slot.bp_commit_read(1, true);
        slot.bp_insert_drain_finalize(0, vec![raw], 1);
        slot
    }

    /// The two queryname kinds merge with their own comparator: `read10` sorts
    /// before `read9` lexicographically and after it naturally. Catches: the lex and
    /// natural arms swapped (both frame records identically, so only the
    /// merge order shows it).
    #[rstest::rstest]
    #[case::lexicographic(SpillKeyKind::QuerynameLex, b"read10")]
    #[case::natural(SpillKeyKind::QuerynameNatural, b"read9")]
    fn queryname_kinds_merge_with_their_own_comparator(
        #[case] kind: SpillKeyKind,
        #[case] first: &[u8],
    ) {
        let slots = vec![queryname_slot(0, b"read9"), queryname_slot(1, b"read10")];
        let mut d = merge_driver_for_kind(
            kind,
            slots,
            MergeCtxErased::Queryname(crate::QnClassTable::new()),
            Arc::new(MergeDemand::new()),
        );
        let mut out = Vec::new();
        assert_eq!(d.try_step_framed(&mut out).unwrap(), MergeStepFramed::Produced);
        let name = fgumi_raw_bam::RawRecordView::new(&out[4..]).read_name().to_vec();
        assert_eq!(name, first, "{kind:?}");
    }

    /// The consolidation driver publishes predictions by SLOT id (review r2 I3):
    /// slots with file ids 0/1 but minted slot ids C/C+1; after the first record
    /// (slot 0 wins, slot 1 is runner-up) the prediction is C+1. Catches: the
    /// driver still publishing `file_id`, which collides after cascading.
    #[test]
    fn consolidation_driver_publishes_its_prediction() {
        const C: u32 = CONSOLIDATION_SLOT_IDS_START;
        let s0 = loaded_slot_raw(0, C, &[(1, b"a"), (2, b"b")]);
        let s1 = loaded_slot_raw(1, C + 1, &[(3, b"c"), (4, b"d")]);
        let demand = Arc::new(MergeDemand::new());
        let mut d = merge_driver_for_kind(
            SpillKeyKind::TemplateK40,
            vec![s0, s1],
            MergeCtxErased::None,
            Arc::clone(&demand),
        );
        let mut out = Vec::new();
        assert_eq!(d.try_step_framed(&mut out).unwrap(), MergeStepFramed::Produced);
        assert_eq!(demand.predicted(), Some(C + 1));
        assert_eq!(demand.frontier(), Some(C));
        let mut want = Vec::new();
        crate::frame_keyed_record_into(&mut want, &template_key(1), b"a").unwrap();
        assert_eq!(out, want, "the winner framed with its key prefix");
    }
}

mod bam_slot_tests {
    //! Body-keyed slots over BAM inputs (`fgumi merge`): the header skip, keys
    //! from the bodies, the per-source order check and the full template key.
    use super::*;
    use crate::keys::RawCoordinateKey;
    use crate::{TemplateKey, TemplateKeySeed};
    use fgumi_raw_bam::{RawRecord, SamBuilder, flags};
    use std::ops::Range;

    fn one_ref_header() -> Header {
        use noodles::sam::header::record::value::{Map, map::ReferenceSequence};
        Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(std::num::NonZeroUsize::new(1_000_000).unwrap()),
            )
            .build()
    }

    fn write_records(path: &Path, header: &Header, records: &[RawRecord]) {
        let mut w = fgumi_bam_io::create_raw_bam_writer(path, header, 1, 1).unwrap();
        for r in records {
            w.write_raw_record(r.as_ref()).unwrap();
        }
        w.finish().unwrap();
    }

    /// `n` coordinate-sorted mapped records (`r{i}` at pos 10·i); returns the
    /// path and each record's body bytes in file order.
    fn write_test_bam(dir: &Path, n: usize) -> (PathBuf, Vec<Vec<u8>>) {
        let recs: Vec<RawRecord> = (0..n)
            .map(|i| {
                let mut b = SamBuilder::new();
                b.read_name(format!("r{i}").as_bytes())
                    .ref_id(0)
                    .pos(i32::try_from(10 * i).unwrap())
                    .mapq(60)
                    .cigar_ops(&[8 << 4])
                    .sequence(b"ACGTACGT")
                    .qualities(&[30; 8]);
                b.build()
            })
            .collect();
        let p = dir.join("t.bam");
        write_records(&p, &one_ref_header(), &recs);
        (p, recs.iter().map(|r| r.as_ref().to_vec()).collect())
    }

    fn template_test_record_at(
        name: &str,
        pos: i32,
        cb: Option<&str>,
        mi: Option<&str>,
    ) -> RawRecord {
        let mut b = SamBuilder::new();
        b.read_name(name.as_bytes())
            .flags(flags::PAIRED | flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(pos)
            .mapq(60)
            .cigar_ops(&[8 << 4])
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8])
            .mate_ref_id(0)
            .mate_pos(pos + 100)
            .template_length(108)
            .add_string_tag(SamTag::MC, b"8M");
        if let Some(cb) = cb {
            b.add_string_tag(SamTag::CB, cb.as_bytes());
        }
        if let Some(mi) = mi {
            b.add_string_tag(SamTag::MI, mi.as_bytes());
        }
        b.build()
    }

    fn template_test_record(cb: Option<&str>, mi: Option<&str>) -> Vec<u8> {
        template_test_record_at("q", 100, cb, mi).as_ref().to_vec()
    }

    /// One R1 per template `p{i}` at pos 100·i with MI `"{i}"` and CB only on
    /// odd i, sorted by the full template key; returns path and bodies in file
    /// order.
    fn write_template_sorted_test_bam(
        dir: &Path,
        name: &str,
        range: Range<u32>,
    ) -> (PathBuf, Vec<Vec<u8>>) {
        let header = one_ref_header();
        let mut recs: Vec<RawRecord> = range
            .map(|i| {
                template_test_record_at(
                    &format!("p{i}"),
                    i32::try_from(100 * i).unwrap(),
                    (i % 2 == 1).then_some("ACGT"),
                    Some(&i.to_string()),
                )
            })
            .collect();
        let ctx = TemplateKeySeed::from_header(&header, Some(SamTag::CB), KeyTypesSpec::Full)
            .provision(recs[0].as_ref());
        recs.sort_by_key(|r| ctx.extract_one(r.as_ref()).0);
        let p = dir.join(format!("{name}.bam"));
        write_records(&p, &header, &recs);
        (p, recs.iter().map(|r| r.as_ref().to_vec()).collect())
    }

    /// Decompressed length of a BAM's header block (magic, text, references):
    /// the bytes a body-keyed slot skips.
    fn bam_header_len_for_test(path: &Path) -> u64 {
        let mut raw = Vec::new();
        noodles_bgzf::io::Reader::new(std::fs::File::open(path).unwrap())
            .read_to_end(&mut raw)
            .unwrap();
        let l_text = u32::from_le_bytes(raw[4..8].try_into().unwrap()) as usize;
        let mut at = 8 + l_text;
        let n_ref = u32::from_le_bytes(raw[at..at + 4].try_into().unwrap());
        at += 4;
        for _ in 0..n_ref {
            let l_name = u32::from_le_bytes(raw[at..at + 4].try_into().unwrap()) as usize;
            at += 4 + l_name + 4;
        }
        at as u64
    }

    fn read_header_for_test(path: &Path) -> Header {
        create_raw_bam_reader(path, 1).unwrap().1
    }

    /// Ingest the whole file as one last slice and serve every frame, so the
    /// merge never stalls.
    fn feed_whole_file(slot: &SortMergeSlot, path: &Path) {
        let bytes = std::fs::read(path).unwrap();
        slot.bp_note_issued(bytes.len() as u64);
        slot.bp_ingest_slice(0, fgumi_bam_io::pread::SliceBufferPool::new(1).lease(bytes), true)
            .unwrap();
        let mut dec = crate::SpillBlockDecompressor::new();
        while let Some(block) = slot.bp_claim_raw(u64::MAX) {
            block.decompress_and_publish(&mut dec).unwrap();
        }
        slot.bp_drain_and_finalize();
    }

    fn open_file_slot(path: &Path, id: u32) -> Arc<SortMergeSlot> {
        let f = std::fs::File::open(path).unwrap();
        let len = f.metadata().unwrap().len();
        let slot =
            open_bam_slot(Arc::new(f), id, bam_header_len_for_test(path), Some(len)).unwrap();
        feed_whole_file(&slot, path);
        slot
    }

    /// Bodies framed embedded-key in one block, body-keyed, at EOF.
    fn loaded_body_slot(file_id: u32, bodies: &[Vec<u8>]) -> Arc<SortMergeSlot> {
        body_slot_in_blocks(file_id, bodies, None)
    }

    /// [`loaded_body_slot`], its framed bytes split into two blocks at
    /// `split_at` when given, so a record can straddle the block boundary.
    fn body_slot_in_blocks(
        file_id: u32,
        bodies: &[Vec<u8>],
        split_at: Option<usize>,
    ) -> Arc<SortMergeSlot> {
        let raw: Vec<u8> = bodies
            .iter()
            .flat_map(|b| {
                let mut v = Vec::new();
                crate::frame_embedded_record_into(&mut v, b).unwrap();
                v
            })
            .collect();
        let slot = Arc::new(
            SortMergeSlot::new(
                file_id,
                Arc::new(tempfile::tempfile().unwrap()),
                0,
                0,
                crate::SpillCodec::Bgzf,
            )
            .body_keyed(),
        );
        let blocks = match split_at {
            Some(at) => vec![raw[..at].to_vec(), raw[at..].to_vec()],
            None => vec![raw],
        };
        slot.bp_commit_read(blocks.len(), true);
        let count = blocks.len();
        slot.bp_insert_drain_finalize(0, blocks, count);
        slot
    }

    fn coordinate_body(tid: i32, pos: i32, tag: u32) -> Vec<u8> {
        let mut body = vec![0u8; 36];
        body[0..4].copy_from_slice(&tid.to_le_bytes());
        body[4..8].copy_from_slice(&pos.to_le_bytes());
        body[32..36].copy_from_slice(&tag.to_le_bytes());
        body
    }

    /// The skipped bytes are the BAM header; the first parsed record is the
    /// BAM's first record. Catches: skip applied to compressed bytes, or off by
    /// one block.
    #[test]
    fn skip_bytes_skips_the_bam_header() {
        let dir = tempfile::tempdir().unwrap();
        let (bam, records) = write_test_bam(dir.path(), 200);
        let slot = open_file_slot(&bam, 0);
        let mut d = MergeDriver::<RawCoordinateKey>::from_slots(vec![slot], vec![], 1, ());
        let first = match d.try_step().unwrap() {
            MergeStep::Produced(b) => b.to_vec(),
            other => panic!("{other:?}"),
        };
        assert_eq!(first, records[0]);
    }

    /// Review B5 (r1): BAM slots carry no key prefix; keys come from the body.
    /// Catches: the 40-byte `TemplateKey` prefix read off the BAM body.
    #[test]
    fn template_bam_slot_parses_records_with_no_key_prefix() {
        let dir = tempfile::tempdir().unwrap();
        let (a, ra) = write_template_sorted_test_bam(dir.path(), "a", 0..50);
        let (b, rb) = write_template_sorted_test_bam(dir.path(), "b", 25..75);
        let header = read_header_for_test(&a);
        let ctx = TemplateKeySeed::from_header(&header, Some(SamTag::CB), KeyTypesSpec::Full)
            .provision(&ra[0]);
        let slots = vec![open_file_slot(&a, 0), open_file_slot(&b, 1)];
        let mut d =
            MergeDriver::<TemplateKey>::from_slots(slots, vec![], 1, Some(Arc::clone(&ctx)));
        let mut got = Vec::new();
        while let MergeStep::Produced(body) = d.try_step().unwrap() {
            got.push(body.to_vec());
        }
        let mut want: Vec<Vec<u8>> = ra.into_iter().chain(rb).collect();
        want.sort_by_key(|body| ctx.extract_one(body).0); // stable: input order on ties
        assert_eq!(got, want);
    }

    /// Catches: the order check missing, or reporting the wrong input/index.
    #[test]
    fn input_order_check_rejects_a_descending_pair() {
        let a = loaded_body_slot(
            0,
            &[coordinate_body(0, 10, 1), coordinate_body(0, 20, 2), coordinate_body(0, 15, 3)],
        );
        let b = loaded_body_slot(1, &[coordinate_body(0, 12, 4)]);
        let mut d = MergeDriver::<RawCoordinateKey>::from_slots(vec![a, b], vec![], 1, ())
            .with_input_order_check(vec!["a.bam".into(), "b.bam".into()], "coordinate");
        let err = loop {
            match d.try_step() {
                Ok(MergeStep::Produced(_)) => {}
                Ok(other) => panic!("{other:?}"),
                Err(e) => break e,
            }
        };
        // a.bam: positions 10, 20, 15 → the third record (1-based 3) sorts before its predecessor.
        assert_eq!(err.to_string(), input_order_violation("a.bam", "coordinate", 3));
    }

    /// A template-coordinate input is in order when its core lanes are: the
    /// name-hash tie-breaker is fgumi's own, and an input sorted by another
    /// tool (samtools breaks ties by name) need not agree with it. Catches: the
    /// order check comparing the full key and rejecting such an input.
    #[test]
    fn template_input_order_check_ignores_the_name_hash_tie_breaker() {
        let header = one_ref_header();
        let mut recs: Vec<Vec<u8>> = ["p0", "p1"]
            .iter()
            .map(|name| template_test_record_at(name, 100, None, None).as_ref().to_vec())
            .collect();
        let ctx = TemplateKeySeed::from_header(&header, Some(SamTag::CB), KeyTypesSpec::Full)
            .provision(&recs[0]);
        // Descending full key: same core lanes, so only the name hash orders them.
        recs.sort_by_key(|body| std::cmp::Reverse(ctx.extract_one(body).0));
        let (first, second) = (ctx.extract_one(&recs[0]).0, ctx.extract_one(&recs[1]).0);
        assert_eq!(second.core_cmp(&first), std::cmp::Ordering::Equal);
        assert_eq!(second.cmp(&first), std::cmp::Ordering::Less);

        let slot = loaded_body_slot(0, &recs);
        let mut d = MergeDriver::<TemplateKey>::from_slots(vec![slot], vec![], 1, Some(ctx))
            .with_input_order_check(vec!["a.bam".into()], "template-coordinate");
        let mut got = Vec::new();
        while let MergeStep::Produced(body) = d.try_step().unwrap() {
            got.push(body.to_vec());
        }
        assert_eq!(got, recs);
    }

    /// A user's record is held to the BAM layout before the record rewriter
    /// walks its aux data, both when the record lies in one block and when it
    /// is assembled across two. Catches: the rewrite running first, where a
    /// short record panics in the rewriter instead of failing cleanly.
    #[rstest::rstest]
    #[case::in_one_block(None)]
    #[case::across_blocks(Some(3))]
    fn a_short_record_is_rejected_before_the_rewrite(#[case] split_at: Option<usize>) {
        use std::sync::atomic::{AtomicBool, Ordering};
        let slot = body_slot_in_blocks(0, &[vec![0u8; 10]], split_at);
        let rewritten = Arc::new(AtomicBool::new(false));
        let seen = Arc::clone(&rewritten);
        let rewrite: RecordRewrite = Arc::new(move |_, _| seen.store(true, Ordering::Relaxed));
        let mut d = MergeDriver::<RawCoordinateKey>::from_slots(vec![slot], vec![], 1, ())
            .with_record_rewrite(&rewrite);
        let err = d.try_step().expect_err("a 10-byte record is not a BAM record");
        assert!(
            format!("{err:#}").contains("block_size 10 is shorter than a BAM record"),
            "{err:#}"
        );
        assert!(!rewritten.load(Ordering::Relaxed), "the rewriter saw an unchecked record");
    }

    /// Review B4 (r1): the merge context yields the full key, CB and MI included.
    #[rstest::rstest]
    #[case::no_cb_no_mi(None, None)]
    #[case::cb_only(Some("AAAA"), None)]
    #[case::mi_only(None, Some("7"))]
    #[case::cb_and_mi(Some("CCCC"), Some("7/A"))]
    fn template_merge_context_extracts_the_full_key(
        #[case] cb: Option<&str>,
        #[case] mi: Option<&str>,
    ) {
        let header = one_ref_header();
        let rec = template_test_record(cb, mi);
        let ctx = TemplateKeySeed::from_header(&header, Some(SamTag::CB), KeyTypesSpec::Full)
            .provision(&template_test_record(None, None));
        let merge_ctx: <TemplateKey as RawSortKey>::MergeCtx = Some(Arc::clone(&ctx));
        assert_eq!(
            TemplateKey::extract_from_record_in(&merge_ctx, &mut (), &rec),
            extract_template_key_inline(&rec, &ctx.lib_lookup, ctx.cell_tag, &ctx.cb_hasher)
        );
    }
}

mod slot_ingest_tests {
    //! The slot's positional handle and slice ingest: `open_spill_slot`'s body
    //! range, slice reorder, the raw stash, EOF finalize, and the FIFO mirror,
    //! observed through the merge's own read path (`slot_try_load_block`).
    use super::*;
    use crate::codec::SpillCodec;
    use crate::merge_slots::{IngestOutcome, SortMergeSlot};
    use crate::spill_block::{SpillBlockCompressor, spill_magic, spill_trailer};
    use crate::spill_block_reader::SpillBlockDecompressor;
    use fgumi_bam_io::pread::{SliceBufferPool, SliceLease};
    use std::io::Write;
    use std::sync::atomic::Ordering;

    /// Write a spill file of `blocks` raw blocks (each `per_block` bytes of a
    /// pattern) in `codec`; returns its dir, path and whole contents.
    fn spill_file(
        codec: SpillCodec,
        blocks: usize,
    ) -> (tempfile::TempDir, std::path::PathBuf, Vec<u8>) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("run.spill");
        let mut bytes = spill_magic(codec).to_vec();
        let mut c = SpillBlockCompressor::new(codec, 1).unwrap();
        for b in 0..blocks {
            let raw: Vec<u8> =
                (0..20_000).map(|i| u8::try_from((i * 31 + b * 7) % 251).unwrap()).collect();
            bytes.extend_from_slice(&c.compress_block(&raw).unwrap());
        }
        bytes.extend_from_slice(spill_trailer(codec));
        std::fs::File::create(&path).unwrap().write_all(&bytes).unwrap();
        (dir, path, bytes)
    }

    fn lease_of(pool: &Arc<SliceBufferPool>, bytes: &[u8]) -> SliceLease {
        pool.lease(bytes.to_vec())
    }

    /// Ingest as the read planner would: issued first, then the slice.
    fn ingest(
        slot: &SortMergeSlot,
        seq: u32,
        lease: SliceLease,
        last: bool,
    ) -> std::io::Result<IngestOutcome> {
        slot.bp_note_issued(lease.len() as u64);
        slot.bp_ingest_slice(seq, lease, last)
    }

    /// Claim every stashed frame and publish it, as decompress workers would.
    fn serve_stash(slot: &SortMergeSlot) {
        let mut dec = SpillBlockDecompressor::new();
        while let Some(b) = slot.bp_claim_raw(0) {
            b.decompress_and_publish(&mut dec).unwrap();
        }
    }

    /// Run `f` on another thread; fail if it does not finish within 5 s.
    fn within_5s(f: impl FnOnce() + Send + 'static) {
        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            f();
            let _ = tx.send(());
        });
        rx.recv_timeout(std::time::Duration::from_secs(5)).expect("did not finish within 5 s");
    }

    #[rstest::rstest]
    #[case::bgzf(SpillCodec::Bgzf, 0)]
    #[case::zstd(SpillCodec::Zstd, 4)]
    fn open_spill_slot_reads_codec_and_body_range(
        #[case] codec: SpillCodec,
        #[case] body_start: u64,
    ) {
        let (_dir, path, bytes) = spill_file(codec, 2);
        let slot = open_spill_slot(&path, 3).unwrap();
        assert_eq!(slot.codec, codec);
        assert_eq!((slot.body_start(), slot.len()), (body_start, bytes.len() as u64));
        let r = slot.reader.lock().unwrap();
        assert!(r.slices.is_empty() && r.parser.carry_len() == 0);
    }

    /// A last slice holding only the BGZF EOF marker, after
    /// every earlier frame was claimed and drained, still finalizes the slot;
    /// the merge's read path reaches `Eof` under a watchdog.
    #[test]
    fn eof_marker_only_slice_reaches_eof() {
        let (_dir, path, bytes) = spill_file(SpillCodec::Bgzf, 2);
        let slot = open_spill_slot(&path, 0).unwrap();
        let pool = SliceBufferPool::new(4);
        let body_end = bytes.len() - spill_trailer(SpillCodec::Bgzf).len();
        let out = ingest(&slot, 0, lease_of(&pool, &bytes[..body_end]), false).unwrap();
        assert_eq!((out.frames, out.hit_eof), (2, false));
        serve_stash(&slot);
        let out = ingest(&slot, 1, lease_of(&pool, &bytes[body_end..]), true).unwrap();
        assert_eq!(out, IngestOutcome { frames: 0, stashed_bytes: 0, hit_eof: true });
        let s2 = Arc::clone(&slot);
        within_5s(move || {
            let mut parser = SlotParserState::new();
            loop {
                match slot_try_load_block(&s2, &mut parser).unwrap() {
                    BlockLoad::Eof => break,
                    _ => std::thread::yield_now(),
                }
            }
        });
        assert!(slot.is_drained());
    }

    /// Slices that land out of order are parsed in file order: the stash's
    /// sequence numbers are dense and its frames equal the sequential reader's.
    #[rstest::rstest]
    #[case::bgzf(SpillCodec::Bgzf)]
    #[case::zstd(SpillCodec::Zstd)]
    fn out_of_order_slices_are_reassembled_by_seq(#[case] codec: SpillCodec) {
        let (_dir, path, bytes) = spill_file(codec, 5);
        let slot = open_spill_slot(&path, 0).unwrap();
        let body = &bytes[usize::try_from(slot.body_start()).unwrap()..];
        let oracle = SpillBlockDecompressor::new().read_raw(&mut &body[..], codec, 64).unwrap();
        let pool = SliceBufferPool::new(8);
        let cut = body.len() / 3 + 7;
        let parts = [&body[..cut], &body[cut..2 * cut], &body[2 * cut..]];
        assert_eq!(ingest(&slot, 2, lease_of(&pool, parts[2]), true).unwrap().frames, 0);
        assert_eq!(ingest(&slot, 1, lease_of(&pool, parts[1]), false).unwrap().frames, 0);
        let out = ingest(&slot, 0, lease_of(&pool, parts[0]), false).unwrap();
        assert_eq!((out.frames, out.hit_eof), (oracle.len(), true));
        let stash = slot.raw_stash.lock().unwrap();
        let seqs: Vec<u64> = stash.iter().map(|b| b.seq).collect();
        assert_eq!(seqs, (0..oracle.len() as u64).collect::<Vec<_>>());
        let frames: Vec<Vec<u8>> = stash.iter().map(|b| b.frame.to_vec()).collect();
        assert_eq!(frames, oracle);
        assert_eq!(slot.stash_len.load(Ordering::Relaxed) as usize, oracle.len());
        assert_eq!(slot.issued_bytes.load(Ordering::Relaxed), 0, "every issued byte was parsed");
    }

    /// A last slice that ends mid-frame fails the slot.
    #[test]
    fn truncated_last_slice_marks_the_slot_failed() {
        let (_dir, path, bytes) = spill_file(SpillCodec::Zstd, 2);
        let slot = open_spill_slot(&path, 0).unwrap();
        let pool = SliceBufferPool::new(2);
        let err = ingest(&slot, 0, lease_of(&pool, &bytes[4..bytes.len() - 9]), true).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::UnexpectedEof);
        assert!(slot.has_error() && slot.queue_eof.load(Ordering::Acquire));
    }

    /// The front escape admits the last block while the FIFO
    /// is full; it waits in `reorder`. Once the consumer frees room, the
    /// drain-only path moves it and sets `queue_eof`. The oracle is the merge's
    /// read path under a watchdog.
    #[test]
    fn last_block_inserted_into_full_fifo_still_reaches_eof() {
        let cap = crate::merge_slots::PHASE2_DECOMP_CAP;
        let slot = Arc::new(SortMergeSlot::for_test(0, SpillCodec::Bgzf));
        slot.bp_stash_frames_for_test(
            (0..=cap).map(|i| vec![u8::try_from(i % 251).unwrap()]).collect(),
            true,
        );
        for _ in 0..cap {
            let b = slot.bp_claim_raw(u64::MAX).unwrap();
            slot.bp_insert_drain_finalize(b.seq, vec![b.frame.to_vec()], 1);
        }
        assert_eq!(slot.fifo_len(), cap, "the FIFO is full");
        let last = slot.bp_claim_raw(u64::MAX).expect("the front is admitted over a full FIFO");
        slot.bp_insert_drain_finalize(last.seq, vec![last.frame.to_vec()], 1);
        assert_eq!(slot.fifo_len(), cap, "the last block waits in reorder");
        assert!(!slot.queue_eof.load(Ordering::Acquire));
        let s2 = Arc::clone(&slot);
        within_5s(move || {
            let mut parser = SlotParserState::new();
            loop {
                // What a decompress pass / self-serve does when it finds no claim.
                s2.bp_drain_and_finalize();
                match slot_try_load_block(&s2, &mut parser).unwrap() {
                    BlockLoad::Eof => break,
                    _ => std::thread::yield_now(),
                }
            }
        });
    }

    /// Borrowed frames claimed and decompressed by several threads deliver
    /// every block in order, and every slice returns to its pool once its
    /// frames are dropped (20 runs).
    #[test]
    fn stash_claims_with_borrowed_frames_release_every_slice() {
        for _ in 0..20 {
            let (_dir, path, bytes) = spill_file(SpillCodec::Zstd, 12);
            let slot = open_spill_slot(&path, 0).unwrap();
            let pool = SliceBufferPool::new(8);
            let body = &bytes[4..];
            let oracle: Vec<Vec<u8>> = {
                let mut dec = SpillBlockDecompressor::new();
                dec.read_raw(&mut &body[..], SpillCodec::Zstd, 64)
                    .unwrap()
                    .iter()
                    .map(|f| dec.decompress_one(SpillCodec::Zstd, f).unwrap())
                    .collect()
            };
            let cut = body.len() / 4 + 3;
            let chunks: Vec<&[u8]> = body.chunks(cut).collect();
            let n = chunks.len();
            std::thread::scope(|sc| {
                for _ in 0..4 {
                    let slot = &slot;
                    sc.spawn(move || {
                        let mut dec = SpillBlockDecompressor::new();
                        while !slot.queue_eof.load(Ordering::Acquire) {
                            if let Some(b) = slot.bp_claim_raw(u64::MAX) {
                                b.decompress_and_publish(&mut dec).unwrap();
                            } else {
                                slot.bp_drain_and_finalize();
                                std::thread::yield_now();
                            }
                        }
                    });
                }
                // Slices land in reverse order; the slot parses them in file order.
                for i in (0..n).rev() {
                    let seq = u32::try_from(i).unwrap();
                    ingest(&slot, seq, lease_of(&pool, chunks[i]), i == n - 1).unwrap();
                }
                let mut parser = SlotParserState::new();
                let mut got = Vec::new();
                loop {
                    match slot_try_load_block(&slot, &mut parser).unwrap() {
                        BlockLoad::Loaded => got.push(std::mem::take(&mut parser.current_buf)),
                        BlockLoad::Eof => break,
                        BlockLoad::WouldBlock => std::thread::yield_now(),
                    }
                }
                assert_eq!(got, oracle);
            });
            assert_eq!(pool.resident_bytes(), 0, "every slice returned to its pool");
            assert_eq!(slot.stash_bytes.load(Ordering::Relaxed), 0, "every charge released");
        }
    }

    /// The stash charges what it holds: a slice in full while any frame
    /// borrowing it is stashed (not the frames' lengths), released with its
    /// last stashed frame; plus the parser's carried partial frame. Every
    /// charge is released once every frame is claimed.
    #[rstest::rstest]
    #[case::bgzf(SpillCodec::Bgzf)]
    #[case::zstd(SpillCodec::Zstd)]
    fn stash_charges_the_slices_its_frames_pin(#[case] codec: SpillCodec) {
        let (_dir, path, bytes) = spill_file(codec, 6);
        let slot = open_spill_slot(&path, 0).unwrap();
        let pool = SliceBufferPool::new(8);
        let body = &bytes[usize::try_from(slot.body_start()).unwrap()..];
        let cut = body.len() / 2 + 5;
        let charged = || slot.stash_bytes.load(Ordering::Relaxed);
        let carry = || slot.reader.lock().unwrap().parser.carry_capacity() as u64;

        let out = ingest(&slot, 0, lease_of(&pool, &body[..cut]), false).unwrap();
        let carry0 = carry();
        assert!(out.frames >= 2 && carry0 > 0, "a straddling frame and several within");
        assert_eq!(charged(), cut as u64 + carry0, "the whole slice and the carry");
        assert_eq!(out.stashed_bytes, i64::try_from(charged()).unwrap());
        for _ in 1..out.frames {
            drop(slot.bp_claim_raw(0).unwrap());
            assert_eq!(charged(), cut as u64 + carry0, "the slice stays charged");
        }
        drop(slot.bp_claim_raw(0).unwrap());
        assert_eq!(charged(), carry0, "the slice is released with its last frame");

        let out = ingest(&slot, 1, lease_of(&pool, &body[cut..]), true).unwrap();
        assert!(out.hit_eof);
        assert_eq!(carry(), 0);
        let stashed: u64 = slot.raw_stash.lock().unwrap().iter().map(|b| b.charge).sum();
        assert_eq!(charged(), stashed, "the carry moved into the owned frame's charge");
        while slot.bp_claim_raw(0).is_some() {}
        assert_eq!(charged(), 0);
    }

    /// The lock-free FIFO-length mirror equals the FIFO's length after a push
    /// (drain) and after the consumer's pop.
    #[test]
    fn fifo_len_mirror_tracks_push_and_pop() {
        let slot = Arc::new(SortMergeSlot::for_test(0, SpillCodec::Bgzf));
        slot.bp_commit_read(3, true);
        slot.bp_insert_drain_finalize(0, vec![vec![1], vec![2], vec![3]], 3);
        assert_eq!((slot.fifo_len_relaxed(), slot.fifo_len()), (3, 3));
        let mut parser = SlotParserState::new();
        assert!(matches!(slot_try_load_block(&slot, &mut parser).unwrap(), BlockLoad::Loaded));
        assert_eq!((slot.fifo_len_relaxed(), slot.fifo_len()), (2, 2));
        slot.fifo_cap.store(1, Ordering::Relaxed);
        assert_eq!(slot.reorder_len_relaxed(), 0);
    }
}
