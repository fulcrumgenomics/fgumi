//! Unit tests for `SpillBlockCompress::compress_event` (the `StepCtx`-free core).

use super::*;
use crate::sort::protocol::MemoryChunkErased;
use fgumi_sort::{InMemoryChunk, RawCoordinateKey, SpillBlockDecompressor, SpillKeyKind};
use rstest::rstest;

fn raw_block(
    ordinal: u64,
    file_id: u32,
    is_last: bool,
    records_ingested_so_far: u64,
    bytes: Vec<u8>,
) -> SpillBlockEvent {
    SpillBlockEvent::Block {
        ordinal,
        file_id,
        key_kind: SpillKeyKind::TemplateK24,
        is_last_in_file: is_last,
        records_ingested_so_far,
        bytes,
    }
}

#[rstest]
#[case(SpillCodec::Zstd)]
#[case(SpillCodec::Bgzf)]
fn block_payload_is_compressed_and_routing_preserved(#[case] codec: SpillCodec) {
    let mut step = SpillBlockCompress::new(codec, 1, 1 << 20);
    let raw = vec![0xABu8; 4096];
    let out = step.compress_event(raw_block(7, 2, true, 123, raw.clone())).unwrap();
    let SpillBlockEvent::Block {
        ordinal,
        file_id,
        key_kind,
        is_last_in_file,
        records_ingested_so_far,
        bytes,
    } = out
    else {
        panic!("expected Block");
    };
    assert_eq!(ordinal, 7, "ordinal must be preserved ({codec:?})");
    assert_eq!(file_id, 2, "file_id must be preserved ({codec:?})");
    assert_eq!(key_kind, SpillKeyKind::TemplateK24, "key_kind must be preserved ({codec:?})");
    assert!(is_last_in_file, "is_last must be preserved ({codec:?})");
    assert_eq!(
        records_ingested_so_far, 123,
        "records_ingested_so_far must pass through compression unchanged ({codec:?})"
    );
    assert_ne!(bytes, raw, "payload must change (be compressed) ({codec:?})");
    assert!(!bytes.is_empty(), "compressed payload non-empty ({codec:?})");

    // Round-trip through the matching decoder (independent oracle): a compressor
    // that silently corrupts still changes the bytes, so `assert_ne!` alone is
    // too weak. `read_raw` parses the block framing (zstd length prefix / BGZF
    // block) and `decompress_one` inverts the codec; the recovered payload must
    // equal the exact input.
    let mut dec = SpillBlockDecompressor::new();
    let mut cursor = std::io::Cursor::new(&bytes[..]);
    let frames = dec.read_raw(&mut cursor, codec, 64).unwrap();
    let mut round = Vec::new();
    for frame in &frames {
        round.extend_from_slice(&dec.decompress_one(codec, frame).unwrap());
    }
    assert_eq!(round, raw, "compressed block must round-trip back to input ({codec:?})");
}

#[test]
fn residual_and_announced_pass_through_unchanged() {
    let mut step = SpillBlockCompress::new(SpillCodec::Zstd, 1, 1 << 20);

    let chunk = MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(vec![(
        RawCoordinateKey { sort_key: 1 },
        vec![9u8; 8],
    )]));
    let residual = SpillBlockEvent::Residual { ordinal: 5, chunk, records_ingested_so_far: 1 };
    let out = step.compress_event(residual).unwrap();
    let SpillBlockEvent::Residual { ordinal, records_ingested_so_far, .. } = out else {
        panic!("expected Residual");
    };
    assert_eq!(ordinal, 5, "residual ordinal must pass through unchanged");
    assert_eq!(
        records_ingested_so_far, 1,
        "residual records_ingested_so_far must pass through unchanged (not reset)"
    );

    let announced = SpillBlockEvent::AllAnnounced {
        ordinal: 6,
        slot_count: 2,
        memory_chunk_count: 1,
        total_records: 3,
    };
    let out = step.compress_event(announced).unwrap();
    assert!(matches!(
        out,
        SpillBlockEvent::AllAnnounced {
            ordinal: 6,
            slot_count: 2,
            memory_chunk_count: 1,
            total_records: 3,
        }
    ));
}

#[test]
fn clone_starts_with_fresh_lazy_compressor() {
    let mut step = SpillBlockCompress::new(SpillCodec::Zstd, 1, 1 << 20);
    // Force the original to build its compressor.
    let _ = step.compress_event(raw_block(0, 0, true, 0, vec![1u8; 16])).unwrap();
    assert!(step.compressor.is_some(), "original built its compressor");
    let fresh = step.clone();
    assert!(fresh.compressor.is_none(), "clone must start with no compressor");
}

#[test]
fn new_worker_copy_is_independent_of_the_template() {
    let mut template = SpillBlockCompress::new(SpillCodec::Bgzf, 3, 4096);
    let _ = template.compress_event(raw_block(0, 0, true, 0, vec![2u8; 32])).unwrap();

    let worker = template.new_worker_copy();
    // Config is inherited...
    assert_eq!(worker.codec, SpillCodec::Bgzf);
    assert_eq!(worker.compression, 3);
    assert_eq!(worker.output_byte_limit, 4096);
    // ...but the compressor is per-worker, so workers cannot share codec state.
    assert!(worker.compressor.is_none(), "each worker builds its own compressor lazily");
}

#[test]
fn profile_advertises_parallel_byordinal_with_a_byte_bounded_queue() {
    let step = SpillBlockCompress::new(SpillCodec::Zstd, 1, 8192);
    let profile = step.profile();
    assert_eq!(profile.name, "SpillBlockCompress");
    // Parallel + ByItemOrdinal is what lets any worker compress any block while
    // `SpillWrite` still receives a dense, in-order stream.
    assert_eq!(profile.kind, StepKind::Parallel);
    assert!(!profile.sticky);
    assert_eq!(profile.branch_ordering, vec![BranchOrdering::ByItemOrdinal]);
    match profile.output_queues.as_slice() {
        [QueueSpec::ByteBounded { limit_bytes }] => assert_eq!(*limit_bytes, 8192),
        other => panic!("expected a single byte-bounded queue, got {other:?}"),
    }
}

#[rstest]
#[case::zstd(SpillCodec::Zstd)]
#[case::bgzf(SpillCodec::Bgzf)]
fn an_empty_block_round_trips_to_empty(#[case] codec: SpillCodec) {
    let mut step = SpillBlockCompress::new(codec, 1, 1 << 20);
    let out = step.compress_event(raw_block(0, 0, true, 0, Vec::new())).unwrap();
    let SpillBlockEvent::Block { bytes, .. } = out else { panic!("expected Block") };

    let mut dec = SpillBlockDecompressor::new();
    let mut cursor = std::io::Cursor::new(&bytes[..]);
    let frames = dec.read_raw(&mut cursor, codec, 64).unwrap();
    let mut round = Vec::new();
    for frame in &frames {
        round.extend_from_slice(&dec.decompress_one(codec, frame).unwrap());
    }
    assert!(round.is_empty(), "an empty block must decompress back to empty ({codec:?})");
}

#[rstest]
#[case::zstd(SpillCodec::Zstd)]
#[case::bgzf(SpillCodec::Bgzf)]
fn consecutive_blocks_reuse_one_compressor_and_stay_independent(#[case] codec: SpillCodec) {
    // The compressor is built once and reused across blocks; each block must still
    // decode standalone, since `SpillWrite` may interleave files.
    let mut step = SpillBlockCompress::new(codec, 1, 1 << 20);
    let payloads = [vec![0x11u8; 512], vec![0x22u8; 1024], vec![0x33u8; 64]];

    for (i, payload) in payloads.iter().enumerate() {
        let out = step
            .compress_event(raw_block(i as u64, 0, i == payloads.len() - 1, 0, payload.clone()))
            .unwrap();
        let SpillBlockEvent::Block { bytes, .. } = out else { panic!("expected Block") };

        let mut dec = SpillBlockDecompressor::new();
        let mut cursor = std::io::Cursor::new(&bytes[..]);
        let frames = dec.read_raw(&mut cursor, codec, 64).unwrap();
        let mut round = Vec::new();
        for frame in &frames {
            round.extend_from_slice(&dec.decompress_one(codec, frame).unwrap());
        }
        assert_eq!(&round, payload, "block {i} must decode standalone ({codec:?})");
    }
    assert!(step.compressor.is_some(), "the compressor is built once and retained");
}

/// Spill compression under `--sort-threads`: the phase-1 admission contract.
#[test]
fn capped_step_follows_the_admission_contract() {
    fgumi_pipeline_core::testing::assert_admission_contract(
        "sort-phase1",
        |cap| SpillBlockCompress::new(SpillCodec::Zstd, 1, 1 << 20).with_phase_cap(cap),
        SpillBlockEvent::AllAnnounced {
            ordinal: 0,
            slot_count: 0,
            memory_chunk_count: 0,
            total_records: 0,
        },
    );
}

/// A coordinate chunk of `n` records with bodies of varied length.
fn range_chunk(n: u64) -> Arc<MemoryChunkErased> {
    let recs = (0..n)
        .map(|k| {
            (
                RawCoordinateKey { sort_key: k },
                vec![u8::try_from(k % 200).unwrap(); 20 + usize::try_from(k * 7 % 60).unwrap()],
            )
        })
        .collect();
    Arc::new(MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(recs)))
}

fn range(chunk: &Arc<MemoryChunkErased>, ordinal: u64, start: u64, end: u64) -> SpillBlockEvent {
    let framed_bytes =
        (start..end).map(|i| chunk.framed_len(usize::try_from(i).unwrap())).sum::<usize>();
    SpillBlockEvent::Range {
        ordinal,
        file_id: 3,
        key_kind: chunk.key_kind(),
        is_last_in_file: false,
        records_ingested_so_far: 99,
        chunk: Arc::clone(chunk),
        start,
        end,
        framed_bytes: u32::try_from(framed_bytes).unwrap(),
        frame_bytes: 0,
    }
}

/// A range compresses to exactly the bytes its framed block compresses to, and
/// keeps its routing.
#[rstest]
#[case(SpillCodec::Zstd)]
#[case(SpillCodec::Bgzf)]
fn range_compresses_byte_identically_to_the_framed_block(#[case] codec: SpillCodec) {
    let chunk = range_chunk(40);
    let mut framed = Vec::new();
    for i in 5..31 {
        chunk.frame_record_into(i, &mut framed).unwrap();
    }
    let mut step = SpillBlockCompress::new(codec, 1, 1 << 20);
    let from_range = step.compress_event(range(&chunk, 7, 5, 31)).unwrap();
    let from_block = SpillBlockCompress::new(codec, 1, 1 << 20)
        .compress_event(SpillBlockEvent::Block {
            ordinal: 7,
            file_id: 3,
            key_kind: SpillKeyKind::Coordinate,
            is_last_in_file: false,
            records_ingested_so_far: 99,
            bytes: framed,
        })
        .unwrap();
    let (
        SpillBlockEvent::Block { ordinal, file_id, records_ingested_so_far, bytes: a, .. },
        SpillBlockEvent::Block { bytes: b, .. },
    ) = (from_range, from_block)
    else {
        panic!("both arms produce compressed blocks");
    };
    assert_eq!((ordinal, file_id, records_ingested_so_far), (7, 3, 99));
    assert_eq!(a, b, "{codec:?}: a range compresses like its framed block");
}

/// One worker frames every range into the same scratch buffer: after the first
/// block its address and capacity never change.
#[test]
fn the_raw_scratch_is_retained_across_blocks() {
    let chunk = range_chunk(60);
    let mut step = SpillBlockCompress::new(SpillCodec::Zstd, 1, 1 << 20);
    step.compress_event(range(&chunk, 0, 0, 20)).unwrap();
    let first = step.raw_scratch();
    assert!(
        first.1 >= fgumi_bgzf::BGZF_MAX_BLOCK_SIZE,
        "the worker keeps a full-block scratch after compressing (capacity {})",
        first.1
    );
    step.compress_event(range(&chunk, 1, 20, 40)).unwrap();
    assert_eq!(step.raw_scratch(), first, "block 2 reused the scratch");
    step.compress_event(range(&chunk, 2, 40, 60)).unwrap();
    assert_eq!(step.raw_scratch(), first, "block 3 reused the scratch");
}

/// A range whose framing disagrees with the planner's `framed_bytes`, or whose
/// end runs past the chunk, is an `io::Error` in every build, as the serial
/// run-final path's `check_run_final_block` is, never a debug-only assert or a
/// panic.
#[rstest]
#[case::framed_bytes_disagree(5, 31, 1)]
#[case::end_past_the_chunk(5, 41, 0)]
#[case::start_after_end(31, 5, 0)]
fn a_range_the_planner_and_framing_disagree_on_is_an_error(
    #[case] start: u64,
    #[case] end: u64,
    #[case] framed_bytes_off_by: u32,
) {
    let chunk = range_chunk(40);
    let SpillBlockEvent::Range { framed_bytes, .. } = range(&chunk, 7, 5, 31) else {
        unreachable!("range() builds a Range");
    };
    let event = SpillBlockEvent::Range {
        ordinal: 7,
        file_id: 3,
        key_kind: chunk.key_kind(),
        is_last_in_file: false,
        records_ingested_so_far: 99,
        chunk: Arc::clone(&chunk),
        start,
        end,
        framed_bytes: framed_bytes + framed_bytes_off_by,
        frame_bytes: 0,
    };
    let Err(err) = SpillBlockCompress::new(SpillCodec::Zstd, 1, 1 << 20).compress_event(event)
    else {
        panic!("a disagreeing range must be rejected");
    };
    assert!(err.to_string().starts_with("SpillBlockCompress: "), "{err}");
}

/// A range spanning several frames compresses to exactly the concatenation of
/// those frames compressed one range each: the planner's cut rule re-cuts it,
/// and each frame stays an independent unit, so the spill bytes do not change.
#[rstest]
#[case(SpillCodec::Zstd)]
#[case(SpillCodec::Bgzf)]
fn a_multi_frame_range_compresses_like_its_frames_one_by_one(#[case] codec: SpillCodec) {
    const FRAME: usize = 300;
    let chunk = range_chunk(60);
    let mut cuts = vec![2usize];
    while *cuts.last().unwrap() < 55 {
        let next =
            crate::sort::spill_gather::plan_next_block(&chunk, *cuts.last().unwrap(), FRAME).0;
        cuts.push(next.min(55));
    }
    assert!(cuts.len() > 3, "the range spans several frames: {cuts:?}");
    let mut whole = range(&chunk, 7, 2, 55);
    if let SpillBlockEvent::Range { frame_bytes, .. } = &mut whole {
        *frame_bytes = u32::try_from(FRAME).unwrap();
    }
    let mut step = SpillBlockCompress::new(codec, 1, 1 << 20);
    let SpillBlockEvent::Block { bytes: multi, ordinal, .. } = step.compress_event(whole).unwrap()
    else {
        panic!("a range compresses to a block");
    };
    let mut one_by_one = Vec::new();
    for pair in cuts.windows(2) {
        let event = range(&chunk, 0, pair[0] as u64, pair[1] as u64);
        let SpillBlockEvent::Block { bytes, .. } = step.compress_event(event).unwrap() else {
            panic!("a range compresses to a block");
        };
        one_by_one.extend_from_slice(&bytes);
    }
    assert_eq!(ordinal, 7);
    assert_eq!(multi, one_by_one, "{codec:?}: a multi-frame range is its frames concatenated");
}

/// With a buffer pool, the compressed block is written into a recycled buffer.
#[test]
fn compressed_blocks_reuse_pooled_buffers() {
    let pool = SpillBufPool::new();
    let recycled = Vec::with_capacity(1 << 20);
    let ptr = recycled.as_ptr();
    pool.give(recycled);
    let mut step =
        SpillBlockCompress::new(SpillCodec::Zstd, 1, 1 << 20).with_buf_pool(Arc::clone(&pool));
    let chunk = range_chunk(40);
    let SpillBlockEvent::Block { bytes, .. } =
        step.compress_event(range(&chunk, 0, 0, 20)).unwrap()
    else {
        panic!("a range compresses to a block");
    };
    assert_eq!(bytes.as_ptr(), ptr, "the block reused the pooled buffer");
    assert!(pool.take().capacity() == 0, "the pool's one buffer was taken");
}

/// The pool hands back what it was given (cleared), ignores empty buffers, and
/// keeps at most `SPILL_BUF_POOL_CAP` buffers.
#[test]
fn the_buffer_pool_is_bounded_and_skips_empty_buffers() {
    let pool = SpillBufPool::new();
    pool.give(Vec::new());
    assert_eq!(pool.take().capacity(), 0, "an unallocated buffer is not kept");
    for _ in 0..SPILL_BUF_POOL_CAP + 10 {
        pool.give(vec![1u8; 16]);
    }
    let mut kept = 0;
    loop {
        let buf = pool.take();
        if buf.capacity() == 0 {
            break;
        }
        assert!(buf.is_empty(), "a recycled buffer comes back cleared");
        kept += 1;
    }
    assert_eq!(kept, SPILL_BUF_POOL_CAP);
}
