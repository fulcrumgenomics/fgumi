//! Unit tests for `SpillCompress::compress_event` (the `StepCtx`-free core).

use super::*;
use crate::sort::protocol::MemoryChunkErased;
use fgumi_sort::{InMemoryChunk, RawCoordinateKey, SpillBlockDecompressor};
use rstest::rstest;

fn raw_block(ordinal: u64, file_id: u32, is_last: bool, bytes: Vec<u8>) -> SpillBlockEvent {
    SpillBlockEvent::Block {
        ordinal,
        file_id,
        is_last_in_file: is_last,
        records_ingested_so_far: 0,
        bytes,
    }
}

#[rstest]
#[case(SpillCodec::Zstd)]
#[case(SpillCodec::Bgzf)]
fn block_payload_is_compressed_and_routing_preserved(#[case] codec: SpillCodec) {
    let mut step = SpillCompress::new(codec, 1, 1 << 20);
    let raw = vec![0xABu8; 4096];
    let out = step.compress_event(raw_block(7, 2, true, raw.clone())).unwrap();
    let SpillBlockEvent::Block { ordinal, file_id, is_last_in_file, bytes, .. } = out else {
        panic!("expected Block");
    };
    assert_eq!(ordinal, 7, "ordinal must be preserved ({codec:?})");
    assert_eq!(file_id, 2, "file_id must be preserved ({codec:?})");
    assert!(is_last_in_file, "is_last must be preserved ({codec:?})");
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
    let mut step = SpillCompress::new(SpillCodec::Zstd, 1, 1 << 20);

    let chunk = MemoryChunkErased::Coordinate(InMemoryChunk::from_owned_records(vec![(
        RawCoordinateKey { sort_key: 1 },
        vec![9u8; 8],
    )]));
    let residual = SpillBlockEvent::Residual { ordinal: 5, chunk, records_ingested_so_far: 1 };
    let out = step.compress_event(residual).unwrap();
    assert!(matches!(out, SpillBlockEvent::Residual { ordinal: 5, .. }));

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
    let mut step = SpillCompress::new(SpillCodec::Zstd, 1, 1 << 20);
    // Force the original to build its compressor.
    let _ = step.compress_event(raw_block(0, 0, true, vec![1u8; 16])).unwrap();
    assert!(step.compressor.is_some(), "original built its compressor");
    let fresh = step.clone();
    assert!(fresh.compressor.is_none(), "clone must start with no compressor");
}
