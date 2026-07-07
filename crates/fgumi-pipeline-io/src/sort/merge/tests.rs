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
    chunks.push(MemoryChunkErased::Coordinate(chunk));

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
    let slot = Arc::new(SortMergeSlot::new(
        0,
        std::io::BufReader::new(tempfile::tempfile().unwrap()),
        fgumi_sort::SpillCodec::Bgzf,
    ));
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
