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
