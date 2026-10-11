//! Opening a spill slot allocates no read buffer: the slot carries a
//! positional handle, and read buffers come from the read path's slice pool
//! on first use. Its own binary, because it installs a counting global
//! allocator (`dhat`, in testing mode).
#![cfg(not(loom))]

use std::io::Write;

#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

#[test]
fn open_spill_slot_allocates_no_read_buffer() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("run.spill");
    let mut file = std::fs::File::create(&path).unwrap();
    file.write_all(&vec![0u8; 1 << 16]).unwrap();
    file.write_all(&fgumi_bgzf::BGZF_EOF).unwrap();
    drop(file);
    let _profiler = dhat::Profiler::builder().testing().build();
    let before = dhat::HeapStats::get().total_bytes;
    let slot = fgumi_sort::open_spill_slot(&path, 0).unwrap();
    let allocated = dhat::HeapStats::get().total_bytes - before;
    assert!(allocated < 8 * 1024, "opening a slot allocated {allocated} bytes (no 8 KiB reader)");
    // The parse-state seam is test support (`test-utils`, on in the workspace
    // test build); the allocation check above stands without it.
    #[cfg(feature = "test-utils")]
    assert!(slot.parse_state_is_empty_for_test());
    drop(slot);
}
