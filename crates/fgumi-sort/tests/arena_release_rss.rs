//! A retired arena pool's last arena returns its pages to the OS as soon as
//! the thread that allocated it collects: at once when that thread frees it
//! (`ArenaPool::release` calls `force_mi_collect`), at the owner's next
//! collect when another thread does. A live pool keeps a returned arena
//! resident for reuse. Its own binary,
//! because the behaviour under test is mimalloc's (the CLI's global
//! allocator): a freed arena-served block only schedules a purge, so without
//! the collect its pages stay resident for mimalloc's purge delay.
//!
//! 256 MiB is below mimalloc's `arena_max_object_size` (2 GiB), so the segment
//! is served from a mimalloc arena on Linux and macOS alike, which is the lazy
//! purge path; nothing is allocated between the free and the measurement, so
//! macOS (whose `phys_footprint` re-counts a reused range) measures it too.
//!
//! RSS is sampled on Linux and macOS only, so the binary is empty elsewhere.
#![cfg(all(not(loom), any(target_os = "linux", target_os = "macos")))]

use std::sync::{Arc, Mutex};

use fgumi_sort::{ArenaPool, force_mi_collect, process_rss_bytes};
use rstest::rstest;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

const MIB: u64 = 1 << 20;

/// The cases measure process RSS, so they must not overlap (`cargo test`
/// runs a binary's tests on concurrent threads; nextest runs each alone).
static SERIAL: Mutex<()> = Mutex::new(());

/// Where the arena's last reference drops, and who collects.
#[derive(Clone, Copy, Debug)]
enum Release {
    /// On the thread that allocated it (`ArenaPool::release` collects there).
    AllocatingThread,
    /// On another thread; the allocating thread then calls
    /// `force_mi_collect`, as the merge (which shares the sort's coordination
    /// thread with `ReadBlocks`, the arena's allocator) does before it builds
    /// a range plan. mimalloc frees a block of another thread's page only
    /// when that thread collects, so `release`'s own collect cannot return it.
    OtherThreadThenOwnerCollects,
    /// On another thread, and nobody collects on the allocating thread
    /// (the dropping thread's own collect, in `ArenaPool::release` for a
    /// retired pool, still runs).
    OtherThread,
}

/// Fill a 256 MiB pooled arena with 240 MiB and release it as `release` says
/// (after retiring the pool when `retire`). A retired arena dropped on
/// another thread keeps its pages until the allocating thread collects:
/// the dropping thread's collect cannot return them, which is why the range
/// merge collects on the coordination thread before it builds its plan.
/// Columns: whether the pool is retired, the release, whether RSS falls by
/// at least 200 MiB at once.
#[rstest]
#[case::retired_release_returns_the_pages_at_once(true, Release::AllocatingThread, true)]
#[case::retired_cross_thread_release_returns_them_at_the_owners_collect(
    true,
    Release::OtherThreadThenOwnerCollects,
    true
)]
#[case::retired_cross_thread_release_waits_for_the_owner(true, Release::OtherThread, false)]
#[case::pooled_release_keeps_the_pages(false, Release::AllocatingThread, false)]
#[case::pooled_cross_thread_release_keeps_the_pages(false, Release::OtherThread, false)]
fn arena_release_rss(#[case] retire: bool, #[case] release: Release, #[case] returns_pages: bool) {
    let _serial = SERIAL.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    let pool = ArenaPool::new(1, 256 << 20);
    let mut arena = pool.try_acquire().expect("a fresh pool has an arena");
    arena.reserve_full_capacity();
    let chunk = vec![0xA5u8; 1 << 20];
    for _ in 0..240 {
        arena.extend_from_slice(&chunk);
    }
    drop(chunk);
    let arena = Arc::new(arena);
    let filled = process_rss_bytes().expect("RSS is sampled on Linux and macOS");
    if retire {
        pool.retire();
    }
    match release {
        Release::AllocatingThread => drop(arena),
        Release::OtherThreadThenOwnerCollects | Release::OtherThread => {
            std::thread::spawn(move || drop(arena)).join().expect("the dropping thread");
        }
    }
    if matches!(release, Release::OtherThreadThenOwnerCollects) {
        force_mi_collect();
    }
    let after = process_rss_bytes().expect("RSS is sampled on Linux and macOS");
    let fell = filled.saturating_sub(after);
    assert_eq!(
        fell >= 200 * MIB,
        returns_pages,
        "RSS {} MiB filled, {} MiB after the release (fell {} MiB)",
        filled / MIB,
        after / MIB,
        fell / MIB
    );
}
