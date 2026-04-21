//! Observation hooks used by gate tests to verify that the pool spawns the
//! expected set of threads and no per-stage worker threads.
//!
//! Always compiled because the overhead is negligible: a single lazily
//! allocated `Mutex<Vec<String>>`. Production builds touch this registry only
//! on thread spawn and thread exit, which is not on any hot path.
//!
//! Pool workers call [`register_thread`] on startup and [`unregister_thread`]
//! on clean exit. The `pipeline::source` and `pipeline::sink` helper
//! threads do the same under their well-known names so the gate tests can
//! enumerate and assert on the live thread set.

use std::sync::Mutex;

static THREAD_NAMES: Mutex<Vec<String>> = Mutex::new(Vec::new());

/// Record `name` as a currently live thread.
///
/// # Panics
/// Panics if the thread-name registry mutex is poisoned, which would mean a
/// previous writer panicked while holding the lock.
pub fn register_thread(name: &str) {
    THREAD_NAMES.lock().unwrap().push(name.to_string());
}

/// Remove the first entry matching `name`.
///
/// # Panics
/// Panics if the thread-name registry mutex is poisoned.
pub fn unregister_thread(name: &str) {
    let mut g = THREAD_NAMES.lock().unwrap();
    if let Some(pos) = g.iter().position(|n| n == name) {
        g.swap_remove(pos);
    }
}

/// Snapshot the current set of live thread names.
///
/// # Panics
/// Panics if the thread-name registry mutex is poisoned.
#[must_use]
pub fn snapshot_threads() -> Vec<String> {
    THREAD_NAMES.lock().unwrap().clone()
}

/// Clear the registry (tests call this in setup to avoid cross-pollution).
///
/// # Panics
/// Panics if the thread-name registry mutex is poisoned.
pub fn clear_threads() {
    THREAD_NAMES.lock().unwrap().clear();
}
