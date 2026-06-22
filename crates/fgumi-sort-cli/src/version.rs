//! Version string for the `@PG` record emitted by `fgumi sort`.
//!
//! This crate is framework-light and standalone, so its built-in default is
//! the Cargo package version. The umbrella `fgumi` binary derives a richer
//! version string (git commit + dirty flag) via `built`; it installs that
//! richer string once at startup through [`set_version_override`] so the
//! `@PG VN` field on sorted BAMs matches the rest of the `fgumi` toolchain.

use std::sync::OnceLock;

/// Process-global version override. Set at most once by the umbrella binary's
/// `main` so the `@PG` record carries the git-augmented version. When unset
/// (e.g. a standalone `fgumi-sort-cli` consumer), the Cargo package version is
/// used.
static VERSION_OVERRIDE: OnceLock<String> = OnceLock::new();

/// Install the version string used in the sort `@PG` record.
///
/// Idempotent: only the first call takes effect (subsequent calls are ignored),
/// matching the "set once at startup" contract. Returns `true` if this call
/// installed the value, or `false` if a value was already set.
pub fn set_version_override(version: String) -> bool {
    VERSION_OVERRIDE.set(version).is_ok()
}

/// Returns the version string used in the sort `@PG` record.
///
/// Prefers the override installed via [`set_version_override`]; otherwise falls
/// back to the Cargo package version.
#[must_use]
pub fn version_string() -> String {
    VERSION_OVERRIDE.get().cloned().unwrap_or_else(|| env!("CARGO_PKG_VERSION").to_string())
}
