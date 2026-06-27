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

#[cfg(test)]
mod tests {
    use super::*;

    /// The override is a process-global `OnceLock`, so its first-write-wins
    /// contract and the `version_string` fallback must be exercised in ONE test
    /// (a second test could observe state left by the first). This test must be
    /// the only thing in the crate's test binary that sets the override.
    #[test]
    fn set_version_override_is_first_write_wins_and_drives_version_string() {
        // Before any override is installed, `version_string` falls back to the
        // Cargo package version. This test is the only thing in the crate's test
        // binary that touches the override, so the fallback is still observable
        // here — assert it so a broken fallback is caught in the same one-shot test.
        assert_eq!(version_string(), env!("CARGO_PKG_VERSION"));

        // First write installs the value and reports success.
        assert!(set_version_override("X".to_string()), "first set must install the value");
        // A second write is ignored (first-write-wins) and reports failure.
        assert!(
            !set_version_override("Y".to_string()),
            "second set must be ignored, leaving the first value"
        );
        // `version_string` reflects the first-installed override, not the Cargo
        // fallback nor the second (ignored) value.
        assert_eq!(version_string(), "X");
    }
}
