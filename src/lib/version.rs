#![allow(clippy::doc_markdown)] // Generated file contains OPT_LEVEL without backticks

use std::sync::LazyLock;

include!(concat!(env!("OUT_DIR"), "/built.rs"));

/// The name of the global allocator used at runtime.
pub const ALLOCATOR: &str = "mimalloc";

/// Version of the software including
/// - Git commit hash
/// - Git dirty info (whether the repo had uncommitted changes)
/// - Cargo package version if no git info found
pub static VERSION: LazyLock<String> = LazyLock::new(|| {
    let prefix = if let Some(s) = GIT_COMMIT_HASH {
        format!("{PKG_VERSION}-{s}")
    } else {
        // This shouldn't happen
        PKG_VERSION.to_string()
    };
    let suffix = match GIT_DIRTY {
        Some(true) => "-dirty",
        _ => "",
    };
    format!("{prefix}{suffix}")
});

/// Multi-line extended version string printed by `fgumi --version`.
///
/// Includes crate version, git commit (short SHA + dirty flag), build timestamp,
/// rust compiler version, target triple, enabled feature flags, and allocator.
/// Invaluable for bug reports since it pins down the exact binary in use.
pub static LONG_VERSION: LazyLock<String> = LazyLock::new(|| {
    let commit = match (GIT_COMMIT_HASH_SHORT, GIT_DIRTY) {
        (Some(sha), Some(true)) => format!("{sha} (dirty)"),
        (Some(sha), Some(false)) => format!("{sha} (clean)"),
        (Some(sha), None) => sha.to_string(),
        (None, _) => "unknown".to_string(),
    };
    let build_date = BUILT_TIME_UTC;
    let features =
        if FEATURES_LOWERCASE_STR.is_empty() { "(none)" } else { FEATURES_LOWERCASE_STR };
    // Clap's `--version` handling prints `{bin_name} {long_version}`, so we start
    // with just the package version (no leading package name) to avoid "fgumi fgumi".
    let rustc = RUSTC_VERSION;
    let target = TARGET;
    let profile = PROFILE;
    let version = PKG_VERSION;
    let allocator = ALLOCATOR;
    format!(
        "{version}\n\
         commit:     {commit}\n\
         built:      {build_date}\n\
         rustc:      {rustc}\n\
         target:     {target}\n\
         profile:    {profile}\n\
         features:   {features}\n\
         allocator:  {allocator}"
    )
});

/// Returns a `&'static str` view of [`LONG_VERSION`] suitable for use as a clap
/// `long_version` attribute. Leaks the string once; negligible memory cost.
#[must_use]
pub fn long_version_str() -> &'static str {
    static LEAKED: LazyLock<&'static str> =
        LazyLock::new(|| Box::leak(LONG_VERSION.clone().into_boxed_str()));
    &LEAKED
}
