#![allow(clippy::doc_markdown)] // Generated file contains OPT_LEVEL without backticks

use std::sync::LazyLock;

include!(concat!(env!("OUT_DIR"), "/built.rs"));

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
