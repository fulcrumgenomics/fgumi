//! Cgroup-aware system resource detection.
//!
//! The authoritative implementations live in `fgumi-cli-common`; this module
//! re-exports them so existing call-site paths (`crate::system::detect_total_memory`,
//! `crate::system::detect_cpu_count`) continue to resolve.
pub use fgumi_cli_common::{detect_cpu_count, detect_total_memory};
