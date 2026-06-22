//! Custom error types for fgumi operations.
//!
//! The authoritative definitions live in `fgumi-cli-common`; this module
//! re-exports them so existing call-site paths (`crate::errors::FgumiError`,
//! `crate::errors::Result`) continue to resolve.
pub use fgumi_cli_common::{FgumiError, Result};
