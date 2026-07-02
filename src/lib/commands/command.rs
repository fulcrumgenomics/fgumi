//! Command trait definition for CLI commands.
//!
//! Re-export shim. The [`Command`] trait that the fgumi CLI commands implement
//! now lives in the `fgumi-cli-common` crate. Re-exporting it here keeps
//! `crate::commands::command::Command` resolving AND ensures there is ONE
//! canonical `Command` trait shared across the crates that define commands.
//! (`fgumi sort` is the exception: `Sort` is a pure args struct and is
//! dispatched directly to `commands::sort::execute_sort_command`, so it does
//! not implement this trait.)
pub use fgumi_cli_common::Command;
