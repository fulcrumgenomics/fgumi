//! Command trait definition for CLI commands.
//!
//! Re-export shim. The [`Command`] trait that all fgumi CLI commands implement
//! now lives in the `fgumi-cli-common` crate. Re-exporting it here keeps
//! `crate::commands::command::Command` resolving AND ensures there is ONE
//! canonical `Command` trait shared with `fgumi-sort-cli` (so the umbrella's
//! `enum_dispatch` `Commands` enum and `fgumi-sort-cli`'s `Sort` impl agree on
//! the same trait).
pub use fgumi_cli_common::Command;
