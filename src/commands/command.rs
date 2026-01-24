//! Command trait definition for CLI commands.
//!
//! This module defines the [`Command`] trait that all fgumi CLI commands implement.
//! The trait uses `enum_dispatch` for efficient dynamic dispatch across command variants.

use anyhow::Result;
use enum_dispatch::enum_dispatch;

/// Trait implemented by all fgumi CLI commands.
///
/// Each command provides an `execute` method that runs the command's main logic.
/// The `command_line` parameter contains the full command invocation for @PG records.
#[enum_dispatch]
pub trait Command {
    #[allow(clippy::missing_errors_doc)]
    fn execute(&self, command_line: &str) -> Result<()>;
}
