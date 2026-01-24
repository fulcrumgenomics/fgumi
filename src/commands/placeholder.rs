//! Placeholder command for project skeleton
//!
//! This is a temporary command that exists only to satisfy enum_dispatch's
//! requirement that the Subcommand enum be non-empty. It will be removed
//! once the first real command (correct) is added in Phase 1.
//!
//! The command is hidden from help output via `#[clap(hide = true)]`.

use anyhow::Result;
use clap::Parser;

use super::command::Command;

/// Placeholder command - will be removed when first real command is added.
///
/// This exists because enum_dispatch requires at least one variant in the
/// Subcommand enum. It's hidden from users via clap's hide attribute.
#[derive(Parser, Debug)]
#[command(about = "Internal placeholder - not for use")]
pub struct Placeholder;

impl Command for Placeholder {
    fn execute(&self, _command_line: &str) -> Result<()> {
        // This should never be called since the command is hidden
        anyhow::bail!(
            "Placeholder command should not be executed. Use --help to see available commands."
        )
    }
}
