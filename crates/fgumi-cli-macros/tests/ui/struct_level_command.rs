//! Struct-level clap configuration is not carried onto the generated companion,
//! so it must fail the build rather than silently apply to only one of the two
//! commands.
use fgumi_cli_macros::multi_options;

#[multi_options("probe", "Probe Options")]
#[derive(clap::Args, Debug, Clone)]
#[command(next_help_heading = "Somewhere Else")]
pub struct StructLevelOpts {
    #[arg(long)]
    pub value: u32,
}

fn main() {}
