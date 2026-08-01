//! `requires` names another argument by its unprefixed id, which dangles once the Multi
//! struct renames fields, so it must fail the build.
use fgumi_cli_macros::multi_options;

#[multi_options("probe", "Probe Options")]
#[derive(clap::Args, Debug, Clone)]
pub struct CrossRefOpts {
    #[arg(long)]
    pub first: u32,
    #[arg(long, requires = "first")]
    pub second: u32,
}

fn main() {}
