//! `validate()` enforces required-ness by observing absence, which a defaulted
//! field and a bare `bool` never exhibit — so the requirement would be silently
//! unenforceable on the re-exposed flag.
use fgumi_cli_macros::multi_options;

#[multi_options("probe", "Probe Options")]
#[derive(clap::Args, Debug, Clone)]
pub struct UnenforceableOpts {
    #[arg(long, default_value_t = 3, required = true)]
    pub defaulted: u32,
    #[arg(long, required = true)]
    pub toggled: bool,
}

fn main() {}
