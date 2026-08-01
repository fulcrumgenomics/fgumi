//! clap's `key(value)` call form arrives as `Meta::List` and slips past every classifier,
//! so it must be rejected for the keys the macro rewrites.
use fgumi_cli_macros::multi_options;

#[multi_options("probe", "Probe Options")]
#[derive(clap::Args, Debug, Clone)]
pub struct CallFormOpts {
    #[arg(long("renamed"))]
    pub value: u32,
}

fn main() {}
