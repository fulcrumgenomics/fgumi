//! Every classifier keys on `#[arg(...)]`, so clap's legacy `#[clap(...)]`
//! spelling would be silently ignored — a `#[clap(skip)]` field would be exposed
//! as a required CLI flag. Reject it instead.
use fgumi_cli_macros::multi_options;

#[multi_options("probe", "Probe Options")]
#[derive(clap::Args, Debug, Clone)]
pub struct LegacyOpts {
    #[clap(skip)]
    pub hidden: u32,
    #[arg(long)]
    pub value: u32,
}

fn main() {}
