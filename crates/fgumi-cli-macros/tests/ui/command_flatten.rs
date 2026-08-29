//! Re-prefixing a clap-flattened struct's fields is out of scope, so
//! `#[command(flatten)]` must fail the build rather than silently drop the nesting.
use fgumi_cli_macros::multi_options;

#[derive(clap::Args, Debug, Clone)]
pub struct Inner {
    #[arg(long)]
    pub inner_value: u32,
}

#[multi_options("probe", "Probe Options")]
#[derive(clap::Args, Debug, Clone)]
pub struct FlattenOpts {
    #[command(flatten)]
    pub inner: Inner,
}

fn main() {}
