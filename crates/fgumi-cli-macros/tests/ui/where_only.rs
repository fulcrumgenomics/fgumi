use fgumi_cli_macros::multi_options;

#[multi_options("probe", "Probe Options")]
#[derive(clap::Args, Debug, Clone)]
pub struct WhereOnlyOpts
where
    u32: Clone,
{
    #[arg(long)]
    pub value: u32,
}

fn main() {}
