//! A positional argument has no flag name to prefix, and clap panics outright when
//! one carries a `long` — which the generated companion always emits. Several
//! stages flattened into one runall command would have ambiguous positionals
//! anyway, so both spellings must fail the build.
use fgumi_cli_macros::multi_options;

#[multi_options("probe", "Probe Options")]
#[derive(clap::Args, Debug, Clone)]
pub struct PositionalOpts {
    pub implicit_positional: u32,
    #[arg(index = 2)]
    pub explicit_positional: u32,
}

fn main() {}
