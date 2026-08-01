//! The macro classifies fields from their literal `#[arg(...)]` attributes and
//! cannot evaluate a cfg predicate, so a clap attribute hidden behind
//! `#[cfg_attr]` would be ignored and the field silently misclassified.
use fgumi_cli_macros::multi_options;

#[multi_options("probe", "Probe Options")]
#[derive(clap::Args, Debug, Clone)]
pub struct ConditionalOpts {
    #[cfg_attr(unix, arg(long, default_value_t = 3))]
    pub knob: u32,
}

fn main() {}
