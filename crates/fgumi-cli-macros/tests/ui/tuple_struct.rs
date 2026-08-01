//! `multi_options` requires named fields; a tuple struct must fail the build.
use fgumi_cli_macros::multi_options;

#[multi_options("probe", "Probe Options")]
pub struct TupleOpts(pub u32);

fn main() {}
