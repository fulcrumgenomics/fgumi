//! `multi_options` generates a companion struct from named fields, so applying it
//! to any other item — an enum here — must fail the build.
use fgumi_cli_macros::multi_options;

#[multi_options("probe", "Probe Options")]
pub enum NotAStruct {
    First,
    Second,
}

fn main() {}
