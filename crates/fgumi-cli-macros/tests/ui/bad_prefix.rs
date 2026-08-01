//! The prefix is spliced into both the flag name (`--<prefix>::<flag>`) and the
//! generated field identifier (`<prefix>_<field>`), so a non-identifier prefix
//! must be rejected at the literal rather than panicking inside `format_ident!`.
use fgumi_cli_macros::multi_options;

#[multi_options("", "Probe Options")]
#[derive(clap::Args, Debug, Clone)]
pub struct EmptyPrefixOpts {
    #[arg(long)]
    pub value: u32,
}

#[multi_options("2fast", "Probe Options")]
#[derive(clap::Args, Debug, Clone)]
pub struct LeadingDigitOpts {
    #[arg(long)]
    pub value: u32,
}

fn main() {}
