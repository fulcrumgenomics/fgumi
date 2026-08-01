//! The companion struct and both conversion impls are emitted without generic
//! parameters, so a generic options struct would expand into code that cannot
//! compile — with every error naming the type parameter rather than the macro.
//!
//! The bounds are chosen so clap's own derive is satisfied; the only diagnostic
//! left is the macro's.
use fgumi_cli_macros::multi_options;

#[multi_options("probe", "Probe Options")]
#[derive(clap::Args, Debug, Clone)]
pub struct GenericOpts<T>
where
    T: Clone + Send + Sync + 'static + std::str::FromStr,
    <T as std::str::FromStr>::Err: std::error::Error + Send + Sync + 'static,
{
    #[arg(long)]
    pub value: T,
}

fn main() {}
