//! The generated companion inherits the annotated struct's visibility, so a
//! private options struct must not produce a `pub` companion that escapes its
//! module.
mod inner {
    use fgumi_cli_macros::multi_options;

    #[multi_options("probe", "Probe Options")]
    #[derive(clap::Args, Debug, Clone)]
    struct PrivateOpts {
        #[arg(long, default_value_t = 1)]
        value: u32,
    }
}

fn main() {
    let _ = inner::MultiPrivateOpts { probe_value: 1 };
}
