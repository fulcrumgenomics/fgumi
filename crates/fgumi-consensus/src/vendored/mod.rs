// Vendored BAM codec from noodles — see bam_codec/mod.rs for details on why this is vendored.
// Allow clippy lints and dead code in vendored modules.
#[allow(dead_code, clippy::all, clippy::pedantic, clippy::nursery)]
pub(crate) mod bam_codec;
