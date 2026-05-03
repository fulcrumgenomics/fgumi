//! Chain runner implementations for `fgumi runall`.
//!
//! PR 1 exposed [`NotImplementedRunner`] — the fallback that matches
//! every `(start_from, stop_after)` shape and bails with "not yet
//! implemented". PR 2 adds [`ExtractChainRunner`] — the
//! `(Extract, Extract)` chain on top of `unified_pipeline::fastq`. PR 3
//! will add `correct.rs` and `extract_correct.rs`. Each new chain
//! runner is registered ahead of `NotImplementedRunner` in
//! `Runall::build_chain_registry`; the fallback keeps handling shapes
//! that haven't been implemented yet.

pub(crate) mod extract;
pub(crate) mod not_implemented;

pub(crate) use extract::ExtractChainRunner;
pub(crate) use not_implemented::NotImplementedRunner;
