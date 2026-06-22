//! Standalone, framework-light `fgumi sort` command.
//!
//! This crate carries the `fgumi sort` CLI command (`Sort`), its options
//! (`SortOptions`, `SortOrderArg`), and the sort-specific finalize hooks and
//! step factory. `Sort::execute` builds the typed-step pipeline directly via
//! [`fgumi_pipeline_core::Pipeline::builder`] rather than going through the
//! umbrella's monolithic `ChainBuilder`, so the dependency graph stays light
//! (`fgumi-pipeline-core`, `fgumi-pipeline-io`, `fgumi-sort`, `fgumi-bam-io`,
//! `fgumi-sam`, `fgumi-cli-common`, `fgumi-cli-macros`) with no pull on
//! `fgumi-consensus`, `fgumi-umi`, or `fgumi-simd-fastq`.

#![deny(unsafe_code)]

pub mod chains;
pub mod sort;
pub mod version;

pub use fgumi_cli_common::Command;
pub use sort::{Sort, SortOptions, SortOrderArg, TMP_DIRS_ENV, parse_cell_tag, resolve_tmp_dirs};
