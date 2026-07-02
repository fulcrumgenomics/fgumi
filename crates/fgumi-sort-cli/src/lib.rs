//! Standalone, framework-light `fgumi sort` command.
//!
//! This crate carries the `fgumi sort` CLI surface (`Sort`, `SortOptions`,
//! `SortOrderArg`), its `--verify` read-and-check path, and the
//! `IndexBamFinalizeHook` BAI indexer. The option types are reused by `runall`
//! (via `MultiSortOptions`); the execution path for standalone `fgumi sort`
//! lives in the umbrella `fgumi` crate (`commands::sort::execute_sort_command`),
//! which routes through the shared chain builder. Keeping this crate
//! framework-light keeps the dependency graph small (`fgumi-pipeline-core`,
//! `fgumi-pipeline-io`, `fgumi-sort`, `fgumi-bam-io`, `fgumi-sam`,
//! `fgumi-cli-common`, `fgumi-cli-macros`) with no pull on `fgumi-consensus`,
//! `fgumi-umi`, or `fgumi-simd-fastq`.

#![deny(unsafe_code)]

pub mod chains;
pub mod sort;
pub mod version;

pub use sort::{Sort, SortOptions, SortOrderArg, TMP_DIRS_ENV, parse_cell_tag, resolve_tmp_dirs};
