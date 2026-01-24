//! CLI command implementations for fgumi.
//!
//! This module contains all the command implementations for the fgumi CLI tool.
//! Each submodule implements a specific command (extract, group, filter, etc.).
//!
//! # Command Categories
//!
//! ## Pre-alignment
//! - [`extract`] - Extract UMIs from FASTQ and create unmapped BAM
//!
//! ## Post-alignment
//! - [`group`] - Group reads by UMI into molecular families
//! - [`correct`] - Correct UMI sequences against a fixed list
//! - [`filter`] - Filter reads based on various criteria
//!
//! ## Consensus Calling
//! - [`simplex`] - Call single-strand molecular consensus
//! - [`duplex`] - Call duplex (double-strand) molecular consensus
//! - [`codec`] - CODEC-style consensus calling
//!
//! ## Utilities
//! - [`clip`] - Clip overlapping read pairs
//! - [`duplex_metrics`] - Calculate duplex sequencing metrics

#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_possible_wrap)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_lossless)]
#![allow(clippy::missing_errors_doc)]
#![allow(clippy::missing_panics_doc)]
#![allow(clippy::unused_self)]
#![allow(clippy::unnecessary_wraps)]
#![allow(clippy::similar_names)]
#![allow(clippy::needless_pass_by_value)]
#![allow(clippy::match_same_arms)]
#![allow(clippy::must_use_candidate)]
#![allow(clippy::items_after_statements)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::fn_params_excessive_bools)]
#![allow(clippy::redundant_else)]
#![allow(clippy::manual_let_else)]
#![allow(clippy::needless_continue)]

pub mod clip;
pub mod codec;
pub mod command;
pub mod common;
#[cfg(feature = "compare")]
pub mod compare;
pub mod consensus_runner;
pub mod correct;
pub mod downsample;
pub mod duplex;
pub mod duplex_metrics;
pub mod extract;
pub mod fastq;
pub mod filter;
pub mod group;
pub mod review;
pub mod simplex;
#[cfg(feature = "simulate")]
pub mod simulate;
pub mod sort;
pub mod zipper;
