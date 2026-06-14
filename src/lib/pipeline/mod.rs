//! The typed-step pipeline framework for `--threads N` mode.
//!
//! All multi-threaded commands run on the **typed-step framework** in
//! [`core`] (execution engine) and [`chains`] (declarative chain
//! construction). A pipeline is a chain of [`core::step::Step`]s, each
//! declaring a [`core::step::StepKind`] (`Serial` / `Parallel` /
//! `Exclusive`); the engine runs all `N` worker threads through a
//! round-robin dispatch loop with `Affinity`-based pinning for I/O
//! sources/sinks. This keeps `--threads N` a strict thread cap with no
//! separate I/O thread pools.
//!
//! A BAM pipeline's canonical step sequence is:
//!
//! ```text
//! Read → Decompress → FindBoundaries → Decode → Group → Process → Serialize → Compress → Write
//! ```
//!
//! # Module structure
//!
//! - [`core`]: typed-step execution engine (worker loop, queues, drain
//!   protocol, reorder buffers).
//! - [`chains`]: declarative chain construction (`build_for`, `ChainBuilder`,
//!   per-command step factories).
//! - [`steps`]: the concrete `Step` implementations (read, decompress,
//!   boundaries, decode, group, process, serialize, compress, write, …).
//!
//! Shared grouping/decoded-record domain types ([`fgumi_bam_io::GroupKey`],
//! [`fgumi_bam_io::DecodedRecord`], [`fgumi_bam_io::Grouper`], …) live in the
//! `fgumi-bam-io` crate, next to the raw-record helpers they operate on.
//!
//! # Adding a new command
//!
//! Implement the per-stage steps as [`core::step::Step`]s under [`steps`],
//! then wire them into a chain via [`chains`] (one `Stage` variant, one bag
//! slot, one validator entry, one `add_<stage>` method, one factory per
//! step). See the project `CLAUDE.md` "Chain-builder façade" notes for the
//! full checklist.

pub mod chains;
pub mod core;
pub mod steps;
