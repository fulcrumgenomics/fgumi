//! Vendored code that will be removed once upstream PRs are merged.
//!
//! This module contains temporary implementations that are waiting on:
//! - noodles-bgzf: `MultithreadedWriter` with position tracking (<https://github.com/zaeleus/noodles/pull/371>)
//!
//! The remaining vendored modules (`bgzf_multithreaded`) should be removed once their
//! respective upstream PRs are merged and released.

pub mod bgzf_multithreaded;

// Re-export MultithreadedWriter types
pub use bgzf_multithreaded::{
    BlockInfo, BlockInfoRx, Builder as MultithreadedWriterBuilder, MultithreadedWriter,
};
