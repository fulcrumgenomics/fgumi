//! Property tests exercising narrow but load-bearing invariants.
//!
//! Each sub-module targets a single property:
//! - [`umi_adjacency`]: Adjacency UMI clustering is deterministic across runs.
//! - [`sort_order`]: The template-coordinate sort key induces a total order.

mod sort_order;
mod umi_adjacency;
