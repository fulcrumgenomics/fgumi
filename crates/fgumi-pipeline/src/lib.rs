//! Pipeline crate: parallel Zone 3 stages, work-stealing scheduler, and pipeline builder.

#![deny(unsafe_code)]

pub mod pipeline;
pub mod stage;
pub mod stages;

pub use stage::{PipelineStage, SequencedItem};
