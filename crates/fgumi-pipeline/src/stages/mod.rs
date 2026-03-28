//! Stage implementations for the runall pipeline.

pub mod aligner;
pub mod compress;
pub mod consensus;
pub mod filter;
pub mod group_assign;
pub mod position_batch;

pub use aligner::{AlignerPreset, AlignerProcess};
pub use compress::{CompressStage, CompressedBatch};
pub use consensus::{CallerFactory, ConsensusStage};
pub use filter::FilterStage;
pub use group_assign::{GroupAndAssignStage, MiGroup, PositionGroupBatch};
pub use position_batch::PositionBatcher;
