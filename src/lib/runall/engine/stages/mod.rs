//! Concrete [`crate::runall::engine::Stage`] implementations for the
//! Sort → `PositionBatch` → `GroupAssign` → Consensus → Filter chain.

pub mod align_and_merge;
pub mod bgzf_compress;
pub mod coalesce;
pub mod consensus;
pub mod correct;
pub mod count;
pub mod extract;
pub mod fastq_decompress;
pub mod fastq_pair;
pub mod fastq_parse;
pub mod filter;
pub mod group_assign;
pub mod mi_assign;
pub mod position_batch;
pub mod sort;
pub mod tofastq;

pub use align_and_merge::AlignAndMerge;
pub use bgzf_compress::BgzfCompress;
pub use consensus::{CallerFactory, ConsensusStage};
pub use correct::{CorrectMetricsAccumulator, CorrectStage, UmiCounts};
pub use count::CountStage;
pub use extract::ExtractStage;
pub use fastq_decompress::BgzfDecompress;
pub use fastq_pair::{FastqBlockMerge, FastqPair};
pub use fastq_parse::{FastqBlockParse, FastqParse};
pub use filter::FilterStage;
pub use group_assign::GroupAssignStage;
pub use mi_assign::MiAssignStage;
pub use position_batch::PositionBatchStage;
pub use sort::SortStage;
pub use tofastq::{ToFastq, bam_records_to_fastq};

// Re-export domain types needed by stage consumers (v2-native).
pub use crate::runall::engine::grouping_types::{MiGroup, MiGroupBatch, PositionGroupBatch};
pub use fgumi_consensus::caller::ConsensusOutput;
