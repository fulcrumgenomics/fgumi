//! Metrics collection and reporting for fgumi operations.
//!
//! This module provides structured metric types organized by command:
//! - [`clip`] - Read clipping metrics
//! - [`group`] - UMI grouping metrics
//! - [`consensus`] - Consensus calling metrics (shared by simplex/duplex/codec)
//! - [`duplex`] - Duplex sequencing QC metrics
//! - [`correct`] - UMI correction metrics
//! - [`writer`] - Metrics file I/O utilities
//!
//! # Traits
//!
//! - [`Metric`] - Core trait for serializable metrics
//! - [`ProcessingMetrics`] - Common interface for input/output metrics

// Re-export core items from fgumi-metrics
pub use fgumi_metrics::{FLOAT_PRECISION, Metric, ProcessingMetrics, format_float};

// Re-export submodules for path compatibility (e.g. fgumi_lib::metrics::consensus::ConsensusMetrics)
pub use fgumi_metrics::clip;
pub use fgumi_metrics::consensus;
pub use fgumi_metrics::correct;
pub use fgumi_metrics::duplex;
pub use fgumi_metrics::group;
pub use fgumi_metrics::writer;

// Re-export commonly used types
pub use clip::{ClippingMetrics, ClippingMetricsCollection, ReadType};
pub use consensus::{ConsensusKvMetric, ConsensusMetrics};
pub use correct::UmiCorrectionMetrics;
pub use duplex::{
    DuplexFamilySizeMetric, DuplexMetricsCollector, DuplexUmiMetric, DuplexYieldMetric,
    FamilySizeMetric, UmiMetric,
};
pub use group::{FamilySizeMetrics, UmiGroupingMetrics};
pub use writer::write_metrics;
