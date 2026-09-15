//! Sink steps (BAM file writers, future FASTQ writers).

/// Terminal metrics-recording sink for the `simplex-metrics`/`duplex-metrics`
/// chain path. Gated with the metrics recording machinery it consumes.
#[cfg(feature = "consensus")]
pub(crate) mod metrics;
pub mod write_bgzf;
pub mod write_raw;
