//! Aligner I/O mode selection and dedicated I/O thread management.
//!
//! Two modes:
//! - [`BatchCoordinated`](AlignerIoMode::BatchCoordinated): uses
//!   [`AlignerBatchState`](super::batch_state::AlignerBatchState) for `-K` flow control.
//!   Two lightweight I/O threads (stdin-writer, stdout-reader) handle blocking pipe I/O.
//! - [`DedicatedThreads`](AlignerIoMode::DedicatedThreads): no batch coordination.
//!   Two I/O threads handle all pipe I/O.
//!
//! In both modes, I/O threads do NOT count against `--threads`.

/// How the pipeline communicates with the aligner subprocess.
#[derive(Debug, Clone)]
pub enum AlignerIoMode {
    /// Batch-coordinated with bwa `-K` flow control.
    BatchCoordinated {
        /// bwa's `-K` value in bases.
        k_bases: u64,
    },
    /// Dedicated I/O threads without batch coordination.
    /// Used for non-bwa aligners or when `-K` is not specified.
    DedicatedThreads,
}

impl AlignerIoMode {
    /// Default `-K` value matching bwa mem's common usage.
    pub const DEFAULT_K_BASES: u64 = 150_000_000;

    /// Create a batch-coordinated mode with the given `-K` value.
    #[must_use]
    pub fn batch_coordinated(k_bases: u64) -> Self {
        Self::BatchCoordinated { k_bases }
    }

    /// Create a dedicated-threads mode (no batch coordination).
    #[must_use]
    pub fn dedicated() -> Self {
        Self::DedicatedThreads
    }

    /// Whether this mode uses batch coordination.
    #[must_use]
    pub fn is_batch_coordinated(&self) -> bool {
        matches!(self, Self::BatchCoordinated { .. })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_batch_coordinated_mode() {
        let mode = AlignerIoMode::batch_coordinated(100_000);
        assert!(mode.is_batch_coordinated());
        match mode {
            AlignerIoMode::BatchCoordinated { k_bases } => assert_eq!(k_bases, 100_000),
            AlignerIoMode::DedicatedThreads => panic!("expected BatchCoordinated"),
        }
    }

    #[test]
    fn test_dedicated_mode() {
        let mode = AlignerIoMode::dedicated();
        assert!(!mode.is_batch_coordinated());
        assert!(matches!(mode, AlignerIoMode::DedicatedThreads));
    }

    #[test]
    fn test_default_k_bases() {
        assert_eq!(AlignerIoMode::DEFAULT_K_BASES, 150_000_000);
        let mode = AlignerIoMode::batch_coordinated(AlignerIoMode::DEFAULT_K_BASES);
        match mode {
            AlignerIoMode::BatchCoordinated { k_bases } => {
                assert_eq!(k_bases, 150_000_000);
            }
            AlignerIoMode::DedicatedThreads => panic!("expected BatchCoordinated"),
        }
    }

    #[test]
    fn test_clone_and_debug() {
        let mode = AlignerIoMode::batch_coordinated(1000);
        let cloned = mode.clone();
        assert!(cloned.is_batch_coordinated());
        // Ensure Debug is implemented (compile-time check, runtime smoke test).
        let debug_str = format!("{mode:?}");
        assert!(debug_str.contains("BatchCoordinated"));
    }
}
