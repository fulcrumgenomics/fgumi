//! Per-thread auto-tuned defaults for the BAM step library, preserving the
//! per-step batch-size and thread-count heuristics from the legacy
//! `PipelineConfig::auto_tuned` (removed in the issue #330 migration).
//!
//! Centralizes the per-step batch sizes and thread-count heuristics so
//! command builders don't have to know the magic numbers. Phase 4
//! migrations should construct via `BamPipelineTuning::auto_tuned(threads)`
//! and pass the resulting struct into each step's constructor.

/// Per-thread auto-tuned defaults for the BAM step library.
///
/// Mirrors legacy `PipelineConfig::auto_tuned`: `blocks_per_read_batch`
/// scales with thread count to reduce downstream queue-op pulse rate;
/// `template_batch_size` defaults to 500 templates/batch (matches
/// legacy's `target_templates_per_batch`).
#[derive(Debug, Clone, Copy)]
pub struct BamPipelineTuning {
    /// Number of pipeline worker threads.
    pub threads: usize,
    /// BGZF blocks per `ReadBgzfBlocks` emission. Matches legacy's
    /// `blocks_per_read_batch`: 16 (1-3 threads), 32 (4-7), 48 (8-15),
    /// 64 (16+).
    pub blocks_per_batch: usize,
    /// Templates per `GroupBam` output batch. Legacy default: 500.
    pub template_batch_size: usize,
    /// Per-step output queue byte budget. 4 MiB per step is a reasonable
    /// laptop-class default; large workloads should raise this via the
    /// `--queue-memory` flag (see `PipelineConfig::queue_memory_total`).
    pub per_step_byte_limit: u64,
    /// BGZF compression level for output (1-12). Default: 1 for the
    /// round-trip preset (fastest); production commands should use 6.
    pub compression_level: u32,
}

impl BamPipelineTuning {
    /// Auto-tuned defaults for `threads` worker threads, mirroring the legacy
    /// `PipelineConfig::auto_tuned` from the now-deprecated unified pipeline.
    /// Deliberately not an intra-doc link: that code is deleted once the
    /// commands migrate off it, and a link would break again then.
    #[must_use]
    pub fn auto_tuned(threads: usize) -> Self {
        let threads = threads.max(1);
        let blocks_per_batch = match threads {
            1..=3 => 16,
            4..=7 => 32,
            8..=15 => 48,
            _ => 64,
        };
        Self {
            threads,
            blocks_per_batch,
            template_batch_size: 500,
            per_step_byte_limit: 4 * 1024 * 1024,
            compression_level: 1,
        }
    }

    /// Override the BGZF compression level (default 1 — fastest).
    /// Production commands typically want 6.
    #[must_use]
    pub fn with_compression_level(mut self, level: u32) -> Self {
        self.compression_level = level;
        self
    }

    /// Override the per-step byte limit. Use this with the global memory
    /// budget from the `--queue-memory` flag.
    #[must_use]
    pub fn with_per_step_byte_limit(mut self, bytes: u64) -> Self {
        self.per_step_byte_limit = bytes;
        self
    }
}

impl Default for BamPipelineTuning {
    fn default() -> Self {
        Self::auto_tuned(4)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn auto_tuned_matches_legacy_block_counts() {
        assert_eq!(BamPipelineTuning::auto_tuned(1).blocks_per_batch, 16);
        assert_eq!(BamPipelineTuning::auto_tuned(3).blocks_per_batch, 16);
        assert_eq!(BamPipelineTuning::auto_tuned(4).blocks_per_batch, 32);
        assert_eq!(BamPipelineTuning::auto_tuned(7).blocks_per_batch, 32);
        assert_eq!(BamPipelineTuning::auto_tuned(8).blocks_per_batch, 48);
        assert_eq!(BamPipelineTuning::auto_tuned(15).blocks_per_batch, 48);
        assert_eq!(BamPipelineTuning::auto_tuned(16).blocks_per_batch, 64);
        assert_eq!(BamPipelineTuning::auto_tuned(64).blocks_per_batch, 64);
    }

    #[test]
    fn template_batch_size_matches_legacy_default() {
        assert_eq!(BamPipelineTuning::auto_tuned(8).template_batch_size, 500);
    }

    #[test]
    fn min_one_thread() {
        assert_eq!(BamPipelineTuning::auto_tuned(0).threads, 1);
    }

    /// The defaults documented on the struct: 4 MiB per step and level 1 (the
    /// round-trip preset). Commands override both, so a silent drift here would
    /// change every chain that does not.
    #[test]
    fn auto_tuned_defaults_match_the_documented_budget_and_level() {
        let tuning = BamPipelineTuning::auto_tuned(8);
        assert_eq!(tuning.per_step_byte_limit, 4 * 1024 * 1024);
        assert_eq!(tuning.compression_level, 1);
    }

    /// The builder-style overrides must change only their own field — a command
    /// setting a production compression level must not also perturb the queue
    /// budget, and vice versa.
    #[test]
    fn the_builder_overrides_change_only_their_own_field() {
        let base = BamPipelineTuning::auto_tuned(8);

        let leveled = base.with_compression_level(6);
        assert_eq!(leveled.compression_level, 6);
        assert_eq!(leveled.per_step_byte_limit, base.per_step_byte_limit);
        assert_eq!(leveled.blocks_per_batch, base.blocks_per_batch);
        assert_eq!(leveled.threads, base.threads);

        let budgeted = base.with_per_step_byte_limit(64 * 1024);
        assert_eq!(budgeted.per_step_byte_limit, 64 * 1024);
        assert_eq!(budgeted.compression_level, base.compression_level);
        assert_eq!(budgeted.blocks_per_batch, base.blocks_per_batch);

        // Chaining applies both.
        let both = base.with_compression_level(6).with_per_step_byte_limit(64 * 1024);
        assert_eq!(both.compression_level, 6);
        assert_eq!(both.per_step_byte_limit, 64 * 1024);
    }

    /// `Default` is documented as `auto_tuned(4)`; pin it so the default queue
    /// budget cannot drift away from the 4-thread tuning without a test failing.
    #[test]
    fn default_is_the_four_thread_tuning() {
        let default = BamPipelineTuning::default();
        let explicit = BamPipelineTuning::auto_tuned(4);
        assert_eq!(default.threads, explicit.threads);
        assert_eq!(default.blocks_per_batch, explicit.blocks_per_batch);
        assert_eq!(default.template_batch_size, explicit.template_batch_size);
        assert_eq!(default.per_step_byte_limit, explicit.per_step_byte_limit);
        assert_eq!(default.compression_level, explicit.compression_level);
    }
}
