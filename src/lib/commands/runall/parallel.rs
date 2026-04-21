//! Filter and group-filter configuration builders for `runall`.

use std::sync::Arc;

use anyhow::{Result, anyhow};

use crate::commands::filter::CollectedFilterMetrics;
use crate::consensus::filter::FilterConfig;
use crate::grouper::GroupFilterConfig;
use crate::per_thread_accumulator::PerThreadAccumulator;

use super::Runall;
use super::plan::StageSpec;
use super::planner::FilterStatsState;

impl Runall {
    /// Build a [`FilterConfig`] from the CLI filter options.
    ///
    /// Uses the same defaults as the standalone `fgumi filter` command for
    /// error rate and no-call fraction thresholds.
    ///
    /// # Errors
    ///
    /// Returns an error if `--filter::min-reads` was not supplied.
    /// `Runall::validate` rejects this case before reaching the planner, so
    /// production call sites surface the error at plan-build time rather than
    /// here. Test code calling this directly must supply `--filter::min-reads`.
    pub fn build_filter_config(&self) -> Result<FilterConfig> {
        let min_reads = self.filter_opts.filter_min_reads.as_deref().unwrap_or(&[]);
        if min_reads.is_empty() {
            return Err(anyhow!(
                "--filter::min-reads is required when the plan includes the \
                 consensus or filter stage"
            ));
        }
        if min_reads.len() > 3 {
            return Err(anyhow!(
                "--filter::min-reads must have 1-3 values, got {}",
                min_reads.len()
            ));
        }
        Ok(self.filter_config_with_min_reads(min_reads))
    }

    /// Plan-time filter config: returns a placeholder for `--explain` where
    /// `--filter::min-reads` is bypassed, otherwise requires it to be set.
    pub(crate) fn filter_config_for_plan(&self, explain: bool) -> Result<FilterConfig> {
        let min_reads = self.filter_opts.filter_min_reads.as_deref().unwrap_or(&[]);
        if explain && min_reads.is_empty() {
            return Ok(self.filter_config_with_min_reads(&[0]));
        }
        self.build_filter_config()
    }

    /// Build the filter [`StageSpec`] plus an optional [`FilterStatsState`]
    /// to finalize. When `--filter::stats` is set, allocates a per-thread
    /// accumulator sized to `threads`; the same `Arc` is handed to the stage
    /// (for writes) and the state (for the post-run fold + TSV write).
    pub(crate) fn filter_stage_spec_with_stats(
        &self,
        explain: bool,
    ) -> Result<(StageSpec, Option<FilterStatsState>)> {
        let config = self.filter_config_for_plan(explain)?;
        let (metrics, state) = match self.filter_opts.filter_stats.as_ref() {
            Some(path) => {
                let acc = PerThreadAccumulator::<CollectedFilterMetrics>::new(self.threads.max(1));
                let state = FilterStatsState { metrics: Arc::clone(&acc), path: path.clone() };
                (Some(acc), Some(state))
            }
            None => (None, None),
        };
        let spec = StageSpec::Filter {
            config,
            filter_by_template: self.filter_opts.filter_filter_by_template,
            metrics,
        };
        Ok((spec, state))
    }

    /// Build a [`FilterConfig`] with an explicit `min_reads` slice. Delegates
    /// to [`FilterConfig::new`], which fans 1-3 values out to
    /// `[total, AB, BA]` thresholds — matching standalone `fgumi filter`'s
    /// handling of `--min-reads`, `--max-read-error-rate`, and
    /// `--max-base-error-rate`.
    pub(crate) fn filter_config_with_min_reads(&self, min_reads: &[usize]) -> FilterConfig {
        FilterConfig::new(
            min_reads,
            &self.filter_opts.filter_max_read_error_rate,
            &self.filter_opts.filter_max_base_error_rate,
            Some(self.filter_opts.filter_min_base_quality),
            self.filter_opts.filter_min_mean_base_quality,
            self.filter_opts.filter_max_no_call_fraction,
        )
    }

    /// Build a [`GroupFilterConfig`] from the CLI group options.
    ///
    /// Controls template filtering before UMI assignment (MAPQ, unmapped,
    /// non-PF, UMI validation), matching the standalone `fgumi group`
    /// command.
    pub fn build_group_filter_config(&self) -> GroupFilterConfig {
        let umi_tag: [u8; 2] = *crate::sam::SamTag::RX;
        GroupFilterConfig {
            umi_tag,
            min_mapq: self.group_opts.group_min_mapq,
            include_non_pf: false,
            min_umi_length: None,
            no_umi: false,
            allow_unmapped: self.group_opts.group_allow_unmapped,
        }
    }
}
