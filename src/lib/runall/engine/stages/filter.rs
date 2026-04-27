//! Filter stage: quality filtering, methylation masking, and
//! template-level emit decision for consensus reads.
//!
//! Parallel pool stage. Consumes the framed BAM bytes produced by
//! [`crate::runall::engine::stages::consensus::ConsensusStage`], applies
//! per-base masking and per-read filters, then walks consecutive
//! same-QNAME records into templates and hands them to the shared
//! [`fgumi_consensus::filter::filter_template_records`] orchestration
//! — the same function the standalone `fgumi filter` command runs.
//!
//! ## Input
//!
//! [`ConsensusOutput`] — a framed BAM-record buffer (each record
//! prefixed with a 4-byte little-endian `block_size`) where consensus
//! R1/R2 pairs share a QNAME and appear consecutively.
//!
//! ## Output
//!
//! [`ConsensusOutput`] — surviving framed records in the same format.
//! Per-base masks have been applied in place; per-template filter
//! failures remove all contributing records when
//! `filter_by_template` is `true`. The output `count` field is not
//! recomputed (set to `0`); downstream consumers read `data` directly.
//!
//! ## Ordering guarantees
//!
//! Input record order is preserved within each surviving template.
//! Templates themselves appear in the same order as they were seen in
//! the input buffer.
//!
//! ## Memory model
//!
//! 1:1 per call, with output capacity preallocated to input length.
//! Scratch span-vector `Vec<(usize, usize)>` grows proportional to
//! record count. No persistent per-worker state beyond the
//! immutable configuration.
//!
//! ## Determinism
//!
//! Byte-identical per batch: filtering is a pure function of
//! `(record bytes, FilterConfig, MethylationRef)`. The reference
//! genome lookups (`methylation_ref`) are also deterministic.

use std::sync::Arc;

use anyhow::{Result, anyhow};
use fgumi_consensus::caller::ConsensusOutput;
use fgumi_consensus::filter::{
    FilterConfig, FilterResult, MethylationDepthThresholds, MethylationTags,
    check_conversion_fraction_raw_with_ref_bases_and_tags, compute_read_stats, filter_duplex_read,
    filter_read, filter_template_records, is_duplex_consensus, mask_bases, mask_duplex_bases,
    mask_methylation_depth_duplex_raw_with_tags, mask_methylation_depth_simplex_raw_with_tags,
    mask_strand_methylation_agreement_raw_with_ref_bases_and_tags, resolve_ref_bases_for_record,
};
use fgumi_raw_bam::fields::{aux_data_slice, read_name};

use crate::commands::filter::CollectedFilterMetrics;
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::runall::engine::stage::{Parallelism, Stage};

/// Pipeline `FilterStage`.
///
/// `Clone` is derived so the work-stealing pool can instantiate per-worker
/// copies — the filter owns no shared mutable state, so cloning is a full
/// independent copy.
#[derive(Clone)]
pub struct FilterStage {
    /// Filter thresholds and per-base masking configuration.
    config: FilterConfig,
    /// Require single-strand agreement for duplex base masking.
    require_ss_agreement: bool,
    /// Methylation depth thresholds for EM-Seq filtering.
    methylation_depth_thresholds: Option<MethylationDepthThresholds>,
    /// Whether to mask `CpG` sites where AB/BA strands disagree on methylation (duplex only).
    require_strand_methylation_agreement: bool,
    /// Minimum bisulfite/enzymatic conversion fraction at non-CpG cytosines.
    min_conversion_fraction: Option<f64>,
    /// Methylation mode (EM-Seq, TAPs, or Disabled).
    methylation_mode: fgumi_consensus::MethylationMode,
    /// Reference + contig names for methylation reference-dependent filters.
    methylation_ref: crate::commands::common::MethylationRef,
    /// When true, apply the per-template filter rule: drop *all* records of a
    /// template if any primary record fails. When false, filter each record
    /// independently. Matches the standalone `fgumi filter --filter-by-template`
    /// flag (default `true` on the CLI).
    filter_by_template: bool,
    /// Per-thread filter-stats accumulator. Workers merge per-batch totals
    /// into their own slot; the planner folds slots and writes the TSV at
    /// pipeline end if `--filter::stats` was supplied.
    metrics: Option<Arc<PerThreadAccumulator<CollectedFilterMetrics>>>,
}

impl FilterStage {
    /// Construct a new `FilterStage` with default (disabled) methylation options.
    ///
    /// For callers that need methylation masking, use [`FilterStage::with_methylation`].
    #[must_use]
    pub fn new(config: FilterConfig, filter_by_template: bool) -> Self {
        Self {
            config,
            require_ss_agreement: false,
            methylation_depth_thresholds: None,
            require_strand_methylation_agreement: false,
            min_conversion_fraction: None,
            methylation_mode: fgumi_consensus::MethylationMode::Disabled,
            methylation_ref: None,
            filter_by_template,
            metrics: None,
        }
    }

    /// Attach a per-thread metrics accumulator. Each worker merges its
    /// per-batch `(total, passed, failed, bases_masked)` into its own slot;
    /// the planner folds slots and writes the `--filter::stats` TSV at
    /// pipeline end.
    #[must_use]
    pub fn with_metrics(
        mut self,
        metrics: Arc<PerThreadAccumulator<CollectedFilterMetrics>>,
    ) -> Self {
        self.metrics = Some(metrics);
        self
    }

    /// Construct a new `FilterStage` with full configuration including methylation.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn with_methylation(
        config: FilterConfig,
        require_ss_agreement: bool,
        methylation_depth_thresholds: Option<MethylationDepthThresholds>,
        require_strand_methylation_agreement: bool,
        min_conversion_fraction: Option<f64>,
        methylation_mode: fgumi_consensus::MethylationMode,
        methylation_ref: crate::commands::common::MethylationRef,
        filter_by_template: bool,
    ) -> Self {
        Self {
            config,
            require_ss_agreement,
            methylation_depth_thresholds,
            require_strand_methylation_agreement,
            min_conversion_fraction,
            methylation_mode,
            methylation_ref,
            filter_by_template,
            metrics: None,
        }
    }

    /// Apply masking and filtering to a single raw BAM record.
    ///
    /// Returns `(bases_masked, pass)`. The record bytes may be modified in
    /// place (base masking). `pass = true` iff the record passes every
    /// per-read filter; the template-level rule is applied by the caller via
    /// [`filter_template_records`].
    #[allow(
        clippy::too_many_lines,
        reason = "linear sequence of masking + filter checks; extracting sub-fns \
                  would obscure the (bases_masked, pass) accumulator flow"
    )]
    fn process_record(&self, record: &mut [u8]) -> Result<(u64, bool)> {
        // Determine whether this is a duplex consensus record (presence of aD/bD tags).
        let aux = aux_data_slice(record);
        let is_duplex = is_duplex_consensus(aux);

        // Pre-compute thresholds once to avoid redundant method calls.
        let duplex_t = self.config.duplex_thresholds();
        let ss_t = self.config.effective_single_strand_thresholds();

        // --- Base masking (accumulate count across all masking steps) ---
        let mut bases_masked: u64 = 0;
        if is_duplex {
            if let Some((cc, ab, ba)) = duplex_t {
                bases_masked += mask_duplex_bases(
                    record,
                    cc,
                    ab,
                    ba,
                    self.config.min_base_quality,
                    self.require_ss_agreement,
                )? as u64;
            }
        } else if let Some(thresholds) = ss_t {
            bases_masked += mask_bases(record, thresholds, self.config.min_base_quality)? as u64;
        }

        // --- Methylation masking (EM-Seq) ---
        let needs_methylation = self.methylation_depth_thresholds.is_some()
            || (self.require_strand_methylation_agreement && is_duplex)
            || self.min_conversion_fraction.is_some();
        let methylation_tags =
            if needs_methylation { Some(MethylationTags::from_record(record)) } else { None };

        if let Some(ref thresholds) = self.methylation_depth_thresholds {
            let tags =
                methylation_tags.as_ref().expect("methylation_tags set when thresholds present");
            bases_masked += if is_duplex {
                mask_methylation_depth_duplex_raw_with_tags(record, thresholds, tags)? as u64
            } else {
                mask_methylation_depth_simplex_raw_with_tags(record, thresholds.duplex, tags)?
                    as u64
            };
        }

        // Resolve reference bases for strand agreement and conversion fraction.
        let needs_ref_bases = (self.require_strand_methylation_agreement && is_duplex)
            || self.min_conversion_fraction.is_some();
        let ref_base_map = if needs_ref_bases {
            self.methylation_ref.as_ref().and_then(|(reference, ref_names)| {
                resolve_ref_bases_for_record(record, reference.as_ref(), ref_names)
            })
        } else {
            None
        };

        // Strand methylation agreement masking (duplex only).
        if self.require_strand_methylation_agreement && is_duplex {
            bases_masked += mask_strand_methylation_agreement_raw_with_ref_bases_and_tags(
                record,
                ref_base_map.as_deref(),
                methylation_tags
                    .as_ref()
                    .expect("methylation_tags set when strand agreement enabled"),
            )? as u64;
        }

        // --- Per-read filter ---
        let filter_result = {
            let aux = aux_data_slice(record);
            if is_duplex {
                if let Some((cc, ab, ba)) = duplex_t {
                    filter_duplex_read(aux, cc, ab, ba)?
                } else {
                    FilterResult::Pass
                }
            } else if let Some(thresholds) = ss_t {
                filter_read(aux, thresholds)?
            } else {
                FilterResult::Pass
            }
        };

        if filter_result != FilterResult::Pass {
            return Ok((bases_masked, false));
        }

        // --- Global per-read checks (no-call fraction, mean base quality) ---
        let (no_call_count, mean_base_qual) = compute_read_stats(record);

        let seq_len = fgumi_raw_bam::fields::l_seq(record) as usize;
        if seq_len > 0 {
            #[expect(
                clippy::cast_precision_loss,
                reason = "precision loss is acceptable for no-call fraction comparison"
            )]
            let no_call_fraction = no_call_count as f64 / seq_len as f64;
            if no_call_fraction > self.config.max_no_call_fraction {
                return Ok((bases_masked, false));
            }
        }

        if let Some(min_mean_qual) = self.config.min_mean_base_quality {
            if mean_base_qual < min_mean_qual {
                return Ok((bases_masked, false));
            }
        }

        // Conversion fraction filter (EM-Seq/TAPs read-level).
        if let Some(min_frac) = self.min_conversion_fraction {
            if !check_conversion_fraction_raw_with_ref_bases_and_tags(
                record,
                min_frac,
                ref_base_map.as_deref(),
                methylation_tags
                    .as_ref()
                    .expect("methylation_tags set when conversion fraction enabled"),
                self.methylation_mode,
            ) {
                return Ok((bases_masked, false));
            }
        }

        Ok((bases_masked, true))
    }
}

impl Stage for FilterStage {
    type Input = ConsensusOutput;
    type Output = ConsensusOutput;

    #[allow(
        clippy::too_many_lines,
        clippy::similar_names,
        reason = "linear two-pass scan; qname_i/qname_j + fs_i/fs_j are group-boundary aliases"
    )]
    #[tracing::instrument(name = "filter", skip_all)]
    fn process(&mut self, input: Self::Input, out: &mut dyn FnMut(Self::Output)) -> Result<()> {
        let data = input.data;
        let mut result: Vec<u8> = Vec::with_capacity(data.len());

        // Walk the framed input and collect (frame_start, frame_end) spans.
        // Pre-allocating frame spans here lets us group consecutive same-QNAME
        // records in a single linear pass (consensus callers emit R1+R2 of a
        // molecule contiguously with a shared QNAME) without hashing.
        let mut spans: Vec<(usize, usize)> = Vec::new();
        let mut offset = 0usize;
        while offset < data.len() {
            if offset + 4 > data.len() {
                return Err(anyhow!(
                    "truncated ConsensusOutput: expected 4-byte block_size at offset {offset}"
                ));
            }
            let block_size = u32::from_le_bytes(
                data[offset..offset + 4].try_into().expect("slice is exactly 4 bytes"),
            ) as usize;
            let frame_start = offset;
            offset += 4;
            if offset + block_size > data.len() {
                return Err(anyhow!(
                    "truncated ConsensusOutput: block_size {block_size} extends past buffer end"
                ));
            }
            spans.push((frame_start, offset + block_size));
            offset += block_size;
        }

        // Per-batch metrics merged into a worker slot once at end of process().
        let mut batch_metrics = CollectedFilterMetrics::default();

        // Group consecutive same-QNAME spans into templates and hand each
        // template to the shared [`filter_template_records`] orchestration.
        let mut i = 0usize;
        while i < spans.len() {
            let (fs_i, fe_i) = spans[i];
            let qname_i = read_name(&data[fs_i + 4..fe_i]);
            let mut j = i + 1;
            while j < spans.len() {
                let (fs_j, fe_j) = spans[j];
                let qname_j = read_name(&data[fs_j + 4..fe_j]);
                if qname_j != qname_i {
                    break;
                }
                j += 1;
            }

            let mut template_records: Vec<fgumi_raw_bam::RawRecord> = spans[i..j]
                .iter()
                .map(|&(fs, fe)| fgumi_raw_bam::RawRecord::from(data[fs + 4..fe].to_vec()))
                .collect();

            let outcome = filter_template_records(
                &mut template_records,
                self.filter_by_template,
                |record| self.process_record(record.as_mut_vec().as_mut_slice()),
            )?;

            batch_metrics.total_bases_masked += outcome.bases_masked;
            for (local_idx, keep) in outcome.keep_mask.iter().enumerate() {
                batch_metrics.total_records += 1;
                if *keep {
                    batch_metrics.passed_records += 1;
                    let rec: &[u8] = template_records[local_idx].as_ref();
                    let block_size = u32::try_from(rec.len()).map_err(|_| {
                        anyhow!("filter output record exceeds u32::MAX bytes (len = {})", rec.len())
                    })?;
                    result.extend_from_slice(&block_size.to_le_bytes());
                    result.extend_from_slice(rec);
                } else {
                    batch_metrics.failed_records += 1;
                }
            }

            i = j;
        }

        if let Some(ref metrics) = self.metrics {
            metrics.with_slot(|slot| slot.merge(batch_metrics));
        }

        // Count isn't recomputed (the upstream stage populates it for BGZF
        // writes and downstream consumers care about `data`, not `count`).
        out(ConsensusOutput { data: result, count: 0 });
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }

    fn output_memory_estimate(&self, output: &Self::Output) -> usize {
        output.data.len()
    }

    fn name(&self) -> &'static str {
        "Filter"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_consensus::filter::FilterThresholds;

    fn minimal_filter_config() -> FilterConfig {
        FilterConfig::for_single_strand(
            FilterThresholds { min_reads: 1, max_read_error_rate: 1.0, max_base_error_rate: 1.0 },
            None,
            None,
            1.0,
        )
    }

    #[test]
    fn test_filter_empty_input() {
        let mut stage = FilterStage::new(minimal_filter_config(), true);
        let input = ConsensusOutput { data: vec![], count: 0 };
        let mut captured = None;
        stage
            .process(input, &mut |v| {
                captured = Some(v);
            })
            .unwrap();
        let out = captured.expect("stage must emit");
        assert_eq!(out.data.len(), 0);
        assert_eq!(stage.name(), "Filter");
        assert_eq!(stage.parallelism(), Parallelism::Parallel);
    }
}
