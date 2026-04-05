//! Pipeline entry point: run from grouped input (MI tags already assigned).

use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{Context, Result};
use fgumi_lib::bam_io::{create_raw_bam_reader, create_raw_bam_writer};
use fgumi_lib::consensus_caller::ConsensusOutput;
use fgumi_lib::logging::{OperationTimer, log_consensus_summary};
use fgumi_lib::mi_group::RawMiGroupIterator;
use fgumi_lib::progress::ProgressTracker;
use fgumi_pipeline::stage::PipelineStage;
use fgumi_pipeline::stages::filter::FilterStage;
use fgumi_raw_bam::RawRecord;
use log::info;

use super::Runall;
use super::metrics::PipelineMetrics;
use crate::commands::consensus_runner::{ConsensusStatsOps, create_unmapped_consensus_header};

impl Runall {
    /// Run the pipeline starting from grouped input (MI tags already assigned).
    ///
    /// Reads records grouped by MI tag, calls consensus on each group, and writes
    /// filtered consensus reads to the output BAM. This is a serial implementation
    /// suitable for correctness validation; parallel execution will be added later.
    pub(super) fn run_from_group(
        &self,
        command_line: &str,
        tmp_output: &Path,
        cancel: &Arc<AtomicBool>,
    ) -> Result<()> {
        let timer = OperationTimer::new("Pipeline (from group)");
        let input = self.single_input()?;

        info!("Starting pipeline (from group, serial mode)");
        info!("Input:  {}", input.display());
        info!("Output: {}", self.output.display());
        info!("MI tag: {}", self.consensus_opts.consensus_mi_tag);
        info!("Min reads: {}", self.filter_opts.filter_min_reads);
        info!("Min input base quality: {}", self.filter_opts.filter_min_base_quality);

        // Open input BAM to read the header
        let (header, read_name_prefix) = {
            let (_reader, header) = create_raw_bam_reader(input, 1)?;
            let prefix =
                self.consensus_opts.consensus_read_name_prefix.clone().unwrap_or_else(|| {
                    fgumi_lib::consensus_caller::make_prefix_from_header(&header)
                });
            (header, prefix)
        };

        // Build output header (unmapped consensus reads)
        let output_header = create_unmapped_consensus_header(
            &header,
            &self.consensus_opts.consensus_read_group_id,
            "Pipeline consensus",
            command_line,
        )?;

        // Create output writer
        let mut writer = create_raw_bam_writer(
            tmp_output,
            &output_header,
            self.threads,
            self.compression_level,
        )?;

        // Create consensus caller and filter stage
        let mut caller = self.build_consensus_caller(read_name_prefix)?;
        let filter_stage = FilterStage::new(
            self.build_filter_config(),
            self.filter_opts.filter_require_ss_agreement,
        );

        // Create raw record iterator over input BAM
        let (mut raw_reader, _) = create_raw_bam_reader(input, self.threads)?;
        let raw_record_iter = std::iter::from_fn(move || {
            let mut record = RawRecord::new();
            match raw_reader.read_record(&mut record) {
                Ok(0) => None,
                Ok(_) => Some(Ok(record.into_inner())),
                Err(e) => Some(Err(e.into())),
            }
        });

        // Group by MI tag and process
        let mi_group_iter =
            RawMiGroupIterator::new(raw_record_iter, &self.consensus_opts.consensus_mi_tag);
        let progress = ProgressTracker::new("Consensus reads written").with_interval(100_000);
        let mut total_consensus_reads: usize = 0;
        let mut total_mi_groups: usize = 0;
        let mut total_input_records: usize = 0;

        let consensus_start = std::time::Instant::now();
        for result in mi_group_iter {
            if cancel.load(Ordering::Relaxed) {
                anyhow::bail!("Pipeline cancelled by signal");
            }
            let (umi, records) = result.context("Failed to read MI group")?;
            total_mi_groups += 1;
            total_input_records += records.len();

            let output: ConsensusOutput = caller
                .consensus_reads(records)
                .with_context(|| format!("Failed to call consensus for MI: {umi}"))?;

            if output.count > 0 {
                let filtered = filter_stage
                    .process(output)
                    .with_context(|| format!("Failed to filter consensus for MI: {umi}"))?;
                if !filtered.is_empty() {
                    // Count records in filtered output (each prefixed with 4-byte block_size).
                    let mut count = 0usize;
                    let mut off = 0;
                    while off + 4 <= filtered.len() {
                        let bs =
                            u32::from_le_bytes(filtered[off..off + 4].try_into().expect("4 bytes"))
                                as usize;
                        off += 4 + bs;
                        count += 1;
                    }
                    total_consensus_reads += count;
                    writer
                        .write_raw_bytes(&filtered)
                        .context("Failed to write filtered consensus reads")?;
                }
            }

            progress.log_if_needed(1);
        }
        let consensus_elapsed = consensus_start.elapsed().as_secs_f64();

        progress.log_final();
        info!(
            "Consensus complete: {} input records → {} MI groups → {} consensus reads",
            total_input_records, total_mi_groups, total_consensus_reads
        );

        // Finish output
        writer.finish().context("Failed to finish output BAM")?;

        // Log summary
        let stats = caller.statistics();
        let consensus_metrics = stats.to_metrics();
        log_consensus_summary(&consensus_metrics);

        info!("Total consensus reads written: {total_consensus_reads}");
        info!(
            "Pipeline complete: {} input records → {} consensus reads",
            total_input_records, total_consensus_reads
        );
        timer.log_completion(total_consensus_reads as u64);

        // Write pipeline stage metrics if requested
        if let Some(ref metrics_path) = self.metrics {
            let mut pm = PipelineMetrics::new();
            pm.add(
                "consensus",
                consensus_elapsed,
                total_input_records as u64,
                total_consensus_reads as u64,
            );
            pm.write_to_file(metrics_path)?;
            info!("Pipeline metrics written to {}", metrics_path.display());
        }

        Ok(())
    }
}
