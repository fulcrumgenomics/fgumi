//! Pipeline entry point: run from mapped BAM (sort mode).

use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use anyhow::{Context, Result};
use fgumi_lib::bam_io::{create_raw_bam_reader, create_raw_bam_writer};
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::sort::sink::{BamWriterSink, SortedRecordSink};
use fgumi_lib::sort::{RawExternalSorter, SortOrder};
use fgumi_pipeline::pipeline::builder::{PipelineBuilder, Stage as PipelineStage};
use log::info;

use super::Runall;
use super::metrics::PipelineMetrics;
use super::options::StopAfter;

impl Runall {
    /// Run the pipeline starting from a mapped BAM (serial sort mode).
    ///
    /// Sorts the input by template-coordinate, detects position-group boundaries,
    /// groups records by queryname into templates, assigns UMIs, calls consensus on
    /// each molecule group, and writes filtered output.
    ///
    /// When `stop_after` permits the full group→consensus→filter chain, delegates to
    /// the parallel [`PipelineBuilder`] (which works with any thread count). Otherwise
    /// falls back to the serial [`ProcessingSink`] path for `--stop-after sort/group`.
    pub(super) fn run_from_sort(
        &self,
        command_line: &str,
        tmp_output: &Path,
        cancel: &Arc<AtomicBool>,
    ) -> Result<()> {
        let input = self.single_input()?;

        info!("Input:  {}", input.display());
        info!("Output: {}", self.output.display());
        info!("Strategy: {:?}", self.group_opts.group_strategy);
        info!("Max edits: {}", self.group_opts.group_max_edits);
        info!("UMI tag: {}", self.group_opts.group_umi_tag);
        info!("Min reads: {}", self.filter_opts.filter_min_reads);
        info!("Min input base quality: {}", self.filter_opts.filter_min_base_quality);

        // Read header from input BAM (needed for both serial and --stop-after sort paths)
        let (_reader, header) = create_raw_bam_reader(input, 1)?;

        // --stop-after sort: write the sorted BAM directly without grouping/consensus.
        // Use the input header (with reference sequences) rather than a consensus header.
        if self.stop_after == Some(StopAfter::Sort) {
            let timer = OperationTimer::new("Pipeline (from sort, stop after sort)");
            let output_header = header.clone();
            let writer = create_raw_bam_writer(
                tmp_output,
                &output_header,
                self.threads,
                self.compression_level,
            )?;
            let mut sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
                .memory_limit(self.sort_opts.sort_memory_limit)
                .threads(self.threads)
                .temp_compression(1);
            if let Some(ref temp_dir) = self.sort_opts.sort_temp_dir {
                sorter = sorter.temp_dir(temp_dir.clone());
            }
            let mut bam_sink = BamWriterSink::new(writer);
            let sort_stats = sorter.sort_template_to_sink(input, &header, &mut bam_sink)?;
            bam_sink.finish().context("Failed to finish sorted BAM writer")?;
            info!(
                "Sort complete (stopped after sort): {} records in, {} records out",
                sort_stats.total_records, sort_stats.output_records
            );
            timer.log_completion(sort_stats.output_records);
            return Ok(());
        }

        // Parallel path: use PipelineBuilder for sort → group → consensus → filter
        if self.should_use_parallel_pipeline() {
            let timer = OperationTimer::new("Pipeline (from sort, parallel)");
            info!("Starting pipeline (from sort, parallel mode, {} threads)", self.threads);
            let config = self.build_pipeline_config(input, tmp_output, PipelineStage::Sort)?;
            let pipeline_start = std::time::Instant::now();
            PipelineBuilder::run(&config).context("Parallel pipeline failed")?;
            let pipeline_elapsed = pipeline_start.elapsed().as_secs_f64();
            info!("Parallel pipeline complete in {pipeline_elapsed:.1}s");
            timer.log_completion(0);

            if let Some(ref metrics_path) = self.metrics {
                let mut pm = PipelineMetrics::new();
                pm.add("sort+group+consensus (parallel)", pipeline_elapsed, 0, 0);
                pm.write_to_file(metrics_path)?;
                info!("Pipeline metrics written to {}", metrics_path.display());
            }
            return Ok(());
        }

        // Serial fallback path
        info!("Starting pipeline (from sort, serial mode)");

        let input = input.clone();
        let header_clone = header.clone();
        let sort_memory_limit = self.sort_opts.sort_memory_limit;
        let sort_temp_dir = self.sort_opts.sort_temp_dir.clone();
        let threads = self.threads;
        self.run_serial_pipeline(
            &header,
            tmp_output,
            move |sink| {
                let mut sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
                    .memory_limit(sort_memory_limit)
                    .threads(threads)
                    .temp_compression(1);
                if let Some(ref temp_dir) = sort_temp_dir {
                    sorter = sorter.temp_dir(temp_dir.clone());
                }
                sorter.sort_template_to_sink(&input, &header_clone, sink)
            },
            cancel,
            "from sort",
            command_line,
        )
    }
}
