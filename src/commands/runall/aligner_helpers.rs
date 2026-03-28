//! Aligner subprocess helpers and shared sort/process logic.

use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use anyhow::{Context, Result, anyhow};
use fgumi_lib::bam_io::{create_bam_reader, create_raw_bam_writer};
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::sort::sink::{BamWriterSink, SortedRecordSink};
use fgumi_lib::sort::{RawExternalSorter, SortOrder};
use fgumi_lib::template::{Template, TemplateIterator};
use fgumi_pipeline::pipeline::builder::{PipelineBuilder, Stage as PipelineStage};
use fgumi_pipeline::stages::aligner::AlignerProcess;
use fgumi_umi::TagInfo;
use log::info;
use noodles::sam::Header;

use super::Runall;
use super::metrics::PipelineMetrics;
use super::options::StopAfter;
use super::run_zipper::zipper_merge_to_channel;
use crate::commands::zipper::build_output_header;

impl Runall {
    pub(super) fn run_fastq_align_to_file(
        &self,
        fgumi_exe: &std::path::Path,
        unmapped_bam: &std::path::Path,
        mapped_sam: &std::path::Path,
        aligner_command: &str,
        reference: &std::path::Path,
    ) -> Result<()> {
        let mut aligner =
            self.spawn_fastq_aligner(fgumi_exe, unmapped_bam, aligner_command, reference)?;
        let mut aligner_stdout = aligner.take_stdout().expect("aligner stdout was piped");
        let mut mapped_sam_file =
            std::fs::File::create(mapped_sam).context("Failed to create temp mapped SAM file")?;
        std::io::copy(&mut aligner_stdout, &mut mapped_sam_file)
            .context("Failed to write aligner output to temp SAM file")?;
        drop(aligner_stdout);
        aligner.wait().context("Aligner subprocess failed")?;
        info!("Alignment complete: {}", mapped_sam.display());
        Ok(())
    }

    /// Spawn `fgumi fastq | aligner` as a subprocess pipeline, returning the
    /// [`AlignerProcess`] handle.
    ///
    /// The caller is responsible for reading stdout and calling [`AlignerProcess::wait`].
    pub(super) fn spawn_fastq_aligner(
        &self,
        fgumi_exe: &std::path::Path,
        unmapped_bam: &std::path::Path,
        aligner_command: &str,
        reference: &std::path::Path,
    ) -> Result<AlignerProcess> {
        let aligner_cmd_str = aligner_command
            .replace("{ref}", &reference.display().to_string())
            .replace("{threads}", &self.aligner_opts.aligner_threads.to_string());

        let fastq_align_cmd = format!(
            "{} fastq -i {} | {}",
            shell_escape(fgumi_exe.display().to_string()),
            shell_escape(unmapped_bam.display().to_string()),
            aligner_cmd_str
        );

        AlignerProcess::spawn(&fastq_align_cmd, 50)
    }

    /// Spawn `fgumi fastq | aligner` and pipe stdout directly into zipper merge + sort.
    ///
    /// Eliminates the temp mapped SAM file by reading the aligner's stdout as a SAM
    /// stream. The aligner process is waited on after the pipeline completes.
    pub(super) fn run_align_zipper_sort(
        &self,
        command_line: &str,
        unmapped_bam: &std::path::Path,
        mut aligner: AlignerProcess,
        tmp_output: &Path,
        cancel: &Arc<AtomicBool>,
        timer: &OperationTimer,
    ) -> Result<()> {
        let aligner_stdout = aligner.take_stdout().expect("aligner stdout was piped");
        let buf_reader = BufReader::with_capacity(256 * 1024, aligner_stdout);

        self.run_zipper_sort_consensus(
            command_line,
            unmapped_bam,
            buf_reader,
            tmp_output,
            cancel,
            timer,
        )?;

        aligner.wait().context("Aligner subprocess failed")?;
        Ok(())
    }

    /// Run zipper merge + sort + consensus + filter from unmapped BAM and a mapped SAM
    /// `Read` source.
    ///
    /// The mapped source can be a file, an aligner's stdout pipe, or any other `Read`
    /// implementation. This allows piping aligner output directly without a temp file.
    ///
    /// Helper shared by `run_from_correct`, `run_from_fastq`, `run_from_align`,
    /// `run_from_extract`.
    pub(super) fn run_zipper_sort_consensus<R: std::io::BufRead + Send + 'static>(
        &self,
        command_line: &str,
        unmapped_bam: &std::path::Path,
        mapped_sam_reader: R,
        tmp_output: &Path,
        cancel: &Arc<AtomicBool>,
        timer: &OperationTimer,
    ) -> Result<()> {
        // reference and dict validated in validate(); unwrap is safe
        let reference = self.reference.as_ref().expect("reference validated in validate()");
        let dict_path = fgumi_lib::reference::find_dict_path(reference)
            .expect("reference dict validated in validate()");

        // Open unmapped BAM
        let (mut unmapped_reader, unmapped_header) = create_bam_reader(unmapped_bam, 1)?;
        let unmapped_header_clone = unmapped_header.clone();
        let (unmapped_tx, unmapped_rx) = std::sync::mpsc::sync_channel::<Result<Template>>(256);
        std::thread::spawn(move || {
            let iter = TemplateIterator::new(
                unmapped_reader
                    .record_bufs(&unmapped_header_clone)
                    .map(|r| r.map_err(anyhow::Error::from)),
            );
            for template in iter {
                if unmapped_tx.send(template).is_err() {
                    break;
                }
            }
        });
        let unmapped_iter = std::iter::from_fn(move || unmapped_rx.recv().ok());

        // Read mapped SAM from the provided source
        let mut sam_reader = noodles::sam::io::Reader::new(mapped_sam_reader);
        let mapped_header = sam_reader.read_header()?;
        let mh = mapped_header.clone();
        let (mapped_tx, mapped_rx) = std::sync::mpsc::sync_channel::<Result<Template>>(256);
        std::thread::spawn(move || {
            let iter = TemplateIterator::new(
                sam_reader.record_bufs(&mh).map(|r| r.map_err(anyhow::Error::from)),
            );
            for template in iter {
                if mapped_tx.send(template).is_err() {
                    break;
                }
            }
        });
        let mapped_iter: Box<dyn Iterator<Item = Result<Template>> + Send> =
            Box::new(std::iter::from_fn(move || mapped_rx.recv().ok()));

        // Build output header
        let output_header = build_output_header(&unmapped_header, &mapped_header, &dict_path)?;
        let output_header = crate::commands::common::add_pg_record(output_header, command_line)?;

        let tag_info = TagInfo::new(
            self.zipper_opts.zipper_tags_to_remove.clone().unwrap_or_default(),
            self.zipper_opts.zipper_tags_to_reverse.clone().unwrap_or_default(),
            self.zipper_opts.zipper_tags_to_revcomp.clone().unwrap_or_default(),
        );

        // Stream merged records through a bounded channel to avoid materializing all
        // records in memory.
        let (merge_tx, merge_rx) = crossbeam_channel::bounded::<Vec<u8>>(10_000);
        let merge_header = output_header.clone();
        let merge_tag_info = tag_info.clone();
        let merge_skip_pa = self.zipper_opts.zipper_skip_pa_tags;

        let producer = std::thread::Builder::new()
            .name("zipper-merge".into())
            .spawn(move || -> Result<u64> {
                zipper_merge_to_channel(
                    unmapped_iter,
                    mapped_iter,
                    &merge_header,
                    &merge_tag_info,
                    merge_skip_pa,
                    merge_tx,
                )
            })
            .context("Failed to spawn zipper merge thread")?;

        // Sort the streamed records and run group → consensus → filter.
        self.sort_and_process_records(
            merge_rx.iter(),
            &output_header,
            command_line,
            tmp_output,
            cancel,
            timer,
            "from zipper (shared)",
        )?;

        // Wait for the producer to finish and propagate any errors.
        let record_count =
            producer.join().map_err(|_| anyhow!("Zipper merge thread panicked"))??;
        info!("Zipper merge complete, {record_count} raw records streamed");

        Ok(())
    }

    /// Sort streamed records, then run group → consensus → filter via the parallel
    /// pipeline builder (if eligible) or the serial [`ProcessingSink`] fallback.
    ///
    /// Shared by `run_from_zipper`, `run_from_extract`, and `run_zipper_sort_consensus`.
    ///
    /// Accepts an iterator of raw BAM record bytes, allowing the caller to stream records
    /// from a channel without materializing them all in memory.
    ///
    /// For the parallel path, the records are written to a temporary unsorted BAM file,
    /// then passed to [`PipelineBuilder::run`] with `Stage::Sort`. For the serial path,
    /// records are sorted and processed inline via the external sorter and
    /// `ProcessingSink`.
    #[expect(clippy::too_many_arguments, reason = "shared helper consolidates duplicated code")]
    pub(super) fn sort_and_process_records(
        &self,
        merged_records: impl Iterator<Item = Vec<u8>>,
        output_header: &Header,
        command_line: &str,
        tmp_output: &Path,
        cancel: &Arc<AtomicBool>,
        timer: &OperationTimer,
        label: &str,
    ) -> Result<()> {
        // --stop-after sort: write the sorted BAM directly without grouping/consensus
        if self.stop_after == Some(StopAfter::Sort) {
            let writer = create_raw_bam_writer(
                tmp_output,
                output_header,
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
            let sort_stats = sorter.sort_from_iter(merged_records, output_header, &mut bam_sink)?;
            bam_sink.finish().context("Failed to finish sorted BAM writer")?;
            info!(
                "Sort complete (stopped after sort): {} records in, {} records out",
                sort_stats.total_records, sort_stats.output_records
            );
            timer.log_completion(sort_stats.output_records);
            return Ok(());
        }

        // Parallel path: sort-accumulate directly from iterator, then merge into Zone 3.
        if self.should_use_parallel_pipeline() {
            info!("Starting pipeline ({label}, parallel mode, {} threads)", self.threads);

            // 1. Sort-accumulate directly from the iterator (no temp BAM).
            let sort_threads = 2.max(self.threads.saturating_sub(2));
            let mut sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
                .memory_limit(self.sort_opts.sort_memory_limit)
                .threads(sort_threads)
                .temp_compression(1);
            if let Some(ref temp_dir) = self.sort_opts.sort_temp_dir {
                sorter = sorter.temp_dir(temp_dir.clone());
            }
            let sorted_chunks = sorter
                .sort_accumulate_from_iter(merged_records, output_header)
                .context("Sort-accumulate from iterator failed")?;
            info!(
                "Sort-accumulate complete: {} records, {} disk chunks",
                sorted_chunks.stats.total_records, sorted_chunks.stats.chunks_written
            );

            // 2. Build Zone 3 pipeline config and run from sorted chunks.
            let config = self.build_pipeline_config_from_header(
                output_header,
                tmp_output, // input path unused by run_from_sorted_chunks
                tmp_output,
                PipelineStage::Sort,
            )?;
            let pipeline_start = std::time::Instant::now();
            PipelineBuilder::run_from_sorted_chunks(&config, output_header, sorted_chunks)
                .context("Parallel pipeline failed")?;
            let pipeline_elapsed = pipeline_start.elapsed().as_secs_f64();
            info!("Parallel pipeline ({label}) complete in {pipeline_elapsed:.1}s");
            timer.log_completion(0);

            if let Some(ref metrics_path) = self.metrics {
                let mut pm = PipelineMetrics::new();
                pm.add(
                    &format!("sort+group+consensus ({label}, parallel)"),
                    pipeline_elapsed,
                    0,
                    0,
                );
                pm.write_to_file(metrics_path)?;
                info!("Pipeline metrics written to {}", metrics_path.display());
            }
            return Ok(());
        }

        // Serial fallback path
        info!("Starting pipeline ({label}, serial mode)");

        let output_header_clone = output_header.clone();
        let sort_memory_limit = self.sort_opts.sort_memory_limit;
        let sort_temp_dir = self.sort_opts.sort_temp_dir.clone();
        let threads = self.threads;
        self.run_serial_pipeline(
            output_header,
            tmp_output,
            move |sink| {
                let mut sorter = RawExternalSorter::new(SortOrder::TemplateCoordinate)
                    .memory_limit(sort_memory_limit)
                    .threads(threads)
                    .temp_compression(1);
                if let Some(ref temp_dir) = sort_temp_dir {
                    sorter = sorter.temp_dir(temp_dir.clone());
                }
                sorter.sort_from_iter(merged_records, &output_header_clone, sink)
            },
            cancel,
            label,
            command_line,
        )
    }
}

// ============================================================================
// Shell helpers
// ============================================================================

/// Wrap a path string in single quotes for safe use in a shell command.
///
/// Replaces any embedded single quotes with `'\''` (end-quote, literal-quote, re-open-quote).
/// This ensures paths with spaces or special characters are passed correctly to `/bin/bash -c`.
pub(super) fn shell_escape(s: String) -> String {
    format!("'{}'", s.replace('\'', "'\\''"))
}
