//! Pipeline entry points: run from extract, correct, fastq, and align stages.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use anyhow::{Context, Result, anyhow};
use fgumi_lib::bam_io::create_raw_bam_writer;
use fgumi_lib::extract::{
    self, ExtractParams, QualityEncoding, ReadGroupMetadata, build_unmapped_bam_header,
    extract_template_and_fastq, make_raw_records,
};
use fgumi_lib::fastq::{FastqSet, ReadSetIterator};
use fgumi_lib::logging::OperationTimer;
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::template::{Template, TemplateIterator};
use fgumi_pipeline::stages::aligner::AlignerProcess;
use fgumi_simd_fastq::SimdFastqReader;
use fgumi_umi::TagInfo;
use log::info;
use noodles::sam::Header;
use read_structure::ReadStructure;

use super::Runall;
use super::aligner_helpers::shell_escape;
use super::options::StopAfter;
use super::run_zipper::zipper_merge_to_channel;
use crate::commands::zipper::build_output_header;

impl Runall {
    /// Run `fgumi fastq | aligner` to a temp file, then copy the result to the output.
    ///
    /// Shared by `run_from_extract`, `run_from_correct`, and `run_from_fastq` when
    /// `--stop-after align` is specified.
    fn run_stop_after_align(
        &self,
        fgumi_exe: &std::path::Path,
        unmapped_bam: &std::path::Path,
        aligner_command: &str,
        reference: &std::path::Path,
        tmp_output: &std::path::Path,
        timer: &OperationTimer,
    ) -> Result<()> {
        let temp_dir = tempfile::tempdir().context("Failed to create temp directory")?;
        let mapped_sam = temp_dir.path().join("mapped.sam");
        self.run_fastq_align_to_file(
            fgumi_exe,
            unmapped_bam,
            &mapped_sam,
            aligner_command,
            reference,
        )?;
        std::fs::copy(&mapped_sam, tmp_output).with_context(|| {
            format!(
                "Failed to copy aligned SAM {} to {}",
                mapped_sam.display(),
                tmp_output.display()
            )
        })?;
        info!("Stopped after align stage");
        timer.log_completion(0);
        Ok(())
    }

    /// Spawn an aligner with the given shell command, capture its output to a temp file,
    /// then copy the result to the output.
    ///
    /// Used by `run_from_align` for `--stop-after align` where the input is a FASTQ
    /// file (not a BAM), so the aligner command uses shell stdin redirection instead of
    /// `fgumi fastq` piping.
    fn run_stop_after_align_with_command(
        &self,
        shell_cmd: &str,
        tmp_output: &std::path::Path,
        timer: &OperationTimer,
    ) -> Result<()> {
        let temp_dir = tempfile::tempdir().context("Failed to create temp directory")?;
        let mapped_sam = temp_dir.path().join("mapped.sam");
        let mut aligner = AlignerProcess::spawn(shell_cmd, 50)?;
        let mut aligner_stdout = aligner.take_stdout().expect("aligner stdout was piped");
        let mut mapped_sam_file =
            std::fs::File::create(&mapped_sam).context("Failed to create temp mapped SAM file")?;
        std::io::copy(&mut aligner_stdout, &mut mapped_sam_file)
            .context("Failed to write aligner output to temp SAM file")?;
        drop(aligner_stdout);
        aligner.wait().context("Aligner subprocess failed")?;
        info!("Alignment complete: {}", mapped_sam.display());
        std::fs::copy(&mapped_sam, tmp_output).with_context(|| {
            format!(
                "Failed to copy aligned SAM {} to {}",
                mapped_sam.display(),
                tmp_output.display()
            )
        })?;
        info!("Stopped after align stage");
        timer.log_completion(0);
        Ok(())
    }

    /// Build an unmapped BAM header for the extract stage.
    ///
    /// Creates a SAM header with sort order, read group (including all metadata
    /// fields), comments, and PG record, matching what `fgumi extract` would
    /// produce. Delegates the core header construction to
    /// [`build_unmapped_bam_header`] in the library.
    fn build_unmapped_header(&self, command_line: &str) -> Result<Header> {
        let sample =
            self.extract_opts.extract_sample.as_ref().expect("--sample validated in validate()");
        let library =
            self.extract_opts.extract_library.as_ref().expect("--library validated in validate()");

        let metadata = ReadGroupMetadata {
            platform: Some(self.extract_opts.extract_platform.clone()),
            platform_unit: self.extract_opts.extract_platform_unit.clone(),
            platform_model: self.extract_opts.extract_platform_model.clone(),
            barcode: self.extract_opts.extract_barcode.clone(),
            sequencing_center: self.extract_opts.extract_sequencing_center.clone(),
            predicted_insert_size: self.extract_opts.extract_predicted_insert_size,
            description: self.extract_opts.extract_description.clone(),
            run_date: self.extract_opts.extract_run_date.clone(),
        };

        let header = build_unmapped_bam_header(
            &self.consensus_opts.consensus_read_group_id,
            sample,
            library,
            &metadata,
        )?;

        // Re-build with comments and @PG record added.
        let mut builder = Header::builder().set_header(
            header.header().expect("build_unmapped_bam_header always sets header record").clone(),
        );
        for (id, rg) in header.read_groups() {
            builder = builder.add_read_group(id.clone(), rg.clone());
        }
        if let Some(ref comments) = self.extract_opts.extract_comment {
            for comment in comments {
                builder = builder.add_comment(comment.clone());
            }
        }
        builder = crate::commands::common::add_pg_to_builder(builder, command_line)?;

        Ok(builder.build())
    }

    /// Build [`ExtractParams`] from the pipeline's CLI options.
    fn build_extract_params(&self) -> ExtractParams {
        ExtractParams {
            read_group_id: self.consensus_opts.consensus_read_group_id.clone(),
            umi_tag: self.group_opts.group_umi_tag.as_bytes().try_into().unwrap_or(*b"RX"),
            cell_tag: self.extract_opts.extract_cell_tag.as_bytes().try_into().unwrap_or(*b"CB"),
            umi_qual_tag: self
                .extract_opts
                .extract_umi_qual_tag
                .as_ref()
                .map(|t| t.as_bytes().try_into().unwrap_or(*b"QX")),
            cell_qual_tag: self
                .extract_opts
                .extract_cell_qual_tag
                .as_ref()
                .map(|t| t.as_bytes().try_into().unwrap_or(*b"CY")),
            single_tag: self
                .extract_opts
                .extract_single_tag
                .as_ref()
                .map(|t| t.as_bytes().try_into().unwrap_or(*b"BX")),
            annotate_read_names: self.extract_opts.extract_annotate_read_names,
            extract_umis_from_read_names: self.extract_opts.extract_extract_umis_from_read_names,
            store_sample_barcode_qualities: self
                .extract_opts
                .extract_store_sample_barcode_qualities,
        }
    }

    /// Parse read structures from the pipeline's extract options, defaulting to `+T`
    /// for each input FASTQ if none are specified.
    fn get_read_structures(&self) -> Result<Vec<ReadStructure>> {
        let rs_strings = self.extract_opts.extract_read_structures.as_deref().unwrap_or_default();
        if rs_strings.is_empty() {
            // Default: one `+T` per input FASTQ (template-only)
            Ok(vec![ReadStructure::from_str("+T")?; self.input.len()])
        } else {
            rs_strings
                .iter()
                .map(|s| ReadStructure::from_str(s).map_err(anyhow::Error::from))
                .collect()
        }
    }

    /// Detect quality encoding by sampling the first N records from the first FASTQ.
    fn detect_quality_encoding(&self) -> Result<QualityEncoding> {
        let first_input =
            self.input.first().ok_or_else(|| anyhow!("No input FASTQ files provided"))?;

        let reader_box: Box<dyn std::io::BufRead + Send> = Box::new(std::io::BufReader::new(
            fgoxide::io::Io::new(5, 1024 * 1024).new_reader(first_input)?,
        ));
        let mut fq_reader = SimdFastqReader::with_capacity(reader_box, 1024 * 1024);

        let mut sample_quals = Vec::new();
        for _i in 0..extract::QUALITY_DETECTION_SAMPLE_SIZE {
            match fq_reader.next() {
                Some(Ok(rec)) => sample_quals.push(rec.quality),
                Some(Err(e)) => return Err(e.into()),
                None => break,
            }
        }

        QualityEncoding::detect(&sample_quals)
    }

    /// Run the in-process extract pipeline: read FASTQs in-process, extract UMIs, pipe FASTQ
    /// text to the aligner, and stream unmapped templates to the zipper merge.
    ///
    /// This eliminates the temp unmapped BAM and the `fgumi extract` + `fgumi fastq`
    /// subprocesses. The aligner starts receiving data immediately (streaming).
    ///
    /// When [`should_use_parallel_pipeline`](Self::should_use_parallel_pipeline) returns
    /// `true`, uses the two-phase pipeline: Phase A (work-stealing pool for Extract,
    /// `ToFastq`, Zipper) → Phase B (Zone 3 from sorted chunks). Otherwise, falls back
    /// to the serial dedicated-thread approach.
    fn run_inline_extract(
        &self,
        command_line: &str,
        tmp_output: &Path,
        cancel: &Arc<AtomicBool>,
        timer: &OperationTimer,
    ) -> Result<()> {
        let aligner_command = self.require_aligner_command()?;
        let (reference, dict_path) = self.require_reference()?;

        // Build extract configuration
        let read_structures = self.get_read_structures()?;
        let extract_params = self.build_extract_params();
        let encoding = self.detect_quality_encoding()?;
        info!("Detected quality encoding: {encoding:?}");

        // Build unmapped header (for zipper merge output header construction)
        let unmapped_header = self.build_unmapped_header(command_line)?;

        // Spawn aligner process (reads interleaved FASTQ from stdin)
        let aligner_cmd_str = aligner_command
            .replace("{ref}", &reference.display().to_string())
            .replace("{threads}", &self.aligner_opts.aligner_threads.to_string());
        let mut aligner = AlignerProcess::spawn(&aligner_cmd_str, 50)?;

        // ---- Two-phase parallel path ----
        if self.should_use_parallel_pipeline() {
            info!("Using two-phase parallel pipeline (Phase A + Phase B)");
            return self.run_inline_extract_two_phase(
                command_line,
                tmp_output,
                cancel,
                timer,
                aligner,
                &unmapped_header,
                &dict_path,
                &read_structures,
                extract_params,
                encoding,
            );
        }

        // ---- Serial fallback path (existing) ----
        let aligner_stdin = aligner.take_stdin().expect("aligner stdin was piped");
        let aligner_stdout = aligner.take_stdout().expect("aligner stdout was piped");

        // Channel for unmapped templates (producer -> zipper merge)
        let (unmapped_tx, unmapped_rx) = std::sync::mpsc::sync_channel::<Result<Template>>(1024);

        // Thread A: FASTQ reader -> extract UMIs -> send (unmapped template, FASTQ bytes)
        let input_paths = self.input.clone();
        let producer = std::thread::Builder::new()
            .name("extract-producer".into())
            .spawn(move || -> Result<u64> {
                extract_producer(
                    &input_paths,
                    &read_structures,
                    &extract_params,
                    encoding,
                    unmapped_tx,
                    aligner_stdin,
                )
            })
            .context("Failed to spawn extract producer thread")?;

        // Read aligner stdout as SAM
        let buf_reader = std::io::BufReader::with_capacity(256 * 1024, aligner_stdout);
        let mut sam_reader = noodles::sam::io::Reader::new(buf_reader);
        let mapped_header = sam_reader.read_header()?;
        let (mapped_tx, mapped_rx) = std::sync::mpsc::sync_channel::<Result<Template>>(256);
        let mapped_header_clone = mapped_header.clone();
        std::thread::Builder::new()
            .name("aligner-reader".into())
            .spawn(move || {
                let iter = TemplateIterator::new(
                    sam_reader
                        .record_bufs(&mapped_header_clone)
                        .map(|r| r.map_err(anyhow::Error::from)),
                );
                for template in iter {
                    if mapped_tx.send(template).is_err() {
                        break;
                    }
                }
            })
            .context("Failed to spawn aligner reader thread")?;

        let unmapped_iter = std::iter::from_fn(move || unmapped_rx.recv().ok());
        let mapped_iter: Box<dyn Iterator<Item = Result<Template>> + Send> =
            Box::new(std::iter::from_fn(move || mapped_rx.recv().ok()));

        // Build output header from unmapped + mapped + reference dict
        let output_header = build_output_header(&unmapped_header, &mapped_header, &dict_path)?;
        let output_header = crate::commands::common::add_pg_record(output_header, command_line)?;

        let tag_info = TagInfo::new(
            self.zipper_opts.zipper_tags_to_remove.clone().unwrap_or_default(),
            self.zipper_opts.zipper_tags_to_reverse.clone().unwrap_or_default(),
            self.zipper_opts.zipper_tags_to_revcomp.clone().unwrap_or_default(),
        );

        // Stream merged records through a bounded channel
        let (merge_tx, merge_rx) = crossbeam_channel::bounded::<Vec<u8>>(10_000);
        let merge_header = output_header.clone();
        let merge_tag_info = tag_info.clone();
        let merge_skip_pa = self.zipper_opts.zipper_skip_pa_tags;

        let zipper_producer = std::thread::Builder::new()
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

        // Sort the streamed records and run group -> consensus -> filter.
        self.sort_and_process_records(
            merge_rx.iter(),
            &output_header,
            command_line,
            tmp_output,
            cancel,
            timer,
            "from extract (inline)",
        )?;

        // Wait for zipper merge to finish
        let record_count =
            zipper_producer.join().map_err(|_| anyhow!("Zipper merge thread panicked"))??;
        info!("Zipper merge complete, {record_count} raw records streamed");

        // Wait for extract producer to finish
        let extract_count =
            producer.join().map_err(|_| anyhow!("Extract producer thread panicked"))??;
        info!("Extract producer complete, {extract_count} records extracted");

        // Wait for aligner process to finish
        aligner.wait().context("Aligner subprocess failed")?;

        Ok(())
    }

    /// Two-phase parallel path for inline extract: Phase A (work-stealing pool) → Phase B
    /// (Zone 3 from sorted chunks).
    ///
    /// Reads the SAM header from aligner stdout on the main thread (resolving the
    /// chicken-and-egg issue), builds the output header, then delegates to
    /// [`PipelineBuilder::run_two_phase`].
    #[expect(clippy::too_many_arguments, reason = "wiring method passes through all dependencies")]
    fn run_inline_extract_two_phase(
        &self,
        command_line: &str,
        tmp_output: &Path,
        cancel: &Arc<AtomicBool>,
        timer: &OperationTimer,
        mut aligner: AlignerProcess,
        unmapped_header: &Header,
        dict_path: &Path,
        read_structures: &[ReadStructure],
        extract_params: ExtractParams,
        encoding: QualityEncoding,
    ) -> Result<()> {
        use fgumi_pipeline::pipeline::builder::{PipelineBuilder, TwoPhaseConfig};
        use fgumi_pipeline::pipeline::phase_a_graph::MergeFn;

        let aligner_stdin = aligner.take_stdin().expect("aligner stdin was piped");
        let aligner_stdout = aligner.take_stdout().expect("aligner stdout was piped");

        let tag_info = Arc::new(TagInfo::new(
            self.zipper_opts.zipper_tags_to_remove.clone().unwrap_or_default(),
            self.zipper_opts.zipper_tags_to_reverse.clone().unwrap_or_default(),
            self.zipper_opts.zipper_tags_to_revcomp.clone().unwrap_or_default(),
        ));

        // Build the FASTQ source iterator from input files.
        let input_paths = self.input.clone();
        let read_structs = read_structures.to_vec();
        let extract_source = build_extract_source(&input_paths, &read_structs)?;

        // Build the merge function (zipper::merge from the binary crate).
        let merge_fn: MergeFn = crate::commands::zipper::merge;

        // Build a closure that constructs the output header from the mapped header.
        // This runs on the stdout-reader thread after it reads the SAM header from
        // the aligner, breaking the deadlock where the main thread waited for the
        // header before starting the pool that feeds the aligner.
        let unmapped_hdr = unmapped_header.clone();
        let dict = dict_path.to_path_buf();
        let cmd = command_line.to_string();
        let build_header_fn = Box::new(
            move |mapped_header: &noodles::sam::Header| -> anyhow::Result<noodles::sam::Header> {
                let output = build_output_header(&unmapped_hdr, mapped_header, &dict)?;
                crate::commands::common::add_pg_record(output, &cmd)
            },
        );

        // Build the Phase B pipeline config with a placeholder header — the real header
        // will be available before Phase B starts (after Phase A completes).
        // `build_pipeline_config_from_header` needs the header for read-name prefix
        // generation, so we build a temporary one from the unmapped header + dict only.
        // However, Phase B gets the real header from the channel, so we need to build
        // the pipeline config lazily too. For now, we use the unmapped header as a
        // placeholder — the consensus read-name prefix is derived from RG entries which
        // are in the unmapped header.
        let pipeline_config = self.build_pipeline_config_from_header(
            unmapped_header,
            tmp_output, // input path unused by run_from_sorted_chunks
            tmp_output,
            fgumi_pipeline::pipeline::builder::Stage::Sort,
        )?;

        let sort_threads = 2.max(self.threads.saturating_sub(2));

        // Build correction config if UMI sequences are provided.
        let correction_config = {
            use fgumi_pipeline::pipeline::phase_a_graph::CorrectionConfig;
            let umis_inline = self.correct_opts.correct_umis.as_deref().unwrap_or_default();
            let umi_files = self.correct_opts.correct_umi_files.as_deref().unwrap_or_default();
            if umis_inline.is_empty() && umi_files.is_empty() {
                None
            } else {
                let (umi_sequences, umi_length) =
                    fgumi_lib::correct::load_umi_sequences(umis_inline, umi_files)?;
                let encoded_umi_set =
                    Arc::new(fgumi_lib::correct::EncodedUmiSet::new(&umi_sequences));
                let umi_tag_bytes: [u8; 2] = self
                    .group_opts
                    .group_umi_tag
                    .as_bytes()
                    .try_into()
                    .map_err(|_| anyhow!("UMI tag must be exactly 2 bytes"))?;
                info!(
                    "UMI correction enabled: {} UMIs, length {}, max_mismatches {}, cache_size {}",
                    umi_sequences.len(),
                    umi_length,
                    self.correct_opts.correct_max_mismatches,
                    self.correct_opts.correct_cache_size,
                );
                Some(CorrectionConfig {
                    encoded_umi_set,
                    umi_length,
                    max_mismatches: self.correct_opts.correct_max_mismatches,
                    min_distance_diff: 2, // matches the standalone `correct` command default
                    umi_tag: umi_tag_bytes,
                    cache_size: self.correct_opts.correct_cache_size,
                })
            }
        };

        // Save a reference to the metrics collector before pipeline_config is consumed.
        let metrics_collector = pipeline_config.metrics_collector.clone();

        let two_phase_config = TwoPhaseConfig {
            extract_source,
            extract_params,
            quality_encoding: encoding,
            aligner_process: aligner,
            aligner_stdin,
            aligner_stdout,
            build_output_header_fn: build_header_fn,
            tag_info,
            skip_pa_tags: self.zipper_opts.zipper_skip_pa_tags,
            merge_fn,
            batch_state: None,
            sort_memory_limit: self.sort_opts.sort_memory_limit,
            sort_temp_dir: self.sort_opts.sort_temp_dir.clone(),
            sort_threads,
            pipeline_config,
            correction_config,
            threads: self.threads,
            cancel: Arc::clone(cancel),
        };

        let pipeline_start = std::time::Instant::now();
        PipelineBuilder::run_two_phase(two_phase_config)
            .context("Two-phase parallel pipeline failed")?;
        let pipeline_elapsed = pipeline_start.elapsed().as_secs_f64();
        info!("Two-phase parallel pipeline complete in {pipeline_elapsed:.1}s");
        timer.log_completion(0);

        self.write_parallel_metrics(metrics_collector.as_ref())?;

        Ok(())
    }

    pub(super) fn run_from_extract(
        &self,
        command_line: &str,
        tmp_output: &Path,
        cancel: &Arc<AtomicBool>,
    ) -> Result<()> {
        let timer = OperationTimer::new("Pipeline (from extract)");

        let sample =
            self.extract_opts.extract_sample.as_ref().expect("--sample validated in validate()");
        let library =
            self.extract_opts.extract_library.as_ref().expect("--library validated in validate()");
        info!("Starting pipeline (from extract)");
        info!("FASTQ inputs: {:?}", self.input);
        info!("Output: {}", self.output.display());
        info!("Sample: {sample}, Library: {library}");

        // --stop-after extract: in-process extract to output BAM (no aligner)
        if self.stop_after == StopAfter::Extract {
            return self.run_stop_after_extract(command_line, tmp_output, &timer);
        }

        // --stop-after fastq: extract to unmapped BAM, then convert to interleaved FASTQ
        if self.stop_after == StopAfter::Fastq {
            return self.run_stop_after_fastq_via_extract(command_line, tmp_output, &timer);
        }

        let aligner_command = self.require_aligner_command()?;
        info!("Aligner: {aligner_command}");

        // --stop-after align: in-process extract -> aligner -> capture SAM output
        if self.stop_after == StopAfter::Align {
            return self.run_stop_after_extract_align(command_line, tmp_output, &timer);
        }

        // Full pipeline: use in-process extract (no temp BAM, streaming)
        info!("Using in-process extract (streaming to aligner)");
        self.run_inline_extract(command_line, tmp_output, cancel, &timer)
    }

    /// In-process extract for `--stop-after extract`: read FASTQs, build unmapped BAM
    /// records, and write directly to the output file. No aligner or zipper involved.
    fn run_stop_after_extract(
        &self,
        command_line: &str,
        tmp_output: &Path,
        timer: &OperationTimer,
    ) -> Result<()> {
        info!("Using in-process extract (--stop-after extract)");

        // Build extract configuration
        let read_structures = self.get_read_structures()?;
        let extract_params = self.build_extract_params();
        let encoding = self.detect_quality_encoding()?;
        info!("Detected quality encoding: {encoding:?}");

        // Build unmapped BAM header and writer
        let header = self.build_unmapped_header(command_line)?;
        let mut writer = create_raw_bam_writer(tmp_output, &header, 1, 6)?;

        let mut fq_iterators = create_fastq_iterators(&self.input, &read_structures)?;

        let progress =
            ProgressTracker::new("Extracted records (stop-after extract)").with_interval(1_000_000);
        let mut record_count: u64 = 0;
        let mut next_read_sets = Vec::with_capacity(fq_iterators.len());

        loop {
            next_read_sets.clear();
            for iter in fq_iterators.iter_mut() {
                if let Some(rec) = iter.next() {
                    next_read_sets.push(rec);
                } else {
                    break;
                }
            }

            if next_read_sets.is_empty() {
                break;
            }

            anyhow::ensure!(
                next_read_sets.len() == fq_iterators.len(),
                "FASTQ sources out of sync: got {} records but expected {}",
                next_read_sets.len(),
                fq_iterators.len()
            );

            // Combine read sets from all input files
            let combined = FastqSet::combine_readsets(std::mem::take(&mut next_read_sets));

            // Build raw BAM records (with block_size prefix) and write directly
            let extracted = make_raw_records(&combined, encoding, &extract_params)?;
            writer
                .write_raw_bytes(&extracted.data)
                .context("Failed to write extracted records to BAM")?;

            record_count += extracted.num_records;
            progress.log_if_needed(extracted.num_records);
        }

        writer.finish().context("Failed to finish output BAM writer")?;
        progress.log_final();
        info!("In-process extract complete: {record_count} records written");
        info!("Stopped after extract stage");
        timer.log_completion(record_count);
        Ok(())
    }

    /// In-process extract for `--stop-after fastq`: extract to a temp unmapped BAM, then
    /// convert to interleaved FASTQ. Produces the same FASTQ that would be piped to the
    /// aligner in the full pipeline.
    fn run_stop_after_fastq_via_extract(
        &self,
        command_line: &str,
        tmp_output: &Path,
        timer: &OperationTimer,
    ) -> Result<()> {
        info!("Using in-process extract + fastq (--stop-after fastq)");

        // Extract to a temp unmapped BAM first
        let temp_dir = tempfile::tempdir().context("Failed to create temp directory")?;
        let temp_bam = temp_dir.path().join("extracted.bam");

        let read_structures = self.get_read_structures()?;
        let extract_params = self.build_extract_params();
        let encoding = self.detect_quality_encoding()?;
        info!("Detected quality encoding: {encoding:?}");

        let header = self.build_unmapped_header(command_line)?;
        let mut writer = create_raw_bam_writer(&temp_bam, &header, 1, 6)?;
        let mut fq_iterators = create_fastq_iterators(&self.input, &read_structures)?;

        let progress =
            ProgressTracker::new("Extracted records (stop-after fastq)").with_interval(1_000_000);
        let mut record_count: u64 = 0;
        let mut next_read_sets = Vec::with_capacity(fq_iterators.len());

        loop {
            next_read_sets.clear();
            for iter in fq_iterators.iter_mut() {
                if let Some(rec) = iter.next() {
                    next_read_sets.push(rec);
                } else {
                    break;
                }
            }

            if next_read_sets.is_empty() {
                break;
            }

            anyhow::ensure!(
                next_read_sets.len() == fq_iterators.len(),
                "FASTQ sources out of sync: got {} records but expected {}",
                next_read_sets.len(),
                fq_iterators.len()
            );

            let combined = FastqSet::combine_readsets(std::mem::take(&mut next_read_sets));
            let extracted = make_raw_records(&combined, encoding, &extract_params)?;
            writer
                .write_raw_bytes(&extracted.data)
                .context("Failed to write extracted records to temp BAM")?;

            record_count += extracted.num_records;
            progress.log_if_needed(extracted.num_records);
        }

        writer.finish().context("Failed to finish temp BAM writer")?;
        progress.log_final();
        info!("Extract complete: {record_count} records written to temp BAM");

        // Convert the temp BAM to interleaved FASTQ
        let written = write_bam_to_fastq_file(&temp_bam, tmp_output)?;
        info!("Stopped after fastq stage ({written} FASTQ records written)");
        timer.log_completion(written);
        Ok(())
    }

    /// In-process extract for `--stop-after align`: read FASTQs, build unmapped records,
    /// pipe FASTQ to aligner, and capture aligner SAM output to the output file.
    fn run_stop_after_extract_align(
        &self,
        _command_line: &str,
        tmp_output: &Path,
        timer: &OperationTimer,
    ) -> Result<()> {
        info!("Using in-process extract + align (--stop-after align)");
        let aligner_command = self.require_aligner_command()?;
        let (reference, _dict_path) = self.require_reference()?;

        // Build extract configuration
        let read_structures = self.get_read_structures()?;
        let extract_params = self.build_extract_params();
        let encoding = self.detect_quality_encoding()?;
        info!("Detected quality encoding: {encoding:?}");

        // Spawn aligner process
        let aligner_cmd_str = aligner_command
            .replace("{ref}", &reference.display().to_string())
            .replace("{threads}", &self.aligner_opts.aligner_threads.to_string());
        let mut aligner = AlignerProcess::spawn(&aligner_cmd_str, 50)?;
        let aligner_stdin = aligner.take_stdin().expect("aligner stdin was piped");
        let mut aligner_stdout = aligner.take_stdout().expect("aligner stdout was piped");

        // Producer thread: read FASTQs, write interleaved FASTQ to aligner stdin
        let input_paths = self.input.clone();
        let producer = std::thread::Builder::new()
            .name("extract-producer-align".into())
            .spawn(move || -> Result<u64> {
                extract_fastq_only_producer(
                    &input_paths,
                    &read_structures,
                    &extract_params,
                    encoding,
                    aligner_stdin,
                )
            })
            .context("Failed to spawn extract producer thread")?;

        // Capture aligner output to the output file
        let mut output_file =
            std::fs::File::create(tmp_output).context("Failed to create output SAM file")?;
        std::io::copy(&mut aligner_stdout, &mut output_file)
            .context("Failed to write aligner output to output file")?;
        drop(aligner_stdout);

        // Wait for producer and aligner
        let extract_count =
            producer.join().map_err(|_| anyhow!("Extract producer thread panicked"))??;
        info!("Extract producer complete: {extract_count} templates processed");

        aligner.wait().context("Aligner subprocess failed")?;

        info!("Stopped after align stage");
        timer.log_completion(0);
        Ok(())
    }

    /// Run the pipeline starting from UMI correction on an unmapped BAM.
    ///
    /// Stages: correct UMIs (in-process) -> fastq -> aligner -> zipper merge -> sort ->
    /// group -> consensus -> filter -> write BAM.
    ///
    /// UMI correction is performed in-process using functions from
    /// [`fgumi_lib::correct`], eliminating the `fgumi correct` subprocess.
    pub(super) fn run_from_correct(
        &self,
        command_line: &str,
        tmp_output: &Path,
        cancel: &Arc<AtomicBool>,
    ) -> Result<()> {
        let timer = OperationTimer::new("Pipeline (from correct)");
        let input = self.single_input()?;

        info!("Starting pipeline (from correct)");
        info!("Input:  {}", input.display());
        info!("Output: {}", self.output.display());

        // Step 1: Correct UMIs in-process
        info!("Step 1/3: Correcting UMIs (in-process)...");
        let umis_inline = self.correct_opts.correct_umis.as_deref().unwrap_or_default();
        let umi_files = self.correct_opts.correct_umi_files.as_deref().unwrap_or_default();
        let (umi_sequences, umi_length) =
            fgumi_lib::correct::load_umi_sequences(umis_inline, umi_files)?;
        let encoded_umi_set = fgumi_lib::correct::EncodedUmiSet::new(&umi_sequences);

        let umi_tag_bytes: [u8; 2] = self
            .group_opts
            .group_umi_tag
            .as_bytes()
            .try_into()
            .map_err(|_| anyhow!("UMI tag must be exactly 2 bytes"))?;
        let max_mismatches = self.correct_opts.correct_max_mismatches;
        // Pipeline uses default min_distance_diff of 2 (matching the correct command default)
        let min_distance_diff = 2;

        let temp_dir = tempfile::tempdir().context("Failed to create temp directory")?;
        let corrected_bam = temp_dir.path().join("corrected.bam");

        {
            let (mut reader, header) = fgumi_lib::bam_io::create_bam_reader(input, 1)?;
            let header_clone = header.clone();
            let mut writer = create_raw_bam_writer(&corrected_bam, &header, 1, 6)?;
            let progress = ProgressTracker::new("Corrected templates").with_interval(1_000_000);
            let mut template_count: u64 = 0;
            let mut cache = Some(lru::LruCache::<Vec<u8>, fgumi_lib::correct::UmiMatch>::new(
                std::num::NonZero::new(100_000).unwrap(),
            ));

            // Iterate templates from the unmapped BAM using TemplateIterator
            let template_iter = TemplateIterator::new(
                reader.record_bufs(&header_clone).map(|r| r.map_err(anyhow::Error::from)),
            );

            let umi_tag = noodles::sam::alignment::record::data::field::Tag::new(
                umi_tag_bytes[0],
                umi_tag_bytes[1],
            );
            let mut encode_buf: Vec<u8> = Vec::new();

            for template_result in template_iter {
                let template = template_result?;

                // Extract UMI from the first record
                let umi_opt = template.records.first().and_then(|rec| {
                    rec.data().get(&umi_tag).map(|v| {
                        if let noodles::sam::alignment::record_buf::data::field::Value::String(s) =
                            v
                        {
                            String::from_utf8_lossy(s).to_string()
                        } else {
                            panic!("UMI tag has non-string value");
                        }
                    })
                });

                if let Some(umi) = umi_opt {
                    let correction = fgumi_lib::correct::compute_template_correction(
                        &umi,
                        umi_length,
                        false, // revcomp: pipeline doesn't expose this option
                        max_mismatches,
                        min_distance_diff,
                        &encoded_umi_set,
                        &mut cache,
                    );

                    // Write corrected records (only if matched)
                    if correction.matched {
                        for record in &template.records {
                            let mut record = record.clone();
                            if correction.needs_correction {
                                if correction.has_mismatches {
                                    record.data_mut().insert(
                                        noodles::sam::alignment::record::data::field::Tag::ORIGINAL_UMI_BARCODE_SEQUENCE,
                                        noodles::sam::alignment::record_buf::data::field::Value::String(
                                            correction.original_umi.clone().into(),
                                        ),
                                    );
                                }
                                if let Some(ref corrected) = correction.corrected_umi {
                                    record.data_mut().insert(
                                        umi_tag,
                                        noodles::sam::alignment::record_buf::data::field::Value::String(
                                            corrected.clone().into(),
                                        ),
                                    );
                                }
                            }
                            encode_buf.clear();
                            fgumi_lib::vendored::bam_codec::encode_record_buf(
                                &mut encode_buf,
                                &header,
                                &record,
                            )
                            .context("Failed to encode corrected record")?;
                            writer
                                .write_raw_record(&encode_buf)
                                .context("Failed to write corrected record")?;
                        }
                    }
                    // Unmatched templates are dropped (pipeline doesn't write rejects)
                } else {
                    // No UMI tag -- write records through unchanged
                    for record in &template.records {
                        encode_buf.clear();
                        fgumi_lib::vendored::bam_codec::encode_record_buf(
                            &mut encode_buf,
                            &header,
                            record,
                        )
                        .context("Failed to encode record")?;
                        writer
                            .write_raw_record(&encode_buf)
                            .context("Failed to write uncorrected record")?;
                    }
                }

                template_count += 1;
                progress.log_if_needed(1);
            }

            writer.finish().context("Failed to finish corrected BAM")?;
            progress.log_final();
            info!("In-process correction complete: {template_count} templates processed");
        }

        // --stop-after correct: copy the corrected unmapped BAM to the output
        if self.stop_after == StopAfter::Correct {
            std::fs::copy(&corrected_bam, tmp_output).with_context(|| {
                format!(
                    "Failed to copy corrected BAM {} to {}",
                    corrected_bam.display(),
                    tmp_output.display()
                )
            })?;
            info!("Stopped after correct stage");
            timer.log_completion(0);
            return Ok(());
        }

        // --stop-after fastq: convert the corrected unmapped BAM to interleaved FASTQ
        if self.stop_after == StopAfter::Fastq {
            let written = write_bam_to_fastq_file(&corrected_bam, tmp_output)?;
            info!("Stopped after fastq stage ({written} FASTQ records written)");
            timer.log_completion(0);
            return Ok(());
        }

        let (reference, _dict_path) = self.require_reference()?;
        let aligner_command = self.require_aligner_command()?;

        // Step 2: FASTQ + aligner (in-process BAM-to-FASTQ conversion)
        info!("Step 2/3: Aligning (in-process BAM to FASTQ | aligner)...");

        if self.stop_after == StopAfter::Align {
            let fgumi_exe =
                std::env::current_exe().context("Failed to determine fgumi binary path")?;
            return self.run_stop_after_align(
                &fgumi_exe,
                &corrected_bam,
                &aligner_command,
                &reference,
                tmp_output,
                &timer,
            );
        }

        // Spawn aligner, feed it FASTQ from the corrected BAM in-process
        let aligner_cmd_str = aligner_command
            .replace("{ref}", &reference.display().to_string())
            .replace("{threads}", &self.aligner_opts.aligner_threads.to_string());
        let mut aligner = AlignerProcess::spawn(&aligner_cmd_str, 50)?;
        let aligner_stdin = aligner.take_stdin().expect("aligner stdin was piped");

        // Spawn a thread to read corrected BAM -> write FASTQ to aligner stdin
        let corrected_bam_clone = corrected_bam.clone();
        let bam_to_fastq = std::thread::Builder::new()
            .name("bam-to-fastq".into())
            .spawn(move || -> Result<u64> {
                bam_to_fastq_writer(&corrected_bam_clone, aligner_stdin)
            })
            .context("Failed to spawn BAM-to-FASTQ thread")?;

        // Step 3: Zipper merge + sort + consensus + filter
        info!("Step 3/3: Zipper merge -> sort -> consensus -> filter...");
        let aligner_stdout = aligner.take_stdout().expect("aligner stdout was piped");
        let buf_reader = std::io::BufReader::with_capacity(256 * 1024, aligner_stdout);
        self.run_zipper_sort_consensus(
            command_line,
            &corrected_bam,
            buf_reader,
            tmp_output,
            cancel,
            &timer,
        )?;

        // Wait for BAM-to-FASTQ thread
        let fastq_count =
            bam_to_fastq.join().map_err(|_| anyhow!("BAM-to-FASTQ thread panicked"))??;
        info!("BAM-to-FASTQ complete: {fastq_count} records written");

        aligner.wait().context("Aligner subprocess failed")?;

        Ok(())
    }

    /// Run the pipeline starting from FASTQ conversion of an unmapped BAM.
    ///
    /// Stages: fastq -> aligner -> zipper merge -> sort -> group -> consensus ->
    /// filter -> write BAM.
    pub(super) fn run_from_fastq(
        &self,
        command_line: &str,
        tmp_output: &Path,
        cancel: &Arc<AtomicBool>,
    ) -> Result<()> {
        let timer = OperationTimer::new("Pipeline (from fastq)");
        let input = self.single_input()?;

        info!("Starting pipeline (from fastq)");
        info!("Input:  {}", input.display());
        info!("Output: {}", self.output.display());

        // --stop-after fastq: convert the unmapped BAM to interleaved FASTQ
        if self.stop_after == StopAfter::Fastq {
            let written = write_bam_to_fastq_file(input, tmp_output)?;
            info!("Stopped after fastq stage ({written} FASTQ records written)");
            timer.log_completion(0);
            return Ok(());
        }

        let (reference, _dict_path) = self.require_reference()?;
        let aligner_command = self.require_aligner_command()?;
        let fgumi_exe = std::env::current_exe().context("Failed to determine fgumi binary path")?;

        // Step 1: FASTQ + aligner
        info!("Step 1/2: Aligning (fastq | aligner)...");

        if self.stop_after == StopAfter::Align {
            return self.run_stop_after_align(
                &fgumi_exe,
                input,
                &aligner_command,
                &reference,
                tmp_output,
                &timer,
            );
        }

        // Pipe aligner stdout directly to zipper merge (no temp SAM file).
        let aligner = self.spawn_fastq_aligner(&fgumi_exe, input, &aligner_command, &reference)?;

        // Step 2: Zipper merge + sort + consensus + filter
        info!("Step 2/2: Zipper merge -> sort -> consensus -> filter...");
        self.run_align_zipper_sort(command_line, input, aligner, tmp_output, cancel, &timer)?;

        Ok(())
    }

    /// Run the pipeline starting from alignment of a FASTQ file.
    ///
    /// Stages: aligner -> zipper merge -> sort -> group -> consensus ->
    /// filter -> write BAM. Requires `--unmapped` for the unmapped BAM used in zipper
    /// merge, and `--input` is the FASTQ file to align.
    pub(super) fn run_from_align(
        &self,
        command_line: &str,
        tmp_output: &Path,
        cancel: &Arc<AtomicBool>,
    ) -> Result<()> {
        let timer = OperationTimer::new("Pipeline (from align)");
        let input = self.single_input()?;
        let unmapped_path = self.unmapped.as_ref().expect("--unmapped validated in validate()");
        let (reference, _dict_path) = self.require_reference()?;
        let aligner_command = self.require_aligner_command()?;

        info!("Starting pipeline (from align)");
        info!("FASTQ input: {}", input.display());
        info!("Unmapped BAM: {}", unmapped_path.display());
        info!("Output: {}", self.output.display());

        // Step 1: Align the FASTQ
        info!("Step 1/2: Aligning FASTQ...");
        let aligner_cmd_str = aligner_command
            .replace("{ref}", &reference.display().to_string())
            .replace("{threads}", &self.aligner_opts.aligner_threads.to_string());

        // Use AlignerProcess (via /bin/bash -c) so quoted arguments are handled correctly.
        // Feed the FASTQ file via shell stdin redirection.
        let shell_cmd =
            format!("{aligner_cmd_str} < {}", shell_escape(input.display().to_string()));

        // --stop-after align: capture aligner output to a temp file for copying
        if self.stop_after == StopAfter::Align {
            return self.run_stop_after_align_with_command(&shell_cmd, tmp_output, &timer);
        }

        // Pipe aligner stdout directly to zipper merge (no temp SAM file).
        let aligner = AlignerProcess::spawn(&shell_cmd, 50)?;

        // Step 2: Zipper merge + sort + consensus + filter
        info!("Step 2/2: Zipper merge -> sort -> consensus -> filter...");
        self.run_align_zipper_sort(
            command_line,
            unmapped_path,
            aligner,
            tmp_output,
            cancel,
            &timer,
        )?;

        Ok(())
    }
}

/// Create FASTQ read-set iterators from input paths and read structures.
///
/// Opens each FASTQ file with buffered I/O and pairs it with the corresponding read
/// structure, returning a vector of [`ReadSetIterator`]s ready for sequential iteration.
/// This is shared by [`extract_producer`], [`extract_fastq_only_producer`], and
/// `run_stop_after_extract` to avoid duplicating the reader setup logic.
fn create_fastq_iterators(
    input_paths: &[PathBuf],
    read_structures: &[ReadStructure],
) -> Result<Vec<ReadSetIterator>> {
    let fq_readers: Vec<Box<dyn std::io::BufRead + Send>> = input_paths
        .iter()
        .map(|p| -> Result<Box<dyn std::io::BufRead + Send>> {
            Ok(Box::new(std::io::BufReader::new(
                fgoxide::io::Io::new(5, 1024 * 1024).new_reader(p)?,
            )))
        })
        .collect::<Result<Vec<_>>>()?;

    let fq_sources: Vec<SimdFastqReader<Box<dyn std::io::BufRead + Send>>> =
        fq_readers.into_iter().map(|fq| SimdFastqReader::with_capacity(fq, 1024 * 1024)).collect();

    let fq_iterators: Vec<ReadSetIterator> = fq_sources
        .into_iter()
        .zip(read_structures.iter())
        .map(|(source, rs)| ReadSetIterator::new(rs.clone(), source, Vec::new()))
        .collect();

    Ok(fq_iterators)
}

/// Extract producer thread: reads FASTQ files, applies read structures, builds
/// unmapped templates and FASTQ text, sends them to the appropriate consumers.
///
/// For each template read from the FASTQ files:
/// 1. Apply read structures to split segments (UMI, template, barcode, etc.)
/// 2. Build an unmapped `Template` with UMI/barcode tags (for zipper merge)
/// 3. Build interleaved FASTQ text (for aligner stdin)
/// 4. Send the unmapped template through `unmapped_tx`
/// 5. Write the FASTQ text to `aligner_stdin`
fn extract_producer(
    input_paths: &[PathBuf],
    read_structures: &[ReadStructure],
    extract_params: &ExtractParams,
    encoding: QualityEncoding,
    unmapped_tx: std::sync::mpsc::SyncSender<Result<Template>>,
    mut aligner_stdin: std::process::ChildStdin,
) -> Result<u64> {
    let progress = ProgressTracker::new("Extracted records (inline)").with_interval(1_000_000);

    let mut fq_iterators = create_fastq_iterators(input_paths, read_structures)?;

    let mut record_count: u64 = 0;
    let mut next_read_sets = Vec::with_capacity(fq_iterators.len());

    loop {
        next_read_sets.clear();
        for iter in fq_iterators.iter_mut() {
            if let Some(rec) = iter.next() {
                next_read_sets.push(rec);
            } else {
                break;
            }
        }

        if next_read_sets.is_empty() {
            break;
        }

        anyhow::ensure!(
            next_read_sets.len() == fq_iterators.len(),
            "FASTQ sources out of sync: got {} records but expected {}",
            next_read_sets.len(),
            fq_iterators.len()
        );

        // Combine read sets from all input files
        let combined = FastqSet::combine_readsets(std::mem::take(&mut next_read_sets));

        // Build unmapped template + FASTQ text in one pass
        let extracted = extract_template_and_fastq(&combined, encoding, extract_params)?;

        // Send unmapped template to zipper merge
        if unmapped_tx.send(Ok(extracted.template)).is_err() {
            // Consumer dropped -- pipeline was cancelled
            anyhow::bail!("Extract producer: unmapped consumer dropped the channel");
        }

        // Write FASTQ text to aligner stdin
        aligner_stdin
            .write_all(&extracted.fastq_bytes)
            .with_context(|| "Extract producer: failed to write FASTQ to aligner stdin")?;

        record_count += 1;
        progress.log_if_needed(1);
    }

    // Close aligner stdin to signal EOF
    drop(aligner_stdin);

    progress.log_final();
    info!("Extract producer finished: {record_count} templates processed");
    Ok(record_count)
}

/// Extract producer thread for `--stop-after align`: reads FASTQ files, applies read
/// structures, and writes only interleaved FASTQ text to the aligner stdin.
///
/// Unlike [`extract_producer`], this does not build unmapped `Template` objects since
/// there is no zipper merge stage. This is more efficient for the stop-after-align path.
fn extract_fastq_only_producer(
    input_paths: &[PathBuf],
    read_structures: &[ReadStructure],
    extract_params: &ExtractParams,
    encoding: QualityEncoding,
    mut aligner_stdin: std::process::ChildStdin,
) -> Result<u64> {
    let progress =
        ProgressTracker::new("Extracted records (stop-after align)").with_interval(1_000_000);

    let mut fq_iterators = create_fastq_iterators(input_paths, read_structures)?;

    let mut record_count: u64 = 0;
    let mut next_read_sets = Vec::with_capacity(fq_iterators.len());

    loop {
        next_read_sets.clear();
        for iter in fq_iterators.iter_mut() {
            if let Some(rec) = iter.next() {
                next_read_sets.push(rec);
            } else {
                break;
            }
        }

        if next_read_sets.is_empty() {
            break;
        }

        anyhow::ensure!(
            next_read_sets.len() == fq_iterators.len(),
            "FASTQ sources out of sync: got {} records but expected {}",
            next_read_sets.len(),
            fq_iterators.len()
        );

        // Combine read sets from all input files
        let combined = FastqSet::combine_readsets(std::mem::take(&mut next_read_sets));

        // Build unmapped template + FASTQ text; we only need the FASTQ bytes
        let extracted = extract_template_and_fastq(&combined, encoding, extract_params)?;

        // Write FASTQ text to aligner stdin
        aligner_stdin
            .write_all(&extracted.fastq_bytes)
            .with_context(|| "Extract producer: failed to write FASTQ to aligner stdin")?;

        record_count += 1;
        progress.log_if_needed(1);
    }

    // Close aligner stdin to signal EOF
    drop(aligner_stdin);

    progress.log_final();
    info!("Extract-align producer finished: {record_count} templates processed");
    Ok(record_count)
}

/// Read records from a BAM file and write them as interleaved FASTQ to the given writer.
///
/// This replaces the subprocess `fgumi fastq -i <bam>` by performing the BAM-to-FASTQ
/// conversion in-process. Excludes secondary and supplementary alignments (flag 0x900),
/// matching the default behavior of `fgumi fastq`.
fn bam_to_fastq_writer<W: Write>(bam_path: &Path, mut writer: W) -> Result<u64> {
    use crate::commands::fastq::write_fastq_record;

    let (mut reader, _header) = fgumi_lib::bam_io::create_bam_reader(bam_path, 1)?;
    let mut buf_writer = std::io::BufWriter::with_capacity(64 * 1024 * 1024, &mut writer);

    let mut total: u64 = 0;
    let mut written: u64 = 0;
    let mut seq_buf: Vec<u8> = Vec::with_capacity(512);
    let mut qual_buf: Vec<u8> = Vec::with_capacity(512);

    for result in reader.records() {
        let record = result.context("Failed to read BAM record for FASTQ conversion")?;
        total += 1;

        let flags = record.flags();
        // Skip secondary and supplementary alignments (matches fgumi fastq defaults)
        if (flags.bits() & 0x900) != 0 {
            continue;
        }

        write_fastq_record(&mut buf_writer, &record, flags, false, &mut seq_buf, &mut qual_buf)?;
        written += 1;
    }

    buf_writer.flush().context("Failed to flush FASTQ output")?;
    drop(buf_writer);
    drop(writer);

    info!("BAM-to-FASTQ: read {total} records, wrote {written} FASTQ records");
    Ok(written)
}

/// Convert a BAM file to interleaved FASTQ, writing the output to a file.
///
/// Used by `--stop-after fastq` to emit the same interleaved FASTQ that would
/// normally be piped to the aligner.
fn write_bam_to_fastq_file(bam_path: &Path, output: &Path) -> Result<u64> {
    let file = std::fs::File::create(output)
        .with_context(|| format!("Failed to create FASTQ output: {}", output.display()))?;
    bam_to_fastq_writer(bam_path, file)
}

/// Build an [`ExtractSource`] from input FASTQ paths and read structures.
///
/// Creates a multi-file FASTQ iterator that yields one [`FastqSet`] per template
/// (combining read sets from all input files). This is the Phase A equivalent of the
/// serial `extract_producer` function.
fn build_extract_source(
    input_paths: &[PathBuf],
    read_structures: &[ReadStructure],
) -> Result<fgumi_pipeline::pipeline::phase_a_graph::ExtractSource> {
    use fgumi_pipeline::pipeline::phase_a_graph::ExtractSource;

    let mut fq_iterators = create_fastq_iterators(input_paths, read_structures)?;

    // Build a boxed iterator that yields combined FastqSets from all input files.
    let iter: Box<dyn Iterator<Item = FastqSet> + Send> = Box::new(std::iter::from_fn(move || {
        let mut next_read_sets = Vec::with_capacity(fq_iterators.len());
        for fq_iter in fq_iterators.iter_mut() {
            if let Some(rec) = fq_iter.next() {
                next_read_sets.push(rec);
            } else {
                return None;
            }
        }
        if next_read_sets.is_empty() {
            return None;
        }
        // If we got fewer read sets than iterators, the files are out of sync.
        // This will be caught by the assertion in the non-parallel path, but
        // for the iterator interface we just stop cleanly.
        if next_read_sets.len() != fq_iterators.len() {
            return None;
        }
        Some(FastqSet::combine_readsets(next_read_sets))
    }));

    Ok(ExtractSource::new(iter))
}
