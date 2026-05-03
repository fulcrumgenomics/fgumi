//! `--start-from extract --stop-after extract` chain runner.
//!
//! `ExtractChainRunner` takes a fully-configured `Runall` invocation,
//! detects compression and quality encoding for the input FASTQs, builds
//! the unmapped BAM header, and runs
//! [`crate::unified_pipeline::run_fastq_pipeline`] with closures that
//! convert each `FastqTemplate` into raw BAM record bytes via
//! [`crate::commands::extract::make_raw_records_static`]. The
//! pipeline's payload type is
//! [`crate::commands::extract::RawExtractedRecords`], whose
//! [`crate::unified_pipeline::MemoryEstimate`] impl drives the
//! memory-bounded backpressure between the process and serialize
//! steps.
//!
//! The runner is gated to the `(Extract, Extract)` chain shape; any
//! other shape falls through to the
//! [`super::not_implemented::NotImplementedRunner`] fallback. When the
//! caller passes `--metrics`, the runner times its `run_fastq_pipeline`
//! call and writes a single `MetricsRow { stage: "extract", ... }` to
//! the supplied TSV path (extract is a 1:1 transform, so
//! `records_in == records_out`).

use std::fs::File;
use std::io::{BufRead, Write};
use std::path::PathBuf;
use std::str::FromStr;
use std::sync::Arc;

use anyhow::{Context, Result};
use noodles::sam::Header;
use read_structure::ReadStructure;

use crate::commands::common::QueueMemoryOptions;
use crate::commands::extract::{
    BUFFER_SIZE, CompressionFormat, ExtractConfig, ExtractHeaderMetadata,
    QUALITY_DETECTION_SAMPLE_SIZE, QualityEncoding, RawExtractedRecords, build_extract_header,
    detect_compression_format, make_raw_records_static, open_fastq_reader,
    strip_read_suffix_extract,
};
use crate::commands::runall::Runall;
use crate::commands::runall::dispatch::{ChainContext, ChainRunner, DispatchContext};
use crate::commands::runall::infra::metrics::{MetricsRow, MetricsWriter};
use crate::commands::runall::options::{StartFrom, StopAfter};
use crate::fastq::FastqSet;
use crate::grouper::FastqTemplate;
use crate::sam::SamTag;
use crate::unified_pipeline::{FastqPipelineConfig, run_fastq_pipeline};
use fgumi_simd_fastq::SimdFastqReader;

/// Stable runner name surfaced via `--explain` and debug logs.
pub(crate) const EXTRACT_CHAIN_RUNNER_NAME: &str = "ExtractChainRunner";

/// Chain runner for `--start-from extract --stop-after extract`.
///
/// Built once per `Runall::execute` call; owns the cloned subset of
/// `Runall` fields needed to construct the unmapped BAM header and the
/// pipeline closures. Cloning keeps the runner `'static` (so it can be
/// pushed into the trait-object registry) while avoiding any wider
/// `Runall` lifetime gymnastics.
///
/// `clippy::struct_excessive_bools` is `#[allow]`-ed because each bool
/// here mirrors a distinct CLI flag on standalone `fgumi extract`;
/// collapsing them into a state machine would obscure the 1:1 mapping
/// without buying us anything (we never branch on combinations).
#[allow(clippy::struct_excessive_bools)]
pub(crate) struct ExtractChainRunner {
    inputs: Vec<PathBuf>,
    threads: usize,
    compression_level: u32,
    read_structures: Vec<String>,
    sample: Option<String>,
    library: Option<String>,
    read_group_id: String,
    barcode: Option<String>,
    platform: String,
    platform_unit: Option<String>,
    platform_model: Option<String>,
    sequencing_center: Option<String>,
    predicted_insert_size: Option<u32>,
    description: Option<String>,
    run_date: Option<String>,
    comments: Vec<String>,
    single_tag: Option<String>,
    annotate_read_names: bool,
    extract_umis_from_read_names: bool,
    store_sample_barcode_qualities: bool,
    store_umi_quals: bool,
    store_cell_quals: bool,
    /// `--extract::clipping-attribute`. Standalone `fgumi extract` declares
    /// this flag but does not read it for FASTQ input (no existing alignment
    /// clipping to adjust). Carried here for CLI parity; we log a warning
    /// when it's set and otherwise ignore it.
    clipping_attribute: Option<String>,
    /// `--queue-memory` / `--queue-memory-per-thread` from the top-level
    /// `Runall`. Cloned (rather than per-field-decomposed) so the runner
    /// can call `QueueMemoryOptions::calculate_memory_limit` and
    /// `QueueMemoryOptions::log_memory_config` exactly the way the
    /// standalone `Extract` command does, keeping memory-budget behavior
    /// in lock-step.
    queue_memory: QueueMemoryOptions,
}

impl ExtractChainRunner {
    /// Build the runner from a parsed `Runall` invocation.
    ///
    /// Builds unconditionally: the registry consults `supports()` to
    /// decide whether to dispatch here, so a `Runall` invocation with
    /// `start_from != Extract` (and therefore unset
    /// `--extract::sample`/`--extract::library`) is allowed to flow
    /// through this constructor. The required-when-running fields are
    /// stored as `Option`s and re-checked inside `run_chain()` so a
    /// missing flag surfaces as a clean `anyhow::Error` rather than a
    /// panic.
    pub(crate) fn from_runall(runall: &Runall) -> Self {
        let opts = &runall.extract_opts;
        Self {
            inputs: runall.input.clone(),
            threads: runall.threads,
            compression_level: runall.compression_level,
            read_structures: opts.extract_read_structures.clone().unwrap_or_default(),
            sample: opts.extract_sample.clone(),
            library: opts.extract_library.clone(),
            read_group_id: runall.consensus_opts.consensus_read_group_id.clone(),
            barcode: opts.extract_barcode.clone(),
            platform: opts.extract_platform.clone(),
            platform_unit: opts.extract_platform_unit.clone(),
            platform_model: opts.extract_platform_model.clone(),
            sequencing_center: opts.extract_sequencing_center.clone(),
            predicted_insert_size: opts.extract_predicted_insert_size,
            description: opts.extract_description.clone(),
            run_date: opts.extract_run_date.clone(),
            comments: opts.extract_comment.clone().unwrap_or_default(),
            single_tag: opts.extract_single_tag.clone(),
            annotate_read_names: opts.extract_annotate_read_names,
            extract_umis_from_read_names: opts.extract_extract_umis_from_read_names,
            store_sample_barcode_qualities: opts.extract_store_sample_barcode_qualities,
            store_umi_quals: opts.extract_store_umi_quals,
            store_cell_quals: opts.extract_store_cell_quals,
            clipping_attribute: opts.extract_clipping_attribute.clone(),
            queue_memory: runall.queue_memory.clone(),
        }
    }

    /// Build the [`FastqPipelineConfig`] for `run_chain`, applying the
    /// queue-memory budget computed from `--queue-memory` /
    /// `--queue-memory-per-thread` (parity with standalone
    /// `fgumi extract`). Extracted for unit-testing the wiring without
    /// running the full pipeline.
    fn build_pipeline_config(&self, all_bgzf: bool) -> Result<FastqPipelineConfig> {
        let pipeline_threads = self.threads.max(1);
        let queue_memory_limit_bytes =
            self.queue_memory.calculate_memory_limit(pipeline_threads)?;
        self.queue_memory.log_memory_config(pipeline_threads, queue_memory_limit_bytes);
        Ok(FastqPipelineConfig::new(pipeline_threads, all_bgzf, self.compression_level)
            .with_queue_memory_limit(queue_memory_limit_bytes))
    }

    /// Resolve the read structures, defaulting to `+T` per input when
    /// `--extract::read-structures` was not supplied. Mirrors the
    /// standalone `Extract::get_read_structures`.
    fn resolve_read_structures(&self) -> Result<Vec<ReadStructure>> {
        if self.read_structures.is_empty() {
            (0..self.inputs.len())
                .map(|_| ReadStructure::from_str("+T").map_err(anyhow::Error::from))
                .collect()
        } else {
            self.read_structures
                .iter()
                .map(|s| ReadStructure::from_str(s).map_err(anyhow::Error::from))
                .collect()
        }
    }

    /// Detect the FASTQ quality encoding by sampling the leading records
    /// of the first input. Mirrors the standalone `Extract` flow.
    fn detect_quality_encoding(&self) -> Result<QualityEncoding> {
        let first =
            self.inputs.first().ok_or_else(|| anyhow::anyhow!("no FASTQ inputs supplied"))?;
        // For sampling, single-thread decompression + no async-reader is
        // sufficient — the main pipeline opens its own readers below.
        let reader = open_fastq_reader(first, 1, false)
            .with_context(|| format!("opening {} for quality detection", first.display()))?;
        let mut fq_reader = SimdFastqReader::with_capacity(reader, BUFFER_SIZE);
        let mut sample_quals: Vec<Vec<u8>> = Vec::with_capacity(QUALITY_DETECTION_SAMPLE_SIZE);
        for _ in 0..QUALITY_DETECTION_SAMPLE_SIZE {
            match fq_reader.next() {
                Some(Ok(rec)) => sample_quals.push(rec.quality),
                Some(Err(e)) => return Err(e.into()),
                None => break,
            }
        }
        QualityEncoding::detect(&sample_quals)
    }

    /// Build the per-template `process_fn` for `run_fastq_pipeline`.
    ///
    /// Captures `read_structures`, `encoding`, and the user's extract
    /// options as cheap-to-clone owned values; the closure produces a
    /// [`RawExtractedRecords`] (raw BAM bytes + record count) per
    /// template. Send + Sync + 'static so it can be installed in the
    /// pipeline's process slot.
    fn build_process_fn(
        read_structures: Vec<ReadStructure>,
        encoding: QualityEncoding,
        cfg: ExtractConfig,
    ) -> impl Fn(FastqTemplate) -> std::io::Result<RawExtractedRecords> + Send + Sync + 'static
    {
        let read_structures = Arc::new(read_structures);
        let cfg = Arc::new(cfg);
        move |template: FastqTemplate| -> std::io::Result<RawExtractedRecords> {
            // Validate paired read names match across streams, mirroring the
            // standalone extract command's synchronized-mode check.
            if template.records.len() >= 2 {
                let base_name = strip_read_suffix_extract(template.records[0].name());
                for (i, record) in template.records.iter().enumerate().skip(1) {
                    let other_base = strip_read_suffix_extract(record.name());
                    if base_name != other_base {
                        return Err(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!(
                                "FASTQ files out of sync: R1 has '{}', R{} has '{}'",
                                String::from_utf8_lossy(base_name),
                                i + 1,
                                String::from_utf8_lossy(other_base),
                            ),
                        ));
                    }
                }
            }

            // Convert each FastqRecord to a FastqSet using its read structure.
            let mut fastq_sets: Vec<FastqSet> = Vec::with_capacity(template.records.len());
            for (record, rs) in template.records.iter().zip(read_structures.iter()) {
                let fastq_set = FastqSet::from_record_with_structure(
                    record.name(),
                    record.sequence(),
                    record.quality(),
                    rs,
                    &[],
                )
                .map_err(std::io::Error::other)?;
                fastq_sets.push(fastq_set);
            }
            let combined = FastqSet::combine_readsets(fastq_sets);
            make_raw_records_static(&combined, encoding, &cfg)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))
        }
    }

    /// Build the identity-shaped `serialize_fn` for `run_fastq_pipeline`.
    ///
    /// `make_raw_records_static` already produces raw BAM record bytes,
    /// so this closure just appends them to the worker's scratch buffer
    /// and reports the record count for upstream accounting.
    fn build_serialize_fn()
    -> impl Fn(RawExtractedRecords, &Header, &mut Vec<u8>) -> std::io::Result<u64>
    + Send
    + Sync
    + 'static {
        move |batch: RawExtractedRecords, _header: &Header, scratch: &mut Vec<u8>| {
            scratch.extend_from_slice(&batch.data);
            Ok(batch.num_records)
        }
    }

    /// Run the chain end-to-end against the prepared
    /// [`ChainContext::tmp_output`] path.
    ///
    /// When `metrics_path` is `Some`, writes a single
    /// [`MetricsRow`] tagged `"extract"` with wall-clock duration and
    /// `records_in == records_out` (extract performs no filtering).
    fn run_chain(
        &self,
        command_line: &str,
        output: &std::path::Path,
        metrics_path: Option<&std::path::Path>,
    ) -> Result<()> {
        anyhow::ensure!(
            !self.inputs.is_empty(),
            "ExtractChainRunner: at least one --input required"
        );
        let sample = self.sample.as_deref().ok_or_else(|| {
            anyhow::anyhow!("--extract::sample is required for the extract chain")
        })?;
        let library = self.library.as_deref().ok_or_else(|| {
            anyhow::anyhow!("--extract::library is required for the extract chain")
        })?;

        // Detect compression format per input; require uniformity so
        // the pipeline configuration (BGZF vs decompressed-readers) is
        // unambiguous.
        let formats: Vec<CompressionFormat> = self
            .inputs
            .iter()
            .map(|p| {
                detect_compression_format(p)
                    .with_context(|| format!("detecting compression format for {}", p.display()))
            })
            .collect::<Result<_>>()?;
        let all_bgzf = formats.iter().all(|f| *f == CompressionFormat::Bgzf);

        // Detect quality encoding from the first FASTQ.
        let encoding = self.detect_quality_encoding()?;
        log::info!("ExtractChainRunner: detected quality encoding {encoding:?}");

        // Resolve read structures and validate count vs inputs.
        let read_structures = self.resolve_read_structures()?;
        anyhow::ensure!(
            read_structures.len() == self.inputs.len(),
            "ExtractChainRunner: --extract::read-structures must have one entry per --input \
             (got {} structures vs {} inputs)",
            read_structures.len(),
            self.inputs.len(),
        );

        // Build the unmapped BAM header (matches `Extract::create_header`).
        let metadata = ExtractHeaderMetadata {
            barcode: self.barcode.as_deref(),
            platform: &self.platform,
            platform_unit: self.platform_unit.as_deref(),
            platform_model: self.platform_model.as_deref(),
            sequencing_center: self.sequencing_center.as_deref(),
            predicted_insert_size: self.predicted_insert_size,
            description: self.description.as_deref(),
            run_date: self.run_date.as_deref(),
        };
        let header = build_extract_header(
            &self.read_group_id,
            sample,
            library,
            &self.comments,
            &metadata,
            command_line,
        )?;

        // `--extract::clipping-attribute` is a no-op for FASTQ input
        // (matches standalone `fgumi extract`, where the flag is parsed
        // but never consumed by the FASTQ→BAM path). Warn the user so
        // they don't think the value took effect.
        if self.clipping_attribute.is_some() {
            log::warn!("--extract::clipping-attribute has no effect on FASTQ input; ignoring");
        }

        // Build the closure config (mirrors `Extract`'s ExtractConfig
        // construction).
        let cfg = ExtractConfig {
            read_group_id: self.read_group_id.clone(),
            store_umi_quals: self.store_umi_quals,
            store_cell_quals: self.store_cell_quals,
            single_tag: self.single_tag.as_ref().map(|t| {
                SamTag::from_str(t).expect("--extract::single-tag length validated by validate.rs")
            }),
            annotate_read_names: self.annotate_read_names,
            extract_umis_from_read_names: self.extract_umis_from_read_names,
            store_sample_barcode_qualities: self.store_sample_barcode_qualities,
        };

        // Build pipeline config. For Gzip / Plain inputs we open the
        // decompressed readers up front; for BGZF inputs we pass `None`
        // and let `run_fastq_pipeline` open the paths itself (so its
        // multi-threaded BGZF reader can engage).
        let pipeline_threads = self.threads.max(1);
        let config = self.build_pipeline_config(all_bgzf)?;

        let decompressed_readers: Option<Vec<Box<dyn BufRead + Send>>> = if all_bgzf {
            None
        } else {
            let decomp_threads = pipeline_threads;
            Some(
                self.inputs
                    .iter()
                    .map(|p| open_fastq_reader(p, decomp_threads, false))
                    .collect::<Result<Vec<_>>>()?,
            )
        };

        let process_fn = Self::build_process_fn(read_structures, encoding, cfg);
        let serialize_fn = Self::build_serialize_fn();

        let writer: Box<dyn Write + Send> = Box::new(std::io::BufWriter::new(
            File::create(output)
                .with_context(|| format!("create output BAM {}", output.display()))?,
        ));

        let start = std::time::Instant::now();
        let records_out = run_fastq_pipeline(
            config,
            &self.inputs,
            decompressed_readers,
            &header,
            writer,
            process_fn,
            serialize_fn,
        )
        .with_context(|| format!("ExtractChainRunner pipeline writing to {}", output.display()))?;
        let wall_time_secs = start.elapsed().as_secs_f64();

        if let Some(path) = metrics_path {
            // The extract stage is a 1:1 transform — every input
            // template becomes one output BAM record per template
            // segment, no filtering. `run_fastq_pipeline` returns the
            // same record count we'd have computed for `records_in`.
            let mut writer = MetricsWriter::create(path)
                .with_context(|| format!("opening --metrics path {} for write", path.display()))?;
            writer
                .write_row(&MetricsRow {
                    stage: "extract".into(),
                    wall_time_secs,
                    records_in: records_out,
                    records_out,
                })
                .with_context(|| format!("writing extract metrics row to {}", path.display()))?;
        }

        Ok(())
    }
}

impl ChainRunner for ExtractChainRunner {
    fn name(&self) -> &'static str {
        EXTRACT_CHAIN_RUNNER_NAME
    }

    fn supports(&self, ctx: &DispatchContext<'_>) -> bool {
        ctx.start_from == StartFrom::Extract && ctx.stop_after == StopAfter::Extract
    }

    fn run(&self, ctx: ChainContext<'_>) -> Result<()> {
        let metrics_path = ctx.dispatch.metrics_path;
        self.run_chain(ctx.command_line, &ctx.tmp_output, metrics_path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::runall::infra::progress::ProgressMode;
    use crate::commands::runall::options::{
        AlignerOptions, CodecOptions, ConsensusMode, ConsensusOptions, CorrectOptions,
        DuplexOptions, ExtractOptions, FilterOptions, GroupOptions, MultiAlignerOptions,
        MultiCodecOptions, MultiConsensusOptions, MultiCorrectOptions, MultiDuplexOptions,
        MultiExtractOptions, MultiFilterOptions, MultiGroupOptions, MultiSortOptions,
        MultiZipperOptions, SortOptions, ZipperOptions,
    };
    use std::path::Path;

    /// Build a minimal `Runall` for `(Extract, Extract)` with
    /// `--extract::sample`/`--extract::library` set. Defaults the read
    /// structures to `["8M+T", "+T"]`; callers that need different
    /// shapes (e.g. cell-barcode segments for the CY-tag test) should
    /// build via [`runall_for_extract_chain_with_opts`].
    fn runall_for_extract_chain(inputs: Vec<PathBuf>, output: PathBuf) -> Runall {
        let extract = ExtractOptions {
            sample: Some("S".to_string()),
            library: Some("L".to_string()),
            read_structures: Some(vec!["8M+T".to_string(), "+T".to_string()]),
            ..ExtractOptions::default()
        };
        runall_for_extract_chain_with_opts(inputs, output, extract)
    }

    /// As [`runall_for_extract_chain`] but with caller-supplied
    /// [`ExtractOptions`] (used by the QX/CY tag tests to flip
    /// `store_umi_quals` / `store_cell_quals` and tune read structures).
    fn runall_for_extract_chain_with_opts(
        inputs: Vec<PathBuf>,
        output: PathBuf,
        extract: ExtractOptions,
    ) -> Runall {
        Runall {
            input: inputs,
            output,
            reference: None,
            start_from: StartFrom::Extract,
            stop_after: StopAfter::Extract,
            consensus_mode: ConsensusMode::Simplex,
            threads: 4,
            progress: ProgressMode::None,
            compression_level: 1,
            metrics: None,
            explain: false,
            consensus_metrics: None,
            intervals: None,
            methylation_mode: None,
            restore_unconverted_bases: false,
            min_methylation_depth: vec![],
            require_strand_methylation_agreement: false,
            min_conversion_fraction: None,
            sort_opts: MultiSortOptions::from(SortOptions::default()),
            group_opts: MultiGroupOptions::from(GroupOptions::default()),
            filter_opts: MultiFilterOptions::from(FilterOptions::default()),
            consensus_opts: MultiConsensusOptions::from(ConsensusOptions::default()),
            extract_opts: MultiExtractOptions::from(extract),
            aligner_opts: MultiAlignerOptions::from(AlignerOptions::default()),
            zipper_opts: MultiZipperOptions::from(ZipperOptions::default()),
            correct_opts: MultiCorrectOptions::from(CorrectOptions::default()),
            duplex_opts: MultiDuplexOptions::from(DuplexOptions::default()),
            codec_opts: MultiCodecOptions::from(CodecOptions::default()),
            queue_memory: QueueMemoryOptions::default(),
        }
    }

    fn dispatch_ctx(start: StartFrom, stop: StopAfter) -> DispatchContext<'static> {
        DispatchContext { start_from: start, stop_after: stop, metrics_path: None }
    }

    fn write_paired_gzip_fastq(r1: &Path, r2: &Path, num_pairs: usize) {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        const DNA: [u8; 4] = [b'A', b'C', b'G', b'T'];

        let f1 = File::create(r1).unwrap();
        let mut e1 = GzEncoder::new(f1, Compression::default());
        let f2 = File::create(r2).unwrap();
        let mut e2 = GzEncoder::new(f2, Compression::default());
        for i in 0..num_pairs {
            let umi: [u8; 8] = std::array::from_fn(|k| DNA[(i + k) % 4]);
            let umi_str = std::str::from_utf8(&umi).unwrap();
            let r1_seq = format!("{umi_str}AAAAAAAAAAAA");
            let r2_seq = "CCCCCCCCCCCC";
            writeln!(e1, "@read{i}/1").unwrap();
            writeln!(e1, "{r1_seq}").unwrap();
            writeln!(e1, "+").unwrap();
            writeln!(e1, "{}", "I".repeat(r1_seq.len())).unwrap();
            writeln!(e2, "@read{i}/2").unwrap();
            writeln!(e2, "{r2_seq}").unwrap();
            writeln!(e2, "+").unwrap();
            writeln!(e2, "{}", "I".repeat(r2_seq.len())).unwrap();
        }
        e1.finish().unwrap();
        e2.finish().unwrap();
    }

    #[test]
    fn supports_extract_to_extract_with_or_without_metrics() {
        let dir = tempfile::TempDir::new().unwrap();
        let inputs = vec![dir.path().join("r1.fq.gz"), dir.path().join("r2.fq.gz")];
        let output = dir.path().join("out.bam");
        // Inputs don't have to exist for `supports` — just for `run`.
        let runall = runall_for_extract_chain(inputs, output);
        let runner = ExtractChainRunner::from_runall(&runall);

        assert!(runner.supports(&dispatch_ctx(StartFrom::Extract, StopAfter::Extract)));
        assert!(!runner.supports(&dispatch_ctx(StartFrom::Extract, StopAfter::Correct)));
        assert!(!runner.supports(&dispatch_ctx(StartFrom::Correct, StopAfter::Correct)));
        assert!(!runner.supports(&dispatch_ctx(StartFrom::Group, StopAfter::Group)));

        let metrics_path = dir.path().join("metrics.tsv");
        let with_metrics = DispatchContext {
            start_from: StartFrom::Extract,
            stop_after: StopAfter::Extract,
            metrics_path: Some(metrics_path.as_path()),
        };
        assert!(
            runner.supports(&with_metrics),
            "extract chain runner now writes its own --metrics row; the registry \
             must dispatch to it even when --metrics is set"
        );
    }

    #[test]
    fn name_is_stable() {
        let dir = tempfile::TempDir::new().unwrap();
        let runall = runall_for_extract_chain(
            vec![dir.path().join("r1.fq.gz"), dir.path().join("r2.fq.gz")],
            dir.path().join("out.bam"),
        );
        let runner = ExtractChainRunner::from_runall(&runall);
        assert_eq!(runner.name(), EXTRACT_CHAIN_RUNNER_NAME);
    }

    #[test]
    fn run_chain_writes_a_valid_bam_for_paired_gzip_input() {
        let dir = tempfile::TempDir::new().unwrap();
        let r1 = dir.path().join("R1.fastq.gz");
        let r2 = dir.path().join("R2.fastq.gz");
        let out_tmp = dir.path().join("extract.bam.tmp");

        write_paired_gzip_fastq(&r1, &r2, 50);

        let runall = runall_for_extract_chain(vec![r1.clone(), r2.clone()], out_tmp.clone());
        let runner = ExtractChainRunner::from_runall(&runall);
        runner.run_chain("fgumi runall", &out_tmp, None).expect("chain runs");

        let mut reader =
            noodles::bam::io::reader::Builder.build_from_path(&out_tmp).expect("open BAM");
        let header = reader.read_header().expect("read header");
        let mut record = noodles::sam::alignment::RecordBuf::default();
        let mut count = 0usize;
        while reader.read_record_buf(&header, &mut record).expect("read") != 0 {
            count += 1;
        }
        // Two records per template (R1 + R2), 50 templates.
        assert_eq!(count, 100);
    }

    /// Read a string-valued aux tag from a BAM record buffer (mirrors
    /// `commands::extract::tests::get_tag_string`). Returns `None` if
    /// the tag is absent or not a string.
    fn get_tag_string(
        record: &noodles::sam::alignment::RecordBuf,
        tag_name: &str,
    ) -> Option<String> {
        use noodles::sam::alignment::record::data::field::Tag;
        use noodles::sam::alignment::record_buf::data::field::Value as RecordBufValue;
        let tag_bytes = tag_name.as_bytes();
        let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);
        record.data().get(&tag).and_then(|value| match value {
            RecordBufValue::String(s) => Some(String::from_utf8_lossy(s.as_ref()).to_string()),
            _ => None,
        })
    }

    /// Count records and report whether any record carries `tag_name`
    /// as a string aux tag.
    fn count_records_and_tag_hits(out: &std::path::Path, tag_name: &str) -> (usize, usize) {
        let mut reader = noodles::bam::io::reader::Builder.build_from_path(out).expect("open BAM");
        let header = reader.read_header().expect("read header");
        let mut record = noodles::sam::alignment::RecordBuf::default();
        let mut total = 0usize;
        let mut hits = 0usize;
        while reader.read_record_buf(&header, &mut record).expect("read") != 0 {
            total += 1;
            if get_tag_string(&record, tag_name).is_some() {
                hits += 1;
            }
        }
        (total, hits)
    }

    #[test]
    fn metrics_tsv_emits_one_row() {
        let dir = tempfile::TempDir::new().unwrap();
        let r1 = dir.path().join("R1.fastq.gz");
        let r2 = dir.path().join("R2.fastq.gz");
        let out_tmp = dir.path().join("extract.bam.tmp");
        let metrics = dir.path().join("metrics.tsv");

        write_paired_gzip_fastq(&r1, &r2, 50);

        let runall = runall_for_extract_chain(vec![r1.clone(), r2.clone()], out_tmp.clone());
        let runner = ExtractChainRunner::from_runall(&runall);
        runner.run_chain("fgumi runall", &out_tmp, Some(metrics.as_path())).expect("chain runs");

        let contents = std::fs::read_to_string(&metrics).expect("read metrics");
        let mut lines = contents.lines();
        assert_eq!(
            lines.next().unwrap(),
            "stage\twall_time_secs\trecords_in\trecords_out",
            "header row"
        );
        let row = lines.next().expect("one extract row");
        let cols: Vec<&str> = row.split('\t').collect();
        assert_eq!(cols.len(), 4, "row has 4 columns: {row}");
        assert_eq!(cols[0], "extract", "stage column");
        let wall: f64 = cols[1].parse().expect("wall_time_secs is a float");
        assert!(wall >= 0.0, "wall_time_secs is non-negative: {wall}");
        // 50 templates * 2 records (R1 + R2).
        assert_eq!(cols[2], "100", "records_in: 50 templates × 2 records");
        assert_eq!(cols[3], "100", "records_out: extract is 1:1, no filtering");
        assert!(lines.next().is_none(), "exactly one data row");
    }

    #[test]
    fn store_umi_quals_writes_qx_tag() {
        let dir = tempfile::TempDir::new().unwrap();
        let r1 = dir.path().join("R1.fastq.gz");
        let r2 = dir.path().join("R2.fastq.gz");
        let out_tmp = dir.path().join("extract.bam.tmp");

        write_paired_gzip_fastq(&r1, &r2, 5);

        let extract = ExtractOptions {
            sample: Some("S".to_string()),
            library: Some("L".to_string()),
            // R1 = 8M+T → 8 bp UMI then template; R2 = +T template.
            read_structures: Some(vec!["8M+T".to_string(), "+T".to_string()]),
            store_umi_quals: true,
            ..ExtractOptions::default()
        };
        let runall = runall_for_extract_chain_with_opts(
            vec![r1.clone(), r2.clone()],
            out_tmp.clone(),
            extract,
        );
        let runner = ExtractChainRunner::from_runall(&runall);
        runner.run_chain("fgumi runall", &out_tmp, None).expect("chain runs");

        let (total, qx_hits) = count_records_and_tag_hits(&out_tmp, "QX");
        assert!(total > 0, "produced records");
        assert!(
            qx_hits > 0,
            "at least one record must carry a QX tag when store_umi_quals=true \
             (got {qx_hits}/{total})"
        );
    }

    #[test]
    fn store_cell_quals_writes_cy_tag() {
        let dir = tempfile::TempDir::new().unwrap();
        let r1 = dir.path().join("R1.fastq.gz");
        let r2 = dir.path().join("R2.fastq.gz");
        let out_tmp = dir.path().join("extract.bam.tmp");

        // Reuse the standard fixture — R1 has 8 bp + 12 bp template.
        // We carve those leading 8 bp into a *cellular* barcode (`C`)
        // segment so the extractor populates the CY tag; R2 stays
        // template-only.
        write_paired_gzip_fastq(&r1, &r2, 5);

        let extract = ExtractOptions {
            sample: Some("S".to_string()),
            library: Some("L".to_string()),
            read_structures: Some(vec!["8C+T".to_string(), "+T".to_string()]),
            store_cell_quals: true,
            ..ExtractOptions::default()
        };
        let runall = runall_for_extract_chain_with_opts(
            vec![r1.clone(), r2.clone()],
            out_tmp.clone(),
            extract,
        );
        let runner = ExtractChainRunner::from_runall(&runall);
        runner.run_chain("fgumi runall", &out_tmp, None).expect("chain runs");

        let (total, cy_hits) = count_records_and_tag_hits(&out_tmp, "CY");
        assert!(total > 0, "produced records");
        assert!(
            cy_hits > 0,
            "at least one record must carry a CY tag when store_cell_quals=true \
             (got {cy_hits}/{total})"
        );
    }

    #[test]
    fn from_runall_propagates_queue_memory_options() {
        // `Runall.queue_memory` carries the `--queue-memory` /
        // `--queue-memory-per-thread` flags. The chain runner must clone
        // those into its own state so it can build a `FastqPipelineConfig`
        // with the same memory budget the standalone `fgumi extract`
        // would compute.
        let dir = tempfile::TempDir::new().unwrap();
        let mut runall = runall_for_extract_chain(
            vec![dir.path().join("r1.fq.gz"), dir.path().join("r2.fq.gz")],
            dir.path().join("out.bam"),
        );
        runall.queue_memory = QueueMemoryOptions {
            queue_memory: "2GB".to_string(),
            queue_memory_per_thread: false,
            queue_memory_limit_mb: None,
        };

        let runner = ExtractChainRunner::from_runall(&runall);
        assert_eq!(runner.queue_memory.queue_memory, "2GB");
        assert!(!runner.queue_memory.queue_memory_per_thread);
    }

    #[test]
    fn build_pipeline_config_applies_queue_memory_limit() {
        // End-to-end check that `--queue-memory` reaches the
        // `FastqPipelineConfig`: 1 GiB total (per-thread = false) → exactly
        // 1 GiB queue-memory limit on the resulting config, regardless of
        // `--threads`.
        let dir = tempfile::TempDir::new().unwrap();
        let mut runall = runall_for_extract_chain(
            vec![dir.path().join("r1.fq.gz"), dir.path().join("r2.fq.gz")],
            dir.path().join("out.bam"),
        );
        runall.threads = 4;
        runall.queue_memory = QueueMemoryOptions {
            queue_memory: "1GiB".to_string(),
            queue_memory_per_thread: false,
            queue_memory_limit_mb: None,
        };

        let runner = ExtractChainRunner::from_runall(&runall);
        let config = runner.build_pipeline_config(true).expect("config builds");
        assert_eq!(
            config.queue_memory_limit,
            1024 * 1024 * 1024,
            "1 GiB total → 1 GiB queue_memory_limit on the FastqPipelineConfig",
        );
    }

    #[test]
    fn build_pipeline_config_per_thread_multiplies_by_threads() {
        // Per-thread mode: 256 MiB × 4 threads → 1 GiB total. Confirms
        // the per-thread axis matches the standalone command's behavior.
        let dir = tempfile::TempDir::new().unwrap();
        let mut runall = runall_for_extract_chain(
            vec![dir.path().join("r1.fq.gz"), dir.path().join("r2.fq.gz")],
            dir.path().join("out.bam"),
        );
        runall.threads = 4;
        runall.queue_memory = QueueMemoryOptions {
            queue_memory: "256MiB".to_string(),
            queue_memory_per_thread: true,
            queue_memory_limit_mb: None,
        };

        let runner = ExtractChainRunner::from_runall(&runall);
        let config = runner.build_pipeline_config(true).expect("config builds");
        assert_eq!(
            config.queue_memory_limit,
            4 * 256 * 1024 * 1024,
            "256 MiB × 4 threads → 1 GiB queue_memory_limit",
        );
    }
}
