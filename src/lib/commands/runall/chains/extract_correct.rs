//! `--start-from extract --stop-after correct` chain runner.
//!
//! Fuses the extract and correct steps into a single
//! [`crate::unified_pipeline::run_fastq_pipeline`] `process_fn` closure: each
//! `FastqTemplate` is converted to BAM bytes via
//! [`crate::commands::extract::make_raw_records_static`] and the per-template
//! UMI correction is applied in-place to those bytes. No intermediate BAM is
//! materialised between the two stages.
//!
//! The closure body is `Fn` (not `FnMut`); per-worker mutable state — the LRU
//! cache used by [`CorrectUmis::compute_template_correction`] — lives in
//! `thread_local!` storage, mirroring the standalone `fgumi correct`
//! pipeline.
//!
//! When the caller passes `--metrics`, the runner times its
//! `run_fastq_pipeline` call and emits two
//! [`crate::commands::runall::infra::metrics::MetricsRow`]s to the supplied
//! TSV: an `extract` row followed by a `correct` row. Per-stage wall time
//! and `records_in` / `records_out` are simplifications until the fused
//! `process_fn` grows per-stage accumulators (PR #329 wires
//! `--correct::metrics`, which exposes that infrastructure):
//!
//! * Both rows share the total wall-clock duration of the fused pipeline.
//! * Both rows have `records_in == records_out` set to the post-correction
//!   record count emitted by `run_fastq_pipeline`. (Without
//!   `--correct::rejects` rejected templates are silently dropped, but the
//!   chain runner does not expose a separate dropped-record counter today;
//!   PR #329 adds that path.)

use std::cell::RefCell;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::Arc;

use anyhow::{Context, Result};
use lru::LruCache;
use noodles::sam::Header;
use read_structure::ReadStructure;

use crate::commands::common::{QueueMemoryOptions, add_pg_record};
use crate::commands::correct::{CorrectUmis, EncodedUmiSet, TemplateCorrection, UmiMatch};
use crate::commands::extract::{
    BUFFER_SIZE, CompressionFormat, ExtractConfig, ExtractHeaderMetadata,
    QUALITY_DETECTION_SAMPLE_SIZE, QualityEncoding, RawExtractedRecords, build_extract_header,
    detect_compression_format, make_raw_records_static, open_fastq_reader,
    strip_read_suffix_extract,
};
use crate::commands::runall::Runall;
use crate::commands::runall::chains::correct::load_umi_sequences;
use crate::commands::runall::dispatch::{ChainContext, ChainRunner, DispatchContext};
use crate::commands::runall::infra::metrics::{MetricsRow, MetricsWriter};
use crate::commands::runall::options::{StartFrom, StopAfter};
use crate::fastq::FastqSet;
use crate::grouper::FastqTemplate;
use crate::sam::SamTag;
use crate::unified_pipeline::{FastqPipelineConfig, run_fastq_pipeline};
use fgumi_raw_bam::{aux_data_slice, find_string_tag, update_string_tag};
use fgumi_simd_fastq::SimdFastqReader;

/// Stable runner name surfaced via `--explain` and debug logs.
pub(crate) const EXTRACT_CORRECT_CHAIN_RUNNER_NAME: &str = "ExtractCorrectChainRunner";

/// Chain runner for `--start-from extract --stop-after correct`.
///
/// Built once per `Runall::execute` call; owns the cloned subset of
/// `Runall` fields needed to construct the unmapped BAM header, the
/// extract closure config, and the correction parameters.
#[allow(clippy::struct_excessive_bools)] // many distinct extract + correct flags
pub(crate) struct ExtractCorrectChainRunner {
    // Extract inputs
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
    // Correction inputs
    umis: Vec<String>,
    umi_files: Vec<PathBuf>,
    max_mismatches: usize,
    min_distance_diff: usize,
    cache_size: usize,
    revcomp: bool,
    dont_store_original_umis: bool,
    // Chain-runner-incompatible flags (gated off; PR #329 will wire them).
    rejects_unset: bool,
    metrics_unset: bool,
    min_corrected_unset: bool,
    /// `--queue-memory` / `--queue-memory-per-thread` from the top-level
    /// `Runall`. Cloned so the runner can build a `FastqPipelineConfig`
    /// with the same memory budget the standalone `fgumi extract` /
    /// `fgumi correct` commands would for the same flag values.
    queue_memory: QueueMemoryOptions,
}

impl ExtractCorrectChainRunner {
    /// Build the runner from a parsed [`Runall`] invocation.
    pub(crate) fn from_runall(runall: &Runall) -> Self {
        let extract_opts = &runall.extract_opts;
        let correct_opts = &runall.correct_opts;
        Self {
            inputs: runall.input.clone(),
            threads: runall.threads,
            compression_level: runall.compression_level,
            read_structures: extract_opts.extract_read_structures.clone().unwrap_or_default(),
            sample: extract_opts.extract_sample.clone(),
            library: extract_opts.extract_library.clone(),
            read_group_id: runall.consensus_opts.consensus_read_group_id.clone(),
            barcode: extract_opts.extract_barcode.clone(),
            platform: extract_opts.extract_platform.clone(),
            platform_unit: extract_opts.extract_platform_unit.clone(),
            platform_model: extract_opts.extract_platform_model.clone(),
            sequencing_center: extract_opts.extract_sequencing_center.clone(),
            predicted_insert_size: extract_opts.extract_predicted_insert_size,
            description: extract_opts.extract_description.clone(),
            run_date: extract_opts.extract_run_date.clone(),
            comments: extract_opts.extract_comment.clone().unwrap_or_default(),
            single_tag: extract_opts.extract_single_tag.clone(),
            annotate_read_names: extract_opts.extract_annotate_read_names,
            extract_umis_from_read_names: extract_opts.extract_extract_umis_from_read_names,
            store_sample_barcode_qualities: extract_opts.extract_store_sample_barcode_qualities,
            umis: correct_opts.correct_umis.clone().unwrap_or_default(),
            umi_files: correct_opts.correct_umi_files.clone().unwrap_or_default(),
            max_mismatches: correct_opts.correct_max_mismatches,
            min_distance_diff: correct_opts.correct_min_distance,
            cache_size: correct_opts.correct_cache_size,
            revcomp: correct_opts.correct_revcomp,
            dont_store_original_umis: correct_opts.correct_dont_store_original_umis,
            rejects_unset: correct_opts.correct_rejects.is_none(),
            metrics_unset: correct_opts.correct_metrics.is_none(),
            min_corrected_unset: correct_opts.correct_min_corrected.is_none(),
            queue_memory: runall.queue_memory.clone(),
        }
    }

    /// Build the [`FastqPipelineConfig`] for `run_chain`, applying the
    /// queue-memory budget computed from `--queue-memory` /
    /// `--queue-memory-per-thread`. Mirrors the helper on
    /// [`super::extract::ExtractChainRunner`] and the standalone
    /// `fgumi extract` flow. Extracted for unit-testing the wiring without
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

    /// Detect the FASTQ quality encoding by sampling the leading records of
    /// the first input. Mirrors the standalone `Extract` flow.
    fn detect_quality_encoding(&self) -> Result<QualityEncoding> {
        let first =
            self.inputs.first().ok_or_else(|| anyhow::anyhow!("no FASTQ inputs supplied"))?;
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

    /// Build the fused per-template `process_fn` for `run_fastq_pipeline`.
    ///
    /// The closure runs once per [`FastqTemplate`]:
    /// 1. extract the template into BAM bytes via [`make_raw_records_static`];
    /// 2. read the template's UMI from the first emitted record (all records in
    ///    one template share the same UMI by construction);
    /// 3. compute the correction once via [`CorrectUmis::compute_template_correction`],
    ///    using a `thread_local!` LRU cache keyed by the original UMI bytes;
    /// 4. either drop the template (UMI not matched) or apply the correction
    ///    in-place to every record's RX (and OX) tag.
    ///
    /// The returned closure is `Send + Sync + 'static` and plugs directly
    /// into the FASTQ pipeline's `process_fn` slot.
    #[allow(clippy::too_many_arguments)]
    fn build_process_fn(
        read_structures: Vec<ReadStructure>,
        encoding: QualityEncoding,
        cfg: ExtractConfig,
        encoded_umi_set: Arc<EncodedUmiSet>,
        umi_length: usize,
        revcomp: bool,
        max_mismatches: usize,
        min_distance_diff: usize,
        dont_store_original_umis: bool,
        cache_size: usize,
    ) -> impl Fn(FastqTemplate) -> io::Result<RawExtractedRecords> + Send + Sync + 'static {
        let read_structures = Arc::new(read_structures);
        let cfg = Arc::new(cfg);
        let umi_tag: [u8; 2] = *SamTag::RX;

        move |template: FastqTemplate| -> io::Result<RawExtractedRecords> {
            // Validate paired read names match across streams (mirrors the
            // standalone extract command's synchronized-mode check).
            if template.records.len() >= 2 {
                let base_name = strip_read_suffix_extract(template.records[0].name());
                for (i, record) in template.records.iter().enumerate().skip(1) {
                    let other_base = strip_read_suffix_extract(record.name());
                    if base_name != other_base {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
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

            // Step 1: extract.
            let mut fastq_sets: Vec<FastqSet> = Vec::with_capacity(template.records.len());
            for (record, rs) in template.records.iter().zip(read_structures.iter()) {
                let fastq_set = FastqSet::from_record_with_structure(
                    record.name(),
                    record.sequence(),
                    record.quality(),
                    rs,
                    &[],
                )
                .map_err(io::Error::other)?;
                fastq_sets.push(fastq_set);
            }
            let combined = FastqSet::combine_readsets(fastq_sets);
            let extracted = make_raw_records_static(&combined, encoding, &cfg)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;

            if extracted.num_records == 0 {
                return Ok(extracted);
            }

            // Step 2: read the UMI from the first record. If any record lacks
            // an RX tag, drop the whole template — matches standalone
            // `fgumi correct`'s behaviour for missing-UMI input with no
            // `--rejects` configured.
            let first_record = match iter_length_prefixed(&extracted.data).next() {
                Some(Ok(slice)) => slice,
                Some(Err(e)) => {
                    return Err(io::Error::other(format!("extract_correct: {e}")));
                }
                None => return Ok(extracted),
            };
            let aux = aux_data_slice(first_record);
            let Some(umi_bytes) = find_string_tag(aux, umi_tag) else {
                return Ok(RawExtractedRecords { data: Vec::new(), num_records: 0 });
            };
            let umi_str = match std::str::from_utf8(umi_bytes) {
                Ok(s) => s.to_owned(),
                Err(e) => {
                    return Err(io::Error::other(format!(
                        "extract_correct: RX tag is not valid UTF-8: {e}"
                    )));
                }
            };

            // Step 3: correction. Per-worker LRU cache lives in TLS so the
            // closure body can stay `Fn` (not `FnMut`).
            thread_local! {
                static CACHE: RefCell<Option<LruCache<Vec<u8>, UmiMatch>>> =
                    const { RefCell::new(None) };
            }
            let correction = CACHE.with(|cell| {
                let mut cache_ref = cell.borrow_mut();
                if cache_ref.is_none() && cache_size > 0 {
                    let cap =
                        NonZero::new(cache_size).expect("cache_size > 0 checked immediately above");
                    *cache_ref = Some(LruCache::new(cap));
                }
                CorrectUmis::compute_template_correction(
                    &umi_str,
                    umi_length,
                    revcomp,
                    max_mismatches,
                    min_distance_diff,
                    &encoded_umi_set,
                    &mut cache_ref,
                )
            });

            if !correction.matched {
                // Drop all records in this template (no `--rejects` support yet).
                return Ok(RawExtractedRecords { data: Vec::new(), num_records: 0 });
            }

            // Step 4: apply correction to every record in the template,
            // length-prefixing each record on the way out.
            let mut kept_data: Vec<u8> = Vec::with_capacity(extracted.data.len());
            let mut kept_count: u64 = 0;
            for record_slice in iter_length_prefixed(&extracted.data) {
                let record_slice =
                    record_slice.map_err(|e| io::Error::other(format!("extract_correct: {e}")))?;
                let mut record: Vec<u8> = record_slice.to_vec();
                apply_correction_to_record_bytes(
                    &mut record,
                    &correction,
                    umi_tag,
                    dont_store_original_umis,
                );
                let block_size = u32::try_from(record.len())
                    .map_err(|_| io::Error::other("extract_correct: record too large for u32"))?;
                kept_data.extend_from_slice(&block_size.to_le_bytes());
                kept_data.extend_from_slice(&record);
                kept_count += 1;
            }

            Ok(RawExtractedRecords { data: kept_data, num_records: kept_count })
        }
    }

    /// Identity `serialize_fn`: the fused process closure already produces
    /// length-prefixed BAM bytes (post-correction); the serializer just
    /// hands those bytes off into the worker scratch buffer.
    fn build_serialize_fn()
    -> impl Fn(RawExtractedRecords, &Header, &mut Vec<u8>) -> io::Result<u64> + Send + Sync + 'static
    {
        move |batch: RawExtractedRecords, _header: &Header, scratch: &mut Vec<u8>| {
            scratch.extend_from_slice(&batch.data);
            Ok(batch.num_records)
        }
    }

    /// Run the chain end-to-end against the prepared
    /// [`ChainContext::tmp_output`] path.
    ///
    /// When `metrics_path` is `Some`, writes two rows: an `extract` row
    /// followed by a `correct` row. Both rows currently share the total
    /// wall-clock duration of the fused pipeline and the same
    /// `records_in == records_out` value (the post-correction count
    /// emitted by `run_fastq_pipeline`). Per-stage timing and a separate
    /// pre-correction record count require the fused `process_fn` to grow
    /// per-stage accumulators — that infrastructure lands with PR #329's
    /// `--correct::metrics` wiring.
    fn run_chain(
        &self,
        command_line: &str,
        output: &Path,
        metrics_path: Option<&Path>,
    ) -> Result<()> {
        anyhow::ensure!(
            !self.inputs.is_empty(),
            "ExtractCorrectChainRunner: at least one --input required"
        );
        anyhow::ensure!(
            !self.umis.is_empty() || !self.umi_files.is_empty(),
            "ExtractCorrectChainRunner: --correct::umis or --correct::umi-files is required"
        );
        let sample = self.sample.as_deref().ok_or_else(|| {
            anyhow::anyhow!("--extract::sample is required for the extract→correct chain")
        })?;
        let library = self.library.as_deref().ok_or_else(|| {
            anyhow::anyhow!("--extract::library is required for the extract→correct chain")
        })?;

        // Detect compression format per input; require uniformity.
        let formats: Vec<CompressionFormat> = self
            .inputs
            .iter()
            .map(|p| {
                detect_compression_format(p)
                    .with_context(|| format!("detecting compression format for {}", p.display()))
            })
            .collect::<Result<_>>()?;
        let all_bgzf = formats.iter().all(|f| *f == CompressionFormat::Bgzf);

        let encoding = self.detect_quality_encoding()?;
        log::info!("ExtractCorrectChainRunner: detected quality encoding {encoding:?}");

        let read_structures = self.resolve_read_structures()?;
        anyhow::ensure!(
            read_structures.len() == self.inputs.len(),
            "ExtractCorrectChainRunner: --extract::read-structures must have one entry per --input \
             (got {} structures vs {} inputs)",
            read_structures.len(),
            self.inputs.len(),
        );

        // Load UMI list for correction.
        let (umi_sequences, umi_length) = load_umi_sequences(&self.umis, &self.umi_files)?;
        let encoded_umi_set = Arc::new(EncodedUmiSet::new(&umi_sequences));

        // Build the unmapped BAM header (matches `Extract::create_header`),
        // then stamp a `@PG` record for runall.
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
        let header = add_pg_record(header, command_line)?;

        // Build the closure config (mirrors `Extract`'s ExtractConfig
        // construction; matches the standalone defaults for fields we don't
        // expose at the runall layer yet).
        let cfg = ExtractConfig {
            read_group_id: self.read_group_id.clone(),
            store_umi_quals: false,
            store_cell_quals: false,
            single_tag: self.single_tag.as_ref().map(|t| {
                SamTag::from_str(t).expect("--extract::single-tag length validated by validate.rs")
            }),
            annotate_read_names: self.annotate_read_names,
            extract_umis_from_read_names: self.extract_umis_from_read_names,
            store_sample_barcode_qualities: self.store_sample_barcode_qualities,
        };

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

        let process_fn = Self::build_process_fn(
            read_structures,
            encoding,
            cfg,
            encoded_umi_set,
            umi_length,
            self.revcomp,
            self.max_mismatches,
            self.min_distance_diff,
            self.dont_store_original_umis,
            self.cache_size,
        );
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
        .with_context(|| {
            format!("ExtractCorrectChainRunner pipeline writing to {}", output.display())
        })?;
        let total_wall = start.elapsed().as_secs_f64();

        if let Some(path) = metrics_path {
            // The fused process_fn does not yet expose per-stage
            // accumulators, so both rows share the total wall time and a
            // single records_in / records_out value (the post-correction
            // record count). PR #329 wires per-stage counters via
            // `--correct::metrics`; once those land, this block should
            // split the timings and populate `extract.records_out` /
            // `correct.records_in` from the pre-correction count.
            let mut writer = MetricsWriter::create(path)
                .with_context(|| format!("opening --metrics path {} for write", path.display()))?;
            writer
                .write_row(&MetricsRow {
                    stage: "extract".into(),
                    wall_time_secs: total_wall,
                    records_in: records_out,
                    records_out,
                })
                .with_context(|| format!("writing extract metrics row to {}", path.display()))?;
            writer
                .write_row(&MetricsRow {
                    stage: "correct".into(),
                    wall_time_secs: total_wall,
                    records_in: records_out,
                    records_out,
                })
                .with_context(|| format!("writing correct metrics row to {}", path.display()))?;
        }

        Ok(())
    }
}

impl ChainRunner for ExtractCorrectChainRunner {
    fn name(&self) -> &'static str {
        EXTRACT_CORRECT_CHAIN_RUNNER_NAME
    }

    /// Match `(Extract, Correct)` when:
    /// * a UMI source is provided (inline `--correct::umis` or
    ///   `--correct::umi-files`); and
    /// * none of the `--correct::rejects`, `--correct::metrics`, or
    ///   `--correct::min-corrected` flags are set (those still gate the
    ///   runner off and fall through to `NotImplementedRunner`; PR #329
    ///   will wire them).
    ///
    /// The top-level `--metrics` flag does NOT gate this runner — it is
    /// honoured directly by `run_chain`, which writes one `extract` row
    /// followed by one `correct` row to the supplied TSV.
    fn supports(&self, ctx: &DispatchContext<'_>) -> bool {
        if !self.rejects_unset || !self.metrics_unset || !self.min_corrected_unset {
            return false;
        }
        if self.umis.is_empty() && self.umi_files.is_empty() {
            return false;
        }
        ctx.start_from == StartFrom::Extract && ctx.stop_after == StopAfter::Correct
    }

    fn run(&self, ctx: ChainContext<'_>) -> Result<()> {
        let metrics_path = ctx.dispatch.metrics_path;
        self.run_chain(ctx.command_line, &ctx.tmp_output, metrics_path)
    }
}

/// Apply UMI correction to a complete BAM record's bytes.
///
/// The closure operates on the raw bytes emitted by
/// [`make_raw_records_static`] (i.e. without the 4-byte block-size prefix);
/// the standalone `CorrectUmis::apply_correction_to_raw` operates on a
/// `&mut RawRecord`. They do the same work — update the RX (and optionally
/// OX) tags in place — but the chain runner avoids the `RawRecord` wrapper
/// because the FASTQ pipeline never reads back as `RawRecord`s. Behaviour
/// matches standalone correct exactly: only update tags when correction is
/// actually needed, store original UMI in OX only when there are real
/// mismatches.
fn apply_correction_to_record_bytes(
    record: &mut Vec<u8>,
    correction: &TemplateCorrection,
    umi_tag: [u8; 2],
    dont_store_original_umis: bool,
) {
    if correction.needs_correction {
        if let Some(ref corrected) = correction.corrected_umi {
            update_string_tag(record, umi_tag, corrected.as_bytes());
        }
        if !dont_store_original_umis && correction.has_mismatches {
            update_string_tag(record, *SamTag::OX, correction.original_umi.as_bytes());
        }
    }
}

/// Iterate over a length-prefixed BAM record stream as produced by
/// [`make_raw_records_static`]. Each item is the record's payload bytes
/// (no 4-byte block-size prefix). Returns an `Err` for truncated input.
fn iter_length_prefixed(data: &[u8]) -> LengthPrefixedIter<'_> {
    LengthPrefixedIter { data, cursor: 0 }
}

struct LengthPrefixedIter<'a> {
    data: &'a [u8],
    cursor: usize,
}

impl<'a> Iterator for LengthPrefixedIter<'a> {
    type Item = io::Result<&'a [u8]>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cursor >= self.data.len() {
            return None;
        }
        if self.cursor + 4 > self.data.len() {
            return Some(Err(io::Error::other("truncated length prefix")));
        }
        let prefix: [u8; 4] =
            self.data[self.cursor..self.cursor + 4].try_into().expect("4-byte slice");
        let block_size = u32::from_le_bytes(prefix) as usize;
        self.cursor += 4;
        let end = self.cursor + block_size;
        if end > self.data.len() {
            return Some(Err(io::Error::other("truncated record body")));
        }
        let slice = &self.data[self.cursor..end];
        self.cursor = end;
        Some(Ok(slice))
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

    fn runall_for_extract_correct(
        inputs: Vec<PathBuf>,
        output: PathBuf,
        umis: Vec<String>,
    ) -> Runall {
        let extract = ExtractOptions {
            sample: Some("S".to_string()),
            library: Some("L".to_string()),
            read_structures: Some(vec!["8M+T".to_string(), "+T".to_string()]),
            ..ExtractOptions::default()
        };
        let correct = CorrectOptions { umis: Some(umis), ..CorrectOptions::default() };
        Runall {
            input: inputs,
            output,
            reference: None,
            start_from: StartFrom::Extract,
            stop_after: StopAfter::Correct,
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
            correct_opts: MultiCorrectOptions::from(correct),
            duplex_opts: MultiDuplexOptions::from(DuplexOptions::default()),
            codec_opts: MultiCodecOptions::from(CodecOptions::default()),
            queue_memory: QueueMemoryOptions::default(),
        }
    }

    fn dispatch_ctx(start: StartFrom, stop: StopAfter) -> DispatchContext<'static> {
        DispatchContext { start_from: start, stop_after: stop, metrics_path: None }
    }

    /// Write a paired-end gzip-compressed FASTQ fixture with a cycling
    /// UMI per record. Mirrors the `write_paired_gzip_fastq*` helpers
    /// used in the extract chain runner tests; kept local so this module
    /// doesn't depend on test code in a sibling file.
    fn write_paired_gzip_fastq_with_cycling_umi(r1: &Path, r2: &Path, num_pairs: usize) {
        use flate2::Compression;
        use flate2::write::GzEncoder;
        const UMI_CYCLE: &[&[u8]] = &[b"AAAAAAAA", b"CCCCCCCC", b"GGGGGGGG", b"TTTTTTTT"];

        let f1 = File::create(r1).unwrap();
        let mut e1 = GzEncoder::new(f1, Compression::default());
        let f2 = File::create(r2).unwrap();
        let mut e2 = GzEncoder::new(f2, Compression::default());
        for i in 0..num_pairs {
            let umi = UMI_CYCLE[i % UMI_CYCLE.len()];
            let umi_str = std::str::from_utf8(umi).unwrap();
            let r1_seq = format!("{umi_str}AAAAAAAA");
            let r2_seq = "CCCCCCCC";
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
    fn supports_only_extract_to_correct_with_umi_source() {
        let dir = tempfile::TempDir::new().unwrap();
        let runall = runall_for_extract_correct(
            vec![dir.path().join("r1.fq.gz"), dir.path().join("r2.fq.gz")],
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        let runner = ExtractCorrectChainRunner::from_runall(&runall);

        assert!(runner.supports(&dispatch_ctx(StartFrom::Extract, StopAfter::Correct)));
        assert!(!runner.supports(&dispatch_ctx(StartFrom::Extract, StopAfter::Extract)));
        assert!(!runner.supports(&dispatch_ctx(StartFrom::Correct, StopAfter::Correct)));
        assert!(!runner.supports(&dispatch_ctx(StartFrom::Group, StopAfter::Group)));

        let metrics_path = dir.path().join("metrics.tsv");
        let with_metrics = DispatchContext {
            start_from: StartFrom::Extract,
            stop_after: StopAfter::Correct,
            metrics_path: Some(metrics_path.as_path()),
        };
        assert!(
            runner.supports(&with_metrics),
            "extract→correct chain runner now writes its own --metrics rows; the registry \
             must dispatch to it even when --metrics is set"
        );
    }

    #[test]
    fn does_not_support_when_correct_rejects_set() {
        let dir = tempfile::TempDir::new().unwrap();
        let mut runall = runall_for_extract_correct(
            vec![dir.path().join("r1.fq.gz"), dir.path().join("r2.fq.gz")],
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        runall.correct_opts.correct_rejects = Some(dir.path().join("rejects.bam"));
        let runner = ExtractCorrectChainRunner::from_runall(&runall);
        assert!(
            !runner.supports(&dispatch_ctx(StartFrom::Extract, StopAfter::Correct)),
            "--correct::rejects must keep gating the runner off (PR #329)"
        );
    }

    #[test]
    fn does_not_support_when_correct_metrics_set() {
        let dir = tempfile::TempDir::new().unwrap();
        let mut runall = runall_for_extract_correct(
            vec![dir.path().join("r1.fq.gz"), dir.path().join("r2.fq.gz")],
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        runall.correct_opts.correct_metrics = Some(dir.path().join("correct-metrics.tsv"));
        let runner = ExtractCorrectChainRunner::from_runall(&runall);
        assert!(
            !runner.supports(&dispatch_ctx(StartFrom::Extract, StopAfter::Correct)),
            "--correct::metrics must keep gating the runner off (PR #329)"
        );
    }

    #[test]
    fn does_not_support_when_correct_min_corrected_set() {
        let dir = tempfile::TempDir::new().unwrap();
        let mut runall = runall_for_extract_correct(
            vec![dir.path().join("r1.fq.gz"), dir.path().join("r2.fq.gz")],
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        runall.correct_opts.correct_min_corrected = Some(0.5);
        let runner = ExtractCorrectChainRunner::from_runall(&runall);
        assert!(
            !runner.supports(&dispatch_ctx(StartFrom::Extract, StopAfter::Correct)),
            "--correct::min-corrected must keep gating the runner off (PR #329)"
        );
    }

    #[test]
    fn from_runall_propagates_queue_memory_options() {
        let dir = tempfile::TempDir::new().unwrap();
        let mut runall = runall_for_extract_correct(
            vec![dir.path().join("r1.fq.gz"), dir.path().join("r2.fq.gz")],
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        runall.queue_memory = QueueMemoryOptions {
            queue_memory: "2GB".to_string(),
            queue_memory_per_thread: false,
            queue_memory_limit_mb: None,
        };

        let runner = ExtractCorrectChainRunner::from_runall(&runall);
        assert_eq!(runner.queue_memory.queue_memory, "2GB");
        assert!(!runner.queue_memory.queue_memory_per_thread);
    }

    #[test]
    fn build_pipeline_config_applies_queue_memory_limit() {
        // 1 GiB total (per-thread = false) → exactly 1 GiB queue-memory
        // limit on the resulting `FastqPipelineConfig`, regardless of
        // `--threads`.
        let dir = tempfile::TempDir::new().unwrap();
        let mut runall = runall_for_extract_correct(
            vec![dir.path().join("r1.fq.gz"), dir.path().join("r2.fq.gz")],
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        runall.threads = 4;
        runall.queue_memory = QueueMemoryOptions {
            queue_memory: "1GiB".to_string(),
            queue_memory_per_thread: false,
            queue_memory_limit_mb: None,
        };

        let runner = ExtractCorrectChainRunner::from_runall(&runall);
        let config = runner.build_pipeline_config(true).expect("config builds");
        assert_eq!(
            config.queue_memory_limit,
            1024 * 1024 * 1024,
            "1 GiB total → 1 GiB queue_memory_limit on the FastqPipelineConfig",
        );
    }

    #[test]
    fn metrics_tsv_emits_two_rows() {
        // Wire `--metrics` end-to-end through `run_chain`. Asserts header,
        // ordering (`extract` then `correct`), and shared timing /
        // record-count values per the simplification documented on
        // `run_chain`.
        let dir = tempfile::TempDir::new().unwrap();
        let r1 = dir.path().join("R1.fastq.gz");
        let r2 = dir.path().join("R2.fastq.gz");
        let out_tmp = dir.path().join("extract_correct.bam.tmp");
        let metrics = dir.path().join("metrics.tsv");

        // 8 records, all UMIs in the allowed list — no template is dropped
        // by correction, so per-row records_in == records_out == 16
        // (8 templates × 2 segments).
        write_paired_gzip_fastq_with_cycling_umi(&r1, &r2, 8);

        let runall = runall_for_extract_correct(
            vec![r1.clone(), r2.clone()],
            out_tmp.clone(),
            vec![
                "AAAAAAAA".to_string(),
                "CCCCCCCC".to_string(),
                "GGGGGGGG".to_string(),
                "TTTTTTTT".to_string(),
            ],
        );
        let runner = ExtractCorrectChainRunner::from_runall(&runall);
        runner
            .run_chain("fgumi runall", &out_tmp, Some(metrics.as_path()))
            .expect("chain runs with --metrics");

        let contents = std::fs::read_to_string(&metrics).expect("read metrics");
        let mut lines = contents.lines();
        assert_eq!(
            lines.next().unwrap(),
            "stage\twall_time_secs\trecords_in\trecords_out",
            "header row"
        );
        let extract_row = lines.next().expect("extract row");
        let correct_row = lines.next().expect("correct row");
        let extract_cols: Vec<&str> = extract_row.split('\t').collect();
        let correct_cols: Vec<&str> = correct_row.split('\t').collect();

        assert_eq!(extract_cols[0], "extract", "row 1 stage");
        assert_eq!(correct_cols[0], "correct", "row 2 stage");
        // Both rows share the total wall time per the documented
        // simplification — assert the timestamps are byte-for-byte equal
        // so a future per-stage split is forced through this test.
        assert_eq!(extract_cols[1], correct_cols[1], "wall_time_secs is shared today");
        assert_eq!(extract_cols[2], "16", "extract.records_in: 8 templates × 2 segments");
        assert_eq!(extract_cols[3], "16", "extract.records_out");
        assert_eq!(correct_cols[2], "16", "correct.records_in");
        assert_eq!(correct_cols[3], "16", "correct.records_out");
        assert!(lines.next().is_none(), "exactly two data rows");
    }

    #[test]
    fn does_not_support_when_no_umi_source() {
        let dir = tempfile::TempDir::new().unwrap();
        let runall = runall_for_extract_correct(
            vec![dir.path().join("r1.fq.gz"), dir.path().join("r2.fq.gz")],
            dir.path().join("out.bam"),
            Vec::new(),
        );
        let runner = ExtractCorrectChainRunner::from_runall(&runall);
        assert!(!runner.supports(&dispatch_ctx(StartFrom::Extract, StopAfter::Correct)));
    }

    #[test]
    fn name_is_stable() {
        let dir = tempfile::TempDir::new().unwrap();
        let runall = runall_for_extract_correct(
            vec![dir.path().join("r1.fq.gz"), dir.path().join("r2.fq.gz")],
            dir.path().join("out.bam"),
            vec!["AAAAAAAA".to_string()],
        );
        let runner = ExtractCorrectChainRunner::from_runall(&runall);
        assert_eq!(runner.name(), EXTRACT_CORRECT_CHAIN_RUNNER_NAME);
    }
}
