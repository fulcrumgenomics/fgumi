//! `ProcessingSink`: bridges sort merge output into UMI grouping + consensus.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{Context, Result};
use fgumi_consensus::filter::{
    FilterConfig, FilterResult, compute_read_stats_raw, filter_duplex_read_raw, filter_read_raw,
    is_duplex_consensus_raw, mask_bases_raw, mask_duplex_bases_raw,
};
use fgumi_lib::bam_io::RawBamWriter;
use fgumi_lib::consensus_caller::{ConsensusCaller, ConsensusOutput};
use fgumi_lib::grouper::{
    FilterMetrics, GroupFilterConfig, extract_umi_from_template, filter_template_raw,
    get_pair_orientation_raw, group_by_queryname,
};
use fgumi_lib::progress::ProgressTracker;
use fgumi_lib::sort::inline_buffer::TemplateKey;
use fgumi_lib::sort::sink::SortedRecordSink;
use fgumi_lib::template::Template;
use fgumi_metrics::downsampling::compute_hash_fraction;
use fgumi_metrics::inline_collector::InlineCollector;
use fgumi_metrics::template_info::{ReadInfoKey, TemplateInfo};
use fgumi_raw_bam::fields::{aux_data_slice, l_seq};
use fgumi_umi::{Umi, UmiAssigner};
use log::info;
use noodles::sam::Header;

use super::options::StopAfter;

pub(super) struct ProcessingSink {
    pub(super) writer: Option<RawBamWriter>,
    pub(super) caller: Box<dyn ConsensusCaller>,
    assigner: Box<dyn UmiAssigner>,
    umi_tag: [u8; 2],
    /// Filter configuration for consensus read filtering/masking.
    filter_config: FilterConfig,
    /// Require single-strand agreement for duplex base masking.
    require_ss_agreement: bool,
    /// Filter configuration for template filtering before UMI assignment
    /// (MAPQ, unmapped, non-PF, UMI validation — matching standalone group command).
    group_filter_config: GroupFilterConfig,
    /// Stage to stop after.
    stop_after: StopAfter,
    /// Cancellation flag set by signal handler (SIGINT/SIGTERM).
    cancel: Arc<AtomicBool>,

    /// Records accumulated for the current position group.
    current_batch: Vec<Vec<u8>>,
    /// `(primary, secondary, cb_hash)` of the current position group.
    current_pos_key: Option<(u64, u64, u64)>,

    /// Running count of consensus reads written.
    pub(super) total_consensus_reads: usize,
    /// Number of position groups processed.
    pub(super) position_groups_processed: u64,
    /// Total input records seen across all position groups.
    pub(super) total_input_records: u64,
    /// Total molecule groups (MI assignments) processed.
    total_mi_groups: u64,
    /// Next global molecule ID offset for unique MI values across position groups.
    next_mi: u64,
    /// Progress tracker for consensus reads (periodic logging).
    progress: ProgressTracker,
    /// Progress tracker for position groups (periodic logging).
    group_progress: ProgressTracker,
    /// Optional inline metrics collector (enabled via --consensus-metrics).
    metrics_collector: Option<InlineCollector>,
    /// Path prefix for metrics output files.
    metrics_prefix: Option<PathBuf>,
    /// BAM header for `ref_name` lookups during metrics collection.
    header: Option<Header>,
}

impl ProcessingSink {
    #[expect(clippy::too_many_arguments, reason = "sink needs all stage configs")]
    pub(super) fn new(
        writer: RawBamWriter,
        caller: Box<dyn ConsensusCaller>,
        assigner: Box<dyn UmiAssigner>,
        umi_tag: [u8; 2],
        filter_config: FilterConfig,
        require_ss_agreement: bool,
        group_filter_config: GroupFilterConfig,
        stop_after: StopAfter,
        cancel: Arc<AtomicBool>,
        metrics_collector: Option<InlineCollector>,
        metrics_prefix: Option<PathBuf>,
        header: Option<Header>,
    ) -> Self {
        Self {
            writer: Some(writer),
            caller,
            assigner,
            umi_tag,
            filter_config,
            require_ss_agreement,
            group_filter_config,
            stop_after,
            cancel,
            current_batch: Vec::new(),
            current_pos_key: None,
            total_consensus_reads: 0,
            position_groups_processed: 0,
            total_input_records: 0,
            total_mi_groups: 0,
            next_mi: 0,
            progress: ProgressTracker::new("Consensus reads written").with_interval(100_000),
            group_progress: ProgressTracker::new("Position groups processed").with_interval(10_000),
            metrics_collector,
            metrics_prefix,
            header,
        }
    }

    /// Process a completed position group: group by queryname, build templates,
    /// assign UMIs, call consensus on each molecule group, and write output.
    fn process_position_group(&mut self) -> Result<()> {
        if self.current_batch.is_empty() {
            return Ok(());
        }

        let records = std::mem::take(&mut self.current_batch);
        self.position_groups_processed += 1;
        self.total_input_records += records.len() as u64;
        self.group_progress.log_if_needed(1);

        // 1. Group records by queryname (consecutive same-name runs).
        let record_groups = group_by_queryname(records);
        if record_groups.is_empty() {
            return Ok(());
        }

        // 2. Build Templates from each name group.
        let mut templates: Vec<Template> = Vec::with_capacity(record_groups.len());
        for group in record_groups {
            templates.push(Template::from_raw_records(group)?);
        }

        // 2b. Filter templates (same criteria as standalone group command:
        //     both-unmapped, MAPQ < min_mapq, non-PF, UMI validation).
        let mut filter_metrics = FilterMetrics::default();
        templates
            .retain(|t| filter_template_raw(t, &self.group_filter_config, &mut filter_metrics));

        if templates.is_empty() {
            return Ok(());
        }

        // 3. Extract UMIs and run UMI assignment, optionally splitting by pair orientation.
        //    Templates with empty or missing UMIs are skipped (retain MoleculeId::None).
        if self.assigner.split_templates_by_pair_orientation() {
            let mut subgroups: HashMap<(bool, bool), Vec<usize>> = HashMap::new();
            for (idx, template) in templates.iter().enumerate() {
                let orientation = get_pair_orientation_raw(template);
                subgroups.entry(orientation).or_default().push(idx);
            }

            for indices in subgroups.values() {
                let (valid_indices, umis): (Vec<usize>, Vec<Umi>) = indices
                    .iter()
                    .filter_map(|&i| {
                        let umi = extract_umi_from_template(&templates[i], self.umi_tag).ok()?;
                        if umi.is_empty() { None } else { Some((i, umi)) }
                    })
                    .unzip();
                if umis.is_empty() {
                    continue;
                }
                let molecule_ids = self.assigner.assign(&umis);
                for (idx, mi) in valid_indices.into_iter().zip(molecule_ids) {
                    templates[idx].mi = mi;
                }
            }
        } else {
            let (valid_indices, umis): (Vec<usize>, Vec<Umi>) = templates
                .iter()
                .enumerate()
                .filter_map(|(i, t)| {
                    let umi = extract_umi_from_template(t, self.umi_tag).ok()?;
                    if umi.is_empty() { None } else { Some((i, umi)) }
                })
                .unzip();
            if !umis.is_empty() {
                let molecule_ids = self.assigner.assign(&umis);
                for (idx, mi) in valid_indices.into_iter().zip(molecule_ids) {
                    templates[idx].mi = mi;
                }
            }
        }

        // 4. Collect inline metrics if enabled.
        if let Some(ref mut collector) = self.metrics_collector {
            let template_infos =
                build_template_infos(&templates, self.umi_tag, self.header.as_ref());
            collector.record_coordinate_group(&template_infos);
        }

        // --stop-after group: write records with MI tags, grouped by MI, and skip
        // consensus.  Records must be consecutive by MI so that `RawMiGroupIterator`
        // can re-read them correctly when resuming with `--start-from group`.
        if self.stop_after == StopAfter::Group {
            use fgumi_lib::sort::bam_fields;
            let mi_tag_bytes = [b'M', b'I'];

            // Compute global MI offset so MI values are unique across position groups.
            let max_local_id =
                templates.iter().filter_map(|t| t.mi.id()).max().map_or(0, |id| id + 1);
            let base_mi = self.next_mi;
            self.next_mi += max_local_id;

            // Collect records by MI so same-MI records are written consecutively.
            let mut mi_buf = String::with_capacity(16);
            let mut mi_groups: HashMap<u64, Vec<Vec<u8>>> = HashMap::new();
            let mut unassigned: Vec<Vec<u8>> = Vec::new();
            for template in templates {
                let mi = template.mi;
                if let Some(id) = mi.id() {
                    mi.write_with_offset(base_mi, &mut mi_buf);
                    let mut raw_records = template.into_raw_records().unwrap_or_default();
                    for record in &mut raw_records {
                        bam_fields::update_string_tag(record, &mi_tag_bytes, mi_buf.as_bytes());
                    }
                    mi_groups.entry(id).or_default().extend(raw_records);
                } else {
                    unassigned.extend(template.into_raw_records().unwrap_or_default());
                }
            }

            let writer = self.writer.as_mut().expect("write called after finish");
            for records in mi_groups.into_values() {
                for record in &records {
                    writer.write_raw_record(record).context("Failed to write grouped record")?;
                    self.total_consensus_reads += 1;
                }
            }
            for record in &unassigned {
                writer.write_raw_record(record).context("Failed to write grouped record")?;
                self.total_consensus_reads += 1;
            }
            self.progress.log_if_needed(1);
            return Ok(());
        }

        // 5. Inject MI tags into raw records and regroup by molecule ID.
        use fgumi_lib::sort::bam_fields;
        let mi_tag_bytes = [b'M', b'I'];
        let mut mi_buf = String::with_capacity(16);
        let mut mi_map: HashMap<u64, Vec<Vec<u8>>> = HashMap::new();
        for template in templates {
            let mi = template.mi;
            if let Some(id) = mi.id() {
                mi.write_with_offset(0, &mut mi_buf);
                let mut raw_records = template.into_raw_records().unwrap_or_default();
                for record in &mut raw_records {
                    bam_fields::update_string_tag(record, &mi_tag_bytes, mi_buf.as_bytes());
                }
                mi_map.entry(id).or_default().extend(raw_records);
            }
            // Templates with MoleculeId::None are dropped (unassigned).
        }
        self.total_mi_groups += mi_map.len() as u64;

        // 6. Call consensus on each molecule group, then filter.
        let duplex_t = self.filter_config.duplex_thresholds();
        let ss_t = self.filter_config.effective_single_strand_thresholds();
        let writer = self.writer.as_mut().expect("write called after finish");

        for (_mi, records) in mi_map {
            let output: ConsensusOutput =
                self.caller.consensus_reads(records).context("Failed to call consensus")?;

            if output.count == 0 {
                continue;
            }

            // Parse consensus output and apply filtering to each record.
            let mut data = output.data;
            let mut offset = 0;
            while offset + 4 <= data.len() {
                let block_size = u32::from_le_bytes(
                    data[offset..offset + 4].try_into().expect("slice is exactly 4 bytes"),
                ) as usize;
                offset += 4;
                if offset + block_size > data.len() {
                    break;
                }

                let record = &mut data[offset..offset + block_size];

                // Apply masking + filtering (same logic as FilterStage).
                let aux = aux_data_slice(record);
                let is_duplex = is_duplex_consensus_raw(aux);

                if is_duplex {
                    if let Some((cc, ab, ba)) = duplex_t {
                        mask_duplex_bases_raw(
                            record,
                            cc,
                            ab,
                            ba,
                            self.filter_config.min_base_quality,
                            self.require_ss_agreement,
                        )?;
                    }
                } else if let Some(thresholds) = ss_t {
                    mask_bases_raw(record, thresholds, self.filter_config.min_base_quality)?;
                }

                let aux = aux_data_slice(record);
                let filter_result = if is_duplex {
                    if let Some((cc, ab, ba)) = duplex_t {
                        filter_duplex_read_raw(aux, cc, ab, ba)?
                    } else {
                        FilterResult::Pass
                    }
                } else if let Some(thresholds) = ss_t {
                    filter_read_raw(aux, thresholds)?
                } else {
                    FilterResult::Pass
                };

                if filter_result != FilterResult::Pass {
                    offset += block_size;
                    continue;
                }

                // Global checks: no-call fraction, mean base quality.
                let (no_call_count, mean_base_qual) = compute_read_stats_raw(record);
                let seq_len = l_seq(record) as usize;
                if seq_len > 0 {
                    #[expect(clippy::cast_precision_loss, reason = "acceptable for fraction")]
                    let no_call_fraction = no_call_count as f64 / seq_len as f64;
                    if no_call_fraction > self.filter_config.max_no_call_fraction {
                        offset += block_size;
                        continue;
                    }
                }
                if let Some(min_mean_qual) = self.filter_config.min_mean_base_quality {
                    if mean_base_qual < min_mean_qual {
                        offset += block_size;
                        continue;
                    }
                }

                // Record passes — write with block_size prefix.
                writer.write_raw_record(record).context("Failed to write filtered record")?;
                self.total_consensus_reads += 1;
                offset += block_size;
            }

            self.progress.log_if_needed(1);
        }

        Ok(())
    }
}

/// Converts templates with assigned MIs into `TemplateInfo` for metrics collection.
fn build_template_infos(
    templates: &[Template],
    umi_tag: [u8; 2],
    header: Option<&Header>,
) -> Vec<TemplateInfo> {
    let mut mi_buf = String::with_capacity(16);
    templates
        .iter()
        .filter(|t| t.mi.id().is_some())
        .map(|t| {
            t.mi.write_with_offset(0, &mut mi_buf);
            let mi = mi_buf.clone();
            let rx = extract_umi_from_template(t, umi_tag).unwrap_or_default();
            let (ref_name, position, end_position, read_info_key) =
                extract_position_info(t, header);
            let read_name = String::from_utf8_lossy(&t.name);
            let hash_fraction = compute_hash_fraction(&read_name);

            TemplateInfo { mi, rx, ref_name, position, end_position, hash_fraction, read_info_key }
        })
        .collect()
}

/// Extracts reference name, insert bounds, and `ReadInfoKey` from a template's
/// raw R1 and R2 records.
///
/// The `ReadInfoKey` uses unclipped 5' positions and strand from both mates,
/// ordered so the earlier-mapping read comes first (matching the standalone
/// metrics path).
#[expect(
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    reason = "genomic positions never exceed i32::MAX; ref_id is validated non-negative"
)]
fn extract_position_info(
    template: &Template,
    header: Option<&Header>,
) -> (Option<String>, Option<i32>, Option<i32>, ReadInfoKey) {
    use fgumi_raw_bam::cigar::{alignment_end_from_raw, alignment_start_from_raw};
    use fgumi_raw_bam::fields::ref_id;

    let r1 = template.raw_r1();
    let r2 = template.raw_r2();

    let ref_name = [r1, r2].into_iter().flatten().find_map(|raw| {
        let rid = ref_id(raw);
        if rid < 0 {
            return None;
        }
        header.and_then(|h| {
            h.reference_sequences().get_index(rid as usize).map(|(name, _)| name.to_string())
        })
    });

    // Compute insert start (min of R1/R2 alignment starts) and
    // end (max of R1/R2 alignment ends), using 1-based coordinates.
    let r1_start = r1.and_then(|raw| alignment_start_from_raw(raw).map(|p| p as i32));
    let r2_start = r2.and_then(|raw| alignment_start_from_raw(raw).map(|p| p as i32));
    let r1_end = r1.and_then(|raw| alignment_end_from_raw(raw).map(|p| p as i32));
    let r2_end = r2.and_then(|raw| alignment_end_from_raw(raw).map(|p| p as i32));

    let position = match (r1_start, r2_start) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) | (None, Some(a)) => Some(a),
        (None, None) => None,
    };

    let end_position = match (r1_end, r2_end) {
        (Some(a), Some(b)) => Some(a.max(b)),
        (Some(a), None) | (None, Some(a)) => Some(a),
        (None, None) => None,
    };

    // Build ReadInfoKey from unclipped 5' positions and strand, matching fgbio's ReadInfo.
    let read_info_key = build_read_info_key(r1, r2);

    (ref_name, position, end_position, read_info_key)
}

/// Builds a [`ReadInfoKey`] from raw R1 and R2 BAM records.
///
/// Uses unclipped 5' positions, reference indices, and strand from both mates.
/// Fields are ordered so the earlier-mapping read comes first (matching fgbio).
/// Returns `ReadInfoKey::default()` if neither mate is mapped.
#[expect(clippy::cast_sign_loss, reason = "ref_id validated non-negative before cast")]
fn build_read_info_key(r1: Option<&[u8]>, r2: Option<&[u8]>) -> ReadInfoKey {
    use fgumi_raw_bam::cigar::unclipped_5prime_from_raw_bam;
    use fgumi_raw_bam::fields::{flags as raw_flags, ref_id};

    /// Extract (`ref_index`, `unclipped_5prime`, `is_reverse`) from a raw record, or `None` if unmapped.
    fn mate_info(raw: &[u8]) -> Option<(usize, i32, bool)> {
        let rid = ref_id(raw);
        if rid < 0 {
            return None;
        }
        let flg = raw_flags(raw);
        let unmapped = (flg & fgumi_raw_bam::fields::flags::UNMAPPED) != 0;
        if unmapped {
            return None;
        }
        let five_prime = unclipped_5prime_from_raw_bam(raw);
        let reverse = (flg & fgumi_raw_bam::fields::flags::REVERSE) != 0;
        Some((rid as usize, five_prime, reverse))
    }

    let r1_info = r1.and_then(mate_info);
    let r2_info = r2.and_then(mate_info);

    match (r1_info, r2_info) {
        (Some((ref1, s1, strand1)), Some((ref2, s2, strand2))) => {
            // Order so earlier-mapping read comes first
            if (Some(ref1), s1) <= (Some(ref2), s2) {
                ReadInfoKey {
                    ref_index1: Some(ref1),
                    start1: s1,
                    strand1,
                    ref_index2: Some(ref2),
                    start2: s2,
                    strand2,
                }
            } else {
                ReadInfoKey {
                    ref_index1: Some(ref2),
                    start1: s2,
                    strand1: strand2,
                    ref_index2: Some(ref1),
                    start2: s1,
                    strand2: strand1,
                }
            }
        }
        (Some((ref1, s1, strand1)), None) => {
            ReadInfoKey { ref_index1: Some(ref1), start1: s1, strand1, ..ReadInfoKey::default() }
        }
        (None, Some((ref2, s2, strand2))) => ReadInfoKey {
            ref_index1: Some(ref2),
            start1: s2,
            strand1: strand2,
            ..ReadInfoKey::default()
        },
        (None, None) => ReadInfoKey::default(),
    }
}

impl SortedRecordSink for ProcessingSink {
    fn emit(&mut self, key: &TemplateKey, record_bytes: Vec<u8>) -> Result<()> {
        if self.cancel.load(Ordering::Relaxed) {
            anyhow::bail!("Pipeline cancelled by signal");
        }

        let pos_key = (key.primary, key.secondary, key.cb_hash);

        if self.current_pos_key != Some(pos_key) {
            // Position changed: process the accumulated batch.
            self.process_position_group()?;
            self.current_pos_key = Some(pos_key);
        }

        self.current_batch.push(record_bytes);
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        // Process the final position group.
        self.process_position_group()?;
        self.progress.log_final();
        self.group_progress.log_final();
        info!(
            "Group/consensus complete: {} input records → {} position groups → {} MI groups → {} output reads",
            self.total_input_records,
            self.position_groups_processed,
            self.total_mi_groups,
            self.total_consensus_reads
        );
        // Write inline metrics if enabled.
        if let (Some(collector), Some(prefix)) =
            (self.metrics_collector.as_ref(), self.metrics_prefix.as_ref())
        {
            collector.write_metrics(prefix)?;
            info!("Consensus metrics written to {}.* ", prefix.display());
        }

        if let Some(writer) = self.writer.take() {
            writer.finish().context("Failed to finish output BAM")?;
        }
        Ok(())
    }
}
