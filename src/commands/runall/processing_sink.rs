//! `ProcessingSink`: bridges sort merge output into UMI grouping + consensus.

use std::collections::HashMap;
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
use fgumi_raw_bam::fields::{aux_data_slice, l_seq};
use fgumi_umi::{Umi, UmiAssigner};
use log::info;

use super::options::StopAfter;

pub(super) struct ProcessingSink {
    pub(super) writer: Option<RawBamWriter>,
    pub(super) caller: Box<dyn ConsensusCaller>,
    assigner: Box<dyn UmiAssigner>,
    umi_tag: [u8; 2],
    /// Filter configuration for consensus read filtering/masking.
    filter_config: FilterConfig,
    /// Filter configuration for template filtering before UMI assignment
    /// (MAPQ, unmapped, non-PF, UMI validation — matching standalone group command).
    group_filter_config: GroupFilterConfig,
    /// Stage to stop after (None = run through filter, the default).
    stop_after: Option<StopAfter>,
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
    /// Progress tracker for consensus reads (periodic logging).
    progress: ProgressTracker,
    /// Progress tracker for position groups (periodic logging).
    group_progress: ProgressTracker,
}

impl ProcessingSink {
    #[expect(clippy::too_many_arguments, reason = "sink needs all stage configs")]
    pub(super) fn new(
        writer: RawBamWriter,
        caller: Box<dyn ConsensusCaller>,
        assigner: Box<dyn UmiAssigner>,
        umi_tag: [u8; 2],
        filter_config: FilterConfig,
        group_filter_config: GroupFilterConfig,
        stop_after: Option<StopAfter>,
        cancel: Arc<AtomicBool>,
    ) -> Self {
        Self {
            writer: Some(writer),
            caller,
            assigner,
            umi_tag,
            filter_config,
            group_filter_config,
            stop_after,
            cancel,
            current_batch: Vec::new(),
            current_pos_key: None,
            total_consensus_reads: 0,
            position_groups_processed: 0,
            total_input_records: 0,
            total_mi_groups: 0,
            progress: ProgressTracker::new("Consensus reads written").with_interval(100_000),
            group_progress: ProgressTracker::new("Position groups processed").with_interval(10_000),
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

        // --stop-after group: write sorted records with MI tags and skip consensus
        if self.stop_after == Some(StopAfter::Group) {
            use fgumi_lib::sort::bam_fields;
            let mi_tag_bytes = [b'M', b'I']; // always write MI tag for grouped output
            let mut mi_buf = String::with_capacity(16);
            let mut scratch = Vec::with_capacity(512);
            let writer = self.writer.as_mut().expect("write called after finish");
            for template in templates {
                let mi = template.mi;
                let has_mi = mi.is_assigned();
                if has_mi {
                    mi.write_with_offset(0, &mut mi_buf);
                }
                for raw in template.into_raw_records().unwrap_or_default() {
                    if has_mi {
                        scratch.clear();
                        scratch.extend_from_slice(&raw);
                        bam_fields::update_string_tag(
                            &mut scratch,
                            &mi_tag_bytes,
                            mi_buf.as_bytes(),
                        );
                        writer
                            .write_raw_record(&scratch)
                            .context("Failed to write grouped record")?;
                    } else {
                        writer.write_raw_record(&raw).context("Failed to write grouped record")?;
                    }
                    self.total_consensus_reads += 1;
                }
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
                            false,
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
        if let Some(writer) = self.writer.take() {
            writer.finish().context("Failed to finish output BAM")?;
        }
        Ok(())
    }
}
