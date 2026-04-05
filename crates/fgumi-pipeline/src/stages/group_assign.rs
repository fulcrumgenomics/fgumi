//! `GroupAndAssign` pipeline stage: template build, UMI assign, and MI regroup.
//!
//! This stage takes a batch of raw BAM records sharing the same genomic position,
//! groups them by queryname into [`Template`]s, runs UMI assignment to assign
//! molecule IDs, then regroups records by molecule ID into [`MiGroup`]s ready
//! for consensus calling.

use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use anyhow::Result;
use fgumi_lib::grouper::{
    FilterMetrics, GroupFilterConfig, extract_umi_from_template, filter_template_raw,
    get_pair_orientation_raw, group_by_queryname,
};
use fgumi_lib::template::Template;
use fgumi_metrics::downsampling::compute_hash_fraction;
use fgumi_metrics::inline_collector::InlineCollector;
use fgumi_metrics::template_info::{ReadInfoKey, TemplateInfo};
use fgumi_raw_bam::update_string_tag;
use fgumi_umi::{Umi, UmiAssigner};
use noodles::sam::Header;

use crate::stage::PipelineStage;

/// A batch of raw BAM records sharing the same genomic position.
///
/// Records within the batch are sorted by template-coordinate key, which includes
/// a name hash. Records with the same queryname are adjacent but name-hash
/// collisions are possible, so actual queryname comparison is required.
pub struct PositionGroupBatch {
    /// Raw BAM record bytes (without `block_size` prefix).
    pub records: Vec<Vec<u8>>,
    /// Position key `(primary, secondary)` from the template-coordinate sort key.
    ///
    /// Note: `cb_hash` (cellular barcode hash) is intentionally excluded from this key.
    /// `cb_hash` is used only for batch boundary detection in [`PositionBatcher`] — records
    /// with different `cb_hash` values are split into separate batches upstream, so within
    /// a single `PositionGroupBatch` all records share the same `cb_hash`.
    pub position_key: (u64, u64),
}

/// A group of raw BAM records sharing the same molecule ID, ready for consensus.
pub struct MiGroup {
    /// Raw BAM record bytes for all templates assigned to this molecule.
    pub records: Vec<Vec<u8>>,
    /// The assigned molecule ID.
    pub mi: u64,
    /// Precomputed sum of record byte lengths, avoiding redundant recomputation.
    pub byte_size: usize,
}

/// Pipeline stage that groups raw BAM records by queryname, assigns UMIs,
/// and regroups by molecule ID.
pub struct GroupAndAssignStage {
    assigner: Arc<dyn UmiAssigner>,
    umi_tag: [u8; 2],
    /// The two-byte BAM tag for the assigned molecule ID (e.g., `b"MI"`).
    assign_tag: [u8; 2],
    /// Template filter config (MAPQ, unmapped, non-PF, UMI validation).
    group_filter_config: GroupFilterConfig,
    /// Global base offset for molecule IDs, incremented atomically across
    /// position groups so that MI values are globally unique.
    next_mi: AtomicU64,
    /// Optional shared inline metrics collector (protected by mutex for thread safety).
    metrics_collector: Option<Arc<Mutex<InlineCollector>>>,
    /// BAM header for `ref_name` lookups during metrics collection.
    header: Option<Arc<Header>>,
}

impl GroupAndAssignStage {
    /// Create a new `GroupAndAssignStage`.
    ///
    /// # Arguments
    ///
    /// * `assigner` - The UMI assignment strategy to use.
    /// * `umi_tag` - The two-byte BAM tag containing the UMI (e.g., `b"RX"`).
    /// * `assign_tag` - The two-byte BAM tag for the assigned molecule ID (e.g., `b"MI"`).
    /// * `group_filter_config` - Template filter config (MAPQ, unmapped, etc.).
    pub fn new(
        assigner: Arc<dyn UmiAssigner>,
        umi_tag: [u8; 2],
        assign_tag: [u8; 2],
        group_filter_config: GroupFilterConfig,
    ) -> Self {
        Self {
            assigner,
            umi_tag,
            assign_tag,
            group_filter_config,
            next_mi: AtomicU64::new(0),
            metrics_collector: None,
            header: None,
        }
    }

    /// Enables inline metrics collection with a shared collector and BAM header.
    ///
    /// The collector is shared across all worker threads (protected by a mutex).
    /// The header is used for reference name lookups during template-to-metrics
    /// conversion.
    #[must_use]
    pub fn with_metrics(
        mut self,
        collector: Arc<Mutex<InlineCollector>>,
        header: Arc<Header>,
    ) -> Self {
        self.metrics_collector = Some(collector);
        self.header = Some(header);
        self
    }

    /// Extract the UMI string from a template, delegating to the shared helper.
    fn extract_umi(&self, template: &Template) -> Result<Umi> {
        extract_umi_from_template(template, self.umi_tag)
    }

    /// Converts templates with assigned MIs into `TemplateInfo` for metrics collection.
    #[expect(
        clippy::cast_possible_wrap,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "genomic positions never exceed i32::MAX; ref_id is validated non-negative"
    )]
    fn build_template_infos(&self, templates: &[Template]) -> Vec<TemplateInfo> {
        use fgumi_raw_bam::cigar::{alignment_end_from_raw, alignment_start_from_raw};
        use fgumi_raw_bam::fields::ref_id;

        let mut mi_buf = String::with_capacity(16);
        templates
            .iter()
            .filter(|t| t.mi.id().is_some())
            .map(|t| {
                t.mi.write_with_offset(0, &mut mi_buf);
                let mi = mi_buf.clone();
                let rx = extract_umi_from_template(t, self.umi_tag).unwrap_or_default();
                let read_name = String::from_utf8_lossy(&t.name);
                let hash_fraction = compute_hash_fraction(&read_name);

                // Extract position info from raw R1 and R2 bytes
                let r1 = t.raw_r1();
                let r2 = t.raw_r2();

                let ref_name = [r1, r2].into_iter().flatten().find_map(|raw| {
                    let rid = ref_id(raw);
                    if rid < 0 {
                        return None;
                    }
                    self.header.as_ref().and_then(|h| {
                        h.reference_sequences()
                            .get_index(rid as usize)
                            .map(|(name, _)| name.to_string())
                    })
                });

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

                // Build ReadInfoKey from unclipped 5' positions and strand
                let read_info_key = Self::build_read_info_key(r1, r2);

                TemplateInfo {
                    mi,
                    rx,
                    ref_name,
                    position,
                    end_position,
                    hash_fraction,
                    read_info_key,
                }
            })
            .collect()
    }

    /// Builds a [`ReadInfoKey`] from raw R1 and R2 BAM records.
    ///
    /// Uses unclipped 5' positions, reference indices, and strand from both mates.
    /// Fields are ordered so the earlier-mapping read comes first (matching fgbio).
    #[expect(clippy::cast_sign_loss, reason = "ref_id validated non-negative before cast")]
    fn build_read_info_key(r1: Option<&[u8]>, r2: Option<&[u8]>) -> ReadInfoKey {
        use fgumi_raw_bam::cigar::unclipped_5prime_from_raw_bam;
        use fgumi_raw_bam::fields::{flags as raw_flags, ref_id};

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
            (Some((ref1, s1, strand1)), None) => ReadInfoKey {
                ref_index1: Some(ref1),
                start1: s1,
                strand1,
                ..ReadInfoKey::default()
            },
            (None, Some((ref2, s2, strand2))) => ReadInfoKey {
                ref_index1: Some(ref2),
                start1: s2,
                strand1: strand2,
                ..ReadInfoKey::default()
            },
            (None, None) => ReadInfoKey::default(),
        }
    }
}

impl PipelineStage for GroupAndAssignStage {
    type Input = PositionGroupBatch;
    type Output = Vec<MiGroup>;

    fn process(&self, input: Self::Input) -> Result<Self::Output> {
        // 1. Group records by queryname (consecutive runs with same name).
        let record_groups = group_by_queryname(input.records);

        if record_groups.is_empty() {
            return Ok(Vec::new());
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
            return Ok(Vec::new());
        }

        // 3. Extract UMIs and run UMI assignment, optionally splitting by pair orientation.
        //    Templates with empty or missing UMIs are skipped (retain MoleculeId::None)
        //    to match the standalone group command's filtering behavior.
        if self.assigner.split_templates_by_pair_orientation() {
            let mut subgroups: HashMap<(bool, bool), Vec<usize>> = HashMap::new();
            for (idx, template) in templates.iter().enumerate() {
                let orientation = get_pair_orientation_raw(template);
                subgroups.entry(orientation).or_default().push(idx);
            }

            for indices in subgroups.values() {
                // Filter to templates with non-empty UMIs
                let (valid_indices, umis): (Vec<usize>, Vec<Umi>) = indices
                    .iter()
                    .filter_map(|&i| {
                        let umi = self.extract_umi(&templates[i]).ok()?;
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
            // Filter to templates with non-empty UMIs
            let (valid_indices, umis): (Vec<usize>, Vec<Umi>) = templates
                .iter()
                .enumerate()
                .filter_map(|(i, t)| {
                    let umi = self.extract_umi(t).ok()?;
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

        // 4b. Collect inline metrics if enabled.
        if let Some(ref collector_arc) = self.metrics_collector {
            let template_infos = self.build_template_infos(&templates);
            match collector_arc.lock() {
                Ok(mut collector) => collector.record_coordinate_group(&template_infos),
                Err(e) => log::warn!("Skipping metrics for position group (mutex poisoned): {e}"),
            }
        }

        // 5. Compute the max local MI id so we can reserve a contiguous range of
        //    global IDs for this position group.
        let max_local_id = templates.iter().filter_map(|t| t.mi.id()).max().map_or(0, |id| id + 1);
        let base_mi = self.next_mi.fetch_add(max_local_id, Ordering::Relaxed);

        // 6. Inject MI tags into raw records and regroup by molecule ID.
        let mut mi_map: HashMap<u64, Vec<Vec<u8>>> = HashMap::new();
        let mut mi_buf = String::new();

        for template in templates {
            let mi = template.mi;
            if !mi.is_assigned() {
                // Templates with MoleculeId::None are dropped (unassigned).
                continue;
            }
            let mi_value = mi.write_with_offset(base_mi, &mut mi_buf);
            let global_id = base_mi + mi.id().unwrap_or(0);
            let mut raw_records = template.into_raw_records().unwrap_or_default();
            for record in &mut raw_records {
                update_string_tag(record, &self.assign_tag, mi_value);
            }
            mi_map.entry(global_id).or_default().extend(raw_records);
        }

        let mi_groups: Vec<MiGroup> = mi_map
            .into_iter()
            .map(|(mi, records)| {
                let byte_size = records.iter().map(Vec::len).sum();
                MiGroup { records, mi, byte_size }
            })
            .collect();

        Ok(mi_groups)
    }

    fn output_memory_estimate(&self, output: &Self::Output) -> usize {
        output.iter().map(|g| g.byte_size).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::fields::read_name;
    use fgumi_umi::IdentityUmiAssigner;

    /// Build a minimal raw BAM record with the specified name, flags, and UMI.
    ///
    /// Creates a mapped paired-end record with a 4-base sequence and a single
    /// `RX:Z:<umi>` auxiliary tag.
    #[allow(clippy::cast_possible_truncation)]
    fn make_raw_record(name: &[u8], flag: u16, umi: &[u8]) -> Vec<u8> {
        let seq_len: usize = 4;
        let l_read_name = (name.len() + 1) as u8; // +1 for NUL
        // Single CIGAR op: 4M
        let cigar_ops: &[u32] = &[(seq_len as u32) << 4];
        let n_cigar_op = cigar_ops.len() as u16;
        let seq_bytes = seq_len.div_ceil(2);
        let total = 32 + l_read_name as usize + cigar_ops.len() * 4 + seq_bytes + seq_len;
        let mut buf = vec![0u8; total];

        // ref_id = 0
        buf[0..4].copy_from_slice(&0i32.to_le_bytes());
        // pos = 100
        buf[4..8].copy_from_slice(&100i32.to_le_bytes());
        // l_read_name
        buf[8] = l_read_name;
        // mapq = 60
        buf[9] = 60;
        // n_cigar_op
        buf[12..14].copy_from_slice(&n_cigar_op.to_le_bytes());
        // flags
        buf[14..16].copy_from_slice(&flag.to_le_bytes());
        // l_seq
        buf[16..20].copy_from_slice(&(seq_len as u32).to_le_bytes());
        // mate_ref_id = 0
        buf[20..24].copy_from_slice(&0i32.to_le_bytes());
        // mate_pos = 200
        buf[24..28].copy_from_slice(&200i32.to_le_bytes());
        // tlen = 150
        buf[28..32].copy_from_slice(&150i32.to_le_bytes());

        // Read name (NUL-terminated)
        let name_start = 32;
        buf[name_start..name_start + name.len()].copy_from_slice(name);
        buf[name_start + name.len()] = 0;

        // CIGAR
        let cigar_start = name_start + l_read_name as usize;
        for (i, &op) in cigar_ops.iter().enumerate() {
            let off = cigar_start + i * 4;
            buf[off..off + 4].copy_from_slice(&op.to_le_bytes());
        }

        // Append RX:Z:<umi>\0 tag
        if !umi.is_empty() {
            buf.extend_from_slice(b"RXZ");
            buf.extend_from_slice(umi);
            buf.push(0);
        }

        buf
    }

    /// SAM flags for paired, first segment (R1)
    const PAIRED_R1: u16 = 0x1 | 0x40; // paired + first_segment
    /// SAM flags for paired, second segment (R2)
    const PAIRED_R2: u16 = 0x1 | 0x80; // paired + last_segment

    #[test]
    fn test_group_by_queryname_groups_consecutive_records() {
        let r1 = make_raw_record(b"read1", PAIRED_R1, b"ACGT");
        let r2 = make_raw_record(b"read1", PAIRED_R2, b"ACGT");
        let r3 = make_raw_record(b"read2", PAIRED_R1, b"TTTT");
        let r4 = make_raw_record(b"read2", PAIRED_R2, b"TTTT");

        let groups = group_by_queryname(vec![r1, r2, r3, r4]);

        assert_eq!(groups.len(), 2);
        assert_eq!(read_name(&groups[0][0]), b"read1");
        assert_eq!(groups[0].len(), 2);
        assert_eq!(read_name(&groups[1][0]), b"read2");
        assert_eq!(groups[1].len(), 2);
    }

    #[test]
    fn test_group_by_queryname_empty_input() {
        let groups = group_by_queryname(vec![]);
        assert!(groups.is_empty());
    }

    #[test]
    fn test_group_by_queryname_single_record() {
        let r = make_raw_record(b"read1", PAIRED_R1, b"ACGT");
        let groups = group_by_queryname(vec![r]);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 1);
    }

    #[test]
    fn test_process_two_templates_same_umi() {
        // Two templates with the same UMI should get the same molecule ID
        let assigner = Arc::new(IdentityUmiAssigner::new());
        let stage = GroupAndAssignStage::new(
            assigner,
            *b"RX",
            *b"MI",
            GroupFilterConfig::with_defaults(*b"RX"),
        );

        let r1a = make_raw_record(b"readA", PAIRED_R1, b"ACGT");
        let r2a = make_raw_record(b"readA", PAIRED_R2, b"ACGT");
        let r1b = make_raw_record(b"readB", PAIRED_R1, b"ACGT");
        let r2b = make_raw_record(b"readB", PAIRED_R2, b"ACGT");

        let batch = PositionGroupBatch { records: vec![r1a, r2a, r1b, r2b], position_key: (0, 0) };

        let result = stage.process(batch).unwrap();
        // IdentityUmiAssigner groups identical UMIs together
        assert_eq!(result.len(), 1, "Same UMI should yield one MI group");
        assert_eq!(result[0].records.len(), 4, "All 4 records in single group");
    }

    #[test]
    fn test_process_two_templates_different_umis() {
        let assigner = Arc::new(IdentityUmiAssigner::new());
        let stage = GroupAndAssignStage::new(
            assigner,
            *b"RX",
            *b"MI",
            GroupFilterConfig::with_defaults(*b"RX"),
        );

        let r1a = make_raw_record(b"readA", PAIRED_R1, b"AAAA");
        let r2a = make_raw_record(b"readA", PAIRED_R2, b"AAAA");
        let r1b = make_raw_record(b"readB", PAIRED_R1, b"CCCC");
        let r2b = make_raw_record(b"readB", PAIRED_R2, b"CCCC");

        let batch = PositionGroupBatch { records: vec![r1a, r2a, r1b, r2b], position_key: (0, 0) };

        let result = stage.process(batch).unwrap();
        assert_eq!(result.len(), 2, "Different UMIs should yield two MI groups");
        for group in &result {
            assert_eq!(group.records.len(), 2, "Each group should have 2 records (R1+R2)");
        }
    }

    #[test]
    fn test_process_empty_batch() {
        let assigner = Arc::new(IdentityUmiAssigner::new());
        let stage = GroupAndAssignStage::new(
            assigner,
            *b"RX",
            *b"MI",
            GroupFilterConfig::with_defaults(*b"RX"),
        );

        let batch = PositionGroupBatch { records: vec![], position_key: (0, 0) };

        let result = stage.process(batch).unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_output_memory_estimate() {
        let assigner = Arc::new(IdentityUmiAssigner::new());
        let stage = GroupAndAssignStage::new(
            assigner,
            *b"RX",
            *b"MI",
            GroupFilterConfig::with_defaults(*b"RX"),
        );

        let groups = vec![
            MiGroup { records: vec![vec![0u8; 100], vec![0u8; 50]], mi: 0, byte_size: 150 },
            MiGroup { records: vec![vec![0u8; 200]], mi: 1, byte_size: 200 },
        ];

        assert_eq!(stage.output_memory_estimate(&groups), 350);
    }
}
