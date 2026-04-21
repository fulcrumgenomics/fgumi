//! `GroupAssignStage`: template build, UMI assign, and MI regroup.
//!
//! Parallel pool stage. Groups records of one position bucket by
//! queryname into [`Template`]s, runs UMI assignment to stamp molecule
//! IDs via a shared [`UmiAssigner`], then regroups records by molecule
//! ID into [`MiGroup`]s ready for consensus calling.
//!
//! ## Input
//!
//! [`PositionGroupBatch`] — concat-byte BAM records sharing a
//! template-coordinate position key, as emitted by
//! [`crate::runall::engine::stages::position_batch::PositionBatchStage`].
//!
//! ## Output
//!
//! [`MiGroupBatch`] — one [`MiGroup`] per assigned molecule, each
//! holding concat-byte BAM records. `ordinal` and `position_key` are
//! propagated from the input batch unchanged.
//!
//! ## Ordering guarantees
//!
//! MI groups are emitted in **ascending MI order** (via a
//! [`BTreeMap`]), and within each group records are sorted by `(MI,
//! QNAME)` — matching the standalone `fgumi group` output. Input
//! records are grouped by consecutive same-QNAME runs (no hashing), so
//! upstream must present reads of a template contiguously (as
//! `PositionBatchStage` does).
//!
//! ## Memory model
//!
//! 1:1 per input batch. Intermediate per-record `Vec<u8>` copies are
//! made from the concat-byte input for the mutation-heavy assigner
//! path; see D6 in the type-unification design spec for why a zero-
//! copy refactor was deferred.
//!
//! ## Determinism
//!
//! Deterministic per batch. The shared [`UmiAssigner`] behind an
//! `Arc<dyn UmiAssigner>` uses an internal atomic counter to mint
//! globally-unique MI values; the counter increments in the order
//! pool workers call `assign()`, so cross-run MI *values* can differ
//! — but the grouping itself (which reads cluster into which MI) is
//! stable. Metrics collection lives in a [`PerThreadAccumulator`], so
//! each worker merges into its own slot without cross-worker contention.

use std::collections::BTreeMap;
use std::sync::Arc;

use anyhow::Result;
use fgumi_metrics::downsampling::compute_hash_fraction;
use fgumi_metrics::inline_collector::InlineCollector;
use fgumi_metrics::template_info::{ReadInfoKey, TemplateInfo};
use fgumi_raw_bam::update_string_tag;
use fgumi_umi::UmiAssigner;
use noodles::sam::Header;

use crate::grouper::{
    FilterMetrics, GroupFilterConfig, assign_umi_groups_impl, extract_umi_from_template,
    filter_template_raw, group_by_queryname,
};
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::runall::engine::grouping_types::{
    MiGroup, MiGroupBatch, PositionGroupBatch, iter_length_prefixed,
};
use crate::runall::engine::stage::{Parallelism, Stage};
use crate::template::Template;

/// Pipeline group-assign stage.
///
/// `Clone` via shared state: `assigner`, `next_mi`, and `metrics_collector`
/// all live behind `Arc` so cloned instances share them. The work-stealing
/// pool instantiates a per-worker stage while preserving global MI
/// uniqueness — all workers increment the same `AtomicU64` counter through
/// their shared `Arc`.
#[derive(Clone)]
pub struct GroupAssignStage {
    assigner: Arc<dyn UmiAssigner>,
    umi_tag: [u8; 2],
    /// Two-byte BAM tag for the assigned molecule ID (e.g., `b"MI"`).
    assign_tag: [u8; 2],
    group_filter_config: GroupFilterConfig,
    metrics_collector: Option<Arc<PerThreadAccumulator<InlineCollector>>>,
    /// Optional header for `ref_name` lookups during metrics collection.
    header: Option<Arc<Header>>,
}

impl GroupAssignStage {
    /// Create a new `GroupAssignStage`.
    #[must_use]
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
            metrics_collector: None,
            header: None,
        }
    }

    /// Enable inline metrics collection.
    #[must_use]
    pub fn with_metrics(
        mut self,
        collector: Arc<PerThreadAccumulator<InlineCollector>>,
        header: Arc<Header>,
    ) -> Self {
        self.metrics_collector = Some(collector);
        self.header = Some(header);
        self
    }

    /// Convert assigned templates into [`TemplateInfo`] for metrics.
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

                let r1 = t.r1();
                let r2 = t.r2();

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

                let read_info_key =
                    Self::build_read_info_key(r1.map(AsRef::as_ref), r2.map(AsRef::as_ref));

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

    #[expect(clippy::cast_sign_loss, reason = "ref_id validated non-negative before cast")]
    fn build_read_info_key(r1: Option<&[u8]>, r2: Option<&[u8]>) -> ReadInfoKey {
        use fgumi_raw_bam::RawRecordView;
        use fgumi_raw_bam::cigar::unclipped_5prime_from_raw_bam;
        use fgumi_raw_bam::fields::ref_id;

        fn mate_info(raw: &[u8]) -> Option<(usize, i32, bool)> {
            let rid = ref_id(raw);
            if rid < 0 {
                return None;
            }
            let flg = RawRecordView::new(raw).flags();
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

    /// Core processing logic: `PositionGroupBatch` (concat bytes) ->
    /// `MiGroupBatch` (one per input, with MIs sorted ascending).
    ///
    /// Parses the input concat bytes into per-record `Vec<u8>` buffers for the
    /// mutation-heavy grouping/assignment body; emits each inner `MiGroup`
    /// with concat-byte `data` and stable MI ordering via `BTreeMap`.
    fn group_assign(&self, input: &PositionGroupBatch) -> Result<MiGroupBatch> {
        // 1. Parse concat bytes into per-record buffers.
        let records: Vec<Vec<u8>> = iter_length_prefixed(&input.data)
            .map(|r| r.map(<[u8]>::to_vec))
            .collect::<Result<Vec<_>>>()?;

        let empty_batch = || MiGroupBatch {
            groups: Vec::new(),
            ordinal: input.ordinal,
            position_key: input.position_key,
        };

        // 2. Group records by queryname (consecutive runs with same name).
        let record_groups = group_by_queryname(records);

        if record_groups.is_empty() {
            return Ok(empty_batch());
        }

        // 3. Build Templates from each name group.
        let mut templates: Vec<Template> = Vec::with_capacity(record_groups.len());
        for group in record_groups {
            let raw_records: Vec<fgumi_raw_bam::RawRecord> =
                group.into_iter().map(fgumi_raw_bam::RawRecord::from).collect();
            templates.push(Template::from_records(raw_records)?);
        }

        // 3b. Filter templates (same criteria as standalone group command).
        let mut filter_metrics = FilterMetrics::default();
        templates
            .retain(|t| filter_template_raw(t, &self.group_filter_config, &mut filter_metrics));

        if templates.is_empty() {
            return Ok(empty_batch());
        }

        // 4. UMI assignment via the shared library orchestration (pair-
        //    orientation split, UMI canonicalization, truncation, and
        //    `assigner.assign()`). The assigner returns globally-unique IDs
        //    directly via its internal atomic counter, so no `base_mi` offset
        //    is applied here.
        assign_umi_groups_impl(
            &mut templates,
            self.assigner.as_ref(),
            self.umi_tag,
            self.group_filter_config.min_umi_length,
            self.group_filter_config.no_umi,
        )?;

        // 5. Inline metrics (if enabled). Each worker merges into its own
        // slot; slots are folded at pipeline end by the planner.
        if let Some(ref collector_arc) = self.metrics_collector {
            let template_infos = self.build_template_infos(&templates);
            collector_arc.with_slot(|collector| collector.record_coordinate_group(&template_infos));
        }

        // 6. Sort templates by (MI, QNAME) so within-MI-group record ordering
        //    matches the standalone `fgumi group` output — commands/group.rs
        //    does the same `.sort_by(|a, b| a.mi.to_vec_index().cmp(&b.mi.to_vec_index())
        //    .then_with(|| a.name.cmp(&b.name)))` before serializing. Matters
        //    downstream: consensus calling averages over reads in input order,
        //    so a different within-group order produces slightly different
        //    consensus qualities even from the same set of records.
        templates.sort_by(|a, b| {
            let a_idx = a.mi.to_vec_index();
            let b_idx = b.mi.to_vec_index();
            a_idx.cmp(&b_idx).then_with(|| a.name.cmp(&b.name))
        });

        // 7. Inject MI tags into raw records and regroup by molecule ID. Use
        //    `BTreeMap` so downstream consumers see MI groups in stable
        //    ascending MI order — required by `Coalesce<MiGroupBatch>` and the
        //    byte-compare gate. MI values come straight from
        //    `MoleculeId::write_with_offset(0, _)` since assigner IDs are
        //    already globally unique.
        let mut mi_map: BTreeMap<u64, Vec<Vec<u8>>> = BTreeMap::new();
        let mut mi_buf = String::new();

        for template in templates {
            let mi = template.mi;
            if !mi.is_assigned() {
                continue;
            }
            let mi_value = mi.write_with_offset(0, &mut mi_buf);
            let global_id = mi.id().unwrap_or(0);
            let mut raw_records = template.into_records();
            for record in &mut raw_records {
                update_string_tag(record.as_mut_vec(), self.assign_tag, mi_value);
            }
            mi_map
                .entry(global_id)
                .or_default()
                .extend(raw_records.into_iter().map(fgumi_raw_bam::RawRecord::into_inner));
        }

        // 8. Serialize each MI group into concat-byte `MiGroup`.
        let mi_groups: Vec<MiGroup> = mi_map
            .into_iter()
            .map(|(mi, records)| {
                let (data, record_count) =
                    crate::runall::engine::grouping_types::serialize_records(&records)?;
                Ok::<_, anyhow::Error>(MiGroup { data, record_count, mi })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(MiGroupBatch {
            groups: mi_groups,
            ordinal: input.ordinal,
            position_key: input.position_key,
        })
    }
}

impl Stage for GroupAssignStage {
    type Input = PositionGroupBatch;
    type Output = MiGroupBatch;

    #[tracing::instrument(name = "group_assign", skip_all)]
    fn process(&mut self, input: Self::Input, out: &mut dyn FnMut(Self::Output)) -> Result<()> {
        let batch = self.group_assign(&input)?;
        out(batch);
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }

    fn output_memory_estimate(&self, output: &Self::Output) -> usize {
        output.groups.iter().map(|g| g.data.len()).sum()
    }

    fn name(&self) -> &'static str {
        "GroupAssign"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::grouping_types::{iter_length_prefixed, serialize_records};
    use fgumi_raw_bam::fields::read_name;
    use fgumi_umi::IdentityUmiAssigner;

    /// Minimal raw BAM record with the given name, flag, and UMI.
    #[allow(clippy::cast_possible_truncation)]
    fn make_raw_record(name: &[u8], flag: u16, umi: &[u8]) -> Vec<u8> {
        let seq_len: usize = 4;
        let l_read_name = (name.len() + 1) as u8; // +1 for NUL
        let cigar_ops: &[u32] = &[(seq_len as u32) << 4];
        let n_cigar_op = cigar_ops.len() as u16;
        let seq_bytes = seq_len.div_ceil(2);
        let total = 32 + l_read_name as usize + cigar_ops.len() * 4 + seq_bytes + seq_len;
        let mut buf = vec![0u8; total];

        buf[0..4].copy_from_slice(&0i32.to_le_bytes()); // ref_id
        buf[4..8].copy_from_slice(&100i32.to_le_bytes()); // pos
        buf[8] = l_read_name;
        buf[9] = 60; // mapq
        buf[12..14].copy_from_slice(&n_cigar_op.to_le_bytes());
        buf[14..16].copy_from_slice(&flag.to_le_bytes());
        buf[16..20].copy_from_slice(&(seq_len as u32).to_le_bytes());
        buf[20..24].copy_from_slice(&0i32.to_le_bytes());
        buf[24..28].copy_from_slice(&200i32.to_le_bytes());
        buf[28..32].copy_from_slice(&150i32.to_le_bytes());

        let name_start = 32;
        buf[name_start..name_start + name.len()].copy_from_slice(name);
        buf[name_start + name.len()] = 0;

        let cigar_start = name_start + l_read_name as usize;
        for (i, &op) in cigar_ops.iter().enumerate() {
            let off = cigar_start + i * 4;
            buf[off..off + 4].copy_from_slice(&op.to_le_bytes());
        }

        if !umi.is_empty() {
            buf.extend_from_slice(b"RXZ");
            buf.extend_from_slice(umi);
            buf.push(0);
        }

        buf
    }

    const PAIRED_R1: u16 = 0x1 | 0x40;
    const PAIRED_R2: u16 = 0x1 | 0x80;

    fn mk_batch(records: &[Vec<u8>], position_key: (u64, u64)) -> PositionGroupBatch {
        let (data, record_count) = serialize_records(records).unwrap();
        PositionGroupBatch { data, record_count, ordinal: 0, position_key }
    }

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
        let assigner = Arc::new(IdentityUmiAssigner::new());
        let stage = GroupAssignStage::new(
            assigner,
            *b"RX",
            *b"MI",
            GroupFilterConfig::with_defaults(*b"RX"),
        );

        let r1a = make_raw_record(b"readA", PAIRED_R1, b"ACGT");
        let r2a = make_raw_record(b"readA", PAIRED_R2, b"ACGT");
        let r1b = make_raw_record(b"readB", PAIRED_R1, b"ACGT");
        let r2b = make_raw_record(b"readB", PAIRED_R2, b"ACGT");

        let batch = mk_batch(&[r1a, r2a, r1b, r2b], (0, 0));

        let result = stage.group_assign(&batch).unwrap();
        assert_eq!(result.groups.len(), 1, "Same UMI should yield one MI group");
        assert_eq!(result.groups[0].record_count, 4, "All 4 records in single group");
        // Verify all 4 records are present in the concat bytes
        let records: Vec<&[u8]> =
            iter_length_prefixed(&result.groups[0].data).map(Result::unwrap).collect();
        assert_eq!(records.len(), 4);
    }

    #[test]
    fn test_process_two_templates_different_umis() {
        let assigner = Arc::new(IdentityUmiAssigner::new());
        let stage = GroupAssignStage::new(
            assigner,
            *b"RX",
            *b"MI",
            GroupFilterConfig::with_defaults(*b"RX"),
        );

        let r1a = make_raw_record(b"readA", PAIRED_R1, b"AAAA");
        let r2a = make_raw_record(b"readA", PAIRED_R2, b"AAAA");
        let r1b = make_raw_record(b"readB", PAIRED_R1, b"CCCC");
        let r2b = make_raw_record(b"readB", PAIRED_R2, b"CCCC");

        let batch = mk_batch(&[r1a, r2a, r1b, r2b], (0, 0));

        let result = stage.group_assign(&batch).unwrap();
        assert_eq!(result.groups.len(), 2, "Different UMIs should yield two MI groups");
        for group in &result.groups {
            assert_eq!(group.record_count, 2, "Each group should have 2 records (R1+R2)");
        }
    }

    #[test]
    fn test_process_empty_batch() {
        let assigner = Arc::new(IdentityUmiAssigner::new());
        let stage = GroupAssignStage::new(
            assigner,
            *b"RX",
            *b"MI",
            GroupFilterConfig::with_defaults(*b"RX"),
        );

        let batch = mk_batch(&[], (0, 0));

        let result = stage.group_assign(&batch).unwrap();
        assert!(result.groups.is_empty());
    }

    #[test]
    fn test_output_memory_estimate() {
        let assigner = Arc::new(IdentityUmiAssigner::new());
        let stage = GroupAssignStage::new(
            assigner,
            *b"RX",
            *b"MI",
            GroupFilterConfig::with_defaults(*b"RX"),
        );

        // Build two MiGroups with known byte counts via serialize_records.
        let g0_records: Vec<Vec<u8>> = vec![vec![0u8; 100], vec![0u8; 50]];
        let (g0_data, g0_count) = serialize_records(&g0_records).unwrap();
        let g1_records: Vec<Vec<u8>> = vec![vec![0u8; 200]];
        let (g1_data, g1_count) = serialize_records(&g1_records).unwrap();

        let mi_batch = MiGroupBatch {
            groups: vec![
                MiGroup { data: g0_data, record_count: g0_count, mi: 0 },
                MiGroup { data: g1_data, record_count: g1_count, mi: 1 },
            ],
            ordinal: 0,
            position_key: (0, 0),
        };

        // Expected: each record gets +4 bytes for the block_size prefix.
        // g0: 100 + 4 + 50 + 4 = 158; g1: 200 + 4 = 204; total = 362.
        assert_eq!(stage.output_memory_estimate(&mi_batch), 362);
        assert_eq!(stage.name(), "GroupAssign");
        assert_eq!(stage.parallelism(), Parallelism::Parallel);
    }
}
