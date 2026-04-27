//! UMI correction stage.
//!
//! Parallel pool stage: one [`CorrectStage`] instance per pool worker
//! (constructed via `stage_factory` so each worker owns a fresh LRU
//! cache). Matches against the shared [`EncodedUmiSet`]; metrics are
//! accumulated into a shared mutex and dumped to TSV by the caller
//! after the pipeline finishes.
//!
//! ## Input
//!
//! [`SerializedBatch`] — length-prefixed BAM records with an `RX` tag
//! carrying the raw UMI. Typically produced by `ExtractStage` or
//! `BamBatchedSource`.
//!
//! ## Output
//!
//! [`SerializedBatch`] — kept records with the corrected UMI written
//! back into the `RX` tag (and original UMI optionally stored); any
//! rejected records are split out into `secondary.data`. The
//! `ordinal` field is preserved unchanged from input to output.
//!
//! ## Ordering guarantees
//!
//! Record order within the batch is preserved: `CorrectStage` walks
//! the input linearly and appends kept records to `primary` and
//! rejected records to `secondary`. Cross-batch ordering is the
//! upstream `ordinal` sequence, preserved verbatim.
//!
//! ## Memory model
//!
//! 1:1 with input plus small LRU overhead. Output buffers are
//! preallocated to `input.primary.data.len()`. Each worker caches up
//! to `cache_size` recent UMI → `UmiMatch` mappings; `cache_size == 0`
//! disables the cache entirely.
//!
//! ## Determinism
//!
//! Byte-identical across runs for a given input: correction is a pure
//! function of `(RX tag, EncodedUmiSet, thresholds)`. The LRU cache
//! only affects performance, not results — keys are uppercase UMI
//! segment bytes, matching the standalone command's cache key.

use std::num::NonZero;
use std::sync::Arc;

use ahash::AHashMap;
use anyhow::{Context, Result};
use lru::LruCache;

use crate::correct::{EncodedUmiSet, UmiMatch};
use crate::per_thread_accumulator::PerThreadAccumulator;
use crate::runall::engine::grouping_types::iter_length_prefixed;
use crate::runall::engine::output_types::SerializedBatch;
use crate::runall::engine::stage::{Parallelism, Stage};

/// Per-slot accumulator held inside a [`PerThreadAccumulator`]. Each worker
/// merges its batch-local totals into the slot it owns; the planner folds
/// slots into final metrics at pipeline end.
#[derive(Default)]
pub struct CorrectMetricsAccumulator {
    /// Per-UMI counts (key is the UMI string, or the unmatched sentinel).
    pub per_umi: AHashMap<String, UmiCounts>,
    /// Total templates (records) seen.
    pub total_records: u64,
    /// Records whose RX tag was missing.
    pub missing_umis: u64,
    /// Records whose RX tag was the wrong length.
    pub wrong_length: u64,
    /// Records whose RX tag was too far from any reference.
    pub mismatched: u64,
}

impl CorrectMetricsAccumulator {
    /// Fold `other` into `self`. Used both per-batch (local → slot) and at
    /// pipeline end (slot → global).
    pub fn merge(&mut self, other: Self) {
        self.total_records += other.total_records;
        self.missing_umis += other.missing_umis;
        self.wrong_length += other.wrong_length;
        self.mismatched += other.mismatched;
        for (umi, counts) in other.per_umi {
            let entry = self.per_umi.entry(umi).or_default();
            entry.total_matches += counts.total_matches;
            entry.perfect_matches += counts.perfect_matches;
            entry.one_mismatch_matches += counts.one_mismatch_matches;
            entry.two_mismatch_matches += counts.two_mismatch_matches;
            entry.other_matches += counts.other_matches;
        }
    }
}

/// Counts for a single (reference or unmatched) UMI bucket.
#[derive(Default, Clone)]
pub struct UmiCounts {
    pub total_matches: u64,
    pub perfect_matches: u64,
    pub one_mismatch_matches: u64,
    pub two_mismatch_matches: u64,
    pub other_matches: u64,
}

/// Parallel stage that corrects UMIs in-place on serialized BAM records.
pub struct CorrectStage {
    encoded_umi_set: Arc<EncodedUmiSet>,
    umi_length: usize,
    max_mismatches: usize,
    min_distance_diff: usize,
    revcomp: bool,
    dont_store_original_umis: bool,
    track_rejects: bool,
    unmatched_umi: String,
    rx_tag: [u8; 2],
    cache: Option<LruCache<Vec<u8>, UmiMatch>>,
    metrics: Arc<PerThreadAccumulator<CorrectMetricsAccumulator>>,
}

impl CorrectStage {
    /// Construct a new [`CorrectStage`]. `cache_size == 0` disables the cache.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        encoded_umi_set: Arc<EncodedUmiSet>,
        umi_length: usize,
        max_mismatches: usize,
        min_distance_diff: usize,
        revcomp: bool,
        dont_store_original_umis: bool,
        track_rejects: bool,
        cache_size: usize,
        metrics: Arc<PerThreadAccumulator<CorrectMetricsAccumulator>>,
    ) -> Self {
        let cache = NonZero::new(cache_size).map(LruCache::new);
        let unmatched_umi = "N".repeat(umi_length);
        Self {
            encoded_umi_set,
            umi_length,
            max_mismatches,
            min_distance_diff,
            revcomp,
            dont_store_original_umis,
            track_rejects,
            unmatched_umi,
            rx_tag: *crate::sam::SamTag::RX,
            cache,
            metrics,
        }
    }
}

impl Stage for CorrectStage {
    type Input = SerializedBatch;
    type Output = SerializedBatch;

    #[allow(clippy::too_many_lines)]
    #[tracing::instrument(name = "correct", skip_all)]
    fn process(&mut self, input: Self::Input, output: &mut dyn FnMut(Self::Output)) -> Result<()> {
        use crate::correct::{
            RejectionReason, apply_correction_to_raw, compute_template_correction,
        };
        use crate::sort::bam_fields;

        let SerializedBatch { primary, secondary: _, ordinal } = input;

        let mut local = CorrectMetricsAccumulator::default();
        let mut kept_data: Vec<u8> = Vec::with_capacity(primary.data.len());
        let mut kept_count: u64 = 0;
        let mut rejected_data: Vec<u8> = Vec::new();
        let mut rejected_count: u64 = 0;

        for record_slice in iter_length_prefixed(&primary.data) {
            let record_slice = record_slice.context("CorrectStage: framed-BAM parse error")?;
            let mut record: Vec<u8> = record_slice.to_vec();
            local.total_records += 1;

            let aux = bam_fields::aux_data_slice(&record);
            let umi_bytes = bam_fields::find_string_tag(aux, self.rx_tag);
            let Some(umi_slice) = umi_bytes else {
                local.missing_umis += 1;
                let entry = local.per_umi.entry(self.unmatched_umi.clone()).or_default();
                entry.total_matches += 1;
                if self.track_rejects {
                    let size = u32::try_from(record.len())
                        .context("CorrectStage: reject record too large for u32 block_size")?;
                    rejected_data.extend_from_slice(&size.to_le_bytes());
                    rejected_data.extend_from_slice(&record);
                    rejected_count += 1;
                }
                continue;
            };
            let umi_string = String::from_utf8_lossy(umi_slice).into_owned();

            let correction = compute_template_correction(
                &umi_string,
                self.umi_length,
                self.revcomp,
                self.max_mismatches,
                self.min_distance_diff,
                &self.encoded_umi_set,
                &mut self.cache,
            );

            if correction.matched {
                for m in &correction.matches {
                    if m.matched {
                        let entry = local.per_umi.entry(m.umi.clone()).or_default();
                        entry.total_matches += 1;
                        match m.mismatches {
                            0 => entry.perfect_matches += 1,
                            1 => entry.one_mismatch_matches += 1,
                            2 => entry.two_mismatch_matches += 1,
                            _ => entry.other_matches += 1,
                        }
                    }
                }

                apply_correction_to_raw(
                    &mut record,
                    &correction,
                    self.rx_tag,
                    self.dont_store_original_umis,
                );

                let new_size = u32::try_from(record.len())
                    .context("CorrectStage: record too large for u32 block_size")?;
                kept_data.extend_from_slice(&new_size.to_le_bytes());
                kept_data.extend_from_slice(&record);
                kept_count += 1;
            } else {
                match correction.rejection_reason {
                    RejectionReason::WrongLength => local.wrong_length += 1,
                    RejectionReason::Mismatched => local.mismatched += 1,
                    RejectionReason::None => {}
                }
                let entry = local.per_umi.entry(self.unmatched_umi.clone()).or_default();
                entry.total_matches += 1;
                if self.track_rejects {
                    let size = u32::try_from(record.len())
                        .context("CorrectStage: reject record too large for u32 block_size")?;
                    rejected_data.extend_from_slice(&size.to_le_bytes());
                    rejected_data.extend_from_slice(&record);
                    rejected_count += 1;
                }
            }
        }

        self.metrics.with_slot(|slot| slot.merge(local));

        let secondary = if self.track_rejects {
            Some(crate::runall::engine::output_types::RawBytes {
                data: rejected_data,
                record_count: rejected_count,
            })
        } else {
            None
        };
        output(SerializedBatch {
            primary: crate::runall::engine::output_types::RawBytes {
                data: kept_data,
                record_count: kept_count,
            },
            secondary,
            ordinal,
        });
        Ok(())
    }

    fn parallelism(&self) -> Parallelism {
        Parallelism::Parallel
    }

    fn output_memory_estimate(&self, out: &Self::Output) -> usize {
        out.primary.data.capacity() + out.secondary.as_ref().map_or(0, |s| s.data.capacity())
    }

    fn name(&self) -> &'static str {
        "CorrectStage"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::runall::engine::output_types::{RawBytes, SerializedBatch};

    /// Build a minimal BAM record with the given RX tag using the fgumi-raw-bam
    /// test builder. Unmapped (tid=-1, pos=-1), flag=4 (unmapped), empty CIGAR,
    /// zero sequence length, aux contains just a single Z-type RX tag.
    fn make_record_with_rx(rx: &[u8]) -> Vec<u8> {
        // Build aux bytes: "RX" tag + 'Z' type + rx_value + null terminator.
        let mut aux = Vec::with_capacity(3 + rx.len() + 1);
        aux.extend_from_slice(b"RX");
        aux.push(b'Z');
        aux.extend_from_slice(rx);
        aux.push(0);
        fgumi_raw_bam::testutil::make_bam_bytes(
            -1,
            -1,
            4, // unmapped
            b"r",
            &[],
            0,
            -1,
            -1,
            &aux,
        )
    }

    fn batch_from_records(records: &[Vec<u8>], ordinal: u64) -> SerializedBatch {
        let mut data = Vec::new();
        for r in records {
            let size = u32::try_from(r.len()).unwrap();
            data.extend_from_slice(&size.to_le_bytes());
            data.extend_from_slice(r);
        }
        SerializedBatch {
            primary: RawBytes { data, record_count: records.len() as u64 },
            secondary: None,
            ordinal,
        }
    }

    #[test]
    fn test_correct_stage_all_match_pass_through() {
        let metrics = PerThreadAccumulator::<CorrectMetricsAccumulator>::new(1);
        let mut stage =
            CorrectStage::new(tiny_umi_set(), 8, 1, 2, false, false, false, 128, metrics.clone());
        let r0 = make_record_with_rx(b"AAAAAAAA");
        let r1 = make_record_with_rx(b"CCCCCCCC");
        let input = batch_from_records(&[r0, r1], 3);

        let mut captured = None;
        stage.process(input, &mut |v| captured = Some(v)).unwrap();
        let out = captured.expect("stage must emit");

        assert_eq!(out.ordinal, 3);
        assert_eq!(out.primary.record_count, 2);
        assert!(out.secondary.is_none());

        let m = metrics.slots()[0].lock();
        assert_eq!(m.total_records, 2);
        assert_eq!(m.per_umi.get("AAAAAAAA").unwrap().perfect_matches, 1);
        assert_eq!(m.per_umi.get("CCCCCCCC").unwrap().perfect_matches, 1);
        assert_eq!(m.mismatched, 0);
    }

    #[test]
    fn test_correct_stage_drops_non_matching() {
        let metrics = PerThreadAccumulator::<CorrectMetricsAccumulator>::new(1);
        let mut stage =
            CorrectStage::new(tiny_umi_set(), 8, 1, 2, false, false, false, 128, metrics.clone());
        let r_match = make_record_with_rx(b"AAAAAAAA");
        let r_reject = make_record_with_rx(b"ACGTACGT"); // >1 mismatch from all refs
        let input = batch_from_records(&[r_match, r_reject], 0);

        let mut captured = None;
        stage.process(input, &mut |v| captured = Some(v)).unwrap();
        let out = captured.unwrap();
        assert_eq!(out.primary.record_count, 1);
        let m = metrics.slots()[0].lock();
        assert_eq!(m.total_records, 2);
        assert_eq!(m.mismatched, 1);
        assert_eq!(m.per_umi.get("NNNNNNNN").unwrap().total_matches, 1);
    }

    #[test]
    fn test_correct_stage_single_mismatch_counts_and_appends_ox() {
        let metrics = PerThreadAccumulator::<CorrectMetricsAccumulator>::new(1);
        let mut stage =
            CorrectStage::new(tiny_umi_set(), 8, 1, 2, false, false, false, 128, metrics.clone());
        let r = make_record_with_rx(b"AAAAAAAC");
        let input = batch_from_records(&[r], 0);

        let mut captured = None;
        stage.process(input, &mut |v| captured = Some(v)).unwrap();
        let out = captured.unwrap();
        assert_eq!(out.primary.record_count, 1);

        let block_size = u32::from_le_bytes(out.primary.data[0..4].try_into().unwrap()) as usize;
        let record = &out.primary.data[4..4 + block_size];
        let aux = crate::sort::bam_fields::aux_data_slice(record);
        let rx = crate::sort::bam_fields::find_string_tag(aux, b"RX").unwrap();
        let ox = crate::sort::bam_fields::find_string_tag(aux, b"OX").unwrap();
        assert_eq!(rx, b"AAAAAAAA");
        assert_eq!(ox, b"AAAAAAAC");

        let m = metrics.slots()[0].lock();
        assert_eq!(m.per_umi.get("AAAAAAAA").unwrap().one_mismatch_matches, 1);
    }

    #[test]
    fn test_correct_stage_dont_store_ox_when_disabled() {
        let metrics = PerThreadAccumulator::<CorrectMetricsAccumulator>::new(1);
        let mut stage =
            CorrectStage::new(tiny_umi_set(), 8, 1, 2, false, true, false, 128, metrics);
        let r = make_record_with_rx(b"AAAAAAAC");
        let input = batch_from_records(&[r], 0);

        let mut captured = None;
        stage.process(input, &mut |v| captured = Some(v)).unwrap();
        let out = captured.unwrap();
        let block_size = u32::from_le_bytes(out.primary.data[0..4].try_into().unwrap()) as usize;
        let record = &out.primary.data[4..4 + block_size];
        let aux = crate::sort::bam_fields::aux_data_slice(record);
        assert!(
            crate::sort::bam_fields::find_string_tag(aux, b"OX").is_none(),
            "OX must not be added when dont_store_original_umis=true"
        );
    }

    #[test]
    fn test_correct_stage_preserves_ordinal() {
        let metrics = PerThreadAccumulator::<CorrectMetricsAccumulator>::new(1);
        let mut stage = CorrectStage::new(tiny_umi_set(), 8, 1, 2, false, false, false, 0, metrics);
        for ord in [0u64, 1, 17, u64::MAX - 1] {
            let r = make_record_with_rx(b"AAAAAAAA");
            let input = batch_from_records(&[r], ord);
            let mut captured = None;
            stage.process(input, &mut |v| captured = Some(v)).unwrap();
            assert_eq!(captured.unwrap().ordinal, ord);
        }
    }

    fn tiny_umi_set() -> Arc<EncodedUmiSet> {
        let seqs: Vec<String> = ["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "TTTTTTTT"]
            .iter()
            .map(|s| (*s).to_string())
            .collect();
        Arc::new(EncodedUmiSet::new(&seqs))
    }

    #[test]
    fn test_correct_stage_new_with_cache() {
        let metrics = PerThreadAccumulator::<CorrectMetricsAccumulator>::new(1);
        let stage =
            CorrectStage::new(tiny_umi_set(), 8, 1, 2, false, false, true, 128, metrics.clone());
        assert_eq!(stage.umi_length, 8);
        assert_eq!(stage.unmatched_umi, "NNNNNNNN");
        assert!(stage.cache.is_some());
        assert!(stage.track_rejects);
        assert_eq!(stage.name(), "CorrectStage");
        assert_eq!(stage.parallelism(), Parallelism::Parallel);
    }

    #[test]
    fn test_correct_stage_new_no_cache() {
        let metrics = PerThreadAccumulator::<CorrectMetricsAccumulator>::new(1);
        let stage = CorrectStage::new(tiny_umi_set(), 8, 1, 2, false, false, false, 0, metrics);
        assert!(stage.cache.is_none());
        assert!(!stage.track_rejects);
    }

    #[test]
    fn test_correct_stage_factory_produces_independent_caches() {
        let metrics = PerThreadAccumulator::<CorrectMetricsAccumulator>::new(1);
        let umi_set = tiny_umi_set();
        let factory = {
            let umi_set = umi_set.clone();
            let metrics = metrics.clone();
            move || {
                CorrectStage::new(
                    umi_set.clone(),
                    8,
                    1,
                    2,
                    false,
                    false,
                    true,
                    128,
                    metrics.clone(),
                )
            }
        };
        let a = factory();
        let b = factory();
        let a_ptr = a.cache.as_ref().unwrap() as *const _;
        let b_ptr = b.cache.as_ref().unwrap() as *const _;
        assert_ne!(a_ptr, b_ptr);
        assert!(Arc::ptr_eq(&a.metrics, &b.metrics));
    }

    #[test]
    fn test_correct_stage_rejects_routed_to_secondary() {
        let metrics = PerThreadAccumulator::<CorrectMetricsAccumulator>::new(1);
        let mut stage = CorrectStage::new(
            tiny_umi_set(),
            8,
            1,
            2,
            false,
            false,
            true, // track_rejects
            0,
            metrics,
        );
        let r_ok = make_record_with_rx(b"AAAAAAAA");
        let r_bad = make_record_with_rx(b"ACGTACGT");
        let input = batch_from_records(&[r_ok, r_bad], 5);

        let mut captured = None;
        stage.process(input, &mut |v| captured = Some(v)).unwrap();
        let out = captured.unwrap();

        assert_eq!(out.ordinal, 5);
        assert_eq!(out.primary.record_count, 1);
        let sec = out.secondary.expect("rejects channel must be present");
        assert_eq!(sec.record_count, 1);
        assert!(!out.primary.data.is_empty());
        assert!(!sec.data.is_empty());
    }

    #[test]
    fn test_correct_stage_track_rejects_off_emits_none() {
        let metrics = PerThreadAccumulator::<CorrectMetricsAccumulator>::new(1);
        let mut stage = CorrectStage::new(
            tiny_umi_set(),
            8,
            1,
            2,
            false,
            false,
            false, // track_rejects off
            0,
            metrics,
        );
        let r_bad = make_record_with_rx(b"ACGTACGT");
        let input = batch_from_records(&[r_bad], 0);
        let mut captured = None;
        stage.process(input, &mut |v| captured = Some(v)).unwrap();
        assert!(captured.unwrap().secondary.is_none());
    }
}
