//! Consensus-read filtering pipeline stage.
//!
//! This stage consumes the serialized output of [`ConsensusStage`] — a buffer of
//! concatenated BAM records with 4-byte little-endian `block_size` prefixes — and
//! applies quality-based masking and per-read filtering to produce a flat `Vec<u8>` of
//! framed BAM record bytes (with `block_size` prefixes) ready for BGZF compression.
//!
//! For each record the stage:
//! 1. Optionally masks low-quality or low-depth bases in place.
//! 2. Checks per-read quality metrics against `FilterConfig` thresholds.
//! 3. Checks global `max_no_call_fraction` and optional `min_mean_base_quality`.
//! 4. Passes or discards the record.
//!
//! Duplex consensus records (detected by the presence of `aD`/`bD` aux tags) are
//! dispatched to the duplex-aware raw filter. Single-strand records use the simpler
//! per-read filter.

use anyhow::{Result, anyhow};
use fgumi_consensus::filter::{
    FilterConfig, FilterResult, compute_read_stats_raw, filter_duplex_read_raw, filter_read_raw,
    is_duplex_consensus_raw, mask_bases_raw, mask_duplex_bases_raw,
};
use fgumi_raw_bam::fields::aux_data_slice;

use crate::stage::PipelineStage;
use fgumi_consensus::caller::ConsensusOutput;

// ============================================================================
// FilterStage
// ============================================================================

/// Pipeline stage that applies consensus-read filters to `ConsensusOutput`.
///
/// The stage iterates over the framed BAM records embedded in a `ConsensusOutput`
/// buffer, applies base masking and per-read quality filters according to `config`,
/// and returns raw bytes for records that pass.
pub struct FilterStage {
    /// Filter thresholds and per-base masking configuration.
    config: FilterConfig,
}

impl FilterStage {
    /// Create a new `FilterStage` with the given filter configuration.
    ///
    /// # Arguments
    ///
    /// * `config` - Filter configuration controlling masking thresholds and filter thresholds.
    #[must_use]
    pub fn new(config: FilterConfig) -> Self {
        Self { config }
    }

    /// Apply masking and filtering to a single raw BAM record.
    ///
    /// Returns `Ok(true)` if the record passes all filters and `Ok(false)` if it should
    /// be discarded. The record bytes may be modified in place (base masking).
    fn process_record(&self, record: &mut [u8]) -> Result<bool> {
        // Determine whether this is a duplex consensus record (presence of aD/bD tags).
        let aux = aux_data_slice(record);
        let is_duplex = is_duplex_consensus_raw(aux);

        // Pre-compute thresholds once to avoid redundant method calls.
        let duplex_t = self.config.duplex_thresholds();
        let ss_t = self.config.effective_single_strand_thresholds();

        // --- Base masking + per-read filter (combined to avoid duplicate dispatch) ---
        let filter_result = if is_duplex {
            if let Some((cc, ab, ba)) = duplex_t {
                mask_duplex_bases_raw(
                    record,
                    cc,
                    ab,
                    ba,
                    self.config.min_base_quality,
                    false, // require_ss_agreement: off — pipeline does not expose this flag yet
                )?;
                let aux = aux_data_slice(record);
                filter_duplex_read_raw(aux, cc, ab, ba)?
            } else {
                FilterResult::Pass
            }
        } else if let Some(thresholds) = ss_t {
            mask_bases_raw(record, thresholds, self.config.min_base_quality)?;
            let aux = aux_data_slice(record);
            filter_read_raw(aux, thresholds)?
        } else {
            FilterResult::Pass
        };

        if filter_result != FilterResult::Pass {
            return Ok(false);
        }

        // --- Global per-read checks (no-call fraction, mean base quality) ---
        let (no_call_count, mean_base_qual) = compute_read_stats_raw(record);

        let seq_len = fgumi_raw_bam::fields::l_seq(record) as usize;
        if seq_len > 0 {
            #[expect(
                clippy::cast_precision_loss,
                reason = "precision loss is acceptable for no-call fraction comparison"
            )]
            let no_call_fraction = no_call_count as f64 / seq_len as f64;
            if no_call_fraction > self.config.max_no_call_fraction {
                return Ok(false);
            }
        }

        if let Some(min_mean_qual) = self.config.min_mean_base_quality {
            if mean_base_qual < min_mean_qual {
                return Ok(false);
            }
        }

        Ok(true)
    }
}

impl PipelineStage for FilterStage {
    /// Input: framed BAM records from consensus calling (each record prefixed with 4-byte
    /// little-endian `block_size`).
    type Input = ConsensusOutput;
    /// Output: framed BAM record bytes (with `block_size` prefixes), ready for
    /// BGZF compression.
    type Output = Vec<u8>;

    /// Iterate over framed records in `input.data`, apply masking and filtering, and collect
    /// the framed bytes (with `block_size` prefix) of passing records into the output buffer.
    ///
    /// # Errors
    ///
    /// Returns an error if a record is too short or aux data cannot be parsed.
    fn process(&self, input: Self::Input) -> Result<Self::Output> {
        let mut result: Vec<u8> = Vec::with_capacity(input.data.len());
        let mut data = input.data;
        let mut offset = 0usize;

        while offset < data.len() {
            // Read 4-byte little-endian block_size prefix.
            if offset + 4 > data.len() {
                return Err(anyhow!(
                    "truncated ConsensusOutput: expected 4-byte block_size at offset {offset}"
                ));
            }
            let block_size = u32::from_le_bytes(
                data[offset..offset + 4].try_into().expect("slice is exactly 4 bytes"),
            ) as usize;
            let frame_start = offset;
            offset += 4;

            if offset + block_size > data.len() {
                return Err(anyhow!(
                    "truncated ConsensusOutput: block_size {block_size} extends past buffer end"
                ));
            }

            // Process the record in-place within the data buffer to avoid an extra copy.
            let record = &mut data[offset..offset + block_size];
            if self.process_record(record)? {
                // Write the 4-byte block_size prefix followed by the record bytes.
                // BAM requires each record to be preceded by its block_size.
                result.extend_from_slice(&data[frame_start..offset + block_size]);
            }
            offset += block_size;
        }

        Ok(result)
    }

    /// Estimate memory usage as the number of raw bytes in the output.
    fn output_memory_estimate(&self, output: &Self::Output) -> usize {
        output.len()
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_consensus::filter::{FilterConfig, FilterThresholds};

    // ---- helpers -----------------------------------------------------------

    /// Build a minimal framed raw BAM record (4-byte `block_size` prefix followed by the
    /// record body). The record has a 4-base sequence with qualities all set to 40, no
    /// aux tags, so it should survive the default lenient filter config.
    #[allow(clippy::cast_possible_truncation)]
    fn make_framed_record(name: &[u8], seq: &[u8], qual: &[u8]) -> Vec<u8> {
        assert_eq!(seq.len(), qual.len());
        let seq_len = seq.len();
        let l_read_name = (name.len() + 1) as u8;
        let n_cigar_op: u16 = 1;
        let seq_bytes = seq_len.div_ceil(2);

        // BAM fixed header (32) + read name + 1 cigar op (4 bytes) + seq + qual
        let body_len = 32 + l_read_name as usize + n_cigar_op as usize * 4 + seq_bytes + seq_len;
        let mut body = vec![0u8; body_len];

        // ref_id = 0, pos = 100
        body[0..4].copy_from_slice(&0i32.to_le_bytes());
        body[4..8].copy_from_slice(&100i32.to_le_bytes());
        body[8] = l_read_name;
        body[9] = 60; // mapq
        body[10..12].copy_from_slice(&0u16.to_le_bytes()); // bin
        body[12..14].copy_from_slice(&n_cigar_op.to_le_bytes());
        body[14..16].copy_from_slice(&0u16.to_le_bytes()); // flags (unmapped=0 is fine for filter)
        body[16..20].copy_from_slice(&(seq_len as u32).to_le_bytes());
        body[20..24].copy_from_slice(&0i32.to_le_bytes()); // mate_ref_id
        body[24..28].copy_from_slice(&200i32.to_le_bytes()); // mate_pos
        body[28..32].copy_from_slice(&150i32.to_le_bytes()); // tlen

        // Read name (NUL-terminated)
        let name_start = 32;
        body[name_start..name_start + name.len()].copy_from_slice(name);
        body[name_start + name.len()] = 0;

        // CIGAR: one op encoding seq_len * M (op type 0)
        let cigar_start = name_start + l_read_name as usize;
        let cigar_op = (seq_len as u32) << 4; // M (op kind 0, encoding = len << 4)
        body[cigar_start..cigar_start + 4].copy_from_slice(&cigar_op.to_le_bytes());

        // Sequence (4-bit packed, 2 bases per byte)
        // Base encoding: A=1, C=2, G=4, T=8, N=15
        let base_enc = |b: u8| -> u8 {
            match b {
                b'A' => 1,
                b'C' => 2,
                b'G' => 4,
                b'T' => 8,
                _ => 15,
            }
        };
        let seq_start = cigar_start + n_cigar_op as usize * 4;
        for (i, &b) in seq.iter().enumerate() {
            let byte_idx = seq_start + i / 2;
            if i % 2 == 0 {
                body[byte_idx] = base_enc(b) << 4;
            } else {
                body[byte_idx] |= base_enc(b);
            }
        }

        // Quality scores (Phred, not +33 encoded in BAM)
        let qual_start = seq_start + seq_bytes;
        body[qual_start..qual_start + seq_len].copy_from_slice(qual);

        // Prepend 4-byte block_size prefix
        let mut framed = Vec::with_capacity(4 + body_len);
        framed.extend_from_slice(&(body_len as u32).to_le_bytes());
        framed.extend_from_slice(&body);
        framed
    }

    /// Build a `ConsensusOutput` from a slice of already-framed record buffers.
    fn make_consensus_output(framed_records: &[Vec<u8>]) -> ConsensusOutput {
        let mut data = Vec::new();
        let mut count = 0usize;
        for rec in framed_records {
            data.extend_from_slice(rec);
            count += 1;
        }
        ConsensusOutput { data, count }
    }

    /// A permissive `FilterConfig` with `min_reads` = 1, error rates = 1.0, no masking.
    fn permissive_config() -> FilterConfig {
        FilterConfig::for_single_strand(
            FilterThresholds { min_reads: 1, max_read_error_rate: 1.0, max_base_error_rate: 1.0 },
            None,
            None,
            1.0,
        )
    }

    // ---- tests -------------------------------------------------------------

    #[test]
    fn test_filter_stage_passes_valid_record() {
        let stage = FilterStage::new(permissive_config());

        let rec = make_framed_record(b"read1", b"ACGT", &[40, 40, 40, 40]);
        let input = make_consensus_output(&[rec]);

        let output = stage.process(input).expect("should succeed");
        // The record should pass; output should be non-empty.
        assert!(!output.is_empty(), "passing record should produce output bytes");
    }

    #[test]
    fn test_filter_stage_empty_input() {
        let stage = FilterStage::new(permissive_config());
        let input = ConsensusOutput { data: vec![], count: 0 };
        let output = stage.process(input).expect("empty input should succeed");
        assert!(output.is_empty());
    }

    #[test]
    fn test_filter_stage_multiple_records_all_pass() {
        let stage = FilterStage::new(permissive_config());

        let rec1 = make_framed_record(b"read1", b"ACGT", &[40, 40, 40, 40]);
        let rec2 = make_framed_record(b"read2", b"TTTT", &[35, 35, 35, 35]);
        let input = make_consensus_output(&[rec1, rec2]);

        let output = stage.process(input).expect("should succeed");
        assert!(!output.is_empty(), "both records should pass");
    }

    #[test]
    fn test_filter_stage_output_memory_estimate() {
        let stage = FilterStage::new(permissive_config());
        let output: Vec<u8> = vec![0u8; 256];
        assert_eq!(stage.output_memory_estimate(&output), 256);
    }

    #[test]
    fn test_filter_stage_truncated_input_returns_error() {
        let stage = FilterStage::new(permissive_config());
        // Only 2 bytes — not enough for a 4-byte block_size prefix.
        let input = ConsensusOutput { data: vec![0u8, 1u8], count: 1 };
        assert!(stage.process(input).is_err(), "truncated input should error");
    }
}
