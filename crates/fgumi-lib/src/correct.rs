//! Core UMI correction logic for matching observed UMIs to a fixed set of known sequences.
//!
//! This module provides the pure-function core of UMI correction, extracted so that both
//! the `correct` command and the `pipeline` command can call it directly without spawning
//! a subprocess.
//!
//! # Overview
//!
//! The main workflow is:
//! 1. Load known UMI sequences with [`load_umi_sequences`]
//! 2. Build an [`EncodedUmiSet`] for efficient matching
//! 3. For each template, extract the UMI with [`extract_and_validate_template_umi_raw`]
//! 4. Compute the correction with [`compute_template_correction`]
//! 5. Apply the correction with [`apply_correction_to_raw`]

use anyhow::{Result, bail};
use log::info;

use crate::bitenc::BitEnc;
use crate::dna::reverse_complement_str;
use crate::sort::bam_fields;

/// Result of matching an observed UMI to an expected UMI.
///
/// Contains information about whether the match was acceptable and the details
/// of the best matching UMI.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UmiMatch {
    /// Whether the match is acceptable based on mismatch and distance thresholds.
    pub matched: bool,

    /// The fixed UMI sequence that was the closest match.
    pub umi: String,

    /// The number of mismatches between the observed and best matching UMI.
    pub mismatches: usize,
}

/// Reason a UMI correction was rejected.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum RejectionReason {
    /// UMI segments have the wrong length.
    WrongLength,
    /// No acceptable match found in the UMI set.
    Mismatched,
    /// Not rejected.
    #[default]
    None,
}

/// Result of UMI correction for a template.
///
/// This struct holds the result of correcting a UMI once for an entire template,
/// which can then be applied to all records in that template.
#[derive(Debug)]
pub struct TemplateCorrection {
    /// Whether the UMI matched successfully.
    pub matched: bool,
    /// The corrected UMI string (if matched).
    pub corrected_umi: Option<String>,
    /// The original UMI string.
    pub original_umi: String,
    /// Whether correction was needed (mismatches > 0 or revcomp).
    pub needs_correction: bool,
    /// Whether there were actual mismatches (not just revcomp).
    pub has_mismatches: bool,
    /// Match details for metrics.
    pub matches: Vec<UmiMatch>,
    /// Rejection reason if not matched.
    pub rejection_reason: RejectionReason,
}

/// Pre-encoded set of known UMI sequences for fast matching.
///
/// Stores both the original byte sequences and their `BitEnc` representations
/// for efficient Hamming distance computation.
#[derive(Clone)]
pub struct EncodedUmiSet {
    /// Original byte sequences for fallback comparison
    bytes: Vec<Vec<u8>>,
    /// Bit-encoded sequences for fast comparison (None if encoding failed)
    encoded: Vec<Option<BitEnc>>,
    /// Original strings for output
    strings: Vec<String>,
}

impl EncodedUmiSet {
    /// Create a new `EncodedUmiSet` from string sequences.
    ///
    /// All sequences are converted to uppercase for case-insensitive matching.
    /// Pre-encodes all sequences that contain only ACGT bases.
    /// Sequences with other characters will use fallback byte comparison.
    #[must_use]
    pub fn new(sequences: &[String]) -> Self {
        // Store uppercase bytes for comparison
        let bytes: Vec<Vec<u8>> =
            sequences.iter().map(|s| s.bytes().map(|b| b.to_ascii_uppercase()).collect()).collect();
        let encoded: Vec<Option<BitEnc>> = bytes.iter().map(|b| BitEnc::from_bytes(b)).collect();
        // Store uppercase strings for output
        let strings: Vec<String> = sequences.iter().map(|s| s.to_uppercase()).collect();

        Self { bytes, encoded, strings }
    }

    /// Get the number of UMIs in the set.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.bytes.len()
    }

    /// Check if the set is empty.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.bytes.is_empty()
    }

    /// Get the original UMI strings (uppercase).
    #[inline]
    #[must_use]
    pub fn strings(&self) -> &[String] {
        &self.strings
    }
}

/// Count mismatches between two byte slices, stopping early once `max_mismatches` is exceeded.
///
/// # Arguments
///
/// * `a` - First sequence
/// * `b` - Second sequence
/// * `max_mismatches` - Stop counting after this many mismatches
///
/// # Returns
///
/// The number of mismatches, capped at `max_mismatches + 1`.
#[must_use]
pub fn count_mismatches_with_max(a: &[u8], b: &[u8], max_mismatches: usize) -> usize {
    let mut mismatches = 0;
    let min_len = a.len().min(b.len());

    for i in 0..min_len {
        if a[i] != b[i] {
            mismatches += 1;
            if mismatches > max_mismatches {
                return mismatches;
            }
        }
    }

    // Count any remaining bases in the longer sequence as mismatches
    mismatches += a.len().abs_diff(b.len());
    mismatches
}

/// Find the best matching UMI using bit-encoded comparison.
///
/// Uses fast XOR + popcount comparison when both observed and expected UMIs
/// can be bit-encoded. Falls back to byte comparison otherwise.
///
/// This function tracks the best match index during iteration and only clones
/// the matching string once at the end, avoiding repeated string allocations.
#[must_use]
pub fn find_best_match_encoded(
    observed: &[u8],
    umi_set: &EncodedUmiSet,
    max_mismatches: usize,
    min_distance_diff: usize,
) -> UmiMatch {
    let mut best_index: Option<usize> = None;
    let mut best_mismatches = usize::MAX;
    let mut second_best_mismatches = usize::MAX;

    // Try to encode the observed UMI
    let observed_encoded = BitEnc::from_bytes(observed);

    match observed_encoded {
        Some(obs_enc) => {
            // Fast path: use bit-encoded comparison
            for (i, fixed_enc) in umi_set.encoded.iter().enumerate() {
                let mismatches = if let Some(enc) = fixed_enc {
                    // Both can be compared with BitEnc
                    enc.hamming_distance(&obs_enc) as usize
                } else {
                    // Fixed UMI couldn't be encoded, fall back to bytes
                    count_mismatches_with_max(observed, &umi_set.bytes[i], second_best_mismatches)
                };

                if mismatches < best_mismatches {
                    second_best_mismatches = best_mismatches;
                    best_mismatches = mismatches;
                    best_index = Some(i);
                } else if mismatches < second_best_mismatches {
                    second_best_mismatches = mismatches;
                }
            }
        }
        None => {
            // Slow path: observed UMI can't be encoded, use byte comparison
            for (i, fixed_umi) in umi_set.bytes.iter().enumerate() {
                let mismatches =
                    count_mismatches_with_max(observed, fixed_umi, second_best_mismatches);

                if mismatches < best_mismatches {
                    second_best_mismatches = best_mismatches;
                    best_mismatches = mismatches;
                    best_index = Some(i);
                } else if mismatches < second_best_mismatches {
                    second_best_mismatches = mismatches;
                }
            }
        }
    }

    // Build result with single clone at end
    let Some(idx) = best_index else {
        return UmiMatch { matched: false, umi: String::new(), mismatches: usize::MAX };
    };

    let matched = if best_mismatches <= max_mismatches {
        let distance_to_second = second_best_mismatches.saturating_sub(best_mismatches);
        distance_to_second >= min_distance_diff
    } else {
        false
    };

    UmiMatch { matched, umi: umi_set.strings[idx].clone(), mismatches: best_mismatches }
}

/// Finds pairs of UMIs within a specified edit distance.
///
/// Identifies all pairs of UMIs that are within `distance` edits of each other,
/// which helps detect potentially ambiguous UMI sets.
///
/// # Arguments
///
/// * `umis` - Slice of UMI sequences to check
/// * `distance` - Maximum edit distance threshold
///
/// # Returns
///
/// A vector of tuples `(umi1, umi2, actual_distance)` for pairs within the threshold.
#[must_use]
pub fn find_umi_pairs_within_distance(
    umis: &[String],
    distance: usize,
) -> Vec<(String, String, usize)> {
    let mut pairs = Vec::new();
    for i in 0..umis.len() {
        for j in (i + 1)..umis.len() {
            let d = count_mismatches_with_max(umis[i].as_bytes(), umis[j].as_bytes(), distance + 1);
            if d <= distance {
                pairs.push((umis[i].clone(), umis[j].clone(), d));
            }
        }
    }
    pairs
}

/// Load UMI sequences from inline strings and/or files.
///
/// Collects UMIs from both inline values and files (one UMI per line), deduplicates them,
/// validates that all UMIs have the same length, and returns the sorted list with its length.
///
/// # Arguments
///
/// * `umis` - Inline UMI sequence strings
/// * `umi_files` - Paths to files containing UMI sequences (one per line)
///
/// # Returns
///
/// A tuple of (sorted UMI sequences, UMI length).
///
/// # Errors
///
/// Returns an error if:
/// - No UMIs are provided
/// - UMI files cannot be read
/// - UMIs have different lengths
pub fn load_umi_sequences(
    umis: &[String],
    umi_files: &[std::path::PathBuf],
) -> Result<(Vec<String>, usize)> {
    let mut umi_set: std::collections::HashSet<String> =
        umis.iter().map(|s| s.to_uppercase()).collect();

    for file in umi_files {
        let content = std::fs::read_to_string(file)?;
        for line in content.lines() {
            let umi = line.trim().to_uppercase();
            if !umi.is_empty() {
                umi_set.insert(umi);
            }
        }
    }

    if umi_set.is_empty() {
        bail!("No UMIs provided.");
    }

    let mut umi_sequences: Vec<String> = umi_set.into_iter().collect();
    umi_sequences.sort_unstable();

    // Check all UMIs have the same length
    let first_len = umi_sequences[0].len();
    if !umi_sequences.iter().all(|u| u.len() == first_len) {
        bail!("All UMIs must have the same length.");
    }

    info!("Loaded {} UMI sequences of length {}", umi_sequences.len(), first_len);
    Ok((umi_sequences, first_len))
}

/// Compute UMI correction for a template (called once per template).
///
/// Splits the UMI by hyphens (for duplex UMIs), optionally reverse-complements each
/// segment, and finds the best match in the `EncodedUmiSet` for each segment.
///
/// # Arguments
///
/// * `umi` - The observed UMI string (may be hyphen-separated for duplex)
/// * `umi_length` - Expected length of each UMI segment
/// * `revcomp` - Whether to reverse-complement UMI segments before matching
/// * `max_mismatches` - Maximum allowed mismatches per segment
/// * `min_distance_diff` - Minimum difference between best and second-best match
/// * `encoded_umi_set` - Pre-encoded set of known UMI sequences
/// * `cache` - Optional LRU cache for memoizing per-segment match results
#[allow(clippy::too_many_arguments)]
pub fn compute_template_correction(
    umi: &str,
    umi_length: usize,
    revcomp: bool,
    max_mismatches: usize,
    min_distance_diff: usize,
    encoded_umi_set: &EncodedUmiSet,
    cache: &mut Option<lru::LruCache<Vec<u8>, UmiMatch>>,
) -> TemplateCorrection {
    let original_umi = umi.to_string();

    // Split and optionally reverse complement
    let sequences: Vec<String> = if revcomp {
        umi.split('-').map(reverse_complement_str).rev().collect()
    } else {
        umi.split('-').map(std::string::ToString::to_string).collect()
    };

    // Check length
    if sequences.iter().any(|s| s.len() != umi_length) {
        return TemplateCorrection {
            matched: false,
            corrected_umi: None,
            original_umi,
            needs_correction: false,
            has_mismatches: false,
            matches: Vec::new(),
            rejection_reason: RejectionReason::WrongLength,
        };
    }

    // Find matches for each UMI segment
    let mut matches = Vec::with_capacity(sequences.len());
    for seq in &sequences {
        // Uppercase using bytes directly - avoids String allocation
        let seq_bytes: Vec<u8> = seq.bytes().map(|b| b.to_ascii_uppercase()).collect();

        let umi_match = if let Some(c) = cache {
            if let Some(cached) = c.get(&seq_bytes[..]) {
                cached.clone()
            } else {
                let result = find_best_match_encoded(
                    &seq_bytes,
                    encoded_umi_set,
                    max_mismatches,
                    min_distance_diff,
                );
                c.put(seq_bytes, result.clone());
                result
            }
        } else {
            find_best_match_encoded(&seq_bytes, encoded_umi_set, max_mismatches, min_distance_diff)
        };

        matches.push(umi_match);
    }

    // Determine if all segments matched
    let all_matched = matches.iter().all(|m| m.matched);
    let has_mismatches = matches.iter().any(|m| m.mismatches > 0);
    let needs_correction = has_mismatches || revcomp;

    if all_matched {
        let corrected_umi: String =
            matches.iter().map(|m| m.umi.clone()).collect::<Vec<_>>().join("-");
        TemplateCorrection {
            matched: true,
            corrected_umi: Some(corrected_umi),
            original_umi,
            needs_correction,
            has_mismatches,
            matches,
            rejection_reason: RejectionReason::None,
        }
    } else {
        TemplateCorrection {
            matched: false,
            corrected_umi: None,
            original_umi,
            needs_correction: false,
            has_mismatches: false,
            matches,
            rejection_reason: RejectionReason::Mismatched,
        }
    }
}

/// Extract and validate UMI from raw-byte records in a template.
///
/// Reads the UMI tag from the first record and validates that all other records
/// in the template have the same UMI value.
///
/// # Arguments
///
/// * `raw_records` - Slice of raw BAM record byte vectors
/// * `umi_tag` - Two-byte UMI tag (e.g., `b"RX"`)
///
/// # Returns
///
/// `Some(umi)` if a UMI was found, `None` if no records have the tag.
///
/// # Panics
///
/// Panics if records have mismatched UMIs, inconsistent UMI presence, or a non-string
/// UMI tag type.
#[must_use]
pub fn extract_and_validate_template_umi_raw(
    raw_records: &[Vec<u8>],
    umi_tag: [u8; 2],
) -> Option<String> {
    if raw_records.is_empty() {
        return None;
    }

    // Guard against truncated records
    if raw_records.iter().any(|r| r.len() < 32) {
        return None;
    }

    let first_aux = bam_fields::aux_data_slice(&raw_records[0]);
    // Work with &[u8] slices to avoid per-record String allocation
    let first_umi_bytes = bam_fields::find_string_tag(first_aux, &umi_tag);

    // If tag exists but is not a Z-type string, panic (matching RecordBuf behavior)
    if first_umi_bytes.is_none() {
        if let Some(tag_type) = bam_fields::find_tag_type(first_aux, &umi_tag) {
            panic!(
                "UMI tag {:?} exists but has non-string type '{}', expected 'Z'",
                std::str::from_utf8(&umi_tag).unwrap_or("??"),
                tag_type as char,
            );
        }
    }

    for raw in &raw_records[1..] {
        let aux = bam_fields::aux_data_slice(raw);
        let current_umi_bytes = bam_fields::find_string_tag(aux, &umi_tag);

        match (first_umi_bytes, current_umi_bytes) {
            (Some(first), Some(current)) if first != current => {
                panic!(
                    "Template has mismatched UMIs: first={:?}, current={:?}",
                    String::from_utf8_lossy(first),
                    String::from_utf8_lossy(current)
                );
            }
            (Some(_), None) | (None, Some(_)) => {
                panic!("Template has inconsistent UMI presence across records");
            }
            _ => {}
        }
    }

    // Only allocate String once at the end
    first_umi_bytes.map(|b| String::from_utf8_lossy(b).into_owned())
}

/// Apply UMI correction to a raw BAM record.
///
/// Updates the UMI tag in-place with the corrected value, and optionally stores
/// the original UMI in the `OX` tag.
///
/// # Arguments
///
/// * `record` - Mutable raw BAM record bytes
/// * `correction` - The computed correction result
/// * `umi_tag` - Two-byte UMI tag (e.g., `b"RX"`)
/// * `dont_store_original_umis` - If true, skip storing the original UMI in `OX`
pub fn apply_correction_to_raw(
    record: &mut Vec<u8>,
    correction: &TemplateCorrection,
    umi_tag: [u8; 2],
    dont_store_original_umis: bool,
) {
    if correction.needs_correction {
        // Write corrected UMI first (in-place update avoids scanning past OX)
        if let Some(ref corrected) = correction.corrected_umi {
            bam_fields::update_string_tag(record, &umi_tag, corrected.as_bytes());
        }

        // Store original UMI if there were actual mismatches
        // Use update_string_tag to avoid duplicate OX tags
        if !dont_store_original_umis && correction.has_mismatches {
            bam_fields::update_string_tag(record, b"OX", correction.original_umi.as_bytes());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FIXED_UMIS: &[&str] = &["AAAAAA", "CCCCCC", "GGGGGG", "TTTTTT"];

    fn get_encoded_umi_set() -> EncodedUmiSet {
        let strings: Vec<String> = FIXED_UMIS.iter().map(|s| (*s).to_string()).collect();
        EncodedUmiSet::new(&strings)
    }

    #[test]
    fn test_find_best_match_perfect() {
        let umi_set = get_encoded_umi_set();

        let hit = find_best_match_encoded(b"AAAAAA", &umi_set, 2, 2);
        assert!(hit.matched);
        assert_eq!(hit.mismatches, 0);
        assert_eq!(hit.umi, "AAAAAA");
    }

    #[test]
    fn test_find_best_match_with_mismatches() {
        let umi_set = get_encoded_umi_set();

        let hit = find_best_match_encoded(b"AAAAAC", &umi_set, 2, 2);
        assert!(hit.matched);
        assert_eq!(hit.mismatches, 1);
        assert_eq!(hit.umi, "AAAAAA");
    }

    #[test]
    fn test_find_best_match_too_many_mismatches() {
        let umi_set = get_encoded_umi_set();

        let hit = find_best_match_encoded(b"ACGTAC", &umi_set, 1, 1);
        assert!(!hit.matched);
    }

    #[test]
    fn test_count_mismatches_identical() {
        assert_eq!(count_mismatches_with_max(b"AAAAAA", b"AAAAAA", 10), 0);
    }

    #[test]
    fn test_count_mismatches_one() {
        assert_eq!(count_mismatches_with_max(b"AAAAAA", b"AAAAAT", 10), 1);
    }

    #[test]
    fn test_count_mismatches_early_exit() {
        assert_eq!(count_mismatches_with_max(b"AAAAAA", b"CCCCCC", 2), 3);
    }

    #[test]
    fn test_find_umi_pairs_within_distance_finds_close_pair() {
        let umis = vec!["AAAA".to_string(), "AAAT".to_string()];
        let pairs = find_umi_pairs_within_distance(&umis, 2);
        assert_eq!(pairs.len(), 1);
    }

    #[test]
    fn test_find_umi_pairs_within_distance_no_close_pair() {
        let umis = vec!["AAAA".to_string(), "CCCC".to_string()];
        let pairs = find_umi_pairs_within_distance(&umis, 1);
        assert!(pairs.is_empty());
    }

    #[test]
    fn test_compute_template_correction_match() {
        let umi_set = get_encoded_umi_set();
        let result = compute_template_correction("AAAAAC", 6, false, 2, 2, &umi_set, &mut None);
        assert!(result.matched);
        assert_eq!(result.corrected_umi, Some("AAAAAA".to_string()));
        assert!(result.has_mismatches);
        assert!(result.needs_correction);
    }

    #[test]
    fn test_compute_template_correction_wrong_length() {
        let umi_set = get_encoded_umi_set();
        let result = compute_template_correction("AAA", 6, false, 2, 2, &umi_set, &mut None);
        assert!(!result.matched);
        assert_eq!(result.rejection_reason, RejectionReason::WrongLength);
    }

    #[test]
    fn test_load_umi_sequences_from_strings() {
        let umis = vec!["AAAA".to_string(), "CCCC".to_string(), "aaaa".to_string()];
        let (seqs, len) = load_umi_sequences(&umis, &[]).unwrap();
        assert_eq!(len, 4);
        // Should deduplicate AAAA and aaaa (both uppercase to AAAA)
        assert_eq!(seqs.len(), 2);
    }

    #[test]
    fn test_load_umi_sequences_empty_fails() {
        let result = load_umi_sequences(&[], &[]);
        assert!(result.is_err());
    }
}
