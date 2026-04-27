//! Template information types for UMI metrics collection.
//!
//! Provides the core data structures for tracking template-level information
//! during metrics computation, including grouping keys and strand metadata.

/// Read name and template information for downsampling.
#[derive(Clone)]
pub struct TemplateInfo {
    /// Molecular identifier tag value (e.g. "1/A").
    pub mi: String,
    /// Raw UMI tag value (e.g. "AAA-TTT").
    pub rx: String,
    /// Reference sequence name, if mapped.
    pub ref_name: Option<String>,
    /// Alignment start position (1-based), if mapped.
    pub position: Option<i32>,
    /// Alignment end position (1-based), if mapped.
    pub end_position: Option<i32>,
    /// Hash fraction for deterministic downsampling (computed once per template).
    pub hash_fraction: f64,
    /// Coordinate+strand grouping key matching fgbio's `ReadInfo` structure.
    /// Used to sub-group templates within a position group by their read pair geometry.
    pub read_info_key: ReadInfoKey,
}

/// Grouping key matching fgbio's `ReadInfo` structure.
///
/// Fields are ordered so the earlier-mapping read comes first.
#[derive(Clone, Default, PartialEq, Eq, Hash)]
pub struct ReadInfoKey {
    /// Reference sequence index for read 1.
    pub ref_index1: Option<usize>,
    /// Unclipped 5' position for read 1.
    pub start1: i32,
    /// `true` if read 1 is reverse-complemented.
    pub strand1: bool,
    /// Reference sequence index for read 2.
    pub ref_index2: Option<usize>,
    /// Unclipped 5' position for read 2.
    pub start2: i32,
    /// `true` if read 2 is reverse-complemented.
    pub strand2: bool,
}

/// Pre-computed metadata for a template within a coordinate group.
pub struct TemplateMetadata<'a> {
    /// Reference to the underlying template info.
    pub template: &'a TemplateInfo,
    /// MI tag value with strand suffix stripped (e.g. "1" from "1/A").
    pub base_umi: &'a str,
    /// `true` if this template belongs to the A strand.
    pub is_a_strand: bool,
    /// `true` if this template belongs to the B strand.
    pub is_b_strand: bool,
}

/// Pre-computes metadata for each template in a coordinate group.
///
/// Parses the MI tag to determine strand assignment and extract the base UMI
/// (MI value without the `/A` or `/B` suffix).
#[must_use]
pub fn compute_template_metadata(group: &[TemplateInfo]) -> Vec<TemplateMetadata<'_>> {
    group
        .iter()
        .map(|t| {
            let (base_umi, is_a, is_b) = if t.mi.ends_with("/A") {
                (&t.mi[..t.mi.len() - 2], true, false)
            } else if t.mi.ends_with("/B") {
                (&t.mi[..t.mi.len() - 2], false, true)
            } else {
                (t.mi.as_str(), false, false)
            };
            TemplateMetadata { template: t, base_umi, is_a_strand: is_a, is_b_strand: is_b }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_template(mi: &str) -> TemplateInfo {
        TemplateInfo {
            mi: mi.to_string(),
            rx: "AAA-TTT".to_string(),
            ref_name: Some("chr1".to_string()),
            position: Some(100),
            end_position: Some(200),
            hash_fraction: 0.5,
            read_info_key: ReadInfoKey::default(),
        }
    }

    #[test]
    fn test_compute_template_metadata_a_strand() {
        let templates = vec![make_template("42/A")];
        let metadata = compute_template_metadata(&templates);
        assert_eq!(metadata.len(), 1);
        assert_eq!(metadata[0].base_umi, "42");
        assert!(metadata[0].is_a_strand);
        assert!(!metadata[0].is_b_strand);
    }

    #[test]
    fn test_compute_template_metadata_b_strand() {
        let templates = vec![make_template("42/B")];
        let metadata = compute_template_metadata(&templates);
        assert_eq!(metadata.len(), 1);
        assert_eq!(metadata[0].base_umi, "42");
        assert!(!metadata[0].is_a_strand);
        assert!(metadata[0].is_b_strand);
    }

    #[test]
    fn test_compute_template_metadata_no_suffix() {
        let templates = vec![make_template("42")];
        let metadata = compute_template_metadata(&templates);
        assert_eq!(metadata.len(), 1);
        assert_eq!(metadata[0].base_umi, "42");
        assert!(!metadata[0].is_a_strand);
        assert!(!metadata[0].is_b_strand);
    }

    #[test]
    fn test_compute_template_metadata_mixed() {
        let templates = vec![make_template("1/A"), make_template("1/B"), make_template("2")];
        let metadata = compute_template_metadata(&templates);
        assert_eq!(metadata.len(), 3);
        assert_eq!(metadata[0].base_umi, "1");
        assert!(metadata[0].is_a_strand);
        assert_eq!(metadata[1].base_umi, "1");
        assert!(metadata[1].is_b_strand);
        assert_eq!(metadata[2].base_umi, "2");
        assert!(!metadata[2].is_a_strand);
        assert!(!metadata[2].is_b_strand);
    }
}
