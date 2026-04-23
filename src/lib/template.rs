//! Template data structure for grouping reads by query name.
//!
//! This module provides the `Template` structure for organizing all SAM/BAM records
//! that share the same query name. In paired-end sequencing, a template typically
//! contains:
//! - Primary R1 and R2 alignments
//! - Supplementary alignments (chimeric/split alignments)
//! - Secondary alignments (alternative mappings)
//!
//! # Structure
//!
//! The `Template` organizes reads into categories:
//! - **Primary reads**: The main alignment for R1 and R2
//! - **Supplementary**: Additional alignments for chimeric reads (same read, different locations)
//! - **Secondary**: Alternative alignments at different confidence levels
//!
//! # Usage
//!
//! Templates can be built incrementally using `Builder`, created from an iterator
//! of records, or iterated over from a BAM file using `TemplateIterator`.
//!
//! # Examples
//!
//! ```rust,ignore
//! // Build a template from records
//! let mut builder = Builder::new();
//! for record in records {
//!     builder.push(record)?;
//! }
//! let template = builder.build()?;
//!
//! // Access primary reads
//! if let Some(r1) = template.r1() {
//!     // Process R1
//! }
//! ```

use crate::sam::SamTag;
use crate::unified_pipeline::MemoryEstimate;
use anyhow::{Result, bail};
use fgumi_raw_bam::{RawRecord, RawRecordView};

// Re-export PairOrientation from sam module for backwards compatibility
pub use crate::sam::PairOrientation;

// MoleculeId is now defined in fgumi-umi crate
pub use fgumi_umi::MoleculeId;

/// A template representing all reads with the same query name.
///
/// In paired-end sequencing, this typically includes R1 and R2, plus any secondary
/// or supplementary alignments. All records with the same name are grouped together.
///
/// This structure mirrors the Scala Template class from fgbio, organizing reads into:
/// - Primary R1 and R2 reads
/// - Supplementary alignments for R1 and R2
/// - Secondary alignments for R1 and R2
#[derive(Debug, Clone)]
pub struct Template {
    /// The query name (QNAME) shared by all records in this template
    pub name: Vec<u8>,
    /// Raw BAM records (without `block_size` prefix).
    pub records: Vec<RawRecord>,
    /// Primary R1 read (first segment, non-secondary, non-supplementary)
    pub r1: Option<(usize, usize)>,
    /// Primary R2 read (second segment, non-secondary, non-supplementary)
    pub r2: Option<(usize, usize)>,
    /// Supplementary alignments for R1
    pub r1_supplementals: Option<(usize, usize)>,
    /// Supplementary alignments for R2
    pub r2_supplementals: Option<(usize, usize)>,
    /// Secondary alignments for R1
    pub r1_secondaries: Option<(usize, usize)>,
    /// Secondary alignments for R2
    pub r2_secondaries: Option<(usize, usize)>,
    /// Assigned molecule ID from UMI grouping (set during group command)
    pub mi: MoleculeId,
}

/// A batch of templates for parallel processing.
pub type TemplateBatch = Vec<Template>;

impl Template {
    /// Creates a new empty template with the given name.
    ///
    /// # Arguments
    ///
    /// * `name` - The query name for this template
    #[must_use]
    pub fn new(name: Vec<u8>) -> Self {
        Template {
            name,
            records: Vec::new(),
            r1: None,
            r2: None,
            r1_supplementals: None,
            r2_supplementals: None,
            r1_secondaries: None,
            r2_secondaries: None,
            mi: MoleculeId::None,
        }
    }

    /// Returns the primary R1 record if present.
    #[must_use]
    pub fn r1(&self) -> Option<&RawRecord> {
        self.r1.map(|_| &self.records[0])
    }

    /// Returns the primary R2 record if present.
    #[must_use]
    pub fn r2(&self) -> Option<&RawRecord> {
        self.r2.map(|(i, _)| &self.records[i])
    }

    /// Returns all records as a slice.
    #[must_use]
    pub fn records(&self) -> &[RawRecord] {
        &self.records
    }

    /// Returns mutable access to all records.
    pub fn records_mut(&mut self) -> &mut [RawRecord] {
        &mut self.records
    }

    /// Consumes self and returns all records as a `Vec`.
    #[must_use]
    pub fn into_records(self) -> Vec<RawRecord> {
        self.records
    }

    /// Returns an iterator over primary (non-secondary, non-supplementary) records.
    ///
    /// Yields the R1 primary record and/or the R2 primary record if present.
    pub fn primary_reads(&self) -> impl Iterator<Item = &RawRecord> {
        let records = &self.records;
        [self.r1.map(|(s, _)| s), self.r2.map(|(s, _)| s)]
            .into_iter()
            .flatten()
            .filter_map(move |start| records.get(start))
    }

    /// Returns the total count of records in this template.
    #[must_use]
    pub fn read_count(&self) -> usize {
        self.records.len()
    }

    /// Builds a `Template` from raw BAM byte records, categorizing by flags.
    ///
    /// The records are categorized using `RawRecordView::flags()` to determine
    /// R1/R2/supplementary/secondary status.
    ///
    /// # Errors
    ///
    /// Returns an error if multiple primary R1s or R2s are found.
    #[allow(clippy::too_many_lines)]
    pub fn from_records(mut raw_records: Vec<RawRecord>) -> Result<Self> {
        use crate::sort::bam_fields;

        if raw_records.is_empty() {
            bail!("No records given to Template::from_records");
        }

        // Guard against truncated records
        for (i, r) in raw_records.iter().enumerate() {
            if r.len() < bam_fields::MIN_BAM_RECORD_LEN {
                bail!(
                    "Raw BAM record {i} too short to parse ({} < {})",
                    r.len(),
                    bam_fields::MIN_BAM_RECORD_LEN
                );
            }
            let l_rn = r[8] as usize;
            if r.len() < 32 + l_rn {
                bail!(
                    "Raw BAM record {i} truncated: l_read_name={l_rn} but only {} bytes after header",
                    r.len() - 32
                );
            }
        }

        // Extract name from first record
        let name = bam_fields::read_name(&raw_records[0]).to_vec();

        // Verify all records share the same QNAME (matching Builder behavior).
        // Done before the fast path so mismatched names cannot slip past into a Template.
        for rec in raw_records.iter().skip(1) {
            if bam_fields::read_name(rec) != name.as_slice() {
                bail!("Template name mismatch in from_raw_records");
            }
        }

        // Fast path for common 2-record paired-end case (no supplementals/secondaries)
        if raw_records.len() == 2 {
            let f0 = RawRecordView::new(&raw_records[0]).flags();
            let f1 = RawRecordView::new(&raw_records[1]).flags();
            let neither_sec_supp =
                (f0 | f1) & (bam_fields::flags::SECONDARY | bam_fields::flags::SUPPLEMENTARY) == 0;
            if neither_sec_supp {
                // Use same R1/R2 logic as general path: R1 = !paired || first_segment
                let is_r1_0 = (f0 & bam_fields::flags::PAIRED) == 0
                    || (f0 & bam_fields::flags::FIRST_SEGMENT) != 0;
                let is_r1_1 = (f1 & bam_fields::flags::PAIRED) == 0
                    || (f1 & bam_fields::flags::FIRST_SEGMENT) != 0;
                // Only use fast path if exactly one R1 and one R2
                if is_r1_0 != is_r1_1 {
                    if !is_r1_0 {
                        // rec[0]=R2, rec[1]=R1 — swap so R1 is at index 0
                        raw_records.swap(0, 1);
                    }
                    return Ok(Template {
                        name,
                        records: raw_records,
                        r1: Some((0, 1)),
                        r2: Some((1, 2)),
                        r1_supplementals: None,
                        r2_supplementals: None,
                        r1_secondaries: None,
                        r2_secondaries: None,
                        mi: MoleculeId::None,
                    });
                }
                // Both R1 or both R2 — fall through to general path for error handling
            }
        }

        // Verify all records share the same QNAME (matching Builder behavior)
        for rec in raw_records.iter().skip(1) {
            if bam_fields::read_name(rec) != name.as_slice() {
                bail!("Template name mismatch in from_records");
            }
        }

        // Categorize records into: r1, r2, r1_supp, r2_supp, r1_sec, r2_sec
        let mut r1_idx: Option<usize> = None;
        let mut r2_idx: Option<usize> = None;
        let mut r1_supp: Vec<usize> = Vec::new();
        let mut r2_supp: Vec<usize> = Vec::new();
        let mut r1_sec: Vec<usize> = Vec::new();
        let mut r2_sec: Vec<usize> = Vec::new();

        for (i, rec) in raw_records.iter().enumerate() {
            let flg = RawRecordView::new(rec).flags();
            let is_secondary = (flg & bam_fields::flags::SECONDARY) != 0;
            let is_supplementary = (flg & bam_fields::flags::SUPPLEMENTARY) != 0;
            let is_paired = (flg & bam_fields::flags::PAIRED) != 0;
            let is_first = (flg & bam_fields::flags::FIRST_SEGMENT) != 0;
            let is_r1 = !is_paired || is_first;

            if is_r1 {
                if is_secondary {
                    r1_sec.push(i);
                } else if is_supplementary {
                    r1_supp.push(i);
                } else if r1_idx.is_some() {
                    bail!(
                        "Multiple non-secondary, non-supplemental R1 records for read '{name:?}'"
                    );
                } else {
                    r1_idx = Some(i);
                }
            } else if is_secondary {
                r2_sec.push(i);
            } else if is_supplementary {
                r2_supp.push(i);
            } else if r2_idx.is_some() {
                bail!("Multiple non-secondary, non-supplemental R2 records for read '{name:?}'");
            } else {
                r2_idx = Some(i);
            }
        }

        // Build ordered Vec: r1, r2, r1_supplementals, r2_supplementals, r1_secondaries, r2_secondaries
        // (supplementals/secondaries reversed to match fgbio's Scala List prepend ordering)
        let mut ordered: Vec<RawRecord> = Vec::with_capacity(raw_records.len());
        let mut take = |idx: usize| -> RawRecord { std::mem::take(&mut raw_records[idx]) };

        // Track indices
        let r1_pair = if let Some(idx) = r1_idx {
            ordered.push(take(idx));
            Some((ordered.len() - 1, ordered.len()))
        } else {
            None
        };

        let r2_pair = if let Some(idx) = r2_idx {
            ordered.push(take(idx));
            Some((ordered.len() - 1, ordered.len()))
        } else {
            None
        };

        let r1_supp_pair = if r1_supp.is_empty() {
            None
        } else {
            r1_supp.reverse();
            let start = ordered.len();
            for idx in &r1_supp {
                ordered.push(take(*idx));
            }
            Some((start, ordered.len()))
        };

        let r2_supp_pair = if r2_supp.is_empty() {
            None
        } else {
            r2_supp.reverse();
            let start = ordered.len();
            for idx in &r2_supp {
                ordered.push(take(*idx));
            }
            Some((start, ordered.len()))
        };

        let r1_sec_pair = if r1_sec.is_empty() {
            None
        } else {
            r1_sec.reverse();
            let start = ordered.len();
            for idx in &r1_sec {
                ordered.push(take(*idx));
            }
            Some((start, ordered.len()))
        };

        let r2_sec_pair = if r2_sec.is_empty() {
            None
        } else {
            r2_sec.reverse();
            let start = ordered.len();
            for idx in &r2_sec {
                ordered.push(take(*idx));
            }
            Some((start, ordered.len()))
        };

        Ok(Template {
            name,
            records: ordered,
            r1: r1_pair,
            r2: r2_pair,
            r1_supplementals: r1_supp_pair,
            r2_supplementals: r2_supp_pair,
            r1_secondaries: r1_sec_pair,
            r2_secondaries: r2_sec_pair,
            mi: MoleculeId::None,
        })
    }

    /// Returns the pair orientation if both R1 and R2 are present and mapped to the same chromosome.
    ///
    /// This matches htsjdk's `SamPairUtil.getPairOrientation()` and fgbio's `Template.pairOrientation`.
    ///
    /// # Returns
    ///
    /// - `Some(PairOrientation::FR)` - Forward-Reverse ("innie") orientation
    /// - `Some(PairOrientation::RF)` - Reverse-Forward ("outie") orientation
    /// - `Some(PairOrientation::Tandem)` - Both reads on the same strand
    /// - `None` - If R1 or R2 is missing, unmapped, or they're on different chromosomes
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// let template = Template::from_records(records)?;
    /// if let Some(orientation) = template.pair_orientation() {
    ///     match orientation {
    ///         PairOrientation::FR => println!("Innie pair"),
    ///         PairOrientation::RF => println!("Outie pair"),
    ///         PairOrientation::Tandem => println!("Tandem pair"),
    ///     }
    /// }
    /// ```
    #[must_use]
    pub fn pair_orientation(&self) -> Option<PairOrientation> {
        use crate::sort::bam_fields;
        let r1 = self.r1()?;
        let r2 = self.r2()?;
        let f1 = RawRecordView::new(r1).flags();
        let f2 = RawRecordView::new(r2).flags();

        // Both must be mapped
        if (f1 & bam_fields::flags::UNMAPPED) != 0 || (f2 & bam_fields::flags::UNMAPPED) != 0 {
            return None;
        }

        // Must be on the same reference
        if bam_fields::ref_id(r1) != bam_fields::ref_id(r2) {
            return None;
        }

        // Use R1 to determine orientation (htsjdk uses the first record passed)
        Some(get_pair_orientation_raw(r1))
    }

    /// Fixes mate information for paired-end reads.
    ///
    /// Modifies `self.records` in place, setting RNEXT, PNEXT, TLEN, mate flags,
    /// MQ, MC, and ms tags for primary pairs and supplementary alignments.
    ///
    /// # Errors
    ///
    /// Returns an error if the records are malformed.
    ///
    /// # Panics
    ///
    /// Panics if internal index invariants are violated (should not happen
    /// with well-formed templates).
    #[allow(clippy::too_many_lines)]
    pub fn fix_mate_info(&mut self) -> Result<()> {
        use crate::sort::bam_fields;

        let rr = &mut self.records;

        if rr.is_empty() {
            return Ok(());
        }

        // Fix mate info for primary R1/R2 pair
        if let (Some((r1_i, _)), Some((r2_i, _))) = (self.r1, self.r2) {
            let r1_is_unmapped =
                (RawRecordView::new(&rr[r1_i]).flags() & bam_fields::flags::UNMAPPED) != 0;
            let r2_is_unmapped =
                (RawRecordView::new(&rr[r2_i]).flags() & bam_fields::flags::UNMAPPED) != 0;

            // Get alignment scores for mate score tags
            let r1_as = bam_fields::find_int_tag(bam_fields::aux_data_slice(&rr[r1_i]), SamTag::AS);
            let r2_as = bam_fields::find_int_tag(bam_fields::aux_data_slice(&rr[r2_i]), SamTag::AS);

            if !r1_is_unmapped && !r2_is_unmapped {
                // Case 1: Both mapped
                self.set_mate_info_both_mapped(r1_i, r2_i);
            } else if r1_is_unmapped && r2_is_unmapped {
                // Case 2: Both unmapped
                self.set_mate_info_both_unmapped(r1_i, r2_i);
            } else {
                // Case 3: One mapped, one unmapped
                self.set_mate_info_one_unmapped(r1_i, r2_i, r1_is_unmapped);
            }

            // Set mate score tags (ms) in all cases
            let rr = &mut self.records;
            if let Some(as_value) = r2_as {
                if let Ok(v) = i32::try_from(as_value) {
                    bam_fields::remove_tag(rr[r1_i].as_mut_vec(), SamTag::MS);
                    bam_fields::append_signed_int_tag(rr[r1_i].as_mut_vec(), SamTag::MS, v);
                }
            }
            if let Some(as_value) = r1_as {
                if let Ok(v) = i32::try_from(as_value) {
                    bam_fields::remove_tag(rr[r2_i].as_mut_vec(), SamTag::MS);
                    bam_fields::append_signed_int_tag(rr[r2_i].as_mut_vec(), SamTag::MS, v);
                }
            }
        }

        // Fix mate info for R1 supplementals (mate is primary R2)
        if let Some((r2_i, _)) = self.r2 {
            let rr = &self.records;
            let r2_ref_id = bam_fields::ref_id(&rr[r2_i]);
            let r2_pos = bam_fields::pos(&rr[r2_i]);
            let r2_flags = RawRecordView::new(&rr[r2_i]).flags();
            let r2_is_reverse = (r2_flags & bam_fields::flags::REVERSE) != 0;
            let r2_is_unmapped = (r2_flags & bam_fields::flags::UNMAPPED) != 0;
            let r2_tlen = bam_fields::template_length(&rr[r2_i]);
            let r2_mapq = bam_fields::mapq(&rr[r2_i]);
            let r2_cigar_str = bam_fields::cigar_to_string_from_raw(&rr[r2_i]);
            let r2_as = bam_fields::find_int_tag(bam_fields::aux_data_slice(&rr[r2_i]), SamTag::AS);

            if let Some((start, end)) = self.r1_supplementals {
                let rr = &mut self.records;
                for rec in &mut rr[start..end] {
                    bam_fields::set_mate_ref_id(rec, r2_ref_id);
                    bam_fields::set_mate_pos(rec, r2_pos);
                    set_mate_flags(rec, r2_is_reverse, r2_is_unmapped);
                    bam_fields::set_template_length(rec, -r2_tlen);

                    let mq_val = if r2_mapq == 255 { 255 } else { i32::from(r2_mapq) };
                    bam_fields::update_int_tag(rec.as_mut_vec(), SamTag::MQ, mq_val);

                    if !r2_cigar_str.is_empty() && r2_cigar_str != "*" && !r2_is_unmapped {
                        bam_fields::update_string_tag(
                            rec.as_mut_vec(),
                            SamTag::MC,
                            r2_cigar_str.as_bytes(),
                        );
                    } else {
                        bam_fields::remove_tag(rec.as_mut_vec(), SamTag::MC);
                    }

                    if let Some(as_value) = r2_as {
                        if let Ok(v) = i32::try_from(as_value) {
                            bam_fields::remove_tag(rec.as_mut_vec(), SamTag::MS);
                            bam_fields::append_signed_int_tag(rec.as_mut_vec(), SamTag::MS, v);
                        }
                    }
                }
            }
        }

        // Fix mate info for R2 supplementals (mate is primary R1)
        if let Some((r1_i, _)) = self.r1 {
            let rr = &self.records;
            let r1_ref_id = bam_fields::ref_id(&rr[r1_i]);
            let r1_pos = bam_fields::pos(&rr[r1_i]);
            let r1_flags = RawRecordView::new(&rr[r1_i]).flags();
            let r1_is_reverse = (r1_flags & bam_fields::flags::REVERSE) != 0;
            let r1_is_unmapped = (r1_flags & bam_fields::flags::UNMAPPED) != 0;
            let r1_tlen = bam_fields::template_length(&rr[r1_i]);
            let r1_mapq = bam_fields::mapq(&rr[r1_i]);
            let r1_cigar_str = bam_fields::cigar_to_string_from_raw(&rr[r1_i]);
            let r1_as = bam_fields::find_int_tag(bam_fields::aux_data_slice(&rr[r1_i]), SamTag::AS);

            if let Some((start, end)) = self.r2_supplementals {
                let rr = &mut self.records;
                for rec in &mut rr[start..end] {
                    bam_fields::set_mate_ref_id(rec, r1_ref_id);
                    bam_fields::set_mate_pos(rec, r1_pos);
                    set_mate_flags(rec, r1_is_reverse, r1_is_unmapped);
                    bam_fields::set_template_length(rec, -r1_tlen);

                    let mq_val = if r1_mapq == 255 { 255 } else { i32::from(r1_mapq) };
                    bam_fields::update_int_tag(rec.as_mut_vec(), SamTag::MQ, mq_val);

                    if !r1_cigar_str.is_empty() && r1_cigar_str != "*" && !r1_is_unmapped {
                        bam_fields::update_string_tag(
                            rec.as_mut_vec(),
                            SamTag::MC,
                            r1_cigar_str.as_bytes(),
                        );
                    } else {
                        bam_fields::remove_tag(rec.as_mut_vec(), SamTag::MC);
                    }

                    if let Some(as_value) = r1_as {
                        if let Ok(v) = i32::try_from(as_value) {
                            bam_fields::remove_tag(rec.as_mut_vec(), SamTag::MS);
                            bam_fields::append_signed_int_tag(rec.as_mut_vec(), SamTag::MS, v);
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Raw-byte: both reads mapped. Sets mate info, TLEN, MQ, MC.
    fn set_mate_info_both_mapped(&mut self, r1_i: usize, r2_i: usize) {
        use crate::sort::bam_fields;

        let rr = &self.records;

        // Get R2's info for R1
        let r2_ref_id = bam_fields::ref_id(&rr[r2_i]);
        let r2_pos = bam_fields::pos(&rr[r2_i]);
        let r2_is_reverse =
            (RawRecordView::new(&rr[r2_i]).flags() & bam_fields::flags::REVERSE) != 0;
        let r2_mapq = bam_fields::mapq(&rr[r2_i]);
        let r2_cigar_str = bam_fields::cigar_to_string_from_raw(&rr[r2_i]);

        // Get R1's info for R2
        let r1_ref_id = bam_fields::ref_id(&rr[r1_i]);
        let r1_pos = bam_fields::pos(&rr[r1_i]);
        let r1_is_reverse =
            (RawRecordView::new(&rr[r1_i]).flags() & bam_fields::flags::REVERSE) != 0;
        let r1_mapq = bam_fields::mapq(&rr[r1_i]);
        let r1_cigar_str = bam_fields::cigar_to_string_from_raw(&rr[r1_i]);

        // Compute insert size before mutating
        let insert_size = compute_insert_size_raw(&rr[r1_i], &rr[r2_i]);

        let rr = &mut self.records;

        // Set mate info on R1 from R2
        bam_fields::set_mate_ref_id(&mut rr[r1_i], r2_ref_id);
        bam_fields::set_mate_pos(&mut rr[r1_i], r2_pos);
        set_mate_flags(&mut rr[r1_i], r2_is_reverse, false);
        let r2_mq_val = if r2_mapq == 255 { 255 } else { i32::from(r2_mapq) };
        bam_fields::update_int_tag(rr[r1_i].as_mut_vec(), SamTag::MQ, r2_mq_val);
        if !r2_cigar_str.is_empty() && r2_cigar_str != "*" {
            bam_fields::update_string_tag(
                rr[r1_i].as_mut_vec(),
                SamTag::MC,
                r2_cigar_str.as_bytes(),
            );
        } else {
            bam_fields::remove_tag(rr[r1_i].as_mut_vec(), SamTag::MC);
        }

        // Set mate info on R2 from R1
        bam_fields::set_mate_ref_id(&mut rr[r2_i], r1_ref_id);
        bam_fields::set_mate_pos(&mut rr[r2_i], r1_pos);
        set_mate_flags(&mut rr[r2_i], r1_is_reverse, false);
        let r1_mq_val = if r1_mapq == 255 { 255 } else { i32::from(r1_mapq) };
        bam_fields::update_int_tag(rr[r2_i].as_mut_vec(), SamTag::MQ, r1_mq_val);
        if !r1_cigar_str.is_empty() && r1_cigar_str != "*" {
            bam_fields::update_string_tag(
                rr[r2_i].as_mut_vec(),
                SamTag::MC,
                r1_cigar_str.as_bytes(),
            );
        } else {
            bam_fields::remove_tag(rr[r2_i].as_mut_vec(), SamTag::MC);
        }

        // Set insert size
        bam_fields::set_template_length(&mut rr[r1_i], insert_size);
        bam_fields::set_template_length(&mut rr[r2_i], -insert_size);
    }

    /// Raw-byte: both reads unmapped. Clears ref/pos, removes MQ/MC, TLEN=0.
    fn set_mate_info_both_unmapped(&mut self, r1_i: usize, r2_i: usize) {
        use crate::sort::bam_fields;

        let rr = &self.records;
        let r1_is_reverse =
            (RawRecordView::new(&rr[r1_i]).flags() & bam_fields::flags::REVERSE) != 0;
        let r2_is_reverse =
            (RawRecordView::new(&rr[r2_i]).flags() & bam_fields::flags::REVERSE) != 0;

        let rr = &mut self.records;

        // R1: set to unmapped coordinates
        bam_fields::set_ref_id(&mut rr[r1_i], -1);
        bam_fields::set_pos(&mut rr[r1_i], -1);
        bam_fields::set_mate_ref_id(&mut rr[r1_i], -1);
        bam_fields::set_mate_pos(&mut rr[r1_i], -1);
        set_mate_flags(&mut rr[r1_i], r2_is_reverse, true);
        bam_fields::remove_tag(rr[r1_i].as_mut_vec(), SamTag::MQ);
        bam_fields::remove_tag(rr[r1_i].as_mut_vec(), SamTag::MC);
        bam_fields::set_template_length(&mut rr[r1_i], 0);

        // R2: set to unmapped coordinates
        bam_fields::set_ref_id(&mut rr[r2_i], -1);
        bam_fields::set_pos(&mut rr[r2_i], -1);
        bam_fields::set_mate_ref_id(&mut rr[r2_i], -1);
        bam_fields::set_mate_pos(&mut rr[r2_i], -1);
        set_mate_flags(&mut rr[r2_i], r1_is_reverse, true);
        bam_fields::remove_tag(rr[r2_i].as_mut_vec(), SamTag::MQ);
        bam_fields::remove_tag(rr[r2_i].as_mut_vec(), SamTag::MC);
        bam_fields::set_template_length(&mut rr[r2_i], 0);
    }

    /// Raw-byte: one mapped, one unmapped. Places unmapped at mapped coords.
    fn set_mate_info_one_unmapped(&mut self, r1_i: usize, r2_i: usize, r1_is_unmapped: bool) {
        use crate::sort::bam_fields;

        let (mapped_i, unmapped_i) = if r1_is_unmapped { (r2_i, r1_i) } else { (r1_i, r2_i) };

        let rr = &self.records;
        let mapped_ref_id = bam_fields::ref_id(&rr[mapped_i]);
        let mapped_pos = bam_fields::pos(&rr[mapped_i]);
        let mapped_flags = RawRecordView::new(&rr[mapped_i]).flags();
        let mapped_is_reverse = (mapped_flags & bam_fields::flags::REVERSE) != 0;
        let mapped_mapq = bam_fields::mapq(&rr[mapped_i]);
        let mapped_cigar_str = bam_fields::cigar_to_string_from_raw(&rr[mapped_i]);

        let unmapped_is_reverse =
            (RawRecordView::new(&rr[unmapped_i]).flags() & bam_fields::flags::REVERSE) != 0;

        let rr = &mut self.records;

        // Place unmapped read at mapped read's coordinates
        bam_fields::set_ref_id(&mut rr[unmapped_i], mapped_ref_id);
        bam_fields::set_pos(&mut rr[unmapped_i], mapped_pos);

        // Set mate info on mapped read (mate is unmapped)
        bam_fields::set_mate_ref_id(&mut rr[mapped_i], mapped_ref_id);
        bam_fields::set_mate_pos(&mut rr[mapped_i], mapped_pos);
        set_mate_flags(&mut rr[mapped_i], unmapped_is_reverse, true);
        bam_fields::remove_tag(rr[mapped_i].as_mut_vec(), SamTag::MQ);
        bam_fields::remove_tag(rr[mapped_i].as_mut_vec(), SamTag::MC);
        bam_fields::set_template_length(&mut rr[mapped_i], 0);

        // Set mate info on unmapped read (mate is mapped)
        bam_fields::set_mate_ref_id(&mut rr[unmapped_i], mapped_ref_id);
        bam_fields::set_mate_pos(&mut rr[unmapped_i], mapped_pos);
        set_mate_flags(&mut rr[unmapped_i], mapped_is_reverse, false);
        let mq_val = if mapped_mapq == 255 { 255 } else { i32::from(mapped_mapq) };
        bam_fields::update_int_tag(rr[unmapped_i].as_mut_vec(), SamTag::MQ, mq_val);
        if !mapped_cigar_str.is_empty() && mapped_cigar_str != "*" {
            bam_fields::update_string_tag(
                rr[unmapped_i].as_mut_vec(),
                SamTag::MC,
                mapped_cigar_str.as_bytes(),
            );
        } else {
            bam_fields::remove_tag(rr[unmapped_i].as_mut_vec(), SamTag::MC);
        }
        bam_fields::set_template_length(&mut rr[unmapped_i], 0);
    }
}

/// Sets mate flags (`MATE_REVERSE`, `MATE_UNMAPPED`) on a raw BAM record.
fn set_mate_flags(record: &mut [u8], mate_is_reverse: bool, mate_is_unmapped: bool) {
    use crate::sort::bam_fields;
    let mut f = RawRecordView::new(record).flags();
    f &= !bam_fields::flags::MATE_REVERSE;
    if mate_is_reverse {
        f |= bam_fields::flags::MATE_REVERSE;
    }
    f &= !bam_fields::flags::MATE_UNMAPPED;
    if mate_is_unmapped {
        f |= bam_fields::flags::MATE_UNMAPPED;
    }
    bam_fields::set_flags(record, f);
}

/// Computes insert size (TLEN) from two raw BAM records.
///
/// Uses 0-based pos from BAM; converts to 1-based for the calculation.
fn compute_insert_size_raw(rec1: &[u8], rec2: &[u8]) -> i32 {
    use crate::sort::bam_fields;

    let f1 = RawRecordView::new(rec1).flags();
    let f2 = RawRecordView::new(rec2).flags();

    // If either read is unmapped, return 0
    if (f1 & bam_fields::flags::UNMAPPED) != 0 || (f2 & bam_fields::flags::UNMAPPED) != 0 {
        return 0;
    }

    // If reads are on different references, return 0
    if bam_fields::ref_id(rec1) != bam_fields::ref_id(rec2) {
        return 0;
    }

    // pos is 0-based in BAM; convert to 1-based for the calculation
    let pos1 = bam_fields::pos(rec1) + 1;
    let pos2 = bam_fields::pos(rec2) + 1;

    // alignment end (1-based inclusive) = pos_1based + ref_len - 1
    let ref_len1 = bam_fields::reference_length_from_raw_bam(rec1);
    let ref_len2 = bam_fields::reference_length_from_raw_bam(rec2);
    let end1 = pos1 + ref_len1 - 1;
    let end2 = pos2 + ref_len2 - 1;

    // 5' position: forward=start, reverse=end
    let first_5prime = if (f1 & bam_fields::flags::REVERSE) != 0 { end1 } else { pos1 };
    let second_5prime = if (f2 & bam_fields::flags::REVERSE) != 0 { end2 } else { pos2 };

    let adjustment = if second_5prime >= first_5prime { 1 } else { -1 };
    second_5prime - first_5prime + adjustment
}

/// Determines the pair orientation for a paired read using raw BAM bytes.
///
/// This matches htsjdk's `SamPairUtil.getPairOrientation()` algorithm exactly.
///
/// # Arguments
/// * `record` - Raw BAM bytes of a paired, mapped read with a mapped mate on the same reference
///
/// # Returns
/// The pair orientation (FR, RF, or Tandem)
#[allow(clippy::cast_possible_wrap)]
fn get_pair_orientation_raw(record: &[u8]) -> PairOrientation {
    use crate::sort::bam_fields;
    let f = RawRecordView::new(record).flags();
    let is_reverse = (f & bam_fields::flags::REVERSE) != 0;
    let mate_reverse = (f & bam_fields::flags::MATE_REVERSE) != 0;

    // Same strand = TANDEM
    if is_reverse == mate_reverse {
        return PairOrientation::Tandem;
    }

    // FR vs RF using htsjdk's logic:
    // positiveStrandFivePrimePos = readIsOnReverseStrand ? mateStart : alignmentStart
    // negativeStrandFivePrimePos = readIsOnReverseStrand ? alignmentEnd : alignmentStart + insertSize
    // FR if positiveStrandFivePrimePos < negativeStrandFivePrimePos
    let alignment_start = bam_fields::pos(record) + 1; // 0-based -> 1-based
    let mate_start = bam_fields::mate_pos(record) + 1;
    let insert_size = bam_fields::template_length(record);

    let (positive_five_prime, negative_five_prime) = if is_reverse {
        let ref_len = bam_fields::reference_length_from_raw_bam(record);
        let end = alignment_start + ref_len - 1;
        (i64::from(mate_start), i64::from(end))
    } else {
        (i64::from(alignment_start), i64::from(alignment_start) + i64::from(insert_size))
    };

    if positive_five_prime < negative_five_prime {
        PairOrientation::FR
    } else {
        PairOrientation::RF
    }
}

// ============================================================================
// TemplateIterator
// ============================================================================

/// Iterator that groups consecutive raw BAM records by query name into `Template` objects.
///
/// Reads `RawRecord`s and groups them into `Template` objects using [`Template::from_records`].
/// The input must be queryname-sorted or queryname-grouped.
///
/// # Requirements
///
/// The input records MUST be queryname-sorted or grouped. Records with the same name must
/// appear consecutively.
pub struct TemplateIterator<R>
where
    R: std::io::Read,
{
    reader: fgumi_raw_bam::RawBamReader<R>,
    /// Lookahead record (peeked from the stream but not yet consumed)
    pending: Option<RawRecord>,
    /// Whether the underlying reader has been exhausted
    exhausted: bool,
}

impl<R: std::io::Read> TemplateIterator<R> {
    /// Creates a new `RawTemplateIterator` from a [`RawBamReader`].
    ///
    /// The reader must already have had its header consumed (e.g., via
    /// [`create_raw_bam_reader`]).
    pub fn new(reader: fgumi_raw_bam::RawBamReader<R>) -> Self {
        TemplateIterator { reader, pending: None, exhausted: false }
    }
}

impl<R: std::io::Read> Iterator for TemplateIterator<R> {
    type Item = Result<Template>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted && self.pending.is_none() {
            return None;
        }

        let mut batch: Vec<RawRecord> = Vec::with_capacity(2);

        // Take any previously-peeked record as the first in this batch
        if let Some(peeked) = self.pending.take() {
            batch.push(peeked);
        }

        // Read records until the query name changes or we exhaust the reader
        loop {
            if self.exhausted {
                break;
            }
            let mut rec = RawRecord::new();
            match self.reader.read_record(&mut rec) {
                Ok(0) => {
                    self.exhausted = true;
                    break;
                }
                Ok(_) => {
                    let name = fgumi_raw_bam::read_name(&rec).to_vec();
                    if batch.is_empty() {
                        batch.push(rec);
                    } else {
                        let first_name = fgumi_raw_bam::read_name(&batch[0]).to_vec();
                        if name == first_name {
                            batch.push(rec);
                        } else {
                            // Different name — save for the next template
                            self.pending = Some(rec);
                            break;
                        }
                    }
                }
                Err(e) => return Some(Err(anyhow::Error::from(e))),
            }
        }

        if batch.is_empty() {
            return None;
        }

        Some(Template::from_records(batch))
    }
}

// ============================================================================
// MemoryEstimate Implementations
// ============================================================================

impl MemoryEstimate for Template {
    fn estimate_heap_size(&self) -> usize {
        // name: Vec<u8>
        let name_size = self.name.capacity();

        // records: Vec<RawRecord>
        let records_size = self.records.iter().map(RawRecord::capacity).sum::<usize>()
            + self.records.capacity() * std::mem::size_of::<RawRecord>();

        name_size + records_size
    }
}

impl MemoryEstimate for TemplateBatch {
    fn estimate_heap_size(&self) -> usize {
        self.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>()
            + self.capacity() * std::mem::size_of::<Template>()
    }
}

#[cfg(test)]
#[allow(clippy::similar_names)]
mod tests {
    use super::*;
    use fgumi_raw_bam::{RawRecord, SamBuilder as RawSamBuilder, flags as raw_flags};

    // SAM flag constants
    const FLAG_PAIRED: u16 = 0x1;
    const FLAG_READ1: u16 = 0x40;
    const FLAG_READ2: u16 = 0x80;
    const FLAG_SECONDARY: u16 = 0x100;
    const FLAG_SUPPLEMENTARY: u16 = 0x800;
    const FLAG_REVERSE: u16 = 0x10;
    const FLAG_MATE_REVERSE: u16 = 0x20;
    const FLAG_UNMAPPED: u16 = 0x4;
    const FLAG_MATE_UNMAPPED: u16 = 0x8;

    /// Create a simple test raw record with given name and flags.
    fn create_test_raw(name: &[u8], flags: u16) -> RawRecord {
        let mut b = RawSamBuilder::new();
        b.read_name(name).sequence(b"ACGT").qualities(&[30; 4]).flags(flags);
        b.build()
    }

    /// Create a mapped test raw record with position and CIGAR.
    fn create_mapped_raw(name: &[u8], flags: u16, pos: usize, mapq: u8) -> RawRecord {
        create_mapped_raw_with_tlen(name, flags, pos, mapq, 0)
    }

    /// Create a mapped test raw record with position, CIGAR, and TLEN.
    fn create_mapped_raw_with_tlen(
        name: &[u8],
        flags: u16,
        pos: usize,
        mapq: u8,
        tlen: i32,
    ) -> RawRecord {
        use fgumi_raw_bam::testutil::encode_op;
        let mut b = RawSamBuilder::new();
        b.read_name(name)
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(flags)
            .ref_id(0)
            .pos(i32::try_from(pos).expect("pos fits i32") - 1)
            .mapq(mapq)
            .template_length(tlen)
            .cigar_ops(&[encode_op(0, 100)]);
        b.build()
    }

    /// Create a mapped test raw record with full fields.
    fn create_mapped_raw_with_flags(
        name: &[u8],
        flags: u16,
        pos: usize,
        mapq: u8,
        tlen: i32,
        mate_ref_id: Option<usize>,
        mate_pos: Option<usize>,
    ) -> RawRecord {
        use fgumi_raw_bam::testutil::encode_op;
        let mut b = RawSamBuilder::new();
        b.read_name(name)
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(flags)
            .ref_id(0)
            .pos(i32::try_from(pos).expect("pos fits i32") - 1)
            .mapq(mapq)
            .template_length(tlen)
            .cigar_ops(&[encode_op(0, 100)]);
        if let Some(mref) = mate_ref_id {
            b.mate_ref_id(i32::try_from(mref).expect("mate_ref_id fits i32"));
        }
        if let Some(mpos) = mate_pos {
            b.mate_pos(i32::try_from(mpos).expect("mate_pos fits i32") - 1);
        }
        b.build()
    }

    #[test]
    fn test_template_new() {
        let template = Template::new(b"read1".to_vec());
        assert_eq!(template.name, b"read1");
        assert_eq!(template.read_count(), 0);
    }

    #[test]
    fn test_from_records_single_r1() -> Result<()> {
        let r1 = create_test_raw(b"read1", 0);
        let template = Template::from_records(vec![r1])?;
        assert_eq!(template.name, b"read1");
        assert_eq!(template.read_count(), 1);
        assert!(template.r1().is_some());
        Ok(())
    }

    #[test]
    fn test_from_records_paired() -> Result<()> {
        let r1 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r2 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2);
        let template = Template::from_records(vec![r1, r2])?;
        assert_eq!(template.read_count(), 2);
        Ok(())
    }

    #[test]
    fn test_from_records_supplementary() -> Result<()> {
        let r1 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1);
        let supp = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY);
        let template = Template::from_records(vec![r1, supp])?;
        assert_eq!(template.read_count(), 2);
        assert!(template.r1().is_some());
        Ok(())
    }

    #[test]
    fn test_from_records_error_name_mismatch() {
        let r1 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1);
        let supp = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY);
        let r2_wrong = create_test_raw(b"read2", FLAG_PAIRED | FLAG_READ2);
        let result = Template::from_records(vec![r1, supp, r2_wrong]);
        assert!(result.is_err());
    }

    #[test]
    fn test_from_records_error_multiple_r1() {
        let r1a = create_test_raw(b"read1", 0);
        let r1b = create_test_raw(b"read1", 0);
        let result = Template::from_records(vec![r1a, r1b]);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Multiple non-secondary"));
    }

    #[test]
    fn test_from_records_error_multiple_r2() {
        let r2a = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2);
        let r2b = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2);
        let result = Template::from_records(vec![r2a, r2b]);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Multiple non-secondary"));
    }

    #[test]
    fn test_from_records_error_empty() {
        let result = Template::from_records(vec![]);
        assert!(result.is_err());
    }

    #[test]
    fn test_template_from_records_basic() -> Result<()> {
        let r1 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r2 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2);
        let template = Template::from_records(vec![r1, r2])?;
        assert_eq!(template.read_count(), 2);
        assert!(template.r1().is_some());
        assert!(template.r2().is_some());
        Ok(())
    }

    #[test]
    fn test_primary_reads_iter() -> Result<()> {
        let r1 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r2 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2);
        let template = Template::from_records(vec![r1, r2])?;
        let primaries: Vec<_> = template.primary_reads().collect();
        assert_eq!(primaries.len(), 2);
        Ok(())
    }

    #[test]
    fn test_into_records_consumes() -> Result<()> {
        let r1 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r2 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2);
        let template = Template::from_records(vec![r1, r2])?;
        assert!(template.r1().is_some());
        assert!(template.r2().is_some());
        let recs = template.into_records();
        assert_eq!(recs.len(), 2);
        Ok(())
    }

    #[test]
    fn test_read_count_with_supplementals() -> Result<()> {
        let r1 = create_test_raw(b"read1", 0);
        let supp = create_test_raw(b"read1", FLAG_SUPPLEMENTARY);
        let sec = create_test_raw(b"read1", FLAG_SECONDARY);
        let template = Template::from_records(vec![r1, supp, sec])?;
        assert_eq!(template.read_count(), 3);
        Ok(())
    }

    #[test]
    fn test_r1_supplementals_range() -> Result<()> {
        let r1 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1);
        let supp1 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY);
        let supp2 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY);
        let template = Template::from_records(vec![r1, supp1, supp2])?;
        // r1_supplementals is a range (start, end) into records
        assert!(template.r1_supplementals.is_some());
        let (start, end) = template.r1_supplementals.unwrap();
        assert_eq!(end - start, 2);
        Ok(())
    }

    #[test]
    fn test_r2_supplementals_range() -> Result<()> {
        let r1 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r2 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2);
        let supp = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY);
        let template = Template::from_records(vec![r1, r2, supp])?;
        assert!(template.r2_supplementals.is_some());
        let (start, end) = template.r2_supplementals.unwrap();
        assert_eq!(end - start, 1);
        Ok(())
    }

    #[test]
    fn test_r1_secondaries_range() -> Result<()> {
        let r1 = create_test_raw(b"read1", 0);
        let sec = create_test_raw(b"read1", FLAG_SECONDARY);
        let template = Template::from_records(vec![r1, sec])?;
        assert!(template.r1_secondaries.is_some());
        let (start, end) = template.r1_secondaries.unwrap();
        assert_eq!(end - start, 1);
        Ok(())
    }

    #[test]
    fn test_r2_secondaries_range() -> Result<()> {
        let r2a = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2);
        let r2b = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SECONDARY);
        let r2c = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SECONDARY);
        let template = Template::from_records(vec![r2a, r2b, r2c])?;
        assert!(template.r2_secondaries.is_some());
        let (start, end) = template.r2_secondaries.unwrap();
        assert_eq!(end - start, 2);
        Ok(())
    }

    /// Helper to create a mapped record with position and CIGAR for `fix_mate_info` tests.
    /// Delegates to the raw helper.
    fn create_mapped_record(name: &[u8], flags: u16, pos: usize, mapq: u8) -> RawRecord {
        create_mapped_raw(name, flags, pos, mapq)
    }

    /// Helper to create a mapped record with position, CIGAR, and TLEN for `fix_mate_info` tests.
    /// Delegates to the raw helper.
    fn create_mapped_record_with_tlen(
        name: &[u8],
        flags: u16,
        pos: usize,
        mapq: u8,
        tlen: i32,
    ) -> RawRecord {
        create_mapped_raw_with_tlen(name, flags, pos, mapq, tlen)
    }

    /// Tests that `fix_mate_info` correctly sets ms tag when AS is stored as `Int8`.
    /// This was a bug where only `Int32` AS values were recognized.
    #[test]
    fn test_fix_mate_info_sets_ms_tag_with_int8_as() -> Result<()> {
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ1)
            .ref_id(0)
            .pos(99)
            .mapq(30)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .add_int_tag(SamTag::AS, 55);
        let r1 = b.build();

        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ2)
            .ref_id(0)
            .pos(199)
            .mapq(40)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .add_int_tag(SamTag::AS, 44);
        let r2 = b.build();

        let mut template = Template::from_records(vec![r1, r2])?;
        template.fix_mate_info()?;

        // R1 should have ms tag with R2's AS value (44)
        assert_eq!(
            template.records()[0].tags().find_int(SamTag::MS),
            Some(44),
            "R1 ms should be 44"
        );
        // R2 should have ms tag with R1's AS value (55)
        assert_eq!(
            template.records()[1].tags().find_int(SamTag::MS),
            Some(55),
            "R2 ms should be 55"
        );

        Ok(())
    }

    /// Tests that `fix_mate_info` correctly sets ms tag when AS is stored as `UInt8`
    #[test]
    fn test_fix_mate_info_sets_ms_tag_with_uint8_as() -> Result<()> {
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ1)
            .ref_id(0)
            .pos(99)
            .mapq(30)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .add_int_tag(SamTag::AS, 77);
        let r1 = b.build();

        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ2)
            .ref_id(0)
            .pos(199)
            .mapq(40)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .add_int_tag(SamTag::AS, 88);
        let r2 = b.build();

        let mut template = Template::from_records(vec![r1, r2])?;
        template.fix_mate_info()?;

        assert_eq!(
            template.records()[0].tags().find_int(SamTag::MS),
            Some(88),
            "R1 ms should be 88"
        );
        assert_eq!(
            template.records()[1].tags().find_int(SamTag::MS),
            Some(77),
            "R2 ms should be 77"
        );

        Ok(())
    }

    /// Tests that `fix_mate_info` correctly sets ms tag when AS is stored as `Int16`
    #[test]
    fn test_fix_mate_info_sets_ms_tag_with_int16_as() -> Result<()> {
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ1)
            .ref_id(0)
            .pos(99)
            .mapq(30)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .add_int_tag(SamTag::AS, 1000);
        let r1 = b.build();

        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ2)
            .ref_id(0)
            .pos(199)
            .mapq(40)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .add_int_tag(SamTag::AS, 2000);
        let r2 = b.build();

        let mut template = Template::from_records(vec![r1, r2])?;
        template.fix_mate_info()?;

        assert_eq!(
            template.records()[0].tags().find_int(SamTag::MS),
            Some(2000),
            "R1 ms should be 2000"
        );
        assert_eq!(
            template.records()[1].tags().find_int(SamTag::MS),
            Some(1000),
            "R2 ms should be 1000"
        );

        Ok(())
    }

    /// Tests that `fix_mate_info` correctly sets ms tag when AS is stored as `Int32`
    #[test]
    fn test_fix_mate_info_sets_ms_tag_with_int32_as() -> Result<()> {
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ1)
            .ref_id(0)
            .pos(99)
            .mapq(30)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .add_int_tag(SamTag::AS, 100_000);
        let r1 = b.build();

        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ2)
            .ref_id(0)
            .pos(199)
            .mapq(40)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .add_int_tag(SamTag::AS, 200_000);
        let r2 = b.build();

        let mut template = Template::from_records(vec![r1, r2])?;
        template.fix_mate_info()?;

        assert_eq!(
            template.records()[0].tags().find_int(SamTag::MS),
            Some(200_000),
            "R1 ms should be 200_000"
        );
        assert_eq!(
            template.records()[1].tags().find_int(SamTag::MS),
            Some(100_000),
            "R2 ms should be 100_000"
        );

        Ok(())
    }

    /// Tests that `fix_mate_info` does not set ms tag when AS is missing
    #[test]
    fn test_fix_mate_info_no_ms_tag_without_as() -> Result<()> {
        let r1 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30);
        let r2 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40);

        let mut template = Template::from_records(vec![r1, r2])?;
        template.fix_mate_info()?;

        // Neither record should have ms tag since AS is missing
        assert!(
            template.records()[0].tags().find_int(SamTag::MS).is_none(),
            "R1 should not have ms tag"
        );
        assert!(
            template.records()[1].tags().find_int(SamTag::MS).is_none(),
            "R2 should not have ms tag"
        );

        Ok(())
    }

    /// Tests that `fix_mate_info` sets ms tag on supplementary alignments
    #[test]
    fn test_fix_mate_info_sets_ms_tag_on_supplementals() -> Result<()> {
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ1)
            .ref_id(0)
            .pos(99)
            .mapq(30)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .add_int_tag(SamTag::AS, 55);
        let r1 = b.build();

        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ2)
            .ref_id(0)
            .pos(199)
            .mapq(40)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .add_int_tag(SamTag::AS, 44);
        let r2 = b.build();

        let r1_supp =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY, 300, 20);

        let mut template = Template::from_records(vec![r1, r2, r1_supp])?;
        template.fix_mate_info()?;

        // After from_records: R1[0], R2[1], R1_supp[2]
        assert_eq!(
            template.records()[2].tags().find_int(SamTag::MS),
            Some(44),
            "R1 supplementary ms should be R2's AS value (44)"
        );

        Ok(())
    }

    /// Tests that `fix_mate_info` sets TLEN on R1 supplementary alignments.
    /// TLEN should be set to negative of mate primary's (R2) TLEN.
    #[test]
    fn test_fix_mate_info_sets_tlen_on_r1_supplementals() -> Result<()> {
        let r1 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30, 200);
        let r2 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40, -200);
        let r1_supp = create_mapped_record_with_tlen(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY,
            300,
            20,
            0,
        );

        let mut template = Template::from_records(vec![r1, r2, r1_supp])?;
        template.fix_mate_info()?;

        // R1(pos=100,forward) and R2(pos=200,forward) → insert_size = 101
        // R1 supplementary TLEN = -(-101) = 101
        assert_eq!(
            template.records()[2].template_length(),
            101,
            "R1 supplementary TLEN should be negative of R2's TLEN"
        );

        Ok(())
    }

    /// Tests that `fix_mate_info` sets TLEN on R2 supplementary alignments.
    /// TLEN should be set to negative of mate primary's (R1) TLEN.
    #[test]
    fn test_fix_mate_info_sets_tlen_on_r2_supplementals() -> Result<()> {
        let r1 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30, 300);
        let r2 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40, -300);
        let r2_supp = create_mapped_record_with_tlen(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY,
            400,
            25,
            0,
        );

        let mut template = Template::from_records(vec![r1, r2, r2_supp])?;
        template.fix_mate_info()?;

        // R1(pos=100,forward) and R2(pos=200,forward) → insert_size = 101
        // R2 supplementary TLEN = -(101) = -101
        assert_eq!(
            template.records()[2].template_length(),
            -101,
            "R2 supplementary TLEN should be negative of R1's TLEN"
        );

        Ok(())
    }

    /// Tests that `fix_mate_info` handles multiple supplementary alignments
    #[test]
    fn test_fix_mate_info_sets_tlen_on_multiple_supplementals() -> Result<()> {
        let r1 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30, 500);
        let r2 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40, -500);
        let r1_supp1 = create_mapped_record_with_tlen(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY,
            300,
            20,
            0,
        );
        let r1_supp2 = create_mapped_record_with_tlen(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY,
            400,
            15,
            0,
        );
        let r2_supp1 = create_mapped_record_with_tlen(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY,
            500,
            25,
            0,
        );
        let r2_supp2 = create_mapped_record_with_tlen(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY,
            600,
            10,
            0,
        );

        // from_records orders: R1[0], R2[1], r1_supps reversed [2,3], r2_supps reversed [4,5]
        let mut template =
            Template::from_records(vec![r1, r2, r1_supp1, r1_supp2, r2_supp1, r2_supp2])?;
        template.fix_mate_info()?;

        // R1(pos=100,forward) and R2(pos=200,forward) → insert_size = 101
        // R1 supplementaries should have TLEN = -(-101) = 101
        assert_eq!(template.records()[2].template_length(), 101, "R1 supp1 TLEN");
        assert_eq!(template.records()[3].template_length(), 101, "R1 supp2 TLEN");

        // R2 supplementaries should have TLEN = -(101) = -101
        assert_eq!(template.records()[4].template_length(), -101, "R2 supp1 TLEN");
        assert_eq!(template.records()[5].template_length(), -101, "R2 supp2 TLEN");

        Ok(())
    }

    /// Tests that `fix_mate_info` recalculates TLEN for supplementaries based on primary positions
    #[test]
    fn test_fix_mate_info_tlen_recalculated() -> Result<()> {
        let r1 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30, 0);
        let r2 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40, 0);
        let r1_supp = create_mapped_record_with_tlen(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY,
            300,
            20,
            999,
        );

        let mut template = Template::from_records(vec![r1, r2, r1_supp])?;
        template.fix_mate_info()?;

        // R1(pos=100,forward) and R2(pos=200,forward) → insert_size = 101
        assert_eq!(
            template.records()[2].template_length(),
            101,
            "R1 supplementary TLEN should be recalculated from positions"
        );

        Ok(())
    }

    /// Tests that `fix_mate_info` doesn't affect supplementaries when mate primary is missing
    #[test]
    fn test_fix_mate_info_no_mate_primary() -> Result<()> {
        let r1 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30, 200);
        let r1_supp = create_mapped_record_with_tlen(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY,
            300,
            20,
            777,
        );

        let mut template = Template::from_records(vec![r1, r1_supp])?;
        template.fix_mate_info()?;

        // No R2 primary, so R1 supplementary TLEN should remain unchanged
        assert_eq!(
            template.records()[1].template_length(),
            777,
            "R1 supplementary TLEN should be unchanged without R2 primary"
        );

        Ok(())
    }

    /// Tests record ordering in Template matches fgbio's ordering.
    /// Order should be: R1, R2, `r1_supps`, `r2_supps`, `r1_secondaries`, `r2_secondaries`
    #[test]
    fn test_template_record_ordering() -> Result<()> {
        let r1 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r2 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2);
        let r1_supp1 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY);
        let r1_supp2 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY);
        let r2_supp = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY);
        let r1_sec = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SECONDARY);
        let r2_sec1 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SECONDARY);
        let r2_sec2 = create_test_raw(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SECONDARY);

        // Pass in scrambled order; from_records sorts into canonical order
        let template = Template::from_records(vec![
            r2_sec1, r1_supp2, r1, r2_supp, r1_sec, r2, r2_sec2, r1_supp1,
        ])?;

        // Verify index ranges
        assert_eq!(template.r1, Some((0, 1)), "R1 should be at index 0");
        assert_eq!(template.r2, Some((1, 2)), "R2 should be at index 1");
        assert_eq!(template.r1_supplementals, Some((2, 4)), "R1 supps should be at indices 2-3");
        assert_eq!(template.r2_supplementals, Some((4, 5)), "R2 supps should be at index 4");
        assert_eq!(template.r1_secondaries, Some((5, 6)), "R1 secs should be at index 5");
        assert_eq!(template.r2_secondaries, Some((6, 8)), "R2 secs should be at indices 6-7");

        // Verify record count
        assert_eq!(template.read_count(), 8);

        // Verify flags match expected ordering
        let recs = template.records();
        assert!(recs[0].is_first_segment() && !recs[0].is_supplementary(), "R1 primary");
        assert!(recs[1].is_last_segment() && !recs[1].is_supplementary(), "R2 primary");
        assert!(recs[2].is_first_segment() && recs[2].is_supplementary(), "R1 supp[0]");
        assert!(recs[3].is_first_segment() && recs[3].is_supplementary(), "R1 supp[1]");
        assert!(recs[4].is_last_segment() && recs[4].is_supplementary(), "R2 supp");
        assert!(recs[5].is_first_segment() && recs[5].is_secondary(), "R1 sec");
        assert!(recs[6].is_last_segment() && recs[6].is_secondary(), "R2 sec[0]");
        assert!(recs[7].is_last_segment() && recs[7].is_secondary(), "R2 sec[1]");

        Ok(())
    }

    /// Tests that records within supplementary/secondary groups are in reverse input order.
    /// This matches fgbio's behavior which uses Scala List prepend (r :: list).
    #[test]
    fn test_record_ordering_within_groups_is_reversed() -> Result<()> {
        let r1 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30);
        let r2 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 30);
        // R1 supplementaries input order: pos 1000, 2000, 3000
        let r1_supp1 =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY, 1000, 20);
        let r1_supp2 =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY, 2000, 20);
        let r1_supp3 =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY, 3000, 20);
        // R2 supplementaries input order: pos 4000, 5000
        let r2_supp1 =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY, 4000, 20);
        let r2_supp2 =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY, 5000, 20);
        // R1 secondaries input order: pos 6000, 7000
        let r1_sec1 =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SECONDARY, 6000, 10);
        let r1_sec2 =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SECONDARY, 7000, 10);

        let template = Template::from_records(vec![
            r1, r2, r1_supp1, r1_supp2, r1_supp3, r2_supp1, r2_supp2, r1_sec1, r1_sec2,
        ])?;

        let recs = template.records();
        // Primaries first
        assert_eq!(recs[0].pos() + 1, 100, "R1 primary");
        assert_eq!(recs[1].pos() + 1, 200, "R2 primary");

        // R1 supplementaries reversed: 3000, 2000, 1000
        assert_eq!(recs[2].pos() + 1, 3000, "R1 supp in reverse order");
        assert_eq!(recs[3].pos() + 1, 2000, "R1 supp in reverse order");
        assert_eq!(recs[4].pos() + 1, 1000, "R1 supp in reverse order");

        // R2 supplementaries reversed: 5000, 4000
        assert_eq!(recs[5].pos() + 1, 5000, "R2 supp in reverse order");
        assert_eq!(recs[6].pos() + 1, 4000, "R2 supp in reverse order");

        // R1 secondaries reversed: 7000, 6000
        assert_eq!(recs[7].pos() + 1, 7000, "R1 sec in reverse order");
        assert_eq!(recs[8].pos() + 1, 6000, "R1 sec in reverse order");

        Ok(())
    }

    /// Helper to create a mapped record with specific flags including reverse strand.
    /// Delegates to the raw helper.
    fn create_mapped_record_with_flags(
        name: &[u8],
        flags: u16,
        pos: usize,
        mapq: u8,
        tlen: i32,
        mate_ref_id: Option<usize>,
        mate_pos: Option<usize>,
    ) -> RawRecord {
        create_mapped_raw_with_flags(name, flags, pos, mapq, tlen, mate_ref_id, mate_pos)
    }

    /// Tests that `fix_mate_info` sets mate info on R1 supplementary alignments correctly.
    /// The mate is R2 primary, so supplementary should get R2's info.
    #[test]
    fn test_fix_mate_info_sets_full_mate_info_on_r1_supplementals() -> Result<()> {
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ1)
            .ref_id(0)
            .pos(99)
            .mapq(30)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .template_length(200)
            .mate_ref_id(0)
            .mate_pos(199)
            .add_int_tag(SamTag::AS, 100);
        let r1 = b.build();

        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE)
            .ref_id(0)
            .pos(199)
            .mapq(40)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .template_length(-200)
            .mate_ref_id(0)
            .mate_pos(99)
            .add_int_tag(SamTag::AS, 150);
        let r2 = b.build();

        // R1 supp has wrong mate info (will be corrected by fix_mate_info)
        let r1_supp = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY | FLAG_MATE_REVERSE,
            500,
            25,
            0,
            Some(0),
            Some(500),
        );

        let mut template = Template::from_records(vec![r1, r2, r1_supp])?;
        template.fix_mate_info()?;

        let supp = &template.records()[2];

        // Mate position should be R2's position (200 = 0-based 199 + 1)
        assert_eq!(supp.mate_pos() + 1, 200, "R1 supp should have mate pos from R2");
        // Mate ref should be 0 (R2's reference)
        assert_eq!(supp.mate_ref_id(), 0, "R1 supp should have mate ref from R2");
        // Mate reverse should be set (R2 is reverse)
        assert!(supp.is_mate_reverse(), "R1 supp should have mate_reverse since R2 is reverse");
        // Mate unmapped should NOT be set
        assert!(!supp.is_mate_unmapped(), "R1 supp should NOT have mate_unmapped");
        // TLEN should be -(-200) = 200
        assert_eq!(supp.template_length(), 200, "R1 supp TLEN should be -(-200) = 200");
        // MQ tag should be R2's mapq (40)
        assert_eq!(supp.tags().find_int(SamTag::MQ), Some(40), "R1 supp MQ should be R2's mapq");
        // MC tag should be present
        assert!(supp.tags().find_string(SamTag::MC).is_some(), "R1 supp should have MC tag");
        // ms tag should be R2's AS value (150)
        assert_eq!(supp.tags().find_int(SamTag::MS), Some(150), "R1 supp ms should be R2's AS");

        Ok(())
    }

    /// Tests that `fix_mate_info` sets mate info on R2 supplementary alignments correctly.
    /// The mate is R1 primary, so supplementary should get R1's info.
    #[test]
    fn test_fix_mate_info_sets_full_mate_info_on_r2_supplementals() -> Result<()> {
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ1)
            .ref_id(0)
            .pos(99)
            .mapq(35)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .template_length(300)
            .mate_ref_id(0)
            .mate_pos(199)
            .add_int_tag(SamTag::AS, 120);
        let r1 = b.build();

        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ2)
            .ref_id(0)
            .pos(199)
            .mapq(45)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .template_length(-300)
            .mate_ref_id(0)
            .mate_pos(99)
            .add_int_tag(SamTag::AS, 180);
        let r2 = b.build();

        let r2_supp = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY,
            600,
            20,
            0,
            Some(0),
            Some(600),
        );

        let mut template = Template::from_records(vec![r1, r2, r2_supp])?;
        template.fix_mate_info()?;

        let supp = &template.records()[2];

        // Mate position should be R1's position (100)
        assert_eq!(supp.mate_pos() + 1, 100, "R2 supp should have mate pos from R1");
        // Mate reverse should NOT be set (R1 is forward)
        assert!(
            !supp.is_mate_reverse(),
            "R2 supp should NOT have mate_reverse since R1 is forward"
        );
        // TLEN = -(101) = -101
        // R1(pos=100,forward) and R2(pos=200,forward) → insert_size = 101
        assert_eq!(supp.template_length(), -101, "R2 supp TLEN should be -101");
        // MQ tag should be R1's mapq (35)
        assert_eq!(supp.tags().find_int(SamTag::MQ), Some(35), "R2 supp MQ should be R1's mapq");
        // MC tag should be present
        assert!(supp.tags().find_string(SamTag::MC).is_some(), "R2 supp should have MC tag");
        // ms tag should be R1's AS value (120)
        assert_eq!(supp.tags().find_int(SamTag::MS), Some(120), "R2 supp ms should be R1's AS");

        Ok(())
    }

    /// Tests that `fix_mate_info` correctly handles unmapped mate.
    /// When mate is unmapped, `mate_unmapped` flag should be set and MC tag should not be set.
    #[test]
    fn test_fix_mate_info_supplemental_with_unmapped_mate() -> Result<()> {
        // R1 primary knows its mate is unmapped
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_UNMAPPED,
            100,
            30,
            0,
            Some(0),
            Some(100),
        );

        // R2 is unmapped (no CIGAR, placed at R1's position)
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ2 | FLAG_UNMAPPED)
            .ref_id(0)
            .pos(99)
            .mapq(0)
            .template_length(0)
            .mate_ref_id(0)
            .mate_pos(99);
        let r2 = b.build();

        // R1 supp originally has incorrect MATE_REVERSE and TLEN
        let r1_supp = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY | FLAG_MATE_REVERSE,
            500,
            25,
            999,
            Some(0),
            Some(500),
        );

        let mut template = Template::from_records(vec![r1, r2, r1_supp])?;
        template.fix_mate_info()?;

        let supp = &template.records()[2];

        assert!(supp.is_mate_unmapped(), "R1 supp should have mate_unmapped since R2 is unmapped");
        assert!(
            !supp.is_mate_reverse(),
            "R1 supp should NOT have mate_reverse when R2 is unmapped"
        );
        assert!(supp.tags().find_string(SamTag::MC).is_none(), "R1 supp should NOT have MC tag");
        assert_eq!(supp.template_length(), 0, "TLEN should be 0 when mate is unmapped");

        Ok(())
    }

    /// Tests that `fix_mate_info` correctly clears `mate_reverse` flag when mate is on forward strand.
    #[test]
    fn test_fix_mate_info_clears_incorrect_mate_reverse_flag() -> Result<()> {
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1,
            100,
            30,
            200,
            Some(0),
            Some(200),
        );
        let r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2, // no FLAG_REVERSE
            200,
            40,
            -200,
            Some(0),
            Some(100),
        );
        // R1 supp incorrectly has MATE_REVERSE
        let r1_supp = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY | FLAG_MATE_REVERSE,
            500,
            25,
            0,
            Some(0),
            Some(500),
        );

        let mut template = Template::from_records(vec![r1, r2, r1_supp])?;

        // Before fix, supplementary has MATE_REVERSE set incorrectly
        assert!(
            template.records()[2].is_mate_reverse(),
            "Before fix: supp incorrectly has mate_reverse"
        );

        template.fix_mate_info()?;

        // After fix, MATE_REVERSE should be cleared since R2 is forward
        assert!(
            !template.records()[2].is_mate_reverse(),
            "After fix: supp should NOT have mate_reverse"
        );

        Ok(())
    }

    /// Tests that `fix_mate_info` handles multiple supplementary alignments.
    #[test]
    fn test_fix_mate_info_multiple_supplementals() -> Result<()> {
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1,
            100,
            30,
            400,
            Some(0),
            Some(300),
        );

        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE)
            .ref_id(0)
            .pos(299)
            .mapq(50)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .template_length(-400)
            .mate_ref_id(0)
            .mate_pos(99)
            .add_int_tag(SamTag::AS, 200);
        let r2 = b.build();

        let r1_supp1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY,
            500,
            20,
            0,
            Some(0),
            Some(500),
        );
        let r1_supp2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY,
            700,
            15,
            0,
            Some(0),
            Some(700),
        );

        let mut template = Template::from_records(vec![r1, r2, r1_supp1, r1_supp2])?;
        template.fix_mate_info()?;

        // R1(pos=100,forward) and R2(pos=300,100M,reverse) → 5' are 100 and 399 → insert_size = 300
        for i in 2..=3 {
            let supp = &template.records()[i];
            assert_eq!(supp.mate_pos() + 1, 300, "Supp {i} should have mate pos 300");
            assert!(supp.is_mate_reverse(), "Supp {i} should have mate_reverse");
            assert_eq!(supp.template_length(), 300, "Supp {i} TLEN should be 300");
        }

        Ok(())
    }

    // ============================================================================
    // Tests for setMateInfo (both mapped case) - ported from htsjdk
    // ============================================================================

    /// Test `fix_mate_info` correctly sets TLEN for normal innie pair
    #[test]
    fn test_fix_mate_info_sets_tlen_normal_innie() -> Result<()> {
        // R1: pos=1, 50M, forward; R2: pos=500, 100M, reverse
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 50])
            .qualities(&[30u8; 50])
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(0)
            .mapq(30)
            .cigar_ops(&[50u32 << 4]);
        let r1 = b.build();

        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE)
            .ref_id(0)
            .pos(499)
            .mapq(30)
            .cigar_ops(&[100u32 << 4]);
        let r2 = b.build();

        let mut template = Template::from_records(vec![r1, r2])?;
        template.fix_mate_info()?;

        // R1 5' = 1 (forward start), R2 5' = 599 (reverse end = 500+100-1)
        // insert_size = 599 - 1 + 1 = 599
        assert_eq!(template.records()[0].template_length(), 599, "R1 TLEN");
        assert_eq!(template.records()[1].template_length(), -599, "R2 TLEN");

        Ok(())
    }

    /// Test `fix_mate_info` correctly sets mate CIGAR (MC tag)
    #[test]
    fn test_fix_mate_info_sets_mate_cigar() -> Result<()> {
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 50])
            .qualities(&[30u8; 50])
            .flags(raw_flags::PAIRED | raw_flags::FIRST_SEGMENT)
            .ref_id(0)
            .pos(0)
            .mapq(30)
            .cigar_ops(&[50u32 << 4]); // 50M
        let r1 = b.build();

        // R2 CIGAR: 25M5I20M
        let cigar = vec![25u32 << 4, (5u32 << 4) | 1, 20u32 << 4];
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 50])
            .qualities(&[30u8; 50])
            .flags(raw_flags::PAIRED | raw_flags::LAST_SEGMENT | raw_flags::REVERSE)
            .ref_id(0)
            .pos(499)
            .mapq(30)
            .cigar_ops(&cigar);
        let r2 = b.build();

        let mut template = Template::from_records(vec![r1, r2])?;
        template.fix_mate_info()?;

        // R1's MC should be R2's CIGAR (25M5I20M)
        let r1_mc = template.records()[0].tags().find_string(SamTag::MC);
        assert!(r1_mc.is_some(), "R1 should have MC tag");
        assert_eq!(r1_mc.unwrap(), b"25M5I20M", "R1 MC should be R2's CIGAR");

        // R2's MC should be R1's CIGAR (50M)
        let r2_mc = template.records()[1].tags().find_string(SamTag::MC);
        assert!(r2_mc.is_some(), "R2 should have MC tag");
        assert_eq!(r2_mc.unwrap(), b"50M", "R2 MC should be R1's CIGAR");

        Ok(())
    }

    /// Test `fix_mate_info` correctly handles both unmapped case
    #[test]
    fn test_fix_mate_info_both_unmapped() -> Result<()> {
        // R1: unmapped, placed at ref 0, has MQ tag
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGT")
            .qualities(&[30u8; 4])
            .flags(FLAG_PAIRED | FLAG_READ1 | FLAG_UNMAPPED)
            .ref_id(0)
            .pos(0)
            .mapq(0)
            .add_int_tag(SamTag::MQ, 30);
        let r1 = b.build();

        // R2: unmapped, placed at ref 0, has MC tag
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGT")
            .qualities(&[30u8; 4])
            .flags(FLAG_PAIRED | FLAG_READ2 | FLAG_UNMAPPED)
            .ref_id(0)
            .pos(0)
            .mapq(0)
            .add_string_tag(SamTag::MC, b"100M");
        let r2 = b.build();

        let mut template = Template::from_records(vec![r1, r2])?;
        template.fix_mate_info()?;

        // Both should have ref_id = -1 (unmapped)
        assert_eq!(template.records()[0].ref_id(), -1, "R1 ref should be cleared");
        assert_eq!(template.records()[1].ref_id(), -1, "R2 ref should be cleared");

        // Both should have mate_unmapped flag set
        assert!(template.records()[0].is_mate_unmapped(), "R1 mate_unmapped");
        assert!(template.records()[1].is_mate_unmapped(), "R2 mate_unmapped");

        // MQ and MC tags should be removed
        assert!(
            template.records()[0].tags().find_int(SamTag::MQ).is_none(),
            "R1 MQ should be removed"
        );
        assert!(
            template.records()[1].tags().find_string(SamTag::MC).is_none(),
            "R2 MC should be removed"
        );

        // TLEN should be 0
        assert_eq!(template.records()[0].template_length(), 0, "R1 TLEN should be 0");
        assert_eq!(template.records()[1].template_length(), 0, "R2 TLEN should be 0");

        Ok(())
    }

    /// Test `fix_mate_info` correctly handles one mapped, one unmapped case
    #[test]
    fn test_fix_mate_info_one_unmapped() -> Result<()> {
        // R1 is mapped
        let r1 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30);

        // R2 is unmapped (no CIGAR, placed at some position)
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(b"ACGT")
            .qualities(&[30u8; 4])
            .flags(FLAG_PAIRED | FLAG_READ2 | FLAG_UNMAPPED)
            .mapq(0);
        let r2 = b.build();

        let mut template = Template::from_records(vec![r1, r2])?;
        template.fix_mate_info()?;

        // Unmapped R2 should be placed at R1's ref_id (0) and position (100, i.e. 0-based 99)
        assert_eq!(template.records()[1].ref_id(), 0, "R2 should be placed at R1's reference");
        assert_eq!(template.records()[1].pos() + 1, 100, "R2 should be placed at R1's position");

        // R1 (mapped) should NOT have MQ tag (mate is unmapped)
        assert!(
            template.records()[0].tags().find_int(SamTag::MQ).is_none(),
            "R1 should NOT have MQ"
        );

        // R2 (unmapped) SHOULD have MQ tag (mate is mapped)
        assert!(template.records()[1].tags().find_int(SamTag::MQ).is_some(), "R2 should have MQ");

        // R1 should have mate_unmapped flag set
        assert!(template.records()[0].is_mate_unmapped(), "R1 mate_unmapped should be set");

        // R2 should NOT have mate_unmapped flag set (R1 is mapped)
        assert!(!template.records()[1].is_mate_unmapped(), "R2 mate_unmapped should NOT be set");

        // TLEN should be 0 for both
        assert_eq!(template.records()[0].template_length(), 0, "R1 TLEN");
        assert_eq!(template.records()[1].template_length(), 0, "R2 TLEN");

        // MC tag: R1 should not have it (mate unmapped), R2 should have it (mate mapped)
        assert!(template.records()[0].tags().find_string(SamTag::MC).is_none(), "R1 NOT have MC");
        assert!(
            template.records()[1].tags().find_string(SamTag::MC).is_some(),
            "R2 should have MC"
        );

        Ok(())
    }

    // ==================== pair_orientation tests ====================

    /// Test FR pair: R1 forward at pos 100, R2 reverse at pos 200 (positive insert size)
    /// This matches htsjdk's test case for getPairOrientation
    #[test]
    fn test_pair_orientation_fr_pair() -> Result<()> {
        // R1: forward strand at pos 100, mate (R2) is reverse at pos 200
        // positive 5' = 100 (R1 start), negative 5' = 100 + 200 = 300
        // 100 < 300 => FR
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_REVERSE,
            100,
            30,
            200, // positive TLEN for FR
            Some(0),
            Some(200),
        );
        let r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE,
            200,
            30,
            -200,
            Some(0),
            Some(100),
        );

        let template = Template::from_records(vec![r1, r2])?;
        assert_eq!(template.pair_orientation(), Some(PairOrientation::FR));
        Ok(())
    }

    /// Test RF pair: R1 forward but positioned after mate's 5' end
    /// This is an "outie" pair
    #[test]
    fn test_pair_orientation_rf_pair() -> Result<()> {
        // R1: forward at pos 300, R2: reverse at pos 100
        // positive 5' = 300 (R1 start), negative 5' = 300 + (-200) = 100
        // 300 >= 100 => RF
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_REVERSE,
            300,
            30,
            -200, // negative TLEN for RF
            Some(0),
            Some(100),
        );
        let r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE,
            100,
            30,
            200,
            Some(0),
            Some(300),
        );

        let template = Template::from_records(vec![r1, r2])?;
        assert_eq!(template.pair_orientation(), Some(PairOrientation::RF));
        Ok(())
    }

    /// Test TANDEM pair: both reads on forward strand
    #[test]
    fn test_pair_orientation_tandem_both_forward() -> Result<()> {
        // Both R1 and R2 on forward strand => TANDEM
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1, // forward, mate also forward (no MATE_REVERSE)
            100,
            30,
            100,
            Some(0),
            Some(200),
        );
        let r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2, // forward
            200,
            30,
            -100,
            Some(0),
            Some(100),
        );

        let template = Template::from_records(vec![r1, r2])?;
        assert_eq!(template.pair_orientation(), Some(PairOrientation::Tandem));
        Ok(())
    }

    /// Test TANDEM pair: both reads on reverse strand
    #[test]
    fn test_pair_orientation_tandem_both_reverse() -> Result<()> {
        // Both R1 and R2 on reverse strand => TANDEM
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_REVERSE | FLAG_MATE_REVERSE,
            100,
            30,
            100,
            Some(0),
            Some(200),
        );
        let r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE | FLAG_MATE_REVERSE,
            200,
            30,
            -100,
            Some(0),
            Some(100),
        );

        let template = Template::from_records(vec![r1, r2])?;
        assert_eq!(template.pair_orientation(), Some(PairOrientation::Tandem));
        Ok(())
    }

    /// Test `pair_orientation` returns None when R1 is missing
    #[test]
    fn test_pair_orientation_no_r1() -> Result<()> {
        let r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE,
            100,
            30,
            0,
            Some(0),
            Some(100),
        );

        let template = Template::from_records(vec![r2])?;
        assert_eq!(template.pair_orientation(), None);
        Ok(())
    }

    /// Test `pair_orientation` returns None when R2 is missing
    #[test]
    fn test_pair_orientation_no_r2() -> Result<()> {
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1,
            100,
            30,
            0,
            Some(0),
            Some(100),
        );

        let template = Template::from_records(vec![r1])?;
        assert_eq!(template.pair_orientation(), None);
        Ok(())
    }

    /// Test `pair_orientation` returns None when R1 is unmapped
    #[test]
    fn test_pair_orientation_r1_unmapped() -> Result<()> {
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_UNMAPPED,
            100,
            30,
            0,
            Some(0),
            Some(100),
        );
        let r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2,
            100,
            30,
            0,
            Some(0),
            Some(100),
        );

        let template = Template::from_records(vec![r1, r2])?;
        assert_eq!(template.pair_orientation(), None);
        Ok(())
    }

    /// Test `pair_orientation` returns None when R2 is unmapped
    #[test]
    fn test_pair_orientation_r2_unmapped() -> Result<()> {
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1,
            100,
            30,
            0,
            Some(0),
            Some(100),
        );
        let r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_UNMAPPED,
            100,
            30,
            0,
            Some(0),
            Some(100),
        );

        let template = Template::from_records(vec![r1, r2])?;
        assert_eq!(template.pair_orientation(), None);
        Ok(())
    }

    /// Test `pair_orientation` returns None when reads are on different chromosomes
    #[test]
    fn test_pair_orientation_different_chromosomes() -> Result<()> {
        // R1 on ref 0, mate_ref_id=1 (different chrom)
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_REVERSE,
            100,
            30,
            200,
            Some(0),
            Some(200),
        );
        // R2 on ref 1
        let mut b = RawSamBuilder::new();
        b.read_name(b"read1")
            .sequence(&[b'A'; 100])
            .qualities(&[30u8; 100])
            .flags(FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE)
            .ref_id(1)
            .pos(199)
            .mapq(30)
            .cigar_ops(&[fgumi_raw_bam::testutil::encode_op(0, 100)])
            .template_length(-200)
            .mate_ref_id(0)
            .mate_pos(99);
        let r2 = b.build();

        let template = Template::from_records(vec![r1, r2])?;
        assert_eq!(template.pair_orientation(), None);
        Ok(())
    }

    /// Test FR pair from reverse strand perspective (R1 is reverse, R2 is forward)
    #[test]
    fn test_pair_orientation_fr_from_reverse_read() -> Result<()> {
        // R1: reverse at pos 200, R2: forward at pos 100
        // For R1 (reverse): positive 5' = mate_start = 100, negative 5' = alignment_end = 299 (200 + 100 - 1)
        // 100 < 299 => FR
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_REVERSE, // R1 reverse, mate forward (no MATE_REVERSE)
            200,
            30,
            -200, // negative because R1 is after R2
            Some(0),
            Some(100),
        );
        let r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_MATE_REVERSE, // R2 forward, mate (R1) reverse
            100,
            30,
            200,
            Some(0),
            Some(200),
        );

        let template = Template::from_records(vec![r1, r2])?;
        assert_eq!(template.pair_orientation(), Some(PairOrientation::FR));
        Ok(())
    }

    // ========================================================================
    // from_records error-path tests
    // ========================================================================

    #[test]
    fn test_from_records_truncated_header() {
        let short = RawRecord::from(vec![0u8; 20]);
        let err = Template::from_records(vec![short]).unwrap_err();
        assert!(err.to_string().contains("too short"), "Error: {err}");
    }

    #[test]
    fn test_from_records_truncated_read_name() {
        let mut buf = vec![0u8; 34];
        buf[8] = 10; // claims 10 bytes for read name, but only 2 available
        let err = Template::from_records(vec![RawRecord::from(buf)]).unwrap_err();
        assert!(err.to_string().contains("truncated"), "Error: {err}");
    }

    // ========================================================================
    // write_with_offset tests
    // ========================================================================

    #[test]
    fn test_write_with_offset_none() {
        let mi = MoleculeId::None;
        let mut buf = String::new();
        let result = mi.write_with_offset(100, &mut buf);
        assert!(result.is_empty());
    }

    #[test]
    fn test_write_with_offset_single() {
        let mi = MoleculeId::Single(5);
        let mut buf = String::new();
        let result = mi.write_with_offset(100, &mut buf);
        assert_eq!(result, b"105");
    }

    #[test]
    fn test_write_with_offset_paired_a() {
        let mi = MoleculeId::PairedA(3);
        let mut buf = String::new();
        let result = mi.write_with_offset(10, &mut buf);
        assert_eq!(result, b"13/A");
    }

    #[test]
    fn test_write_with_offset_paired_b() {
        let mi = MoleculeId::PairedB(3);
        let mut buf = String::new();
        let result = mi.write_with_offset(10, &mut buf);
        assert_eq!(result, b"13/B");
    }

    #[test]
    fn test_write_with_offset_reuses_buffer() {
        let mut buf = String::new();
        let mi1 = MoleculeId::Single(1);
        let _ = mi1.write_with_offset(0, &mut buf);
        assert_eq!(buf, "1");

        // Reuse buffer — should clear and overwrite
        let mi2 = MoleculeId::PairedA(99);
        let result = mi2.write_with_offset(0, &mut buf);
        assert_eq!(result, b"99/A");
        assert_eq!(buf, "99/A");
    }

    #[test]
    fn test_write_with_offset_matches_to_string() {
        // Verify write_with_offset produces identical output to to_string_with_offset
        let mut buf = String::new();
        for mi in [
            MoleculeId::None,
            MoleculeId::Single(0),
            MoleculeId::Single(42),
            MoleculeId::PairedA(7),
            MoleculeId::PairedB(7),
        ] {
            let expected = mi.to_string_with_offset(100);
            let result = mi.write_with_offset(100, &mut buf);
            assert_eq!(result, expected.as_bytes(), "Mismatch for {mi:?}");
        }
    }
}
