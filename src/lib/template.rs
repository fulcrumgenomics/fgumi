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

use crate::sam::{record_utils, to_smallest_signed_int};
use crate::unified_pipeline::MemoryEstimate;
use anyhow::{Result, anyhow, bail};
use noodles::sam::alignment::record::Cigar;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
use std::fmt::Write as _;

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
    /// The records (empty in raw-byte mode)
    pub records: Vec<RecordBuf>,
    /// Raw BAM record bytes (without `block_size` prefix). Present only in raw-byte mode.
    pub raw_records: Option<Vec<Vec<u8>>>,
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

/// A builder for constructing `Template` instances.
///
/// The builder allows incremental construction of a template by adding records
/// one at a time. It automatically categorizes each record into the appropriate
/// slot (primary R1/R2, supplementary, or secondary) based on SAM flags.
///
/// # Examples
///
/// ```rust,ignore
/// let mut builder = Builder::new();
/// builder.push(r1_record)?;
/// builder.push(r2_record)?;
/// builder.push(r1_supplementary)?;
/// let template = builder.build()?;
/// ```
#[derive(Debug, Default)]
pub struct Builder {
    /// The query name (QNAME) shared by all records in this template
    name: Option<Vec<u8>>,
    /// Primary R1 read (first segment, non-secondary, non-supplementary)
    r1: Option<RecordBuf>,
    /// Primary R2 read (second segment, non-secondary, non-supplementary)
    r2: Option<RecordBuf>,
    /// Supplementary alignments for R1
    r1_supplementals: Vec<RecordBuf>,
    /// Supplementary alignments for R2
    r2_supplementals: Vec<RecordBuf>,
    /// Secondary alignments for R1
    r1_secondaries: Vec<RecordBuf>,
    /// Secondary alignments for R2
    r2_secondaries: Vec<RecordBuf>,
}

impl Builder {
    /// Creates a new empty template builder.
    ///
    /// # Returns
    ///
    /// A new `Builder` instance ready to accept records
    fn new() -> Self {
        Builder {
            name: None,
            r1: None,
            r2: None,
            r1_supplementals: Vec::new(),
            r2_supplementals: Vec::new(),
            r1_secondaries: Vec::new(),
            r2_secondaries: Vec::new(),
        }
    }

    /// Adds a record to this builder, organizing it into the appropriate category.
    ///
    /// # Arguments
    ///
    /// * `record` - The record to add
    ///
    /// # Errors
    ///
    /// Returns an error if multiple primary R1s or R2s are encountered
    ///
    /// # Panics
    ///
    /// Panics if template name is unexpectedly missing when reporting duplicate records.
    pub fn push(&mut self, record: RecordBuf) -> Result<&mut Builder> {
        // Get record name as bytes for comparison (no allocation yet)
        let record_name_bytes: &[u8] = record.name().map_or(&[], <_ as AsRef<[u8]>>::as_ref);

        if let Some(cur_name) = &self.name {
            // Compare bytes directly without allocating
            if cur_name.as_slice() != record_name_bytes {
                bail!(
                    "Template name mismatch: found '{record_name_bytes:?}', expected '{cur_name:?}'"
                );
            }
        } else {
            // Only allocate Vec when setting the name for the first time
            self.name = Some(record_name_bytes.to_vec());
        }

        let flags = record.flags();
        let is_r1 = !flags.is_segmented() || flags.is_first_segment();

        if is_r1 {
            if flags.is_secondary() {
                self.r1_secondaries.push(record);
            } else if flags.is_supplementary() {
                self.r1_supplementals.push(record);
            } else if self.r1.is_some() {
                return Err(anyhow!(
                    "Multiple non-secondary, non-supplemental R1 records for read '{:?}'",
                    self.name.as_ref().unwrap()
                ));
            } else {
                self.r1 = Some(record);
            }
        } else if flags.is_secondary() {
            self.r2_secondaries.push(record);
        } else if flags.is_supplementary() {
            self.r2_supplementals.push(record);
        } else if self.r2.is_some() {
            return Err(anyhow!(
                "Multiple non-secondary, non-supplemental R2 records for read '{:?}'",
                self.name.as_ref().unwrap()
            ));
        } else {
            self.r2 = Some(record);
        }

        Ok(self)
    }

    /// Returns the number of records currently in the builder.
    ///
    /// # Returns
    ///
    /// Total count of records across all categories
    #[must_use]
    pub fn len(&self) -> usize {
        let mut count = 0;
        count += usize::from(self.r1.is_some());
        count += usize::from(self.r2.is_some());
        count += self.r1_supplementals.len();
        count += self.r2_supplementals.len();
        count += self.r1_secondaries.len();
        count += self.r2_secondaries.len();
        count
    }

    /// Checks if the builder contains no records.
    ///
    /// # Returns
    ///
    /// `true` if the builder is empty, `false` otherwise
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Builds a new template from this builder.
    ///
    /// The builder may be re-used.
    ///
    /// # Errors
    ///
    /// If the template could not be built.
    pub fn build(&mut self) -> Result<Template> {
        let r1_end: usize = usize::from(self.r1.is_some());
        let r2_end: usize = r1_end + usize::from(self.r2.is_some());
        let r1_supplementals_end: usize = r2_end + self.r1_supplementals.len();
        let r2_supplementals_end: usize = r1_supplementals_end + self.r2_supplementals.len();
        let r1_secondaries_end: usize = r2_supplementals_end + self.r1_secondaries.len();
        let r2_secondaries_end: usize = r1_secondaries_end + self.r2_secondaries.len();

        let r1: Option<(usize, usize)> = if self.r1.is_some() { Some((0, r1_end)) } else { None };
        let r2: Option<(usize, usize)> =
            if self.r2.is_some() { Some((r1_end, r2_end)) } else { None };
        let r1_supplementals: Option<(usize, usize)> = if self.r1_supplementals.is_empty() {
            None
        } else {
            Some((r2_end, r1_supplementals_end))
        };
        let r2_supplementals: Option<(usize, usize)> = if self.r2_supplementals.is_empty() {
            None
        } else {
            Some((r1_supplementals_end, r2_supplementals_end))
        };
        let r1_secondaries: Option<(usize, usize)> = if self.r1_secondaries.is_empty() {
            None
        } else {
            Some((r2_supplementals_end, r1_secondaries_end))
        };
        let r2_secondaries: Option<(usize, usize)> = if self.r2_secondaries.is_empty() {
            None
        } else {
            Some((r1_secondaries_end, r2_secondaries_end))
        };

        let Some(name) = self.name.take() else {
            return Err(anyhow!("No records given to template builder"));
        };

        // Pre-allocate records Vec with exact capacity to avoid reallocations
        let mut records: Vec<RecordBuf> = Vec::with_capacity(self.len());
        if let Some(rec) = self.r1.take() {
            records.push(rec);
        }
        if let Some(rec) = self.r2.take() {
            records.push(rec);
        }
        // Reverse the supplementary and secondary vectors to match fgbio's behavior.
        // fgbio uses Scala List prepend (r :: list) which results in reverse input order.
        self.r1_supplementals.reverse();
        self.r2_supplementals.reverse();
        self.r1_secondaries.reverse();
        self.r2_secondaries.reverse();
        records.append(&mut self.r1_supplementals);
        records.append(&mut self.r2_supplementals);
        records.append(&mut self.r1_secondaries);
        records.append(&mut self.r2_secondaries);

        Ok(Template {
            name,
            records,
            raw_records: None,
            r1,
            r2,
            r1_supplementals,
            r2_supplementals,
            r1_secondaries,
            r2_secondaries,
            mi: MoleculeId::None,
        })
    }
}

impl Template {
    /// Returns the primary R1 read if present.
    ///
    /// The primary R1 is the main (non-secondary, non-supplementary) alignment
    /// for the first read in a pair (or the only read for single-end data).
    ///
    /// # Returns
    ///
    /// Reference to the primary R1 record, or `None` if not present.
    /// Returns `None` in raw-byte mode (use `raw_r1()` instead).
    #[must_use]
    pub fn r1(&self) -> Option<&RecordBuf> {
        if self.records.is_empty() {
            return None;
        }
        self.r1.map(|_| &self.records[0])
    }

    /// Returns the primary R2 read if present.
    /// Returns `None` in raw-byte mode (use `raw_r2()` instead).
    ///
    /// The primary R2 is the main (non-secondary, non-supplementary) alignment
    /// for the second read in a paired-end template.
    ///
    /// # Returns
    ///
    /// Reference to the primary R2 record, or `None` if not present
    #[must_use]
    pub fn r2(&self) -> Option<&RecordBuf> {
        if self.records.is_empty() {
            return None;
        }
        self.r2.map(|(i, _)| &self.records[i])
    }

    /// Returns supplementary alignments for R1.
    /// Returns empty in raw-byte mode.
    ///
    /// Supplementary alignments represent additional mapping locations for chimeric
    /// reads (reads that map to multiple locations due to structural variants or
    /// fusion events).
    ///
    /// # Returns
    ///
    /// Vector of references to R1 supplementary records
    #[must_use]
    pub fn r1_supplementals(&self) -> Vec<&RecordBuf> {
        if let Some((start, end)) = self.r1_supplementals {
            if end <= self.records.len() {
                return self.records[start..end].iter().collect();
            }
        }
        Vec::new()
    }

    /// Returns supplementary alignments for R2.
    /// Returns empty in raw-byte mode.
    ///
    /// Supplementary alignments represent additional mapping locations for chimeric
    /// reads (reads that map to multiple locations due to structural variants or
    /// fusion events).
    ///
    /// # Returns
    ///
    /// Vector of references to R2 supplementary records
    #[must_use]
    pub fn r2_supplementals(&self) -> Vec<&RecordBuf> {
        if let Some((start, end)) = self.r2_supplementals {
            if end <= self.records.len() {
                return self.records[start..end].iter().collect();
            }
        }
        Vec::new()
    }

    /// Returns secondary alignments for R1.
    /// Returns empty in raw-byte mode.
    ///
    /// Secondary alignments represent alternative mapping locations for reads with
    /// ambiguous alignments. These differ from supplementary alignments in that they
    /// represent the full read aligned to multiple locations (not chimeric reads).
    ///
    /// # Returns
    ///
    /// Vector of references to R1 secondary records
    #[must_use]
    pub fn r1_secondaries(&self) -> Vec<&RecordBuf> {
        if let Some((start, end)) = self.r1_secondaries {
            if end <= self.records.len() {
                return self.records[start..end].iter().collect();
            }
        }
        Vec::new()
    }

    /// Returns secondary alignments for R2.
    /// Returns empty in raw-byte mode.
    ///
    /// Secondary alignments represent alternative mapping locations for reads with
    /// ambiguous alignments. These differ from supplementary alignments in that they
    /// represent the full read aligned to multiple locations (not chimeric reads).
    ///
    /// # Returns
    ///
    /// Vector of references to R2 secondary records
    #[must_use]
    pub fn r2_secondaries(&self) -> Vec<&RecordBuf> {
        if let Some((start, end)) = self.r2_secondaries {
            if end <= self.records.len() {
                return self.records[start..end].iter().collect();
            }
        }
        Vec::new()
    }

    /// Creates a new empty template with the given name
    ///
    /// # Arguments
    ///
    /// * `name` - The query name for this template
    ///
    /// # Returns
    ///
    /// An empty template with no records
    #[must_use]
    pub fn new(name: Vec<u8>) -> Self {
        Template {
            name,
            records: Vec::new(),
            raw_records: None,
            r1: None,
            r2: None,
            r1_supplementals: None,
            r2_supplementals: None,
            r1_secondaries: None,
            r2_secondaries: None,
            mi: MoleculeId::None,
        }
    }

    /// Creates a Template from an iterator of records with the same query name
    ///
    /// Automatically organizes records into r1/r2, primary/supplementary/secondary categories.
    /// Assumes all records have the same query name.
    ///
    /// # Arguments
    ///
    /// * `records` - Iterator of records to organize into a template
    ///
    /// # Returns
    ///
    /// A Template with records organized by type, or an error if multiple primary reads found
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Multiple non-secondary, non-supplemental R1s are found
    /// - Multiple non-secondary, non-supplemental R2s are found
    /// - Iterator is empty (no name available)
    ///
    /// # Panics
    ///
    /// - If not all the records have the same query name
    pub fn from_records<I>(records: I) -> Result<Self>
    where
        I: IntoIterator<Item = RecordBuf>,
    {
        let mut builder = Builder::new();

        for record in records {
            builder.push(record)?;
        }

        builder.build()
    }

    /// Returns an iterator over primary (non-secondary, non-supplementary) reads in this template
    ///
    /// Returns both R1 and R2 primary reads if present.
    ///
    /// # Returns
    ///
    /// Iterator over primary read records
    pub fn primary_reads(&self) -> impl Iterator<Item = &RecordBuf> {
        let records = &self.records;
        self.r1.iter().chain(self.r2.iter()).filter_map(move |(start, _)| records.get(*start))
    }

    /// Consumes the template and returns owned primary reads.
    ///
    /// This method takes ownership of the template and returns the primary R1 and R2
    /// reads without cloning. Use this when you need owned records and won't access
    /// the template again.
    ///
    /// # Returns
    ///
    /// Tuple of `(Option<RecordBuf>, Option<RecordBuf>)` containing owned R1 and R2 records
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let template = Template::from_records(records)?;
    /// let (r1, r2) = template.into_primary_reads();
    /// if let Some(r1) = r1 {
    ///     writer.write_record(r1)?; // No clone needed
    /// }
    /// ```
    #[must_use]
    pub fn into_primary_reads(mut self) -> (Option<RecordBuf>, Option<RecordBuf>) {
        if self.records.is_empty() {
            return (None, None);
        }
        // Extract r2 first (higher index) to avoid shifting r1's position
        let r2 = self.r2.and_then(|(i, _)| self.records.get_mut(i).map(std::mem::take));
        let r1 = self.r1.and_then(|_| self.records.get_mut(0).map(std::mem::take));
        (r1, r2)
    }

    /// Returns an iterator over all supplementary and secondary reads
    ///
    /// # Returns
    ///
    /// Iterator over non-primary reads
    pub fn all_supplementary_and_secondary(&self) -> impl Iterator<Item = &RecordBuf> {
        let records = &self.records;
        let iter = self
            .r1_supplementals
            .iter()
            .chain(self.r2_supplementals.iter())
            .chain(self.r1_secondaries.iter())
            .chain(self.r2_secondaries.iter());
        iter.flat_map(move |(start, end)| {
            let s = (*start).min(records.len());
            let e = (*end).min(records.len());
            records[s..e].iter()
        })
    }

    /// Returns an iterator over all reads in this template
    ///
    /// Includes primary, supplementary, and secondary alignments.
    ///
    /// # Returns
    ///
    /// Iterator over all records
    pub fn all_reads(&self) -> impl Iterator<Item = &RecordBuf> {
        self.primary_reads().chain(self.all_supplementary_and_secondary())
    }

    /// Consumes the template and returns all records as a vector.
    ///
    /// This is useful when you need to take ownership of all records
    /// for further processing.
    #[must_use]
    pub fn into_records(self) -> Vec<RecordBuf> {
        self.records
    }

    /// Returns an iterator over all R1 reads (primary, supplementary, and secondary)
    ///
    /// # Returns
    ///
    /// Iterator over all R1 records
    pub fn all_r1s(&self) -> impl Iterator<Item = &RecordBuf> {
        let records = &self.records;
        let iter =
            self.r1.iter().chain(self.r1_secondaries.iter()).chain(self.r1_supplementals.iter());
        iter.flat_map(move |(start, end)| {
            let s = (*start).min(records.len());
            let e = (*end).min(records.len());
            records[s..e].iter()
        })
    }

    /// Returns an iterator over all R2 reads (primary, supplementary, and secondary)
    ///
    /// # Returns
    ///
    /// Iterator over all R2 records
    pub fn all_r2s(&self) -> impl Iterator<Item = &RecordBuf> {
        let records = &self.records;
        let iter =
            self.r2.iter().chain(self.r2_secondaries.iter()).chain(self.r2_supplementals.iter());
        iter.flat_map(move |(start, end)| {
            let s = (*start).min(records.len());
            let e = (*end).min(records.len());
            records[s..e].iter()
        })
    }

    /// Returns the total count of records in this template
    ///
    /// # Returns
    ///
    /// Total number of records (primary, supplementary, and secondary)
    #[must_use]
    pub fn read_count(&self) -> usize {
        if let Some(ref rr) = self.raw_records { rr.len() } else { self.records.len() }
    }

    /// Returns true if this template is in raw-byte mode.
    #[inline]
    #[must_use]
    pub fn is_raw_byte_mode(&self) -> bool {
        self.raw_records.is_some()
    }

    /// Returns the raw R1 bytes if in raw-byte mode.
    #[must_use]
    pub fn raw_r1(&self) -> Option<&[u8]> {
        let rr = self.raw_records.as_ref()?;
        self.r1.map(|_| rr[0].as_slice())
    }

    /// Returns the raw R2 bytes if in raw-byte mode.
    #[must_use]
    pub fn raw_r2(&self) -> Option<&[u8]> {
        let rr = self.raw_records.as_ref()?;
        self.r2.map(|(i, _)| rr[i].as_slice())
    }

    /// Returns all raw records if in raw-byte mode.
    #[must_use]
    pub fn all_raw_records(&self) -> Option<&[Vec<u8>]> {
        self.raw_records.as_deref()
    }

    /// Consumes self and returns the raw records if in raw-byte mode.
    #[must_use]
    pub fn into_raw_records(self) -> Option<Vec<Vec<u8>>> {
        self.raw_records
    }

    /// Returns mutable access to all raw records if in raw-byte mode.
    pub fn all_raw_records_mut(&mut self) -> Option<&mut [Vec<u8>]> {
        self.raw_records.as_deref_mut()
    }

    /// Build a Template from raw BAM byte records, categorizing by flags.
    ///
    /// The records are categorized using `bam_fields::flags()` to determine
    /// R1/R2/supplementary/secondary status, with the same index-pair scheme
    /// as `Builder::build()`.
    ///
    /// # Errors
    ///
    /// Returns an error if multiple primary R1s or R2s are found.
    #[allow(clippy::too_many_lines)]
    pub fn from_raw_records(mut raw_records: Vec<Vec<u8>>) -> Result<Self> {
        use crate::sort::bam_fields;

        if raw_records.is_empty() {
            bail!("No records given to from_raw_records");
        }

        // Guard against truncated records
        for (i, r) in raw_records.iter().enumerate() {
            if r.len() < bam_fields::MIN_BAM_HEADER_LEN {
                bail!(
                    "Raw BAM record {i} too short to parse ({} < {})",
                    r.len(),
                    bam_fields::MIN_BAM_HEADER_LEN
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

        // Fast path for common 2-record paired-end case (no supplementals/secondaries)
        if raw_records.len() == 2 {
            let f0 = bam_fields::flags(&raw_records[0]);
            let f1 = bam_fields::flags(&raw_records[1]);
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
                        records: Vec::new(),
                        raw_records: Some(raw_records),
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
                bail!("Template name mismatch in from_raw_records");
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
            let flg = bam_fields::flags(rec);
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
        // (matching Builder::build() ordering, with supplementals/secondaries reversed)
        let mut ordered: Vec<Vec<u8>> = Vec::with_capacity(raw_records.len());
        let mut take = |idx: usize| -> Vec<u8> { std::mem::take(&mut raw_records[idx]) };

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
            records: Vec::new(),
            raw_records: Some(ordered),
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
        let r1 = self.r1()?;
        let r2 = self.r2()?;

        // Both must be mapped
        if r1.flags().is_unmapped() || r2.flags().is_unmapped() {
            return None;
        }

        // Must be on the same reference
        if r1.reference_sequence_id() != r2.reference_sequence_id() {
            return None;
        }

        // Use R1 to determine orientation (htsjdk uses the first record passed)
        Some(get_pair_orientation(r1))
    }

    /// Fixes mate information for paired-end reads.
    ///
    /// This method follows htsjdk's `SamPairUtil.setMateInfo` behavior exactly:
    ///
    /// **For primary R1/R2 pairs:**
    /// - Case 1 (both mapped): Sets mate ref, mate pos, mate strand flag, mate unmapped=false,
    ///   MQ tag, MC tag (if valid CIGAR), and computes TLEN
    /// - Case 2 (both unmapped): Clears ref/pos to unmapped values, sets mate unmapped=true,
    ///   removes MQ and MC tags, sets TLEN=0
    /// - Case 3 (one mapped, one unmapped): Places unmapped read at mapped read's coordinates,
    ///   sets appropriate mate info, MQ only on unmapped read, TLEN=0
    ///
    /// **For supplementary alignments:**
    /// - Sets mate information based on the mate's primary alignment
    /// - Sets ms (mate score) tag if AS (alignment score) is available
    ///
    /// # Errors
    ///
    /// Returns an error if the CIGAR operations cannot be parsed (malformed CIGAR data)
    #[allow(clippy::too_many_lines)]
    pub fn fix_mate_info(&mut self) -> Result<()> {
        if self.records.is_empty() {
            return Ok(());
        }

        let mapq_tag = Tag::new(b'M', b'Q');
        let mate_cigar_tag = Tag::new(b'M', b'C');
        let mate_score_tag = Tag::new(b'm', b's');
        let align_score_tag = Tag::new(b'A', b'S');

        // Fix mate info for primary R1/R2 pair
        // Following htsjdk's SamPairUtil.setMateInfo exactly
        if let (Some((r1_i, _)), Some((r2_i, _))) = (self.r1, self.r2) {
            let r1_is_unmapped = self.records[r1_i].flags().is_unmapped();
            let r2_is_unmapped = self.records[r2_i].flags().is_unmapped();

            // Get alignment scores for mate score tags (used in all cases)
            let r1_as = self.records[r1_i].data().get(&align_score_tag).and_then(extract_int_value);
            let r2_as = self.records[r2_i].data().get(&align_score_tag).and_then(extract_int_value);

            if !r1_is_unmapped && !r2_is_unmapped {
                // Case 1: Both reads are mapped
                self.set_mate_info_both_mapped(r1_i, r2_i, mapq_tag, mate_cigar_tag)?;
            } else if r1_is_unmapped && r2_is_unmapped {
                // Case 2: Both reads are unmapped
                self.set_mate_info_both_unmapped(r1_i, r2_i, mapq_tag, mate_cigar_tag);
            } else {
                // Case 3: One mapped, one unmapped
                self.set_mate_info_one_unmapped(
                    r1_i,
                    r2_i,
                    r1_is_unmapped,
                    mapq_tag,
                    mate_cigar_tag,
                )?;
            }

            // Set mate score tags (ms) in all cases
            // Use Int8 to match fgbio's tag type encoding
            if let Some(as_value) = r2_as {
                self.records[r1_i]
                    .data_mut()
                    .insert(mate_score_tag, to_smallest_signed_int(as_value));
            }
            if let Some(as_value) = r1_as {
                self.records[r2_i]
                    .data_mut()
                    .insert(mate_score_tag, to_smallest_signed_int(as_value));
            }
        }

        // Fix mate info for R1 supplementals (mate is primary R2)
        // Following htsjdk's setMateInformationOnSupplementalAlignment:
        // - Mate reference index (RNEXT)
        // - Mate alignment start (PNEXT)
        // - Mate reverse strand flag (0x20)
        // - Mate unmapped flag (0x8)
        // - TLEN is set to negative of mate primary's TLEN
        // - MQ tag is set to mate primary's mapping quality
        // - MC tag is set to mate primary's CIGAR
        // - ms tag is set to mate primary's AS value
        if let Some((r2_i, _)) = self.r2 {
            let r2_ref_id = self.records[r2_i].reference_sequence_id();
            let r2_pos = self.records[r2_i].alignment_start();
            let r2_is_reverse = self.records[r2_i].flags().is_reverse_complemented();
            let r2_is_unmapped = self.records[r2_i].flags().is_unmapped();
            let r2_tlen = self.records[r2_i].template_length();
            let r2_mapq = self.records[r2_i].mapping_quality();
            let r2_cigar_str = cigar_to_string(self.records[r2_i].cigar())?;
            let r2_as = self.records[r2_i].data().get(&align_score_tag).and_then(extract_int_value);

            if let Some((start, end)) = self.r1_supplementals {
                for i in start..end {
                    // Set mate reference and position
                    *self.records[i].mate_reference_sequence_id_mut() = r2_ref_id;
                    *self.records[i].mate_alignment_start_mut() = r2_pos;

                    // Set mate flags (reverse strand and unmapped)
                    set_mate_flags(&mut self.records[i], r2_is_reverse, r2_is_unmapped);

                    // Set TLEN to negative of mate primary's TLEN
                    *self.records[i].template_length_mut() = -r2_tlen;

                    // Set MQ tag - use Int8 to match fgbio's tag type encoding
                    let r2_mapq_value = r2_mapq.map_or(255, |m| i32::from(u8::from(m)));
                    self.records[i]
                        .data_mut()
                        .insert(mapq_tag, to_smallest_signed_int(r2_mapq_value));

                    // Set MC tag if CIGAR is available
                    if !r2_cigar_str.is_empty() && r2_cigar_str != "*" && !r2_is_unmapped {
                        self.records[i]
                            .data_mut()
                            .insert(mate_cigar_tag, BufValue::from(r2_cigar_str.clone()));
                    }

                    // Set ms tag if AS is available - use Int8 to match fgbio
                    if let Some(as_value) = r2_as {
                        self.records[i]
                            .data_mut()
                            .insert(mate_score_tag, to_smallest_signed_int(as_value));
                    }
                }
            }
        }

        // Fix mate info for R2 supplementals (mate is primary R1)
        if let Some((r1_i, _)) = self.r1 {
            let r1_ref_id = self.records[r1_i].reference_sequence_id();
            let r1_pos = self.records[r1_i].alignment_start();
            let r1_is_reverse = self.records[r1_i].flags().is_reverse_complemented();
            let r1_is_unmapped = self.records[r1_i].flags().is_unmapped();
            let r1_tlen = self.records[r1_i].template_length();
            let r1_mapq = self.records[r1_i].mapping_quality();
            let r1_cigar_str = cigar_to_string(self.records[r1_i].cigar())?;
            let r1_as = self.records[r1_i].data().get(&align_score_tag).and_then(extract_int_value);

            if let Some((start, end)) = self.r2_supplementals {
                for i in start..end {
                    // Set mate reference and position
                    *self.records[i].mate_reference_sequence_id_mut() = r1_ref_id;
                    *self.records[i].mate_alignment_start_mut() = r1_pos;

                    // Set mate flags (reverse strand and unmapped)
                    set_mate_flags(&mut self.records[i], r1_is_reverse, r1_is_unmapped);

                    // Set TLEN to negative of mate primary's TLEN
                    *self.records[i].template_length_mut() = -r1_tlen;

                    // Set MQ tag - use Int8 to match fgbio's tag type encoding
                    let r1_mapq_value = r1_mapq.map_or(255, |m| i32::from(u8::from(m)));
                    self.records[i]
                        .data_mut()
                        .insert(mapq_tag, to_smallest_signed_int(r1_mapq_value));

                    // Set MC tag if CIGAR is available
                    if !r1_cigar_str.is_empty() && r1_cigar_str != "*" && !r1_is_unmapped {
                        self.records[i]
                            .data_mut()
                            .insert(mate_cigar_tag, BufValue::from(r1_cigar_str.clone()));
                    }

                    // Set ms tag if AS is available - use Int8 to match fgbio
                    if let Some(as_value) = r1_as {
                        self.records[i]
                            .data_mut()
                            .insert(mate_score_tag, to_smallest_signed_int(as_value));
                    }
                }
            }
        }

        Ok(())
    }

    /// Sets mate information when both reads are mapped.
    /// Following htsjdk's SamPairUtil.setMateInfo case 1.
    fn set_mate_info_both_mapped(
        &mut self,
        r1_i: usize,
        r2_i: usize,
        mapq_tag: Tag,
        mate_cigar_tag: Tag,
    ) -> Result<()> {
        // Get R2's info for R1
        let r2_ref_id = self.records[r2_i].reference_sequence_id();
        let r2_pos = self.records[r2_i].alignment_start();
        let r2_is_reverse = self.records[r2_i].flags().is_reverse_complemented();
        let r2_mapq = self.records[r2_i].mapping_quality();
        let r2_cigar_str = cigar_to_string(self.records[r2_i].cigar())?;

        // Get R1's info for R2
        let r1_ref_id = self.records[r1_i].reference_sequence_id();
        let r1_pos = self.records[r1_i].alignment_start();
        let r1_is_reverse = self.records[r1_i].flags().is_reverse_complemented();
        let r1_mapq = self.records[r1_i].mapping_quality();
        let r1_cigar_str = cigar_to_string(self.records[r1_i].cigar())?;

        // Set mate info on R1 from R2
        *self.records[r1_i].mate_reference_sequence_id_mut() = r2_ref_id;
        *self.records[r1_i].mate_alignment_start_mut() = r2_pos;
        set_mate_flags(&mut self.records[r1_i], r2_is_reverse, false);
        // Use Int8 for MQ to match fgbio's tag type encoding
        let r2_mapq_value = r2_mapq.map_or(255, |m| i32::from(u8::from(m)));
        self.records[r1_i].data_mut().insert(mapq_tag, to_smallest_signed_int(r2_mapq_value));
        if !r2_cigar_str.is_empty() && r2_cigar_str != "*" {
            self.records[r1_i].data_mut().insert(mate_cigar_tag, BufValue::from(r2_cigar_str));
        }

        // Set mate info on R2 from R1
        *self.records[r2_i].mate_reference_sequence_id_mut() = r1_ref_id;
        *self.records[r2_i].mate_alignment_start_mut() = r1_pos;
        set_mate_flags(&mut self.records[r2_i], r1_is_reverse, false);
        // Use Int8 for MQ to match fgbio's tag type encoding
        let r1_mapq_value = r1_mapq.map_or(255, |m| i32::from(u8::from(m)));
        self.records[r2_i].data_mut().insert(mapq_tag, to_smallest_signed_int(r1_mapq_value));
        if !r1_cigar_str.is_empty() && r1_cigar_str != "*" {
            self.records[r2_i].data_mut().insert(mate_cigar_tag, BufValue::from(r1_cigar_str));
        }

        // Compute and set insert size (TLEN)
        let insert_size = compute_insert_size(&self.records[r1_i], &self.records[r2_i]);
        *self.records[r1_i].template_length_mut() = insert_size;
        *self.records[r2_i].template_length_mut() = -insert_size;

        Ok(())
    }

    /// Sets mate information when both reads are unmapped.
    /// Following htsjdk's SamPairUtil.setMateInfo case 2.
    fn set_mate_info_both_unmapped(
        &mut self,
        r1_i: usize,
        r2_i: usize,
        mapq_tag: Tag,
        mate_cigar_tag: Tag,
    ) {
        // Get strand info before modifying
        let r1_is_reverse = self.records[r1_i].flags().is_reverse_complemented();
        let r2_is_reverse = self.records[r2_i].flags().is_reverse_complemented();

        // R1: set to unmapped coordinates
        *self.records[r1_i].reference_sequence_id_mut() = None;
        *self.records[r1_i].alignment_start_mut() = None;
        *self.records[r1_i].mate_reference_sequence_id_mut() = None;
        *self.records[r1_i].mate_alignment_start_mut() = None;
        set_mate_flags(&mut self.records[r1_i], r2_is_reverse, true);
        self.records[r1_i].data_mut().remove(&mapq_tag);
        self.records[r1_i].data_mut().remove(&mate_cigar_tag);
        *self.records[r1_i].template_length_mut() = 0;

        // R2: set to unmapped coordinates
        *self.records[r2_i].reference_sequence_id_mut() = None;
        *self.records[r2_i].alignment_start_mut() = None;
        *self.records[r2_i].mate_reference_sequence_id_mut() = None;
        *self.records[r2_i].mate_alignment_start_mut() = None;
        set_mate_flags(&mut self.records[r2_i], r1_is_reverse, true);
        self.records[r2_i].data_mut().remove(&mapq_tag);
        self.records[r2_i].data_mut().remove(&mate_cigar_tag);
        *self.records[r2_i].template_length_mut() = 0;
    }

    /// Sets mate information when one read is mapped and one is unmapped.
    /// Following htsjdk's SamPairUtil.setMateInfo case 3.
    fn set_mate_info_one_unmapped(
        &mut self,
        r1_i: usize,
        r2_i: usize,
        r1_is_unmapped: bool,
        mapq_tag: Tag,
        mate_cigar_tag: Tag,
    ) -> Result<()> {
        let (mapped_i, unmapped_i) = if r1_is_unmapped { (r2_i, r1_i) } else { (r1_i, r2_i) };

        // Get mapped read's info
        let mapped_ref_id = self.records[mapped_i].reference_sequence_id();
        let mapped_pos = self.records[mapped_i].alignment_start();
        let mapped_is_reverse = self.records[mapped_i].flags().is_reverse_complemented();
        let mapped_mapq = self.records[mapped_i].mapping_quality();
        let mapped_cigar_str = cigar_to_string(self.records[mapped_i].cigar())?;

        // Get unmapped read's strand (for mate strand flag on mapped read)
        let unmapped_is_reverse = self.records[unmapped_i].flags().is_reverse_complemented();

        // Place unmapped read at mapped read's coordinates (htsjdk behavior)
        *self.records[unmapped_i].reference_sequence_id_mut() = mapped_ref_id;
        *self.records[unmapped_i].alignment_start_mut() = mapped_pos;

        // Set mate info on mapped read (mate is unmapped)
        // After placing unmapped at mapped's coords, unmapped's ref/pos = mapped's ref/pos
        *self.records[mapped_i].mate_reference_sequence_id_mut() = mapped_ref_id;
        *self.records[mapped_i].mate_alignment_start_mut() = mapped_pos;
        set_mate_flags(&mut self.records[mapped_i], unmapped_is_reverse, true);
        self.records[mapped_i].data_mut().remove(&mapq_tag);
        self.records[mapped_i].data_mut().remove(&mate_cigar_tag);
        *self.records[mapped_i].template_length_mut() = 0;

        // Set mate info on unmapped read (mate is mapped)
        *self.records[unmapped_i].mate_reference_sequence_id_mut() = mapped_ref_id;
        *self.records[unmapped_i].mate_alignment_start_mut() = mapped_pos;
        set_mate_flags(&mut self.records[unmapped_i], mapped_is_reverse, false);
        // Use Int8 for MQ to match fgbio's tag type encoding
        let mapped_mapq_value = mapped_mapq.map_or(255, |m| i32::from(u8::from(m)));
        self.records[unmapped_i]
            .data_mut()
            .insert(mapq_tag, to_smallest_signed_int(mapped_mapq_value));
        if !mapped_cigar_str.is_empty() && mapped_cigar_str != "*" {
            self.records[unmapped_i]
                .data_mut()
                .insert(mate_cigar_tag, BufValue::from(mapped_cigar_str));
        }
        *self.records[unmapped_i].template_length_mut() = 0;

        Ok(())
    }
}

/// Computes the insert size (TLEN) for two reads.
/// Following htsjdk's SamPairUtil.computeInsertSize exactly.
///
/// The insert size is computed using 5' positions:
/// - For forward strand reads: alignment start (leftmost position)
/// - For reverse strand reads: alignment end (rightmost position)
///
/// Returns 0 if either read is unmapped or if reads are on different references.
fn compute_insert_size(rec1: &RecordBuf, rec2: &RecordBuf) -> i32 {
    // If either read is unmapped, return 0
    if rec1.flags().is_unmapped() || rec2.flags().is_unmapped() {
        return 0;
    }

    // If reads are on different references, return 0
    if rec1.reference_sequence_id() != rec2.reference_sequence_id() {
        return 0;
    }

    // Get alignment start positions
    let pos1 = rec1.alignment_start().map_or(0, |p| i32::try_from(usize::from(p)).unwrap_or(0));
    let pos2 = rec2.alignment_start().map_or(0, |p| i32::try_from(usize::from(p)).unwrap_or(0));

    // Calculate alignment ends (1-based, inclusive)
    let end1 = pos1 + alignment_length(rec1.cigar()) - 1;
    let end2 = pos2 + alignment_length(rec2.cigar()) - 1;

    // Use 5' positions: alignment end for reverse strand, alignment start for forward
    // This matches htsjdk's computeInsertSize exactly
    let first_5prime = if rec1.flags().is_reverse_complemented() { end1 } else { pos1 };
    let second_5prime = if rec2.flags().is_reverse_complemented() { end2 } else { pos2 };

    // Compute insert size with adjustment
    let adjustment = if second_5prime >= first_5prime { 1 } else { -1 };
    second_5prime - first_5prime + adjustment
}

/// Calculates the alignment length from a CIGAR (reference-consuming operations).
fn alignment_length(cigar: &noodles::sam::alignment::record_buf::Cigar) -> i32 {
    use noodles::sam::alignment::record::cigar::op::Kind;
    let mut len = 0i32;
    for op in cigar.iter().flatten() {
        match op.kind() {
            Kind::Match
            | Kind::Deletion
            | Kind::Skip
            | Kind::SequenceMatch
            | Kind::SequenceMismatch => {
                len += i32::try_from(op.len()).unwrap_or(0);
            }
            _ => {}
        }
    }
    len
}

/// Determines the pair orientation for a paired read.
///
/// This matches htsjdk's `SamPairUtil.getPairOrientation()` algorithm exactly.
///
/// # Arguments
/// * `record` - A paired, mapped read with a mapped mate on the same reference
///
/// # Returns
/// The pair orientation (FR, RF, or Tandem)
///
/// # Panics
/// The caller must ensure the record is paired, mapped, has a mapped mate,
/// and is on the same reference as its mate.
#[allow(clippy::cast_possible_wrap)]
fn get_pair_orientation(record: &RecordBuf) -> PairOrientation {
    let is_reverse = record.flags().is_reverse_complemented();
    let mate_reverse = record.flags().is_mate_reverse_complemented();

    // Same strand = TANDEM
    if is_reverse == mate_reverse {
        return PairOrientation::Tandem;
    }

    // Now determine if FR or RF using htsjdk's logic:
    // positiveStrandFivePrimePos = readIsOnReverseStrand ? mateStart : alignmentStart
    // negativeStrandFivePrimePos = readIsOnReverseStrand ? alignmentEnd : alignmentStart + insertSize
    // FR if positiveStrandFivePrimePos < negativeStrandFivePrimePos

    let alignment_start = record.alignment_start().map_or(0, usize::from);
    let mate_start = record.mate_alignment_start().map_or(0, usize::from);
    let insert_size = record.template_length();

    let (positive_five_prime, negative_five_prime) = if is_reverse {
        // This read is on reverse strand, mate is on positive strand
        // positiveStrandFivePrimePos = mateStart
        // negativeStrandFivePrimePos = alignmentEnd (this read's end)
        let end = record_utils::alignment_end(record).unwrap_or(alignment_start);
        (mate_start as i64, end as i64)
    } else {
        // This read is on positive strand, mate is on reverse strand
        // positiveStrandFivePrimePos = alignmentStart (this read's start)
        // negativeStrandFivePrimePos = alignmentStart + insertSize
        (alignment_start as i64, alignment_start as i64 + i64::from(insert_size))
    };

    // FR if positive strand 5' < negative strand 5'
    if positive_five_prime < negative_five_prime {
        PairOrientation::FR
    } else {
        PairOrientation::RF
    }
}

/// Sets the mate flags (`MATE_REVERSE_COMPLEMENTED` and `MATE_UNMAPPED`) on a record
/// based on the mate's flags.
///
/// This follows htsjdk's `setMateInformationOnSupplementalAlignment` behavior:
/// - Sets `MATE_REVERSE_COMPLEMENTED` (0x20) if mate is reverse complemented
/// - Sets `MATE_UNMAPPED` (0x8) if mate is unmapped
fn set_mate_flags(record: &mut RecordBuf, mate_is_reverse: bool, mate_is_unmapped: bool) {
    use noodles::sam::alignment::record::Flags;

    let mut flags = record.flags();

    // Clear and then conditionally set MATE_REVERSE_COMPLEMENTED
    flags.remove(Flags::MATE_REVERSE_COMPLEMENTED);
    if mate_is_reverse {
        flags.insert(Flags::MATE_REVERSE_COMPLEMENTED);
    }

    // Clear and then conditionally set MATE_UNMAPPED
    flags.remove(Flags::MATE_UNMAPPED);
    if mate_is_unmapped {
        flags.insert(Flags::MATE_UNMAPPED);
    }

    *record.flags_mut() = flags;
}

/// Extracts an i32 value from any integer `BufValue` variant.
///
/// SAM tags can store integers in various sizes (Int8, Int16, Int32, `UInt8`, `UInt16`, `UInt32`).
/// This function extracts the value as an i32 from any of these variants.
///
/// # Arguments
///
/// * `value` - The `BufValue` to extract from
///
/// # Returns
///
/// `Some(i32)` if the value is any integer type, `None` otherwise.
fn extract_int_value(value: &BufValue) -> Option<i32> {
    match value {
        BufValue::Int8(i) => Some(i32::from(*i)),
        BufValue::Int16(i) => Some(i32::from(*i)),
        BufValue::Int32(i) => Some(*i),
        BufValue::UInt8(i) => Some(i32::from(*i)),
        BufValue::UInt16(i) => Some(i32::from(*i)),
        BufValue::UInt32(i) => i32::try_from(*i).ok(),
        _ => None,
    }
}

/// Converts a CIGAR operation kind to its SAM character representation.
///
/// Maps CIGAR operation types to their single-character SAM format codes:
/// - M (Match/Mismatch)
/// - I (Insertion)
/// - D (Deletion)
/// - N (Skipped region)
/// - S (Soft clip)
/// - H (Hard clip)
/// - P (Padding)
/// - = (Sequence match)
/// - X (Sequence mismatch)
///
/// # Arguments
///
/// * `kind` - The CIGAR operation kind to convert
///
/// # Returns
///
/// Single character representing the operation in SAM format
///
/// # Panics
///
/// This function does not panic. The match is exhaustive over all CIGAR operation kinds.
fn kind_to_char(kind: noodles::sam::alignment::record::cigar::op::Kind) -> char {
    use noodles::sam::alignment::record::cigar::op::Kind;
    match kind {
        Kind::Match => 'M',
        Kind::Insertion => 'I',
        Kind::Deletion => 'D',
        Kind::Skip => 'N',
        Kind::SoftClip => 'S',
        Kind::HardClip => 'H',
        Kind::Pad => 'P',
        Kind::SequenceMatch => '=',
        Kind::SequenceMismatch => 'X',
    }
}

/// Converts a CIGAR to its string representation
///
/// # Arguments
///
/// * `cigar` - The CIGAR to convert
///
/// # Returns
///
/// String representation of the CIGAR (e.g., "100M5I3D") or "*" if empty
///
/// # Errors
///
/// Returns an error if the CIGAR operations cannot be parsed (malformed CIGAR data)
fn cigar_to_string(cigar: &noodles::sam::alignment::record_buf::Cigar) -> Result<String> {
    if cigar.is_empty() {
        return Ok(String::from("*"));
    }

    // Pre-allocate: estimate ~4 chars per op (e.g., "10M" = 3 chars, "100M" = 4 chars)
    let mut result = String::with_capacity(cigar.as_ref().len() * 4);
    for op in cigar.iter() {
        let op = op?;
        let _ = write!(result, "{}{}", op.len(), kind_to_char(op.kind()));
    }
    Ok(result)
}

/// Iterator that groups consecutive records by query name into `Template` objects.
///
/// This iterator consumes a stream of SAM/BAM records and yields complete `Template`
/// objects. It assumes the input is queryname-sorted or queryname-grouped, meaning all
/// records with the same query name appear consecutively in the stream.
///
/// # Requirements
///
/// The input records MUST be queryname-sorted or grouped. If records with the same
/// name are not consecutive, they will be placed in separate templates.
///
/// # Examples
///
/// ```rust,ignore
/// let record_iter = bam_reader.records();
/// let template_iter = TemplateIterator::new(record_iter);
///
/// for template in template_iter {
///     let template = template?;
///     // Process template
/// }
/// ```
pub struct TemplateIterator<I>
where
    I: Iterator<Item = Result<RecordBuf>>,
{
    record_iter: I,
    builder: Builder,
}

impl<I> TemplateIterator<I>
where
    I: Iterator<Item = Result<RecordBuf>>,
{
    /// Creates a new `TemplateIterator`
    ///
    /// # Arguments
    ///
    /// * `record_iter` - Iterator over records to group into templates
    ///
    /// # Returns
    ///
    /// A new `TemplateIterator` ready to produce Templates
    pub fn new(record_iter: I) -> Self {
        TemplateIterator { record_iter, builder: Builder::new() }
    }
}

impl<I> Iterator for TemplateIterator<I>
where
    I: Iterator<Item = Result<RecordBuf>>,
{
    type Item = Result<Template>;

    /// Returns the next Template from the record iterator
    ///
    /// Reads records until the query name changes, then returns the completed
    /// template. On EOF, returns any remaining template, then None.
    ///
    /// # Returns
    ///
    /// * `Some(Ok(Template))` - Next template successfully read
    /// * `Some(Err(_))` - Error from underlying iterator
    /// * `None` - No more templates (EOF reached)
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.record_iter.next() {
                None => {
                    // EOF - return any remaining template
                    if self.builder.is_empty() {
                        return None;
                    }
                    return Some(self.builder.build());
                }
                Some(Err(e)) => {
                    // Error from underlying iterator
                    return Some(Err(e));
                }
                Some(Ok(record)) => {
                    let name: Vec<u8> = if let Some(n) = record.name() {
                        Vec::from(<_ as AsRef<[u8]>>::as_ref(n))
                    } else {
                        Vec::new()
                    };

                    if self.builder.name.iter().any(|n| *n == name) || self.builder.is_empty() {
                        if let Err(e) = self.builder.push(record) {
                            return Some(Err(e));
                        }
                    } else {
                        let template: std::result::Result<Template, anyhow::Error> =
                            self.builder.build();
                        if let Err(e) = self.builder.push(record) {
                            return Some(Err(e));
                        }
                        return Some(template);
                    }
                }
            }
        }
    }
}

// ============================================================================
// MemoryEstimate Implementations
// ============================================================================

impl MemoryEstimate for Template {
    fn estimate_heap_size(&self) -> usize {
        // name: Vec<u8>
        let name_size = self.name.capacity();

        // records: Vec<RecordBuf>
        // Each RecordBuf contains:
        // - name: Option<BString> (Vec<u8>)
        // - sequence: Sequence (Vec<u8>)
        // - quality_scores: QualityScores (Vec<u8>)
        // - cigar: Cigar (Vec<Op>)
        // - data: Data (IndexMap of tags to values)
        let records_size: usize = self.records.iter().map(estimate_record_buf_heap_size).sum();
        let records_vec_overhead = self.records.capacity() * std::mem::size_of::<RecordBuf>();

        let raw_records_size = self.raw_records.as_ref().map_or(0, |rr| {
            rr.iter().map(Vec::capacity).sum::<usize>()
                + rr.capacity() * std::mem::size_of::<Vec<u8>>()
        });

        name_size + records_size + records_vec_overhead + raw_records_size
    }
}

impl MemoryEstimate for TemplateBatch {
    fn estimate_heap_size(&self) -> usize {
        self.iter().map(MemoryEstimate::estimate_heap_size).sum::<usize>()
            + self.capacity() * std::mem::size_of::<Template>()
    }
}

/// Estimate heap size of a `RecordBuf`.
///
/// This function estimates the heap memory used by a noodles `RecordBuf`,
/// including its name, sequence, quality scores, CIGAR, and data fields.
///
/// noodles `RecordBuf` layout:
/// - name: `Option<BString>` (`Vec<u8>`) - heap allocated
/// - sequence: `Sequence` (`Vec<u8>`, 1 byte per base, unpacked ASCII)
/// - `quality_scores`: `QualityScores` (`Vec<u8>`)
/// - cigar: `Cigar` (`Vec<Op>`, each `Op` is 4 bytes)
/// - data: Data (`IndexMap`<Tag, Value>) - significant overhead from hash table
#[must_use]
pub fn estimate_record_buf_heap_size(record: &RecordBuf) -> usize {
    // Name: Option<BString> which is a Vec<u8>
    let name_size = record.name().map_or(0, |n| n.len());

    // Sequence: noodles stores bases as unpacked ASCII (1 byte per base in Vec<u8>)
    // .len() returns number of bases, which equals the heap allocation in bytes
    let seq_len = record.sequence().len();

    // Quality scores: Vec<u8>, one byte per base
    let qual_len = record.quality_scores().len();

    // CIGAR: Vec<Op> where each Op is 4 bytes (u32)
    let cigar_ops = record.cigar().as_ref().len();
    let cigar_size = cigar_ops * 4;

    // Data fields: IndexMap<Tag, Value>
    // IndexMap stores entries in a Vec<Bucket<(K, V)>> plus a hash table Vec<usize>.
    // - Each entry: Tag (2 bytes) + Value enum (~40 bytes on stack for noodles Value)
    //   + padding = ~48 bytes per entry in the entries vec
    // - Hash table: capacity * 8 bytes (indices) + capacity * 8 bytes (hashes)
    //   With ~87.5% load factor, capacity ~= count * 1.15
    // - Values with heap allocation (String, Array): add ~24 bytes (Vec header) + data
    //   Typical BAM tags: RX/QX (string, ~20 bytes), MI (int, no heap), CB/CR (string, ~16 bytes)
    //   Estimate ~32 bytes average heap per string-type field, ~50% of fields are strings
    let data_fields = record.data().iter().count();
    let entry_capacity = (data_fields * 115) / 100 + 1; // ~1.15x for load factor
    let entries_size = data_fields * 48; // entry vec storage
    let hash_table_size = entry_capacity * 16; // hash + index arrays
    let value_heap_size = data_fields * 16; // average heap per field (50% strings × 32 bytes)
    let data_size = entries_size + hash_table_size + value_heap_size;

    // Note: RecordBuf inline struct overhead (flags, position, etc.) is accounted for
    // by the caller via `capacity * size_of::<RecordBuf>()`, so we only count heap allocations here.
    name_size + seq_len + qual_len + cigar_size + data_size
}

#[cfg(test)]
#[allow(clippy::similar_names)]
mod tests {
    use super::*;
    use crate::sam::builder::RecordBuilder;
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::record::cigar::Op;
    use noodles::sam::alignment::record::cigar::op::Kind;
    use noodles::sam::alignment::record_buf::Cigar as CigarBuf;

    fn create_test_record(name: &[u8], flags: u16) -> RecordBuf {
        RecordBuilder::new()
            .name(std::str::from_utf8(name).unwrap())
            .flags(Flags::from(flags))
            .build()
    }

    // SAM flag constants
    const FLAG_PAIRED: u16 = 0x1;
    const FLAG_READ1: u16 = 0x40;
    const FLAG_READ2: u16 = 0x80;
    const FLAG_SECONDARY: u16 = 0x100;
    const FLAG_SUPPLEMENTARY: u16 = 0x800;

    #[test]
    fn test_template_new() {
        let template = Template::new(b"read1".to_vec());
        assert_eq!(template.name, b"read1");
        assert_eq!(template.read_count(), 0);
    }

    #[test]
    fn test_builder_empty() {
        let builder = Builder::new();
        assert!(builder.is_empty());
        assert_eq!(builder.len(), 0);
    }

    #[test]
    fn test_builder_push_r1_only() -> Result<()> {
        let mut builder = Builder::new();
        let r1 = create_test_record(b"read1", 0);
        builder.push(r1)?;
        assert_eq!(builder.len(), 1);
        assert!(!builder.is_empty());
        Ok(())
    }

    #[test]
    fn test_builder_push_paired_end() -> Result<()> {
        let mut builder = Builder::new();
        let r1 = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r2 = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2);
        builder.push(r1)?;
        builder.push(r2)?;
        assert_eq!(builder.len(), 2);
        Ok(())
    }

    #[test]
    fn test_builder_push_supplementary() -> Result<()> {
        let mut builder = Builder::new();
        let supp = create_test_record(b"read1", FLAG_SUPPLEMENTARY);
        builder.push(supp)?;
        assert_eq!(builder.len(), 1);
        Ok(())
    }

    #[test]
    fn test_builder_push_secondary() -> Result<()> {
        let mut builder = Builder::new();
        let sec = create_test_record(b"read1", FLAG_SECONDARY);
        builder.push(sec)?;
        assert_eq!(builder.len(), 1);
        Ok(())
    }

    #[test]
    fn test_builder_error_on_name_mismatch() {
        let mut builder = Builder::new();
        let r1 = create_test_record(b"read1", 0);
        let r2 = create_test_record(b"read2", 0);
        builder.push(r1).unwrap();
        let result = builder.push(r2);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("mismatch"));
    }

    #[test]
    fn test_builder_error_on_multiple_r1() {
        let mut builder = Builder::new();
        let r1a = create_test_record(b"read1", 0);
        let r1b = create_test_record(b"read1", 0);
        builder.push(r1a).unwrap();
        let result = builder.push(r1b);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Multiple non-secondary"));
    }

    #[test]
    fn test_builder_error_on_multiple_r2() {
        let mut builder = Builder::new();
        let r2a = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2);
        let r2b = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2);
        builder.push(r2a).unwrap();
        let result = builder.push(r2b);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Multiple non-secondary"));
    }

    #[test]
    fn test_builder_build_empty_error() {
        let mut builder = Builder::new();
        let result = builder.build();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("No records"));
    }

    #[test]
    fn test_builder_build_success() -> Result<()> {
        let mut builder = Builder::new();
        let r1 = create_test_record(b"read1", 0);
        builder.push(r1)?;
        let template = builder.build()?;
        assert_eq!(template.name, b"read1");
        assert_eq!(template.read_count(), 1);
        assert!(template.r1().is_some());
        Ok(())
    }

    #[test]
    fn test_template_from_records() -> Result<()> {
        let records = vec![
            create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1),
            create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2),
        ];
        let template = Template::from_records(records)?;
        assert_eq!(template.read_count(), 2);
        assert!(template.r1().is_some());
        assert!(template.r2().is_some());
        Ok(())
    }

    #[test]
    fn test_template_primary_reads() -> Result<()> {
        let mut builder = Builder::new();
        builder.push(create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1))?;
        builder.push(create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2))?;
        let template = builder.build()?;
        let primaries: Vec<_> = template.primary_reads().collect();
        assert_eq!(primaries.len(), 2);
        Ok(())
    }

    #[test]
    fn test_template_into_primary_reads() -> Result<()> {
        let records = vec![
            create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1),
            create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2),
        ];
        let template = Template::from_records(records)?;

        // Verify we have both reads before consuming
        assert!(template.r1().is_some());
        assert!(template.r2().is_some());

        // Consume the template and get owned reads
        let (r1, r2) = template.into_primary_reads();

        // Both should be Some
        assert!(r1.is_some());
        assert!(r2.is_some());

        // Verify the records have the expected flags
        let r1 = r1.unwrap();
        let r2 = r2.unwrap();
        assert!(r1.flags().is_first_segment());
        assert!(r2.flags().is_last_segment());

        Ok(())
    }

    #[test]
    fn test_template_into_primary_reads_r1_only() -> Result<()> {
        let records = vec![create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1)];
        let template = Template::from_records(records)?;

        let (r1, r2) = template.into_primary_reads();

        assert!(r1.is_some());
        assert!(r2.is_none());

        Ok(())
    }

    #[test]
    fn test_template_all_reads() -> Result<()> {
        let mut builder = Builder::new();
        builder.push(create_test_record(b"read1", 0))?;
        builder.push(create_test_record(b"read1", FLAG_SUPPLEMENTARY))?;
        builder.push(create_test_record(b"read1", FLAG_SECONDARY))?;
        let template = builder.build()?;
        let all: Vec<_> = template.all_reads().collect();
        assert_eq!(all.len(), 3);
        Ok(())
    }

    #[test]
    fn test_template_all_r1s() -> Result<()> {
        let mut builder = Builder::new();
        builder.push(create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1))?;
        builder
            .push(create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY))?;
        let template = builder.build()?;
        let all_r1s: Vec<_> = template.all_r1s().collect();
        assert_eq!(all_r1s.len(), 2);
        Ok(())
    }

    #[test]
    fn test_template_all_r2s() -> Result<()> {
        let mut builder = Builder::new();
        builder.push(create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2))?;
        builder.push(create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SECONDARY))?;
        let template = builder.build()?;
        let all_r2s: Vec<_> = template.all_r2s().collect();
        assert_eq!(all_r2s.len(), 2);
        Ok(())
    }

    #[test]
    fn test_template_supplementals_and_secondaries() -> Result<()> {
        let mut builder = Builder::new();
        builder.push(create_test_record(b"read1", 0))?;
        builder.push(create_test_record(b"read1", FLAG_SUPPLEMENTARY))?;
        builder.push(create_test_record(b"read1", FLAG_SECONDARY))?;
        let template = builder.build()?;
        let non_primary: Vec<_> = template.all_supplementary_and_secondary().collect();
        assert_eq!(non_primary.len(), 2);
        Ok(())
    }

    #[test]
    fn test_cigar_to_string() {
        use noodles::sam::alignment::record_buf::Cigar;
        let cigar = Cigar::default();
        assert_eq!(cigar_to_string(&cigar).unwrap(), "*");
    }

    #[test]
    fn test_cigar_to_string_with_ops() {
        let ops =
            vec![Op::new(Kind::Match, 10), Op::new(Kind::Deletion, 2), Op::new(Kind::Match, 5)];
        let cigar = CigarBuf::from(ops);
        let result = cigar_to_string(&cigar).unwrap();
        assert_eq!(result, "10M2D5M");
    }

    #[test]
    fn test_template_builder_reuse() -> Result<()> {
        let mut builder = Builder::new();
        builder.push(create_test_record(b"read1", 0))?;
        let template1 = builder.build()?;
        assert_eq!(template1.name, b"read1");
        builder.push(create_test_record(b"read2", 0))?;
        let template2 = builder.build()?;
        assert_eq!(template2.name, b"read2");
        Ok(())
    }

    #[test]
    fn test_template_r1_supplementals() -> Result<()> {
        let mut builder = Builder::new();
        builder.push(create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1))?;
        builder
            .push(create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY))?;
        builder
            .push(create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY))?;
        let template = builder.build()?;
        let supps = template.r1_supplementals();
        assert_eq!(supps.len(), 2);
        Ok(())
    }

    #[test]
    fn test_template_r2_supplementals() -> Result<()> {
        let mut builder = Builder::new();
        builder.push(create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2))?;
        builder
            .push(create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY))?;
        let template = builder.build()?;
        let supps = template.r2_supplementals();
        assert_eq!(supps.len(), 1);
        Ok(())
    }

    #[test]
    fn test_template_r1_secondaries() -> Result<()> {
        let mut builder = Builder::new();
        builder.push(create_test_record(b"read1", 0))?;
        builder.push(create_test_record(b"read1", FLAG_SECONDARY))?;
        let template = builder.build()?;
        let secs = template.r1_secondaries();
        assert_eq!(secs.len(), 1);
        Ok(())
    }

    #[test]
    fn test_template_r2_secondaries() -> Result<()> {
        let mut builder = Builder::new();
        builder.push(create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2))?;
        builder.push(create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SECONDARY))?;
        builder.push(create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SECONDARY))?;
        let template = builder.build()?;
        let secs = template.r2_secondaries();
        assert_eq!(secs.len(), 2);
        Ok(())
    }

    // Tests for extract_int_value helper function
    #[test]
    fn test_extract_int_value_int8() {
        let value = BufValue::Int8(42);
        assert_eq!(extract_int_value(&value), Some(42));
    }

    #[test]
    fn test_extract_int_value_int16() {
        let value = BufValue::Int16(1234);
        assert_eq!(extract_int_value(&value), Some(1234));
    }

    #[test]
    fn test_extract_int_value_int32() {
        let value = BufValue::Int32(123_456);
        assert_eq!(extract_int_value(&value), Some(123_456));
    }

    #[test]
    fn test_extract_int_value_uint8() {
        let value = BufValue::UInt8(255);
        assert_eq!(extract_int_value(&value), Some(255));
    }

    #[test]
    fn test_extract_int_value_uint16() {
        let value = BufValue::UInt16(65535);
        assert_eq!(extract_int_value(&value), Some(65535));
    }

    #[test]
    fn test_extract_int_value_uint32() {
        let value = BufValue::UInt32(100_000);
        assert_eq!(extract_int_value(&value), Some(100_000));
    }

    #[test]
    fn test_extract_int_value_uint32_overflow() {
        // Value too large to fit in i32
        let value = BufValue::UInt32(u32::MAX);
        assert_eq!(extract_int_value(&value), None);
    }

    #[test]
    fn test_extract_int_value_string_returns_none() {
        let value = BufValue::String("test".into());
        assert_eq!(extract_int_value(&value), None);
    }

    #[test]
    fn test_extract_int_value_negative() {
        let value = BufValue::Int8(-42);
        assert_eq!(extract_int_value(&value), Some(-42));

        let value = BufValue::Int16(-1234);
        assert_eq!(extract_int_value(&value), Some(-1234));

        let value = BufValue::Int32(-123_456);
        assert_eq!(extract_int_value(&value), Some(-123_456));
    }

    /// Helper to create a mapped record with position and CIGAR for `fix_mate_info` tests
    fn create_mapped_record(name: &[u8], flags: u16, pos: usize, mapq: u8) -> RecordBuf {
        create_mapped_record_with_tlen(name, flags, pos, mapq, 0)
    }

    /// Helper to create a mapped record with position, CIGAR, and TLEN for `fix_mate_info` tests
    fn create_mapped_record_with_tlen(
        name: &[u8],
        flags: u16,
        pos: usize,
        mapq: u8,
        tlen: i32,
    ) -> RecordBuf {
        let is_read1 = (flags & FLAG_READ1) != 0;
        let is_paired = (flags & FLAG_PAIRED) != 0;

        let mut builder = RecordBuilder::new()
            .name(std::str::from_utf8(name).unwrap())
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .cigar("100M")
            .reference_sequence_id(0)
            .alignment_start(pos)
            .mapping_quality(mapq)
            .template_length(tlen);

        if is_paired {
            builder = builder.first_segment(is_read1);
        }

        // Handle secondary/supplementary flags - directly set flags for non-standard combinations
        let mut record = builder.build();
        // Apply any additional flags not covered by builder methods
        *record.flags_mut() = Flags::from(flags);

        record
    }

    /// Tests that `fix_mate_info` correctly sets ms tag when AS is stored as `Int8`.
    /// This was a bug where only `Int32` AS values were recognized.
    #[test]
    fn test_fix_mate_info_sets_ms_tag_with_int8_as() -> Result<()> {
        let as_tag = Tag::new(b'A', b'S');
        let ms_tag = Tag::new(b'm', b's');

        let mut r1 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30);
        r1.data_mut().insert(as_tag, BufValue::Int8(55)); // AS as Int8 (small value)

        let mut r2 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40);
        r2.data_mut().insert(as_tag, BufValue::Int8(44)); // AS as Int8 (small value)

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // R1 should have ms tag with R2's AS value (44)
        let r1_ms = template.records[0].data().get(&ms_tag);
        assert!(r1_ms.is_some(), "R1 should have ms tag");
        assert_eq!(extract_int_value(r1_ms.unwrap()), Some(44));

        // R2 should have ms tag with R1's AS value (55)
        let r2_ms = template.records[1].data().get(&ms_tag);
        assert!(r2_ms.is_some(), "R2 should have ms tag");
        assert_eq!(extract_int_value(r2_ms.unwrap()), Some(55));

        Ok(())
    }

    /// Tests that `fix_mate_info` correctly sets ms tag when AS is stored as `UInt8`
    #[test]
    fn test_fix_mate_info_sets_ms_tag_with_uint8_as() -> Result<()> {
        let as_tag = Tag::new(b'A', b'S');
        let ms_tag = Tag::new(b'm', b's');

        let mut r1 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30);
        r1.data_mut().insert(as_tag, BufValue::UInt8(77)); // AS as UInt8

        let mut r2 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40);
        r2.data_mut().insert(as_tag, BufValue::UInt8(88)); // AS as UInt8

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // R1 should have ms tag with R2's AS value (88)
        let r1_ms = template.records[0].data().get(&ms_tag);
        assert!(r1_ms.is_some(), "R1 should have ms tag");
        assert_eq!(extract_int_value(r1_ms.unwrap()), Some(88));

        // R2 should have ms tag with R1's AS value (77)
        let r2_ms = template.records[1].data().get(&ms_tag);
        assert!(r2_ms.is_some(), "R2 should have ms tag");
        assert_eq!(extract_int_value(r2_ms.unwrap()), Some(77));

        Ok(())
    }

    /// Tests that `fix_mate_info` correctly sets ms tag when AS is stored as `Int16`
    #[test]
    fn test_fix_mate_info_sets_ms_tag_with_int16_as() -> Result<()> {
        let as_tag = Tag::new(b'A', b'S');
        let ms_tag = Tag::new(b'm', b's');

        let mut r1 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30);
        r1.data_mut().insert(as_tag, BufValue::Int16(1000)); // AS as Int16

        let mut r2 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40);
        r2.data_mut().insert(as_tag, BufValue::Int16(2000)); // AS as Int16

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // R1 should have ms tag with R2's AS value (2000)
        let r1_ms = template.records[0].data().get(&ms_tag);
        assert!(r1_ms.is_some(), "R1 should have ms tag");
        assert_eq!(extract_int_value(r1_ms.unwrap()), Some(2000));

        // R2 should have ms tag with R1's AS value (1000)
        let r2_ms = template.records[1].data().get(&ms_tag);
        assert!(r2_ms.is_some(), "R2 should have ms tag");
        assert_eq!(extract_int_value(r2_ms.unwrap()), Some(1000));

        Ok(())
    }

    /// Tests that `fix_mate_info` correctly sets ms tag when AS is stored as `Int32`
    #[test]
    fn test_fix_mate_info_sets_ms_tag_with_int32_as() -> Result<()> {
        let as_tag = Tag::new(b'A', b'S');
        let ms_tag = Tag::new(b'm', b's');

        let mut r1 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30);
        r1.data_mut().insert(as_tag, BufValue::Int32(100_000)); // AS as Int32

        let mut r2 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40);
        r2.data_mut().insert(as_tag, BufValue::Int32(200_000)); // AS as Int32

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // R1 should have ms tag with R2's AS value
        let r1_ms = template.records[0].data().get(&ms_tag);
        assert!(r1_ms.is_some(), "R1 should have ms tag");
        assert_eq!(extract_int_value(r1_ms.unwrap()), Some(200_000));

        // R2 should have ms tag with R1's AS value
        let r2_ms = template.records[1].data().get(&ms_tag);
        assert!(r2_ms.is_some(), "R2 should have ms tag");
        assert_eq!(extract_int_value(r2_ms.unwrap()), Some(100_000));

        Ok(())
    }

    /// Tests that `fix_mate_info` does not set ms tag when AS is missing
    #[test]
    fn test_fix_mate_info_no_ms_tag_without_as() -> Result<()> {
        let ms_tag = Tag::new(b'm', b's');

        let r1 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30);
        let r2 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40);

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // Neither record should have ms tag since AS is missing
        assert!(template.records[0].data().get(&ms_tag).is_none(), "R1 should not have ms tag");
        assert!(template.records[1].data().get(&ms_tag).is_none(), "R2 should not have ms tag");

        Ok(())
    }

    /// Tests that `fix_mate_info` sets ms tag on supplementary alignments
    #[test]
    fn test_fix_mate_info_sets_ms_tag_on_supplementals() -> Result<()> {
        let as_tag = Tag::new(b'A', b'S');
        let ms_tag = Tag::new(b'm', b's');

        let mut r1 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30);
        r1.data_mut().insert(as_tag, BufValue::Int8(55));

        let mut r2 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40);
        r2.data_mut().insert(as_tag, BufValue::Int8(44));

        let mut r1_supp =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY, 300, 20);
        r1_supp.data_mut().insert(as_tag, BufValue::Int8(33));

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        builder.push(r1_supp)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // R1 supplementary should have ms tag with R2's AS value (44)
        let r1_supp_ms = template.records[2].data().get(&ms_tag);
        assert!(r1_supp_ms.is_some(), "R1 supplementary should have ms tag");
        assert_eq!(extract_int_value(r1_supp_ms.unwrap()), Some(44));

        Ok(())
    }

    /// Tests that `fix_mate_info` sets TLEN on R1 supplementary alignments.
    /// TLEN should be set to negative of mate primary's (R2) TLEN.
    #[test]
    fn test_fix_mate_info_sets_tlen_on_r1_supplementals() -> Result<()> {
        // R1 primary with TLEN=200
        let r1 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30, 200);
        // R2 primary with TLEN=-200
        let r2 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40, -200);
        // R1 supplementary with original TLEN=0
        let r1_supp = create_mapped_record_with_tlen(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY,
            300,
            20,
            0,
        );

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        builder.push(r1_supp)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // R1 supplementary TLEN should be -(-101) = 101 (negative of R2's recalculated TLEN)
        // R1(pos=100,forward) and R2(pos=200,forward) → 5' positions are 100 and 200 → insert_size = 101
        assert_eq!(
            template.records[2].template_length(),
            101,
            "R1 supplementary TLEN should be negative of R2's TLEN"
        );

        Ok(())
    }

    /// Tests that `fix_mate_info` sets TLEN on R2 supplementary alignments.
    /// TLEN should be set to negative of mate primary's (R1) TLEN.
    #[test]
    fn test_fix_mate_info_sets_tlen_on_r2_supplementals() -> Result<()> {
        // R1 primary with TLEN=300
        let r1 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30, 300);
        // R2 primary with TLEN=-300
        let r2 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40, -300);
        // R2 supplementary with original TLEN=0
        let r2_supp = create_mapped_record_with_tlen(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY,
            400,
            25,
            0,
        );

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        builder.push(r2_supp)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // R2 supplementary TLEN should be -(101) = -101 (negative of R1's recalculated TLEN)
        // R1(pos=100,forward) and R2(pos=200,forward) → 5' positions are 100 and 200 → insert_size = 101
        assert_eq!(
            template.records[2].template_length(),
            -101,
            "R2 supplementary TLEN should be negative of R1's TLEN"
        );

        Ok(())
    }

    /// Tests that `fix_mate_info` handles multiple supplementary alignments
    #[test]
    fn test_fix_mate_info_sets_tlen_on_multiple_supplementals() -> Result<()> {
        // R1 primary with TLEN=500
        let r1 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30, 500);
        // R2 primary with TLEN=-500
        let r2 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40, -500);
        // Two R1 supplementaries
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
        // Two R2 supplementaries
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

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        builder.push(r1_supp1)?;
        builder.push(r1_supp2)?;
        builder.push(r2_supp1)?;
        builder.push(r2_supp2)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // Order: R1[0], R2[1], R1_supp1[2], R1_supp2[3], R2_supp1[4], R2_supp2[5]
        // R1(pos=100,forward) and R2(pos=200,forward) → 5' positions are 100 and 200 → insert_size = 101

        // R1 supplementaries should have TLEN = -(-101) = 101
        assert_eq!(template.records[2].template_length(), 101, "R1 supp1 TLEN");
        assert_eq!(template.records[3].template_length(), 101, "R1 supp2 TLEN");

        // R2 supplementaries should have TLEN = -(101) = -101
        assert_eq!(template.records[4].template_length(), -101, "R2 supp1 TLEN");
        assert_eq!(template.records[5].template_length(), -101, "R2 supp2 TLEN");

        Ok(())
    }

    /// Tests that `fix_mate_info` recalculates TLEN for supplementaries based on primary positions
    #[test]
    fn test_fix_mate_info_tlen_recalculated() -> Result<()> {
        // Primary reads with TLEN=0 (initial value, will be recalculated)
        // Note: even if original TLEN is 0, fix_mate_info recalculates based on positions
        let r1 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30, 0);
        let r2 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 40, 0);
        let r1_supp = create_mapped_record_with_tlen(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY,
            300,
            20,
            999, // Non-zero original TLEN (will be recalculated)
        );

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        builder.push(r1_supp)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // R1(pos=100,forward) and R2(pos=200,forward) → 5' positions are 100 and 200 → insert_size = 101
        // R1 supplementary TLEN = -R2.TLEN = -(-101) = 101
        assert_eq!(
            template.records[2].template_length(),
            101,
            "R1 supplementary TLEN should be recalculated from positions"
        );

        Ok(())
    }

    /// Tests that `fix_mate_info` doesn't affect supplementaries when mate primary is missing
    #[test]
    fn test_fix_mate_info_no_mate_primary() -> Result<()> {
        // Only R1 primary (no R2 primary)
        let r1 = create_mapped_record_with_tlen(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30, 200);
        // R1 supplementary
        let r1_supp = create_mapped_record_with_tlen(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY,
            300,
            20,
            777, // Original TLEN
        );

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r1_supp)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // R1 supplementary TLEN should remain unchanged since there's no R2 primary
        assert_eq!(
            template.records[1].template_length(),
            777,
            "R1 supplementary TLEN should be unchanged without R2 primary"
        );

        Ok(())
    }

    /// Tests record ordering in Template matches fgbio's ordering
    /// Order should be: R1, R2, `r1_supps`, `r2_supps`, `r1_secondaries`, `r2_secondaries`
    #[test]
    fn test_template_record_ordering() -> Result<()> {
        let r1 = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r2 = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2);
        let r1_supp1 = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY);
        let r1_supp2 = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY);
        let r2_supp = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY);
        let r1_sec = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SECONDARY);
        let r2_sec1 = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SECONDARY);
        let r2_sec2 = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SECONDARY);

        // Push in scrambled order
        let mut builder = Builder::new();
        builder.push(r2_sec1)?;
        builder.push(r1_supp2)?;
        builder.push(r1)?;
        builder.push(r2_supp)?;
        builder.push(r1_sec)?;
        builder.push(r2)?;
        builder.push(r2_sec2)?;
        builder.push(r1_supp1)?;
        let template = builder.build()?;

        // Verify indices
        assert_eq!(template.r1, Some((0, 1)), "R1 should be at index 0");
        assert_eq!(template.r2, Some((1, 2)), "R2 should be at index 1");
        assert_eq!(template.r1_supplementals, Some((2, 4)), "R1 supps should be at indices 2-3");
        assert_eq!(template.r2_supplementals, Some((4, 5)), "R2 supps should be at index 4");
        assert_eq!(template.r1_secondaries, Some((5, 6)), "R1 secs should be at index 5");
        assert_eq!(template.r2_secondaries, Some((6, 8)), "R2 secs should be at indices 6-7");

        // Verify all_reads() returns in correct order
        let all_reads: Vec<_> = template.all_reads().collect();
        assert_eq!(all_reads.len(), 8);

        // Verify flags match expected ordering
        assert!(
            all_reads[0].flags().is_first_segment() && !all_reads[0].flags().is_supplementary()
        );
        assert!(all_reads[1].flags().is_last_segment() && !all_reads[1].flags().is_supplementary());
        assert!(all_reads[2].flags().is_first_segment() && all_reads[2].flags().is_supplementary());
        assert!(all_reads[3].flags().is_first_segment() && all_reads[3].flags().is_supplementary());
        assert!(all_reads[4].flags().is_last_segment() && all_reads[4].flags().is_supplementary());
        assert!(all_reads[5].flags().is_first_segment() && all_reads[5].flags().is_secondary());
        assert!(all_reads[6].flags().is_last_segment() && all_reads[6].flags().is_secondary());
        assert!(all_reads[7].flags().is_last_segment() && all_reads[7].flags().is_secondary());

        Ok(())
    }

    /// Tests that records within supplementary/secondary groups are in reverse input order.
    /// This matches fgbio's behavior which uses Scala List prepend (r :: list).
    #[test]
    fn test_record_ordering_within_groups_is_reversed() -> Result<()> {
        // Create records with distinct positions so we can verify ordering
        let r1 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30);
        let r2 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ2, 200, 30);

        // R1 supplementaries with positions 1000, 2000, 3000 (input order)
        let r1_supp1 =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY, 1000, 20);
        let r1_supp2 =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY, 2000, 20);
        let r1_supp3 =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY, 3000, 20);

        // R2 supplementaries with positions 4000, 5000 (input order)
        let r2_supp1 =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY, 4000, 20);
        let r2_supp2 =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY, 5000, 20);

        // R1 secondaries with positions 6000, 7000 (input order)
        let r1_sec1 =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SECONDARY, 6000, 10);
        let r1_sec2 =
            create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_SECONDARY, 7000, 10);

        // Push in input order
        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        builder.push(r1_supp1)?;
        builder.push(r1_supp2)?;
        builder.push(r1_supp3)?;
        builder.push(r2_supp1)?;
        builder.push(r2_supp2)?;
        builder.push(r1_sec1)?;
        builder.push(r1_sec2)?;

        let template = builder.build()?;

        // Helper to extract position as usize
        let get_pos = |rec: &RecordBuf| -> usize { rec.alignment_start().map_or(0, |p| p.get()) };

        // Verify R1 and R2 primaries are first
        assert_eq!(get_pos(&template.records[0]), 100, "R1 primary");
        assert_eq!(get_pos(&template.records[1]), 200, "R2 primary");

        // Verify R1 supplementaries are in REVERSE input order: 3000, 2000, 1000
        assert_eq!(get_pos(&template.records[2]), 3000, "R1 supp should be in reverse order");
        assert_eq!(get_pos(&template.records[3]), 2000, "R1 supp should be in reverse order");
        assert_eq!(get_pos(&template.records[4]), 1000, "R1 supp should be in reverse order");

        // Verify R2 supplementaries are in REVERSE input order: 5000, 4000
        assert_eq!(get_pos(&template.records[5]), 5000, "R2 supp should be in reverse order");
        assert_eq!(get_pos(&template.records[6]), 4000, "R2 supp should be in reverse order");

        // Verify R1 secondaries are in REVERSE input order: 7000, 6000
        assert_eq!(get_pos(&template.records[7]), 7000, "R1 sec should be in reverse order");
        assert_eq!(get_pos(&template.records[8]), 6000, "R1 sec should be in reverse order");

        Ok(())
    }

    /// Helper to create a mapped record with specific flags including reverse strand
    fn create_mapped_record_with_flags(
        name: &[u8],
        flags: u16,
        pos: usize,
        mapq: u8,
        tlen: i32,
        mate_ref_id: Option<usize>,
        mate_pos: Option<usize>,
    ) -> RecordBuf {
        let is_read1 = (flags & FLAG_READ1) != 0;
        let is_paired = (flags & FLAG_PAIRED) != 0;

        let mut builder = RecordBuilder::new()
            .name(std::str::from_utf8(name).unwrap())
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .cigar("100M")
            .reference_sequence_id(0)
            .alignment_start(pos)
            .mapping_quality(mapq)
            .template_length(tlen);

        if is_paired {
            builder = builder.first_segment(is_read1);
        }

        if let Some(mate_ref) = mate_ref_id {
            builder = builder.mate_reference_sequence_id(mate_ref);
        }
        if let Some(mate_p) = mate_pos {
            builder = builder.mate_alignment_start(mate_p);
        }

        // Build and then override with exact flags
        let mut record = builder.build();
        *record.flags_mut() = Flags::from(flags);

        record
    }

    // Flag constants for reverse strand
    const FLAG_REVERSE: u16 = 0x10;
    const FLAG_MATE_REVERSE: u16 = 0x20;
    const FLAG_UNMAPPED: u16 = 0x4;
    const FLAG_MATE_UNMAPPED: u16 = 0x8;

    /// Tests that `fix_mate_info` sets mate info on R1 supplementary alignments correctly.
    /// The mate is R2 primary, so supplementary should get R2's info.
    #[test]
    fn test_fix_mate_info_sets_full_mate_info_on_r1_supplementals() -> Result<()> {
        let mq_tag = Tag::new(b'M', b'Q');
        let mc_tag = Tag::new(b'M', b'C');
        let as_tag = Tag::new(b'A', b'S');
        let ms_tag = Tag::new(b'm', b's');

        // R1 primary at pos=100, forward strand, TLEN=200
        let mut r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1,
            100,
            30,
            200,
            Some(0),
            Some(200),
        );
        r1.data_mut().insert(as_tag, BufValue::Int32(100));

        // R2 primary at pos=200, reverse strand (0x10), TLEN=-200
        let mut r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE,
            200,
            40,
            -200,
            Some(0),
            Some(100),
        );
        r2.data_mut().insert(as_tag, BufValue::Int32(150));

        // R1 supplementary at pos=500, originally has wrong mate info (will be corrected)
        let r1_supp = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY | FLAG_MATE_REVERSE, // has wrong mate_reverse
            500,
            25,
            0, // wrong TLEN
            Some(0),
            Some(500), // wrong mate pos (points to self)
        );

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        builder.push(r1_supp)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // Verify R1 supplementary (at index 2) has correct mate info from R2 primary
        let supp = &template.records[2];

        // Mate position should be R2's position (200)
        assert_eq!(
            supp.mate_alignment_start().map(|p| p.get()),
            Some(200),
            "R1 supp should have mate pos from R2"
        );

        // Mate reference should be R2's reference (0)
        assert_eq!(
            supp.mate_reference_sequence_id(),
            Some(0),
            "R1 supp should have mate ref from R2"
        );

        // Mate reverse flag should be set (R2 is reverse)
        assert!(
            supp.flags().is_mate_reverse_complemented(),
            "R1 supp should have mate_reverse set since R2 is reverse"
        );

        // Mate unmapped flag should NOT be set (R2 is mapped)
        assert!(
            !supp.flags().is_mate_unmapped(),
            "R1 supp should NOT have mate_unmapped since R2 is mapped"
        );

        // TLEN should be negative of R2's TLEN
        assert_eq!(supp.template_length(), 200, "R1 supp TLEN should be -(-200) = 200");

        // MQ tag should be R2's mapping quality (40)
        let mq_value = supp.data().get(&mq_tag).and_then(extract_int_value);
        assert_eq!(mq_value, Some(40), "R1 supp MQ should be R2's mapq");

        // MC tag should be R2's CIGAR (100M)
        let mc_value = supp.data().get(&mc_tag);
        assert!(mc_value.is_some(), "R1 supp should have MC tag");

        // ms tag should be R2's AS value (150)
        let ms_value = supp.data().get(&ms_tag).and_then(extract_int_value);
        assert_eq!(ms_value, Some(150), "R1 supp ms should be R2's AS value");

        Ok(())
    }

    /// Tests that `fix_mate_info` sets mate info on R2 supplementary alignments correctly.
    /// The mate is R1 primary, so supplementary should get R1's info.
    #[test]
    fn test_fix_mate_info_sets_full_mate_info_on_r2_supplementals() -> Result<()> {
        let mq_tag = Tag::new(b'M', b'Q');
        let mc_tag = Tag::new(b'M', b'C');
        let as_tag = Tag::new(b'A', b'S');
        let ms_tag = Tag::new(b'm', b's');

        // R1 primary at pos=100, forward strand, TLEN=300
        let mut r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1,
            100,
            35,
            300,
            Some(0),
            Some(200),
        );
        r1.data_mut().insert(as_tag, BufValue::Int32(120));

        // R2 primary at pos=200, forward strand, TLEN=-300
        let mut r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2,
            200,
            45,
            -300,
            Some(0),
            Some(100),
        );
        r2.data_mut().insert(as_tag, BufValue::Int32(180));

        // R2 supplementary at pos=600, originally has wrong mate info
        let r2_supp = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_SUPPLEMENTARY,
            600,
            20,
            0, // wrong TLEN
            Some(0),
            Some(600), // wrong mate pos
        );

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        builder.push(r2_supp)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // Verify R2 supplementary (at index 2) has correct mate info from R1 primary
        let supp = &template.records[2];

        // Mate position should be R1's position (100)
        assert_eq!(
            supp.mate_alignment_start().map(|p| p.get()),
            Some(100),
            "R2 supp should have mate pos from R1"
        );

        // Mate reverse flag should NOT be set (R1 is forward)
        assert!(
            !supp.flags().is_mate_reverse_complemented(),
            "R2 supp should NOT have mate_reverse since R1 is forward"
        );

        // TLEN should be negative of R1's recalculated TLEN
        // R1(pos=100,forward) and R2(pos=200,forward) → 5' positions are 100 and 200 → insert_size = 101
        assert_eq!(supp.template_length(), -101, "R2 supp TLEN should be -(101) = -101");

        // MQ tag should be R1's mapping quality (35)
        let mq_value = supp.data().get(&mq_tag).and_then(extract_int_value);
        assert_eq!(mq_value, Some(35), "R2 supp MQ should be R1's mapq");

        // MC tag should be present
        assert!(supp.data().get(&mc_tag).is_some(), "R2 supp should have MC tag");

        // ms tag should be R1's AS value (120)
        let ms_value = supp.data().get(&ms_tag).and_then(extract_int_value);
        assert_eq!(ms_value, Some(120), "R2 supp ms should be R1's AS value");

        Ok(())
    }

    /// Tests that `fix_mate_info` correctly handles unmapped mate.
    /// When mate is unmapped, `mate_unmapped` flag should be set and MC tag should not be set.
    #[test]
    fn test_fix_mate_info_supplemental_with_unmapped_mate() -> Result<()> {
        let mc_tag = Tag::new(b'M', b'C');

        // R1 primary at pos=100, forward strand
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_UNMAPPED, // R1 knows mate is unmapped
            100,
            30,
            0, // TLEN=0 when mate unmapped
            Some(0),
            Some(100), // unmapped mate placed at R1's position
        );

        // R2 is unmapped but placed at R1's position (standard convention)
        let mut r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_UNMAPPED,
            100, // same position as R1 (convention for unmapped mate)
            0,   // mapq 0 for unmapped
            0,   // TLEN=0
            Some(0),
            Some(100),
        );
        // Clear CIGAR for unmapped read
        *r2.cigar_mut() = CigarBuf::default();

        // R1 supplementary, originally has incorrect mate_reverse flag set
        let r1_supp = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY | FLAG_MATE_REVERSE, // wrong
            500,
            25,
            999, // wrong TLEN
            Some(0),
            Some(500), // wrong mate pos
        );

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        builder.push(r1_supp)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // Verify R1 supplementary has correct flags for unmapped mate
        let supp = &template.records[2];

        // Mate unmapped flag should be set
        assert!(
            supp.flags().is_mate_unmapped(),
            "R1 supp should have mate_unmapped since R2 is unmapped"
        );

        // Mate reverse flag should NOT be set for unmapped mate
        assert!(
            !supp.flags().is_mate_reverse_complemented(),
            "R1 supp should NOT have mate_reverse when R2 is unmapped (unmapped reads have no orientation)"
        );

        // MC tag should NOT be set for unmapped mate
        assert!(
            supp.data().get(&mc_tag).is_none(),
            "R1 supp should NOT have MC tag when mate is unmapped"
        );

        // TLEN should be 0 (negative of mate's 0)
        assert_eq!(supp.template_length(), 0, "TLEN should be 0 when mate is unmapped");

        Ok(())
    }

    /// Tests that `fix_mate_info` correctly clears `mate_reverse` flag when mate is on forward strand.
    /// This tests the scenario where supplementary had incorrect `mate_reverse` flag.
    #[test]
    fn test_fix_mate_info_clears_incorrect_mate_reverse_flag() -> Result<()> {
        // R1 primary at pos=100, forward strand
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1,
            100,
            30,
            200,
            Some(0),
            Some(200),
        );

        // R2 primary at pos=200, forward strand (NOT reverse)
        let r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2, // no FLAG_REVERSE
            200,
            40,
            -200,
            Some(0),
            Some(100),
        );

        // R1 supplementary incorrectly has MATE_REVERSE set
        let r1_supp = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_SUPPLEMENTARY | FLAG_MATE_REVERSE, // incorrectly set
            500,
            25,
            0,
            Some(0),
            Some(500),
        );

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        builder.push(r1_supp)?;
        let mut template = builder.build()?;

        // Before fix, supplementary has MATE_REVERSE set incorrectly
        assert!(
            template.records[2].flags().is_mate_reverse_complemented(),
            "Before fix: supp incorrectly has mate_reverse"
        );

        template.fix_mate_info()?;

        // After fix, MATE_REVERSE should be cleared since R2 is forward
        assert!(
            !template.records[2].flags().is_mate_reverse_complemented(),
            "After fix: supp should NOT have mate_reverse since R2 is forward"
        );

        Ok(())
    }

    /// Tests that `fix_mate_info` handles multiple supplementary alignments.
    #[test]
    fn test_fix_mate_info_multiple_supplementals() -> Result<()> {
        let as_tag = Tag::new(b'A', b'S');

        // R1 primary, forward
        let r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1,
            100,
            30,
            400,
            Some(0),
            Some(300),
        );

        // R2 primary, reverse
        let mut r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE,
            300,
            50,
            -400,
            Some(0),
            Some(100),
        );
        r2.data_mut().insert(as_tag, BufValue::Int32(200));

        // Two R1 supplementaries
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

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        builder.push(r1_supp1)?;
        builder.push(r1_supp2)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // Both supplementaries should have correct mate info from R2
        // Note: after build(), supplementaries are at indices 2 and 3 (reversed order)
        // R1(pos=100,forward) and R2(pos=300,100M,reverse) → 5' positions are 100 and 399 → insert_size = 300
        for i in 2..=3 {
            let supp = &template.records[i];

            assert_eq!(
                supp.mate_alignment_start().map(|p| p.get()),
                Some(300),
                "Supp {i} should have mate pos 300"
            );
            assert!(
                supp.flags().is_mate_reverse_complemented(),
                "Supp {i} should have mate_reverse"
            );
            assert_eq!(supp.template_length(), 300, "Supp {i} TLEN should be 300");
        }

        Ok(())
    }

    // ============================================================================
    // Tests for compute_insert_size (ported from htsjdk SamPairUtilTest)
    // ============================================================================

    /// Helper to create a mapped record with specific strand orientation and CIGAR for insert size tests
    fn create_insert_size_record(
        name: &[u8],
        flags: u16,
        pos: usize,
        read_length: usize,
        ref_id: Option<usize>,
    ) -> RecordBuf {
        let is_read1 = (flags & FLAG_READ1) != 0;
        let is_paired = (flags & FLAG_PAIRED) != 0;

        let mut builder = RecordBuilder::new()
            .name(std::str::from_utf8(name).unwrap())
            .sequence(&"A".repeat(read_length))
            .qualities(&vec![30u8; read_length])
            .cigar(&format!("{read_length}M"))
            .alignment_start(pos)
            .mapping_quality(30);

        if let Some(rid) = ref_id {
            builder = builder.reference_sequence_id(rid);
        }

        if is_paired {
            builder = builder.first_segment(is_read1);
        }

        // Build and then override with exact flags
        let mut record = builder.build();
        *record.flags_mut() = Flags::from(flags);

        record
    }

    /// Test `compute_insert_size` for "normal innie" FR pair (htsjdk test case)
    /// R1: pos=1, length=100, forward
    /// R2: pos=500, length=100, reverse
    /// Expected: 5' positions are 1 (forward start) and 599 (reverse end)
    /// Insert size = 599 - 1 + 1 = 599
    #[test]
    fn test_compute_insert_size_normal_innie() {
        let r1 = create_insert_size_record(b"read1", FLAG_PAIRED | FLAG_READ1, 1, 100, Some(0));
        let r2 = create_insert_size_record(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE,
            500,
            100,
            Some(0),
        );

        let insert_size = compute_insert_size(&r1, &r2);
        // R1 5' = 1 (forward, start), R2 5' = 599 (reverse, end = 500+100-1)
        // second_5prime (599) >= first_5prime (1), so adjustment = +1
        // 599 - 1 + 1 = 599
        assert_eq!(insert_size, 599);
    }

    /// Test `compute_insert_size` for overlapping innie
    /// R1: pos=1, length=100, forward (ends at 100)
    /// R2: pos=50, length=100, reverse (ends at 149)
    #[test]
    fn test_compute_insert_size_overlapping_innie() {
        let r1 = create_insert_size_record(b"read1", FLAG_PAIRED | FLAG_READ1, 1, 100, Some(0));
        let r2 = create_insert_size_record(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE,
            50,
            100,
            Some(0),
        );

        let insert_size = compute_insert_size(&r1, &r2);
        // R1 5' = 1, R2 5' = 149 (50+100-1)
        // 149 - 1 + 1 = 149
        assert_eq!(insert_size, 149);
    }

    /// Test `compute_insert_size` for completely overlapping innie
    /// R1: pos=1, length=100, forward
    /// R2: pos=1, length=100, reverse
    #[test]
    fn test_compute_insert_size_completely_overlapping_innie() {
        let r1 = create_insert_size_record(b"read1", FLAG_PAIRED | FLAG_READ1, 1, 100, Some(0));
        let r2 = create_insert_size_record(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE,
            1,
            100,
            Some(0),
        );

        let insert_size = compute_insert_size(&r1, &r2);
        // R1 5' = 1, R2 5' = 100 (1+100-1)
        // 100 - 1 + 1 = 100
        assert_eq!(insert_size, 100);
    }

    /// Test `compute_insert_size` for outie pair (RF orientation)
    /// R1: pos=1, length=100, reverse
    /// R2: pos=500, length=100, forward
    #[test]
    fn test_compute_insert_size_normal_outie() {
        let r1 = create_insert_size_record(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_REVERSE,
            1,
            100,
            Some(0),
        );
        let r2 = create_insert_size_record(b"read1", FLAG_PAIRED | FLAG_READ2, 500, 100, Some(0));

        let insert_size = compute_insert_size(&r1, &r2);
        // R1 5' = 100 (reverse, end = 1+100-1), R2 5' = 500 (forward, start)
        // 500 - 100 + 1 = 401
        assert_eq!(insert_size, 401);
    }

    /// Test `compute_insert_size` for forward tandem
    /// R1: pos=1, length=100, reverse
    /// R2: pos=500, length=100, reverse
    #[test]
    fn test_compute_insert_size_forward_tandem() {
        let r1 = create_insert_size_record(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_REVERSE,
            1,
            100,
            Some(0),
        );
        let r2 = create_insert_size_record(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE,
            500,
            100,
            Some(0),
        );

        let insert_size = compute_insert_size(&r1, &r2);
        // R1 5' = 100 (reverse), R2 5' = 599 (reverse)
        // 599 - 100 + 1 = 500
        assert_eq!(insert_size, 500);
    }

    /// Test `compute_insert_size` for reverse tandem
    /// R1: pos=1, length=100, forward
    /// R2: pos=500, length=100, forward
    #[test]
    fn test_compute_insert_size_reverse_tandem() {
        let r1 = create_insert_size_record(b"read1", FLAG_PAIRED | FLAG_READ1, 1, 100, Some(0));
        let r2 = create_insert_size_record(b"read1", FLAG_PAIRED | FLAG_READ2, 500, 100, Some(0));

        let insert_size = compute_insert_size(&r1, &r2);
        // R1 5' = 1 (forward), R2 5' = 500 (forward)
        // 500 - 1 + 1 = 500
        assert_eq!(insert_size, 500);
    }

    /// Test `compute_insert_size` when second read is before first (negative insert size)
    /// R1: pos=500, length=100, forward
    /// R2: pos=1, length=100, reverse
    #[test]
    fn test_compute_insert_size_negative() {
        let r1 = create_insert_size_record(b"read1", FLAG_PAIRED | FLAG_READ1, 500, 100, Some(0));
        let r2 = create_insert_size_record(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE,
            1,
            100,
            Some(0),
        );

        let insert_size = compute_insert_size(&r1, &r2);
        // R1 5' = 500 (forward), R2 5' = 100 (reverse)
        // second_5prime (100) < first_5prime (500), so adjustment = -1
        // 100 - 500 + (-1) = -401
        assert_eq!(insert_size, -401);
    }

    /// Test `compute_insert_size` returns 0 when reads are on different references
    #[test]
    fn test_compute_insert_size_different_references() {
        let r1 = create_insert_size_record(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 100, Some(0));
        let r2 = create_insert_size_record(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE,
            200,
            100,
            Some(1), // different reference
        );

        let insert_size = compute_insert_size(&r1, &r2);
        assert_eq!(insert_size, 0, "Insert size should be 0 for different references");
    }

    /// Test `compute_insert_size` returns 0 when R1 is unmapped
    #[test]
    fn test_compute_insert_size_r1_unmapped() {
        let r1 = create_insert_size_record(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_UNMAPPED,
            100,
            100,
            Some(0),
        );
        let r2 = create_insert_size_record(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE,
            200,
            100,
            Some(0),
        );

        let insert_size = compute_insert_size(&r1, &r2);
        assert_eq!(insert_size, 0, "Insert size should be 0 when R1 is unmapped");
    }

    /// Test `compute_insert_size` returns 0 when R2 is unmapped
    #[test]
    fn test_compute_insert_size_r2_unmapped() {
        let r1 = create_insert_size_record(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 100, Some(0));
        let r2 = create_insert_size_record(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_UNMAPPED,
            200,
            100,
            Some(0),
        );

        let insert_size = compute_insert_size(&r1, &r2);
        assert_eq!(insert_size, 0, "Insert size should be 0 when R2 is unmapped");
    }

    /// Test `compute_insert_size` returns 0 when both reads are unmapped
    #[test]
    fn test_compute_insert_size_both_unmapped() {
        let r1 = create_insert_size_record(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_UNMAPPED,
            100,
            100,
            Some(0),
        );
        let r2 = create_insert_size_record(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_UNMAPPED,
            200,
            100,
            Some(0),
        );

        let insert_size = compute_insert_size(&r1, &r2);
        assert_eq!(insert_size, 0, "Insert size should be 0 when both reads are unmapped");
    }

    /// Test `compute_insert_size` with complex CIGAR (insertions don't affect reference length)
    #[test]
    fn test_compute_insert_size_with_indels() {
        // CIGAR: 50M10I40M = 90M on reference (insertions don't consume reference)
        let r1 = RecordBuilder::new()
            .name("read1")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("50M10I40M")
            .first_segment(true)
            .build();

        // CIGAR: 50M5D50M = 105 on reference (deletions consume reference)
        let mut r2 = RecordBuilder::new()
            .name("read1")
            .sequence(&"A".repeat(100))
            .qualities(&[30u8; 100])
            .reference_sequence_id(0)
            .alignment_start(200)
            .cigar("50M5D50M")
            .first_segment(false)
            .reverse_complement(true)
            .build();
        // Ensure FLAG_REVERSE is set
        *r2.flags_mut() = Flags::from(FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE);

        let insert_size = compute_insert_size(&r1, &r2);
        // R1: alignment length = 50 + 40 = 90, end = 90, 5' = 1 (forward)
        // R2: alignment length = 50 + 5 + 50 = 105, end = 304, 5' = 304 (reverse)
        // 304 - 1 + 1 = 304
        assert_eq!(insert_size, 304);
    }

    // ============================================================================
    // Tests for setMateInfo (both mapped case) - ported from htsjdk
    // ============================================================================

    /// Test `fix_mate_info` correctly sets TLEN for normal innie pair
    #[test]
    fn test_fix_mate_info_sets_tlen_normal_innie() -> Result<()> {
        let r1 = create_insert_size_record(b"read1", FLAG_PAIRED | FLAG_READ1, 1, 100, Some(0));
        let r2 = create_insert_size_record(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE,
            500,
            100,
            Some(0),
        );

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // Expected insert size: 599
        assert_eq!(template.records[0].template_length(), 599, "R1 TLEN");
        assert_eq!(template.records[1].template_length(), -599, "R2 TLEN");

        Ok(())
    }

    /// Test `fix_mate_info` correctly sets mate CIGAR (MC tag)
    #[test]
    fn test_fix_mate_info_sets_mate_cigar() -> Result<()> {
        let mc_tag = Tag::new(b'M', b'C');

        let r1 = RecordBuilder::new()
            .name("read1")
            .sequence(&"A".repeat(50))
            .qualities(&[30u8; 50])
            .reference_sequence_id(0)
            .alignment_start(1)
            .cigar("50M")
            .first_segment(true)
            .build();

        // Complex CIGAR for R2: 25M5I20M
        let mut r2 = RecordBuilder::new()
            .name("read1")
            .sequence(&"A".repeat(50))
            .qualities(&[30u8; 50])
            .reference_sequence_id(0)
            .alignment_start(500)
            .cigar("25M5I20M")
            .first_segment(false)
            .reverse_complement(true)
            .build();
        // Ensure FLAG_REVERSE is set
        *r2.flags_mut() = Flags::from(FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE);

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // R1's MC should be R2's CIGAR
        let r1_mc = template.records[0].data().get(&mc_tag);
        assert!(r1_mc.is_some(), "R1 should have MC tag");
        if let Some(BufValue::String(mc)) = r1_mc {
            assert_eq!(&mc[..], b"25M5I20M", "R1 MC should be R2's CIGAR");
        } else {
            panic!("MC tag should be a string");
        }

        // R2's MC should be R1's CIGAR
        let r2_mc = template.records[1].data().get(&mc_tag);
        assert!(r2_mc.is_some(), "R2 should have MC tag");
        if let Some(BufValue::String(mc)) = r2_mc {
            assert_eq!(&mc[..], b"50M", "R2 MC should be R1's CIGAR");
        } else {
            panic!("MC tag should be a string");
        }

        Ok(())
    }

    /// Test `fix_mate_info` correctly handles both unmapped case
    #[test]
    fn test_fix_mate_info_both_unmapped() -> Result<()> {
        let mq_tag = Tag::new(b'M', b'Q');
        let mc_tag = Tag::new(b'M', b'C');

        let mut r1 = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ1 | FLAG_UNMAPPED);
        *r1.reference_sequence_id_mut() = Some(0); // Originally had a reference
        *r1.data_mut() = [(mq_tag, BufValue::Int32(30))].into_iter().collect(); // Has MQ tag

        let mut r2 = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_UNMAPPED);
        *r2.reference_sequence_id_mut() = Some(0);
        *r2.data_mut() = [(mc_tag, BufValue::String("100M".into()))].into_iter().collect();

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // Both should have reference cleared (None)
        assert!(template.records[0].reference_sequence_id().is_none(), "R1 ref should be cleared");
        assert!(template.records[1].reference_sequence_id().is_none(), "R2 ref should be cleared");

        // Both should have mate_unmapped flag set
        assert!(template.records[0].flags().is_mate_unmapped(), "R1 mate_unmapped");
        assert!(template.records[1].flags().is_mate_unmapped(), "R2 mate_unmapped");

        // MQ and MC tags should be removed
        assert!(template.records[0].data().get(&mq_tag).is_none(), "R1 MQ should be removed");
        assert!(template.records[1].data().get(&mc_tag).is_none(), "R2 MC should be removed");

        // TLEN should be 0
        assert_eq!(template.records[0].template_length(), 0, "R1 TLEN should be 0");
        assert_eq!(template.records[1].template_length(), 0, "R2 TLEN should be 0");

        Ok(())
    }

    /// Test `fix_mate_info` correctly handles one mapped, one unmapped case
    #[test]
    fn test_fix_mate_info_one_unmapped() -> Result<()> {
        let mq_tag = Tag::new(b'M', b'Q');
        let mc_tag = Tag::new(b'M', b'C');

        // R1 is mapped
        let r1 = create_mapped_record(b"read1", FLAG_PAIRED | FLAG_READ1, 100, 30);

        // R2 is unmapped
        let r2 = create_test_record(b"read1", FLAG_PAIRED | FLAG_READ2 | FLAG_UNMAPPED);

        let mut builder = Builder::new();
        builder.push(r1)?;
        builder.push(r2)?;
        let mut template = builder.build()?;

        template.fix_mate_info()?;

        // Unmapped R2 should be placed at R1's position
        assert_eq!(
            template.records[1].reference_sequence_id(),
            Some(0),
            "R2 should be placed at R1's reference"
        );
        assert_eq!(
            template.records[1].alignment_start().map(|p| p.get()),
            Some(100),
            "R2 should be placed at R1's position"
        );

        // R1 (mapped) should NOT have MQ tag (mate is unmapped)
        assert!(
            template.records[0].data().get(&mq_tag).is_none(),
            "R1 should NOT have MQ when mate is unmapped"
        );

        // R2 (unmapped) SHOULD have MQ tag (mate is mapped)
        assert!(
            template.records[1].data().get(&mq_tag).is_some(),
            "R2 should have MQ when mate is mapped"
        );

        // R1 should have mate_unmapped flag set
        assert!(template.records[0].flags().is_mate_unmapped(), "R1 mate_unmapped should be set");

        // R2 should NOT have mate_unmapped flag set (R1 is mapped)
        assert!(
            !template.records[1].flags().is_mate_unmapped(),
            "R2 mate_unmapped should NOT be set"
        );

        // TLEN should be 0 for both
        assert_eq!(template.records[0].template_length(), 0, "R1 TLEN");
        assert_eq!(template.records[1].template_length(), 0, "R2 TLEN");

        // MC tag: R1 should not have it (mate unmapped), R2 should have it (mate mapped)
        assert!(
            template.records[0].data().get(&mc_tag).is_none(),
            "R1 should NOT have MC when mate is unmapped"
        );
        assert!(
            template.records[1].data().get(&mc_tag).is_some(),
            "R2 should have MC when mate is mapped"
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
        let mut r1 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | FLAG_MATE_REVERSE,
            100,
            30,
            200,
            Some(1), // mate on different chromosome
            Some(200),
        );
        // Override R1's reference to be chr 0
        *r1.reference_sequence_id_mut() = Some(0);

        let mut r2 = create_mapped_record_with_flags(
            b"read1",
            FLAG_PAIRED | FLAG_READ2 | FLAG_REVERSE,
            200,
            30,
            -200,
            Some(0),
            Some(100),
        );
        // R2 is on chr 1
        *r2.reference_sequence_id_mut() = Some(1);

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
    // Tests for RecordBuf accessor panic in raw-byte mode
    // ========================================================================

    /// Helper to create a minimal raw BAM record for testing raw-byte mode.
    #[allow(clippy::cast_possible_truncation)]
    fn make_minimal_raw_bam(name: &[u8], flags: u16) -> Vec<u8> {
        let l_read_name = (name.len() + 1) as u8; // +1 for null terminator
        let total = 32 + l_read_name as usize; // minimal: header + name only
        let mut buf = vec![0u8; total];

        buf[8] = l_read_name;
        buf[14..16].copy_from_slice(&flags.to_le_bytes());

        let name_start = 32;
        buf[name_start..name_start + name.len()].copy_from_slice(name);
        buf[name_start + name.len()] = 0; // null terminator

        buf
    }

    #[test]
    fn test_from_raw_records_creates_raw_mode_template() {
        let raw = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ1);
        let template = Template::from_raw_records(vec![raw]).unwrap();

        // Verify it's in raw-byte mode
        assert!(template.is_raw_byte_mode());
        // Verify records vec is empty (all data is in raw_records)
        assert!(template.records.is_empty());
        // Verify raw accessors work
        assert!(template.raw_r1().is_some());
    }

    #[test]
    fn test_r1_accessor_returns_none_in_raw_mode() {
        // After fix: calling r1() on a raw-mode Template returns None instead of panicking
        let raw = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ1);
        let template = Template::from_raw_records(vec![raw]).unwrap();

        // r1() should return None in raw mode (use raw_r1() instead)
        assert!(template.r1().is_none());
        // raw_r1() should still work
        assert!(template.raw_r1().is_some());
    }

    #[test]
    fn test_r2_accessor_returns_none_in_raw_mode() {
        // After fix: calling r2() on a raw-mode Template returns None instead of panicking
        let r1 = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r2 = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ2);
        let template = Template::from_raw_records(vec![r1, r2]).unwrap();

        assert!(template.is_raw_byte_mode());

        // r2() should return None in raw mode (use raw_r2() instead)
        assert!(template.r2().is_none());
        // raw_r2() should still work
        assert!(template.raw_r2().is_some());
    }

    // ========================================================================
    // from_raw_records fast path tests
    // ========================================================================

    #[test]
    fn test_from_raw_records_fast_path_normal_order() {
        // R1 first, R2 second — fast path with no swap
        let r1 = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r2 = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ2);
        let template = Template::from_raw_records(vec![r1, r2]).unwrap();

        assert!(template.is_raw_byte_mode());
        assert!(template.raw_r1().is_some());
        assert!(template.raw_r2().is_some());
        // Verify R1 is at index 0 (has FIRST_SEGMENT flag)
        let r1_flags = crate::sort::bam_fields::flags(template.raw_r1().unwrap());
        assert_ne!(r1_flags & FLAG_READ1, 0);
    }

    #[test]
    fn test_from_raw_records_fast_path_swap() {
        // R2 first, R1 second — fast path should swap them
        let r1 = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r2 = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ2);
        let template = Template::from_raw_records(vec![r2, r1]).unwrap();

        assert!(template.is_raw_byte_mode());
        // After swap, R1 should be at index 0
        let r1_flags = crate::sort::bam_fields::flags(template.raw_r1().unwrap());
        assert_ne!(r1_flags & FLAG_READ1, 0);
        let r2_flags = crate::sort::bam_fields::flags(template.raw_r2().unwrap());
        assert_ne!(r2_flags & FLAG_READ2, 0);
    }

    #[test]
    fn test_from_raw_records_fast_path_fallthrough_both_r1() {
        // Both records are R1 — fast path should fall through to general path (error)
        let r1a = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r1b = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ1);
        let result = Template::from_raw_records(vec![r1a, r1b]);
        assert!(result.is_err());
    }

    #[test]
    fn test_from_raw_records_fast_path_with_secondary() {
        // 2 records but one is secondary — should skip fast path
        let r1 = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ1);
        let sec = make_minimal_raw_bam(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | crate::sort::bam_fields::flags::SECONDARY,
        );
        let template = Template::from_raw_records(vec![r1, sec]).unwrap();
        assert!(template.is_raw_byte_mode());
        assert!(template.raw_r1().is_some());
    }

    #[test]
    fn test_from_raw_records_with_supplementary() {
        // R1 primary + R1 supplementary + R2 primary — general path
        let r1 = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r1_supp = make_minimal_raw_bam(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | crate::sort::bam_fields::flags::SUPPLEMENTARY,
        );
        let r2 = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ2);
        let template = Template::from_raw_records(vec![r1, r1_supp, r2]).unwrap();
        assert!(template.is_raw_byte_mode());
        assert!(template.raw_r1().is_some());
        assert!(template.raw_r2().is_some());
    }

    #[test]
    fn test_from_raw_records_single_unpaired() {
        // Single unpaired record — should be treated as R1
        let r = make_minimal_raw_bam(b"read1", 0); // no PAIRED flag
        let template = Template::from_raw_records(vec![r]).unwrap();
        assert!(template.is_raw_byte_mode());
        assert!(template.raw_r1().is_some());
        assert!(template.raw_r2().is_none());
    }

    #[test]
    fn test_from_raw_records_empty() {
        let result = Template::from_raw_records(vec![]);
        assert!(result.is_err());
    }

    #[test]
    fn test_from_raw_records_name_mismatch() {
        // Two records with different names — should error in general path
        let r1 = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r1_supp = make_minimal_raw_bam(
            b"read1",
            FLAG_PAIRED | FLAG_READ1 | crate::sort::bam_fields::flags::SUPPLEMENTARY,
        );
        let r2_wrong = make_minimal_raw_bam(b"read2", FLAG_PAIRED | FLAG_READ2);
        let result = Template::from_raw_records(vec![r1, r1_supp, r2_wrong]);
        assert!(result.is_err());
    }

    #[test]
    fn test_from_raw_records_all_raw_records_mut() {
        let r1 = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r2 = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ2);
        let mut template = Template::from_raw_records(vec![r1, r2]).unwrap();
        // Should be able to get mutable access
        let recs = template.all_raw_records_mut().unwrap();
        assert_eq!(recs.len(), 2);
    }

    #[test]
    fn test_from_raw_records_into_raw_records() {
        let r1 = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ1);
        let r2 = make_minimal_raw_bam(b"read1", FLAG_PAIRED | FLAG_READ2);
        let template = Template::from_raw_records(vec![r1, r2]).unwrap();
        let recs = template.into_raw_records().unwrap();
        assert_eq!(recs.len(), 2);
    }

    #[test]
    fn test_from_raw_records_truncated_header() {
        // Record too short (< 32 bytes)
        let short = vec![0u8; 20];
        let err = Template::from_raw_records(vec![short]).unwrap_err();
        assert!(err.to_string().contains("too short"), "Error: {err}");
    }

    #[test]
    fn test_from_raw_records_truncated_read_name() {
        // Record has 32 bytes but l_read_name claims more bytes than available
        let mut buf = vec![0u8; 34]; // 32 header + 2 bytes
        buf[8] = 10; // l_read_name=10, but only 2 bytes after header
        let err = Template::from_raw_records(vec![buf]).unwrap_err();
        assert!(err.to_string().contains("truncated"), "Error: {err}");
    }

    #[test]
    fn test_from_raw_records_valid_l_read_name() {
        // Record with l_read_name that exactly fits
        let rec = make_minimal_raw_bam(b"test", FLAG_PAIRED | FLAG_READ1);
        assert!(Template::from_raw_records(vec![rec]).is_ok());
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
