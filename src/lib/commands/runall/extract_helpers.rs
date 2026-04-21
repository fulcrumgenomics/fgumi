//! Extract-path config builders consumed by the `runall` planner.

use std::str::FromStr;

use crate::extract::{
    self, ExtractParams, QualityEncoding, ReadGroupMetadata, build_unmapped_bam_header,
};
use anyhow::{Result, anyhow};
use fgumi_simd_fastq::SimdFastqReader;
use noodles::sam::Header;
use read_structure::ReadStructure;

use super::Runall;

impl Runall {
    /// Build an unmapped BAM header for the extract stage.
    ///
    /// Creates a SAM header with sort order, read group (including all metadata
    /// fields), comments, and PG record, matching what `fgumi extract` would
    /// produce. Delegates the core header construction to
    /// [`build_unmapped_bam_header`] in the library.
    pub(crate) fn build_unmapped_header(&self, command_line: &str) -> Result<Header> {
        let sample =
            self.extract_opts.extract_sample.as_ref().expect("--sample validated in validate()");
        let library =
            self.extract_opts.extract_library.as_ref().expect("--library validated in validate()");

        let metadata = ReadGroupMetadata {
            platform: Some(self.extract_opts.extract_platform.clone()),
            platform_unit: self.extract_opts.extract_platform_unit.clone(),
            platform_model: self.extract_opts.extract_platform_model.clone(),
            barcode: self.extract_opts.extract_barcode.clone(),
            sequencing_center: self.extract_opts.extract_sequencing_center.clone(),
            predicted_insert_size: self.extract_opts.extract_predicted_insert_size,
            description: self.extract_opts.extract_description.clone(),
            run_date: self.extract_opts.extract_run_date.clone(),
        };

        let header = build_unmapped_bam_header(
            &self.consensus_opts.consensus_read_group_id,
            sample,
            library,
            &metadata,
        )?;

        let mut builder = Header::builder().set_header(
            header.header().expect("build_unmapped_bam_header always sets header record").clone(),
        );
        for (id, rg) in header.read_groups() {
            builder = builder.add_read_group(id.clone(), rg.clone());
        }
        if let Some(ref comments) = self.extract_opts.extract_comment {
            for comment in comments {
                builder = builder.add_comment(comment.clone());
            }
        }
        builder = crate::commands::common::add_pg_to_builder(builder, command_line)?;

        Ok(builder.build())
    }

    /// Build [`ExtractParams`] from the pipeline's CLI options.
    pub(crate) fn build_extract_params(&self) -> ExtractParams {
        ExtractParams {
            read_group_id: self.consensus_opts.consensus_read_group_id.clone(),
            umi_tag: *crate::sam::SamTag::RX,
            cell_tag: *crate::sam::SamTag::CB,
            umi_qual_tag: None,
            cell_qual_tag: None,
            single_tag: self
                .extract_opts
                .extract_single_tag
                .as_ref()
                .map(|t| t.as_bytes().try_into().unwrap_or(*b"BX")),
            annotate_read_names: self.extract_opts.extract_annotate_read_names,
            extract_umis_from_read_names: self.extract_opts.extract_extract_umis_from_read_names,
            store_sample_barcode_qualities: self
                .extract_opts
                .extract_store_sample_barcode_qualities,
        }
    }

    /// Parse read structures from the pipeline's extract options, defaulting to `+T`
    /// for each input FASTQ if none are specified.
    pub(crate) fn get_read_structures(&self) -> Result<Vec<ReadStructure>> {
        let rs_strings = self.extract_opts.extract_read_structures.as_deref().unwrap_or_default();
        if rs_strings.is_empty() {
            Ok(vec![ReadStructure::from_str("+T")?; self.input.len()])
        } else {
            rs_strings
                .iter()
                .map(|s| ReadStructure::from_str(s).map_err(anyhow::Error::from))
                .collect()
        }
    }

    /// Detect quality encoding by sampling the first N records from the first FASTQ.
    pub(crate) fn detect_quality_encoding(&self) -> Result<QualityEncoding> {
        let first_input =
            self.input.first().ok_or_else(|| anyhow!("No input FASTQ files provided"))?;

        let reader_box: Box<dyn std::io::BufRead + Send> = Box::new(std::io::BufReader::new(
            fgoxide::io::Io::new(5, 1024 * 1024).new_reader(first_input)?,
        ));
        let mut fq_reader = SimdFastqReader::with_capacity(reader_box, 1024 * 1024);

        let mut sample_quals = Vec::new();
        for _i in 0..extract::QUALITY_DETECTION_SAMPLE_SIZE {
            match fq_reader.next() {
                Some(Ok(rec)) => sample_quals.push(rec.quality),
                Some(Err(e)) => return Err(e.into()),
                None => break,
            }
        }

        QualityEncoding::detect(&sample_quals)
    }
}
