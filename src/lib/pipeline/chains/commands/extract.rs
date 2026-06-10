//! Chain builder support for [`crate::pipeline::chains::Stage::Extract`].
//!
//! Provides the extract-specific finalize hook, the `build_fastq_header`
//! helper that synthesizes an unmapped-BAM `@HD`/`@RG`/`@CO` header from
//! [`ExtractOptions`], and the `build_extract_chain` delegate.
//!
//! The header builder reproduces the same `@HD`, `@RG`, and `@CO` records
//! that [`Extract::create_header`] emits, minus the `@PG` record —
//! [`ChainBuilder::new`] adds that uniformly for all chains.
//!
//! [`ExtractOptions`]: crate::commands::extract::ExtractOptions
//! [`Extract::create_header`]: crate::commands::extract::Extract
//! [`ChainBuilder::new`]: crate::pipeline::chains::builder::ChainBuilder

use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Result;
use log::info;
use noodles::sam::header::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::map::builder::Builder;
use noodles::sam::header::record::value::map::header::group_order;
use noodles::sam::header::record::value::map::header::sort_order;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use noodles::sam::header::record::value::{
    Map as HeaderRecordMap,
    map::{Header as HeaderRecord, Tag as HeaderTag},
};

use crate::commands::extract::ExtractOptions;
use crate::logging::OperationTimer;
use crate::pipeline::chains::FinalizeHook;

// ─────────────────────────────────────────────────────────────────────────────
// ExtractFinalizeHook
// ─────────────────────────────────────────────────────────────────────────────

/// Post-pipeline finalize hook for the extract stage.
///
/// Logs the records-emitted count and calls `timer.log_completion` so the
/// per-stage wall-time line appears in the log.
pub(crate) struct ExtractFinalizeHook {
    /// Counter atomically incremented by `ExtractStep` as it emits
    /// `BamTemplateBatch`es downstream.
    pub(crate) records_emitted: Arc<AtomicU64>,
    pub(crate) timer: OperationTimer,
}

impl FinalizeHook for ExtractFinalizeHook {
    fn finalize(self: Box<Self>) -> Result<()> {
        let ExtractFinalizeHook { records_emitted, timer } = *self;
        let n = records_emitted.load(Ordering::Relaxed);
        info!("Extract: completed ({n} records emitted)");
        timer.log_completion(n);
        Ok(())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// build_fastq_header — reproduce Extract::create_header sans @PG
// ─────────────────────────────────────────────────────────────────────────────

/// Conditionally add a read-group tag value to a [`Builder<ReadGroup>`].
///
/// If `value` is `Some`, inserts the tag with the value. Otherwise returns
/// the builder unchanged. Mirrors [`Extract::add_to_read_group`].
fn add_to_read_group(
    rg: Builder<ReadGroup>,
    tag: noodles::sam::header::record::value::map::tag::Other<rg_tag::Standard>,
    value: Option<&String>,
) -> Builder<ReadGroup> {
    if let Some(v) = value { rg.insert(tag, v.clone()) } else { rg }
}

/// Build a SAM header for an unmapped BAM produced from FASTQ input.
///
/// The returned header contains:
/// - `@HD` with `VN:1.6`, `SO:unsorted`, `GO:query`.
/// - `@RG` with the read-group fields from `extract_opts` (SM, LB, PL, PU).
/// - `@CO` comment lines from `extract_opts.comments`.
///
/// The header does **not** include a `@PG` record — [`ChainBuilder::new`]
/// injects that uniformly for every chain.
///
/// This reproduces the header structure of [`Extract::create_header`] so the
/// two code paths produce byte-identical headers (modulo `@PG`).
///
/// [`ChainBuilder::new`]: crate::pipeline::chains::builder::ChainBuilder
/// [`Extract::create_header`]: crate::commands::extract::Extract
pub(crate) fn build_fastq_header(extract_opts: &ExtractOptions) -> Result<Header> {
    let mut header = Header::builder();

    // @HD: sort order = unsorted, group order = query.
    let HeaderTag::Other(so_tag) = HeaderTag::from([b'S', b'O']) else { unreachable!() };
    let HeaderTag::Other(go_tag) = HeaderTag::from([b'G', b'O']) else { unreachable!() };
    let map = HeaderRecordMap::<HeaderRecord>::builder()
        .insert(so_tag, sort_order::UNSORTED)
        .insert(go_tag, group_order::QUERY)
        .build()?;
    header = header.set_header(map);

    // @CO comment lines.
    for comment in &extract_opts.comments {
        header = header.add_comment(comment.clone());
    }

    // @RG read group. SM and LB are required; all other tags are optional.
    let mut rg = Map::<ReadGroup>::builder();
    rg = add_to_read_group(rg, rg_tag::SAMPLE, Some(&extract_opts.sample));
    rg = add_to_read_group(rg, rg_tag::LIBRARY, Some(&extract_opts.library));
    rg = add_to_read_group(rg, rg_tag::BARCODE, extract_opts.barcode.as_ref());
    rg = add_to_read_group(rg, rg_tag::PLATFORM, extract_opts.platform.as_ref());
    rg = add_to_read_group(rg, rg_tag::PLATFORM_UNIT, extract_opts.platform_unit.as_ref());
    rg = add_to_read_group(rg, rg_tag::PLATFORM_MODEL, extract_opts.platform_model.as_ref());
    rg = add_to_read_group(rg, rg_tag::SEQUENCING_CENTER, extract_opts.sequencing_center.as_ref());
    rg = add_to_read_group(
        rg,
        rg_tag::PREDICTED_MEDIAN_INSERT_SIZE,
        extract_opts.predicted_insert_size.map(|i| i.to_string()).as_ref(),
    );
    rg = add_to_read_group(rg, rg_tag::DESCRIPTION, extract_opts.description.as_ref());
    rg = add_to_read_group(rg, rg_tag::PRODUCED_AT, extract_opts.run_date.as_ref());

    header = header.add_read_group(extract_opts.read_group_id.clone(), rg.build()?);

    Ok(header.build())
}

// ─────────────────────────────────────────────────────────────────────────────
// build_extract_chain — 10-line delegate
// ─────────────────────────────────────────────────────────────────────────────

/// Build a [`BuiltPipeline`] for an extract-only chain.
///
/// `spec.stages == [Stage::Extract]`. Delegates to [`ChainBuilder`].
///
/// # Errors
///
/// Returns validation errors (missing extract options, invalid source type)
/// or any underlying pipeline construction error.
///
/// [`BuiltPipeline`]: crate::pipeline::chains::BuiltPipeline
/// [`ChainBuilder`]: crate::pipeline::chains::builder::ChainBuilder
#[allow(clippy::needless_pass_by_value)]
pub fn build_extract_chain(
    spec: crate::pipeline::chains::ChainSpec,
) -> Result<crate::pipeline::chains::BuiltPipeline> {
    use crate::pipeline::chains::{Stage, builder::StagePosition};
    let mut chain = crate::pipeline::chains::builder::ChainBuilder::new(&spec)?;
    chain.add_source()?;
    chain.add_stage(Stage::Extract, StagePosition::Terminal)?;
    chain.add_sink()?;
    chain.build()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::extract::QualityEncoding;
    use bstr::BString;

    /// Build a minimal `ExtractOptions` for tests.
    fn test_extract_opts() -> ExtractOptions {
        ExtractOptions {
            sample: "TestSample".to_string(),
            library: "TestLibrary".to_string(),
            platform: Some("illumina".to_string()),
            platform_unit: Some("HXXXXX.1.ATCACG".to_string()),
            read_group_id: "A".to_string(),
            comments: vec!["comment one".to_string(), "comment two".to_string()],
            barcode: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            run_date: None,
            quality_encoding: QualityEncoding::Standard,
            store_umi_quals: false,
            store_cell_quals: false,
            single_tag: None,
            annotate_read_names: false,
            extract_umis_from_read_names: false,
            store_sample_barcode_qualities: false,
            async_reader: false,
        }
    }

    #[test]
    fn build_fastq_header_has_expected_hd() {
        let opts = test_extract_opts();
        let header = build_fastq_header(&opts).expect("build_fastq_header");
        let hd = header.header().expect("@HD should be present");
        // sort order = unsorted
        let HeaderTag::Other(so_tag) = HeaderTag::from([b'S', b'O']) else { unreachable!() };
        assert_eq!(
            hd.other_fields().get(&so_tag).map(std::string::ToString::to_string),
            Some("unsorted".to_string())
        );
        // group order = query
        let HeaderTag::Other(go_tag) = HeaderTag::from([b'G', b'O']) else { unreachable!() };
        assert_eq!(
            hd.other_fields().get(&go_tag).map(std::string::ToString::to_string),
            Some("query".to_string())
        );
    }

    #[test]
    fn build_fastq_header_has_expected_rg() {
        let opts = test_extract_opts();
        let header = build_fastq_header(&opts).expect("build_fastq_header");
        let rg_id = BString::from("A");
        let rg = header.read_groups().get(&rg_id).expect("@RG ID:A should be present");
        assert_eq!(
            rg.other_fields().get(&rg_tag::SAMPLE).map(std::string::ToString::to_string),
            Some("TestSample".to_string())
        );
        assert_eq!(
            rg.other_fields().get(&rg_tag::LIBRARY).map(std::string::ToString::to_string),
            Some("TestLibrary".to_string())
        );
        assert_eq!(
            rg.other_fields().get(&rg_tag::PLATFORM).map(std::string::ToString::to_string),
            Some("illumina".to_string())
        );
        assert_eq!(
            rg.other_fields().get(&rg_tag::PLATFORM_UNIT).map(std::string::ToString::to_string),
            Some("HXXXXX.1.ATCACG".to_string())
        );
    }

    #[test]
    fn build_fastq_header_has_expected_comments() {
        let opts = test_extract_opts();
        let header = build_fastq_header(&opts).expect("build_fastq_header");
        let comments: Vec<String> =
            header.comments().iter().map(std::string::ToString::to_string).collect();
        assert_eq!(comments, vec!["comment one", "comment two"]);
    }

    #[test]
    fn build_fastq_header_has_no_pg() {
        let opts = test_extract_opts();
        let header = build_fastq_header(&opts).expect("build_fastq_header");
        assert!(header.programs().as_ref().is_empty(), "@PG should not be present");
    }

    #[test]
    fn build_fastq_header_optional_fields_omitted() {
        let mut opts = test_extract_opts();
        opts.platform = None;
        opts.platform_unit = None;
        let header = build_fastq_header(&opts).expect("build_fastq_header");
        let rg_id = BString::from("A");
        let rg = header.read_groups().get(&rg_id).expect("@RG ID:A");
        assert!(rg.other_fields().get(&rg_tag::PLATFORM).is_none());
        assert!(rg.other_fields().get(&rg_tag::PLATFORM_UNIT).is_none());
    }
}
