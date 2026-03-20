use anyhow::Result;
use std::path::Path;

use crate::generate_metrics;
use crate::generate_summary;
use crate::generate_tools;

/// Generate all documentation markdown files.
///
/// Creates auto-generated tool reference pages, metric reference pages, and the
/// SUMMARY.md table of contents in the mdBook source directory.
pub fn generate() -> Result<()> {
    let docs_src = Path::new("docs/src");

    // Generate tool reference pages
    let tool_pages = generate_tools::generate(docs_src)?;
    eprintln!("Generated {} tool pages", tool_pages.len());

    // Generate metric reference pages
    let metric_pages = generate_metrics::generate(docs_src)?;
    eprintln!("Generated {} metric pages", metric_pages.len());

    // Generate SUMMARY.md
    generate_summary::generate(docs_src, &tool_pages, &metric_pages)?;
    eprintln!("Generated SUMMARY.md");

    Ok(())
}
