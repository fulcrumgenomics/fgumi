use anyhow::Result;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::fs;
use std::path::Path;

use crate::generate_metrics::MetricPage;
use crate::generate_tools::ToolPage;

/// Generate `SUMMARY.md` for mdBook from static guide pages and dynamic tool/metric pages.
pub fn generate(
    docs_src: &Path,
    tool_pages: &[ToolPage],
    metric_pages: &[MetricPage],
) -> Result<()> {
    let mut md = String::new();

    md.push_str("# Summary\n\n");
    md.push_str("[Introduction](index.md)\n\n");

    // Guide section (static hand-written pages)
    md.push_str("# User Guide\n\n");
    let guide_pages = [
        ("Getting Started", "guide/getting-started.md"),
        ("Read Structures", "guide/read-structures.md"),
        ("UMI Grouping", "guide/umi-grouping.md"),
        ("Consensus Calling", "guide/consensus-calling.md"),
        ("Duplex Consensus Calling", "guide/duplex-consensus-calling.md"),
        ("Tracking Reads", "guide/tracking-reads.md"),
        ("Best Practices", "guide/best-practices.md"),
        ("Performance Tuning", "guide/performance-tuning.md"),
        ("Working with Metrics", "guide/working-with-metrics.md"),
        ("Migration from fgbio", "guide/migration-from-fgbio.md"),
    ];

    for (title, path) in &guide_pages {
        if docs_src.join(path).exists() {
            let _ = writeln!(md, "- [{title}]({path})");
        }
    }
    md.push('\n');

    // Tool reference section (auto-generated from clap)
    md.push_str("# Tool Reference\n\n");
    md.push_str("- [Overview](tools/README.md)\n");

    let mut by_category: BTreeMap<String, Vec<&ToolPage>> = BTreeMap::new();
    for page in tool_pages {
        by_category.entry(page.category.clone()).or_default().push(page);
    }

    for (category, tools) in &by_category {
        let _ = writeln!(md, "  - [{category}]()");
        for tool in tools {
            let _ = writeln!(md, "    - [{}]({})", tool.name, tool.path);
        }
    }
    md.push('\n');

    // Metrics reference section (auto-generated from syn)
    md.push_str("# Metrics Reference\n\n");
    md.push_str("- [Overview](metrics/README.md)\n");
    for page in metric_pages {
        let _ = writeln!(md, "  - [{}]({})", page.name, page.path);
    }
    md.push('\n');

    fs::write(docs_src.join("SUMMARY.md"), md)?;
    Ok(())
}
