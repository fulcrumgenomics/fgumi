use anyhow::Result;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::fs;
use std::path::Path;

use crate::generate_metrics::MetricPage;
use crate::generate_tools::ToolPage;

/// Generate `SUMMARY.md` for mdBook from static guide pages and dynamic tool/metric pages.
#[allow(clippy::too_many_lines)]
pub fn generate(
    docs_src: &Path,
    tool_pages: &[ToolPage],
    metric_pages: &[MetricPage],
) -> Result<()> {
    let mut md = String::new();

    md.push_str("# Summary\n\n");
    md.push_str("[Home](index.md)\n\n");

    // ── User Guide ────────────────────────────────────────────────────────────
    md.push_str("# User Guide\n\n");

    // Top-level entry point
    let entry = ("Getting Started", "guide/getting-started.md");
    if docs_src.join(entry.1).exists() {
        let _ = writeln!(md, "- [{}]({})", entry.0, entry.1);
    }

    // Core Concepts sub-group
    let core_concepts = [
        ("Read Structures", "guide/read-structures.md"),
        ("UMI Grouping", "guide/umi-grouping.md"),
        ("Tracking Reads", "guide/tracking-reads.md"),
    ];
    let core_exists = core_concepts.iter().any(|(_, p)| docs_src.join(p).exists());
    if core_exists {
        md.push_str("- [Core Concepts]()\n");
        for (title, path) in &core_concepts {
            if docs_src.join(path).exists() {
                let _ = writeln!(md, "  - [{title}]({path})");
            }
        }
    }

    // Consensus sub-group
    let consensus = [
        ("Consensus Calling", "guide/consensus-calling.md"),
        ("Duplex Consensus Calling", "guide/duplex-consensus-calling.md"),
    ];
    let consensus_exists = consensus.iter().any(|(_, p)| docs_src.join(p).exists());
    if consensus_exists {
        md.push_str("- [Consensus Calling]()\n");
        for (title, path) in &consensus {
            if docs_src.join(path).exists() {
                let _ = writeln!(md, "  - [{title}]({path})");
            }
        }
    }

    // Advanced Topics sub-group
    let advanced = [
        ("Best Practices", "guide/best-practices.md"),
        ("Performance Tuning", "guide/performance-tuning.md"),
        ("Working with Metrics", "guide/working-with-metrics.md"),
        ("Migration from fgbio", "guide/migration-from-fgbio.md"),
    ];
    let advanced_exists = advanced.iter().any(|(_, p)| docs_src.join(p).exists());
    if advanced_exists {
        md.push_str("- [Advanced Topics]()\n");
        for (title, path) in &advanced {
            if docs_src.join(path).exists() {
                let _ = writeln!(md, "  - [{title}]({path})");
            }
        }
    }

    md.push('\n');

    // ── Tool Reference ────────────────────────────────────────────────────────
    md.push_str("# Tool Reference\n\n");
    md.push_str("- [Index](tools/README.md)\n");

    let mut by_category: BTreeMap<String, Vec<&ToolPage>> = BTreeMap::new();
    for page in tool_pages {
        by_category.entry(page.category.clone()).or_default().push(page);
    }

    // Pipeline order: merge GROUP + DEDUP into one sidebar section
    let pipeline: &[(&str, &[&str])] = &[
        ("UMI Extraction", &["UMI EXTRACTION"]),
        ("Alignment", &["ALIGNMENT"]),
        ("Grouping & Deduplication", &["GROUP", "DEDUP"]),
        ("Consensus Calling", &["CONSENSUS"]),
        ("Post-Consensus", &["POST-CONSENSUS"]),
        ("Utilities", &["UTILITIES"]),
    ];

    for (display_name, source_categories) in pipeline {
        let tools: Vec<&&ToolPage> = source_categories
            .iter()
            .filter_map(|cat| by_category.get(*cat))
            .flat_map(|v| v.iter())
            .collect();
        if !tools.is_empty() {
            let _ = writeln!(md, "- [{display_name}]()");
            for tool in tools {
                let _ = writeln!(md, "  - [{}]({})", tool.name, tool.path);
            }
        }
    }

    // Any categories not covered by the pipeline order
    let covered: std::collections::HashSet<&str> =
        pipeline.iter().flat_map(|(_, cats)| cats.iter().copied()).collect();
    for (category, tools) in &by_category {
        if !covered.contains(category.as_str()) {
            let _ = writeln!(md, "- [{category}]()");
            for tool in tools {
                let _ = writeln!(md, "  - [{}]({})", tool.name, tool.path);
            }
        }
    }

    md.push('\n');

    // ── Metrics Reference ─────────────────────────────────────────────────────
    md.push_str("# Metrics Reference\n\n");
    md.push_str("- [Index](metrics/README.md)\n");

    // Group metrics by prefix into logical sections
    let metric_groups: &[(&str, &[&str])] = &[
        ("UMI", &["Umi"]),
        ("Grouping", &["FamilySize", "PositionGroup"]),
        ("Duplex", &["Duplex"]),
        ("Simplex", &["Simplex"]),
        ("Consensus", &["Consensus"]),
        ("Post-Processing", &["Clipping"]),
    ];

    let mut assigned: std::collections::HashSet<&str> = std::collections::HashSet::new();

    for (group_name, prefixes) in metric_groups {
        let matches: Vec<&MetricPage> = metric_pages
            .iter()
            .filter(|p| prefixes.iter().any(|pfx| p.name.starts_with(pfx)))
            .collect();
        if !matches.is_empty() {
            let _ = writeln!(md, "- [{group_name}]()");
            for page in &matches {
                let _ = writeln!(md, "  - [{}]({})", page.name, page.path);
                assigned.insert(page.name.as_str());
            }
        }
    }

    // Any metrics not matched by a group
    let unmatched: Vec<&MetricPage> =
        metric_pages.iter().filter(|p| !assigned.contains(p.name.as_str())).collect();
    if !unmatched.is_empty() {
        md.push_str("- [Other]()\n");
        for page in unmatched {
            let _ = writeln!(md, "  - [{}]({})", page.name, page.path);
        }
    }

    md.push('\n');

    fs::write(docs_src.join("SUMMARY.md"), md)?;
    Ok(())
}
