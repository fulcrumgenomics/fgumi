use anyhow::Result;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::fs;
use std::path::Path;

use crate::generate_metrics::MetricPage;
use crate::generate_tools::ToolPage;

/// Emit one collapsible User Guide sub-group: a bracketed section header
/// followed by a `guide/<page>.md` link for each entry whose page exists on
/// disk. The whole group (header included) is omitted when none of its pages
/// exist, so a stripped-down docs tree never renders an empty section.
fn push_group(md: &mut String, docs_src: &Path, header: &str, entries: &[(&str, &str)]) {
    if entries.iter().any(|(_, path)| docs_src.join(path).exists()) {
        let _ = writeln!(md, "- [{header}]()");
        for (title, path) in entries {
            if docs_src.join(path).exists() {
                let _ = writeln!(md, "  - [{title}]({path})");
            }
        }
    }
}

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

    // Running Pipelines — fgumi runall vs individual commands. A standalone
    // top-level leaf, listed right after Getting Started.
    let running = ("Running Pipelines", "guide/running-pipelines.md");
    if docs_src.join(running.1).exists() {
        let _ = writeln!(md, "- [{}]({})", running.0, running.1);
    }

    // Sub-groups: each renders a collapsible sidebar section, omitted entirely
    // when none of its pages exist. Standalone leaf pages (above) are listed
    // before these so the sidebar's plain links and expandable section headers
    // form two visually distinct zones.
    push_group(
        &mut md,
        docs_src,
        "Core Concepts",
        &[
            ("Read Structures", "guide/read-structures.md"),
            ("UMI Grouping", "guide/umi-grouping.md"),
            ("Tracking Reads", "guide/tracking-reads.md"),
        ],
    );
    push_group(
        &mut md,
        docs_src,
        "Consensus Calling",
        &[
            ("Consensus Calling", "guide/consensus-calling.md"),
            ("Duplex Consensus Calling", "guide/duplex-consensus-calling.md"),
        ],
    );
    push_group(&mut md, docs_src, "Methylation", &[("Pipeline Guide", "guide/methylation.md")]);
    push_group(
        &mut md,
        docs_src,
        "NanoSeq (Duplex-Seq)",
        &[("Pipeline Guide", "guide/nanoseq.md")],
    );
    push_group(
        &mut md,
        docs_src,
        "Advanced Topics",
        &[
            ("Best Practices", "guide/best-practices.md"),
            ("Performance Tuning", "guide/performance-tuning.md"),
            ("Working with Metrics", "guide/working-with-metrics.md"),
            ("Migration from fgbio", "guide/migration-from-fgbio.md"),
        ],
    );
    push_group(
        &mut md,
        docs_src,
        "Reference",
        &[("Glossary", "guide/glossary.md"), ("Troubleshooting", "guide/troubleshooting.md")],
    );

    md.push('\n');

    // ── Tool Reference ────────────────────────────────────────────────────────
    md.push_str("# Tool Reference\n\n");
    md.push_str("- [Index](tools/index.md)\n");

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
    md.push_str("- [Index](metrics/index.md)\n");

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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    /// Generate `SUMMARY.md` for a docs tree that either has or lacks the
    /// `NanoSeq` guide page, and return the markdown it wrote.
    fn summary_with_nanoseq(present: bool) -> String {
        let tmp = TempDir::new().expect("failed to create temp dir");
        let docs_src = tmp.path();
        let guide_dir = docs_src.join("guide");
        fs::create_dir_all(&guide_dir).expect("failed to create guide dir");
        if present {
            fs::write(guide_dir.join("nanoseq.md"), "# NanoSeq\n").expect("failed to write guide");
        }
        generate(docs_src, &[], &[]).expect("generate failed");
        fs::read_to_string(docs_src.join("SUMMARY.md")).expect("failed to read SUMMARY.md")
    }

    #[test]
    fn nanoseq_group_present_when_guide_exists() {
        let md = summary_with_nanoseq(true);
        assert!(md.contains("- [NanoSeq (Duplex-Seq)]()"), "missing NanoSeq group:\n{md}");
        assert!(
            md.contains("  - [Pipeline Guide](guide/nanoseq.md)"),
            "missing NanoSeq guide link:\n{md}"
        );
    }

    #[test]
    fn nanoseq_group_absent_when_guide_missing() {
        let md = summary_with_nanoseq(false);
        assert!(
            !md.contains("NanoSeq"),
            "NanoSeq group should be omitted when guide is missing:\n{md}"
        );
        assert!(
            !md.contains("guide/nanoseq.md"),
            "NanoSeq link should be omitted when guide is missing:\n{md}"
        );
    }

    /// Generate `SUMMARY.md` for a docs tree that optionally contains the
    /// Running Pipelines leaf and the Reference-group guide pages, and return
    /// the markdown it wrote.
    fn summary_with_pages(pages: &[&str]) -> String {
        let tmp = TempDir::new().expect("failed to create temp dir");
        let docs_src = tmp.path();
        let guide_dir = docs_src.join("guide");
        fs::create_dir_all(&guide_dir).expect("failed to create guide dir");
        for page in pages {
            fs::write(guide_dir.join(page), "# Page\n").expect("failed to write guide");
        }
        generate(docs_src, &[], &[]).expect("generate failed");
        fs::read_to_string(docs_src.join("SUMMARY.md")).expect("failed to read SUMMARY.md")
    }

    #[test]
    fn running_pipelines_leaf_present_when_guide_exists() {
        let md = summary_with_pages(&["running-pipelines.md"]);
        assert!(
            md.contains("- [Running Pipelines](guide/running-pipelines.md)"),
            "missing Running Pipelines leaf:\n{md}"
        );
    }

    #[test]
    fn running_pipelines_leaf_absent_when_guide_missing() {
        let md = summary_with_pages(&[]);
        assert!(
            !md.contains("guide/running-pipelines.md"),
            "Running Pipelines link should be omitted when guide is missing:\n{md}"
        );
    }

    #[test]
    fn reference_group_present_when_a_guide_exists() {
        let md = summary_with_pages(&["glossary.md", "troubleshooting.md"]);
        assert!(md.contains("- [Reference]()"), "missing Reference group:\n{md}");
        assert!(md.contains("  - [Glossary](guide/glossary.md)"), "missing Glossary link:\n{md}");
        assert!(
            md.contains("  - [Troubleshooting](guide/troubleshooting.md)"),
            "missing Troubleshooting link:\n{md}"
        );
    }

    #[test]
    fn reference_group_present_with_only_one_guide() {
        // Pins the `.any()` group guard and the per-entry existence filter: with
        // only Glossary present the group still renders, but the missing
        // Troubleshooting link is omitted (an `.all()` group guard, or an
        // unconditional per-entry emit, would fail this).
        let md = summary_with_pages(&["glossary.md"]);
        assert!(
            md.contains("- [Reference]()"),
            "Reference group should render with one guide:\n{md}"
        );
        assert!(md.contains("  - [Glossary](guide/glossary.md)"), "missing Glossary link:\n{md}");
        assert!(
            !md.contains("guide/troubleshooting.md"),
            "Troubleshooting link should be omitted when its guide is missing:\n{md}"
        );
    }

    #[test]
    fn reference_group_absent_when_guides_missing() {
        let md = summary_with_pages(&[]);
        assert!(
            !md.contains("- [Reference]()"),
            "Reference group should be omitted when its guides are missing:\n{md}"
        );
    }
}
