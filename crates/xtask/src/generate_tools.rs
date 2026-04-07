use anyhow::Result;
use clap::CommandFactory;
use std::collections::BTreeMap;
use std::fmt::Write;
use std::fs;
use std::path::Path;

/// A parsed tool page with its metadata.
pub struct ToolPage {
    /// Command name (e.g., "extract").
    pub name: String,
    /// Category (e.g., "UMI Extraction").
    pub category: String,
    /// Relative path within `docs/src` (e.g., "tools/extract.md").
    pub path: String,
}

/// Generate tool reference markdown pages from clap command introspection.
///
/// Uses `clap::CommandFactory` to walk the command tree from the fgumi
/// binary, extracting each subcommand's name, description, and arguments.
pub fn generate(docs_src: &Path) -> Result<Vec<ToolPage>> {
    let tools_dir = docs_src.join("tools");
    if tools_dir.exists() {
        fs::remove_dir_all(&tools_dir)?;
    }
    fs::create_dir_all(&tools_dir)?;

    let commands = collect_commands();

    let mut pages = Vec::new();
    let mut by_category: BTreeMap<String, Vec<(String, String, String)>> = BTreeMap::new();

    for (name, category, description, long_about, cmd) in &commands {
        let markdown = render_tool_page(name, category, description, long_about.as_deref(), cmd);
        let filename = format!("{name}.md");
        let path = tools_dir.join(&filename);
        fs::write(&path, &markdown)?;

        let rel_path = format!("tools/{filename}");
        pages.push(ToolPage {
            name: name.clone(),
            category: category.clone(),
            path: rel_path.clone(),
        });

        by_category.entry(category.clone()).or_default().push((
            name.clone(),
            rel_path,
            description.clone(),
        ));
    }

    let index = render_tools_index(&by_category);
    fs::write(tools_dir.join("README.md"), index)?;

    Ok(pages)
}

type CommandInfo = (String, String, String, Option<String>, clap::Command);

/// Collect all commands by building each command's clap parser individually.
fn collect_commands() -> Vec<CommandInfo> {
    #[cfg(feature = "compare")]
    use fgumi_lib::commands::compare;
    #[cfg(feature = "simulate")]
    use fgumi_lib::commands::simulate;
    use fgumi_lib::commands::{
        clip, codec, correct, dedup, downsample, duplex, duplex_metrics, extract, fastq, filter,
        group, merge, review, simplex, simplex_metrics, sort, zipper,
    };

    // Each command type with its CLI name. The name override is needed because
    // clap derives names from struct names, which may not match the subcommand
    // names used in the `Subcommand` enum.
    type CommandBuilder = (&'static str, fn() -> clap::Command);
    let parsers: Vec<CommandBuilder> = vec![
        ("extract", <extract::Extract as CommandFactory>::command),
        ("correct", <correct::CorrectUmis as CommandFactory>::command),
        ("fastq", <fastq::Fastq as CommandFactory>::command),
        ("zipper", <zipper::Zipper as CommandFactory>::command),
        ("sort", <sort::Sort as CommandFactory>::command),
        ("group", <group::GroupReadsByUmi as CommandFactory>::command),
        ("dedup", <dedup::MarkDuplicates as CommandFactory>::command),
        ("simplex", <simplex::Simplex as CommandFactory>::command),
        ("duplex", <duplex::Duplex as CommandFactory>::command),
        ("codec", <codec::Codec as CommandFactory>::command),
        ("filter", <filter::Filter as CommandFactory>::command),
        ("clip", <clip::Clip as CommandFactory>::command),
        ("duplex-metrics", <duplex_metrics::DuplexMetrics as CommandFactory>::command),
        ("review", <review::Review as CommandFactory>::command),
        ("downsample", <downsample::Downsample as CommandFactory>::command),
        ("merge", <merge::Merge as CommandFactory>::command),
        ("simplex-metrics", <simplex_metrics::SimplexMetrics as CommandFactory>::command),
        #[cfg(feature = "compare")]
        ("compare", <compare::Compare as CommandFactory>::command),
        #[cfg(feature = "simulate")]
        ("simulate", <simulate::Simulate as CommandFactory>::command),
    ];

    parsers
        .into_iter()
        .map(|(cli_name, build_cmd)| {
            let cmd = build_cmd();
            let name = cli_name.to_string();
            let about_raw = cmd.get_about().map(ToString::to_string).unwrap_or_default();
            let (category, description) = parse_about(&about_raw);
            let long_about = cmd.get_long_about().map(ToString::to_string);
            (name, category, description, long_about, cmd)
        })
        .collect()
}

/// Parse the ANSI-encoded `about` string into (category, description).
///
/// Format: `\x1b[38;5;NNNm[CATEGORY]\x1b[0m \x1b[36mDescription\x1b[0m`
fn parse_about(about: &str) -> (String, String) {
    let stripped = strip_ansi(about);
    if stripped.starts_with('[') {
        if let Some(end_bracket) = stripped.find(']') {
            let category = stripped[1..end_bracket].to_string();
            let description = stripped[end_bracket + 1..].trim().to_string();
            return (category, description);
        }
    }
    ("Uncategorized".to_string(), stripped)
}

/// Strip ANSI escape sequences from a string.
fn strip_ansi(s: &str) -> String {
    let mut result = String::with_capacity(s.len());
    let mut in_escape = false;
    for c in s.chars() {
        if c == '\x1b' {
            in_escape = true;
        } else if in_escape {
            if c.is_ascii_alphabetic() {
                in_escape = false;
            }
        } else {
            result.push(c);
        }
    }
    result
}

/// Escape HTML angle brackets in text, preserving markdown code spans and blocks.
fn escape_html_in_text(s: &str) -> String {
    let mut result = String::with_capacity(s.len());
    let mut in_code_block = false;
    for line in s.lines() {
        if line.starts_with("```") {
            in_code_block = !in_code_block;
            result.push_str(line);
        } else if in_code_block {
            result.push_str(line);
        } else {
            let mut in_code_span = false;
            for c in line.chars() {
                if c == '`' {
                    in_code_span = !in_code_span;
                    result.push(c);
                } else if !in_code_span && c == '<' {
                    result.push_str("&lt;");
                } else if !in_code_span && c == '>' {
                    result.push_str("&gt;");
                } else {
                    result.push(c);
                }
            }
        }
        result.push('\n');
    }
    result
}

/// Render a single tool's markdown page.
fn render_tool_page(
    name: &str,
    category: &str,
    description: &str,
    long_about: Option<&str>,
    cmd: &clap::Command,
) -> String {
    let mut md = String::new();

    let _ = writeln!(md, "# {name}\n");
    let _ = writeln!(md, "**Category:** {category}\n");
    let _ = writeln!(md, "{description}\n");

    if let Some(long) = long_about {
        md.push_str("## Description\n\n");
        md.push_str(&escape_html_in_text(&strip_ansi(long)));
        md.push('\n');
    }

    // Collect arguments (skip help and version which are auto-added)
    let args: Vec<_> =
        cmd.get_arguments().filter(|a| a.get_id() != "help" && a.get_id() != "version").collect();

    if !args.is_empty() {
        md.push_str("## Arguments\n\n");
        md.push_str("| Flag | Description | Default |\n");
        md.push_str("|------|-------------|--------|\n");

        for arg in &args {
            let flag = format_flag(arg);
            let help = arg
                .get_help()
                .map(|s| strip_ansi(&s.to_string()))
                .unwrap_or_default()
                .replace('|', "\\|")
                .replace('\n', " ")
                .replace('<', "&lt;")
                .replace('>', "&gt;");

            let default = format_default(arg);
            let _ = writeln!(md, "| `{flag}` | {help} | {default} |");
        }
        md.push('\n');
    }

    // Check for subcommands (e.g., compare has bams/metrics)
    let subs: Vec<_> = cmd.get_subcommands().collect();
    if !subs.is_empty() {
        md.push_str("## Subcommands\n\n");
        for sub in subs {
            let sub_name = sub.get_name();
            let sub_about = sub.get_about().map(|s| strip_ansi(&s.to_string())).unwrap_or_default();
            let _ = writeln!(md, "### {sub_name}\n");
            let _ = writeln!(md, "{sub_about}\n");

            let sub_args: Vec<_> = sub
                .get_arguments()
                .filter(|a| a.get_id() != "help" && a.get_id() != "version")
                .collect();

            if !sub_args.is_empty() {
                md.push_str("| Flag | Description | Default |\n");
                md.push_str("|------|-------------|--------|\n");
                for arg in &sub_args {
                    let flag = format_flag(arg);
                    let help = arg
                        .get_help()
                        .map(|s| strip_ansi(&s.to_string()))
                        .unwrap_or_default()
                        .replace('|', "\\|")
                        .replace('\n', " ")
                        .replace('<', "&lt;")
                        .replace('>', "&gt;");
                    let default = format_default(arg);
                    let _ = writeln!(md, "| `{flag}` | {help} | {default} |");
                }
                md.push('\n');
            }
        }
    }

    md
}

/// Format an argument's flag representation (e.g., `-i, --input <FILE>`).
fn format_flag(arg: &clap::Arg) -> String {
    let mut parts = Vec::new();
    if let Some(short) = arg.get_short() {
        parts.push(format!("-{short}"));
    }
    if let Some(long) = arg.get_long() {
        parts.push(format!("--{long}"));
    }
    if parts.is_empty() {
        parts.push(arg.get_id().to_string());
    }

    let flag = parts.join(", ");
    let value_names: Vec<_> = arg.get_value_names().unwrap_or_default().to_vec();
    if value_names.is_empty() {
        if arg.get_action().takes_values() {
            format!("{flag} <{}>", arg.get_id().to_string().to_uppercase())
        } else {
            flag
        }
    } else {
        let names = value_names.iter().map(|n| format!("<{n}>")).collect::<Vec<_>>().join(" ");
        format!("{flag} {names}")
    }
}

/// Format an argument's default value.
fn format_default(arg: &clap::Arg) -> String {
    let defaults: Vec<_> =
        arg.get_default_values().iter().map(|v| v.to_string_lossy().to_string()).collect();
    if defaults.is_empty() {
        if arg.is_required_set() { "**required**".to_string() } else { String::new() }
    } else {
        format!("`{}`", defaults.join(", "))
    }
}

/// Render the tools index page grouped by category.
fn render_tools_index(by_category: &BTreeMap<String, Vec<(String, String, String)>>) -> String {
    let mut md = String::new();
    md.push_str("# Tool Reference\n\n");
    md.push_str("Auto-generated from fgumi command definitions.\n\n");

    for (category, tools) in by_category {
        let _ = writeln!(md, "## {category}\n");
        md.push_str("| Command | Description |\n");
        md.push_str("|---------|-------------|\n");
        for (name, path, description) in tools {
            let display_path = path.strip_prefix("tools/").unwrap_or(path);
            let _ = writeln!(md, "| [`{name}`]({display_path}) | {description} |");
        }
        md.push('\n');
    }

    md
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strip_ansi() {
        let input = "\x1b[38;5;30m[UMI EXTRACTION]\x1b[0m \x1b[36mExtract UMIs\x1b[0m";
        assert_eq!(strip_ansi(input), "[UMI EXTRACTION] Extract UMIs");
    }

    #[test]
    fn test_parse_about() {
        let input = "\x1b[38;5;30m[UMI EXTRACTION]\x1b[0m \x1b[36mExtract UMIs from FASTQ\x1b[0m";
        let (cat, desc) = parse_about(input);
        assert_eq!(cat, "UMI EXTRACTION");
        assert_eq!(desc, "Extract UMIs from FASTQ");
    }

    #[test]
    fn test_parse_about_no_ansi() {
        let (cat, desc) = parse_about("[CONSENSUS] Call reads");
        assert_eq!(cat, "CONSENSUS");
        assert_eq!(desc, "Call reads");
    }

    #[test]
    fn test_collect_commands_finds_all() {
        let commands = collect_commands();
        let names: Vec<_> = commands.iter().map(|(n, ..)| n.as_str()).collect();
        assert!(names.contains(&"extract"));
        assert!(names.contains(&"group"));
        assert!(names.contains(&"simplex"));
        assert!(names.contains(&"duplex"));
        assert!(names.contains(&"filter"));
        assert!(names.contains(&"merge"));
        assert!(names.contains(&"simplex-metrics"));
        assert!(names.len() >= 17);
    }

    #[test]
    fn test_render_tool_page_contains_args() {
        let commands = collect_commands();
        let (name, category, description, long_about, cmd) =
            commands.iter().find(|(n, ..)| n == "extract").unwrap();
        let md = render_tool_page(name, category, description, long_about.as_deref(), cmd);
        assert!(md.contains("# extract"));
        assert!(md.contains("## Arguments"));
        assert!(md.contains("--input"));
    }
}
