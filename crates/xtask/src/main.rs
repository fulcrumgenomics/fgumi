#![deny(unsafe_code)]
use anyhow::{Context, Result};
use clap::Parser;
use std::path::PathBuf;

mod check_tag_literals;
mod generate_docs;
mod generate_metrics;
mod generate_summary;
mod generate_tools;

#[derive(Parser, Debug)]
#[command(name = "xtask", about = "Build tasks for fgumi")]
enum Xtask {
    /// Generate documentation markdown and serve with mdbook
    Docs,
    /// Generate documentation markdown and build with mdbook
    DocsBuild,
    /// Generate documentation markdown only (no mdbook)
    GenerateDocs,
    /// Scan workspace for bare 2-byte SAM-tag byte literals that should use `SamTag` constants
    CheckTagLiterals,
}

fn main() -> Result<()> {
    let args = Xtask::parse();
    match args {
        Xtask::GenerateDocs => generate_docs::generate()?,
        Xtask::Docs => {
            generate_docs::generate()?;
            run_mdbook(&["serve", "docs"])?;
        }
        Xtask::DocsBuild => {
            generate_docs::generate()?;
            run_mdbook(&["build", "docs"])?;
        }
        Xtask::CheckTagLiterals => check_tag_literals_cmd()?,
    }
    Ok(())
}

fn check_tag_literals_cmd() -> Result<()> {
    let root = workspace_root()?;
    let findings = check_tag_literals::scan_workspace(&root)?;
    if findings.is_empty() {
        println!("check-tag-literals: OK (no bare SAM-tag byte literals found)");
    } else {
        for f in &findings {
            eprintln!(
                "{}:{}  {} (payload=b\"{}{}\")",
                f.path.display(),
                f.line,
                f.snippet,
                f.payload[0] as char,
                f.payload[1] as char,
            );
        }
        anyhow::bail!(
            "{} bare SAM-tag byte literal(s) found outside the allowlist; \
             replace with SamTag::XX constants or add to PAYLOAD_ALLOWLIST with justification",
            findings.len()
        );
    }
    Ok(())
}

fn workspace_root() -> Result<PathBuf> {
    // CARGO_MANIFEST_DIR points to crates/xtask; go up two levels to the workspace root.
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .with_context(|| "failed to resolve workspace root")
}

fn run_mdbook(args: &[&str]) -> Result<()> {
    let status = std::process::Command::new("mdbook").args(args).status()?;
    if !status.success() {
        anyhow::bail!("mdbook failed with exit code: {status}");
    }
    Ok(())
}
