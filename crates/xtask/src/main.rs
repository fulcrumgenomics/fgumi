#![deny(unsafe_code)]
use anyhow::Result;
use clap::Parser;

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
    }
    Ok(())
}

fn run_mdbook(args: &[&str]) -> Result<()> {
    let status = std::process::Command::new("mdbook").args(args).status()?;
    if !status.success() {
        anyhow::bail!("mdbook failed with exit code: {status}");
    }
    Ok(())
}
