//! R guards and assertion shared by the duplex- and simplex-metrics PDF tests.

use std::fs;
use std::path::PathBuf;

use clap::Parser;
use fgumi_lib::commands::command::Command;
use fgumi_lib::commands::shared_metrics::{is_r_available, is_r_pdf_reproducible};
use tempfile::TempDir;

/// Runs `subcommand` twice on the BAM that `write_input_bam` writes and asserts both runs'
/// `<prefix>.<pdf_suffix>` PDFs are byte-identical; skips without R 4.5.0+, ggplot2 and scales.
pub fn assert_pdf_byte_identical_across_runs<C: Parser + Command>(
    subcommand: &str,
    pdf_suffix: &str,
    write_input_bam: impl FnOnce(&PathBuf),
) {
    if !is_r_available() {
        eprintln!("Skipping: R with ggplot2 and scales not available");
        return;
    }
    if !is_r_pdf_reproducible() {
        eprintln!("Skipping: R lacks pdf(timestamp, producer), added in R 4.5.0");
        return;
    }
    let temp_dir = TempDir::new().expect("Failed to create temp dir");
    let input_bam = temp_dir.path().join("input.bam");
    write_input_bam(&input_bam);

    let render = |run: &str| {
        let output_prefix = temp_dir.path().join(run);
        let cmd = C::try_parse_from([
            subcommand,
            "--input",
            input_bam.to_str().unwrap(),
            "--output",
            output_prefix.to_str().unwrap(),
        ])
        .unwrap_or_else(|e| panic!("failed to parse {subcommand} args: {e}"));
        cmd.execute(&format!("fgumi {subcommand}"))
            .unwrap_or_else(|e| panic!("{subcommand} command failed: {e}"));
        fs::read(format!("{}.{pdf_suffix}", output_prefix.display())).expect("PDF not created")
    };
    assert!(render("first") == render("second"), "PDF bytes differ between identical runs");
}
