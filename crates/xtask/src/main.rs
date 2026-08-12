#![deny(unsafe_code)]
use anyhow::{Context, Result};
use clap::Parser;
use std::path::{Path, PathBuf};

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
    /// Validate that the crates.io publish list names every publishable workspace crate, and
    /// only those, in dependency order
    CheckPublishOrder,
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
        Xtask::CheckPublishOrder => check_publish_order_cmd()?,
    }
    Ok(())
}

/// Run `scripts/publish-crates.sh --check`.
///
/// This delegates to the script rather than reimplementing the check in Rust so
/// there is exactly one definition of the publish list and its ordering rules.
/// The release workflow and the `publish-order` CI job already invoke the script
/// directly; a second, independent implementation here could drift from it, and a
/// local gate that disagrees with CI is worse than no local gate.
fn check_publish_order_cmd() -> Result<()> {
    check_publish_order_in(&workspace_root()?)
}

/// Run `<root>/scripts/publish-crates.sh --check`, failing if the script does.
///
/// `root` is a parameter rather than resolved inside, so the delegation itself --
/// which script runs, that it is passed `--check`, and that a non-zero exit
/// propagates as an error rather than a silent pass -- is testable against a stub.
/// Testing it against the real workspace could not cover the failure path: the
/// publish list there is, by construction, always correct.
fn check_publish_order_in(root: &Path) -> Result<()> {
    let script = root.join("scripts/publish-crates.sh");
    // Invoked through `bash` rather than executed directly so the check does not
    // depend on the file's executable bit surviving a checkout or export.
    let status = std::process::Command::new("bash")
        .arg(&script)
        .arg("--check")
        .current_dir(root)
        .status()
        .with_context(|| format!("failed to run {}", script.display()))?;
    if !status.success() {
        anyhow::bail!("publish-crates.sh --check failed with exit code: {status}");
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Build a fake workspace root whose `publish-crates.sh` exits with `exit_code`
    /// and records the arguments it was handed in `args.txt`.
    fn fake_root(exit_code: i32) -> tempfile::TempDir {
        let dir = tempfile::tempdir().expect("temp dir");
        let scripts = dir.path().join("scripts");
        std::fs::create_dir_all(&scripts).expect("scripts dir");
        let mut f = std::fs::File::create(scripts.join("publish-crates.sh")).expect("stub script");
        // Deliberately left non-executable: the production path runs the script
        // through `bash`, so it must not depend on the executable bit.
        writeln!(f, "#!/usr/bin/env bash").expect("write stub");
        writeln!(f, "printf '%s\\n' \"$@\" > \"$(dirname \"$0\")/../args.txt\"")
            .expect("write stub");
        writeln!(f, "exit {exit_code}").expect("write stub");
        dir
    }

    /// A passing script is a passing check, and it is invoked with `--check`.
    ///
    /// The argument assertion is the load-bearing half: `publish-crates.sh` with no
    /// arguments prints usage and exits non-zero, and with `--list` it prints the
    /// crates and exits zero — so a check that dropped the flag would still pass
    /// here while validating nothing.
    #[test]
    fn publish_order_check_passes_and_forwards_the_check_flag() {
        let root = fake_root(0);
        check_publish_order_in(root.path()).expect("a script that exits 0 must pass the check");
        let args = std::fs::read_to_string(root.path().join("args.txt")).expect("args recorded");
        assert_eq!(args.trim(), "--check", "the script must be invoked with --check");
    }

    /// A failing script fails the check, rather than being swallowed.
    ///
    /// This is the whole point of the gate: `publish-crates.sh --check` exits 1 when
    /// the publish list is missing a crate, carries a stale entry, or is misordered.
    #[test]
    fn publish_order_check_propagates_script_failure() {
        let root = fake_root(1);
        let err = check_publish_order_in(root.path())
            .expect_err("a script that exits non-zero must fail the check");
        assert!(
            err.to_string().contains("publish-crates.sh --check failed"),
            "the error must name the script that failed, got: {err}"
        );
    }

    /// The real workspace passes its own gate, through the real path resolution.
    ///
    /// The stub tests above pin the delegation's behaviour but hand it a root, so they
    /// never exercise `workspace_root()` or the `scripts/publish-crates.sh` location
    /// relative to it. A wrong join there would leave every stub test green while the
    /// gate found no script in the repo it is meant to check. This is also a standing
    /// assertion that `CRATES` is currently correct, which is what the gate exists to
    /// say.
    #[test]
    fn publish_order_check_passes_on_this_workspace() {
        check_publish_order_cmd()
            .expect("the workspace's own publish list must satisfy the check it ships");
    }

    /// A missing script is an error, not a silent pass.
    ///
    /// `bash` exits 127 for a script it cannot open, so this lands on the same
    /// non-zero path as a real failure — pinned so a future rewrite that resolves
    /// the script differently cannot turn "not found" into "nothing to check".
    #[test]
    fn publish_order_check_fails_when_the_script_is_missing() {
        let root = tempfile::tempdir().expect("temp dir");
        assert!(
            check_publish_order_in(root.path()).is_err(),
            "a missing publish-crates.sh must fail the check"
        );
    }
}
