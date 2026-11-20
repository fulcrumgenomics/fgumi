//! Index ownership and failed finalization must not remove unrelated files or
//! leave an index advertising incomplete output.

use fgumi_htslib::{CramWriteOptions, CramWriter};

const HEADER: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n";

#[cfg(unix)]
#[rstest::rstest]
#[case::dangling_index_symlink(false)]
#[case::index_is_a_directory(true)]
fn failed_reference_setup_preserves_an_existing_index(#[case] index_directory: bool) {
    let dir = tempfile::tempdir().unwrap();
    let reference = dir.path().join("reference.fa");
    std::fs::write(&reference, b">chr1\nACGT\n").unwrap();
    // Both paths prevent reference-index setup even when tests run as root.
    let fai = dir.path().join("reference.fa.fai");
    if index_directory {
        std::fs::create_dir(&fai).unwrap();
    } else {
        std::os::unix::fs::symlink(dir.path().join("missing/reference.fa.fai"), &fai).unwrap();
    }
    let output = dir.path().join("output.cram");
    let index = dir.path().join("output.cram.crai");
    std::fs::write(&index, b"existing index").unwrap();
    let opts =
        CramWriteOptions { reference: Some(&reference), crai: Some(&index), ..Default::default() };

    let Err(error) = CramWriter::create(&output, HEADER, &opts) else {
        panic!("an unusable reference index must fail setup");
    };
    assert!(error.to_string().contains("setting the CRAM reference"), "{error}");
    assert_eq!(std::fs::read(index).unwrap(), b"existing index");
}

#[cfg(unix)]
#[test]
fn failed_final_flush_removes_the_owned_index() {
    const CHILD: &str = "FGUMI_TEST_FAILED_FINAL_FLUSH";
    if std::env::var_os(CHILD).is_some() {
        let dir = tempfile::tempdir().unwrap();
        let output = dir.path().join("output.cram");
        let index = dir.path().join("output.cram.crai");
        let opts = CramWriteOptions { no_ref: true, crai: Some(&index), ..Default::default() };
        let writer = CramWriter::create(&output, HEADER, &opts).unwrap();
        assert!(index.exists(), "the writer created its index");
        let (jobs, finisher) = writer.seal().unwrap();
        assert!(jobs.is_empty());
        rustix::process::setrlimit(
            rustix::process::Resource::Fsize,
            rustix::process::Rlimit { current: Some(0), maximum: Some(0) },
        )
        .unwrap();
        let error = finisher.finish().unwrap_err();
        assert!(error.to_string().contains("closing the CRAM output"), "{error}");
        assert!(!index.exists(), "a failed close must remove the writer-owned index");
        return;
    }

    // htslib buffers this header-only output until close. A zero file-size
    // limit forces the real final flush to fail, independently of disk space.
    let output = std::process::Command::new("sh")
        .args(["-c", "trap '' XFSZ; exec \"$1\" --exact failed_final_flush_removes_the_owned_index --nocapture", "sh"])
        .arg(std::env::current_exe().unwrap())
        .env(CHILD, "1")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "final-flush child failed: {}{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}
