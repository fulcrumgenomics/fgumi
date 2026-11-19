//! Index ownership and failed finalization must not remove unrelated files or
//! leave an index advertising incomplete output.

use fgumi_htslib::{CramWriteOptions, CramWriter};

const HEADER: &str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n";

#[cfg(unix)]
#[rstest::rstest]
#[case::dangling_index_symlink(false)]
#[case::index_is_a_directory(true)]
fn failed_reference_setup_preserves_existing_output_and_index(#[case] index_directory: bool) {
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
    std::fs::write(&output, b"existing output").unwrap();
    std::fs::write(&index, b"existing index").unwrap();
    let opts =
        CramWriteOptions { reference: Some(&reference), crai: Some(&index), ..Default::default() };

    let Err(error) = CramWriter::create(&output, HEADER, &opts) else {
        panic!("an unusable reference index must fail setup");
    };
    assert!(error.to_string().contains("the CRAM reference index"), "{error}");
    // The reference is checked before the destination is created (truncated).
    assert_eq!(std::fs::read(output).unwrap(), b"existing output");
    assert_eq!(std::fs::read(index).unwrap(), b"existing index");
}

/// `data` as one BGZF block followed by the BGZF EOF block.
fn bgzf(data: &[u8]) -> Vec<u8> {
    use std::io::Write;
    let mut deflate =
        flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::default());
    deflate.write_all(data).unwrap();
    let cdata = deflate.finish().unwrap();
    let mut crc = flate2::Crc::new();
    crc.update(data);
    let bsize = u16::try_from(18 + cdata.len() + 8 - 1).unwrap();
    let mut out = vec![0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0];
    out.extend_from_slice(&bsize.to_le_bytes());
    out.extend_from_slice(&cdata);
    out.extend_from_slice(&crc.sum().to_le_bytes());
    out.extend_from_slice(&u32::try_from(data.len()).unwrap().to_le_bytes());
    // The 28-byte BGZF EOF block.
    out.extend_from_slice(&[
        0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0, 0x1b, 0, 3, 0, 0, 0, 0, 0,
        0, 0, 0, 0,
    ]);
    out
}

/// A BGZF-compressed reference whose `.fai` is fine but whose `.gzi` is
/// missing or unreadable fails before the destination is created: htslib
/// would load the `.gzi` only after truncating it. The same reference with
/// both sidecars intact sets up, so the `.gzi` is the only fault.
#[rstest::rstest]
#[case::missing_gzi(None)]
#[case::empty_gzi(Some(&b""[..]))]
#[case::truncated_gzi(Some(&[1u8, 0, 0][..]))]
fn bgzf_reference_with_bad_gzi_preserves_existing_output(#[case] gzi_bytes: Option<&[u8]>) {
    let dir = tempfile::tempdir().unwrap();
    let reference = dir.path().join("reference.fa.gz");
    std::fs::write(&reference, bgzf(b">chr1\nACGTACGTAC\n")).unwrap();
    let fai = dir.path().join("reference.fa.gz.fai");
    let gzi = dir.path().join("reference.fa.gz.gzi");
    let opts = CramWriteOptions { reference: Some(&reference), ..Default::default() };
    let header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10\n";

    // Setup builds both sidecars and succeeds.
    let writer = CramWriter::create(&dir.path().join("first.cram"), header, &opts).unwrap();
    let (jobs, finisher) = writer.seal().unwrap();
    assert!(jobs.is_empty());
    finisher.finish().unwrap();
    assert!(fai.exists() && gzi.exists(), "setup indexed the BGZF reference");

    match gzi_bytes {
        None => std::fs::remove_file(&gzi).unwrap(),
        Some(bytes) => std::fs::write(&gzi, bytes).unwrap(),
    }
    let output = dir.path().join("output.cram");
    std::fs::write(&output, b"existing output").unwrap();

    let Err(error) = CramWriter::create(&output, header, &opts) else {
        panic!("a BGZF reference with a bad .gzi must fail setup");
    };
    assert!(error.to_string().contains("reference.fa.gz.gzi"), "{error}");
    assert_eq!(std::fs::read(output).unwrap(), b"existing output");
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
