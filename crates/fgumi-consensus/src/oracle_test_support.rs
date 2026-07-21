//! Shared test-only helpers for the offline-pinned fgbio-oracle tests in
//! [`codec_caller`](crate::codec_caller) and
//! [`duplex_caller`](crate::duplex_caller).
//!
//! Each oracle test drives fgumi's own caller on a fixture and asserts the
//! emitted depth/error tags equal values CAPTURED FROM A REAL fgbio RUN on the
//! exact same input BAM, rather than re-derived from the implementation's own
//! formula. This closes the CodeRabbit finding that the hand-authored
//! saturation tests "repeat the implementation's expected formulas, so a shared
//! misunderstanding of fgbio's cap-before-sum behavior would still pass".
//!
//! Equivalence: one [`SamBuilder`] produces the fixture and writes ONE BAM. That
//! exact file is both (1) what the offline fgbio run consumes and (2) what the
//! fgumi test reads back into `RawRecord`s. See the `#[ignore]`d `regen_*` tests
//! in each caller module to re-emit the fixture BAM for a fresh fgbio capture.

use std::path::Path;

use fgumi_raw_bam::{ParsedBamRecord, RawRecord, SamTag};
use fgumi_sam::builder::SamBuilder;

/// Read a BAM at `path` back into a `Vec<RawRecord>` in file order — the same
/// records the offline fgbio run sees (byte-identical input equivalence).
pub(crate) fn raw_records_from_bam(path: &Path) -> Vec<RawRecord> {
    let (mut reader, _header) =
        fgumi_bam_io::create_raw_bam_reader(path, 1).expect("open fixture BAM");
    let mut records = Vec::new();
    let mut rec = RawRecord::new();
    loop {
        let n = reader.read_record(&mut rec).expect("read raw record");
        if n == 0 {
            break;
        }
        records.push(rec.clone());
    }
    records
}

/// Write `builder` to a temp BAM and read the records straight back. The
/// permanent oracle tests use this so the fgumi caller consumes exactly the
/// bytes the offline fgbio run would (a temp file, never overwriting a
/// maintainer's `FGUMI_ORACLE_BAM_OUT`).
pub(crate) fn records_from_builder(builder: &SamBuilder) -> Vec<RawRecord> {
    let tmp = tempfile::NamedTempFile::new().expect("temp fixture BAM");
    builder.write_bam(tmp.path()).expect("write fixture BAM");
    raw_records_from_bam(tmp.path())
}

/// Re-emit a fixture BAM for a fresh fgbio capture: writes to
/// `FGUMI_ORACLE_BAM_OUT` if set, otherwise a temp file whose path is logged.
/// Used only by the `#[ignore]`d `regen_*` tests.
pub(crate) fn regen_write(builder: &SamBuilder) {
    if let Some(out) = std::env::var_os("FGUMI_ORACLE_BAM_OUT") {
        let path = Path::new(&out);
        builder.write_bam(path).expect("write fixture BAM to FGUMI_ORACLE_BAM_OUT");
        eprintln!("wrote fixture BAM to {}", path.display());
    } else {
        let tmp = tempfile::NamedTempFile::new().expect("temp fixture BAM");
        let (_file, path) = tmp.keep().expect("persist temp fixture BAM");
        builder.write_bam(&path).expect("write fixture BAM");
        eprintln!("wrote fixture BAM to {}", path.display());
    }
}

/// The full set of fgbio-captured (pinned) consensus tag values for one
/// consensus record. Bundling them keeps each `#[rstest]` case table readable
/// (expected values live in the table) and each case's arg count small.
pub(crate) struct ExpectedOracleTags {
    pub(crate) cd: i32,
    pub(crate) cm: i32,
    pub(crate) ce: f32,
    pub(crate) ad: i32,
    pub(crate) am: i32,
    pub(crate) ae: f32,
    pub(crate) bd: i32,
    pub(crate) bm: i32,
    pub(crate) be: f32,
    pub(crate) ad_bases: Vec<i16>,
    pub(crate) bd_bases: Vec<i16>,
    pub(crate) ae_bases: Vec<i16>,
    pub(crate) be_bases: Vec<i16>,
}

impl ExpectedOracleTags {
    /// Assert every emitted depth/error tag on `rec` equals the pinned value.
    pub(crate) fn assert_matches(&self, rec: &ParsedBamRecord) {
        assert_eq!(rec.get_int_tag(SamTag::CD), Some(self.cd), "cD (fgbio-pinned)");
        assert_eq!(rec.get_int_tag(SamTag::CM), Some(self.cm), "cM (fgbio-pinned)");
        assert_eq!(rec.get_int_tag(SamTag::AD), Some(self.ad), "aD (fgbio-pinned)");
        assert_eq!(rec.get_int_tag(SamTag::AM), Some(self.am), "aM (fgbio-pinned)");
        assert_eq!(rec.get_int_tag(SamTag::BD), Some(self.bd), "bD (fgbio-pinned)");
        assert_eq!(rec.get_int_tag(SamTag::BM), Some(self.bm), "bM (fgbio-pinned)");

        let ce = rec.get_float_tag(SamTag::CE).expect("cE present");
        assert!((ce - self.ce).abs() < 1e-6, "cE (fgbio-pinned): expected {}, got {ce}", self.ce);
        let ae = rec.get_float_tag(SamTag::AE).expect("aE present");
        assert!((ae - self.ae).abs() < 1e-6, "aE (fgbio-pinned): expected {}, got {ae}", self.ae);
        let be = rec.get_float_tag(SamTag::BE).expect("bE present");
        assert!((be - self.be).abs() < 1e-6, "bE (fgbio-pinned): expected {}, got {be}", self.be);

        assert_eq!(
            rec.get_i16_array_tag(SamTag::AD_BASES),
            Some(self.ad_bases.clone()),
            "ad (fgbio-pinned)"
        );
        assert_eq!(
            rec.get_i16_array_tag(SamTag::BD_BASES),
            Some(self.bd_bases.clone()),
            "bd (fgbio-pinned)"
        );
        assert_eq!(
            rec.get_i16_array_tag(SamTag::AE_BASES),
            Some(self.ae_bases.clone()),
            "ae (fgbio-pinned)"
        );
        assert_eq!(
            rec.get_i16_array_tag(SamTag::BE_BASES),
            Some(self.be_bases.clone()),
            "be (fgbio-pinned)"
        );
    }
}
