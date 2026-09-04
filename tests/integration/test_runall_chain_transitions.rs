//! In-process `build_for`+run tests for each multi-stage transition `runall`
//! depends on. These validate the chain plumbing (type-erased hand-offs)
//! BEFORE the CLI is wired, per spec §9.
use std::fs;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use fgumi_lib::commands::common::{
    CompressionOptions, QueueMemoryOptions, SchedulerOptions, ThreadingOptions,
};
use fgumi_lib::commands::correct::CorrectOptions;
use fgumi_lib::commands::extract::ExtractRunallOptions;
use fgumi_lib::pipeline::chains::{
    ChainSpec, SinkSpec, SourceSpec, Stage, StageOptionsBag, build_for,
};
use tempfile::TempDir;

/// Write a gzip-compressed FASTQ from `(name, seq, qual)` records.
fn write_gzip_fastq(path: &Path, records: &[(&str, &str, &str)]) {
    let file = fs::File::create(path).expect("create gzip fastq");
    let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    for (name, seq, qual) in records {
        writeln!(encoder, "@{name}\n{seq}\n+\n{qual}").expect("write fastq record");
    }
    encoder.finish().expect("finish gzip fastq");
}

/// Build a small paired gzip FASTQ (`r1.fq.gz`, `r2.fq.gz`) with a 4 bp UMI on
/// R1 only (read structures `4M+T` / `+T`) — 2 families x 3 read pairs each.
/// The R1 UMI for each family is one of the correct-step's own whitelist
/// entries (`ACGT`, `TGCA`), so every extracted record is an exact match and
/// the correct step keeps it (a mismatched/ambiguous UMI would be dropped,
/// which would make a non-empty-output assertion meaningless). Returns
/// `(r1_path, r2_path)`.
fn write_tiny_umi_fastqs(dir: &Path) -> (PathBuf, PathBuf) {
    let r1 = dir.join("r1.fq.gz");
    let r2 = dir.join("r2.fq.gz");
    // (r1 umi (4bp, matches the correct-step whitelist), r1 template (12bp),
    // r2 template (12bp)).
    let families =
        [("ACGT", "ACGTACGTACGT", "GGTTAACCGGTT"), ("TGCA", "TGCATGCATGCA", "CCAATTGGCCAA")];
    let r1_qual = "I".repeat(4 + 12); // 4 (UMI) + 12 (template)
    let r2_qual = "I".repeat(12);
    let mut r1_owned: Vec<(String, String)> = Vec::new();
    let mut r2_owned: Vec<(String, String)> = Vec::new();
    for (fi, (umi, r1_tmpl, r2_tmpl)) in families.iter().enumerate() {
        for i in 0..3 {
            let name = format!("fam{fi}_{i}");
            r1_owned.push((name.clone(), format!("{umi}{r1_tmpl}"))); // 4 + 12 = 16
            r2_owned.push((name, (*r2_tmpl).to_string())); // 12
        }
    }
    let r1_slices: Vec<(&str, &str, &str)> =
        r1_owned.iter().map(|(n, s)| (n.as_str(), s.as_str(), r1_qual.as_str())).collect();
    let r2_slices: Vec<(&str, &str, &str)> =
        r2_owned.iter().map(|(n, s)| (n.as_str(), s.as_str(), r2_qual.as_str())).collect();
    write_gzip_fastq(&r1, &r1_slices);
    write_gzip_fastq(&r2, &r2_slices);
    (r1, r2)
}

/// Extract→Correct is the one `runall` transition that is unwired at the
/// chain-builder level: `add_correct` used to unconditionally prepend a
/// `GroupByQueryname` step (input `DecodedRecordBatch`) onto whatever tail
/// preceded it. Extract's output tail is `BamTemplateBatch`, so feeding an
/// extract-fed chain into correct was a type-erased runtime panic, not a
/// build-time error. This builds `[Stage::Extract, Stage::Correct]` directly
/// (bypassing the not-yet-existing `RunAll` command) and asserts it runs to
/// completion and produces a non-empty corrected unmapped BAM.
#[test]
fn extract_to_correct_chain_builds_and_runs() {
    let tmp = TempDir::new().unwrap();
    let (r1, r2) = write_tiny_umi_fastqs(tmp.path());
    let out = tmp.path().join("corrected.bam");

    let extract = ExtractRunallOptions {
        inputs: vec![r1.clone(), r2.clone()],
        read_structures: vec!["4M+T".parse().unwrap(), "+T".parse().unwrap()],
        sample: "s1".into(),
        library: "lib1".into(),
        ..ExtractRunallOptions::default()
    };
    let correct = CorrectOptions {
        umis: vec!["ACGT".into(), "TGCA".into()], // 4 bp whitelist, not 2 bp
        min_distance_diff: 1,
        ..CorrectOptions::default()
    };

    let bag = StageOptionsBag {
        extract: Some(extract.to_extract_options()),
        correct: Some(correct),
        ..Default::default()
    };

    let spec = ChainSpec {
        stages: vec![Stage::Extract, Stage::Correct],
        source: SourceSpec::fastqs(
            vec![r1, r2],
            vec!["4M+T".parse().unwrap(), "+T".parse().unwrap()],
        )
        .unwrap(),
        sink: SinkSpec::Bam(out.clone()),
        stage_opts: bag,
        threading: ThreadingOptions { threads: None },
        compression: CompressionOptions::default(),
        scheduler: SchedulerOptions::default(),
        queue_memory: QueueMemoryOptions::default(),
        async_reader: false,
        read_streams: fgumi_bam_io::ReadStreams::Fixed(1),
        verify_crc: false,
        command_line: "test".into(),
    };

    build_for(spec).expect("build").run().expect("run without panic");
    assert!(fs::metadata(&out).unwrap().len() > 0, "corrected BAM is non-empty");
}
