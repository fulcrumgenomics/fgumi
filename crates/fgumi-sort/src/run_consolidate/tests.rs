use super::*;
use crate::external::GenericKeyedChunkReader;
use rstest::rstest;
use tempfile::TempDir;

/// A template-coordinate key (40-byte prefix, not embedded in the body) whose
/// order is fully determined by `primary`; `name_hash` is fixed so equal
/// `primary` values produce equal keys and the tie-break is observable.
#[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
fn template_key(primary: u64) -> TemplateKey {
    TemplateKey::new(
        0,
        primary as i32,
        false,
        i32::MAX,
        i32::MAX,
        false,
        0,
        0,
        (0, false),
        0,
        false,
    )
}

/// A minimal BAM record body carrying `tid`/`pos` where the coordinate key is
/// extracted from, plus a trailing `tag` that makes each record distinguishable
/// even when its key ties with another's.
fn coordinate_body(tid: i32, pos: i32, tag: u32) -> Vec<u8> {
    let mut body = vec![0u8; 36];
    body[0..4].copy_from_slice(&tid.to_le_bytes());
    body[4..8].copy_from_slice(&pos.to_le_bytes());
    body[32..36].copy_from_slice(&tag.to_le_bytes());
    body
}

/// Write `records` as one spill run with the spill kernels, using an odd raw
/// block size so records routinely straddle compressed-unit boundaries.
fn write_run<K: RawSortKey>(
    path: &Path,
    codec: SpillCodec,
    records: &[(K, Vec<u8>)],
    block_size: usize,
) {
    let mut writer = SpillRunWriter::create(path, codec, 1, block_size).expect("create run");
    for (key, body) in records {
        writer.write_record(key, body).expect("write record");
    }
    writer.finish().expect("finish run");
}

/// Drive a merger to completion one record at a time, so every record is
/// produced by a separate, resumed `step` call.
fn merge_one_at_a_time(merger: &mut dyn RunMergerDyn) -> u64 {
    loop {
        match merger.step(1, u64::MAX).expect("merge step") {
            RunMergeProgress::Working => {}
            RunMergeProgress::Done { records } => return records,
        }
    }
}

/// Read a run back with an independent reader: the retired engine's keyed
/// chunk reader, which the arena spill format is byte-compatible with.
fn read_back_template(path: &Path) -> Vec<(TemplateKey, Vec<u8>)> {
    let mut reader = GenericKeyedChunkReader::<TemplateKey>::open(path, None).expect("open");
    let mut body = Vec::new();
    let mut out = Vec::new();
    while let Some(key) = reader.next_record(&mut body).expect("read record") {
        out.push((key, body.clone()));
    }
    out
}

/// The records a correct merge must produce: the concatenation of the inputs in
/// run order, stably sorted by key — exactly what a single k-way merge over the
/// unconsolidated runs emits.
fn stable_merge_expectation<K: RawSortKey>(runs: &[Vec<(K, Vec<u8>)>]) -> Vec<(K, Vec<u8>)> {
    let mut all: Vec<(K, Vec<u8>)> = runs.iter().flatten().cloned().collect();
    all.sort_by(|a, b| a.0.cmp(&b.0));
    all
}

#[rstest]
#[case::zstd(SpillCodec::Zstd)]
#[case::bgzf(SpillCodec::Bgzf)]
fn merge_is_the_stable_merge_of_its_inputs_in_run_order(#[case] codec: SpillCodec) {
    let dir = TempDir::new().unwrap();
    // Three runs whose keys overlap and repeat across runs (keys 0..40 step 3
    // shifted per run, so many keys appear in more than one run). Each body
    // names its run and index, so a tie resolved in the wrong order shows up
    // as a body mismatch even though the keys agree.
    let runs: Vec<Vec<(TemplateKey, Vec<u8>)>> = (0..3u8)
        .map(|run| {
            (0..40u64)
                .map(|i| {
                    (
                        template_key(i / 2 + u64::from(run)),
                        vec![run, u8::try_from(i).unwrap(), 0xAB],
                    )
                })
                .collect()
        })
        .collect();
    let mut inputs = Vec::new();
    for (i, run) in runs.iter().enumerate() {
        let path = dir.path().join(format!("chunk_{i}.keyed"));
        write_run(&path, codec, run, 23);
        inputs.push(path);
    }
    let output = dir.path().join("merged.keyed");
    let spec =
        RunMergeSpec { inputs: &inputs, output: &output, codec, compression: 1, block_size: 29 };
    let mut merger = new_run_merger(SpillKeyKind::TemplateK40, &spec).expect("open merger");

    let written = merge_one_at_a_time(merger.as_mut());

    let expected = stable_merge_expectation(&runs);
    assert_eq!(written, expected.len() as u64);
    assert_eq!(read_back_template(&output), expected);
}

#[test]
fn merge_of_embedded_keys_re_extracts_them_from_the_body() {
    let dir = TempDir::new().unwrap();
    // Coordinate keys live inside the record body, so the file carries no key
    // prefix; the merger must recover them from the body to order the output.
    let runs: Vec<Vec<(RawCoordinateKey, Vec<u8>)>> = (0..4u32)
        .map(|run| {
            (0..25u32)
                .map(|i| {
                    let body =
                        coordinate_body(0, (i * 2 + run % 2).try_into().unwrap(), run * 100 + i);
                    (RawCoordinateKey::extract_from_record(&body), body)
                })
                .collect()
        })
        .collect();
    let mut inputs = Vec::new();
    for (i, run) in runs.iter().enumerate() {
        let path = dir.path().join(format!("chunk_{i}.keyed"));
        write_run(&path, SpillCodec::Zstd, run, 50);
        inputs.push(path);
    }
    let output = dir.path().join("merged.keyed");
    let spec = RunMergeSpec {
        inputs: &inputs,
        output: &output,
        codec: SpillCodec::Zstd,
        compression: 1,
        block_size: 64 * 1024,
    };
    let mut merger = new_run_merger(SpillKeyKind::Coordinate, &spec).expect("open merger");
    assert_eq!(
        merger.step(usize::MAX, u64::MAX).expect("merge"),
        RunMergeProgress::Done { records: 100 }
    );

    let mut reader = SpillRunReader::open(&output).expect("open output");
    let mut dec = SpillBlockDecompressor::new();
    let (mut key_buf, mut body) = (Vec::new(), Vec::new());
    let mut got = Vec::new();
    while let Some(key) =
        reader.next_record::<RawCoordinateKey>(&mut dec, 0, &mut key_buf, &mut body).unwrap()
    {
        got.push((key, body.clone()));
    }
    assert_eq!(got, stable_merge_expectation(&runs));
}

#[test]
fn empty_inputs_are_skipped_and_an_all_empty_merge_writes_an_empty_run() {
    let dir = TempDir::new().unwrap();
    let empty_a = dir.path().join("a.keyed");
    let full = dir.path().join("b.keyed");
    let empty_c = dir.path().join("c.keyed");
    write_run::<TemplateKey>(&empty_a, SpillCodec::Bgzf, &[], 64);
    let records = vec![(template_key(1), vec![1]), (template_key(2), vec![2])];
    write_run(&full, SpillCodec::Bgzf, &records, 64);
    write_run::<TemplateKey>(&empty_c, SpillCodec::Bgzf, &[], 64);

    let output = dir.path().join("merged.keyed");
    let inputs = vec![empty_a.clone(), full, empty_c.clone()];
    let spec = RunMergeSpec {
        inputs: &inputs,
        output: &output,
        codec: SpillCodec::Bgzf,
        compression: 1,
        block_size: 64,
    };
    let mut merger = new_run_merger(SpillKeyKind::TemplateK40, &spec).unwrap();
    assert_eq!(merge_one_at_a_time(merger.as_mut()), 2);
    assert_eq!(read_back_template(&output), records);

    let all_empty = dir.path().join("empty.keyed");
    let inputs = vec![empty_a, empty_c];
    let spec = RunMergeSpec {
        inputs: &inputs,
        output: &all_empty,
        codec: SpillCodec::Bgzf,
        compression: 1,
        block_size: 64,
    };
    let mut merger = new_run_merger(SpillKeyKind::TemplateK40, &spec).unwrap();
    assert_eq!(merger.step(10, u64::MAX).unwrap(), RunMergeProgress::Done { records: 0 });
    assert_eq!(
        merger.step(10, u64::MAX).unwrap(),
        RunMergeProgress::Done { records: 0 },
        "Done is sticky"
    );
    assert!(read_back_template(&all_empty).is_empty());
}

#[test]
fn a_run_truncated_mid_record_is_an_unexpected_eof() {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("chunk.keyed");
    // Write the framing by hand into one uncompressed-stream zstd block that
    // stops inside the second record's body.
    let mut raw = Vec::new();
    frame_keyed_record_into(&mut raw, &template_key(1), &[1, 2, 3]).unwrap();
    frame_keyed_record_into(&mut raw, &template_key(2), &[4, 5, 6, 7]).unwrap();
    raw.truncate(raw.len() - 2);
    let mut compressor = SpillBlockCompressor::new(SpillCodec::Zstd, 1).unwrap();
    let mut file = File::create(&path).unwrap();
    file.write_all(spill_magic(SpillCodec::Zstd)).unwrap();
    file.write_all(&compressor.compress_block(&raw).unwrap()).unwrap();
    drop(file);

    let output = dir.path().join("merged.keyed");
    let inputs = vec![path];
    let spec = RunMergeSpec {
        inputs: &inputs,
        output: &output,
        codec: SpillCodec::Zstd,
        compression: 1,
        block_size: 64,
    };
    let mut merger =
        new_run_merger(SpillKeyKind::TemplateK40, &spec).expect("first record is intact");
    let err = merger.step(usize::MAX, u64::MAX).expect_err("second record is truncated");
    assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
}

#[test]
fn the_output_path_must_not_already_exist() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("chunk.keyed");
    write_run(&input, SpillCodec::Zstd, &[(template_key(1), vec![1])], 64);
    let inputs = vec![input.clone()];
    // Pointing the output at a live input must fail rather than truncate it.
    let spec = RunMergeSpec {
        inputs: &inputs,
        output: &input,
        codec: SpillCodec::Zstd,
        compression: 1,
        block_size: 64,
    };
    let Err(err) = new_run_merger(SpillKeyKind::TemplateK40, &spec) else {
        panic!("creating the output over an existing file must fail");
    };
    assert_eq!(err.kind(), io::ErrorKind::AlreadyExists);
}

/// The byte budget bounds a call independently of the record budget, so a call
/// over long reads stays short even when the record budget is huge.
#[test]
fn the_byte_budget_bounds_a_step_independently_of_the_record_budget() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("chunk.keyed");
    let records: Vec<(TemplateKey, Vec<u8>)> =
        (0..10u64).map(|i| (template_key(i), vec![7u8; 100])).collect();
    write_run(&input, SpillCodec::Zstd, &records, 64);
    let output = dir.path().join("merged.keyed");
    let inputs = vec![input];
    let spec = RunMergeSpec {
        inputs: &inputs,
        output: &output,
        codec: SpillCodec::Zstd,
        compression: 1,
        block_size: 64,
    };
    let mut merger = new_run_merger(SpillKeyKind::TemplateK40, &spec).unwrap();
    // 250 bytes of 100-byte bodies: the third record crosses the budget and
    // ends the call, whatever the record budget allows.
    assert_eq!(merger.step(usize::MAX, 250).unwrap(), RunMergeProgress::Working);
    assert_eq!(merger.records_written(), 3);
    // A budget smaller than one record still makes progress.
    assert_eq!(merger.step(usize::MAX, 1).unwrap(), RunMergeProgress::Working);
    assert_eq!(merger.records_written(), 4);
    while merger.step(usize::MAX, 250).unwrap() == RunMergeProgress::Working {}
    assert_eq!(read_back_template(&output), records);
}
