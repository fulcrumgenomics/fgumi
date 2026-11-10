//! Range tasks against the serial merge oracle.

use super::*;
use crate::inline::TemplateKey32;
use crate::keys::RawCoordinateKey;
use crate::range_merge::plan::{MAX_RANGES, RangeSizing};
use crate::range_merge::test_support::{frame_bodies, serial_reference, write_indexed_run};
use crate::spill_block::{SpillBlockCompressor, SpillKeyKind, spill_magic};
use crate::{CbKey32, SpillCodec, TemplateKey, TemplateKey24, TertKey32};
use proptest::prelude::*;
use proptest::test_runner::{Config as ProptestConfig, TestRunner};
use rstest::rstest;
use std::fmt::Debug;
use std::os::unix::fs::FileExt;
use std::path::Path;

/// Proptest cases per rstest case: 25 rstest cases x 64 keeps the property
/// well under a minute in a debug build.
const PROPTEST_CASES: u32 = 64;

/// A splitmix64 stream: the generated data is a pure function of the seed.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `0..n` (0 when `n == 0`).
    fn below(&mut self, n: u64) -> u64 {
        if n == 0 { 0 } else { self.next() % n }
    }

    fn between(&mut self, lo: usize, hi: usize) -> usize {
        lo + usize::try_from(self.below((hi - lo + 1) as u64)).unwrap()
    }
}

/// A BAM-shaped body of `len >= 32` bytes: `tid`, `pos` and `flag` where the
/// coordinate key reads them, and `id` in the mate fields so every body is
/// distinct.
fn bam_body(tid: i32, pos: i32, flag: u16, id: u64, len: usize) -> Vec<u8> {
    assert!(len >= fgumi_raw_bam::MIN_BAM_RECORD_LEN);
    #[allow(clippy::cast_possible_truncation)]
    let mut b = vec![id as u8; len];
    b[0..4].copy_from_slice(&tid.to_le_bytes());
    b[4..8].copy_from_slice(&pos.to_le_bytes());
    b[14..16].copy_from_slice(&flag.to_le_bytes());
    b[20..28].copy_from_slice(&id.to_le_bytes());
    b
}

/// A merge key type the property runs over.
trait Lane: RawSortKey + Default + Copy + Send + Sync + Debug + 'static {
    /// The spill kind whose runs this key merges.
    const KIND: SpillKeyKind;

    /// A record with key ordinal `ord` (equal ordinals give equal keys;
    /// `unmapped` sorts last) and a `body_len`-byte body tagged with `id`.
    fn record(ord: u64, unmapped: bool, id: u64, body_len: usize) -> (Self, Vec<u8>);
}

impl Lane for RawCoordinateKey {
    const KIND: SpillKeyKind = SpillKeyKind::Coordinate;

    #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
    fn record(ord: u64, unmapped: bool, id: u64, body_len: usize) -> (Self, Vec<u8>) {
        let body = if unmapped {
            bam_body(-1, -1, 0, id, body_len)
        } else {
            let tid = ((ord >> 32) & 0x0FFF_FFFF) as i32;
            let pos = (ord & 0x3FFF_FFFF) as i32;
            let flag = if ord >> 30 & 1 == 1 { fgumi_raw_bam::flags::REVERSE } else { 0 };
            bam_body(tid, pos, flag, id, body_len)
        };
        (RawCoordinateKey::extract_from_record(&body), body)
    }
}

fn key32(ord: u64, unmapped: bool) -> TemplateKey32 {
    if unmapped {
        return TemplateKey32 {
            primary: u64::MAX,
            secondary: u64::MAX,
            opt: u64::MAX,
            name_hash_upper: u64::MAX,
        };
    }
    TemplateKey32 {
        primary: ord >> 16,
        secondary: ord & 0xFFFF,
        opt: ord.rotate_left(3),
        name_hash_upper: ord.rotate_left(7),
    }
}

impl Lane for TemplateKey24 {
    const KIND: SpillKeyKind = SpillKeyKind::TemplateK24;

    fn record(ord: u64, unmapped: bool, id: u64, body_len: usize) -> (Self, Vec<u8>) {
        let k = key32(ord, unmapped);
        let key = TemplateKey24 {
            primary: k.primary,
            secondary: k.secondary,
            name_hash_upper: k.name_hash_upper,
        };
        (key, bam_body(0, 0, 0, id, body_len))
    }
}

impl Lane for CbKey32 {
    const KIND: SpillKeyKind = SpillKeyKind::TemplateCb32;

    fn record(ord: u64, unmapped: bool, id: u64, body_len: usize) -> (Self, Vec<u8>) {
        (CbKey32(key32(ord, unmapped)), bam_body(0, 0, 0, id, body_len))
    }
}

impl Lane for TertKey32 {
    const KIND: SpillKeyKind = SpillKeyKind::TemplateTert32;

    fn record(ord: u64, unmapped: bool, id: u64, body_len: usize) -> (Self, Vec<u8>) {
        (TertKey32(key32(ord, unmapped)), bam_body(0, 0, 0, id, body_len))
    }
}

impl Lane for TemplateKey {
    const KIND: SpillKeyKind = SpillKeyKind::TemplateK40;

    fn record(ord: u64, unmapped: bool, id: u64, body_len: usize) -> (Self, Vec<u8>) {
        let k = key32(ord, unmapped);
        let key = TemplateKey {
            primary: k.primary,
            secondary: k.secondary,
            cb_hash: k.opt,
            tertiary: k.opt.rotate_left(11),
            name_hash_upper: k.name_hash_upper,
        };
        (key, bam_body(0, 0, 0, id, body_len))
    }
}

/// Input shapes, after the prototype's selftest.
#[derive(Clone, Copy, Debug)]
enum Shape {
    HeavyTies,
    AllEqualKeys,
    AllUnmapped,
    WideKeys,
    SingleRun,
    MoreRangesThanRecords,
    OneRecord,
    BigRecordsByteCap,
    ManyRuns,
    EqualKeysCrossFrames,
}

const SHAPES: [Shape; 10] = [
    Shape::HeavyTies,
    Shape::AllEqualKeys,
    Shape::AllUnmapped,
    Shape::WideKeys,
    Shape::SingleRun,
    Shape::MoreRangesThanRecords,
    Shape::OneRecord,
    Shape::BigRecordsByteCap,
    Shape::ManyRuns,
    Shape::EqualKeysCrossFrames,
];

/// What a shape generates.
struct Params {
    records: usize,
    runs: usize,
    keyspace: u64,
    unmapped_pct: u64,
    body: (usize, usize),
    frame_bytes: usize,
}

impl Shape {
    fn params(self) -> Params {
        let p = |records, runs, keyspace, unmapped_pct, body, frame_bytes| Params {
            records,
            runs,
            keyspace,
            unmapped_pct,
            body,
            frame_bytes,
        };
        match self {
            Shape::HeavyTies => p(600, 7, 50, 5, (32, 300), 2048),
            Shape::AllEqualKeys => p(300, 5, 1, 0, (32, 200), 1024),
            Shape::AllUnmapped => p(200, 4, 1000, 100, (32, 200), 1024),
            Shape::WideKeys => p(600, 13, u64::MAX, 2, (32, 600), 4096),
            Shape::SingleRun => p(400, 1, 500, 3, (32, 300), 2048),
            Shape::MoreRangesThanRecords => p(40, 3, 10, 10, (32, 100), 512),
            Shape::OneRecord => p(1, 1, 10, 0, (32, 100), 512),
            Shape::BigRecordsByteCap => p(30, 3, 100, 1, (64 * 1024 + 1, 100 * 1024), 64 * 1024),
            Shape::ManyRuns => p(800, 40, 2000, 5, (32, 200), 1024),
            Shape::EqualKeysCrossFrames => p(300, 2, 3, 0, (32, 200), 512),
        }
    }
}

/// Records of each run, each sorted by key.
type Generated<K> = Vec<Vec<(K, Vec<u8>)>>;

fn generate<K: Lane>(params: &Params, seed: u64) -> Generated<K> {
    let mut rng = Rng(seed);
    let sources = params.runs;
    let mut all: Vec<Vec<(K, Vec<u8>)>> = vec![Vec::new(); sources];
    for id in 0..params.records as u64 {
        let ord = if params.keyspace == u64::MAX { rng.next() } else { rng.below(params.keyspace) };
        let unmapped = rng.below(100) < params.unmapped_pct;
        let len = rng.between(params.body.0, params.body.1);
        let source = usize::try_from(rng.below(sources as u64)).unwrap();
        all[source].push(K::record(ord, unmapped, id, len));
    }
    for source in &mut all {
        source.sort_by(|a, b| a.0.cmp(&b.0));
    }
    all
}

/// Raw bytes the plan counts: framed run records (with their key).
fn raw_bytes<K: Lane>(runs: &[Vec<(K, Vec<u8>)>]) -> u64 {
    let key_bytes = K::KIND.serialized_key_bytes();
    let run: usize = runs.iter().flatten().map(|(_, b)| key_bytes + 4 + b.len()).sum();
    run as u64
}

/// Write `runs` as indexed spill runs under `dir` and plan ranges over them.
fn plan_of<K: Lane>(
    dir: &Path,
    runs: &[Vec<(K, Vec<u8>)>],
    frame_bytes: usize,
    target_range_bytes: u64,
) -> Arc<dyn RangeMerge> {
    let slots = runs
        .iter()
        .enumerate()
        .map(|(i, recs)| {
            write_indexed_run(dir, u32::try_from(i).unwrap(), K::KIND, recs, frame_bytes)
        })
        .collect();
    let total = runs.iter().map(Vec::len).sum::<usize>() as u64;
    let sizing = RangeSizing { target_range_bytes, max_ranges: MAX_RANGES };
    Arc::new(RangePlan::<K>::new(slots, total, sizing).expect("plan"))
}

/// A [`RecordSink`] that frames each body as the range output does.
#[derive(Default)]
struct CollectSink {
    framed: Vec<u8>,
    lens: Vec<u32>,
}

impl RecordSink for CollectSink {
    fn push_record(&mut self, body: &[u8]) -> io::Result<()> {
        let len = u32::try_from(body.len()).expect("a body fits u32");
        self.framed.extend_from_slice(&len.to_le_bytes());
        self.framed.extend_from_slice(body);
        self.lens.push(len);
        Ok(())
    }
}

/// A task merges through `merge` or `merge_into`, never both: switching is an
/// error that ends the task.
#[rstest]
#[case::buffers_then_sink(false)]
#[case::sink_then_buffers(true)]
fn switching_merge_targets_errors(#[case] sink_first: bool) {
    let runs = generate::<RawCoordinateKey>(&Shape::HeavyTies.params(), 7);
    let dir = tempfile::tempdir().expect("tempdir");
    let plan = plan_of(dir.path(), &runs, Shape::HeavyTies.params().frame_bytes, u64::MAX);
    let mut task = Arc::clone(&plan).task(0);
    task.read().unwrap();
    let mut sink = CollectSink::default();
    if sink_first {
        assert!(task.merge_into(&mut sink, 1).unwrap().is_none());
        let err = task.merge(1).expect_err("switching targets");
        assert!(err.to_string().contains("switched"), "{err}");
    } else {
        assert!(task.merge(1).unwrap().is_none());
        let err = task.merge_into(&mut sink, 1).expect_err("switching targets");
        assert!(err.to_string().contains("switched"), "{err}");
    }
}

fn check_ranges<K: Lane>(shape: Shape, seed: u64, target_ranges: u64) -> Result<(), TestCaseError> {
    let params = shape.params();
    let runs = generate::<K>(&params, seed);
    let dir = tempfile::tempdir().expect("tempdir");
    let target = (raw_bytes(&runs) / target_ranges).max(1);
    let plan = plan_of(dir.path(), &runs, params.frame_bytes, target);

    let mut framed = Vec::new();
    let mut lens = Vec::new();
    let mut records = 0u64;
    let (mut kept_run, mut decoded) = (0u64, 0u64);
    // Every third range (j % 3 == 1) merges into a sink, a few records per
    // call, as the merge thread's fused ranges do; the others build their
    // output in buffers. A range is handed the previous range's output
    // buffers, stale contents included, whenever that range built them, as
    // the merge recycles them: so both sink ranges (j = 1, 4, ...) and buffer
    // ranges (j = 3, 6, ...) run on recycled buffers.
    let mut spare: Option<(Vec<u8>, Vec<u32>)> = None;
    let mut buffer_ranges_reusing = 0usize;
    for j in 0..plan.range_count() {
        let sink_range = j % 3 == 1;
        let mut task = Arc::clone(&plan).task(j);
        if let Some((f, l)) = spare.take() {
            buffer_ranges_reusing += usize::from(!sink_range);
            task.reuse_output_buffers(f, l);
        }
        task.read().map_err(|e| TestCaseError::fail(format!("range {j} read: {e}")))?;
        let out = if sink_range {
            let mut sink = CollectSink::default();
            let stats = loop {
                match task.merge_into(&mut sink, 7) {
                    Ok(Some(stats)) => break stats,
                    Ok(None) => {}
                    Err(e) => return Err(TestCaseError::fail(format!("range {j} merge: {e}"))),
                }
            };
            RangeOut { framed: sink.framed, lens: sink.lens, stats }
        } else {
            task.merge(u64::MAX)
                .map_err(|e| TestCaseError::fail(format!("range {j} merge: {e}")))?
                .expect("an unbounded merge finishes the range")
        };
        prop_assert_eq!(out.stats.records, out.lens.len() as u64);
        records += out.stats.records;
        kept_run += out.stats.kept_run_bytes;
        decoded += out.stats.decoded_bytes;
        framed.extend_from_slice(&out.framed);
        lens.extend_from_slice(&out.lens);
        if !sink_range {
            spare = Some((out.framed, out.lens));
        }
    }
    prop_assert!(
        plan.range_count() < 4 || buffer_ranges_reusing > 0,
        "no buffer range ran on a recycled buffer over {} ranges",
        plan.range_count()
    );

    // Shapes with dozens of frames cut several ranges whenever more than one
    // is asked for.
    if target_ranges > 1
        && matches!(
            shape,
            Shape::HeavyTies | Shape::AllEqualKeys | Shape::WideKeys | Shape::ManyRuns
        )
    {
        prop_assert!(plan.range_count() > 1, "{:?}: one range for {} wanted", shape, target_ranges);
    }
    let want = serial_reference(&runs);
    let want_lens: Vec<u32> = want.iter().map(|b| u32::try_from(b.len()).unwrap()).collect();
    prop_assert_eq!(records, plan.total_records());
    prop_assert_eq!(lens, want_lens);
    // Records are counted as kept run bytes in their decoded spill form
    // (`[key][u32 len][body]`), and every kept record was decoded from a
    // frame at least as large. One range decodes
    // every frame once and keeps every record: amplification exactly 1.
    let key_size = crate::merge_driver::slot_key_size::<K>().expect("a fixed-width key");
    let run_decoded: usize = runs.iter().flatten().map(|(_, b)| key_size + 4 + b.len()).sum();
    prop_assert_eq!(kept_run, run_decoded as u64);
    prop_assert!(decoded >= kept_run, "decoded {} < kept run {}", decoded, kept_run);
    if plan.range_count() == 1 {
        prop_assert_eq!(decoded, kept_run);
    }
    prop_assert!(
        framed == frame_bodies(&want),
        "{shape:?} seed {seed}: {} ranges do not concatenate to the serial merge",
        plan.range_count()
    );
    Ok(())
}

fn run_property<K: Lane>(target_ranges: u64) {
    let config = ProptestConfig {
        cases: PROPTEST_CASES,
        failure_persistence: None,
        ..ProptestConfig::default()
    };
    let mut runner = TestRunner::new(config);
    let strategy = (prop::sample::select(SHAPES.to_vec()), any::<u64>());
    runner.run(&strategy, |(shape, seed)| check_ranges::<K>(shape, seed, target_ranges)).unwrap();
}

#[rstest]
#[case::coordinate(SpillKeyKind::Coordinate)]
#[case::k24(SpillKeyKind::TemplateK24)]
#[case::cb32(SpillKeyKind::TemplateCb32)]
#[case::tert32(SpillKeyKind::TemplateTert32)]
#[case::k40(SpillKeyKind::TemplateK40)]
fn ranges_concatenate_to_the_serial_merge(
    #[case] kind: SpillKeyKind,
    #[values(1, 2, 7, 64, 1000)] target_ranges: u64,
) {
    match kind {
        SpillKeyKind::Coordinate => run_property::<RawCoordinateKey>(target_ranges),
        SpillKeyKind::TemplateK24 => run_property::<TemplateKey24>(target_ranges),
        SpillKeyKind::TemplateCb32 => run_property::<CbKey32>(target_ranges),
        SpillKeyKind::TemplateTert32 => run_property::<TertKey32>(target_ranges),
        SpillKeyKind::TemplateK40 => run_property::<TemplateKey>(target_ranges),
        SpillKeyKind::QuerynameLex | SpillKeyKind::QuerynameNatural => {
            unreachable!("queryname runs are not indexed")
        }
    }
}

/// Coordinate records `(key ordinal, body length)` with ids from `first_id`.
fn coordinate_records(spec: &[(u64, usize)], first_id: u64) -> Vec<(RawCoordinateKey, Vec<u8>)> {
    let mut recs: Vec<_> = spec
        .iter()
        .zip(first_id..)
        .map(|(&(ord, len), id)| RawCoordinateKey::record(ord, false, id, len))
        .collect();
    recs.sort_by(|a, b| a.0.cmp(&b.0));
    recs
}

#[derive(Clone, Copy, Debug)]
enum Corrupt {
    BadZstd,
    FrameLenMismatch,
    TruncatedFile,
    RecordCrossesFrame,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Phase {
    Read,
    Merge,
}

/// A run whose two frames split a record: frame 0 holds record A and the
/// first bytes of record B, frame 1 the rest of B. The index claims each frame
/// starts with A's key.
fn write_split_record_run(path: &Path) -> Arc<SortMergeSlot> {
    let recs = coordinate_records(&[(5, 40), (6, 40)], 0);
    let mut raw = Vec::new();
    for (k, body) in &recs {
        crate::frame_keyed_record_into(&mut raw, k, body).unwrap();
    }
    let cut = 44 + 10;
    let mut compressor = SpillBlockCompressor::new(SpillCodec::Zstd, 1).unwrap();
    let mut file = spill_magic(SpillCodec::Zstd).to_vec();
    let mut index = RunIndex::new(8);
    let mut first_key = Vec::new();
    recs[0].0.write_to(&mut first_key).unwrap();
    for (part, records) in [(&raw[..cut], 2), (&raw[cut..], 1)] {
        let mut unit = Vec::new();
        compressor.compress_block_into(part, &mut unit).unwrap();
        let comp = u32::try_from(unit.len()).unwrap();
        index.push_frame(comp, u32::try_from(part.len()).unwrap(), records, &first_key);
        file.extend_from_slice(&unit);
    }
    std::fs::write(path, &file).unwrap();
    index.rebase(spill_magic(SpillCodec::Zstd).len() as u64);
    let slot = crate::open_spill_slot(path, 0, crate::FramesEnd::AtEof).unwrap();
    slot.attach_run_index(Arc::new(index));
    slot
}

#[rstest]
#[case::bad_zstd(Corrupt::BadZstd, Phase::Merge, "zstd")]
#[case::frame_len_mismatch(Corrupt::FrameLenMismatch, Phase::Merge, "does not match the index's")]
#[case::truncated_file(Corrupt::TruncatedFile, Phase::Read, "truncated")]
#[case::record_crosses_frame(Corrupt::RecordCrossesFrame, Phase::Merge, "not record-aligned")]
fn corrupt_frame_errors(#[case] corrupt: Corrupt, #[case] fails_in: Phase, #[case] needle: &str) {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("run-0.zsp");
    let (slot, records) = if let Corrupt::RecordCrossesFrame = corrupt {
        (write_split_record_run(&path), 3)
    } else {
        let spec: Vec<(u64, usize)> = (0..50).map(|i| (i, 60)).collect();
        let recs = coordinate_records(&spec, 0);
        (write_indexed_run(dir.path(), 0, SpillKeyKind::Coordinate, &recs, 512), 50)
    };
    let first_comp_len = u32::try_from(slot.run_index().unwrap().comp_len(0)).unwrap();
    let sizing = RangeSizing { target_range_bytes: u64::MAX, max_ranges: 1 };
    let plan: Arc<dyn RangeMerge> =
        Arc::new(RangePlan::<RawCoordinateKey>::new(vec![slot], records, sizing).unwrap());

    // Corrupt the file after planning, as a disk fault would.
    let file = std::fs::OpenOptions::new().write(true).open(&path).unwrap();
    let body_start = spill_magic(SpillCodec::Zstd).len() as u64;
    match corrupt {
        // The zstd magic number of frame 0.
        Corrupt::BadZstd => file.write_all_at(&[0u8; 4], body_start + 4).unwrap(),
        // Frame 0's prefix claims 4 bytes more than its index allows.
        Corrupt::FrameLenMismatch => {
            file.write_all_at(&first_comp_len.to_le_bytes(), body_start).unwrap();
        }
        Corrupt::TruncatedFile => {
            let len = file.metadata().unwrap().len();
            file.set_len(len - 10).unwrap();
        }
        Corrupt::RecordCrossesFrame => {}
    }

    let mut task = plan.task(0);
    let read = task.read();
    let result = match read {
        Err(e) => Err((Phase::Read, e)),
        Ok(()) => task.merge(u64::MAX).map_err(|e| (Phase::Merge, e)),
    };
    let (phase, err) = result.expect_err("a corrupt run must fail the range");
    assert_eq!(phase, fails_in, "{err}");
    assert!(err.to_string().contains(needle), "{err}");
}

/// Write `recs` as an indexed run, then replace its slot with one over the
/// same file whose index is `edit(index)` sparsified by `factor`.
fn sparse_run(
    dir: &Path,
    file_id: u32,
    recs: &[(RawCoordinateKey, Vec<u8>)],
    factor: u32,
    edit: impl FnOnce(&RunIndex) -> RunIndex,
) -> Arc<SortMergeSlot> {
    let dense = write_indexed_run(dir, file_id, SpillKeyKind::Coordinate, recs, 512);
    let mut index = edit(dense.run_index().unwrap());
    index.sparsify(factor).unwrap();
    let slot = crate::open_spill_slot(
        &dir.join(format!("run-{file_id}.zsp")),
        file_id,
        crate::FramesEnd::AtEof,
    )
    .unwrap();
    slot.attach_run_index(Arc::new(index));
    slot
}

/// Two interleaved runs (even and odd key ordinals, 7 records a frame) whose
/// indexes are sparsified by `factor`: every range merge, including those
/// whose bound falls inside an entry, concatenates to the serial merge.
/// Run 0 gets a short last entry (its frame count is not a multiple of
/// `factor`).
#[rstest]
#[case::by_two(2, 100, 80, 4)]
#[case::by_three_short_last_group(3, 98, 60, 5)]
#[case::by_four_one_range(4, 100, 100, 1)]
#[case::by_one_is_dense(1, 100, 80, 4)]
fn sparse_index_ranges_match_the_serial_merge(
    #[case] factor: u32,
    #[case] n0: u64,
    #[case] n1: u64,
    #[case] ranges: u64,
) {
    let spec = |parity: u64, n: u64| -> Vec<(u64, usize)> {
        (0..n).map(|i| (2 * i + parity, 60)).collect()
    };
    let recs = [coordinate_records(&spec(0, n0), 0), coordinate_records(&spec(1, n1), 1000)];
    let dir = tempfile::tempdir().unwrap();
    let slots: Vec<_> = recs
        .iter()
        .enumerate()
        .map(|(i, r)| sparse_run(dir.path(), u32::try_from(i).unwrap(), r, factor, RunIndex::clone))
        .collect();
    let index0 = Arc::clone(slots[0].run_index().unwrap());
    assert_eq!(index0.stride(), factor);
    assert!(factor == 1 || index0.covered_frames() % u64::from(factor) != 0);

    let total = n0 + n1;
    let target = (total * 72 / ranges).max(1);
    let sizing = RangeSizing { target_range_bytes: target, max_ranges: MAX_RANGES };
    let plan = Arc::new(RangePlan::<RawCoordinateKey>::new(slots, total, sizing).unwrap());
    assert_eq!(plan.range_count() > 1, ranges > 1);

    // Some range bound falls strictly inside one of run 0's entries.
    if factor > 1 && ranges > 1 {
        let inside = (0..plan.range_count()).filter_map(|j| plan.bounds(j).1).any(|hi| {
            (0..index0.len() - 1).any(|i| {
                index0.key_at::<RawCoordinateKey>(i).unwrap() < hi.key
                    && hi.key < index0.key_at::<RawCoordinateKey>(i + 1).unwrap()
            })
        });
        assert!(inside, "no range bound fell inside an entry");
    }

    let span_raw: Vec<u64> = (0..plan.range_count()).map(|j| plan.span_raw_bytes(j)).collect();
    let dyn_plan: Arc<dyn RangeMerge> = plan;
    let (mut framed, mut records) = (Vec::new(), 0u64);
    for (j, &span_raw) in span_raw.iter().enumerate() {
        let mut task = Arc::clone(&dyn_plan).task(j);
        task.read().unwrap();
        let out = task.merge(u64::MAX).unwrap().expect("done");
        // The output was reserved at the span's raw bytes and never grew.
        assert_eq!(out.framed.capacity() as u64, span_raw, "range {j} output capacity");
        records += out.stats.records;
        framed.extend_from_slice(&out.framed);
    }
    assert_eq!(records, total);
    let want = frame_bodies(&serial_reference(&recs));
    assert!(framed == want, "sparse ranges differ from the serial merge");
}

#[derive(Clone, Copy, Debug)]
enum SparseCorrupt {
    /// The second frame of entry 0 claims to run past the entry's end.
    PrefixPastEntryEnd,
    /// Entry 0's index raw length is one more than its frames decode to.
    RawLenDiffers,
}

#[rstest]
#[case::prefix_past_entry_end(SparseCorrupt::PrefixPastEntryEnd, "does not match the index")]
#[case::raw_len_differs(SparseCorrupt::RawLenDiffers, "the index says")]
fn sparse_entry_tiling_errors(#[case] corrupt: SparseCorrupt, #[case] needle: &str) {
    let dir = tempfile::tempdir().unwrap();
    let spec: Vec<(u64, usize)> = (0..60).map(|i| (i, 60)).collect();
    let recs = coordinate_records(&spec, 0);
    let edit = |dense: &RunIndex| {
        let mut index = RunIndex::new(dense.key_width());
        for i in 0..dense.len() {
            let bump = u32::from(matches!(corrupt, SparseCorrupt::RawLenDiffers) && i == 0);
            index.push_frame(
                u32::try_from(dense.comp_len(i)).unwrap(),
                dense.raw_len(i) + bump,
                dense.records(i),
                dense.key_bytes(i),
            );
        }
        index.rebase(dense.offset(0));
        index
    };
    let slot = sparse_run(dir.path(), 0, &recs, 2, edit);
    let index = Arc::clone(slot.run_index().unwrap());
    assert!(index.len() > 1 && index.stride() == 2);
    let sizing = RangeSizing { target_range_bytes: u64::MAX, max_ranges: 1 };
    let plan: Arc<dyn RangeMerge> =
        Arc::new(RangePlan::<RawCoordinateKey>::new(vec![slot], 60, sizing).unwrap());

    if let SparseCorrupt::PrefixPastEntryEnd = corrupt {
        // Entry 0's frames: [prefix][zstd] twice. Read the first prefix to find
        // the second, and make that one claim the whole entry's length.
        let path = dir.path().join("run-0.zsp");
        let file = std::fs::OpenOptions::new().read(true).write(true).open(&path).unwrap();
        let mut prefix = [0u8; 4];
        file.read_exact_at(&mut prefix, index.offset(0)).unwrap();
        let second = index.offset(0) + 4 + u64::from(u32::from_le_bytes(prefix));
        assert!(second < index.end(0), "entry 0 holds two frames");
        let claim = u32::try_from(index.comp_len(0)).unwrap();
        file.write_all_at(&claim.to_le_bytes(), second).unwrap();
    }

    let mut task = plan.task(0);
    task.read().unwrap();
    let err = task.merge(u64::MAX).expect_err("a mistiled entry must fail the range");
    assert!(err.to_string().contains("corrupt spill frame 0"), "{err}");
    assert!(err.to_string().contains(needle), "{err}");
}

#[rstest]
#[case::one(1, 999)]
#[case::sixteen(16, 62)]
#[case::all_but_one(999, 1)]
#[case::exact(1000, 0)]
#[case::unbounded(u64::MAX, 0)]
fn merge_yields_in_slices(#[case] max_records: u64, #[case] expected_nones: usize) {
    let params = Shape::HeavyTies.params();
    let params = Params { records: 1000, ..params };
    let runs = generate::<RawCoordinateKey>(&params, 7);
    let dir = tempfile::tempdir().unwrap();
    let plan = plan_of(dir.path(), &runs, params.frame_bytes, u64::MAX);
    assert_eq!(plan.range_count(), 1);

    let mut whole = Arc::clone(&plan).task(0);
    whole.read().unwrap();
    let whole = whole.merge(u64::MAX).unwrap().expect("done");

    let mut task = Arc::clone(&plan).task(0);
    task.read().unwrap();
    let mut nones = 0;
    let sliced = loop {
        match task.merge(max_records).unwrap() {
            Some(out) => break out,
            None => nones += 1,
        }
    };
    assert_eq!(nones, expected_nones);
    assert_eq!(sliced.stats.records, 1000);
    assert_eq!(sliced.lens, whole.lens);
    assert!(sliced.framed == whole.framed, "a sliced merge differs from the whole merge");
}

/// Merge every range of `plan` and return the concatenated output and each
/// range's record count, checking each range's span against the plan's
/// range-size bound for `sizing`.
fn merge_all_bounded(
    plan: &Arc<RangePlan<RawCoordinateKey>>,
    sizing: RangeSizing,
) -> Result<(Vec<u8>, Vec<u64>), TestCaseError> {
    let limit = plan.range_bound(sizing);
    let dyn_plan: Arc<dyn RangeMerge> = Arc::clone(plan) as Arc<dyn RangeMerge>;
    let (mut framed, mut records) = (Vec::new(), Vec::new());
    for j in 0..dyn_plan.range_count() {
        let raw = plan.span_raw_bytes(j);
        prop_assert!(raw <= limit, "range {} spans {} raw bytes > {}", j, raw, limit);
        let mut task = Arc::clone(&dyn_plan).task(j);
        task.read().map_err(|e| TestCaseError::fail(format!("range {j} read: {e}")))?;
        let out = task
            .merge(u64::MAX)
            .map_err(|e| TestCaseError::fail(format!("range {j} merge: {e}")))?
            .expect("an unbounded merge finishes the range");
        records.push(out.stats.records);
        framed.extend_from_slice(&out.framed);
    }
    Ok((framed, records))
}

/// Coordinate inputs dominated by one key.
#[derive(Clone, Copy, Debug)]
enum TieShape {
    /// Every record at one mapped position.
    AllEqualKeys,
    /// A few mapped records, then an unmapped tail (one key) of most records.
    UnmappedTail,
    /// Scattered positions with a deep pileup at one position in the middle.
    Pileup,
    /// Scattered positions, a pileup and an unmapped tail.
    Mix,
}

impl TieShape {
    /// Record `i`'s key ordinal and whether it is unmapped.
    fn key(self, rng: &mut Rng) -> (u64, bool) {
        let scattered = rng.below(10_000);
        let roll = rng.below(100);
        match self {
            TieShape::AllEqualKeys => (500, false),
            TieShape::UnmappedTail => (scattered, roll < 80),
            TieShape::Pileup => (if roll < 70 { 5_000 } else { scattered }, false),
            TieShape::Mix => (if roll < 35 { 5_000 } else { scattered }, roll >= 65),
        }
    }
}

/// A group of equal keys far larger than the range target is cut into many
/// ranges (it was one range when bounds fell only between distinct keys),
/// each within `ceil(total / R) + k x W` raw bytes, and the ranges still
/// concatenate to the serial merge.
#[rstest]
#[case::all_equal_keys(TieShape::AllEqualKeys, 4, 1)]
#[case::all_equal_keys_sparse(TieShape::AllEqualKeys, 4, 3)]
#[case::unmapped_tail(TieShape::UnmappedTail, 5, 1)]
#[case::unmapped_tail_sparse(TieShape::UnmappedTail, 5, 2)]
#[case::pileup(TieShape::Pileup, 3, 1)]
#[case::pileup_sparse(TieShape::Pileup, 3, 4)]
#[case::mix(TieShape::Mix, 6, 1)]
#[case::mix_sparse(TieShape::Mix, 6, 2)]
fn equal_key_groups_split_into_bounded_ranges(
    #[case] shape: TieShape,
    #[case] runs: usize,
    #[case] stride: u32,
) {
    const RECORDS: u64 = 4000;
    const RANGES: u64 = 16;
    let mut rng = Rng(0xB1);
    let mut recs: Vec<Vec<(RawCoordinateKey, Vec<u8>)>> = vec![Vec::new(); runs];
    for id in 0..RECORDS {
        let (ord, unmapped) = shape.key(&mut rng);
        let len = rng.between(32, 200);
        let run = usize::try_from(rng.below(runs as u64)).unwrap();
        recs[run].push(RawCoordinateKey::record(ord, unmapped, id, len));
    }
    for run in &mut recs {
        run.sort_by(|a, b| a.0.cmp(&b.0));
    }
    let dir = tempfile::tempdir().unwrap();
    let slots: Vec<_> = recs
        .iter()
        .enumerate()
        .map(|(i, r)| sparse_run(dir.path(), u32::try_from(i).unwrap(), r, stride, RunIndex::clone))
        .collect();
    let target = raw_bytes(&recs) / RANGES;
    let sizing = RangeSizing { target_range_bytes: target, max_ranges: MAX_RANGES };
    let plan = Arc::new(RangePlan::<RawCoordinateKey>::new(slots, RECORDS, sizing).unwrap());
    // Bounds land inside the dominant key's group.
    assert!(plan.range_count() as u64 >= RANGES - 1, "{} ranges", plan.range_count());

    let (framed, records) = merge_all_bounded(&plan, sizing).unwrap();
    assert_eq!(records.iter().sum::<u64>(), RECORDS);
    let largest = *records.iter().max().unwrap();
    assert!(largest * 4 <= RECORDS, "largest range holds {largest} of {RECORDS} records");
    assert!(framed == frame_bodies(&serial_reference(&recs)), "{shape:?}: ranges differ");
}

/// Per-run records `(key ordinal, unmapped, body length)`, the index stride
/// and the number of ranges to aim for.
type TieInputs = (Vec<Vec<(u64, bool, usize)>>, u32, u64);

/// Keys from a tiny domain (three positions and unmapped), so nearly every
/// range bound falls inside a group of equal keys.
fn tie_inputs() -> impl Strategy<Value = TieInputs> {
    (
        prop::collection::vec(
            prop::collection::vec((0u64..3, prop::bool::weighted(0.3), 32usize..150), 0..150),
            1..6,
        ),
        1u32..=4,
        1u64..40,
    )
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 64, failure_persistence: None, ..ProptestConfig::default() })]

    /// Heavily duplicated keys, random run counts and sizes, strides 1 to 4:
    /// ranges concatenate to the serial merge and stay within the bound.
    #[test]
    fn heavy_duplicate_ranges_match_serial_and_stay_bounded(
        (spec, stride, ranges) in tie_inputs()
    ) {
        let mut id = 0u64;
        let recs: Vec<Vec<(RawCoordinateKey, Vec<u8>)>> = spec
            .iter()
            .map(|run| {
                let mut recs: Vec<_> = run
                    .iter()
                    .map(|&(ord, unmapped, len)| {
                        id += 1;
                        RawCoordinateKey::record(ord, unmapped, id, len)
                    })
                    .collect();
                recs.sort_by(|a, b| a.0.cmp(&b.0));
                recs
            })
            .collect();
        let dir = tempfile::tempdir().unwrap();
        let slots: Vec<_> = recs
            .iter()
            .enumerate()
            .map(|(i, r)| {
                sparse_run(dir.path(), u32::try_from(i).unwrap(), r, stride, RunIndex::clone)
            })
            .collect();
        let total = recs.iter().map(Vec::len).sum::<usize>() as u64;
        let target = (raw_bytes(&recs) / ranges).max(1);
        let sizing = RangeSizing { target_range_bytes: target, max_ranges: MAX_RANGES };
        let plan = Arc::new(RangePlan::<RawCoordinateKey>::new(slots, total, sizing).unwrap());
        let (framed, records) = merge_all_bounded(&plan, sizing)?;
        prop_assert_eq!(records.iter().sum::<u64>(), total);
        prop_assert!(framed == frame_bodies(&serial_reference(&recs)), "ranges differ from serial");
    }
}
