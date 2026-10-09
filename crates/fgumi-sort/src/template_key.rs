//! Template-coordinate sort key: the per-sort seed, the provisioned context
//! (variant + dropped-lane baseline, immutable and `Arc`-shared), and the one
//! extraction implementation used by both the per-record accumulator and the
//! batched extraction over arena byte ranges. A batch is a
//! `(byte_start, byte_end, n_records)` range of `[u32 LE block_size][body]*`
//! frames in the shared arena, re-walked here instead of carried as an extents
//! array. Its refs come back as [`ExtractedRefs`], which the ref container in
//! [`template_arena`](crate::template_arena) appends; this module depends on
//! nothing in that one.

use std::io;
use std::sync::Arc;

use fgumi_raw_bam::{MIN_BAM_RECORD_LEN, SamTag};
use noodles::sam::Header;

use crate::arena_pool::PooledSegmentedBuf;
use crate::external::{
    DroppedLaneViolation, KeyTypesSpec, LibraryLookup, TemplateKeyVariant, cb_hasher,
    dropped_lane_error, extract_template_key_inline, select_template_variant, verify_dropped_lanes,
    violation_read_name,
};
use crate::inline::{
    CbKey32, TemplateKey, TemplateKey24, TemplateKey40, TemplateLaneKey, TemplateRecordRef,
    TertKey32,
};
use crate::prefetch::{KEY_PREFETCH_DISTANCE, prefetch_read_l1};

/// Run `$body` against the inner collection of whichever arm `$self` holds
/// (works for both [`ExtractedRefs`] and
/// [`TemplateArenaRefs`](crate::template_arena::TemplateArenaRefs)). `$body` is
/// monomorphized per arm, so a fifth lane is one edit here, not one per method.
macro_rules! with_lane_refs {
    ($ty:ident, $self:expr, $refs:ident => $body:expr) => {
        match $self {
            $ty::K24($refs) => $body,
            $ty::Cb32($refs) => $body,
            $ty::Tert32($refs) => $body,
            $ty::K40($refs) => $body,
        }
    };
}

/// The [`TemplateKeyVariant`] an arm of `$ty` (either lane-ref enum) encodes.
macro_rules! lane_variant {
    ($ty:ident, $self:expr) => {
        match $self {
            $ty::K24(_) => $crate::template_key::variant_of(false, false),
            $ty::Cb32(_) => $crate::template_key::variant_of(true, false),
            $ty::Tert32(_) => $crate::template_key::variant_of(false, true),
            $ty::K40(_) => $crate::template_key::variant_of(true, true),
        }
    };
}

pub(crate) use {lane_variant, with_lane_refs};

pub(crate) const fn variant_of(cb: bool, tertiary: bool) -> TemplateKeyVariant {
    TemplateKeyVariant { cb, tertiary }
}

/// One batch's refs as extracted by
/// [`TemplateKeyContext::extract_batch_erased`], in record order. Same four arms
/// as [`TemplateArenaRefs`](crate::template_arena::TemplateArenaRefs), which
/// appends them.
#[derive(Debug)]
pub enum ExtractedRefs {
    /// 24-byte lane.
    K24(Box<[TemplateRecordRef<TemplateKey24>]>),
    /// 32-byte lane carrying `cb_hash`.
    Cb32(Box<[TemplateRecordRef<CbKey32>]>),
    /// 32-byte lane carrying the tertiary word.
    Tert32(Box<[TemplateRecordRef<TertKey32>]>),
    /// Full 40-byte key.
    K40(Box<[TemplateRecordRef<TemplateKey40>]>),
}

impl ExtractedRefs {
    /// Number of refs in the batch.
    #[must_use]
    pub fn len(&self) -> usize {
        with_lane_refs!(Self, self, v => v.len())
    }

    /// `true` iff the batch is empty (a run's empty terminal batch).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Logical bytes the refs occupy (`len × ref_width`). This is what the refs
    /// hold, not what a container they are appended to has allocated: a growing
    /// `Vec` can hold up to about twice its length in capacity.
    #[must_use]
    pub fn byte_len(&self) -> usize {
        self.len() * self.variant().ref_width()
    }

    /// The variant this batch's arm encodes.
    #[must_use]
    pub fn variant(&self) -> TemplateKeyVariant {
        lane_variant!(Self, self)
    }
}

/// Per-sort configuration the template key needs before the first record is
/// seen: library lookup, cell-barcode tag + fixed-seed hasher, `--key-types`.
#[derive(Clone)]
pub struct TemplateKeySeed {
    lib_lookup: LibraryLookup,
    cell_tag: Option<SamTag>,
    /// Fixed seed — feeds the sort key. Never reseed (byte identity).
    cb_hasher: ahash::RandomState,
    key_types: KeyTypesSpec,
    header_library_varies: bool,
}

impl TemplateKeySeed {
    /// Derive the seed from the BAM header: the library lookup and CB hasher are
    /// built exactly as
    /// [`RawExternalSorter::into_template_chunk_sorter`](crate::RawExternalSorter::into_template_chunk_sorter)
    /// builds them, so both paths key records identically.
    #[must_use]
    pub fn from_header(header: &Header, cell_tag: Option<SamTag>, key_types: KeyTypesSpec) -> Self {
        let lib_lookup = LibraryLookup::from_header(header);
        let header_library_varies = lib_lookup.distinct_header_ordinals() > 1;
        Self { lib_lookup, cell_tag, cb_hasher: cb_hasher(), key_types, header_library_varies }
    }

    /// Extract the sort's FIRST record's full key, choose the lane variant from
    /// it (`select_template_variant`, the same choice the owned sorter makes on
    /// its first push), and freeze both into a shared context. Called once per
    /// sort, with the sort's first record in input order: the variant and the
    /// dropped-lane baseline are properties of the sort, so every worker and
    /// every batch must see the same ones.
    ///
    /// # Panics
    ///
    /// Panics if `first_body` is shorter than the 32 fixed BAM bytes
    /// ([`MIN_BAM_RECORD_LEN`]); the key extractor indexes those fields. Every
    /// front checks this first: the arena boundary scan rejects such a frame,
    /// [`TemplateArenaAccumulator::push`](crate::TemplateArenaAccumulator::push)
    /// returns an error, and `extract_batch` returns `InvalidData`.
    #[must_use]
    pub fn provision(&self, first_body: &[u8]) -> Arc<TemplateKeyContext> {
        let first_key = extract_template_key_inline(
            first_body,
            &self.lib_lookup,
            self.cell_tag,
            &self.cb_hasher,
        );
        let variant =
            select_template_variant(Some(&first_key), self.key_types, self.header_library_varies);
        Arc::new(TemplateKeyContext {
            lib_lookup: self.lib_lookup.clone(),
            cell_tag: self.cell_tag,
            cb_hasher: self.cb_hasher.clone(),
            first_key,
            variant,
        })
    }
}

/// Everything a key extraction needs, immutable after provisioning and shared
/// by `Arc` with every batch.
pub struct TemplateKeyContext {
    pub(crate) lib_lookup: LibraryLookup,
    pub(crate) cell_tag: Option<SamTag>,
    pub(crate) cb_hasher: ahash::RandomState,
    /// The first record's full key: the dropped-lane verify baseline.
    pub first_key: TemplateKey,
    /// The lanes the narrowed key retains.
    pub variant: TemplateKeyVariant,
}

/// One batch's refs and its lowest-indexed violation. A batch whose
/// `violation` is `Some` must fail the run, never be sealed (see
/// [`TemplateKeyContext::extract_batch`]).
#[derive(Debug)]
pub struct ExtractedBatch<K: TemplateLaneKey> {
    /// Refs in record order, offsets global to the arena.
    pub refs: Box<[TemplateRecordRef<K>]>,
    /// The lowest-indexed dropped-lane violation in the batch, if any.
    pub violation: Option<KeyViolation>,
}

/// A dropped-lane violation tagged with the record that carried it.
/// `record_index` is run-relative; batches are consumed in input order, so the
/// first violation seen is the first in the input.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KeyViolation {
    /// Run-relative index of the offending record.
    pub record_index: u64,
    /// Which lane disagreed with the first record.
    pub violation: DroppedLaneViolation,
    /// The offending record's read name (lossy UTF-8).
    pub name: String,
}

impl KeyViolation {
    /// The exact `io::Error` the serial path produces:
    /// `InvalidData` carrying `format!("{:#}", dropped_lane_error(name, v))`.
    #[must_use]
    pub fn into_io_error(self) -> io::Error {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{:#}", dropped_lane_error(&self.name, self.violation)),
        )
    }
}

/// How many refs `extract_batch` reserves for a batch claiming `n_records`
/// frames in `span_len` bytes. `n_records` is not trusted: the span cannot hold
/// more than `span_len / (4 + MIN_BAM_RECORD_LEN)` frames, and a corrupt count
/// must surface as `InvalidData` ("fewer frames than `n_records`"), not as an
/// allocation abort.
fn ref_reservation(n_records: u32, span_len: usize) -> usize {
    (n_records as usize).min(span_len / (4 + MIN_BAM_RECORD_LEN))
}

/// An `InvalidData` error naming the malformed batch by its first record.
fn malformed(first_record_index: u64, what: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("template key batch (first record {first_record_index}): {what}"),
    )
}

impl TemplateKeyContext {
    /// The single extraction implementation: full key + dropped-lane verify
    /// against `first_key`. Both the per-record accumulator and the batched
    /// path call this, so they cannot drift.
    ///
    /// # Panics
    ///
    /// Panics if `body` is shorter than [`MIN_BAM_RECORD_LEN`], as
    /// [`TemplateKeySeed::provision`] does; both callers check first.
    #[inline]
    #[must_use]
    pub fn extract_one(&self, body: &[u8]) -> (TemplateKey, Option<DroppedLaneViolation>) {
        let full =
            extract_template_key_inline(body, &self.lib_lookup, self.cell_tag, &self.cb_hasher);
        let violation = verify_dropped_lanes(&self.first_key, &full, self.variant);
        (full, violation)
    }

    /// Re-walk `[byte_start, byte_end)` of `arena` as exactly `n_records`
    /// `[u32 LE block_size][body]` frames; extract every key; keep the LOWEST
    /// indexed violation; build `K`-lane refs whose `offset` is the body's
    /// global arena offset. Prefetches `KEY_PREFETCH_DISTANCE` ahead within the
    /// batch's span, so the last few records of a batch get no lead
    /// (negligible at 4096 records).
    ///
    /// Every record gets a ref, a violating one included: the violation is
    /// deferred because batches run out of order, and the caller must fail the
    /// run with the LOWEST-indexed violation across all of them. A batch whose
    /// `violation` is `Some` is therefore never sealed — the caller turns the
    /// run's first violation into its error with [`KeyViolation::into_io_error`],
    /// exactly as the per-record path rejects the record it cannot key.
    ///
    /// # Errors
    ///
    /// `InvalidData` naming the batch, and the frame by its record index, byte
    /// offset and `block_size`, if a frame is shorter than the 32 fixed BAM
    /// bytes, runs past `byte_end`, fewer than `n_records` frames are present,
    /// or bytes remain after the last frame — impossible under the boundary
    /// step's cut rule, so never silently truncated.
    ///
    /// # Panics
    ///
    /// Panics if `[byte_start, byte_end)` does not lie within the written bytes
    /// of one arena segment (the arena front reserves one segment per run, and
    /// the boundary step cuts batches only from bytes already inflated into it).
    pub fn extract_batch<K: TemplateLaneKey>(
        &self,
        arena: &PooledSegmentedBuf,
        byte_start: u64,
        byte_end: u64,
        n_records: u32,
        first_record_index: u64,
    ) -> io::Result<ExtractedBatch<K>> {
        let start = usize::try_from(byte_start)
            .map_err(|_| malformed(first_record_index, "start overflows usize"))?;
        let end = usize::try_from(byte_end)
            .map_err(|_| malformed(first_record_index, "end overflows usize"))?;
        if end < start {
            return Err(malformed(first_record_index, "end precedes start"));
        }
        let span = arena.slice(start, end - start);
        let mut refs: Vec<TemplateRecordRef<K>> =
            Vec::with_capacity(ref_reservation(n_records, span.len()));
        let mut violation: Option<KeyViolation> = None;
        let mut cur = 0usize;
        for i in 0..u64::from(n_records) {
            let record = first_record_index + i;
            let at = byte_start + cur as u64;
            if cur + 4 > span.len() {
                return Err(malformed(
                    first_record_index,
                    &format!(
                        "fewer frames than n_records ({n_records}): record {record} has no \
                         block_size at byte {at}"
                    ),
                ));
            }
            let bs = u32::from_le_bytes([span[cur], span[cur + 1], span[cur + 2], span[cur + 3]]);
            if (bs as usize) < MIN_BAM_RECORD_LEN {
                return Err(malformed(
                    first_record_index,
                    &format!(
                        "record {record} at byte {at}: block_size {bs} is shorter than a BAM \
                         record"
                    ),
                ));
            }
            let body_start = cur + 4;
            let body_end = body_start + bs as usize;
            if body_end > span.len() {
                return Err(malformed(
                    first_record_index,
                    &format!(
                        "record {record} at byte {at}: block_size {bs} runs past byte_end \
                         {byte_end}"
                    ),
                ));
            }
            if let Some(ahead) = span.get(cur + KEY_PREFETCH_DISTANCE) {
                prefetch_read_l1(ahead);
            }
            let body = &span[body_start..body_end];
            let (full, lane) = self.extract_one(body);
            if violation.is_none()
                && let Some(v) = lane
            {
                violation = Some(KeyViolation {
                    record_index: record,
                    violation: v,
                    name: violation_read_name(body),
                });
            }
            refs.push(TemplateRecordRef {
                key: K::from_full(&full),
                offset: byte_start + body_start as u64,
                len: bs,
                padding: 0,
            });
            cur = body_end;
        }
        if cur != span.len() {
            return Err(malformed(
                first_record_index,
                &format!(
                    "{} bytes remain after the last frame (byte {})",
                    span.len() - cur,
                    byte_start + cur as u64
                ),
            ));
        }
        Ok(ExtractedBatch { refs: refs.into_boxed_slice(), violation })
    }

    /// Variant-erased dispatch over the four lane widths: the caller holds one
    /// context per sort and need not name the lane type the first record chose.
    /// The violation carries the same contract as [`Self::extract_batch`]'s: a
    /// batch that reports one is never sealed.
    ///
    /// # Errors
    ///
    /// As [`Self::extract_batch`].
    pub fn extract_batch_erased(
        &self,
        arena: &PooledSegmentedBuf,
        byte_start: u64,
        byte_end: u64,
        n_records: u32,
        first_record_index: u64,
    ) -> io::Result<(ExtractedRefs, Option<KeyViolation>)> {
        macro_rules! run {
            ($k:ty, $arm:ident) => {{
                let b = self.extract_batch::<$k>(
                    arena,
                    byte_start,
                    byte_end,
                    n_records,
                    first_record_index,
                )?;
                (ExtractedRefs::$arm(b.refs), b.violation)
            }};
        }
        Ok(match (self.variant.cb, self.variant.tertiary) {
            (false, false) => run!(TemplateKey24, K24),
            (true, false) => run!(CbKey32, Cb32),
            (false, true) => run!(TertKey32, Tert32),
            (true, true) => run!(TemplateKey40, K40),
        })
    }
}

#[cfg(test)]
mod tests {
    use fgumi_raw_bam::SamTag;
    use fgumi_raw_bam::testutil::make_bam_bytes;
    use noodles::sam::Header;

    use super::*;
    use crate::arena_pool::PooledSegmentedBuf;
    use crate::external::{KeyTypesSpec, extract_template_key_inline};
    use crate::inline::{TemplateKey24, TemplateKey40};
    use crate::segmented_buf::SegmentedBuf;

    /// `CB:Z:<value>` aux bytes.
    fn cb_aux(value: &[u8]) -> Vec<u8> {
        let mut aux = b"CBZ".to_vec();
        aux.extend_from_slice(value);
        aux.push(0);
        aux
    }

    /// `MI:Z:<value>` aux bytes (feeds the tertiary lane).
    fn mi_aux(value: &[u8]) -> Vec<u8> {
        let mut aux = b"MIZ".to_vec();
        aux.extend_from_slice(value);
        aux.push(0);
        aux
    }

    /// A mapped, paired record body (names chosen so `len + 1` is a multiple of 4,
    /// per `make_bam_bytes`'s alignment note).
    fn rec(pos: i32, name: &[u8], aux: &[u8]) -> Vec<u8> {
        make_bam_bytes(0, pos, 0x3, name, &[], 20, 0, pos + 50, aux)
    }

    /// One-segment arena holding `base` bytes of padding, then the records as
    /// `[u32 LE len][body]*` frames. Returns the arena and `(byte_start, byte_end)`.
    fn arena_with(records: &[Vec<u8>], base: usize) -> (PooledSegmentedBuf, u64, u64) {
        let mut frames = Vec::new();
        for r in records {
            frames.extend_from_slice(&u32::try_from(r.len()).unwrap().to_le_bytes());
            frames.extend_from_slice(r);
        }
        let mut buf = SegmentedBuf::with_capacity(0, base + frames.len() + 64);
        if base > 0 {
            let off = buf.extend_from_slice(&vec![0u8; base]);
            assert_eq!(off, 0);
        }
        let start = buf.extend_from_slice(&frames);
        assert_eq!(start, base, "frames must be contiguous right after the padding");
        (PooledSegmentedBuf::unpooled(buf), start as u64, (start + frames.len()) as u64)
    }

    fn seed(cell_tag: Option<SamTag>, key_types: KeyTypesSpec) -> TemplateKeySeed {
        TemplateKeySeed::from_header(&Header::default(), cell_tag, key_types)
    }

    /// For every lane width, each ref's key must equal the narrowing of the key
    /// `extract_template_key_inline` computes directly (an oracle outside the
    /// code under test), and its offset/len must point at the record's body.
    /// Records carry CB and MI so the Cb32/Tert32/K40 lanes hold non-zero values;
    /// every record has the same CB/MI as the first, so no lane is dropped.
    #[rstest::rstest]
    #[case::k24(KeyTypesSpec::None)]
    #[case::cb32(KeyTypesSpec::Explicit { cb: true, tertiary: false })]
    #[case::tert32(KeyTypesSpec::Explicit { cb: false, tertiary: true })]
    #[case::k40(KeyTypesSpec::Full)]
    fn batch_keys_equal_serial_keys(#[case] spec: KeyTypesSpec) {
        let mut aux = cb_aux(b"ACGT");
        aux.extend_from_slice(&mi_aux(b"7"));
        let records: Vec<Vec<u8>> = [5, 3, 9, 3, 1, 8, 2, 7]
            .iter()
            .enumerate()
            .map(|(i, &p)| rec(100 + p, format!("k{i}a").as_bytes(), &aux))
            .collect();
        let ctx = seed(Some(SamTag::CB), spec).provision(&records[0]);
        let (arena, start, end) = arena_with(&records, 0);
        macro_rules! check {
            ($k:ty) => {{
                let batch = ctx.extract_batch::<$k>(&arena, start, end, 8, 0).unwrap();
                assert!(batch.violation.is_none());
                assert_eq!(batch.refs.len(), records.len());
                let mut off = start + 4;
                for (r, body) in batch.refs.iter().zip(&records) {
                    let direct = extract_template_key_inline(
                        body,
                        &ctx.lib_lookup,
                        ctx.cell_tag,
                        &ctx.cb_hasher,
                    );
                    assert_eq!(r.key, <$k as TemplateLaneKey>::from_full(&direct));
                    assert_eq!((r.offset, r.len as usize), (off, body.len()));
                    off += body.len() as u64 + 4;
                }
            }};
        }
        match (ctx.variant.cb, ctx.variant.tertiary) {
            (false, false) => check!(TemplateKey24),
            (true, false) => check!(crate::inline::CbKey32),
            (false, true) => check!(crate::inline::TertKey32),
            (true, true) => check!(TemplateKey40),
        }
        for r in &records {
            let direct =
                extract_template_key_inline(r, &ctx.lib_lookup, ctx.cell_tag, &ctx.cb_hasher);
            assert_eq!(ctx.extract_one(r).0, direct, "extract_one is exactly the inline extractor");
        }
    }

    /// Two offenders; the earlier one is reported, at its run-relative index.
    #[test]
    fn batch_reports_the_lowest_indexed_violation() {
        let records = vec![
            rec(100, b"r0a", &cb_aux(b"AAAA")),
            rec(101, b"r1a", &cb_aux(b"AAAA")),
            rec(102, b"r2a", &cb_aux(b"CCCC")),
            rec(103, b"r3a", &cb_aux(b"AAAA")),
            rec(104, b"r4a", &cb_aux(b"GGGG")),
        ];
        // `None` drops the CB lane, so differing barcodes violate.
        let ctx = seed(Some(SamTag::CB), KeyTypesSpec::None).provision(&records[0]);
        let (arena, start, end) = arena_with(&records, 0);
        let batch = ctx.extract_batch::<TemplateKey24>(&arena, start, end, 5, 1000).unwrap();
        let v = batch.violation.expect("the dropped CB lane must be reported");
        assert_eq!(v.record_index, 1002, "r2a is record 2 of a batch starting at 1000");
        assert_eq!(v.name, "r2a");
        assert_eq!(batch.refs.len(), 5, "every record still gets a ref");
    }

    /// A frame whose `l_read_name` runs past its body still keys (the extractor
    /// bounds-checks the name) and, when it drops a lane, is reported under a
    /// placeholder name — the violation branch never slices past the body.
    #[test]
    fn a_violation_whose_read_name_overruns_its_frame_is_reported_not_a_panic() {
        let mut bad = rec(101, b"r1a", &[]);
        bad[8] = 200;
        assert!(32 + 199 > bad.len(), "the name must overrun the body");
        let records = vec![rec(100, b"r0a", &cb_aux(b"AAAA")), bad];
        // `None` drops the CB lane; the bad record has no CB, so it violates.
        let ctx = seed(Some(SamTag::CB), KeyTypesSpec::None).provision(&records[0]);
        let (arena, start, end) = arena_with(&records, 0);
        let batch = ctx.extract_batch::<TemplateKey24>(&arena, start, end, 2, 0).unwrap();
        let v = batch.violation.expect("the missing CB is a dropped-lane violation");
        assert_eq!((v.record_index, v.name.as_str()), (1, "<malformed read name>"));
    }

    /// Real batches start at or after `FRONT_REGION` (8 MiB); ref offsets are global
    /// arena offsets of the BODY (length prefix excluded).
    #[test]
    fn extract_batch_offsets_are_global_arena_offsets() {
        const BASE: usize = 8 * 1024 * 1024;
        let records: Vec<Vec<u8>> =
            (0..4).map(|i| rec(200 + i, format!("s{i}a").as_bytes(), &[])).collect();
        let ctx = seed(None, KeyTypesSpec::Full).provision(&records[0]);
        let (arena, start, end) = arena_with(&records, BASE);
        assert_eq!(start, BASE as u64);
        let batch = ctx.extract_batch::<TemplateKey40>(&arena, start, end, 4, 0).unwrap();
        assert_eq!(batch.refs.len(), records.len());
        for (r, original) in batch.refs.iter().zip(&records) {
            assert!(r.offset >= BASE as u64 + 4);
            assert_eq!(
                arena.slice(usize::try_from(r.offset).unwrap(), r.len as usize),
                &original[..]
            );
        }
    }

    /// Malformed batches are an `InvalidData` error naming the batch and the
    /// offending frame, never a panic or a short ref array. Each case pins the
    /// guard that fired by its message.
    #[rstest::rstest]
    #[case::more_records_than_frames(
        4,
        0,
        "fewer frames than n_records (4): record 80 has no block_size"
    )]
    #[case::frame_runs_past_end(3, 7, "record 79 at byte ")]
    #[case::trailing_bytes_after_last_frame(2, 0, "bytes remain after the last frame")]
    fn extract_batch_rejects_malformed_frames(
        #[case] n_records: u32,
        #[case] trim_end: u64,
        #[case] want: &str,
    ) {
        let records: Vec<Vec<u8>> =
            (0..3).map(|i| rec(300 + i, format!("m{i}a").as_bytes(), &[])).collect();
        let ctx = seed(None, KeyTypesSpec::Full).provision(&records[0]);
        let (arena, start, end) = arena_with(&records, 0);
        let Err(err) =
            ctx.extract_batch::<TemplateKey40>(&arena, start, end - trim_end, n_records, 77)
        else {
            panic!("malformed batch must be rejected");
        };
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        let msg = err.to_string();
        assert!(msg.starts_with("template key batch (first record 77): "), "{msg}");
        assert!(msg.contains(want), "{msg}");
        if trim_end > 0 {
            assert!(msg.contains(&format!("runs past byte_end {}", end - trim_end)), "{msg}");
        }
    }

    /// A corrupt `n_records` near `u32::MAX` is an `InvalidData` error, not an
    /// allocation abort: the ref reservation is capped by what the span holds.
    #[test]
    fn extract_batch_rejects_a_corrupt_record_count_without_reserving_it() {
        let records: Vec<Vec<u8>> =
            (0..3).map(|i| rec(300 + i, format!("c{i}a").as_bytes(), &[])).collect();
        let ctx = seed(None, KeyTypesSpec::Full).provision(&records[0]);
        let (arena, start, end) = arena_with(&records, 0);
        let Err(err) = ctx.extract_batch::<TemplateKey40>(&arena, start, end, u32::MAX - 1, 9)
        else {
            panic!("a corrupt record count must be rejected");
        };
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("fewer frames than n_records"), "{err}");
    }

    /// The reservation never exceeds the frames the span can hold, and is the
    /// claimed count when that fits.
    #[test]
    fn ref_reservation_is_capped_by_the_span() {
        let frame = 4 + MIN_BAM_RECORD_LEN;
        assert_eq!(ref_reservation(u32::MAX - 1, 3 * frame), 3);
        assert_eq!(ref_reservation(u32::MAX, 0), 0);
        assert_eq!(ref_reservation(2, 3 * frame), 2);
    }

    /// A frame shorter than the 32 fixed BAM bytes is rejected before the key
    /// extractor indexes its fixed fields; the message locates the frame.
    #[test]
    fn extract_batch_rejects_a_frame_shorter_than_a_bam_record() {
        let first = rec(1, b"f0a", &[]);
        let ctx = seed(None, KeyTypesSpec::Full).provision(&first);
        let (arena, start, end) = arena_with(&[vec![0u8; 12]], 0);
        let Err(err) = ctx.extract_batch::<TemplateKey40>(&arena, start, end, 1, 5) else {
            panic!("a 12-byte frame must be rejected");
        };
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        assert!(
            err.to_string().contains(&format!(
                "record 5 at byte {start}: block_size 12 is shorter than a BAM record"
            )),
            "{err}"
        );
    }

    /// The erased dispatch picks the arm matching the provisioned variant, for
    /// all four lane widths.
    #[rstest::rstest]
    #[case::k24(KeyTypesSpec::None, TemplateKeyVariant { cb: false, tertiary: false })]
    #[case::cb32(
        KeyTypesSpec::Explicit { cb: true, tertiary: false },
        TemplateKeyVariant { cb: true, tertiary: false }
    )]
    #[case::tert32(
        KeyTypesSpec::Explicit { cb: false, tertiary: true },
        TemplateKeyVariant { cb: false, tertiary: true }
    )]
    #[case::k40(KeyTypesSpec::Full, TemplateKeyVariant { cb: true, tertiary: true })]
    fn erased_dispatch_matches_the_variant(
        #[case] spec: KeyTypesSpec,
        #[case] want: TemplateKeyVariant,
    ) {
        let records: Vec<Vec<u8>> =
            (0..5).map(|i| rec(400 + i, format!("e{i}a").as_bytes(), &[])).collect();
        let ctx = seed(None, spec).provision(&records[0]);
        assert_eq!(ctx.variant, want);
        let (arena, start, end) = arena_with(&records, 0);
        let (refs, violation) = ctx.extract_batch_erased(&arena, start, end, 5, 0).unwrap();
        assert!(violation.is_none());
        assert_eq!(refs.variant(), want);
        assert_eq!(refs.len(), 5);
        assert_eq!(refs.byte_len(), 5 * want.ref_width());
        let arm_ok = matches!(
            (&refs, want.cb, want.tertiary),
            (ExtractedRefs::K24(_), false, false)
                | (ExtractedRefs::Cb32(_), true, false)
                | (ExtractedRefs::Tert32(_), false, true)
                | (ExtractedRefs::K40(_), true, true)
        );
        assert!(arm_ok);
    }

    /// The batched violation's error is byte-identical to the one the per-record
    /// accumulator returns for the same offending record (`TemplateStrategy`
    /// wraps that one with `{:#}`, as `into_io_error` does).
    #[test]
    fn violation_message_matches_the_serial_path() {
        let records =
            vec![rec(100, b"r0a", &cb_aux(b"AAAA")), rec(101, b"offender", &cb_aux(b"CCCC"))];
        let mut acc = crate::TemplateArenaAccumulator::from_header(
            &Header::default(),
            Some(SamTag::CB),
            KeyTypesSpec::None,
        );
        acc.push(&records[0], 4, u32::try_from(records[0].len()).unwrap()).unwrap();
        let serial = acc
            .push(&records[1], 200, u32::try_from(records[1].len()).unwrap())
            .expect_err("the accumulator rejects the differing CB");

        let ctx = seed(Some(SamTag::CB), KeyTypesSpec::None).provision(&records[0]);
        let (arena, start, end) = arena_with(&records, 0);
        let (_, violation) = ctx.extract_batch_erased(&arena, start, end, 2, 0).unwrap();
        let batched = violation.expect("the batch reports the differing CB").into_io_error();
        assert_eq!(batched.kind(), std::io::ErrorKind::InvalidData);
        assert_eq!(batched.to_string(), format!("{serial:#}"));
    }

    /// An empty batch (a run's empty terminal batch) is valid: no refs, no error.
    #[test]
    fn empty_batch_is_valid() {
        let first = rec(1, b"z0a", &[]);
        let ctx = seed(None, KeyTypesSpec::Auto).provision(&first);
        let (arena, start, _) = arena_with(&[first], 0);
        let b = ctx.extract_batch::<TemplateKey24>(&arena, start, start, 0, 9).unwrap();
        assert!(b.refs.is_empty() && b.violation.is_none());
    }
}
