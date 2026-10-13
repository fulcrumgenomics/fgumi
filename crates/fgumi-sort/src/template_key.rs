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
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;
use std::collections::{HashMap, HashSet};

use crate::arena_pool::PooledSegmentedBuf;
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

    /// The lanes `--key-types` fixes before any record is seen; `None` under
    /// `auto`, where the first record chooses them.
    #[must_use]
    pub fn known_variant(&self) -> Option<TemplateKeyVariant> {
        (self.key_types != KeyTypesSpec::Auto)
            .then(|| select_template_variant(None, self.key_types, self.header_library_varies))
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
        self.provision_inner(first_body)
    }

    /// The full-width context a merge of independently sorted inputs extracts
    /// its keys with (`fgumi merge`): every lane, CB included, whatever this
    /// seed's `--key-types`. No first record is needed: a full key drops no
    /// lane, so there is no dropped-lane baseline to check against.
    #[must_use]
    pub fn merge_context(&self) -> Arc<TemplateKeyContext> {
        Arc::new(TemplateKeyContext {
            lib_lookup: self.lib_lookup.clone(),
            cell_tag: self.cell_tag,
            cb_hasher: self.cb_hasher.clone(),
            first_key: TemplateKey::default(),
            variant: select_template_variant(None, KeyTypesSpec::Full, self.header_library_varies),
        })
    }

    fn provision_inner(&self, first_body: &[u8]) -> Arc<TemplateKeyContext> {
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
    /// against `first_key`. The per-record accumulator, the batched path and a
    /// body-keyed merge slot (`fgumi merge`) call this, so they cannot drift.
    ///
    /// # Panics
    ///
    /// Panics if `body` is shorter than [`MIN_BAM_RECORD_LEN`], as
    /// [`TemplateKeySeed::provision`] does; every caller checks first (the
    /// merge slot parser checks the whole record layout).
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
        // Every `first_record_index + i` below is at most the range's last
        // index. Byte offsets need no such check: `byte_start + cur` never
        // exceeds `byte_end`, because `cur` never exceeds `span.len()`.
        if first_record_index.checked_add(u64::from(n_records.saturating_sub(1))).is_none() {
            return Err(malformed(first_record_index, "record index range overflows u64"));
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

// ============================================================================
// Library Lookup for Template-Coordinate Sort
// ============================================================================

/// Deterministic hasher for cell barcode hashing in template-coordinate sort.
///
/// Delegates to the shared [`fgumi_dna::deterministic_state`] so the fixed seed table lives in
/// exactly one place; hash values remain reproducible across runs.
#[must_use]
pub fn cb_hasher() -> ahash::RandomState {
    fgumi_dna::deterministic_state()
}

/// Maps read group ID -> library ordinal for O(1) comparison.
///
/// Pre-computes ordinals by sorting library names alphabetically.
/// Empty/unknown library sorts first (ordinal 0).
#[derive(Clone)]
pub struct LibraryLookup {
    /// RG ID -> library ordinal
    rg_to_ordinal: HashMap<Vec<u8>, u32>,
    /// Deterministic hasher for read name hashing, constructed once for reuse.
    hasher: ahash::RandomState,
}

impl LibraryLookup {
    /// Build lookup from BAM header.
    #[must_use]
    #[allow(clippy::cast_possible_truncation)]
    pub fn from_header(header: &Header) -> Self {
        // Collect all unique library names from read groups
        let mut libraries: Vec<String> = header
            .read_groups()
            .iter()
            .filter_map(|(_, rg)| {
                rg.other_fields().get(&rg_tag::LIBRARY).map(std::string::ToString::to_string)
            })
            .collect();

        // Sort alphabetically and deduplicate
        libraries.sort();
        libraries.dedup();

        // Build library name -> ordinal mapping
        // Empty string gets ordinal 0, then libraries in sorted order
        let mut lib_to_ordinal: HashMap<String, u32> = HashMap::new();
        lib_to_ordinal.insert(String::new(), 0);
        for (i, lib) in libraries.iter().enumerate() {
            lib_to_ordinal.insert(lib.clone(), (i + 1) as u32);
        }

        // Build RG ID -> ordinal mapping
        let rg_to_ordinal: HashMap<Vec<u8>, u32> = header
            .read_groups()
            .iter()
            .map(|(id, rg)| {
                let lib = rg
                    .other_fields()
                    .get(&rg_tag::LIBRARY)
                    .map(std::string::ToString::to_string)
                    .unwrap_or_default();
                let ordinal = *lib_to_ordinal.get(&lib).unwrap_or(&0);
                (id.to_vec(), ordinal)
            })
            .collect();

        // Arbitrary fixed seeds — chosen for uniqueness, not cryptographic strength.
        let hasher = ahash::RandomState::with_seeds(
            0x517c_c1b7_2722_0a95,
            0x1234_5678_90ab_cdef,
            0xfedc_ba98_7654_3210,
            0x0123_4567_89ab_cdef,
        );

        Self { rg_to_ordinal, hasher }
    }

    /// Hash a read name deterministically.
    #[inline]
    #[must_use]
    pub fn hash_name(&self, name: &[u8]) -> u64 {
        self.hasher.hash_one(name)
    }

    /// Get library ordinal for a record (from RG tag in aux data).
    ///
    /// Only used in tests; production code uses [`ordinal_from_rg`](Self::ordinal_from_rg)
    /// with pre-extracted RG bytes from the single-pass aux scan.
    #[cfg(test)]
    #[must_use]
    pub fn get_ordinal(&self, bam: &[u8]) -> u32 {
        fgumi_raw_bam::RawRecordView::new(bam)
            .tags()
            .find_string(SamTag::RG)
            .and_then(|rg| self.rg_to_ordinal.get(rg))
            .copied()
            .unwrap_or(0)
    }

    /// Get library ordinal from pre-extracted RG tag bytes.
    #[inline]
    #[must_use]
    pub fn ordinal_from_rg(&self, rg: Option<&[u8]>) -> u32 {
        rg.and_then(|rg| self.rg_to_ordinal.get(rg)).copied().unwrap_or(0)
    }

    /// Number of distinct library ordinals realized by the header's read groups.
    ///
    /// Used to decide whether the template-key `tertiary` library lane can vary:
    /// `<= 1` means every header read group maps to the same library ordinal, so
    /// the lane is provisionally constant (decode-time verify still backstops the
    /// case where a record carries an RG absent from the header → ordinal 0).
    /// Returns 0 when the header declares no read groups.
    #[must_use]
    pub fn distinct_header_ordinals(&self) -> usize {
        self.rg_to_ordinal.values().copied().collect::<HashSet<u32>>().len()
    }
}

/// Extract a packed `TemplateKey` directly from BAM record bytes.
///
/// This function computes the template-coordinate sort key inline, avoiding
/// heap allocations for the read name by using a hash instead.
///
/// When `cell_tag` is `Some`, the CB (cellular barcode) tag value is hashed
/// and included in the sort key between neg2 and MI, matching fgbio's order.
#[must_use]
pub fn extract_template_key_inline(
    bam_bytes: &[u8],
    lib_lookup: &LibraryLookup,
    cell_tag: Option<SamTag>,
    cb_hasher: &ahash::RandomState,
) -> TemplateKey {
    use fgumi_raw_bam;
    use fgumi_raw_bam::{flags, mate_unclipped_5prime, unclipped_5prime_raw};

    // Single-pass extraction of all aux tags (MI, RG, cell barcode, MC)
    let aux = fgumi_raw_bam::extract_template_aux_tags(bam_bytes, cell_tag);
    let mi = aux.mi;
    let library = lib_lookup.ordinal_from_rg(aux.rg);
    let cb_hash = aux.cell.map_or(0u64, |cb_bytes| cb_hasher.hash_one(cb_bytes));

    // Extract fields from raw bytes
    let v = fgumi_raw_bam::RawRecordView::new(bam_bytes);
    let tid = v.ref_id();
    let pos = v.pos();
    let l_read_name = v.l_read_name() as usize;
    let flag = v.flags();
    let mate_tid = v.mate_ref_id();
    let mate_pos = v.mate_pos();

    // Extract flags
    let is_unmapped = (flag & flags::UNMAPPED) != 0;
    let mate_unmapped = (flag & flags::MATE_UNMAPPED) != 0;
    let is_reverse = (flag & flags::REVERSE) != 0;
    let mate_reverse = (flag & flags::MATE_REVERSE) != 0;
    let is_paired = (flag & flags::PAIRED) != 0;

    // Hash read name (exclude null terminator)
    let name_len = l_read_name.saturating_sub(1);
    let name = if name_len > 0 && 32 + name_len <= bam_bytes.len() {
        &bam_bytes[32..32 + name_len]
    } else {
        &[]
    };
    let name_hash = lib_lookup.hash_name(name);

    // Secondary/supplementary reads cannot reconstruct their template coordinate
    // from their own record (they carry their own and their mate's position, but
    // not their own primary's). When `fgumi zipper` has stamped the exact
    // template coordinate into the `tc` tag, key on it so the read sorts at its
    // primary pair's coordinate — a placement per-record keying (samtools/fgbio)
    // cannot achieve. Falls through to the position-based computation when absent.
    let is_secondary = (flag & flags::SECONDARY) != 0;
    let is_supplementary = (flag & flags::SUPPLEMENTARY) != 0;
    if is_secondary || is_supplementary {
        let aux_slice = fgumi_raw_bam::aux_data_slice(bam_bytes);
        if let Some([tid1, pos1, neg1, tid2, pos2, neg2]) =
            fgumi_raw_bam::read_tc_template_coordinate(aux_slice)
        {
            // `tc` is canonicalized by (tid, pos) only (see zipper's
            // `add_template_coordinate_tags_raw`), whereas the per-record path
            // additionally tiebreaks on strand for equal positions. In the
            // degenerate case where both primaries share an identical tid and
            // unclipped-5' position but differ in strand, the supplementary's
            // secondary-lane strand order can differ from its primaries'. This is
            // harmless: the primary lane (tid1, tid2, pos1) still matches, so the
            // read co-locates with its template; only a same-coordinate tiebreak
            // differs.
            //
            // `TemplateKey::new`'s final argument is the `is_upper` bit, packed
            // into `name_hash_upper` as the *last* tiebreak (after position,
            // strand, CB, library, MI, and name). For a mapped read it answers
            // "is this read's own primary the upper (lane-2) end?", derived from
            // this read's position vs. its mate. A secondary/supplementary read
            // cannot recover that: `tc` is canonicalized by (tid, pos) and does
            // not record which segment (R1/R2) landed on which lane, and the
            // read's own alignment position is exactly the misleading coordinate
            // this branch exists to avoid. We use `is_read2` because it mirrors
            // the primary's `is_upper` for the standard FR pair (R1 is the lower
            // lane -> is_upper=false; R2 is the upper lane -> is_upper=true) and
            // is fully deterministic; see
            // `test_extract_template_key_secondary_supplementary_matches_primary_via_tc`,
            // which asserts the tc-keyed key is byte-identical to the primary's.
            // Any residual disagreement (an exotic pair whose R2 is the lower
            // lane) only reorders same-name records at an identical coordinate,
            // so it cannot affect co-location or dedup grouping.
            let is_read2 = (flag & 0x80) != 0;
            return TemplateKey::new(
                tid1,
                pos1,
                neg1 != 0,
                tid2,
                pos2,
                neg2 != 0,
                cb_hash,
                library,
                mi,
                name_hash,
                is_read2,
            );
        }
    }

    // Handle unmapped reads
    if is_unmapped {
        if is_paired && !mate_unmapped {
            // Unmapped read with mapped mate - use mate's position as primary key
            let mate_unclipped =
                aux.mc.map_or(mate_pos, |mc| mate_unclipped_5prime(mate_pos, mate_reverse, mc));

            return TemplateKey::new(
                mate_tid,
                mate_unclipped,
                mate_reverse,
                i32::MAX,
                i32::MAX,
                false,
                cb_hash,
                library,
                mi,
                name_hash,
                true, // Unmapped read is always "upper" relative to mapped mate
            );
        }

        // Completely unmapped - sort to end. Still carry the read's library and
        // MI lanes so a fully-unmapped read realizes the same tertiary as its
        // mapped, same-library peers (otherwise the dropped-lane verify treats a
        // single-library file as if the library varied; see #375).
        let is_read2 = (flag & 0x80) != 0; // is_last_segment flag
        return TemplateKey::unmapped(name_hash, cb_hash, library, mi, is_read2);
    }

    // Calculate unclipped 5' position for this read (zero-alloc: reads cigar directly)
    let this_pos = unclipped_5prime_raw(bam_bytes, pos, is_reverse);

    // Calculate mate's unclipped 5' position
    let mate_unclipped = if is_paired && !mate_unmapped {
        aux.mc.map_or(mate_pos, |mc| mate_unclipped_5prime(mate_pos, mate_reverse, mc))
    } else {
        mate_pos
    };

    // Determine canonical ordering
    let (tid1, tid2, pos1, pos2, neg1, neg2, is_upper) = if is_paired && !mate_unmapped {
        // Samtools logic: is_upper if pos > mate_pos, or (pos == mate_pos && this read is reverse)
        let is_upper = (tid, this_pos) > (mate_tid, mate_unclipped)
            || ((tid, this_pos) == (mate_tid, mate_unclipped) && is_reverse);

        if is_upper {
            // Swap: mate's position comes first
            (mate_tid, tid, mate_unclipped, this_pos, mate_reverse, is_reverse, true)
        } else {
            // No swap: this read's position comes first
            (tid, mate_tid, this_pos, mate_unclipped, is_reverse, mate_reverse, false)
        }
    } else {
        // Unpaired or mate unmapped - use MAX for tid2/pos2
        (tid, i32::MAX, this_pos, i32::MAX, is_reverse, false, false)
    };

    TemplateKey::new(tid1, pos1, neg1, tid2, pos2, neg2, cb_hash, library, mi, name_hash, is_upper)
}

// ============================================================================
// Template-key variant selection
// ============================================================================

/// Which optional lanes a chosen template key retains.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct TemplateKeyVariant {
    /// Retain the `cb_hash` lane.
    pub cb: bool,
    /// Retain the tertiary (library|mi) lane.
    pub tertiary: bool,
}

impl TemplateKeyVariant {
    /// The widest variant (every optional lane): what a template ref is charged
    /// at before the sort's first record chooses its lanes.
    pub const WIDEST: Self = Self { cb: true, tertiary: true };

    /// Number of u64 lanes: 3 + cb + tertiary.
    #[must_use]
    pub fn lanes(self) -> usize {
        3 + usize::from(self.cb) + usize::from(self.tertiary)
    }

    /// Bytes one accumulated `TemplateRecordRef<K>` occupies for this variant:
    /// the narrowed key (8 bytes per lane) plus the 16-byte
    /// `offset`/`len`/`padding` tail. Derived from [`Self::lanes`], which owns the
    /// lane count; pinned against `size_of` by `ref_width_matches_the_ref_size`.
    /// A logical size: records × `ref_width` is what the refs hold, while a
    /// growing `Vec` of them can have up to about twice that allocated.
    #[must_use]
    pub fn ref_width(self) -> usize {
        8 * self.lanes() + 16
    }
}

/// Parsed `--key-types` spec controlling which optional sort-key lanes are kept.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub enum KeyTypesSpec {
    /// Auto-detect from the first record (default).
    #[default]
    Auto,
    /// Force all lanes (N=5).
    Full,
    /// Force no optional lanes (N=3).
    None,
    /// Force a specific optional set.
    Explicit {
        /// Retain the `cb_hash` lane.
        cb: bool,
        /// Retain the tertiary (library|mi) lane.
        tertiary: bool,
    },
}

/// Mask for the MI component of the `tertiary` key lane (low 48 bits).
/// `tertiary` packs `library_ordinal << 48 | (mi_value << 1 | !mi_suffix)`, so the
/// high 16 bits are the library ordinal and the low 48 bits encode MI. A nonzero
/// low-48 region means the first record carries an MI tag.
const TERTIARY_MI_MASK: u64 = (1u64 << 48) - 1;

/// Choose the template-key variant for this sort.
///
/// `first_key` is the full key extracted from the first record (`None` for an
/// empty input — then the variant is irrelevant; default to lite). `spec` is the
/// parsed `--key-types` override. `header_library_varies` is `true` when the
/// header realizes more than one distinct library ordinal (see
/// [`LibraryLookup::distinct_header_ordinals`]).
///
/// Selection only *provisions* the variant; the decode-time verify (fill loop)
/// still *guarantees* the dropped lanes are constant. Under `Auto`, the tertiary
/// library lane is kept only when it can actually vary: the first record carries
/// an MI tag, or the header realizes more than one library ordinal. A single,
/// constant library ordinal in the high bits is droppable — the header informs
/// provisioning, while verify backstops the rare case of a record whose RG is
/// absent from the header (which realizes ordinal 0).
#[must_use]
pub fn select_template_variant(
    first_key: Option<&TemplateKey>,
    spec: KeyTypesSpec,
    header_library_varies: bool,
) -> TemplateKeyVariant {
    match spec {
        KeyTypesSpec::Full => TemplateKeyVariant { cb: true, tertiary: true },
        KeyTypesSpec::None => TemplateKeyVariant { cb: false, tertiary: false },
        KeyTypesSpec::Explicit { cb, tertiary } => TemplateKeyVariant { cb, tertiary },
        KeyTypesSpec::Auto => match first_key {
            None => TemplateKeyVariant { cb: false, tertiary: false },
            Some(k) => {
                // Keep the tertiary lane only if it can actually vary: MI present on
                // the first record (low-48 bits), or the header realizes more than
                // one library ordinal. A single (constant) library ordinal in the
                // high bits is droppable — decode-time verify backstops the rare
                // case of a record whose RG is absent from the header (ordinal 0).
                let mi_present = (k.tertiary & TERTIARY_MI_MASK) != 0;
                TemplateKeyVariant {
                    cb: k.cb_hash != 0,
                    tertiary: mi_present || header_library_varies,
                }
            }
        },
    }
}

// ============================================================================
// Dropped-lane violation detection
// ============================================================================

/// A dropped-lane violation discovered during decode-time verify.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DroppedLaneViolation {
    /// A record carried a CB tag absent from the first record.
    Cb,
    /// A record carried an MI value differing from the first record.
    Mi,
    /// A record's realized library ordinal differed from the first record.
    Library,
}

impl DroppedLaneViolation {
    /// The `--key-types` token that re-includes this lane.
    #[must_use]
    pub fn key_types_token(self) -> &'static str {
        match self {
            DroppedLaneViolation::Cb => "cb",
            DroppedLaneViolation::Mi => "mi",
            DroppedLaneViolation::Library => "library",
        }
    }
}

/// Verify that the lanes the chosen `variant` dropped are constant for `cur`
/// relative to the sort's `first` full key. Returns the first violation found,
/// or `None` if all dropped lanes match.
///
/// Tertiary sub-fields are decoded only to attribute the violation: the high 16
/// bits are the library ordinal; the low 48 bits are `(mi_value << 1) | !mi_suffix`
/// (`mi_value` occupies bits 47..1; `!mi_suffix` is bit 0).
#[must_use]
pub fn verify_dropped_lanes(
    first: &TemplateKey,
    cur: &TemplateKey,
    variant: TemplateKeyVariant,
) -> Option<DroppedLaneViolation> {
    if !variant.cb && cur.cb_hash != first.cb_hash {
        return Some(DroppedLaneViolation::Cb);
    }
    if !variant.tertiary && cur.tertiary != first.tertiary {
        // Attribute: library is the high 16 bits, mi the low 48.
        let lib_changed = (cur.tertiary >> 48) != (first.tertiary >> 48);
        return Some(if lib_changed {
            DroppedLaneViolation::Library
        } else {
            DroppedLaneViolation::Mi
        });
    }
    None
}

/// The read name of a record that violated a dropped lane, for its error
/// message: the bytes [`RawRecordView::read_name`](fgumi_raw_bam::RawRecordView::read_name)
/// returns when `l_read_name` fits inside `body`, else `<malformed read name>`.
/// The key extractor tolerates a name that overruns its body (it bounds-checks
/// the name), so such a record can reach the violation branch, which must not
/// slice past the body while formatting the error.
pub(crate) fn violation_read_name(body: &[u8]) -> String {
    let l = body.get(8).map_or(0, |&b| usize::from(b));
    let name: &[u8] = if l > 1 {
        match body.get(32..32 + l - 1) {
            Some(name) => name,
            None => return "<malformed read name>".to_owned(),
        }
    } else {
        &[]
    };
    String::from_utf8_lossy(name).into_owned()
}

/// Build the actionable error message for a dropped-lane violation.
#[must_use]
pub(crate) fn dropped_lane_error(name: &str, v: DroppedLaneViolation) -> anyhow::Error {
    let field = match v {
        DroppedLaneViolation::Cb => "CB",
        DroppedLaneViolation::Mi => "MI",
        DroppedLaneViolation::Library => "library",
    };
    anyhow::anyhow!(
        "record {name} carries a {field} value absent from the input's first record; \
         re-run with --key-types {}",
        v.key_types_token(),
    )
}

#[cfg(test)]
mod tests {
    use fgumi_raw_bam::SamTag;
    use fgumi_raw_bam::testutil::make_bam_bytes;
    use noodles::sam::Header;

    use super::*;
    use crate::arena_pool::PooledSegmentedBuf;
    use crate::inline::{TemplateKey24, TemplateKey40};
    use crate::segmented_buf::SegmentedBuf;
    use crate::template_key::{KeyTypesSpec, extract_template_key_inline};

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

    /// A corrupt `first_record_index` whose record range overflows `u64` is an
    /// `InvalidData` error before any frame is read, not a debug-build panic
    /// or a wrapped index in a release build's error messages.
    #[test]
    fn extract_batch_rejects_a_record_index_range_that_overflows() {
        let records: Vec<Vec<u8>> =
            (0..3).map(|i| rec(300 + i, format!("o{i}a").as_bytes(), &[])).collect();
        let ctx = seed(None, KeyTypesSpec::Full).provision(&records[0]);
        let (arena, start, end) = arena_with(&records, 0);
        let Err(err) = ctx.extract_batch::<TemplateKey40>(&arena, start, end, 3, u64::MAX - 1)
        else {
            panic!("an overflowing record index range must be rejected");
        };
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("record index range overflows u64"), "{err}");
        // The last index of the range is exactly `u64::MAX`: still valid.
        assert!(ctx.extract_batch::<TemplateKey40>(&arena, start, end, 3, u64::MAX - 2).is_ok());
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

#[cfg(test)]
mod moved_key_tests {
    use super::*;
    use bstr::BString;
    use fgumi_raw_bam::flags;
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReadGroup;
    use rstest::rstest;

    // ========================================================================
    // LibraryLookup tests
    // ========================================================================

    #[test]
    fn test_library_lookup_empty_header() {
        let header = Header::builder().build();
        let lookup = LibraryLookup::from_header(&header);
        assert!(lookup.rg_to_ordinal.is_empty());
    }

    #[test]
    fn test_library_lookup_single_rg() {
        let rg = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let header = Header::builder().add_read_group(BString::from("rg1"), rg).build();

        let lookup = LibraryLookup::from_header(&header);
        assert_eq!(lookup.rg_to_ordinal.len(), 1);
        // LibA is the only library, so it gets ordinal 1 (0 is reserved for empty/unknown)
        assert_eq!(
            *lookup.rg_to_ordinal.get(b"rg1".as_slice()).expect("rg1 should be in ordinal map"),
            1
        );
    }

    #[test]
    fn test_library_lookup_multiple_libraries() {
        let rg_a = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibC"))
            .build()
            .expect("valid");
        let rg_b = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let rg_c = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibB"))
            .build()
            .expect("valid");

        let header = Header::builder()
            .add_read_group(BString::from("rg1"), rg_a)
            .add_read_group(BString::from("rg2"), rg_b)
            .add_read_group(BString::from("rg3"), rg_c)
            .build();

        let lookup = LibraryLookup::from_header(&header);
        assert_eq!(lookup.rg_to_ordinal.len(), 3);

        // Libraries sorted alphabetically: LibA=1, LibB=2, LibC=3
        let rg2 = *lookup.rg_to_ordinal.get(b"rg2".as_slice()).expect("rg2");
        let rg3 = *lookup.rg_to_ordinal.get(b"rg3".as_slice()).expect("rg3");
        let rg1 = *lookup.rg_to_ordinal.get(b"rg1".as_slice()).expect("rg1");
        assert_eq!(rg2, 1); // LibA
        assert_eq!(rg3, 2); // LibB
        assert_eq!(rg1, 3); // LibC
    }

    #[test]
    fn test_distinct_header_ordinals() {
        // (i) no read groups -> 0.
        let empty = LibraryLookup::from_header(&Header::builder().build());
        assert_eq!(empty.distinct_header_ordinals(), 0);

        // (ii) one RG with an LB -> 1.
        let rg_a = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let one = LibraryLookup::from_header(
            &Header::builder().add_read_group(BString::from("rg1"), rg_a).build(),
        );
        assert_eq!(one.distinct_header_ordinals(), 1);

        // (iii) two RGs with two different LBs -> 2.
        let rg_b = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let rg_c = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibB"))
            .build()
            .expect("valid");
        let two = LibraryLookup::from_header(
            &Header::builder()
                .add_read_group(BString::from("rg1"), rg_b)
                .add_read_group(BString::from("rg2"), rg_c)
                .build(),
        );
        assert_eq!(two.distinct_header_ordinals(), 2);

        // (iv) two RGs with the SAME LB -> 1.
        let rg_d = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let rg_e = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let same = LibraryLookup::from_header(
            &Header::builder()
                .add_read_group(BString::from("rg1"), rg_d)
                .add_read_group(BString::from("rg2"), rg_e)
                .build(),
        );
        assert_eq!(same.distinct_header_ordinals(), 1);

        // (v) one RG with an LB + one RG WITHOUT an LB -> 2 (ordinals {1, 0}).
        let rg_f = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let rg_no_lb = Map::<ReadGroup>::builder().build().expect("valid");
        let mixed = LibraryLookup::from_header(
            &Header::builder()
                .add_read_group(BString::from("rg1"), rg_f)
                .add_read_group(BString::from("rg2"), rg_no_lb)
                .build(),
        );
        assert_eq!(mixed.distinct_header_ordinals(), 2);
    }

    #[test]
    fn test_library_lookup_unknown_rg_returns_zero() {
        let rg = Map::<ReadGroup>::builder()
            .insert(rg_tag::LIBRARY, String::from("LibA"))
            .build()
            .expect("valid");
        let header = Header::builder().add_read_group(BString::from("rg1"), rg).build();

        let lookup = LibraryLookup::from_header(&header);
        // A BAM record with no RG tag or unknown RG should return ordinal 0
        // We can test get_ordinal with a minimal BAM record that has no aux data
        let mut bam = vec![0u8; 36];
        bam[8] = 4; // l_read_name = 4
        bam[32..36].copy_from_slice(b"rea\0");
        assert_eq!(lookup.get_ordinal(&bam), 0);
    }

    // ========================================================================
    // find_string_tag_in_record tests (via fgumi_raw_bam)
    // ========================================================================

    /// Helper to build minimal BAM bytes with aux data appended.
    /// Fixed 32-byte header + read name "rea\0" (`l_read_name=4`) + no cigar + no seq + aux.
    fn build_bam_with_aux(aux_data: &[u8]) -> Vec<u8> {
        let l_read_name: u8 = 4; // "rea" + null
        let mut bam = vec![0u8; 32];
        bam[8] = l_read_name;
        // n_cigar_op = 0 (bytes 12-13 already zero)
        // l_seq = 0 (bytes 16-19 already zero)
        // read name
        bam.extend_from_slice(b"rea\0");
        // no cigar, no seq, no qual
        // append aux data
        bam.extend_from_slice(aux_data);
        bam
    }

    #[rstest::rstest]
    #[case::present(b"RGZgroup1\0".as_slice(), Some(b"group1".as_slice()))]
    #[case::absent(b"".as_slice(), None)]
    fn test_find_rg_tag(#[case] aux_data: &[u8], #[case] expected: Option<&[u8]>) {
        let bam = build_bam_with_aux(aux_data);
        assert_eq!(
            fgumi_raw_bam::RawRecordView::new(&bam).tags().find_string(SamTag::RG),
            expected
        );
    }

    #[test]
    fn test_find_rg_tag_after_other_tags() {
        // XY:i:42 followed by RG:Z:mygroup — aux built dynamically for the integer tag
        let mut aux = Vec::new();
        aux.extend_from_slice(b"XYi");
        aux.extend_from_slice(&42i32.to_le_bytes());
        aux.extend_from_slice(b"RGZmygroup\0");
        let bam = build_bam_with_aux(&aux);
        assert_eq!(
            fgumi_raw_bam::RawRecordView::new(&bam).tags().find_string(SamTag::RG),
            Some(b"mygroup".as_slice())
        );
    }

    // ========================================================================
    // extract_template_key_inline cell_tag tests
    // ========================================================================

    fn test_cb_hasher() -> ahash::RandomState {
        cb_hasher()
    }

    /// Build minimal BAM bytes with mapped read at (tid, pos) on the forward strand,
    /// with optional aux data appended.
    #[allow(clippy::cast_possible_truncation)]
    fn build_mapped_bam(tid: i32, pos: i32, name: &[u8], aux: &[u8]) -> Vec<u8> {
        let l_read_name = (name.len() + 1) as u8; // +1 for null terminator
        let mut bam = vec![0u8; 32];
        // ref_id at offset 0..4
        bam[0..4].copy_from_slice(&tid.to_le_bytes());
        // pos at offset 4..8
        bam[4..8].copy_from_slice(&pos.to_le_bytes());
        // l_read_name at offset 8
        bam[8] = l_read_name;
        // flags at offset 14..16: paired + proper pair = 0x03
        bam[14..16].copy_from_slice(&3u16.to_le_bytes());
        // mate_ref_id at offset 20..24: same tid
        bam[20..24].copy_from_slice(&tid.to_le_bytes());
        // mate_pos at offset 24..28: same pos
        bam[24..28].copy_from_slice(&pos.to_le_bytes());
        // read name
        bam.extend_from_slice(name);
        bam.push(0); // null terminator
        // no cigar, no seq, no qual
        // aux data
        bam.extend_from_slice(aux);
        bam
    }

    /// Build CB:Z: aux tag bytes.
    fn cb_aux(value: &[u8]) -> Vec<u8> {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"CBZ");
        aux.extend_from_slice(value);
        aux.push(0); // null terminator
        aux
    }

    #[test]
    fn test_extract_template_key_cb_present_has_nonzero_hash() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        let aux = cb_aux(b"ACGTACGT");
        let bam = build_mapped_bam(0, 100, b"read1", &aux);

        let key =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        assert_ne!(key.cb_hash, 0, "CB present should produce non-zero cb_hash");
    }

    #[test]
    fn test_extract_template_key_cb_absent_has_zero_hash() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        // No CB tag in aux data
        let bam = build_mapped_bam(0, 100, b"read1", &[]);

        let key =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        assert_eq!(key.cb_hash, 0, "missing CB tag should produce cb_hash=0");
    }

    #[test]
    fn test_extract_template_key_cell_tag_none_has_zero_hash() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        let aux = cb_aux(b"ACGTACGT");
        let bam = build_mapped_bam(0, 100, b"read1", &aux);

        let key = extract_template_key_inline(&bam, &lib_lookup, None, &test_cb_hasher());
        assert_eq!(key.cb_hash, 0, "cell_tag=None should produce cb_hash=0");
    }

    #[test]
    fn test_extract_template_key_different_cb_values_differ() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);

        let aux1 = cb_aux(b"ACGTACGT");
        let bam1 = build_mapped_bam(0, 100, b"read1", &aux1);
        let key1 =
            extract_template_key_inline(&bam1, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());

        let aux2 = cb_aux(b"TGCATGCA");
        let bam2 = build_mapped_bam(0, 100, b"read1", &aux2);
        let key2 =
            extract_template_key_inline(&bam2, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());

        assert_ne!(
            key1.cb_hash, key2.cb_hash,
            "different CB values should produce different hashes"
        );
    }

    #[test]
    fn test_extract_template_key_cb_hash_is_deterministic() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        let aux = cb_aux(b"ACGTACGT");
        let bam = build_mapped_bam(0, 100, b"read1", &aux);

        let key1 =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        let key2 =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        assert_eq!(key1.cb_hash, key2.cb_hash, "same input should produce same cb_hash");
    }

    /// Encode a `tc:B:i` aux tag carrying the 6-element template-coordinate
    /// array `[tid1, pos1, neg1, tid2, pos2, neg2]` that `fgumi zipper` stamps
    /// onto secondary/supplementary reads.
    fn tc_aux(vals: &[i32; 6]) -> Vec<u8> {
        let mut aux = Vec::new();
        aux.extend_from_slice(b"tcB"); // tag "tc", array type 'B'
        aux.push(b'i'); // subtype: int32
        #[allow(clippy::cast_possible_truncation)]
        aux.extend_from_slice(&(vals.len() as u32).to_le_bytes());
        for v in vals {
            aux.extend_from_slice(&v.to_le_bytes());
        }
        aux
    }

    /// Build minimal BAM bytes for a mapped secondary/supplementary alignment:
    /// paired, mate mapped, at `own_pos`, with its mate at `mate_pos`, carrying the
    /// caller-supplied `sec_supp_flag` (`flags::SECONDARY` or `flags::SUPPLEMENTARY`),
    /// plus optional aux. Both flags route through the same `tc` keying branch.
    #[allow(clippy::cast_possible_truncation)]
    fn build_secondary_or_supplementary_bam(
        own_pos: i32,
        mate_pos: i32,
        sec_supp_flag: u16,
        name: &[u8],
        aux: &[u8],
    ) -> Vec<u8> {
        let l_read_name = (name.len() + 1) as u8;
        let mut bam = vec![0u8; 32];
        bam[0..4].copy_from_slice(&0i32.to_le_bytes()); // ref_id = chr0
        bam[4..8].copy_from_slice(&own_pos.to_le_bytes());
        bam[8] = l_read_name;
        // flags: PAIRED | (SECONDARY or SUPPLEMENTARY) -- mate is mapped.
        bam[14..16].copy_from_slice(&(flags::PAIRED | sec_supp_flag).to_le_bytes());
        bam[20..24].copy_from_slice(&0i32.to_le_bytes()); // mate ref_id = chr0
        bam[24..28].copy_from_slice(&mate_pos.to_le_bytes());
        bam.extend_from_slice(name);
        bam.push(0);
        bam.extend_from_slice(aux);
        bam
    }

    /// Build minimal BAM bytes for a mapped primary read on `chr0` with the
    /// caller-supplied `flag`, aligned at `own_pos` with its mate at `mate_pos`.
    /// No CIGAR, so the unclipped 5' position equals `own_pos`/`mate_pos`. Used to
    /// key the two primaries of a pair through the per-record (mapped) path so
    /// their `TemplateKey`s can be compared against the tc-keyed sec/supp keys.
    #[allow(clippy::cast_possible_truncation)]
    fn build_primary_bam(own_pos: i32, mate_pos: i32, flag: u16, name: &[u8]) -> Vec<u8> {
        let l_read_name = (name.len() + 1) as u8;
        let mut bam = vec![0u8; 32];
        bam[0..4].copy_from_slice(&0i32.to_le_bytes()); // ref_id = chr0
        bam[4..8].copy_from_slice(&own_pos.to_le_bytes());
        bam[8] = l_read_name;
        bam[14..16].copy_from_slice(&flag.to_le_bytes());
        bam[20..24].copy_from_slice(&0i32.to_le_bytes()); // mate ref_id = chr0
        bam[24..28].copy_from_slice(&mate_pos.to_le_bytes());
        bam.extend_from_slice(name);
        bam.push(0);
        bam
    }

    /// A secondary/supplementary alignment carries its own position and its mate's
    /// position, but **not** its own primary's position. Per-record keying (samtools
    /// / fgbio) can therefore only place it by the mate coordinate, which can drift a
    /// fragment-length from the template's true coordinate. Because `fgumi zipper`
    /// pre-computes the exact template coordinate into the `tc` tag, fgumi's
    /// template-coordinate sort can place the read **exactly** at its primary pair's
    /// coordinate — a placement samtools/fgbio cannot achieve. This is the novel
    /// fgumi capability; both SECONDARY and SUPPLEMENTARY take the tc branch.
    #[rstest]
    #[case::secondary(flags::SECONDARY)]
    #[case::supplementary(flags::SUPPLEMENTARY)]
    fn test_extract_template_key_secondary_supplementary_uses_tc_tag(#[case] sec_supp_flag: u16) {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);

        // Sec/supp of a left-anchored R1: the primary is at pos 1000, its mate (R2
        // primary) at 1400, but this alignment is at 9000. The exact template
        // coordinate (pos1) is 1000, carried in `tc`.
        let tc = tc_aux(&[0, 1000, 0, 0, 1400, 1]);
        let supp = build_secondary_or_supplementary_bam(9000, 1400, sec_supp_flag, b"tmpl", &tc);
        let supp_key = extract_template_key_inline(&supp, &lib_lookup, None, &test_cb_hasher());

        // A different template whose true coordinate (pos1) is 1200 — between the
        // read's exact coordinate (1000, from tc) and the mate-only approximation
        // (1400).
        let filler = build_mapped_bam(0, 1200, b"fillr", &[]);
        let filler_key = extract_template_key_inline(&filler, &lib_lookup, None, &test_cb_hasher());

        assert!(
            supp_key < filler_key,
            "sec/supp must sort at its template coordinate (pos1=1000 via tc), before \
             the pos1=1200 template; per-record keying misplaces it at the mate \
             coordinate (1400), after the 1200 template"
        );
    }

    /// Without a `tc` tag, a secondary/supplementary must fall back to the per-record
    /// (own-position / mate) computation — the prior behavior — so the tc branch
    /// never changes keying for reads zipper did not stamp.
    #[rstest]
    #[case::secondary(flags::SECONDARY)]
    #[case::supplementary(flags::SUPPLEMENTARY)]
    fn test_extract_template_key_secondary_supplementary_without_tc_falls_back(
        #[case] sec_supp_flag: u16,
    ) {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);

        // Same geometry as the tc test, but no `tc` tag.
        let supp = build_secondary_or_supplementary_bam(9000, 1400, sec_supp_flag, b"tmpl", &[]);
        let supp_key = extract_template_key_inline(&supp, &lib_lookup, None, &test_cb_hasher());

        // Fallback keys off the mate coordinate (1400), so it sorts AFTER the
        // pos1=1200 template — the unchanged, pre-tc behavior.
        let filler = build_mapped_bam(0, 1200, b"fillr", &[]);
        let filler_key = extract_template_key_inline(&filler, &lib_lookup, None, &test_cb_hasher());

        assert!(
            supp_key > filler_key,
            "without tc, a sec/supp keeps the per-record (mate=1400) placement"
        );
    }

    /// Tie-case: with `(tid1, pos1, tid2, pos2)` held identical between a primary
    /// pair and its tc-stamped sec/supp reads, the tc branch must place the
    /// sec/supp read at the *exact same* `TemplateKey` as its primary — not merely
    /// nearby. This is the strong property the reviewer asked for: it pins down the
    /// full key (position lanes, strand lanes, CB, library, MI, name, **and** the
    /// `is_upper`/`is_read2` tiebreak), verifying tc ordering matches the mapped
    /// path rather than only asserting a coarse `<`/`>` placement.
    ///
    /// Geometry is a standard FR pair on `chr0`: R1 forward at unclipped-5' 1000
    /// (the lower/lane-1 end), R2 reverse at 1400 (the upper/lane-2 end). The two
    /// primaries key through the per-record path; the two sec/supp reads key
    /// through the tc branch. `is_read2` mirrors the primary `is_upper` here (R1 ->
    /// false, R2 -> true), so each sec/supp key is byte-identical to its primary's.
    #[rstest]
    #[case::secondary(flags::SECONDARY)]
    #[case::supplementary(flags::SUPPLEMENTARY)]
    fn test_extract_template_key_secondary_supplementary_matches_primary_via_tc(
        #[case] sec_supp_flag: u16,
    ) {
        const FIRST_SEGMENT: u16 = 0x40;
        const SECOND_SEGMENT: u16 = 0x80;
        const REVERSE: u16 = 0x10;
        const MATE_REVERSE: u16 = 0x20;

        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        let key =
            |bam: &[u8]| extract_template_key_inline(bam, &lib_lookup, None, &test_cb_hasher());

        // Primaries through the per-record path. R1 (forward, lower) is lane 1;
        // R2 (reverse, upper) is lane 2. Both canonicalize to the same coordinate.
        let r1_primary =
            build_primary_bam(1000, 1400, flags::PAIRED | FIRST_SEGMENT | MATE_REVERSE, b"tmpl");
        let r2_primary =
            build_primary_bam(1400, 1000, flags::PAIRED | SECOND_SEGMENT | REVERSE, b"tmpl");
        let r1_primary_key = key(&r1_primary);
        let r2_primary_key = key(&r2_primary);

        // The two primaries co-locate but differ only in the is_upper tiebreak.
        assert_eq!(
            r1_primary_key.core_cmp(&r2_primary_key),
            std::cmp::Ordering::Equal,
            "the pair's two primaries share (tid1, pos1, tid2, pos2)"
        );
        assert_ne!(
            r1_primary_key, r2_primary_key,
            "primaries differ only in the is_upper tiebreak (R1 lower, R2 upper)"
        );

        // Sec/supp reads through the tc branch, stamped with the pair's exact
        // canonical template coordinate, but aligned far away (own pos 9000).
        let tc = tc_aux(&[0, 1000, 0, 0, 1400, 1]);
        let r1_ss = build_secondary_or_supplementary_bam(
            9000,
            1400,
            sec_supp_flag | FIRST_SEGMENT,
            b"tmpl",
            &tc,
        );
        let r2_ss = build_secondary_or_supplementary_bam(
            9000,
            1000,
            sec_supp_flag | SECOND_SEGMENT,
            b"tmpl",
            &tc,
        );

        assert_eq!(
            key(&r1_ss),
            r1_primary_key,
            "R1 sec/supp keyed via tc must be byte-identical to the R1 primary key"
        );
        assert_eq!(
            key(&r2_ss),
            r2_primary_key,
            "R2 sec/supp keyed via tc must be byte-identical to the R2 primary key"
        );
    }

    #[test]
    fn test_extract_template_key_unmapped_with_cb() {
        let header = Header::builder().build();
        let lib_lookup = LibraryLookup::from_header(&header);
        let aux = cb_aux(b"ACGTACGT");

        // Build unmapped read (flag = PAIRED | UNMAPPED | MATE_UNMAPPED = 0x0D)
        let mut bam = vec![0u8; 32];
        bam[8] = 6; // l_read_name = 6 ("read1" + null)
        bam[14..16].copy_from_slice(&0x000Du16.to_le_bytes()); // flags
        // ref_id = -1 (unmapped)
        bam[0..4].copy_from_slice(&(-1i32).to_le_bytes());
        bam[4..8].copy_from_slice(&(-1i32).to_le_bytes()); // pos = -1
        bam[20..24].copy_from_slice(&(-1i32).to_le_bytes()); // mate_ref_id = -1
        bam[24..28].copy_from_slice(&(-1i32).to_le_bytes()); // mate_pos = -1
        bam.extend_from_slice(b"read1\0");
        bam.extend_from_slice(&aux);

        let key =
            extract_template_key_inline(&bam, &lib_lookup, Some(SamTag::CB), &test_cb_hasher());
        assert_ne!(key.cb_hash, 0, "unmapped read with CB should have non-zero cb_hash");
        assert_eq!(key.primary, u64::MAX, "unmapped both-mates should have MAX primary");
    }

    // ========================================================================
    // verify_dropped_lanes / dropped_lane_error tests
    // ========================================================================

    fn full_key(cb: u64, tertiary: u64) -> TemplateKey {
        TemplateKey { primary: 0, secondary: 0, cb_hash: cb, tertiary, name_hash_upper: 0 }
    }

    const LITE: TemplateKeyVariant = TemplateKeyVariant { cb: false, tertiary: false };

    #[test]
    fn verify_detects_cb_appearing() {
        let v = verify_dropped_lanes(&full_key(0, 0), &full_key(123, 0), LITE);
        assert_eq!(v, Some(DroppedLaneViolation::Cb));
        assert_eq!(v.unwrap().key_types_token(), "cb");
    }

    #[test]
    fn verify_detects_mi_change() {
        let first = full_key(0, 0);
        let cur = full_key(0, 0b10); // mi_value bit set, library still 0
        assert_eq!(verify_dropped_lanes(&first, &cur, LITE), Some(DroppedLaneViolation::Mi));
    }

    #[test]
    fn verify_detects_library_change() {
        let first = full_key(0, 0);
        let cur = full_key(0, 1u64 << 48); // library ordinal 1
        assert_eq!(verify_dropped_lanes(&first, &cur, LITE), Some(DroppedLaneViolation::Library));
    }

    #[test]
    fn verify_passes_when_dropped_lanes_constant() {
        let variant = TemplateKeyVariant { cb: true, tertiary: false };
        assert_eq!(verify_dropped_lanes(&full_key(1, 5), &full_key(2, 5), variant), None);
    }

    #[test]
    fn verify_error_message_names_field_and_token() {
        let e = dropped_lane_error("read42", DroppedLaneViolation::Library);
        let msg = e.to_string();
        assert!(
            msg.contains("read42")
                && msg.contains("library")
                && msg.contains("--key-types library"),
            "{msg}"
        );
    }

    // ========================================================================
    // select_template_variant tests
    // ========================================================================

    fn key(cb: u64, tertiary: u64) -> TemplateKey {
        TemplateKey { primary: 1, secondary: 2, cb_hash: cb, tertiary, name_hash_upper: 3 }
    }

    #[test]
    fn auto_lite_when_no_cb_no_tertiary() {
        let v = select_template_variant(Some(&key(0, 0)), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: false });
        assert_eq!(v.lanes(), 3);
    }

    #[test]
    fn auto_cb_only_when_first_record_has_cb() {
        let v = select_template_variant(Some(&key(42, 0)), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: true, tertiary: false });
        assert_eq!(v.lanes(), 4);
    }

    #[test]
    fn auto_tertiary_only_when_first_record_has_mi_or_library() {
        // 0xABCD is in the low-48 bits => MI present => keep tertiary, even though
        // the header library does not vary.
        let v = select_template_variant(Some(&key(0, 0xABCD)), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: true });
        assert_eq!(v.lanes(), 4);
    }

    #[test]
    fn auto_full_when_both_present() {
        let v = select_template_variant(Some(&key(7, 9)), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: true, tertiary: true });
        assert_eq!(v.lanes(), 5);
    }

    // tertiary high bits only = a single library ordinal, no MI. With a non-varying
    // header, the constant library lane is dropped -> lite.
    #[test]
    fn auto_drops_constant_single_library_lane() {
        let k = key(0, 1u64 << 48); // library ordinal 1, no MI
        let v = select_template_variant(Some(&k), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: false });
        assert_eq!(v.lanes(), 3);
    }

    #[test]
    fn auto_keeps_tertiary_when_header_library_varies() {
        let k = key(0, 1u64 << 48); // library ordinal 1, no MI
        let v = select_template_variant(Some(&k), KeyTypesSpec::Auto, true); // >1 library
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: true });
        assert_eq!(v.lanes(), 4);
    }

    #[test]
    fn auto_keeps_tertiary_for_mi_even_with_constant_library() {
        // MI present (low bits) -> keep tertiary even though header library is constant.
        let k = key(0, (1u64 << 48) | 0b10); // library ordinal 1 + MI value bit
        let v = select_template_variant(Some(&k), KeyTypesSpec::Auto, false);
        assert_eq!(v, TemplateKeyVariant { cb: false, tertiary: true });
    }

    #[test]
    fn none_forces_lite() {
        assert_eq!(
            select_template_variant(Some(&key(99, 99)), KeyTypesSpec::None, true).lanes(),
            3
        );
    }

    #[test]
    fn full_forces_all_lanes() {
        assert_eq!(select_template_variant(Some(&key(0, 0)), KeyTypesSpec::Full, false).lanes(), 5);
    }

    #[test]
    fn explicit_spec_passthrough() {
        assert_eq!(
            select_template_variant(
                None,
                KeyTypesSpec::Explicit { cb: true, tertiary: false },
                false
            )
            .lanes(),
            4
        );
    }

    #[test]
    fn empty_input_defaults_to_lite() {
        assert_eq!(select_template_variant(None, KeyTypesSpec::Auto, false).lanes(), 3);
    }

    /// `ref_width` is the bytes one accumulated template ref occupies for the
    /// variant; memory accounting multiplies records by it, so it must
    /// equal the real `size_of` for every lane, not a hand-kept table.
    #[rstest]
    #[case::k24(false, false, std::mem::size_of::<TemplateRecordRef<TemplateKey24>>())]
    #[case::cb32(true, false, std::mem::size_of::<TemplateRecordRef<CbKey32>>())]
    #[case::tert32(false, true, std::mem::size_of::<TemplateRecordRef<TertKey32>>())]
    #[case::k40(true, true, std::mem::size_of::<TemplateRecordRef<TemplateKey>>())]
    fn ref_width_matches_the_ref_size(
        #[case] cb: bool,
        #[case] tertiary: bool,
        #[case] size: usize,
    ) {
        assert_eq!(TemplateKeyVariant { cb, tertiary }.ref_width(), size);
    }
}
