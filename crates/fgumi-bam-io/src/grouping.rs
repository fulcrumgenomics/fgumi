//! Shared grouping and decoded-record domain types.
//!
//! Used by both the group command path (`grouper`, `mi_group`,
//! `commands::{group,dedup}` in the main crate) and the typed-step pipeline
//! (`pipeline::steps::group`, `steps::parse::decode`).
//!
//! These are BAM-record domain types — pre-computed grouping keys, the
//! decoded-record representation, the batching-weight and grouper traits.
//! They live here in `fgumi-bam-io`, next to [`crate::MemoryEstimate`] and
//! the [`fgumi_raw_bam`] raw-record helpers they operate on, rather than in
//! the pipeline crate that merely consumes them.

use std::io;
use std::sync::Arc;

use noodles::sam::alignment::record::data::field::Tag;

use crate::library::LibraryIndex;
use fgumi_raw_bam::{RawRecord, RawRecordView};

pub use crate::mem_estimate::MemoryEstimate;

// ============================================================================
// GroupKey - Pre-computed grouping key for fast comparison in Group step
// ============================================================================

/// Pre-computed grouping key for fast comparison in Group step.
///
/// All fields are integers/hashes for O(1) comparison. This is computed during
/// the parallel Decode step so the serial Group step only does integer comparisons.
///
/// For paired-end reads, positions are normalized so the lower position comes first.
/// For single-end reads, the mate fields use `UNKNOWN_*` sentinel values.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GroupKey {
    // Position info (normalized: lower position first)
    /// Reference sequence index for position 1 (lower).
    pub ref_id1: i32,
    /// Unclipped 5' position for position 1.
    pub pos1: i32,
    /// Strand for position 1 (0=forward, 1=reverse).
    pub strand1: u8,
    /// Reference sequence index for position 2 (higher or mate).
    pub ref_id2: i32,
    /// Unclipped 5' position for position 2.
    pub pos2: i32,
    /// Strand for position 2.
    pub strand2: u8,

    // Grouping metadata
    /// Library index (pre-computed from RG tag via header lookup).
    pub library_idx: u16,
    /// Hash of cell barcode (0 if none).
    pub cell_hash: u64,

    // For name-based grouping within position groups
    /// Hash of QNAME for fast name comparison.
    pub name_hash: u64,
}

impl GroupKey {
    /// Sentinel value for unknown reference ID (unpaired reads).
    pub const UNKNOWN_REF: i32 = i32::MAX;
    /// Sentinel value for unknown position (unpaired reads).
    pub const UNKNOWN_POS: i32 = i32::MAX;
    /// Sentinel value for unknown strand (unpaired reads).
    pub const UNKNOWN_STRAND: u8 = u8::MAX;

    /// Create a `GroupKey` for a paired-end read with mate info.
    ///
    /// Positions are automatically normalized so the lower position comes first.
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn paired(
        ref_id: i32,
        pos: i32,
        strand: u8,
        mate_ref_id: i32,
        mate_pos: i32,
        mate_strand: u8,
        library_idx: u16,
        cell_hash: u64,
        name_hash: u64,
    ) -> Self {
        // Normalize: put lower position first (matching ReadInfo behavior)
        let (ref_id1, pos1, strand1, ref_id2, pos2, strand2) =
            if (ref_id, pos, strand) <= (mate_ref_id, mate_pos, mate_strand) {
                (ref_id, pos, strand, mate_ref_id, mate_pos, mate_strand)
            } else {
                (mate_ref_id, mate_pos, mate_strand, ref_id, pos, strand)
            };

        Self { ref_id1, pos1, strand1, ref_id2, pos2, strand2, library_idx, cell_hash, name_hash }
    }

    /// Create a `GroupKey` for a single-end/unpaired read.
    #[must_use]
    pub fn single(
        ref_id: i32,
        pos: i32,
        strand: u8,
        library_idx: u16,
        cell_hash: u64,
        name_hash: u64,
    ) -> Self {
        Self {
            ref_id1: ref_id,
            pos1: pos,
            strand1: strand,
            ref_id2: Self::UNKNOWN_REF,
            pos2: Self::UNKNOWN_POS,
            strand2: Self::UNKNOWN_STRAND,
            library_idx,
            cell_hash,
            name_hash,
        }
    }

    /// Returns the position-only key for grouping by genomic position.
    ///
    /// This is used by `RecordPositionGrouper` to determine if records belong to
    /// the same position group (ignoring name).
    #[must_use]
    pub fn position_key(&self) -> (i32, i32, u8, i32, i32, u8, u16, u64) {
        (
            self.ref_id1,
            self.pos1,
            self.strand1,
            self.ref_id2,
            self.pos2,
            self.strand2,
            self.library_idx,
            self.cell_hash,
        )
    }

    /// Whether this key carries a mate position, i.e. it was built by
    /// [`Self::paired`] rather than [`Self::single`].
    ///
    /// `strand2` is set to [`Self::UNKNOWN_STRAND`] by the single-end
    /// constructor, so it is the discriminator between the two shapes.
    #[must_use]
    pub fn has_mate_position(&self) -> bool {
        self.strand2 != Self::UNKNOWN_STRAND
    }
}

impl PartialOrd for GroupKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for GroupKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.position_key()
            .cmp(&other.position_key())
            .then_with(|| self.name_hash.cmp(&other.name_hash))
    }
}

impl Default for GroupKey {
    fn default() -> Self {
        Self {
            ref_id1: Self::UNKNOWN_REF,
            pos1: Self::UNKNOWN_POS,
            strand1: Self::UNKNOWN_STRAND,
            ref_id2: Self::UNKNOWN_REF,
            pos2: Self::UNKNOWN_POS,
            strand2: Self::UNKNOWN_STRAND,
            library_idx: 0,
            cell_hash: 0,
            name_hash: 0,
        }
    }
}

// ============================================================================
// DecodedRecord - Record with pre-computed grouping key
// ============================================================================

/// A decoded BAM record with its pre-computed grouping key.
///
/// This is the output of the Decode step and input to the Group step.
/// The key is computed during the parallel Decode step so that the
/// serial Group step only needs to do fast integer comparisons.
///
/// # Cached-UMI invariant
///
/// `umi_value_offset` / `umi_value_len` cache the position of the UMI value
/// *within* `data`. The invariant is: **the cache is only valid while `data` is
/// unmodified.** Any mutation of `data` that could shift or overwrite the UMI
/// value bytes invalidates the cache. To enforce this cheaply, the sole
/// mutable-bytes accessor ([`Self::raw_bytes_mut`]) resets the cache to
/// [`Self::UMI_OFFSET_UNCACHED`] on hand-out, so a `cached_umi()` read after a
/// mutation falls back to a re-scan instead of slicing stale bytes.
#[derive(Debug)]
pub struct DecodedRecord {
    /// Pre-computed grouping key.
    pub key: GroupKey,
    /// Raw BAM record bytes.
    data: RawRecord,
    /// Cached record-relative offset of the UMI tag's value bytes (i.e. the
    /// first byte after the 2-byte tag header and the 1-byte type byte). The
    /// slice `data[umi_value_offset..umi_value_offset + umi_value_len]` yields
    /// the UMI bytes without the trailing NUL.
    ///
    /// Set to [`Self::UMI_OFFSET_UNCACHED`] when no UMI position was cached
    /// during decode (UMI tag missing, not Z-typed, or caching disabled).
    umi_value_offset: u32,
    /// Cached UMI value length in bytes, paired with `umi_value_offset`.
    umi_value_len: u16,
}

impl DecodedRecord {
    /// Sentinel value for `umi_value_offset` indicating no cached UMI position.
    /// Chosen as `u32::MAX` so it can never collide with a real BAM record
    /// offset (BAM records are bounded well under 4 GiB).
    pub const UMI_OFFSET_UNCACHED: u32 = u32::MAX;

    /// Create a decoded record from raw bytes, skipping noodles decode.
    ///
    /// Accepts anything that converts `Into<RawRecord>` (e.g. a bare `Vec<u8>` or
    /// an already-constructed `RawRecord`).
    #[must_use]
    pub fn from_raw_bytes(raw: impl Into<RawRecord>, key: GroupKey) -> Self {
        Self {
            key,
            data: raw.into(),
            umi_value_offset: Self::UMI_OFFSET_UNCACHED,
            umi_value_len: 0,
        }
    }

    /// Attach a cached UMI value position to this decoded record.
    ///
    /// `umi_value_offset` is the record-relative offset of the first UMI value
    /// byte (after the tag header), and `umi_value_len` is the value length
    /// (excluding any trailing NUL).
    pub fn set_cached_umi(&mut self, umi_value_offset: u32, umi_value_len: u16) {
        self.umi_value_offset = umi_value_offset;
        self.umi_value_len = umi_value_len;
    }

    /// Returns the cached UMI bytes (without trailing NUL) if a position was
    /// recorded during decode and still falls within the raw record bytes.
    /// Returns `None` if no cache is set or the cached position is out of range.
    #[must_use]
    pub fn cached_umi(&self) -> Option<&[u8]> {
        if self.umi_value_offset == Self::UMI_OFFSET_UNCACHED {
            return None;
        }
        let start = self.umi_value_offset as usize;
        let end = start.checked_add(self.umi_value_len as usize)?;
        let value = self.data.as_ref().get(start..end)?;
        // The cache only ever holds a Z-typed UMI value (NUL-free by the BAM
        // spec; the trailing NUL is excluded by `set_cached_umi`). If a mutation
        // shifted the record bytes such that this offset now lands on a different
        // region, the slice would typically contain an interior NUL — a cheap,
        // zero-release-cost canary for a stale cache that survived bounds checks.
        debug_assert!(
            !value.contains(&0),
            "stale cached UMI: offset {start} sliced bytes containing an interior NUL \
             (the record was likely mutated without invalidating the UMI cache)"
        );
        Some(value)
    }

    /// Returns the cached UMI offset ([`Self::UMI_OFFSET_UNCACHED`] when absent)
    /// and length.
    #[must_use]
    pub fn cached_umi_position(&self) -> (u32, u16) {
        (self.umi_value_offset, self.umi_value_len)
    }

    /// Returns the cached UMI `(offset, len)` position, or `None` when no
    /// position was recorded during decode (the [`Self::UMI_OFFSET_UNCACHED`]
    /// sentinel). Keeps the sentinel check at the owning layer so callers can
    /// branch on a plain `Option` instead of comparing against the sentinel.
    #[must_use]
    pub fn cached_umi_position_opt(&self) -> Option<(u32, u16)> {
        if self.umi_value_offset == Self::UMI_OFFSET_UNCACHED {
            None
        } else {
            Some((self.umi_value_offset, self.umi_value_len))
        }
    }

    /// Returns a reference to the raw bytes.
    #[must_use]
    pub fn raw_bytes(&self) -> &[u8] {
        self.data.as_ref()
    }

    /// Immutable access to the underlying [`RawRecord`], for read-only
    /// consumers (e.g. the FASTQ-encode step) that need the typed accessors
    /// (`flags()`, sequence, qualities, tags) rather than the raw byte slice.
    #[must_use]
    pub fn record(&self) -> &RawRecord {
        &self.data
    }

    /// Mutable access to the underlying [`RawRecord`]. Used by mid-chain
    /// `Parallel` Process steps that need to mutate record bytes (e.g.,
    /// MQ bumping, tag rewriting) without rebuilding the `DecodedRecord`
    /// or invalidating its pre-computed `GroupKey`. Caller must not
    /// change the record's identity (qname, library bytes, cell barcode)
    /// — those feed `key`, which we don't recompute here.
    ///
    /// Handing out mutable bytes resets the cached UMI position to
    /// [`Self::UMI_OFFSET_UNCACHED`] (see the type-level cached-UMI invariant):
    /// a length-changing edit before the UMI tag could shift the value while
    /// leaving the old offset in-bounds, so a later `cached_umi()` would slice
    /// stale-but-valid bytes. Clearing forces a re-scan instead. This is a
    /// single field store and the cache is only repopulated by the decode step.
    pub fn raw_bytes_mut(&mut self) -> &mut RawRecord {
        self.umi_value_offset = Self::UMI_OFFSET_UNCACHED;
        self.umi_value_len = 0;
        &mut self.data
    }

    /// Takes the [`RawRecord`] out.
    #[must_use]
    pub fn into_raw_bytes(self) -> RawRecord {
        self.data
    }
}

impl MemoryEstimate for DecodedRecord {
    fn estimate_heap_size(&self) -> usize {
        // RawRecord::capacity() returns the inner Vec<u8> capacity.
        self.data.capacity()
    }
}
// Vec<DecodedRecord>, Vec<RecordBuf>, Vec<u8>, RecordBuf, () — all provided
// by the blanket/foreign impls in fgumi_bam_io::mem_estimate.

// ============================================================================
// GroupKeyConfig - Configuration for computing `GroupKey` during Decode
// ============================================================================

/// Configuration for computing `GroupKey` during the Decode step.
///
/// When this is provided to the pipeline, the Decode step will compute
/// full `GroupKey` values for each record. This moves expensive computations
/// (CIGAR parsing, tag extraction) from the serial Group step to the parallel
/// Decode step.
#[derive(Debug, Clone)]
pub struct GroupKeyConfig {
    /// Library index for fast RG → library lookup.
    pub library_index: Arc<LibraryIndex>,
    /// Tag used for cell barcode extraction. None skips cell extraction.
    pub cell_tag: Option<Tag>,
    /// How much of the [`GroupKey`] the Decode step computes per record. See
    /// [`KeyMode`]: `Full` (position + RG/CB/MC + name hash), `NameHashOnly`
    /// (name hash alone, for queryname grouping), or `None` (no key at all —
    /// `GroupKey::default()`, for stages that never read the key, e.g. the
    /// BAM→FASTQ encode).
    pub key_mode: KeyMode,
    /// UMI tag (raw 2-byte form) whose value position should be cached on each
    /// [`DecodedRecord`] during decode. `None` disables caching — downstream
    /// code must fall back to scanning aux data. This is orthogonal to
    /// `key_mode`: the UMI cache scan is gated solely on `umi_tag` being set
    /// (in practice only the Group stage sets it).
    pub umi_tag: Option<[u8; 2]>,
}

/// How much of the [`GroupKey`] the Decode step computes for each record.
///
/// Each variant does strictly less work than the one above it (declaration
/// order `Full` → `NameHashOnly` → `None`); a stage picks the cheapest level
/// whose output it actually reads. Changing the level only
/// changes the (discarded) key, never the record bytes, so output is unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum KeyMode {
    /// Full key: the CIGAR 5′-position walk, the RG/CB/MC aux-tag extraction
    /// pass, and the read-name hash. Needed by group/dedup/consensus.
    #[default]
    Full,
    /// Read-name hash only; the rest of the [`GroupKey`] stays at default.
    /// Skips the CIGAR walk and the aux-tag pass. For queryname-grouping stages
    /// (e.g. `correct`) that read only [`GroupKey::name_hash`]. See
    /// [`name_hash_key`].
    NameHashOnly,
    /// No key at all: the Decode step emits `GroupKey::default()` and computes
    /// nothing. For stages that never read any `GroupKey` field — e.g. the
    /// terminal BAM→FASTQ encode, which reads only the raw record. Skips even
    /// the per-record name hash `NameHashOnly` still pays.
    None,
}

impl GroupKeyConfig {
    /// Create a new `GroupKeyConfig` that computes the full position/cell key.
    #[must_use]
    pub fn new(library_index: LibraryIndex, cell_tag: Tag) -> Self {
        Self {
            library_index: Arc::new(library_index),
            cell_tag: Some(cell_tag),
            key_mode: KeyMode::Full,
            umi_tag: None,
        }
    }

    /// Create a `GroupKeyConfig` without cell barcode extraction.
    #[must_use]
    pub fn new_raw_no_cell(library_index: LibraryIndex) -> Self {
        Self {
            library_index: Arc::new(library_index),
            cell_tag: None,
            key_mode: KeyMode::Full,
            umi_tag: None,
        }
    }

    /// Create a `GroupKeyConfig` that computes only the read-name hash.
    ///
    /// For queryname-grouping stages (e.g. `correct`) whose grouper reads only
    /// [`GroupKey::name_hash`]; skips the CIGAR position walk and the aux-tag
    /// extraction pass. `library_index` is retained (unused for the key) so
    /// the config shape is uniform across the pipeline.
    #[must_use]
    pub fn name_hash_only(library_index: LibraryIndex) -> Self {
        Self {
            library_index: Arc::new(library_index),
            cell_tag: None,
            key_mode: KeyMode::NameHashOnly,
            umi_tag: None,
        }
    }

    /// Create a `GroupKeyConfig` that computes no key at all: the Decode step
    /// emits `GroupKey::default()` for every record. For stages that never read
    /// any `GroupKey` field — e.g. the terminal BAM→FASTQ encode, which reads
    /// only the raw record (name/flags/SEQ/QUAL, and the UMI tag straight off
    /// the record). Cheapest of all: skips even the per-record name hash.
    /// `library_index` is retained (unused) so the config shape is uniform.
    #[must_use]
    pub fn no_key(library_index: LibraryIndex) -> Self {
        Self {
            library_index: Arc::new(library_index),
            cell_tag: None,
            key_mode: KeyMode::None,
            umi_tag: None,
        }
    }

    /// Enable UMI position caching for the given raw 2-byte UMI tag (e.g. `*b"RX"`).
    ///
    /// The Decode step will record the UMI value's record-relative position on
    /// each [`DecodedRecord`] so downstream UMI lookups can slice it directly
    /// via [`DecodedRecord::cached_umi`] without re-scanning aux data.
    #[must_use]
    pub fn with_umi_tag(mut self, umi_tag: [u8; 2]) -> Self {
        self.umi_tag = Some(umi_tag);
        self
    }
}

impl Default for GroupKeyConfig {
    fn default() -> Self {
        Self {
            library_index: Arc::new(LibraryIndex::default()),
            cell_tag: Some(Tag::from([b'C', b'B'])), // Default cell barcode tag (CB)
            key_mode: KeyMode::Full,
            umi_tag: None,
        }
    }
}

/// Hash a record's read name, mapping an empty name to `hash_name(None)` so the
/// raw path matches the noodles path (which sees `None` for an empty name).
///
/// Shared by [`name_hash_key`] and [`compute_group_key_from_raw`]: both must
/// produce the same hash for the same record, and a second copy of this block
/// could drift out of parity silently.
fn raw_name_hash(raw: &[u8]) -> u64 {
    let name = fgumi_raw_bam::read_name(raw);
    if name.is_empty() {
        LibraryIndex::hash_name(None)
    } else {
        LibraryIndex::hash_name(Some(name))
    }
}

/// Compute a [`GroupKey`] containing only the read-name hash, leaving all
/// position/strand/library/cell fields at their default.
///
/// Reproduces exactly the `name_hash` that [`compute_group_key_from_raw`]
/// computes (empty name → `hash_name(None)`), so a queryname grouper sees an
/// identical hash. Skips the CIGAR 5′-position walk and the RG/CB/MC aux-tag
/// extraction pass.
///
/// # Panics
///
/// Panics if `raw` is not a validated BAM record payload (same contract as
/// [`compute_group_key_from_raw`]).
#[must_use]
pub fn name_hash_key(raw: &[u8]) -> GroupKey {
    GroupKey { name_hash: raw_name_hash(raw), ..GroupKey::default() }
}

/// Compute a record's [`GroupKey`] at the requested [`KeyMode`]: `None` →
/// `GroupKey::default()` (no work), `NameHashOnly` → [`name_hash_key`], `Full`
/// → [`compute_group_key_from_raw`]. The single place the three levels are
/// mapped, so every decode consumer (BAM and SAM) stays in parity.
///
/// # Panics
///
/// For `NameHashOnly`/`Full`, panics if `raw` is not a validated BAM record
/// payload (same contract as the underlying key functions). `None` never
/// inspects `raw`.
#[must_use]
pub fn key_for_mode(
    mode: KeyMode,
    raw: &[u8],
    library_index: &LibraryIndex,
    cell_tag: Option<Tag>,
) -> GroupKey {
    key_and_umi_for_mode(mode, raw, library_index, cell_tag, None).0
}

/// Like [`key_for_mode`], but also captures the UMI tag's value position in the
/// same aux-data scan that builds the key (only `KeyMode::Full` reads aux data).
///
/// Returns the key plus the record-relative `(offset, len)` of the UMI value
/// bytes when `umi_tag` is set and found during the key's aux scan; `None`
/// otherwise (UMI not requested, tag absent, or a key mode / degenerate path
/// that does not read aux data). For `None`, callers that need the UMI cached
/// fall back to a standalone, crate-external scan (`cache_umi_position` in the
/// pipeline decode step), preserving behaviour while avoiding a second aux pass
/// on the common `Full` path. See `compute_group_key_and_umi_from_raw` for the
/// exact tag-resolution semantics (which match how this scan already resolves
/// RG/cell/MC, not the standalone scan's, on spec-violating duplicate tags).
///
/// # Panics
///
/// For `NameHashOnly`/`Full`, panics if `raw` is not a validated BAM record
/// payload (same contract as [`key_for_mode`] and the underlying key functions).
/// `None` never inspects `raw`.
#[must_use]
pub fn key_and_umi_for_mode(
    mode: KeyMode,
    raw: &[u8],
    library_index: &LibraryIndex,
    cell_tag: Option<Tag>,
    umi_tag: Option<[u8; 2]>,
) -> (GroupKey, Option<(u32, u16)>) {
    match mode {
        KeyMode::None => (GroupKey::default(), None),
        KeyMode::NameHashOnly => (name_hash_key(raw), None),
        KeyMode::Full => compute_group_key_and_umi_from_raw(raw, library_index, cell_tag, umi_tag),
    }
}

/// Compute a `GroupKey` directly from raw BAM bytes, matching `compute_group_key()` exactly.
///
/// Uses 1-based coordinate helpers to produce identical keys to the noodles path.
///
/// # Panics
///
/// Panics if `raw` is not a validated BAM record payload (as produced by the BAM
/// reader / raw-record pipeline). Callers must not pass arbitrary external bytes;
/// raw-field accessors will panic on malformed or truncated input.
#[must_use]
pub fn compute_group_key_from_raw(
    raw: &[u8],
    library_index: &LibraryIndex,
    cell_tag: Option<noodles::sam::alignment::record::data::field::Tag>,
) -> GroupKey {
    compute_group_key_and_umi_from_raw(raw, library_index, cell_tag, None).0
}

/// Return the aux-data offset (via [`fgumi_raw_bam::aux_data_offset_from_record`])
/// alongside the aux-data slice, sharing the single header decode. Mirrors
/// [`fgumi_raw_bam::aux_data_slice`]'s bounds guard so the two agree on when aux
/// data is absent (truncated / out-of-range records yield `(None, &[])`).
#[inline]
fn aux_slice_with_offset(raw: &[u8]) -> (Option<usize>, &[u8]) {
    match fgumi_raw_bam::aux_data_offset_from_record(raw) {
        Some(offset) if offset <= raw.len() => (Some(offset), &raw[offset..]),
        _ => (None, &[]),
    }
}

/// Convert an aux-relative UMI value position (as returned by
/// [`fgumi_raw_bam::extract_aux_string_tags`]) into the record-relative
/// `(offset, len)` expected by [`crate::DecodedRecord::set_cached_umi`]. Uses the
/// same offset arithmetic and width guards (`u32`-narrow the aux offset,
/// `checked_add` the aux-relative offset) as the standalone `cache_umi_position`
/// scan, so a folded capture yields the same cached position for any spec-legal
/// record. The *tag lookup* upstream of this conversion is
/// `extract_aux_string_tags`, whose duplicate/mistyped-tag semantics differ from
/// the standalone scan's — see `compute_group_key_and_umi_from_raw`.
#[inline]
fn record_relative_umi(
    aux_offset: Option<usize>,
    umi_position: Option<(u32, u16)>,
) -> Option<(u32, u16)> {
    let (aux_rel_off, len) = umi_position?;
    let base = u32::try_from(aux_offset?).ok()?;
    let offset = base.checked_add(aux_rel_off)?;
    Some((offset, len))
}

/// Like [`compute_group_key_from_raw`], but also captures the UMI tag's value
/// position in the SAME aux-data scan used to build the key, avoiding a second
/// pass over aux data during decode.
///
/// Returns the key plus the record-relative `(offset, len)` of the UMI value
/// bytes (excluding the trailing NUL) when `umi_tag` is set and the record is
/// keyed via a path that scans aux data (primary reads, and secondary /
/// supplementary reads carrying a `tc` tag). Returns `None` for the UMI position
/// when `umi_tag` is `None`, the tag is absent, or the key falls back to a
/// name-only path that never reads aux data — callers cache the UMI via a
/// standalone scan (`cache_umi_position`) in those cases.
///
/// # Tag-resolution semantics
///
/// The UMI is resolved by the SAME [`fgumi_raw_bam::extract_aux_string_tags`]
/// pass that resolves RG/cell/MC for the key, so the cached UMI is consistent
/// with the key's own tag resolution. For **spec-legal** records — each aux tag
/// present at most once (SAM §1.5) — this yields the identical position the
/// standalone `find_string_tag_position` scan (used by the `cache_umi_position`
/// fallback) would. The two deliberately differ only on **malformed** records
/// with a duplicated or type-shadowed UMI tag: `extract_aux_string_tags` skips a
/// non-`Z` entry and takes a later `Z` copy, whereas `find_tag_position` resolves
/// the first id match and rejects it if non-`Z`. This is the exact,
/// already-accepted divergence documented for the `MC` tag on `validate_mc_tag`
/// (see `src/lib/grouper.rs`) — the relaxation makes the cached value agree with
/// the value the grouping key actually uses, rather than losing it.
///
/// # Panics
///
/// Panics if `raw` is not a validated BAM record payload (as produced by the BAM
/// reader / raw-record pipeline). Callers must not pass arbitrary external bytes;
/// raw-field accessors will panic on malformed or truncated input.
#[must_use]
pub(crate) fn compute_group_key_and_umi_from_raw(
    raw: &[u8],
    library_index: &LibraryIndex,
    cell_tag: Option<noodles::sam::alignment::record::data::field::Tag>,
    umi_tag: Option<[u8; 2]>,
) -> (GroupKey, Option<(u32, u16)>) {
    // Extract name hash (match noodles path: empty name → None → hash 0)
    let name_hash = raw_name_hash(raw);

    // Check secondary/supplementary
    let flg = RawRecordView::new(raw).flags();
    let is_secondary = (flg & fgumi_raw_bam::flags::SECONDARY) != 0;
    let is_supplementary = (flg & fgumi_raw_bam::flags::SUPPLEMENTARY) != 0;
    if is_secondary || is_supplementary {
        // A secondary/supplementary read cannot compute its own template
        // coordinate (it lacks its own primary's position). When `fgumi zipper`
        // has stamped the exact coordinate into `tc`, key on it so the read
        // groups into the same position group as its primary — instead of an
        // UNKNOWN position that relies on the read sorting adjacent to its
        // primary. Library/cell are extracted the same way as primaries so the
        // position key matches (they share their template's RG/CB).
        let (aux_offset, aux_data) = aux_slice_with_offset(raw);
        if let Some([tid1, pos1, neg1, tid2, pos2, neg2]) =
            fgumi_raw_bam::read_tc_template_coordinate(aux_data)
        {
            let cell_tag_bytes = cell_tag.map_or([0u8; 2], |t| [t.as_ref()[0], t.as_ref()[1]]);
            let aux_tags =
                fgumi_raw_bam::extract_aux_string_tags(aux_data, cell_tag_bytes, umi_tag);
            let library_idx =
                aux_tags.rg.map_or(0, |rg| library_index.get(LibraryIndex::hash_rg(rg)));
            let cell_hash = aux_tags.cell.map_or(0, |cb| LibraryIndex::hash_cell_barcode(Some(cb)));
            let umi = record_relative_umi(aux_offset, aux_tags.umi_position);
            return (
                GroupKey::paired(
                    tid1,
                    pos1,
                    u8::from(neg1 != 0),
                    tid2,
                    pos2,
                    u8::from(neg2 != 0),
                    library_idx,
                    cell_hash,
                    name_hash,
                ),
                umi,
            );
        }
        return (GroupKey { name_hash, ..GroupKey::default() }, None);
    }

    // Own position (1-based, matching noodles) — zero-allocation CIGAR iteration
    let reverse = (flg & fgumi_raw_bam::flags::REVERSE) != 0;
    let own_pos = fgumi_raw_bam::unclipped_5prime_from_raw_bam(raw);

    // A mapped record with an empty or truncated CIGAR has no computable
    // unclipped 5' position, so `unclipped_5prime_from_raw_bam` returns the
    // `i32::MAX` sentinel. Fall back to the name-only key rather than letting the
    // sentinel reach a position slot — otherwise distinct templates sharing
    // ref/strand/library/cell would collide on `i32::MAX` (`position_key`
    // excludes `name_hash`). Matches the secondary/supplementary fallback above.
    if own_pos == i32::MAX {
        return (GroupKey { name_hash, ..GroupKey::default() }, None);
    }

    let own_ref_id = fgumi_raw_bam::ref_id(raw);
    let strand = u8::from(reverse);

    // Single-pass aux tag extraction (RG, cell barcode, MC, and — when requested
    // — the UMI value position, folded into this one scan that builds the key).
    let (aux_offset, aux_data) = aux_slice_with_offset(raw);
    let cell_tag_bytes = cell_tag.map_or([0u8; 2], |t| [t.as_ref()[0], t.as_ref()[1]]);
    let aux_tags = fgumi_raw_bam::extract_aux_string_tags(aux_data, cell_tag_bytes, umi_tag);
    let umi = record_relative_umi(aux_offset, aux_tags.umi_position);

    let library_idx = if let Some(rg) = aux_tags.rg {
        let rg_hash = LibraryIndex::hash_rg(rg);
        library_index.get(rg_hash)
    } else {
        0
    };

    let cell_hash =
        if let Some(cb) = aux_tags.cell { LibraryIndex::hash_cell_barcode(Some(cb)) } else { 0 };

    // Check if paired
    let is_paired = (flg & fgumi_raw_bam::flags::PAIRED) != 0;
    if !is_paired {
        return (
            GroupKey::single(own_ref_id, own_pos, strand, library_idx, cell_hash, name_hash),
            umi,
        );
    }

    // Mate info — guard against MATE_UNMAPPED (matching noodles path)
    let mate_unmapped = (flg & fgumi_raw_bam::flags::MATE_UNMAPPED) != 0;
    let mate_reverse = (flg & fgumi_raw_bam::flags::MATE_REVERSE) != 0;
    let mate_strand = u8::from(mate_reverse);
    let raw_mate_ref_id = fgumi_raw_bam::mate_ref_id(raw);
    let raw_mate_pos = fgumi_raw_bam::mate_pos(raw);

    // Get mate unclipped 5' position via MC tag (skip if mate is unmapped)
    let mate_pos_result = if mate_unmapped {
        None
    } else {
        aux_tags
            .mc
            .map(|mc| fgumi_raw_bam::mate_unclipped_5prime_1based(raw_mate_pos, mate_reverse, mc))
    };

    let key = match mate_pos_result {
        Some(mp) => GroupKey::paired(
            own_ref_id,
            own_pos,
            strand,
            raw_mate_ref_id,
            mp,
            mate_strand,
            library_idx,
            cell_hash,
            name_hash,
        ),
        None => {
            // No MC tag — fall back to single-end behavior
            GroupKey::single(own_ref_id, own_pos, strand, library_idx, cell_hash, name_hash)
        }
    };
    (key, umi)
}

/// Groups a stream of in-order [`DecodedRecord`]s into completed groups.
///
/// Implementors maintain partial groups across `add_records` calls and emit
/// completed ones; `finish` flushes any trailing partial group at EOF. Used
/// by the pipeline's Group step and the standalone grouping commands.
pub trait Grouper: Send {
    /// The type of group produced by this grouper.
    type Group: Send;

    /// Add decoded records to the grouper.
    ///
    /// Records are guaranteed to be in order (from template-coordinate sorted BAM).
    /// The grouper maintains partial groups waiting for more records.
    ///
    /// Each `DecodedRecord` contains the record plus a pre-computed `GroupKey`
    /// for fast comparison (position, name hash, library, etc.).
    ///
    /// Returns completed groups (may be empty if more records are needed).
    ///
    /// # Errors
    ///
    /// Returns an I/O error if grouping logic encounters invalid data.
    fn add_records(&mut self, records: Vec<DecodedRecord>) -> io::Result<Vec<Self::Group>>;

    /// Signal that no more input will arrive (EOF).
    ///
    /// Returns any remaining partial group.
    ///
    /// # Errors
    ///
    /// Returns an I/O error if finalizing the grouper fails.
    fn finish(&mut self) -> io::Result<Option<Self::Group>>;

    /// Returns true if the grouper has a partial group.
    fn has_pending(&self) -> bool;
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::SamBuilder;
    use fgumi_raw_bam::SamTag;
    use fgumi_raw_bam::flags;
    use rstest::rstest;

    // ========================================================================
    // compute_group_key_from_raw — primary (fully-populated) path
    // ========================================================================

    /// Read group `RG1` resolves to library `libA`; `RG2` to `libB`.
    fn library_index_with_two_groups() -> LibraryIndex {
        use noodles::sam::header::record::value::Map;
        use noodles::sam::header::record::value::map::ReadGroup;
        use noodles::sam::header::record::value::map::read_group::tag as rg_tag;

        let mut header = noodles::sam::Header::builder();
        for (id, library) in [("RG1", "libA"), ("RG2", "libB")] {
            let rg = Map::<ReadGroup>::builder()
                .insert(rg_tag::LIBRARY, String::from(library))
                .build()
                .expect("read group builds");
            header = header.add_read_group(bstr::BString::from(id), rg);
        }
        LibraryIndex::from_header(&header.build()).expect("builds")
    }

    /// Build one mapped mate of a pair carrying `RG`, `CB` and `MC`.
    fn paired_record_with_tags(
        pos: i32,
        mate_pos: i32,
        reverse: bool,
        mate_reverse: bool,
    ) -> fgumi_raw_bam::RawRecord {
        let mut flag = flags::PAIRED;
        if reverse {
            flag |= flags::REVERSE;
        }
        if mate_reverse {
            flag |= flags::MATE_REVERSE;
        }
        let mut b = SamBuilder::new();
        b.ref_id(0)
            .pos(pos)
            .flags(flag)
            .mate_ref_id(0)
            .mate_pos(mate_pos)
            .read_name(b"pair_full")
            .cigar_ops(&[cigar_m(50)])
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30])
            .add_string_tag(SamTag::RG, b"RG1")
            .add_string_tag(SamTag::CB, b"ACGTACGTA")
            .add_string_tag(SamTag::MC, b"50M");
        b.build()
    }

    /// The fully-populated branch: paired, mate mapped, `MC` present, `RG`
    /// resolving through a non-default `LibraryIndex`, and a cell barcode read via
    /// `cell_tag`. Every other test in this module lands on a fallback (secondary,
    /// supplementary, or paired-without-`MC`) with `LibraryIndex::default()` and
    /// `cell_tag: None`, so `library_idx`, `cell_hash` and the mate slot are
    /// otherwise never exercised — and all three participate in key equality.
    #[test]
    fn paired_with_mc_populates_mate_library_and_cell_fields() {
        let lib = library_index_with_two_groups();
        let cb = noodles::sam::alignment::record::data::field::Tag::CELL_BARCODE_ID;

        let rec = paired_record_with_tags(1000, 1200, false, true);
        let key = compute_group_key_from_raw(rec.as_ref(), &lib, Some(cb));

        // Mate slot is populated, not left at the single-end sentinels.
        assert_ne!(key.ref_id2, GroupKey::UNKNOWN_REF, "MC must populate the mate slot");
        assert_ne!(key.strand2, GroupKey::UNKNOWN_STRAND);

        // Pin the mate coordinate: it must be the MC-derived unclipped 5' position,
        // not the raw `mate_pos`. `!= UNKNOWN_POS` alone would accept a wrong
        // CIGAR walk.
        let expected_mate_pos = fgumi_raw_bam::mate_unclipped_5prime_1based(
            fgumi_raw_bam::mate_pos(rec.as_ref()),
            true,
            b"50M",
        );
        assert_eq!(key.pos2, expected_mate_pos, "mate slot must hold the MC-derived position");

        // RG resolved through the index — not the unknown bucket.
        assert_eq!(key.library_idx, lib.get(LibraryIndex::hash_rg(b"RG1")));
        assert_ne!(key.library_idx, 0, "RG1 must resolve to a real library");

        // Cell barcode hashed because `cell_tag` named CB.
        assert_eq!(key.cell_hash, LibraryIndex::hash_cell_barcode(Some(b"ACGTACGTA")));
        assert_ne!(key.cell_hash, 0);

        assert_eq!(key.name_hash, LibraryIndex::hash_name(Some(b"pair_full")));

        // Own slot holds this record's own unclipped 5' position.
        assert_eq!(key.pos1, fgumi_raw_bam::unclipped_5prime_from_raw_bam(rec.as_ref()));
        assert_eq!(key.ref_id1, fgumi_raw_bam::ref_id(rec.as_ref()));
    }

    /// Both mates of the same template must normalize to one key — this is what
    /// makes them group together. The mate sees own/mate swapped and the strands
    /// exchanged.
    #[test]
    fn paired_with_mc_normalizes_both_mates_to_the_same_key() {
        let lib = library_index_with_two_groups();
        let cb = noodles::sam::alignment::record::data::field::Tag::CELL_BARCODE_ID;

        let forward = paired_record_with_tags(1000, 1200, false, true);
        let reverse = paired_record_with_tags(1200, 1000, true, false);

        let key_forward = compute_group_key_from_raw(forward.as_ref(), &lib, Some(cb));
        let key_reverse = compute_group_key_from_raw(reverse.as_ref(), &lib, Some(cb));

        assert_eq!(key_forward, key_reverse, "both mates of a template must share one key");
    }

    /// A secondary/supplementary read cannot derive its own template coordinate,
    /// so `fgumi zipper` stamps the primary's into the `tc` aux array. Keying on
    /// it puts the read in the SAME position group as its primary; without it the
    /// read falls back to a name-only key and relies on sorting adjacency.
    #[test]
    fn secondary_with_tc_tag_keys_on_the_stamped_template_coordinate() {
        let lib = library_index_with_two_groups();

        let mut b = SamBuilder::new();
        b.ref_id(0)
            .pos(5000)
            .flags(flags::PAIRED | flags::SECONDARY)
            .read_name(b"sec_with_tc")
            .cigar_ops(&[cigar_m(50)])
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30])
            .add_string_tag(SamTag::RG, b"RG1")
            // tc = [tid1, pos1, neg1, tid2, pos2, neg2] — the primary's coordinate.
            .add_array_i32(SamTag::TC, &[0, 1001, 0, 0, 1249, 1]);
        let rec = b.build();

        let key = compute_group_key_from_raw(rec.as_ref(), &lib, None);

        // Both slots come from `tc`, NOT from this record's own pos (5000).
        assert_eq!((key.ref_id1, key.pos1, key.strand1), (0, 1001, 0));
        assert_eq!((key.ref_id2, key.pos2, key.strand2), (0, 1249, 1));
        assert_ne!(key.pos1, 5001, "the record's own position must not be used");

        // Library still resolves, so the key matches its primary's.
        assert_eq!(key.library_idx, lib.get(LibraryIndex::hash_rg(b"RG1")));
        assert_eq!(key.name_hash, LibraryIndex::hash_name(Some(b"sec_with_tc")));
    }

    /// Without a `tc` tag the secondary path still falls back to the name-only
    /// key — including for a mapped record with an EMPTY cigar, where
    /// `unclipped_5prime_from_raw_bam` returns `i32::MAX`. The fallback must not
    /// leak that sentinel into a position slot.
    #[test]
    fn secondary_without_tc_falls_back_to_name_only_even_with_an_empty_cigar() {
        let lib = library_index_with_two_groups();

        let mut b = SamBuilder::new();
        b.ref_id(0)
            .pos(5000)
            .flags(flags::SUPPLEMENTARY)
            .read_name(b"sup_no_tc")
            .cigar_ops(&[])
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30]);
        let rec = b.build();

        let key = compute_group_key_from_raw(rec.as_ref(), &lib, None);

        assert_eq!(key, GroupKey { name_hash: key.name_hash, ..GroupKey::default() });
        assert_eq!(key.name_hash, LibraryIndex::hash_name(Some(b"sup_no_tc")));
        assert_eq!(key.pos1, GroupKey::default().pos1, "no i32::MAX sentinel may leak in");
    }

    /// The PRIMARY path has the same sentinel hazard: a mapped, unpaired primary
    /// with an empty CIGAR has no computable unclipped 5' position, so it must
    /// fall back to the name-only key too. Without the guard the `i32::MAX`
    /// sentinel would reach `pos1` and merge distinct templates sharing
    /// ref/strand/library/cell (`position_key` excludes `name_hash`).
    #[test]
    fn primary_with_empty_cigar_falls_back_to_name_only() {
        let lib = library_index_with_two_groups();

        // Default flags (0) => mapped, unpaired, primary.
        let mut b = SamBuilder::new();
        b.ref_id(0)
            .pos(5000)
            .read_name(b"primary_no_cigar")
            .cigar_ops(&[])
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30]);
        let rec = b.build();

        let key = compute_group_key_from_raw(rec.as_ref(), &lib, None);

        assert_eq!(key, GroupKey { name_hash: key.name_hash, ..GroupKey::default() });
        assert_eq!(key.name_hash, LibraryIndex::hash_name(Some(b"primary_no_cigar")));
        assert_eq!(key.pos1, GroupKey::default().pos1, "no i32::MAX sentinel may leak in");
    }

    /// `library_idx` is part of key equality: the same position read under a
    /// different read group must NOT group with it.
    #[test]
    fn differing_library_splits_the_key() {
        let lib = library_index_with_two_groups();

        let mut b = SamBuilder::new();
        b.ref_id(0)
            .pos(1000)
            .flags(flags::PAIRED)
            .mate_ref_id(0)
            .mate_pos(1200)
            .read_name(b"pair_full")
            .cigar_ops(&[cigar_m(50)])
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30])
            .add_string_tag(SamTag::RG, b"RG2")
            .add_string_tag(SamTag::MC, b"50M");
        let other_library = b.build();

        let key_a = compute_group_key_from_raw(
            paired_record_with_tags(1000, 1200, false, true).as_ref(),
            &lib,
            None,
        );
        let key_b = compute_group_key_from_raw(other_library.as_ref(), &lib, None);

        assert_ne!(key_a.library_idx, key_b.library_idx);
        assert_ne!(key_a, key_b, "a different library must not group with the first");
    }

    /// CIGAR op `(len << 4) | op_code`; op 0 = `M`.
    fn cigar_m(len: u32) -> u32 {
        len << 4
    }

    // ========================================================================
    // GroupKey::paired normalization
    // ========================================================================

    #[test]
    fn paired_normalizes_swapped_own_and_mate_to_equal_keys() {
        // Two reads of the same template: one sees (own=A, mate=B), the other
        // sees (own=B, mate=A). After normalization both must yield the same
        // GroupKey, so they group together.
        let a = GroupKey::paired(0, 100, 0, 0, 200, 1, 7, 99, 42);
        let b = GroupKey::paired(0, 200, 1, 0, 100, 0, 7, 99, 42);
        assert_eq!(a, b, "swapped own/mate positions must normalize to equal keys");

        // The lower (ref_id, pos, strand) tuple is placed in slot 1.
        assert_eq!((a.ref_id1, a.pos1, a.strand1), (0, 100, 0));
        assert_eq!((a.ref_id2, a.pos2, a.strand2), (0, 200, 1));
    }

    #[test]
    fn paired_keys_order_by_normalized_position_then_name() {
        // Position 1 < position 2 orders the keys; equal positions fall back to
        // name_hash ordering.
        let lower = GroupKey::paired(0, 100, 0, 0, 200, 0, 0, 0, 1);
        let higher = GroupKey::paired(0, 150, 0, 0, 200, 0, 0, 0, 1);
        assert!(lower < higher, "lower position-1 must sort first");

        let same_pos_low_name = GroupKey::paired(0, 100, 0, 0, 200, 0, 0, 0, 1);
        let same_pos_high_name = GroupKey::paired(0, 100, 0, 0, 200, 0, 0, 0, 2);
        assert!(same_pos_low_name < same_pos_high_name, "equal positions must order by name_hash",);
    }

    // ========================================================================
    // compute_group_key_from_raw: secondary / supplementary fallback
    // ========================================================================

    #[test]
    fn secondary_and_supplementary_records_yield_name_hash_only_key() {
        let lib = LibraryIndex::default();
        let expected_name_hash = LibraryIndex::hash_name(Some(b"rec1"));

        for flag in [flags::SECONDARY, flags::SUPPLEMENTARY] {
            // Give the record a real mapped position so we can prove the
            // position fields are *not* derived for secondary/supplementary.
            let mut b = SamBuilder::new();
            b.ref_id(0)
                .pos(500)
                .flags(flag)
                .read_name(b"rec1")
                .cigar_ops(&[cigar_m(50)])
                .sequence(b"ACGT")
                .qualities(&[30, 30, 30, 30]);
            let rec = b.build();

            let key = compute_group_key_from_raw(rec.as_ref(), &lib, None);

            // Only the name hash is set; every position field stays at default.
            assert_eq!(key, GroupKey { name_hash: expected_name_hash, ..GroupKey::default() });
            assert_eq!(key.ref_id1, GroupKey::UNKNOWN_REF);
            assert_eq!(key.pos1, GroupKey::UNKNOWN_POS);
            assert_eq!(key.strand1, GroupKey::UNKNOWN_STRAND);
        }
    }

    // ========================================================================
    // compute_group_key_from_raw: paired record missing MC falls back to single
    // ========================================================================

    #[test]
    fn paired_without_mc_tag_falls_back_to_single_semantics() {
        let lib = LibraryIndex::default();

        // Paired + mate mapped, but no MC tag: cannot compute the mate's
        // unclipped 5' position, so the key must use single-end semantics
        // (mate fields left at the UNKNOWN sentinels).
        let mut b = SamBuilder::new();
        b.ref_id(0)
            .pos(1000)
            .flags(flags::PAIRED)
            .mate_ref_id(0)
            .mate_pos(1200)
            .read_name(b"pair_no_mc")
            .cigar_ops(&[cigar_m(50)])
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30]);
        let rec = b.build();

        let key = compute_group_key_from_raw(rec.as_ref(), &lib, None);

        // Mate fields fall back to the single-end sentinels.
        assert_eq!(key.ref_id2, GroupKey::UNKNOWN_REF);
        assert_eq!(key.pos2, GroupKey::UNKNOWN_POS);
        assert_eq!(key.strand2, GroupKey::UNKNOWN_STRAND);

        // The own position is still populated; the key matches the single-end
        // key built from the same fields.
        let expected = GroupKey::single(
            key.ref_id1,
            key.pos1,
            key.strand1,
            key.library_idx,
            key.cell_hash,
            LibraryIndex::hash_name(Some(b"pair_no_mc")),
        );
        assert_eq!(key, expected);
    }

    // ========================================================================
    // compute_group_key_and_umi_from_raw: in-line UMI capture parity
    // ========================================================================

    /// The record-relative UMI position a standalone `find_string_tag_position`
    /// scan would produce — i.e. exactly what the old two-pass
    /// `cache_umi_position` cached — for cross-checking the folded capture.
    fn standalone_umi_position(raw: &[u8], tag: [u8; 2]) -> Option<(u32, u16)> {
        let aux_offset = fgumi_raw_bam::aux_data_offset_from_record(raw)?;
        let aux = raw.get(aux_offset..)?;
        let (rel_off, len) = fgumi_raw_bam::find_string_tag_position(aux, tag)?;
        Some((u32::try_from(aux_offset).ok()? + rel_off, len))
    }

    #[rstest]
    // Primary paired read with RX + MC: the common consensus/group/dedup shape.
    #[case::paired_with_rx(
        flags::PAIRED,
        true,
        Some(&b"ACGTACGT"[..]),
    )]
    // Primary single-end read carrying RX.
    #[case::single_with_rx(0, false, Some(&b"TTGGCCAA"[..]))]
    // Primary read with no RX tag: capture must be `None` and fall through to
    // the standalone scan (which also finds nothing).
    #[case::primary_without_rx(flags::PAIRED, true, None)]
    fn compute_group_key_and_umi_captures_umi_matching_standalone_scan(
        #[case] flag: u16,
        #[case] add_mc: bool,
        #[case] umi: Option<&[u8]>,
    ) {
        let lib = LibraryIndex::default();
        let mut b = SamBuilder::new();
        b.ref_id(0)
            .pos(1000)
            .flags(flag)
            .read_name(b"rec")
            .cigar_ops(&[cigar_m(50)])
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30]);
        if flag & flags::PAIRED != 0 {
            b.mate_ref_id(0).mate_pos(1200);
        }
        // RG is present in every real grouping input; include it so the scan
        // has to walk past a leading tag before reaching RX.
        b.add_string_tag(SamTag::RG, b"RG1");
        if add_mc {
            b.add_string_tag(SamTag::MC, b"50M");
        }
        if let Some(u) = umi {
            b.add_string_tag(SamTag::RX, u);
        }
        let rec = b.build();

        let (key, umi_pos) =
            compute_group_key_and_umi_from_raw(rec.as_ref(), &lib, None, Some(*SamTag::RX));

        // The key is identical to the non-UMI-capturing entry point.
        assert_eq!(key, compute_group_key_from_raw(rec.as_ref(), &lib, None));
        // The captured position is byte-for-byte what the old standalone scan
        // (`cache_umi_position` → `find_string_tag_position`) would have cached.
        assert_eq!(umi_pos, standalone_umi_position(rec.as_ref(), *SamTag::RX));
        // And, when a UMI is present, it slices back to the original bytes.
        if let Some(expected) = umi {
            let (off, len) = umi_pos.expect("UMI captured");
            let (off, len) = (off as usize, len as usize);
            assert_eq!(&rec.as_ref()[off..off + len], expected);
        } else {
            assert_eq!(umi_pos, None);
        }
    }

    #[test]
    fn compute_group_key_and_umi_returns_none_when_umi_tag_not_requested() {
        let lib = LibraryIndex::default();
        let mut b = SamBuilder::new();
        b.ref_id(0)
            .pos(1000)
            .flags(0)
            .read_name(b"rec")
            .cigar_ops(&[cigar_m(50)])
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30]);
        b.add_string_tag(SamTag::RX, b"ACGTACGT");
        let rec = b.build();

        // `umi_tag = None` must never capture a position even though RX is present.
        let (_key, umi_pos) = compute_group_key_and_umi_from_raw(rec.as_ref(), &lib, None, None);
        assert_eq!(umi_pos, None);
    }

    /// The `tc`-stamped secondary/supplementary branch also folds the UMI capture
    /// into its aux scan (grouping.rs, the secondary branch), a path the primary
    /// `#[case]`s above do not exercise. Assert the captured position matches the
    /// standalone scan and slices back to the RX bytes.
    #[rstest]
    #[case::secondary(flags::PAIRED | flags::SECONDARY)]
    #[case::supplementary(flags::PAIRED | flags::SUPPLEMENTARY)]
    fn compute_group_key_and_umi_captures_umi_on_the_tc_secondary_path(#[case] flag: u16) {
        let lib = library_index_with_two_groups();
        let mut b = SamBuilder::new();
        b.ref_id(0)
            .pos(5000)
            .flags(flag)
            .read_name(b"sec_with_tc")
            .cigar_ops(&[cigar_m(50)])
            .sequence(b"ACGT")
            .qualities(&[30, 30, 30, 30])
            .add_string_tag(SamTag::RG, b"RG1")
            .add_array_i32(SamTag::TC, &[0, 1001, 0, 0, 1249, 1])
            .add_string_tag(SamTag::RX, b"GGTTCCAA");
        let rec = b.build();

        let (key, umi_pos) =
            compute_group_key_and_umi_from_raw(rec.as_ref(), &lib, None, Some(*SamTag::RX));

        // Key is unchanged from the non-UMI-capturing entry point (still keyed on tc).
        assert_eq!(key, compute_group_key_from_raw(rec.as_ref(), &lib, None));
        // The folded capture equals the standalone scan and slices to the RX bytes.
        assert_eq!(umi_pos, standalone_umi_position(rec.as_ref(), *SamTag::RX));
        let (off, len) = umi_pos.expect("UMI captured on the tc-secondary path");
        let (off, len) = (off as usize, len as usize);
        assert_eq!(&rec.as_ref()[off..off + len], b"GGTTCCAA");
    }

    // ========================================================================
    // name_hash_key parity with compute_group_key_from_raw
    // ========================================================================

    #[test]
    fn name_hash_key_matches_compute_group_key_name_hash() {
        let lib = LibraryIndex::default();

        // A representative spread of records: mapped single-end, secondary,
        // paired-without-MC, and an empty-named record (exercises the
        // hash_name(None) branch shared by both functions).
        let mut single = SamBuilder::new();
        single
            .ref_id(0)
            .pos(10)
            .read_name(b"alpha")
            .cigar_ops(&[cigar_m(8)])
            .sequence(b"ACGTACGT")
            .qualities(&[30; 8]);

        let mut secondary = SamBuilder::new();
        secondary.ref_id(0).pos(20).flags(flags::SECONDARY).read_name(b"beta");

        let mut paired = SamBuilder::new();
        paired
            .ref_id(0)
            .pos(30)
            .flags(flags::PAIRED)
            .mate_ref_id(0)
            .mate_pos(60)
            .read_name(b"gamma");

        let mut empty_name = SamBuilder::new();
        empty_name.ref_id(0).pos(40).read_name(b"");

        for mut builder in [single, secondary, paired, empty_name] {
            let rec = builder.build();
            let raw = rec.as_ref();
            assert_eq!(
                name_hash_key(raw).name_hash,
                compute_group_key_from_raw(raw, &lib, None).name_hash,
                "name_hash parity mismatch",
            );
        }
    }

    // ========================================================================
    // key_for_mode dispatch
    // ========================================================================

    /// `key_for_mode` maps each [`KeyMode`] to the matching key function, and the
    /// `None` arm returns `GroupKey::default()` without inspecting `raw` — the
    /// documented contract that lets the fastq encode decode with no key at all.
    #[test]
    fn key_for_mode_dispatches_each_level() {
        let lib = LibraryIndex::default();

        // `None` must not touch `raw`: an empty (unvalidated) slice that would
        // panic in `name_hash_key`/`compute_group_key_from_raw` is accepted here.
        assert_eq!(
            key_for_mode(KeyMode::None, &[], &lib, None),
            GroupKey::default(),
            "KeyMode::None must return the default key without reading raw",
        );

        let mut b = SamBuilder::new();
        b.ref_id(0)
            .pos(10)
            .read_name(b"delta")
            .cigar_ops(&[cigar_m(4)])
            .sequence(b"ACGT")
            .qualities(&[30; 4]);
        let rec = b.build();
        let raw = rec.as_ref();

        assert_eq!(
            key_for_mode(KeyMode::NameHashOnly, raw, &lib, None),
            name_hash_key(raw),
            "KeyMode::NameHashOnly must match name_hash_key",
        );
        assert_eq!(
            key_for_mode(KeyMode::Full, raw, &lib, None),
            compute_group_key_from_raw(raw, &lib, None),
            "KeyMode::Full must match compute_group_key_from_raw",
        );
    }

    // ========================================================================
    // DecodedRecord cached-UMI invalidation (S6-002)
    // ========================================================================

    /// Build a `DecodedRecord` carrying an RX UMI tag with its value position
    /// cached, exactly as the Decode step's `cache_umi_position` would.
    fn decoded_with_cached_umi(umi: &[u8]) -> DecodedRecord {
        use fgumi_raw_bam::SamTag;
        let mut b = SamBuilder::new();
        b.read_name(b"read1").sequence(b"ACGT").qualities(&[30; 4]);
        b.add_string_tag(SamTag::RX, umi);
        let rec = b.build();
        let bytes = rec.as_ref();
        let aux_offset =
            fgumi_raw_bam::aux_data_offset_from_record(bytes).expect("aux offset present");
        let aux = &bytes[aux_offset..];
        let (off_in_aux, len) =
            fgumi_raw_bam::find_string_tag_position(aux, *SamTag::RX).expect("RX tag present");
        let offset = u32::try_from(aux_offset).expect("aux offset fits u32") + off_in_aux;
        let mut d = DecodedRecord::from_raw_bytes(rec, GroupKey::default());
        d.set_cached_umi(offset, len);
        d
    }

    #[test]
    fn cached_umi_returns_value_then_raw_bytes_mut_invalidates() {
        // The cache resolves to the UMI value up front; handing out mutable bytes
        // must reset the cache to the sentinel so a later read re-scans rather
        // than slicing stale-but-in-range bytes (see the cached-UMI invariant).
        let mut decoded = decoded_with_cached_umi(b"ACGTACGT");
        assert_eq!(decoded.cached_umi(), Some(b"ACGTACGT".as_ref()), "cache populated up front");
        assert_ne!(decoded.cached_umi_position().0, DecodedRecord::UMI_OFFSET_UNCACHED);

        let _ = decoded.raw_bytes_mut();

        assert_eq!(
            decoded.cached_umi_position().0,
            DecodedRecord::UMI_OFFSET_UNCACHED,
            "raw_bytes_mut resets the cache to the uncached sentinel",
        );
        assert!(decoded.cached_umi().is_none(), "cached_umi returns None after mutation");
    }
}
