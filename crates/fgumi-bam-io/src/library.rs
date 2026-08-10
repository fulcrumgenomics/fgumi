//! Library lookup tables built from SAM `@RG` headers.
//!
//! [`LibraryLookup`] maps read-group IDs to library-name strings;
//! [`LibraryIndex`] is the hash-based hot-path variant returning numeric
//! library indices (used by `grouping::compute_group_key_from_raw`).
//!
//! Relocated here from the main crate's `read_info` module so the grouping
//! domain types (which depend on `LibraryIndex`) can live in this crate
//! alongside [`fgumi_raw_bam::RawRecord`]-adjacent helpers and `MemoryEstimate`.

use std::collections::HashMap;
use std::sync::Arc;

use bstr::ByteSlice;
use noodles::sam::header::Header;
use noodles::sam::header::record::value::map::read_group::tag as rg_tag;

/// A lookup table mapping read group IDs to library names.
///
/// This is built from the SAM header's `@RG` lines and used to resolve the library
/// name (`LB` field) from a record's `RG` tag. This matches fgbio's behavior where
/// grouping uses the library name, not the read group ID.
///
/// Uses `Arc<str>` for library names to avoid cloning strings for every read.
///
/// # Note: `LibraryLookup` vs `LibraryIndex`
///
/// Both `LibraryLookup` and [`LibraryIndex`] exist for different use cases:
/// - `LibraryLookup`: String-based lookup returning library names. Used where
///   the actual library-name string is needed (e.g. the main crate's
///   `ReadInfo::from`).
/// - [`LibraryIndex`]: Hash-based lookup returning numeric indices. Used by
///   [`compute_group_key_from_raw`](crate::compute_group_key_from_raw) in the
///   hot path where only equality comparison matters, avoiding string
///   allocations.
pub type LibraryLookup = Arc<HashMap<String, Arc<str>>>;

/// Shared "unknown" library string to avoid repeated allocations.
static UNKNOWN_LIBRARY: std::sync::LazyLock<Arc<str>> =
    std::sync::LazyLock::new(|| Arc::from("unknown"));

/// Returns the shared "unknown" library name (`Arc<str>` of `"unknown"`),
/// used as the fallback when a read group has no `LB` field.
#[must_use]
pub fn unknown_library() -> Arc<str> {
    Arc::clone(&UNKNOWN_LIBRARY)
}

/// Builds a library lookup table from a SAM header.
///
/// Iterates through all `@RG` lines in the header and creates a mapping from
/// read group ID to library name. If a read group has no `LB` field, it maps
/// to "unknown" (matching fgbio's behavior).
///
/// # Arguments
///
/// * `header` - The SAM header containing `@RG` lines
///
/// # Returns
///
/// An `Arc<HashMap>` mapping read group IDs to library names
#[must_use]
pub fn build_library_lookup(header: &Header) -> LibraryLookup {
    let mut lookup = HashMap::new();

    for (id, rg) in header.read_groups() {
        // Get the LB field from the read group's other_fields
        let library: Arc<str> = rg
            .other_fields()
            .get(&rg_tag::LIBRARY)
            .map_or_else(|| Arc::clone(&UNKNOWN_LIBRARY), |s| Arc::from(s.to_string()));
        lookup.insert(id.to_string(), library);
    }

    Arc::new(lookup)
}

// ============================================================================
// LibraryIndex - Fast RG to library index mapping for GroupKey computation
// ============================================================================

/// Fast lookup from `RG` tag value to library index.
///
/// This provides `O(1)` library lookup during Decode using string hashing,
/// replacing the `O(n)` string comparison in the original `LibraryLookup`.
#[derive(Debug, Clone)]
pub struct LibraryIndex {
    /// Map from `RG` string hash to library index.
    lookup: ahash::AHashMap<u64, u16>,
    /// Library names for each index (for output/debugging).
    names: Vec<Arc<str>>,
    /// Unknown library index (always 0).
    unknown_idx: u16,
}

impl LibraryIndex {
    /// Build a library index from a SAM header.
    ///
    /// Each unique library name gets a sequential index starting from 0.
    /// Index 0 is reserved for "unknown" library.
    ///
    /// # Panics
    ///
    /// Panics if the header contains more than 65,535 distinct libraries.
    #[must_use]
    pub fn from_header(header: &Header) -> Self {
        let mut lookup = ahash::AHashMap::new();
        let mut names = vec![Arc::clone(&UNKNOWN_LIBRARY)]; // Index 0 = unknown
        let mut library_to_idx: ahash::AHashMap<Arc<str>, u16> = ahash::AHashMap::new();
        library_to_idx.insert(Arc::clone(&UNKNOWN_LIBRARY), 0);

        for (id, rg) in header.read_groups() {
            // Get library name from LB field
            let library: Arc<str> = rg
                .other_fields()
                .get(&rg_tag::LIBRARY)
                .map_or_else(|| Arc::clone(&UNKNOWN_LIBRARY), |s| Arc::from(s.to_string()));

            // Get or create library index
            let lib_idx = *library_to_idx.entry(library.clone()).or_insert_with(|| {
                let idx: u16 =
                    names.len().try_into().expect("too many distinct libraries for u16 index");
                names.push(library);
                idx
            });

            // Hash the RG string and map to library index
            let rg_hash = Self::hash_rg(id.as_bytes());
            lookup.insert(rg_hash, lib_idx);
        }

        Self { lookup, names, unknown_idx: 0 }
    }

    /// Get the library index for a read group hash.
    ///
    /// Returns 0 (unknown) if the `RG` hash is not found.
    #[must_use]
    pub fn get(&self, rg_hash: u64) -> u16 {
        *self.lookup.get(&rg_hash).unwrap_or(&self.unknown_idx)
    }

    /// Get the library name for an index.
    #[must_use]
    pub fn library_name(&self, idx: u16) -> &Arc<str> {
        self.names.get(idx as usize).unwrap_or(&self.names[0])
    }

    /// Hash a byte slice using `AHash`. Returns 0 for `None`.
    ///
    /// This is the single hashing implementation used by all `hash_*` methods.
    #[must_use]
    pub fn hash_bytes(bytes: Option<&[u8]>) -> u64 {
        use ahash::AHasher;
        use std::hash::{Hash, Hasher};
        match bytes {
            Some(b) => {
                let mut hasher = AHasher::default();
                b.hash(&mut hasher);
                hasher.finish()
            }
            None => 0,
        }
    }

    /// Hash an `RG` tag value for lookup.
    #[must_use]
    pub fn hash_rg(rg_bytes: &[u8]) -> u64 {
        Self::hash_bytes(Some(rg_bytes))
    }

    /// Hash a cell barcode for `GroupKey`.
    #[must_use]
    pub fn hash_cell_barcode(cell_bytes: Option<&[u8]>) -> u64 {
        Self::hash_bytes(cell_bytes)
    }

    /// Hash a read name for `GroupKey`.
    #[must_use]
    pub fn hash_name(name_bytes: Option<&[u8]>) -> u64 {
        Self::hash_bytes(name_bytes)
    }
}

impl Default for LibraryIndex {
    fn default() -> Self {
        Self {
            lookup: ahash::AHashMap::new(),
            names: vec![Arc::clone(&UNKNOWN_LIBRARY)],
            unknown_idx: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::header::record::value::Map;
    use noodles::sam::header::record::value::map::ReadGroup;
    use noodles::sam::header::record::value::map::read_group::tag as rg_tag;

    /// Header with three read groups: two sharing a library, one with no `LB`.
    fn header_with_read_groups() -> Header {
        let mut header = Header::builder();
        // RG1 and RG3 deliberately share a library so the dedup path is exercised.
        for (id, library) in [("RG1", "libA"), ("RG2", "libB"), ("RG3", "libA")] {
            let rg = Map::<ReadGroup>::builder()
                .insert(rg_tag::LIBRARY, String::from(library))
                .build()
                .expect("read group builds");
            header = header.add_read_group(bstr::BString::from(id), rg);
        }
        // A fourth group with no LB at all, which must fall back to "unknown".
        header = header.add_read_group(
            bstr::BString::from("RG4"),
            Map::<ReadGroup>::builder().build().expect("read group builds"),
        );
        header.build()
    }

    /// `library_idx` is the field `GroupKey` equality depends on, so read groups
    /// sharing a library MUST collapse to one index and distinct libraries MUST
    /// NOT — the first would merge unrelated UMI groups, the second would split
    /// one group in half.
    #[test]
    fn from_header_assigns_one_index_per_distinct_library() {
        let index = LibraryIndex::from_header(&header_with_read_groups());

        let idx1 = index.get(LibraryIndex::hash_rg(b"RG1"));
        let idx2 = index.get(LibraryIndex::hash_rg(b"RG2"));
        let idx3 = index.get(LibraryIndex::hash_rg(b"RG3"));

        assert_eq!(idx1, idx3, "read groups sharing a library share an index");
        assert_ne!(idx1, idx2, "distinct libraries get distinct indices");
        assert_ne!(idx1, 0, "a resolved library is never the unknown index");
        assert_ne!(idx2, 0, "a resolved library is never the unknown index");

        assert_eq!(index.library_name(idx1).as_ref(), "libA");
        assert_eq!(index.library_name(idx2).as_ref(), "libB");
    }

    /// A read group with no `LB` maps to the shared "unknown" library at index 0,
    /// which is also the miss value — so an absent `LB` and an absent read group
    /// group together rather than each forming a singleton.
    #[test]
    fn from_header_maps_missing_library_and_unknown_read_group_to_index_zero() {
        let index = LibraryIndex::from_header(&header_with_read_groups());

        assert_eq!(index.get(LibraryIndex::hash_rg(b"RG4")), 0, "no LB field → unknown");
        assert_eq!(index.get(LibraryIndex::hash_rg(b"ABSENT")), 0, "unseen RG → unknown");
        assert_eq!(index.library_name(0).as_ref(), "unknown");
    }

    /// Out-of-range indices fall back to "unknown" rather than panicking.
    #[test]
    fn library_name_saturates_to_unknown_for_an_out_of_range_index() {
        let index = LibraryIndex::from_header(&header_with_read_groups());
        assert_eq!(index.library_name(u16::MAX).as_ref(), "unknown");
    }

    /// The string-keyed lookup carries the same LB resolution as the hashed one,
    /// including the "unknown" fallback.
    #[test]
    fn build_library_lookup_maps_each_read_group_to_its_library() {
        let lookup = build_library_lookup(&header_with_read_groups());

        assert_eq!(lookup.len(), 4);
        assert_eq!(lookup.get("RG1").expect("RG1 present").as_ref(), "libA");
        assert_eq!(lookup.get("RG2").expect("RG2 present").as_ref(), "libB");
        assert_eq!(lookup.get("RG3").expect("RG3 present").as_ref(), "libA");
        assert_eq!(lookup.get("RG4").expect("RG4 present").as_ref(), "unknown");
    }

    /// An empty header yields an index that resolves everything to unknown.
    #[test]
    fn from_header_with_no_read_groups_resolves_everything_to_unknown() {
        let index = LibraryIndex::from_header(&Header::default());
        assert_eq!(index.get(LibraryIndex::hash_rg(b"RG1")), 0);
        assert!(build_library_lookup(&Header::default()).is_empty());
    }

    /// `hash_bytes(None)` is the documented zero sentinel, and the typed wrappers
    /// agree with it — `GroupKey`'s `cell_hash`/`name_hash` rely on that.
    #[test]
    fn hash_helpers_agree_and_treat_none_as_zero() {
        assert_eq!(LibraryIndex::hash_bytes(None), 0);
        assert_eq!(LibraryIndex::hash_cell_barcode(None), 0);
        assert_eq!(LibraryIndex::hash_name(None), 0);
        assert_eq!(
            LibraryIndex::hash_name(Some(b"read1")),
            LibraryIndex::hash_bytes(Some(b"read1"))
        );
        assert_ne!(LibraryIndex::hash_name(Some(b"read1")), 0);
    }
}
