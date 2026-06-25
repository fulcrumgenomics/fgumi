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
/// This is built from the SAM header's @RG lines and used to resolve the library
/// name (LB field) from a record's RG tag. This matches fgbio's behavior where
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
/// Iterates through all @RG lines in the header and creates a mapping from
/// read group ID to library name. If a read group has no LB field, it maps
/// to "unknown" (matching fgbio's behavior).
///
/// # Arguments
///
/// * `header` - The SAM header containing @RG lines
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

/// Fast lookup from RG tag value to library index.
///
/// This provides O(1) library lookup during Decode using string hashing,
/// replacing the O(n) string comparison in the original `LibraryLookup`.
#[derive(Debug, Clone)]
pub struct LibraryIndex {
    /// Map from RG string hash to library index.
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
    /// Returns 0 (unknown) if the RG hash is not found.
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

    /// Hash an RG tag value for lookup.
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
