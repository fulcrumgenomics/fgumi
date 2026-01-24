//! Consensus sequence data structure.
//!
//! This module provides [`ConsensusSequence`], a wrapper type that enforces
//! invariants for parallel vectors representing consensus base calls.

/// Consensus sequence with per-base annotations.
///
/// This struct maintains parallel vectors for bases, quality scores, depths, and errors.
/// All vectors are guaranteed to have the same length through the API.
///
/// # Layout
///
/// Uses a Structure-of-Arrays layout for cache-friendly access:
/// - Sequential access to all bases is cache-efficient
/// - Sequential access to all qualities is cache-efficient
/// - etc.
///
/// # Example
///
/// ```
/// use fgumi_lib::consensus::ConsensusSequence;
///
/// let mut seq = ConsensusSequence::new();
/// seq.push(b'A', 30, 10, 0);
/// seq.push(b'C', 25, 8, 1);
///
/// assert_eq!(seq.len(), 2);
/// assert_eq!(seq.bases(), &[b'A', b'C']);
/// assert_eq!(seq.quals(), &[30, 25]);
/// ```
#[derive(Debug, Clone, Default)]
pub struct ConsensusSequence {
    /// Consensus sequence bases
    bases: Vec<u8>,
    /// Consensus quality scores (Phred-scaled)
    quals: Vec<u8>,
    /// Per-base read depth
    depths: Vec<u16>,
    /// Per-base error count
    errors: Vec<u16>,
}

impl ConsensusSequence {
    /// Creates a new empty consensus sequence.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates a new consensus sequence with the given capacity.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            bases: Vec::with_capacity(capacity),
            quals: Vec::with_capacity(capacity),
            depths: Vec::with_capacity(capacity),
            errors: Vec::with_capacity(capacity),
        }
    }

    /// Creates a consensus sequence from parallel vectors.
    ///
    /// # Panics
    ///
    /// Panics if the vectors have different lengths.
    #[must_use]
    pub fn from_vecs(bases: Vec<u8>, quals: Vec<u8>, depths: Vec<u16>, errors: Vec<u16>) -> Self {
        let len = bases.len();
        assert!(
            quals.len() == len && depths.len() == len && errors.len() == len,
            "All vectors must have the same length: bases={}, quals={}, depths={}, errors={}",
            len,
            quals.len(),
            depths.len(),
            errors.len()
        );
        Self { bases, quals, depths, errors }
    }

    /// Returns the length of the consensus sequence.
    #[must_use]
    pub fn len(&self) -> usize {
        self.bases.len()
    }

    /// Returns true if the consensus sequence is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }

    /// Appends a single base with its annotations.
    pub fn push(&mut self, base: u8, qual: u8, depth: u16, error: u16) {
        self.bases.push(base);
        self.quals.push(qual);
        self.depths.push(depth);
        self.errors.push(error);
    }

    /// Extends this sequence with another.
    pub fn extend(&mut self, other: &ConsensusSequence) {
        self.bases.extend_from_slice(&other.bases);
        self.quals.extend_from_slice(&other.quals);
        self.depths.extend_from_slice(&other.depths);
        self.errors.extend_from_slice(&other.errors);
    }

    /// Clears all vectors.
    pub fn clear(&mut self) {
        self.bases.clear();
        self.quals.clear();
        self.depths.clear();
        self.errors.clear();
    }

    /// Returns the consensus bases.
    #[must_use]
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }

    /// Returns mutable access to the consensus bases.
    ///
    /// # Note
    ///
    /// The caller should not change the length (slices prevent this anyway).
    #[must_use]
    pub fn bases_mut(&mut self) -> &mut [u8] {
        &mut self.bases
    }

    /// Returns the quality scores.
    #[must_use]
    pub fn quals(&self) -> &[u8] {
        &self.quals
    }

    /// Returns mutable access to the quality scores.
    ///
    /// # Note
    ///
    /// The caller should not change the length (slices prevent this anyway).
    #[must_use]
    pub fn quals_mut(&mut self) -> &mut [u8] {
        &mut self.quals
    }

    /// Returns the per-base depths.
    #[must_use]
    pub fn depths(&self) -> &[u16] {
        &self.depths
    }

    /// Returns mutable access to the per-base depths.
    ///
    /// # Note
    ///
    /// The caller should not change the length (slices prevent this anyway).
    #[must_use]
    pub fn depths_mut(&mut self) -> &mut [u16] {
        &mut self.depths
    }

    /// Returns the per-base error counts.
    #[must_use]
    pub fn errors(&self) -> &[u16] {
        &self.errors
    }

    /// Returns mutable access to the per-base error counts.
    ///
    /// # Note
    ///
    /// The caller should not change the length (slices prevent this anyway).
    #[must_use]
    pub fn errors_mut(&mut self) -> &mut [u16] {
        &mut self.errors
    }

    /// Returns the maximum depth across all positions.
    #[must_use]
    pub fn max_depth(&self) -> u16 {
        self.depths.iter().copied().max().unwrap_or(0)
    }

    /// Returns the minimum depth across all positions.
    #[must_use]
    pub fn min_depth(&self) -> u16 {
        self.depths.iter().copied().min().unwrap_or(0)
    }

    /// Returns the error rate (total errors / total depth).
    #[must_use]
    pub fn error_rate(&self) -> f32 {
        let total_depth: u32 = self.depths.iter().map(|&d| u32::from(d)).sum();
        if total_depth == 0 {
            return 0.0;
        }
        let total_errors: u32 = self.errors.iter().map(|&e| u32::from(e)).sum();
        total_errors as f32 / total_depth as f32
    }

    /// Decomposes this sequence into its component vectors.
    ///
    /// This is useful when the caller needs ownership of the vectors.
    #[must_use]
    pub fn into_vecs(self) -> (Vec<u8>, Vec<u8>, Vec<u16>, Vec<u16>) {
        (self.bases, self.quals, self.depths, self.errors)
    }

    /// Pads the sequence to a new length by adding values to the left or right.
    ///
    /// # Arguments
    /// * `new_length` - Target length (must be >= current length)
    /// * `left` - If true, pad on the left side; otherwise pad on the right
    /// * `base` - Base to use for padding (default: 'n')
    /// * `qual` - Quality to use for padding (default: 2)
    ///
    /// # Panics
    /// Panics if `new_length` < current length
    #[must_use]
    pub fn padded(&self, new_length: usize, left: bool, base: u8, qual: u8) -> Self {
        let current_len = self.bases.len();
        assert!(
            new_length >= current_len,
            "new_length ({new_length}) must be >= current length ({current_len})"
        );

        if new_length == current_len {
            return self.clone();
        }

        let pad_len = new_length - current_len;

        // Pre-allocate with exact capacity (avoids intermediate allocations)
        let mut bases = Vec::with_capacity(new_length);
        let mut quals = Vec::with_capacity(new_length);
        let mut depths = Vec::with_capacity(new_length);
        let mut errors = Vec::with_capacity(new_length);

        if left {
            // Padding on left: [padding][original]
            bases.resize(pad_len, base);
            quals.resize(pad_len, qual);
            depths.resize(pad_len, 0);
            errors.resize(pad_len, 0);

            bases.extend_from_slice(&self.bases);
            quals.extend_from_slice(&self.quals);
            depths.extend_from_slice(&self.depths);
            errors.extend_from_slice(&self.errors);
        } else {
            // Padding on right: [original][padding]
            bases.extend_from_slice(&self.bases);
            quals.extend_from_slice(&self.quals);
            depths.extend_from_slice(&self.depths);
            errors.extend_from_slice(&self.errors);

            bases.resize(new_length, base);
            quals.resize(new_length, qual);
            depths.resize(new_length, 0);
            errors.resize(new_length, 0);
        }

        Self { bases, quals, depths, errors }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_is_empty() {
        let seq = ConsensusSequence::new();
        assert!(seq.is_empty());
        assert_eq!(seq.len(), 0);
    }

    #[test]
    fn test_push_and_accessors() {
        let mut seq = ConsensusSequence::new();
        seq.push(b'A', 30, 10, 0);
        seq.push(b'C', 25, 8, 1);
        seq.push(b'G', 20, 5, 2);

        assert_eq!(seq.len(), 3);
        assert!(!seq.is_empty());
        assert_eq!(seq.bases(), b"ACG");
        assert_eq!(seq.quals(), &[30, 25, 20]);
        assert_eq!(seq.depths(), &[10, 8, 5]);
        assert_eq!(seq.errors(), &[0, 1, 2]);
    }

    #[test]
    fn test_from_vecs() {
        let seq =
            ConsensusSequence::from_vecs(vec![b'A', b'T'], vec![30, 25], vec![10, 8], vec![0, 1]);

        assert_eq!(seq.len(), 2);
        assert_eq!(seq.bases(), b"AT");
    }

    #[test]
    #[should_panic(expected = "All vectors must have the same length")]
    fn test_from_vecs_mismatched_lengths() {
        let _ = ConsensusSequence::from_vecs(
            vec![b'A', b'T'],
            vec![30], // Wrong length
            vec![10, 8],
            vec![0, 1],
        );
    }

    #[test]
    fn test_extend() {
        let mut seq1 = ConsensusSequence::from_vecs(vec![b'A'], vec![30], vec![10], vec![0]);
        let seq2 =
            ConsensusSequence::from_vecs(vec![b'T', b'C'], vec![25, 20], vec![8, 5], vec![1, 0]);

        seq1.extend(&seq2);
        assert_eq!(seq1.len(), 3);
        assert_eq!(seq1.bases(), b"ATC");
    }

    #[test]
    fn test_max_min_depth() {
        let seq = ConsensusSequence::from_vecs(
            vec![b'A', b'T', b'C'],
            vec![30, 25, 20],
            vec![10, 5, 15],
            vec![0, 0, 0],
        );

        assert_eq!(seq.max_depth(), 15);
        assert_eq!(seq.min_depth(), 5);
    }

    #[test]
    fn test_error_rate() {
        let seq = ConsensusSequence::from_vecs(
            vec![b'A', b'T', b'C', b'G'],
            vec![30, 30, 30, 30],
            vec![10, 10, 10, 10], // Total depth = 40
            vec![0, 1, 1, 2],     // Total errors = 4
        );

        // Error rate = 4/40 = 0.1
        assert!((seq.error_rate() - 0.1).abs() < 0.001);
    }

    #[test]
    fn test_error_rate_zero_depth() {
        let seq = ConsensusSequence::new();
        assert!((seq.error_rate() - 0.0).abs() < f32::EPSILON);
    }

    #[test]
    fn test_into_vecs() {
        let seq =
            ConsensusSequence::from_vecs(vec![b'A', b'T'], vec![30, 25], vec![10, 8], vec![0, 1]);

        let (bases, quals, depths, errors) = seq.into_vecs();
        assert_eq!(bases, vec![b'A', b'T']);
        assert_eq!(quals, vec![30, 25]);
        assert_eq!(depths, vec![10, 8]);
        assert_eq!(errors, vec![0, 1]);
    }

    #[test]
    fn test_padded_right() {
        let seq =
            ConsensusSequence::from_vecs(vec![b'A', b'C'], vec![30, 25], vec![10, 8], vec![0, 1]);

        let padded = seq.padded(4, false, b'N', 2);
        assert_eq!(padded.len(), 4);
        assert_eq!(padded.bases(), b"ACNN");
        assert_eq!(padded.quals(), &[30, 25, 2, 2]);
        assert_eq!(padded.depths(), &[10, 8, 0, 0]);
    }

    #[test]
    fn test_padded_left() {
        let seq =
            ConsensusSequence::from_vecs(vec![b'A', b'C'], vec![30, 25], vec![10, 8], vec![0, 1]);

        let padded = seq.padded(4, true, b'N', 2);
        assert_eq!(padded.len(), 4);
        assert_eq!(padded.bases(), b"NNAC");
        assert_eq!(padded.quals(), &[2, 2, 30, 25]);
    }

    #[test]
    fn test_padded_same_length() {
        let seq =
            ConsensusSequence::from_vecs(vec![b'A', b'C'], vec![30, 25], vec![10, 8], vec![0, 1]);

        let padded = seq.padded(2, false, b'N', 2);
        assert_eq!(padded.len(), 2);
        assert_eq!(padded.bases(), seq.bases());
    }

    #[test]
    fn test_clear() {
        let mut seq =
            ConsensusSequence::from_vecs(vec![b'A', b'C'], vec![30, 25], vec![10, 8], vec![0, 1]);

        seq.clear();
        assert!(seq.is_empty());
        assert_eq!(seq.len(), 0);
    }

    #[test]
    fn test_with_capacity() {
        let seq = ConsensusSequence::with_capacity(100);
        assert!(seq.is_empty());
        // Capacity is an implementation detail, but we can verify the struct is usable
    }
}
