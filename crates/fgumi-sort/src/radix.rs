//! Radix sort for coordinate-based keys.
//!
//! This module provides O(n) radix sort for coordinate sorting, which is
//! faster than comparison-based O(n log n) sorts for large arrays.
//!
//! # Algorithm
//!
//! Uses LSD (Least Significant Digit) radix sort with 8-bit radix (256 buckets).
//! Adaptive: only sorts the bytes actually needed based on max values in dataset.

use std::cmp::Ordering;

/// Threshold below which we use insertion sort instead of radix sort.
const RADIX_THRESHOLD: usize = 256;

// ============================================================================
// Samtools-style Adaptive Radix Sort
// ============================================================================

/// Adaptive radix sort for coordinate keys, following samtools approach.
///
/// Key optimizations from samtools:
/// 1. Only sort the bytes actually needed (based on max tid/pos in dataset)
/// 2. Pack pos+reverse+tid into little-endian bytes for LSD radix sort
/// 3. Unmapped reads (tid=-1) sort to end
///
/// # Arguments
/// * `entries` - Slice of (`packed_key`, `record_index`) pairs to sort in place
/// * `nref` - Number of reference sequences (for mapping unmapped tid)
#[allow(clippy::uninit_vec, unsafe_code)]
pub fn radix_sort_coordinate_adaptive<T: Clone>(
    entries: &mut [(u64, T)],
    max_tid: u32,
    max_pos: u64,
) {
    let n = entries.len();
    if n < RADIX_THRESHOLD {
        insertion_sort_by_key(entries, |(k, _)| *k);
        return;
    }

    // Calculate bytes needed for pos and tid (like samtools)
    let pos_bytes = bytes_needed_u64(max_pos);
    let tid_bytes = bytes_needed_u32(max_tid);
    let total_bytes = pos_bytes + tid_bytes;

    if total_bytes == 0 {
        return; // All same value, already sorted
    }

    // Allocate auxiliary buffer
    let mut aux: Vec<(u64, T)> = Vec::with_capacity(n);
    unsafe {
        aux.set_len(n);
    }

    let mut src = entries as *mut [(u64, T)];
    let mut dst = aux.as_mut_slice() as *mut [(u64, T)];

    // LSD radix sort - byte by byte from least significant
    for byte_idx in 0..total_bytes {
        let src_slice = unsafe { &*src };
        let dst_slice = unsafe { &mut *dst };

        // Count occurrences of each byte value
        let mut counts = [0usize; 256];
        for (key, _) in src_slice {
            let byte = ((key >> (byte_idx * 8)) & 0xFF) as usize;
            counts[byte] += 1;
        }

        // Convert to cumulative offsets
        let mut total = 0;
        for count in &mut counts {
            let c = *count;
            *count = total;
            total += c;
        }

        // Scatter elements to destination
        for item in src_slice {
            let byte = ((item.0 >> (byte_idx * 8)) & 0xFF) as usize;
            let dest_idx = counts[byte];
            counts[byte] += 1;
            dst_slice[dest_idx] = item.clone();
        }

        // Swap src and dst
        std::mem::swap(&mut src, &mut dst);
    }

    // If odd number of passes, copy back to original buffer
    if total_bytes % 2 == 1 {
        let src_slice = unsafe { &*src };
        entries.clone_from_slice(src_slice);
    }
}

/// Calculate number of bytes needed to represent a u64 value.
#[inline]
pub(crate) fn bytes_needed_u64(val: u64) -> usize {
    if val == 0 {
        return 0;
    }
    ((64 - val.leading_zeros()) as usize).div_ceil(8)
}

/// Calculate number of bytes needed to represent a u32 value.
#[inline]
fn bytes_needed_u32(val: u32) -> usize {
    if val == 0 {
        return 0;
    }
    ((32 - val.leading_zeros()) as usize).div_ceil(8)
}

// ============================================================================
// Fixed-Array Heap (ks_heapadjust style)
// ============================================================================

/// In-place heap adjustment (sift-down), following samtools `ks_heapadjust`.
///
/// This is more efficient than Rust's `BinaryHeap` for the merge phase because:
/// 1. No allocation on push/pop
/// 2. Single sift-down operation per element
/// 3. Works with fixed-size array
///
/// # Arguments
/// * `heap` - The heap array (max-heap by default, use reverse comparator for min-heap)
/// * `i` - Index to sift down from
/// * `n` - Heap size (may be less than array length)
/// * `lt` - Less-than comparator (for max-heap, swap if child > parent)
#[inline]
#[allow(unsafe_code)]
pub fn heap_sift_down<T, F>(heap: &mut [T], mut i: usize, n: usize, lt: &F)
where
    F: Fn(&T, &T) -> bool,
{
    let tmp = unsafe { std::ptr::read(&raw const heap[i]) };

    loop {
        let left = 2 * i + 1;
        if left >= n {
            break;
        }

        // Find larger child
        let right = left + 1;
        let mut child = left;
        if right < n && lt(&heap[left], &heap[right]) {
            child = right;
        }

        // If tmp >= largest child, we're done
        if !lt(&tmp, &heap[child]) {
            break;
        }

        // Move child up
        unsafe {
            std::ptr::copy_nonoverlapping(&raw const heap[child], &raw mut heap[i], 1);
        }
        i = child;
    }

    unsafe {
        std::ptr::write(&raw mut heap[i], tmp);
    }
}

/// Build a heap from an unsorted array (heapify).
#[inline]
pub fn heap_make<T, F>(heap: &mut [T], lt: &F)
where
    F: Fn(&T, &T) -> bool,
{
    let n = heap.len();
    if n <= 1 {
        return;
    }

    // Start from last non-leaf node and sift down
    for i in (0..n / 2).rev() {
        heap_sift_down(heap, i, n, lt);
    }
}

/// Pop the top element and restore heap property.
///
/// Returns the new heap size (n - 1). The popped element is moved to heap[n-1].
#[inline]
pub fn heap_pop<T, F>(heap: &mut [T], n: usize, lt: &F) -> usize
where
    F: Fn(&T, &T) -> bool,
{
    if n == 0 {
        return 0;
    }
    if n == 1 {
        return 0;
    }

    // Swap top with last element
    heap.swap(0, n - 1);

    // Sift down the new top
    let new_n = n - 1;
    if new_n > 0 {
        heap_sift_down(heap, 0, new_n, lt);
    }

    new_n
}

/// Replace the top element and restore heap property.
///
/// More efficient than pop + push when the heap size stays the same.
#[inline]
pub fn heap_replace_top<T, F>(heap: &mut [T], new_value: T, n: usize, lt: &F)
where
    F: Fn(&T, &T) -> bool,
{
    heap[0] = new_value;
    heap_sift_down(heap, 0, n, lt);
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Binary insertion sort for small arrays.
///
/// Uses binary search to find insertion point, reducing comparisons from
/// O(n²) to O(n log n) while maintaining O(n²) moves.
#[inline]
pub fn insertion_sort_by_key<T, K: Ord, F: Fn(&T) -> K>(arr: &mut [T], key_fn: F) {
    for i in 1..arr.len() {
        // Binary search for insertion point
        let key = key_fn(&arr[i]);
        let insert_pos = arr[..i].partition_point(|x| key_fn(x) <= key);

        // Rotate to insert element at correct position
        if insert_pos < i {
            arr[insert_pos..=i].rotate_right(1);
        }
    }
}

/// Binary insertion sort with comparison function.
#[inline]
pub fn binary_insertion_sort<T, F>(arr: &mut [T], compare: F)
where
    F: Fn(&T, &T) -> Ordering,
{
    for i in 1..arr.len() {
        // Binary search for insertion point
        let mut lo = 0;
        let mut hi = i;

        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if compare(&arr[mid], &arr[i]) == Ordering::Greater {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }

        // Rotate to insert element at correct position
        if lo < i {
            arr[lo..=i].rotate_right(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Packed coordinate key for testing radix sort.
    ///
    /// Layout: `[tid:16][pos:32][reverse:1][padding:15]`
    #[derive(Clone, Copy, Eq, PartialEq, Debug)]
    #[repr(transparent)]
    struct PackedCoordinateKey(u64);

    impl PackedCoordinateKey {
        #[allow(clippy::cast_sign_loss)]
        fn new(tid: i32, pos: i32, reverse: bool) -> Self {
            let tid_bits = if tid < 0 { 0xFFFF_u64 } else { (tid as u64) & 0xFFFF };
            let pos_bits = if pos < 0 { 0xFFFF_FFFF_u64 } else { (pos as u64) & 0xFFFF_FFFF };
            let reverse_bit = u64::from(reverse);
            Self((tid_bits << 48) | (pos_bits << 16) | (reverse_bit << 15))
        }
    }

    impl PartialOrd for PackedCoordinateKey {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }

    impl Ord for PackedCoordinateKey {
        fn cmp(&self, other: &Self) -> Ordering {
            self.0.cmp(&other.0)
        }
    }

    /// Pack coordinate fields for radix sorting (samtools style, test-only).
    #[allow(clippy::cast_sign_loss)]
    fn pack_coordinate_for_radix(tid: i32, pos: i32, reverse: bool, nref: u32) -> u64 {
        let tid_val = if tid < 0 { nref } else { tid as u32 };
        let pos_val = (((pos + 1) as u64) << 1) | u64::from(reverse);
        pos_val | (u64::from(tid_val) << 40)
    }

    #[test]
    fn test_packed_coordinate_key() {
        let k1 = PackedCoordinateKey::new(0, 100, false);
        let k2 = PackedCoordinateKey::new(0, 200, false);
        let k3 = PackedCoordinateKey::new(1, 100, false);
        let k4 = PackedCoordinateKey::new(0, 100, true);

        assert!(k1 < k2); // Same tid, pos1 < pos2
        assert!(k1 < k3); // tid1 < tid2
        assert!(k1 < k4); // Same tid+pos, forward < reverse
    }

    #[test]
    fn test_packed_coordinate_key_unmapped() {
        let mapped = PackedCoordinateKey::new(0, 100, false);
        let unmapped = PackedCoordinateKey::new(-1, -1, false);

        assert!(mapped < unmapped); // Unmapped sorts last
    }

    #[test]
    fn test_insertion_sort_by_key_packed() {
        let mut entries: Vec<(u64, i32)> = vec![(5, 50), (3, 30), (8, 80), (1, 10), (4, 40)];
        let mut expected = entries.clone();
        expected.sort_by_key(|(k, _)| *k);

        insertion_sort_by_key(&mut entries, |(k, _)| *k);

        assert_eq!(entries, expected);
    }

    #[test]
    fn test_radix_sort_coordinate_adaptive_large() {
        let mut entries: Vec<(u64, usize)> = (0..1000).rev().map(|i| (i as u64, i)).collect();
        let mut expected = entries.clone();
        expected.sort_by_key(|(k, _)| *k);

        // max_tid=100, max_pos large enough to cover all keys
        radix_sort_coordinate_adaptive(&mut entries, 100, 1000);

        assert_eq!(entries, expected);
    }

    #[test]
    fn test_insertion_sort() {
        let mut arr = vec![5, 3, 8, 1, 4, 2, 7, 6];
        binary_insertion_sort(&mut arr, std::cmp::Ord::cmp);
        assert_eq!(arr, vec![1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn test_insertion_sort_by_key() {
        let mut arr: Vec<(i32, &str)> = vec![(5, "five"), (3, "three"), (8, "eight"), (1, "one")];
        insertion_sort_by_key(&mut arr, |(k, _)| *k);
        assert_eq!(arr[0].0, 1);
        assert_eq!(arr[1].0, 3);
        assert_eq!(arr[2].0, 5);
        assert_eq!(arr[3].0, 8);
    }

    #[test]
    fn test_heap_operations() {
        // Test min-heap (use > as lt for min-heap)
        let lt = |a: &i32, b: &i32| *a > *b;

        let mut heap = vec![5, 3, 8, 1, 4, 2, 7, 6];
        heap_make(&mut heap, &lt);

        // Pop elements should come out in sorted order
        let mut sorted = Vec::new();
        let mut n = heap.len();
        while n > 0 {
            sorted.push(heap[0]);
            n = heap_pop(&mut heap, n, &lt);
        }

        assert_eq!(sorted, vec![1, 2, 3, 4, 5, 6, 7, 8]);
    }

    #[test]
    fn test_pack_coordinate_for_radix() {
        let nref = 100u32;

        // Mapped read
        let k1 = pack_coordinate_for_radix(0, 100, false, nref);
        let k2 = pack_coordinate_for_radix(0, 200, false, nref);
        let k3 = pack_coordinate_for_radix(1, 100, false, nref);
        let k4 = pack_coordinate_for_radix(0, 100, true, nref);

        assert!(k1 < k2); // Same tid, pos1 < pos2
        assert!(k1 < k3); // tid1 < tid2
        assert!(k1 < k4); // Same tid+pos, forward < reverse

        // Unmapped sorts last
        let unmapped = pack_coordinate_for_radix(-1, -1, false, nref);
        assert!(k1 < unmapped);
        assert!(k3 < unmapped);
    }

    #[test]
    fn test_bytes_needed() {
        assert_eq!(bytes_needed_u64(0), 0);
        assert_eq!(bytes_needed_u64(255), 1);
        assert_eq!(bytes_needed_u64(256), 2);
        assert_eq!(bytes_needed_u64(65535), 2);
        assert_eq!(bytes_needed_u64(65536), 3);
        assert_eq!(bytes_needed_u64(u64::MAX), 8);

        assert_eq!(bytes_needed_u32(0), 0);
        assert_eq!(bytes_needed_u32(255), 1);
        assert_eq!(bytes_needed_u32(256), 2);
        assert_eq!(bytes_needed_u32(u32::MAX), 4);
    }

    #[test]
    fn test_radix_sort_adaptive() {
        let mut entries: Vec<(u64, usize)> = vec![
            (pack_coordinate_for_radix(1, 500, false, 100), 0),
            (pack_coordinate_for_radix(0, 100, true, 100), 1),
            (pack_coordinate_for_radix(0, 100, false, 100), 2),
            (pack_coordinate_for_radix(2, 0, false, 100), 3),
            (pack_coordinate_for_radix(-1, -1, false, 100), 4), // unmapped
        ];

        // Max tid=2, max pos=(500+1)<<1=1002
        radix_sort_coordinate_adaptive(&mut entries, 100, 1002);

        // Expected order: (0,100,false), (0,100,true), (1,500,false), (2,0,false), unmapped
        assert_eq!(entries[0].1, 2); // tid=0, pos=100, rev=false
        assert_eq!(entries[1].1, 1); // tid=0, pos=100, rev=true
        assert_eq!(entries[2].1, 0); // tid=1, pos=500, rev=false
        assert_eq!(entries[3].1, 3); // tid=2, pos=0, rev=false
        assert_eq!(entries[4].1, 4); // unmapped
    }
}
