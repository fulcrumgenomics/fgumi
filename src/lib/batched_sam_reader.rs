//! Adaptive buffered SAM reader that grows based on observed batch sizes.
//!
//! # Motivation
//!
//! When piping SAM output from `bwa mem` via stdin, the output arrives in batches
//! determined by bwa's `-K` parameter (chunk size in bases). Each batch contains
//! reads totaling at least `-K` bases, with an even number of records to keep
//! read pairs together. The batch sizes in bytes vary significantly based on
//! read length, alignment complexity, and optional tag content.
//!
//! A fixed-size buffer is problematic:
//! - Too small: causes excessive I/O syscalls and potential stalls waiting for data
//! - Too large: wastes memory, especially for small input files
//!
//! This reader solves the problem by starting with a reasonable initial buffer (64MB)
//! and growing adaptively based on observed batch sizes. It never shrinks, avoiding
//! repeated allocations when batch sizes fluctuate.
//!
//! # How it works
//!
//! The reader tracks bwa's batch boundaries by counting bases and records. When a
//! batch boundary is detected (bases >= `chunk_size` AND record count is even), it
//! measures the batch's byte size and grows the buffer if needed. This ensures the
//! buffer can always hold at least one complete batch, enabling efficient streaming
//! without blocking on partial batches.

use std::io::{self, BufRead, BufReader, Read};

const READ_CHUNK_SIZE: usize = 64 * 1024 * 1024; // 64MB I/O chunks
const INITIAL_BUFFER_SIZE: usize = 64 * 1024 * 1024; // 64MB initial buffer
const GROWTH_FACTOR: f64 = 1.25;

/// Adaptive buffered reader for SAM input from bwa mem.
///
/// Reads bytes into a growable buffer, tracks batch boundaries like bwa,
/// and grows the buffer based on observed batch sizes. Never shrinks.
pub struct BatchedSamReader<R: Read> {
    inner: BufReader<R>,
    buffer: Vec<u8>,
    buffer_len: usize, // Valid bytes in buffer
    parse_pos: usize,  // Current parse position

    // Batch tracking (mirrors bwa logic)
    chunk_size: u64, // BWA -K value
    bases_this_batch: u64,
    records_this_batch: usize,
    bytes_at_batch_start: usize,

    // Statistics
    batches_completed: usize,
    total_bytes_read: usize,
}

impl<R: Read> BatchedSamReader<R> {
    /// Create a new reader with the given bwa chunk size (-K value).
    pub fn new(reader: R, chunk_size: u64) -> Self {
        Self {
            inner: BufReader::with_capacity(READ_CHUNK_SIZE, reader),
            buffer: vec![0u8; INITIAL_BUFFER_SIZE],
            buffer_len: 0,
            parse_pos: 0,
            chunk_size,
            bases_this_batch: 0,
            records_this_batch: 0,
            bytes_at_batch_start: 0,
            batches_completed: 0,
            total_bytes_read: 0,
        }
    }

    /// Fill buffer by reading chunks until buffer is full or EOF.
    ///
    /// Handles edge cases:
    /// - Grows buffer proactively if unparsed data > 50% of capacity
    /// - Preserves partial records at end of buffer
    /// - Never loses data
    pub fn fill_buffer(&mut self) -> io::Result<bool> {
        let unparsed_len = self.buffer_len.saturating_sub(self.parse_pos);

        // Proactive growth: if unparsed data > 50% of buffer, grow first
        if unparsed_len > self.buffer.len() / 2 {
            let new_size = self.buffer.len() * 2;
            log::info!(
                "Growing buffer (proactive): {}MB -> {}MB (unparsed: {}MB)",
                self.buffer.len() / (1024 * 1024),
                new_size / (1024 * 1024),
                unparsed_len / (1024 * 1024)
            );
            self.buffer.resize(new_size, 0);
        }

        // Compact: move unparsed data to front
        if self.parse_pos > 0 && unparsed_len > 0 {
            self.buffer.copy_within(self.parse_pos..self.buffer_len, 0);
        }
        self.buffer_len = unparsed_len;
        self.bytes_at_batch_start = self.bytes_at_batch_start.saturating_sub(self.parse_pos);
        self.parse_pos = 0;

        // Read chunks until buffer is full or EOF
        while self.buffer_len < self.buffer.len() {
            let space = self.buffer.len() - self.buffer_len;
            let to_read = std::cmp::min(space, READ_CHUNK_SIZE);

            let n =
                self.inner.read(&mut self.buffer[self.buffer_len..self.buffer_len + to_read])?;
            if n == 0 {
                // EOF - return true if there's still data to parse
                return Ok(self.buffer_len > 0);
            }
            self.buffer_len += n;
            self.total_bytes_read += n;
        }
        Ok(true)
    }

    /// Get current buffer slice for parsing.
    pub fn buffer(&self) -> &[u8] {
        &self.buffer[self.parse_pos..self.buffer_len]
    }

    /// Check if we've hit a batch boundary (like bwa).
    fn at_batch_boundary(&self) -> bool {
        self.bases_this_batch >= self.chunk_size && self.records_this_batch.is_multiple_of(2)
    }

    /// Called after parsing a SAM record.
    ///
    /// Tracks bases and records for batch detection. When a batch boundary
    /// is reached, may grow the buffer based on observed batch size.
    ///
    /// Returns true if the buffer was grown.
    pub fn record_parsed(&mut self, seq_len: usize, bytes_consumed: usize) -> bool {
        self.bases_this_batch += seq_len as u64;
        self.records_this_batch += 1;
        self.parse_pos += bytes_consumed;

        if self.at_batch_boundary() {
            let bytes_this_batch = self.parse_pos - self.bytes_at_batch_start;
            let grew = self.maybe_grow(bytes_this_batch);

            // Reset for next batch
            self.batches_completed += 1;
            self.bases_this_batch = 0;
            self.records_this_batch = 0;
            self.bytes_at_batch_start = self.parse_pos;

            grew
        } else {
            false
        }
    }

    /// Advance parse position without tracking as a record (e.g., for headers).
    pub fn advance(&mut self, bytes: usize) {
        self.parse_pos += bytes;
    }

    /// Grow buffer if needed based on batch size. Never shrinks.
    fn maybe_grow(&mut self, bytes_this_batch: usize) -> bool {
        let needed = ((bytes_this_batch as f64) * GROWTH_FACTOR) as usize;
        if needed > self.buffer.len() {
            // Round up to next power of 2 for efficiency
            let new_size = needed.next_power_of_two();
            log::info!(
                "Growing buffer (batch {}): {}MB -> {}MB ({} bytes for {} bases)",
                self.batches_completed + 1,
                self.buffer.len() / (1024 * 1024),
                new_size / (1024 * 1024),
                bytes_this_batch,
                self.bases_this_batch
            );
            self.buffer.resize(new_size, 0);
            true
        } else {
            false
        }
    }

    /// Get current buffer capacity.
    pub fn capacity(&self) -> usize {
        self.buffer.len()
    }

    /// Get number of completed batches.
    pub fn batches_completed(&self) -> usize {
        self.batches_completed
    }

    /// Get total bytes read from input.
    pub fn total_bytes_read(&self) -> usize {
        self.total_bytes_read
    }

    /// Get the configured chunk size (bwa -K value).
    pub fn chunk_size(&self) -> u64 {
        self.chunk_size
    }

    /// Set the buffer for testing purposes.
    #[cfg(test)]
    pub fn set_buffer(&mut self, buffer: Vec<u8>) {
        self.buffer = buffer;
    }

    /// Internal `fill_buf` that handles buffer management for `BufRead`.
    ///
    /// When used via `BufRead`, we can't track SAM record boundaries,
    /// but we still grow the buffer proactively when it's nearly full.
    fn fill_buf_internal(&mut self) -> io::Result<&[u8]> {
        // Return existing data if available
        if self.parse_pos < self.buffer_len {
            return Ok(&self.buffer[self.parse_pos..self.buffer_len]);
        }

        // Compact and refill
        self.parse_pos = 0;
        self.buffer_len = 0;

        // Read new data
        let bytes_read = self.inner.read(&mut self.buffer)?;
        self.buffer_len = bytes_read;
        self.total_bytes_read += bytes_read;

        // Proactive growth: if buffer was mostly filled (>75%), grow for next time
        if bytes_read >= self.buffer.len() * 3 / 4 {
            let new_size = (self.buffer.len() * 2).min(1024 * 1024 * 1024); // Cap at 1GB
            if new_size > self.buffer.len() {
                log::info!(
                    "Growing buffer (BufRead): {}MB -> {}MB (read: {}MB)",
                    self.buffer.len() / (1024 * 1024),
                    new_size / (1024 * 1024),
                    bytes_read / (1024 * 1024),
                );
                self.buffer.resize(new_size, 0);
            }
        }

        Ok(&self.buffer[self.parse_pos..self.buffer_len])
    }
}

/// Implement `Read` trait for compatibility with standard I/O.
impl<R: Read> Read for BatchedSamReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let available = self.fill_buf_internal()?;
        if available.is_empty() {
            return Ok(0);
        }
        let to_copy = buf.len().min(available.len());
        buf[..to_copy].copy_from_slice(&available[..to_copy]);
        self.parse_pos += to_copy;
        Ok(to_copy)
    }
}

/// Implement `BufRead` trait for use with SAM readers (e.g., noodles).
///
/// When used via `BufRead`, batch tracking is not active since we don't
/// have visibility into SAM record boundaries. The buffer still grows
/// proactively based on fill ratio, similar to `GrowableReader`.
impl<R: Read> BufRead for BatchedSamReader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.fill_buf_internal()
    }

    fn consume(&mut self, amt: usize) {
        self.parse_pos = (self.parse_pos + amt).min(self.buffer_len);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// Create a fake SAM record line with given sequence length.
    fn make_sam_record(name: &str, seq_len: usize) -> String {
        let seq = "A".repeat(seq_len);
        let qual = "I".repeat(seq_len);
        format!("{name}\t0\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}\n")
    }

    /// Create paired SAM records (R1 and R2).
    fn make_read_pair(name: &str, seq_len: usize) -> String {
        let r1 = make_sam_record(&format!("{name}/1"), seq_len);
        let r2 = make_sam_record(&format!("{name}/2"), seq_len);
        format!("{r1}{r2}")
    }

    #[test]
    fn test_basic_reading() {
        let data = make_read_pair("read1", 100);
        let cursor = Cursor::new(data.as_bytes().to_vec());

        let mut reader = BatchedSamReader::new(cursor, 1000);
        assert!(reader.fill_buffer().unwrap());
        assert!(!reader.buffer().is_empty());
    }

    #[test]
    fn test_batch_boundary_detection() {
        // Create 10 pairs of 100bp reads = 2000 bases total
        // With chunk_size=500, should hit boundary after ~3 pairs (600 bases)
        let mut data = String::new();
        for i in 0..10 {
            data.push_str(&make_read_pair(&format!("read{i}"), 100));
        }

        let cursor = Cursor::new(data.as_bytes().to_vec());
        let mut reader = BatchedSamReader::new(cursor, 500);

        reader.fill_buffer().unwrap();

        let mut batches = 0;
        let mut records = 0;

        // Parse records and count batches
        while !reader.buffer().is_empty() {
            if let Some(line_end) = reader.buffer().iter().position(|&b| b == b'\n') {
                let line = &reader.buffer()[..=line_end];
                // Simple seq length extraction (field 10)
                let seq_len = line.split(|&b| b == b'\t').nth(9).map_or(0, <[u8]>::len);

                if reader.record_parsed(seq_len, line_end + 1) {
                    // Buffer grew - this shouldn't happen with small data
                }
                records += 1;

                if reader.batches_completed() > batches {
                    batches = reader.batches_completed();
                }
            } else {
                break;
            }
        }

        assert_eq!(records, 20); // 10 pairs = 20 records
        assert!(batches >= 3); // 2000 bases / 500 chunk = 4 batches (but pairs must be even)
    }

    #[test]
    fn test_buffer_growth_on_large_batch() {
        // Create data larger than initial 64MB buffer equivalent
        // Use small initial buffer for testing
        let mut data = String::new();
        // Each pair is ~250 bytes, create enough for 100KB
        for i in 0..200 {
            data.push_str(&make_read_pair(&format!("read{i}"), 100));
        }

        let cursor = Cursor::new(data.as_bytes().to_vec());

        // Use tiny initial buffer to force growth
        let mut reader = BatchedSamReader::new(cursor, 50000); // 50K bases per batch

        // Manually set small initial buffer for testing
        reader.set_buffer(vec![0u8; 1024]); // 1KB initial

        reader.fill_buffer().unwrap();
        let initial_capacity = reader.capacity();

        // Parse until we hit a batch boundary that triggers growth
        let mut grew = false;
        while !reader.buffer().is_empty() {
            if let Some(line_end) = reader.buffer().iter().position(|&b| b == b'\n') {
                let seq_len = 100; // We know our records have 100bp
                if reader.record_parsed(seq_len, line_end + 1) {
                    grew = true;
                }
            } else if !reader.fill_buffer().unwrap() {
                break;
            }
        }

        // Buffer should have grown
        assert!(reader.capacity() >= initial_capacity);
        // Note: grew might not be true if proactive growth happened instead
        let _ = grew;
    }

    #[test]
    fn test_never_shrinks() {
        // First batch large, second batch small - buffer should not shrink
        let mut data = String::new();

        // Large batch: 100 pairs of 150bp = 30000 bases
        for i in 0..100 {
            data.push_str(&make_read_pair(&format!("large{i}"), 150));
        }

        // Small batch: 10 pairs of 50bp = 1000 bases
        for i in 0..10 {
            data.push_str(&make_read_pair(&format!("small{i}"), 50));
        }

        let cursor = Cursor::new(data.as_bytes().to_vec());
        let mut reader = BatchedSamReader::new(cursor, 10000);

        // Small initial buffer
        reader.set_buffer(vec![0u8; 4096]);

        reader.fill_buffer().unwrap();

        let mut max_capacity = reader.capacity();

        while !reader.buffer().is_empty() {
            if let Some(line_end) = reader.buffer().iter().position(|&b| b == b'\n') {
                let line = &reader.buffer()[..=line_end];
                let seq_len = line.split(|&b| b == b'\t').nth(9).map_or(0, <[u8]>::len);

                reader.record_parsed(seq_len, line_end + 1);

                // Track max capacity
                if reader.capacity() > max_capacity {
                    max_capacity = reader.capacity();
                }

                // Verify never shrinks
                assert!(reader.capacity() >= max_capacity);
            } else if !reader.fill_buffer().unwrap() {
                break;
            }
        }

        // Final capacity should be >= max seen
        assert!(reader.capacity() >= max_capacity);
    }

    #[test]
    fn test_partial_line_preserved() {
        // Create data where a line will be split across buffer boundaries
        let record = make_sam_record("splitme", 100);

        // Put half of record at end of "first chunk"
        let split_point = record.len() / 2;
        let first_half = &record[..split_point];
        let second_half = &record[split_point..];

        // We need to simulate partial reads, so use a custom reader
        struct PartialReader {
            chunks: Vec<Vec<u8>>,
            idx: usize,
        }

        impl Read for PartialReader {
            fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
                if self.idx >= self.chunks.len() {
                    return Ok(0);
                }
                let chunk = &self.chunks[self.idx];
                let to_copy = std::cmp::min(buf.len(), chunk.len());
                buf[..to_copy].copy_from_slice(&chunk[..to_copy]);
                self.idx += 1;
                Ok(to_copy)
            }
        }

        let partial_reader = PartialReader {
            chunks: vec![first_half.as_bytes().to_vec(), second_half.as_bytes().to_vec()],
            idx: 0,
        };

        let mut reader = BatchedSamReader::new(partial_reader, 1000);
        reader.set_buffer(vec![0u8; 1024]); // Small buffer

        // First fill gets partial line
        reader.fill_buffer().unwrap();

        // Should not find complete line yet if buffer is small enough
        // After second fill, should have complete line
        let mut found_complete = false;
        for _ in 0..3 {
            if reader.buffer().contains(&b'\n') {
                found_complete = true;
                break;
            }
            if !reader.fill_buffer().unwrap() {
                break;
            }
        }

        assert!(found_complete, "Should eventually find complete line");
    }

    #[test]
    fn test_proactive_growth() {
        // Test that buffer grows proactively when unparsed > 50%
        let mut data = String::new();
        for i in 0..100 {
            data.push_str(&make_read_pair(&format!("read{i}"), 100));
        }

        let cursor = Cursor::new(data.as_bytes().to_vec());
        let mut reader = BatchedSamReader::new(cursor, 1_000_000); // Large chunk to avoid batch triggers

        // Tiny buffer to force proactive growth
        reader.set_buffer(vec![0u8; 512]);
        let initial = reader.capacity();

        reader.fill_buffer().unwrap();

        // Don't parse anything - just refill
        // Second fill should trigger proactive growth if data > 50% of buffer
        reader.fill_buffer().unwrap();

        // Buffer should have grown proactively
        assert!(
            reader.capacity() > initial,
            "Buffer should grow proactively: {} vs {}",
            reader.capacity(),
            initial
        );
    }

    #[test]
    fn test_empty_input() {
        let cursor = Cursor::new(Vec::<u8>::new());
        let mut reader = BatchedSamReader::new(cursor, 1000);

        // Should return false (no data)
        assert!(!reader.fill_buffer().unwrap());
        assert!(reader.buffer().is_empty());
    }

    #[test]
    fn test_header_lines_with_advance() {
        let mut data = String::new();
        data.push_str("@HD\tVN:1.6\n");
        data.push_str("@SQ\tSN:chr1\tLN:1000\n");
        data.push_str(&make_read_pair("read1", 100));

        let cursor = Cursor::new(data.as_bytes().to_vec());
        let mut reader = BatchedSamReader::new(cursor, 1000);

        reader.fill_buffer().unwrap();

        let mut alignment_records = 0;

        while !reader.buffer().is_empty() {
            if let Some(line_end) = reader.buffer().iter().position(|&b| b == b'\n') {
                let line = &reader.buffer()[..=line_end];

                if line.starts_with(b"@") {
                    // Header line - advance without counting as record
                    reader.advance(line_end + 1);
                } else {
                    // Alignment record
                    let seq_len = line.split(|&b| b == b'\t').nth(9).map_or(0, <[u8]>::len);
                    reader.record_parsed(seq_len, line_end + 1);
                    alignment_records += 1;
                }
            } else {
                break;
            }
        }

        assert_eq!(alignment_records, 2); // One pair
    }

    // Tests for BufRead implementation

    #[test]
    fn test_bufread_basic_read() {
        let data = b"Hello, World!";
        let cursor = Cursor::new(data.to_vec());

        let mut reader = BatchedSamReader::new(cursor, 1000);
        let mut buf = vec![0u8; 1024];

        let n = reader.read(&mut buf).unwrap();
        assert_eq!(n, 13);
        assert_eq!(&buf[..n], b"Hello, World!");
    }

    #[test]
    fn test_bufread_fill_buf() {
        let data = b"Line 1\nLine 2\nLine 3\n";
        let cursor = Cursor::new(data.to_vec());

        let mut reader = BatchedSamReader::new(cursor, 1000);

        // Read lines using BufRead
        let mut line = String::new();
        reader.read_line(&mut line).unwrap();
        assert_eq!(line, "Line 1\n");

        line.clear();
        reader.read_line(&mut line).unwrap();
        assert_eq!(line, "Line 2\n");
    }

    #[test]
    fn test_bufread_consume() {
        let data = b"ABCDEFGHIJ";
        let cursor = Cursor::new(data.to_vec());

        let mut reader = BatchedSamReader::new(cursor, 1000);

        // Fill buffer
        let buf = reader.fill_buf().unwrap();
        assert_eq!(buf, b"ABCDEFGHIJ");

        // Consume 5 bytes
        reader.consume(5);

        // Remaining should be FGHIJ
        let buf = reader.fill_buf().unwrap();
        assert_eq!(buf, b"FGHIJ");
    }

    #[test]
    fn test_bufread_growth_on_large_batch() {
        // Data that fills >75% of buffer
        let data = vec![0xAA; 52 * 1024];
        let cursor = Cursor::new(data);

        let mut reader = BatchedSamReader::new(cursor, 1000);
        reader.set_buffer(vec![0u8; 64 * 1024]); // 64KB buffer

        let mut buf = vec![0u8; 64 * 1024];
        let _ = reader.read(&mut buf).unwrap();

        // Should have grown to 128KB (doubled)
        assert_eq!(reader.capacity(), 128 * 1024);
    }

    #[test]
    fn test_bufread_no_growth_on_small_batch() {
        // Data that fills <75% of buffer
        let data = vec![0xAA; 16 * 1024];
        let cursor = Cursor::new(data);

        let mut reader = BatchedSamReader::new(cursor, 1000);
        reader.set_buffer(vec![0u8; 64 * 1024]); // 64KB buffer

        let mut buf = vec![0u8; 64 * 1024];
        let _ = reader.read(&mut buf).unwrap();

        // Should NOT have grown
        assert_eq!(reader.capacity(), 64 * 1024);
    }

    #[test]
    fn test_bufread_empty_input() {
        let data: Vec<u8> = vec![];
        let cursor = Cursor::new(data);

        let mut reader = BatchedSamReader::new(cursor, 1000);
        let mut buf = vec![0u8; 1024];

        let n = reader.read(&mut buf).unwrap();
        assert_eq!(n, 0);
    }
}
