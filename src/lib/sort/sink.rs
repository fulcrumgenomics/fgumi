//! Sink abstraction for consuming sorted records from the merge loop.
//!
//! The [`SortedRecordSink`] trait decouples the sort/merge phase from the output destination,
//! allowing sorted records to flow into a BAM writer, a downstream pipeline stage, or a
//! test collector without changing the merge logic.

use anyhow::Result;

use super::inline_buffer::TemplateKey;
use crate::bam_io::RawBamWriter;

/// Consumes sorted records one at a time from the sort merge loop.
///
/// Implementations receive each record in sorted order along with its sort key.
/// After all records have been emitted, [`finish`](SortedRecordSink::finish) is called
/// to finalize any underlying resources (e.g. flush and close a BAM writer).
pub trait SortedRecordSink {
    /// Write a single sorted record.
    ///
    /// `key` is the sort key that was used to order this record.
    /// `record_bytes` contains the raw BAM record bytes (without the 4-byte `block_size` prefix).
    ///
    /// # Errors
    ///
    /// Returns an error if the record cannot be written to the underlying destination.
    fn emit(&mut self, key: &TemplateKey, record_bytes: Vec<u8>) -> Result<()>;

    /// Finalize the sink, flushing any buffered data and releasing resources.
    ///
    /// # Errors
    ///
    /// Returns an error if finalization fails (e.g. flushing or closing the output).
    fn finish(&mut self) -> Result<()>;
}

/// Default sink that writes sorted records directly to a BAM file.
///
/// Wraps a [`RawBamWriter`] and forwards each emitted record as raw BAM bytes.
/// The writer is consumed on [`finish`](SortedRecordSink::finish) to ensure the
/// BGZF EOF block is written.
pub struct BamWriterSink {
    writer: Option<RawBamWriter>,
    records_written: u64,
}

impl BamWriterSink {
    /// Create a new sink backed by the given BAM writer.
    ///
    /// The writer should already have its BAM header written before records are emitted.
    #[must_use]
    pub fn new(writer: RawBamWriter) -> Self {
        Self { writer: Some(writer), records_written: 0 }
    }

    /// Returns the number of records written so far.
    #[must_use]
    pub fn records_written(&self) -> u64 {
        self.records_written
    }
}

impl SortedRecordSink for BamWriterSink {
    fn emit(&mut self, _key: &TemplateKey, record_bytes: Vec<u8>) -> Result<()> {
        self.writer.as_mut().expect("emit called after finish").write_raw_record(&record_bytes)?;
        self.records_written += 1;
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        if let Some(writer) = self.writer.take() {
            writer.finish()?;
        }
        Ok(())
    }
}

/// Test sink that collects emitted records for verification.
///
/// Available as `pub(crate)` so that other modules (e.g. `raw.rs`) can reuse the
/// same sink in their tests without duplicating the definition.
#[cfg(test)]
pub(crate) struct CollectingSink {
    pub(crate) records: Vec<(TemplateKey, Vec<u8>)>,
    pub(crate) finished: bool,
}

#[cfg(test)]
impl CollectingSink {
    pub(crate) fn new() -> Self {
        Self { records: Vec::new(), finished: false }
    }
}

#[cfg(test)]
impl SortedRecordSink for CollectingSink {
    fn emit(&mut self, key: &TemplateKey, record_bytes: Vec<u8>) -> Result<()> {
        self.records.push((*key, record_bytes));
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        self.finished = true;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_collecting_sink_stores_records() {
        let mut sink = CollectingSink::new();
        let key = TemplateKey::default();

        sink.emit(&key, vec![1, 2, 3]).expect("emit should succeed");
        sink.emit(&key, vec![4, 5, 6]).expect("emit should succeed");
        sink.finish().expect("finish should succeed");

        assert_eq!(sink.records.len(), 2);
        assert_eq!(sink.records[0].1, vec![1, 2, 3]);
        assert_eq!(sink.records[1].1, vec![4, 5, 6]);
        assert_eq!(sink.records[0].0, TemplateKey::default());
        assert!(sink.finished);
    }

    #[test]
    fn test_collecting_sink_finish_is_idempotent() {
        let mut sink = CollectingSink::new();
        sink.finish().expect("first finish should succeed");
        sink.finish().expect("second finish should succeed");
        assert!(sink.finished);
    }
}
