//! Sequential BAM-file sink with reorder buffer.
//!
//! Consumes [`CompressedBatch`] items in arbitrary ordinal order and writes
//! them to a BAM file in strict ordinal order. Internally uses a
//! [`ReorderBuffer`] keyed by the batch's `ordinal`.
//!
//! On startup, writes the BAM header(s) via [`write_bam_header`]. On each
//! popped-in-order batch, appends `primary.data` (already a concatenation
//! of complete BGZF blocks) and, if present, `secondary.data` to the
//! respective file. On clean finish, writes a [`BGZF_EOF`] marker per file
//! and asserts the reorder buffer is empty — returning a pipeline error
//! that names the pending ordinals if not. Under cancellation, skips the
//! completion-validation check so that the real error propagates.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use fgumi_bgzf::BGZF_EOF;
use noodles::sam::Header;

use crate::bam_io::write_bam_header;
use crate::runall::engine::backoff::Backoff;
use crate::runall::engine::cancel::CancelToken;
use crate::runall::engine::output_types::CompressedBatch;
use crate::runall::engine::reorder::ReorderBuffer;
use crate::runall::engine::sink::{InputQueue, Sink};

/// Sequential sink that reorders [`CompressedBatch`] items by their
/// `ordinal` and appends their compressed bytes to one or two BAM files.
pub struct BamFileWrite {
    primary_path: PathBuf,
    primary_header: Header,
    secondary_path: Option<PathBuf>,
    secondary_header: Option<Header>,
    /// Unused in current logic but reserved for backpressure/queue tuning
    /// if the engine later allows sinks to publish their own capacity hints.
    _queue_capacity: usize,
}

impl BamFileWrite {
    /// Construct a new sink.
    ///
    /// # Arguments
    ///
    /// * `primary_path` — path to the primary output BAM.
    /// * `primary_header` — SAM header to serialize at the start of the primary file.
    /// * `secondary_path` — optional path to a secondary output BAM
    ///   (e.g., filter rejects). Must be `Some` iff `secondary_header` is `Some`.
    /// * `secondary_header` — optional SAM header for the secondary file.
    /// * `queue_capacity` — currently reserved for future backpressure hints.
    #[must_use]
    pub fn new(
        primary_path: PathBuf,
        primary_header: Header,
        secondary_path: Option<PathBuf>,
        secondary_header: Option<Header>,
        queue_capacity: usize,
    ) -> Self {
        debug_assert_eq!(
            secondary_path.is_some(),
            secondary_header.is_some(),
            "secondary_path and secondary_header must be Some/None together",
        );
        Self {
            primary_path,
            primary_header,
            secondary_path,
            secondary_header,
            _queue_capacity: queue_capacity,
        }
    }
}

impl Sink for BamFileWrite {
    type Input = CompressedBatch;

    fn run(
        self: Box<Self>,
        input: Box<dyn InputQueue<Self::Input>>,
        cancel: CancelToken,
    ) -> Result<()> {
        let this = *self;

        // Open files and write headers once.
        let mut primary = BufWriter::new(
            File::create(&this.primary_path)
                .with_context(|| format!("open primary BAM {}", this.primary_path.display()))?,
        );
        write_header_as_bgzf(&mut primary, &this.primary_header)
            .context("write primary BAM header")?;

        let mut secondary: Option<(PathBuf, BufWriter<File>)> =
            match (this.secondary_path.as_ref(), this.secondary_header.as_ref()) {
                (Some(path), Some(header)) => {
                    let mut w = BufWriter::new(
                        File::create(path)
                            .with_context(|| format!("open secondary BAM {}", path.display()))?,
                    );
                    write_header_as_bgzf(&mut w, header).context("write secondary BAM header")?;
                    Some((path.clone(), w))
                }
                _ => None,
            };

        // Reorder buffer keyed by ordinal.
        let mut reorder: ReorderBuffer<CompressedBatch> = ReorderBuffer::new();
        let mut backoff = Backoff::new();

        loop {
            if cancel.is_cancelled() {
                // Skip completion validation; propagate promptly.
                return Ok(());
            }

            if let Some(item) = input.pop() {
                let batch = item.item;
                let ordinal = batch.ordinal;
                reorder.push(ordinal, batch);
                while let Some(ready) = reorder.pop_ready() {
                    let primary_count = ready.primary.record_count;
                    primary
                        .write_all(&ready.primary.data)
                        .context("append compressed primary batch")?;
                    if let Some((ref path, ref mut w)) = secondary {
                        if let Some(sec) = ready.secondary {
                            w.write_all(&sec.data).with_context(|| {
                                format!("append compressed secondary batch to {}", path.display())
                            })?;
                        }
                    }
                    if primary_count > 0 {
                        crate::progress::records_written(primary_count);
                    }
                }
                backoff.reset();
            } else if input.is_drained() {
                break;
            } else {
                backoff.snooze();
            }
        }

        // Drain any remaining in-order items (shouldn't be any if queue is
        // drained and ordinals are dense, but pop_ready is safe to call on
        // an empty buffer).
        while let Some(ready) = reorder.pop_ready() {
            let primary_count = ready.primary.record_count;
            primary
                .write_all(&ready.primary.data)
                .context("drain: append compressed primary batch")?;
            if let Some((_, ref mut w)) = secondary {
                if let Some(sec) = ready.secondary {
                    w.write_all(&sec.data).context("drain: append compressed secondary batch")?;
                }
            }
            if primary_count > 0 {
                crate::progress::records_written(primary_count);
            }
        }

        // Completion validation: the reorder buffer MUST be empty. If it is
        // not, an upstream stage either lost a batch or produced duplicates
        // — either is a bug worth reporting.
        if !reorder.is_empty() {
            anyhow::bail!(
                "BamFileWrite: {} batches still pending in reorder buffer \
                 (next_expected={}); upstream stage produced non-dense ordinals",
                reorder.len(),
                reorder.next_expected(),
            );
        }

        // BGZF EOF marker + file close.
        primary.write_all(&BGZF_EOF).context("write primary BGZF EOF")?;
        primary.flush().context("flush primary BAM")?;
        if let Some((path, mut w)) = secondary {
            w.write_all(&BGZF_EOF)
                .with_context(|| format!("write secondary BGZF EOF for {}", path.display()))?;
            w.flush().with_context(|| format!("flush secondary BAM {}", path.display()))?;
        }

        Ok(())
    }

    fn name(&self) -> &'static str {
        "BamFileWrite"
    }
}

/// Write a BAM header and then a single terminating BGZF flush boundary.
///
/// `write_bam_header` emits raw (uncompressed) BAM header bytes. For them to
/// form a valid BAM file when concatenated with downstream pre-compressed
/// blocks, the header must itself be written as one or more BGZF blocks.
/// We do this by wrapping the header bytes in an [`InlineBgzfCompressor`]
/// and flushing, producing a compact header-only BGZF stream.
fn write_header_as_bgzf(writer: &mut impl Write, header: &Header) -> Result<()> {
    use fgumi_bgzf::InlineBgzfCompressor;
    let mut header_buf: Vec<u8> = Vec::new();
    write_bam_header(&mut header_buf, header).context("serialize BAM header")?;

    let mut compressor = InlineBgzfCompressor::new(1);
    compressor.write_all(&header_buf).context("compress BAM header")?;
    compressor.flush().context("flush BAM header compressor")?;
    for block in compressor.take_blocks() {
        writer.write_all(&block.data).context("write header BGZF block")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::runall::engine::memory::MemoryTracker;
    use crate::runall::engine::output_types::CompressedBytes;
    use crate::runall::engine::queue::StageQueue;
    use crate::runall::engine::stage::SequencedItem;
    use fgumi_bgzf::InlineBgzfCompressor;
    use std::sync::Arc;

    fn empty_header() -> Header {
        use noodles::sam::header::record::value::Map;
        Header::builder()
            .set_header(Map::<noodles::sam::header::record::value::map::Header>::new(
                noodles::sam::header::record::value::map::header::Version::new(1, 6),
            ))
            .build()
    }

    /// Produce a one-batch `CompressedBatch` whose primary `data` contains the
    /// BGZF-compressed payload `bytes`. Makes assertions possible without
    /// needing real serialized BAM records.
    fn compressed_batch(bytes: &[u8], ordinal: u64) -> CompressedBatch {
        let mut c = InlineBgzfCompressor::new(1);
        c.write_all(bytes).unwrap();
        c.flush().unwrap();
        let mut data = Vec::new();
        for b in c.take_blocks() {
            data.extend_from_slice(&b.data);
        }
        CompressedBatch {
            primary: CompressedBytes { data, record_count: 1 },
            secondary: None,
            ordinal,
        }
    }

    fn new_queue() -> Arc<StageQueue<CompressedBatch>> {
        let tracker = Arc::new(MemoryTracker::new(1_000_000_000));
        Arc::new(StageQueue::new("test_sink_in", 64, 100_000_000, tracker))
    }

    #[test]
    fn test_out_of_order_ordinals_write_in_order() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let q = new_queue();
        // Push ordinals 2, 0, 1. Write contents are ordinal-tagged for recognition.
        q.push(SequencedItem::new(2, compressed_batch(b"BATCH-2", 2), 16)).unwrap();
        q.push(SequencedItem::new(0, compressed_batch(b"BATCH-0", 0), 16)).unwrap();
        q.push(SequencedItem::new(1, compressed_batch(b"BATCH-1", 1), 16)).unwrap();
        q.close();

        let sink =
            Box::new(BamFileWrite::new(tmp.path().to_path_buf(), empty_header(), None, None, 16));
        sink.run(Box::new(q), CancelToken::new()).expect("sink must succeed");

        // Open the written file and confirm the decompressed payload ordering.
        // Header BGZF block(s) come first, then BATCH-0, BATCH-1, BATCH-2, EOF.
        let bytes = std::fs::read(tmp.path()).unwrap();
        let joined = decompress_all_bgzf(&bytes);
        // After the BAM header bytes, the three batch payloads appear in
        // strict ordinal order.
        let marker0 = joined
            .windows(b"BATCH-0".len())
            .position(|w| w == b"BATCH-0")
            .expect("BATCH-0 must be present");
        let marker1 = joined
            .windows(b"BATCH-1".len())
            .position(|w| w == b"BATCH-1")
            .expect("BATCH-1 must be present");
        let marker2 = joined
            .windows(b"BATCH-2".len())
            .position(|w| w == b"BATCH-2")
            .expect("BATCH-2 must be present");
        assert!(marker0 < marker1 && marker1 < marker2, "ordinal ordering violated");
    }

    #[test]
    fn test_completion_validation_detects_gap() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let q = new_queue();
        // Ordinals 0, 2 — ordinal 1 is missing. The sink should fail to drain.
        q.push(SequencedItem::new(0, compressed_batch(b"BATCH-0", 0), 16)).unwrap();
        q.push(SequencedItem::new(2, compressed_batch(b"BATCH-2", 2), 16)).unwrap();
        q.close();

        let sink =
            Box::new(BamFileWrite::new(tmp.path().to_path_buf(), empty_header(), None, None, 16));
        let err = sink
            .run(Box::new(q), CancelToken::new())
            .expect_err("sink must fail on non-dense ordinals");
        let msg = err.to_string();
        assert!(
            msg.contains("pending") && msg.contains("next_expected=1"),
            "error must name the pending ordinal: {msg}"
        );
    }

    #[test]
    fn test_empty_input_produces_valid_empty_bam() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let q = new_queue();
        q.close();

        let sink =
            Box::new(BamFileWrite::new(tmp.path().to_path_buf(), empty_header(), None, None, 16));
        sink.run(Box::new(q), CancelToken::new()).expect("sink must succeed");

        let bytes = std::fs::read(tmp.path()).unwrap();
        // Must contain at least the BGZF-compressed header and the EOF marker.
        assert!(bytes.len() >= BGZF_EOF.len(), "empty BAM must include BGZF EOF marker");
        assert!(bytes.ends_with(&BGZF_EOF), "empty BAM must end with BGZF EOF");
    }

    #[test]
    fn test_header_written_exactly_once() {
        // Header is serialized at sink startup, before the input loop. We
        // verify this indirectly: with 3 batches, there is still exactly one
        // BAM magic ("BAM\x01") marker in the decompressed stream.
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let q = new_queue();
        for i in 0..3u64 {
            let payload = format!("PAYLOAD-{i}");
            q.push(SequencedItem::new(i, compressed_batch(payload.as_bytes(), i), 32)).unwrap();
        }
        q.close();

        let sink =
            Box::new(BamFileWrite::new(tmp.path().to_path_buf(), empty_header(), None, None, 16));
        sink.run(Box::new(q), CancelToken::new()).expect("sink must succeed");

        let bytes = std::fs::read(tmp.path()).unwrap();
        let joined = decompress_all_bgzf(&bytes);
        let bam_magic = b"BAM\x01";
        let magic_count = joined.windows(bam_magic.len()).filter(|w| *w == bam_magic).count();
        assert_eq!(magic_count, 1, "BAM magic must appear exactly once");
    }

    #[test]
    fn test_cancellation_skips_completion_validation() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let q = new_queue();
        // Leave the queue open and non-drained so the sink would normally
        // spin; set cancel before the first pop.
        let cancel = CancelToken::new();
        cancel.cancel();

        let sink =
            Box::new(BamFileWrite::new(tmp.path().to_path_buf(), empty_header(), None, None, 16));
        // With cancel already set, the sink returns Ok(()) without writing
        // EOF and without running the pending-batch check (nothing to check).
        sink.run(Box::new(q), cancel).expect("cancelled sink must return Ok");
    }

    /// Pull every BGZF block from `bytes` and return the concatenated
    /// decompressed payload. Useful for asserting on the logical contents of
    /// a produced BAM file without reconstructing full records.
    ///
    /// Note: the real `fgumi_bgzf::decompress_block` requires a mutable
    /// `libdeflater::Decompressor`; we create one here for the test helper.
    fn decompress_all_bgzf(bytes: &[u8]) -> Vec<u8> {
        use libdeflater::Decompressor;
        // Strip the trailing EOF marker if present — it decompresses to zero
        // bytes anyway but keeps the helper tolerant of variants.
        let effective: &[u8] =
            if bytes.ends_with(&BGZF_EOF) { &bytes[..bytes.len() - BGZF_EOF.len()] } else { bytes };
        let mut reader = std::io::Cursor::new(effective);
        let mut decompressor = Decompressor::new();
        let mut out = Vec::new();
        loop {
            let blocks = fgumi_bgzf::read_raw_blocks(&mut reader, 1)
                .expect("read_raw_blocks succeeds on well-formed input");
            if blocks.is_empty() {
                break;
            }
            for block in &blocks {
                let payload = fgumi_bgzf::decompress_block(block, &mut decompressor)
                    .expect("decompress_block succeeds on well-formed blocks");
                out.extend_from_slice(&payload);
            }
        }
        out
    }
}
