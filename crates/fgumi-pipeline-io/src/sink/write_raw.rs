//! `WriteRawFile` sink step: writes a byte stream verbatim to a file or stdout,
//! with **no** container header and **no** BGZF EOF marker.
//!
//! Unlike [`super::write_bgzf::WriteBgzfFile`] (which is BAM-specific — it emits
//! a BAM header on open and a BGZF EOF on drain), this sink just concatenates
//! the `bytes` of each block it receives. It backs FASTQ output: the chain's
//! FASTQ-encode step produces `DecompressedBlock`s of FASTQ text, which either
//! go straight here (plain output / stdout) or through `BgzfCompress` first
//! (`.gz`/`.bgz` output, producing `BgzfBlock`s — still just bytes to write).
//!
//! `Serial + Affinity::Writer + sticky`, matching `WriteBgzfFile`: exactly one
//! shared instance drains the (reorder-ordered) block stream to the sink.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use parking_lot::Mutex;

use crate::types::{BgzfBlock, DecompressedBlock};
use fgumi_pipeline_core::{
    item::HeapSize,
    step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// A pipeline block whose payload is a run of bytes to write verbatim.
pub trait RawBytesBlock: Send + HeapSize + 'static {
    /// The bytes to write for this block.
    fn bytes(&self) -> &[u8];
}

impl RawBytesBlock for DecompressedBlock {
    fn bytes(&self) -> &[u8] {
        &self.bytes
    }
}

impl RawBytesBlock for BgzfBlock {
    fn bytes(&self) -> &[u8] {
        &self.bytes
    }
}

/// `Serial + sticky` sink that writes each block's bytes verbatim to a file or
/// stdout. Generic over the block type so it serves both the plain
/// (`DecompressedBlock`) and BGZF-compressed (`BgzfBlock`) FASTQ tails.
pub struct WriteRawFile<B> {
    state: Mutex<Option<BufWriter<Box<dyn Write + Send>>>>,
    /// Bytes appended once, after the last block, on a clean drain. Empty for
    /// plain output; the 28-byte BGZF EOF marker for BGZF output so the `.gz`
    /// stream is a complete, non-truncated BGZF file.
    trailer: &'static [u8],
    _marker: std::marker::PhantomData<fn(B)>,
}

impl<B> WriteRawFile<B> {
    /// Open `path` for writing (`-` selects stdout), appending `trailer` once
    /// after the final block on clean completion. No header is written.
    ///
    /// # Errors
    ///
    /// Returns I/O errors from opening the file.
    pub fn new<P: AsRef<Path>>(path: P, trailer: &'static [u8]) -> io::Result<Self> {
        let inner: Box<dyn Write + Send> = if path.as_ref().as_os_str() == "-" {
            Box::new(io::stdout())
        } else {
            Box::new(File::create(path.as_ref())?)
        };
        Ok(Self {
            state: Mutex::new(Some(BufWriter::with_capacity(256 * 1024, inner))),
            trailer,
            _marker: std::marker::PhantomData,
        })
    }
}

impl<B: RawBytesBlock> Step for WriteRawFile<B> {
    type Input = B;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: "WriteRawFile",
            kind: StepKind::Serial,
            sticky: true,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn affinity(&self) -> Affinity {
        Affinity::Writer
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let mut guard = self.state.lock();
        let Some(out) = guard.as_mut() else {
            return Ok(StepOutcome::Finished);
        };

        if let Some(block) = ctx.input.pop() {
            out.write_all(block.bytes())?;
            return Ok(StepOutcome::Progress);
        }

        if ctx.input.is_drained() {
            // Clean end-of-stream: append the trailer (e.g. the BGZF EOF marker)
            // exactly once, then flush and retire the sink.
            if !self.trailer.is_empty() {
                out.write_all(self.trailer)?;
            }
            out.flush()?;
            let _ = guard.take();
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

impl<B> Drop for WriteRawFile<B> {
    /// Cleanup-only drop: flush buffered bytes but deliberately do **not** write
    /// `trailer`. The trailer (e.g. the BGZF EOF marker for `.gz` output) is
    /// written exclusively by the drained-finish path in `try_run`, which then
    /// takes the state so this drop is a no-op for a cleanly-finished sink. If
    /// state is still present here the stream was aborted mid-way, and stamping
    /// the BGZF EOF onto a truncated `.gz` would make it look complete and hide
    /// the truncation from readers — so we withhold it, exactly as
    /// `WriteBgzfFile::drop` does.
    fn drop(&mut self) {
        if let Some(mut out) = self.state.lock().take() {
            let _ = out.flush();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn raw_bytes_block_exposes_payload() {
        let plain = DecompressedBlock { batch_serial: 0, bytes: b"ACGT".to_vec() };
        assert_eq!(RawBytesBlock::bytes(&plain), b"ACGT");
        let bgzf = BgzfBlock { batch_serial: 0, bytes: b"\x1f\x8b".to_vec(), uncompressed_size: 4 };
        assert_eq!(RawBytesBlock::bytes(&bgzf), b"\x1f\x8b");
    }

    #[test]
    fn profile_is_serial_writer_sink() {
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let step = WriteRawFile::<DecompressedBlock>::new(&path, b"").unwrap();
        let profile = step.profile();
        assert_eq!(profile.name, "WriteRawFile");
        assert_eq!(profile.kind, StepKind::Serial);
        assert!(profile.sticky);
        assert_eq!(step.affinity(), Affinity::Writer);
    }

    /// A sink dropped before the drained-finish path (an aborted stream) must
    /// NOT append its trailer — stamping the BGZF EOF onto a truncated `.gz`
    /// would hide the truncation. Mirrors
    /// `WriteBgzfFile::drop_before_finish_does_not_append_eof_marker`.
    #[test]
    fn drop_before_finish_omits_trailer() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().to_path_buf();
        let trailer = b"\x1f\x8bTRAILER";
        let step = WriteRawFile::<BgzfBlock>::new(&path, trailer).unwrap();
        // Write some payload bytes directly, then drop without draining (abort).
        {
            let mut guard = step.state.lock();
            guard.as_mut().unwrap().write_all(b"PAYLOAD").unwrap();
        }
        drop(step);

        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(bytes, b"PAYLOAD", "aborted stream must contain payload only, no trailer");
        assert!(!bytes.ends_with(trailer), "aborted stream must not end with the trailer");
    }

    #[test]
    fn dash_path_selects_stdout_without_error() {
        // `-` must construct a stdout-backed sink without touching the filesystem.
        let step = WriteRawFile::<DecompressedBlock>::new("-", b"").unwrap();
        assert!(step.state.lock().is_some());
    }
}
