//! `WriteBgzfFile` sink step. `Serial` + `Affinity::Writer`. Receives
//! pre-compressed `BgzfBlock`s from `BgzfCompress` and writes them
//! directly to disk.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use fgumi_bgzf::{BGZF_EOF, InlineBgzfCompressor};
use noodles::sam::Header;
use parking_lot::Mutex;

use crate::types::BgzfBlock;
use fgumi_pipeline_core::{
    header::HeaderHandle,
    step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile},
};

/// `Serial + sticky` BAM sink that consumes pre-compressed `BgzfBlock`s.
pub struct WriteBgzfFile {
    state: Mutex<Option<WriterState>>,
    name: &'static str,
}

struct WriterState {
    out: BufWriter<File>,
    pending_header: Option<PendingHeader>,
}

struct PendingHeader {
    handle: HeaderHandle,
    compression_level: u32,
}

impl WriteBgzfFile {
    /// Open `path`, BGZF-compress and write the BAM header bytes, return
    /// the sink ready to receive `BgzfBlock`s.
    ///
    /// # Errors
    ///
    /// Returns I/O errors from path open or header write.
    pub fn new<P: AsRef<Path>>(
        path: P,
        header: &Header,
        compression_level: u32,
    ) -> io::Result<Self> {
        let file = File::create(path.as_ref())?;
        let mut out = BufWriter::with_capacity(256 * 1024, file);

        let mut header_bytes = Vec::new();
        fgumi_bam_io::write_bam_header(&mut header_bytes, header)
            .map_err(|e| io::Error::other(format!("write_bam_header: {e}")))?;

        let mut hc = InlineBgzfCompressor::new(compression_level);
        hc.write_all(&header_bytes)?;
        hc.flush()?;
        hc.write_blocks_to(&mut out)?;

        Ok(Self {
            state: Mutex::new(Some(WriterState { out, pending_header: None })),
            name: "WriteBgzfFile",
        })
    }

    /// Open `path` and return the sink with the BAM header write
    /// deferred until an upstream step resolves `handle`.
    ///
    /// # Errors
    ///
    /// Returns I/O errors from path open. Header-write errors are
    /// surfaced from `try_run` once the handle resolves.
    pub fn new_with_handle<P: AsRef<Path>>(
        path: P,
        handle: HeaderHandle,
        compression_level: u32,
    ) -> io::Result<Self> {
        let file = File::create(path.as_ref())?;
        let out = BufWriter::with_capacity(256 * 1024, file);
        Ok(Self {
            state: Mutex::new(Some(WriterState {
                out,
                pending_header: Some(PendingHeader { handle, compression_level }),
            })),
            name: "WriteBgzfFile",
        })
    }

    fn try_write_pending_header(state: &mut WriterState) -> io::Result<bool> {
        let Some(pending) = state.pending_header.as_ref() else {
            return Ok(true);
        };
        let header_clone = match pending.handle.try_get() {
            None => return Ok(false),
            Some(Err(e)) => return Err(e),
            Some(Ok(h)) => h.clone(),
        };
        let level = pending.compression_level;

        let mut header_bytes = Vec::new();
        fgumi_bam_io::write_bam_header(&mut header_bytes, &header_clone)
            .map_err(|e| io::Error::other(format!("write_bam_header: {e}")))?;
        let mut hc = InlineBgzfCompressor::new(level);
        hc.write_all(&header_bytes)?;
        hc.flush()?;
        hc.write_blocks_to(&mut state.out)?;

        state.pending_header = None;
        Ok(true)
    }
}

impl Step for WriteBgzfFile {
    type Input = BgzfBlock;
    type Outputs = ();

    fn profile(&self) -> StepProfile {
        StepProfile {
            name: self.name,
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
        let Some(state) = guard.as_mut() else {
            return Ok(StepOutcome::Finished);
        };

        let header_ready = Self::try_write_pending_header(state)?;
        if header_ready {
            if let Some(block) = ctx.input.pop() {
                state.out.write_all(&block.bytes)?;
                return Ok(StepOutcome::Progress);
            }
        }

        if ctx.input.is_drained() {
            if !header_ready {
                return Err(io::Error::other(
                    "WriteBgzfFile: input drained before HeaderHandle was resolved",
                ));
            }
            state.out.write_all(&BGZF_EOF)?;
            state.out.flush()?;
            let _ = guard.take();
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

impl Drop for WriteBgzfFile {
    /// Cleanup-only drop. The BGZF EOF marker is written **exclusively** by the
    /// drained-finish path in `try_run` (which then takes the state so this drop
    /// is a no-op for a normally-finished sink). If state is still present here
    /// the sink was dropped before that path ran — i.e. an aborted/partial
    /// stream — so we deliberately do **not** append `BGZF_EOF`: stamping the
    /// EOF marker onto a truncated BAM would make it look like a complete stream
    /// and hide the truncation from downstream readers. We only flush whatever
    /// bytes were already buffered so the on-disk file reflects what was written
    /// (and stays detectably truncated). A still-pending header means nothing
    /// valid was written, so leave the file empty.
    fn drop(&mut self) {
        let mut guard = self.state.lock();
        if let Some(mut state) = guard.take() {
            if state.pending_header.is_some() {
                return;
            }
            let _ = state.out.flush();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn empty_header() -> Header {
        Header::default()
    }

    #[test]
    fn profile_advertises_serial_writer_sink() {
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let header = empty_header();
        let step = WriteBgzfFile::new(&path, &header, 1).unwrap();
        let profile = step.profile();
        assert_eq!(profile.name, "WriteBgzfFile");
        assert_eq!(profile.kind, StepKind::Serial);
        assert!(profile.sticky);
        assert_eq!(step.affinity(), Affinity::Writer);
        assert_eq!(profile.output_queues.len(), 0);
        assert_eq!(profile.branch_ordering.len(), 0);
    }

    #[test]
    fn header_only_round_trip() {
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let header = empty_header();
        let step = WriteBgzfFile::new(&path, &header, 1).unwrap();
        let mut guard = step.state.lock();
        let mut state = guard.take().expect("state present");
        state.out.write_all(&BGZF_EOF).unwrap();
        state.out.flush().unwrap();
        drop(guard);

        let bytes = std::fs::read(&path).unwrap();
        assert!(bytes.len() >= 28, "BGZF EOF + header should be at least 28 bytes");
        assert_eq!(&bytes[0..2], &[0x1f, 0x8b], "BGZF/gzip magic at start");
        let tail = &bytes[bytes.len() - 28..];
        assert_eq!(tail, &BGZF_EOF, "file ends with BGZF EOF marker");
    }

    #[test]
    fn new_with_handle_defers_header_until_resolved() {
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let handle = HeaderHandle::new();
        let step = WriteBgzfFile::new_with_handle(&path, handle.clone(), 1).unwrap();

        let bytes_before = std::fs::read(&path).unwrap();
        assert_eq!(bytes_before.len(), 0, "no bytes written until header resolves");

        {
            let mut guard = step.state.lock();
            let state = guard.as_mut().expect("state present");
            assert!(state.pending_header.is_some(), "handle still pending");
            let wrote = WriteBgzfFile::try_write_pending_header(state).unwrap();
            assert!(!wrote, "unresolved handle should yield without writing");
            assert!(state.pending_header.is_some(), "still pending after no-op probe");
        }
        let bytes_mid = std::fs::read(&path).unwrap();
        assert_eq!(bytes_mid.len(), 0, "still nothing on disk after no-op probe");

        handle.set(empty_header()).expect("first set");
        {
            let mut guard = step.state.lock();
            let state = guard.as_mut().expect("state present");
            let wrote = WriteBgzfFile::try_write_pending_header(state).unwrap();
            assert!(wrote, "resolved handle should write");
            assert!(state.pending_header.is_none(), "pending slot cleared");
            state.out.flush().unwrap();
        }
        let bytes_after = std::fs::read(&path).unwrap();
        assert!(bytes_after.len() >= 2, "header BGZF block emitted");
        assert_eq!(&bytes_after[0..2], &[0x1f, 0x8b], "BGZF/gzip magic at start");
    }

    #[test]
    fn new_with_handle_propagates_poison() {
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let handle = HeaderHandle::new();
        let step = WriteBgzfFile::new_with_handle(&path, handle.clone(), 1).unwrap();

        handle.poison(io::Error::new(io::ErrorKind::BrokenPipe, "aligner died")).unwrap();
        let mut guard = step.state.lock();
        let state = guard.as_mut().expect("state present");
        let err = WriteBgzfFile::try_write_pending_header(state).expect_err("poison");
        assert_eq!(err.kind(), io::ErrorKind::BrokenPipe);
        assert_eq!(err.to_string(), "aligner died");
    }

    #[test]
    fn drop_before_finish_does_not_append_eof_marker() {
        // A sink dropped before the drained-finish path in `try_run` (e.g. a
        // pipeline abort) must NOT append the BGZF EOF marker. Appending it
        // would stamp a "complete stream" signature onto a truncated BAM,
        // hiding the truncation from downstream readers. See the `Drop` doc.
        let path = tempfile::NamedTempFile::new().unwrap();
        let path_buf = path.path().to_path_buf();
        let header = empty_header();
        // `new` eagerly writes the header (pending_header is None), so the
        // only thing standing between this state and a valid EOF marker is
        // the `try_run` drained-finish path, which we never reach.
        let step = WriteBgzfFile::new(&path_buf, &header, 1).unwrap();
        drop(step);

        let bytes = std::fs::read(&path_buf).unwrap();
        assert!(bytes.len() >= 28, "header bytes should be on disk");
        let tail = &bytes[bytes.len() - 28..];
        assert_ne!(tail, &BGZF_EOF, "aborted output must not end with a valid BGZF EOF marker");
    }

    #[test]
    fn drop_with_unresolved_handle_leaves_empty_file() {
        let path = tempfile::NamedTempFile::new().unwrap();
        let path_buf = path.path().to_path_buf();
        let handle = HeaderHandle::new();
        let step = WriteBgzfFile::new_with_handle(&path_buf, handle, 1).unwrap();
        drop(step);

        let bytes = std::fs::read(&path_buf).unwrap();
        assert_eq!(bytes.len(), 0, "Drop with unresolved handle must skip EOF — see Drop doc");
    }
}
