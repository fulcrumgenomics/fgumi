//! `WriteBgzfFile` sink step. `Serial` + `Affinity::Writer`. Receives
//! pre-compressed `BgzfBlock`s from `BgzfCompress` and writes them
//! directly to disk.
//!
//! Worker `N - 1` (the last worker) is the only worker that ever
//! attempts the sink's mutex; other workers `Skip` this step in
//! dispatch. Mirrors legacy's "T(N-1) owns Write" pattern in
//! `pipeline/scheduler/mod.rs:282-286`. Worker N-1 still
//! round-robins through every other step; the priority loop in the
//! driver brings it back to the sink as upstream produces work.
//!
//! ### Header
//!
//! Two construction modes, eager and lazy.
//!
//! **Eager** ([`WriteBgzfFile::new`]) writes the SAM header bytes at
//! construction time: the bytes are framed
//! (`BAM_MAGIC` + `l_text` + text + `n_ref` + per-ref records) by
//! [`fgumi_bam_io::write_bam_header`], BGZF-compressed via an
//! `InlineBgzfCompressor`, and emitted to the output file before the
//! first `try_run`. Header bytes flush as their own BGZF block(s),
//! ensuring they're a clean prefix to the data blocks that follow.
//!
//! **Lazy** ([`WriteBgzfFile::new_with_handle`]) opens the file but
//! parks the header write until an upstream step resolves the
//! [`HeaderHandle`]. On every `try_run`, the writer first probes the
//! handle via [`HeaderHandle::try_get`]: if `None`, it returns
//! `StepOutcome::NoProgress`; if `Some(Ok(header))`, it writes the
//! header bytes the same way the eager path does and then proceeds to
//! drain `ctx.input`; if `Some(Err(e))`, it propagates the error.
//! Used by `AlignAndMergeStep` to inject the aligner's `@PG` (plus any
//! runtime `@RG`/`@CO`) into the merged output header.
//!
//! ### Trailer
//!
//! Once `ctx.input` is drained, `try_run` emits the 28-byte BGZF EOF marker,
//! flushes/closes the file, and returns `StepOutcome::Finished` (the writer
//! state is taken so a re-dispatch is a no-op). [`Drop`] remains the
//! best-effort net for shutdown paths that never reach a clean drain (a
//! sibling-step error or a panic): see the `Drop` impl.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use fgumi_bgzf::{BGZF_EOF, InlineBgzfCompressor};
use noodles::sam::Header;
use parking_lot::Mutex;

use crate::pipeline::core::header::HeaderHandle;
use crate::pipeline::core::step::{Affinity, Step, StepCtx, StepKind, StepOutcome, StepProfile};
use crate::pipeline::steps::types::BgzfBlock;

/// `Exclusive + sticky` BAM sink that consumes pre-compressed `BgzfBlock`s.
///
/// State is held behind a `Mutex<Option<...>>` because the `Step` trait's
/// uniform `Clone` bound applies. The runtime never `Clone`s `Exclusive`
/// steps; clones panic. The owning worker accesses `&mut self`.
pub struct WriteBgzfFile {
    state: Mutex<Option<WriterState>>,
    name: &'static str,
}

struct WriterState {
    /// Buffered file. In the eager-header path the header has already
    /// been written before `try_run` ever runs and `pending_header` is
    /// `None`; in the lazy path `pending_header` carries the
    /// `HeaderHandle` and the compression level used to frame the
    /// header BGZF block(s), and `try_run` consumes it as soon as an
    /// upstream step resolves the handle.
    out: BufWriter<File>,
    pending_header: Option<PendingHeader>,
}

struct PendingHeader {
    handle: HeaderHandle,
    compression_level: u32,
}

impl WriteBgzfFile {
    /// Open `path`, BGZF-compress and write the BAM header bytes, return
    /// the sink ready to receive `BgzfBlock`s. The header must be exactly
    /// what should appear in the output (caller has already applied any
    /// `@PG`-record updates).
    ///
    /// `compression_level` is used **only for the header bytes**; data
    /// `BgzfBlock`s arrive already compressed at whatever level the
    /// upstream `BgzfCompress` step was configured with.
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
        // Match the legacy sort writer's 256 KiB output BufWriter
        // (`crates/fgumi-sort/src/external.rs::sort_*` paths). The
        // default 8 KiB buffer overflows on every BGZF block write
        // (~64 KiB compressed) — 32× too small to amortize syscalls.
        let mut out = BufWriter::with_capacity(256 * 1024, file);

        // BAM header bytes (BAM_MAGIC + l_text + text + n_ref + refs).
        let mut header_bytes = Vec::new();
        fgumi_bam_io::write_bam_header(&mut header_bytes, header)
            .map_err(|e| io::Error::other(format!("write_bam_header: {e}")))?;

        // BGZF-compress the header bytes via a one-shot inline compressor
        // and emit them to the output file. This lets the header occupy
        // its own BGZF block(s); the data blocks that follow are appended
        // verbatim from upstream.
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
    /// **deferred** until an upstream step resolves `handle`.
    ///
    /// On every `try_run`, the writer probes `handle.try_get()` until
    /// it observes a resolved value. While unresolved, `try_run`
    /// returns `StepOutcome::NoProgress` without popping from
    /// `ctx.input`, so block-buffered backpressure builds upstream
    /// naturally. Once resolved, the BAM header bytes are framed and
    /// BGZF-compressed (using `compression_level`) and emitted before
    /// any data block — identical on-disk layout to the eager path.
    ///
    /// `compression_level` is used **only for the header bytes**;
    /// data `BgzfBlock`s arrive already compressed at whatever level
    /// the upstream `BgzfCompress` step was configured with.
    ///
    /// # Errors
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

    /// Consume the pending header handle if present and resolved.
    ///
    /// Returns:
    /// - `Ok(true)` if the header has already been written or was
    ///   written in this call — caller may proceed to drain input.
    /// - `Ok(false)` if the handle is still unresolved — caller
    ///   should yield without popping from input.
    /// - `Err` if the handle is poisoned or the underlying I/O fails.
    fn try_write_pending_header(state: &mut WriterState) -> io::Result<bool> {
        let Some(pending) = state.pending_header.as_ref() else {
            return Ok(true);
        };
        // Header clone is a one-time cost on the resolve transition.
        // Cloning lets us drop the borrow on `pending_header` before
        // mutating `state.out`, avoiding a borrow-checker conflict.
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
            // Worker N-1 (the Affinity::Writer target) drives the sink
            // sticky — drains the writer queue in tight bursts before
            // yielding to round-robin. Without this the writer is
            // visited only once per priority loop pass, which lets the
            // upstream BgzfCompress queue backpressure unevenly and
            // produces large run-to-run wall-time variance on the
            // 2.4M-record benchmark.
            sticky: true,
            output_queues: vec![],
            branch_ordering: vec![],
        }
    }

    fn affinity(&self) -> Affinity {
        // Worker N-1 is the dedicated writer; other workers Skip this step.
        // Matches legacy `scheduler/mod.rs:282-286` where T(N-1) owns Write.
        Affinity::Writer
    }

    fn try_run(&mut self, ctx: &mut StepCtx<'_, Self>) -> io::Result<StepOutcome> {
        let mut guard = self.state.lock();
        // State already taken (trailer written on a prior pass) — the writer
        // is closed; a re-dispatch is an idempotent no-op.
        let Some(state) = guard.as_mut() else {
            return Ok(StepOutcome::Finished);
        };

        // Lazy-header path: resolve + write the header before popping input,
        // otherwise we'd consume a block we can't yet flush in the correct
        // on-disk order. `header_ready` is false while the handle is unresolved.
        let header_ready = Self::try_write_pending_header(state)?;
        if header_ready {
            if let Some(block) = ctx.input.pop() {
                state.out.write_all(&block.bytes)?;
                return Ok(StepOutcome::Progress);
            }
        }

        // No block written this call. If upstream is drained, emit the trailer
        // and finish.
        if ctx.input.is_drained() {
            // The header must be resolved by now (drain happens after every
            // upstream is drained). If it isn't, the step that owns it never
            // set it (typically a crashed sibling) — surface explicitly rather
            // than emit a truncated header-less BAM.
            if !header_ready {
                return Err(io::Error::other(
                    "WriteBgzfFile: input drained before HeaderHandle was resolved",
                ));
            }
            state.out.write_all(&BGZF_EOF)?;
            state.out.flush()?;
            // Close the writer so `Drop` (and any re-dispatch) is a no-op.
            let _ = guard.take();
            return Ok(StepOutcome::Finished);
        }
        Ok(StepOutcome::NoProgress)
    }
}

impl Drop for WriteBgzfFile {
    /// Best-effort BGZF EOF emission for shutdown paths that bypass the clean
    /// `try_run` completion (where the trailer is normally written) — e.g., a
    /// pipeline error in another step that triggers cancellation before this
    /// sink observes its input drained, or a panic anywhere in the worker pool.
    /// Without this, a partial run leaves a header-only BAM that downstream
    /// tools (samtools, noodles) flag as truncated.
    ///
    /// **Lazy-header path:** if the writer was constructed with
    /// `new_with_handle` and the handle is still unresolved at Drop
    /// time, the BAM header bytes were never written. Emitting just
    /// the 28-byte BGZF EOF marker would produce an unrecognisable
    /// artifact (no `BAM_MAGIC`, no `@SQ` table). We deliberately
    /// skip EOF in that case: the file is left at 0 bytes — clearly
    /// "nothing was written" rather than a corrupt fragment.
    fn drop(&mut self) {
        let mut guard = self.state.lock();
        if let Some(mut state) = guard.take() {
            if state.pending_header.is_some() {
                return;
            }
            // Ignore errors: we're already on a shutdown path and there's
            // no useful place to surface them.
            let _ = state.out.write_all(&BGZF_EOF);
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
        // Write the trailer (BGZF EOF) directly, as `try_run`'s drained
        // completion would. Manual cleanup since constructing a real StepCtx in
        // unit tests is awkward; framework integration is tested in tests.rs.
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

        // Before the handle resolves the file should be empty.
        let bytes_before = std::fs::read(&path).unwrap();
        assert_eq!(bytes_before.len(), 0, "no bytes written until header resolves");

        // First try_write_pending_header probe must observe pending.
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

        // Resolve the handle and probe again — header bytes flush.
        handle.set(empty_header()).expect("first set");
        {
            let mut guard = step.state.lock();
            let state = guard.as_mut().expect("state present");
            let wrote = WriteBgzfFile::try_write_pending_header(state).unwrap();
            assert!(wrote, "resolved handle should write");
            assert!(state.pending_header.is_none(), "pending slot cleared");
            // BufWriter buffers up to 256 KiB; force a flush so the
            // on-disk inspection below sees the bytes.
            state.out.flush().unwrap();
        }
        let bytes_after = std::fs::read(&path).unwrap();
        assert!(bytes_after.len() >= 2, "header BGZF block emitted");
        assert_eq!(&bytes_after[0..2], &[0x1f, 0x8b], "BGZF/gzip magic at start");

        // Subsequent probe is a no-op.
        {
            let mut guard = step.state.lock();
            let state = guard.as_mut().expect("state present");
            let wrote = WriteBgzfFile::try_write_pending_header(state).unwrap();
            assert!(wrote, "already-resolved fast path returns true");
        }
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
    fn lazy_path_byte_identical_to_eager_via_delayed_set() {
        // Genuine parity test for the deferred path: probe returns
        // `None`, then we `set`, then we probe again, then EOF. This
        // exercises every branch of `try_write_pending_header` —
        // unlike `new_with_handle_static_header_round_trip` which
        // fast-paths through a pre-set handle.
        let header = empty_header();

        let lazy_path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let handle = HeaderHandle::new();
        let lazy_step = WriteBgzfFile::new_with_handle(&lazy_path, handle.clone(), 1).unwrap();
        {
            let mut guard = lazy_step.state.lock();
            let state = guard.as_mut().expect("state present");
            let wrote = WriteBgzfFile::try_write_pending_header(state).unwrap();
            assert!(!wrote, "first probe with unresolved handle");
        }
        handle.set(header.clone()).expect("first set");
        {
            let mut guard = lazy_step.state.lock();
            let state = guard.as_mut().expect("state present");
            let wrote = WriteBgzfFile::try_write_pending_header(state).unwrap();
            assert!(wrote, "second probe writes header");
            state.out.write_all(&BGZF_EOF).unwrap();
            state.out.flush().unwrap();
            let _ = guard.take();
        }
        let lazy_bytes = std::fs::read(&lazy_path).unwrap();

        let eager_path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let eager_step = WriteBgzfFile::new(&eager_path, &header, 1).unwrap();
        {
            let mut guard = eager_step.state.lock();
            let state = guard.as_mut().expect("state present");
            state.out.write_all(&BGZF_EOF).unwrap();
            state.out.flush().unwrap();
            let _ = guard.take();
        }
        let eager_bytes = std::fs::read(&eager_path).unwrap();

        assert_eq!(
            lazy_bytes, eager_bytes,
            "lazy + delayed-set must match eager bytes after a None→resolve→Ok probe sequence"
        );
    }

    #[test]
    fn unresolved_handle_yields_false_from_helper() {
        // The drained-before-resolve path in `try_run` is:
        // `if !try_write_pending_header(state)? { return Err(...) }`.
        // We can't easily construct a `StepCtx` in a unit test, so
        // we verify the load-bearing precondition directly: the
        // helper returns Ok(false) (no panic, no error) for an
        // unresolved handle, which is what that `!` branch keys
        // off of. Integration coverage of the full error surface
        // lives in `pipeline/steps/tests.rs`.
        let path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let handle = HeaderHandle::new();
        let step = WriteBgzfFile::new_with_handle(&path, handle, 1).unwrap();

        let mut guard = step.state.lock();
        let state = guard.as_mut().expect("state present");
        let wrote = WriteBgzfFile::try_write_pending_header(state).unwrap();
        assert!(!wrote, "unresolved handle yields false (drain converts this to Err)");
        assert!(state.pending_header.is_some(), "slot still pending");
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

    #[test]
    fn new_with_handle_static_header_round_trip() {
        // Lazy path with a handle that was pre-set via `from_header`
        // must produce a byte-identical file to the eager path. Note:
        // this is a unit test of the internal helper; it does not
        // drive `Step::try_run`. The genuine
        // None→resolve→Ok lazy-probe path is exercised by
        // `lazy_path_byte_identical_to_eager_via_delayed_set`.
        let header = empty_header();

        let lazy_path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let lazy_step = WriteBgzfFile::new_with_handle(
            &lazy_path,
            HeaderHandle::from_header(header.clone()),
            1,
        )
        .unwrap();
        {
            let mut guard = lazy_step.state.lock();
            let state = guard.as_mut().expect("state present");
            WriteBgzfFile::try_write_pending_header(state).unwrap();
            state.out.write_all(&BGZF_EOF).unwrap();
            state.out.flush().unwrap();
            let _ = guard.take();
        }
        let lazy_bytes = std::fs::read(&lazy_path).unwrap();

        let eager_path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        let eager_step = WriteBgzfFile::new(&eager_path, &header, 1).unwrap();
        {
            let mut guard = eager_step.state.lock();
            let state = guard.as_mut().expect("state present");
            state.out.write_all(&BGZF_EOF).unwrap();
            state.out.flush().unwrap();
            let _ = guard.take();
        }
        let eager_bytes = std::fs::read(&eager_path).unwrap();

        assert_eq!(lazy_bytes, eager_bytes, "lazy + pre-set handle must match eager bytes");
    }
}
