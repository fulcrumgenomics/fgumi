//! CRAM reading and writing through a vendored htslib whose CPU-heavy steps
//! run on the caller's threads.
//!
//! htslib parallelizes CRAM with its own pthread pool. This crate instead
//! exposes the three-stage shape htslib uses internally so the stages can run
//! on fgumi's pipeline workers:
//!
//! * **Read:** [`CramReader::next_slice`] (serial) yields [`SliceJob`]s;
//!   [`SliceJob::decode`] (any thread, concurrently) returns the slice's
//!   records as framed BAM bytes (`block_size` prefix included).
//! * **Write:** [`CramWriter::write_raw`] (serial) packs records into
//!   containers; full containers come back as [`EncodeJob`]s from
//!   [`CramWriter::drain_jobs`]; [`EncodeJob::encode`] (any thread,
//!   concurrently) compresses one; [`EncodedContainer::write`] (serial, in job
//!   order) writes it. [`CramWriter::seal`] flushes the final container and
//!   returns a [`CramFinisher`] that closes the file (and CRAI) once every
//!   container is written.
//!
//! The write path relies on `vendor/container-executor.patch`, which adds a
//! container-executor hook to htslib's CRAM encoder; the read path re-states
//! htslib's single-threaded slice loop in `csrc/fgumi_hts.c`.
//!
//! CRAM input older than version 2.1 is refused: it has no EOF container, so a
//! truncated file could not be told from a complete one.
//!
//! CRAM output bytes are not reproducible from run to run, even for identical
//! input and settings: htslib picks each block's codec adaptively, from trial
//! compressions. The decoded records are the same every time.
//!
//! All `unsafe` in fgumi's CRAM support lives in this crate.
#![allow(unsafe_code)]
// Every `unsafe` block states its invariant (CLAUDE.md, "Approved FFI crate").
#![warn(clippy::undocumented_unsafe_blocks)]

// Linked for the C libraries htslib calls (zlib, libdeflate, bzip2, lzma).
use bzip2_sys as _;
use libdeflate_sys as _;
use libz_sys as _;
use lzma_sys as _;

use std::ffi::{CString, c_char, c_int, c_void};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::ptr::NonNull;
use std::sync::{Arc, Mutex};

type Ptr = *mut c_void;
type ReadCb = extern "C" fn(ctx: Ptr, buf: *mut c_void, n: usize) -> isize;
type CloseCb = extern "C" fn(ctx: Ptr);
type ExecutorCb = extern "C" fn(arg: Ptr, fd: Ptr, c: Ptr) -> c_int;

unsafe extern "C" {
    fn hts_open(fname: *const c_char, mode: *const c_char) -> Ptr;
    fn hts_close(fp: Ptr) -> c_int;
    fn hts_set_fai_filename(fp: Ptr, fn_aux: *const c_char) -> c_int;
    fn sam_hdr_read(fp: Ptr) -> Ptr;
    fn sam_hdr_write(fp: Ptr, h: Ptr) -> c_int;
    fn sam_hdr_parse(len: usize, text: *const c_char) -> Ptr;
    fn sam_hdr_destroy(h: Ptr);
    fn sam_idx_init(fp: Ptr, h: Ptr, min_shift: c_int, fnidx: *const c_char) -> c_int;
    fn bam_init1() -> Ptr;
    fn bam_destroy1(b: Ptr);
    fn cram_free_container(c: Ptr);
    fn cram_free_slice(s: Ptr);

    fn fgumi_hts_open_callback_reader(
        ctx: Ptr,
        rd: ReadCb,
        cl: CloseCb,
        name: *const c_char,
        cram_version: *mut c_int,
    ) -> Ptr;
    fn fgumi_hts_cram_fd(fp: Ptr) -> Ptr;
    fn fgumi_cram_set_shared_ref(fd: Ptr);
    fn fgumi_hdr_text(h: Ptr, len: *mut usize) -> *const c_char;
    fn fgumi_cram_next_slice(fd: Ptr, cur: Ptr, c_out: *mut Ptr, s_out: *mut Ptr) -> c_int;
    fn fgumi_cram_decode_slice(
        fd: Ptr,
        c: Ptr,
        s: Ptr,
        out: *mut *mut u8,
        out_len: *mut usize,
    ) -> i64;
    fn fgumi_cram_slice_bytes(s: Ptr) -> usize;
    fn fgumi_cram_slice_decoded_size(c: Ptr, s: Ptr, out: *mut u64) -> c_int;
    fn fgumi_cram_abandon(fd: Ptr);
    fn fgumi_free(p: Ptr);
    fn fgumi_cram_writer_options(
        fp: Ptr,
        no_ref: c_int,
        use_bzip2: c_int,
        use_lzma: c_int,
    ) -> c_int;
    fn fgumi_cram_set_file_id(fd: Ptr, id: *const c_char);
    fn fgumi_hts_write_raw(fp: Ptr, h: Ptr, scratch: Ptr, rec: *const u8, len: usize) -> c_int;

    // The container-executor API (vendor/container-executor.patch, proposed
    // upstream as samtools/htslib's `cram_set_container_executor`).
    fn cram_set_container_executor(fd: Ptr, f: Option<ExecutorCb>, arg: Ptr) -> c_int;
    fn cram_encode_container(fd: Ptr, c: Ptr) -> c_int;
    fn cram_flush_encoded_container(fd: Ptr, c: Ptr) -> c_int;
    fn cram_flush(fd: Ptr) -> c_int;
}

/// `s` as a C string; `what` names it in the error if it holds a NUL byte.
fn cstr(s: &str, what: &str) -> io::Result<CString> {
    CString::new(s).map_err(|_| {
        io::Error::new(io::ErrorKind::InvalidInput, format!("{what} contains a NUL byte"))
    })
}

fn path_cstr(p: &Path) -> io::Result<CString> {
    use std::os::unix::ffi::OsStrExt;
    CString::new(p.as_os_str().as_bytes()).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("path {} contains a NUL byte", p.display()),
        )
    })
}

/// A parsed SAM header owned by Rust (destroyed on drop).
struct OwnedHeader(Ptr);

impl Drop for OwnedHeader {
    fn drop(&mut self) {
        // SAFETY: the header is ours (from `sam_hdr_parse`); htslib keeps its
        // own copy of any header it was given, so nothing else points here.
        unsafe { sam_hdr_destroy(self.0) };
    }
}

fn hts_err(what: &str) -> io::Error {
    io::Error::other(format!("htslib: {what} failed (see htslib messages above)"))
}

/// Header text of an htslib header handle.
fn header_text(h: Ptr) -> io::Result<String> {
    let mut len = 0usize;
    // SAFETY: `h` is a live header; htslib returns a pointer to `len` bytes it
    // owns, valid until the header is modified or destroyed (not before we copy).
    let p = unsafe { fgumi_hdr_text(h, &raw mut len) };
    if p.is_null() {
        return Ok(String::new());
    }
    // SAFETY: as above.
    let bytes = unsafe { std::slice::from_raw_parts(p.cast::<u8>(), len) };
    String::from_utf8(bytes.to_vec())
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "CRAM header is not UTF-8"))
}

// ---------------------------------------------------------------------------
// Reading

/// The open file. Shared by the reader and every in-flight [`SliceJob`]; the
/// file is closed when the last of them drops.
struct ReaderInner {
    fp: Ptr,
    fd: Ptr,
}

// SAFETY: htslib's CRAM decoder is designed for slice decodes to run on pool
// threads concurrently with the main thread reading further containers (its
// shared state is guarded by `ref_lock`/`metrics_lock`/`bam_list_lock`). The
// serial operations (`fgumi_cram_next_slice`) are only reachable through
// `&mut CramReader`, so they never run concurrently with each other.
unsafe impl Send for ReaderInner {}
// SAFETY: see above.
unsafe impl Sync for ReaderInner {}

impl Drop for ReaderInner {
    fn drop(&mut self) {
        // SAFETY: `fp` was opened by htslib and every container/slice that
        // referenced it has been freed (they hold an `Arc` to this).
        unsafe { hts_close(self.fp) };
    }
}

/// One container: freed once the reader has moved past it and every slice job
/// that references it has finished.
struct Container {
    ptr: Ptr,
    reader: Arc<ReaderInner>,
}
// SAFETY: slice decodes share the parsed compression header read-only. The
// vendored decoder guards its mutable `c->unsorted` state with `ref_lock`;
// serial slice iteration uses separate fields. The container is freed once.
unsafe impl Send for Container {}
// SAFETY: see above.
unsafe impl Sync for Container {}
impl Drop for Container {
    fn drop(&mut self) {
        // SAFETY: last reference; no slice of it is still decoding. `ptr` is
        // the non-null container `fgumi_cram_next_slice` returned. The file
        // (`reader`) is released after this body, so it outlives the container.
        unsafe { cram_free_container(self.ptr) };
    }
}

/// A CRAM reader. [`CramReader::next_slice`] is the serial step.
pub struct CramReader {
    inner: Arc<ReaderInner>,
    current: Option<Arc<Container>>,
    header_text: String,
    serial: u64,
    finished: bool,
}

/// One slice to decode: `Send`, decodable on any thread.
pub struct SliceJob {
    slice: Ptr,
    container: Arc<Container>,
    serial: u64,
    compressed_bytes: u64,
}
// SAFETY: a slice is owned by exactly one job and decoded by whichever thread
// holds it (htslib's own pool does the same).
unsafe impl Send for SliceJob {}

impl Drop for SliceJob {
    fn drop(&mut self) {
        // SAFETY: the job owns its slice (non-null, from
        // `fgumi_cram_next_slice`). The container (an `Arc` field, which in
        // turn holds the file) is dropped after this body, so the slice is
        // freed before its container, and the container before the file.
        unsafe { cram_free_slice(self.slice) };
    }
}

/// A decoded slice: its records as framed BAM (each with its `block_size`).
pub struct DecodedSlice {
    pub serial: u64,
    pub records: u64,
    pub bytes: Vec<u8>,
}

impl SliceJob {
    /// 0-based position of this slice in the file (consecutive).
    #[must_use]
    pub fn serial(&self) -> u64 {
        self.serial
    }

    /// Compressed bytes this undecoded slice holds (for queue accounting).
    #[must_use]
    pub fn compressed_bytes(&self) -> u64 {
        self.compressed_bytes
    }

    /// Practical metadata-based estimate for admitting a slice before decode.
    ///
    /// Charges the container's bases, uncompressed blocks, record overhead,
    /// repeated tag labels and scratch/copy slack. The regenerated-MD reserve
    /// is capped at four times the bases and blocks, so ordinary mapped slices
    /// can still decode concurrently. Exceptional repeated long deletions can
    /// exceed it. This is not a hard allocation or RSS bound: codecs and shared
    /// reference caches can also allocate independently.
    ///
    /// # Errors
    ///
    /// Returns `InvalidData` for negative sizes or checked-arithmetic overflow.
    pub fn decoded_size_estimate(&self) -> io::Result<u64> {
        let mut bytes = 0;
        // SAFETY: the job holds its container and owns its slice; the estimator
        // reads immutable header metadata and writes only this out-parameter.
        if unsafe { fgumi_cram_slice_decoded_size(self.container.ptr, self.slice, &raw mut bytes) }
            != 0
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid CRAM decoded size estimate",
            ));
        }
        Ok(bytes)
    }

    /// Decode the slice to framed BAM records.
    ///
    /// # Errors
    ///
    /// Returns an error if htslib fails to decode the slice (corrupt data, or a
    /// missing / mismatched reference).
    pub fn decode(self) -> io::Result<DecodedSlice> {
        let mut out: *mut u8 = std::ptr::null_mut();
        let mut len = 0usize;
        // SAFETY: `fd`, the container and the slice are live (held by `self`);
        // concurrent decodes of distinct slices are what htslib's pool does.
        let n = unsafe {
            fgumi_cram_decode_slice(
                self.container.reader.fd,
                self.container.ptr,
                self.slice,
                &raw mut out,
                &raw mut len,
            )
        };
        if n < 0 {
            return Err(hts_err("CRAM slice decode"));
        }
        let bytes = if out.is_null() {
            Vec::new()
        } else {
            // SAFETY: `out` is a malloc'd buffer of `len` initialized bytes that
            // we own; copied, then freed with the matching allocator.
            let v = unsafe { std::slice::from_raw_parts(out, len) }.to_vec();
            // SAFETY: `out` came from the glue's `malloc`/`realloc` and is not
            // used again; `fgumi_free` is the matching `free`.
            unsafe { fgumi_free(out.cast()) };
            v
        };
        Ok(DecodedSlice { serial: self.serial, records: n.unsigned_abs(), bytes })
    }
}

extern "C" fn read_cb(ctx: Ptr, buf: *mut c_void, n: usize) -> isize {
    // SAFETY: `ctx` is the `Box<Box<dyn Read + Send>>` leaked in
    // `CramReader::from_reader`; htslib calls back only from the thread that
    // is driving the reader, never concurrently.
    let reader = unsafe { &mut *ctx.cast::<Box<dyn Read + Send>>() };
    // SAFETY: htslib passes a writable buffer of `n` bytes. Initialize it
    // before creating a slice: safe `Read` implementations may inspect it.
    unsafe { std::ptr::write_bytes(buf.cast::<u8>(), 0, n) };
    // SAFETY: the exclusive native buffer now contains `n` initialized bytes.
    let buf = unsafe { std::slice::from_raw_parts_mut(buf.cast::<u8>(), n) };
    loop {
        match reader.read(buf) {
            Ok(k) => match isize::try_from(k) {
                Ok(count) if k <= n => return count,
                _ => {
                    log::error!("CRAM input read returned an invalid count: {k} for {n} bytes");
                    return -1;
                }
            },
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => {
                log::error!("CRAM input read failed: {e}");
                return -1;
            }
        }
    }
}

extern "C" fn close_cb(ctx: Ptr) {
    // SAFETY: reclaims the box leaked in `from_reader`; htslib calls this once.
    drop(unsafe { Box::from_raw(ctx.cast::<Box<dyn Read + Send>>()) });
}

impl CramReader {
    /// Open a CRAM stream from any reader (file, stdin, FIFO, a stream whose
    /// first bytes were already sniffed and replayed). `name` is used in
    /// htslib messages. `reference` is the FASTA used to decode; without it
    /// htslib falls back to `REF_PATH`/`REF_CACHE` and the `@SQ UR`/`M5` tags.
    ///
    /// # Errors
    ///
    /// Returns an error if the stream is not CRAM, is CRAM older than version
    /// 2.1 (which has no EOF container, so truncation cannot be detected), or
    /// its header is unreadable.
    pub fn from_reader(
        reader: Box<dyn Read + Send>,
        name: &str,
        reference: Option<&Path>,
    ) -> io::Result<Self> {
        let name_c = cstr(name, "the CRAM input name")?;
        let ctx = Box::into_raw(Box::new(reader)).cast::<c_void>();
        let mut version: c_int = 0;
        // SAFETY: on failure htslib closes the backend, which frees `ctx`;
        // `version` is a live out-parameter.
        let fp = unsafe {
            fgumi_hts_open_callback_reader(
                ctx,
                read_cb,
                close_cb,
                name_c.as_ptr(),
                &raw mut version,
            )
        };
        if fp.is_null() {
            if version > 0 && version < (2 << 8 | 1) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "{name}: CRAM version {}.{} is not supported: versions before 2.1 \
                         have no EOF container, so a truncated file cannot be told from a \
                         complete one; re-encode it with `samtools view -C` (or convert it \
                         to BAM)",
                        version >> 8,
                        version & 0xff
                    ),
                ));
            }
            return Err(hts_err("opening CRAM input"));
        }
        // SAFETY: `fp` is a live file; null when it is not CRAM.
        let fd = unsafe { fgumi_hts_cram_fd(fp) };
        if fd.is_null() {
            // SAFETY: `fp` is ours and unused; closing it frees the context.
            unsafe { hts_close(fp) };
            return Err(io::Error::new(io::ErrorKind::InvalidData, "input is not CRAM"));
        }
        // SAFETY: `fd` is the live CRAM handle of `fp`, before any slice is
        // read or decoded. Slices decode concurrently, so the reference must
        // be shared whole (as htslib does under its own pool).
        unsafe { fgumi_cram_set_shared_ref(fd) };
        let inner = Arc::new(ReaderInner { fp, fd });
        if let Some(r) = reference {
            let r = path_cstr(r)?;
            // SAFETY: live file, NUL-terminated path.
            if unsafe { hts_set_fai_filename(inner.fp, r.as_ptr()) } != 0 {
                return Err(hts_err("setting the CRAM reference"));
            }
        }
        // SAFETY: live file; we own the returned copy and destroy it.
        let h = unsafe { sam_hdr_read(inner.fp) };
        if h.is_null() {
            return Err(hts_err("reading the CRAM header"));
        }
        let text = header_text(h);
        // SAFETY: `h` is ours.
        unsafe { sam_hdr_destroy(h) };
        Ok(Self { inner, current: None, header_text: text?, serial: 0, finished: false })
    }

    /// Open a CRAM file by path.
    ///
    /// # Errors
    ///
    /// As [`CramReader::from_reader`], plus failing to open `path`.
    // A test / diagnostic helper (fgumi itself opens through `from_reader`).
    // Hidden rather than `cfg(test)`: integration tests link the crate built
    // without `cfg(test)`.
    #[doc(hidden)]
    pub fn open(path: &Path, reference: Option<&Path>) -> io::Result<Self> {
        let f = std::fs::File::open(path)?;
        Self::from_reader(Box::new(f), &path.display().to_string(), reference)
    }

    /// The SAM header text.
    #[must_use]
    pub fn header_text(&self) -> &str {
        &self.header_text
    }

    /// The next slice to decode, or `None` at end of file. Serial.
    ///
    /// # Errors
    ///
    /// Returns an error if a container or slice cannot be read, including a
    /// stream that ends without its EOF container or has stray bytes after
    /// it (both `InvalidData`).
    pub fn next_slice(&mut self) -> io::Result<Option<SliceJob>> {
        if self.finished {
            return Ok(None);
        }
        let cur = self.current.as_ref().map_or(std::ptr::null_mut(), |c| c.ptr);
        let mut c_out: Ptr = std::ptr::null_mut();
        let mut s_out: Ptr = std::ptr::null_mut();
        // SAFETY: serial (`&mut self`); `cur` is live (we hold an `Arc`).
        let r =
            unsafe { fgumi_cram_next_slice(self.inner.fd, cur, &raw mut c_out, &raw mut s_out) };
        match r {
            0 => {
                self.current = None;
                self.finished = true;
                Ok(None)
            }
            1 => {
                let container = match &self.current {
                    Some(c) if c.ptr == c_out => Arc::clone(c),
                    _ => {
                        // A newly read container: we own it now.
                        let c = Arc::new(Container { ptr: c_out, reader: Arc::clone(&self.inner) });
                        self.current = Some(Arc::clone(&c));
                        c
                    }
                };
                let serial = self.serial;
                self.serial += 1;
                // SAFETY: `s_out` is the live slice just read.
                let compressed_bytes = unsafe { fgumi_cram_slice_bytes(s_out) } as u64;
                Ok(Some(SliceJob { slice: s_out, container, serial, compressed_bytes }))
            }
            // Not `UnexpectedEof`: readers layered on this (e.g. a BGZF
            // reader) treat that kind as a clean end of stream.
            -2 => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CRAM stream ends without its EOF container (truncated)",
            )),
            -3 => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "CRAM stream has stray bytes after its EOF container",
            )),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "htslib: reading the next CRAM container or slice failed: the stream is \
                 corrupt, cut short inside a container, or unreadable (see htslib messages \
                 above)",
            )),
        }
    }
}

#[cfg(test)]
mod reader_callback_tests {
    use super::*;
    use rstest::rstest;

    struct InspectBuffer;
    impl Read for InspectBuffer {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            if buf.iter().any(|&byte| byte != 0) {
                return Err(io::Error::other("Read received nonzero scratch bytes"));
            }
            buf[0] = 42;
            Ok(1)
        }
    }

    #[test]
    fn callback_initializes_buffer_before_read() {
        let mut reader: Box<dyn Read + Send> = Box::new(InspectBuffer);
        let mut bytes = [0x55u8; 4];
        let count = read_cb((&raw mut reader).cast(), bytes.as_mut_ptr().cast(), bytes.len());
        assert_eq!(count, 1);
        assert_eq!(bytes[0], 42);
    }

    struct CountReader(usize);
    impl Read for CountReader {
        fn read(&mut self, _buf: &mut [u8]) -> io::Result<usize> {
            Ok(self.0)
        }
    }

    #[rstest]
    #[case::short_read(1, 1)]
    #[case::exact_read(4, 4)]
    #[case::eof(0, 0)]
    #[case::oversized_read(5, -1)]
    #[case::unrepresentable_count(usize::MAX, -1)]
    fn callback_validates_read_count(#[case] count: usize, #[case] expected: isize) {
        let mut reader: Box<dyn Read + Send> = Box::new(CountReader(count));
        let mut bytes = [0u8; 4];
        assert_eq!(
            read_cb((&raw mut reader).cast(), bytes.as_mut_ptr().cast(), bytes.len()),
            expected
        );
    }
}

// ---------------------------------------------------------------------------
// Writing

/// The open output file, shared by the writer, every in-flight job, and the
/// finisher.
struct WriterInner {
    fp: Ptr,
    fd: Ptr,
    /// Serializes `cram_flush_encoded_container` and checks job order.
    write: Mutex<WriteState>,
    /// The CRAI successfully opened by this writer, removed on failed completion.
    crai: Option<PathBuf>,
}
struct WriteState {
    next: u64,
    closed: bool,
}
// SAFETY: htslib runs container encodes on pool threads concurrently with the
// producer packing further records, and writes containers from one thread at
// a time; this type reproduces exactly that: encodes are free-threaded,
// writes go through `write` (a mutex), and record packing is reachable only
// through `&mut CramWriter`. The writer and the producer touch disjoint
// `cram_fd` state (the writer: `fp`, `idxfp`, read-only `version`/`vv`).
unsafe impl Send for WriterInner {}
// SAFETY: see above.
unsafe impl Sync for WriterInner {}

impl WriterInner {
    /// Remove only the index this writer successfully opened.
    fn remove_index(&self) {
        if let Some(crai) = &self.crai {
            let _ = std::fs::remove_file(crai);
        }
    }
}

impl Drop for WriterInner {
    fn drop(&mut self) {
        let st = self.write.get_mut().unwrap_or_else(std::sync::PoisonError::into_inner);
        if !st.closed {
            // An abandoned write (an error path). Drop the open container and
            // mark the file incomplete, so htslib writes no EOF block: a
            // truncated CRAM must not read as complete (as the BAM sink never
            // stamps an EOF marker on an abandoned stream). The CRAI htslib
            // finishes on close would index an incomplete file, so remove it.
            // SAFETY: last reference to the file; no job, container or hook
            // can reach it any more.
            unsafe {
                fgumi_cram_abandon(self.fd);
                hts_close(self.fp);
            }
            self.remove_index();
        }
    }
}

/// Containers handed over by htslib's container executor, in file order.
#[derive(Default)]
struct DispatchSink {
    containers: Vec<Ptr>,
}

extern "C" fn dispatch_cb(ctx: Ptr, _fd: Ptr, c: Ptr) -> c_int {
    // SAFETY: `ctx` is the writer's `DispatchSink`, owned through a raw
    // pointer (`CramWriter::sink`) for its whole life and never through a
    // reference that outlives one access, so this pointer's provenance stays
    // valid. htslib calls the executor only from inside
    // `sam_write1`/`cram_flush`, which run on the thread holding
    // `&mut CramWriter`, so no other access to the sink is live.
    let sink = unsafe { &mut *ctx.cast::<DispatchSink>() };
    sink.containers.push(c);
    0
}

/// A write-side container, owned by exactly one [`EncodeJob`] or
/// [`EncodedContainer`] at a time, with the file it belongs to. Freed on drop
/// unless [`OwnedContainer::take`] handed it to htslib.
struct OwnedContainer {
    ptr: Ptr,
    inner: Arc<WriterInner>,
    serial: u64,
}

// SAFETY: one owner per container, and htslib's own pool encodes containers
// on whichever thread picks them up; the file it points into is `Send + Sync`
// (see `WriterInner`).
unsafe impl Send for OwnedContainer {}

impl OwnedContainer {
    /// Give up ownership of the container (to an htslib call that frees it).
    fn take(&mut self) -> Ptr {
        std::mem::replace(&mut self.ptr, std::ptr::null_mut())
    }
}

impl Drop for OwnedContainer {
    fn drop(&mut self) {
        if !self.ptr.is_null() {
            // SAFETY: never written (a write takes the pointer first), so the
            // container is still ours and freed exactly once, here.
            unsafe { cram_free_container(self.ptr) };
        }
    }
}

/// A container to encode: `Send`, encodable on any thread.
pub struct EncodeJob {
    container: OwnedContainer,
    approx_bytes: u64,
}

/// An encoded container awaiting its in-order write.
pub struct EncodedContainer {
    container: OwnedContainer,
}

impl EncodeJob {
    /// The container's position in write order (0-based, consecutive).
    #[must_use]
    pub fn serial(&self) -> u64 {
        self.container.serial
    }

    /// Approximate raw (BAM-encoded) bytes of the container's records, for
    /// queue accounting.
    #[must_use]
    pub fn approx_bytes(&self) -> u64 {
        self.approx_bytes
    }

    /// Compress the container. Runs concurrently with other encodes and with
    /// the producer.
    ///
    /// # Errors
    ///
    /// Returns an error if htslib fails to encode (e.g. a reference mismatch).
    pub fn encode(self) -> io::Result<EncodedContainer> {
        let container = self.container;
        // SAFETY: we own the container; `fd` is live (held by `inner`). On
        // failure the container is still ours, and dropping it frees it.
        if unsafe { cram_encode_container(container.inner.fd, container.ptr) } != 0 {
            return Err(hts_err("CRAM container encode"));
        }
        Ok(EncodedContainer { container })
    }
}

impl EncodedContainer {
    /// The container's position in write order (0-based, consecutive).
    #[must_use]
    pub fn serial(&self) -> u64 {
        self.container.serial
    }

    /// Write the container to the file (and its CRAI entries), then free it.
    ///
    /// # Errors
    ///
    /// Returns an error on an out-of-order write or an I/O failure.
    pub fn write(mut self) -> io::Result<()> {
        let inner = Arc::clone(&self.container.inner);
        let mut st = inner.write.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        if st.closed || st.next != self.container.serial {
            return Err(io::Error::other(format!(
                "CRAM container {} written out of order (expected {})",
                self.container.serial, st.next
            )));
        }
        let ptr = self.container.take();
        // SAFETY: serialized by the mutex, in dispatch order; frees `ptr`.
        let rc = unsafe { cram_flush_encoded_container(inner.fd, ptr) };
        st.next += 1;
        if rc != 0 {
            return Err(hts_err("CRAM container write"));
        }
        Ok(())
    }
}

/// Closes the output once every container has been written.
pub struct CramFinisher {
    inner: Arc<WriterInner>,
    total: u64,
}

impl CramFinisher {
    /// Number of containers the writer dispatched in total.
    // A test / diagnostic helper, hidden rather than `cfg(test)`: integration
    // tests link the crate built without `cfg(test)`.
    #[doc(hidden)]
    #[must_use]
    pub fn total_containers(&self) -> u64 {
        self.total
    }

    /// Write the EOF container, finish the CRAI (if any) and close the file.
    ///
    /// # Errors
    ///
    /// Returns an error if containers are still unwritten or close fails.
    pub fn finish(self) -> io::Result<()> {
        let mut st = self.inner.write.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        if st.next != self.total {
            return Err(io::Error::other(format!(
                "CRAM close with {} of {} containers written",
                st.next, self.total
            )));
        }
        st.closed = true;
        // SAFETY: every container is written and the writer (the executor's
        // context) is gone; nothing else uses the file now. For CRAM,
        // `hts_close` writes the EOF block and finishes the CRAI itself.
        if unsafe { hts_close(self.inner.fp) } != 0 {
            self.inner.remove_index();
            return Err(hts_err("closing the CRAM output"));
        }
        Ok(())
    }
}

/// A CRAM writer whose container encodes run on the caller's threads.
pub struct CramWriter {
    inner: Arc<WriterInner>,
    hdr: OwnedHeader,
    scratch: Ptr,
    /// The executor's context, from `Box::into_raw`: accessed only through
    /// this pointer (never a `Box`/reference that would invalidate the copy
    /// htslib holds) and freed in `Drop` after the executor is detached.
    sink: NonNull<DispatchSink>,
    dispatched: u64,
    bytes_in_open: u64,
}
// SAFETY: the writer is used by one thread at a time (`&mut self`); moving it
// between threads is fine (htslib holds no thread affinity).
unsafe impl Send for CramWriter {}

impl Drop for CramWriter {
    fn drop(&mut self) {
        // The executor's context is `self.sink`, which dies with us while the
        // file (`WriterInner`) may live on in jobs or the finisher: detach it
        // first.
        // SAFETY: serial (`&mut self`); htslib only calls the executor from
        // `sam_write1`/`cram_flush`, which need `&mut CramWriter`.
        unsafe { cram_set_container_executor(self.inner.fd, None, std::ptr::null_mut()) };
        // SAFETY: the sink came from `Box::into_raw` in `create` and is freed
        // only here; the executor (htslib's copy of the pointer) was detached
        // above, so nothing can reach it any more.
        let sink = unsafe { Box::from_raw(self.sink.as_ptr()) };
        // Discard containers never taken as jobs.
        for c in sink.containers {
            // SAFETY: owned by us, never dispatched to a job.
            unsafe { cram_free_container(c) };
        }
        // SAFETY: ours (null if `create` failed to allocate it; accepted).
        unsafe { bam_destroy1(self.scratch) };
    }
}

/// The iterator [`CramWriter::drain_jobs`] returns.
struct DrainJobs<'a> {
    containers: std::vec::Drain<'a, Ptr>,
    inner: &'a Arc<WriterInner>,
    dispatched: &'a mut u64,
    approx_bytes: u64,
}

impl Iterator for DrainJobs<'_> {
    type Item = EncodeJob;

    fn next(&mut self) -> Option<EncodeJob> {
        let ptr = self.containers.next()?;
        let serial = *self.dispatched;
        *self.dispatched += 1;
        let container = OwnedContainer { ptr, inner: Arc::clone(self.inner), serial };
        Some(EncodeJob { container, approx_bytes: self.approx_bytes })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.containers.size_hint()
    }
}

impl ExactSizeIterator for DrainJobs<'_> {}

impl Drop for DrainJobs<'_> {
    fn drop(&mut self) {
        // Free the containers not taken (each job's drop frees its container).
        self.by_ref().for_each(drop);
    }
}

/// How to create a CRAM output (see [`CramWriter::create`]).
#[derive(Debug, Clone, Copy, Default)]
pub struct CramWriteOptions<'a> {
    /// Reference FASTA to encode mapped reads against. Without one, htslib
    /// looks the reference up itself (`REF_PATH`/`REF_CACHE`, `@SQ UR:`) and,
    /// failing that, embeds the reference bases a container needs.
    pub reference: Option<&'a Path>,
    /// Store bases verbatim instead of against any reference (htslib's
    /// `no_ref`): valid for any reads and what unmapped data needs anyway, but
    /// mapped reads compress worse. Ignored when `reference` is given.
    pub no_ref: bool,
    /// Write a CRAI here as the file is written (coordinate-sorted output only).
    pub crai: Option<&'a Path>,
    /// The file-definition ID (at most 20 bytes are kept); htslib otherwise
    /// uses the path the file was opened under.
    pub file_id: Option<&'a str>,
    /// Allow htslib's bzip2 block codec (slower, sometimes smaller).
    pub use_bzip2: bool,
    /// Allow htslib's lzma block codec (slower, sometimes smaller).
    pub use_lzma: bool,
}

impl CramWriter {
    /// Create `path` (`-` for stdout) as CRAM with `header_text`, as `opts`
    /// says (reference, CRAI, file ID, codecs).
    ///
    /// # Errors
    ///
    /// Returns an error if the file, reference or header cannot be set up.
    pub fn create(path: &Path, header_text: &str, opts: &CramWriteOptions<'_>) -> io::Result<Self> {
        // Validate every input before `hts_open` creates (truncates) `path`,
        // so a bad reference, header or name never destroys an existing file.
        let p = path_cstr(path)?;
        let mode = cstr("wc", "the htslib mode")?;
        let text = cstr(header_text, "the CRAM output header text")?;
        let file_id = opts.file_id.map(|id| cstr(id, "the CRAM file ID")).transpose()?;
        let crai = opts.crai.map(path_cstr).transpose()?;
        let reference = match opts.reference {
            Some(r) => {
                // Existence and readability only: htslib builds a missing
                // `.fai` itself.
                std::fs::File::open(r).map_err(|e| {
                    io::Error::new(
                        e.kind(),
                        format!("cannot read the CRAM reference {}: {e}", r.display()),
                    )
                })?;
                Some(path_cstr(r)?)
            }
            None => None,
        };
        // SAFETY: NUL-terminated text of the given length.
        let hdr = unsafe { sam_hdr_parse(header_text.len(), text.as_ptr()) };
        if hdr.is_null() {
            return Err(hts_err("parsing the output header"));
        }
        let hdr = OwnedHeader(hdr);
        // SAFETY: NUL-terminated strings.
        let fp = unsafe { hts_open(p.as_ptr(), mode.as_ptr()) };
        if fp.is_null() {
            return Err(io::Error::other(format!("cannot create CRAM output {}", path.display())));
        }
        // SAFETY: `fp` is a live file.
        let fd = unsafe { fgumi_hts_cram_fd(fp) };
        if fd.is_null() {
            // SAFETY: `fp` is ours and unused.
            unsafe { hts_close(fp) };
            return Err(io::Error::other("htslib did not open the output as CRAM"));
        }
        // SAFETY: live CRAM file. Containers encode concurrently, so the
        // reference must be shared whole (as htslib does under its own pool).
        unsafe { fgumi_cram_set_shared_ref(fd) };
        // From here on `WriterInner`'s drop abandons the file on any error.
        let inner = Arc::new(WriterInner {
            fp,
            fd,
            write: Mutex::new(WriteState { next: 0, closed: false }),
            crai: None,
        });
        // SAFETY: allocation only; checked below, freed in Drop.
        let scratch = unsafe { bam_init1() };
        // SAFETY: `Box::into_raw` never returns null.
        let sink = unsafe { NonNull::new_unchecked(Box::into_raw(Box::<DispatchSink>::default())) };
        // From here on `CramWriter`'s drop frees the sink, scratch and header.
        let mut w = Self { inner, hdr, scratch, sink, dispatched: 0, bytes_in_open: 0 };
        if w.scratch.is_null() {
            return Err(io::Error::other("htslib could not allocate a record"));
        }
        // SAFETY: live file, before the header is written.
        let rc = unsafe {
            fgumi_cram_writer_options(
                fp,
                c_int::from(opts.no_ref && opts.reference.is_none()),
                c_int::from(opts.use_bzip2),
                c_int::from(opts.use_lzma),
            )
        };
        if rc != 0 {
            return Err(hts_err("setting CRAM output options"));
        }
        if let Some(id) = &file_id {
            // SAFETY: live file, before the header (which writes the file
            // definition); NUL-terminated id.
            unsafe { fgumi_cram_set_file_id(fd, id.as_ptr()) };
        }
        if let Some(r) = &reference {
            // SAFETY: live file, NUL-terminated path.
            if unsafe { hts_set_fai_filename(fp, r.as_ptr()) } != 0 {
                return Err(hts_err("setting the CRAM reference"));
            }
        }
        // Install the executor first, as htslib documents: before the header.
        let ctx = w.sink.as_ptr().cast::<c_void>();
        // SAFETY: the sink is a raw heap allocation that lives (at a stable
        // address) until `CramWriter::drop`, which detaches the executor before
        // freeing it; it is only ever accessed through this same pointer.
        if unsafe { cram_set_container_executor(w.inner.fd, Some(dispatch_cb), ctx) } != 0 {
            return Err(hts_err("installing the CRAM container executor"));
        }
        // SAFETY: live file and header.
        if unsafe { sam_hdr_write(w.inner.fp, w.hdr.0) } != 0 {
            return Err(hts_err("writing the CRAM header"));
        }
        if let Some(c) = &crai {
            // SAFETY: live file/header; must precede the first record.
            if unsafe { sam_idx_init(w.inner.fp, w.hdr.0, 0, c.as_ptr()) } != 0 {
                return Err(hts_err("starting the CRAI index"));
            }
            // No jobs can hold the file before `create` returns.
            Arc::get_mut(&mut w.inner)
                .ok_or_else(|| io::Error::other("CRAM writer shared before index setup completed"))?
                .crai = opts.crai.map(Path::to_path_buf);
        }
        Ok(w)
    }

    /// Add one raw BAM record (WITHOUT its 4-byte `block_size` prefix).
    ///
    /// # Errors
    ///
    /// Returns an error if the record is malformed (including a reference id
    /// outside the header, or a CIGAR whose query length is not the sequence
    /// length) or htslib rejects it.
    pub fn write_raw(&mut self, record: &[u8]) -> io::Result<()> {
        // SAFETY: serial (`&mut self`); `record` is a valid slice.
        let rc = unsafe {
            fgumi_hts_write_raw(
                self.inner.fp,
                self.hdr.0,
                self.scratch,
                record.as_ptr(),
                record.len(),
            )
        };
        if rc != 0 {
            let why = match rc {
                -5 => "a reference id outside the header",
                -6 => "a CIGAR whose query length differs from the sequence length",
                -4 => "a malformed BAM record",
                _ => "a record htslib rejected",
            };
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("cannot write CRAM: {why}"),
            ));
        }
        self.bytes_in_open += record.len() as u64;
        Ok(())
    }

    /// Containers filled since the last call, as jobs, in write order.
    ///
    /// Dropping the iterator early still turns every remaining container into
    /// a job and drops it, so none leaks; those containers are then never
    /// written, and [`CramFinisher::finish`] fails.
    pub fn drain_jobs(&mut self) -> impl Iterator<Item = EncodeJob> + '_ {
        // SAFETY: the sink lives until `Drop` and is reached only through
        // `self.sink`'s pointer; `&mut self` (held by the returned iterator)
        // means no FFI call (and so no `dispatch_cb`) is in progress while it
        // lives, so this is the only live access.
        let sink = unsafe { &mut *self.sink.as_ptr() };
        // Approximate per-container size: the raw bytes packed since the last
        // drain, split evenly.
        let n = sink.containers.len() as u64;
        let per = if n == 0 { 0 } else { std::mem::take(&mut self.bytes_in_open) / n };
        DrainJobs {
            containers: sink.containers.drain(..),
            inner: &self.inner,
            dispatched: &mut self.dispatched,
            approx_bytes: per.max(1),
        }
    }

    /// Flush the final partial container. Returns the remaining jobs and the
    /// finisher that closes the file once every job is written.
    ///
    /// # Errors
    ///
    /// Returns an error if htslib fails to flush.
    pub fn seal(mut self) -> io::Result<(Vec<EncodeJob>, CramFinisher)> {
        // SAFETY: serial; with an executor installed, `cram_flush` hands the
        // final partial container to our sink instead of writing it.
        if unsafe { cram_flush(self.inner.fd) } != 0 {
            return Err(hts_err("flushing the final CRAM container"));
        }
        let jobs = self.drain_jobs().collect();
        let finisher = CramFinisher { inner: Arc::clone(&self.inner), total: self.dispatched };
        Ok((jobs, finisher))
    }
}
