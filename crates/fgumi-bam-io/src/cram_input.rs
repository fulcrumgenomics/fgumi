//! CRAM input, decoded through htslib at the normalization boundary (feature
//! `cram`).
//!
//! A CRAM stream is presented as the same BGZF BAM byte stream every fgumi
//! reader consumes (see [`crate::sam_input`]): the header is re-encoded as a
//! BAM header and each CRAM slice is decoded straight to framed BAM records,
//! then BGZF-framed uncompressed. Nothing downstream knows the input was CRAM.
//!
//! Pipeline (chain) inputs do not come through here: they are read slice by
//! slice and decoded on the pipeline's own workers (`ReadCramSlices` /
//! `DecodeCramSlices` in `fgumi-pipeline-io`). This path serves readers
//! outside the chain layer. Decoding is the expensive part, so slices decode in
//! parallel on a pool this module owns; a reader thread pulls containers and
//! slices off the stream serially (htslib's own structure, with fgumi's threads
//! in place of htslib's pool). Output order is the file order.
//!
//! The decode pool is process-global and built lazily by the first CRAM input
//! opened here. Its width is the command's `--threads` when `main` passes it
//! through [`set_cram_decode_threads`], and otherwise half the CPUs (at least
//! one). A later call to [`set_cram_decode_threads`] does not resize a pool
//! that already exists. These threads are in addition to the command's own
//! workers.
//!
//! The reference comes from [`set_cram_reference`] (which `main` calls with
//! `fgumi --reference`, or the subcommand's own reference), else the
//! `FGUMI_REFERENCE` environment variable, else htslib's own lookup
//! (`REF_PATH`/`REF_CACHE`, and `@SQ UR:` paths). [`cram_reference_configured`]
//! answers whether any of these is set up, for readers and writers alike.

use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, OnceLock, RwLock};

/// Estimated decoded bytes retained by one standalone input, independent of
/// decode width. One oversized slice may run alone to preserve progress. This
/// bounds admission charges, not htslib's reference cache or exact process RSS.
const DECODE_BYTE_LIMIT: u64 = 256 << 20;

struct DecodeBudget {
    limit: u64,
    state: Mutex<DecodeBudgetState>,
    changed: Condvar,
}

#[derive(Default)]
struct DecodeBudgetState {
    used: u64,
    cancelled: bool,
}

impl DecodeBudget {
    fn new(limit: u64) -> Arc<Self> {
        Arc::new(Self { limit, state: Mutex::default(), changed: Condvar::new() })
    }

    /// Only the serial reader waits here; rayon workers never wait on capacity.
    fn acquire(self: &Arc<Self>, bytes: u64) -> Option<DecodeLease> {
        let mut state = self.state.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        while !state.cancelled && state.used != 0 && bytes > self.limit.saturating_sub(state.used) {
            state = self.changed.wait(state).unwrap_or_else(std::sync::PoisonError::into_inner);
        }
        if state.cancelled {
            return None;
        }
        state.used += bytes;
        Some(DecodeLease { budget: Arc::clone(self), bytes })
    }

    fn cancel(&self) {
        let mut state = self.state.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        state.cancelled = true;
        self.changed.notify_all();
    }
}

struct DecodeLease {
    budget: Arc<DecodeBudget>,
    bytes: u64,
}

impl Drop for DecodeLease {
    fn drop(&mut self) {
        let mut state = self.budget.state.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
        state.used -= self.bytes;
        self.budget.changed.notify_all();
    }
}

/// The lease follows the bytes through worker results and the active chunk.
struct DecodedChunk {
    bytes: Vec<u8>,
    lease: DecodeLease,
}

use anyhow::{Context, Result};
use crossbeam_channel::{Receiver, bounded};
use fgumi_htslib::CramReader;
use noodles::sam;

static REFERENCE: RwLock<Option<PathBuf>> = RwLock::new(None);
static DECODE_THREADS: AtomicUsize = AtomicUsize::new(0);
static POOL: OnceLock<Arc<rayon::ThreadPool>> = OnceLock::new();

/// Set the reference FASTA used to decode CRAM inputs in this process.
pub fn set_cram_reference(path: &Path) {
    *REFERENCE.write().unwrap_or_else(std::sync::PoisonError::into_inner) =
        Some(path.to_path_buf());
}

/// The reference to decode CRAM with, if one is configured.
#[must_use]
pub fn cram_reference() -> Option<PathBuf> {
    REFERENCE
        .read()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .clone()
        .or_else(|| std::env::var_os("FGUMI_REFERENCE").map(PathBuf::from))
}

/// Whether htslib's own reference lookup is set up (`REF_PATH` or `REF_CACHE`
/// in the environment).
fn htslib_reference_lookup_configured() -> bool {
    std::env::var_os("REF_PATH").is_some() || std::env::var_os("REF_CACHE").is_some()
}

/// Whether any CRAM reference source is configured: [`cram_reference`] (the
/// `--reference` flag, a subcommand's own reference, or `FGUMI_REFERENCE`), or
/// htslib's `REF_PATH`/`REF_CACHE` lookup.
///
/// This is the one statement of the lookup chain: CRAM input uses it to decide
/// whether a decode failure is likely a missing reference, and CRAM output to
/// decide whether to write reference-free.
#[must_use]
pub fn cram_reference_configured() -> bool {
    cram_reference().is_some() || htslib_reference_lookup_configured()
}

/// How to give a CRAM reference, for messages shown when none is configured
/// (see [`cram_reference_configured`]).
///
/// `--reference` is a top-level flag, so it goes before the subcommand.
pub const NO_REFERENCE_HINT: &str = "no reference configured: pass `fgumi --reference <FASTA> \
     <command> ...` (before the subcommand), or the command's own --ref/-r where it has one, \
     or set FGUMI_REFERENCE, REF_PATH or REF_CACHE";

/// What a CRAM input error names: the input, and, for a slice-decode failure
/// with no reference configured, how to give one (the usual cause).
///
/// Cheap to clone (the label is shared), so it can move into decode workers.
#[derive(Clone, Debug)]
pub struct CramErrorContext {
    label: Arc<str>,
    has_reference: bool,
}

impl CramErrorContext {
    /// Context for input `label`, decoded with or without a configured
    /// reference.
    #[must_use]
    pub fn new(label: &str, has_reference: bool) -> Self {
        Self { label: Arc::from(label), has_reference }
    }

    /// Context for input `label`, with the reference state taken from
    /// [`cram_reference_configured`].
    #[must_use]
    pub fn from_global(label: &str) -> Self {
        Self::new(label, cram_reference_configured())
    }

    /// Wrap a failure to read the next container or slice off the stream:
    /// truncation, an I/O error, or a malformed container. `e` already says
    /// what went wrong; this only names the input. No reference is involved,
    /// so no reference hint.
    #[must_use]
    pub fn wrap_read(&self, e: &io::Error) -> io::Error {
        io::Error::new(e.kind(), format!("CRAM input {}: {e}", self.label))
    }

    /// Wrap a slice-decode failure, adding [`NO_REFERENCE_HINT`] when no
    /// reference is configured.
    #[must_use]
    pub fn wrap_decode(&self, e: &io::Error) -> io::Error {
        let msg = if self.has_reference {
            format!("CRAM input {}: {e}", self.label)
        } else {
            format!("CRAM input {}: {e} ({NO_REFERENCE_HINT})", self.label)
        };
        io::Error::new(e.kind(), msg)
    }
}

/// Set how many threads decode CRAM slices (default: half the CPUs, at least
/// one). Takes effect for the first CRAM input opened in the process.
pub fn set_cram_decode_threads(n: usize) {
    DECODE_THREADS.store(n.max(1), Ordering::Relaxed);
}

fn pool() -> Result<Arc<rayon::ThreadPool>> {
    if let Some(pool) = POOL.get() {
        return Ok(Arc::clone(pool));
    }
    let n = match DECODE_THREADS.load(Ordering::Relaxed) {
        0 => std::thread::available_parallelism().map_or(2, |n| (n.get() / 2).max(1)),
        n => n,
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .thread_name(|i| format!("fgumi-cram-{i}"))
        .build()
        .context("Failed to start the CRAM decode threads")?;
    // A racing initializer may have won; either pool is equivalent.
    Ok(Arc::clone(POOL.get_or_init(|| Arc::new(pool))))
}

/// Open a CRAM stream (whose first bytes may already be replayed) with the
/// configured reference, and parse its header.
///
/// # Errors
///
/// Returns an error if the stream is not CRAM or its header cannot be read.
pub fn open_cram_reader(
    input: Box<dyn Read + Send>,
    path: &Path,
) -> Result<(CramReader, sam::Header)> {
    let reference = cram_reference();
    let reader = CramReader::from_reader(input, &path.display().to_string(), reference.as_deref())
        .with_context(|| format!("Failed to open CRAM input: {}", path.display()))?;
    let header: sam::Header = reader
        .header_text()
        .parse()
        .with_context(|| format!("Failed to parse the CRAM header: {}", path.display()))?;
    Ok((reader, header))
}

/// Raw (un-BGZF-framed) BAM byte stream decoded from CRAM: the BAM header,
/// then each slice's records in file order.
pub(crate) struct CramToBamBytes {
    /// Pending chunk and read position within it.
    chunk: Vec<u8>,
    pos: usize,
    /// One receiver per slice, in file order.
    slices: Receiver<Receiver<io::Result<DecodedChunk>>>,
    chunk_lease: Option<DecodeLease>,
    decode_budget: Arc<DecodeBudget>,
    /// The serial reader thread (exits on EOF, error, or a dropped consumer);
    /// joined at end of stream so a panic there is not mistaken for EOF.
    /// `None` once joined: the stream has ended.
    reader: Option<std::thread::JoinHandle<()>>,
}

impl CramToBamBytes {
    /// Start decoding `input` (whose first bytes may already be replayed).
    pub(crate) fn open(input: Box<dyn Read + Send>, path: &Path) -> Result<Self> {
        let (mut reader, header) = open_cram_reader(input, path)?;
        let mut first = Vec::new();
        crate::writer::write_bam_header(&mut first, &header)?;

        let context = CramErrorContext::from_global(&path.display().to_string());
        let pool = pool()?;
        // Count bounds receiver overhead; byte admission independently covers
        // queued/running jobs, completed results and the partly consumed chunk.
        let depth = pool.current_num_threads() * 2;
        let decode_budget = DecodeBudget::new(DECODE_BYTE_LIMIT);
        let reader_budget = Arc::clone(&decode_budget);
        let (tx, rx) = bounded::<Receiver<io::Result<DecodedChunk>>>(depth);
        let handle = std::thread::Builder::new()
            .name("fgumi-cram-read".into())
            .spawn(move || {
                loop {
                    let (otx, orx) = bounded(1);
                    match reader.next_slice() {
                        Ok(Some(job)) => {
                            let estimate = match job.decoded_size_estimate() {
                                Ok(estimate) => estimate,
                                Err(error) => {
                                    let _ = otx.send(Err(context.wrap_read(&error)));
                                    let _ = tx.send(orx);
                                    return;
                                }
                            };
                            let Some(lease) = reader_budget.acquire(estimate) else {
                                return;
                            };
                            if tx.send(orx).is_err() {
                                return; // consumer gone
                            }
                            let context = context.clone();
                            pool.spawn(move || {
                                let _ = otx.send(
                                    job.decode()
                                        .map(|d| DecodedChunk { bytes: d.bytes, lease })
                                        .map_err(|e| context.wrap_decode(&e)),
                                );
                            });
                        }
                        Ok(None) => return,
                        Err(e) => {
                            let _ = otx.send(Err(context.wrap_read(&e)));
                            let _ = tx.send(orx);
                            return;
                        }
                    }
                }
            })
            .context("Failed to spawn the CRAM reader thread")?;
        Ok(Self {
            chunk: first,
            pos: 0,
            slices: rx,
            chunk_lease: None,
            decode_budget,
            reader: Some(handle),
        })
    }
}

impl Read for CramToBamBytes {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        while self.pos == self.chunk.len() {
            // Free the exhausted allocation before releasing its admission.
            self.chunk = Vec::new();
            self.pos = 0;
            self.chunk_lease = None;
            if self.reader.is_none() {
                return Ok(0);
            }
            if let Ok(slice) = self.slices.recv() {
                let chunk =
                    slice.recv().map_err(|_| io::Error::other("CRAM decode worker lost"))??;
                self.chunk = chunk.bytes;
                self.chunk_lease = Some(chunk.lease);
            } else {
                // The reader thread hung up: end of file -- unless it died, in
                // which case the tail was never read and that is an error.
                if let Some(handle) = self.reader.take()
                    && handle.join().is_err()
                {
                    return Err(io::Error::other("CRAM reader thread panicked"));
                }
                self.chunk.clear();
            }
            self.pos = 0;
        }
        let n = buf.len().min(self.chunk.len() - self.pos);
        buf[..n].copy_from_slice(&self.chunk[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

impl Drop for CramToBamBytes {
    fn drop(&mut self) {
        // Wake a reader waiting for capacity before dropping receiver results.
        // Do not join: next_slice may be blocked on an external input stream.
        self.decode_budget.cancel();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case::within_allowance(8, 8)]
    #[case::oversized_alone(32, 8)]
    fn decoded_budget_waits_until_prior_lease_releases(#[case] first: u64, #[case] second: u64) {
        let budget = DecodeBudget::new(8);
        let lease = budget.acquire(first).unwrap();
        let (tx, rx) = bounded(1);
        let worker_budget = Arc::clone(&budget);
        let worker = std::thread::spawn(move || tx.send(worker_budget.acquire(second)).unwrap());
        assert!(rx.recv_timeout(std::time::Duration::from_millis(20)).is_err());
        drop(lease);
        let second_lease = rx.recv_timeout(std::time::Duration::from_secs(2)).unwrap().unwrap();
        assert_eq!(budget.state.lock().unwrap().used, second);
        drop(second_lease);
        worker.join().unwrap();
        assert_eq!(budget.state.lock().unwrap().used, 0);
    }

    #[test]
    fn dropping_stream_cancels_blocked_decode_admission() {
        let budget = DecodeBudget::new(1);
        let lease = budget.acquire(2).unwrap();
        let (tx, rx) = bounded(1);
        let worker_budget = Arc::clone(&budget);
        let worker =
            std::thread::spawn(move || tx.send(worker_budget.acquire(2).is_none()).unwrap());
        let (_slices_tx, slices) = bounded(1);
        let stream = CramToBamBytes {
            chunk: vec![0],
            pos: 0,
            slices,
            chunk_lease: Some(lease),
            decode_budget: Arc::clone(&budget),
            reader: None,
        };
        drop(stream);
        assert!(rx.recv_timeout(std::time::Duration::from_secs(2)).unwrap());
        worker.join().unwrap();
        assert_eq!(budget.state.lock().unwrap().used, 0);
        assert!(budget.acquire(1).is_none(), "cancellation is permanent");
    }

    #[rstest]
    #[case::worker_error(false)]
    #[case::worker_disconnected(true)]
    fn decoded_stream_can_retry_after_result_error(#[case] disconnected: bool) {
        let budget = DecodeBudget::new(8);
        let (result_tx, result_rx) = bounded(1);
        if !disconnected {
            result_tx.send(Err(io::Error::other("injected decode error"))).unwrap();
        }
        drop(result_tx);
        let (slices_tx, slices) = bounded(1);
        slices_tx.send(result_rx).unwrap();
        drop(slices_tx);
        let mut stream = CramToBamBytes {
            chunk: vec![1, 2, 3],
            pos: 0,
            slices,
            chunk_lease: None,
            decode_budget: budget,
            reader: Some(std::thread::spawn(|| {})),
        };
        assert_eq!(stream.read(&mut [0; 4]).unwrap(), 3);
        let error = stream.read(&mut [0; 4]).unwrap_err();
        assert!(error.to_string().contains(if disconnected {
            "worker lost"
        } else {
            "injected decode error"
        }));
        assert_eq!(stream.read(&mut [0; 4]).unwrap(), 0);
        assert_eq!(stream.read(&mut [0; 4]).unwrap(), 0);
    }

    #[test]
    fn decoded_stream_repeated_eof_after_nonempty_chunk_is_stable() {
        let budget = DecodeBudget::new(8);
        let (result_tx, result_rx) = bounded(1);
        result_tx
            .send(Ok(DecodedChunk { bytes: vec![1, 2, 3], lease: budget.acquire(8).unwrap() }))
            .unwrap();
        let (slices_tx, slices) = bounded(1);
        slices_tx.send(result_rx).unwrap();
        drop(slices_tx);
        let mut stream = CramToBamBytes {
            chunk: Vec::new(),
            pos: 0,
            slices,
            chunk_lease: None,
            decode_budget: budget,
            reader: Some(std::thread::spawn(|| {})),
        };
        let mut actual = Vec::new();
        stream.read_to_end(&mut actual).unwrap();
        assert_eq!(actual, [1, 2, 3]);
        for _ in 0..3 {
            assert_eq!(stream.read(&mut [0; 4]).unwrap(), 0);
        }
    }

    #[test]
    fn decoded_chunk_holds_lease_across_result_queue_and_partial_read() {
        let budget = DecodeBudget::new(8);
        let lease = budget.acquire(8).unwrap();
        let (result_tx, result_rx) = bounded(1);
        result_tx.send(Ok(DecodedChunk { bytes: vec![1, 2, 3], lease })).unwrap();
        let (slices_tx, slices) = bounded(1);
        slices_tx.send(result_rx).unwrap();
        drop(slices_tx);
        let mut stream = CramToBamBytes {
            chunk: Vec::new(),
            pos: 0,
            slices,
            chunk_lease: None,
            decode_budget: Arc::clone(&budget),
            reader: Some(std::thread::spawn(|| {})),
        };
        let mut byte = [0];
        assert_eq!(stream.read(&mut byte).unwrap(), 1);
        assert_eq!(byte, [1]);
        assert_eq!(budget.state.lock().unwrap().used, 8);
        drop(stream);
        assert_eq!(budget.state.lock().unwrap().used, 0);
    }

    /// The reference hint goes only on a decode failure with no reference
    /// configured; a read failure (truncation, I/O) never carries it.
    #[rstest]
    #[case::decode_without_reference(true, false, true)]
    #[case::decode_with_reference(true, true, false)]
    #[case::read_without_reference(false, false, false)]
    #[case::read_with_reference(false, true, false)]
    fn reference_hint_only_on_decode_failures_without_a_reference(
        #[case] decode: bool,
        #[case] has_reference: bool,
        #[case] expect_hint: bool,
    ) {
        let context = CramErrorContext::new("in.cram", has_reference);
        let cause = io::Error::other("htslib failed");
        let err = if decode { context.wrap_decode(&cause) } else { context.wrap_read(&cause) };
        let msg = err.to_string();
        assert!(msg.starts_with("CRAM input in.cram: htslib failed"), "unexpected message: {msg}");
        assert_eq!(msg.contains(NO_REFERENCE_HINT), expect_hint, "hint placement: {msg}");
    }

    /// The hint names the top-level flag before the subcommand, and every
    /// fallback the lookup chain honors.
    #[rstest]
    #[case::flag_before_subcommand("fgumi --reference <FASTA> <command>")]
    #[case::subcommand_ref("--ref/-r")]
    #[case::fgumi_reference_env("FGUMI_REFERENCE")]
    #[case::ref_path_env("REF_PATH")]
    #[case::ref_cache_env("REF_CACHE")]
    fn hint_names_every_reference_source(#[case] needle: &str) {
        assert!(NO_REFERENCE_HINT.contains(needle), "hint lacks {needle}: {NO_REFERENCE_HINT}");
    }
}
