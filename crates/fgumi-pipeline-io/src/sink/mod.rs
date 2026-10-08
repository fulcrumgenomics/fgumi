pub mod write_bgzf;
pub mod write_raw;

/// Test doubles shared by the sink tests.
#[cfg(test)]
pub(crate) mod test_support {
    use std::io::{self, Write};
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Arc, Mutex};

    use fgumi_bam_io::OutputSink;

    /// The error a failing [`CloseProbe`] close returns.
    pub(crate) const CLOSE_ERROR: &str = "input/output error at close";

    /// An in-memory [`OutputSink`] that records its bytes, counts closes, and
    /// can fail its close, modelling write-back (or an NFS flush) failing only
    /// when the file is closed. Clones share state.
    #[derive(Clone, Default)]
    pub(crate) struct CloseProbe {
        written: Arc<Mutex<Vec<u8>>>,
        closes: Arc<AtomicUsize>,
        fail_close: bool,
    }

    impl CloseProbe {
        /// A probe whose close fails with [`CLOSE_ERROR`].
        pub(crate) fn failing_close() -> Self {
            Self { fail_close: true, ..Self::default() }
        }

        pub(crate) fn bytes(&self) -> Vec<u8> {
            self.written.lock().unwrap().clone()
        }

        pub(crate) fn closes(&self) -> usize {
            self.closes.load(Ordering::SeqCst)
        }
    }

    impl Write for CloseProbe {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.written.lock().unwrap().extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    impl OutputSink for CloseProbe {
        fn close(self: Box<Self>) -> io::Result<()> {
            self.closes.fetch_add(1, Ordering::SeqCst);
            if self.fail_close {
                return Err(io::Error::other(CLOSE_ERROR));
            }
            Ok(())
        }
    }
}
