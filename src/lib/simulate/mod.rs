//! Simulation utilities for generating synthetic sequencing data.
//!
//! This module provides utilities for generating synthetic FASTQ and BAM files
//! for benchmarking and testing the fgumi pipeline.
//!
//! # Modules
//!
//! - [`rng`] - Seeded random number generator utilities
//! - [`quality`] - Position-dependent quality score models
//! - [`insert_size`] - Insert size distribution models
//! - [`family_size`] - Family size distribution models
//! - [`strand_bias`] - A/B strand ratio models for duplex
//! - [`fastq_writer`] - Gzipped FASTQ file writer

pub mod family_size;
pub mod fastq_writer;
pub mod insert_size;
pub mod parallel_gzip_writer;
pub mod quality;
pub mod rng;
pub mod strand_bias;

pub use family_size::FamilySizeDistribution;
pub use fastq_writer::FastqWriter;
pub use insert_size::InsertSizeModel;
pub use quality::{PositionQualityModel, ReadPairQualityBias};
pub use rng::create_rng;
pub use strand_bias::StrandBiasModel;

/// Flush and close a buffered simulation output (see
/// [`fgumi_bam_io::OutputSink::close`]), naming `path` in any error.
///
/// # Errors
///
/// Returns an error if flushing or closing the output fails.
pub(crate) fn close_output<W: fgumi_bam_io::OutputSink>(
    writer: std::io::BufWriter<W>,
    path: &std::path::Path,
) -> anyhow::Result<()> {
    use anyhow::Context;

    fgumi_bam_io::close_buffered(writer)
        .with_context(|| format!("Failed to close {}", path.display()))
}

/// Test doubles shared by the simulation tests.
#[cfg(test)]
pub(crate) mod test_support {
    use std::io::{self, Read, Write};
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::{Arc, Mutex};

    use fgumi_bam_io::OutputSink;

    /// The error a failing [`CloseProbe`] close returns.
    pub(crate) const CLOSE_ERROR: &str = "input/output error at close";

    /// An in-memory [`OutputSink`] that records its bytes, counts closes, and
    /// can fail its close, modelling an error reported at close (such as an
    /// NFS flush failing). Clones share state.
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
            self.written.lock().expect("probe lock").clone()
        }

        pub(crate) fn closes(&self) -> usize {
            self.closes.load(Ordering::SeqCst)
        }

        /// The bytes written, decompressed as a (possibly concatenated) gzip
        /// stream.
        pub(crate) fn decoded(&self) -> Vec<u8> {
            let mut decoded = Vec::new();
            flate2::read::MultiGzDecoder::new(self.bytes().as_slice())
                .read_to_end(&mut decoded)
                .expect("the probe must hold a decodable gzip stream");
            decoded
        }
    }

    impl Write for CloseProbe {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.written.lock().expect("probe lock").extend_from_slice(buf);
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

#[cfg(test)]
mod tests {
    use super::close_output;
    use super::test_support::{CLOSE_ERROR, CloseProbe};
    use std::io::{BufWriter, Write};
    use std::path::Path;

    #[test]
    fn close_output_flushes_and_closes_once() {
        let probe = CloseProbe::default();
        let mut writer = BufWriter::new(probe.clone());
        writer.write_all(b"truth\n").unwrap();
        close_output(writer, Path::new("truth.tsv")).unwrap();
        assert_eq!(probe.bytes(), b"truth\n");
        assert_eq!(probe.closes(), 1);
    }

    /// A close failure surfaces naming the output, since a simulate command
    /// writes several files.
    #[test]
    fn close_output_names_the_path_in_close_errors() {
        let writer = BufWriter::new(CloseProbe::failing_close());
        let err = close_output(writer, Path::new("out/truth.tsv")).expect_err("must fail");
        let msg = format!("{err:#}");
        assert!(msg.contains("out/truth.tsv") && msg.contains(CLOSE_ERROR), "{msg}");
    }
}
