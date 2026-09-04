//! FASTQ-writing helpers shared by the `runall` integration tests.
//!
//! `runall`'s extract-fusion tests (`test_runall_command.rs`,
//! `test_runall_chain_transitions.rs`) both need to stage small gzip-compressed
//! FASTQ fixtures; this used to be duplicated verbatim in both files.

use std::fs;
use std::io::Write as _;
use std::path::Path;

/// Write a gzip-compressed FASTQ from `(name, seq, qual)` records.
pub fn write_gzip_fastq(path: &Path, records: &[(&str, &str, &str)]) {
    let file = fs::File::create(path).expect("create gzip fastq");
    let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    for (name, seq, qual) in records {
        writeln!(encoder, "@{name}\n{seq}\n+\n{qual}").expect("write fastq record");
    }
    encoder.finish().expect("finish gzip fastq");
}
