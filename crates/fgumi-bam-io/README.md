# fgumi-bam-io

BAM I/O factories and pipeline primitives for the [fgumi](https://github.com/fulcrumgenomics/fgumi) family of crates.

This crate provides:

- BAM and raw-BAM reader/writer factories built on top of `noodles` and libdeflate-backed BGZF
- BAI sidecar writing (`.bam.bai` convention)
- `ProgressTracker` — thread-safe progress reporting for fgumi's pipelines
- `ReorderBuffer` — maintains sequential output order across parallel workers
- the `MemoryEstimate` trait used for fgumi's pipeline backpressure
- a vendored, position-tracking variant of noodles' multithreaded BGZF writer

## Usage

Open a BAM reader and create a writer with the same header (single-threaded; pass a higher
thread count for multithreaded BGZF compression):

```rust,no_run
use fgumi_bam_io::reader::create_bam_reader;
use fgumi_bam_io::writer::create_bam_writer;
use std::path::Path;

// Open a reader (1 decompression thread) and a writer that reuses the header
// (1 compression thread, fast compression level 1).
let (mut reader, header) = create_bam_reader(Path::new("input.bam"), 1).unwrap();
let mut writer = create_bam_writer(Path::new("output.bam"), &header, 1, 1).unwrap();
// ... read records from `reader` and write them to `writer` ...
```

## Status

Pre-1.0. APIs may evolve. Pin a specific minor version.

## License

MIT — see the workspace `LICENSE` file.
