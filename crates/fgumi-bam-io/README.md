# fgumi-bam-io

BAM I/O factories and pipeline primitives for the [fgumi](https://github.com/fulcrumgenomics/fgumi) family of crates.

This crate provides:

- BAM and raw-BAM reader/writer factories built on top of `noodles` and libdeflate-backed BGZF
- BAI index writing
- The `ProgressTracker` used by fgumi's pipelines
- The `ReorderBuffer` used to maintain output order through parallel workers
- The `MemoryEstimate` trait used by fgumi's pipeline backpressure
- A vendored, position-tracking variant of noodles' multithreaded BGZF writer

## Status

Pre-1.0. APIs may evolve. Pin a specific minor version.

## License

MIT — see the workspace `LICENSE` file.
