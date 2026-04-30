# fgumi-sort

High-performance external merge-sort engine for BAM files.

This crate powers `fgumi sort` and is published so that downstream consumers can build
their own BAM sort tools without taking on the rest of fgumi.

## Sort orders

- `coordinate` — standard genomic coordinate sort (tid → pos → strand)
- `queryname` (lexicographic or natural) — read-name ordering
- `template-coordinate` — paired-end template position, used for UMI grouping

## Performance

External merge-sort with parallel radix sort over packed sort keys, K-way merge via a
loser tree, async I/O prefetch, libdeflate-backed BGZF, and free-space-aware temp-dir
round-robin. Faster than `samtools sort` on the workloads it targets.

## Status

Pre-1.0. APIs may evolve. Pin a specific minor version.

## License

MIT — see the workspace `LICENSE` file.
