# fgumi-sort

High-performance external merge-sort engine for BAM files.

This crate powers `fgumi sort` and is published so that downstream consumers can build
their own BAM sort tools without taking on the rest of fgumi.

## Sort orders

- `coordinate` — standard genomic coordinate sort (tid → pos → strand)
- `queryname` (lexicographic or natural) — read-name ordering
- `template-coordinate` — paired-end template position, used for UMI grouping

## Performance

External merge-sort built from:

- parallel radix sort over packed sort keys for in-memory chunks
- K-way merge via a loser tree to combine spill files
- async I/O prefetch and spill compression (BGZF legacy / Zstd default)
- free-space-aware, round-robin spill across multiple temp directories

Designed to outperform `samtools sort` on the large-BAM sort workloads it
targets; actual speedup depends on the data, sort order, thread count, and
hardware.

## Usage

### Sort a BAM file

```rust,no_run
use fgumi_sort::{RawExternalSorter, SortOrder};
use fgumi_sort::keys::QuerynameComparator;
use std::path::Path;

// Coordinate sort (for IGV, variant callers, etc.)
RawExternalSorter::new(SortOrder::Coordinate)
    .threads(4)
    .output_compression(6)
    .sort(Path::new("input.bam"), Path::new("sorted.bam"))
    .expect("sort failed");

// Queryname sort — natural ordering (read1 < read2 < read10)
RawExternalSorter::new(SortOrder::Queryname(QuerynameComparator::Natural))
    .threads(4)
    .sort(Path::new("input.bam"), Path::new("by-name.bam"))
    .expect("sort failed");

// Template-coordinate sort (required before fgumi group)
RawExternalSorter::new(SortOrder::TemplateCoordinate)
    .threads(4)
    .sort(Path::new("input.bam"), Path::new("template-coord.bam"))
    .expect("sort failed");
```

### Verify an existing sort order

```rust,no_run
use fgumi_sort::{RawBamRecordReader, SortOrder, extract_coordinate_key_inline, verify_sort_order};
use fgumi_sort::keys::QuerynameComparator;
use std::fs::File;

// Check that a file is coordinate-sorted
let file = File::open("candidate.bam").unwrap();
let mut reader = RawBamRecordReader::new(file).unwrap();
reader.skip_header().unwrap();

let nref: u32 = 25; // number of reference sequences in the header
let (total, violations, first_violation) = verify_sort_order(
    reader,
    |bam| extract_coordinate_key_inline(bam, nref),
    |key, prev| key < prev,
)
.unwrap();

if violations == 0 {
    println!("File is coordinate-sorted ({total} records)");
} else {
    let (rec_num, name) = first_violation.unwrap();
    println!("{violations} violations; first at record {rec_num} ({name})");
}
```

### Memory and temp-dir tuning

```rust,no_run
use fgumi_sort::{RawExternalSorter, SortOrder};
use std::path::{Path, PathBuf};

RawExternalSorter::new(SortOrder::Coordinate)
    .threads(8)
    // Limit in-memory sort buffer to ~4 GiB (default: 80% of available RAM)
    .max_memory_bytes(4 * 1024 * 1024 * 1024)
    // Write spill files to a fast scratch volume
    .tmp_dirs(vec![PathBuf::from("/scratch/sort-spill")])
    .output_compression(6)
    .sort(Path::new("large.bam"), Path::new("large_sorted.bam"))
    .expect("sort failed");
```

## Status

Pre-1.0. APIs may evolve. Pin a specific minor version.

## License

MIT — see the workspace `LICENSE` file.
