# fgumi-sort

Engine components for sorting and merging BAM records.

This crate holds the pieces `fgumi sort`, `fgumi merge` and `fgumi simulate` are
built from: sort keys and their inline extractors, the in-memory chunk sorters,
the spill-file block format, the k-way merge driver and its loser tree, the merge
slots the spill read planner feeds, and a sort-order verifier. The sort itself is
a chain-native pipeline that lives in the `fgumi` crate; this crate no longer
ships a standalone sorter.

## Sort orders

- `coordinate` — standard genomic coordinate sort (tid → pos → strand)
- `queryname` (lexicographic or natural) — read-name ordering
- `template-coordinate` — paired-end template position, used for UMI grouping

## Usage

### Sort a BAM file

Use the `fgumi sort` command, or, from Rust, `fgumi_lib::sort::sort_bam` (in the
`fgumi` crate), which builds the same chain the command runs. The full option set
is `fgumi_lib::commands::sort::Sort`.

```rust,ignore
use fgumi_lib::commands::sort::SortOrderArg;
use fgumi_lib::sort::{KeyTypesSpec, SortConfig, sort_bam};
use std::path::PathBuf;

let stats = sort_bam(&SortConfig {
    input: PathBuf::from("input.bam"),
    output: PathBuf::from("sorted.bam"),
    order: SortOrderArg::Coordinate,
    threads: 8,
    // A 4 GiB in-memory budget for the whole run (None: 768 MiB per thread).
    max_memory: Some(4 * 1024 * 1024 * 1024),
    // Spill to a fast scratch volume (empty: the system temp directory).
    tmp_dirs: vec![PathBuf::from("/scratch/sort-spill")],
    compression: 6,
    write_index: true,
    key_types: KeyTypesSpec::default(),
})
.expect("sort failed");
println!("sorted {} records", stats.total_records);
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

## Status

Pre-1.0. APIs may evolve. Pin a specific minor version.

## License

MIT — see the workspace `LICENSE` file.
