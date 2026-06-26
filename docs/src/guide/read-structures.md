# Read Structures

## Overview

A *Read Structure* is a string that describes how the bases in a sequencing read should be allocated into logical segments. It serves a similar purpose to `--use-bases-mask` in Illumina's bcl-convert, but provides additional capabilities. (If you are familiar with bcl-convert, read structures are the more flexible equivalent.)

A read structure is a sequence of `<length><operator>` segments, where `<length>` is a positive integer or, in the final segment, `+`. For example, `8M+T` means **8** bases of UMI (`M`), then **all remaining** bases (`+`) of template (`T`). Only the last segment may use `+`, to mean "whatever bases remain."

You supply one read structure per FASTQ, primarily to `fgumi extract`, to tell it where the UMI, template, sample-barcode, cell-barcode, and skip bases sit in each read.

## Operators

Five kinds of operator are supported:

| Operator | Name | Meaning |
|----------|------|---------|
| `T` | Template | Reads of template (e.g. genomic DNA, RNA) |
| `B` | Sample Barcode | Index sequence for sample identification |
| `M` | Molecular Barcode | UMI sequence for identifying the source molecule |
| `C` | Cell Barcode | Index sequence for identifying the cell (single-cell) |
| `S` | Skip | Bases to skip or ignore (e.g. monotemplate from library prep) |

## Rules

- Any number of segments >= 1 is valid
- The length of each segment must be a positive integer >= 1, or `+`
- Only the last segment in a read structure may use `+` for its length
- Adjacent segments may use the same operator (e.g. `6B6B+T` is valid if two sample indices are ligated separately)

## Examples

### Simple paired-end (2x150bp, no indices)

Per-read structures: `+T`, `+T`

### Paired-end with 8bp sample index

Per-read structures: `+T`, `8B`, `+T`

### Paired-end with inline 6bp UMI in R1

Per-read structures: `6M+T`, `8B`, `+T`

The first 6 bases of R1 are the UMI, followed by template.

### Duplex sequencing with dual barcoding and UMI + monotemplate

Per-read structures: `10M5S+T`, `8B`, `8B`, `10M5S+T`

Both R1 and R2 start with a 10bp UMI followed by 5bp of monotemplate (skipped), then template.

### Single-cell with cell barcodes and UMI

Per-read structures: `5C30S5C3S8M+T`, `8B`, `+T`

R1 contains two cell barcodes separated by linker sequences, then a UMI, then template.

## Troubleshooting

`fgumi extract` rejects a read structure whose fixed-length segments add up to more bases than a read contains, or whose number of read structures does not match the number of input FASTQs. If you hit an error:

- Count the fixed-length bases in each structure and confirm they fit within the read length; use `+` for the final segment so it absorbs whatever remains.
- Supply exactly one read structure per FASTQ, in the same order as `--inputs`.
- If all UMI families come out as singletons afterward, the `M` segment is likely the wrong length or in the wrong read — see [Troubleshooting](troubleshooting.md).

## Technical Details

Most users never need this section — the operators, rules, and examples above cover real read structures. It is included for completeness. fgumi uses the [`read-structure`](https://crates.io/crates/read-structure) crate for parsing and validation; the full grammar is:

```text
<read-structure>     ::= <fixed-structure> <segment>
<fixed-structure>    ::= "" | <fixed-length> <operator> <fixed-structure>
<segment>            ::= <fixed-length> <operator> | <variable-length> <operator>
<operator>           ::= "T" | "B" | "M" | "C" | "S"
<fixed-length>       ::= <non-zero-digit>{<digit>}
<variable-length>    ::= "+"
<any-length>         ::= <fixed-length> | <variable-length>
<non-zero-digit>     ::= "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9"
<digit>              ::= "0" | <non-zero-digit>
```

## See Also

- [Getting Started](getting-started.md) — read structures in the context of `fgumi extract`
- [UMI Grouping](umi-grouping.md) — what happens to the extracted UMIs next
