# Read Structures

## Overview

A *Read Structure* is a string that describes how the bases in a sequencing run should be allocated into logical reads. It serves a similar purpose to the `--use-bases-mask` in Illumina's bcl-convert, but provides additional capabilities.

A Read Structure is a sequence of `<number><operator>` pairs (called *segments*). The last segment may use `+` instead of a number to mean "whatever bases remain." fgumi uses the [`read-structure`](https://crates.io/crates/read-structure) crate for parsing and validation.

Read structures are used primarily in `fgumi extract` to specify where UMI bases, template bases, and other sequences are located in each FASTQ read.

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

## Formal Grammar

```
<read-structure>     ::= <fixed-structure> <variable-segment>
<fixed-structure>    ::= "" | <fixed-length> <operator> <fixed-structure>
<variable-segment>   ::= "" | <variable-length> <operator>
<segment>            ::= <any-length><operator>
<operator>           ::= "T" | "B" | "M" | "C" | "S"
<fixed-length>       ::= <non-zero-digit>{<digit>}
<variable-length>    ::= "+"
<any-length>         ::= <fixed-length> | <variable-length>
<non-zero-digit>     ::= "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9"
<digit>              ::= "0" | <non-zero-digit>
```
