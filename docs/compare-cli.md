# fgumi compare CLI Reference

Developer tools for comparing BAM files and metrics outputs for testing and validation.

**Requires:** `cargo build --features compare`

## Commands Overview

| Command | Description | Primary Use Case |
|---------|-------------|------------------|
| `fgumi compare bams` | Compare two BAM files | Validate BAM output equality |
| `fgumi compare metrics` | Compare two TSV metrics files | Compare fgumi/fgbio metrics |

---

## fgumi compare bams

Compare two BAM files for equality of core SAM fields and tag values.

### Usage

```bash
fgumi compare bams <BAM1> <BAM2> [OPTIONS]
```

### Required Arguments

| Option | Type | Description |
|--------|------|-------------|
| `<BAM1>` | PATH | First BAM file |
| `<BAM2>` | PATH | Second BAM file |

### Comparison Modes

| Mode | Description | Use Case |
|------|-------------|----------|
| `full` (default) | MI grouping check + full content comparison | Same-tool, same-order `group` output |
| `content` | Pure record-by-record comparison of all fields | `extract`, `zipper`, `correct`, `dedup`, `filter`, `clip`, `simplex`, `duplex`, `codec` output (the latter three via the saturation-aware `ExactConsensus` predicate) |
| `grouping` | MI equivalence only, or (with `--ignore-order`) key-join content + MI bijection | Cross-tool `group` output, where MI values/order may differ (see `--command group`) |

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-c, --command` | STRING | (none) | Preset `--mode`/`--ignore-order` for a specific fgumi pipeline stage (see [Command Presets](#command-presets---command) below). An explicit `--mode`/`--ignore-order` overrides the preset. |
| `--mode` | STRING | `full` | Comparison mode: `full`, `content`, or `grouping` |
| `-m, --max-diffs` | INT | 10 | Maximum differences to report |
| `-q, --quiet` | FLAG | false | Only exit code indicates result (0=equal, 1=different) |
| `--ignore-order` | FLAG | false | Ignore record order (for parallel consensus output) |
| `--buffer-size` | INT | 1000 | Initial buffer size for `--ignore-order` mode |
| `-t, --threads` | INT | 1 | Threads for BGZF decompression and comparison |
| `--batch-size` | INT | 10000 | Records per batch for parallel processing |
| `--sort-memory` | STRING | `512M` | Total (not per-thread; not multiplied by `--threads`) memory budget for the internal queryname-canonicalization sort. Consumed only by `--command group`'s key-join engine; ignored by every other mode/preset. |
| `--sort-tmp-dir` | PATH | (platform default) | Temporary directory for the `--command group` key-join canonicalization sort's spill files. Repeatable (may be passed more than once); same semantics as `fgumi sort -T`/`--tmp-dir`. |

### Command Presets (`--command`)

`--command <stage>` (short form `-c`) applies the canonical `--mode`/`--ignore-order`
defaults for comparing BAM output from a specific fgumi pipeline stage, so you don't
need to remember which mode/flags a given command's output requires. This is the
primary, recommended way to compare a pipeline stage's output — especially for
cross-tool comparisons (e.g. fgumi vs. fgbio) where MI values or record order may
legitimately differ. An explicit `--mode` or `--ignore-order` overrides the preset.

Presets: `extract`, `zipper`, `sort`, `correct`, `dedup`, `group`, `simplex`, `duplex`,
`codec`, `filter`. See [Recommended Settings by Command](#recommended-settings-by-command)
below for what each resolves to.

```bash
fgumi compare bams a.bam b.bam --command extract
fgumi compare bams a.bam b.bam --command group
fgumi compare bams a.bam b.bam --command simplex
```

### Comparison Fields

The tool compares:
- **Core SAM fields**: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL
- **Tag values**: Order-independent comparison (tags may appear in different order)

### Recommended Settings by Command

These are the same defaults each `--command <stage>` preset applies automatically
(see [Command Presets](#command-presets---command) above); use `--command` rather
than setting `--mode`/`--ignore-order` by hand where a preset exists.

| fgumi Command | --mode | --ignore-order | Notes |
|---------------|--------|----------------|-------|
| `extract` | content | false | No MI tags; deterministic |
| `zipper` | content | false | Preserves MI tags unchanged |
| `group` | grouping (key-join engine) | true | Cross-tool: compares via key-join, pairing records by identity after canonicalizing both inputs to queryname order. Content is checked under `ExactMinusMi` (every field except MI) plus a separate fgumi-MI/fgbio-MI bijection check, since MI values or record order may legitimately differ between tools. |
| `simplex` | content | false | Saturation-aware exact comparison (`ExactConsensus`): tolerates the accepted `cD`/`cM`/`cE` depth-saturation divergence, otherwise exact |
| `duplex` | content | false | Saturation-aware exact comparison (`ExactConsensus`): tolerates the accepted depth-saturation divergence (`cD`/`cM`/`cE` plus duplex `aD`/`aM`/`bD`/`bM`/`aE`/`bE`), otherwise exact |
| `codec` | content | false | Saturation-aware exact comparison (`ExactConsensus`), same carve-out as `simplex`/`duplex` |
| `filter` | content | false | Passes through MI/depth tags unchanged; uses the same saturation-aware `ExactConsensus` predicate as consensus output |
| `clip` | content | false | Does not modify MI tags |
| `correct` | content | false | Modifies RX tag only, not MI |
| `downsample` | content | false | Deterministic with seed |
| `review` | content | false | Preserves MI tags |

### Examples

```bash
# Compare extract output (no MI tags)
fgumi compare bams extracted1.bam extracted2.bam --mode content

# Compare group output (cross-tool: key-join + MI bijection check)
fgumi compare bams grouped1.bam grouped2.bam --command group

# Compare consensus output (saturation-aware exact content comparison)
fgumi compare bams consensus1.bam consensus2.bam --command simplex
fgumi compare bams consensus1.bam consensus2.bam --command duplex
fgumi compare bams consensus1.bam consensus2.bam --command codec

# Compare with multi-threading for large files
fgumi compare bams large1.bam large2.bam --mode content --threads 8

# Quiet mode for CI pipelines
fgumi compare bams bam1.bam bam2.bam --quiet && echo "Files match"
```

### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Files are equal |
| 1 | Files are different |

---

## fgumi compare metrics

Compare two TSV metrics files by outer-joining their data rows on a key column
set, so rows may legitimately appear in different order (or one side may be
missing a row entirely -- e.g. fgumi#498's dense-vs-sparse family-size rows)
without the comparison cascading into spurious diffs on every following line.
Only float *representation* is tolerated (fgbio's `DecimalFormat("0.######")`
vs fgumi's full precision); everything else -- including the row set and any
non-key value -- is compared faithfully.

### Usage

```bash
fgumi compare metrics <FILE1> <FILE2> [OPTIONS]
```

### Required Arguments

| Option | Type | Description |
|--------|------|-------------|
| `<FILE1>` | PATH | First TSV file |
| `<FILE2>` | PATH | Second TSV file |

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--precision` | INT | 6 | Round floats to this many decimal places (-1 to disable) |
| `--rel-tol` | FLOAT | 1e-9 | Relative tolerance for float comparison |
| `--abs-tol` | FLOAT | 1e-9 | Absolute tolerance for float comparison |
| `--key-columns` | STRING | (first column) | Comma-separated column name(s) or 0-based index(es) to join rows on, e.g. `ab_size,ba_size` |
| `-m, --max-diffs` | INT | 20 | Maximum differences to display |
| `-q, --quiet` | FLAG | false | Only exit code indicates result |
| `-v, --verbose` | FLAG | false | Print success message when files match |

### Comparison Behavior

- Both files' header rows must declare the **same column set** (order does not
  matter). A column-set difference is itself reported as `DIFFER` -- row values
  cannot be meaningfully aligned without it.
- Data rows are joined on the key column(s) (default: the first column) as an
  **outer join**: a key present in both files has its remaining columns
  compared; a key present in only one file is reported individually (e.g.
  `row 3 only in file1`) instead of cascading into every following row.
- Within a matched key: **integers** are compared exactly, **floats** are
  optionally rounded then compared with tolerance, **strings** are compared
  exactly.
- A key that appears more than once in a file (a duplicate key) is matched
  against the other file's same-key rows as a multiset, not by on-disk order.

### Float Comparison

Floats are considered equal if either condition is met:
- `|a - b| <= abs_tol`
- `|a - b| <= rel_tol * max(|a|, |b|)`

### Examples

```bash
# Basic comparison with default precision
fgumi compare metrics file1.txt file2.txt

# Compare with 6 decimal precision
fgumi compare metrics file1.txt file2.txt --precision 6

# Looser tolerance for cross-platform comparison
fgumi compare metrics file1.txt file2.txt --rel-tol 1e-6 --abs-tol 1e-9

# Multi-column key, e.g. duplex_family_sizes keyed on (ab_size, ba_size)
fgumi compare metrics file1.txt file2.txt --key-columns ab_size,ba_size

# Quiet mode for CI
fgumi compare metrics fgumi_metrics.txt fgbio_metrics.txt --quiet
```

### Use Cases

1. **Comparing fgumi vs fgbio output**: Different tools may produce slightly different floating-point representations, and rows may legitimately appear in a different order or set (tracked fgbio-parity gaps, e.g. fgumi#498)
2. **Cross-platform validation**: Float formatting may differ between platforms
3. **Regression testing**: Verify metrics haven't changed between versions

### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Files are equal (within tolerance) |
| 1 | Files are different |

---

## Workflow Examples

### Validating a Full Pipeline

```bash
# Run pipeline twice with different implementations
fgumi extract ... --output extract1.bam
fgbio ExtractUmisFromBam ... --output extract2.bam

# Compare outputs
fgumi compare bams extract1.bam extract2.bam --mode content

# Compare downstream metrics
fgumi duplex-metrics --input duplex1.bam --output metrics1
fgumi duplex-metrics --input duplex2.bam --output metrics2
fgumi compare metrics metrics1.family_sizes.txt metrics2.family_sizes.txt
```

### CI/CD Integration

```bash
#!/bin/bash
set -e

# Run tests and compare outputs
fgumi group --input test.bam --output actual.bam --strategy adjacency
fgumi compare bams expected.bam actual.bam --mode full --quiet

echo "BAM comparison passed"
```
