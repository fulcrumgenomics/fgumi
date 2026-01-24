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
| `full` (default) | MI grouping check + full content comparison | `group` output |
| `content` | Pure record-by-record comparison of all fields | `extract`, `filter`, `clip`, `correct` output |
| `grouping` | MI equivalence only (ignores other content) | `simplex`, `duplex`, `codec` output |

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--mode` | STRING | `full` | Comparison mode: `full`, `content`, or `grouping` |
| `-m, --max-diffs` | INT | 10 | Maximum differences to report |
| `-q, --quiet` | FLAG | false | Only exit code indicates result (0=equal, 1=different) |
| `--ignore-order` | FLAG | false | Ignore record order (for parallel consensus output) |
| `--buffer-size` | INT | 1000 | Initial buffer size for `--ignore-order` mode |
| `-t, --threads` | INT | 1 | Threads for BGZF decompression and comparison |
| `--batch-size` | INT | 10000 | Records per batch for parallel processing |

### Comparison Fields

The tool compares:
- **Core SAM fields**: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL
- **Tag values**: Order-independent comparison (tags may appear in different order)

### Recommended Settings by Command

| fgumi Command | --mode | --ignore-order | Notes |
|---------------|--------|----------------|-------|
| `extract` | content | false | No MI tags; deterministic |
| `zipper` | content | false | Preserves MI tags unchanged |
| `group` | full | false | Has MI tags; deterministic |
| `simplex` | grouping | true | Non-deterministic with `--threads` |
| `duplex` | grouping | true | Non-deterministic with `--threads` |
| `codec` | grouping | true | Non-deterministic with `--threads` |
| `filter` | content | false | Passes through MI tags unchanged |
| `clip` | content | false | Does not modify MI tags |
| `correct` | content | false | Modifies RX tag only, not MI |
| `downsample` | content | false | Deterministic with seed |
| `review` | content | false | Preserves MI tags |

### Examples

```bash
# Compare extract output (no MI tags)
fgumi compare bams extracted1.bam extracted2.bam --mode content

# Compare group output (has MI tags)
fgumi compare bams grouped1.bam grouped2.bam --mode full

# Compare consensus output (MI grouping, non-deterministic order)
fgumi compare bams consensus1.bam consensus2.bam --mode grouping --ignore-order

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

Compare two TSV metrics files with optional float precision rounding.

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
| `-m, --max-diffs` | INT | 20 | Maximum differences to display |
| `-q, --quiet` | FLAG | false | Only exit code indicates result |
| `-v, --verbose` | FLAG | false | Print success message when files match |

### Comparison Behavior

- **Integers**: Compared exactly
- **Floats**: Optionally rounded, then compared with tolerance
- **Strings**: Compared exactly

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

# Quiet mode for CI
fgumi compare metrics fgumi_metrics.txt fgbio_metrics.txt --quiet
```

### Use Cases

1. **Comparing fgumi vs fgbio output**: Different tools may produce slightly different floating-point representations
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
