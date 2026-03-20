# Working with Metrics

fgumi commands produce structured metrics files for quality control and analysis. This guide covers the file formats, terminology, and how to work with the outputs.

## Commands that Produce Metrics

| Command | Metrics Output | Flag |
|---------|---------------|------|
| `filter` | Filtering pass/fail statistics | `--stats` |
| `simplex` | Consensus calling statistics | `--stats` |
| `duplex` | Consensus calling statistics | `--stats` |
| `codec` | Consensus calling statistics | `--stats` |
| `dedup` | Deduplication metrics and family size histogram | `--metrics`, `--family-size-histogram` |
| `duplex-metrics` | Comprehensive duplex QC metrics | `--output` (prefix) |
| `group` | Family size histogram and grouping metrics | `--family-size-histogram`, `--grouping-metrics` |

See the [Metrics Reference](../metrics/README.md) for field-level documentation of each metric type.

## File Formats

Most metrics files are tab-separated values (TSV) with a header row. There are two formats:

### Horizontal TSV (Most Commands)

A header row followed by a single data row. Used by `dedup`, `codec`, `duplex-metrics`, and `group`.

```text
total_templates	unique_templates	duplicate_templates	duplicate_rate
25000	18750	6250	0.25
```

### Vertical Key-Value (Simplex/Duplex)

The `simplex` and `duplex` commands use a three-column format with one metric per row:

```text
key	value	description
raw_reads_considered	50000	Total raw reads considered from input file
raw_reads_used	41800	Total count of raw reads used in consensus reads
consensus_reads_emitted	12000	Total number of consensus reads (R1+R2=2) emitted
```

This format is compatible with fgbio's `CallMolecularConsensusReads` output.

### Filter Stats (Special Case)

The `filter --stats` output uses a two-column key-value format **without a header row**:

```text
total_reads	10000
passed_reads	8542
pass_rate	0.8542
```

## Duplex Metrics Terminology

The `duplex-metrics` command uses specific terminology for family types:

| Prefix | Name | Definition |
|--------|------|------------|
| **CS** | Coordinate-Strand | Families defined by genome coordinates and strand only (no UMI information) |
| **SS** | Single-Stranded | Families defined by coordinates, strand, and UMI. Two SS families from the same molecule (e.g., 50/A and 50/B) are counted separately |
| **DS** | Double-Stranded | Collapsed across SS families from the same molecule. SS families from opposite strands become one DS family |

The duplex-metrics output files include:

| File | Description |
|------|-------------|
| `<prefix>.family_sizes.txt` | Family size distribution by type (CS/SS/DS) |
| `<prefix>.duplex_family_sizes.txt` | Duplex family sizes by A→B and B→A strand counts |
| `<prefix>.duplex_yield_metrics.txt` | Summary QC metrics at subsampling levels (5%–100%) |
| `<prefix>.umi_counts.txt` | UMI observation frequencies |
| `<prefix>.duplex_umi_counts.txt` | Duplex UMI pair frequencies (optional, `--duplex-umi-counts`) |
| `<prefix>.duplex_qc.pdf` | QC plots (requires R with ggplot2) |

## Comparing Metrics

Use `fgumi compare metrics` to compare metrics files between runs:

```bash
fgumi compare metrics file1.txt file2.txt --precision 6 --rel-tol 1e-6
```

This is useful for validating that pipeline changes produce equivalent results. See the [compare documentation](https://github.com/fulcrumgenomics/fgumi/blob/main/docs/compare-cli.md) for details.
