# UMI Grouping

## Overview

`fgumi group` assigns reads that appear to come from the same original molecule to the same group by writing a shared Molecular Identifier (`MI`) tag. Grouping relies on template-coordinate sort order.

This page describes:
1. How reads and templates are filtered before grouping
2. How mapping coordinates and UMIs identify reads from the same molecule
3. Template-coordinate sort order
4. Cell barcode support
5. Metrics output

## Filtering Reads and Templates

A *read* is a single sequenced strand. A *template* is all reads sharing the same query name (typically a read pair).

| Concept | Definition | Example |
|---------|------------|---------|
| **Read** | A single sequenced strand (R1 or R2) | `@read123/1` |
| **Template** | The full fragment, represented by both reads in a pair | `@read123` includes both /1 and /2 |

Reads and templates are filtered before grouping to prevent splitting reads from a single molecule into separate groups.

**Individual reads** are filtered if:
- Flagged as secondary (unless `--include-secondary`)
- Flagged as supplementary (unless `--include-supplementary`)

**All reads for a template** are filtered if:
- All reads for the template are unmapped (unless `--allow-unmapped`)
- Any non-secondary, non-supplementary read has mapping quality < `--min-map-q`
- Any UMI sequence contains one or more `N` bases
- `--min-umi-length` is specified and the UMI does not meet the length requirement

### Grouping Unmapped Reads

By default, templates where all reads are unmapped are excluded from grouping. Pass `--allow-unmapped`
to include them. This is useful for workflows where some templates genuinely fail to align
(e.g. cell-free DNA fragments that fall outside the target region) but should still be counted
and may share UMIs with mapped templates from the same molecule:

```bash
fgumi group \
  --input sorted.bam \
  --output grouped.bam \
  --strategy adjacency \
  --allow-unmapped
```

## Grouping Strategies

Grouping is performed by one of four strategies:

### identity

Groups only reads with identical UMI sequences. It is the simplest and fastest strategy, but a single sequencing error in a UMI splits one molecule into separate groups — so use it only for data exploration, not production grouping.

### edit

Reads are clustered into groups such that each read within a group has at least one other read in the group with <= `--edits` differences, and there are no inter-group pairings with <= `--edits` differences. Effective when there are small numbers of reads per UMI, but breaks down at very high UMI coverage.

### adjacency

A version of the directed adjacency method described in [umi_tools](http://dx.doi.org/10.1101/051755). It merges UMIs that differ by a few bases, but only when the read counts form a *count gradient* — a high-count UMI absorbs a near-identical low-count one (the low-count UMI is treated as a sequencing error of the abundant parent), while two similarly abundant UMIs are kept separate (they are likely distinct molecules). **Recommended for most simplex and CODEC workflows.**

### paired

Similar to adjacency but for duplex sequencing where each template has two UMIs (one from each strand). Expects UMI sequences stored in a single tag separated by a hyphen (e.g. `ACGT-CCGG`). Allows one UMI to be absent (e.g. `ACGT-` or `-ACGT`).

The molecular IDs produced have structure: `{base}/{A|B}`. For example, UMI pairs `AAAA-GGGG` and `GGGG-AAAA` map to `1/A` and `1/B` respectively. See [Tracking Reads](tracking-reads.md) for details. **Recommended for duplex workflows.**

The `edit`, `adjacency`, and `paired` strategies use the `--edits` parameter to control matching of non-identical UMIs.

## Cell Barcode Support

When processing data with cell barcodes (e.g. single-cell sequencing), reads at the same genomic
position are partitioned by cell barcode *before* UMI assignment. This ensures that reads from
different cells are never grouped together, even if they share a UMI and mapping position.

The cell barcode is read from the standard `CB` tag. No correction or
error-handling is performed on cell barcodes — they must be corrected upstream before grouping.

Cell barcodes are detected automatically across the entire pipeline — no additional flags are
needed. The consensus callers validate that all source reads in a group share the same cell
barcode and propagate it to the output consensus read.

## Metrics Output

`fgumi group` can emit three types of metrics files. They can be specified individually or all at
once with the `--metrics` prefix flag.

### Using `--metrics` (recommended)

The `-M`/`--metrics` flag writes all three metrics files under a single prefix in one step:

```bash
fgumi group \
  --input sorted.bam \
  --output grouped.bam \
  --strategy adjacency \
  --metrics my_sample
```

This produces:
- `my_sample.family_sizes.txt` — histogram of UMI family sizes
- `my_sample.grouping_metrics.txt` — overall grouping statistics
- `my_sample.position_group_sizes.txt` — histogram of UMI families per genomic position

### Using individual flags

The three metrics can also be written to explicit paths:

```bash
fgumi group \
  --input sorted.bam \
  --output grouped.bam \
  --strategy adjacency \
  --family-size-histogram family_sizes.txt \
  --grouping-metrics grouping_metrics.txt
```

Note: `position_group_sizes.txt` is only available via `--metrics`. The individual flags
`--family-size-histogram` and `--grouping-metrics` can be used alongside `--metrics`.

### Family sizes

The `family_sizes.txt` file is a histogram of how many templates belong to each UMI family (on paired-end data a read pair counts as one template). A large
fraction of singleton families usually indicates insufficient sequencing depth relative to library
complexity, or a UMI extraction error (the wrong bases ending up in `RX`).

### Grouping metrics

The `grouping_metrics.txt` file contains summary statistics about the grouping run, including
total reads, accepted reads, discarded reads by reason, and UMI assignment counts.

### Position group sizes

The `position_group_sizes.txt` file is a histogram of how many distinct UMI families were
observed at each unique genomic position (coordinate + strand). A distribution skewed toward
large position groups may indicate high on-target duplication or UMI exhaustion.

## Template-Coordinate Sort Order

`fgumi group` requires its input to be template-coordinate sorted. The header must advertise
`SO:unsorted`, `GO:query`, and `SS:template-coordinate`; without `SS:template-coordinate` the
input is treated as queryname-grouped (e.g. FASTQ-order output from `fgumi extract`) and
rejected with an actionable error pointing back here. `fgumi group` does **not** sort
internally — pre-sort with:

```bash
fgumi sort --order template-coordinate --input aligned.bam --output sorted.bam
```

The streaming grouper relies on records that share a position key being consecutive in the
input, which is what template-coordinate sort guarantees. Any other ordering (queryname,
coordinate, FASTQ-order) would split each true molecule across many small groups and assign
distinct `MI` values to reads that should share one.

For single-cell data, the `CB` cell barcode tag is automatically incorporated in the sort key,
keeping templates from different cells at the same locus separate:

```bash
fgumi sort --order template-coordinate --input aligned.bam --output sorted.bam
```

Template-coordinate order sorts reads by:
1. The earlier unclipped 5' coordinate of the read pair
2. The higher unclipped 5' coordinate of the read pair
3. Strand orientation
4. The cellular barcode (CB tag, if present)
5. The molecular identifier (MI tag, if present)
6. Read name
7. Library (from read group)
8. Whether R1 has the lower coordinates of the pair

Reads grouped by `fgumi group` with the same MI will share the same outer start/stop coordinates.
Because 5' coordinates are strand-aware, reads from opposite strands with the same UMI and
position will not be grouped together (they belong to different strands of the same duplex
molecule).

## See Also

- [Consensus Calling](consensus-calling.md) and [Duplex Consensus Calling](duplex-consensus-calling.md) — the next stage after grouping
- [Tracking Reads](tracking-reads.md) — the `MI` tag and `/A`·`/B` strand labels grouping assigns
- [Working with Metrics](working-with-metrics.md) — interpreting the grouping metrics
- [Best Practices](best-practices.md) — grouping in the recommended workflow
