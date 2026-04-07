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

Only reads with identical UMI sequences are grouped together. This is simpler and faster than other strategies, but should usually be avoided because sequencing errors in the UMI will split reads from the same molecule into separate groups. Useful for data exploration.

### edit

Reads are clustered into groups such that each read within a group has at least one other read in the group with <= `--edits` differences, and there are no inter-group pairings with <= `--edits` differences. Effective when there are small numbers of reads per UMI, but breaks down at very high UMI coverage.

### adjacency

A version of the directed adjacency method described in [umi_tools](http://dx.doi.org/10.1101/051755) that allows for errors between UMIs but only when there is a count gradient. **Recommended for most simplex and CODEC workflows.**

### paired

Similar to adjacency but for duplex sequencing where each template has two UMIs (one from each strand). Expects UMI sequences stored in a single tag separated by a hyphen (e.g. `ACGT-CCGG`). Allows one UMI to be absent (e.g. `ACGT-` or `-ACGT`).

The molecular IDs produced have structure: `{base}/{A|B}`. For example, UMI pairs `AAAA-GGGG` and `GGGG-AAAA` map to `1/A` and `1/B` respectively. See [Tracking Reads](tracking-reads.md) for details. **Recommended for duplex workflows.**

The `edit`, `adjacency`, and `paired` strategies use the `--edits` parameter to control matching of non-identical UMIs.

## Cell Barcode Support

When processing data with cell barcodes (e.g. single-cell sequencing), reads at the same genomic
position are partitioned by cell barcode *before* UMI assignment. This ensures that reads from
different cells are never grouped together, even if they share a UMI and mapping position.

The cell barcode is read from the tag specified by `--cell-tag` (default: `CB`). No correction or
error-handling is performed on cell barcodes — they must be corrected upstream before grouping.

Pass `--cell-tag` consistently through the entire pipeline:

```bash
fgumi sort   --order template-coordinate --cell-tag CB ...
fgumi group  --cell-tag CB ...
fgumi simplex --cell-tag CB ...
```

The consensus callers validate that all source reads in a group share the same cell barcode and
propagate it to the output consensus read.

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

The `family_sizes.txt` file is a histogram of how many reads belong to each UMI family. A large
fraction of singleton families may indicate UMI collisions, over-sequencing, or UMI extraction
errors.

### Grouping metrics

The `grouping_metrics.txt` file contains summary statistics about the grouping run, including
total reads, accepted reads, discarded reads by reason, and UMI assignment counts.

### Position group sizes

The `position_group_sizes.txt` file is a histogram of how many distinct UMI families were
observed at each unique genomic position (coordinate + strand). A distribution skewed toward
large position groups may indicate high on-target duplication or UMI exhaustion.

## Template-Coordinate Sort Order

If the input is not sorted in template-coordinate order, `fgumi group` will internally sort the
reads. To avoid this overhead, pre-sort with:

```bash
fgumi sort --order template-coordinate --input aligned.bam --output sorted.bam
```

For single-cell data, include `--cell-tag CB` so that the sort key incorporates the cell barcode,
keeping templates from different cells at the same locus separate:

```bash
fgumi sort --order template-coordinate --cell-tag CB --input aligned.bam --output sorted.bam
```

Template-coordinate order sorts reads by:
1. The earlier unclipped 5' coordinate of the read pair
2. The higher unclipped 5' coordinate of the read pair
3. Strand orientation
4. The cellular barcode (CB tag, if `--cell-tag` is set)
5. The molecular identifier (MI tag, if present)
6. Read name
7. Library (from read group)
8. Whether R1 has the lower coordinates of the pair

Reads grouped by `fgumi group` with the same MI will share the same outer start/stop coordinates.
Because 5' coordinates are strand-aware, reads from opposite strands with the same UMI and
position will not be grouped together (they belong to different strands of the same duplex
molecule).

See also: [Consensus Calling](consensus-calling.md), [Duplex Consensus Calling](duplex-consensus-calling.md), [Best Practices](best-practices.md)
