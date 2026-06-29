# Duplex Consensus Calling

## Overview

Duplex consensus calling takes reads generated from both strands of a double-stranded source molecule and produces consensus reads with extremely low error rates. This is the process used in duplex sequencing methods such as those described by [Kennedy et al](https://www.ncbi.nlm.nih.gov/pubmed/25299156), where UMIs are attached to each end of the source molecule.

The mathematical model is similar to [single-strand consensus calling](consensus-calling.md), but the mechanics differ because reads from both strands must be combined.

Duplex consensus calling is run *after* grouping reads with `fgumi group --strategy paired`.

## Process

Starting from a group of reads identified as originating from the same double-stranded molecule, the two strands are labeled A and B. The process proceeds through these steps:

1. **Split** reads into four sub-groups: A1 (strand A, read 1), A2, B1, B2
2. **Unmap** and revert to sequencing order
3. **Quality trim** (optional, recommended)
4. **Mask** remaining low-quality bases to `N`
5. **Trim to insert length** to avoid reading into adapters
6. **Filter by CIGAR** to ensure reads are in phase
7. **Call four single-strand consensus reads** (one each for A1, A2, B1, B2)
8. **Call two duplex consensus reads** by combining A1+B2 and A2+B1

### Splitting Reads into Groups

Reads are split by their strand sub-family (`A` or `B`, assigned by `fgumi group --strategy paired` from sequencing order — see [Tracking Reads](tracking-reads.md)) and by whether they are sequencing read 1 or 2. R1s from the `A` sub-family correspond to R2s from the `B` sub-family, and vice versa, which is why the duplex consensuses are formed from A1+B2 and A2+B1.

### Quality Trimming

Pass `--trim` to end-trim low-quality bases before consensus calling. This is highly recommended: it reduces consensus disagreements and no-calls (`N`s). Trimming uses the same running-sum algorithm as BWA.

### Masking Low-Quality Bases

Bases below the minimum quality threshold are converted to `N`s so they are not used in consensus calling. If quality trimming is disabled, reads are truncated to remove contiguous trailing `N`s.

### Trimming to Insert Length

Reads longer than the insert length read into adapter sequence. For duplex data, A1 and B2 reads may read into *different* adapter sequences. Calling consensus across different adapters produces many disagreements and no-calls, potentially causing consensus reads to be erroneously filtered. fgumi therefore trims reads to the insert length before consensus calling.

### CIGAR Filtering

Without multiple alignment, length errors (indels) in raw reads cause reads to be out of phase with each other. For example:

```text
1: ACGTGACTGACTAGCTTTTTTT-AGACTAGCTACTACT
2: ACGTGACTGACTAGCTTTTTTT-AGACTAGCTACTACT
3: ACGTGACTGACTAGCTTTTTTTT-GACTAGCTACTACT
```

Read 3 has an extra T, causing many disagreements with reads 1 and 2.

To handle this, reads are grouped by compatible CIGAR alignments, and only the largest group is used for consensus. This is performed independently on A1+B2 and B1+A2 reads.

### Calling Single-Strand Consensus Reads

Four single-strand consensus reads are generated (A1, A2, B1, B2) using the standard [consensus calling model](consensus-calling.md).

### Calling Duplex Consensus Reads

The final duplex R1 and R2 are produced by merging the appropriate A and B reads base-by-base:

1. **Bases agree**: quality = Q(A) + Q(B)
2. **Bases disagree, different qualities**: base = higher quality base, quality = Q(higher) - Q(lower)
3. **Bases disagree, same quality**: base is arbitrarily from A, quality = 2 (minimum Phred score)

## The min-reads Parameter

### For Simplex Consensus

`fgumi simplex` accepts a single `--min-reads` value.

### For Duplex Consensus

`fgumi duplex` and `fgumi filter` accept one, two, or three `--min-reads` values. When you supply fewer than three, the last is repeated — `80 40` becomes `80 40 40`, and `10` becomes `10 10 10`.

| Position | Applies to | Meaning |
|----------|-----------|---------|
| 1st | the final **duplex** read | minimum total raw reads across both strands |
| 2nd | the **better-supported** single-strand consensus | minimum reads for the strand with *more* support |
| 3rd | the **less-supported** single-strand consensus | minimum reads for the strand with *less* support |

The 2nd value must be at least the 3rd: it applies to the strand that has *more* reads, so it makes no sense to demand more reads from the weaker strand than the stronger one.

**Example:** `--min-reads 7 3 1` requires at least 7 total raw reads supporting the duplex consensus, at least 3 raw reads for the better-supported single-strand consensus, and at least 1 raw read for the other.

## See Also

- [Consensus Calling](consensus-calling.md) — the underlying single-strand model
- [Tracking Reads](tracking-reads.md) — how `/A`·`/B` strand labels and the per-strand tags are assigned
- [UMI Grouping](umi-grouping.md) — the `paired` strategy that duplex requires
- [Working with Metrics](working-with-metrics.md) — duplex QC metrics
