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

Reads are split by strand of origin (A or B) and whether they are sequencing read 1 or 2. R1s from strand A correspond to R2s from strand B, and vice versa.

### Quality Trimming

Reads can be end-trimmed to remove low-quality bases. This is highly recommended as it reduces disagreements in the consensus and fewer no-calls (`N`s). Trimming uses the same running-sum algorithm as BWA.

### Masking Low-Quality Bases

Bases below the minimum quality threshold are converted to `N`s so they are not used in consensus calling. If quality trimming is disabled, reads are truncated to remove contiguous trailing `N`s.

### Trimming to Insert Length

Reads longer than the insert length read into adapter sequence. For duplex data, A1 and B2 reads may read into *different* adapter sequences. Calling consensus across different adapters produces many disagreements and no-calls, potentially causing consensus reads to be erroneously filtered. Reads are therefore trimmed to insert length before consensus calling.

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

`fgumi simplex` and `fgumi filter` accept a single `--min-reads` value.

### For Duplex Consensus

`fgumi duplex` and `fgumi filter` accept one, two, or three `--min-reads` values. If fewer than three values are supplied, the last is repeated (e.g. `80 40` becomes `80 40 40`, `10` becomes `10 10 10`).

The values control:
1. **First value**: minimum total raw reads across both single-strand consensuses for the final duplex read
2. **Second value**: minimum reads for the single-strand consensus with *more* support
3. **Third value**: minimum reads for the single-strand consensus with *less* support

If values two and three differ, the more stringent value must come first.

**Example:** `--min-reads 7 3 1` requires:
- At least 7 total raw reads supporting the duplex consensus
- At least 3 raw reads for the better-supported single-strand consensus
- At least 1 raw read for the other single-strand consensus
