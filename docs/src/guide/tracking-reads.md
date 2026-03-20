# Tracking Reads through Grouping and Duplex Consensus Calling

This guide describes conventions for tracking reads from raw data through grouping and duplex consensus calling. It covers how molecular identifiers relate to strand assignment and how consensus tags encode single-strand and duplex information.

## Top and Bottom Strand for Raw Reads

`fgumi group` assigns the same molecular ID to raw reads from the same source molecule, with trailing `/A` and `/B` to indicate which strand they belong to (top or bottom, AB or BA).

**Convention:** The `/A` raw reads are those where the 5' unclipped position of read one (of the pair) is less than or equal to the 5' unclipped position of read two. The 5' unclipped position is relative to sequencing order, not the reference genome strand.

For example:

```
x: R1----------------->    <-------------------R2
y: R2----------------->    <-------------------R1
z: R1----------------->
     <-----------------R2
```

- `x` gets `/A` (R1's 5' end is at or before R2's 5' end)
- `y` gets `/B` (R1's 5' end is after R2's 5' end in sequencing order)
- `z` gets `/A` (even though fully overlapped, R1's 5' end is earlier)

## Single-Strand Reads Relative to Duplex Consensus

`fgumi duplex` writes single-strand information into SAM tags for each duplex consensus read. Which single-strand consensus goes into the "AB" vs "BA" tags is determined as follows:

1. **Both strands present:** Information for raw reads with `/A` in their molecular ID goes into "AB" tags; `/B` reads go into "BA" tags.
2. **Only one strand present:** The "AB" tags contain the single-strand consensus that was generated. The "BA" tags contain only per-read tags (no consensus data).

The duplex consensus sequence has the same strand orientation as the "AB" single-strand consensus.

## Consensus Tags

SAM tags used for single-strand and duplex consensus reads:

| Value | AB Tag | BA Tag | Final Tag |
|-------|--------|--------|-----------|
| Per-read depth | `aD` | `bD` | `cD` |
| Per-read min depth | `aM` | `bM` | `cM` |
| Per-read error rate | `aE` | `bE` | `cE` |
| Per-base depth | `ad` | `bd` | `cd` |
| Per-base error count | `ae` | `be` | `ce` |
| Per-base bases | `ac` | `bc` | (bases) |
| Per-base quals | `aq` | `bq` | (quals) |

**Convention:** The second letter in the tag is lowercase for per-base values and uppercase for per-read values.
