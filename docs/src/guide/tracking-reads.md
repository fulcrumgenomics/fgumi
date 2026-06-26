# Tracking Reads through Grouping and Duplex Consensus Calling

This guide describes conventions for tracking reads from raw data through grouping and duplex consensus calling. It covers how molecular identifiers relate to strand assignment and how consensus tags encode single-strand and duplex information.

You need this page when you filter or interpret duplex consensus reads by strand, query the per-strand consensus tags in a downstream tool, or debug why reads landed in the wrong molecule. For routine pipelines, the defaults handle all of this automatically.

## Strand Labels (`/A` and `/B`) for Raw Reads

`fgumi group --strategy paired` assigns the same molecular ID base to raw reads from one source molecule, with a trailing `/A` or `/B` marking which of the two strands a read pair came from.

**Convention:** a read pair is labeled `/A` when read 1's 5′ end comes at or before read 2's 5′ end *in sequencing order* (ignoring soft-clipping and independent of the reference strand); otherwise it is `/B`.

For example:

```text
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

### Reading a duplex consensus read's tags

For a duplex consensus read built from molecule `1` (raw reads tagged `1/A` and `1/B`):

- `cD` is the final duplex depth; `aD` and `bD` are the depths of the `/A` and `/B` single-strand consensuses that were combined.
- A high `cD` with a small value for one of the single-strand consensuses means the duplex rests on one well-covered strand and one thinly-covered strand — exactly the case the per-strand `--min-reads` thresholds control (their order is by support, not by `/A` vs `/B`). See [Duplex Consensus Calling](duplex-consensus-calling.md).
- The per-base `ad`/`ae` (and `bd`/`be`) arrays let `fgumi filter` mask individual low-support or high-error bases.

## See Also

- [Duplex Consensus Calling](duplex-consensus-calling.md) — how the `/A`·`/B` strands become a duplex consensus
- [UMI Grouping](umi-grouping.md) — where the `MI` tags and strand labels are assigned
- [Working with Metrics](working-with-metrics.md) — interpreting the resulting QC metrics
