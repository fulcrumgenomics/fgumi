# Troubleshooting

Common errors and unexpected results, with their causes and fixes.

## `fgumi extract` rejects the read structures

**Cause.** The number of read structures does not match the number of `--inputs` FASTQs, or a structure's fixed-length segments add up to more bases than a read contains.

**Fix.** Supply exactly one read structure per FASTQ, in the same order as `--inputs`, and use `+` for the final segment so it absorbs whatever bases remain (e.g. `8M+T` rather than `8M142T` when read length varies). See [Read Structures](read-structures.md) for the syntax and worked examples.

## "Input BAM must be template-coordinate sorted"

```text
Input BAM must be template-coordinate sorted (header must advertise
SO:unsorted, GO:query, and SS:template-coordinate).
```

**Cause.** `fgumi group`, `fgumi dedup`, and `fgumi downsample` require template-coordinate order, and the input is in some other order.

**Fix.** Sort with fgumi (not samtools — see below):

```bash
fgumi sort -i input.bam -o sorted.bam --order template-coordinate
```

## `samtools sort` output gives wrong grouping or dedup results

**Cause.** `fgumi zipper` writes a `tc` tag on secondary/supplementary reads so they sort and deduplicate correctly. `samtools sort --template-coordinate` ignores the `tc` tag, so its output orders those reads incorrectly for `fgumi group`/`dedup`.

**Fix.** Always sort with `fgumi sort --order template-coordinate` after `fgumi zipper`. Reserve `samtools sort` for plain coordinate or queryname sorting of BAMs that will not feed `group`, `dedup`, or `downsample`.

## Every UMI family has size 1 (all singletons)

**Cause.** Reads from the same molecule are not being grouped together. Common reasons:

- The UMI length in `--read-structures` is wrong, so the wrong bases are stored in `RX`.
- UMI diversity is high relative to depth (genuinely few duplicates).
- The `identity` grouping strategy is splitting families on sequencing errors in the UMI.

**Fix.** Verify the read structure against your library design ([Read Structures](read-structures.md)); switch from `identity` to `adjacency` so single-base UMI errors are tolerated ([UMI Grouping](umi-grouping.md)); and inspect the family-size histogram from `fgumi group --metrics` ([Working with Metrics](working-with-metrics.md)).

## Low consensus yield

**Cause.** `--min-reads` is set higher than most families can satisfy, or families are too small.

**Fix.** Lower `--min-reads` (it is always safe to call with `--min-reads 1` and filter later), confirm UMI extraction, and run `fgumi simplex-metrics` / `fgumi duplex-metrics` on the grouped BAM to see how yield scales with depth.

## Excessive duplex strand disagreement / many no-calls

**Cause.** The two strands disagree at many positions, often from low input base quality or reads extending past the insert.

**Fix.** Raise `--min-input-base-quality`, enable `--trim` to remove low-quality read ends, and tune the per-strand `--min-reads` triple `[duplex, AB, BA]`. See [Duplex Consensus Calling](duplex-consensus-calling.md).

## Out-of-memory / process killed

**Cause.** `--max-memory` bounds only the pipeline queue, and it is *per thread* by default, so total RSS scales with `--threads`. A high thread count can exceed host RAM.

**Fix.** On a fixed-RAM host, use `--max-memory auto` (self-throttles to the host) or set a fixed total budget with `--max-memory <N> --memory-per-thread false`; reduce `--threads`. See the [Performance Tuning](performance-tuning.md) memory section for the full model.

## `fgumi runall` produces wrong or empty output (input doesn't match `--start-from`)

**Cause.** `runall` cannot verify that your input file is actually in the state `--start-from` claims. If the input is in a different state — e.g. a raw unsorted BAM passed to `--start-from group`, or an ungrouped BAM passed to `--start-from consensus` — the pipeline runs on whatever it finds and silently produces wrong (or empty) results rather than erroring.

**Fix.** Make the input match the declared start stage. Check the BAM's sort order in its header (`samtools view -H in.bam | grep @HD`) and confirm the expected tags are present (`RX` after extract, `MI` after group). The [stage model table](running-pipelines.md#the-stage-model) lists what each start stage expects. When in doubt, start from an earlier stage so `runall` does the upstream work itself.

## `--stop-after align` is rejected (runall)

**Cause.** Raw aligner output before the zipper-merge has lost every original tag, so `align` is not a valid stop point.

**Fix.** Use `--stop-after zipper` as the earliest stop after `--start-from align`. See [Running Pipelines](running-pipelines.md).

## `fgumi group` rejects `--strategy paired` with another flag

```text
Paired strategy cannot be used with --min-umi-length
--no-umi cannot be used with --strategy paired
```

**Cause.** The `paired` strategy (for duplex data) is incompatible with `--min-umi-length` and `--no-umi`.

**Fix.** Drop the conflicting flag. Use `paired` only for duplex workflows; use `adjacency` for simplex and CODEC.
