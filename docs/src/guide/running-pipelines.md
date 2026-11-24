# Running Pipelines

Most fgumi work is a *pipeline*: extract UMIs, align, sort, group, call consensus, filter — or, when you want duplicate-marked reads rather than consensus reads, extract, align, sort, dedup. You can run each stage as its own command — writing a BAM to disk between every step — or fuse a run of stages into a single streaming process with `fgumi runall`.

## When to use `runall`

> **Rule of thumb:** run a single stage with its own command (`fgumi sort`, `fgumi group`, `fgumi simplex`, …). Chain **two or more** stages with `fgumi runall`.

`runall` streams records directly from one stage to the next, so every intermediate BAM a sequential run would write is elided. Its output is **record-identical** to running the equivalent standalone commands in turn — the same aligned records in the same order, and the same `@HD` `SO`/`GO`/`SS`, `@SQ`, and `@RG` header fields. The output is *not* byte-for-byte identical: the fused chain writes a single `@PG` record where the staged run writes one per command, and the `@HD` `VN` and `@CO` fields are not compared. Same algorithm, just without the round-trips to disk.

Reach for individual commands when you want to inspect or reuse an intermediate BAM (for example, to try several filter settings against one consensus BAM), or when you need a metrics file that `runall` does not emit, such as `simplex-metrics` or `duplex-metrics` (see [Caveats](#caveats)).

`runall` processes one input. If you have several BAMs from separate sequencing runs, combine them first with `fgumi merge --order template-coordinate`, then point `runall` at the merged file with `--start-from group`. Use `--start-from sort` only for raw aligned input that still needs sorting.

## The stage model

`runall` understands a fixed, linearly ordered set of stages, with two alternative branches after `sort`:

```text
extract < correct < align < zipper < sort < group < consensus < filter
                                            sort < dedup
```

A chain ends at either `dedup` or the group → consensus → filter branch, never both.

| Stage | What it does | Reads |
|-------|--------------|-------|
| `extract` | UMIs from FASTQ → unmapped BAM | FASTQ (`--extract::inputs`) — the only stage that reads FASTQ rather than a BAM |
| `correct` | Correct UMIs against a fixed whitelist | unmapped BAM (`--input`) |
| `align` | Run the aligner subprocess and zipper-merge the result | unmapped BAM (`--input`) + `--ref` + an aligner (`--aligner::preset` or `--aligner::command`) |
| `zipper` | Merge a mapped BAM with its unmapped BAM, restoring tags | mapped BAM (`--input`) + `--unmapped` + `--ref` (with a companion `.dict`) |
| `sort` | Sort into template-coordinate order | aligned BAM (`--input`) |
| `group` | Group reads by UMI, assigning MI tags | template-coordinate-sorted BAM (`--input`) |
| `consensus` | Call consensus (algorithm set by `--consensus`) | grouped, MI-tagged BAM (`--input`) |
| `filter` | Filter/mask consensus reads | consensus BAM (`--input`) |
| `dedup` | Mark (or remove) duplicates by UMI, instead of group → consensus → filter | template-coordinate-sorted BAM (`--input`) |

`--realign` adds a `realign` stage and a coordinate `sort` after the last of these; see [Re-aligning consensus reads](#re-aligning-consensus-reads-with---realign).

## Selecting the slice to run

Pick the slice with `--start-from` and `--stop-after` (both required). `--start-from` declares the state of your input; `--stop-after` names the last stage to run. They must satisfy `start-from ≤ stop-after` in the order above, and `--stop-after` cannot be `align` — it is start-only, so the valid stop stages are `extract`, `correct`, `zipper`, `sort`, `group`, `consensus`, `filter`, and `dedup`. `--stop-after dedup` accepts `--start-from` `extract`, `correct`, `align`, `zipper`, `sort`, or `dedup`; combining `dedup` with `group`, `consensus`, or `filter` in either flag is an error. Both are required because `runall` cannot infer either one: the same BAM could be valid input to several stages, and you may want to stop anywhere along the chain.

> **Important — match the input to `--start-from`.** `runall` validates the stage *range* but **cannot** verify that your input file is actually in the state `--start-from` claims. If you pass, for example, a raw unsorted BAM to `--start-from group`, the run may fail in a later stage or silently produce incorrect output, depending on what the mismatched input trips. Confirm your input matches the declared start stage (the [stage model](#the-stage-model) table lists what each stage expects).

When the slice reaches the `consensus` stage, choose the algorithm with `--consensus`:

| `--consensus` | Consensus caller | Required grouping strategy |
|---------------|------------------|----------------------------|
| `simplex` | single-strand ([Consensus Calling](consensus-calling.md)) | `adjacency` (or `identity`/`edit`) |
| `duplex` | double-strand ([Duplex Consensus Calling](duplex-consensus-calling.md)) | **`paired`** |
| `codec` | CODEC duplex from single read pairs | `adjacency` (or `identity`/`edit`) |

The algorithm is chosen by `--consensus`, **not** by the stage name. `--stop-after consensus` produces the same BAM as the matching `fgumi simplex`/`duplex`/`codec` command. **`--consensus duplex` requires `--group::strategy paired`** (so MIs carry the `/A`·`/B` strand suffixes duplex needs); `simplex` and `codec` require a non-paired strategy. Mismatching these is a hard error.

### Automatically inserted stages

You only name the *first* and *last* stages; `runall` fills in any stages in between that the chain requires. The one case that surprises people: the `align` and `zipper` stages emit reads in queryname order, but `group` needs template-coordinate order, so `runall` always inserts a `sort` step between them. For example, `--start-from extract --stop-after consensus` expands to extract → align → sort → group → consensus, `--start-from zipper --stop-after group` expands to zipper → sort → group, and `--start-from align --stop-after dedup` expands to align → sort → dedup (`dedup` needs template-coordinate order too). You never write the inserted `sort` yourself, but you can still tune it with `--sort::*` flags.

## Per-stage options

Tune any stage with prefixed flags of the form `--<stage>::<flag>`. These mirror the standalone command's options exactly:

| Standalone | `runall` |
|------------|----------|
| `fgumi sort --max-memory 4G` | `--sort::max-memory 4G` |
| `fgumi group --strategy adjacency` | `--group::strategy adjacency` |
| `fgumi duplex --min-reads 1,1,0` | `--duplex::min-reads 1,1,0` |
| `fgumi extract --read-structures +T +T` | `--extract::read-structures +T +T` |
| `fgumi dedup --remove-duplicates` | `--dedup::remove-duplicates` |

Cross-cutting options stay at the top level: `--threads`, `--compression-level`, `--max-memory`, `--ref`, `--rejects`, `--stats`, `--methylation-mode`. Run `fgumi runall --help` for the full, grouped list.

## Examples

Each example chains two or more stages, so each one is a job for `runall` rather than a sequence of standalone commands.

**Sort → group → simplex consensus**, replacing three commands and two intermediate BAMs:

```bash
fgumi runall \
  --start-from sort --stop-after consensus --consensus simplex \
  --input unsorted.bam --output consensus.bam \
  --group::strategy adjacency --group::edits 1 \
  --simplex::min-reads 1 \
  --threads 8
```

**Group → duplex consensus** (duplex requires the `paired` strategy):

```bash
fgumi runall \
  --start-from group --stop-after consensus --consensus duplex \
  --input sorted.bam --output duplex.bam \
  --group::strategy paired --group::edits 1 \
  --duplex::min-reads 1,1,0 \
  --threads 8
```

**Consensus → filter, fused** — no intermediate consensus BAM is written:

```bash
fgumi runall \
  --start-from group --stop-after filter --consensus simplex \
  --input sorted.bam --output filtered.bam \
  --group::strategy adjacency \
  --simplex::min-reads 1 \
  --filter::min-reads 1 --filter::max-read-error-rate 0.025 \
  --threads 8
```

**FASTQ → consensus**, the whole pipeline in one command (`runall` inserts `align`, `sort`, and `group`):

```bash
fgumi runall \
  --start-from extract --stop-after consensus --consensus simplex \
  --extract::inputs R1.fq.gz R2.fq.gz \
  --extract::read-structures 8M+T +T \
  --extract::sample MySample --extract::library MyLibrary \
  --ref ref.fa --aligner::preset bwa-mem3 --output consensus.bam \
  --group::strategy adjacency \
  --simplex::min-reads 1 \
  --threads 8
```

**FASTQ → duplicate-marked BAM**, for UMI data where consensus calling is not wanted (for example as a UMI-aware replacement for Picard `MarkDuplicates`; choose an aligner that suits the data, such as a splice-aware one through `--aligner::command` for RNA-seq). `runall` inserts `align` and `sort`, and the sorted reads stream straight into `dedup`:

```bash
fgumi runall \
  --start-from extract --stop-after dedup \
  --extract::inputs R1.fq.gz R2.fq.gz \
  --extract::read-structures 8M+T +T \
  --extract::sample MySample --extract::library MyLibrary \
  --ref ref.fa --aligner::preset bwa-mem3 --output deduped.bam \
  --dedup::strategy adjacency \
  --dedup::metrics deduped.metrics.txt \
  --threads 8
```

`--stop-after dedup` produces the same BAM as `fgumi dedup` on the sorted reads, and its metrics outputs (`--dedup::metrics`, `--dedup::family-size-histogram`, `--dedup::duplication-ladder`) are the same files. Add `--dedup::remove-duplicates` to drop duplicates instead of marking them. To deduplicate a BAM that is already template-coordinate sorted, use `--start-from dedup --stop-after dedup`.

## Re-aligning consensus reads with `--realign`

Consensus reads come out unmapped. `--realign` aligns them at the end of the same run: it appends a `realign` stage and a coordinate `sort` after the last stage, so a run ending at `consensus` or `filter` writes aligned, coordinate-sorted consensus reads (it is rejected for any other `--stop-after`). For example, `--start-from group --stop-after filter --realign` runs group → consensus → filter → realign → sort (coordinate), and `--start-from filter --stop-after filter --realign` filters an existing consensus BAM, then re-aligns and sorts the reads it keeps. Filter behaves as it does without `--realign`: it keeps the same reads, and `--filter::rejects` and `--filter::stats` are written as usual (the rejects are unaligned consensus reads).

`realign` is the standalone `fgumi realign` (see [Best Practices](best-practices.md#step-2a2-re-alignment)): for each end of a consensus read it keeps the aligner's primary if it lies inside the region the molecule's raw reads aligned to, and otherwise promotes the first in-window alternative from `XA`, or, failing that, swaps the primary with an in-window supplementary piece of the same read. The two ends of a pair the aligner flagged proper move together or not at all: no swap is made for one end unless its mate is also outside its own window with an in-window alternative of its own. That region (the `rw` tag) is recorded by group, so `--realign` implies `--group::realign-window` (and rejects an explicit `--group::realign-window false`); with `--start-from consensus` or `--start-from filter` the input must already come from `fgumi group --realign-window`.

- **Aligner.** `realign` uses the same aligner as the `align` stage: `--aligner::preset` or `--aligner::command`, `--aligner-bin`, `--ref` (with its index and a `.dict`), and `--methylation-mode`. These are required even when the run does not align raw reads. One exception: under `--methylation-mode` (`em-seq` or `taps`), duplex consensus reads are re-aligned with the plain preset (no `--meth`), since a duplex read's sequence is the molecule's own (its methylation is in `MM`/`ML`); simplex consensus reads keep their conversions and re-align conversion-aware with the `bwa-mem3` and `bwa-mem3-inproc` presets (`--meth`, or `--meth=taps`); the `bwa` preset re-aligns `taps` plain. With `--start-from filter`, pass `--consensus duplex` when the input is duplex consensus, since runall cannot otherwise tell; the duplex re-align then does not use `--methylation-mode`, so drop it unless one of filter's methylation options needs it. This applies to presets only: an `--aligner::command` runs as written for realign too, so for duplex input do not put `--meth` in it (runall warns if you do).
- **`--realign::xa-max-hits <primary>[,<alt>]`** sets bwa's `-h` for the realign stage only: only alternatives listed in `XA` can be promoted, so a larger value widens the candidates. bwa omits `XA` entirely when more than `-h` hits score near the best, so the most repetitive reads get no candidates at all. It requires `--realign`. Unset, the aligner keeps its own default (bwa's `5,200`) and an `--aligner::command` is run as written. Set, with `--aligner::command`, `-h` goes right after the command's first `mem` word, or at the end of the command (with a warning) if there is none; a command that already passes `-h` keeps its own (with a warning).
- **Per-base tags.** realign reverses (and reverse-complements) the per-base consensus tags of reverse-strand reads itself, as `--zipper::tags-to-reverse Consensus --zipper::tags-to-revcomp Consensus` would; `--zipper::tags-to-*` apply only to the raw-read merge of the `align` or `zipper` stage (runall warns when they are set with `--realign`). Do not pass `--reverse-per-base-tags` to a later `fgumi filter`.
- **Sorting.** The trailing sort always sorts by coordinate and reuses the `--sort::*` memory, temporary-directory and thread options. `--sort::order` does not apply to it; runall warns that the order is ignored only when you pass one.
- **Filter runs first.** With `--stop-after filter`, filter sees unaligned consensus reads and keeps, rejects and counts them as it does without `--realign`. Filter's methylation options follow the same rule as without `--realign`: accepted after `duplex`, rejected after `simplex` (see Methylation limits below). With `--start-from filter`, `--consensus` names the input's consensus type: `--consensus duplex` accepts them, and `--consensus simplex` (or `codec`) rejects them, since on single-strand consensus input they check nothing (the reads are not yet aligned). Without `--consensus` runall cannot tell simplex from duplex input, so it accepts them and warns; run `fgumi filter` after `fgumi realign` instead for single-strand input.
- **Memory.** A run with `--realign` holds more at once than one without. On a chain that also sorts raw reads there are two sorts, and each may use up to `--sort::max-memory` (the raw-read sort's merge can overlap the coordinate sort's buffer), so budget up to twice that. The reference is read on demand: realign reads it only to rescore the alignments it promotes, which are few, so it reads just those regions through the FASTA's `.fai` index rather than holding the FASTA (on hg38 with in-process bwa-mem3, about 3.3 GB lower peak memory than holding it, with identical output). A compressed (gzip or BGZF) FASTA, one without a `.fai` (`samtools faidx`), or one that is not a regular file (a pipe, such as `<(cat ref.fa)`), cannot be read on demand, so realign loads it into memory (decompressing a compressed one) with a warning; `--realign::in-memory-reference` loads it into memory regardless. When a methylation-aware consensus caller already holds the same `--ref` in memory, realign shares that copy. The aligner index: with `--aligner::preset bwa-mem3-inproc` one index is shared with the `align` stage, except for duplex consensus under `--methylation-mode`, where the `align` stage loads the `bwa-mem3 index --meth` dual index and `realign` (which aligns duplex reads plain) the plain index, so both are resident for the whole run. A subprocess preset (`bwa-mem3`, `bwa`) or `--aligner::command` on a chain that also aligns raw reads runs two aligner processes for the whole run, each holding its own copy of the index. Budget for these on top of the usual run.

## `runall` vs. individual commands

| | `runall` | Individual commands |
|---|----------|---------------------|
| Intermediate BAMs | none (in-memory) | written to disk between stages |
| Output | record-identical to the command chain (header `@PG`/`@HD VN` differ) | the same records |
| Per-stage metrics/histograms | the stages' own, opt-in (`--correct::metrics`, `--group::metrics`, `--dedup::metrics`, or `--all-metrics`); not `simplex-metrics`/`duplex-metrics` | available (e.g. `fgumi group --metrics`) |
| Inspect/reuse an intermediate | no | yes — the BAM is on disk |
| Best for | production runs, multi-stage chains | single stages, debugging, metrics, experimentation |

For the recommended end-to-end pipelines and parameter choices, see [Best Practices](best-practices.md).

## Caveats

- **Only the stages' own metrics.** A fused chain writes the metrics a stage exposes as a `--<stage>::` option (`--correct::metrics`, `--group::metrics`, dedup's `--dedup::metrics`, `--dedup::family-size-histogram` and `--dedup::duplication-ladder`), each only when asked for. `--all-metrics` derives the standard metrics and histogram paths, but dedup's duplication ladder remains opt-in via `--dedup::duplication-ladder`. It does not run the separate QC commands `simplex-metrics` and `duplex-metrics`; run those on the output when you need them — see [Working with Metrics](working-with-metrics.md).
- **`dedup` replaces the consensus branch.** A chain ends at either `dedup` or group → consensus → filter, so `--start-from group|consensus|filter --stop-after dedup` and `--start-from dedup --stop-after group|consensus|filter` are rejected, as is `--realign` with `--stop-after dedup` (there are no consensus reads to re-align). Dedup is configured only by `--dedup::*`: `--group::*` options do not apply to it, and `--consensus` is ignored (with a warning).
- **No `--stop-after align`.** Raw aligner output before the zipper-merge has lost every original tag, so the earliest stop after `--start-from align` is `zipper`.
- **Input state is not validated** against `--start-from` — see the warning under [Selecting the slice to run](#selecting-the-slice-to-run).
- **UMI rejects need `--correct::rejects` when `correct` is fused mid-chain.** On a chain that runs past `correct`, pass `--correct::rejects <path>` to capture correct's UMI rejects. The top-level `--rejects` does not capture them there: on a chain that reaches a consensus stage it is that stage's rejects file, and on one that does not it is unused (the run warns). On a `correct` self-pair (`--stop-after correct`), the top-level `--rejects` takes precedence over `--correct::rejects`.
- **Strategy constraints.** `--consensus duplex` requires `--group::strategy paired`; `--consensus simplex` and `--consensus codec` require a non-paired strategy.
- **Methylation limits.** On a chain that includes the `align` stage (`--start-from align`, or an `extract`/`correct` chain, which inserts it), `--methylation-mode` makes alignment conversion-aware: the `bwa-mem3` and `bwa-mem3-inproc` presets align with `bwa-mem3 mem --meth --meth-scoring genomic` under `em-seq` and `bwa-mem3 mem --meth=taps --meth-scoring neutral` under `taps`, and both need the dual index built by `bwa-mem3 index --meth <ref>`. From bwa-mem3 0.18.0 that command also builds the plain index and refuses to replace an existing one, so beside a plain index run `bwa-mem3 index --meth --force <ref>`; earlier releases do not accept `--force`, and do not need it. The `bwa-mem3` preset also accepts a methylation-only index (`bwa-mem3 index --meth --meth-only <ref>`, which skips the plain `<ref>.bwt.2bit.64`); `bwa-mem3-inproc` needs the full dual index, since it also loads the plain one. Duplex consensus reads re-align without `--meth` (`--realign`, `fgumi realign`), so that step needs the plain `<ref>.bwt.2bit.64`, and a `--meth-only` index is rejected there at startup. `--meth=taps` places TAPs reads better than plain alignment and writes `XM` calls with the TAPs meaning (a bare `--meth` would invert them); see the [Methylation](methylation.md) guide. The `bwa` preset has no conversion-aware mode: it is rejected under `em-seq` and aligns `taps` plain. With `--aligner::command`, the command runs as written, so it must pass the mode itself (`bwa-mem3 mem --meth …` or `bwameth.py …` for EM-seq, `bwa-mem3 mem --meth=taps …` for TAPs, with the `=`). `--methylation-mode` is unsupported with `--consensus codec`. Filter's methylation options (`--filter::min-methylation-depth`, `--filter::require-strand-methylation-agreement`, `--filter::min-conversion-fraction`) need aligned simplex reads. runall cannot re-align between the consensus caller and filter (`--realign` re-aligns after filter), so it rejects them on a chain where filter follows `simplex`, with or without `--realign`; re-align the consensus reads (for example with `fgumi realign`) and run `--start-from filter` instead. Duplex reads carry their methylation sites in their own sequence, so after `duplex` the options are accepted. On a chain that reaches a consensus stage they also need `--methylation-mode` (and `--ref`): without it the consensus caller writes no `cu`/`ct` counts, so the options are rejected. See [Methylation](methylation.md).
- **`zipper` start needs a dictionary.** `--start-from zipper` requires `--ref` with a companion `.dict` file (create one with `samtools dict`).

## See Also

- [Getting Started](getting-started.md) — the step-by-step tutorial these stages come from
- [Best Practices](best-practices.md) — recommended end-to-end pipelines and parameters
- [Troubleshooting](troubleshooting.md) — fixes for common `runall` and sort-order errors
- [Glossary](glossary.md) — UMI, template-coordinate, strand labels, and BAM tags
