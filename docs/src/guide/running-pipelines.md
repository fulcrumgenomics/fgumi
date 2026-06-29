# Running Pipelines

Most fgumi work is a *pipeline*: extract UMIs, align, sort, group, call consensus, filter. You can run each stage as its own command — writing a BAM to disk between every step — or fuse a run of stages into a single streaming process with `fgumi runall`.

## When to use `runall`

> **Rule of thumb:** run a single stage with its own command (`fgumi sort`, `fgumi group`, `fgumi simplex`, …). Chain **two or more** stages with `fgumi runall`.

`runall` streams records directly from one stage to the next, so every intermediate BAM a sequential run would write is elided. Its output is **byte-for-byte identical** to running the equivalent standalone commands in turn — the same algorithm, just without the round-trips to disk.

Reach for individual commands when you want to inspect or reuse an intermediate BAM (for example, to try several filter settings against one consensus BAM), or when you need a per-stage metrics file that `runall` does not emit (see [Caveats](#caveats)).

`runall` processes one input. If you have several BAMs from separate sequencing runs, combine them first with `fgumi merge --order template-coordinate`, then point `runall` at the merged file with `--start-from group`. Use `--start-from sort` only for raw aligned input that still needs sorting.

## The stage model

`runall` understands a fixed, linearly ordered set of stages:

```text
extract < correct < align < zipper < sort < group < consensus < filter
```

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

## Selecting the slice to run

Pick the slice with `--start-from` and `--stop-after` (both required). `--start-from` declares the state of your input; `--stop-after` names the last stage to run. They must satisfy `start-from ≤ stop-after` in the order above, and `--stop-after` cannot be `align` — it is start-only, so the valid stop stages are `extract`, `correct`, `zipper`, `sort`, `group`, `consensus`, and `filter`. Both are required because `runall` cannot infer either one: the same BAM could be valid input to several stages, and you may want to stop anywhere along the chain.

> **Important — match the input to `--start-from`.** `runall` validates the stage *range* but **cannot** verify that your input file is actually in the state `--start-from` claims. If you pass, for example, a raw unsorted BAM to `--start-from group`, the run may fail in a later stage or silently produce incorrect output, depending on what the mismatched input trips. Confirm your input matches the declared start stage (the [stage model](#the-stage-model) table lists what each stage expects).

When the slice reaches the `consensus` stage, choose the algorithm with `--consensus`:

| `--consensus` | Consensus caller | Required grouping strategy |
|---------------|------------------|----------------------------|
| `simplex` | single-strand ([Consensus Calling](consensus-calling.md)) | `adjacency` (or `identity`/`edit`) |
| `duplex` | double-strand ([Duplex Consensus Calling](duplex-consensus-calling.md)) | **`paired`** |
| `codec` | CODEC duplex from single read pairs | `adjacency` (or `identity`/`edit`) |

The algorithm is chosen by `--consensus`, **not** by the stage name. `--stop-after consensus` produces the same BAM as the matching `fgumi simplex`/`duplex`/`codec` command. **`--consensus duplex` requires `--group::strategy paired`** (so MIs carry the `/A`·`/B` strand suffixes duplex needs); `simplex` and `codec` require a non-paired strategy. Mismatching these is a hard error.

### Automatically inserted stages

You only name the *first* and *last* stages; `runall` fills in any stages in between that the chain requires. The one case that surprises people: the `align` and `zipper` stages emit reads in queryname order, but `group` needs template-coordinate order, so `runall` always inserts a `sort` step between them. For example, `--start-from extract --stop-after consensus` expands to extract → align → sort → group → consensus, and `--start-from zipper --stop-after group` expands to zipper → sort → group. You never write the inserted `sort` yourself, but you can still tune it with `--sort::*` flags.

## Per-stage options

Tune any stage with prefixed flags of the form `--<stage>::<flag>`. These mirror the standalone command's options exactly:

| Standalone | `runall` |
|------------|----------|
| `fgumi sort --max-memory 4G` | `--sort::max-memory 4G` |
| `fgumi group --strategy adjacency` | `--group::strategy adjacency` |
| `fgumi duplex --min-reads 1,1,0` | `--duplex::min-reads 1,1,0` |
| `fgumi extract --read-structures +T +T` | `--extract::read-structures +T +T` |

Cross-cutting options stay at the top level: `--threads`, `--compression-level`, `--max-memory`, `--ref`, `--rejects`, `--stats`, `--raw-tag` (default `RX`), `--assign-tag` (default `MI`), `--methylation-mode`. Run `fgumi runall --help` for the full, grouped list.

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
  --ref ref.fa --aligner::preset bwa --output consensus.bam \
  --group::strategy adjacency \
  --simplex::min-reads 1 \
  --threads 8
```

## `runall` vs. individual commands

| | `runall` | Individual commands |
|---|----------|---------------------|
| Intermediate BAMs | none (in-memory) | written to disk between stages |
| Output | identical to the command chain | the same bytes |
| Per-stage metrics/histograms | not emitted | available (e.g. `fgumi group --metrics`) |
| Inspect/reuse an intermediate | no | yes — the BAM is on disk |
| Best for | production runs, multi-stage chains | single stages, debugging, metrics, experimentation |

For the recommended end-to-end pipelines and parameter choices, see [Best Practices](best-practices.md).

## Caveats

- **No per-stage metrics.** Fused chains do not emit the metrics or histograms that standalone stages do (`fgumi group --metrics`, `simplex-metrics`, `duplex-metrics`). Run the standalone stage when you need them — see [Working with Metrics](working-with-metrics.md).
- **No `--stop-after align`.** Raw aligner output before the zipper-merge has lost every original tag, so the earliest stop after `--start-from align` is `zipper`.
- **Input state is not validated** against `--start-from` — see the warning under [Selecting the slice to run](#selecting-the-slice-to-run).
- **UMI rejects are discarded when `correct` is fused mid-chain.** To capture them, run `fgumi correct --rejects` as a separate step.
- **Strategy constraints.** `--consensus duplex` requires `--group::strategy paired`; `--consensus simplex` and `--consensus codec` require a non-paired strategy.
- **Methylation limits.** `--methylation-mode` is not supported with `--start-from align` or with `--consensus codec`; use it with `simplex` or `duplex`. See [Methylation](methylation.md).
- **`zipper` start needs a dictionary.** `--start-from zipper` requires `--ref` with a companion `.dict` file (create one with `samtools dict`).

## See Also

- [Getting Started](getting-started.md) — the step-by-step tutorial these stages come from
- [Best Practices](best-practices.md) — recommended end-to-end pipelines and parameters
- [Troubleshooting](troubleshooting.md) — fixes for common `runall` and sort-order errors
- [Glossary](glossary.md) — UMI, template-coordinate, strand labels, and BAM tags
