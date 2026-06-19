# Design: `fgumi simulate aligner`

Date: 2026-06-18
Branch: `nh/feat-simulate-aligner-fr` (off `feat-runall`)

## Problem

fgumi-benchmarks' `runall_end_to_end` rule fakes the aligner by replaying a
pre-captured aligned BAM (so the benchmark measures fgumi's pipeline, not bwa's
compute). Today that replay is `workflow/scripts/replay-aligner.sh`, which reads
stdin and writes stdout in **lockstep on a single thread** — and that deadlocks
fgumi's AlignAndMerge (verified: hangs indefinitely; real bwa completes). Replace
it with a small, tested fgumi subcommand.

## Design

`fgumi simulate aligner --replay-bam <bam> <ref>` — a fake aligner subprocess
invoked via `--aligner::command`. Two independent threads:

- **Drainer:** read stdin to EOF, discard.
- **Emitter:** copy the replay BAM (header + all records, verbatim, in order) to
  stdout as BGZF, using fgumi's own BAM reader/writer (no samtools).

The threads never wait on each other, so it cannot deadlock. A real streaming
aligner reads and writes concurrently; two threads capture that. Backpressure
comes for free from the OS pipes plus fgumi's own in-flight gate.

This is deliberately the opposite of the shell script's one-thread lockstep —
that coupling was the entire bug.

## Details

- Lives at `src/lib/commands/simulate/aligner.rs`, behind the existing
  `simulate` feature gate, registered in `simulate/mod.rs` like the other
  `simulate` subcommands.
- Emit the header before/independently of draining (the emitter is its own
  thread, so this is automatic).
- Stream records **verbatim including secondary/supplementary** — fgumi's
  reader regroups by queryname; filtering would break its template counts.
- Accept and ignore the positional `{ref}` (and any `{threads}`) token —
  fgumi's command template requires `{ref}`.
- Treat `BrokenPipe` on stdout as a clean exit (fgumi closed the read end).

## Testing

Core is a function over generic `Read`/`Write`, tested in-process,
deterministic, sub-second:
- Liveness: a reader that withholds reads while stuffing stdin still completes
  (the case the shell script deadlocks on).
- Verbatim passthrough incl. non-primary records, header first.
- `BrokenPipe` → clean exit.

## Benchmark integration (fgumi-benchmarks)

Swap `runall_end_to_end`'s aligner command to `fgumi simulate aligner
--replay-bam {aligned_bam} {ref}`; delete `replay-aligner.sh`. Pre-captured
`shared/*.aligned.bam` fixtures stay (the verify rule needs the real bwa
alignments). Build the Docker image with `--features simulate`.

## Out of scope (add only if a real need appears)

Per-aligner `--model` presets, `-K`/chunk and thread knobs, emit-rate throttle.
The two-thread copy is all the benchmark needs; an optional `--rate`/`--chunk`
knob on the emitter can simulate a slower aligner later without changing the
structure.
