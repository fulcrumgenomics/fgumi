#!/usr/bin/env bash
# Reproduce the agilent-hs2 benchmark reported in the performance chapter.
#
# Runs the full step-by-step pipeline (extract -> bwa -> zipper -> sort ->
# group -> simplex -> filter) and the equivalent `fgumi runall` invocation,
# wrapping each with /usr/bin/time so wall-clock can be compared stage-by-
# stage and end-to-end.
#
# Requirements:
#   * `bwa` and `samtools` on PATH
#   * `fgumi` binary (pass the path as the second argument)
#   * Agilent HS2 paired FASTQ files (R1/R2) — override via FASTQ_DIR
#   * hs38DH.fa reference with a `bwa index` and samtools .dict alongside —
#     override via REF
#
# Usage:
#   ./scripts/bench-repro.sh <work-dir> <fgumi-bin>
#
# Example:
#   ./scripts/bench-repro.sh /tmp/fgumi-bench ./target/release/fgumi

set -euo pipefail

WORK=${1:?work dir required}
FGUMI=${2:?fgumi binary path required}
FASTQ_DIR=${FASTQ_DIR:-/Volumes/scratch-00001/fgumi-benchmarks/data/raw/vendor}
REF=${REF:-/Volumes/scratch-00001/fgumi-benchmarks/data/resources/hs38DH.fa}
THREADS=${THREADS:-8}
BWA_CHUNK=${BWA_CHUNK:-150000000}

mkdir -p "$WORK"
cd "$WORK"

R1="$FASTQ_DIR/agilent-hs2_1.fastq.gz"
R2="$FASTQ_DIR/agilent-hs2_2.fastq.gz"

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "error: FASTQ files not found at $R1 / $R2" >&2
    echo "       override FASTQ_DIR to point at an agilent-hs2 FASTQ directory" >&2
    exit 1
fi
if [[ ! -f "$REF" ]]; then
    echo "error: reference FASTA not found at $REF" >&2
    echo "       override REF to point at hs38DH.fa (with bwa index + .dict)" >&2
    exit 1
fi

echo "=== v1 chain (step-by-step) ==="

echo "--- extract ---"
/usr/bin/time -p "$FGUMI" extract \
    --inputs "$R1" "$R2" \
    --output v1.unmapped.bam \
    --sample agilent-hs2 \
    --library lib1 \
    --read-structures 3M2S+T 3M2S+T

echo "--- fastq | bwa mem ---"
/usr/bin/time -p bash -c \
    "'$FGUMI' fastq --input v1.unmapped.bam | bwa mem -p -K $BWA_CHUNK -t $THREADS '$REF' - > v1.mapped.sam"

echo "--- zipper ---"
/usr/bin/time -p "$FGUMI" zipper \
    --input v1.mapped.sam \
    --unmapped v1.unmapped.bam \
    --reference "$REF" \
    --output v1.zippered.bam

echo "--- sort (template-coordinate) ---"
/usr/bin/time -p "$FGUMI" sort \
    --order template-coordinate \
    --input v1.zippered.bam \
    --output v1.sorted.bam

echo "--- group ---"
/usr/bin/time -p "$FGUMI" group \
    --input v1.sorted.bam \
    --output v1.group.bam \
    --strategy adjacency

echo "--- simplex ---"
/usr/bin/time -p "$FGUMI" simplex \
    --input v1.group.bam \
    --output v1.simplex.bam \
    --min-reads 1

echo "--- filter ---"
/usr/bin/time -p "$FGUMI" filter \
    --input v1.simplex.bam \
    --output v1.filter.bam \
    --min-reads 1 \
    --min-base-quality 2

echo ""
echo "=== v2 runall ==="
/usr/bin/time -p "$FGUMI" runall \
    --start-from extract \
    --input "$R1" "$R2" \
    --reference "$REF" \
    --output v2.final.bam \
    --threads "$THREADS" \
    --aligner::preset bwa-mem \
    --aligner::threads "$THREADS" \
    --aligner::chunk-size "$BWA_CHUNK" \
    --extract::sample agilent-hs2 \
    --extract::library lib1 \
    --extract::read-structures 3M2S+T 3M2S+T

echo ""
echo "Done. Compare v1.filter.bam vs v2.final.bam."
