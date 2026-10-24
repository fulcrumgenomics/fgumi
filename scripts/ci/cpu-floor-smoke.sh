#!/usr/bin/env bash
# CPU-floor smoke test for an `aligner-bwa-mem3` build of fgumi on x86_64.
#
# AVX2 is the x86_64 floor for this build: the linked bwa-mem3 compiles its
# non-kernel code with -mavx2, and only the SIMD kernels are multi-tier,
# dispatched at runtime. #1063 was an AVX-512 kernel symbol leaking into the
# AVX2 path, which crashed every non-AVX-512 x86 CPU with SIGILL. This runs the
# binary under the qemu-user Haswell model (AVX2, no AVX-512) so that class of
# leak fails CI regardless of the runner's real CPU: `--version`, `extract` and
# an in-process `runall` alignment of error-bearing reads (reaching extension
# SW, `ksw_global2` and mate rescue) must all work there.
#
# Usage: cpu-floor-smoke.sh <fgumi-bin> <workdir>
# Requires `bwa-mem3`, `qemu-x86_64`, `samtools`, `python3` and `mawk` on PATH.
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "usage: $0 <fgumi-bin> <workdir>" >&2
  exit 2
fi
for tool in bwa-mem3 qemu-x86_64 samtools python3 mawk; do
  command -v "$tool" > /dev/null || { echo "FAIL: '$tool' is not on PATH" >&2; exit 1; }
done
fgumi="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
workdir="$2"
mkdir -p "$workdir"
cd "$workdir"

# The versioned name: the bare `Haswell-noTSX` is an alias "configured by
# machine type", which qemu-user has no machine for.
cpu="Haswell-v2" # Haswell without TSX: AVX2, no AVX-512
n_pairs=200

# `qemu-x86_64 -cpu help` prints the model list but exits non-zero, which
# `set -e` would turn into a silent abort; capture the list regardless.
cpu_models="$(qemu-x86_64 -cpu help || true)"
if ! grep -qw -- "$cpu" <<< "$cpu_models"; then
  echo "$cpu_models"
  echo "FAIL: qemu-x86_64 has no CPU model '$cpu'" >&2
  exit 1
fi

# A 20 kb random reference (the same generator as check.yml's e2e-parity job),
# indexed natively: indexing does not exercise the fgumi binary under test.
python3 -c "import random; random.seed(42); print('>ref'); print(''.join(random.choice('ACGT') for _ in range(20000)))" > ref.fa
bwa-mem3 index ref.fa > index.log 2>&1 || { echo "FAIL: bwa-mem3 index" >&2; cat index.log >&2; exit 1; }
samtools faidx ref.fa
samtools dict ref.fa -o ref.dict

# Read pairs sampled from the reference (R1: an 8 bp UMI then 92 bp forward;
# R2: the reverse complement 300 bp downstream), with errors chosen to reach
# the SIMD paths #1063 crashed in rather than bwa-mem3's exact-match shortcuts:
# two substitutions in every read and a mid-read 2 bp indel in half of them
# force extension SW and `ksw_global2`, and in every tenth pair R2 carries a
# mismatch every 8 bp (no 19-mer seed), so it can only align by mate rescue.
python3 - "$n_pairs" <<'EOF'
import random, sys
n = int(sys.argv[1])
seq = "".join(l.strip() for l in open("ref.fa") if not l.startswith(">"))
comp = str.maketrans("ACGT", "TGCA")
random.seed(7)

def substitute(read, positions):
    read = list(read)
    for i in positions:
        read[i] = random.choice([b for b in "ACGT" if b != read[i]])
    return "".join(read)

def mutate(read, length):
    read = substitute(read, random.sample(range(len(read)), 2))
    if random.random() < 0.5:
        i = random.randrange(30, 70)
        if random.random() < 0.5:
            read = read[:i] + read[i + 2:]
        else:
            read = read[:i] + "".join(random.choice("ACGT") for _ in range(2)) + read[i:]
    return read[:length]

with open("R1.fq", "w") as r1, open("R2.fq", "w") as r2:
    for i in range(n):
        p = random.randrange(0, len(seq) - 420)
        umi = "".join(random.choice("ACGT") for _ in range(8))
        b1 = umi + mutate(seq[p:p + 94], 92)
        fwd2 = seq[p + 300:p + 402]
        if i % 10 == 0:
            b2 = substitute(fwd2[:100], range(4, 100, 8)).translate(comp)[::-1]
        else:
            b2 = mutate(fwd2, 100).translate(comp)[::-1]
        r1.write(f"@r{i}\n{b1}\n+\n{'I' * len(b1)}\n")
        r2.write(f"@r{i}\n{b2}\n+\n{'I' * len(b2)}\n")
EOF

# Runs fgumi under the CPU model, keeping stderr in <label>.err; returns
# fgumi's exit status without tripping `set -e`.
run_on() {
  local label="$1"
  shift
  local status=0
  qemu-x86_64 -cpu "$cpu" "$fgumi" "$@" > "$label.out" 2> "$label.err" || status=$?
  echo "$status"
}

fail() {
  echo "FAIL: $1" >&2
  [ -n "${2:-}" ] && [ -f "$2" ] && tail -20 "$2" >&2
  exit 1
}

status=$(run_on version --version)
[ "$status" -eq 0 ] || fail "fgumi --version on $cpu exited $status" version.err

status=$(run_on extract extract --inputs R1.fq R2.fq \
  --read-structures 8M+T +T --sample s --library l --output u.bam)
[ "$status" -eq 0 ] || fail "fgumi extract on $cpu exited $status" extract.err

status=$(run_on inproc runall --start-from extract --stop-after zipper \
  --extract::inputs R1.fq R2.fq --extract::read-structures 8M+T +T \
  --extract::sample s --extract::library l \
  --ref ref.fa --aligner::preset bwa-mem3-inproc --output h.bam)
[ "$status" -eq 0 ] || fail "runall bwa-mem3-inproc on $cpu exited $status" inproc.err
[ -s h.bam ] || fail "runall bwa-mem3-inproc on $cpu wrote no h.bam" inproc.err
# Mapped primaries (unmapped templates are written too, so a record count alone
# cannot tell a kernel that aligns nothing from one that works) and proper
# pairs, which need the rescued mates.
mapped=$(samtools view -c -F 0x904 h.bam)
[ "$mapped" -eq $((2 * n_pairs)) ] \
  || fail "expected all $((2 * n_pairs)) primary records mapped in h.bam, got $mapped"
proper=$(samtools view -c -f 0x2 -F 0x904 h.bam)
[ "$proper" -ge $((2 * n_pairs * 95 / 100)) ] \
  || fail "expected >= 95% of $((2 * n_pairs)) records properly paired in h.bam, got $proper"
# The flag checks would also pass if every indel were soft-clipped away, so
# require the fixture's 2 bp indels to reach a primary CIGAR (a lower bound, not
# the exact count, which depends on the aligner).
indel=$(samtools view -F 0x904 h.bam | mawk '$6 ~ /(^|[MIDNSHP=X])2[ID]/ {n++} END {print n+0}')
[ "$indel" -ge 1 ] \
  || fail "expected at least one primary alignment with a 2 bp insertion or deletion CIGAR in h.bam, got $indel"

echo "cpu-floor smoke: OK ($cpu maps $mapped records in process, $proper properly paired, $indel with a 2 bp indel)"
