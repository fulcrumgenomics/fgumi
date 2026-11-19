#!/usr/bin/env bash
#
# Bump fgumi's optional `bwa-mem3-rs` dependency to the latest crates.io
# release, and keep every pin that tracks the bwa-mem3 it vendors in step.
#
# A `bwa-mem3-rs` bump is not just a Cargo.toml edit. The crate's
# `bwa-mem3-sys` vendors a specific fg-labs/bwa-mem3 commit
# (`vendor/COMMIT`), and CI pins that same bwa-mem3 in two places:
#   - `e2e-parity` builds the reference CLI at the vendored sha (and fails
#     fast if the sha does not match the locked crate);
#   - `cpu-floor` installs the matching bwa-mem3 from bioconda.
# Comments in check.yml and the parity test also name the crate version, the
# vendored version and its short sha.
#
# Usage: scripts/ci/bwa-mem3-rs-bump.sh [--check] [--version X.Y.Z] [--summary FILE]
#   --check         Print current/latest and exit; change nothing.
#   --version       Bump to this version instead of the latest on crates.io.
#   --summary FILE  Write a Markdown summary of the bump (for a PR body).
#
# When $GITHUB_OUTPUT is set, writes `current`, `latest` and `needs_bump`.
#
# Requires: bash, curl, jq, mawk, python3, GNU-compatible `sort -V`; cargo and
# git only when bumping (not for --check).

set -euo pipefail

readonly CRATE="bwa-mem3-rs"
readonly SYS_CRATE="bwa-mem3-sys"
readonly INDEX_URL="https://index.crates.io/bw/a-/${CRATE}"
readonly UPSTREAM="https://github.com/fg-labs/bwa-mem3"
readonly CHECK_YML=".github/workflows/check.yml"
readonly PARITY_TEST="tests/align_inproc_parity.rs"
# Cargo.toml target table the dependency lives under (see its comment there).
readonly TARGET_CFG='cfg(any(target_arch = "x86_64", target_arch = "aarch64"))'

die() { echo "FATAL: $*" >&2; exit 1; }
is_semver() { [[ "$1" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; }
is_sha() { [[ "$1" =~ ^[0-9a-f]{40}$ ]]; }
# True when $1 sorts strictly before $2.
version_lt() {
  [ "$1" != "$2" ] && [ "$(printf '%s\n%s\n' "$1" "$2" | sort -V | head -1)" = "$1" ]
}
out() { if [ -n "${GITHUB_OUTPUT:-}" ]; then echo "$1=$2" >> "$GITHUB_OUTPUT"; fi; }

check_only=false
target=""
summary=""
while [ $# -gt 0 ]; do
  case "$1" in
    --check) check_only=true ;;
    --version) target="${2:?--version needs a value}"; shift ;;
    --summary) summary="${2:?--summary needs a value}"; shift ;;
    -h|--help) sed -n '2,/^$/s/^# \{0,1\}//p' "$0"; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
  shift
done

cd "$(dirname "$0")/../.."

requested=false
if [ -n "$target" ]; then
  is_semver "$target" || die "--version must be a stable X.Y.Z version, got '$target'"
  requested=true
fi

# The version requirement in Cargo.toml, e.g. `0.4.4` from
# `bwa-mem3-rs = { version = "0.4.4", optional = true }`.
current="$(sed -nE "s/^${CRATE} = \{ version = \"([^\"]+)\".*/\1/p" Cargo.toml)"
if [ -z "$current" ] || [ "$(wc -l <<< "$current")" -ne 1 ] || ! is_semver "$current"; then
  die "could not find exactly one '${CRATE} = { version = \"X.Y.Z\", ... }' line in Cargo.toml"
fi

if [ "$requested" = false ]; then
  # The sparse index has one JSON line per published version. Take the newest
  # stable (no pre-release suffix), non-yanked one.
  index="$(curl -fsS --retry 3 "$INDEX_URL")" || die "could not fetch ${INDEX_URL}"
  target="$(jq -r 'select(.yanked | not) | .vers' <<< "$index" \
    | { grep -E '^[0-9]+\.[0-9]+\.[0-9]+$' || true; } | sort -V | tail -1)"
  [ -n "$target" ] || die "no stable, non-yanked ${CRATE} release listed at ${INDEX_URL}"
fi

needs_bump=false
if version_lt "$current" "$target"; then
  needs_bump=true
fi
echo "${CRATE}: current=${current} target=${target} needs_bump=${needs_bump}"
if [ "$requested" = true ] && [ "$needs_bump" = false ]; then
  echo "Requested ${target} is not newer than the current ${current}; nothing to do (downgrades are not supported)."
fi
out current "$current"
out latest "$target"
out needs_bump "$needs_bump"

if [ "$check_only" = true ] || [ "$needs_bump" = false ]; then
  exit 0
fi

# The bwa-mem3 commit vendored by the locked bwa-mem3-sys (the same lookup the
# e2e-parity job does). `cargo metadata` downloads the crate if needed.
vendored_sha() {
  local meta manifest sha
  meta="$(cargo metadata --format-version 1 --locked --features aligner-bwa-mem3)" \
    || die "cargo metadata failed while resolving ${SYS_CRATE}"
  manifest="$(jq -r --arg n "$SYS_CRATE" \
    '[.packages[] | select(.name == $n) | .manifest_path][0] // empty' <<< "$meta")"
  [ -n "$manifest" ] || die "${SYS_CRATE} is not in the locked dependency graph"
  sha="$(tr -d '[:space:]' < "$(dirname "$manifest")/vendor/COMMIT")" \
    || die "could not read ${SYS_CRATE}'s vendor/COMMIT"
  is_sha "$sha" || die "${SYS_CRATE}'s vendor/COMMIT is not a 40-hex sha: '$sha'"
  echo "$sha"
}

# The upstream release version (`X.Y.Z`, no leading `v`) whose commit is $1,
# or empty when the sha is not a stable release tag. Annotated tags carry the
# commit on their peeled `^{}` line; pre-release tags (`v0.1.0-pre`) are
# ignored.
tag_for_sha() {
  local refs
  refs="$(git ls-remote --tags "$UPSTREAM")" || die "could not list the tags of ${UPSTREAM}"
  mawk -v sha="$1" '
    $1 == sha && $2 ~ /^refs\/tags\/v[0-9]+\.[0-9]+\.[0-9]+(\^\{\})?$/ {
      t = $2; sub(/^refs\/tags\/v/, "", t); sub(/\^\{\}$/, "", t); print t; exit
    }' <<< "$refs"
}

# Replace every match of the Python regex $2 in file $1 with the literal $3.
# Prints the number of replacements made (0 when nothing matched).
subst() {
  python3 - "$@" <<'PY'
import re, sys
path, pattern, new = sys.argv[1:4]
with open(path) as fh:
    text = fh.read()
text, n = re.subn(pattern, lambda _m: new, text)
if n:
    with open(path, "w") as fh:
        fh.write(text)
print(n)
PY
}

# Sum of `subst` over both pinned files.
subst_both() {
  local total=0 f
  for f in "$CHECK_YML" "$PARITY_TEST"; do
    total=$((total + $(subst "$f" "$1" "$2")))
  done
  echo "$total"
}

old_sha="$(vendored_sha)"
old_ver="$(tag_for_sha "$old_sha")"

echo "Bumping ${CRATE} ${current} -> ${target}"
# `cargo add` edits the existing entry in place (keeping `optional` and the
# surrounding comments). `--precise` then pins exactly the requested version;
# bwa-mem3-sys follows because bwa-mem3-rs requires `^` its own version.
cargo add "${CRATE}@${target}" --optional --target "$TARGET_CFG" -p fgumi
cargo update -p "$CRATE" --precise "$target"
grep -qE "^${CRATE} = \{ version = \"${target//./\\.}\"" Cargo.toml \
  || die "cargo add did not leave '${CRATE} = { version = \"${target}\" ...' in Cargo.toml"

new_sha="$(vendored_sha)"
new_ver="$(tag_for_sha "$new_sha")"
echo "vendored bwa-mem3: ${old_ver:-untagged} (${old_sha}) -> ${new_ver:-untagged} (${new_sha})"

updated=()
notes=()
cur_re="${current//./\\.}"

# Comments that name the crate version, e.g. "`bwa-mem3-rs` 0.4.4". The
# lookahead keeps 0.4.4 from matching inside 0.4.40.
n="$(subst_both "(?<=\`${CRATE}\` )${cur_re}(?![0-9.])" "$target")"
if [ "$n" -gt 0 ]; then
  updated+=("${n} \`${CRATE}\` version mention(s) in \`${CHECK_YML}\` / \`${PARITY_TEST}\`")
else
  notes+=("No \`${CRATE}\` ${current} mention found in \`${CHECK_YML}\` or \`${PARITY_TEST}\`; their version comments may already be stale.")
fi

# Check the actual reference pin even for crate-only releases, where the
# vendored commit is unchanged. A matching pin needs no edit; only advance
# the previously locked pin automatically, leaving unrelated pins for review.
reference_sha="$(sed -nE 's/^[[:space:]]*sha="([0-9a-f]{40})"[[:space:]]*$/\1/p' "$CHECK_YML")"
if [ "$reference_sha" = "$new_sha" ]; then
  : # Already in sync with the lockfile.
elif [ "$reference_sha" = "$old_sha" ]; then
  n="$(subst "$CHECK_YML" "(?m)(?<=sha=\")${old_sha}(?=\"[ \\t]*$)" "$new_sha")"
  updated+=("\`e2e-parity\` reference-CLI sha in \`${CHECK_YML}\`")
else
  notes+=("\`${CHECK_YML}\` does not have a unique \`e2e-parity\` sha matching the locked vendored commit. Set it to \`${new_sha}\` by hand.")
fi

if [ "$old_sha" != "$new_sha" ]; then
  # Short-sha mentions in comments (not followed by more hex, so the full sha
  # is never touched here).
  n="$(subst_both "\\b${old_sha:0:8}(?![0-9a-f])" "${new_sha:0:8}")"
  if [ "$n" -gt 0 ]; then
    updated+=("${n} short-sha mention(s) in comments")
  fi

  if [ -z "$new_ver" ] || [ -z "$old_ver" ]; then
    notes+=("The old (\`${old_sha:0:8}\`) or new (\`${new_sha:0:8}\`) vendored commit is not a stable fg-labs/bwa-mem3 release tag, so the \`vendored v<version>\` comments were left alone. Update them by hand.")
  else
    # Anchored on "vendored v" so action version comments (`# vX.Y.Z`) and
    # longer versions sharing the prefix are left alone.
    n="$(subst_both "(?<=vendored v)${old_ver//./\\.}(?![0-9.])" "$new_ver")"
    if [ "$n" -gt 0 ]; then
      updated+=("${n} \`vendored v<version>\` comment(s)")
    else
      notes+=("No \`vendored v${old_ver}\` comment found to update.")
    fi
  fi

  # Anything still naming the old commit was missed above.
  for f in "$CHECK_YML" "$PARITY_TEST"; do
    if grep -q "${old_sha:0:8}" "$f"; then
      notes+=("\`${f}\` still mentions the old vendored commit \`${old_sha:0:8}\`.")
    fi
  done
fi

# Check cpu-floor's actual pin on every bump, including crate-only releases.
# A matching pin needs no lookup or report; bioconda can lag upstream releases.
bioconda_pin="$(sed -nE 's/^[[:space:]]*create-args: bwa-mem3=([0-9]+\.[0-9]+\.[0-9]+) .*/\1/p' "$CHECK_YML")"
if [ -z "$new_ver" ]; then
  notes+=("The bioconda version cannot be resolved because the new vendored commit (\`${new_sha:0:8}\`) is not a stable fg-labs/bwa-mem3 release tag. The \`cpu-floor\` pin was left alone; check it by hand.")
elif [ -z "$bioconda_pin" ] || [ "$(wc -l <<< "$bioconda_pin")" -ne 1 ]; then
  notes+=("Could not find a unique \`create-args: bwa-mem3=X.Y.Z\` pin in \`${CHECK_YML}\`; set \`cpu-floor\` to \`bwa-mem3=${new_ver}\` by hand.")
elif [ "$bioconda_pin" != "$new_ver" ]; then
  if curl -fsS --retry 3 https://api.anaconda.org/package/bioconda/bwa-mem3 \
      | jq -e --arg v "$new_ver" '.versions | index($v)' > /dev/null; then
    n="$(subst "$CHECK_YML" "(?<=create-args: bwa-mem3=)${bioconda_pin//./\\.}(?= )" "$new_ver")"
    updated+=("\`cpu-floor\` bioconda pin (\`bwa-mem3=${new_ver}\`)")
  else
    notes+=("bioconda does not have \`bwa-mem3=${new_ver}\` yet (or the lookup failed), so \`cpu-floor\` still installs \`bwa-mem3=${bioconda_pin}\`, which does not match the vendored version. Bump it in \`${CHECK_YML}\` once the bioconda package is available.")
  fi
fi

for n in "${notes[@]+"${notes[@]}"}"; do echo "NOTE: $n" >&2; done

if [ -n "$summary" ]; then
  {
    echo "Bumps \`${CRATE}\` from \`${current}\` to \`${target}\` ([crates.io](https://crates.io/crates/${CRATE}/${target}), [releases](https://github.com/fg-labs/bwa-mem3-rs/releases))."
    echo
    if [ "$old_sha" = "$new_sha" ]; then
      echo "The vendored bwa-mem3 is unchanged (\`${old_ver:-${old_sha:0:8}}\`)."
    else
      echo "The vendored bwa-mem3 changes from \`${old_ver:-untagged}\` (\`${old_sha:0:8}\`) to \`${new_ver:-untagged}\` (\`${new_sha:0:8}\`) ([compare](${UPSTREAM}/compare/${old_sha}...${new_sha}))."
    fi
    if [ "${#updated[@]}" -gt 0 ]; then
      echo
      echo "Also updated:"
      for u in "${updated[@]}"; do echo "- $u"; done
    fi
    if [ "${#notes[@]}" -gt 0 ]; then
      echo
      echo "**Needs attention:**"
      for n in "${notes[@]}"; do echo "- $n"; done
    fi
    echo
    echo "Opened as a draft by the \`bwa-mem3-rs-bump\` workflow, which only edits pins. If the crate's API changed, CI fails until the adaptation is pushed onto this branch. Check \`e2e-parity\` (byte-parity against the vendored CLI) before marking it ready."
  } > "$summary"
fi
