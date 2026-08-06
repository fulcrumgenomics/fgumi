#!/usr/bin/env bash
# Single source of truth for the crates.io publish order, plus the checks that
# keep it honest.
#
# This lives in a script rather than inline in `publish.yml` for one reason: the
# checks are only useful *before* a release. Inline, they ran exclusively in the
# publish job, which triggers on push to `main` — so a mis-ordered or incomplete
# list was discovered mid-publish, on main, after some crates had already gone
# out. `cargo publish` cannot be undone; a version can only be yanked. Every
# check here is a pure metadata/list check needing no registry access, so CI can
# run it on every PR and fail the PR instead of the release.
#
# Usage:
#   publish-crates.sh --list    print the crates in publish order, one per line
#   publish-crates.sh --check   validate the list against the workspace; exit 1 on error
set -euo pipefail

# Workspace crates in topological (dependency) order: leaf crates first, root
# crate last. Update when adding a workspace crate or changing inter-crate
# dependencies — `--check` fails the build if you forget.
CRATES=(
  fgumi-dna
  fgumi-bgzf
  fgumi-simd-fastq
  fgumi-tag
  fgumi-raw-bam
  fgumi-bam-io
  fgumi-umi
  fgumi-sam
  fgumi-metrics
  fgumi-sort
  fgumi-consensus
  fgumi
)

list_crates() {
  printf '%s\n' "${CRATES[@]}"
}

check_crates() {
  local meta ws_json workspace_crates listed_crates missing extra order_errors
  meta=$(cargo metadata --no-deps --format-version 1)

  # Crates marked `publish = false` (serialized as `publish: []`) are excluded —
  # e.g. internal xtask tooling.
  workspace_crates=$(echo "$meta" | jq -r '.packages[] | select(.publish != []) | .name' | sort)
  listed_crates=$(printf '%s\n' "${CRATES[@]}" | sort)

  missing=$(comm -23 <(echo "$workspace_crates") <(echo "$listed_crates"))
  if [ -n "$missing" ]; then
    echo "ERROR: workspace crates missing from the publish list:"
    echo "$missing"
    echo "Add them to CRATES in ${BASH_SOURCE[0]} in the correct dependency order."
    return 1
  fi

  extra=$(comm -13 <(echo "$workspace_crates") <(echo "$listed_crates"))
  if [ -n "$extra" ]; then
    echo "ERROR: the publish list contains entries not in the workspace:"
    echo "$extra"
    echo "Remove or rename stale entries before publishing."
    return 1
  fi

  # Every workspace dependency of a crate must appear EARLIER in the list.
  # `cargo publish` resolves each crate's dependencies against the crates.io
  # index, so a dependency published later does not yet exist when its dependent
  # is published — which fails mid-loop and leaves a partial release.
  ws_json=$(echo "$meta" | jq '[.packages[] | select(.publish != []) | .name]')
  index_of() {
    local target="$1" i
    for i in "${!CRATES[@]}"; do
      if [ "${CRATES[$i]}" = "$target" ]; then echo "$i"; return; fi
    done
    echo "-1"
  }
  order_errors=""
  for i in "${!CRATES[@]}"; do
    local crate deps j
    crate="${CRATES[$i]}"
    # Dev-dependencies with no version requirement (path-only, `req == "*"`) are
    # stripped from the manifest cargo uploads, so they impose no ordering
    # constraint. This is what makes the root crate's self dev-dependency —
    # `fgumi = { path = ".", features = ["test-utils"] }`, the only way to enable
    # `test-utils` for the integration test targets — publishable: without this
    # filter it reads as "fgumi depends on fgumi", which no ordering satisfies.
    # Dev-dependencies that DO carry a version survive packaging and must still
    # be published first.
    deps=$(echo "$meta" | jq -r --arg c "$crate" --argjson ws "$ws_json" '
      .packages[] | select(.name == $c) | .dependencies[]
      | select(.kind != "dev" or .req != "*") | .name
      | select(. as $d | $ws | index($d))' | sort -u)
    for dep in $deps; do
      j=$(index_of "$dep")
      if [ "$j" -ge "$i" ]; then
        order_errors="${order_errors}  ${crate} depends on ${dep}, which must be listed earlier\n"
      fi
    done
  done
  if [ -n "$order_errors" ]; then
    echo "ERROR: the publish list is not in valid dependency order:"
    printf "%b" "$order_errors"
    echo "Reorder CRATES so every dependency precedes its dependents."
    return 1
  fi

  echo "publish list OK: ${#CRATES[@]} crates, complete and in dependency order"
}

usage() {
  echo "usage: $(basename "$0") --list|--check"
}

case "${1:-}" in
  --list) list_crates ;;
  --check) check_crates ;;
  -h | --help) usage ;;
  *)
    usage >&2
    exit 2
    ;;
esac
