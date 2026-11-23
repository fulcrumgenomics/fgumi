#!/usr/bin/env bash
# Exercise start-server.sh against a stub sccache: a server that starts must
# turn the wrapper on, and one that does not must leave it off without failing.
set -euo pipefail

script="$(dirname "$0")/start-server.sh"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
mkdir "$tmp/bin"
# pkill stub: never touch real processes.
printf '#!/bin/sh\nexit 0\n' > "$tmp/bin/pkill"
chmod +x "$tmp/bin/pkill"

check() { # <stub exit code> <expected GITHUB_ENV>
  printf '#!/bin/sh\nexit %s\n' "$1" > "$tmp/sccache"
  chmod +x "$tmp/sccache"
  : > "$tmp/env"
  PATH="$tmp/bin:$PATH" SCCACHE_PATH="$tmp/sccache" GITHUB_ENV="$tmp/env" bash "$script" > /dev/null
  if [[ "$(cat "$tmp/env")" != "$2" ]]; then
    echo "server exit $1: expected GITHUB_ENV [$2], got [$(cat "$tmp/env")]" >&2
    exit 1
  fi
}

check 0 $'SCCACHE_GHA_ENABLED=true\nRUSTC_WRAPPER=sccache'
check 2 'RUSTC_WRAPPER='
echo "start-server.sh: both outcomes behave"
