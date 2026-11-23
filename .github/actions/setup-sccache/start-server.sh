#!/usr/bin/env bash
# Start the sccache server against the GitHub Actions cache and route rustc
# through it, or, when it cannot start, leave rustc unwrapped so the build runs
# uncached. See action.yml for why.
set -uo pipefail

# Waits up to sccache's startup timeout (10s by default). With no idle timeout
# the server lives for the whole job, so a long stretch without compiles cannot
# make a later rustc start a new one that might not come up.
if SCCACHE_GHA_ENABLED=true SCCACHE_IDLE_TIMEOUT=0 "$SCCACHE_PATH" --start-server; then
  echo "SCCACHE_GHA_ENABLED=true" >> "$GITHUB_ENV"
  echo "RUSTC_WRAPPER=sccache" >> "$GITHUB_ENV"
  exit 0
fi

# Defensive: make sure no half-started server survives into later steps.
pkill -x sccache 2> /dev/null || true
# Empty, cargo uses no wrapper, even if a workflow `env:` set one.
echo "RUSTC_WRAPPER=" >> "$GITHUB_ENV"
echo "::warning title=sccache unavailable::The sccache server did not start (see the log above); building without a compilation cache. Re-running the job retries the cache."
