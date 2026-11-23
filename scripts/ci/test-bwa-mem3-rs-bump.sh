#!/usr/bin/env bash
# Offline regression tests for the updater and its workflow's publishing steps.
set -euo pipefail

root="$(cd "$(dirname "$0")/../.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
readonly root tmp
real_git="$(command -v git)"
export REAL_GIT="$real_git"
export BRANCH_PREFIX=deps/bwa-mem3-rs- MARKER=bwa-mem3-rs-bump:pushed-sha
export BRANCH=deps/bwa-mem3-rs-0.5.0 CURRENT=0.4.4 LATEST=0.5.0 APP_SLUG=test

# Execute the actual workflow step, not a copy of the publishing logic.
step() {
  mawk -v name="$1" '
    /^      - name:/ { found = index($0, name) > 0; run = 0 }
    found && /^        run: \|/ { run = 1; next }
    run && /^          / { sub(/^          /, ""); print; next }
    run && NF { exit }
  ' "$root/.github/workflows/bwa-mem3-rs-bump.yml"
}
step 'Commit, push, and open the draft PR' > "$tmp/publish.sh"
step 'Supersede older bump PRs' > "$tmp/cleanup.sh"
step 'Check for an existing PR' > "$tmp/existing.sh"

mkdir "$tmp/bin"
cat > "$tmp/bin/gh" <<'SH'
#!/usr/bin/env bash
set -euo pipefail
case "$1 $2" in
  'api /users/test%5Bbot%5D') printf '123\n' ;;
  'pr list')
    printf 'list\n' >> "$CASE/list.log"
    query='.'
    while [ "$#" -gt 0 ]; do
      case "$1" in
        --jq|-q) query="$2"; shift ;;
      esac
      shift
    done
    jq -r "$query" "$CASE/prs.json" ;;
  'pr create'|'pr close'|'pr comment')
    printf '%s\n' "$*" >> "$CASE/gh.log"
    jq -cn --args '$ARGS.positional' -- "$@" >> "$CASE/gh-argv.jsonl"
    if [ "$2" = create ]; then
      while [ "$#" -gt 0 ]; do
        if [ "$1" = --body-file ]; then cp "$2" "$CASE/created-body.md"; break; fi
        shift
      done
    fi
    if [[ "$*" = *--delete-branch* ]]; then
      ref="$(jq -r '.[0].headRefName' "$CASE/prs.json")"
      "$REAL_GIT" push -q origin ":refs/heads/$ref"
    fi ;;
  *) printf 'Unexpected gh invocation: %s\n' "$*" >&2; exit 1 ;;
esac
SH
cat > "$tmp/bin/git" <<'SH'
#!/usr/bin/env bash
set -euo pipefail
# Simulate a push that lands after the PR snapshot / branch lookup.
if [ "$1" = push ] && [ -n "${RACE_SHA:-}" ]; then
  "$REAL_GIT" --git-dir="$CASE/remote.git" update-ref "refs/heads/$RACE_BRANCH" "$RACE_SHA"
fi
if [ "$1" = push ] && [ -n "${RACE_DELETE_BRANCH:-}" ]; then
  "$REAL_GIT" --git-dir="$CASE/remote.git" update-ref -d "refs/heads/$RACE_DELETE_BRANCH"
fi
if [ "$1" = ls-remote ] && [ "${LOOKUP_FAIL:-false}" = true ]; then
  printf 'Remote lookup unavailable\n' >&2
  exit 73
fi
exec "$REAL_GIT" "$@"
SH
chmod +x "$tmp/bin/gh" "$tmp/bin/git"
export PATH="$tmp/bin:$PATH"

setup_repo() {
  export CASE="$tmp/$1" RUNNER_TEMP="$tmp/$1" GITHUB_OUTPUT="$tmp/$1/output"
  unset RACE_SHA RACE_BRANCH RACE_DELETE_BRANCH LOOKUP_FAIL
  mkdir -p "$CASE/repo/.github/workflows" "$CASE/repo/tests"
  "$real_git" init -q --bare "$CASE/remote.git"
  "$real_git" init -q "$CASE/repo"
  cd "$CASE/repo"
  git config user.name Test
  git config user.email test@example.invalid
  git config commit.gpgsign false
  git config core.hooksPath /dev/null
  printf 'old\n' > Cargo.toml
  touch Cargo.lock .github/workflows/check.yml tests/align_inproc_parity.rs
  git add Cargo.toml Cargo.lock .github/workflows/check.yml tests/align_inproc_parity.rs
  git commit -qm base
  git remote add origin "$CASE/remote.git"
  git push -q origin HEAD:refs/heads/main
  base="$(git rev-parse HEAD)"
  printf 'adaptation\n' > Cargo.toml
  git commit -qam adaptation
  adaptation="$(git rev-parse HEAD)"
  git push -q origin "HEAD:refs/heads/$BRANCH"
  git switch -q --detach "$base"
  printf 'bumped\n' > Cargo.toml
  printf 'Summary\n' > "$CASE/summary.md"
  printf '[]\n' > "$CASE/prs.json"
  : > "$CASE/gh.log"
}

failures=0
scenarios=0
check() {
  if "$@"; then return; fi
  printf 'FAIL: %s: %s\n' "${test_name}" "$*" >&2
  failures=$((failures + 1))
}

# The existing-PR policy is separate from safe branch publication.
for spec in 'open OPEN false true true' 'merged MERGED true true true' \
    'closed_nightly CLOSED false true true' 'closed_retry CLOSED true true false' \
    'new none false true false' 'dry_run OPEN false false false' \
    'fork_open OPEN false true false'; do
  read -r test_name state RETRY PUBLISH expected_skip <<< "$spec"
  scenarios=$((scenarios + 1))
  export RETRY PUBLISH
  export CASE="$tmp/policy-$test_name" GITHUB_OUTPUT="$tmp/policy-$test_name/output"
  export GITHUB_STEP_SUMMARY="$tmp/policy-$test_name/summary"
  mkdir -p "$CASE"
  fork=false
  [ "$test_name" != fork_open ] || fork=true
  jq -n --arg state "$state" --argjson fork "$fork" \
    'if $state == "none" then [] else [{number:1,state:$state,isCrossRepository:$fork}] end' \
    > "$CASE/prs.json"
  bash "$tmp/existing.sh" > "$CASE/run.log" 2>&1
  check grep -Fxq "branch=$BRANCH" "$GITHUB_OUTPUT"
  check test "$(grep -c '^skip=true$' "$GITHUB_OUTPUT" || true)" -eq "$([ "$expected_skip" = true ] && printf 1 || printf 0)"
  if [ "$PUBLISH" = true ]; then
    check test "$(wc -l < "$CASE/list.log")" -eq 1
  else
    check test ! -e "$CASE/list.log"
  fi
done

for test_name in existing_branch concurrent_branch concurrent_fastforward new_branch; do
  scenarios=$((scenarios + 1))
  setup_repo "$test_name"
  if [ "$test_name" != existing_branch ]; then
    git push -q origin ":refs/heads/$BRANCH"
  fi
  if [ "$test_name" = concurrent_branch ]; then
    export RACE_SHA="$adaptation" RACE_BRANCH="$BRANCH"
  fi
  if [ "$test_name" = concurrent_fastforward ]; then
    export RACE_SHA="$base" RACE_BRANCH="$BRANCH"
  fi
  rc=0
  bash "$tmp/publish.sh" > "$CASE/run.log" 2>&1 || rc=$?
  remote="$(git ls-remote --heads origin "refs/heads/$BRANCH" | cut -f1)"
  if [ "$test_name" = new_branch ]; then
    check test "$rc" -eq 0
    check test "$remote" = "$(git rev-parse HEAD)"
    # shellcheck disable=SC2016
    check jq -se --arg branch "$BRANCH" '
      length == 1 and (.[0] as $a | $a[0:2] == ["pr","create"]
        and $a[2] == "--draft" and $a[3:5] == ["--base","main"]
        and $a[5:7] == ["--head",$branch])' "$CASE/gh-argv.jsonl"
    check grep -Fxq "<!-- ${MARKER}=${remote} -->" "$CASE/created-body.md"
    check grep -Fxq 'Summary' "$CASE/created-body.md"
    check test "$(cat "$GITHUB_OUTPUT")" = 'pushed=true'
  else
    check test "$rc" -ne 0
    expected_remote="$adaptation"
    [ "$test_name" != concurrent_fastforward ] || expected_remote="$base"
    check test "$remote" = "$expected_remote"
    check test ! -s "$CASE/gh.log"
  fi
done

for test_name in unchanged_old_branch concurrent_old_push concurrent_old_delete \
    unchanged_delete_denied cleanup_lookup_failed changed_marker missing_marker; do
  scenarios=$((scenarios + 1))
  setup_repo "$test_name"
  old_branch=deps/bwa-mem3-rs-0.4.5
  git push -q origin "$base:refs/heads/$old_branch"
  marker_sha="$base"
  case "$test_name" in
    changed_marker) marker_sha="$adaptation" ;;
    missing_marker) marker_sha= ;;
  esac
  jq -n --arg sha "$base" --arg ref "$old_branch" --arg marker "$MARKER" --arg marker_sha "$marker_sha" \
    '[{number: 1, isCrossRepository:false, headRefName: $ref, headRefOid: $sha, body: ("<!-- " + $marker + "=" + $marker_sha + " -->")}]' \
    > "$CASE/prs.json"
  if [ "$test_name" = concurrent_old_push ]; then
    # The head advances after gh pr list but before deletion.
    export RACE_SHA="$adaptation" RACE_BRANCH="$old_branch"
    git --git-dir="$CASE/remote.git" update-ref "refs/heads/$old_branch" "$adaptation"
  fi
  if [ "$test_name" = concurrent_old_delete ]; then
    export RACE_DELETE_BRANCH="$old_branch"
  fi
  case "$test_name" in
    unchanged_delete_denied|cleanup_lookup_failed)
      git --git-dir="$CASE/remote.git" config receive.denyDeletes true ;;
  esac
  if [ "$test_name" = cleanup_lookup_failed ]; then export LOOKUP_FAIL=true; fi
  rc=0
  bash "$tmp/cleanup.sh" > "$CASE/run.log" 2>&1 || rc=$?
  remote="$("$real_git" ls-remote --heads origin "refs/heads/$old_branch" | cut -f1)"
  if [ "$test_name" = unchanged_delete_denied ] || [ "$test_name" = cleanup_lookup_failed ]; then
    check test "$rc" -ne 0
    check test "$remote" = "$base"
    check test ! -s "$CASE/gh.log"
    if [ "$test_name" = cleanup_lookup_failed ]; then
      check test "$rc" -eq 73
      check grep -Ei 'error.*(read|look up|lookup|determine).*head' "$CASE/run.log"
    else
      check grep -Fq "$old_branch" "$CASE/run.log"
      check grep -Fq "$base" "$CASE/run.log"
      check grep -Ei 'error.*(delet|remov).*unchanged' "$CASE/run.log"
    fi
  elif [ "$test_name" != unchanged_old_branch ]; then
    check test "$rc" -eq 0
    expected_remote="$base"
    case "$test_name" in
      concurrent_old_push) expected_remote="$adaptation" ;;
      concurrent_old_delete) expected_remote= ;;
    esac
    check test "$remote" = "$expected_remote"
    check grep -q 'pr comment' "$CASE/gh.log"
    check test "$(grep -c 'pr close' "$CASE/gh.log" || true)" -eq 0
    case "$test_name" in
      concurrent_old_push) check grep -Fq "$adaptation" "$CASE/gh.log" ;;
      concurrent_old_delete) check grep -Fq absent "$CASE/gh.log" ;;
    esac
  else
    check test "$rc" -eq 0
    check test -z "$remote"
    check grep -q 'pr close' "$CASE/gh.log"
    check test "$(grep -c -- '--delete-branch' "$CASE/gh.log" || true)" -eq 0
  fi
done

# A mixed snapshot must touch only an older, same-repository bump.
for test_name in mixed_cleanup fork_cleanup; do
  scenarios=$((scenarios + 1))
  setup_repo "$test_name"
  old_branch=deps/bwa-mem3-rs-0.4.5
  newer_branch=deps/bwa-mem3-rs-0.6.0
  fork=false
  old_head="$base"
  if [ "$test_name" = fork_cleanup ]; then fork=true; old_head="$adaptation"; fi
  git push -q origin "$old_head:refs/heads/$old_branch" "$base:refs/heads/$newer_branch"
  # A fork can copy an edited upstream SHA and marker; neither proves ownership.
  jq -n --arg old "$old_branch" --arg current "$BRANCH" --arg newer "$newer_branch" \
    --arg sha "$old_head" --arg newer_sha "$base" --arg current_sha "$adaptation" --arg marker "$MARKER" --argjson fork "$fork" '
      [{number:1,headRefName:$old,headRefOid:$sha,isCrossRepository:$fork,
        body:("<!-- " + $marker + "=" + $sha + " -->")},
       {number:2,headRefName:$current,headRefOid:$current_sha,isCrossRepository:false,
        body:("<!-- " + $marker + "=" + $current_sha + " -->")},
       {number:3,headRefName:$newer,headRefOid:$newer_sha,isCrossRepository:false,
        body:("<!-- " + $marker + "=" + $newer_sha + " -->")}]' > "$CASE/prs.json"
  bash "$tmp/cleanup.sh" > "$CASE/run.log" 2>&1
  git ls-remote --heads origin > "$CASE/refs"
  if [ "$fork" = true ]; then
    printf '%s\trefs/heads/%s\n' "$old_head" "$old_branch" > "$CASE/expected-refs"
    check test ! -s "$CASE/gh.log"
  else
    : > "$CASE/expected-refs"
    # shellcheck disable=SC2016
    check jq -se --arg latest "$LATEST" --arg branch "$BRANCH" '
      . == [["pr","close","1","--comment",
        ("Superseded by the bump to `" + $latest + "` (branch `" + $branch + "`).")]]' \
      "$CASE/gh-argv.jsonl"
  fi
  printf '%s\trefs/heads/%s\n%s\trefs/heads/%s\n%s\trefs/heads/main\n' \
    "$adaptation" "$BRANCH" "$base" "$newer_branch" "$base" >> "$CASE/expected-refs"
  sort "$CASE/expected-refs" > "$CASE/expected-sorted"
  sort "$CASE/refs" > "$CASE/actual-sorted"
  check cmp "$CASE/expected-sorted" "$CASE/actual-sorted"
done

# Cargo and the remote release service are external boundaries. Fake only
# those; execute the complete bump script against generated input files.
mkdir "$tmp/pin-bin"
cat > "$tmp/pin-bin/cargo" <<'SH'
#!/usr/bin/env bash
set -euo pipefail
case "$1" in
  metadata)
    sha="$OLD_SHA"
    [ ! -f updated ] || sha="$NEW_SHA"
    printf '%s\n' "$sha" > vendor/COMMIT
    jq -n --arg path "$PWD/Cargo.toml" '{packages:[{name:"bwa-mem3-sys",manifest_path:$path}]}' ;;
  add) printf 'bwa-mem3-rs = { version = "0.5.0", optional = true }\n' > Cargo.toml ;;
  update) touch updated ;;
  *) exit 1 ;;
esac
SH
cat > "$tmp/pin-bin/git" <<'SH'
#!/usr/bin/env bash
set -euo pipefail
[ "$1 $2" = 'ls-remote --tags' ] || exit 1
if [ "${TAG_STYLE:-lightweight}" = annotated ]; then
  # Both direct and peeled prerelease refs match the new commit before the
  # stable tag; only the stable tag's peeled ref identifies its commit.
  printf '%s\t%s\n' \
    aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa refs/tags/v0.15.0 \
    "$OLD_SHA" 'refs/tags/v0.15.0^{}' \
    "$NEW_SHA" refs/tags/v0.16.0-rc.1 \
    cccccccccccccccccccccccccccccccccccccccc refs/tags/v0.16.0-rc.2 \
    "$NEW_SHA" 'refs/tags/v0.16.0-rc.2^{}' \
    bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb refs/tags/v0.16.0 \
    "$NEW_SHA" 'refs/tags/v0.16.0^{}'
  exit 0
fi
if [ -n "${OLD_TAG-0.15.0}" ]; then
  printf '%s\trefs/tags/v%s\n' "$OLD_SHA" "${OLD_TAG-0.15.0}"
fi
if [ -n "${NEW_TAG-0.16.0}" ]; then
  printf '%s\trefs/tags/v%s\n' "$NEW_SHA" "${NEW_TAG-0.16.0}"
fi
SH
cat > "$tmp/pin-bin/curl" <<'SH'
#!/usr/bin/env bash
set -euo pipefail
[ "${*: -1}" = https://api.anaconda.org/package/bioconda/bwa-mem3 ] || exit 1
# One file per (version, platform, label), as api.anaconda.org reports them;
# `.versions` lists a version once any platform has it.
files() {
  local label="$1" platforms="$2" sep='' v p
  printf '{"versions":["0.14.0","0.15.0","0.16.0"],"files":['
  for v in 0.14.0 0.15.0 0.16.0; do
    for p in $platforms; do
      [ "$v" = 0.14.0 ] && p=linux-64
      printf '%s{"version":"%s","attrs":{"subdir":"%s"},"labels":["%s"]}' "$sep" "$v" "$p" "$label"
      sep=,
      [ "$v" = 0.14.0 ] && break
    done
  done
  printf ']}\n'
}
case "${BIOCONDA_LOOKUP:-available}" in
  available) files main 'linux-64 linux-aarch64 osx-64 osx-arm64' ;;
  unavailable) printf '{"versions":["0.14.0"],"files":[{"version":"0.14.0","attrs":{"subdir":"linux-64"},"labels":["main"]}]}\n' ;;
  other_platform) files main 'linux-aarch64 osx-64 osx-arm64' ;;
  not_main) files dev 'linux-64 linux-aarch64 osx-64 osx-arm64' ;;
  failed) exit 22 ;;
  *) exit 1 ;;
esac
SH
chmod +x "$tmp/pin-bin/"*
export OLD_SHA=1111111111111111111111111111111111111111
for test_name in crate_only_stale crate_only_matching changed_old_pin changed_new_pin changed_stale changed_annotated_tag; do
  scenarios=$((scenarios + 1))
  export CASE="$tmp/$test_name" NEW_SHA="$OLD_SHA" TAG_STYLE=lightweight
  if [ "$test_name" = changed_annotated_tag ]; then TAG_STYLE=annotated; fi
  unset RACE_SHA RACE_BRANCH
  pin="$OLD_SHA"
  case "$test_name" in
    crate_only_stale) pin=3333333333333333333333333333333333333333 ;;
    changed_*) NEW_SHA=2222222222222222222222222222222222222222 ;;
  esac
  case "$test_name" in
    changed_new_pin) pin="$NEW_SHA" ;;
    changed_stale) pin=3333333333333333333333333333333333333333 ;;
  esac
  mkdir -p "$CASE/scripts/ci" "$CASE/.github/workflows" "$CASE/tests" "$CASE/vendor"
  cp "$root/scripts/ci/bwa-mem3-rs-bump.sh" "$CASE/scripts/ci/"
  printf 'bwa-mem3-rs = { version = "0.4.4", optional = true }\n' > "$CASE/Cargo.toml"
  # Backticks are literal Markdown in the generated fixture.
  # shellcheck disable=SC2016
  printf '# `bwa-mem3-rs` 0.4.4, vendored v0.15.0 (`11111111`)\nsha="%s"\ncreate-args: bwa-mem3=0.15.0 samtools\n' \
    "$pin" > "$CASE/.github/workflows/check.yml"
  # shellcheck disable=SC2016
  printf '// `bwa-mem3-rs` 0.4.4, vendored v0.15.0 (`11111111`)\n' \
    > "$CASE/tests/align_inproc_parity.rs"
  PATH="$tmp/pin-bin:$PATH" bash "$CASE/scripts/ci/bwa-mem3-rs-bump.sh" \
    --version 0.5.0 --summary "$CASE/summary.md" > "$CASE/run.log" 2>&1
  expected_version=0.15.0
  expected_short_sha=11111111
  case "$test_name" in
    changed_*) expected_version=0.16.0; expected_short_sha=22222222 ;;
  esac
  check grep -Fxq 'bwa-mem3-rs = { version = "0.5.0", optional = true }' "$CASE/Cargo.toml"
  check grep -Fxq "# \`bwa-mem3-rs\` 0.5.0, vendored v${expected_version} (\`${expected_short_sha}\`)" \
    "$CASE/.github/workflows/check.yml"
  check grep -Fxq "// \`bwa-mem3-rs\` 0.5.0, vendored v${expected_version} (\`${expected_short_sha}\`)" \
    "$CASE/tests/align_inproc_parity.rs"
  check grep -Fxq "create-args: bwa-mem3=${expected_version} samtools" "$CASE/.github/workflows/check.yml"
  if [ "$test_name" = changed_annotated_tag ]; then
    # shellcheck disable=SC2016
    check grep -Fq 'from `0.15.0` (`11111111`) to `0.16.0` (`22222222`)' "$CASE/summary.md"
    check test "$(grep -Ec 'untagged|rc\.|Needs attention' "$CASE/summary.md" || true)" -eq 0
  fi
  case "$test_name" in
    *stale)
      check grep -q 'Needs attention' "$CASE/summary.md"
      check grep -q "Set it to \`$NEW_SHA\`" "$CASE/summary.md"
      check grep -q "$pin" "$CASE/.github/workflows/check.yml" ;;
    *)
      check test "$(grep -c 'Needs attention' "$CASE/summary.md" || true)" -eq 0
      check grep -q "sha=\"$NEW_SHA\"" "$CASE/.github/workflows/check.yml" ;;
  esac
done
unset TAG_STYLE

# Compare the actual bioconda pin, independently of SHA changes and old tags.
# Same-SHA fixtures expose only v0.15.0, avoiding ambiguous tags on one commit.
while read -r test_name sha_change tags bioconda_pin BIOCONDA_LOOKUP expected_pin expected_comment attention; do
  scenarios=$((scenarios + 1))
  export CASE="$tmp/$test_name" NEW_SHA="$OLD_SHA" OLD_TAG=0.15.0 NEW_TAG='' BIOCONDA_LOOKUP
  new_version=0.15.0
  if [ "$sha_change" = changed ]; then
    NEW_SHA=2222222222222222222222222222222222222222
    NEW_TAG=0.16.0
    new_version=0.16.0
  fi
  case "$tags" in
    old_unknown) OLD_TAG= ;;
    new_unknown) NEW_TAG=; new_version= ;;
    both_unknown) OLD_TAG=; NEW_TAG=; new_version= ;;
  esac
  mkdir -p "$CASE/scripts/ci" "$CASE/.github/workflows" "$CASE/tests" "$CASE/vendor"
  cp "$root/scripts/ci/bwa-mem3-rs-bump.sh" "$CASE/scripts/ci/"
  printf 'bwa-mem3-rs = { version = "0.4.4", optional = true }\n' > "$CASE/Cargo.toml"
  # shellcheck disable=SC2016
  printf '# `bwa-mem3-rs` 0.4.4, vendored v0.15.0 (`11111111`)\nsha="%s"\n' \
    "$OLD_SHA" > "$CASE/.github/workflows/check.yml"
  if [ "$bioconda_pin" != missing ]; then
    printf 'create-args: bwa-mem3=%s samtools\n' "$bioconda_pin" >> "$CASE/.github/workflows/check.yml"
  fi
  # shellcheck disable=SC2016
  printf '// `bwa-mem3-rs` 0.4.4, vendored v0.15.0 (`11111111`)\n' \
    > "$CASE/tests/align_inproc_parity.rs"
  PATH="$tmp/pin-bin:$PATH" bash "$CASE/scripts/ci/bwa-mem3-rs-bump.sh" \
    --version 0.5.0 --summary "$CASE/summary.md" > "$CASE/run.log" 2>&1
  check grep -Fxq 'bwa-mem3-rs = { version = "0.5.0", optional = true }' "$CASE/Cargo.toml"
  check grep -Fxq "sha=\"$NEW_SHA\"" "$CASE/.github/workflows/check.yml"
  for file in .github/workflows/check.yml tests/align_inproc_parity.rs; do
    check grep -Fq "vendored v${expected_comment} (\`${NEW_SHA:0:8}\`)" "$CASE/$file"
  done
  if [ "$expected_pin" = missing ]; then
    check test "$(grep -c 'create-args:' "$CASE/.github/workflows/check.yml" || true)" -eq 0
  else
    check grep -Fxq "create-args: bwa-mem3=${expected_pin} samtools" "$CASE/.github/workflows/check.yml"
  fi
  if [ "$bioconda_pin" != "$expected_pin" ]; then
    check grep -Fq "bioconda pin (\`bwa-mem3=${expected_pin}\`)" "$CASE/summary.md"
  else
    check test "$(grep -c 'bioconda pin (' "$CASE/summary.md" || true)" -eq 0
  fi
  for report in "$CASE/summary.md" "$CASE/run.log"; do
    case "$attention" in
      none)
        if [ "$bioconda_pin" = "$expected_pin" ]; then
          check test "$(grep -Ec 'cpu-floor|bioconda' "$report" || true)" -eq 0
        else
          check test "$(grep -Ec 'bioconda.*(yet|left alone|by hand)|cpu-floor.*still' "$report" || true)" -eq 0
        fi ;;
      mismatch)
        check grep -Eq 'cpu-floor|bioconda' "$report"
        check grep -Fq "bwa-mem3=${bioconda_pin}" "$report"
        check grep -Fq "bwa-mem3=${new_version}" "$report" ;;
      unknown)
        check grep -Eiq 'bioconda.*(not.*stable|unknown|resolve)' "$report"
        check test "$(grep -c 'bioconda pin (' "$report" || true)" -eq 0 ;;
      missing)
        check grep -Fq 'Could not find' "$report"
        check grep -Fq "bwa-mem3=${new_version}" "$report" ;;
    esac
    if [ "$tags" = old_unknown ]; then
      check test "$(grep -Ec 'bioconda.*left alone|bioconda.*by hand' "$report" || true)" -eq 0
    fi
  done
done <<'CASES'
bioconda_same_stale_available same known 0.14.0 available 0.15.0 0.15.0 none
bioconda_same_stale_unavailable same known 0.14.0 unavailable 0.14.0 0.15.0 mismatch
bioconda_same_stale_failed same known 0.14.0 failed 0.14.0 0.15.0 mismatch
bioconda_same_stale_other_platform same known 0.14.0 other_platform 0.14.0 0.15.0 mismatch
bioconda_same_matching same known 0.15.0 available 0.15.0 0.15.0 none
bioconda_same_matching_unavailable same known 0.15.0 unavailable 0.15.0 0.15.0 none
bioconda_same_matching_failed same known 0.15.0 failed 0.15.0 0.15.0 none
bioconda_changed_stale_available changed known 0.14.0 available 0.16.0 0.16.0 none
bioconda_changed_stale_unavailable changed known 0.14.0 unavailable 0.14.0 0.16.0 mismatch
bioconda_changed_stale_failed changed known 0.14.0 failed 0.14.0 0.16.0 mismatch
bioconda_changed_stale_other_platform changed known 0.14.0 other_platform 0.14.0 0.16.0 mismatch
bioconda_changed_stale_not_main changed known 0.14.0 not_main 0.14.0 0.16.0 mismatch
bioconda_changed_matching changed known 0.16.0 available 0.16.0 0.16.0 none
bioconda_changed_matching_unavailable changed known 0.16.0 unavailable 0.16.0 0.16.0 none
bioconda_changed_matching_failed changed known 0.16.0 failed 0.16.0 0.16.0 none
bioconda_old_unknown_available changed old_unknown 0.14.0 available 0.16.0 0.15.0 none
bioconda_old_unknown_matching changed old_unknown 0.16.0 available 0.16.0 0.15.0 none
bioconda_old_unknown_unavailable changed old_unknown 0.14.0 unavailable 0.14.0 0.15.0 mismatch
bioconda_old_unknown_failed changed old_unknown 0.14.0 failed 0.14.0 0.15.0 mismatch
bioconda_new_unknown changed new_unknown 0.14.0 available 0.14.0 0.15.0 unknown
bioconda_both_unknown changed both_unknown 0.14.0 available 0.14.0 0.15.0 unknown
bioconda_same_unknown same both_unknown 0.14.0 available 0.14.0 0.15.0 unknown
bioconda_same_missing same known missing available missing 0.15.0 missing
bioconda_changed_missing changed known missing available missing 0.16.0 missing
bioconda_missing_unavailable changed known missing unavailable missing 0.16.0 missing
bioconda_missing_failed same known missing failed missing 0.15.0 missing
CASES
unset OLD_TAG NEW_TAG BIOCONDA_LOOKUP

# Check mode must select only stable, non-yanked releases and never edit inputs.
mkdir "$tmp/detect-bin"
cat > "$tmp/detect-bin/curl" <<'SH'
#!/usr/bin/env bash
set -euo pipefail
[ "$CURL_FAIL" = false ] || exit 22
cat "$CASE/index.jsonl"
SH
chmod +x "$tmp/detect-bin/curl"
for spec in 'detect_numeric auto false 0.4.10 true 0' \
    'detect_equal 0.4.4 false 0.4.4 false 0' \
    'detect_older 0.4.3 false 0.4.3 false 0' \
    'detect_failed auto true none none 1' \
    'detect_empty auto false none none 1'; do
  read -r test_name requested CURL_FAIL expected_latest expected_bump expected_rc <<< "$spec"
  scenarios=$((scenarios + 1))
  export CASE="$tmp/$test_name" CURL_FAIL GITHUB_OUTPUT="$tmp/$test_name/output"
  mkdir -p "$CASE/scripts/ci" "$CASE/.github/workflows" "$CASE/tests"
  cp "$root/scripts/ci/bwa-mem3-rs-bump.sh" "$CASE/scripts/ci/"
  printf 'bwa-mem3-rs = { version = "0.4.4", optional = true }\n' > "$CASE/Cargo.toml"
  printf 'locked\n' > "$CASE/Cargo.lock"
  printf 'workflow\n' > "$CASE/.github/workflows/check.yml"
  printf 'parity\n' > "$CASE/tests/align_inproc_parity.rs"
  : > "$CASE/index.jsonl"
  if [ "$test_name" != detect_empty ]; then
    printf '%s\n' \
      '{"vers":"0.4.9","yanked":false}' \
      '{"vers":"0.4.10","yanked":false}' \
      '{"vers":"0.5.0-rc.1","yanked":false}' \
      '{"vers":"0.5.0","yanked":true}' > "$CASE/index.jsonl"
  fi
  for file in Cargo.toml Cargo.lock .github/workflows/check.yml tests/align_inproc_parity.rs; do
    cp "$CASE/$file" "$CASE/$(basename "$file").before"
  done
  args=(--check)
  [ "$requested" = auto ] || args+=(--version "$requested")
  rc=0
  PATH="$tmp/detect-bin:$PATH" bash "$CASE/scripts/ci/bwa-mem3-rs-bump.sh" "${args[@]}" \
    > "$CASE/run.log" 2>&1 || rc=$?
  if [ "$expected_rc" -eq 0 ]; then
    check test "$rc" -eq 0
    printf 'current=0.4.4\nlatest=%s\nneeds_bump=%s\n' "$expected_latest" "$expected_bump" \
      > "$CASE/expected-output"
    check cmp "$CASE/expected-output" "$GITHUB_OUTPUT"
  else
    check test "$rc" -ne 0
    check test ! -s "$GITHUB_OUTPUT"
    check grep -q 'FATAL:' "$CASE/run.log"
  fi
  for file in Cargo.toml Cargo.lock .github/workflows/check.yml tests/align_inproc_parity.rs; do
    check cmp "$CASE/$(basename "$file").before" "$CASE/$file"
  done
done

if [ "$failures" -ne 0 ]; then
  printf '%s regression assertion(s) failed\n' "$failures" >&2
  exit 1
fi
printf 'All %s updater regression scenarios passed.\n' "$scenarios"
