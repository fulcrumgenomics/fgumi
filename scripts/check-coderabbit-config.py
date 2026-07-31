#!/usr/bin/env python3
"""Validate `.coderabbit.yaml` against the constraints that fail it silently.

CodeRabbit does not report configuration errors on the pull request. When the
file violates its schema, the bot discards the whole file and falls back to the
organization-UI defaults, so the repository silently loses every path
instruction and profile setting with no visible signal. That has happened here
before: the original `tone_instructions` was 293 characters against a 250
character cap, and the config was inert until it was noticed by hand.

This script turns each of those silent failures into a build failure:

* the file parses as YAML at all;
* the documented length caps are respected (a folded scalar is measured after
  folding, which is the length that actually counts);
* every `path_instructions` glob matches at least one tracked file, so an
  instruction cannot go on referring to a module that was renamed or removed;
* no backtick-quoted identifier was split across lines by YAML folding, which
  silently rewrites `foo_bar` into `foo_ bar`, or `Foo::bar` into `Foo:: bar`,
  in the text the reviewer reads.

Run it directly, with no arguments, from anywhere in the repository.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover - environment guard, not logic
    sys.exit("error: PyYAML is required; install it with `pip install pyyaml`")

# Caps from schema.v2.json. Keep in sync if the upstream schema changes; the
# JSON path for each is noted so it can be re-checked quickly.
TONE_INSTRUCTIONS_MAX = 250  # .properties.tone_instructions.maxLength
PATH_INSTRUCTIONS_MAX = 20000  # ...path_instructions.items.properties.instructions
LABELING_INSTRUCTIONS_MAX = 3000  # ...labeling_instructions.items.properties.instructions

# Matches a backtick-quoted span where a trailing `_` or `::` ran into a space,
# which is what YAML folding does to an identifier split across two lines. Both
# separators are needed: `_` catches the snake_case names, and `::` catches the
# fully-qualified paths (`MemoryTracker::is_below_drain_threshold` is 39
# characters, long enough that a rewrap breaks it, and `::` is where an editor
# breaks it). Word characters are required on both sides so that legitimate
# spaced syntax such as `case _ => ()` is not mistaken for a broken name.
BROKEN_IDENTIFIER = re.compile(r"`[^`]*\w(?:_|::) \w[^`]*`")


def repo_root() -> Path:
    """Return the repository root, so the script runs from any subdirectory."""
    result = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        capture_output=True,
        text=True,
        check=True,
    )
    return Path(result.stdout.strip())


def tracked_files(root: Path) -> list[str]:
    """Return every file tracked by git, as repository-relative paths."""
    result = subprocess.run(
        ["git", "ls-files"], cwd=root, capture_output=True, text=True, check=True
    )
    return [line for line in result.stdout.splitlines() if line]


def expand_braces(pattern: str) -> list[str]:
    """Expand `{a,b}` alternations into one pattern per alternative.

    CodeRabbit matches paths with minimatch, which supports brace expansion.
    The recursion resolves innermost-first, so it handles several groups in one
    pattern and nested groups alike. A nested group can yield the same
    alternative twice, which is harmless: `matches` only asks whether any
    alternative hits.
    """
    match = re.search(r"\{([^{}]*)\}", pattern)
    if not match:
        return [pattern]
    expanded = []
    for alternative in match.group(1).split(","):
        substituted = pattern[: match.start()] + alternative + pattern[match.end() :]
        expanded.extend(expand_braces(substituted))
    return expanded


def glob_to_regex(pattern: str) -> re.Pattern[str]:
    """Translate a minimatch-style glob into an anchored regex.

    `**` crosses directory separators (and may match zero of them), `*` and `?`
    do not.
    """
    out = ["^"]
    index = 0
    while index < len(pattern):
        char = pattern[index]
        if pattern.startswith("**/", index):
            out.append("(?:[^/]+/)*")
            index += 3
        elif pattern.startswith("**", index):
            out.append(".*")
            index += 2
        elif char == "*":
            out.append("[^/]*")
            index += 1
        elif char == "?":
            out.append("[^/]")
            index += 1
        else:
            out.append(re.escape(char))
            index += 1
    out.append("$")
    return re.compile("".join(out))


def matches(pattern: str, paths: list[str]) -> list[str]:
    """Return the tracked paths matching `pattern` under minimatch semantics."""
    regexes = [glob_to_regex(alternative) for alternative in expand_braces(pattern)]
    return [path for path in paths if any(regex.match(path) for regex in regexes)]


def check_length(label: str, text: str, cap: int, errors: list[str]) -> None:
    """Record an error when `text` exceeds `cap` characters."""
    if len(text) > cap:
        errors.append(
            f"{label} is {len(text)} characters, over the {cap}-character cap. "
            f"Exceeding it invalidates the entire file, not just this field."
        )


def main() -> int:
    root = repo_root()
    config_path = root / ".coderabbit.yaml"
    if not config_path.is_file():
        print(f"error: {config_path} not found", file=sys.stderr)
        return 1

    try:
        config = yaml.safe_load(config_path.read_text())
    except yaml.YAMLError as exc:
        print(f"error: .coderabbit.yaml is not valid YAML: {exc}", file=sys.stderr)
        return 1

    errors: list[str] = []
    reviews = config.get("reviews", {})

    check_length(
        "tone_instructions", config.get("tone_instructions", ""), TONE_INSTRUCTIONS_MAX, errors
    )

    paths = tracked_files(root)
    path_instructions = reviews.get("path_instructions", [])
    for entry in path_instructions:
        glob = entry["path"]
        check_length(
            f"path_instructions[{glob}]",
            entry["instructions"],
            PATH_INSTRUCTIONS_MAX,
            errors,
        )
        if not matches(glob, paths):
            errors.append(
                f"path_instructions glob {glob!r} matches no tracked file. "
                f"A stale glob silently stops aiming the reviewer at anything; "
                f"update it, or drop it until the code it names exists."
            )

    for entry in reviews.get("labeling_instructions", []):
        check_length(
            f"labeling_instructions[{entry['label']}]",
            entry["instructions"],
            LABELING_INSTRUCTIONS_MAX,
            errors,
        )

    # Folded scalars join their lines with a space, so an identifier wrapped
    # across two lines reaches the reviewer with a space inside it.
    for entry in path_instructions:
        for broken in BROKEN_IDENTIFIER.findall(entry["instructions"]):
            errors.append(
                f"path_instructions[{entry['path']}] contains {broken!r}, an identifier "
                f"split across lines by YAML folding. Rewrap so it stays on one line."
            )

    # `path_filters` is deliberately not checked for dead patterns: entries such
    # as `!**/target/**` are meant to match nothing tracked, and exist to exclude
    # build output that only appears in a working tree.

    if errors:
        print("`.coderabbit.yaml` validation failed:\n", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        print(
            "\nCodeRabbit reports none of these on the PR -- it discards the file "
            "and uses org defaults instead.",
            file=sys.stderr,
        )
        return 1

    print(
        f"check-coderabbit-config: OK "
        f"({len(path_instructions)} path instructions, all globs match tracked files; "
        f"tone_instructions {len(config.get('tone_instructions', ''))}/{TONE_INSTRUCTIONS_MAX})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
