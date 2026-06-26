# Contributing to fgumi

For a full developer guide — project layout, the development lifecycle, and the release process — see [docs/DEVELOPING.md](docs/DEVELOPING.md). [CLAUDE.md](CLAUDE.md) documents the code architecture and key design patterns.

## Development Setup

### Prerequisites

- Rust 1.87.0 or newer (the project uses edition 2024)
- [cargo-nextest](https://nexte.st/) for running tests (`cargo install cargo-nextest --locked`)

The `cargo ci-*` commands below are project-specific aliases defined in [`.cargo/config.toml`](.cargo/config.toml); they wrap the formatting, lint, and test invocations CI runs.

### Install Git Hooks

We use pre-commit hooks to ensure code quality. Install them after cloning:

```bash
./scripts/install-hooks.sh
```

This installs hooks that run before each commit:
- `cargo ci-fmt` - Check code formatting
- `cargo ci-lint` - Run clippy lints

### Running Checks Manually

```bash
# Format check (fails if formatting differs)
cargo ci-fmt

# Lint check (fails on any warnings)
cargo ci-lint

# Run all tests
cargo ci-test
```

### Pre-Commit Hook Options

**Run tests in pre-commit hook:**
```bash
FGUMI_PRECOMMIT_TEST=1 git commit -m "message"
```

**Bypass hooks (use sparingly):**
```bash
git commit --no-verify -m "message"
```

## Code Style

- Run `cargo fmt` before committing
- Fix all clippy warnings
- Add backticks around identifiers in doc comments (e.g., `` `read_name` ``)

## Commit Messages

Use the [Conventional Commits](https://www.conventionalcommits.org/) format: `<type>[(scope)][!]: <description>`. Common types are `feat`, `fix`, `docs`, `refactor`, `perf`, `test`, `build`, `ci`, and `chore`.

```text
feat(group): add paired adjacency strategy
fix(codec): handle empty consensus families
```

## Branch Naming

Name branches `<issue-number>/<user>/<type>-<description>`, e.g. `42/jdidion/fix-fibonacci-calculation`. Never commit directly to `main` or `dev`; use a feature branch.

## Testing

All new features should include tests. Run the full test suite with:

```bash
cargo ci-test
```

## Pull Requests

1. Ensure all CI checks pass (`cargo ci-fmt`, `cargo ci-lint`, `cargo ci-test`)
2. Keep PRs focused and reasonably sized (250-1000 LOC ideal)
3. Include tests for new functionality
