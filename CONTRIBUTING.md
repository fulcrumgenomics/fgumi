# Contributing to fgumi

## Development Setup

### Prerequisites

- Rust (stable toolchain)
- cargo-nextest (for running tests)

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

## Testing

All new features should include tests. Run the full test suite with:

```bash
cargo ci-test
```

## Pull Requests

1. Ensure all CI checks pass (`cargo ci-fmt`, `cargo ci-lint`, `cargo ci-test`)
2. Keep PRs focused and reasonably sized (250-1000 LOC ideal)
3. Include tests for new functionality
