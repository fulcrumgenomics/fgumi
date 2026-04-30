# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Refactor

- Extracted the sort engine into the new `fgumi-sort` crate. The main `fgumi` binary now consumes it as a workspace dependency; behavior is unchanged.
