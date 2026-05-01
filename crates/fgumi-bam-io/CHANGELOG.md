# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Refactor

- Extracted the BAM-pipeline I/O layer into the new `fgumi-bam-io` crate. The main `fgumi` binary now consumes it as a workspace dependency; behavior is unchanged.
