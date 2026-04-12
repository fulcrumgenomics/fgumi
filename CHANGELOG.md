# Changelog

All notable changes to this project will be documented in this file.

## [0.1.3] - 2026-04-11

### Bug Fixes

- Use monotonic counter for chunk file naming during consolidation ([#178](https://github.com/fulcrumgenomics/fgumi/pull/178))
- Exclude BAM bin field from core field comparison ([#208](https://github.com/fulcrumgenomics/fgumi/pull/208))
- Remove unused --sort-order flag from codec and simplex ([#209](https://github.com/fulcrumgenomics/fgumi/pull/209))
- Allow explicit true/false values on all boolean flags ([#210](https://github.com/fulcrumgenomics/fgumi/pull/210))
- Use semantic integer comparison for BAM tag values ([#214](https://github.com/fulcrumgenomics/fgumi/pull/214))
- Dedup --no-umi OOM on production WES data ([#231](https://github.com/fulcrumgenomics/fgumi/pull/231))
- Replace panic!() with graceful error handling in production code ([#223](https://github.com/fulcrumgenomics/fgumi/pull/223))
- Correct Zenodo DOI badge link in README ([#250](https://github.com/fulcrumgenomics/fgumi/pull/250))
- Abort pipeline on deadlock detection when recovery is disabled ([#252](https://github.com/fulcrumgenomics/fgumi/pull/252))
- Prevent pipeline deadlock when held items block reorder buffer progress ([#251](https://github.com/fulcrumgenomics/fgumi/pull/251))
- Prevent OOM and group flush timeout in unified pipeline ([#253](https://github.com/fulcrumgenomics/fgumi/pull/253))

### Documentation

- Improve cell barcode documentation for group and dedup commands ([#177](https://github.com/fulcrumgenomics/fgumi/pull/177))
- Update bam_codec comments to reflect noodles#364 closure ([#182](https://github.com/fulcrumgenomics/fgumi/pull/182))
- Update fulcrum genomics logo with light/dark theme support ([#194](https://github.com/fulcrumgenomics/fgumi/pull/194))
- Add trait overview table and extension guide for unified pipeline ([#222](https://github.com/fulcrumgenomics/fgumi/pull/222))
- Add mdBook documentation site with FG branding ([#243](https://github.com/fulcrumgenomics/fgumi/pull/243))
- Describe fgumi as a research preview instead of alpha ([#248](https://github.com/fulcrumgenomics/fgumi/pull/248))

### Features

- Add CB (cellular barcode) to template-coordinate sort key ([#160](https://github.com/fulcrumgenomics/fgumi/pull/160))
- Accept BAM input for mapped reads ([#181](https://github.com/fulcrumgenomics/fgumi/pull/181)) ([#183](https://github.com/fulcrumgenomics/fgumi/pull/183))
- Add SIMD-accelerated FASTQ parsing via fgumi-simd-fastq crate ([#180](https://github.com/fulcrumgenomics/fgumi/pull/180))
- Add fgumi merge command with loser tree ([#186](https://github.com/fulcrumgenomics/fgumi/pull/186))
- Add --includelist option to simulate fastq-reads ([#198](https://github.com/fulcrumgenomics/fgumi/pull/198))
- Add simplex-metrics command for simplex sequencing QC ([#195](https://github.com/fulcrumgenomics/fgumi/pull/195))
- Add position group size metrics and --metrics prefix output ([#232](https://github.com/fulcrumgenomics/fgumi/pull/232))
- Accept yes/no/y/n/t/f as boolean flag values ([#235](https://github.com/fulcrumgenomics/fgumi/pull/235))
- Add --max-memory=auto with system memory detection ([#236](https://github.com/fulcrumgenomics/fgumi/pull/236))

### Miscellaneous Tasks

- Pin GitHub Actions to full-length commit SHAs ([#215](https://github.com/fulcrumgenomics/fgumi/pull/215))
- Add cargo-audit security job to CI workflow ([#221](https://github.com/fulcrumgenomics/fgumi/pull/221))

### Performance

- Improve multi-threaded sort scaling ([#187](https://github.com/fulcrumgenomics/fgumi/pull/187))
- Replace LSD with MSD hybrid radix sort for template-coordinate ([#191](https://github.com/fulcrumgenomics/fgumi/pull/191))
- Replace RecordBuf with raw byte comparison for compare-bams ([#197](https://github.com/fulcrumgenomics/fgumi/pull/197))
- Reuse buffers in merge phase to reduce allocations ([#219](https://github.com/fulcrumgenomics/fgumi/pull/219))
- Optimize sort pipeline with LoserTree merge, EMBEDDED_IN_RECORD, and queryname specifiers ([#217](https://github.com/fulcrumgenomics/fgumi/pull/217))
- Eliminate RecordBuf re-encode bottleneck with raw-byte merge ([#228](https://github.com/fulcrumgenomics/fgumi/pull/228))
- Implement N+2 worker-pool model for parallel sort I/O ([#242](https://github.com/fulcrumgenomics/fgumi/pull/242))
- Redesign phase 2 work stealing and bound rayon to --threads ([#247](https://github.com/fulcrumgenomics/fgumi/pull/247))
- Add memory probes and force mimalloc arena collection ([#249](https://github.com/fulcrumgenomics/fgumi/pull/249))
- Bounded BGZF FASTQ pipeline with memory backpressure ([#244](https://github.com/fulcrumgenomics/fgumi/pull/244))

### Refactor

- Introduce ChunkNamer to centralize temp file naming ([#179](https://github.com/fulcrumgenomics/fgumi/pull/179))
- Remove vendored BAM decoder, use noodles Reader API ([#185](https://github.com/fulcrumgenomics/fgumi/pull/185))
- Use ClippingMode enum instead of String for --clipping-mode ([#218](https://github.com/fulcrumgenomics/fgumi/pull/218))
- Eliminate RecordBuf duplication by converting to raw-byte processing ([#229](https://github.com/fulcrumgenomics/fgumi/pull/229))
- Replace all unwrap() calls in src/commands/ with expect() ([#224](https://github.com/fulcrumgenomics/fgumi/pull/224))
- Replace all unwrap() calls in src/lib/ with expect() ([#226](https://github.com/fulcrumgenomics/fgumi/pull/226))
- Remove customizable tag options in favor of SAM spec standard tags ([#240](https://github.com/fulcrumgenomics/fgumi/pull/240))
- Remove dead rayon code paths superseded by N+2 pool ([#246](https://github.com/fulcrumgenomics/fgumi/pull/246))
- Introduce SamTag newtype for two-character BAM tag fields ([#241](https://github.com/fulcrumgenomics/fgumi/pull/241))

### Revert

- Restore --max-memory default to 768M ([#245](https://github.com/fulcrumgenomics/fgumi/pull/245))

### Testing

- Add end-to-end regression tests using simulate and compare ([#227](https://github.com/fulcrumgenomics/fgumi/pull/227))

<!-- generated by git-cliff -->

## [0.1.2] - 2026-03-06

### Bug Fixes

- Clamp compute_position fallback and validate sort-order ([#156](https://github.com/fulcrumgenomics/fgumi/pull/156))
- Expand ci-fmt and ci-lint to all workspace crates ([#162](https://github.com/fulcrumgenomics/fgumi/pull/162))

### Documentation

- Add noodles PR links to vendored/raw-record comments ([#159](https://github.com/fulcrumgenomics/fgumi/pull/159))

### Features

- Add --allow-unmapped flag for grouping unmapped reads by UMI ([#39](https://github.com/fulcrumgenomics/fgumi/pull/39))

### Miscellaneous Tasks

- Lockstep workspace versions and fix publish workflow ([#132](https://github.com/fulcrumgenomics/fgumi/pull/132))
- Bump noodles in the all-cargo-deps group ([#131](https://github.com/fulcrumgenomics/fgumi/pull/131))

### Refactor

- Improve code reuse, quality, and efficiency in fgumi-metrics ([#130](https://github.com/fulcrumgenomics/fgumi/pull/130))
- Remove dead code and unify duplicated logic in fgumi-consensus ([#139](https://github.com/fulcrumgenomics/fgumi/pull/139))
- Consolidate duplicated hash logic and remove dead code ([#140](https://github.com/fulcrumgenomics/fgumi/pull/140))
- Remove dead code and improve efficiency in consensus commands ([#154](https://github.com/fulcrumgenomics/fgumi/pull/154))
- Extract complement_base_preserve_case and remove unnecessary binding ([#142](https://github.com/fulcrumgenomics/fgumi/pull/142))
- Deduplicate core lib modules and remove redundant state ([#151](https://github.com/fulcrumgenomics/fgumi/pull/151))
- Extract generic verify helper and deduplicate test code ([#152](https://github.com/fulcrumgenomics/fgumi/pull/152))
- Deduplicate simulate utilities and clip header logic ([#155](https://github.com/fulcrumgenomics/fgumi/pull/155))
- Deduplicate utilities and hoist per-record hasher construction ([#149](https://github.com/fulcrumgenomics/fgumi/pull/149))
- Deduplicate scheduler constants and shared functions ([#150](https://github.com/fulcrumgenomics/fgumi/pull/150))
- Replace TagFamilySizeMetric with FamilySizeMetrics and remove dead parameter ([#153](https://github.com/fulcrumgenomics/fgumi/pull/153))
- Remove dead TemplateCoordinateKey and compare_template_coordinate_raw ([#161](https://github.com/fulcrumgenomics/fgumi/pull/161))

<!-- generated by git-cliff -->

## [0.1.1] - 2026-02-27

### Bug Fixes

- Write BGZF EOF explicitly in pipeline instead of EofWriter ([#126](https://github.com/fulcrumgenomics/fgumi/pull/126))
- Collapse input read group attributes into consensus output header ([#128](https://github.com/fulcrumgenomics/fgumi/pull/128))
- Replace no-op writer.finish(&header) with proper BAM writer finalization ([#127](https://github.com/fulcrumgenomics/fgumi/pull/127))

### Documentation

- Add Zenodo DOI badge to README ([#119](https://github.com/fulcrumgenomics/fgumi/pull/119))
- Add Bioconda badge to README ([#123](https://github.com/fulcrumgenomics/fgumi/pull/123))

### Features

- Make --ref optional to skip reference loading for unmapped reads ([#122](https://github.com/fulcrumgenomics/fgumi/pull/122))
- Add --no-umi mode for position-only grouping ([#56](https://github.com/fulcrumgenomics/fgumi/pull/56))

### Miscellaneous Tasks

- Remove temporary publish=false overrides from release-plz config ([#116](https://github.com/fulcrumgenomics/fgumi/pull/116))
- Bump the all-cargo-deps group with 6 updates ([#120](https://github.com/fulcrumgenomics/fgumi/pull/120))

### Refactor

- Move RawRecord/RawBamReader to fgumi-raw-bam crate ([#124](https://github.com/fulcrumgenomics/fgumi/pull/124))

### Testing

- Add integration tests for review and zipper commands ([#117](https://github.com/fulcrumgenomics/fgumi/pull/117))

<!-- generated by git-cliff -->
