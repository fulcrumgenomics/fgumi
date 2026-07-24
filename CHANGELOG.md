# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

## [0.5.0] - 2026-07-24

### Bug Fixes

- Write tool/metric index as index.md to fix broken readthedocs link ([#456](https://github.com/fulcrumgenomics/fgumi/pull/456))
- Make sequential and parallel assigners agree on UMI case ([#455](https://github.com/fulcrumgenomics/fgumi/pull/455))
- Pair duplex strands sharing an unclipped 5' coordinate ([#485](https://github.com/fulcrumgenomics/fgumi/pull/485))
- Correct consensus reverse/revcomp tag sets to match fgbio ConsensusTags ([#488](https://github.com/fulcrumgenomics/fgumi/pull/488))
- Score duplicates by Picard SUM_OF_BASE_QUALITIES (DEDUP-02)
- Match fgbio FilterConsensusReads parity (FILT-01/02/03/04)
- Match fgbio strict read-name UMI extraction (EXT-01/03/04) ([#489](https://github.com/fulcrumgenomics/fgumi/pull/489))
- Unmap reads when clipping consumes the whole alignment ([#492](https://github.com/fulcrumgenomics/fgumi/pull/492))
- Require both aD and bD to treat a read as duplex ([#495](https://github.com/fulcrumgenomics/fgumi/pull/495))
- IUPAC-aware, case-preserving complement via shared 256-byte LUT ([#491](https://github.com/fulcrumgenomics/fgumi/pull/491))
- Match fgbio CorrectUmis metrics and log level (COR-01/02/03/04) ([#503](https://github.com/fulcrumgenomics/fgumi/pull/503))
- Count empty molecule-end halves, orient duplex UMIs by R1 strand, guard malformed UMIs (DXM-01/02/03, SIM-01) ([#506](https://github.com/fulcrumgenomics/fgumi/pull/506))
- Match fgbio ReviewConsensusVariants output and MI handling (REV-01/02/03/04/05) ([#500](https://github.com/fulcrumgenomics/fgumi/pull/500))
- Accept non-terminal + and reject over-long reads in read structures (R2-RS-01, R2-RS-02) ([#496](https://github.com/fulcrumgenomics/fgumi/pull/496))
- Emit review rows in fgbio coordinate + MI/read-number order ([#501](https://github.com/fulcrumgenomics/fgumi/pull/501))
- Clip chimeric templates, soft-only mate window, drop --sort-order (R2-CLIP-02/03/04) ([#502](https://github.com/fulcrumgenomics/fgumi/pull/502))
- Narrow read-name sync suffix strip to fgbio parity ([#512](https://github.com/fulcrumgenomics/fgumi/pull/512))
- Panic in aOrNotB/ln_a_minus_b when a < b (R2-NUM-01) ([#507](https://github.com/fulcrumgenomics/fgumi/pull/507))
- Sort-order + error-rate guards and duplex rejection rows (CONS-01, DUP-01, R2-MET-05, R2-UCC-01) ([#508](https://github.com/fulcrumgenomics/fgumi/pull/508))
- Align TSV metric output with fgbio Metric format ([#498](https://github.com/fulcrumgenomics/fgumi/pull/498)) ([#504](https://github.com/fulcrumgenomics/fgumi/pull/504))
- Write SS sub-sort tag as <sort-order>:<sub-sort> for fgbio/samtools parity ([#514](https://github.com/fulcrumgenomics/fgumi/pull/514))
- Require query-grouped input (FILT3-02, CLIP3-05) ([#517](https://github.com/fulcrumgenomics/fgumi/pull/517))
- Finish fgbio float format + reads-threshold guards (DXM3-04, DXM3-06, SIMM3-02) ([#520](https://github.com/fulcrumgenomics/fgumi/pull/520))
- Emit fgbio KV metrics format instead of the wide table ([#513](https://github.com/fulcrumgenomics/fgumi/pull/513))
- Reject UMIs of differing length in all assigners (GRP-01) ([#510](https://github.com/fulcrumgenomics/fgumi/pull/510))
- Match fgbio ReviewConsensusVariants parity (REV3-01..15) ([#522](https://github.com/fulcrumgenomics/fgumi/pull/522))
- Count masked bases only for retained primary reads (FILT-05) ([#523](https://github.com/fulcrumgenomics/fgumi/pull/523))
- Emit @HD SO:unsorted to match fgbio ([#526](https://github.com/fulcrumgenomics/fgumi/pull/526))
- [**breaking**] Apply fgbio pre-group filter to simplex/codec, add --allow-unmapped ([#509](https://github.com/fulcrumgenomics/fgumi/pull/509))
- Synthesize @HD when input lacks one ([#527](https://github.com/fulcrumgenomics/fgumi/pull/527))
- Validate input sort order to prevent silent corruption (MERGE3-01) ([#519](https://github.com/fulcrumgenomics/fgumi/pull/519))
- Place secondary/supplementary reads at exact template coordinate via tc ([#529](https://github.com/fulcrumgenomics/fgumi/pull/529))
- Write SS sub-sort tag via SortOrder accessors (R2-HDR-01) ([#531](https://github.com/fulcrumgenomics/fgumi/pull/531))
- Grouping-key correctness — strand tie-break, unsuffixed MI, library/cell partition, reject duplex-to-simplex (DXM3-02/03/05/07, SIMM3-01) ([#532](https://github.com/fulcrumgenomics/fgumi/pull/532))
- Add reverse-orientation edges to the parallel paired assigner (GRP3-01) ([#525](https://github.com/fulcrumgenomics/fgumi/pull/525))
- Upgrade clipping on all template reads, not just the primary pair ([#570](https://github.com/fulcrumgenomics/fgumi/pull/570))
- Strand-normalize overlap clipping, count existing fixed clipping, honor SoftWithMask upgrade (CLIP3-01/02/03) ([#528](https://github.com/fulcrumgenomics/fgumi/pull/528))
- Detect quality encoding from all input FASTQs, not just the first (EXT3-01) ([#550](https://github.com/fulcrumgenomics/fgumi/pull/550))
- Honest missing-quality default + singleton desync warning (FASTQ3-02/03) ([#560](https://github.com/fulcrumgenomics/fgumi/pull/560))
- [**breaking**] Clamp codec/duplex scalar depth+error tags to fgbio's Short ceiling ([#552](https://github.com/fulcrumgenomics/fgumi/pull/552))
- Accept CRLF and unterminated final records (W9b: EXT3-06/07) ([#554](https://github.com/fulcrumgenomics/fgumi/pull/554))
- Align reject-reason labels with fgbio wording (CODEC3-08) ([#566](https://github.com/fulcrumgenomics/fgumi/pull/566))
- Soft-only overlap-clip boundary, MC-independent codec clip, fgbio-faithful quality masking ([#533](https://github.com/fulcrumgenomics/fgumi/pull/533))
- Emit spec sub-sort spelling `lexicographical` in @HD SS (SORT3-10) ([#567](https://github.com/fulcrumgenomics/fgumi/pull/567))
- Template-coordinate order via canonical fgumi-sort + hermetic tests ([#576](https://github.com/fulcrumgenomics/fgumi/pull/576))
- Recompute the BAM bin field after raw POS/CIGAR mutations (clip, zipper) ([#591](https://github.com/fulcrumgenomics/fgumi/pull/591))
- Count final FASTQ batch before publishing read_done ([#598](https://github.com/fulcrumgenomics/fgumi/pull/598))
- Reject degenerate CLI values instead of silently misbehaving ([#601](https://github.com/fulcrumgenomics/fgumi/pull/601))
- Reserve Phase-1 input serials under the input lock ([#600](https://github.com/fulcrumgenomics/fgumi/pull/600))
- Emit faithful duplex consensus tags ([#603](https://github.com/fulcrumgenomics/fgumi/pull/603))
- Do not count a failed UMI-assignment group as accepted ([#610](https://github.com/fulcrumgenomics/fgumi/pull/610))
- Write the BAI sidecar to the samtools path, not an extension-replaced one ([#607](https://github.com/fulcrumgenomics/fgumi/pull/607))
- Emit unmapped consensus header on the single-threaded path ([#602](https://github.com/fulcrumgenomics/fgumi/pull/602))
- Accept EB/EiB memory sizes instead of rejecting them as scientific notation ([#615](https://github.com/fulcrumgenomics/fgumi/pull/615))
- Name the tag-skipping flag after the tag it actually writes ([#617](https://github.com/fulcrumgenomics/fgumi/pull/617))
- Give parallel UMI assigners a persistent molecule-id counter ([#630](https://github.com/fulcrumgenomics/fgumi/pull/630))
- Fail the run when UMI assignment fails ([#631](https://github.com/fulcrumgenomics/fgumi/pull/631))
- [**breaking**] Group --no-umi orientation-agnostically like Picard MarkDuplicates ([#649](https://github.com/fulcrumgenomics/fgumi/pull/649))
- Position-bin placed-but-unmapped reads in the BAI ([#651](https://github.com/fulcrumgenomics/fgumi/pull/651))
- Error on an over-long read name instead of panicking ([#629](https://github.com/fulcrumgenomics/fgumi/pull/629))

### Documentation

- Correct --write-index threading note ([#464](https://github.com/fulcrumgenomics/fgumi/pull/464))
- Document MI tag value divergence between fgumi group and fgbio GroupReadsByUmi ([#499](https://github.com/fulcrumgenomics/fgumi/pull/499))
- Clear fgbio-parity doc tail across commands (W11) ([#558](https://github.com/fulcrumgenomics/fgumi/pull/558))
- Document that dedup filters templates like fgbio GroupReadsByUmi (DEDUP3-01) ([#565](https://github.com/fulcrumgenomics/fgumi/pull/565))
- Fix stale doctest crate paths in consensus/umi and gate doctests in CI ([#573](https://github.com/fulcrumgenomics/fgumi/pull/573))
- Fix broken intra-doc links and gate rustdoc in CI ([#574](https://github.com/fulcrumgenomics/fgumi/pull/574))
- Document EXT3-03 encoding-failure parity with fgbio ([#587](https://github.com/fulcrumgenomics/fgumi/pull/587))
- Fix private intra-doc link breaking the docs build ([#597](https://github.com/fulcrumgenomics/fgumi/pull/597))
- Point Zenodo DOI badge at the concept DOI ([#596](https://github.com/fulcrumgenomics/fgumi/pull/596))
- Correct four statements that are false of the current behavior ([#619](https://github.com/fulcrumgenomics/fgumi/pull/619))
- Document the CODEC model and its quality-masking options ([#616](https://github.com/fulcrumgenomics/fgumi/pull/616))
- Document that --no-umi still splits by strand of origin ([#648](https://github.com/fulcrumgenomics/fgumi/pull/648))

### Features

- Strip old-style /1 and /2 read-number suffixes from QNAME ([#486](https://github.com/fulcrumgenomics/fgumi/pull/486))
- Body error injection, template-coordinate sort fix, and correctness guards ([#541](https://github.com/fulcrumgenomics/fgumi/pull/541))
- Add --include-unmapped to pass through no-mapped-read templates ([#589](https://github.com/fulcrumgenomics/fgumi/pull/589))
- [**breaking**] Harden fgumi compare into a sound, faithful fgbio-parity oracle ([#530](https://github.com/fulcrumgenomics/fgumi/pull/530))
- Add clip preset and lock review-.txt comparison (CMP3-02/03) ([#548](https://github.com/fulcrumgenomics/fgumi/pull/548))
- UMI-in-read-name output and real BGZF for .gz paths ([#605](https://github.com/fulcrumgenomics/fgumi/pull/605))
- Add per-phase --sort-threads / --merge-threads ([#608](https://github.com/fulcrumgenomics/fgumi/pull/608))
- Accept query-grouped input with --allow-unmapped ([#620](https://github.com/fulcrumgenomics/fgumi/pull/620))
- Add --max-temp-files to tune the spill-file consolidation limit ([#643](https://github.com/fulcrumgenomics/fgumi/pull/643))
- Add -t/--target {umi,barcode} to correct RX or BC ([#624](https://github.com/fulcrumgenomics/fgumi/pull/624))

### Miscellaneous Tasks

- Auto-review stacked PRs on non-default base branches ([#453](https://github.com/fulcrumgenomics/fgumi/pull/453))
- Fail the docs build on dead internal links via mdbook-linkcheck2 ([#457](https://github.com/fulcrumgenomics/fgumi/pull/457))
- Run the whole workspace in tests + fix stale fgumi-sam revcomp assertions ([#569](https://github.com/fulcrumgenomics/fgumi/pull/569))
- [**breaking**] Honest MSRV in lockstep with the toolchain + let-chain adoption ([#575](https://github.com/fulcrumgenomics/fgumi/pull/575))
- Gate sort correctness against samtools ([#577](https://github.com/fulcrumgenomics/fgumi/pull/577))
- Compile-gate all features + nightly stress-tests ([#578](https://github.com/fulcrumgenomics/fgumi/pull/578))
- Add a Miri gate for the raw-pointer comparator unsafe ([#581](https://github.com/fulcrumgenomics/fgumi/pull/581))

### Performance

- Stream consensus-read extraction via multi-interval query ([#584](https://github.com/fulcrumgenomics/fgumi/pull/584))
- Verify sort order in a single streaming pass ([#568](https://github.com/fulcrumgenomics/fgumi/pull/568))
- Stream records into the sort instead of an intermediate BAM ([#593](https://github.com/fulcrumgenomics/fgumi/pull/593))
- Bound sort-verify run comparison to O(order divergence) ([#594](https://github.com/fulcrumgenomics/fgumi/pull/594))
- Borrow reference bases in regenerate_alignment_tags ([#604](https://github.com/fulcrumgenomics/fgumi/pull/604))
- Narrow radix passes, reuse the stored queryname NUL, and borrow record bytes on ingest ([#606](https://github.com/fulcrumgenomics/fgumi/pull/606))
- Remove three per-block and per-read allocations ([#614](https://github.com/fulcrumgenomics/fgumi/pull/614))
- Derive the radix bound inside the first counting pass ([#622](https://github.com/fulcrumgenomics/fgumi/pull/622))
- Pool-integrate and optimize the coordinate --write-index merge ([#647](https://github.com/fulcrumgenomics/fgumi/pull/647))
- Validate mate CIGARs from the decoded key, not a second aux walk ([#639](https://github.com/fulcrumgenomics/fgumi/pull/639))

### Refactor

- Consolidate read-name suffix stripping and assert cross-parser parity ([#572](https://github.com/fulcrumgenomics/fgumi/pull/572))
- Unify per-template clip/mate-repair across threading paths ([#582](https://github.com/fulcrumgenomics/fgumi/pull/582))
- Guard FASTQ boundaries completion on empty input queue ([#599](https://github.com/fulcrumgenomics/fgumi/pull/599))

### Testing

- Read stats as fgbio KV metrics in pre-group-filter test ([#585](https://github.com/fulcrumgenomics/fgumi/pull/585))
- Gate command output against its @HD sort-order claim ([#580](https://github.com/fulcrumgenomics/fgumi/pull/580))
- Skip float-parseable garbage in AF-fallback proptest ([#592](https://github.com/fulcrumgenomics/fgumi/pull/592))

<!-- generated by git-cliff -->

## [0.4.0] - 2026-06-20

### Miscellaneous Tasks

- Ignore RUSTSEC-2024-0436 (paste unmaintained, tracked in #26) ([#444](https://github.com/fulcrumgenomics/fgumi/pull/444))

### Refactor

- [**breaking**] Remove deprecated --queue-memory* aliases ([#445](https://github.com/fulcrumgenomics/fgumi/pull/445))

<!-- generated by git-cliff -->

## [0.3.1] - 2026-06-17

### Bug Fixes

- Read process stdin directly in open_bgzf_reader ([#418](https://github.com/fulcrumgenomics/fgumi/pull/418))
- Emit per-base arrays and float error rates in consensus-reads ([#423](https://github.com/fulcrumgenomics/fgumi/pull/423))

### Features

- Support reading BAM input from stdin ([#415](https://github.com/fulcrumgenomics/fgumi/pull/415))
- Support reading BAM input from stdin for correct, codec, duplex, filter ([#419](https://github.com/fulcrumgenomics/fgumi/pull/419))

### Miscellaneous Tasks

- Add CodeRabbit config (assertive profile + path instructions) ([#392](https://github.com/fulcrumgenomics/fgumi/pull/392))

### Refactor

- Route group/dedup/clip stdin gate through BamIoOptions::validate ([#428](https://github.com/fulcrumgenomics/fgumi/pull/428))

### Styling

- Clear clippy --all-features --all-targets warnings ([#427](https://github.com/fulcrumgenomics/fgumi/pull/427))

<!-- generated by git-cliff -->

## [0.3.0] - 2026-06-10

### Bug Fixes

- Make root crate build cleanly with `--no-default-features` ([#314](https://github.com/fulcrumgenomics/fgumi/pull/314))
- Deterministic MoleculeId numbering via opt-in MI Assign stage ([#319](https://github.com/fulcrumgenomics/fgumi/pull/319))
- Address remaining review feedback from PR #323 ([#325](https://github.com/fulcrumgenomics/fgumi/pull/325))
- [**breaking**] Require SS:template-coordinate ([#337](https://github.com/fulcrumgenomics/fgumi/pull/337))
- Drop dense (max_ab × max_ba) grid that OOMs on cfDNA hot spots ([#359](https://github.com/fulcrumgenomics/fgumi/pull/359))
- Accept --compression-level 0 for uncompressed BAM ([#360](https://github.com/fulcrumgenomics/fgumi/pull/360)) ([#361](https://github.com/fulcrumgenomics/fgumi/pull/361))
- Use 128-bit composite hash for read-name keys ([#369](https://github.com/fulcrumgenomics/fgumi/pull/369))

### Documentation

- Use absolute URLs for logo images ([#324](https://github.com/fulcrumgenomics/fgumi/pull/324))
- Recommend uncompressed BAM input for streaming pipelines ([#362](https://github.com/fulcrumgenomics/fgumi/pull/362))

### Features

- Zstd as the default temp-spill codec ([#341](https://github.com/fulcrumgenomics/fgumi/pull/341))
- --parallel-group-min-templates with per-strategy auto thresholds ([#371](https://github.com/fulcrumgenomics/fgumi/pull/371))
- Variable-width template-coordinate sort key, zero-copy merge, and slimmer arena ([#375](https://github.com/fulcrumgenomics/fgumi/pull/375))
- Add host-aware --max-memory to pipeline commands ([#381](https://github.com/fulcrumgenomics/fgumi/pull/381))

### Miscellaneous Tasks

- Add publish dry-run to catch release-blocking bugs pre-merge ([#312](https://github.com/fulcrumgenomics/fgumi/pull/312))
- Add codecov.yml with project coverage threshold ([#374](https://github.com/fulcrumgenomics/fgumi/pull/374))
- Bump codecov-action to v6.0.2 to fix upload GPG verification ([#383](https://github.com/fulcrumgenomics/fgumi/pull/383))

### Performance

- Cache UMI tag position in DecodedRecord to avoid duplicate aux scan ([#343](https://github.com/fulcrumgenomics/fgumi/pull/343))
- Single-pass aux scan captures UMI position ([#334](https://github.com/fulcrumgenomics/fgumi/pull/334)) ([#351](https://github.com/fulcrumgenomics/fgumi/pull/351))
- Bypass libdeflater for deflate stored blocks ([#366](https://github.com/fulcrumgenomics/fgumi/pull/366))
- Pin flate2 zlib-rs backend to restore fast FASTQ decompression ([#387](https://github.com/fulcrumgenomics/fgumi/pull/387))

### Refactor

- Validate and consolidate SAM tag identifiers via SamTag ([#316](https://github.com/fulcrumgenomics/fgumi/pull/316))
- Extract fgumi-bam-io and fgumi-sort workspace crates ([#323](https://github.com/fulcrumgenomics/fgumi/pull/323))
- Migrate hand-rolled rejects writers to first-class secondary_output ([#332](https://github.com/fulcrumgenomics/fgumi/pull/332))
- Typed CodecConsensusError replaces substring matching ([#342](https://github.com/fulcrumgenomics/fgumi/pull/342))
- Replace dead Binomial::new Err arm with expect ([#363](https://github.com/fulcrumgenomics/fgumi/pull/363))
- Delegate ConsensusOutput heap-size estimate to owning type ([#368](https://github.com/fulcrumgenomics/fgumi/pull/368))
- Demote algorithm-internal info! logs to debug ([#367](https://github.com/fulcrumgenomics/fgumi/pull/367))

### Testing

- In-process Command::execute for Phase 1 commands ([#352](https://github.com/fulcrumgenomics/fgumi/pull/352)) ([#353](https://github.com/fulcrumgenomics/fgumi/pull/353))
- In-process Command::execute for Phase 2 commands ([#352](https://github.com/fulcrumgenomics/fgumi/pull/352)) ([#354](https://github.com/fulcrumgenomics/fgumi/pull/354))
- In-process Command::execute for Phase 4 commands ([#352](https://github.com/fulcrumgenomics/fgumi/pull/352)) ([#356](https://github.com/fulcrumgenomics/fgumi/pull/356))
- In-process follow-ups from Opus/Sonnet review of #352 ([#357](https://github.com/fulcrumgenomics/fgumi/pull/357))

<!-- generated by git-cliff -->

### Documentation

- Recommend uncompressed BAM (`bwa-mem3 mem --bam=0`) as the preferred `fgumi zipper` input format for streaming pipelines, with SAM kept as a fallback for aligners that can't emit BAM. Replaced the misleading "BAM input is discouraged" `warn!` on the stdin-BAM path with a neutral `info!`, and dropped the equivalent warning on the file-BAM path. Updated `zipper` CLI help/long_about and the getting-started, best-practices, performance-tuning, and migration-from-fgbio guides accordingly.

### Bug Fixes

- [**breaking**] `fgumi group`, `fgumi dedup`, and `fgumi downsample` now strictly require template-coordinate sorted input — the header must advertise `SO:unsorted`, `GO:query`, **and** `SS:template-coordinate`. The previous implementation accepted any `SO:unsorted GO:query` header even without `SS:template-coordinate`, which caused `fgumi extract` output (FASTQ-order, no `SS`) and other queryname-grouped BAMs to be silently treated as template-coordinate sorted. Because the streaming `RecordPositionGrouper` requires records sharing a position key to be consecutive, this could split a single true molecule across multiple position groups and assign distinct `MI` tags to reads that should share one. The `--allow-unmapped` queryname-sort shortcut has also been removed for the same reason. To migrate: insert `fgumi sort --order template-coordinate` between extract (or any non-TC source) and `fgumi group`/`dedup`/`downsample`.

### Refactor

- Extracted the sort engine into the new `fgumi-sort` crate and the BAM-pipeline I/O layer into the new `fgumi-bam-io` crate. The main `fgumi` binary now consumes both as workspace dependencies; behavior is unchanged.

## [0.2.0] - 2026-04-22

### Bug Fixes

- Add fgumi-simd-fastq to publish list and verify completeness ([#257](https://github.com/fulcrumgenomics/fgumi/pull/257))
- Take last `:`-separated field as UMI from read name ([#264](https://github.com/fulcrumgenomics/fgumi/pull/264))
- [**breaking**] Rename sort-key tag from pa to tc to avoid bwa-mem pa:f collision ([#270](https://github.com/fulcrumgenomics/fgumi/pull/270))
- Emit consecutive MI integers 0..N-1 ([#273](https://github.com/fulcrumgenomics/fgumi/pull/273))
- Cap metric memory with per-thread accumulators ([#285](https://github.com/fulcrumgenomics/fgumi/pull/285)) ([#287](https://github.com/fulcrumgenomics/fgumi/pull/287))
- Suppress deadlock false positive when all queues empty ([#297](https://github.com/fulcrumgenomics/fgumi/pull/297))
- Stream rejects to disk in simplex/duplex/codec/correct ([#293](https://github.com/fulcrumgenomics/fgumi/pull/293))
- Apply --unmapped-fraction to mapped-reads output ([#304](https://github.com/fulcrumgenomics/fgumi/pull/304))
- Use htsjdk-compatible Murmur3 for downsample selection ([#306](https://github.com/fulcrumgenomics/fgumi/pull/306))
- Preserve /A /B strand suffix in paired-UMI MI keys ([#308](https://github.com/fulcrumgenomics/fgumi/pull/308))

### Documentation

- Add EM-Seq user guide ([#171](https://github.com/fulcrumgenomics/fgumi/pull/171))
- Add TAPs guide and update EM-Seq guide for --methylation-mode ([#174](https://github.com/fulcrumgenomics/fgumi/pull/174))
- Correct insert complexity and add micro-bench ([#265](https://github.com/fulcrumgenomics/fgumi/pull/265))
- Fix broken doctests in commands modules ([#300](https://github.com/fulcrumgenomics/fgumi/pull/300))

### Features

- Add EM-Seq methylation-aware consensus calling ([#168](https://github.com/fulcrumgenomics/fgumi/pull/168))
- Add EM-Seq methylation filters and performance improvements ([#169](https://github.com/fulcrumgenomics/fgumi/pull/169))
- Add --restore-unconverted-bases flag for EM-seq ([#170](https://github.com/fulcrumgenomics/fgumi/pull/170))
- Add --methylation-mode for EM-Seq/TAPs consensus ([#172](https://github.com/fulcrumgenomics/fgumi/pull/172))
- Add --methylation-mode to filter for conversion fraction ([#173](https://github.com/fulcrumgenomics/fgumi/pull/173))
- Add --methylation-mode to simulate subcommands ([#263](https://github.com/fulcrumgenomics/fgumi/pull/263))
- Auto-detect BAM input on stdin via BGZF magic ([#267](https://github.com/fulcrumgenomics/fgumi/pull/267))
- Replace --mode with orthogonal flags + add --command preset ([#276](https://github.com/fulcrumgenomics/fgumi/pull/276))
- Add --async-reader flag for prefetch I/O on input BAM ([#261](https://github.com/fulcrumgenomics/fgumi/pull/261))
- Distribute spill files across multiple temp dirs ([#284](https://github.com/fulcrumgenomics/fgumi/pull/284))
- Additive API surface ([#288](https://github.com/fulcrumgenomics/fgumi/pull/288))
- Restore --command preset for per-stage compare defaults ([#307](https://github.com/fulcrumgenomics/fgumi/pull/307))

### Performance

- Speed up --restore-unconverted-bases hot path ([#271](https://github.com/fulcrumgenomics/fgumi/pull/271))
- Async prefetch reader + POSIX_FADV_SEQUENTIAL for BAM/FASTQ inputs ([#258](https://github.com/fulcrumgenomics/fgumi/pull/258))
- Restore unconverted bases on raw BAM bytes ([#274](https://github.com/fulcrumgenomics/fgumi/pull/274))
- Follow-up speedups + SIMD nibble2base ([#282](https://github.com/fulcrumgenomics/fgumi/pull/282))
- Migrate RecordBuf hot paths to raw-byte accessors ([#305](https://github.com/fulcrumgenomics/fgumi/pull/305))

### Refactor

- Migrate fgumi to new raw-bam view/editor API ([#279](https://github.com/fulcrumgenomics/fgumi/pull/279))
- [**breaking**] Narrow re-exports + partial free-fn demotion ([#280](https://github.com/fulcrumgenomics/fgumi/pull/280))
- Use RawRecord methods in restore_unconverted_bases_in_raw_record ([#286](https://github.com/fulcrumgenomics/fgumi/pull/286))
- PerThreadAccumulator for all per-batch metric collection ([#290](https://github.com/fulcrumgenomics/fgumi/pull/290))
- Migrate internals to RawRecord; drop RecordBuf variants ([#291](https://github.com/fulcrumgenomics/fgumi/pull/291))
- Migrate commands, storage, and pipeline to raw-byte records ([#292](https://github.com/fulcrumgenomics/fgumi/pull/292))
- [**breaking**] Collapse Template/DecodedRecord to raw-only; drop dead RecordBuf helpers ([#294](https://github.com/fulcrumgenomics/fgumi/pull/294))
- Migrate workspace tests to fgumi_raw_bam SamBuilder ([#295](https://github.com/fulcrumgenomics/fgumi/pull/295))
- [**breaking**] Delete vendored bam_codec; migrate simulate; propagate all remaining 48-commit work ([#296](https://github.com/fulcrumgenomics/fgumi/pull/296))

### Testing

- Benchmarks + proptests + edge cases ([#281](https://github.com/fulcrumgenomics/fgumi/pull/281))

<!-- generated by git-cliff -->

### Bug Fixes

- Cap `fgumi group` metric memory by merging per-position-group counts into per-thread accumulators instead of a `SegQueue` that grew one `AHashMap` per position group ([#285](https://github.com/fulcrumgenomics/fgumi/issues/285)).
- Replace unbounded `SegQueue<CollectedXxxMetrics>` buffering in `filter`, `clip`, `correct`, `simplex`, `duplex`, and `codec` with a shared `PerThreadAccumulator` helper, capping retained metric memory at `O(threads × distinct keys)` across all pipeline commands. Rejects buffering in the consensus commands is unchanged and tracked separately.

### Breaking Changes

- Rename the template-coordinate sort-key tag written by `fgumi zipper` from `pa` to `tc` to avoid collision with bwa-mem's `pa:f` (primary-alignment score fraction) ([#268](https://github.com/fulcrumgenomics/fgumi/issues/268)). The `--skip-pa-tags` flag on `fgumi zipper` is renamed to `--skip-tc-tags`, and the `missing_pa_tag` dedup metric field is renamed to `missing_tc_tag`. Existing fgumi-zippered BAMs must be re-zippered (or the tag must be manually renamed) before `fgumi dedup` will accept them.

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
