# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- Polars: `number`, `segment`, `numbering_method` and `segmentation_method` paired the numbering with the
  sequence from its first residue instead of from `query_start`, so every residue after a leader or other
  leading flank was shifted by the length of that flank (#53).
- `segment` dropped the residues the aligner leaves out before and after the domain instead of putting
  them in `prefix` and `postfix`, in Rust, Python, Polars and JavaScript. The regions now always join
  back into the input sequence, as documented (#58).
- Python: an `Annotator` used from a thread other than the one that created it raised a `PanicException`
  ("unsendable, but sent to another thread"). It can now be shared across threads.
- The documentation's web tool upper-cased the sequence before numbering it, so it showed residues the
  user hadn't entered; it now numbers the sequence as given, like every other interface. It still drops
  whitespace from the text box, so a wrapped or spaced sequence can be pasted.

### Added
- Rust: `Annotator::domains(sequence)` finds every variable domain in a sequence, ordered by position.
  It keeps the best alignment's domain, then searches the residues before and after it the same way, each
  on its own, until a part is shorter than `MIN_SEQUENCE_LENGTH` or aligns with too little confidence.
  Each domain numbers like `Annotator::number` on the residues it was found in, so a domain missing its
  N-terminal residues after another domain is found as it would be at the start of a sequence. A domain
  shorter than `MIN_SEQUENCE_LENGTH` is not reported, but the residues around it are still searched.
- Rust: `Domain::number(scheme)` numbers a found domain under any scheme from the alignment that found it,
  without aligning again. The first domain of a single-domain sequence numbers exactly like
  `Annotator::number`, including AHo's light-chain position 149. A `Domain`'s span, chain and confidence
  are read through methods (`query_start()`, `chain()`, ...), so they can't drift from its alignment.
- Every interface numbers and segments every variable domain in a sequence, such as both domains of an
  scFv: `number_domains` and `segment_domains` in Rust (`Annotator`), Python (`Annotator`) and Polars,
  `numberDomains` and `segmentDomains` in JavaScript, and `--all-domains` on the CLI's `number` and
  `segment`. Each domain's result is what `number` or `segment` returns for that domain, in sequence
  order. No domain gives an empty list (no CLI record); an invalid sequence gives a single result with
  `error` set, as `number` and `segment` do.
- `segment_domains` puts every residue in exactly one domain's regions: a domain's `prefix` holds the
  residues since the previous domain (or the start of the sequence), and only the last domain has the
  residues after it as its `postfix`, so all domains' regions in order rebuild the sequence.
- CLI: `immunum segment` splits sequences into FR/CDR regions, one record per sequence (or per domain
  with `--all-domains`) with a column or field per region. It takes the same options as `number`. With
  `--all-domains`, both commands add a 0-based `domain` column (TSV) or field (JSON).
- Polars: every expression (`number`, `number_domains`, `segment`, `segment_domains`) takes either
  `chains`, `scheme` and `min_confidence`, or a prebuilt `Annotator` as `annotator=`.
- Rust: `NumberedRecord` and the new `SegmentedRecord` have a `domain` field and `in_domain`, and
  `OutputFormat::write_header` and `write_record` take an `all_domains` switch, with
  `write_segment_header` and `write_segment_record` for segments.
- Rust: `NumberingResult::residues(sequence)` pairs each numbered position with its residue, and
  `NumberingResult::segment(sequence)` splits a numbering into FR/CDR regions, flanks included, without
  numbering again. `sequence` is the whole sequence that was numbered or searched; one too short for the
  numbering is an error. Python, JavaScript, Polars and the CLI now all use these.
- Rust: `numbering::SEGMENT_NAMES` lists the segments in sequence order by the names every interface
  uses, and `SegmentResult::regions()` pairs each with its residues. Python, JavaScript, Polars and the
  documentation's web tool take their region names and order from these instead of their own lists.
- CLI: JSON and JSONL records carry `query_start` and `query_end`, as Python results do (`queryStart` and
  `queryEnd` in JavaScript).
- JavaScript: `schemeSupportsChain(scheme, chain)` tells whether a scheme numbers a chain, by the rule
  the `Annotator` constructor applies; Rust: `Scheme::supports(chain)`. The documentation's web tool
  uses it to enable its chain checkboxes instead of keeping its own copy of the rule.

### Changed
- **Breaking:** JavaScript uses camelCase throughout, as JavaScript and TypeScript code expects: `number`
  results have `queryStart` and `queryEnd` instead of `query_start` and `query_end`, and the
  constructor's third parameter is declared as `minConfidence`.
- **Breaking:** Polars `number` and `numbering_method` return the same struct, with the fields Python's
  `Annotator.number` returns: `chain`, `scheme`, `confidence`, `numbering`,
  `query_start`, `query_end` and `error`. `numbering` is a list of `{position, residue}` structs; explode
  and unnest it for one row per residue. `number` returned `positions` and `residues` lists instead of
  `numbering`, and had no `confidence` (#38); neither had `query_start` or `query_end`.
- Chain and scheme names are parsed once, in Rust, for every interface. Python drops its own alias
  tables, so Python, JavaScript, Polars and the CLI accept the same names and report the same error
  message for an unknown one. The chain groups `ig`, `tcr` and `all`, which only the CLI accepted, now
  work everywhere; a chain named twice (e.g. `["ig", "H"]`) is used once.
- `min_confidence` outside `[0, 1]` is rejected by every interface. Only Python checked it; JavaScript,
  Polars and the CLI accepted any value.
- Polars `number` and `segment` report an unknown chain or scheme when the expression runs, as a
  `ComputeError`, instead of raising `ValueError` when the expression is built.
- Rust: `Chain` and `Scheme` parse errors are `immunum::Error` (`InvalidChain`/`InvalidScheme`) naming
  the accepted values, instead of `strum::ParseError`. New: `Chain::parse_names` and
  `Annotator::from_names`, and `Error::InvalidMinConfidence`.
- Rust: `Annotator` is `Send + Sync`: its alignment buffer is per thread, so one annotator can serve many
  threads. A thread keeps at most about 650 KB of alignment buffer between calls, what a 1,000-residue
  sequence needs; a longer sequence's buffer is freed when its call returns. `Annotator::number` and
  `Annotator::segment` are built on the pieces above; apart from the flanks fixed above, their results
  are unchanged.

### Deprecated
- Polars `numbering_method` and `segmentation_method`. Use `number(expr, annotator=annotator)` and
  `segment(expr, annotator=annotator)`, which return the same.

### Removed
- Rust: `Chain::parse_chain_spec`. Use `Chain::parse_names(spec.split(','))`.

## [1.3.3] - 2026-10-01

### Fixed
- A sequence ending exactly one residue past a domain's last framework position no longer gets that
  residue numbered as part of the domain. A heavy chain + one residue numbered it as an insertion after
  the last position (IMGT `128A`, Kabat `113A`), and a kappa or lambda chain + one residue numbered it at
  IMGT 128, a position no light chain has. In both cases it also extended `query_end` and FR4. Two or more
  trailing residues were already left out. Numbering now matches ANARCI, including AHo's position 149
  after a light chain.
- The documentation release build failed because `task docs-build -- --strict` forwarded `--strict`
  through the `build-wasm-web` dependency into `cargo build`. The wasm-pack task no longer takes
  `CLI_ARGS`.

### Changed
- The TRB validation fixture no longer expects a residue at IMGT 128. The residue after a TRB J-region is
  the first residue of the constant region; TRB germline J genes and ANARCI end at 127.

## [1.3.2] - 2026-09-29

### Added
- The documentation footer now shows the version it was built from, linking to the matching
  GitHub release.

### Changed
- The release workflow builds the documentation site and publishes it to S3

## [1.3.1] - 2026-09-02

### Added
- `regions_for(scheme, chain)` in Python, returning a scheme's FR/CDR boundaries as
  `{region: (start, end)}` with both ends inclusive and keyed by the lowercase region names
  `segment()` uses. These are the same tables numbering assigns residues to, now readable without
  numbering a sequence. Raises `ValueError` for an unknown scheme or chain, and for a pair immunum
  cannot number (Kabat, Chothia, Martin and AHo on TCR chains).
- Rust: `RegionDefinition::spans()` reads the stored region ends back as inclusive `(start, end)`
  pairs, and `Scheme::validate_chain()` is the scheme/chain check `Annotator::new` previously did
  inline.

## [1.3.0] - 2026-08-28

### Added
- Chothia, Martin (extended Chothia) and AHo numbering schemes for antibody chains (IGH, IGK, IGL),
  each derived from the internal IMGT numbering like Kabat. Select them as `chothia` / `c`,
  `martin` / `m`, `aho` / `a` from the CLI (`-s`), Python (`Annotator(..., scheme=...)`, including the
  Polars plugin), WASM and Rust (`Scheme::Chothia`, `Scheme::Martin`, `Scheme::Aho`). Per-residue
  agreement against ANARCI-numbered reference sets is 99.98–100% for all three (see `BENCHMARKS.toml`);
  TCR chains remain IMGT-only.
- Validation fixtures for the new schemes: `fixtures/validation/ab_{H,K,L}_{chothia,martin,aho}.csv`.

### Changed
- **Segmentation breaking change for Kabat:** FR/CDR boundaries are now scheme *and* chain specific
  instead of one shared table, so `segment()` returns different region strings under Kabat.
  Heavy chains use CDR-H1 31–35 (plus 35A/35B), CDR-H2 50–65, CDR-H3 95–102; light chains use
  CDR-L1 24–34, CDR-L2 50–56, CDR-L3 89–97 with FR4 ending at 107. The previous table applied
  CDR1 26–35, CDR2 51–57, CDR3 93–100 to both chains. Position *numbers* are unchanged; only the
  region a residue is assigned to changes. IMGT segmentation is unaffected.
- The CDR-center insertion penalty dropped from -3.0 to -2.0. Ultralong CDR3 domains (e.g. bovine
  VHH-like heavy chains absorbing ~48 residues at IMGT 111/112) are no longer suffix-clipped before
  FR4, so `query_end` and the numbering now cover the full domain. Raised AHo heavy per-residue
  accuracy from 99.87% to 100% and IMGT TRB from 99.79% to 99.80%.
- AHo light chains now emit the trailing position 149 when a residue follows AHo 148, matching
  ANARCI's `number_aho` tail rule. Adds one residue to the numbered range for affected sequences.
- **Rust API breaking change:** `numbering::segment` takes an extra `chain: Chain` argument, and
  `Scheme` gained three variants, so exhaustive `match` arms on `Scheme` in downstream code need
  updating.

## [1.2.0] - 2026-08-04

### Fixed
- Kabat numbering no longer panics when a CDR aligns to a single residue. Heavy CDR1, heavy CDR3 and light CDR2 each had a `deletion_order` shorter than `base_len - 1`, so deleting every base position but one indexed past the end of the slice. The orders now follow ANARCI, which leaves position 23 for heavy CDR1 and 94 for heavy CDR3. Affected roughly 1 in 20 000 truncated but productive VHH reads, all of which numbered correctly under IMGT.

## [1.1.2] - 2026-05-21

### Added
- `py.typed` marker so downstream type checkers (mypy, pyright) pick up inline annotations and the bundled `_internal.pyi` stub per PEP 561.

## [1.1.1] - 2026-05-18

### Changed
- Increased maximum allowed input sequence length to 10000.

### Fixed
- Updated packages to address security vulnerabilities.

## [1.1.0] - 2026-04-09

### Added
- Interactive WASM numbering tool and demo page with UI for confidence settings.
- Sequence validation with length constraints and character checks.
- Error handling in numbering and segmentation processes.
- Logo in README and docs.

### Changed
- **WASM breaking change:** the `Numbering` output type is now `Map<string, string>` instead of `Record<string, string>`. Access residues with `numbering.get("112A")` instead of `numbering["112A"]`, and iterate with `for (const [pos, aa] of numbering)` or `numbering.entries()` instead of `Object.entries(numbering)`. This preserves insertion-code ordering (e.g. `"111A"` / `"112A"` stay between `"111"` and `"112"`).
- `number()` result now exposes `query_start` / `query_end` as inclusive 0-indexed integers marking the aligned region of the input sequence (both `null` on error).

### Fixed
- Clearer error message on numbering failure.
- macOS CI build.

## [1.0.0] - Prior release

See the [GitHub releases page](https://github.com/ENPICOM/immunum-rs/releases) for history prior to 1.1.0.
