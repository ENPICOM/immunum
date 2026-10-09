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
- CLI: an argument that wasn't an existing file was read as a sequence, so a mistyped file name gave one
  invalid-sequence record and exit status 0. An argument is now read as a sequence only when it is all
  letters; otherwise a missing file is an error (exit status 1).
- CLI: the `--min-confidence` help said sequences below it are skipped; they get an error record.
- Python: unpickling a corrupt `Annotator` raised `PanicException`; it raises `ValueError`.

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
  order. The list is never empty: a sequence without a domain or an invalid one gives a single result
  with `error` set, as `number` and `segment` do, and the CLI writes one error record for it.
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
  uses, built from `Region::name()` and the new `numbering::PREFIX` and `numbering::POSTFIX`, and
  `SegmentResult::regions()` pairs each with its residues. Python, JavaScript, Polars and the
  documentation's web tool take their region names and order from these instead of their own lists.
- CLI: JSON and JSONL records carry `query_start` and `query_end`, as Python results do (`queryStart` and
  `queryEnd` in JavaScript).
- Python and JavaScript expose the same scheme lookups: `scheme_supports_chain(scheme, chain)`
  (`schemeSupportsChain` in JavaScript) tells whether a scheme numbers a chain, by the rule the
  `Annotator` constructor applies, and JavaScript gains `regionsFor(scheme, chain)`, which Python has as
  `regions_for`. Both take names and report unknown ones as every interface does; they are
  `immunum::scheme_supports_chain` and `numbering::region_spans` in Rust, with `Scheme::supports(chain)`
  and `Region::name()`. The documentation's web tool uses `schemeSupportsChain` to enable its chain
  checkboxes instead of keeping its own copy of the rule.
- Rust: `per_domain` turns what `number_domains` and `segment_domains` return into one result per
  domain, or the error as the single result. Python, JavaScript, Polars and the CLI all use it, so they
  report an invalid sequence the same way.
- Every result that failed carries `error_kind` next to `error` (`errorKind` in JavaScript), so failed
  rows can be filtered without matching messages: `invalid_sequence`, `low_confidence` or, from the
  `*_domains` methods only, `domain_too_short` (the best alignment is confident but shorter than a
  domain must be, so there is no domain to report). It is a field of Python's `NumberingResult` and
  `SegmenationResult`, of the Polars structs and of CLI JSON records, and a column of CLI TSV output.
- Python: `immunum.Error`, a `ValueError` subclass raised for every setup mistake, with a `kind`:
  `invalid_chain`, `invalid_scheme`, `unsupported_chain` or `invalid_min_confidence`. JavaScript throws
  `Error` objects with the same `kind` (`ImmunumError` in the TypeScript types).
- CLI: when any record has an error, a one-line summary goes to stderr.

### Changed
- **Breaking:** error messages drop their category prefix, which `error_kind` and `kind` now carry, and
  read the same in every interface: `Invalid sequence: sequence length 4 is below minimum 30` is now
  `sequence length 4 is below minimum 30`, `Low confidence: 0.0416 < threshold 0.5000` is now
  `alignment confidence 0.0416 is below min_confidence 0.5000`, and `Invalid chain type: unknown chain
  'IGX' …` is now `unknown chain 'IGX' …`. A scheme asked for a chain it doesn't number names the chain.
- **Breaking:** JavaScript throws `Error` objects with `message` and `kind` instead of plain strings.
- **Breaking:** CLI TSV output has an `error_kind` column after `error`, and JSON records an
  `error_kind` field.
- **Breaking:** Rust errors are split by how they're reported. `immunum::Error` is a setup mistake,
  raised by every interface: `InvalidChain`, `InvalidScheme`, the new `UnsupportedChain` (a scheme asked
  for a chain it has no rules for, which was `InvalidScheme`), `InvalidMinConfidence`, `InvalidPosition`
  and the new `WrongSequence` (a numbering paired with another sequence, which was `InvalidSequence`).
  `immunum::SequenceError` is one sequence that couldn't be numbered, returned by every interface:
  `InvalidSequence`, `LowConfidence` and the new `DomainTooShort`. `Annotator::number`, `segment`,
  `domains`, `number_domains` and `segment_domains` return `Result<_, SequenceError>`, and `domains` is
  an error instead of an empty list when there is no domain. Both enums are `#[non_exhaustive]` and have
  `kind()`. `immunum::Result` takes the error type as a defaulted second parameter. `AlignmentError`,
  `ConsensusParseError`, `PositionMappingError` and `Io` are gone: none could happen through the
  library. `ScoringMatrix::load` returns the matrix itself, and the validation functions return a boxed
  `ValidationError`.
- **Breaking:** Rust `NumberedRecord` and `SegmentedRecord` hold `result: Result<_, SequenceError>` and
  are built with `new`, instead of `success`/`failure` and separate `result`/`error` options.
- **Breaking:** JavaScript uses camelCase throughout, as JavaScript and TypeScript code expects: `number`
  results have `queryStart` and `queryEnd` instead of `query_start` and `query_end`, and the
  constructor's third parameter is declared as `minConfidence`. Every other field name is the same as in
  Python, Polars and the CLI.
- **Breaking:** JavaScript `segment` and `segmentDomains` results that failed set every region to `null`,
  as Python, Polars and the CLI do, instead of leaving the region fields out.
- **Breaking:** Polars `number` and `numbering_method` return the same struct, with the fields Python's
  `Annotator.number` returns: `chain`, `scheme`, `confidence`, `numbering`,
  `query_start`, `query_end` and `error`. `numbering` is a list of `{position, residue}` structs; explode
  and unnest it for one row per residue. `number` returned `positions` and `residues` lists instead of
  `numbering`, and had no `confidence` (#38); neither had `query_start` or `query_end`.
- Chain and scheme names are parsed once, in Rust, for every interface. Python drops its own alias
  tables, so Python, JavaScript, Polars and the CLI accept the same names and report the same error
  message for an unknown one. The chain groups `ig`, `tcr` and `all`, which only the CLI accepted, now
  work everywhere; a chain named twice (e.g. `["ig", "H"]`) is used once. An unknown chain or scheme
  name gets one message everywhere, listing every name the parser accepts (taken from the parser
  itself, the name results report first), plus the chain groups where groups are accepted.
- `min_confidence` outside `[0, 1]` is rejected by every interface. Only Python checked it; JavaScript,
  Polars and the CLI accepted any value.
- Polars builds the annotator when the expression is built, from `chains`, `scheme` and
  `min_confidence` or from `annotator=`, so an unknown name raises `immunum.Error` (a `ValueError`) right
  there. Both forms run through the same plugin function.
- Rust: `Chain` and `Scheme` parse errors are `immunum::Error` naming the accepted values, instead of
  `strum::ParseError`. New: `Chain::parse_names` and `Annotator::from_names`.
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

See the [GitHub releases page](https://github.com/ENPICOM/immunum/releases) for history prior to 1.1.0.
