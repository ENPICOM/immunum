# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
