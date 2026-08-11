# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
