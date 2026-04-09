# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
