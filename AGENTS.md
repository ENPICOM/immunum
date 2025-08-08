# Repository Guidelines

## Project Structure & Module Organization
- `src/`: Rust library and CLI (`lib.rs`, `main.rs`, modules like `annotator.rs`, `sequence_io.rs`).
- `benches/`: Criterion benchmarks (`performance_benchmarks.rs`).
- `tests/`: Python API tests (`test_number_sequence.py`) and WASM tests in `tests/wasm/*.test.js`.
- `examples/`, `fixtures/`: Example data and test FASTA files.
- `wasm_build/`: Output folder for WASM artifacts (gitignored).

## Build, Test, and Development Commands
- Rust
  - Build: `cargo build` (release: `cargo build --release`).
  - Test: `cargo test`.
  - Lint: `cargo fmt --all --check` and `cargo clippy --workspace --all-targets --all-features -- -D warnings`.
  - Run CLI: `cargo run --bin immunum-cli -- -i sequences.fasta -s imgt -c igh igk`.
- Python (via maturin + uv)
  - Env: `uv venv && source .venv/bin/activate && uv sync`.
  - Develop install: `uv run maturin develop --features python`.
  - Tests: `uv run pytest`.
- WASM/Node
  - Build: `npm run build:wasm` (uses `wasm-pack`).
  - Tests: `npm test` (Vitest; runs after building WASM).

## Coding Style & Naming Conventions
- Rust: Edition 2021; format with `rustfmt`; keep modules snake_case and types/enums CamelCase (e.g., `Annotator`, `Scheme`, `Chain`). Ensure Clippy passes with no warnings.
- Python: Type-friendly API; keep files snake_case; lint/format with `ruff` (`uv run ruff check` and `uv run ruff format`).
- JS/WASM: Tests live under `tests/wasm/*.test.js`; use ES modules.

## Testing Guidelines
- Rust: Unit/integration tests via `cargo test`. Benchmarks via `cargo bench` (Criterion HTML in `target/criterion`).
- Python: Pytest files under `tests/` named `test_*.py`; ensure the extension is built (maturin develop) before running.
- WASM: Vitest specs under `tests/wasm/`; run with `npm test`.

## Commit & Pull Request Guidelines
- Commits: Keep messages clear and scoped (no strict convention enforced in history). Group logical changes; prefer imperative subject lines.
- PRs: Describe purpose and scope, include usage notes or screenshots for CLI changes, link related issues, and note performance impacts. Ensure CI parity: run `cargo fmt`, `cargo clippy`, `cargo test`, and relevant Python/WASM tests locally.

## Security & Configuration Tips
- Large inputs: prefer streaming file input (`-i <FASTA/FASTQ>`); gzipped files are supported.
- Reproducibility: record CLI flags in PRs and benchmarks; keep artifacts out of VCS (`wasm_build/`, `target/`).
