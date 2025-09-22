# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Immunum is a high-performance Rust library for numbering antibody and T-cell receptor sequences using standard immunological numbering schemes (IMGT and Kabat). The project provides:
- A CLI tool (`immunum-cli`)
- Python bindings via PyO3
- WebAssembly bindings for web applications

## Build Commands

### Rust Development
```bash
# Build the CLI tool
cargo build

# Run all Rust tests
cargo test

# Run with debug output for specific tests
RUST_BACKTRACE=1 cargo test test_name -- --nocapture

# Run linting and checks
cargo clippy
cargo check

# Clean build artifacts
cargo clean
```

### Python Development
```bash
# Set up Python environment (uses uv package manager)
uv venv && source .venv/bin/activate
uv sync

# Build and install Python package in development mode
uv run maturin develop --features python

# Build release wheel
uv run maturin build --features python --release

# Run Python tests
uv run pytest
```

### WebAssembly Development
```bash
# Build WASM module
wasm-pack build --target web --out-dir wasm_build --features wasm --no-default-features

# Run WASM tests
node test-wasm-node.js
npm test
```

## Architecture

### Core Architecture
The codebase follows a modular Rust architecture with the main entry point being the `Annotator` struct in `src/annotator.rs`. Key architectural components:

- **Annotator**: Main API entry point that handles sequence numbering with pre-built numbering schemes
- **NumberingScheme**: Core data structure containing reference sequences and scoring matrices for specific chain types
- **Sequence Processing**: Handles FASTA/FASTQ file parsing with support for gzip compression
- **Alignment Engine**: Custom Needleman-Wunsch implementation optimized for immunoglobulin sequences

### Key Modules
- `annotator.rs`: Main API with parallel processing support via Rayon
- `schemes.rs`: Numbering scheme data and initialization
- `needleman_wunsch.rs`: Sequence alignment implementation
- `prefiltering.rs`: Performance optimization for multi-chain analysis
- `sequence.rs`: File format handling (FASTA/FASTQ/gzip)
- `result.rs`: Output formatting and annotation results
- `python_bindings.rs`/`wasm_bindings.rs`: Language bindings

### Features System
The project uses Cargo features for conditional compilation:
- `python`: Enables PyO3 Python bindings
- `wasm`: Enables WebAssembly bindings
- Default features include both `python` and `wasm`

### Parallel Processing
The library extensively uses Rayon for parallel processing in:
- Multi-sequence processing (`number_sequences`)
- File processing (`number_file`)
- CLI tool processing

### Testing Strategy
- Rust unit tests cover core functionality
- Integration tests in `tests/` directory for validation
- Python binding tests via pytest
- WASM binding tests via Node.js
- Performance benchmarks via `benchmark.py`

## Development Commands

### Running Single Tests
```bash
# Run specific Rust test
cargo test test_name

# Run Python tests with verbose output
uv run pytest -v

# Run validation tests
cargo test abpdseq_validation_test
```

### Debugging
```bash
# Run with Rust backtrace
RUST_BACKTRACE=1 cargo run -- [args]

# Test with debug output
cargo test test_name -- --nocapture
```

### Multi-threading
The CLI supports parallel processing via the `--threads` option. Set thread count for development testing:
```bash
./target/debug/immunum-cli --threads 4 [other args]
```

## Package Managers
- **Rust**: Standard `cargo` for Rust dependencies
- **Python**: `uv` (modern Python package manager, faster than pip)
- **JavaScript**: `npm` for WASM testing dependencies