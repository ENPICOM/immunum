# immunum

High-performance antibody and TCR sequence numbering in Rust.

## Overview

`immunum` is a pure-Rust library for numbering antibody and T-cell receptor (TCR) variable domain sequences. It uses Needleman-Wunsch semi-global alignment against position-specific scoring matrices (PSSM) built from consensus sequences, with BLOSUM62-based substitution scores.

### Supported chains

| Antibody | TCR |
|----------|-----|
| IGH (heavy) | TRA (alpha) |
| IGK (kappa) | TRB (beta) |
| IGL (lambda) | TRD (delta) |
| | TRG (gamma) |

### Numbering schemes

- **IMGT** — all 7 chain types
- **Kabat** — antibody chains (IGH, IGK, IGL)

Chain type is automatically detected by aligning against all loaded chains and selecting the best match.

## Usage

```rust
use immunum::{Annotator, Chain, Scheme};

let annotator = Annotator::new(
    &[Chain::IGH, Chain::IGK, Chain::IGL],
    Scheme::IMGT,
).unwrap();

let sequence = "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN";

let alignment = annotator.annotate(sequence).unwrap();

println!("Chain: {}", alignment.chain);        // IGH
println!("Confidence: {:.2}", alignment.confidence);

let numbering = alignment.numbering(Scheme::IMGT);
for (aa, pos) in sequence.chars().zip(numbering.iter()) {
    println!("{} -> {}", aa, pos);
}
```

## Development

The project follows a test-driven development workflow.

```bash
# Run tests
cargo test

# Lint
cargo fmt
cargo clippy

# Run validation benchmarks (writes results to BENCHMARKS.toml)
cargo run --release --bin benchmark
```

### Project structure

```
src/
├── lib.rs           # Public API
├── annotator.rs     # Sequence annotation and chain detection
├── alignment.rs     # Needleman-Wunsch semi-global alignment
├── numbering/       # Scheme-specific numbering logic
│   ├── imgt.rs
│   └── kabat.rs
├── scoring.rs       # PSSM and scoring matrices
├── types.rs         # Core domain types (Chain, Scheme, Position)
├── validation.rs    # Validation utilities
└── error.rs         # Error types
resources/
└── consensus/       # Consensus sequence CSVs (compiled into scoring matrices)
fixtures/
└── validation/      # ANARCI-numbered reference datasets
scripts/             # Python tooling for generating consensus data
```

### Design decisions

- **Semi-global alignment** forces full query consumption, preventing long CDR3 regions from being treated as trailing gaps.
- **Anchor positions** at highly conserved FR residues receive 3× gap penalties to stabilize alignment.
- **FR regions** use alignment-based numbering; **CDR regions** use scheme-specific insertion rules.
- Scoring matrices are generated at compile time from consensus data via `build.rs`.
