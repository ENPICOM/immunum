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

## CLI

```bash
immunum number [OPTIONS] [INPUT] [OUTPUT]
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-s, --scheme` | Numbering scheme: `imgt` (`i`), `kabat` (`k`) | `imgt` |
| `-c, --chain` | Chain filter: `h`,`k`,`l`,`a`,`b`,`g`,`d` or groups: `ig`, `tcr`, `all`. Accepts any form (`h`, `heavy`, `igh`), case-insensitive. | `ig` |
| `-f, --format` | Output format: `tsv`, `json`, `jsonl` | `tsv` |

### Input

Accepts a raw sequence, a FASTA file, or stdin (auto-detected):

```bash
immunum number EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS
immunum number sequences.fasta
cat sequences.fasta | immunum number
immunum number - < sequences.fasta
```

### Output

Writes to stdout by default, or to a file if a second positional argument is given:

```bash
immunum number sequences.fasta results.tsv
immunum number -f json sequences.fasta results.json
```

### Examples

```bash
# Kabat scheme, JSON output
immunum number -s kabat -f json EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS

# All chains (antibody + TCR), JSONL output
immunum number -c all -f jsonl sequences.fasta

# TCR sequences only, save to file
immunum number -c tcr tcr_sequences.fasta output.tsv

# Extract sequences from a TSV column and pipe in (see fixtures/ig.tsv)
tail -n +2 fixtures/ig.tsv | cut -f2 | immunum number
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="sequence") c=i} NR>1{print $c}' fixtures/ig.tsv | immunum number

# Filter TSV output to CDR3 positions (111-128 in IMGT)
immunum number sequences.fasta | awk -F'\t' '$4 >= 111 && $4 <= 128'

# Filter to heavy chain results only
immunum number -c all sequences.fasta | awk -F'\t' 'NR==1 || $2=="H"'

# Extract CDR3 sequences with jq
immunum number -f json sequences.fasta | jq '[.[] | {id: .sequence_id, numbering}]'
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
├── main.rs          # CLI binary (immunum number ...)
├── lib.rs           # Public API
├── annotator.rs     # Sequence annotation and chain detection
├── alignment.rs     # Needleman-Wunsch semi-global alignment
├── io.rs            # Input parsing (FASTA, raw) and output formatting (TSV, JSON, JSONL)
├── numbering.rs     # Numbering module entry point
├── numbering/
│   ├── imgt.rs      # IMGT numbering rules
│   └── kabat.rs     # Kabat numbering rules
├── scoring.rs       # PSSM and scoring matrices
├── types.rs         # Core domain types (Chain, Scheme, Position)
├── validation.rs    # Validation utilities
├── error.rs         # Error types
└── bin/
    ├── benchmark.rs       # Validation metrics report
    ├── debug_validation.rs # Alignment mismatch visualization
    └── speed_benchmark.rs  # Performance benchmarks
resources/
└── consensus/       # Consensus sequence CSVs (compiled into scoring matrices)
fixtures/
├── validation/      # ANARCI-numbered reference datasets
├── ig.fasta         # Example antibody sequences
└── ig.tsv           # Example TSV input
scripts/             # Python tooling for generating consensus data
```

### Design decisions

- **Semi-global alignment** forces full query consumption, preventing long CDR3 regions from being treated as trailing gaps.
- **Anchor positions** at highly conserved FR residues receive 3× gap penalties to stabilize alignment.
- **FR regions** use alignment-based numbering; **CDR regions** use scheme-specific insertion rules.
- Scoring matrices are generated at compile time from consensus data via `build.rs`.
