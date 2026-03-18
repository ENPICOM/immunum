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

let result = annotator.number(sequence).unwrap();

println!("Chain: {}", result.chain);        // IGH
println!("Confidence: {:.2}", result.confidence);

for (aa, pos) in sequence.chars().zip(result.positions.iter()) {
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

To orchestrate a project between cargo and python, we use [`task`](http://taskfile.dev).
You can install it with:

```bash
uv tool install go-task-bin
```

And then run `task` or `task --list-all` to get the full list of available tasks.

By default, `dev` profile will be used in all but `behcnmark-*` taks, but you can change it
via providing `PROFILE=release` to your task.

Also, by default, `task` caches results, but you can ignore it by running `task my-task -f`.

### Building local environment

```bash
# build a dev environment
task build

# build a dev environment with --release flag
task build PROFILE=release
```

### Testing

```bash
task test-rust    # test only rust code
task test-python  # test only python code
task test         # test all code
```

### Linting

```bash
task format  # formats python and rust code
task lint    # runs linting for python and rust
```

### Benchmarking

There are multiple benchmarks in the repository. For full list, see `task | grep behcmark`:

```bash
$ task | grep benchmark
* benchmark-accuracy:           Accuracy benchmark across all fixtures (1k sequences, 7 rounds each)
* benchmark-cli:                Behcmark correctness of the CLI tool
* benchmark-comparison:         Speed + correctness benchmark: immunum vs antpack vs anarci (1k IGH sequences)
* benchmark-scaling:            Scaling benchmark: sizes 100..10M (10x steps), 1 round, H/imgt. Pass CLI_ARGS to filter tools, e.g. -- --tools immunum
* benchmark-speed:              Speed benchmark across dataset sizes (100 to 1M sequences, 7 rounds, H/imgt)
* benchmark-speed-polars:       Speed benchmark for immunum polars across all chain/scheme fixtures
```

## Project structure
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
immunum/
└── _internal.pyi    # python stub file for pyo3
└── polars.py        # polars extension module
└── python.py        # python module
```

### Design decisions

- **Semi-global alignment** forces full query consumption, preventing long CDR3 regions from being treated as trailing gaps.
- **Anchor positions** at highly conserved FR residues receive 3× gap penalties to stabilize alignment.
- **FR regions** use alignment-based numbering; **CDR regions** use scheme-specific insertion rules.
- Scoring matrices are generated at compile time from consensus data via `build.rs`.
