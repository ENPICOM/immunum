![Immunum Logo](https://raw.githubusercontent.com/ENPICOM/immunum/master/docs/assets/immunum_logotype.svg)

Immunum is a high-performance antibody and TCR sequence numbering tool for Rust, Python, R, Polars and JS/TS.

Try it in your browser: [interactive demo](https://immunum.enpicom.com/demo/).

[![Crates.io](https://img.shields.io/crates/v/immunum)](https://crates.io/crates/immunum)
[![PyPI](https://img.shields.io/pypi/v/immunum)](https://pypi.org/project/immunum/)
[![npm](https://img.shields.io/npm/v/immunum)](https://www.npmjs.com/package/immunum)
[![License: MIT](https://img.shields.io/crates/l/immunum)](LICENSE)
[![CI](https://img.shields.io/github/actions/workflow/status/ENPICOM/immunum/ci.yml?label=CI)](https://github.com/ENPICOM/immunum/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-immunum.enpicom.com-blue)](https://immunum.enpicom.com)

## Overview

`immunum` is a library for numbering antibody and T-cell receptor (TCR) variable domain sequences. It uses Needleman-Wunsch semi-global alignment against position-specific scoring matrices built from consensus sequences, with BLOSUM62-based substitution scores.

Available as:

- **Rust crate** — core library and CLI
- **Python package** — with a [Polars](https://pola.rs) plugin for vectorized batch processing
- **R package** — with rayon-parallel batch processing, distributed via r-universe
- **npm package** — for Node.js and browsers

### Supported chains

| Antibody     | TCR         |
| ------------ | ----------- |
| IGH (heavy)  | TRA (alpha) |
| IGK (kappa)  | TRB (beta)  |
| IGL (lambda) | TRD (delta) |
|              | TRG (gamma) |

Chain codes: `H` (IGH), `K` (IGK), `L` (IGL), `A` (TRA), `B` (TRB), `D` (TRD), `G` (TRG).

Chain type is automatically detected by aligning against all loaded chains and selecting the best match.

### Numbering schemes

- **IMGT** — all 7 chain types
- **Kabat** — antibody chains (IGH, IGK, IGL)

## Table of Contents

- [Python](#python)
  - [Installation](#installation)
  - [Numbering](#numbering)
  - [Segmentation](#segmentation)
  - [Polars plugin](#polars-plugin)
- [R](#r)
  - [Installation](#installation-1)
  - [Numbering](#numbering-1)
  - [Segmentation](#segmentation-1)
- [JavaScript / npm](#javascript--npm)
  - [Installation](#installation-2)
  - [Usage](#usage)
- [Rust](#rust)
  - [Installation](#installation-3)
  - [Usage](#usage-1)
- [CLI](#cli)
  - [Options](#options)
  - [Input](#input)
  - [Output](#output)
  - [Examples](#examples)
- [Development](#development)
- [Project structure](#project-structure)

## Python

### Installation

```bash
pip install immunum
```

### Numbering

```python
from immunum import Annotator

annotator = Annotator(chains=["H", "K", "L"], scheme="imgt")

sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"

result = annotator.number(sequence)
print(result.chain)       # H
print(result.confidence)  # 0.78
print(result.numbering)   # {"1": "Q", "2": "V", "3": "Q", ...}
```

### Segmentation

`segment` splits the sequence into FR/CDR regions:

```python
from immunum import Annotator

annotator = Annotator(chains=["H", "K", "L"], scheme="imgt")

sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"

result = annotator.segment(sequence)
assert result.fr1 == 'QVQLVQSGAEVKRPGSSVTVSCKAS'
assert result.cdr1 == 'GGSFSTYA'
assert result.fr2 == 'LSWVRQAPGRGLEWMGG'
assert result.cdr2 == 'VIPLLTIT'
assert result.fr3 == 'NYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYC'
assert result.cdr3 == 'AREGTTGKPIGAFAH'
assert result.fr4 == 'WGQGTLVTVSS'
```

### Polars plugin

For batch processing, `immunum.polars` registers elementwise Polars expressions:

```python
import polars as pl
import immunum.polars as imp

df = pl.DataFrame({"sequence": [
    "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS",
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
]})

# Add a struct column with chain, scheme, confidence, numbering
result = df.with_columns(
    imp.number(pl.col("sequence"), chains=["H", "K", "L"], scheme="imgt").alias("numbered")
)

# Add a struct column with FR/CDR segments
result = df.with_columns(
    imp.segment(pl.col("sequence"), chains=["H", "K", "L"], scheme="imgt").alias("segmented")
)
```

The `number` expression returns a struct with fields `chain`, `scheme`, `confidence`, and `numbering` (a struct of position→residue). The `segment` expression returns a struct with fields `fr1`, `cdr1`, `fr2`, `cdr2`, `fr3`, `cdr3`, `fr4`, `prefix`, `postfix`.

## R

### Installation

```r
# install.packages("remotes")
remotes::install_github("ENPICOM/immunum", subdir = "r-immunum", build = FALSE)
```

<details>
<summary>Building from source</summary>

Building from source requires a [Rust toolchain](https://rustup.rs). On Windows you also need [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

Check your setup with:

```r
# install.packages("rextendr")
rextendr::rust_sitrep()
```

Install Rust if needed:

```bash
# macOS / Linux
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Windows: download and run rustup-init.exe from https://rustup.rs
```

</details>

### Numbering

```r
library(immunum)

ann <- Annotator$new(chains = c("H", "K", "L"), scheme = "imgt")

sequence <- "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"

result <- ann$number(sequence)
result$chain       # "H"
result$confidence  # 0.78
result$numbering   # named character vector: c("1"="Q", "2"="V", ...)
result$error       # NULL on success, error message on failure
```

### Segmentation

```r
result <- ann$segment(sequence)
result$fr1   # "QVQLVQSGAEVKRPGSSVTVSCKAS"
result$cdr1  # "GGSFSTYA"
result$cdr3  # "AREGTTGKPIGAFAH"
result$fr4   # "WGQGTLVTVSS"
```

## JavaScript / npm

### Installation

```bash
npm install immunum
```

### Usage

```js
const { Annotator } = require("immunum");

const annotator = new Annotator(["H", "K", "L"], "imgt");

const sequence =
  "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";

const result = annotator.number(sequence);
console.log(result.chain);      // "H"
console.log(result.confidence); // 0.97
console.log(result.numbering);  // { "1": "Q", "2": "V", ... }

const segments = annotator.segment(sequence);
console.log(segments.cdr3); // "AREGTTGKPIGAFAH"

annotator.free(); // or use `using annotator = new Annotator(...)` with explicit resource management
```

## Rust

### Installation

Add to `Cargo.toml`:

```toml
[dependencies]
immunum = "1.1"
```

### Usage

```rust
use immunum::{Annotator, Chain, Scheme};

let annotator = Annotator::new(
    &[Chain::IGH, Chain::IGK, Chain::IGL],
    Scheme::IMGT,
    None, // uses default min_confidence of 0.5
).unwrap();

let sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";

let result = annotator.number(sequence).unwrap();
println!("Chain: {}", result.chain);        // IGH
println!("Confidence: {:.2}", result.confidence);
for (aa, pos) in sequence.chars().zip(result.positions.iter()) {
    println!("{} -> {}", aa, pos);
}

let segments = annotator.segment(sequence).unwrap();
println!("CDR3: {}", segments.cdr3);
```

## CLI

```bash
immunum number [OPTIONS] [INPUT] [OUTPUT]
```

### Options

| Flag           | Description                                                                                                                        | Default |
| -------------- | ---------------------------------------------------------------------------------------------------------------------------------- | ------- |
| `-s, --scheme` | Numbering scheme: `imgt` (`i`), `kabat` (`k`)                                                                                      | `imgt`  |
| `-c, --chain`  | Chain filter: `h`,`k`,`l`,`a`,`b`,`g`,`d` or groups: `ig`, `tcr`, `all`. Accepts any form (`h`, `heavy`, `igh`), case-insensitive. | `ig`    |
| `-f, --format` | Output format: `tsv`, `json`, `jsonl`                                                                                              | `tsv`   |

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

By default, `dev` profile will be used in all but `benchmark-*` tasks, but you can change it
via providing `PROFILE=release` to your task.

Also, by default, `task` caches results, but you can ignore it by running `task my-task -f`.

### Building local environment

```bash
# build a dev environment
task build-local

# build a dev environment with --release flag
task build-local PROFILE=release
```

### Testing

```bash
task test-rust    # test only rust code
task test-python  # test only python code
task test         # test all code
task r:test       # test R package
```

### Linting

```bash
task format  # formats python and rust code
task lint    # runs linting for python and rust
```

### Benchmarking

There are multiple benchmarks in the repository. For full list, see `task | grep benchmark`:

```bash
$ task | grep benchmark
* benchmark-accuracy:           Accuracy benchmark across all fixtures (1k sequences, 7 rounds each)
* benchmark-cli:                Benchmark correctness of the CLI tool
* benchmark-comparison:         Speed + correctness benchmark: immunum vs antpack vs anarci (1k IGH sequences)
* benchmark-scaling:            Scaling benchmark: sizes 100..10M (10x steps), 1 round, H/imgt. Pass CLI_ARGS to filter tools, e.g. -- --tools immunum
* benchmark-speed:              Speed benchmark across dataset sizes (100 to 1M sequences, 7 rounds, H/imgt)
* benchmark-speed-polars:       Speed benchmark for immunum polars across all chain/scheme fixtures
```

R vs Python scaling benchmark:

```bash
task r:bench          # R vs Python polars scaling benchmark with chart
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
├── _internal.pyi    # python stub file for pyo3
├── polars.py        # polars extension module
└── python.py        # python module
r-immunum/           # R package (extendr bindings)
├── R/               # R source files (Annotator, polars wrappers, normalization)
├── src/extendr/     # Rust shim crate (path dep to parent immunum crate)
└── tests/testthat/  # testthat tests (annotator, polars, validation, cross-language)
```

### Design decisions

- **Semi-global alignment** forces full query consumption, preventing long CDR3 regions from being treated as trailing gaps.
- **Anchor positions** at highly conserved FR residues receive 3× gap penalties to stabilize alignment.
- **FR regions** use alignment-based numbering; **CDR regions** use scheme-specific insertion rules.
- Scoring matrices are generated at compile time from consensus data via `build.rs`.
