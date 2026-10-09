![Immunum Logo](https://raw.githubusercontent.com/ENPICOM/immunum/master/docs/assets/immunum_logotype.svg)

Immunum is a high-performance antibody and TCR sequence numbering tool for Rust, Python, Polars and JS/TS.

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
- **Chothia** — antibody chains (IGH, IGK, IGL)
- **Martin** (extended Chothia) — antibody chains (IGH, IGK, IGL)
- **AHo** — antibody chains (IGH, IGK, IGL)

Every scheme is derived from the internal IMGT numbering. Region (FR/CDR) boundaries follow each
scheme's own definition and differ between heavy and light chains.

### Errors

immunum fails in two ways, and every interface reports each the same way:

- **A setup mistake is raised** before any sequence is numbered: Python raises `immunum.Error` (a `ValueError`) and JavaScript throws an `Error`, both with a `kind`; Polars raises when the expression is built; the CLI prints `error: …` and exits with status 1.
- **A sequence that can't be numbered is returned**, so a batch never stops for it: its result has `error` (the message) and `error_kind` (`errorKind` in JavaScript) set, and every other field empty.

| Kind                     | Raised or returned | When                                                         |
| ------------------------ | ------------------ | ------------------------------------------------------------ |
| `invalid_chain`          | raised             | an unknown chain name, or no chains                          |
| `invalid_scheme`         | raised             | an unknown scheme name                                       |
| `unsupported_chain`      | raised             | a scheme other than IMGT asked to number a TCR chain         |
| `invalid_min_confidence` | raised             | `min_confidence` outside `[0, 1]`                            |
| `invalid_sequence`       | returned           | a sequence too short, too long, or holding a non-letter      |
| `low_confidence`         | returned           | no alignment reaches `min_confidence`                        |

```python
import immunum

try:
    immunum.Annotator(chains=["IGX"], scheme="imgt")
except immunum.Error as e:
    print(e.kind)  # invalid_chain

result = immunum.Annotator(chains=["ig"], scheme="imgt").number("AAAA")
print(result.error_kind)  # invalid_sequence
print(result.error)       # sequence length 4 is below minimum 30
```

## Table of Contents

- [Python](#python)
  - [Installation](#installation)
  - [Numbering](#numbering)
  - [Segmentation](#segmentation)
  - [Multiple domains](#multiple-domains)
  - [Polars plugin](#polars-plugin)
- [JavaScript / npm](#javascript--npm)
  - [Installation](#installation-1)
  - [Usage](#usage)
- [Rust](#rust)
  - [Installation](#installation-2)
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

### Multiple domains

`number` and `segment` work on the best-scoring domain. `number_domains` and `segment_domains` do the same for every variable domain in a sequence, such as both domains of an scFv, in sequence order. Each result is what `number` or `segment` returns for that domain:

```python
from immunum import Annotator

heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
kappa = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
linker = "GGGGSGGGGSGGGGS"

annotator = Annotator(chains=["ig"], scheme="imgt")
domains = annotator.number_domains(heavy + linker + kappa)
assert [d.chain for d in domains] == ["H", "K"]

heavy_regions, kappa_regions = annotator.segment_domains(heavy + linker + kappa)
assert kappa_regions.prefix == linker
assert kappa_regions.cdr3 == "QQHYTTPPT"
```

In `segment_domains`, every residue lands in exactly one domain's regions: a domain's `prefix` holds the residues since the previous domain (or the start of the sequence), and only the last domain has the residues after it as its `postfix`, so all domains' regions in order rebuild the sequence.

The lists are empty when no domain aligns with enough confidence. An invalid sequence gives a single result with `error` set, as `number` and `segment` would.

A domain that lacks its first IMGT positions (a light chain starting at position 2, say) and directly follows other residues, such as a linker, can have the residue just before it numbered as its first position. IMGT position 1 is so variable that the sequence alone can't tell a linker residue from the domain's own first residue.

### Polars plugin

For batch processing, `immunum.polars` registers elementwise Polars expressions:

```python
import polars as pl
import immunum.polars as imp

df = pl.DataFrame({"sequence": [
    "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS",
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
]})

# Add a struct column with the fields Annotator.number returns
result = df.with_columns(
    imp.number(pl.col("sequence"), chains=["H", "K", "L"], scheme="imgt").alias("numbered")
)

# Add a struct column with FR/CDR segments
result = df.with_columns(
    imp.segment(pl.col("sequence"), chains=["H", "K", "L"], scheme="imgt").alias("segmented")
)
```

The `number` expression returns a struct with the fields `Annotator.number` returns: `chain`, `scheme`, `confidence`, `numbering` (a list of `{position, residue}` structs), `query_start`, `query_end`, `error` and `error_kind`. Explode `numbering` and unnest it for one row per residue. The `segment` expression returns a struct with fields `prefix`, `fr1`, `cdr1`, `fr2`, `cdr2`, `fr3`, `cdr3`, `fr4`, `postfix`, `error` and `error_kind`. `number_domains` and `segment_domains` return a list of those structs per sequence, one per domain.

Every expression takes either `chains`, `scheme` and `min_confidence`, or a prebuilt `Annotator` as `annotator=`. Either way the annotator is built when the expression is, so an unknown name raises `immunum.Error` right there.

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

// Every domain, e.g. both domains of an scFv; each is what number() / segment() returns for it
const domains = annotator.numberDomains(sequence);
const domainSegments = annotator.segmentDomains(sequence);

annotator.free(); // or use `using annotator = new Annotator(...)` with explicit resource management
```

`schemeSupportsChain(scheme, chain)` tells whether a scheme numbers a chain before you build an annotator, e.g. to offer only valid choices in a UI. `regionsFor(scheme, chain)` returns the scheme's FR/CDR boundaries for the chain, both ends inclusive, without numbering a sequence:

```js
const { regionsFor, schemeSupportsChain } = require("immunum");

schemeSupportsChain("kabat", "H"); // true
schemeSupportsChain("kabat", "B"); // false: only IMGT numbers TCR chains
regionsFor("kabat", "H").cdr1;    // [31, 35]
```

A setup mistake throws an `Error` with a `kind`; a sequence that can't be numbered comes back with `error` and `errorKind` set (see [Errors](#errors)):

```js
try {
  new Annotator(["IGX"], "imgt");
} catch (err) {
  console.log(err.kind);    // "invalid_chain"
  console.log(err.message); // "unknown chain 'IGX' (options: ...)"
}

const failed = new Annotator(["ig"], "imgt").number("AAAA");
console.log(failed.errorKind); // "invalid_sequence"
console.log(failed.error);     // "sequence length 4 is below minimum 30"
```

## Rust

### Installation

Add to `Cargo.toml`:

```toml
[dependencies]
immunum = "0.9"
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

// Every domain, e.g. both domains of an scFv, numbered under the annotator's scheme
for domain in annotator.number_domains(sequence).unwrap() {
    println!("{} at {}..={}", domain.chain, domain.query_start, domain.query_end);
}
for segments in annotator.segment_domains(sequence).unwrap() {
    println!("CDR3: {}", segments.cdr3);
}
```

Setting up returns `immunum::Error`; numbering one sequence returns `immunum::SequenceError`, so a sequence's result can only fail with `InvalidSequence` or `LowConfidence`. Both have `kind()`, the code every other interface reports.

## CLI

```bash
immunum number [OPTIONS] [INPUT] [OUTPUT]
immunum segment [OPTIONS] [INPUT] [OUTPUT]
```

`number` writes one record per numbered residue (TSV) or per sequence (JSON). `segment` writes one record per sequence, with a column or field per region: `prefix`, `fr1`, `cdr1`, `fr2`, `cdr2`, `fr3`, `cdr3`, `fr4`, `postfix`, `error` and `error_kind`. A sequence that can't be numbered gets a record with `error` and `error_kind` set, and the run carries on; a summary of how many failed goes to stderr.

### Options

Both commands take the same options.

| Flag           | Description                                                                                                                        | Default |
| -------------- | ---------------------------------------------------------------------------------------------------------------------------------- | ------- |
| `-s, --scheme` | Numbering scheme: `imgt` (`i`), `kabat` (`k`), `chothia` (`c`), `martin` (`m`), `aho` (`a`)                                          | `imgt`  |
| `-c, --chain`  | Chain filter: `h`,`k`,`l`,`a`,`b`,`g`,`d` or groups: `ig`, `tcr`, `all`. Accepts any form (`h`, `heavy`, `igh`), case-insensitive. | `ig`    |
| `-f, --format` | Output format: `tsv`, `json`, `jsonl`                                                                                              | `tsv`   |
| `--min-confidence` | Minimum alignment confidence, in [0, 1]; a sequence below it gets a record with `error_kind` `low_confidence` | `0.5`   |
| `--all-domains` | Every domain in each sequence (e.g. both domains of an scFv): one record per domain, with a 0-based `domain` column/field. With `segment`, every residue lands in exactly one domain's regions | off     |

### Input

Accepts a raw sequence, a FASTA file, or stdin (auto-detected). An argument is read as a sequence only when it is all letters, so a mistyped file name is an error rather than a sequence:

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

# Both domains of each scFv, one JSONL record per domain
immunum number --all-domains -f jsonl scfvs.fasta

# FR/CDR regions of each sequence, one TSV row per sequence
immunum segment sequences.fasta

# CDR3 of every domain of each scFv
immunum segment --all-domains scfvs.fasta | cut -f1,2,9

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
│   ├── kabat.rs     # Kabat numbering rules
│   ├── chothia.rs   # Chothia numbering rules
│   ├── martin.rs    # Martin (extended Chothia) numbering rules
│   └── aho.rs       # AHo numbering rules
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
```

### Design decisions

- **Semi-global alignment** forces full query consumption, preventing long CDR3 regions from being treated as trailing gaps.
- **Anchor positions** at highly conserved FR residues receive 3× gap penalties to stabilize alignment.
- **FR regions** use alignment-based numbering; **CDR regions** use scheme-specific insertion rules.
- Scoring matrices are generated at compile time from consensus data via `build.rs`.
