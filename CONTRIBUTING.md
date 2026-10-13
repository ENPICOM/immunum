# Contributing

Bug reports, questions and pull requests are welcome on
[GitHub](https://github.com/ENPICOM/immunum/issues). For a numbering you think is wrong, include
the sequence, the chains and scheme you asked for, and the numbering you expected.

## Development

To orchestrate a project between cargo and python, we use [`task`](https://taskfile.dev).
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
task test-python  # test only python code, including the Python examples in README.md and docs/
task test-wasm    # test only the JavaScript bindings
task test         # test all code
task ci           # everything CI runs, including a strict docs build
```

### Linting

```bash
task format  # formats python code and checks rust formatting
task lint    # runs linting for python and rust
```

### Documentation

The site at [immunum.enpicom.com](https://immunum.enpicom.com) is built with MkDocs Material from
`docs/` and `mkdocs.yml`:

```bash
task docs        # serve it locally, rebuilding on changes
task docs-build  # strict build into site/, which fails on broken links
```

- The home and getting-started pages include sections of `README.md`, marked with
  `<!-- --8<-- [start:name] -->`; edit those in the README. The changelog page includes
  `CHANGELOG.md` and this page includes `CONTRIBUTING.md`.
- `scripts/mkdocs_pages.py` fills in the TypeScript declarations on the JavaScript API page from
  the web WASM build, and the accuracy table on the benchmarks page from `BENCHMARKS.toml`.
- Python code blocks in `README.md` and `docs/` run as tests (`pytest-markdown-docs`), so their
  output comments must be real.

### Benchmarking

The validation tools in `src/bin` read the fixtures in a checkout, so they're behind the `dev-tools` feature and not installed with the CLI. `benchmark` writes the per-scheme accuracy in `BENCHMARKS.toml`; `debug-validation <CHAIN> [SCHEME] [HEADER]` shows the alignments behind each mismatch:

```bash
cargo run --release --features dev-tools --bin benchmark
cargo run --release --features dev-tools --bin debug-validation -- H kabat
```

There are multiple benchmarks in the repository. For full list, see `task | grep benchmark`:

```bash
$ task | grep benchmark
* bencher-python:               Run Python benchmarks with bencher. Alerts are disabled, used only for tracking performance
* bencher-rust:                 Run Criterion benchmarks with bencher
* benchmark:                    Run all benchmarks and produce plots for them
* benchmark-accuracy:           Accuracy benchmark across all fixtures (1k sequences, 7 rounds each)
* benchmark-cli:                Speed benchmark using Criterion (benches/speed_benchmark.rs)
* benchmark-comparison:         Speed + correctness benchmark: immunum vs antpack vs anarci (1k IGH sequences)
* benchmark-speed:              Speed benchmark: sizes 100..1M, 3 rounds, H/imgt. Pass CLI_ARGS to filter tools, e.g. -- --tools immunum
* benchmark-speed-polars:       Speed benchmark for immunum polars across all chain/scheme fixtures
* plots:                        Produce plots from benchmark results
```

## Project structure

```
src/
├── main.rs          # CLI binary: runs cli.rs
├── cli.rs           # The immunum command (number / segment), also installed by the Python package
├── lib.rs           # Public API
├── annotator.rs     # Sequence annotation, chain detection and domain search
├── alignment.rs     # Needleman-Wunsch semi-global alignment
├── error.rs         # Error (raised) and SequenceError (returned)
├── io.rs            # Input parsing (FASTA, raw) and output formatting (TSV, JSON, JSONL)
├── numbering.rs     # Numbering module entry point
├── numbering/
│   ├── imgt.rs      # IMGT numbering rules
│   ├── kabat.rs     # Kabat numbering rules
│   ├── chothia.rs   # Chothia numbering rules
│   ├── martin.rs    # Martin (extended Chothia) numbering rules
│   └── aho.rs       # AHo numbering rules
├── scoring.rs       # PSSM and scoring matrices
├── types.rs         # Core domain types (Chain, Scheme, Position, Region)
├── python.rs        # Python bindings (pyo3)
├── polars.rs        # Polars plugin expressions
├── wasm.rs          # JavaScript bindings (wasm-bindgen) and TypeScript types
├── validation.rs    # Validation utilities
└── bin/
    ├── benchmark.rs       # Validation metrics report
    └── debug_validation.rs # Alignment mismatch visualization
benches/             # Criterion and Python speed/correctness benchmarks
resources/
├── consensus/         # Consensus sequence CSVs (compiled into scoring matrices)
└── benchmark_results/ # Benchmark CSVs the plots on the benchmarks page are made from
fixtures/
├── validation/      # ANARCI-numbered reference datasets
├── ig.fasta         # Example antibody sequences
└── ig.tsv           # Example TSV input
docs/                # Documentation site pages, the demo's script and styles
scripts/             # Python tooling for consensus data, benchmark plots and the docs build
tests/               # CLI, Python, Polars and JS tests, and the shared error_cases.json
immunum/
├── __init__.py      # Python module
├── _internal.pyi    # Python stub file for pyo3
└── polars.py        # Polars extension module
```

### Design decisions

- **Semi-global alignment** forces full query consumption, preventing long CDR3 regions from being treated as trailing gaps.
- **Anchor positions** at highly conserved FR residues receive 3× gap penalties to stabilize alignment.
- **FR regions** use alignment-based numbering; **CDR regions** use scheme-specific insertion rules.
- Scoring matrices are generated at compile time from consensus data via `build.rs`.
