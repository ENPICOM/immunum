![Immunum Logo](https://raw.githubusercontent.com/ENPICOM/immunum/master/docs/assets/immunum_logotype.svg)

<!-- The documentation site's pages include the sections between the --8<-- markers (see docs/index.md and docs/getting-started.md), so edit them here. Links in them are absolute so they work on GitHub, PyPI, crates.io, npm and the site. -->

<!-- --8<-- [start:intro] -->
Immunum is a high-performance antibody and T-cell receptor sequence numbering tool for Rust, Python, Polars,
JavaScript and the command line. [Try it in your browser](https://immunum.enpicom.com/demo/).

Documentation: [immunum.enpicom.com](https://immunum.enpicom.com).

[![Crates.io](https://img.shields.io/crates/v/immunum)](https://crates.io/crates/immunum)
[![PyPI](https://img.shields.io/pypi/v/immunum)](https://pypi.org/project/immunum/)
[![npm](https://img.shields.io/npm/v/immunum)](https://www.npmjs.com/package/immunum)
[![License: MIT](https://img.shields.io/crates/l/immunum)](https://github.com/ENPICOM/immunum/blob/main/LICENSE)
[![CI](https://img.shields.io/github/actions/workflow/status/ENPICOM/immunum/ci.yml?label=CI)](https://github.com/ENPICOM/immunum/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-immunum.enpicom.com-blue)](https://immunum.enpicom.com)
<!-- --8<-- [end:intro] -->

<!-- --8<-- [start:overview] -->
## Overview

`immunum` numbers antibody and T-cell receptor variable domain sequences and splits them
into framework and complementarity-determining regions. It uses Needleman-Wunsch
semi-global alignment against position-specific scoring matrices built from consensus sequences,
with BLOSUM62-based substitution scores.

| Antibody     | TCR         |
| ------------ | ----------- |
| IGH (heavy)  | TRA (alpha) |
| IGK (kappa)  | TRB (beta)  |
| IGL (lambda) | TRD (delta) |
|              | TRG (gamma) |

The chain is detected by aligning against every chain you allow and keeping the best match.
Positions follow one of five schemes, all derived from the internal IMGT numbering:

- **IMGT**: all 7 chains
- **Kabat**, **Chothia**, **Martin** (extended Chothia) and **AHo**: antibody chains (IGH, IGK, IGL)

How the schemes and their CDR definitions differ is explained in
[Concepts](https://immunum.enpicom.com/concepts/).
<!-- --8<-- [end:overview] -->

<!-- --8<-- [start:install] -->
## Installation

### Python

```bash
pip install immunum
pip install "immunum[polars]"   # with Polars, for immunum.polars
```

Python 3.10 or newer. Prebuilt wheels cover Linux (x86_64 and aarch64, manylinux), macOS (Apple
Silicon and Intel) and Windows (x64); elsewhere `pip` builds from source, which needs
[Rust](https://rustup.rs). `immunum.polars` needs Polars, which plain `pip install immunum` doesn't
install.

### JavaScript

```bash
npm install immunum
```

The npm package is built for Node.js. It ships its TypeScript types.

### Rust

```bash
cargo add immunum
cargo add immunum --no-default-features   # the library without the CLI's dependencies
```

Rust 1.89 or newer. API reference: [docs.rs/immunum](https://docs.rs/immunum).

### Command line

The Python package installs the `immunum` command, so `uvx` runs it without installing anything
first, and `pip install immunum` puts it on your path:

```bash
uvx immunum number sequences.fasta
```

Or build it with Cargo (Rust 1.89 or newer, from [rustup](https://rustup.rs)):

```bash
cargo install immunum
```

There are no standalone prebuilt binaries. `immunum --version` prints the version you're running.
<!-- --8<-- [end:install] -->

<!-- --8<-- [start:quickstart] -->
## Quickstart

### Python

```python
from immunum import Annotator

annotator = Annotator(chains=["ig"], scheme="imgt")
sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"

result = annotator.number(sequence)
print(result.chain)           # H
print(result.numbering["1"])  # Q

segments = annotator.segment(sequence)
print(segments.cdr3)          # AREGTTGKPIGAFAH
```

More in the [Python guide](https://immunum.enpicom.com/usage/python/).

### Polars

```python
import polars as pl
import immunum.polars as imp

df = pl.DataFrame({"sequence": [
    "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS",
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
]})
cdr3 = df.select(imp.segment("sequence", chains=["ig"], scheme="imgt").struct.field("cdr3"))
print(cdr3["cdr3"].to_list())  # ['AREGTTGKPIGAFAH', 'QQHYTTPPT']
```

More in the [Polars guide](https://immunum.enpicom.com/usage/polars/).

### JavaScript

```js
const { Annotator } = require("immunum");

const annotator = new Annotator(["ig"], "imgt");
const result = annotator.number(
  "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS",
);
console.log(result.chain);              // "H"
console.log(result.numbering.get("1")); // "Q"
annotator.free();
```

More in the [JavaScript guide](https://immunum.enpicom.com/usage/javascript/).

### Rust

```rust
use immunum::{Annotator, Chain, Scheme};

let annotator = Annotator::new(&[Chain::IGH, Chain::IGK, Chain::IGL], Scheme::IMGT, None)?;
let sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";

let segments = annotator.segment(sequence)?;
println!("CDR3: {}", segments.cdr3); // AREGTTGKPIGAFAH
```

More in the [Rust guide](https://immunum.enpicom.com/usage/rust/) and on
[docs.rs](https://docs.rs/immunum).

### Command line

```bash
$ immunum number sequences.fasta | head -3
sequence_id	chain	scheme	confidence	position	residue	error	error_kind
seq_1	H	IMGT	0.7845	1	Q
seq_1	H	IMGT	0.7845	2	V
$ immunum segment -s kabat -f jsonl sequences.fasta
```

More in the [command-line guide](https://immunum.enpicom.com/usage/cli/).

### Errors

A setup mistake, such as an unknown chain name, is raised; a sequence that can't be numbered comes
back as a result with `error` and `error_kind` set, so a batch never stops for one bad sequence. See
[Errors](https://immunum.enpicom.com/usage/errors/).
<!-- --8<-- [end:quickstart] -->

## Documentation

- [Getting started](https://immunum.enpicom.com/getting-started/): installing, and which interface to use
- [Concepts and FAQ](https://immunum.enpicom.com/concepts/): chains, schemes, confidence, prefix and postfix
- [Usage guides](https://immunum.enpicom.com/usage/) and the [API reference](https://immunum.enpicom.com/api/)
- [Upgrading to 2.0](https://immunum.enpicom.com/usage/upgrading/) and the [changelog](CHANGELOG.md)
- [Benchmarks](https://immunum.enpicom.com/benchmarks/)

## Citing immunum

If you use immunum in your research, please cite it: [CITATION.cff](CITATION.cff) has the details,
and GitHub's "Cite this repository" button formats them.

## Contributing

Development setup, tests and the project layout are in [CONTRIBUTING.md](CONTRIBUTING.md).

## License

MIT, see [LICENSE](LICENSE).
