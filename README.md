# Immunum

A high-performance CLI tool, Python library, and WebAssembly module for numbering antibody and T-cell receptor sequences using standard immunological numbering schemes.

## Overview

Immunum provides a command-line interface, Python bindings and WebAssembly module for numbering immunoglobulin (antibody) and T-cell receptor sequences according to established numbering schemes like IMGT and Kabat. The tool is built in Rust for performance and provides Python bindings via PyO3 for easy integration into bioinformatics workflows. It also provides a WebAssembly module for use in web applications.

## Features

- **Multiple Numbering Schemes**: Support for IMGT and Kabat numbering schemes
- **Chain Type Support**: Handle all major immunoglobulin and T-cell receptor chain types:
  - Immunoglobulin Heavy (IGH), Kappa (IGK), and Lambda (IGL) chains
  - T-cell receptor Alpha (TRA), Beta (TRB), Gamma (TRG), and Delta (TRD) chains
- **Flexible Input**: Accept both file paths (FASTA/FASTQ) and direct sequence strings
- **File Format Support**: Parse FASTA and FASTQ files, including gzip-compressed files
- **High Performance**: Rust implementation for fast processing of large sequence datasets
- **Output formats**: #TODO

## Installation

### From Source

#### Prerequisites
- Rust (1.70+)
- Python (3.11+)

#### Build and Install
```bash
# Clone the repository
git clone https://github.com/ENPICOM/immunum.git && cd immunum

# Build the CLI tool
cargo build
```

## Usage

### Command Line Interface

The CLI tool is named `immunum-cli` and provides a simple interface for sequence numbering:

```bash
# Number a sequence string using IMGT scheme for heavy chain
immunum-cli -i "QVQLVQSGAEVKKPGASVKVSCKAS" -s imgt -c igh

# Process a FASTA file with multiple chain types
immunum-cli -i sequences.fasta -s kabat -c igh igk igl -o results.txt

# Process gzipped FASTQ file
immunum-cli -i sequences.fastq.gz -s imgt -c tra trb
```

#### CLI Options

- `-i, --input <INPUT>`: Input file path (FASTA/FASTQ) or sequence string
- `-s, --scheme <SCHEME>`: Numbering scheme (`imgt` or `kabat`) [default: imgt]
- `-c, --chains <CHAINS>`: Chain types to process (multiple values allowed)
- `-o, --output <OUTPUT>`: Output file path (optional, defaults to stdout)

#### Supported Chain Types

| Chain Type | Full Name | Aliases |
|------------|-----------|---------|
| `igh` | Immunoglobulin Heavy | `heavy`, `h` |
| `igk` | Immunoglobulin Kappa | `kappa`, `k` |
| `igl` | Immunoglobulin Lambda | `lambda`, `l` |
| `tra` | T-cell Receptor Alpha | `alpha`, `a` |
| `trb` | T-cell Receptor Beta | `beta`, `b` |
| `trg` | T-cell Receptor Gamma | `gamma`, `g` |
| `trd` | T-cell Receptor Delta | `delta`, `d` |

#### Supported Numbering Schemes

| Scheme | Description | Aliases |
|--------|-------------|---------|
| `imgt` | IMGT numbering scheme | `i` |
| `kabat` | Kabat numbering scheme | `k` |

### Python Library

Build and install the Python package:
```bash
# Install the virtual environment
uv venv && uv sync && source .venv/bin/activate

# Install the Python package in this environment
maturin develop --features python

# OR build the Python package
maturin build --features python
# And install the Python package in your own project
uv add target/wheels/immunum-*.whl
```

Import and use the `immunum` module in your Python scripts:

```python
import immunum

# Number a single sequence
sequence = "QVQLVQSGAEVKKPGASVKVSCKAS"
result = immunum.number_sequence(sequence, "imgt", ["igh"])
print(result)

# Process multiple sequences with different parameters
sequences = [
    "QVQLVQSGAEVKKPGASVKVSCKAS",
    "DIQMTQSPSSLSASVGDRVTITC"
]

for seq in sequences:
    # Use case-insensitive scheme and chain names
    result = immunum.number_sequence(seq, "IMGT", ["IGH", "IGK"])
    print(result)

# Handle errors gracefully
try:
    result = immunum.number_sequence("INVALID", "invalid_scheme", ["igh"])
except ValueError as e:
    print(f"Error: {e}")
```

### WASM Module

The project also provides a WASM module that can be used in JS environments. You can build the WASM module using the following command:

```bash
cd immunum && wasm-pack build --target web --out-dir ../wasm_build --features wasm --no-default-features; cd -
```
Check `test-wasm-node.js` for an example of how to use the WASM module in a Node.js environment.


### Development
The project uses `cargo` for Rust package management and `uv` for Python package management.

### Running Tests

The Rust test should cover all the functionality of the library. You can run the tests using the following command:

```bash
cargo test
```

We also have integration tests for the python and WASM bindings, these are intended to be small and only test the bindings themselves. You can run the tests using the following command:

```bash
uv run pytest
node test-wasm-node.js
```

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add tests for new functionality
5. Run the test suite (`cargo test && pytest`)
6. Commit your changes (`git commit -m 'Add amazing feature'`)
7. Push to the branch (`git push origin feature/amazing-feature`)
8. Open a Pull Request

