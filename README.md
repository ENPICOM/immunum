# Antinum

A high-performance CLI tool and Python library for numbering antibody and T-cell receptor sequences using standard immunological numbering schemes.

## Overview

Antinum provides both a command-line interface and Python bindings for numbering immunoglobulin (antibody) and T-cell receptor sequences according to established numbering schemes like IMGT and Kabat. The tool is built in Rust for performance and provides Python bindings via PyO3 for easy integration into bioinformatics workflows.

## Features

- **Multiple Numbering Schemes**: Support for IMGT and Kabat numbering schemes
- **Chain Type Support**: Handle all major immunoglobulin and T-cell receptor chain types:
  - Immunoglobulin Heavy (IGH), Kappa (IGK), and Lambda (IGL) chains
  - T-cell receptor Alpha (TRA), Beta (TRB), Gamma (TRG), and Delta (TRD) chains
- **Flexible Input**: Accept both file paths (FASTA/FASTQ) and direct sequence strings
- **File Format Support**: Parse FASTA and FASTQ files, including gzip-compressed files
- **Python Integration**: Use as a Python library with comprehensive error handling
- **High Performance**: Rust implementation for fast processing of large sequence datasets

## Installation

### From Source

#### Prerequisites
- Rust (1.70+)
- Python (3.11+)
- Maturin for building Python bindings

#### Build and Install
```bash
# Clone the repository
git clone https://github.com/your-username/antinum.git
cd antinum

# Build the CLI tool
cargo build --release

# Build and install Python package
maturin develop
```

## Usage

### Command Line Interface

The CLI tool is named `immunumber` and provides a simple interface for sequence numbering:

```bash
# Number a sequence string using IMGT scheme for heavy chain
immunumber -i "QVQLVQSGAEVKKPGASVKVSCKAS" -s imgt -c igh

# Process a FASTA file with multiple chain types
immunumber -i sequences.fasta -s kabat -c igh igk igl -o results.txt

# Process gzipped FASTQ file
immunumber -i sequences.fastq.gz -s imgt -c tra trb
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

Import and use the `antinum` module in your Python scripts:

```python
import antinum

# Number a single sequence
sequence = "QVQLVQSGAEVKKPGASVKVSCKAS"
result = antinum.number_sequence(sequence, "imgt", ["igh"])
print(result)

# Process multiple sequences with different parameters
sequences = [
    "QVQLVQSGAEVKKPGASVKVSCKAS",
    "DIQMTQSPSSLSASVGDRVTITC"
]

for seq in sequences:
    # Use case-insensitive scheme and chain names
    result = antinum.number_sequence(seq, "IMGT", ["IGH", "IGK"])
    print(result)

# Handle errors gracefully
try:
    result = antinum.number_sequence("INVALID", "invalid_scheme", ["igh"])
except ValueError as e:
    print(f"Error: {e}")
```
### Development
The project uses `cargo` for Rust package management and `uv` for Python package management.

### Running Tests

#### Rust Tests
```bash
cargo test
```

#### Python Tests
```bash
uv run pytest tests/
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

