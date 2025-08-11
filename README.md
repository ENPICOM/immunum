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

#### Building from Source

To build and install the Python library from source:

```bash
# Clone the repository
git clone https://github.com/ENPICOM/immunum.git && cd immunum

# Install UV (Python package manager) if not already installed
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create and activate virtual environment
uv venv && source .venv/bin/activate

# Install dependencies
uv sync

# Build and install the Python package in development mode
uv run maturin develop --features python
```

#### Installing in Your Python Environment

To install the library in your own Python project:

```bash
# Option 1: Build wheel and install
uv run maturin build --features python --release
pip install target/wheels/immunum-*.whl

# Option 2: Install directly in development mode (for contributors)
uv run maturin develop --features python
```

#### Python API Usage

The Python API provides a clean, object-oriented interface for sequence numbering:

```python
import immunum

# Create an annotator for IMGT scheme with heavy chain support
annotator = immunum.Annotator(
    scheme=immunum.Scheme.IMGT,
    chains=[immunum.Chain.IGH]
)

# Number a single sequence
sequence = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS"
result = annotator.number_sequence(sequence)

print(f"Chain: {result.chain}")
print(f"Scheme: {result.scheme}")
print(f"Identity: {result.identity:.2f}")
print(f"Numbers: {result.numbers[:10]}...")  # First 10 positions

# Access specific regions
cdr_sequences = result.get_cdr_sequences()
print(f"CDR1: {cdr_sequences.get('CDR1', 'Not found')}")
print(f"CDR2: {cdr_sequences.get('CDR2', 'Not found')}")
print(f"CDR3: {cdr_sequences.get('CDR3', 'Not found')}")

# Multi-chain annotator for light and heavy chains
multi_annotator = immunum.Annotator(
    scheme=immunum.Scheme.IMGT,
    chains=[immunum.Chain.IGH, immunum.Chain.IGK, immunum.Chain.IGL]
)

# Process multiple sequences
sequences = [
    "QVQLVQSGAEVKKPGASVKVSCKAS",  # Heavy chain fragment
    "DIQMTQSPSSLSASVGDRVTITC"     # Light chain fragment
]

results = multi_annotator.number_sequences(sequences)
for i, result in enumerate(results):
    print(f"Sequence {i+1}: Chain {result.chain}, Identity {result.identity:.2f}")

# Use parallel processing for better performance with many sequences
parallel_results = multi_annotator.number_sequences(sequences, parallel=True)
print(f"Processed {len(parallel_results)} sequences in parallel")

# Process sequences from a FASTA file
file_results = multi_annotator.number_file("sequences.fasta")
for seq_name, result in file_results:
    print(f"{seq_name}: Chain {result.chain}, Identity {result.identity:.2f}")

# Use parallel processing for better performance with large files
parallel_file_results = multi_annotator.number_file("large_sequences.fasta", parallel=True)
print(f"Processed {len(parallel_file_results)} sequences in parallel")

# Custom scoring parameters
custom_params = immunum.ScoringParams(
    gap_pen_cp=60.0,  # Gap penalty for conserved positions
    gap_pen_fr=30.0,  # Gap penalty for framework regions
    cdr_increase=0.7  # CDR increase factor
)

custom_annotator = immunum.Annotator(
    scheme=immunum.Scheme.KABAT,
    chains=[immunum.Chain.IGH],
    scoring_params=custom_params
)

# Prefiltering is enabled by default for better performance
# You can disable it if needed
no_prefilter_annotator = immunum.Annotator(
    scheme=immunum.Scheme.IMGT,
    chains=[immunum.Chain.IGH, immunum.Chain.IGK, immunum.Chain.IGL, 
            immunum.Chain.TRA, immunum.Chain.TRB],
    disable_prefiltering=True
)

# Error handling
try:
    result = annotator.number_sequence("INVALID_SHORT_SEQUENCE")
except RuntimeError as e:
    print(f"Numbering failed: {e}")
```

#### Key Classes and Methods

- **`Annotator`**: Main class for sequence numbering
  - `number_sequence(sequence)`: Number a single sequence
  - `number_sequences(sequences, parallel=False)`: Number multiple sequences with optional parallel processing
  - `number_file(file_path, parallel=False)`: Process FASTA/FASTQ files with optional parallel processing

- **`AnnotationResult`**: Contains numbering results
  - Properties: `sequence`, `numbers`, `scheme`, `chain`, `identity`, `regions`
  - Methods: `get_region_sequence()`, `get_cdr_sequences()`, `get_framework_sequences()`

- **`ScoringParams`**: Customizable scoring parameters for alignment
  - All parameters have getters/setters for fine-tuning alignment behavior

- **Enums**: `Scheme.IMGT`/`Scheme.KABAT`, `Chain.IGH`/`Chain.IGK`/`Chain.IGL`/etc.

#### Recent API Changes

**v0.2.0+**: The Python API has been streamlined and enhanced:
- **Parallel Processing**: Added `rayon`-based multithreading support for `number_sequences()` and `number_file()` methods
- **Streamlined API**: Removed `PyScoringParams` wrapper - now uses `ScoringParams` directly with automatic property access
- **Consistent naming**: Unified class naming between Python and WASM APIs
- **Enhanced Properties**: All `ScoringParams` properties support direct getting/setting via PyO3 attributes
- **Better Type Support**: Improved type annotations in `immunum.pyi` for better IDE support

#### Performance Tips

- **Use parallel processing** (`parallel=True`) when processing multiple sequences or large files
- **Prefiltering is enabled by default** for optimal performance with multiple chain types
- **Batch processing** is more efficient than processing sequences individually
- **File processing** is optimized for both FASTA and FASTQ formats, including gzip compression

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

