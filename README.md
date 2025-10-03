# Immunum

A high-performance CLI tool, Python library, and WebAssembly module for numbering antibody and T-cell receptor sequences using standard immunological numbering schemes.

## Overview

Immunum provides a command-line interface, Python bindings and WebAssembly module for numbering immunoglobulin (antibody) and T-cell receptor sequences according to established numbering schemes like IMGT and Kabat. The tool is built in Rust for performance and provides Python bindings via PyO3 for easy integration into bioinformatics workflows. It also provides a WebAssembly module for use in web applications.

## Features

- **Multiple Numbering Schemes**: Support for IMGT and Kabat numbering schemes
- **Flexible CDR Definitions**: IMGT, Kabat, and North CDR boundary definitions
- **Chain Type Support**: Handle all major immunoglobulin and T-cell receptor chain types:
  - Immunoglobulin Heavy (IGH), Kappa (IGK), and Lambda (IGL) chains
  - T-cell receptor Alpha (TRA), Beta (TRB), Gamma (TRG), and Delta (TRD) chains  ( Not tested yet )
- **Multi-Chain Analysis**: Detect multiple chain types in single sequences (paired analysis)
- **High Performance**:
  - Rust implementation
  - Parallel processing support for large datasets
  - K-mer prefiltering for fast multi-chain analysis
  - Custom thread pool configuration
- **Flexible Input/Output**:
  - File formats: FASTA, FASTQ (including gzip compression)
  - Direct sequence strings
  - JSON output format with detailed region information
- **Quality Control**: Configurable confidence thresholds and alignment parameters
- **Multiple Interfaces**: CLI tool, Python library, and WebAssembly module

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

The CLI tool is named `immunum-cli` and provides a simple interface for sequence numbering with JSON output:

```bash
# Number a sequence string using IMGT scheme for heavy chain
immunum-cli -i "QVQLVQSGAEVKKPGASVKVSCKAS" -s imgt -c igh -o results.json

# Process a FASTA file with multiple chain types
immunum-cli -i sequences.fasta -s kabat -c igh igk igl -o results.json

# Process gzipped FASTQ file with custom CDR definitions
immunum-cli -i sequences.fastq.gz -s imgt -c tra trb --cdr-definitions north -o results.json

# Use parallel processing with custom confidence threshold
immunum-cli -i large_sequences.fasta -s imgt -c igh igk igl -t 8 -m 0.8 -o results.json

# Enable verbose output to see processing details
immunum-cli -i sequences.fasta -s kabat -c igh -v -o results.json

# Allow multiple chains per sequence (paired chain analysis)
immunum-cli -i paired_sequences.fasta -s imgt -c igh igk igl --max-chains 2 -o results.json

# Disable prefiltering for more exhaustive search
immunum-cli -i sequences.fasta -s imgt -c igh igk igl --disable-prefiltering -o results.json
```

#### CLI Options

**Required:**
- `-i, --input <INPUT>`: Input file path (FASTA/FASTQ) or sequence string
- `-o, --output <OUTPUT>`: Output file path for JSON results

**Numbering Scheme:**
- `-s, --scheme <SCHEME>`: Numbering scheme (`imgt` or `kabat`) [default: imgt]
- `--cdr-definitions <CDR_DEFINITIONS>`: CDR definition scheme (`imgt`, `kabat`, or `north`) [default: matches scheme]

**Chain Selection:**
- `-c, --chains <CHAINS>...`: Chain types to process (multiple values allowed) [default: igh igk igl]
- `--max-chains <MAX_CHAINS>`: Maximum chains per sequence (0=unlimited, 1=single) [default: 1]

**Performance & Quality:**
- `-t, --threads <THREADS>`: Number of threads for parallel processing [default: CPU cores]
- `-m, --min-confidence <MIN_CONFIDENCE>`: Minimum alignment confidence [default: 0.7]
- `--min-kmer-overlap <MIN_KMER_OVERLAP>`: K-mer overlap threshold for prefiltering [default: 0.2]
- `--disable-prefiltering`: Disable prefiltering for more exhaustive search

**Output:**
- `-v, --verbose`: Show detailed progress for each sequence processed

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

#### CDR Definition Schemes

| Scheme | Description | Best For |
|--------|-------------|----------|
| `imgt` | IMGT CDR definitions (default for IMGT scheme) | Standardized analysis |
| `kabat` | Kabat CDR definitions (default for Kabat scheme) | Traditional analysis |
| `north` | North CDR definitions | Alternative CDR boundaries |

#### Output Format

The CLI outputs results in JSON format with the following structure:

```json
[
  {
    "sequence_name": "sequence_1",
    "chains": [
      {
        "sequence": "QVQLVQSGAEVKKPGASVKVSCKAS...",
        "numbers": ["1", "2", "3", "4", "5", "..."],
        "scheme": "IMGT",
        "chain": "IGH",
        "identity": 0.95,
        "start": 0,
        "end": 119,
        "regions": {
          "fr1": {"start": 1, "end": 26, "sequence": "QVQLVQSGAEVKKPGASVKVSCKAS"},
          "cdr1": {"start": 27, "end": 38, "sequence": "GYTFTSYYMH"},
          "fr2": {"start": 39, "end": 55, "sequence": "WVRQAPGQGLEWMG"},
          "cdr2": {"start": 56, "end": 65, "sequence": "IINPSGGSTSYAQKFQ"},
          "fr3": {"start": 66, "end": 104, "sequence": "GRVTMTRDTSTSTVYMELSSLRSEDTAVYYCAR"},
          "cdr3": {"start": 105, "end": 117, "sequence": "WGGRGSYAMDYW"},
          "fr4": {"start": 118, "end": 128, "sequence": "GQGTLVTVSS"}
        }
      }
    ]
  }
]
```

**JSON Fields:**
- `sequence_name`: Original sequence identifier from FASTA header or auto-generated
- `chains`: Array of detected chains in the sequence
- `sequence`: Full input sequence
- `numbers`: Position numbers according to the numbering scheme
- `scheme`: Applied numbering scheme (IMGT/KABAT)
- `chain`: Detected chain type (IGH/IGK/IGL/TRA/TRB/TRG/TRD)
- `identity`: Alignment confidence score (0.0-1.0)
- `start`/`end`: Indices of numbered region within input sequence
- `regions`: CDR and framework regions with positions and sequences

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

The Python API provides a simple interface for sequence numbering with automatic parallel processing:

```python
import immunum

# Create an annotator for IMGT scheme with heavy chain support
annotator = immunum.Annotator(
    scheme=immunum.Scheme.IMGT,
    chains=[immunum.Chain.IGH]
)

# Process multiple sequences (returns list of chain results for each sequence)
sequences = [
    "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS",
    "DIQMTQSPSSLSASVGDRVTITC"
]

# Each sequence returns a list of detected chains: [(numbers, confidence, chain), ...]
results = annotator.number_sequences(sequences)

for i, sequence_chains in enumerate(results):
    print(f"Sequence {i+1}: Found {len(sequence_chains)} chains")
    for numbers, confidence, chain in sequence_chains:
        print(f"  Chain {chain}: {confidence:.2f} confidence")
        print(f"  Numbers: {numbers[:10]}...")  # First 10 positions

# Multi-chain annotator for comprehensive analysis
multi_annotator = immunum.Annotator(
    scheme=immunum.Scheme.IMGT,
    chains=[immunum.Chain.IGH, immunum.Chain.IGK, immunum.Chain.IGL],
    threads=8  # Use 8 threads for parallel processing
)

# Process with higher chain limit for paired sequences
paired_results = multi_annotator.number_sequences(sequences, max_chains=2)
for i, sequence_chains in enumerate(paired_results):
    print(f"Paired sequence {i+1}:")
    for numbers, confidence, chain in sequence_chains:
        print(f"  {chain}: {confidence:.2f}")

# Process sequences with IDs
sequences_with_ids = [
    ("heavy_chain_1", "QVQLVQSGAEVKKPGASVKVSCKAS"),
    ("light_chain_1", "DIQMTQSPSSLSASVGDRVTITC")
]

id_results = multi_annotator.number_sequences(sequences_with_ids)

# Custom parameters for fine-tuning
custom_annotator = immunum.Annotator(
    scheme=immunum.Scheme.KABAT,
    chains=[immunum.Chain.IGH, immunum.Chain.IGK, immunum.Chain.IGL],
    disable_prefiltering=False,  # Use prefiltering for speed
    threads=4,                   # Use 4 threads
    min_confidence=0.8,          # Higher confidence threshold
    min_kmer_overlap=0.3         # Higher k-mer overlap requirement
)


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

