## Development Guidelines

### Python
- Always execute python commands with `uv` and the `.venv` environment

### Rust
- **Test-Driven Development**: Write tests before implementation
  - Don't create nonsensical tests
  - Use `cargo test` to run all tests
  - Use `cargo test --lib <module>::tests::<test_name>` for specific tests

- **Code Quality**: Run before committing
  - `cargo fmt` - Auto-format code
  - `cargo clippy` - Lint and catch common mistakes
  - `cargo test` - Verify all tests pass

- **Validation Benchmarks**: Track accuracy metrics
  - `cargo run --quiet --release --bin benchmark > BENCHMARKS.md` - Update validation metrics
  - Always use `--release` for accurate execution time measurements
  - Commit updated BENCHMARKS.md when making improvements to alignment/scoring
  - See BENCHMARKS.md for current accuracy across all chains

- **Best Practices**:
  - Only create the code needed in the moment
  - Keep the code concise
  - Use `cargo build --release` for optimized builds

## Project Structure

### Core Library (`src/`)
- **lib.rs**: Library entry point, exports public API
- **types.rs**: Core data types (Chain, Position, Region, Scheme)
- **error.rs**: Error types using thiserror
- **scoring.rs**: Position-specific scoring matrices loaded from build-time generated JSON
- **alignment.rs**: Needleman-Wunsch semi-global alignment with IMGT-aware region numbering
- **imgt.rs**: IMGT numbering rules for CDR1/2/3 with palindromic insertions
- **annotator.rs**: High-level API for sequence annotation with auto chain detection
- **consensus.rs**: Consensus sequence parsing and representation
- **validation.rs**: Validation framework comparing numbering against test datasets

### Binaries (`src/bin/`, `src/main.rs`)
- **main.rs**: Example CLI demonstrating annotator usage
- **debug_validation.rs**: Debug tool for visualizing alignment mismatches with expected vs actual positions
- **benchmark.rs**: Generate validation metrics report (outputs to BENCHMARKS.md)

### Build System
- **build.rs**: Compile-time generation of scoring matrices from TSV consensus files with region-aware gap/insertion penalties
  - Penalty constants at top of file for easy tuning
  - Uses BLOSUM62 for amino acid substitution scores

### Key Design Decisions
- Semi-global alignment forces full query consumption to prevent long CDR3s being treated as trailing gaps
- FR regions use alignment-based numbering, CDR regions use scheme-specific numbering rules
- Highly conserved FR positions have 3x gap penalties to enforce anchor points

## Data Files

### Fixtures (`fixtures/validation/`)
Validation datasets with expected IMGT numbering for testing:
- **ab_H_imgt.csv, ab_K_imgt.csv, ab_L_imgt.csv**: Antibody heavy/kappa/lambda chains (2473/1491/383 sequences)
- **tcr_A_imgt.csv, tcr_B_imgt.csv, tcr_G_imgt.csv, tcr_D_imgt.csv**: TCR alpha/beta/gamma/delta chains (865/934/25/23 sequences)

### Resources (`resources/consensus/`)
Consensus sequences with amino acid frequencies, conservation scores, and penalties:
- **ab_imgt_{H,K,L}_consensus.tsv**: IMGT antibody consensus sequences
- **tcr_imgt_{A,B,D,G}_consensus.tsv**: IMGT TCR consensus sequences
- Format: position, amino acids, frequencies, occupancy, gaps, penalties, conservation class, region, allows_insertion

## Python Scripts (`scripts/`)

### Data Processing
- **number_test_sequences.py**: Process FASTA files with AntPack/ANARCI, generate numbered CSV validation files
- **generate_consensus.py**: Generate consensus TSV files from numbered sequences with conservation metrics

### Analysis
- **analyze_sequences.py**: Analyze numbered sequences, generate statistics and plots for CDR/FR regions
- **imgt_kabat_conversion.py**: Convert between IMGT and Kabat CDR3 numbering schemes using sequence indices

### Utilities (`scripts/utils/`)
- **fasta.py**: FASTA file parsing
- **numbering.py**: Numbering utilities and position definitions
- **csv_output.py**: CSV file writing utilities
