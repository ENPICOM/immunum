## Phase 1: Core Data Structure Changes

### 1.1 Replace ndarray with native Rust types

- __File__: `immunum/src/numbering_scheme_type.rs`

  - Replace `Array2<f64>` with `Vec<Vec<f64>>` or `ScoringMatrix` struct
  - Update all matrix access patterns from `matrix[[row, col]]` to `matrix[row][col]`
  - Implement helper methods for matrix operations

### 1.2 Create new ScoringMatrix type

- __New file__: `immunum/src/scoring_matrix.rs`

  - Define `ScoringMatrix` struct with `data: Vec<f64>`, `rows: usize`, `cols: usize`
  - Implement indexing methods: `get(row, col)`, `shape()`, etc.
  - Add serialization support if needed (replace .npy functionality)

### 1.3 Embed consensus data at compile time

- __File__: `immunum/src/constants.rs`

  - Add static consensus data using `include_str!` for each .txt file
  - Create const functions to parse embedded consensus strings
  - Generate static hashmaps for consensus amino acids

## Phase 2: Scoring Matrix Calculation Refactor

### 2.1 Move matrix calculation to constructor

- __File__: `immunum/src/consensus_scoring.rs`

  - Remove file I/O functions: `read_consensus_file()`, `read_scoring_matrix()`
  - Remove `write_scoring_matrix()` and `write_all_scoring_matrices()`
  - Create `calculate_scoring_matrix(consensus, params)` function
  - Keep BLOSUM62 lookup and scoring logic

### 2.2 Update NumberingScheme construction

- __File__: `immunum/src/numbering_scheme_type.rs`

  - Add constructor: `NumberingScheme::new(scheme_type, chain_type, params)`
  - Calculate scoring matrix in constructor using embedded consensus data
  - Remove `scoring_matrix` field initialization from file reading

### 2.3 Update scheme factory functions

- __File__: `immunum/src/schemes.rs`

  - Replace current functions with new constructors
  - Remove file reading calls
  - Add parameter-based variants: `imgt_heavy_with_params(params)`

## Phase 3: API Updates

### 3.1 Core Rust API

- __File__: `immunum/src/lib.rs`

  - Export new constructor functions
  - Update public API documentation
  - Add convenience methods for common use cases

### 3.2 Python bindings

- __File__: `immunum/src/python_bindings.rs`

  - Add `#[pyclass]` for `NumberingScheme` and `ScoringParams`
  - Implement `__new__` methods for Python constructors
  - Add batch processing methods: `number_sequences(sequences: Vec<String>)`
  - Update existing Python API to use new constructors

### 3.3 WASM bindings

- __File__: `immunum/src/wasm_bindings.rs`

  - Add `#[wasm_bindgen]` for new constructor methods
  - Implement JavaScript-friendly parameter passing
  - Add batch processing methods for arrays
  - Update memory management for scheme objects

## Phase 4: Dependency and Build Changes

### 4.1 Update Cargo.toml

- __File__: `immunum/Cargo.toml`

  - Remove `ndarray = "0.16.0"`
  - Remove `ndarray-npy = "0.9.1"`
  - Verify no other dependencies need ndarray

### 4.2 Update workspace Cargo.toml if needed

- __File__: `Cargo.toml`
  - Check if any workspace-level changes are needed

## Phase 5: Testing and Validation

### 5.1 Update unit tests

- __Files__: All test modules in `src/` files

  - Update tests that use `Array2<f64>` to use new matrix type
  - Update tests that read .npy files
  - Add tests for new constructor patterns
  - Verify scoring matrix calculations match previous results

### 5.2 Integration tests

- __File__: `tests/test_number_sequence.py`

  - Update Python tests to use new API
  - Add tests for custom parameters
  - Test batch processing methods

### 5.3 WASM tests

- __File__: `test-wasm-node.js`

  - Update JavaScript tests for new API
  - Test parameter passing and batch methods

## Phase 6: Documentation and Cleanup

### 6.1 Remove obsolete files

- Delete all `.npy` files from `immunum/resources/consensus/`
- Keep `.txt` files as they're now embedded at compile time
- Update `.gitignore` if needed

### 6.2 Update documentation

- __Files__: `docs/PYTHON_README.md`, `docs/WASM_README.md`

  - Document new API patterns
  - Add examples of parameter customization
  - Explain performance characteristics

### 6.3 Update examples

- __Files__: `immunum/examples/`

  - Update any examples to use new API
  - Add examples showing parameter customization

## Phase 7: Performance Validation

### 7.1 Benchmarking

- Create benchmarks comparing old vs new performance
- Verify matrix calculation is only done once per scheme
- Test batch processing performance

### 7.2 Memory usage validation

- Verify memory usage is reasonable
- Test with multiple scheme instances
- Validate WASM memory management
