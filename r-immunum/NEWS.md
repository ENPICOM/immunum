# immunum 1.0.0

* Initial r-universe release. R bindings for the `immunum` Rust crate via
  `extendr`.
* `immunum_version()` returns the version of the linked Rust core. The
  shim crate uses a path dependency to the repo-root `immunum` crate,
  so the R package always builds against the exact same Rust source at
  the same commit.
* `Annotator` R6 class with `$new()` / `$number()` / `$segment()` methods,
  mirroring the Python wrapper one-to-one.
* Vectorized polars batch path: `polars_number()`, `polars_segment()` and
  the matching `*_method()` variants that reuse a pre-built `Annotator`,
  routed through r-polars `$map_batches()` over a rayon-parallel Rust
  primitive.
* Validation support: `validation_fixtures()` and `benchmark_threshold()`
  provide the manifest and accuracy thresholds from `BENCHMARKS.toml`.
  Validation fixtures live in the repo tree at `fixtures/validation/`
  (same as the Python tests) and are not shipped inside the package.
