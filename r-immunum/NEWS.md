# immunum 1.1.0

* R bindings for the `immunum` Rust crate via `extendr`. Depends on the
  published `immunum` 1.1.0 release on crates.io — the package is
  self-contained and does not require the parent repository at install
  time.
* `immunum_version()` returns the version of the linked Rust core.
* `Annotator` R6 class with `$new()` / `$number()` / `$segment()` methods,
  mirroring the Python wrapper one-to-one.
* Validation support: `validation_fixtures()` and `benchmark_threshold()`
  provide the manifest and accuracy thresholds from `BENCHMARKS.toml`.
  Validation fixtures live in the repo tree at `fixtures/validation/`
  (same as the Python tests) and are not shipped inside the package.
