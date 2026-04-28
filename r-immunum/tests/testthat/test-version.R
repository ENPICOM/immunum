test_that("immunum_version returns the linked Rust crate version", {
  v <- immunum_version()
  expect_type(v, "character")
  expect_length(v, 1L)

  # Must look like a SemVer (major.minor.patch with optional pre-release).
  expect_match(v, "^[0-9]+\\.[0-9]+\\.[0-9]+(-[A-Za-z0-9.-]+)?$")
})

test_that("R package version matches the linked Rust crate version", {
  # The R DESCRIPTION version and the Rust shim crate version are kept in
  # lockstep. Both must match the upstream immunum crate version.
  pkg_version <- as.character(utils::packageVersion("immunum"))
  rust_version <- immunum_version()
  expect_equal(pkg_version, rust_version)
})
