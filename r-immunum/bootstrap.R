# bootstrap.R — vendor the parent immunum Rust crate into the R package tree.
#
# r-universe runs this script automatically before `R CMD build` when it
# detects a bootstrap.R in the package root. This is the same pattern used
# by xgboost (https://github.com/dmlc/xgboost/blob/master/R-package/bootstrap.R).
#
# Why: `R CMD build` copies only the package subdir into the source tarball.
# The Cargo.toml path dependency `path = "../../.."` would break at install
# time because the parent crate is not in the tarball. This script copies the
# necessary Rust source into the package and rewrites the path before the
# tarball is created.

VENDOR_DIR <- file.path("src", "extendr", "immunum-source")

if (dir.exists(VENDOR_DIR)) {
  message("* bootstrap.R: immunum source already vendored")
  q(status = 0)
}

# The parent crate should be available (r-universe clones the full repo)
PARENT <- file.path("..", "..", "..")
CARGO_TOML <- file.path(PARENT, "Cargo.toml")

if (!file.exists(CARGO_TOML)) {
  # Try one level up (subdir = "r-immunum" means parent is "..")
  PARENT <- ".."
  CARGO_TOML <- file.path(PARENT, "Cargo.toml")
}

if (!file.exists(CARGO_TOML)) {
  stop("bootstrap.R: cannot find parent immunum crate")
}

message("* bootstrap.R: vendoring immunum source from ", normalizePath(PARENT))
dir.create(VENDOR_DIR, recursive = TRUE)

# Copy only what cargo needs: Cargo.toml, build.rs, src/, resources/consensus/
file.copy(file.path(PARENT, "Cargo.toml"), VENDOR_DIR)
file.copy(file.path(PARENT, "build.rs"), VENDOR_DIR)
file.copy(file.path(PARENT, "src"), VENDOR_DIR, recursive = TRUE)
dir.create(file.path(VENDOR_DIR, "resources"), showWarnings = FALSE)
file.copy(file.path(PARENT, "resources", "consensus"),
          file.path(VENDOR_DIR, "resources"), recursive = TRUE)

# Rewrite Cargo.toml path dependency
cargo_toml <- file.path("src", "extendr", "Cargo.toml")
lines <- readLines(cargo_toml)
lines <- gsub(
  'path = "../../.."',
  'path = "immunum-source"',
  lines,
  fixed = TRUE
)
writeLines(lines, cargo_toml)

message("* bootstrap.R: done")
