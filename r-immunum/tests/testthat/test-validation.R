.fixtures_dir <- function() {
  # From source tree: tests/testthat -> r-immunum -> repo root -> fixtures
  src <- file.path(testthat::test_path(), "..", "..", "..", "fixtures", "validation")
  if (dir.exists(src)) return(normalizePath(src))
  # From installed tests: try IMMUNUM_FIXTURES env var or repo working directory
  env <- Sys.getenv("IMMUNUM_FIXTURES", "")
  if (nzchar(env) && dir.exists(env)) return(normalizePath(env))
  normalizePath(src, mustWork = FALSE)
}

skip_if_no_fixtures <- function() {
  d <- .fixtures_dir()
  testthat::skip_if(
    !dir.exists(d),
    "validation fixtures not available (not in repo tree)"
  )
}

.fixture_path <- function(stem) {
  file.path(.fixtures_dir(), paste0(stem, ".csv"))
}

# ── Manifest / threshold accessors ──────────────────────────────────────────

test_that("validation_fixtures returns the manifest", {
  fx <- validation_fixtures()
  expect_s3_class(fx, "data.frame")
  expect_true(all(c("stem", "scheme", "benchmark", "chains") %in% names(fx)))
  expect_gte(nrow(fx), 10L)
  expect_type(fx$chains, "list")
  expect_type(fx$chains[[1]], "character")
})

test_that("benchmark_threshold returns numeric percentages", {
  expect_equal(benchmark_threshold("imgt.H"), 99.88)
  expect_equal(benchmark_threshold("kabat.K"), 99.66)
  expect_equal(benchmark_threshold("imgt.D"), 100.0)
})

test_that("benchmark_threshold rejects unknown keys", {
  expect_error(benchmark_threshold("nope.X"), "Unknown benchmark key")
  expect_error(benchmark_threshold(NA_character_), "non-NA")
})

# ── Fixture files exist in repo ─────────────────────────────────────────────

test_that("all manifest fixtures exist in the repo tree", {
  skip_if_no_fixtures()
  fx <- validation_fixtures()
  for (stem in fx$stem) {
    path <- .fixture_path(stem)
    expect_true(file.exists(path),
                info = sprintf("missing fixture %s", stem))
  }
})

# ── Accuracy validation ────────────────────────────────────────────────────

.compare_fixture <- function(path, chains, scheme) {
  pl <- polars::pl
  df <- pl$read_csv(path, infer_schema_length = 0L)
  meta <- c("header", "sequence", "species")
  position_cols <- setdiff(df$columns, meta)

  numbered <- polars_number(pl$col("sequence"),
                            chains = chains, scheme = scheme,
                            min_confidence = 0.0)$alias("numbered")
  result <- df$with_columns(numbered)$unnest("numbered")
  rdf <- as.data.frame(result)

  n <- nrow(rdf)
  mismatches <- 0L
  for (i in seq_len(n)) {
    expected <- character(0)
    for (p in position_cols) {
      aa <- rdf[[p]][i]
      if (!is.null(aa) && !is.na(aa) && nzchar(aa)) {
        expected[p] <- aa
      }
    }
    pos_vec <- rdf$positions[[i]]
    res_vec <- rdf$residues[[i]]
    if (is.null(pos_vec) || length(pos_vec) == 0L) {
      mismatches <- mismatches + 1L
      next
    }
    got <- stats::setNames(res_vec, pos_vec)
    if (!identical(sort(names(got)), sort(names(expected))) ||
        !identical(unname(got[names(expected)]), unname(expected))) {
      mismatches <- mismatches + 1L
    }
  }

  perfect <- n - mismatches
  list(
    mismatches  = mismatches,
    total       = n,
    perfect     = perfect,
    perfect_pct = if (n > 0L) 100 * perfect / n else 0
  )
}

test_that("fixtures match BENCHMARKS.toml thresholds", {
  skip_if_no_fixtures()
  testthat::skip_if_not_installed("polars")
  fx <- validation_fixtures()
  for (i in seq_len(nrow(fx))) {
    row <- fx[i, ]
    path <- .fixture_path(row$stem)
    out <- .compare_fixture(path, row$chains[[1]], row$scheme)
    threshold <- benchmark_threshold(row$benchmark)
    expect_gte(
      round(out$perfect_pct, 2),
      threshold,
      label = sprintf(
        "%s (%d/%d perfect, %.2f%%, threshold %.2f%%)",
        row$stem, out$perfect, out$total, out$perfect_pct, threshold
      )
    )
  }
})
