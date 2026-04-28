.fixtures_dir <- function() {
  src <- file.path(testthat::test_path(), "..", "..", "..", "fixtures", "validation")
  if (dir.exists(src)) return(normalizePath(src))
  env <- Sys.getenv("IMMUNUM_FIXTURES", "")
  if (nzchar(env) && dir.exists(env)) return(normalizePath(env))
  normalizePath(src, mustWork = FALSE)
}

skip_if_no_fixtures <- function() {
  testthat::skip_if(
    !dir.exists(.fixtures_dir()),
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
    expect_true(file.exists(.fixture_path(stem)),
                info = sprintf("missing fixture %s", stem))
  }
})

# ── Accuracy validation ────────────────────────────────────────────────────

.compare_fixture <- function(path, chains, scheme) {
  df   <- utils::read.csv(path, colClasses = "character",
                         na.strings = c("", "NA"), check.names = FALSE)
  meta <- c("header", "sequence", "species")
  pos_cols <- setdiff(names(df), meta)

  rows <- immunum:::numbering_batch(df$sequence, chains, scheme,
                                    min_confidence = 0.0)

  n <- nrow(df)
  mismatches <- 0L
  for (i in seq_len(n)) {
    # Expected: non-empty cells in position columns
    mask <- !is.na(df[i, pos_cols]) & nzchar(df[i, pos_cols])
    exp_pos <- pos_cols[mask]
    exp_res <- unlist(df[i, pos_cols[mask]], use.names = FALSE)

    got_pos <- rows$positions[[i]]
    got_res <- rows$residues[[i]]

    if (is.null(got_pos) || length(got_pos) == 0L) {
      mismatches <- mismatches + 1L
      next
    }

    ord <- sort(exp_pos)
    got <- stats::setNames(got_res, got_pos)

    if (!identical(sort(got_pos), ord) ||
        !identical(unname(got[ord]), unname(exp_res[match(ord, exp_pos)]))) {
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
  fx <- validation_fixtures()
  for (i in seq_len(nrow(fx))) {
    row       <- fx[i, ]
    out       <- .compare_fixture(.fixture_path(row$stem), row$chains[[1]], row$scheme)
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
