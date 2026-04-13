skip_if_polars_missing <- function() {
  testthat::skip_if_not(
    requireNamespace("immunum", quietly = TRUE),
    "immunum native library not loaded"
  )
  testthat::skip_if_not_installed("polars")
}

# ── polars_number() ─────────────────────────────────────────────────────────

test_that("polars_number returns a polars expression", {
  skip_if_polars_missing()
  expr <- polars_number(polars::pl$col("sequence"),
                        chains = "IGH", scheme = "IMGT")
  expect_s3_class(expr, "polars_expr")
})

test_that("polars_number on a polars DataFrame produces a struct column", {
  skip_if_polars_missing()
  df <- polars::pl$DataFrame(sequence = IGH_SEQ)
  result <- df$select(
    polars_number(polars::pl$col("sequence"),
                  chains = "IGH", scheme = "IMGT")$alias("numbered")
  )
  expect_true("numbered" %in% result$columns)
  expect_equal(result$height, 1L)
})

test_that("polars_number struct exposes chain/scheme/positions/residues/error", {
  skip_if_polars_missing()
  df <- polars::pl$DataFrame(sequence = IGH_SEQ)
  result <- df$select(
    polars_number(polars::pl$col("sequence"),
                  chains = "IGH", scheme = "IMGT")$alias("numbered")
  )$unnest("numbered")
  cols <- result$columns
  expect_true("chain" %in% cols)
  expect_true("scheme" %in% cols)
  expect_true("positions" %in% cols)
  expect_true("residues" %in% cols)
  expect_true("error" %in% cols)
})

test_that("polars_number error field is null on success", {
  skip_if_polars_missing()
  df <- polars::pl$DataFrame(sequence = IGH_SEQ)
  result <- df$select(
    polars_number(polars::pl$col("sequence"),
                  chains = "IGH", scheme = "IMGT")$alias("numbered")
  )$unnest("numbered")
  expect_true(is.na(result$get_column("error")$to_r_vector()))
})

test_that("polars_number error field is set on failure", {
  skip_if_polars_missing()
  df <- polars::pl$DataFrame(sequence = "AAAAAAAAAAAAAAAA")
  result <- df$select(
    polars_number(polars::pl$col("sequence"),
                  chains = "IGH", scheme = "IMGT")$alias("numbered")
  )$unnest("numbered")
  expect_false(is.na(result$get_column("error")$to_r_vector()))
  expect_true(is.na(result$get_column("chain")$to_r_vector()))
})

test_that("polars_number handles a multi-row batch", {
  skip_if_polars_missing()
  df <- polars::pl$DataFrame(sequence = c(IGH_SEQ, IGL_SEQ))
  result <- df$select(
    polars_number(polars::pl$col("sequence"),
                  chains = ALL_CHAINS_R, scheme = "IMGT")$alias("numbered")
  )
  expect_equal(result$height, 2L)
})

test_that("polars_number accepts a column name shortcut", {
  skip_if_polars_missing()
  df <- polars::pl$DataFrame(sequence = IGH_SEQ)
  result <- df$select(
    polars_number("sequence", chains = "IGH", scheme = "IMGT")$alias("numbered")
  )
  expect_equal(result$height, 1L)
})

test_that("polars_number validates min_confidence range", {
  skip_if_polars_missing()
  expect_error(
    polars_number("sequence", chains = "IGH", scheme = "IMGT",
                  min_confidence = -0.1),
    "must be a single numeric in"
  )
  expect_error(
    polars_number("sequence", chains = "IGH", scheme = "IMGT",
                  min_confidence = 1.5),
    "must be a single numeric in"
  )
})

# ── polars_segment() ────────────────────────────────────────────────────────

test_that("polars_segment returns a polars expression", {
  skip_if_polars_missing()
  expr <- polars_segment(polars::pl$col("sequence"),
                         chains = "IGH", scheme = "IMGT")
  expect_s3_class(expr, "polars_expr")
})

test_that("polars_segment on a polars DataFrame produces a struct column", {
  skip_if_polars_missing()
  df <- polars::pl$DataFrame(sequence = IGH_SEQ)
  result <- df$select(
    polars_segment(polars::pl$col("sequence"),
                   chains = "IGH", scheme = "IMGT")$alias("segmented")
  )
  expect_true("segmented" %in% result$columns)
  expect_equal(result$height, 1L)
})

test_that("polars_segment struct exposes FR/CDR/prefix/postfix/error fields", {
  skip_if_polars_missing()
  df <- polars::pl$DataFrame(sequence = IGH_SEQ)
  result <- df$select(
    polars_segment(polars::pl$col("sequence"),
                   chains = "IGH", scheme = "IMGT")$alias("segmented")
  )$unnest("segmented")
  expected <- c("fr1", "fr2", "fr3", "fr4",
                "cdr1", "cdr2", "cdr3",
                "prefix", "postfix", "error")
  expect_true(all(expected %in% result$columns))
})

test_that("polars_segment error field is null on success", {
  skip_if_polars_missing()
  df <- polars::pl$DataFrame(sequence = IGH_SEQ)
  result <- df$select(
    polars_segment(polars::pl$col("sequence"),
                   chains = "IGH", scheme = "IMGT")$alias("segmented")
  )$unnest("segmented")
  expect_true(is.na(result$get_column("error")$to_r_vector()))
})

test_that("polars_segment error field is set on failure", {
  skip_if_polars_missing()
  df <- polars::pl$DataFrame(sequence = "AAAAAAAAAAAAAAAA")
  result <- df$select(
    polars_segment(polars::pl$col("sequence"),
                   chains = "IGH", scheme = "IMGT")$alias("segmented")
  )$unnest("segmented")
  expect_false(is.na(result$get_column("error")$to_r_vector()))
  expect_true(is.na(result$get_column("fr1")$to_r_vector()))
})

test_that("polars_segment handles a multi-row batch", {
  skip_if_polars_missing()
  df <- polars::pl$DataFrame(sequence = c(IGH_SEQ, IGL_SEQ))
  result <- df$select(
    polars_segment(polars::pl$col("sequence"),
                   chains = AB_CHAINS_R, scheme = "IMGT")$alias("segmented")
  )
  expect_equal(result$height, 2L)
})

# ── polars_numbering_method() / polars_segmentation_method() ────────────────

test_that("polars_numbering_method returns a polars expression", {
  skip_if_polars_missing()
  ann <- Annotator$new(chains = "IGH", scheme = "IMGT")
  expr <- polars_numbering_method(polars::pl$col("sequence"), annotator = ann)
  expect_s3_class(expr, "polars_expr")
})

test_that("polars_numbering_method on a polars DataFrame produces a struct column", {
  skip_if_polars_missing()
  ann <- Annotator$new(chains = "IGH", scheme = "IMGT")
  df <- polars::pl$DataFrame(sequence = IGH_SEQ)
  result <- df$select(
    polars_numbering_method(polars::pl$col("sequence"),
                            annotator = ann)$alias("numbered")
  )
  expect_true("numbered" %in% result$columns)
  expect_equal(result$height, 1L)
})

test_that("polars_numbering_method struct exposes chain/scheme/confidence/numbering/error", {
  skip_if_polars_missing()
  ann <- Annotator$new(chains = "IGH", scheme = "IMGT")
  df <- polars::pl$DataFrame(sequence = IGH_SEQ)
  result <- df$select(
    polars_numbering_method(polars::pl$col("sequence"),
                            annotator = ann)$alias("numbered")
  )
  unnested <- result$unnest("numbered")
  for (field in c("chain", "scheme", "confidence", "error")) {
    expect_true(field %in% unnested$columns,
                info = paste("missing field:", field))
  }
})

test_that("polars_segmentation_method returns a polars expression", {
  skip_if_polars_missing()
  ann <- Annotator$new(chains = "IGH", scheme = "IMGT")
  expr <- polars_segmentation_method(polars::pl$col("sequence"), annotator = ann)
  expect_s3_class(expr, "polars_expr")
})

test_that("polars_segmentation_method on a polars DataFrame produces a struct column", {
  skip_if_polars_missing()
  ann <- Annotator$new(chains = "IGH", scheme = "IMGT")
  df <- polars::pl$DataFrame(sequence = IGH_SEQ)
  result <- df$select(
    polars_segmentation_method(polars::pl$col("sequence"),
                               annotator = ann)$alias("segmented")
  )
  expect_true("segmented" %in% result$columns)
  expect_equal(result$height, 1L)
})

test_that("polars_*_method requires an Annotator R6 instance", {
  skip_if_polars_missing()
  expect_error(
    polars_numbering_method(polars::pl$col("sequence"), annotator = "nope"),
    "must be an `Annotator` R6 instance"
  )
})
