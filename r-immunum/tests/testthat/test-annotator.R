skip_if_not_loaded <- function() {
  testthat::skip_if_not(
    requireNamespace("immunum", quietly = TRUE),
    "immunum native library not loaded"
  )
}

# ── Construction ────────────────────────────────────────────────────────────

test_that("Annotator constructs across the parametrized init matrix", {
  skip_if_not_loaded()
  for (case in INIT_CASES) {
    ann <- Annotator$new(chains = case$chains, scheme = case$scheme)
    expect_s3_class(ann, "Annotator")
  }
})

test_that("Annotator can number the seed sequence for each init case", {
  skip_if_not_loaded()
  for (case in INIT_CASES) {
    ann <- Annotator$new(chains = case$chains, scheme = case$scheme)
    result <- ann$number(case$seq)
    expect_type(result$chain, "character")
    expect_type(result$scheme, "character")
    expect_type(result$confidence, "double")
    expect_true(result$confidence >= 0 && result$confidence <= 1)
  }
})

test_that("Kabat scheme rejects TCR chains", {
  skip_if_not_loaded()
  bad_combos <- list(
    list(chains = "TRA",          scheme = "Kabat"),
    list(chains = "TRB",          scheme = "Kabat"),
    list(chains = c("IGH", "TRA"), scheme = "Kabat"),
    list(chains = ALL_CHAINS_R,   scheme = "Kabat")
  )
  for (combo in bad_combos) {
    expect_error(
      Annotator$new(chains = combo$chains, scheme = combo$scheme),
      "Kabat"
    )
  }
})

test_that("invalid construction args raise", {
  skip_if_not_loaded()
  expect_error(Annotator$new(chains = "INVALID", scheme = "IMGT"))
  expect_error(Annotator$new(chains = "IGH",     scheme = "INVALID"))
  expect_error(Annotator$new(chains = character(), scheme = "IMGT"))
})

test_that("min_confidence is range-checked", {
  skip_if_not_loaded()
  expect_error(
    Annotator$new(chains = "IGH", scheme = "IMGT", min_confidence = -0.1),
    "must be a single numeric in"
  )
  expect_error(
    Annotator$new(chains = "IGH", scheme = "IMGT", min_confidence = 1.1),
    "must be a single numeric in"
  )
  expect_error(
    Annotator$new(chains = "IGH", scheme = "IMGT", min_confidence = c(0.5, 0.6)),
    "must be a single numeric in"
  )
  expect_no_error(Annotator$new(chains = "IGH", scheme = "IMGT", min_confidence = NULL))
  expect_no_error(Annotator$new(chains = "IGH", scheme = "IMGT", min_confidence = 0))
  expect_no_error(Annotator$new(chains = "IGH", scheme = "IMGT", min_confidence = 1))
})

# ── Numbering ───────────────────────────────────────────────────────────────

test_that("number detects an IGH sequence as the H chain", {
  skip_if_not_loaded()
  ann <- Annotator$new(chains = ALL_CHAINS_R, scheme = "IMGT")
  result <- ann$number(IGH_SEQ)
  expect_equal(result$chain, "H")
  expect_equal(result$scheme, "IMGT")
  expect_null(result$error)
})

test_that("number with a single-chain annotator returns that chain", {
  skip_if_not_loaded()
  ann <- Annotator$new(chains = "IGH", scheme = "IMGT")
  result <- ann$number(IGH_SEQ)
  expect_equal(result$chain, "H")
})

test_that("number on an empty sequence returns error", {
  skip_if_not_loaded()
  ann <- Annotator$new(chains = ALL_CHAINS_R, scheme = "IMGT")
  result <- ann$number("")
  expect_false(is.null(result$error))
  expect_null(result$chain)
  expect_null(result$query_start)
  expect_null(result$query_end)
})

test_that("number on an invalid sequence returns error", {
  skip_if_not_loaded()
  ann <- Annotator$new(chains = ALL_CHAINS_R, scheme = "IMGT")
  result <- ann$number("AAAAAAAAAAAAAAAA")
  expect_false(is.null(result$error))
  expect_null(result$chain)
})

test_that("query_start and query_end are returned on success", {
  skip_if_not_loaded()
  ann <- Annotator$new(chains = "IGH", scheme = "IMGT")
  result <- ann$number(IGH_SEQ)
  expect_type(result$query_start, "integer")
  expect_type(result$query_end, "integer")
  expect_true(result$query_start >= 1L)
  expect_true(result$query_end >= result$query_start)
  expect_true(result$query_end <= nchar(IGH_SEQ))
  expect_null(result$error)
})

test_that("number rejects bad sequence input shapes", {
  skip_if_not_loaded()
  ann <- Annotator$new(chains = ALL_CHAINS_R, scheme = "IMGT")
  expect_error(ann$number(NA_character_), "single non-NA")
  expect_error(ann$number(c("AAA", "BBB")), "single non-NA")
  expect_error(ann$number(123), "single non-NA")
})

test_that("confidence is a numeric in [0, 1] for every init case", {
  skip_if_not_loaded()
  for (case in INIT_CASES) {
    ann <- Annotator$new(chains = case$chains, scheme = case$scheme)
    result <- ann$number(case$seq)
    expect_type(result$confidence, "double")
    expect_true(result$confidence >= 0 && result$confidence <= 1)
  }
})

test_that("numbering is a named character vector with one entry per residue", {
  skip_if_not_loaded()
  ann <- Annotator$new(chains = "IGH", scheme = "IMGT")
  result <- ann$number(IGH_SEQ)
  expect_type(result$numbering, "character")
  expect_true(length(result$numbering) > 0L)
  expect_false(is.null(names(result$numbering)))
  expect_true(all(nchar(result$numbering) == 1L))
  expect_true(all(nzchar(names(result$numbering))))
})

# ── Segmentation ────────────────────────────────────────────────────────────

test_that("segmentation returns the nine expected fields", {
  skip_if_not_loaded()
  expected_fields <- c("prefix", "fr1", "cdr1", "fr2", "cdr2",
                       "fr3", "cdr3", "fr4", "postfix")
  for (case in INIT_CASES) {
    ann <- Annotator$new(chains = case$chains, scheme = case$scheme)
    seg <- ann$segment(case$seq)
    expect_setequal(names(seg), expected_fields)
    for (field in expected_fields) {
      expect_type(seg[[field]], "character")
      expect_length(seg[[field]], 1L)
      expect_false(is.na(seg[[field]]))
    }
  }
})

test_that("segmentation on an invalid sequence returns error", {
  skip_if_not_loaded()
  ann <- Annotator$new(chains = "IGH", scheme = "IMGT")
  seg <- ann$segment("AAAAAAAAAAAAAAAA")
  expect_false(is.null(seg$error))
  expect_null(seg$fr1)
})

test_that("known IGH sequence segments to canonical FR1/CDR1/CDR3/FR4", {
  skip_if_not_loaded()
  ann <- Annotator$new(chains = "IGH", scheme = "IMGT")
  seg <- ann$segment(IGH_SEQ)
  expect_equal(seg$fr1, "QVQLVQSGAEVKRPGSSVTVSCKAS")
  expect_equal(seg$cdr1, "GGSFSTYA")
  expect_equal(seg$cdr3, "AREGTTGKPIGAFAH")
  expect_equal(seg$fr4, "WGQGTLVTVSS")
  expect_equal(seg$prefix, "")
  expect_equal(seg$postfix, "")
})
