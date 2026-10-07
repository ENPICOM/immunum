skip_if_not_loaded <- function() {
  testthat::skip_if_not(
    requireNamespace("immunum", quietly = TRUE),
    "immunum native library not loaded"
  )
}

BATCH_SEQS <- c(IGH_SEQ, IGL_SEQ, TRB_SEQ)

# ── numbering_batch columns ──────────────────────────────────────────────────

test_that("numbering_batch returns expected column names", {
  skip_if_not_loaded()
  out <- immunum:::numbering_batch(BATCH_SEQS, ALL_CHAINS_R, "IMGT")
  expect_named(out, c("chain", "scheme", "confidence", "positions", "residues", "error"),
               ignore.order = TRUE)
})

test_that("numbering_batch returns one row per input sequence", {
  skip_if_not_loaded()
  out <- immunum:::numbering_batch(BATCH_SEQS, ALL_CHAINS_R, "IMGT")
  expect_length(out$chain,      length(BATCH_SEQS))
  expect_length(out$confidence, length(BATCH_SEQS))
  expect_length(out$positions,  length(BATCH_SEQS))
  expect_length(out$residues,   length(BATCH_SEQS))
  expect_length(out$error,      length(BATCH_SEQS))
})

test_that("numbering_batch results match single-sequence Annotator$number()", {
  skip_if_not_loaded()
  seqs <- c(IGH_SEQ, IGL_SEQ)
  ann  <- Annotator$new(chains = ALL_CHAINS_R, scheme = "IMGT")
  bat  <- immunum:::numbering_batch(seqs, ALL_CHAINS_R, "IMGT")

  for (i in seq_along(seqs)) {
    single <- ann$number(seqs[[i]])
    expect_equal(bat$chain[[i]],      single$chain)
    expect_equal(bat$confidence[[i]], single$confidence, tolerance = 1e-6)
    expect_equal(bat$positions[[i]],  names(single$numbering))
    expect_equal(bat$residues[[i]],   unname(single$numbering))
    expect_true(is.na(bat$error[[i]]))
  }
})

test_that("numbering_batch error field is set for bad sequences", {
  skip_if_not_loaded()
  seqs <- c(IGH_SEQ, "AAAAAAAAAAAAAAAA", IGL_SEQ)
  out  <- immunum:::numbering_batch(seqs, ALL_CHAINS_R, "IMGT")
  expect_true(is.na(out$error[[1]]))
  expect_false(is.na(out$error[[2]]))
  expect_true(is.na(out$chain[[2]]))
  expect_true(is.na(out$error[[3]]))
})

test_that("numbering_batch treats NA input as missing (no error field set)", {
  skip_if_not_loaded()
  seqs <- c(IGH_SEQ, NA_character_, IGL_SEQ)
  out  <- immunum:::numbering_batch(seqs, ALL_CHAINS_R, "IMGT")
  expect_true(is.na(out$chain[[2]]))
  expect_true(is.na(out$error[[2]]))
})

test_that("numbering_batch works with Kabat scheme", {
  skip_if_not_loaded()
  out <- immunum:::numbering_batch(c(IGH_SEQ, IGL_SEQ), AB_CHAINS_R, "Kabat")
  expect_true(is.na(out$error[[1]]))
  expect_true(is.na(out$error[[2]]))
  expect_equal(out$scheme[[1]], "Kabat")
})

test_that("numbering_batch accepts chain/scheme aliases", {
  skip_if_not_loaded()
  out1 <- immunum:::numbering_batch(c(IGH_SEQ), c("IGH"), "IMGT")
  out2 <- immunum:::numbering_batch(c(IGH_SEQ), c("H"),   "imgt")
  out3 <- immunum:::numbering_batch(c(IGH_SEQ), c("heavy"), "i")
  expect_equal(out1$chain, out2$chain)
  expect_equal(out1$chain, out3$chain)
})

# ── segmentation_batch columns ────────────────────────────────────────────────

test_that("segmentation_batch returns expected column names", {
  skip_if_not_loaded()
  out <- immunum:::segmentation_batch(BATCH_SEQS, ALL_CHAINS_R, "IMGT")
  expect_named(out, c("prefix", "fr1", "cdr1", "fr2", "cdr2",
                      "fr3", "cdr3", "fr4", "postfix", "error"),
               ignore.order = TRUE)
})

test_that("segmentation_batch results match single-sequence Annotator$segment()", {
  skip_if_not_loaded()
  seqs   <- c(IGH_SEQ, IGL_SEQ)
  ann    <- Annotator$new(chains = ALL_CHAINS_R, scheme = "IMGT")
  bat    <- immunum:::segmentation_batch(seqs, ALL_CHAINS_R, "IMGT")
  fields <- c("prefix", "fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4", "postfix")

  for (i in seq_along(seqs)) {
    single <- ann$segment(seqs[[i]])
    for (f in fields) {
      expect_equal(bat[[f]][[i]], single[[f]],
                   info = sprintf("row %d field %s", i, f))
    }
    expect_true(is.na(bat$error[[i]]))
  }
})

test_that("segmentation_batch error field is set for bad sequences", {
  skip_if_not_loaded()
  seqs <- c(IGH_SEQ, "AAAAAAAAAAAAAAAA")
  out  <- immunum:::segmentation_batch(seqs, ALL_CHAINS_R, "IMGT")
  expect_true(is.na(out$error[[1]]))
  expect_false(is.na(out$error[[2]]))
  expect_true(is.na(out$fr1[[2]]))
})

# ── _with variants ────────────────────────────────────────────────────────────

test_that("numbering_batch_with returns the same result as numbering_batch", {
  skip_if_not_loaded()
  ann  <- Annotator$new(chains = ALL_CHAINS_R, scheme = "IMGT")
  seqs <- c(IGH_SEQ, IGL_SEQ)
  ref  <- immunum:::numbering_batch(seqs, ALL_CHAINS_R, "IMGT")
  got  <- immunum:::numbering_batch_with(seqs, ann)

  expect_equal(got$chain,      ref$chain)
  expect_equal(got$scheme,     ref$scheme)
  expect_equal(got$confidence, ref$confidence, tolerance = 1e-6)
  expect_equal(got$positions,  ref$positions)
  expect_equal(got$residues,   ref$residues)
})

test_that("segmentation_batch_with returns the same result as segmentation_batch", {
  skip_if_not_loaded()
  ann    <- Annotator$new(chains = ALL_CHAINS_R, scheme = "IMGT")
  seqs   <- c(IGH_SEQ, IGL_SEQ)
  ref    <- immunum:::segmentation_batch(seqs, ALL_CHAINS_R, "IMGT")
  got    <- immunum:::segmentation_batch_with(seqs, ann)
  fields <- c("prefix", "fr1", "cdr1", "fr2", "cdr2",
              "fr3", "cdr3", "fr4", "postfix", "error")

  for (f in fields) {
    expect_equal(got[[f]], ref[[f]], info = paste("field:", f))
  }
})

test_that("_with variants reject non-Annotator inputs", {
  skip_if_not_loaded()
  expect_error(immunum:::numbering_batch_with(IGH_SEQ, "nope"),
               "must be an Annotator R6 instance")
  expect_error(immunum:::segmentation_batch_with(IGH_SEQ, list()),
               "must be an Annotator R6 instance")
})
