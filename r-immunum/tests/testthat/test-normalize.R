.normalize_chains <- function(...) immunum:::normalize_chains(...)
.normalize_scheme <- function(...) immunum:::normalize_scheme(...)

test_that("short codes resolve to canonical chain codes", {
  expect_equal(.normalize_chains("H"), "IGH")
  expect_equal(.normalize_chains("K"), "IGK")
  expect_equal(.normalize_chains("L"), "IGL")
  expect_equal(.normalize_chains("A"), "TRA")
  expect_equal(.normalize_chains("B"), "TRB")
  expect_equal(.normalize_chains("G"), "TRG")
  expect_equal(.normalize_chains("D"), "TRD")
})

test_that("named aliases resolve to canonical chain codes", {
  expect_equal(.normalize_chains("heavy"),  "IGH")
  expect_equal(.normalize_chains("kappa"),  "IGK")
  expect_equal(.normalize_chains("lambda"), "IGL")
  expect_equal(.normalize_chains("alpha"),  "TRA")
  expect_equal(.normalize_chains("beta"),   "TRB")
  expect_equal(.normalize_chains("gamma"),  "TRG")
  expect_equal(.normalize_chains("delta"),  "TRD")
})

test_that("canonical codes pass through unchanged", {
  expect_equal(.normalize_chains("IGH"), "IGH")
  expect_equal(.normalize_chains(ALL_CHAINS_R), ALL_CHAINS_R)
})

test_that("normalization is case-insensitive", {
  expect_equal(.normalize_chains("igh"), "IGH")
  expect_equal(.normalize_chains("Heavy"), "IGH")
  expect_equal(.normalize_chains("LAMBDA"), "IGL")
  expect_equal(.normalize_chains(c("h", "K", "Lambda")), c("IGH", "IGK", "IGL"))
})

test_that("unknown chain raises informative error", {
  expect_error(.normalize_chains("INVALID"), "Unknown chain")
  expect_error(.normalize_chains("Z"),       "Unknown chain")
  expect_error(.normalize_chains("IGX"),     "Unknown chain")
  expect_error(.normalize_chains(c("H", "INVALID")), "Unknown chain")
})

test_that("invalid input shape raises", {
  expect_error(.normalize_chains(character()), "cannot be empty")
  expect_error(.normalize_chains(NA_character_), "cannot contain NA")
  expect_error(.normalize_chains(123), "must be a character")
})

test_that("scheme aliases resolve to canonical case", {
  expect_equal(.normalize_scheme("IMGT"),  "IMGT")
  expect_equal(.normalize_scheme("imgt"),  "IMGT")
  expect_equal(.normalize_scheme("i"),     "IMGT")
  expect_equal(.normalize_scheme("Kabat"), "Kabat")
  expect_equal(.normalize_scheme("kabat"), "Kabat")
  expect_equal(.normalize_scheme("k"),     "Kabat")
  expect_equal(.normalize_scheme("ImGt"),  "IMGT")
})

test_that("unknown scheme raises informative error", {
  expect_error(.normalize_scheme("INVALID"), "Unknown scheme")
  expect_error(.normalize_scheme("xyz"),     "Unknown scheme")
})

test_that("invalid scheme input shape raises", {
  expect_error(.normalize_scheme(character()), "single non-NA character")
  expect_error(.normalize_scheme(c("IMGT", "Kabat")), "single non-NA character")
  expect_error(.normalize_scheme(NA_character_), "single non-NA character")
})

# ── End-to-end alias equivalence ─────────────────────────────────────────────
# Mirrors TestNormalization::test_alias_produces_identical_result from
# test_python.py: using an alias must produce byte-identical chain, scheme,
# and numbering output as the canonical form.

test_that("alias annotators produce identical numbering to canonical", {
  skip_if_not(
    requireNamespace("immunum", quietly = TRUE),
    "immunum native library not loaded"
  )
  cases <- list(
    list(alias = "H",       canonical = "IGH",       scheme_a = "IMGT",  scheme_c = "IMGT",  seq = IGH_SEQ),
    list(alias = "K",       canonical = "IGK",       scheme_a = "IMGT",  scheme_c = "IMGT",  seq = IGL_SEQ),
    list(alias = "L",       canonical = "IGL",       scheme_a = "IMGT",  scheme_c = "IMGT",  seq = IGL_SEQ),
    list(alias = "A",       canonical = "TRA",       scheme_a = "IMGT",  scheme_c = "IMGT",  seq = TRA_SEQ),
    list(alias = "B",       canonical = "TRB",       scheme_a = "IMGT",  scheme_c = "IMGT",  seq = TRB_SEQ),
    list(alias = "G",       canonical = "TRG",       scheme_a = "IMGT",  scheme_c = "IMGT",  seq = TRG_SEQ),
    list(alias = "D",       canonical = "TRD",       scheme_a = "IMGT",  scheme_c = "IMGT",  seq = TRD_SEQ),
    list(alias = "heavy",   canonical = "IGH",       scheme_a = "IMGT",  scheme_c = "IMGT",  seq = IGH_SEQ),
    list(alias = "kappa",   canonical = "IGK",       scheme_a = "IMGT",  scheme_c = "IMGT",  seq = IGL_SEQ),
    list(alias = "lambda",  canonical = "IGL",       scheme_a = "IMGT",  scheme_c = "IMGT",  seq = IGL_SEQ),
    list(alias = "alpha",   canonical = "TRA",       scheme_a = "IMGT",  scheme_c = "IMGT",  seq = TRA_SEQ),
    list(alias = "beta",    canonical = "TRB",       scheme_a = "IMGT",  scheme_c = "IMGT",  seq = TRB_SEQ),
    list(alias = "gamma",   canonical = "TRG",       scheme_a = "IMGT",  scheme_c = "IMGT",  seq = TRG_SEQ),
    list(alias = "delta",   canonical = "TRD",       scheme_a = "IMGT",  scheme_c = "IMGT",  seq = TRD_SEQ),
    list(alias = "igh",     canonical = "IGH",       scheme_a = "imgt",  scheme_c = "IMGT",  seq = IGH_SEQ),
    list(alias = "IGH",     canonical = "IGH",       scheme_a = "i",     scheme_c = "IMGT",  seq = IGH_SEQ),
    list(alias = AB_CHAINS_R, canonical = AB_CHAINS_R, scheme_a = "k",   scheme_c = "Kabat", seq = IGL_SEQ)
  )
  for (case in cases) {
    alias_r    <- Annotator$new(chains = case$alias,    scheme = case$scheme_a)$number(case$seq)
    canon_r    <- Annotator$new(chains = case$canonical, scheme = case$scheme_c)$number(case$seq)
    expect_equal(alias_r$chain,    canon_r$chain,
                 info = paste("chain mismatch for alias", paste(case$alias, collapse = "+")))
    expect_equal(alias_r$scheme,   canon_r$scheme,
                 info = paste("scheme mismatch for alias", paste(case$alias, collapse = "+")))
    expect_equal(alias_r$numbering, canon_r$numbering,
                 info = paste("numbering mismatch for alias", paste(case$alias, collapse = "+")))
  }
})
