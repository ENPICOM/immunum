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
