skip_if_no_python_immunum <- function() {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if(
    !nzchar(Sys.which("python3")) && !nzchar(Sys.which("python")),
    "No Python installation found"
  )
  old <- Sys.getenv("RETICULATE_PYTHON_FALLBACK", NA)
  Sys.setenv(RETICULATE_PYTHON_FALLBACK = "false")
  on.exit({
    if (is.na(old)) Sys.unsetenv("RETICULATE_PYTHON_FALLBACK")
    else Sys.setenv(RETICULATE_PYTHON_FALLBACK = old)
  })
  tryCatch(
    reticulate::import("immunum"),
    error = function(e) {
      testthat::skip("Python immunum not importable")
    }
  )
}

# ── Single-sequence parity ──────────────────────────────────────────────────

test_that("R and Python Annotator$number() produce identical output", {
  skip_if_no_python_immunum()
  imp <- reticulate::import("immunum")

  py_ann <- imp$Annotator(chains = list("H", "K", "L"), scheme = "imgt")
  r_ann  <- Annotator$new(chains = c("H", "K", "L"), scheme = "imgt")

  for (seq in c(IGH_SEQ, IGL_SEQ)) {
    py_res <- py_ann$number(seq)
    r_res  <- r_ann$number(seq)

    expect_identical(r_res$chain, py_res$chain,
                     info = sprintf("chain mismatch for %.20s...", seq))
    expect_identical(r_res$scheme, py_res$scheme,
                     info = sprintf("scheme mismatch for %.20s...", seq))
    expect_equal(r_res$confidence, py_res$confidence, tolerance = 1e-6,
                 info = sprintf("confidence mismatch for %.20s...", seq))

    py_numb <- unlist(reticulate::py_to_r(py_res$numbering))
    r_numb  <- r_res$numbering
    expect_identical(sort(names(r_numb)), sort(names(py_numb)),
                     info = sprintf("position set mismatch for %.20s...", seq))
    expect_identical(unname(r_numb[sort(names(r_numb))]),
                     unname(py_numb[sort(names(py_numb))]),
                     info = sprintf("residue values mismatch for %.20s...", seq))
  }
})

test_that("R and Python Annotator$segment() produce identical output", {
  skip_if_no_python_immunum()
  imp <- reticulate::import("immunum")

  py_ann <- imp$Annotator(chains = list("H", "K", "L"), scheme = "imgt")
  r_ann  <- Annotator$new(chains = c("H", "K", "L"), scheme = "imgt")

  regions <- c("prefix", "fr1", "cdr1", "fr2", "cdr2",
               "fr3", "cdr3", "fr4", "postfix")

  for (seq in c(IGH_SEQ, IGL_SEQ)) {
    py_res <- py_ann$segment(seq)
    r_res  <- r_ann$segment(seq)

    for (region in regions) {
      py_val <- py_res[[region]]
      r_val  <- r_res[[region]]
      if (is.null(py_val) || identical(py_val, "")) py_val <- NA_character_
      if (is.null(r_val)  || identical(r_val, ""))  r_val  <- NA_character_
      expect_identical(r_val, py_val,
                       info = sprintf("%s mismatch for %.20s...", region, seq))
    }
  }
})

# ── Polars batch parity ────────────────────────────────────────────────────

test_that("R and Python polars_number produce identical struct output", {
  skip_if_no_python_immunum()
  testthat::skip_if_not_installed("polars")

  reticulate::py_run_string("
import polars as pl
import immunum.polars as imp

seqs = [
    'QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS',
    'SALTQPPAVSGTPGQRVTISCSGSDIGRRSVNWYQQFPGTAPKLLIYSNDQRPSVVPDRFSGSKSGTSASLAISGLQSEDEAEYYCAAWDDSLAVFGGGTQLTVGQPKA',
]

df = pl.DataFrame({'sequence': seqs})
result = df.select(
    imp.number(pl.col('sequence'), chains=['H', 'K', 'L'], scheme='IMGT',
               min_confidence=0.0).alias('numbered')
).unnest('numbered')

py_chains = result.get_column('chain').to_list()
py_positions = result.get_column('positions').to_list()
py_residues = result.get_column('residues').to_list()
")

  py_chains    <- reticulate::py$py_chains
  py_positions <- reticulate::py$py_positions
  py_residues  <- reticulate::py$py_residues

  seqs <- c(IGH_SEQ, IGL_SEQ)
  df <- polars::pl$DataFrame(sequence = seqs)
  result <- df$select(
    polars_number(polars::pl$col("sequence"),
                  chains = c("H", "K", "L"), scheme = "IMGT",
                  min_confidence = 0.0)$alias("numbered")
  )$unnest("numbered")

  r_chains    <- result$get_column("chain")$to_r_vector()
  r_positions <- as.list(result$get_column("positions"))
  r_residues  <- as.list(result$get_column("residues"))

  expect_identical(r_chains, py_chains, info = "chain detection differs")

  for (i in seq_along(seqs)) {
    expect_identical(r_positions[[i]], py_positions[[i]],
                     info = sprintf("row %d position lists differ", i))
    expect_identical(r_residues[[i]], py_residues[[i]],
                     info = sprintf("row %d residue lists differ", i))
  }
})
