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
