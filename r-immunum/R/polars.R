# Polars batch path -- mirrors immunum/polars.py.
#
# Uses map_batches() to dispatch to rayon-parallel Rust batch functions.
# Helpers live in R/polars-utils.R.

#' Number sequences inside a polars expression
#'
#' Mirrors `immunum.polars.number()` from the Python wrapper. Returns a
#' struct with fields `chain`, `scheme`, `confidence`, `positions`, `residues`.
#' Failed sequences produce nulls.
#'
#' @param expr A polars expression or column name (string).
#' @param chains Character vector of chain identifiers (case-insensitive).
#' @param scheme `"IMGT"` or `"Kabat"`, case-insensitive.
#' @param min_confidence Optional numeric in `[0, 1]`. Default `NULL` uses 0.5.
#' @return A `polars` expression returning a struct column.
#' @examples
#' \dontrun{
#' df <- polars::pl$DataFrame(sequence = c("QVQLVQSGAEV...", "DIQMTQSPSSL..."))
#' df$with_columns(
#'   polars_number("sequence", chains = c("H", "K", "L"), scheme = "imgt")$alias("numbered")
#' )
#' }
#' @export
polars_number <- function(expr, chains, scheme, min_confidence = NULL) {
  .require_polars()
  canon_chains <- normalize_chains(chains)
  canon_scheme <- normalize_scheme(scheme)
  conf <- .check_min_confidence(min_confidence)
  ret <- .numbering_dtype()

  expr <- .as_polars_expr(expr)
  expr$map_batches(
    function(s) {
      seqs <- s$to_r_vector()
      cols <- .Call(wrap__numbering_batch, seqs, canon_chains, canon_scheme, conf)
      .numbering_to_struct_series(cols, s$name)
    },
    return_dtype = ret
  )
}

#' Segment sequences inside a polars expression
#'
#' Mirrors `immunum.polars.segment()` from the Python wrapper. Returns a
#' struct with fields `prefix`, `fr1`, `cdr1`, `fr2`, `cdr2`, `fr3`,
#' `cdr3`, `fr4`, `postfix`. Failed sequences produce nulls.
#'
#' @inheritParams polars_number
#' @return A `polars` expression returning a struct column.
#' @examples
#' \dontrun{
#' df <- polars::pl$DataFrame(sequence = c("QVQLVQSGAEV...", "DIQMTQSPSSL..."))
#' df$with_columns(
#'   polars_segment("sequence", chains = c("H", "K", "L"), scheme = "imgt")$alias("segmented")
#' )
#' }
#' @export
polars_segment <- function(expr, chains, scheme, min_confidence = NULL) {
  .require_polars()
  canon_chains <- normalize_chains(chains)
  canon_scheme <- normalize_scheme(scheme)
  conf <- .check_min_confidence(min_confidence)
  ret <- .segmentation_dtype()

  expr <- .as_polars_expr(expr)
  expr$map_batches(
    function(s) {
      seqs <- s$to_r_vector()
      cols <- .Call(wrap__segmentation_batch, seqs, canon_chains, canon_scheme, conf)
      .segmentation_to_struct_series(cols, s$name)
    },
    return_dtype = ret
  )
}

#' Number sequences using a pre-built Annotator (polars)
#'
#' Mirrors `immunum.polars.numbering_method()`. Reuses an [Annotator]
#' instance to skip per-batch construction.
#'
#' @param expr A polars expression or column name (string).
#' @param annotator A pre-built [Annotator] R6 instance.
#' @return A `polars` expression returning a struct column.
#' @examples
#' \dontrun{
#' ann <- Annotator$new(chains = c("H", "K", "L"), scheme = "imgt")
#' df <- polars::pl$DataFrame(sequence = c("QVQLVQSGAEV...", "DIQMTQSPSSL..."))
#' df$with_columns(
#'   polars_numbering_method("sequence", annotator = ann)$alias("numbered")
#' )
#' }
#' @export
polars_numbering_method <- function(expr, annotator) {
  .require_polars()
  inner <- .annotator_inner(annotator)
  ret <- .numbering_dtype()

  expr <- .as_polars_expr(expr)
  expr$map_batches(
    function(s) {
      seqs <- s$to_r_vector()
      cols <- .Call(wrap__numbering_batch_with, seqs, inner)
      .numbering_to_struct_series(cols, s$name)
    },
    return_dtype = ret
  )
}

#' Segment sequences using a pre-built Annotator (polars)
#'
#' Mirrors `immunum.polars.segmentation_method()`. Reuses an [Annotator]
#' instance to skip per-batch construction.
#'
#' @inheritParams polars_numbering_method
#' @return A `polars` expression returning a struct column.
#' @examples
#' \dontrun{
#' ann <- Annotator$new(chains = c("H", "K", "L"), scheme = "imgt")
#' df <- polars::pl$DataFrame(sequence = c("QVQLVQSGAEV...", "DIQMTQSPSSL..."))
#' df$with_columns(
#'   polars_segmentation_method("sequence", annotator = ann)$alias("segmented")
#' )
#' }
#' @export
polars_segmentation_method <- function(expr, annotator) {
  .require_polars()
  inner <- .annotator_inner(annotator)
  ret <- .segmentation_dtype()

  expr <- .as_polars_expr(expr)
  expr$map_batches(
    function(s) {
      seqs <- s$to_r_vector()
      cols <- .Call(wrap__segmentation_batch_with, seqs, inner)
      .segmentation_to_struct_series(cols, s$name)
    },
    return_dtype = ret
  )
}
