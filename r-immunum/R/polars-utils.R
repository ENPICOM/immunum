# Internal helpers for the polars batch path.

.require_polars <- function() {
  if (!requireNamespace("polars", quietly = TRUE)) {
    stop(
      "The 'polars' package is required for immunum's polars batch path. ",
      "Install it from r-universe: ",
      "install.packages('polars', repos = c('https://community.r-multiverse.org', getOption('repos')))",
      call. = FALSE
    )
  }
}

.pl <- function() polars::pl

.annotator_inner <- function(annotator) {
  if (!inherits(annotator, "Annotator")) {
    stop("`annotator` must be an `Annotator` R6 instance", call. = FALSE)
  }
  annotator$.__enclos_env__$private$.inner
}

.check_min_confidence <- function(min_confidence) {
  if (is.null(min_confidence)) return(NULL)
  if (!is.numeric(min_confidence) || length(min_confidence) != 1L ||
      is.na(min_confidence) || min_confidence < 0 || min_confidence > 1) {
    stop(
      sprintf("`min_confidence` must be a single numeric in [0, 1], got %s",
              format(min_confidence)),
      call. = FALSE
    )
  }
  as.numeric(min_confidence)
}

.numbering_dtype <- function() {
  pl <- .pl()
  pl$Struct(
    chain      = pl$String,
    scheme     = pl$String,
    confidence = pl$Float64,
    positions  = pl$List(pl$String),
    residues   = pl$List(pl$String),
    error      = pl$String
  )
}

.segmentation_dtype <- function() {
  pl <- .pl()
  pl$Struct(
    prefix  = pl$String,
    fr1     = pl$String,
    cdr1    = pl$String,
    fr2     = pl$String,
    cdr2    = pl$String,
    fr3     = pl$String,
    cdr3    = pl$String,
    fr4     = pl$String,
    postfix = pl$String,
    error   = pl$String
  )
}

# Uses nanoarrow for List(String) columns (~25x faster than direct R -> polars).
.numbering_to_struct_series <- function(cols, name) {
  pl <- .pl()
  if (requireNamespace("nanoarrow", quietly = TRUE)) {
    s_chain  <- polars::as_polars_series(cols$chain,      name = "chain")
    s_scheme <- polars::as_polars_series(cols$scheme,     name = "scheme")
    s_conf   <- polars::as_polars_series(cols$confidence, name = "confidence")
    s_pos    <- polars::as_polars_series(
      nanoarrow::as_nanoarrow_array(cols$positions,
                                    schema = nanoarrow::na_list(nanoarrow::na_string())),
      name = "positions"
    )
    s_res    <- polars::as_polars_series(
      nanoarrow::as_nanoarrow_array(cols$residues,
                                    schema = nanoarrow::na_list(nanoarrow::na_string())),
      name = "residues"
    )
    s_err    <- polars::as_polars_series(cols$error,      name = "error")
    pl$DataFrame(s_chain, s_scheme, s_conf, s_pos, s_res, s_err)$to_struct(name)
  } else {
    do.call(pl$DataFrame, cols)$to_struct(name)
  }
}

.segmentation_to_struct_series <- function(cols, name) {
  pl <- .pl()
  do.call(pl$DataFrame, cols)$to_struct(name)
}

.as_polars_expr <- function(expr) {
  if (is.character(expr) && length(expr) == 1L && !is.na(expr)) {
    return(.pl()$col(expr))
  }
  polars::as_polars_expr(expr)
}
