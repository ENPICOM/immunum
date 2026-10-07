# immunum -- mirrors immunum/__init__.py.
#
# Contains chain/scheme normalization and the Annotator R6 class.

#' immunum: Fast Antibody and T-Cell Receptor Sequence Numbering
#'
#' R bindings for the `immunum` Rust crate. High-performance numbering
#' of antibody and T-cell receptor variable-domain sequences using IMGT
#' and Kabat schemes.
#'
#' @section Main entry points:
#' - [Annotator]: R6 class for numbering and segmenting sequences.
#' - [immunum_version()]: Version of the linked Rust core.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib immunum, .registration = TRUE
## usethis namespace: end
NULL

# ── Chain / scheme alias normalization ──────────────────────────────────────

.CHAIN_ALIASES <- c(
  igh = "IGH",   h = "IGH",   heavy  = "IGH",
  igk = "IGK",   k = "IGK",   kappa  = "IGK",
  igl = "IGL",   l = "IGL",   lambda = "IGL",
  tra = "TRA",   a = "TRA",   alpha  = "TRA",
  trb = "TRB",   b = "TRB",   beta   = "TRB",
  trg = "TRG",   g = "TRG",   gamma  = "TRG",
  trd = "TRD",   d = "TRD",   delta  = "TRD"
)

.SCHEME_ALIASES <- c(
  imgt  = "IMGT",
  i     = "IMGT",
  kabat = "Kabat",
  k     = "Kabat"
)

normalize_chains <- function(chains) {
  if (!is.character(chains)) {
    stop("`chains` must be a character vector", call. = FALSE)
  }
  if (length(chains) == 0L) {
    stop("`chains` cannot be empty", call. = FALSE)
  }
  if (anyNA(chains)) {
    stop("`chains` cannot contain NA", call. = FALSE)
  }

  out <- .CHAIN_ALIASES[tolower(chains)]
  bad <- is.na(out)
  if (any(bad)) {
    valid <- sort(unique(unname(.CHAIN_ALIASES)))
    stop(
      sprintf(
        "Unknown chain(s): %s. Valid chains: %s",
        paste(shQuote(chains[bad]), collapse = ", "),
        paste(valid, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  unname(out)
}

normalize_scheme <- function(scheme) {
  if (!is.character(scheme) || length(scheme) != 1L || is.na(scheme)) {
    stop("`scheme` must be a single non-NA character string", call. = FALSE)
  }
  out <- .SCHEME_ALIASES[tolower(scheme)]
  if (is.na(out)) {
    valid <- sort(unique(unname(.SCHEME_ALIASES)))
    stop(
      sprintf(
        "Unknown scheme: %s. Valid schemes: %s",
        shQuote(scheme),
        paste(valid, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  unname(out)
}

check_sequence <- function(sequence) {
  if (!is.character(sequence) || length(sequence) != 1L || is.na(sequence)) {
    stop("`sequence` must be a single non-NA character string", call. = FALSE)
  }
  invisible(TRUE)
}

# ── Annotator R6 class ──────────────────────────────────────────────────────

#' @rawNamespace export(Annotator)
#' @importFrom R6 R6Class
Annotator <- R6::R6Class(
  "Annotator",
  cloneable = FALSE,
  public = list(
    initialize = function(chains, scheme, min_confidence = NULL) {
      canon_chains <- normalize_chains(chains)
      canon_scheme <- normalize_scheme(scheme)

      tcr <- c("TRA", "TRB", "TRG", "TRD")
      if (canon_scheme == "Kabat" && any(canon_chains %in% tcr)) {
        stop(
          "Kabat scheme only supported for antibody chains (IGH, IGK, IGL)",
          call. = FALSE
        )
      }

      if (!is.null(min_confidence)) {
        if (!is.numeric(min_confidence) || length(min_confidence) != 1L ||
            is.na(min_confidence) || min_confidence < 0 || min_confidence > 1) {
          stop(
            sprintf(
              "`min_confidence` must be a single numeric in [0, 1], got %s",
              format(min_confidence)
            ),
            call. = FALSE
          )
        }
        min_confidence <- as.numeric(min_confidence)
      }

      private$.inner <- .Call(
        wrap__Annotator__new,
        canon_chains,
        canon_scheme,
        min_confidence
      )
      private$.chains <- canon_chains
      private$.scheme <- canon_scheme
      invisible(self)
    },

    number = function(sequence) {
      check_sequence(sequence)
      raw <- .Call(wrap__Annotator__number, private$.inner, sequence)
      if (!is.null(raw$error)) {
        return(list(
          chain       = NULL,
          scheme      = NULL,
          confidence  = NULL,
          numbering   = NULL,
          query_start = NULL,
          query_end   = NULL,
          error       = raw$error
        ))
      }
      list(
        chain       = raw$chain,
        scheme      = raw$scheme,
        confidence  = raw$confidence,
        numbering   = stats::setNames(raw$residues, raw$positions),
        query_start = raw$query_start,
        query_end   = raw$query_end,
        error       = NULL
      )
    },

    segment = function(sequence) {
      check_sequence(sequence)
      raw <- .Call(wrap__Annotator__segment, private$.inner, sequence)
      if (!is.null(raw$error)) {
        return(list(
          fr1     = NULL,
          cdr1    = NULL,
          fr2     = NULL,
          cdr2    = NULL,
          fr3     = NULL,
          cdr3    = NULL,
          fr4     = NULL,
          prefix  = NULL,
          postfix = NULL,
          error   = raw$error
        ))
      }
      raw$error <- NULL
      raw
    },

    print = function(...) {
      cat("<immunum::Annotator>\n")
      cat("  chains:        ", paste(private$.chains, collapse = ", "), "\n", sep = "")
      cat("  scheme:        ", private$.scheme, "\n", sep = "")
      invisible(self)
    }
  ),
  private = list(
    .inner  = NULL,
    .chains = NULL,
    .scheme = NULL
  )
)
