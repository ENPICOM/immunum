# Validation fixture manifest and benchmark thresholds.
#
# Fixtures live at fixtures/validation/ in the repo root (not shipped
# in the installed package). Thresholds are baked in from BENCHMARKS.toml.

.VALIDATION_FIXTURES <- data.frame(
  stem      = c(
    "ab_H_imgt",  "ab_K_imgt",  "ab_L_imgt",
    "ab_H_kabat", "ab_K_kabat", "ab_L_kabat",
    "tcr_A_imgt", "tcr_B_imgt", "tcr_G_imgt", "tcr_D_imgt"
  ),
  scheme    = c(
    "IMGT",  "IMGT",  "IMGT",
    "Kabat", "Kabat", "Kabat",
    "IMGT",  "IMGT",  "IMGT", "IMGT"
  ),
  benchmark = c(
    "imgt.H",  "imgt.K",  "imgt.L",
    "kabat.H", "kabat.K", "kabat.L",
    "imgt.A",  "imgt.B",  "imgt.G", "imgt.D"
  ),
  stringsAsFactors = FALSE
)
.VALIDATION_FIXTURES$chains <- list(
  "IGH", "IGK", "IGL",
  "IGH", "IGK", "IGL",
  "TRA", "TRB", "TRG", "TRD"
)

# Perfect_pct values from BENCHMARKS.toml. Update via tools/sync-benchmarks.R.
.BENCHMARK_THRESHOLDS <- c(
  `imgt.A` = 88.21,
  `imgt.B` = 97.54,
  `imgt.D` = 100,
  `imgt.G` = 96,
  `imgt.H` = 99.88,
  `imgt.K` = 99.73,
  `imgt.L` = 99.46,
  `kabat.H` = 99.88,
  `kabat.K` = 99.66,
  `kabat.L` = 99.46
)

#' List the available validation fixtures
#'
#' Returns the manifest of upstream validation fixtures.
#'
#' @return A `data.frame` with columns `stem`, `chains` (list-column),
#'   `scheme`, and `benchmark`.
#' @examples
#' validation_fixtures()
#' @export
validation_fixtures <- function() {
  .VALIDATION_FIXTURES
}

#' Look up the benchmark perfect-percentage threshold
#'
#' @param benchmark_key Dotted key, e.g. `"imgt.H"` or `"kabat.K"`.
#' @return A single numeric in `[0, 100]`.
#' @examples
#' benchmark_threshold("imgt.H")
#' @export
benchmark_threshold <- function(benchmark_key) {
  if (!is.character(benchmark_key) || length(benchmark_key) != 1L ||
      is.na(benchmark_key)) {
    stop("`benchmark_key` must be a single non-NA character string",
         call. = FALSE)
  }
  out <- .BENCHMARK_THRESHOLDS[benchmark_key]
  if (is.na(out)) {
    stop(
      sprintf(
        "Unknown benchmark key: %s. Known: %s",
        shQuote(benchmark_key),
        paste(names(.BENCHMARK_THRESHOLDS), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  unname(out)
}
