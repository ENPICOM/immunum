# Regenerate the .BENCHMARK_THRESHOLDS constant in R/validation.R from
# the source-of-truth BENCHMARKS.toml at the repo root.
#
# The R package can't read BENCHMARKS.toml at runtime without pulling
# in a TOML parser as a hard dependency. Instead the perfect-pct
# values are baked in as a static R constant; this script keeps them
# fresh after upstream re-runs the benchmark suite.
#
# Run from the r-immunum/ directory after BENCHMARKS.toml changes:
#   Rscript --vanilla tools/sync-benchmarks.R
#
# It rewrites the `.BENCHMARK_THRESHOLDS <- c(...)` block in
# R/validation.R in place. Diff afterwards before committing.

stopifnot(file.exists("../BENCHMARKS.toml"),
          file.exists("R/validation.R"))

lines <- readLines("../BENCHMARKS.toml")

# Tiny single-purpose parser for the section headers + perfect_pct
# lines we care about. BENCHMARKS.toml is a flat dict-of-dicts so we
# don't need a real TOML parser here.
section <- NULL
thresholds <- list()
for (line in lines) {
  line <- trimws(line)
  if (!nzchar(line) || startsWith(line, "#")) next
  m <- regmatches(line, regexec("^\\[(.+)\\]$", line))[[1]]
  if (length(m) == 2L) {
    section <- m[[2]]
    next
  }
  if (!is.null(section)) {
    kv <- regmatches(line, regexec("^perfect_pct\\s*=\\s*([0-9.]+)$", line))[[1]]
    if (length(kv) == 2L) {
      thresholds[[section]] <- as.numeric(kv[[2]])
    }
  }
}

if (length(thresholds) == 0L) {
  stop("No perfect_pct entries found in BENCHMARKS.toml -- aborting")
}

ordered_keys <- sort(names(thresholds))

# Pretty-print the constant body so it diffs cleanly. Backtick-quote
# the dotted keys so they're valid R names.
body_lines <- vapply(seq_along(ordered_keys), function(i) {
  k <- ordered_keys[[i]]
  v <- thresholds[[k]]
  comma <- if (i == length(ordered_keys)) "" else ","
  sprintf("  `%s` = %s%s", k, format(v, nsmall = 0L), comma)
}, character(1))

new_block <- c(".BENCHMARK_THRESHOLDS <- c(", body_lines, ")")

# Splice into R/validation.R between the existing markers.
src <- readLines("R/validation.R")
start <- grep("^\\.BENCHMARK_THRESHOLDS <- c\\($", src)
if (length(start) != 1L) {
  stop("Could not find a unique `.BENCHMARK_THRESHOLDS <- c(` line in R/validation.R")
}
end <- start - 1L + which(startsWith(src[start:length(src)], ")"))[1]
if (is.na(end) || end <= start) {
  stop("Could not find the closing `)` for .BENCHMARK_THRESHOLDS")
}

updated <- c(src[seq_len(start - 1L)], new_block, src[(end + 1L):length(src)])
writeLines(updated, "R/validation.R")

cat("Updated R/validation.R with",
    length(thresholds), "thresholds from BENCHMARKS.toml\n")
for (k in ordered_keys) {
  cat(sprintf("  %-10s %s\n", k, format(thresholds[[k]])))
}
