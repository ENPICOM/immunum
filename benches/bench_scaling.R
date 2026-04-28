# Scaling benchmark: R-immunum vs Python-immunum across batch sizes.
# Produces a comparison chart at docs/assets/benchmark_r_vs_python.svg

library(immunum)
library(polars)
library(ggplot2)

FIXTURES <- normalizePath(file.path("..", "fixtures", "validation"), mustWork = TRUE)
fixture <- file.path(FIXTURES, "ab_H_imgt.csv")
SIZES <- c(100L, 500L, 1000L, 5000L, 10000L, 50000L)
ROUNDS <- 3L
SEED <- 42L

df_full <- pl$read_csv(fixture, infer_schema_length = 0L)

rows <- list()

for (size in SIZES) {
  cat(sprintf("\n=== size = %d ===\n", size))
  df <- df_full$sample(n = size, with_replacement = TRUE, seed = SEED)

  # R polars batch
  times_r <- numeric(ROUNDS)
  for (r in seq_len(ROUNDS)) {
    t0 <- proc.time()["elapsed"]
    df$select(
      polars_number(pl$col("sequence"),
                    chains = "IGH", scheme = "IMGT",
                    min_confidence = 0.0)$alias("numbered")
    )
    times_r[r] <- proc.time()["elapsed"] - t0
  }
  med_r <- median(times_r)
  cat(sprintf("  R polars:      %.3fs\n", med_r))
  rows <- c(rows, list(data.frame(size = size, tool = "R (polars batch)", time_s = med_r)))

  # Python polars batch
  reticulate::py_run_string(sprintf("
import polars
import immunum.polars as imp
import time

df = polars.read_csv('%s', infer_schema=False).sample(n=%d, with_replacement=True, seed=%d)
times = []
for _ in range(%d):
    t0 = time.perf_counter()
    df.select(imp.number(polars.col('sequence'), chains=['IGH'], scheme='IMGT', min_confidence=0.0).alias('numbered'))
    times.append(time.perf_counter() - t0)

median_s = sorted(times)[len(times) // 2]
", gsub("\\\\", "/", fixture), size, SEED, ROUNDS))
  med_py <- reticulate::py$median_s
  cat(sprintf("  Python polars: %.3fs\n", med_py))
  rows <- c(rows, list(data.frame(size = size, tool = "Python (polars batch)", time_s = med_py)))
}

results <- do.call(rbind, rows)
cat("\n=== Results ===\n")
print(results)

# Write CSV
csv_path <- file.path("..", "resources", "benchmark_results", "results_r_vs_python.csv")
write.csv(results, csv_path, row.names = FALSE)
cat(sprintf("\nCSV saved to %s\n", csv_path))

# Produce chart
p <- ggplot(results, aes(x = size, y = time_s, color = tool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  labs(
    title = "immunum: R vs Python polars batch numbering",
    subtitle = "IGH / IMGT, median of 3 rounds",
    x = "Batch size",
    y = "Time (seconds, log scale)",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

svg_path <- file.path("..", "docs", "assets", "benchmark_r_vs_python.svg")
ggsave(svg_path, p, width = 8, height = 5)
cat(sprintf("Chart saved to %s\n", svg_path))
