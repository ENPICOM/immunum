# Benchmark immunum R package throughput.
# Usage:
#   Rscript bench/bench-annotate.R
#   RAYON_NUM_THREADS=1 Rscript bench/bench-annotate.R   # single-thread baseline
# No dependencies beyond the installed immunum package.

library(immunum)

# Time a closure over `reps` iterations, return median of `rounds` such blocks.
.bench_s <- function(f, reps = 10L, rounds = 5L, warmup = 3L) {
  for (i in seq_len(warmup)) f()
  times <- numeric(rounds)
  for (r in seq_len(rounds)) {
    times[[r]] <- system.time(for (i in seq_len(reps)) f())[["elapsed"]]
  }
  median(times) / reps
}

# ── Setup ──────────────────────────────────────────────────────────────────────

nthreads <- Sys.getenv("RAYON_NUM_THREADS", unset = "default")
cat(sprintf("immunum %s  |  R %s  |  RAYON_NUM_THREADS=%s  |  platform=%s\n",
            immunum_version(), getRversion(), nthreads, .Platform$OS.type))

SEQS <- list(
  IGH = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS",
  IGL = "SALTQPPAVSGTPGQRVTISCSGSDIGRRSVNWYQQFPGTAPKLLIYSNDQRPSVVPDRFSGSKSGTSASLAISGLQSEDEAEYYCAAWDDSLAVFGGGTQLTVGQPKA",
  TRB = "GVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVTE"
)
CHAIN_SETS <- list(
  IG  = c("IGH", "IGK", "IGL"),
  TCR = c("TRA", "TRB", "TRG", "TRD"),
  ALL = c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")
)
make_batch <- function(seq, n) rep(seq, length.out = n)
reps_for   <- function(n) max(2L, min(30L, as.integer(6000L / n)))

# ── Single sequence ────────────────────────────────────────────────────────────

cat("\n── Single sequence (Annotator$number) ──────────────────────────────────\n")
ann_all <- Annotator$new(chains = CHAIN_SETS$ALL, scheme = "IMGT")
for (nm in names(SEQS)) {
  seq <- SEQS[[nm]]
  t   <- .bench_s(function() ann_all$number(seq), reps = 500L, rounds = 5L)
  cat(sprintf("  %-8s  %.3f ms/seq\n", nm, t * 1000))
}

# ── Batch numbering ────────────────────────────────────────────────────────────

cat("\n── Batch numbering (numbering_batch) ──────────────────────────────────\n")
cat(sprintf("  %-8s  %6s  %10s  %12s\n", "chain_set", "n", "median_s", "seqs/s"))
for (cs in names(CHAIN_SETS)) {
  chains <- CHAIN_SETS[[cs]]
  for (n in c(100L, 1000L, 10000L)) {
    seqs <- make_batch(SEQS$IGH, n)
    t    <- .bench_s(function() immunum:::numbering_batch(seqs, chains, "IMGT"),
                     reps = reps_for(n), rounds = 3L)
    cat(sprintf("  %-8s  %6d  %10.4f  %12.0f\n", cs, n, t, n / t))
  }
}

# ── Batch segmentation ─────────────────────────────────────────────────────────

cat("\n── Batch segmentation (segmentation_batch) ─────────────────────────────\n")
cat(sprintf("  %6s  %10s  %12s\n", "n", "median_s", "seqs/s"))
for (n in c(100L, 1000L, 10000L)) {
  seqs <- make_batch(SEQS$IGH, n)
  t    <- .bench_s(function() immunum:::segmentation_batch(seqs, CHAIN_SETS$IG, "IMGT"),
                   reps = reps_for(n), rounds = 3L)
  cat(sprintf("  %6d  %10.4f  %12.0f\n", n, t, n / t))
}

# ── Bottleneck breakdown ───────────────────────────────────────────────────────
#
# Three potential costs inside numbering_batch:
#   A) R-side: normalize_chains / normalize_scheme
#   B) Rust: build InnerAnnotator (parse chains, load scoring matrices)
#   C) Rust: parallel alignment (the actual work)
#
# numbering_batch     = A + B + C
# numbering_batch_with = A(zero) + B(zero) + C   (annotator pre-built)
# So: B = numbering_batch - numbering_batch_with  at the same n

cat("\n── Bottleneck breakdown ─────────────────────────────────────────────────\n")

# A: pure R normalization (no Rust call at all)
t_r <- .bench_s(function() {
  immunum:::normalize_chains(CHAIN_SETS$IG)
  immunum:::normalize_scheme("IMGT")
}, reps = 10000L, rounds = 3L)
cat(sprintf("  A  R normalization only:          %.4f ms\n", t_r * 1000))

# B+C vs C: annotator construction cost
ann_ig  <- Annotator$new(chains = CHAIN_SETS$IG, scheme = "IMGT")
seqs_1k <- make_batch(SEQS$IGH, 1000L)
t_fresh <- .bench_s(function() immunum:::numbering_batch(seqs_1k, CHAIN_SETS$IG, "IMGT"),
                    reps = 10L, rounds = 5L)
t_with  <- .bench_s(function() immunum:::numbering_batch_with(seqs_1k, ann_ig),
                    reps = 10L, rounds = 5L)
cat(sprintf("  B  Annotator construction (1k):   %.4f ms  (= fresh - with)\n",
            (t_fresh - t_with) * 1000))
cat(sprintf("  C  Alignment work only (1k):      %.4f ms  (numbering_batch_with)\n",
            t_with * 1000))
cat(sprintf("     Per-sequence alignment time:   %.4f µs\n", t_with / 1000 * 1e6))

# Throughput curve — shows where rayon parallelism saturates
cat("\n── Throughput curve (IG chains, IMGT, numbering_batch_with) ────────────\n")
cat(sprintf("  %8s  %10s  %12s\n", "n", "µs/seq", "seqs/s"))
ann_ig <- Annotator$new(chains = CHAIN_SETS$IG, scheme = "IMGT")
for (n in c(1L, 5L, 20L, 100L, 500L, 2000L, 10000L)) {
  seqs <- make_batch(SEQS$IGH, n)
  reps <- max(5L, min(200L, as.integer(10000L / n)))
  t    <- .bench_s(function() immunum:::numbering_batch_with(seqs, ann_ig),
                   reps = reps, rounds = 3L)
  cat(sprintf("  %8d  %10.2f  %12.0f\n", n, t / n * 1e6, n / t))
}
cat("  (plateau = alignment-bound; steep = thread/overhead-bound)\n")

cat("\nDone.\n")
