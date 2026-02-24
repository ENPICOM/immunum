# Validation Benchmarks

This file tracks accuracy metrics across all supported chains. Metrics are generated from validation datasets and updated via `cargo run --release --bin benchmark`.

**Last Updated**: 2026-02-24

**Execution Time**: 1.90s

## IMGT Summary

| Chain | Total Sequences | Perfect Sequences | Perfect % | Overall Accuracy |
|-------|-----------------|-------------------|-----------|------------------|
| IGH   |            2452 |              2449 |    99.88% |           99.99% |
| IGK   |            1485 |              1480 |    99.66% |           99.98% |
| IGL   |             371 |               369 |    99.46% |           99.99% |
| TRA   |             865 |               757 |    87.51% |           99.29% |
| TRB   |             934 |               904 |    96.79% |           99.70% |
| TRG   |              25 |                25 |   100.00% |          100.00% |
| TRD   |              23 |                23 |   100.00% |          100.00% |

## Kabat Summary

| Chain | Total Sequences | Perfect Sequences | Perfect % | Overall Accuracy |
|-------|-----------------|-------------------|-----------|------------------|
| IGH   |            2452 |              2449 |    99.88% |           99.99% |
| IGK   |            1485 |              1478 |    99.53% |           99.98% |
| IGL   |             371 |               369 |    99.46% |           99.99% |

## Metrics

- **Perfect %**: share of sequences where every residue position matches the reference exactly.
- **Overall Accuracy**: fraction of individual residue positions that match across all sequences.

The test suite enforces minimum thresholds of ≥99% for both metrics.

## Workflow

To update these metrics:
```bash
# Run all validation tests
cargo test

# Generate updated benchmark report (release mode for accurate timing)
cargo run --quiet --release --bin benchmark > BENCHMARKS.md
```
