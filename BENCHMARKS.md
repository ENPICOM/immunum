# Validation Benchmarks

This file tracks accuracy metrics across all supported chains. Metrics are generated from validation datasets and updated via `cargo run --release --bin benchmark`.

**Last Updated**: 2026-02-20

**Execution Time**: 1.66s

## IMGT Summary

| Chain | Total Sequences | Perfect Sequences | Perfect % | Overall Accuracy | Correct Positions | Total Positions |
|-------|-----------------|-------------------|-----------|------------------|-------------------|-----------------|
| IGH   |            2452 |              2449 |    99.88% |           99.99% |            295875 |          295901 |
| IGK   |            1485 |              1480 |    99.66% |           99.98% |            160921 |          160947 |
| IGL   |             371 |               369 |    99.46% |           99.99% |             40288 |           40293 |
| TRA   |             865 |               761 |    87.98% |           99.36% |             94829 |           95438 |
| TRB   |             934 |               904 |    96.79% |           99.70% |            105008 |          105325 |
| TRG   |              25 |                25 |   100.00% |          100.00% |              2799 |            2799 |
| TRD   |              23 |                23 |   100.00% |          100.00% |              2622 |            2622 |

## Kabat Summary

| Chain | Total Sequences | Perfect Sequences | Perfect % | Overall Accuracy | Correct Positions | Total Positions |
|-------|-----------------|-------------------|-----------|------------------|-------------------|-----------------|
| IGH   |            2452 |              2449 |    99.88% |           99.99% |            295876 |          295901 |
| IGK   |            1485 |              1478 |    99.53% |           99.98% |            160920 |          160947 |
| IGL   |             371 |               369 |    99.46% |           99.99% |             40288 |           40293 |

## Quality Thresholds

The test suite enforces these minimum thresholds:
- **Perfect sequence accuracy**: ≥99% (sequences with 100% correct positions)
- **Overall position accuracy**: ≥99% (all positions across all sequences)

## Workflow

To update these metrics:
```bash
# Run all validation tests
cargo test

# Generate updated benchmark report (release mode for accurate timing)
cargo run --quiet --release --bin benchmark > BENCHMARKS.md
```
