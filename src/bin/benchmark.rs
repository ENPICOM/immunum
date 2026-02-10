//! Generate validation benchmark report

use chrono::Local;
use immunum2::types::Chain;
use immunum2::validation::{validate_chain, ChainMetrics};
use std::path::PathBuf;

fn main() {
    let start_time = Local::now();

    let chains = vec![
        (Chain::IGH, "fixtures/validation/ab_H_imgt.csv"),
        (Chain::IGK, "fixtures/validation/ab_K_imgt.csv"),
        (Chain::IGL, "fixtures/validation/ab_L_imgt.csv"),
        (Chain::TRA, "fixtures/validation/tcr_A_imgt.csv"),
        (Chain::TRB, "fixtures/validation/tcr_B_imgt.csv"),
        (Chain::TRG, "fixtures/validation/tcr_G_imgt.csv"),
        (Chain::TRD, "fixtures/validation/tcr_D_imgt.csv"),
    ];

    let mut all_metrics = Vec::new();

    // Collect metrics for each chain
    for (chain, csv_path) in &chains {
        let path = PathBuf::from(csv_path);
        if !path.exists() {
            eprintln!("Warning: Skipping {} - file not found: {}", chain, csv_path);
            continue;
        }

        match validate_chain(*chain, csv_path) {
            Ok(metrics) => all_metrics.push(metrics),
            Err(e) => {
                eprintln!("Error validating {}: {}", chain, e);
                continue;
            }
        }
    }

    let end_time = Local::now();
    let elapsed = end_time.signed_duration_since(start_time);

    // Generate report
    print_benchmark_report(&all_metrics, &end_time, elapsed);
}

fn print_benchmark_report(
    metrics: &[ChainMetrics],
    timestamp: &chrono::DateTime<Local>,
    elapsed: chrono::Duration,
) {
    let date = timestamp.format("%Y-%m-%d").to_string();

    println!("# Validation Benchmarks");
    println!();
    println!("This file tracks accuracy metrics across all supported chains. Metrics are generated from validation datasets and updated via `cargo run --release --bin benchmark`.");
    println!();
    println!("**Last Updated**: {}", date);
    println!();
    println!(
        "**Execution Time**: {:.2}s",
        elapsed.num_milliseconds() as f64 / 1000.0
    );
    println!();

    // Overall summary table
    println!("## Overall Summary");
    println!();
    println!("| Chain | Total Sequences | Perfect Sequences | Perfect % | Overall Accuracy | Correct Positions | Total Positions |");
    println!("|-------|-----------------|-------------------|-----------|------------------|-------------------|-----------------|");

    for m in metrics {
        println!(
            "| {:5} | {:15} | {:17} | {:8.2}% | {:15.2}% | {:17} | {:15} |",
            m.chain.to_string(),
            m.total_sequences,
            m.perfect_sequences,
            m.perfect_percentage(),
            m.overall_accuracy(),
            m.correct_positions,
            m.total_positions
        );
    }

    println!();
    println!("## Quality Thresholds");
    println!();
    println!("The test suite enforces these minimum thresholds:");
    println!("- **Perfect sequence accuracy**: ≥99% (sequences with 100% correct positions)");
    println!("- **Overall position accuracy**: ≥99% (all positions across all sequences)");
    println!();

    println!("## Workflow");
    println!();
    println!("To update these metrics:");
    println!("```bash");
    println!("# Run all validation tests");
    println!("cargo test");
    println!();
    println!("# Generate updated benchmark report (release mode for accurate timing)");
    println!("cargo run --quiet --release --bin benchmark > BENCHMARKS.md");
    println!("```");
}
