//! Generate validation benchmark report

use chrono::Local;
use immunum2::types::{Chain, Scheme};
use immunum2::validation::{validate_chain_with_scheme, ChainMetrics};
use std::path::PathBuf;

fn main() {
    let start_time = Local::now();

    let imgt_chains = vec![
        (Chain::IGH, "fixtures/validation/ab_H_imgt.csv"),
        (Chain::IGK, "fixtures/validation/ab_K_imgt.csv"),
        (Chain::IGL, "fixtures/validation/ab_L_imgt.csv"),
        (Chain::TRA, "fixtures/validation/tcr_A_imgt.csv"),
        (Chain::TRB, "fixtures/validation/tcr_B_imgt.csv"),
        (Chain::TRG, "fixtures/validation/tcr_G_imgt.csv"),
        (Chain::TRD, "fixtures/validation/tcr_D_imgt.csv"),
    ];

    let kabat_chains = vec![
        (Chain::IGH, "fixtures/validation/ab_H_kabat.csv"),
        (Chain::IGK, "fixtures/validation/ab_K_kabat.csv"),
        (Chain::IGL, "fixtures/validation/ab_L_kabat.csv"),
    ];

    let imgt_metrics = collect_metrics(&imgt_chains, Scheme::IMGT);
    let kabat_metrics = collect_metrics(&kabat_chains, Scheme::Kabat);

    let end_time = Local::now();
    let elapsed = end_time.signed_duration_since(start_time);

    print_benchmark_report(&imgt_metrics, &kabat_metrics, &end_time, elapsed);
}

fn collect_metrics(chains: &[(Chain, &str)], scheme: Scheme) -> Vec<ChainMetrics> {
    let mut all_metrics = Vec::new();
    for (chain, csv_path) in chains {
        let path = PathBuf::from(csv_path);
        if !path.exists() {
            eprintln!("Warning: Skipping {} - file not found: {}", chain, csv_path);
            continue;
        }
        match validate_chain_with_scheme(*chain, csv_path, scheme, None) {
            Ok(metrics) => all_metrics.push(metrics),
            Err(e) => eprintln!("Error validating {}: {}", chain, e),
        }
    }
    all_metrics
}

fn print_benchmark_report(
    imgt_metrics: &[ChainMetrics],
    kabat_metrics: &[ChainMetrics],
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

    println!("## IMGT Summary");
    println!();
    print_metrics_table(imgt_metrics);

    println!("## Kabat Summary");
    println!();
    print_metrics_table(kabat_metrics);

    println!("## Metrics");
    println!();
    println!("- **Perfect %**: share of sequences where every residue position matches the reference exactly.");
    println!("- **Overall Accuracy**: fraction of individual residue positions that match across all sequences.");
    println!();
    println!("The test suite enforces minimum thresholds of ≥99% for both metrics.");
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

fn print_metrics_table(metrics: &[ChainMetrics]) {
    println!("| Chain | Total Sequences | Perfect Sequences | Perfect % | Overall Accuracy |");
    println!("|-------|-----------------|-------------------|-----------|------------------|");
    for m in metrics {
        println!(
            "| {:5} | {:15} | {:17} | {:8.2}% | {:15.2}% |",
            m.chain.to_string(),
            m.total_sequences,
            m.perfect_sequences,
            m.perfect_percentage(),
            m.overall_accuracy(),
        );
    }
    println!();
}
