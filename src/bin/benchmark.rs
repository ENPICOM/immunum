//! Generate validation benchmark report

use chrono::Local;
use immunum::types::{Chain, Scheme};
use immunum::validation::{validate_chain_with_scheme, ChainMetrics};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

const BENCHMARKS_FILE: &str = "BENCHMARKS.toml";

#[derive(Debug, Serialize, Deserialize)]
struct BenchmarkReport {
    last_updated: String,
    execution_time_secs: f64,
    #[serde(default)]
    imgt: BTreeMap<String, ChainEntry>,
    #[serde(default)]
    kabat: BTreeMap<String, ChainEntry>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct ChainEntry {
    total_sequences: usize,
    perfect_sequences: usize,
    perfect_pct: f64,
    overall_accuracy: f64,
}

impl From<&ChainMetrics> for ChainEntry {
    fn from(m: &ChainMetrics) -> Self {
        Self {
            total_sequences: m.total_sequences,
            perfect_sequences: m.perfect_sequences,
            perfect_pct: round2(m.perfect_percentage()),
            overall_accuracy: round2(m.overall_accuracy()),
        }
    }
}

fn round2(v: f64) -> f64 {
    (v * 100.0).round() / 100.0
}

fn main() {
    let start_time = Local::now();

    let previous: Option<BenchmarkReport> = fs::read_to_string(BENCHMARKS_FILE)
        .ok()
        .and_then(|c| toml::from_str(&c).ok());

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
    let elapsed_secs = elapsed.num_milliseconds() as f64 / 1000.0;

    let report = BenchmarkReport {
        last_updated: end_time.format("%Y-%m-%d").to_string(),
        execution_time_secs: round2(elapsed_secs),
        imgt: metrics_to_map(&imgt_metrics),
        kabat: metrics_to_map(&kabat_metrics),
    };

    // Print comparison to terminal
    print_comparison(&imgt_metrics, "IMGT", previous.as_ref().map(|p| &p.imgt));
    print_comparison(&kabat_metrics, "Kabat", previous.as_ref().map(|p| &p.kabat));

    // Write TOML
    let toml_str = toml::to_string_pretty(&report).expect("Failed to serialize TOML");
    fs::write(BENCHMARKS_FILE, &toml_str).expect("Failed to write BENCHMARKS.toml");

    println!("Wrote {} ({:.2}s)", BENCHMARKS_FILE, elapsed_secs);
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

fn metrics_to_map(metrics: &[ChainMetrics]) -> BTreeMap<String, ChainEntry> {
    metrics
        .iter()
        .map(|m| (m.chain.to_string(), ChainEntry::from(m)))
        .collect()
}

fn format_delta(new: f64, old: f64) -> String {
    let diff = new - old;
    if diff.abs() < 0.005 {
        String::new()
    } else if diff > 0.0 {
        format!("(+{:.2}%)", diff)
    } else {
        format!("({:.2}%)", diff)
    }
}

fn print_comparison(
    metrics: &[ChainMetrics],
    scheme: &str,
    previous: Option<&BTreeMap<String, ChainEntry>>,
) {
    println!("{} Summary:", scheme);
    println!(
        "  {:5}  {:>10}  {:>10}  {:>12}  {:>12}",
        "Chain", "Perfect", "Accuracy", "Perfect Δ", "Accuracy Δ"
    );

    for m in metrics {
        let chain_name = m.chain.to_string();
        let perfect = m.perfect_percentage();
        let accuracy = m.overall_accuracy();

        let (p_delta, a_delta) = match previous.and_then(|p| p.get(&chain_name)) {
            Some(old) => (
                format_delta(perfect, old.perfect_pct),
                format_delta(accuracy, old.overall_accuracy),
            ),
            None => ("(new)".to_string(), "(new)".to_string()),
        };

        println!(
            "  {:5}  {:>9.2}%  {:>9.2}%  {:>12}  {:>12}",
            chain_name, perfect, accuracy, p_delta, a_delta
        );
    }
    println!();
}
