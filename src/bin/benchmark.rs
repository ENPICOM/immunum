#!/usr/bin/env rustc
//! Immunum Performance Benchmark Tool
//!
//! This standalone tool benchmarks the performance of the Immunum library,
//! tracking sequential vs parallel performance, with and without prefiltering.

use clap::Parser;
use immunum::{
    annotator::Annotator,
    sequence::SequenceRecord,
    types::{CdrDefinitions, Chain, Scheme},
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Write;
use std::time::{Duration, Instant};

#[derive(Parser)]
#[command(name = "immunum-benchmark")]
#[command(about = "Performance benchmark tool for Immunum")]
#[command(version = "0.1.0")]
struct Args {
    /// Number of test sequences to generate for benchmarking
    #[arg(short = 'n', long = "sequences", default_value_t = 1000)]
    num_sequences: usize,

    /// Number of benchmark runs to average
    #[arg(short = 'r', long = "runs", default_value_t = 3)]
    runs: usize,

    /// Output file for benchmark results (JSON format)
    #[arg(short, long, default_value = "benchmark_results.json")]
    output: String,

    /// Number of threads for parallel tests
    #[arg(short = 't', long = "threads", default_value_t = num_cpus::get())]
    threads: usize,

    /// Include memory usage profiling
    #[arg(long)]
    profile_memory: bool,

    /// Test only specific modes (comma-separated: sequential,parallel,prefilter,no-prefilter)
    #[arg(long)]
    modes: Option<String>,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

#[derive(Debug, Serialize, Deserialize)]
struct BenchmarkResult {
    mode: String,
    sequences_processed: usize,
    total_chains_found: usize,
    processing_time_ms: u64,
    sequences_per_second: f64,
    memory_peak_mb: Option<u64>,
    threads_used: usize,
    prefiltering_enabled: bool,
}

#[derive(Debug, Serialize, Deserialize)]
struct BenchmarkSummary {
    timestamp: String,
    system_info: SystemInfo,
    test_parameters: TestParameters,
    results: Vec<BenchmarkResult>,
    performance_comparison: PerformanceComparison,
}

#[derive(Debug, Serialize, Deserialize)]
struct SystemInfo {
    cpu_cores: usize,
    rust_version: String,
    immunum_version: String,
}

#[derive(Debug, Serialize, Deserialize)]
struct TestParameters {
    num_sequences: usize,
    runs_averaged: usize,
    sequence_types: Vec<String>,
}

#[derive(Debug, Serialize, Deserialize)]
struct PerformanceComparison {
    sequential_vs_parallel_speedup: f64,
    prefilter_vs_no_prefilter_speedup: f64,
    memory_efficiency_rating: String,
}

// Sample sequences for testing (mix of heavy and light chains)
const SAMPLE_SEQUENCES: &[&str] = &[
    // Heavy chain sequences
    "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS",
    "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDGGYCSGGSCYFDYWGQGTLVTVSS",
    "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVTITADESTSTAYMELSSLRSEDTAVYYCARSHYGMDVWGQGTTVTVSS",
    // Light chain sequences (kappa)
    "DIQMTQSPSSLSASVGDRVTITCRASQGIRNDLGWYQQKPGKAPKRLIYDASSLESGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQSYSTPWTFGQGTKVEIK",
    "DIVMTQSHKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKALIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELK",
    "DIQMTQSPSSLSASVGDRVTITCRASQSISSWLAWYQQKPGKAPKLLIYDASSLESGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQYYSTPLTFGQGTKVEIK",
    // Light chain sequences (lambda)
    "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSKRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYAGSNNLVFGAGTKLTVL",
    "SYELTQPPSVSVSPGQTARITCSGDALPKKYAYWYQQKPGQSPVLVIYDSKRPSGIPERFSGSNSGNTATLISRVEAGDEADYYCQSADSSGTYVFGTGTKVTVL",
];

fn main() {
    let args = Args::parse();

    println!("🚀 Immunum Performance Benchmark");
    println!("================================");
    println!("Sequences: {}", args.num_sequences);
    println!("Runs: {}", args.runs);
    println!("Threads: {}", args.threads);
    println!();

    // Generate test sequences
    let sequences = generate_test_sequences(args.num_sequences);

    // Determine which modes to test
    let modes = parse_modes(&args.modes);

    // Run benchmarks
    let mut results = Vec::new();

    for mode in modes {
        if args.verbose {
            println!("🧪 Running benchmark: {}", mode);
        }

        let result = run_benchmark_mode(&mode, &sequences, args.runs, args.threads, args.verbose);
        results.push(result);
    }

    // Calculate performance comparisons
    let comparison = calculate_performance_comparison(&results);

    // Create summary
    let summary = BenchmarkSummary {
        timestamp: chrono::Utc::now().to_rfc3339(),
        system_info: SystemInfo {
            cpu_cores: num_cpus::get(),
            rust_version: rustc_version::version().unwrap().to_string(),
            immunum_version: env!("CARGO_PKG_VERSION").to_string(),
        },
        test_parameters: TestParameters {
            num_sequences: args.num_sequences,
            runs_averaged: args.runs,
            sequence_types: vec![
                "Heavy".to_string(),
                "Kappa".to_string(),
                "Lambda".to_string(),
            ],
        },
        results,
        performance_comparison: comparison,
    };

    // Write results to file
    let mut file = File::create(&args.output).expect("Failed to create output file");
    let json = serde_json::to_string_pretty(&summary).expect("Failed to serialize results");
    file.write_all(json.as_bytes())
        .expect("Failed to write results");

    // Print summary
    print_summary(&summary);
    println!("\n📊 Results written to: {}", args.output);
}

fn generate_test_sequences(count: usize) -> Vec<SequenceRecord> {
    let mut sequences = Vec::with_capacity(count);

    for i in 0..count {
        let sequence_str = SAMPLE_SEQUENCES[i % SAMPLE_SEQUENCES.len()];
        let sequence_record = SequenceRecord {
            name: format!("test_sequence_{}", i + 1).into_bytes(),
            sequence: sequence_str.as_bytes().to_vec(),
        };
        sequences.push(sequence_record);
    }

    sequences
}

fn parse_modes(modes_str: &Option<String>) -> Vec<String> {
    match modes_str {
        Some(modes) => modes.split(',').map(|s| s.trim().to_string()).collect(),
        None => vec![
            "sequential".to_string(),
            "parallel".to_string(),
            "prefilter".to_string(),
            "no-prefilter".to_string(),
        ],
    }
}

fn run_benchmark_mode(
    mode: &str,
    sequences: &[SequenceRecord],
    runs: usize,
    threads: usize,
    verbose: bool,
) -> BenchmarkResult {
    let mut total_time = Duration::new(0, 0);
    let mut total_chains = 0;

    for run in 0..runs {
        if verbose {
            println!("  Run {} of {}", run + 1, runs);
        }

        let (duration, chains_found) = run_single_benchmark(mode, sequences, threads);
        total_time += duration;
        total_chains += chains_found;
    }

    let avg_time_ms = total_time.as_millis() / runs as u128;
    let sequences_per_second = sequences.len() as f64 / (avg_time_ms as f64 / 1000.0);

    BenchmarkResult {
        mode: mode.to_string(),
        sequences_processed: sequences.len(),
        total_chains_found: total_chains / runs,
        processing_time_ms: avg_time_ms as u64,
        sequences_per_second,
        memory_peak_mb: None,
        threads_used: if mode == "sequential" { 1 } else { threads },
        prefiltering_enabled: mode != "no-prefilter",
    }
}

fn run_single_benchmark(
    mode: &str,
    sequences: &[SequenceRecord],
    threads: usize,
) -> (Duration, usize) {
    // Configure annotator based on mode
    let disable_prefiltering = mode == "no-prefilter";
    let use_parallel = mode == "parallel" || mode == "prefilter";

    let annotator = Annotator::new(
        Scheme::IMGT,
        vec![Chain::IGH, Chain::IGK, Chain::IGL],
        Some(CdrDefinitions::IMGT),
        disable_prefiltering,
        Some(0.7),
        None,
    )
    .expect("Failed to create annotator");

    // Run benchmark
    let start = Instant::now();

    let results: Vec<_> = if use_parallel {
        // Use a custom thread pool to avoid conflicts with global pool
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("Failed to build thread pool");

        pool.install(|| {
            sequences
                .par_iter()
                .map(|seq| annotator.number_sequence(seq, Some(1)))
                .collect()
        })
    } else {
        sequences
            .iter()
            .map(|seq| annotator.number_sequence(seq, Some(1)))
            .collect()
    };

    let duration = start.elapsed();

    // Count total chains found
    let total_chains = results
        .into_iter()
        .map(|result| match result {
            Ok(chains) => chains.len(),
            Err(_) => 0,
        })
        .sum();

    (duration, total_chains)
}

fn calculate_performance_comparison(results: &[BenchmarkResult]) -> PerformanceComparison {
    let sequential_speed = results
        .iter()
        .find(|r| r.mode == "sequential")
        .map(|r| r.sequences_per_second)
        .unwrap_or(0.0);

    let parallel_speed = results
        .iter()
        .find(|r| r.mode == "parallel")
        .map(|r| r.sequences_per_second)
        .unwrap_or(0.0);

    let prefilter_speed = results
        .iter()
        .find(|r| r.mode == "prefilter")
        .map(|r| r.sequences_per_second)
        .unwrap_or(0.0);

    let no_prefilter_speed = results
        .iter()
        .find(|r| r.mode == "no-prefilter")
        .map(|r| r.sequences_per_second)
        .unwrap_or(0.0);

    let sequential_vs_parallel_speedup = if sequential_speed > 0.0 {
        parallel_speed / sequential_speed
    } else {
        1.0
    };

    let prefilter_vs_no_prefilter_speedup = if no_prefilter_speed > 0.0 {
        prefilter_speed / no_prefilter_speed
    } else {
        1.0
    };

    PerformanceComparison {
        sequential_vs_parallel_speedup,
        prefilter_vs_no_prefilter_speedup,
        memory_efficiency_rating: "Good".to_string(), // TODO: Calculate based on actual memory usage
    }
}

fn print_summary(summary: &BenchmarkSummary) {
    println!("\n📈 Benchmark Results Summary");
    println!("============================");
    println!("System: {} CPU cores", summary.system_info.cpu_cores);
    println!("Rust version: {}", summary.system_info.rust_version);
    println!("Immunum version: {}", summary.system_info.immunum_version);
    println!();

    println!("Performance Results:");
    for result in &summary.results {
        println!(
            "  {}: {:.1} seq/s ({} ms total, {} chains found)",
            result.mode,
            result.sequences_per_second,
            result.processing_time_ms,
            result.total_chains_found
        );
    }

    println!();
    println!("Performance Comparison:");
    println!(
        "  Sequential vs Parallel: {:.2}x speedup",
        summary
            .performance_comparison
            .sequential_vs_parallel_speedup
    );
    println!(
        "  Prefilter vs No-prefilter: {:.2}x speedup",
        summary
            .performance_comparison
            .prefilter_vs_no_prefilter_speedup
    );
}
