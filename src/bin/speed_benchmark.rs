//! Speed benchmark for annotation throughput

use immunum::types::Chain;
use immunum::{read_input, Annotator, Scheme};
use std::time::Instant;

const DEFAULT_FASTA: &str = "fixtures/SRR_HKL_10k.fasta";
const ROUNDS: usize = 5;

fn main() {
    let fasta = std::env::args()
        .nth(1)
        .unwrap_or_else(|| DEFAULT_FASTA.to_string());

    let records = read_input(Some(&fasta)).expect("Failed to read input");
    let sequences: Vec<String> = records.into_iter().map(|r| r.sequence).collect();

    let chains = [Chain::IGH, Chain::IGK, Chain::IGL];
    let mut annotator = Annotator::new(&chains, Scheme::IMGT).expect("Failed to create annotator");

    // Warmup round: populate buffers and CPU caches
    let mut success = 0usize;
    let mut failed = 0usize;
    print!("  warmup...");
    for seq in &sequences {
        match annotator.number(seq) {
            Ok(_) => success += 1,
            Err(_) => failed += 1,
        }
    }
    println!(" done");

    // Benchmark over multiple rounds
    let mut timings = Vec::with_capacity(ROUNDS);

    for round in 0..ROUNDS {
        let start = Instant::now();

        for seq in &sequences {
            let _ = annotator.number(seq);
        }

        let elapsed = start.elapsed();
        let per_sec = sequences.len() as f64 / elapsed.as_secs_f64();
        println!(
            "  round {}: {:.2}s ({:.0} seq/s)",
            round + 1,
            elapsed.as_secs_f64(),
            per_sec
        );
        timings.push(elapsed.as_secs_f64());
    }

    let mean = timings.iter().sum::<f64>() / timings.len() as f64;
    let min = timings.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = timings.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mean_per_sec = sequences.len() as f64 / mean;

    println!();
    println!("{} sequences, {} rounds", sequences.len(), ROUNDS);
    println!("  mean: {:.2}s ({:.0} seq/s)", mean, mean_per_sec);
    println!("  min:  {:.2}s  max: {:.2}s", min, max);
    println!("  ok: {}  failed: {}", success, failed);
}
