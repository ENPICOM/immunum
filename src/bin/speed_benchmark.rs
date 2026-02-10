//! Speed benchmark for annotation throughput

use immunum2::types::Chain;
use immunum2::validation::load_validation_csv;
use immunum2::{Annotator, Scheme};
use std::path::Path;
use std::time::Instant;

fn main() {
    let csv_path = "fixtures/validation/ab_H_imgt.csv";
    let path = Path::new(csv_path);

    let sequences: Vec<String> = load_validation_csv(path)
        .expect("Failed to load validation CSV")
        .into_iter()
        .map(|e| e.sequence)
        .collect();

    let annotator =
        Annotator::new(&[Chain::IGH], Scheme::IMGT).expect("Failed to create annotator");

    // Warm up cache by annotating 100 sequences first (to load scoring matrices, etc.)
    for seq in sequences.iter().take(100) {
        let _ = annotator.annotate(seq);
    }
    
    // Benchmark
    let start = Instant::now();

    for seq in &sequences {
        // annotate and number the output alignment with imgt scheme
        match annotator.annotate(seq) {
            Ok(result) => {
                let _ = result.numbering(Scheme::IMGT);
            }
            Err(e) => eprintln!("Failed to annotate sequence: {}", e),
        }

    }

    let elapsed = start.elapsed();
    let per_sec = sequences.len() as f64 / elapsed.as_secs_f64();

    println!(
        "{} sequences in {:.2}ms ({:.0} seq/s)",
        sequences.len(),
        elapsed.as_secs_f64() * 1000.0,
        per_sec,
    );
}
