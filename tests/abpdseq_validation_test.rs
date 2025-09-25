use immunum::annotator::Annotator;
use immunum::sequence::SequenceRecord;
use immunum::types::{Chain, Scheme};
use serde_json::Value;
use std::collections::HashMap;
use std::fs;
use std::time::Instant;

/// Structure to hold expected numbering result
#[derive(Debug, Clone)]
struct ExpectedResult {
    chain_type: String,
    numbers: Vec<String>,
}

/// Structure to hold validation statistics
#[derive(Debug, Default)]
struct ValidationStats {
    total_sequences: usize,
    successful_sequence_matches: usize,
    total_paired: usize,
    successful_paired_matches: usize,
    numbering_errors: usize,
    execution_time_ms: u128,
}

impl ValidationStats {
    fn success_rate(&self) -> f64 {
        if self.total_sequences == 0 {
            0.0
        } else {
            (self.successful_sequence_matches as f64 / self.total_sequences as f64) * 100.0
        }
    }

    fn print_summary(&self) {
        println!("\n=== AbPdSeq Validation Results ===");
        println!("Total sequences: {}", self.total_sequences);
        println!("  of which paired sequences: {}", self.total_paired);
        println!(
            "Successful sequence matches: {} ({:.2}%)",
            self.successful_sequence_matches,
            self.success_rate()
        );
        println!(
            "  of which paired matches: {} ({:.2}%)",
            self.successful_paired_matches,
            (self.successful_paired_matches as f64 / self.total_paired as f64) * 100.0
        );
        println!("Numbering errors: {}", self.numbering_errors);
        println!("Execution time: {} ms", self.execution_time_ms);
        println!("=====================================");
    }
}

/// Load expected results from JSON file
fn load_expected_results(
    file_path: &str,
) -> Result<HashMap<String, Vec<ExpectedResult>>, Box<dyn std::error::Error>> {
    let content = fs::read_to_string(file_path)?;
    let json: Value = serde_json::from_str(&content)?;
    let mut results = HashMap::new();

    if let Value::Object(map) = json {
        for (sequence_id, data) in map {
            if let Value::Array(entries) = data {
                let mut expected_results = Vec::new();

                for entry in entries {
                    if let Value::Object(obj) = entry {
                        let chain_type = obj
                            .get("chain_type")
                            .and_then(|v| v.as_str())
                            .unwrap_or("unknown")
                            .to_string();

                        let numbers = obj
                            .get("numbers")
                            .and_then(|v| v.as_array())
                            .map(|arr| {
                                arr.iter()
                                    .filter_map(|v| v.as_str().map(|s| s.to_string()))
                                    .collect()
                            })
                            .unwrap_or_default();

                        expected_results.push(ExpectedResult {
                            chain_type,
                            numbers,
                        });
                    }
                }

                results.insert(sequence_id, expected_results);
            }
        }
    }

    Ok(results)
}

/// Run validation using parallel batch processing
fn run_validation_parallel(
    annotator: &Annotator,
    sequences: &[SequenceRecord],
    expected_results: &HashMap<String, Vec<ExpectedResult>>,
    max_count: Option<usize>,
) -> ValidationStats {
    let mut stats = ValidationStats::default();

    let sequences_to_process = if let Some(max) = max_count {
        &sequences[..sequences.len().min(max)]
    } else {
        sequences
    };

    stats.total_sequences = sequences_to_process.len();
    stats.total_paired = expected_results.values().filter(|v| v.len() == 2).count();

    println!(
        "Processing {} sequences in parallel...",
        sequences_to_process.len()
    );

    // Process all sequences in parallel
    let start_time = Instant::now();
    let results = annotator.number(sequences_to_process.to_vec(), Some(2));
    stats.execution_time_ms = start_time.elapsed().as_millis();

    // Validate each result
    for (header, result) in results.into_iter() {
        match result {
            Ok(annotation_result) => {
                if let Some(expected) = expected_results.get(&header) {
                    let mut matches = Vec::new();
                    // Find matching expected result based on chain type
                    let exp_chain_1 = expected.get(0).expect("");
                    let exp_chain_2 = expected.get(1);

                    // Did we find chain1?
                    let own_chain_1_opt = annotation_result
                        .iter()
                        .find(|chain| chain.chain.to_short() == exp_chain_1.chain_type);

                    if let Some(own_chain_1) = own_chain_1_opt {
                        if own_chain_1.numbers == exp_chain_1.numbers {
                            matches.push(own_chain_1)
                        } else {
                            stats.numbering_errors += 1
                        }
                    }

                    // Did we find chain2 (if there is one)
                    if let Some(exp_chain_2) = exp_chain_2 {
                        let own_chain_2_opt = annotation_result
                            .iter()
                            .find(|chain| chain.chain.to_short() == exp_chain_2.chain_type);

                        if let Some(own_chain_2) = own_chain_2_opt {
                            if own_chain_2.numbers == exp_chain_2.numbers {
                                matches.push(own_chain_2)
                            } else {
                                stats.numbering_errors += 1
                            }
                        }
                    }
                    match (matches.len(), expected.len()) {
                        (1, 1) => stats.successful_sequence_matches += 1,
                        (2, 2) => {
                            stats.successful_paired_matches += 1;
                            stats.successful_sequence_matches += 1;
                        }
                        (_, _) => (),
                    }
                }
            }
            Err(e) => {
                if stats.numbering_errors < 5 {
                    println!("Error numbering sequence '{}': {}", header, e);
                }
                stats.numbering_errors += 1;
            }
        }
    }

    stats
}

#[cfg(test)]
mod tests {
    use immunum::sequence::SequenceStream;

    use super::*;

    #[test]
    fn test_abpdseq_validation_subset() {
        // Test with first 50 sequences for quick validation
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            true, // Enable prefiltering to help with chain detection
        )
        .expect("Failed to create annotator");

        let sequence_stream = SequenceStream::new("fixtures/abpdseq_agreed.fasta")
            .expect("Failed to create sequence stream");
        let sequences: Result<Vec<_>, _> = sequence_stream.collect();
        let sequences = match sequences {
            Ok(sequences) => sequences,
            Err(e) => {
                eprintln!("Error collecting sequences: {e}");
                std::process::exit(1);
            }
        };

        let expected_results = load_expected_results("fixtures/abpdseq_agreed_output.json")
            .expect("Failed to load expected results");

        println!(
            "Loaded {} sequences and {} expected results",
            sequences.len(),
            expected_results.len()
        );

        let stats = run_validation_parallel(&annotator, &sequences, &expected_results, Some(50));
        stats.print_summary();

        assert!(
            stats.success_rate() >= 99.0,
            "Success rate too low: {:.2}%",
            stats.success_rate()
        );
    }

    #[test]
    #[ignore] // Ignore by default due to long runtime - run with `cargo test --ignored`
    fn test_abpdseq_validation_full() {
        // Full validation - will take much longer
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            false,
        )
        .expect("Failed to create annotator");

        let sequence_stream = SequenceStream::new("fixtures/abpdseq_agreed.fasta")
            .expect("Failed to create sequence stream");
        let sequences: Result<Vec<_>, _> = sequence_stream.collect();
        let sequences = match sequences {
            Ok(sequences) => sequences,
            Err(e) => {
                eprintln!("Error collecting sequences: {e}");
                std::process::exit(1);
            }
        };

        let expected_results = load_expected_results("fixtures/abpdseq_agreed_output.json")
            .expect("Failed to load expected results");

        let stats = run_validation_parallel(&annotator, &sequences, &expected_results, None);
        stats.print_summary();

        // For full validation, we expect higher accuracy
        assert!(
            stats.success_rate() >= 99.0,
            "Success rate too low for full validation: {:.2}%",
            stats.success_rate()
        );
    }
}
