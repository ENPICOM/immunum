use immunum::annotator::Annotator;
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
    successful_matches: usize,
    failed_matches: usize,
    parsing_errors: usize,
    numbering_errors: usize,
    execution_time_ms: u128,
}

impl ValidationStats {
    fn success_rate(&self) -> f64 {
        if self.total_sequences == 0 {
            0.0
        } else {
            (self.successful_matches as f64 / self.total_sequences as f64) * 100.0
        }
    }

    fn print_summary(&self) {
        println!("\n=== AbPdSeq Validation Results ===");
        println!("Total sequences: {}", self.total_sequences);
        println!("Successful matches: {}", self.successful_matches);
        println!("Failed matches: {}", self.failed_matches);
        println!("Parsing errors: {}", self.parsing_errors);
        println!("Numbering errors: {}", self.numbering_errors);
        println!("Success rate: {:.2}%", self.success_rate());
        println!("Execution time: {} ms", self.execution_time_ms);
        println!("=====================================");
    }
}

/// Parse FASTA file and extract sequences with their metadata
fn parse_fasta_file(
    file_path: &str,
) -> Result<Vec<(String, String, Chain)>, Box<dyn std::error::Error>> {
    let content = fs::read_to_string(file_path)?;
    let mut sequences = Vec::new();
    let mut current_header = String::new();
    let mut current_sequence = String::new();

    for line in content.lines() {
        if line.starts_with('>') {
            // Process previous sequence if exists
            if !current_header.is_empty() && !current_sequence.is_empty() {
                let chain = parse_chain_from_header(&current_header)?;
                sequences.push((current_header.clone(), current_sequence.clone(), chain));
            }

            // Start new sequence
            current_header = line[1..].to_string(); // Remove '>'
            current_sequence.clear();
        } else {
            current_sequence.push_str(line.trim());
        }
    }

    // Process last sequence
    if !current_header.is_empty() && !current_sequence.is_empty() {
        let chain = parse_chain_from_header(&current_header)?;
        sequences.push((current_header, current_sequence, chain));
    }

    Ok(sequences)
}

/// Parse chain type from FASTA header
/// For light chains, we return IGK as default but the annotator will determine the actual type
fn parse_chain_from_header(header: &str) -> Result<Chain, Box<dyn std::error::Error>> {
    // Header format: "3tje_L|Light|L" or similar
    if header.contains("Heavy") || header.ends_with("|H") {
        Ok(Chain::IGH)
    } else if header.contains("Light") || header.ends_with("|L") || header.ends_with("|M") {
        // For light chains, we return IGK as default, but the annotator will determine
        // whether it's actually IGK (Kappa) or IGL (Lambda)
        Ok(Chain::IGK)
    } else {
        Err(format!("Unable to determine chain type from header: {}", header).into())
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
    sequences: &[(String, String, Chain)],
    expected_results: &HashMap<String, Vec<ExpectedResult>>,
    max_count: Option<usize>,
) -> ValidationStats {
    let start_time = Instant::now();
    let mut stats = ValidationStats::default();

    let sequences_to_process = if let Some(max) = max_count {
        &sequences[..sequences.len().min(max)]
    } else {
        sequences
    };

    stats.total_sequences = sequences_to_process.len();

    // Extract just the sequence strings for batch processing
    let sequence_strings: Vec<String> = sequences_to_process
        .iter()
        .map(|(_, seq, _)| seq.clone())
        .collect();

    println!(
        "Processing {} sequences in parallel...",
        sequence_strings.len()
    );

    // Process all sequences in parallel
    let results = annotator.number_sequences(&sequence_strings, true);
    stats.execution_time_ms = start_time.elapsed().as_millis();
    
    // Validate each result
    for (i, result) in results.into_iter().enumerate() {
        let (sequence_id, _, _) = &sequences_to_process[i];

        match result {
            Ok(annotation_result) => {
                if let Some(expected) = expected_results.get(sequence_id) {
                    // Convert our result to comparable format
                    let our_numbers = annotation_result.numbers;
                    let our_chain_type = match annotation_result.chain {
                        Chain::IGH => "H",
                        Chain::IGK => "K",
                        Chain::IGL => "L",
                        _ => "unknown",
                    };

                    // Find matching expected result based on chain type
                    let mut found_match = false;
                    for exp in expected {
                        if exp.chain_type == our_chain_type {
                            if our_numbers == exp.numbers {
                                stats.successful_matches += 1;
                                found_match = true;
                                break;
                            } else {
                                // Print detailed mismatch for debugging (only first few)
                                if stats.failed_matches < 5 {
                                    println!(
                                                "MISMATCH for {}: Expected {} numbers, got {}. First 10 expected: {:?}, First 10 actual: {:?}",
                                                sequence_id,
                                                exp.numbers.len(),
                                                our_numbers.len(),
                                                exp.numbers.iter().take(10).collect::<Vec<_>>(),
                                                our_numbers.iter().take(10).collect::<Vec<_>>()
                                            );
                                }
                                stats.failed_matches += 1;
                                found_match = true;
                                break;
                            }
                        }
                    }

                    if !found_match {
                        if stats.parsing_errors < 5 {
                            println!(
                                "No expected result found for chain type '{}' in sequence '{}'",
                                our_chain_type, sequence_id
                            );
                        }
                        stats.parsing_errors += 1;
                    }
                } else {
                    if stats.parsing_errors < 5 {
                        println!("No expected result found for sequence: {}", sequence_id);
                    }
                    stats.parsing_errors += 1;
                }
            }
            Err(e) => {
                if stats.numbering_errors < 5 {
                    println!("Error numbering sequence '{}': {}", sequence_id, e);
                }
                stats.numbering_errors += 1;
            }
        }
    }


    stats
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_abpdseq_validation_subset() {
        // Test with first 10 sequences for quick validation
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            None,
            Some(true), // Enable prefiltering to help with chain detection
        )
        .expect("Failed to create annotator");

        let sequences =
            parse_fasta_file("fixtures/abpdseq_agreed.fasta").expect("Failed to parse FASTA file");

        let expected_results = load_expected_results("fixtures/abpdseq_agreed_output.json")
            .expect("Failed to load expected results");

        println!(
            "Loaded {} sequences and {} expected results",
            sequences.len(),
            expected_results.len()
        );

        let stats = run_validation_parallel(&annotator, &sequences, &expected_results, Some(500));
        stats.print_summary();

        // Assert that at least some sequences match (allowing for some expected differences)
        // This is a sanity check rather than requiring 100% accuracy
        assert!(
            stats.successful_matches > 0,
            "No sequences matched - there might be a systematic issue"
        );
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
            None,
            Some(true),
        )
        .expect("Failed to create annotator");

        let sequences =
            parse_fasta_file("fixtures/abpdseq_agreed.fasta").expect("Failed to parse FASTA file");

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
