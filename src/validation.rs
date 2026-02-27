//! Validation against numbered test data

use crate::annotator::Annotator;
use crate::error::{Error, Result};
use crate::types::{Chain, Position, Scheme};
use std::collections::HashMap;
use std::fs;
use std::path::Path;

/// A single validation entry with expected numbering
#[derive(Debug, Clone)]
pub struct ValidationEntry {
    pub header: String,
    pub sequence: String,
    pub species: String,
    pub expected_positions: Vec<(Position, char)>,
}

/// Load validation data from CSV file
pub fn load_validation_csv(path: &Path) -> Result<Vec<ValidationEntry>> {
    let content = fs::read_to_string(path)
        .map_err(|e| Error::ConsensusParseError(format!("Failed to read file: {}", e)))?;

    let mut entries = Vec::new();
    let mut lines = content.lines();

    // Parse header to get position labels
    let header_line = lines
        .next()
        .ok_or_else(|| Error::ConsensusParseError("Empty CSV file".to_string()))?;
    let headers: Vec<&str> = header_line.split(',').collect();

    // Skip first three columns (header, sequence, species)
    let position_headers: Vec<Position> = headers[3..]
        .iter()
        .filter_map(|h| h.parse::<Position>().ok())
        .collect();

    // Parse data lines
    for line in lines {
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 3 {
            continue;
        }

        let header = fields[0].to_string();
        let sequence = fields[1].to_string();
        let species = fields[2].to_string();

        // Parse expected positions (skip empty fields)
        let mut expected_positions = Vec::new();
        for (i, field) in fields[3..].iter().enumerate() {
            if !field.is_empty() && i < position_headers.len() {
                let pos = position_headers[i];
                let aa = field.chars().next().unwrap();
                expected_positions.push((pos, aa));
            }
        }

        entries.push(ValidationEntry {
            header,
            sequence,
            species,
            expected_positions,
        });
    }

    Ok(entries)
}

/// Validation result for a single sequence
#[derive(Debug, Clone)]
pub struct ValidationResult {
    pub header: String,
    pub sequence: String,
    pub detected_chain: Chain,
    pub numbering: Vec<Position>,
    pub total_positions: usize,
    pub correct_positions: usize,
    pub incorrect_positions: usize,
    pub missing_positions: usize,
    pub extra_positions: usize,
    pub mismatches: Vec<(usize, Position, Position)>, // (seq_idx, expected, got)
    pub alignment_confidence: f32,
}

impl ValidationResult {
    pub fn accuracy(&self) -> f32 {
        if self.total_positions == 0 {
            return 0.0;
        }
        self.correct_positions as f32 / self.total_positions as f32
    }

    pub fn is_perfect(&self) -> bool {
        self.accuracy() == 1.0
    }
}

/// Aggregated metrics for a chain
#[derive(Debug, Clone)]
pub struct ChainMetrics {
    pub chain: Chain,
    pub scheme: Scheme,
    pub csv_path: String,
    pub total_sequences: usize,
    pub perfect_sequences: usize,
    pub total_positions: usize,
    pub correct_positions: usize,
}

impl ChainMetrics {
    pub fn new(chain: Chain, scheme: Scheme, csv_path: String) -> Self {
        Self {
            chain,
            scheme,
            csv_path,
            total_sequences: 0,
            perfect_sequences: 0,
            total_positions: 0,
            correct_positions: 0,
        }
    }

    pub fn add_result(&mut self, result: &ValidationResult) {
        self.total_sequences += 1;
        self.total_positions += result.total_positions;
        self.correct_positions += result.correct_positions;
        if result.is_perfect() {
            self.perfect_sequences += 1;
        }
    }

    pub fn perfect_percentage(&self) -> f64 {
        if self.total_sequences == 0 {
            return 0.0;
        }
        (self.perfect_sequences as f64 / self.total_sequences as f64) * 100.0
    }

    pub fn overall_accuracy(&self) -> f64 {
        if self.total_positions == 0 {
            return 0.0;
        }
        (self.correct_positions as f64 / self.total_positions as f64) * 100.0
    }
}

/// Validate all sequences for a given chain with IMGT scheme (default)
pub fn validate_chain(chain: Chain, csv_path: &str) -> Result<ChainMetrics> {
    validate_chain_with_scheme(chain, csv_path, Scheme::IMGT, None)
}

/// Validate all sequences for a given chain with a specific scheme and collect metrics
/// Optionally filter by species (e.g., "human", "mouse")
pub fn validate_chain_with_scheme(
    chain: Chain,
    csv_path: &str,
    scheme: Scheme,
    species_filter: Option<&str>,
) -> Result<ChainMetrics> {
    let path = std::path::PathBuf::from(csv_path);
    let entries = load_validation_csv(&path)?;
    let annotator = Annotator::new(&[chain], scheme)?;

    let mut metrics = ChainMetrics::new(chain, scheme, csv_path.to_string());

    for entry in entries.iter() {
        // Apply species filter if provided
        if let Some(species) = species_filter {
            if !entry.species.eq_ignore_ascii_case(species) {
                continue;
            }
        }

        let result = validate_entry(entry, &annotator, scheme)?;
        metrics.add_result(&result);
    }

    Ok(metrics)
}

/// Validate a single entry against the annotator
pub fn validate_entry(
    entry: &ValidationEntry,
    annotator: &Annotator,
    scheme: Scheme,
) -> Result<ValidationResult> {
    // Annotate the sequence
    let result = annotator.annotate(&entry.sequence)?;
    let numbering = result.numbering(scheme);

    // The expected_positions are already aligned to the sequence (no gaps)
    // Each entry is (Position, amino_acid) for each residue
    // Numbering should match 1:1 with the sequence length

    if numbering.len() != entry.sequence.len() {
        return Err(Error::AlignmentError(format!(
            "Numbering length {} doesn't match sequence length {}. Alignment: {:?}",
            numbering.len(),
            entry.sequence.len(),
            result.alignment,
        )));
    }

    let total_positions = entry.expected_positions.len();
    let mut correct_positions = 0;
    let mut incorrect_positions = 0;
    let mut mismatches = Vec::new();

    // Build a map from sequence index to expected position
    let mut expected_by_idx: HashMap<usize, &Position> = HashMap::new();

    for (seq_idx, (expected_pos, _aa)) in entry.expected_positions.iter().enumerate() {
        expected_by_idx.insert(seq_idx, expected_pos);
    }

    // Compare positions
    for (idx, actual_pos) in numbering.iter().enumerate() {
        if let Some(&expected_pos) = expected_by_idx.get(&idx) {
            if actual_pos == expected_pos {
                correct_positions += 1;
            } else {
                incorrect_positions += 1;
                mismatches.push((idx, *expected_pos, *actual_pos));
            }
        }
    }

    let missing_positions = total_positions.saturating_sub(correct_positions + incorrect_positions);
    let extra_positions = 0; // We now ensure lengths match

    Ok(ValidationResult {
        header: entry.header.clone(),
        sequence: entry.sequence.clone(),
        detected_chain: result.chain,
        numbering,
        total_positions,
        correct_positions,
        incorrect_positions,
        missing_positions,
        extra_positions,
        mismatches,
        alignment_confidence: result.confidence,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_load_validation_csv() {
        let path = PathBuf::from("fixtures/validation/ab_H_imgt.csv");
        if !path.exists() {
            return; // Skip if file doesn't exist
        }

        let entries = load_validation_csv(&path).unwrap();
        assert!(!entries.is_empty());

        // Check first entry
        let first = &entries[0];
        assert!(!first.header.is_empty());
        assert!(!first.sequence.is_empty());
        assert!(!first.expected_positions.is_empty());
    }

    /// Helper function to validate all sequences for a given chain with IMGT scheme
    fn validate_chain_sequences(
        chain: Chain,
        csv_path: &str,
        min_perfect_pct: f64,
        min_overall_accuracy: f64,
    ) {
        validate_chain_sequences_filtered(
            chain,
            csv_path,
            Scheme::IMGT,
            None,
            &[],
            min_perfect_pct,
            min_overall_accuracy,
        );
    }

    /// Helper function to validate all sequences for a given chain with a specific scheme
    fn validate_chain_sequences_with_scheme(
        chain: Chain,
        csv_path: &str,
        scheme: Scheme,
        min_perfect_pct: f64,
        min_overall_accuracy: f64,
    ) {
        validate_chain_sequences_filtered(
            chain,
            csv_path,
            scheme,
            None,
            &[],
            min_perfect_pct,
            min_overall_accuracy,
        );
    }

    /// Helper function to validate sequences with species filter and header exclusions
    fn validate_chain_sequences_filtered(
        chain: Chain,
        csv_path: &str,
        scheme: Scheme,
        species_filter: Option<&str>,
        exclude_headers: &[&str],
        min_perfect_pct: f64,
        min_overall_accuracy: f64,
    ) {
        let path = PathBuf::from(csv_path);
        if !path.exists() {
            return; // Skip if file doesn't exist
        }

        let entries = load_validation_csv(&path).unwrap();
        let annotator = Annotator::new(&[chain], scheme).unwrap();

        let mut total_sequences = 0;
        let mut perfect_sequences = 0;
        let mut total_positions = 0;
        let mut correct_positions = 0;
        let mut failures = Vec::new();

        for entry in entries.iter() {
            // Skip excluded sequences
            if exclude_headers.iter().any(|h| entry.header.contains(h)) {
                continue;
            }

            // Apply species filter if provided
            if let Some(species) = species_filter {
                if !entry.species.eq_ignore_ascii_case(species) {
                    continue;
                }
            }

            let result = validate_entry(entry, &annotator, scheme).unwrap();

            total_sequences += 1;
            total_positions += result.total_positions;
            correct_positions += result.correct_positions;

            if result.accuracy() == 1.0 {
                perfect_sequences += 1;
            } else {
                failures.push(result.clone());
            }

            assert_eq!(result.detected_chain, chain);
        }

        let overall_accuracy = (correct_positions as f64) / (total_positions as f64);
        let perfect_percentage = (perfect_sequences as f64) / (total_sequences as f64);

        if !failures.is_empty() {
            println!("\n{} sequences with mismatches:", failures.len());
            for failure in &failures {
                println!(
                    "\n  {}: {:.1}% accuracy ({}/{} correct)",
                    failure.header,
                    failure.accuracy() * 100.0,
                    failure.correct_positions,
                    failure.total_positions
                );

                // Show first 10 mismatches for this sequence
                let mismatches: Vec<_> = failure.mismatches.iter().take(10).collect();
                for (seq_idx, expected, got) in mismatches {
                    let aa = failure.sequence.chars().nth(*seq_idx).unwrap_or('?');
                    println!(
                        "    Seq pos {} ({}): expected {}, got {}",
                        seq_idx, aa, expected, got
                    );
                }
                if failure.mismatches.len() > 10 {
                    println!(
                        "    ... and {} more mismatches",
                        failure.mismatches.len() - 10
                    );
                }
            }
        }

        println!("\n{}", "=".repeat(80));

        println!("\n{}", "=".repeat(80));
        println!("{} Validation Summary", chain);
        println!("{}", "=".repeat(80));
        println!("Total sequences: {}", total_sequences);
        println!(
            "Perfect sequences: {} ({:.1}%)",
            perfect_sequences,
            perfect_percentage * 100.0
        );
        println!(
            "Overall accuracy: {}/{} ({:.2}%)",
            correct_positions,
            total_positions,
            overall_accuracy * 100.0
        );

        assert!(
            perfect_percentage >= min_perfect_pct,
            "{} Sequence accuracy {:.2}% is below {:.0}% threshold",
            chain,
            perfect_percentage * 100.0,
            min_perfect_pct * 100.0
        );

        assert!(
            overall_accuracy >= min_overall_accuracy,
            "{} overall accuracy {:.2}% is below {:.0}% threshold",
            chain,
            overall_accuracy * 100.0,
            min_overall_accuracy * 100.0
        );
    }

    #[test]
    fn test_validate_igh_sequences() {
        validate_chain_sequences(Chain::IGH, "fixtures/validation/ab_H_imgt.csv", 0.99, 0.99);
    }

    #[test]
    fn test_validate_igk_sequences() {
        validate_chain_sequences(Chain::IGK, "fixtures/validation/ab_K_imgt.csv", 0.99, 0.99);
    }

    #[test]
    fn test_validate_igl_sequences() {
        validate_chain_sequences(Chain::IGL, "fixtures/validation/ab_L_imgt.csv", 0.99, 0.99);
    }

    #[test]
    fn test_validate_tra_sequences() {
        validate_chain_sequences(Chain::TRA, "fixtures/validation/tcr_A_imgt.csv", 0.99, 0.99);
    }

    #[test]
    fn test_validate_trb_sequences() {
        validate_chain_sequences(Chain::TRB, "fixtures/validation/tcr_B_imgt.csv", 0.99, 0.99);
    }

    #[test]
    fn test_validate_trg_sequences() {
        validate_chain_sequences(Chain::TRG, "fixtures/validation/tcr_G_imgt.csv", 0.99, 0.99);
    }

    #[test]
    fn test_validate_trd_sequences() {
        validate_chain_sequences(Chain::TRD, "fixtures/validation/tcr_D_imgt.csv", 0.99, 0.99);
    }

    // Kabat validation tests
    #[test]
    fn test_validate_igh_kabat_sequences() {
        validate_chain_sequences_with_scheme(
            Chain::IGH,
            "fixtures/validation/ab_H_kabat.csv",
            Scheme::Kabat,
            0.99,
            0.99,
        );
    }

    #[test]
    fn test_validate_igk_kabat_sequences() {
        validate_chain_sequences_with_scheme(
            Chain::IGK,
            "fixtures/validation/ab_K_kabat.csv",
            Scheme::Kabat,
            0.99,
            0.99,
        );
    }

    #[test]
    fn test_validate_igl_kabat_sequences() {
        validate_chain_sequences_with_scheme(
            Chain::IGL,
            "fixtures/validation/ab_L_kabat.csv",
            Scheme::Kabat,
            0.99,
            0.99,
        );
    }
}
