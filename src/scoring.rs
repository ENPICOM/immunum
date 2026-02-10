//! Scoring matrices for sequence alignment

use crate::error::{Error, Result};
use crate::types::Chain;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A position-specific scoring matrix for a chain type
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScoringMatrix {
    /// The chain type this matrix is for
    pub chain: String,
    /// Scoring data for each position
    pub positions: Vec<PositionScores>,
}

/// Scores and gap penalties for a single position
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PositionScores {
    /// IMGT position number
    pub position: u32,
    /// Scores for each amino acid at this position
    pub scores: HashMap<char, f32>,
    /// Penalty for gap in query sequence (skipping this consensus position)
    pub gap_penalty: f32,
    /// Penalty for gap in consensus (insertion in query at this position)
    pub insertion_penalty: f32,
}

impl ScoringMatrix {
    /// Load a scoring matrix for a specific chain type
    pub fn load(chain: Chain) -> Result<Self> {
        let json = match chain {
            Chain::IGH => include_str!(concat!(env!("OUT_DIR"), "/matrices/IGH.json")),
            Chain::IGK => include_str!(concat!(env!("OUT_DIR"), "/matrices/IGK.json")),
            Chain::IGL => include_str!(concat!(env!("OUT_DIR"), "/matrices/IGL.json")),
            Chain::TRA => include_str!(concat!(env!("OUT_DIR"), "/matrices/TRA.json")),
            Chain::TRB => include_str!(concat!(env!("OUT_DIR"), "/matrices/TRB.json")),
            Chain::TRG => include_str!(concat!(env!("OUT_DIR"), "/matrices/TRG.json")),
            Chain::TRD => include_str!(concat!(env!("OUT_DIR"), "/matrices/TRD.json")),
        };

        serde_json::from_str(json).map_err(|e| {
            Error::ConsensusParseError(format!("Failed to parse scoring matrix: {}", e))
        })
    }

    /// Get the score for an amino acid at a specific position
    pub fn get_score(&self, position: u32, amino_acid: char) -> Option<f32> {
        self.positions
            .iter()
            .find(|p| p.position == position)
            .and_then(|p| p.scores.get(&amino_acid).copied())
    }

    /// Get gap penalties for a specific position
    pub fn get_gap_penalties(&self, position: u32) -> Option<(f32, f32)> {
        self.positions
            .iter()
            .find(|p| p.position == position)
            .map(|p| (p.gap_penalty, p.insertion_penalty))
    }

    /// Get position scores by position number
    pub fn get_position(&self, position: u32) -> Option<&PositionScores> {
        self.positions.iter().find(|p| p.position == position)
    }

    /// Get the number of positions in this matrix
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Check if the matrix is empty
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_igh_matrix() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        assert_eq!(matrix.chain, "IGH");
        assert!(!matrix.is_empty());
        assert!(matrix.len() > 100);
    }

    #[test]
    fn test_load_all_matrices() {
        // Test that all chain matrices can be loaded
        for chain in [
            Chain::IGH,
            Chain::IGK,
            Chain::IGL,
            Chain::TRA,
            Chain::TRB,
            Chain::TRG,
            Chain::TRD,
        ] {
            let matrix = ScoringMatrix::load(chain).unwrap();
            assert!(!matrix.is_empty());
        }
    }

    #[test]
    fn test_get_score() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

        // Position 1 should have scores
        let score = matrix.get_score(1, 'E');
        assert!(score.is_some());

        // Should return None for non-existent position
        let score = matrix.get_score(9999, 'A');
        assert!(score.is_none());
    }

    #[test]
    fn test_gap_penalties() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

        // All positions should have gap penalties
        for pos in &matrix.positions {
            assert!(
                pos.gap_penalty < 0.0,
                "Gap in query penalty should be negative"
            );
            assert!(
                pos.insertion_penalty < 0.0,
                "Gap in consensus penalty should be negative"
            );
        }
    }

    #[test]
    fn test_scores_reasonable() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

        // Check that scores are in reasonable range (BLOSUM62 range is roughly -4 to 11)
        for pos in &matrix.positions {
            for (aa, score) in &pos.scores {
                assert!(
                    *score >= -10.0 && *score <= 15.0,
                    "Score for {} at position {} is out of range: {}",
                    aa,
                    pos.position,
                    score
                );
            }
        }
    }
}
