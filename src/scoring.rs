//! Scoring matrices for sequence alignment

use crate::error::{Error, Result};
use crate::types::Chain;
use serde::{Deserialize, Serialize};

/// Default score for unknown amino acids
const DEFAULT_SCORE: f32 = -4.0;

/// A position-specific scoring matrix for a chain type
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScoringMatrix {
    /// Scoring data for each position
    pub positions: Vec<PositionScores>,
}

/// Scores and gap penalties for a single position
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PositionScores {
    /// IMGT position number
    pub position: u8,
    /// Scores indexed by (amino_acid_byte - b'A'), 26 slots for A-Z
    /// Non-amino-acid letters default to DEFAULT_SCORE
    pub scores: [f32; 26],
    /// Penalty for gap in query sequence (skipping this consensus position)
    pub gap_penalty: f32,
    /// Penalty for gap in consensus (insertion in query at this position)
    pub insertion_penalty: f32,
}

impl PositionScores {
    /// Get score for an amino acid byte (must be uppercase A-Z)
    #[inline(always)]
    pub fn score_for(&self, aa: u8) -> f32 {
        let idx = aa.wrapping_sub(b'A') as usize;
        if idx < 26 {
            // SAFETY: idx is guaranteed < 26, which is the array length
            unsafe { *self.scores.get_unchecked(idx) }
        } else {
            DEFAULT_SCORE
        }
    }
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_igh_matrix() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        assert!(!matrix.positions.is_empty());
        assert!(matrix.positions.len() > 100);
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
            assert!(!matrix.positions.is_empty());
        }
    }

    #[test]
    fn test_gap_penalties() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

        // All positions should have gap penalties
        for pos in &matrix.positions {
            assert!(
                pos.gap_penalty <= 0.0,
                "Gap in query penalty should be non-positive,"
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
            for (i, &score) in pos.scores.iter().enumerate() {
                let aa = (b'A' + i as u8) as char;
                assert!(
                    score >= -10.0 && score <= 15.0,
                    "Score for {} at position {} is out of range: {}",
                    aa,
                    pos.position,
                    score
                );
            }
        }
    }
}
