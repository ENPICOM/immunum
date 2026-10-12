//! Scoring matrices for sequence alignment

#[cfg(feature = "python")]
use pyo3::prelude::*;

use crate::types::Chain;
use serde::{Deserialize, Serialize};

/// Default score for unknown amino acids
const DEFAULT_SCORE: f32 = -4.0;

/// A position-specific scoring matrix for a chain type
#[cfg_attr(feature = "python", pyclass(get_all))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScoringMatrix {
    /// Scoring data for each position
    pub positions: Vec<PositionScores>,
}

/// Scores and gap penalties for a single position
#[cfg_attr(feature = "python", pyclass(get_all))]
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
    /// Best possible score at this position (max of all amino acid scores)
    pub max_score: f32,
    /// Whether this position counts toward the confidence score (occupancy > 0.5)
    pub counts_for_confidence: bool,
}

impl PositionScores {
    /// Get score for an amino acid byte (must be uppercase A-Z)
    #[inline(always)]
    pub fn score_for(&self, aa: u8) -> f32 {
        let idx = aa.wrapping_sub(b'A') as usize;
        if idx < 26 {
            self.scores[idx]
        } else {
            DEFAULT_SCORE
        }
    }
}

impl ScoringMatrix {
    /// The built-in scoring matrix for a chain type
    pub fn load(chain: Chain) -> Self {
        let json = match chain {
            Chain::IGH => include_str!(concat!(env!("OUT_DIR"), "/matrices/IGH.json")),
            Chain::IGK => include_str!(concat!(env!("OUT_DIR"), "/matrices/IGK.json")),
            Chain::IGL => include_str!(concat!(env!("OUT_DIR"), "/matrices/IGL.json")),
            Chain::TRA => include_str!(concat!(env!("OUT_DIR"), "/matrices/TRA.json")),
            Chain::TRB => include_str!(concat!(env!("OUT_DIR"), "/matrices/TRB.json")),
            Chain::TRG => include_str!(concat!(env!("OUT_DIR"), "/matrices/TRG.json")),
            Chain::TRD => include_str!(concat!(env!("OUT_DIR"), "/matrices/TRD.json")),
        };
        // The matrices are generated at build time; one that doesn't parse is a build bug, which
        // `test_load_all_matrices` catches
        serde_json::from_str(json).expect("built-in scoring matrix parses")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_igh_matrix() {
        let matrix = ScoringMatrix::load(Chain::IGH);
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
            let matrix = ScoringMatrix::load(chain);
            assert!(!matrix.positions.is_empty());
        }
    }

    #[test]
    fn test_gap_penalties() {
        let matrix = ScoringMatrix::load(Chain::IGH);

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
        let matrix = ScoringMatrix::load(Chain::IGH);

        // Check that scores are in reasonable range (BLOSUM62 range is roughly -4 to 11)
        for pos in &matrix.positions {
            for (i, &score) in pos.scores.iter().enumerate() {
                let aa = (b'A' + i as u8) as char;
                assert!(
                    (-10.0..=15.0).contains(&score),
                    "Score for {} at position {} is out of range: {}",
                    aa,
                    pos.position,
                    score
                );
            }
        }
    }
}
