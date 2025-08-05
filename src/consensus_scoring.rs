use crate::constants::{ScoringParams, ENCODED_RESIDUES_MAP};
use crate::constants::{ACCEPTED_RESIDUES, BLOSUM62};
use crate::scoring_matrix::ScoringMatrix;
use std::collections::HashMap;

/// Calculate scoring matrix from consensus data and scoring parameters
pub fn calculate_scoring_matrix(
    consensus: &HashMap<u32, Vec<u8>>,
    scoring_params: &ScoringParams,
    scheme_gap_penalty_fn: impl Fn(u32, &ScoringParams) -> (f64, f64),
) -> ScoringMatrix {
    let number_accepted_residues = ACCEPTED_RESIDUES.len();
    let consensus_length = consensus.len();

    let mut matrix = ScoringMatrix::zeros(consensus_length, number_accepted_residues + 2);

    // Fill in scores
    for consensus_position in 0..consensus_length {
        for residue_index in 0..number_accepted_residues {
            let residue: u8 = ACCEPTED_RESIDUES[residue_index];
            let score =
                best_score_consensus((consensus_position + 1) as u32, residue, consensus) as f64;
            matrix.set(consensus_position, residue_index, score);
        }

        // Fill in gap penalties
        let (query_gap, consensus_gap) =
            scheme_gap_penalty_fn((consensus_position + 1) as u32, scoring_params);

        // Query gap
        matrix.set(consensus_position, number_accepted_residues, query_gap);
        // Consensus gap
        matrix.set(
            consensus_position,
            number_accepted_residues + 1,
            consensus_gap,
        );
    }

    matrix
}

/// Finds highest score according to blosum62
fn best_score_consensus(position: u32, residue: u8, consensus: &HashMap<u32, Vec<u8>>) -> i32 {
    // Initialize to minimal score
    let to_check_residues: &Vec<u8> = consensus
        .get(&position)
        .expect("Position outside of consensus");

    // when any is allowed, best score is perfect match
    if to_check_residues.iter().all(|&c| c == b'-') {
        return blosum_lookup(&residue, &residue);
    }
    // iterate over residues to get max score
    blosum_lookup(
        &residue,
        to_check_residues
            .iter()
            .max_by_key(|&c| blosum_lookup(&residue, c))
            .unwrap(),
    )
}

/// Utility to easily find blosum62 score for two residues
fn blosum_lookup(residue1: &u8, residue2: &u8) -> i32 {
    if residue1 < residue2 {
        let lookup: &[u8; 2] = &[*residue1, *residue2];
        *BLOSUM62.get(lookup).unwrap()
    } else {
        let lookup: &[u8; 2] = &[*residue2, *residue1];
        *BLOSUM62.get(lookup).unwrap()
    }
}

/// Encodes sequence to indexes for the scoring matrix
pub(crate) fn encode_sequence(input: &[u8]) -> Vec<u8> {
    // Create a lookup table (once, could be static)
    input
        .iter()
        .map(|&residue| *ENCODED_RESIDUES_MAP.get(&residue).unwrap_or(&128))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::get_scoring_params;

    #[test]
    fn test_lookup() {
        // working
        assert_eq!(blosum_lookup(&b'B', &b'N'), 3);
        assert_eq!(blosum_lookup(&b'W', &b'Y'), 2);
        assert_eq!(blosum_lookup(&b'X', &b'E'), -1);
        // same either way of inputting arguments
        assert_eq!(blosum_lookup(&b'C', &b'S'), blosum_lookup(&b'S', &b'C'));
    }

    #[test]
    fn amino_acid_encoding() {
        assert_eq!(
            encode_sequence("ABCDEFGHIKLMNPQRSTVWXYZ".as_bytes()),
            (0..23).collect::<Vec<u8>>()
        )
    }
}
