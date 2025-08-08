use crate::consensus_scoring::encode_sequence;
use crate::constants::{
    traceback_directions, CONSENSUS_GAP_COLUMN, GAP_PEN_END, GAP_PEN_START, MATCH_CP_MULTIPLIER,
    QUERY_GAP_COLUMN,
};
use crate::numbering_scheme_type::NumberingScheme;

/// Optimized matrix pool for reusing memory between alignments
pub struct MatrixPool {
    dynamic_matrix: Vec<f64>,
    traceback_matrix: Vec<u8>,
    current_capacity: (usize, usize), // (consensus_len, query_len)
}

impl MatrixPool {
    pub fn new() -> Self {
        Self {
            dynamic_matrix: Vec::new(),
            traceback_matrix: Vec::new(),
            current_capacity: (0, 0),
        }
    }

    /// Get matrices with at least the required size, reusing memory when possible
    pub fn get_matrices(
        &mut self,
        consensus_len: usize,
        query_len: usize,
    ) -> (&mut [f64], &mut [u8]) {
        let required_size = (consensus_len + 1) * (query_len + 1);

        // Resize if needed
        if self.dynamic_matrix.len() < required_size {
            self.dynamic_matrix.resize(required_size, 0.0);
            self.traceback_matrix.resize(required_size, 0);
            self.current_capacity = (consensus_len, query_len);
        }

        // Clear only the portion we'll use
        let used_size = (consensus_len + 1) * (query_len + 1);
        self.dynamic_matrix[..used_size].fill(0.0);
        self.traceback_matrix[..used_size].fill(0);

        (
            &mut self.dynamic_matrix[..used_size],
            &mut self.traceback_matrix[..used_size],
        )
    }
}

impl Default for MatrixPool {
    fn default() -> Self {
        Self::new()
    }
}

/// Optimized Needleman-Wunsch alignment with flat matrix and memory reuse
pub fn needleman_wunsch_consensus(
    query_sequence: &[u8],
    scheme: &NumberingScheme,
    matrix_pool: &mut MatrixPool,
    early_termination_threshold: f64,
) -> (Vec<String>, f64) {
    let encoded_query = encode_sequence(query_sequence);
    // The Vec is sized max_position + 1, but positions are 1-indexed, so effective positions are len - 1
    let num_positions_consensus = scheme.consensus_amino_acids.len() - 1;
    let len_query_sequence = query_sequence.len();

    // Get reusable matrices with flat layout for better cache performance
    let (dynamic_matrix, traceback_matrix) =
        matrix_pool.get_matrices(num_positions_consensus, len_query_sequence);

    let cols = len_query_sequence + 1;

    // Initialize first row and column - using flat indexing
    for i in 0..=len_query_sequence {
        let idx = i; // First row: row=0, col=i
        dynamic_matrix[idx] = -(GAP_PEN_START * i as f64);
        traceback_matrix[idx] = traceback_directions::FROM_LEFT;
    }

    for i in 0..=num_positions_consensus {
        let idx = i * cols; // First column: row=i, col=0
        dynamic_matrix[idx] = -(GAP_PEN_START * i as f64);
        traceback_matrix[idx] = traceback_directions::FROM_TOP;
    }

    // Track best score seen so far for early termination
    let mut current_best_score = f64::NEG_INFINITY;
    let mut rows_without_improvement = 0;
    const MAX_ROWS_WITHOUT_IMPROVEMENT: usize = 20;

    // Main dynamic programming loop with flat matrix indexing
    for consensus_position in 1..=num_positions_consensus {
        let mut row_max_score = f64::NEG_INFINITY;

        for query_position in 1..=len_query_sequence {
            let current_idx = consensus_position * cols + query_position;
            let top_idx = (consensus_position - 1) * cols + query_position;
            let left_idx = consensus_position * cols + (query_position - 1);
            let diag_idx = (consensus_position - 1) * cols + (query_position - 1);

            // Get gap penalties - avoid repeated matrix lookups
            let (query_gap_penalty, consensus_gap_penalty) = if consensus_position
                == num_positions_consensus
                || query_position == len_query_sequence
            {
                (GAP_PEN_END, GAP_PEN_END)
            } else {
                (
                    scheme.scoring_matrix[[consensus_position - 1, QUERY_GAP_COLUMN]],
                    scheme.scoring_matrix[[consensus_position - 1, CONSENSUS_GAP_COLUMN]],
                )
            };

            // Calculate scores for three possible moves
            let top_value = dynamic_matrix[top_idx] - consensus_gap_penalty;
            let left_value = dynamic_matrix[left_idx] - query_gap_penalty;
            let mut match_value = dynamic_matrix[diag_idx];

            // Get match score
            let mut best_score = scheme.scoring_matrix[[
                consensus_position - 1,
                encoded_query[query_position - 1] as usize,
            ]];

            // Apply conserved position multiplier - using pre-computed set
            if scheme
                .conserved_positions
                .contains(&(consensus_position as u32))
            {
                best_score *= MATCH_CP_MULTIPLIER;
            }

            match_value += best_score;

            // Find maximum and set traceback direction
            let mut max_score = match_value;
            let mut transfer = traceback_directions::FROM_DIAG;

            if top_value > max_score {
                max_score = top_value;
                transfer = traceback_directions::FROM_TOP;
            }
            if left_value > max_score {
                max_score = left_value;
                transfer = traceback_directions::FROM_LEFT;
            }

            dynamic_matrix[current_idx] = max_score;
            traceback_matrix[current_idx] = transfer;

            // Check for perfect match - using pre-computed data
            let seq_char = query_sequence[query_position - 1];
            if transfer == traceback_directions::FROM_DIAG
                && (consensus_position) < scheme.consensus_amino_acids.len()
                && scheme.consensus_amino_acids[consensus_position].contains(&seq_char)
                && scheme
                    .restricted_sites
                    .contains(&(consensus_position as u32))
            {
                traceback_matrix[current_idx] = traceback_directions::PERFECT_MATCH;
            }

            // Track row maximum for early termination
            if max_score > row_max_score {
                row_max_score = max_score;
            }
        }

        // Early termination check - if we haven't improved in many rows and score is poor
        if row_max_score > current_best_score {
            current_best_score = row_max_score;
            rows_without_improvement = 0;
        } else {
            rows_without_improvement += 1;

            // Early exit if score is very poor and no improvement for many rows
            if rows_without_improvement >= MAX_ROWS_WITHOUT_IMPROVEMENT
                && current_best_score < early_termination_threshold
            {
                // Fill remaining matrix with default values for traceback
                for remaining_consensus in consensus_position + 1..=num_positions_consensus {
                    for remaining_query in 1..=len_query_sequence {
                        let idx = remaining_consensus * cols + remaining_query;
                        dynamic_matrix[idx] = f64::NEG_INFINITY;
                        traceback_matrix[idx] = traceback_directions::FROM_TOP;
                    }
                }
                break;
            }
        }
    }

    // Traceback using flat matrix
    let (numbering, matches) = traceback_alignment_flat(
        len_query_sequence,
        num_positions_consensus,
        traceback_matrix,
        cols,
    );

    let identity = matches as f64 / (scheme.restricted_sites.len() as f64).max(1.0);
    (numbering, identity)
}

/// Optimized traceback with flat matrix indexing
fn traceback_alignment_flat(
    len_query_sequence: usize,
    num_positions_consensus: usize,
    traceback_matrix: &[u8],
    cols: usize,
) -> (Vec<String>, u32) {
    let mut numbering = Vec::new();
    let mut matches = 0u32;

    let mut consensus_pos = num_positions_consensus;
    let mut query_pos = len_query_sequence;

    while consensus_pos > 0 || query_pos > 0 {
        let current_idx = consensus_pos * cols + query_pos;
        let direction = traceback_matrix[current_idx];

        if direction == traceback_directions::FROM_LEFT {
            // Insertion in query (gap in consensus) - add dash for missing consensus position
            numbering.insert(0, "-".to_string());
            query_pos = query_pos.saturating_sub(1);
        } else if direction == traceback_directions::FROM_TOP {
            // Deletion in query (gap in query) - skip, no position added to numbering
            consensus_pos = consensus_pos.saturating_sub(1);
        } else {
            // FROM_DIAG or PERFECT_MATCH - add position number
            numbering.insert(0, consensus_pos.to_string());
            if direction == traceback_directions::PERFECT_MATCH {
                matches += 1;
            }
            consensus_pos = consensus_pos.saturating_sub(1);
            query_pos = query_pos.saturating_sub(1);
        }
    }

    (numbering, matches)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        schemes::get_scheme,
        types::{Chain, Scheme},
    };

    #[test]
    fn test_needleman_wunsch_consensus() {
        let scheme = get_scheme(Scheme::IMGT, Chain::IGH, None);
        let mut pool = MatrixPool::new();

        let test_sequence = b"QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS";

        // Test the function
        let (numbering, identity) = needleman_wunsch_consensus(
            test_sequence,
            &scheme,
            &mut pool,
            -100.0, // Disable early termination
        );

        // Basic sanity checks
        assert!(identity > 0.0, "Identity should be positive");
        assert!(!numbering.is_empty(), "Numbering should not be empty");
        assert!(
            numbering.len() <= test_sequence.len(),
            "Numbering should not be longer than sequence"
        );
    }

    #[test]
    fn test_matrix_pool_reuse() {
        let mut pool = MatrixPool::new();

        // First allocation
        let (mat1, _) = pool.get_matrices(100, 200);
        let initial_capacity = mat1.len();

        // Second allocation with smaller size - should reuse
        let (mat2, _) = pool.get_matrices(50, 100);
        let expected_size = (50 + 1) * (100 + 1);
        assert_eq!(mat2.len(), expected_size); // Should return the slice size for smaller matrices

        // Third allocation with larger size - should expand
        let (mat3, _) = pool.get_matrices(150, 300);
        assert!(mat3.len() >= initial_capacity);
    }
}
