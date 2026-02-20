//! Sequence alignment using Needleman-Wunsch algorithm

use crate::error::{Error, Result};
use crate::scoring::PositionScores;

/// Direction in the alignment traceback matrix
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Direction {
    Match,
    GapInQuery,
    GapInConsensus,
}

/// Represents a position in the alignment result
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignedPosition {
    /// Gap in query (query has no residue here, consensus does)
    QueryGap,
    /// Aligned to consensus position N
    Aligned(u8),
    /// Insertion in query (query has residue, consensus doesn't)
    Insertion,
}

/// Cell in the dynamic programming matrix
#[derive(Debug, Clone)]
struct Cell {
    score: f32,
    direction: Direction,
}

/// Result of aligning a sequence to a consensus
#[derive(Debug, Clone)]
pub struct Alignment {
    /// Alignment score
    pub score: f32,
    /// Aligned positions (merged representation of query and consensus alignment)
    pub positions: Vec<AlignedPosition>,
    /// Start position in consensus
    pub start_pos: u8,
    /// End position in consensus
    pub end_pos: u8,
}

/// Align a query sequence to a scoring matrix using Needleman-Wunsch
pub fn align(query: &str, positions: &[PositionScores]) -> Result<Alignment> {
    let query = query.to_uppercase();
    let query_bytes: Vec<u8> = query.bytes().collect();
    let query_len = query_bytes.len();
    let cons_len = positions.len();

    if query_len == 0 {
        return Err(Error::InvalidSequence("empty sequence".to_string()));
    }

    // Initialize DP matrix
    let mut dp = vec![
        vec![
            Cell {
                score: 0.0,
                direction: Direction::Match,
            };
            cons_len + 1
        ];
        query_len + 1
    ];

    // Initialize first row and column for semi-global alignment
    // Allow query to start anywhere without penalty (gaps at start of query are free)
    for j in 1..=cons_len {
        dp[0][j] = Cell {
            score: 0.0, // No penalty for skipping consensus positions at start
            direction: Direction::GapInQuery,
        };
    }

    // Allow consensus to start anywhere without penalty (gaps at start of consensus are free)
    for row in dp.iter_mut().take(query_len + 1).skip(1) {
        row[0] = Cell {
            score: 0.0, // No penalty for query starting before consensus
            direction: Direction::GapInConsensus,
        };
    }

    // Fill DP matrix
    for i in 1..=query_len {
        let query_aa = query_bytes[i - 1] as char;

        for j in 1..=cons_len {
            let cons_pos = &positions[j - 1];

            // Match/mismatch score
            let match_score = cons_pos.scores.get(&query_aa).copied().unwrap_or(-4.0);
            let from_match = dp[i - 1][j - 1].score + match_score;

            // Gap in query (skip consensus position)
            let from_gap_query = dp[i][j - 1].score + cons_pos.gap_penalty;

            // Gap in consensus (insertion in query sequence)
            let from_gap_cons = dp[i - 1][j].score + cons_pos.insertion_penalty;

            // Choose best
            let (best_score, best_dir) =
                if from_match >= from_gap_query && from_match >= from_gap_cons {
                    (from_match, Direction::Match)
                } else if from_gap_query >= from_gap_cons {
                    (from_gap_query, Direction::GapInQuery)
                } else {
                    (from_gap_cons, Direction::GapInConsensus)
                };

            dp[i][j] = Cell {
                score: best_score,
                direction: best_dir,
            };
        }
    }

    // Find best ending position for semi-global alignment
    // Force full query consumption to prevent long CDR3s from being treated as trailing sequence
    let mut best_i = query_len;
    let mut best_j = cons_len;
    let mut best_score = dp[query_len][cons_len].score;

    // Check last row: query fully used, allow trailing gaps in consensus
    for j in 0..cons_len {
        if dp[query_len][j].score > best_score {
            best_score = dp[query_len][j].score;
            best_i = query_len;
            best_j = j;
        }
    }

    // NOTE: We do NOT allow early termination of query (partial query alignment)
    // This prevents long CDR3 regions from being incorrectly treated as unaligned trailing sequence.
    // The query must align through its full length, but the consensus can have trailing gaps.

    let aligned_positions = traceback(&dp, &query_bytes, positions, best_i, best_j);

    // Find start and end positions from aligned consensus positions
    let start_pos = aligned_positions
        .iter()
        .find_map(|p| match p {
            AlignedPosition::Aligned(n) => Some(*n),
            _ => None,
        })
        .unwrap_or(1);
    let end_pos = aligned_positions
        .iter()
        .rev()
        .find_map(|p| match p {
            AlignedPosition::Aligned(n) => Some(*n),
            _ => None,
        })
        .unwrap_or(128);

    Ok(Alignment {
        score: best_score,
        positions: aligned_positions,
        start_pos,
        end_pos,
    })
}

fn traceback(
    dp: &[Vec<Cell>],
    _query: &[u8],
    positions: &[PositionScores],
    query_len: usize,
    cons_len: usize,
) -> Vec<AlignedPosition> {
    let mut aligned_positions = Vec::new();

    let mut i = query_len;
    let mut j = cons_len;

    // For semi-global alignment, handle unaligned suffix of query
    let query_total_len = _query.len();
    if i < query_total_len {
        // Add remaining query residues as insertions at the end
        for _ in i..query_total_len {
            aligned_positions.push(AlignedPosition::Insertion);
        }
    }

    while i > 0 && j > 0 {
        match dp[i][j].direction {
            Direction::Match => {
                aligned_positions.push(AlignedPosition::Aligned(positions[j - 1].position));
                i -= 1;
                j -= 1;
            }
            Direction::GapInQuery => {
                aligned_positions.push(AlignedPosition::QueryGap);
                j -= 1;
            }
            Direction::GapInConsensus => {
                aligned_positions.push(AlignedPosition::Insertion);
                i -= 1;
            }
        }
    }

    // Handle unaligned prefix of query
    while i > 0 {
        aligned_positions.push(AlignedPosition::Insertion);
        i -= 1;
    }

    // Reverse (we built it backwards)
    aligned_positions.reverse();

    aligned_positions
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::imgt::get_imgt_numbering;
    use crate::scoring::ScoringMatrix;
    use crate::Chain;

    #[test]
    fn test_align_simple() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

        // A simple sequence that should align
        let sequence =
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN";

        let result = align(sequence, &matrix.positions).unwrap();

        // Alignment should produce a score (may be negative for partial sequences)
        assert!(!result.positions.is_empty());
        assert!(result.start_pos > 0);
        assert!(result.end_pos > 0);
    }

    #[test]
    fn test_numbering() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let sequence =
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN";

        let result = align(sequence, &matrix.positions).unwrap();
        let numbering = get_imgt_numbering(&result);

        assert_eq!(numbering.len(), sequence.len());
        assert!(numbering[0].number > 0);
    }

    #[test]
    fn test_empty_sequence() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let result = align("", &matrix.positions);
        assert!(result.is_err());
    }

    #[test]
    fn test_numbering_length_matches_sequence() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

        // Test with different sequences
        let sequences = vec![
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN",
            "QVQLQQSGAELMKPGASVKISCKATGYTFSSYWIEWVKQRPGHGLEWIGEILPGSGSTNY",
            "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR",
        ];

        for seq in sequences {
            let result = align(seq, &matrix.positions).unwrap();
            let numbering = get_imgt_numbering(&result);

            assert_eq!(
                numbering.len(),
                seq.len(),
                "Numbering length {} doesn't match sequence length {} for sequence starting with {}",
                numbering.len(),
                seq.len(),
                &seq[..20.min(seq.len())]
            );
        }
    }

    #[test]
    fn test_partial_sequence_start() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

        // A sequence starting from middle of framework should align without heavy penalty
        let partial_seq = "GLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN";
        let result = align(partial_seq, &matrix.positions).unwrap();

        // Should have a reasonable score (not heavily penalized for missing start)
        assert!(
            result.score > -100.0,
            "Partial sequence should not be heavily penalized"
        );
        // Allow some gaps but not excessive (e.g., <10% of sequence length)
        let gap_count = result
            .positions
            .iter()
            .filter(|p| matches!(p, AlignedPosition::QueryGap))
            .count();
        assert!(
            gap_count < partial_seq.len() / 10,
            "Query should have minimal gaps, found {} gaps",
            gap_count
        );
    }

    #[test]
    fn test_partial_sequence_end() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

        // A sequence ending in middle should align without heavy penalty
        let partial_seq = "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVS";
        let result = align(partial_seq, &matrix.positions).unwrap();

        // Should have a reasonable score
        assert!(
            result.score > -100.0,
            "Partial sequence should not be heavily penalized"
        );
    }

    #[test]
    fn test_fragment_alignment() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

        // A small fragment from the middle should align
        let fragment = "GLEWVSAISKSGGSTYY";
        let result = align(fragment, &matrix.positions).unwrap();

        // Should align without extreme penalties for missing ends
        assert!(
            result.score > -50.0,
            "Fragment should align with reasonable score"
        );
        assert!(
            result.start_pos > 1,
            "Fragment should start after position 1"
        );
        assert!(
            result.end_pos < 128,
            "Fragment should end before last position"
        );
    }
}
