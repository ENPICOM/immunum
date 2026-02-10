//! Sequence alignment using Needleman-Wunsch algorithm

use crate::error::{Error, Result};
use crate::imgt;
use crate::scoring::ScoringMatrix;
use crate::types::{Position, Region};

/// Direction in the alignment traceback matrix
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Direction {
    Match,
    GapInQuery,
    GapInConsensus,
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
    /// Query sequence (input)
    pub query: String,
    /// Aligned query sequence (with gaps)
    pub aligned_query: String,
    /// Aligned consensus positions
    pub aligned_positions: Vec<Option<u32>>,
    /// Start position in consensus
    pub start_pos: u32,
    /// End position in consensus
    pub end_pos: u32,
}

impl Alignment {
    /// Get IMGT-specific numbering: FR regions from alignment, CDR regions from IMGT rules
    pub fn get_imgt_numbering(&self) -> Vec<Position> {
        // First, get alignment-based numbering to identify regions
        let alignment_positions = self.get_alignment_positions();

        // Group consecutive positions by region
        let mut regions: Vec<(Region, Vec<usize>)> = Vec::new();

        for (seq_idx, pos) in alignment_positions.iter().enumerate() {
            let region = imgt::position_to_region(pos.number);

            match regions.last_mut() {
                Some((last_region, indices)) if *last_region == region => {
                    indices.push(seq_idx);
                }
                _ => {
                    regions.push((region, vec![seq_idx]));
                }
            }
        }

        // Build final numbering: use IMGT rules for CDR, alignment for FR
        let mut final_positions = vec![Position::new(1); alignment_positions.len()];

        for (region, indices) in regions {
            let length = indices.len();
            let region_positions = match region {
                Region::CDR1 => imgt::cdr1_numbering(length),
                Region::CDR2 => imgt::cdr2_numbering(length),
                Region::CDR3 => imgt::cdr3_numbering(length),
                _ => {
                    // For FR regions, use alignment-based numbering
                    indices
                        .iter()
                        .map(|&i| alignment_positions[i].clone())
                        .collect()
                }
            };

            for (i, &seq_idx) in indices.iter().enumerate() {
                final_positions[seq_idx] = region_positions[i].clone();
            }
        }

        final_positions
    }

    /// Get Kabat-specific numbering: FR regions mapped from IMGT, CDR regions from Kabat rules
    pub fn get_kabat_numbering(&self, chain: crate::types::Chain) -> Vec<Position> {
        use crate::kabat;
        use crate::types::Chain;

        // Get IMGT alignment-based numbering
        let alignment_positions = self.get_alignment_positions();

        // Convert to Kabat numbering using chain-specific mapping
        match chain {
            Chain::IGH => kabat::imgt_to_kabat_heavy(&alignment_positions),
            Chain::IGK | Chain::IGL => kabat::imgt_to_kabat_light(&alignment_positions),
            _ => panic!("Kabat numbering only supported for antibody chains"),
        }
    }

    /// Get the assigned positions for each residue from alignment only
    fn get_alignment_positions(&self) -> Vec<Position> {
        let mut positions = Vec::new();
        let mut last_consensus_pos: Option<u32> = None;
        let mut insertion_count = 0u32;

        for (i, &cons_pos) in self.aligned_positions.iter().enumerate() {
            // Skip gaps in the query
            if self.aligned_query.chars().nth(i) == Some('-') {
                continue;
            }

            match cons_pos {
                Some(pos) if Some(pos) != last_consensus_pos => {
                    // New consensus position
                    positions.push(Position::new(pos));
                    last_consensus_pos = Some(pos);
                    insertion_count = 0;
                }
                _ => {
                    // Insertion: either gap in consensus or same consensus position repeated
                    let base_pos = last_consensus_pos.unwrap_or(1);
                    insertion_count += 1;
                    positions.push(Position::with_insertion(
                        base_pos,
                        insertion_to_letter(insertion_count),
                    ));
                }
            }
        }

        positions
    }
}

/// Convert insertion count to IMGT letter (1->A, 2->B, ..., 26->Z)
fn insertion_to_letter(count: u32) -> char {
    if count == 0 || count > 26 {
        'A' // Fallback
    } else {
        (b'A' + (count - 1) as u8) as char
    }
}

/// Align a query sequence to a scoring matrix using Needleman-Wunsch
pub fn align(query: &str, matrix: &ScoringMatrix) -> Result<Alignment> {
    let query = query.to_uppercase();
    let query_bytes: Vec<u8> = query.bytes().collect();
    let query_len = query_bytes.len();
    let cons_len = matrix.positions.len();

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
            let cons_pos = &matrix.positions[j - 1];

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

    let (aligned_query, aligned_positions) = traceback(&dp, &query_bytes, matrix, best_i, best_j);

    // Find start and end positions
    let start_pos = aligned_positions.iter().find_map(|&p| p).unwrap_or(1);
    let end_pos = aligned_positions
        .iter()
        .rev()
        .find_map(|&p| p)
        .unwrap_or(128);

    Ok(Alignment {
        score: best_score,
        query: query.clone(),
        aligned_query,
        aligned_positions,
        start_pos,
        end_pos,
    })
}

fn traceback(
    dp: &[Vec<Cell>],
    query: &[u8],
    matrix: &ScoringMatrix,
    query_len: usize,
    cons_len: usize,
) -> (String, Vec<Option<u32>>) {
    let mut aligned_query = String::new();
    let mut aligned_positions = Vec::new();

    let mut i = query_len;
    let mut j = cons_len;

    // For semi-global alignment, handle unaligned suffix of query
    let query_total_len = query.len();
    if i < query_total_len {
        // Add remaining query residues as insertions at the end
        for k in (i..query_total_len).rev() {
            aligned_query.push(query[k] as char);
            aligned_positions.push(None);
        }
    }

    while i > 0 && j > 0 {
        match dp[i][j].direction {
            Direction::Match => {
                aligned_query.push(query[i - 1] as char);
                aligned_positions.push(Some(matrix.positions[j - 1].position));
                i -= 1;
                j -= 1;
            }
            Direction::GapInQuery => {
                aligned_query.push('-');
                aligned_positions.push(Some(matrix.positions[j - 1].position));
                j -= 1;
            }
            Direction::GapInConsensus => {
                aligned_query.push(query[i - 1] as char);
                aligned_positions.push(None);
                i -= 1;
            }
        }
    }

    // Handle unaligned prefix of query
    while i > 0 {
        aligned_query.push(query[i - 1] as char);
        aligned_positions.push(None);
        i -= 1;
    }

    // Reverse (we built it backwards)
    aligned_query = aligned_query.chars().rev().collect();
    aligned_positions.reverse();

    (aligned_query, aligned_positions)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Chain;

    #[test]
    fn test_align_simple() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

        // A simple sequence that should align
        let sequence =
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN";

        let result = align(sequence, &matrix).unwrap();

        // Alignment should produce a score (may be negative for partial sequences)
        assert!(!result.aligned_query.is_empty());
        assert_eq!(result.aligned_positions.len(), result.aligned_query.len());
        assert!(result.start_pos > 0);
        assert!(result.end_pos > 0);
    }

    #[test]
    fn test_numbering() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let sequence =
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN";

        let result = align(sequence, &matrix).unwrap();
        let numbering = result.get_imgt_numbering();

        assert_eq!(numbering.len(), sequence.len());
        assert!(numbering[0].number > 0);
    }

    #[test]
    fn test_empty_sequence() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let result = align("", &matrix);
        assert!(result.is_err());
    }

    #[test]
    fn test_insertion_letter_conversion() {
        assert_eq!(insertion_to_letter(1), 'A');
        assert_eq!(insertion_to_letter(2), 'B');
        assert_eq!(insertion_to_letter(26), 'Z');
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
            let result = align(seq, &matrix).unwrap();
            let numbering = result.get_imgt_numbering();

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
        let result = align(partial_seq, &matrix).unwrap();

        // Should have a reasonable score (not heavily penalized for missing start)
        assert!(
            result.score > -100.0,
            "Partial sequence should not be heavily penalized"
        );
        // Allow some gaps but not excessive (e.g., <10% of sequence length)
        let gap_count = result.aligned_query.matches('-').count();
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
        let result = align(partial_seq, &matrix).unwrap();

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
        let result = align(fragment, &matrix).unwrap();

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
