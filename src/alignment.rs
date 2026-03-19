//! Sequence alignment using Needleman-Wunsch algorithm

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};
use crate::scoring::PositionScores;

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Direction in the alignment traceback matrix
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
enum Direction {
    Match = 0,
    GapInQuery = 1,
    GapInConsensus = 2,
}

impl Direction {
    #[inline(always)]
    fn as_u8(self) -> u8 {
        self as u8
    }

    #[inline(always)]
    fn from_u8(v: u8) -> Self {
        match v {
            0 => Direction::Match,
            1 => Direction::GapInQuery,
            _ => Direction::GapInConsensus,
        }
    }
}

/// Represents a position in the alignment result.
/// The vector of these is always query-indexed (length == query length).
/// Consensus positions skipped by the aligner are not represented here;
/// the covered consensus range is captured by `Alignment::start_pos`/`end_pos`.
#[cfg_attr(feature = "python", pyclass(get_all))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AlignedPosition {
    /// Aligned to consensus position N
    Aligned(u8),
    /// Insertion in query (query has residue, consensus doesn't)
    Insertion(),
}

/// Result of aligning a sequence to a consensus
#[cfg_attr(feature = "python", pyclass(get_all))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Alignment {
    /// Alignment score
    pub score: f32,
    /// Aligned positions (merged representation of query and consensus alignment)
    pub positions: Vec<AlignedPosition>,
    /// Start position in consensus
    pub start_pos: u8,
    /// End position in consensus
    pub end_pos: u8,
    /// Sum of scores at confidence-relevant positions (occupancy > 0.5)
    pub confidence_score: f32,
    /// Sum of max possible scores at confidence-relevant positions
    pub max_confidence_score: f32,
}

/// Reusable buffer for alignment DP matrix to avoid repeated allocation
pub struct AlignBuffer {
    dp_scores: Vec<f32>,
    dp_traceback: Vec<u8>,
}

impl Default for AlignBuffer {
    fn default() -> Self {
        Self::new()
    }
}

impl AlignBuffer {
    pub fn new() -> Self {
        Self {
            dp_scores: Vec::new(),
            dp_traceback: Vec::new(),
        }
    }

    fn ensure_capacity(&mut self, total: usize) {
        if self.dp_scores.len() < total {
            self.dp_scores.resize(total, 0.0);
            self.dp_traceback.resize(total, 0);
        }
    }
}

/// Performs semi-global Needleman-Wunsch alignment of a query amino acid sequence
/// against position-specific scoring matrices derived from consensus sequences.
///
/// Leading/trailing gaps in the consensus are free (semi-global), but the full query
/// is consumed to prevent long CDR3 loops from being treated as trailing gaps.
/// Returns the optimal [`Alignment`] with scored positions, or an error if the sequence is empty.
pub fn align(
    query: &str,
    positions: &[PositionScores],
    align_buffer: Option<&mut AlignBuffer>,
) -> Result<Alignment> {
    let query_bytes = query.as_bytes();
    let query_len = query_bytes.len();
    let cons_len = positions.len();
    let stride = cons_len + 1;

    if query_len == 0 {
        return Err(Error::InvalidSequence("empty sequence".to_string()));
    }

    // Use provided buffer or allocate a local one
    let mut local_buf;
    let buf = match align_buffer {
        Some(buf) => buf,
        None => {
            local_buf = AlignBuffer::new();
            &mut local_buf
        }
    };

    // Reuse buffer, growing only if needed
    let total = (query_len + 1) * stride;
    buf.ensure_capacity(total);
    let dp_scores = &mut buf.dp_scores[..total];
    let dp_traceback = &mut buf.dp_traceback[..total];

    // Initialize first row: all zeros for scores, gaps in query are free (semi-global)
    dp_scores[0] = 0.0;
    dp_traceback[0] = Direction::Match.as_u8();
    for j in 1..stride {
        dp_scores[j] = 0.0;
        dp_traceback[j] = Direction::GapInQuery.as_u8();
    }

    // Initialize first column: gaps in consensus are free (semi-global)
    for i in 1..=query_len {
        dp_scores[i * stride] = 0.0;
        dp_traceback[i * stride] = Direction::GapInConsensus.as_u8();
    }

    // Fill DP matrix
    for i in 1..=query_len {
        let aa = query_bytes[i - 1].to_ascii_uppercase();
        let curr_row = i * stride;
        let prev_row = curr_row - stride;

        for j in 1..=cons_len {
            let cons_pos = &positions[j - 1];

            // Match/mismatch score using direct array index
            let match_score = cons_pos.score_for(aa);
            let from_match = dp_scores[prev_row + j - 1] + match_score;

            // Gap in query (skip consensus position)
            let from_gap_query = dp_scores[curr_row + j - 1] + cons_pos.gap_penalty;

            // Gap in consensus (insertion in query sequence)
            let from_gap_cons = dp_scores[prev_row + j] + cons_pos.insertion_penalty;

            // Choose best
            let (best_score, best_dir) =
                if from_match >= from_gap_query && from_match >= from_gap_cons {
                    (from_match, Direction::Match.as_u8())
                } else if from_gap_query >= from_gap_cons {
                    (from_gap_query, Direction::GapInQuery.as_u8())
                } else {
                    (from_gap_cons, Direction::GapInConsensus.as_u8())
                };

            dp_scores[curr_row + j] = best_score;
            dp_traceback[curr_row + j] = best_dir;
        }
    }

    // Find best ending position for semi-global alignment
    // Force full query consumption to prevent long CDR3s from being treated as trailing sequence
    let last_row = query_len * stride;
    let mut best_j = cons_len;
    let mut best_score = dp_scores[last_row + cons_len];

    // Check last row: query fully used, allow trailing gaps in consensus
    for j in 0..cons_len {
        if dp_scores[last_row + j] > best_score {
            best_score = dp_scores[last_row + j];
            best_j = j;
        }
    }

    // NOTE: We do NOT allow early termination of query (partial query alignment)
    // This prevents long CDR3 regions from being incorrectly treated as unaligned trailing sequence.
    // The query must align through its full length, but the consensus can have trailing gaps.

    let (aligned_positions, confidence_score, max_confidence_score, start_pos, end_pos) =
        build_traceback(
            dp_traceback,
            stride,
            query_bytes,
            positions,
            query_len,
            best_j,
        );

    Ok(Alignment {
        score: best_score,
        positions: aligned_positions,
        start_pos,
        end_pos,
        confidence_score,
        max_confidence_score,
    })
}

fn build_traceback(
    dp_traceback: &[u8],
    stride: usize,
    query_bytes: &[u8],
    positions: &[PositionScores],
    query_end: usize,
    cons_end: usize,
) -> (Vec<AlignedPosition>, f32, f32, u8, u8) {
    // Exactly one entry per query residue
    let mut aligned_positions = Vec::with_capacity(query_end);

    let mut confidence_score = 0.0f32;
    let mut max_confidence_score = 0.0f32;

    // Tracked during backward traversal: first Aligned seen = end_pos, last = start_pos
    let mut start_pos: u8 = 1;
    let mut end_pos: u8 = 128;
    let mut found_aligned = false;

    let mut i = query_end;
    let mut j = cons_end;

    while i > 0 && j > 0 {
        match Direction::from_u8(dp_traceback[i * stride + j]) {
            Direction::Match => {
                let pos = &positions[j - 1];
                aligned_positions.push(AlignedPosition::Aligned(pos.position));

                if pos.counts_for_confidence {
                    let aa = query_bytes[i - 1].to_ascii_uppercase();
                    confidence_score += pos.score_for(aa);
                    max_confidence_score += pos.max_score;
                }

                // Going backwards: first Aligned = end_pos, every Aligned updates start_pos
                if !found_aligned {
                    end_pos = pos.position;
                    found_aligned = true;
                }
                start_pos = pos.position;

                i -= 1;
                j -= 1;
            }

            Direction::GapInQuery => {
                let pos = &positions[j - 1];

                if pos.counts_for_confidence {
                    confidence_score += pos.gap_penalty;
                    max_confidence_score += pos.max_score;
                }

                // No query residue consumed — do not push to aligned_positions
                j -= 1;
            }

            Direction::GapInConsensus => {
                aligned_positions.push(AlignedPosition::Insertion());
                i -= 1;
            }
        }
    }

    while i > 0 {
        aligned_positions.push(AlignedPosition::Insertion());
        i -= 1;
    }

    aligned_positions.reverse();

    (
        aligned_positions,
        confidence_score,
        max_confidence_score,
        start_pos,
        end_pos,
    )
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::ScoringMatrix;
    use crate::Chain;

    fn test_align(query: &str, positions: &[PositionScores]) -> Result<Alignment> {
        align(query, positions, None)
    }
    #[test]
    fn test_align_simple() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

        // A simple sequence that should align
        let sequence =
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN";

        let result = test_align(sequence, &matrix.positions).unwrap();

        // Alignment should produce a score (may be negative for partial sequences)
        assert!(!result.positions.is_empty());
        assert!(result.start_pos > 0);
        assert!(result.end_pos > 0);
    }

    #[test]
    fn test_empty_sequence() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let result = test_align("", &matrix.positions);
        assert!(result.is_err());
    }

    #[test]
    fn test_partial_sequence_start() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

        // A sequence starting from middle of framework should align without heavy penalty
        let partial_seq = "GLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN";
        let result = test_align(partial_seq, &matrix.positions).unwrap();

        // Should have a reasonable score (not heavily penalized for missing start)
        assert!(
            result.score > -100.0,
            "Partial sequence should not be heavily penalized"
        );
        // Positions are query-indexed so length must equal query length
        assert_eq!(
            result.positions.len(),
            partial_seq.len(),
            "Positions length should equal query length"
        );
    }

    #[test]
    fn test_partial_sequence_end() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

        // A sequence ending in middle should align without heavy penalty
        let partial_seq = "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVS";
        let result = test_align(partial_seq, &matrix.positions).unwrap();

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
        let result = test_align(fragment, &matrix.positions).unwrap();

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
