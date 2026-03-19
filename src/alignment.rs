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
/// the covered consensus range is captured by `Alignment::cons_start`/`cons_end`.
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
    /// Aligned positions (length == query length; flanking residues are Insertion())
    pub positions: Vec<AlignedPosition>,
    /// First aligned consensus position (IMGT position number)
    pub cons_start: u8,
    /// Last aligned consensus position (IMGT position number)
    pub cons_end: u8,
    /// Sum of scores at confidence-relevant positions (occupancy > 0.5)
    pub confidence_score: f32,
    /// Sum of max possible scores at confidence-relevant positions
    pub max_confidence_score: f32,
    /// 0-based index of the first antibody residue in the query (0 when no prefix)
    pub query_start: usize,
    /// 0-based index of the last antibody residue in the query (query.len()-1 when no suffix)
    pub query_end: usize,
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
/// Leading/trailing gaps in the consensus and query are free (semi-global).
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

    // Track per-row maximum score and its column index during the fill.
    // This enables true semi-global alignment with free end gaps on all four ends:
    // - Free leading query/consensus gaps: via first row/column = 0 initialization
    // - Free trailing consensus gaps: by picking best j per row (not forcing j = cons_len)
    // - Free trailing query gaps: by picking the row i with the highest per-row max score
    //   (if trailing query residues degrade the score, best_i < query_len)
    let mut row_max_score = vec![f32::NEG_INFINITY; query_len + 1];
    let mut row_best_j = vec![0usize; query_len + 1];

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

            if best_score > row_max_score[i] {
                row_max_score[i] = best_score;
                row_best_j[i] = j;
            }
        }
    }

    // Default: full query consumption (standard semi-global behavior).
    // Override only if early termination gives meaningfully better score (clips suffix).
    // The threshold prevents clipping valid antibody endings that score slightly negative.
    const SUFFIX_CLIP_THRESHOLD: f32 = 50.0;
    let mut best_i = query_len;
    let mut best_j = row_best_j[query_len];
    for i in 1..query_len {
        if row_max_score[i] > row_max_score[query_len] + SUFFIX_CLIP_THRESHOLD
            && row_max_score[i] > row_max_score[best_i]
        {
            best_i = i;
            best_j = row_best_j[i];
        }
    }
    let best_score = row_max_score[best_i];

    let (
        aligned_positions,
        confidence_score,
        max_confidence_score,
        cons_start,
        cons_end,
        query_start,
        query_end,
    ) = build_traceback(
        dp_traceback,
        stride,
        query_bytes,
        positions,
        query_len,
        best_i,
        best_j,
    );

    Ok(Alignment {
        score: best_score,
        positions: aligned_positions,
        cons_start,
        cons_end,
        confidence_score,
        max_confidence_score,
        query_start,
        query_end,
    })
}

fn build_traceback(
    dp_traceback: &[u8],
    stride: usize,
    query_bytes: &[u8],
    positions: &[PositionScores],
    query_len: usize,
    traceback_end_i: usize,
    cons_end: usize,
) -> (Vec<AlignedPosition>, f32, f32, u8, u8, usize, usize) {
    // positions.len() == query_len always; trailing suffix residues become Insertion()
    let mut aligned_positions = Vec::with_capacity(query_len);

    let mut confidence_score = 0.0f32;
    let mut max_confidence_score = 0.0f32;

    // Tracked during backward traversal: first Aligned seen = cons_end, last = cons_start
    let mut cons_start: u8 = 1;
    let mut cons_end_pos: u8 = 128;
    let mut found_aligned = false;

    let mut i = traceback_end_i;
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

                // Going backwards: first Aligned = cons_end, every Aligned updates cons_start
                if !found_aligned {
                    cons_end_pos = pos.position;
                    found_aligned = true;
                }
                cons_start = pos.position;

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

    // i > 0 && j == 0: remaining query residues are the leading prefix (free via first column)
    let query_start = i;
    while i > 0 {
        aligned_positions.push(AlignedPosition::Insertion());
        i -= 1;
    }

    aligned_positions.reverse();

    // Append trailing Insertion() for suffix residues beyond traceback_end_i
    let query_end = traceback_end_i.saturating_sub(1);
    for _ in traceback_end_i..query_len {
        aligned_positions.push(AlignedPosition::Insertion());
    }

    (
        aligned_positions,
        confidence_score,
        max_confidence_score,
        cons_start,
        cons_end_pos,
        query_start,
        query_end,
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
        assert!(result.cons_start > 0);
        assert!(result.cons_end > 0);
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
            result.cons_start > 1,
            "Fragment should start after position 1"
        );
        assert!(
            result.cons_end < 128,
            "Fragment should end before last position"
        );
    }

    // Full IGH sequence (FR1 through FR4) from the task description
    const FULL_IGH: &str = "EVQLVESGGGLVQPGGSLRLSCAASGFNVSYSSIHWVRQAPGKGLEWVAYIYPSSGYTSYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARSYSTKLAMDYWGQGTLVTVSS";

    #[test]
    fn test_normal_antibody_no_flanking() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let result = test_align(FULL_IGH, &matrix.positions).unwrap();

        assert_eq!(result.query_start, 0, "No prefix: query_start should be 0");
        assert_eq!(
            result.query_end,
            FULL_IGH.len() - 1,
            "No suffix: query_end should be last index"
        );
        assert_eq!(result.positions.len(), FULL_IGH.len());
    }

    #[test]
    fn test_trailing_suffix() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let suffix = "AAAAAAA";
        let sequence = format!("{FULL_IGH}{suffix}");
        let result = test_align(&sequence, &matrix.positions).unwrap();

        assert_eq!(
            result.query_end,
            FULL_IGH.len() - 1,
            "Trailing suffix should be excluded: query_end should point to last antibody residue"
        );
        assert_eq!(result.query_start, 0);
        assert_eq!(result.positions.len(), sequence.len());
        // Suffix positions should be Insertion()
        for pos in &result.positions[FULL_IGH.len()..] {
            assert_eq!(*pos, AlignedPosition::Insertion());
        }
    }

    #[test]
    fn test_leading_prefix() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let prefix = "AAAAAAA";
        let sequence = format!("{prefix}{FULL_IGH}");
        let result = test_align(&sequence, &matrix.positions).unwrap();

        assert_eq!(
            result.query_start,
            prefix.len(),
            "Leading prefix should be identified: query_start should point past the prefix"
        );
        assert_eq!(result.query_end, sequence.len() - 1);
        assert_eq!(result.positions.len(), sequence.len());
    }

    #[test]
    fn test_both_flanking() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let prefix = "AAAAAAA";
        let suffix = "AAAAAAA";
        let sequence = format!("{prefix}{FULL_IGH}{suffix}");
        let result = test_align(&sequence, &matrix.positions).unwrap();

        assert_eq!(result.query_start, prefix.len());
        assert_eq!(result.query_end, prefix.len() + FULL_IGH.len() - 1);
        assert_eq!(result.positions.len(), sequence.len());
    }
}
