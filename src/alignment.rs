//! Sequence alignment using Needleman-Wunsch algorithm
use std::ops::Range;

use serde::{Deserialize, Serialize};

use crate::{annotator::translate_codon, scoring::PositionScores};

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
    /// Aligned positions
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

/// Result of tracing back through the DP matrix: everything in an [`Alignment`] except the score
struct Traceback {
    /// Aligned positions
    positions: Vec<AlignedPosition>,
    /// Sum of scores at confidence-relevant positions
    confidence_score: f32,
    /// Sum of max possible scores at confidence-relevant positions
    max_confidence_score: f32,
    /// First aligned consensus position (IMGT position number)
    cons_start: u8,
    /// Last aligned consensus position (IMGT position number)
    cons_end: u8,
    /// 0-based index of the first antibody residue in the query
    query_start: usize,
    /// 0-based index of the last antibody residue in the query
    query_end: usize,
}

impl Traceback {
    fn into_alignment(self, score: f32) -> Alignment {
        Alignment {
            score,
            positions: self.positions,
            cons_start: self.cons_start,
            cons_end: self.cons_end,
            confidence_score: self.confidence_score,
            max_confidence_score: self.max_confidence_score,
            query_start: self.query_start,
            query_end: self.query_end,
        }
    }
}

/// Result of aligning nt query to consensus
#[derive(Debug, Clone)]
pub struct AlignmentNT {
    pub codon_spans: Vec<Range<usize>>,
    pub translated: String,
}

struct TracebackNT {
    pub codon_spans: Vec<Range<usize>>,
    pub translated: Vec<u8>,
}

impl TracebackNT {
    fn into_alignment(self) -> AlignmentNT {
        let translated = self.translated.iter().map(|c| *c as char).collect();
        AlignmentNT {
            codon_spans: self.codon_spans,
            translated,
        }
    }
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
            // Exactly `total`, so the memory held is what the longest query needed. The old contents
            // aren't copied: `align` writes every cell before reading it.
            self.dp_scores = vec![0.0; total];
            self.dp_traceback = vec![0; total];
        }
    }

    /// Free the matrices if they take more than `cells` cells
    pub(crate) fn release_above(&mut self, cells: usize) {
        if self.dp_scores.capacity() > cells {
            *self = Self::new();
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
) -> Alignment {
    let query_bytes = query.as_bytes();
    let query_len = query_bytes.len();
    let cons_len = positions.len();
    let stride = cons_len + 1;

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
    // Pinned with J_MOTIF_BONUS in build.rs: at 35 the bovine ultralong CDR3 (4k3e,
    // `test_ultralong_cdr3_is_not_clipped`) is clipped before its FR4, and at 48 the lambda chain in
    // the Python tests numbered as kappa keeps the constant region's VGQPKA as FR4.
    const SUFFIX_CLIP_THRESHOLD: f32 = 42.0;
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

    let mut best_score = row_max_score[best_i];

    // Drop trailing residues inserted after the last position or matched to an unoccupied one
    // (gap_penalty == 0): a single such residue costs < SUFFIX_CLIP_THRESHOLD and would otherwise
    // be numbered. Take the score from where the walk stops so it matches the kept alignment.
    {
        let (mut i, mut j) = (best_i, best_j);
        while i > 0 && j > 0 {
            let unoccupied = positions[j - 1].gap_penalty == 0.0;
            match Direction::from_u8(dp_traceback[i * stride + j]) {
                Direction::GapInConsensus if j == cons_len => i -= 1,
                Direction::Match if unoccupied => {
                    i -= 1;
                    j -= 1;
                }
                Direction::GapInQuery if unoccupied => j -= 1,
                _ => break,
            }
        }
        if i != best_i {
            best_i = i;
            best_j = j;
            best_score = dp_scores[i * stride + j];
        }
    }

    let traceback = build_traceback(
        dp_traceback,
        stride,
        query_bytes,
        positions,
        query_len,
        best_i,
        best_j,
    );

    traceback.into_alignment(best_score)
}

fn build_traceback(
    dp_traceback: &[u8],
    stride: usize,
    query_bytes: &[u8],
    positions: &[PositionScores],
    query_len: usize,
    traceback_end_i: usize,
    cons_end: usize,
) -> Traceback {
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

    Traceback {
        positions: aligned_positions,
        confidence_score,
        max_confidence_score,
        cons_start,
        cons_end: cons_end_pos,
        query_start,
        query_end,
    }
}

/// Performs semi-global Needleman-Wunsch alignment of a query nucleotide sequence
/// against position-specific scoring matrices derived from consensus sequences.
/// Nucleotides translated on the fly.
///
/// Leading/trailing gaps in the consensus and query are free (semi-global).
/// Returns the optimal [`Alignment`] with scored positions, or an error if the sequence is empty.
pub fn align_nt(
    query: &str,
    positions: &[PositionScores],
    align_buffer: Option<&mut AlignBuffer>,
) -> (Alignment, AlignmentNT) {
    let query_bytes = query.as_bytes();
    let query_len = query_bytes.len();
    let cons_len = positions.len();
    let stride = cons_len + 1;

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
    let mut row_best_column = vec![0usize; query_len + 1];

    // Fill DP matrix
    for row in 1..=query_len {
        let curr_row = row * stride;

        let codon_3nt = if row >= 3 {
            Some(translate_codon(&query_bytes[row - 3..row]))
        } else {
            None
        };

        for column in 1..=cons_len {
            let cons_pos = &positions[column - 1];

            let from_match = if let Some(residue) = codon_3nt {
                let row_3_back = curr_row - 3 * stride;
                dp_scores[row_3_back + column - 1] + cons_pos.score_for(residue)
            } else {
                f32::NEG_INFINITY
            };

            // 3nt gap in query
            let from_gap_query = dp_scores[curr_row + column - 1] + cons_pos.gap_penalty;

            // 3nt insertion in query sequence
            let from_gap_cons = if row >= 3 {
                let row_3_back = curr_row - 3 * stride;
                dp_scores[row_3_back + column] + cons_pos.insertion_penalty
            } else {
                f32::NEG_INFINITY
            };

            // Choose best
            let (best_score, best_dir) =
                if from_match >= from_gap_query && from_match >= from_gap_cons {
                    (from_match, Direction::Match.as_u8())
                } else if from_gap_query >= from_gap_cons {
                    (from_gap_query, Direction::GapInQuery.as_u8())
                } else {
                    (from_gap_cons, Direction::GapInConsensus.as_u8())
                };

            dp_scores[curr_row + column] = best_score;
            dp_traceback[curr_row + column] = best_dir;

            if best_score > row_max_score[row] {
                row_max_score[row] = best_score;
                row_best_column[row] = column;
            }
        }
    }

    // Default: full or partial query consumption (standard semi-global behavior).
    // Override only if early termination gives meaningfully better score (clips suffix).
    // The threshold prevents clipping valid antibody endings that score slightly negative.
    // Pinned with J_MOTIF_BONUS in build.rs: at 35 the bovine ultralong CDR3 (4k3e,
    // `test_ultralong_cdr3_is_not_clipped`) is clipped before its FR4, and at 48 the lambda chain in
    // the Python tests numbered as kappa keeps the constant region's VGQPKA as FR4.
    const SUFFIX_CLIP_THRESHOLD: f32 = 42.0;

    // When we are clipping the nt sequence we do not know which of the last 3 nts
    // will have the best score. We consider 3 positions for default end.
    let mut best_row = query_len;
    for row in query_len - 2..query_len {
        if row_max_score[row] > row_max_score[best_row] {
            best_row = row;
        }
    }
    let mut best_column = row_best_column[best_row];
    let default_score = row_max_score[best_row];

    for i in 1..query_len - 2 {
        if row_max_score[i] > default_score + SUFFIX_CLIP_THRESHOLD
            && row_max_score[i] > row_max_score[best_row]
        {
            best_row = i;
            best_column = row_best_column[i];
        }
    }

    let mut best_score = row_max_score[best_row];

    // Drop trailing residues inserted after the last position or matched to an unoccupied one
    // (gap_penalty == 0): a single such residue costs < SUFFIX_CLIP_THRESHOLD and would otherwise
    // be numbered. Take the score from where the walk stops so it matches the kept alignment.
    {
        let (mut i, mut j) = (best_row, best_column);
        while i > 0 && j > 0 {
            let unoccupied = positions[j - 1].gap_penalty == 0.0;
            match Direction::from_u8(dp_traceback[i * stride + j]) {
                Direction::GapInConsensus if j == cons_len => i -= 3,
                Direction::Match if unoccupied => {
                    i -= 3;
                    j -= 1;
                }
                Direction::GapInQuery if unoccupied => j -= 1,
                _ => break,
            }
        }
        if i != best_row {
            best_row = i;
            best_column = j;
            best_score = dp_scores[i * stride + j];
        }
    }

    let (traceback, traceback_nt) = build_traceback_nt(
        dp_traceback,
        stride,
        query_bytes,
        positions,
        query_len,
        best_row,
        best_column,
    );

    (
        traceback.into_alignment(best_score),
        traceback_nt.into_alignment(),
    )
}

fn build_traceback_nt(
    dp_traceback: &[u8],
    stride: usize,
    query_bytes: &[u8],
    positions: &[PositionScores],
    query_len: usize,
    traceback_end_row: usize,
    cons_end: usize,
) -> (Traceback, TracebackNT) {
    let mut aligned_positions: Vec<AlignedPosition> = Vec::with_capacity(query_len / 3);
    let mut codon_spans: Vec<Range<usize>> = Vec::with_capacity(query_len / 3);
    let mut translated: Vec<u8> = Vec::with_capacity(query_len / 3);

    let mut confidence_score = 0.0f32;
    let mut max_confidence_score = 0.0f32;

    // Tracked during backward traversal: first Aligned seen = cons_end, last = cons_start
    let mut cons_start: u8 = 1;
    let mut cons_end_pos: u8 = 128;
    let mut found_aligned = false;

    let mut row = traceback_end_row;
    let mut column = cons_end;

    while row > 0 && column > 0 {
        match Direction::from_u8(dp_traceback[row * stride + column]) {
            Direction::Match => {
                let pos = &positions[column - 1];
                aligned_positions.push(AlignedPosition::Aligned(pos.position));
                let codon_span = row - 3..row;
                let aa = translate_codon(&query_bytes[codon_span.clone()]);
                codon_spans.push(codon_span);
                translated.push(aa);

                if pos.counts_for_confidence {
                    confidence_score += pos.score_for(aa);
                    max_confidence_score += pos.max_score;
                }

                // Going backwards: first Aligned = cons_end, every Aligned updates cons_start
                if !found_aligned {
                    cons_end_pos = pos.position;
                    found_aligned = true;
                }
                cons_start = pos.position;

                row -= 3;
                column -= 1;
            }

            Direction::GapInQuery => {
                let pos = &positions[column - 1];

                if pos.counts_for_confidence {
                    confidence_score += pos.gap_penalty;
                    max_confidence_score += pos.max_score;
                }

                // No query residue consumed — do not push to aligned_positions
                column -= 1;
            }

            Direction::GapInConsensus => {
                aligned_positions.push(AlignedPosition::Insertion());

                let codon_span = row - 3..row;
                let aa = translate_codon(&query_bytes[codon_span.clone()]);
                codon_spans.push(codon_span);
                translated.push(aa);

                row -= 3;
            }
            _ => {
                panic!("This should not happen!");
            }
        }
    }

    let query_start = row;

    aligned_positions.reverse();
    codon_spans.reverse();
    translated.reverse();

    let query_end = traceback_end_row.saturating_sub(1);
    (
        Traceback {
            positions: aligned_positions,
            confidence_score,
            max_confidence_score,
            cons_start,
            cons_end: cons_end_pos,
            query_start,
            query_end,
        },
        TracebackNT {
            codon_spans,
            translated,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::ScoringMatrix;
    use crate::Chain;

    fn test_align(query: &str, positions: &[PositionScores]) -> Alignment {
        align(query, positions, None)
    }
    #[test]
    fn test_align_simple() {
        let matrix = ScoringMatrix::load(Chain::IGH);

        // A simple sequence that should align
        let sequence =
            "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";

        let result = test_align(sequence, &matrix.positions);

        // Alignment should produce a score (may be negative for partial sequences)
        assert!(!result.positions.is_empty());
        assert!(result.cons_start > 0);
        assert!(result.cons_end > 0);
    }

    #[test]
    fn test_empty_sequence() {
        let matrix = ScoringMatrix::load(Chain::IGH);
        let result = test_align("", &matrix.positions);
        // Empty sequence should align with negative infinity score and no positions
        assert_eq!(result.score, f32::NEG_INFINITY);
        assert!(result.positions.is_empty());
    }

    #[test]
    fn test_partial_sequence_start() {
        let matrix = ScoringMatrix::load(Chain::IGH);

        // A sequence starting from middle of framework should align without heavy penalty
        let partial_seq = "GLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN";
        let result = test_align(partial_seq, &matrix.positions);

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
        let matrix = ScoringMatrix::load(Chain::IGH);

        // A sequence ending in middle should align without heavy penalty
        let partial_seq = "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVS";
        let result = test_align(partial_seq, &matrix.positions);

        // Should have a reasonable score
        assert!(
            result.score > -100.0,
            "Partial sequence should not be heavily penalized"
        );
    }

    #[test]
    fn test_fragment_alignment() {
        let matrix = ScoringMatrix::load(Chain::IGH);

        // A small fragment from the middle should align
        let fragment = "GLEWVSAISKSGGSTYY";
        let result = test_align(fragment, &matrix.positions);

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
        let matrix = ScoringMatrix::load(Chain::IGH);
        let result = test_align(FULL_IGH, &matrix.positions);

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
        let matrix = ScoringMatrix::load(Chain::IGH);
        let suffix = "AAAAAAA";
        let sequence = format!("{FULL_IGH}{suffix}");
        let result = test_align(&sequence, &matrix.positions);

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

    // PDB 1EFQ chain A. IGK consensus position 128 is unoccupied: no kappa has a residue there.
    const KAPPA: &str = "DIVMTQSPDSLAVSLGERATINCKSSQSVLYSSNSKNYLAWYQDKPGQPPKLLIYWASTRESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQYYSTPYSFGQGTKLEIK";

    #[test]
    fn test_one_trailing_residue_is_not_inserted_after_the_last_position() {
        let matrix = ScoringMatrix::load(Chain::IGH);
        for residue in ["A", "G", "R", "W"] {
            let sequence = format!("{FULL_IGH}{residue}");
            let result = test_align(&sequence, &matrix.positions);

            assert_eq!(result.query_end, FULL_IGH.len() - 1, "{residue}");
            assert_eq!(result.positions.len(), sequence.len());
            assert_eq!(
                result.positions[FULL_IGH.len()],
                AlignedPosition::Insertion()
            );
            assert_eq!(result.cons_end, 128);
        }
    }

    #[test]
    fn test_one_trailing_residue_is_not_placed_at_an_unoccupied_position() {
        let matrix = ScoringMatrix::load(Chain::IGK);
        for residue in ["A", "G", "R", "W"] {
            let sequence = format!("{KAPPA}{residue}");
            let result = test_align(&sequence, &matrix.positions);

            assert_eq!(result.query_end, KAPPA.len() - 1, "{residue}");
            assert_eq!(result.positions[KAPPA.len()], AlignedPosition::Insertion());
            assert_eq!(result.cons_end, 127);
        }
    }

    #[test]
    fn test_one_trailing_residue_aligns_the_domain_like_none() {
        for (chain, domain) in [(Chain::IGH, FULL_IGH), (Chain::IGK, KAPPA)] {
            let matrix = ScoringMatrix::load(chain);
            let bare = test_align(domain, &matrix.positions);
            let with_tail = test_align(&format!("{domain}A"), &matrix.positions);

            assert_eq!(with_tail.positions[..domain.len()], bare.positions[..]);
            assert_eq!(
                (
                    with_tail.score,
                    with_tail.confidence_score,
                    with_tail.max_confidence_score
                ),
                (bare.score, bare.confidence_score, bare.max_confidence_score)
            );
        }
    }

    #[test]
    fn test_leading_prefix() {
        let matrix = ScoringMatrix::load(Chain::IGH);
        let prefix = "AAAAAAA";
        let sequence = format!("{prefix}{FULL_IGH}");
        let result = test_align(&sequence, &matrix.positions);

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
        let matrix = ScoringMatrix::load(Chain::IGH);
        let prefix = "AAAAAAA";
        let suffix = "AAAAAAA";
        let sequence = format!("{prefix}{FULL_IGH}{suffix}");
        let result = test_align(&sequence, &matrix.positions);

        assert_eq!(result.query_start, prefix.len());
        assert_eq!(result.query_end, prefix.len() + FULL_IGH.len() - 1);
        assert_eq!(result.positions.len(), sequence.len());
    }

    // Bovine ultralong-CDR3 heavy chain (PDB 4k3e chain H): ~48 residues are absorbed as
    // insertions at IMGT 111/112. Terminating the alignment early is free, so the CDR-center
    // insertion penalty has to stay cheap enough that reaching FR4 still wins.
    const ULTRALONG_IGH: &str = "VQLRESGPSLVQPSQTLSLTCTASGFSKAVGWVRQAPGKALEWLGSIDTGGSTGYNPGLKSRLSITKDNSKSQVSLSVSSVTTEDSATYYCTTVHQETRKTCSDGYIAVDSCGRGQSDGCVNDCNSCYYGWRNCRRQPAIHSYEFHVDAWGRGLLVTVSS";

    #[test]
    fn test_ultralong_cdr3_is_not_clipped() {
        let matrix = ScoringMatrix::load(Chain::IGH);
        let result = test_align(ULTRALONG_IGH, &matrix.positions);

        assert_eq!(result.query_start, 0);
        assert_eq!(
            result.query_end,
            ULTRALONG_IGH.len() - 1,
            "Ultralong CDR3 must not be suffix-clipped: the C-terminal FR4 belongs to the domain"
        );
        assert_eq!(
            result.cons_end, 128,
            "Alignment should reach the end of FR4"
        );
    }

    /// A consensus position that scores `aa` at 10 and every other residue at -5
    fn nt_test_position(position: u8, aa: u8) -> PositionScores {
        let mut scores = [-5.0f32; 26];
        scores[(aa - b'A') as usize] = 10.0;
        PositionScores {
            position,
            scores,
            gap_penalty: -10.0,
            insertion_penalty: -10.0,
            max_score: 10.0,
            counts_for_confidence: true,
        }
    }

    #[test]
    fn test_align_nt_codon_spans_and_translation() {
        // Consensus QVQLV at IMGT positions 1..=5
        let positions: Vec<PositionScores> = b"QVQLV"
            .iter()
            .enumerate()
            .map(|(i, &aa)| nt_test_position(i as u8 + 1, aa))
            .collect();

        // 2 flanking bases, the codons CAG GTG CAG CTG GTG (QVQLV), 1 flanking base
        let query = "GGCAGGTGCAGCTGGTGT";

        let (alignment, alignment_nt) = align_nt(query, &positions, None);

        assert_eq!(alignment.cons_start, 1);
        assert_eq!(alignment.cons_end, 5);
        assert_eq!(alignment.query_start, 2);
        assert_eq!(alignment.query_end, 16);
        assert_eq!(alignment.score, 50.0);

        // Half-open base ranges of each aligned codon in the query
        assert_eq!(
            alignment_nt.codon_spans,
            vec![2..5, 5..8, 8..11, 11..14, 14..17]
        );
        for span in &alignment_nt.codon_spans {
            assert_eq!(span.len(), 3);
        }
        assert_eq!(alignment_nt.translated, "QVQLV");
    }

    // A codon the consensus has no position for (W between QVQ and LV) is an Insertion, and it
    // still has its bases and its residue: positions, codon spans and translation stay one-to-one
    #[test]
    fn test_align_nt_inserted_codon_has_span_and_residue() {
        // Consensus QVQLV at IMGT positions 1..=5
        let positions: Vec<PositionScores> = b"QVQLV"
            .iter()
            .enumerate()
            .map(|(i, &aa)| nt_test_position(i as u8 + 1, aa))
            .collect();

        // CAG GTG CAG TGG CTG GTG: QVQ, an extra W, LV
        let query = "CAGGTGCAGTGGCTGGTG";

        let (alignment, alignment_nt) = align_nt(query, &positions, None);

        // 5 matches (+50) and one inserted codon (-10)
        assert_eq!(alignment.score, 40.0);
        assert_eq!(
            alignment.positions,
            vec![
                AlignedPosition::Aligned(1),
                AlignedPosition::Aligned(2),
                AlignedPosition::Aligned(3),
                AlignedPosition::Insertion(),
                AlignedPosition::Aligned(4),
                AlignedPosition::Aligned(5),
            ]
        );
        assert_eq!(
            alignment_nt.codon_spans,
            vec![0..3, 3..6, 6..9, 9..12, 12..15, 15..18]
        );
        assert_eq!(alignment_nt.translated, "QVQWLV");
        assert_eq!(alignment_nt.codon_spans.len(), alignment.positions.len());
    }
}
