//! Sequence alignment using Needleman-Wunsch algorithm

use serde::{Deserialize, Serialize};

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
    /// Evidence for the chain: at each aligned position that counts for confidence, the score of
    /// the residue there minus what a random residue scores there on average
    pub confidence_score: f32,
    /// The same evidence for a typical sequence of the chain covering the whole domain: at every
    /// position that counts for confidence, its typical score minus its random score. Confidence
    /// is `confidence_score / max_confidence_score`, clamped to 0–1
    pub max_confidence_score: f32,
    /// 0-based index of the first antibody residue in the query (0 when no prefix)
    pub query_start: usize,
    /// 0-based index of the last antibody residue in the query (query.len()-1 when no suffix)
    pub query_end: usize,
}

/// Most scoring matrices aligned in one pass, one per lane: eight `f32` lanes fill a 256-bit vector
const MAX_LANES: usize = 8;

/// The score row for a byte outside A-Z, after the rows for A to Z
const OTHER_LETTER: usize = 26;

/// Scoring matrices packed to align a query against all of them in one pass.
///
/// Filling the alignment matrix is latency-bound: each cell waits for the one to its left. Matrices
/// placed in separate lanes of the same cells don't depend on each other, so aligning against up to
/// [`MAX_LANES`] of them takes about as long as against one. Each lane does exactly the arithmetic a
/// lone alignment does, so its result is the same.
///
/// Only the scores and penalties are packed; the matrices stay with whoever built this, which passes
/// them to [`align_best`] alongside it.
#[derive(Debug, Clone)]
pub(crate) struct PackedScoringMatrices {
    groups: Vec<LaneGroup>,
}

impl PackedScoringMatrices {
    pub(crate) fn new<'a>(matrices: impl IntoIterator<Item = &'a [PositionScores]>) -> Self {
        let matrices: Vec<&[PositionScores]> = matrices.into_iter().collect();
        Self {
            groups: matrices.chunks(MAX_LANES).map(LaneGroup::new).collect(),
        }
    }
}

/// Up to [`MAX_LANES`] matrices, column by column, one lane per matrix. A shorter matrix is padded to
/// the longest with scores and penalties of -∞, which make every padding cell -∞: a cell depends
/// only on cells to its left and above it, so padding never reaches a real cell, and it never beats
/// a lane's row maximum, which keeps that and the traceback within the lane's own columns.
#[derive(Debug, Clone)]
struct LaneGroup {
    /// 1, 4 or [`MAX_LANES`]; the lanes after the first `filled` are padding
    lanes: usize,
    /// How many of the lanes hold a matrix
    filled: usize,
    columns: usize,
    /// `[(letter * columns + column) * lanes + lane]`, the letter being `byte - b'A'` or
    /// [`OTHER_LETTER`]
    scores: Vec<f32>,
    /// `[column * lanes + lane]`
    gap_penalties: Vec<f32>,
    /// `[column * lanes + lane]`
    insertion_penalties: Vec<f32>,
}

impl LaneGroup {
    fn new(matrices: &[&[PositionScores]]) -> Self {
        let lanes = match matrices.len() {
            1 => 1,
            2..=4 => 4,
            _ => MAX_LANES,
        };
        let columns = matrices.iter().map(|m| m.len()).max().unwrap_or(0);
        let padding = f32::NEG_INFINITY;
        let mut scores = vec![padding; (OTHER_LETTER + 1) * columns * lanes];
        let mut gap_penalties = vec![padding; columns * lanes];
        let mut insertion_penalties = vec![padding; columns * lanes];
        for (lane, positions) in matrices.iter().enumerate() {
            for (column, position) in positions.iter().enumerate() {
                for letter in 0..=OTHER_LETTER {
                    scores[(letter * columns + column) * lanes + lane] =
                        position.score_for(b'A' + letter as u8);
                }
                gap_penalties[column * lanes + lane] = position.gap_penalty;
                insertion_penalties[column * lanes + lane] = position.insertion_penalty;
            }
        }
        Self {
            lanes,
            filled: matrices.len(),
            columns,
            scores,
            gap_penalties,
            insertion_penalties,
        }
    }
}

/// Reusable buffer for the alignment matrices, to avoid allocating them for every alignment
pub struct AlignBuffer {
    /// `[(row * (columns + 1) + column) * lanes + lane]`
    dp_scores: Vec<f32>,
    /// [`Direction`]s, laid out as `dp_scores`
    dp_traceback: Vec<u8>,
    /// Each row's best score, `[row * lanes + lane]`
    row_max: Vec<f32>,
    /// The column of each row's best score, laid out as `row_max`
    row_best_column: Vec<u32>,
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
            row_max: Vec::new(),
            row_best_column: Vec::new(),
        }
    }

    fn ensure_capacity(&mut self, cells: usize, rows: usize) {
        // Exactly what is asked, so the memory held is what the longest query needed. The old
        // contents aren't copied: a fill writes every cell before reading it.
        if self.dp_scores.len() < cells {
            self.dp_scores = vec![0.0; cells];
            self.dp_traceback = vec![0; cells];
        }
        if self.row_max.len() < rows {
            self.row_max = vec![0.0; rows];
            self.row_best_column = vec![0; rows];
        }
    }

    /// Free the matrices if they take more than `cells` cells, counting every lane
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
/// Returns the optimal [`Alignment`] with scored positions.
pub fn align(
    query: &str,
    positions: &[PositionScores],
    align_buffer: Option<&mut AlignBuffer>,
) -> Alignment {
    let mut local_buf;
    let buf = match align_buffer {
        Some(buf) => buf,
        None => {
            local_buf = AlignBuffer::new();
            &mut local_buf
        }
    };
    align_best(
        query,
        &PackedScoringMatrices::new([positions]),
        |_| positions,
        buf,
    )
    .1
}

/// The best alignment of `query` against the packed matrices, with the index of its matrix: the
/// highest confidence, or the first matrix of those that tie for it. Confidence is comparable
/// across matrices where the raw score isn't: a chain whose reference residues vary more scores
/// lower even on its own sequences. Each matrix's alignment is the one [`align`] makes against it
/// alone. `matrix(k)` is the `k`-th matrix `packed` was built from: the packed copy has the scores
/// and penalties the fill needs, the matrix itself the details of each position that finding the
/// end and the traceback read.
pub(crate) fn align_best<'a>(
    query: &str,
    packed: &PackedScoringMatrices,
    matrix: impl Fn(usize) -> &'a [PositionScores],
    buf: &mut AlignBuffer,
) -> (usize, Alignment) {
    let query = query.as_bytes();
    let mut best: Option<(usize, Alignment, f32)> = None;
    let mut first_matrix = 0;
    for group in &packed.groups {
        match group.lanes {
            1 => fill::<1>(query, group, buf),
            4 => fill::<4>(query, group, buf),
            _ => fill::<MAX_LANES>(query, group, buf),
        }
        // The group's best lane, the first on a tie
        let mut group_best: Option<(usize, AlignmentEnd, f32)> = None;
        for lane in 0..group.filled {
            let positions = matrix(first_matrix + lane);
            let view = Lane::new(buf, group, lane);
            let end = alignment_end(&view, positions, query.len());
            let confidence = confidence_at(&view, positions, query, &end);
            if group_best.as_ref().is_none_or(|(_, _, c)| confidence > *c) {
                group_best = Some((lane, end, confidence));
            }
        }
        let (lane, end, confidence) = group_best.expect("a group has a matrix");
        // Built now, while the matrices still hold this group
        if best.as_ref().is_none_or(|(_, _, c)| confidence > *c) {
            let alignment = traced_alignment(
                &Lane::new(buf, group, lane),
                matrix(first_matrix + lane),
                query,
                end,
            );
            best = Some((first_matrix + lane, alignment, confidence));
        }
        first_matrix += group.filled;
    }
    let (index, alignment, _) = best.expect("there is a matrix");
    (index, alignment)
}

/// Fill the score and traceback matrices of every lane in `group`, with each row's best cell
// Indexing every per-lane array by lane is what lets the lane loop vectorise
#[allow(clippy::needless_range_loop)]
fn fill<const L: usize>(query: &[u8], group: &LaneGroup, buf: &mut AlignBuffer) {
    let columns = group.columns;
    let stride = columns + 1;
    let rows = query.len() + 1;
    buf.ensure_capacity(rows * stride * L, rows * L);
    let scores = buf.dp_scores[..rows * stride * L].as_chunks_mut::<L>().0;
    let traceback = buf.dp_traceback[..rows * stride * L].as_chunks_mut::<L>().0;
    let row_max = buf.row_max[..rows * L].as_chunks_mut::<L>().0;
    let row_best_column = buf.row_best_column[..rows * L].as_chunks_mut::<L>().0;
    let letter_scores = group.scores.as_chunks::<L>().0;
    let gap_penalties = &group.gap_penalties.as_chunks::<L>().0[..columns];
    let insertion_penalties = &group.insertion_penalties.as_chunks::<L>().0[..columns];

    // Leading gaps in the query and the consensus are free: the first row and column score 0. Their
    // traceback is never read, as the traceback stops on reaching either.
    scores[..stride].fill([0.0; L]);
    row_max[0] = [f32::NEG_INFINITY; L];
    row_best_column[0] = [0; L];

    for (i, &residue) in query.iter().enumerate() {
        let letter = (residue.to_ascii_uppercase().wrapping_sub(b'A') as usize).min(OTHER_LETTER);
        let match_scores = &letter_scores[letter * columns..][..columns];
        let (above, below) = scores.split_at_mut((i + 1) * stride);
        let prev = &above[i * stride..][..stride];
        let curr = &mut below[..stride];
        let directions = &mut traceback[(i + 1) * stride..][..stride];
        curr[0] = [0.0; L];

        let mut left = [0.0f32; L];
        let mut best = [f32::NEG_INFINITY; L];
        let mut best_column = [0u32; L];
        for j in 0..columns {
            let mut cell = [0.0f32; L];
            let mut direction = [0u8; L];
            for l in 0..L {
                let from_match = prev[j][l] + match_scores[j][l];
                // Gap in query (skip consensus position)
                let from_gap_query = left[l] + gap_penalties[j][l];
                // Gap in consensus (insertion in query sequence)
                let from_gap_cons = prev[j + 1][l] + insertion_penalties[j][l];
                // The best of the three, a tie going to the match, then the gap in the query, then
                // the gap in the consensus. Each cell waits on the cell to its left; the better of
                // the match and the gap in the consensus needs only the row above, so it's chosen
                // first, off that wait. Selects rather than branches, so the lane loop vectorises.
                let match_first = from_match >= from_gap_cons;
                let from_above = if match_first {
                    from_match
                } else {
                    from_gap_cons
                };
                let take_above = if match_first {
                    from_above >= from_gap_query
                } else {
                    from_above > from_gap_query
                };
                cell[l] = if take_above {
                    from_above
                } else {
                    from_gap_query
                };
                direction[l] = if !take_above {
                    Direction::GapInQuery.as_u8()
                } else if match_first {
                    Direction::Match.as_u8()
                } else {
                    Direction::GapInConsensus.as_u8()
                };
                let better = cell[l] > best[l];
                best[l] = if better { cell[l] } else { best[l] };
                best_column[l] = if better { j as u32 + 1 } else { best_column[l] };
            }
            curr[j + 1] = cell;
            directions[j + 1] = direction;
            left = cell;
        }
        row_max[i + 1] = best;
        row_best_column[i + 1] = best_column;
    }
}

/// One lane of the filled matrices
struct Lane<'a> {
    buf: &'a AlignBuffer,
    stride: usize,
    lanes: usize,
    lane: usize,
}

impl<'a> Lane<'a> {
    fn new(buf: &'a AlignBuffer, group: &LaneGroup, lane: usize) -> Self {
        Self {
            buf,
            stride: group.columns + 1,
            lanes: group.lanes,
            lane,
        }
    }

    fn score(&self, i: usize, j: usize) -> f32 {
        self.buf.dp_scores[(i * self.stride + j) * self.lanes + self.lane]
    }

    fn direction(&self, i: usize, j: usize) -> Direction {
        Direction::from_u8(self.buf.dp_traceback[(i * self.stride + j) * self.lanes + self.lane])
    }

    fn row_max(&self, i: usize) -> f32 {
        self.buf.row_max[i * self.lanes + self.lane]
    }

    fn row_best_column(&self, i: usize) -> usize {
        self.buf.row_best_column[i * self.lanes + self.lane] as usize
    }
}

/// Where an alignment ends, and its score
struct AlignmentEnd {
    score: f32,
    i: usize,
    j: usize,
}

fn alignment_end(lane: &Lane, positions: &[PositionScores], query_len: usize) -> AlignmentEnd {
    let cons_len = positions.len();

    // Each row's best cell enables true semi-global alignment with free end gaps on all four ends:
    // - Free leading query/consensus gaps: via first row/column = 0 initialization
    // - Free trailing consensus gaps: by picking best j per row (not forcing j = cons_len)
    // - Free trailing query gaps: by picking the row i with the highest per-row max score
    //   (if trailing query residues degrade the score, best_i < query_len)
    //
    // Default: full query consumption (standard semi-global behavior).
    // Override only if early termination gives meaningfully better score (clips suffix).
    // The threshold prevents clipping valid antibody endings that score slightly negative.
    const SUFFIX_CLIP_THRESHOLD: f32 = 50.0;
    let mut best_i = query_len;
    let mut best_j = lane.row_best_column(query_len);
    for i in 1..query_len {
        if lane.row_max(i) > lane.row_max(query_len) + SUFFIX_CLIP_THRESHOLD
            && lane.row_max(i) > lane.row_max(best_i)
        {
            best_i = i;
            best_j = lane.row_best_column(i);
        }
    }
    let mut score = lane.row_max(best_i);

    // Drop trailing residues inserted after the last position or matched to an unoccupied one
    // (gap_penalty == 0): a single such residue costs < SUFFIX_CLIP_THRESHOLD and would otherwise
    // be numbered. Take the score from where the walk stops so it matches the kept alignment.
    let (mut i, mut j) = (best_i, best_j);
    while i > 0 && j > 0 {
        let unoccupied = positions[j - 1].gap_penalty == 0.0;
        match lane.direction(i, j) {
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
        score = lane.score(i, j);
    }
    AlignmentEnd {
        score,
        i: best_i,
        j: best_j,
    }
}

/// Evidence a residue gives for the chain at a position: its score minus what a random residue
/// scores there on average, or nothing at a position that doesn't count for confidence
fn evidence(pos: &PositionScores, residue: u8) -> f32 {
    if pos.counts_for_confidence {
        pos.score_for(residue.to_ascii_uppercase()) - pos.random_score
    } else {
        0.0
    }
}

/// The evidence a typical sequence of the chain gives over the whole domain
fn typical_evidence(positions: &[PositionScores]) -> f32 {
    positions
        .iter()
        .filter(|pos| pos.counts_for_confidence)
        .map(|pos| pos.typical_score - pos.random_score)
        .sum()
}

/// The confidence of the alignment ending at `end`, unclamped, without building it: positions the
/// alignment doesn't reach, or skips, add no evidence and count as random
fn confidence_at(
    lane: &Lane,
    positions: &[PositionScores],
    query: &[u8],
    end: &AlignmentEnd,
) -> f32 {
    let typical = typical_evidence(positions);
    if typical <= 0.0 {
        return 0.0;
    }
    let (mut i, mut j, mut found) = (end.i, end.j, 0.0f32);
    while i > 0 && j > 0 {
        match lane.direction(i, j) {
            Direction::Match => {
                found += evidence(&positions[j - 1], query[i - 1]);
                i -= 1;
                j -= 1;
            }
            Direction::GapInQuery => j -= 1,
            Direction::GapInConsensus => i -= 1,
        }
    }
    found / typical
}

fn traced_alignment(
    lane: &Lane,
    positions: &[PositionScores],
    query_bytes: &[u8],
    end: AlignmentEnd,
) -> Alignment {
    let query_len = query_bytes.len();
    // positions.len() == query_len always; trailing suffix residues become Insertion()
    let mut aligned_positions = Vec::with_capacity(query_len);

    // Positions the alignment doesn't reach, or skips, add no evidence: they count as random
    let mut confidence_score = 0.0f32;
    let max_confidence_score = typical_evidence(positions);

    // Tracked during backward traversal: first Aligned seen = cons_end, last = cons_start
    let mut cons_start: u8 = 1;
    let mut cons_end: u8 = 128;
    let mut found_aligned = false;

    let mut i = end.i;
    let mut j = end.j;

    while i > 0 && j > 0 {
        match lane.direction(i, j) {
            Direction::Match => {
                let pos = &positions[j - 1];
                aligned_positions.push(AlignedPosition::Aligned(pos.position));

                confidence_score += evidence(pos, query_bytes[i - 1]);

                // Going backwards: first Aligned = cons_end, every Aligned updates cons_start
                if !found_aligned {
                    cons_end = pos.position;
                    found_aligned = true;
                }
                cons_start = pos.position;

                i -= 1;
                j -= 1;
            }

            Direction::GapInQuery => {
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

    // Append trailing Insertion() for suffix residues beyond the alignment's end
    let query_end = end.i.saturating_sub(1);
    for _ in end.i..query_len {
        aligned_positions.push(AlignedPosition::Insertion());
    }

    Alignment {
        score: end.score,
        positions: aligned_positions,
        cons_start,
        cons_end,
        confidence_score,
        max_confidence_score,
        query_start,
        query_end,
    }
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
    fn aligning_against_several_matrices_at_once_matches_aligning_against_each() {
        use Chain::*;
        // Every lane width, a group with padding lanes and consensuses of different lengths, and
        // more matrices than lanes; the repeated IGH and IGK tie with the first ones, which must win
        let sets: [&[Chain]; 5] = [
            &[IGH],
            &[IGK, IGL],
            &[IGH, IGK, IGL],
            &[IGH, IGK, IGL, TRA, TRB, TRG, TRD],
            &[TRA, IGH, IGK, IGL, TRB, TRG, TRD, IGH, IGK],
        ];
        let mut queries = Vec::new();
        for fixture in ["ab_H", "ab_K", "ab_L", "tcr_A", "tcr_B", "tcr_G", "tcr_D"] {
            let path = std::path::PathBuf::from(format!("fixtures/validation/{fixture}_imgt.csv"));
            for entry in crate::load_validation_csv(&path)
                .unwrap()
                .into_iter()
                .take(25)
            {
                // Whole, and cut short at both ends so alignments clip and start inside the domain
                let s = entry.sequence;
                queries.push(s[s.len() / 4..s.len() * 3 / 4].to_string());
                queries.push(s);
            }
        }
        let mut buf = AlignBuffer::new();
        for chains in sets {
            let matrices: Vec<ScoringMatrix> =
                chains.iter().map(|&c| ScoringMatrix::load(c)).collect();
            let packed = PackedScoringMatrices::new(matrices.iter().map(|m| &m.positions[..]));
            for query in &queries {
                let alone: Vec<Alignment> = matrices
                    .iter()
                    .map(|m| align(query, &m.positions, None))
                    .collect();
                // The highest confidence wins, the first matrix on a tie
                let confidence = |a: &Alignment| a.confidence_score / a.max_confidence_score;
                let mut expected = 0;
                for (k, alignment) in alone.iter().enumerate() {
                    if confidence(alignment) > confidence(&alone[expected]) {
                        expected = k;
                    }
                }
                let (index, alignment) =
                    align_best(query, &packed, |k| &matrices[k].positions[..], &mut buf);
                assert_eq!(
                    index, expected,
                    "{chains:?} chose a different matrix for {query}"
                );
                assert_eq!(
                    serde_json::to_value(&alignment).unwrap(),
                    serde_json::to_value(&alone[expected]).unwrap(),
                    "{chains:?} aligned {query} differently"
                );
            }
        }
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
}
