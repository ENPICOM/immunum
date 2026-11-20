//! Sequence alignment using Needleman-Wunsch algorithm
use std::{cell::RefCell, collections::HashSet, ops::Range, vec};

use serde::{Deserialize, Serialize};

use crate::{annotator::translate_codon, scoring::PositionScores};

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// How many nucleotides in a row are considered a homopoly run. Homopoly runs
/// are considered hotspots for sequencing errors. The penality for insertion/deletion
/// of nucleotide is decreased in such hotspots.
const MIN_NT_FOR_HOMOPOLY_RUN: usize = 4;

/// Direction in the alignment traceback matrix
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
enum Direction {
    Match = 0,
    GapInQuery = 1,
    GapInConsensus = 2,
    MatchPartial2nt = 3,
    MatchPartial4nt = 4,
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
            2 => Direction::GapInConsensus,
            3 => Direction::MatchPartial2nt,
            4 => Direction::MatchPartial4nt,
            _ => panic!("Should never happen"),
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
            _ => {
                panic!("This should NOT happen!");
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

/// Builds a mask over the query sequence marking nucleotide that are part of a homopoly run.
fn build_homopoly_mask(query: &str) -> Vec<bool> {
    let query_len = query.len();
    if query_len == 0 {
        return vec![false; 0];
    }
    let mut in_homopoly_run: Vec<bool> = vec![false; query_len];
    let mut index_start = 0;
    let mut last_base: char = query.bytes().nth(0).unwrap() as char;
    for (i, base) in query.chars().enumerate().skip(1) {
        if base != last_base {
            if i - index_start >= MIN_NT_FOR_HOMOPOLY_RUN {
                in_homopoly_run[index_start..i].fill(true);
            }
            index_start = i;
        }

        last_base = base;
    }
    if query_len - index_start >= MIN_NT_FOR_HOMOPOLY_RUN {
        in_homopoly_run[index_start..query_len].fill(true);
    }

    in_homopoly_run
}

/// Tries to infer missing nt based on homopoly mask.
/// Returns unknown residue X if it cannot infer it.
fn translate_quasicodon_2nt(
    query_bytes: &[u8],
    codon_start: usize,
    homopoly_run_mask: &Vec<bool>,
) -> u8 {
    let first = query_bytes[codon_start];
    let second = query_bytes[codon_start + 1];
    let first_homo = homopoly_run_mask[codon_start];
    let second_homo = homopoly_run_mask[codon_start + 1];
    // Both are in homopoly run
    if first_homo && second_homo {
        // We can only infer nt if both of the nts are in the same run
        if first == second {
            return translate_codon(&[first, second, second]);
        }
    } else if first_homo {
        return translate_codon(&[first, first, second]);
    } else if second_homo {
        return translate_codon(&[first, second, second]);
    }
    'X' as u8
}

/// Triest to infer codon by deleting a nt in homopoly run.
///
/// It builds every possible translation of resulting codon.
/// Returns unknown residue X if it cannot infer it.
fn translate_quasicodon_4nt(
    query_bytes: &[u8],
    codon_start: usize,
    homopoly_run_mask: &Vec<bool>,
) -> u8 {
    thread_local! {
        static RESIDUES: RefCell<HashSet<u8>> = RefCell::new(HashSet::new());
    }
    RESIDUES.with(|c| {
        let mut infered_residues = c.borrow_mut();
        infered_residues.clear();
        for i in 0..4 {
            // If not in homopoly run we do not consider deletion of nt
            if !homopoly_run_mask[codon_start + i] {
                continue;
            }

            let range1: Range<usize> = codon_start..codon_start + i;
            let range2 = codon_start + i + 1..codon_start + 4;
            println!("Range1: {:?} Range2: {:?}", range1, range2);
            let result: [u8; 3] = [&query_bytes[range1], &query_bytes[range2]]
                .concat()
                .try_into()
                .unwrap();

            let aa = translate_codon(&result);
            infered_residues.insert(aa);
        }

        //Only one residue can be infered from all options
        if infered_residues.len() == 1 {
            *infered_residues.iter().next().unwrap()
        } else {
            'X' as u8
        }
    })
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

    let homopoly_run_mask = build_homopoly_mask(query);

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

            //2nt partial match
            let from_2nt_codon = if row >= 2 {
                let aa = translate_quasicodon_2nt(query_bytes, row - 2, &homopoly_run_mask);
                let row_2_back = curr_row - 2 * stride;
                dp_scores[row_2_back + column - 1]
                    + cons_pos.nt_insertion_penalty
                    + cons_pos.score_for(aa)
            } else {
                f32::NEG_INFINITY
            };

            //4nt partial match
            let from_4nt_codon = if row >= 4 {
                let aa = translate_quasicodon_4nt(query_bytes, row - 4, &homopoly_run_mask);
                let row_4_back = curr_row - 4 * stride;
                dp_scores[row_4_back + column - 1]
                    + cons_pos.nt_deletion_penalty
                    + cons_pos.score_for(aa)
            } else {
                f32::NEG_INFINITY
            };

            //TODO: This is getting ugly. Find a better way.
            // Choose best
            let (best_score, best_dir) = if from_match >= from_gap_query
                && from_match >= from_gap_cons
                && from_match >= from_2nt_codon
                && from_match >= from_4nt_codon
            {
                (from_match, Direction::Match.as_u8())
            } else if from_gap_query >= from_gap_cons
                && from_gap_query >= from_2nt_codon
                && from_gap_query >= from_4nt_codon
            {
                (from_gap_query, Direction::GapInQuery.as_u8())
            } else if from_gap_cons >= from_2nt_codon && from_gap_cons >= from_4nt_codon {
                (from_gap_cons, Direction::GapInConsensus.as_u8())
            } else if from_2nt_codon > from_4nt_codon {
                (from_2nt_codon, Direction::MatchPartial2nt.as_u8())
            } else {
                (from_4nt_codon, Direction::MatchPartial4nt.as_u8())
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
                Direction::MatchPartial2nt if unoccupied => {
                    i -= 2;
                    j -= 1;
                }
                Direction::MatchPartial4nt if unoccupied => {
                    i -= 4;
                    j -= 1;
                }
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
        &homopoly_run_mask,
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
    homopoly_run_mask: &Vec<bool>,
) -> (Traceback, TracebackNT) {
    let mut aligned_positions: Vec<AlignedPosition> = Vec::with_capacity(query_len / 3 + 1);
    let mut codon_spans: Vec<Range<usize>> = Vec::with_capacity(query_len / 3 + 1);
    let mut translated: Vec<u8> = Vec::with_capacity(query_len / 3 + 1);

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
            Direction::MatchPartial2nt => {
                let pos: &PositionScores = &positions[column - 1];
                aligned_positions.push(AlignedPosition::Aligned(pos.position));
                let codon_span = row - 2..row;
                let aa = translate_quasicodon_2nt(query_bytes, row - 2, &homopoly_run_mask);
                codon_spans.push(codon_span);
                translated.push(aa);

                //TODO: what do we do with confidence?
                // if pos.counts_for_confidence {
                //     confidence_score += pos.score_for(aa);
                //     max_confidence_score += pos.max_score;
                // }

                // Going backwards: first Aligned = cons_end, every Aligned updates cons_start
                if !found_aligned {
                    cons_end_pos = pos.position;
                    found_aligned = true;
                }
                cons_start = pos.position;

                row -= 2;
                column -= 1;
            }
            Direction::MatchPartial4nt => {
                let pos: &PositionScores = &positions[column - 1];
                aligned_positions.push(AlignedPosition::Aligned(pos.position));
                let codon_span = row - 4..row;
                let aa = translate_quasicodon_4nt(query_bytes, row - 4, &homopoly_run_mask);
                codon_spans.push(codon_span);
                translated.push(aa);

                //TODO: what do we do with confidence?
                // if pos.counts_for_confidence {
                //     confidence_score += pos.score_for(aa);
                //     max_confidence_score += pos.max_score;
                // }

                // Going backwards: first Aligned = cons_end, every Aligned updates cons_start
                if !found_aligned {
                    cons_end_pos = pos.position;
                    found_aligned = true;
                }
                cons_start = pos.position;

                row -= 4;
                column -= 1;
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

    /// A consensus position that scores `aa` at 10 and every other residue at -5. Its NT penalties
    /// are the hotspot ones if `in_homopolymer_hotspot`.
    fn nt_test_position(position: u8, aa: u8, in_homopolymer_hotspot: bool) -> PositionScores {
        /// The amino acid gap and insertion penalty of every test position: the cost of skipping or
        /// inserting a whole codon
        const NT_TEST_CODON_GAP_PENALTY: f32 = -10.0;
        /// The cost of a two-base codon. Three of them cost more than skipping a codon, so a frame is
        /// only repaired where a base was lost, not to dodge a mismatch.
        const NT_TEST_INSERTION_PENALTY: f32 = -4.0;
        /// The cost of a two-base codon where a homopolymer run is expected: slightly less, so the
        /// repair lands there when several splits translate the same
        const NT_TEST_HOTSPOT_INSERTION_PENALTY: f32 = -3.5;

        let mut scores = [-5.0f32; 26];
        scores[(aa - b'A') as usize] = 10.0;
        let nt_penalty = if in_homopolymer_hotspot {
            NT_TEST_HOTSPOT_INSERTION_PENALTY
        } else {
            NT_TEST_INSERTION_PENALTY
        };
        PositionScores {
            position,
            scores,
            gap_penalty: NT_TEST_CODON_GAP_PENALTY,
            insertion_penalty: NT_TEST_CODON_GAP_PENALTY,
            max_score: 10.0,
            counts_for_confidence: true,
            nt_deletion_penalty: nt_penalty,
            nt_insertion_penalty: nt_penalty,
        }
    }

    /// One consensus position per residue of `consensus`, at IMGT positions 1, 2, ... The
    /// positions in `homopolymer_hotspots` get the lower NT penalties of a homopolymer hotspot.
    fn nt_test_positions(consensus: &[u8], homopolymer_hotspots: &[u8]) -> Vec<PositionScores> {
        consensus
            .iter()
            .zip(1u8..)
            .map(|(&aa, pos)| nt_test_position(pos, aa, homopolymer_hotspots.contains(&pos)))
            .collect()
    }

    #[test]
    fn test_align_nt_codon_spans_and_translation() {
        // Consensus QVQLV at IMGT positions 1..=5
        let positions = nt_test_positions(b"QVQLV", &[]);

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
        let positions = nt_test_positions(b"QVQLV", &[]);

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

    /// Aligns `read`, which has one frameshift (a base lost or gained), against `consensus` and
    /// checks that every position is numbered, that exactly one codon is a quasi-codon (not three
    /// bases long), and that it is one of `quasicodon_spans`. `homopolymer_hotspots` are the IMGT
    /// positions where the repair is expected; they get the lower NT penalties of a homopolymer
    /// hotspot.
    fn assert_frameshift_repaired(
        consensus: &[u8],
        homopolymer_hotspots: &[u8],
        read: &str,
        expected_translation: &str,
        quasicodon_spans: &[Range<usize>],
    ) -> AlignmentNT {
        let positions = nt_test_positions(consensus, homopolymer_hotspots);
        let (alignment, alignment_nt) = align_nt(read, &positions, None);

        let expected_positions: Vec<AlignedPosition> = (1..=consensus.len() as u8)
            .map(AlignedPosition::Aligned)
            .collect();
        assert_eq!(alignment.positions, expected_positions, "{read}");
        assert_eq!(alignment_nt.translated, expected_translation, "{read}");
        assert_eq!(alignment_nt.codon_spans.len(), consensus.len());

        let quasicodons: Vec<&Range<usize>> = alignment_nt
            .codon_spans
            .iter()
            .filter(|span| span.len() != 3)
            .collect();
        assert_eq!(quasicodons.len(), 1, "{:?}", alignment_nt.codon_spans);
        assert!(
            quasicodon_spans.contains(quasicodons[0]),
            "quasi-codon {:?}, expected one of {quasicodon_spans:?}",
            quasicodons[0]
        );
        // The spans tile the read without gaps or overlaps
        for pair in alignment_nt.codon_spans.windows(2) {
            assert_eq!(pair[0].end, pair[1].start);
        }
        alignment_nt
    }

    // Each test below lost one base from a run of at least 5, so the read still has a run of at
    // least MIN_NT_FOR_HOMOPOLY_RUN. Several splits of the read repair it; the hotspot position
    // picks the one whose two-base codon XY is the case under test.

    // XY with only X in the run is completed as XXY. True CAG GTG AAA AAC TTC CAG GTG (QVKNFQV)
    // has AAAAA; the read has AAAA at 6..10. AC at 9..11 becomes AAC (Asn).
    #[test]
    fn test_align_nt_deletion_first_base_in_run_completes_xxy() {
        let alignment_nt = assert_frameshift_repaired(
            b"QVKNFQV",
            &[4],
            "CAGGTGAAAACTTCCAGGTG",
            "QVKNFQV",
            &[9..11],
        );
        assert_eq!(
            alignment_nt.codon_spans,
            vec![0..3, 3..6, 6..9, 9..11, 11..14, 14..17, 17..20]
        );
    }

    // XY with only Y in the run is completed as XYY. True CAG GTG AGG GGG TTC CAG GTG (QVRGFQV)
    // has GGGGG; the read has GGGG at 7..11. AG at 6..8 becomes AGG (Arg).
    #[test]
    fn test_align_nt_deletion_second_base_in_run_completes_xyy() {
        let alignment_nt = assert_frameshift_repaired(
            b"QVRGFQV",
            &[3],
            "CAGGTGAGGGGTTCCAGGTG",
            "QVRGFQV",
            &[6..8],
        );
        assert_eq!(
            alignment_nt.codon_spans,
            vec![0..3, 3..6, 6..8, 8..11, 11..14, 14..17, 17..20]
        );
    }

    // XY with both in the same run: XXY and XYY are the same codon. True CAG GTG GGG GGG GTC TTC
    // CAG (QVGGVFQ) has eight Gs; the read has seven at 5..12. GG at 9..11 becomes GGG (Gly).
    #[test]
    fn test_align_nt_deletion_inside_run_completes_xxx() {
        let alignment_nt = assert_frameshift_repaired(
            b"QVGGVFQ",
            &[4],
            "CAGGTGGGGGGGTCTTCCAG",
            "QVGGVFQ",
            &[9..11],
        );
        assert_eq!(
            alignment_nt.codon_spans,
            vec![0..3, 3..6, 6..9, 9..11, 11..14, 14..17, 17..20]
        );
    }

    // XY with X and Y in two different runs: it is not known which run lost the base, so the
    // codon is X. True CAG GTG AAA ACG GGG TTC CAG (QVKTGFQ) lost the C between AAAA and GGGG;
    // the read has AAAAGGGG at 6..14. AG at 9..11 is X.
    #[test]
    fn test_align_nt_deletion_between_two_runs_is_x() {
        let alignment_nt = assert_frameshift_repaired(
            b"QVKTGFQ",
            &[4],
            "CAGGTGAAAAGGGGTTCCAG",
            "QVKXGFQ",
            &[9..11],
        );
        assert_eq!(
            alignment_nt.codon_spans,
            vec![0..3, 3..6, 6..9, 9..11, 11..14, 14..17, 17..20]
        );
    }

    // No run. True CAG GTG ACG TTC CAG GTG ATC (QVTFQVI); the read lost the C of ACG and has no
    // run longer than 2. The frame is repaired at AG (6..8), but the missing base can't be
    // inferred, so that codon is X. There is no run, so there is no hotspot.
    #[test]
    fn test_align_nt_deletion_outside_run_is_x() {
        let alignment_nt =
            assert_frameshift_repaired(b"QVTFQVI", &[], "CAGGTGAGTTCCAGGTGATC", "QVXFQVI", &[6..8]);
        assert_eq!(
            alignment_nt.codon_spans,
            vec![0..3, 3..6, 6..8, 8..11, 11..14, 14..17, 17..20]
        );
    }

    // Each test below gained one base, so the read has a four-base codon. Only a base in a run of
    // at least MIN_NT_FOR_HOMOPOLY_RUN may be dropped from it: if every such drop gives the same
    // residue that is the residue, otherwise X. Several splits of the read repair it; the hotspot
    // position picks the one whose four-base codon is the case under test.

    // All four bases in the run. True CAG GTG ACC GGG TTC CAG GTG (QVTGFQV) has GGG; the read
    // gained a G and has GGGG at 9..13. Dropping any G gives GGG (Gly).
    #[test]
    fn test_align_nt_insertion_inside_run_drops_run_base() {
        let alignment_nt = assert_frameshift_repaired(
            b"QVTGFQV",
            &[4],
            "CAGGTGACCGGGGTTCCAGGTG",
            "QVTGFQV",
            &[9..13],
        );
        assert_eq!(
            alignment_nt.codon_spans,
            vec![0..3, 3..6, 6..9, 9..13, 13..16, 16..19, 19..22]
        );
    }

    // The run ends inside the codon. True CAG GTG ACG GGA TTC CAG GTG (QVTGFQV) has GGG across
    // ACG|GGA; the read gained a G and has GGGG at 8..12. GGGA at 9..13: dropping any of its
    // three Gs gives GGA (Gly), and A is not in the run, so it is never dropped.
    #[test]
    fn test_align_nt_insertion_run_ends_in_codon_drops_run_base() {
        let alignment_nt = assert_frameshift_repaired(
            b"QVTGFQV",
            &[4],
            "CAGGTGACGGGGATTCCAGGTG",
            "QVTGFQV",
            &[9..13],
        );
        assert_eq!(
            alignment_nt.codon_spans,
            vec![0..3, 3..6, 6..9, 9..13, 13..16, 16..19, 19..22]
        );
    }

    // The run starts inside the codon. True CAG GTG ACC AGG GTC TTC CAG (QVTRVFQ) has GGG across
    // AGG|GTC; the read gained a G and has GGGG at 10..14. AGGG at 9..13: dropping any of its
    // three Gs gives AGG (Arg), and A is not in the run, so it is never dropped.
    #[test]
    fn test_align_nt_insertion_run_starts_in_codon_drops_run_base() {
        let alignment_nt = assert_frameshift_repaired(
            b"QVTRVFQ",
            &[4],
            "CAGGTGACCAGGGGTCTTCCAG",
            "QVTRVFQ",
            &[9..13],
        );
        assert_eq!(
            alignment_nt.codon_spans,
            vec![0..3, 3..6, 6..9, 9..13, 13..16, 16..19, 19..22]
        );
    }

    // Two runs meet inside the codon: it is not known which run gained the base, so the codon is
    // X. The read has AAAAGGGG at 7..15; AAGG at 9..13 gives AGG (Arg) when an A is dropped and
    // AAG (Lys) when a G is. Position 4 is Ser in the consensus, so neither guess matches and the
    // split that repairs the frame elsewhere doesn't outscore X.
    #[test]
    fn test_align_nt_insertion_between_two_runs_is_x() {
        let alignment_nt = assert_frameshift_repaired(
            b"QVQSGTV",
            &[4],
            "CAGGTGCAAAAGGGGTACCGTG",
            "QVQXGTV",
            &[9..13],
        );
        assert_eq!(
            alignment_nt.codon_spans,
            vec![0..3, 3..6, 6..9, 9..13, 13..16, 16..19, 19..22]
        );
    }

    // No run. True CAG GTG ACC TTC CAG GTG ATC (QVTFQVI); the read gained a G inside ACC and has
    // no run longer than 2. The frame is repaired at AGCC (6..10), but no base may be dropped, so
    // that codon is X. There is no run, so there is no hotspot.
    #[test]
    fn test_align_nt_insertion_outside_run_is_x() {
        let alignment_nt = assert_frameshift_repaired(
            b"QVTFQVI",
            &[],
            "CAGGTGAGCCTTCCAGGTGATC",
            "QVXFQVI",
            &[6..10],
        );
        assert_eq!(
            alignment_nt.codon_spans,
            vec![0..3, 3..6, 6..10, 10..13, 13..16, 16..19, 19..22]
        );
    }

    // Every base of a run of at least MIN_NT_FOR_HOMOPOLY_RUN identical bases is marked, and
    // nothing else. Each segment below is one run, and no two neighbouring segments share a base.
    #[test]
    fn test_build_homopoly_mask() {
        assert_eq!(
            MIN_NT_FOR_HOMOPOLY_RUN, 4,
            "the segments below assume runs of 4"
        );
        let segments: [(&str, bool); 12] = [
            ("AAAA", true), // run of exactly 4 at the very start
            ("C", false),   // single base between runs
            ("GGG", false), // one short of a run
            ("T", false),
            ("CCCCC", true), // longer than the minimum
            ("GGGG", true),  // run directly after another run, no base between
            ("TT", false),
            ("ACAC", false),    // alternating, never two alike
            ("GGGGGGGG", true), // long run
            ("C", false),
            ("AAA", false), // short run directly before the final run
            ("TTTT", true), // run of exactly 4 at the very end
        ];

        let query: String = segments.iter().map(|(run, _)| *run).collect();
        let expected: Vec<bool> = segments
            .iter()
            .flat_map(|(run, marked)| std::iter::repeat_n(*marked, run.len()))
            .collect();

        assert_eq!(build_homopoly_mask(&query), expected, "{query}");
    }

    // Sequences too short to hold the start and end cases separately
    #[test]
    fn test_build_homopoly_mask_whole_sequence() {
        assert_eq!(build_homopoly_mask(""), Vec::<bool>::new());
        assert_eq!(build_homopoly_mask("A"), vec![false]);
        assert_eq!(build_homopoly_mask("AAA"), vec![false; 3]);
        // The whole sequence is one run: it starts and ends at the sequence's ends
        assert_eq!(build_homopoly_mask("AAAA"), vec![true; 4]);
        assert_eq!(build_homopoly_mask("AAAAAA"), vec![true; 6]);
    }
}
