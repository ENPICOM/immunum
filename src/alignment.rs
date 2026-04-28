//! Sequence alignment using Needleman-Wunsch algorithm

use serde::{Deserialize, Serialize};

use crate::scoring::{PositionScores, ScoringMatrix};
use crate::types::Region;

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

/// Result of aligning a query window to one contiguous consensus region.
///
/// Coordinates are absolute indices into the full query.
#[derive(Debug, Clone)]
pub struct RegionAlignment {
    pub score: f32,
    /// First absolute query index covered by this region (inclusive).
    pub query_start: usize,
    /// Last absolute query index covered by this region (inclusive).
    /// Equals `query_start - 1` when the region matched zero query residues.
    pub query_end: usize,
    /// One [`AlignedPosition`] per residue in `query[query_start..=query_end]`.
    pub positions: Vec<AlignedPosition>,
    pub confidence_score: f32,
    pub max_confidence_score: f32,
    /// True iff at least one query residue was matched against a consensus position.
    pub has_alignment: bool,
}

/// Per-corner free-gap parameters for a region alignment.
#[derive(Debug, Clone, Copy)]
pub struct RegionAlignParams {
    /// Allow alignment to start at any position in the query window for free.
    pub free_query_prefix: bool,
    /// Allow alignment to end at any position in the query window for free.
    pub free_query_suffix: bool,
    /// Allow alignment to start at any consensus position for free.
    pub free_cons_prefix: bool,
    /// Allow alignment to end at any consensus position for free.
    pub free_cons_suffix: bool,
    /// Score margin a smaller-i endpoint must exceed full-query consumption by
    /// before we clip trailing query residues. Use 0 for inner FRs (window
    /// deliberately extends past region end), ~50 for FR4 (full query is the
    /// natural answer; only clip for clearly non-AB suffixes).
    pub suffix_clip_threshold: f32,
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

    Alignment {
        score: best_score,
        positions: aligned_positions,
        cons_start,
        cons_end,
        confidence_score,
        max_confidence_score,
        query_start,
        query_end,
    }
}

/// Align a query window to a contiguous slice of consensus positions (one region).
///
/// Free-gap behavior at each corner is controlled by [`RegionAlignParams`]. Returns a
/// [`RegionAlignment`] whose `positions` covers exactly `query[query_start..=query_end]`
/// (no leading/trailing `Insertion()` padding outside that range).
pub fn align_region(
    query: &[u8],
    window_start: usize,
    window_end: usize,
    region_positions: &[PositionScores],
    params: RegionAlignParams,
    buf: &mut AlignBuffer,
) -> RegionAlignment {
    debug_assert!(window_end >= window_start);
    debug_assert!(window_end <= query.len());

    let qlen = window_end - window_start;
    let cons_len = region_positions.len();
    let stride = cons_len + 1;

    if qlen == 0 || cons_len == 0 {
        return RegionAlignment {
            score: 0.0,
            query_start: window_start,
            query_end: window_start.saturating_sub(1),
            positions: Vec::new(),
            confidence_score: 0.0,
            max_confidence_score: 0.0,
            has_alignment: false,
        };
    }

    let total = (qlen + 1) * stride;
    buf.ensure_capacity(total);
    let dp_scores = &mut buf.dp_scores[..total];
    let dp_traceback = &mut buf.dp_traceback[..total];

    // First row (i=0): aligned 0 query residues, j consensus positions skipped.
    dp_scores[0] = 0.0;
    dp_traceback[0] = Direction::Match.as_u8();
    if params.free_cons_prefix {
        for j in 1..stride {
            dp_scores[j] = 0.0;
            dp_traceback[j] = Direction::GapInQuery.as_u8();
        }
    } else {
        let mut acc = 0.0;
        for j in 1..stride {
            acc += region_positions[j - 1].gap_penalty;
            dp_scores[j] = acc;
            dp_traceback[j] = Direction::GapInQuery.as_u8();
        }
    }

    // First column (j=0): i query residues as insertions in consensus.
    if params.free_query_prefix {
        for i in 1..=qlen {
            dp_scores[i * stride] = 0.0;
            dp_traceback[i * stride] = Direction::GapInConsensus.as_u8();
        }
    } else {
        // Charge the first region position's insertion penalty per skipped query residue.
        let pen = region_positions[0].insertion_penalty;
        let mut acc = 0.0;
        for i in 1..=qlen {
            acc += pen;
            dp_scores[i * stride] = acc;
            dp_traceback[i * stride] = Direction::GapInConsensus.as_u8();
        }
    }

    // DP fill (same recurrence as `align`).
    for i in 1..=qlen {
        let aa = query[window_start + i - 1].to_ascii_uppercase();
        let curr_row = i * stride;
        let prev_row = curr_row - stride;

        for j in 1..=cons_len {
            let cons_pos = &region_positions[j - 1];

            let from_match = dp_scores[prev_row + j - 1] + cons_pos.score_for(aa);
            let from_gap_query = dp_scores[curr_row + j - 1] + cons_pos.gap_penalty;
            let from_gap_cons = dp_scores[prev_row + j] + cons_pos.insertion_penalty;

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

    // Pick the best endpoint per the free-suffix params.
    let (best_i, best_j) = pick_endpoint(dp_scores, stride, qlen, cons_len, params);
    let best_score = dp_scores[best_i * stride + best_j];

    // Traceback.
    let (positions, confidence_score, max_confidence_score, query_start_in_window, has_alignment) =
        build_region_traceback(
            dp_traceback,
            stride,
            query,
            window_start,
            region_positions,
            best_i,
            best_j,
        );

    let query_start_abs = window_start + query_start_in_window;
    let query_end_abs = if best_i == 0 {
        // No query residue consumed (degenerate / unaligned).
        query_start_abs.saturating_sub(1)
    } else {
        window_start + best_i - 1
    };

    RegionAlignment {
        score: best_score,
        query_start: query_start_abs,
        query_end: query_end_abs,
        positions,
        confidence_score,
        max_confidence_score,
        has_alignment,
    }
}

// =============================================================================
// Sequential FR-anchored alignment
// =============================================================================

/// Maximum slack for the variable-length region between FR1 and FR2 (CDR1).
/// Empirical max in the validation set is around 12; 16 leaves room.
pub const CDR1_MAX_LEN: usize = 20;
/// Maximum slack for the variable-length region between FR2 and FR3 (CDR2).
pub const CDR2_MAX_LEN: usize = 20;
/// Maximum slack for the variable-length region between FR3 and FR4 (CDR3).
/// Tested against 30+; 40 covers all observed cases.
pub const CDR3_MAX_LEN: usize = 40;

/// Per-FR confidence floor. If FR1's alignment scores below this, the chain is
/// rejected immediately without aligning subsequent FRs. Set well below the
/// `DEFAULT_MIN_CONFIDENCE` (0.5) used for the full alignment so genuinely weak
/// antibody/TCR sequences still progress to FR2-FR4 (their full-alignment
/// confidence may pass the overall threshold even if FR1 alone is borderline);
/// the primary purpose of this floor is to skip the remaining FRs for clearly
/// non-immunoglobulin chains during multi-chain selection.
pub const FR1_CONFIDENCE_THRESHOLD: f32 = 0.15;

fn region_slice(matrix: &ScoringMatrix, region: Region) -> &[PositionScores] {
    let (start, end) = region.imgt_range();
    &matrix.positions[(start as usize - 1)..(end as usize)]
}

/// FR3 alignment extends 1 position into CDR3 to anchor IMGT 105 — needed because
/// Kabat heavy treats IMGT 105 as FR3-last while IMGT scheme treats it as CDR3-first.
fn fr3_extended_slice(matrix: &ScoringMatrix) -> &[PositionScores] {
    // IMGT 66..=105 → indices 65..=104 → slice 65..105
    &matrix.positions[65..105]
}

/// FR4 alignment extends 1 position upstream into CDR3 to anchor IMGT 117 — needed
/// because Kabat (both light and heavy) treats IMGT 117 as FR4-first while IMGT
/// scheme treats it as CDR3-last.
fn fr4_extended_slice(matrix: &ScoringMatrix) -> &[PositionScores] {
    // IMGT 117..=128 → indices 116..=127 → slice 116..128
    &matrix.positions[116..128]
}

/// Inner FR params: window deliberately extends past the region end into the
/// adjacent CDR; eager max-clipping is correct (threshold = 0).
const INNER_FR_PARAMS: RegionAlignParams = RegionAlignParams {
    free_query_prefix: true,
    free_query_suffix: true,
    free_cons_prefix: false,
    free_cons_suffix: false,
    suffix_clip_threshold: 0.0,
};

/// FR1 params: as inner FR but free_cons_prefix=true to tolerate truncated N-term.
const FR1_PARAMS: RegionAlignParams = RegionAlignParams {
    free_query_prefix: true,
    free_query_suffix: true,
    free_cons_prefix: true,
    free_cons_suffix: false,
    suffix_clip_threshold: 0.0,
};

/// FR4 params: full query consumption is the natural answer; only clip for
/// clearly non-AB suffixes (threshold = 50, mirroring the legacy `align()`).
/// Free cons suffix tolerates truncated C-term.
const FR4_PARAMS: RegionAlignParams = RegionAlignParams {
    free_query_prefix: true,
    free_query_suffix: true,
    free_cons_prefix: false,
    free_cons_suffix: true,
    suffix_clip_threshold: 50.0,
};

/// Align FR1 only — used for cheap chain triage. Returns the FR1
/// [`RegionAlignment`] which can be compared across chains before deciding which
/// chain to run the remaining FR2-FR4 alignments for.
pub fn align_fr1(query: &str, matrix: &ScoringMatrix, buf: &mut AlignBuffer) -> RegionAlignment {
    let qbytes = query.as_bytes();
    let qlen = qbytes.len();
    if qlen == 0 {
        return RegionAlignment {
            score: f32::NEG_INFINITY,
            query_start: 0,
            query_end: 0,
            positions: Vec::new(),
            confidence_score: 0.0,
            max_confidence_score: 0.0,
            has_alignment: false,
        };
    }
    let fr1_slice = region_slice(matrix, Region::FR1);
    let fr1_window_end = (fr1_slice.len() + CDR1_MAX_LEN).min(qlen);
    align_region(qbytes, 0, fr1_window_end, fr1_slice, FR1_PARAMS, buf)
}

/// Continue from a previously-computed FR1 alignment to produce the full
/// [`Alignment`] across all four FRs (with CDRs inferred from gaps).
///
/// Returns `None` if FR1 fails the per-FR confidence threshold, if any FR alignment
/// is degenerate, or if the resulting CDR lengths are out of bounds.
pub fn align_remaining_frs(
    query: &str,
    matrix: &ScoringMatrix,
    fr1: &RegionAlignment,
    buf: &mut AlignBuffer,
) -> Option<Alignment> {
    let qbytes = query.as_bytes();
    let qlen = qbytes.len();
    if qlen == 0 || !fr1.has_alignment {
        return None;
    }

    let fr1_conf = if fr1.max_confidence_score > 0.0 {
        (fr1.confidence_score / fr1.max_confidence_score).clamp(0.0, 1.0)
    } else {
        0.0
    };
    if fr1_conf < FR1_CONFIDENCE_THRESHOLD {
        return None;
    }

    let fr2_slice = region_slice(matrix, Region::FR2);
    let fr3_slice = fr3_extended_slice(matrix);
    let fr4_slice = fr4_extended_slice(matrix);

    let fr2_window_start = fr1.query_end + 1;
    if fr2_window_start >= qlen {
        return None;
    }
    let fr2_window_end = (fr2_window_start + fr2_slice.len() + CDR1_MAX_LEN).min(qlen);
    let fr2 = align_region(
        qbytes,
        fr2_window_start,
        fr2_window_end,
        fr2_slice,
        INNER_FR_PARAMS,
        buf,
    );
    if !fr2.has_alignment || fr2.query_start <= fr1.query_end {
        return None;
    }

    let fr3_window_start = fr2.query_end + 1;
    if fr3_window_start >= qlen {
        return None;
    }
    let fr3_window_end = (fr3_window_start + fr3_slice.len() + CDR2_MAX_LEN).min(qlen);
    let fr3 = align_region(
        qbytes,
        fr3_window_start,
        fr3_window_end,
        fr3_slice,
        INNER_FR_PARAMS,
        buf,
    );
    if !fr3.has_alignment || fr3.query_start <= fr2.query_end {
        return None;
    }

    let fr4_window_start = fr3.query_end + 1;
    if fr4_window_start >= qlen {
        return None;
    }
    let fr4 = align_region(qbytes, fr4_window_start, qlen, fr4_slice, FR4_PARAMS, buf);
    if !fr4.has_alignment || fr4.query_start <= fr3.query_end {
        return None;
    }

    let cdr1_len = fr2.query_start - fr1.query_end - 1;
    let cdr2_len = fr3.query_start - fr2.query_end - 1;
    let cdr3_len = fr4.query_start - fr3.query_end - 1;
    if cdr1_len > CDR1_MAX_LEN || cdr2_len > CDR2_MAX_LEN || cdr3_len > CDR3_MAX_LEN {
        return None;
    }

    Some(stitch_alignment(qlen, fr1, &fr2, &fr3, &fr4))
}

/// Stitch FR alignments into a single query-indexed `Alignment`.
///
/// CDR residues between FRs are filled with `Aligned(<CDR first IMGT pos>)` so that
/// downstream `apply_numbering` rules pick them up by IMGT range. Residues outside
/// the antibody region (before FR1, after FR4) are `Insertion()`.
fn stitch_alignment(
    qlen: usize,
    fr1: &RegionAlignment,
    fr2: &RegionAlignment,
    fr3: &RegionAlignment,
    fr4: &RegionAlignment,
) -> Alignment {
    let mut positions: Vec<AlignedPosition> = Vec::with_capacity(qlen);

    // Leading prefix (before FR1).
    for _ in 0..fr1.query_start {
        positions.push(AlignedPosition::Insertion());
    }

    // FR1 residues.
    positions.extend_from_slice(&fr1.positions);

    // CDR1: between fr1.query_end and fr2.query_start.
    for _ in (fr1.query_end + 1)..fr2.query_start {
        positions.push(AlignedPosition::Aligned(Region::CDR1.imgt_range().0));
    }

    // FR2.
    positions.extend_from_slice(&fr2.positions);

    // CDR2.
    for _ in (fr2.query_end + 1)..fr3.query_start {
        positions.push(AlignedPosition::Aligned(Region::CDR2.imgt_range().0));
    }

    // FR3.
    positions.extend_from_slice(&fr3.positions);

    // CDR3 middle (between FR3-ext anchor at IMGT 105 and FR4-ext anchor at IMGT 117).
    // Use IMGT 110 — a center CDR3 position contained by all rule sets:
    // IMGT (105..=117), Kabat heavy CDR3 (106..=117), Kabat light CDR3 (105..=116).
    const CDR3_MIDDLE_PLACEHOLDER: u8 = 110;
    for _ in (fr3.query_end + 1)..fr4.query_start {
        positions.push(AlignedPosition::Aligned(CDR3_MIDDLE_PLACEHOLDER));
    }

    // FR4.
    positions.extend_from_slice(&fr4.positions);

    // Trailing suffix (after FR4).
    for _ in (fr4.query_end + 1)..qlen {
        positions.push(AlignedPosition::Insertion());
    }

    debug_assert_eq!(positions.len(), qlen);

    let cons_start = first_aligned_imgt(&positions).unwrap_or(1);
    let cons_end = last_aligned_imgt(&positions).unwrap_or(128);

    Alignment {
        score: fr1.score + fr2.score + fr3.score + fr4.score,
        positions,
        cons_start,
        cons_end,
        confidence_score: fr1.confidence_score
            + fr2.confidence_score
            + fr3.confidence_score
            + fr4.confidence_score,
        max_confidence_score: fr1.max_confidence_score
            + fr2.max_confidence_score
            + fr3.max_confidence_score
            + fr4.max_confidence_score,
        query_start: fr1.query_start,
        query_end: fr4.query_end,
    }
}

fn first_aligned_imgt(positions: &[AlignedPosition]) -> Option<u8> {
    positions.iter().find_map(|p| match p {
        AlignedPosition::Aligned(n) => Some(*n),
        _ => None,
    })
}

fn last_aligned_imgt(positions: &[AlignedPosition]) -> Option<u8> {
    positions.iter().rev().find_map(|p| match p {
        AlignedPosition::Aligned(n) => Some(*n),
        _ => None,
    })
}

fn pick_endpoint(
    dp_scores: &[f32],
    stride: usize,
    qlen: usize,
    cons_len: usize,
    params: RegionAlignParams,
) -> (usize, usize) {
    let cell = |i: usize, j: usize| dp_scores[i * stride + j];

    // Per-row: best j and its score, respecting free_cons_suffix.
    let row_info = |i: usize| -> (f32, usize) {
        if params.free_cons_suffix {
            let mut best_j = cons_len;
            let mut best = cell(i, cons_len);
            for j in 0..cons_len {
                let s = cell(i, j);
                if s > best {
                    best = s;
                    best_j = j;
                }
            }
            (best, best_j)
        } else {
            (cell(i, cons_len), cons_len)
        }
    };

    let (mut best_score, mut best_j) = row_info(qlen);
    let mut best_i = qlen;

    if params.free_query_suffix {
        let baseline = best_score;
        for i in 1..qlen {
            let (s, j) = row_info(i);
            if s > baseline + params.suffix_clip_threshold && s > best_score {
                best_score = s;
                best_i = i;
                best_j = j;
            }
        }
    }

    (best_i, best_j)
}

#[allow(clippy::too_many_arguments)]
fn build_region_traceback(
    dp_traceback: &[u8],
    stride: usize,
    query: &[u8],
    window_start: usize,
    region_positions: &[PositionScores],
    end_i: usize,
    end_j: usize,
) -> (Vec<AlignedPosition>, f32, f32, usize, bool) {
    let mut positions = Vec::new();
    let mut confidence_score = 0.0f32;
    let mut max_confidence_score = 0.0f32;
    let mut has_alignment = false;

    let mut i = end_i;
    let mut j = end_j;

    while i > 0 && j > 0 {
        match Direction::from_u8(dp_traceback[i * stride + j]) {
            Direction::Match => {
                let pos = &region_positions[j - 1];
                positions.push(AlignedPosition::Aligned(pos.position));
                if pos.counts_for_confidence {
                    let aa = query[window_start + i - 1].to_ascii_uppercase();
                    confidence_score += pos.score_for(aa);
                    max_confidence_score += pos.max_score;
                }
                has_alignment = true;
                i -= 1;
                j -= 1;
            }
            Direction::GapInQuery => {
                let pos = &region_positions[j - 1];
                if pos.counts_for_confidence {
                    confidence_score += pos.gap_penalty;
                    max_confidence_score += pos.max_score;
                }
                j -= 1;
            }
            Direction::GapInConsensus => {
                positions.push(AlignedPosition::Insertion());
                i -= 1;
            }
        }
    }

    // Remaining query residues (if any) are the leading window prefix — kept out
    // of `positions` since RegionAlignment only covers query_start..=query_end.
    let query_start_in_window = i;

    positions.reverse();

    (
        positions,
        confidence_score,
        max_confidence_score,
        query_start_in_window,
        has_alignment,
    )
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

    fn test_align(query: &str, positions: &[PositionScores]) -> Alignment {
        align(query, positions, None)
    }
    #[test]
    fn test_align_simple() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

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
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let result = test_align("", &matrix.positions);
        // Empty sequence should align with negative infinity score and no positions
        assert_eq!(result.score, f32::NEG_INFINITY);
        assert!(result.positions.is_empty());
    }

    #[test]
    fn test_partial_sequence_start() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

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
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

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
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();

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
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
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
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
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

    #[test]
    fn test_leading_prefix() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
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

    fn all_free() -> RegionAlignParams {
        RegionAlignParams {
            free_query_prefix: true,
            free_query_suffix: true,
            free_cons_prefix: true,
            free_cons_suffix: true,
            suffix_clip_threshold: 50.0,
        }
    }

    #[test]
    fn test_align_region_full_matches_align() {
        // When called over the entire consensus with every corner free, align_region
        // should reach the same score as align() — the recurrence is identical.
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let mut buf = AlignBuffer::new();

        let full = align(FULL_IGH, &matrix.positions, Some(&mut buf));
        let region = align_region(
            FULL_IGH.as_bytes(),
            0,
            FULL_IGH.len(),
            &matrix.positions,
            all_free(),
            &mut buf,
        );

        assert!(region.has_alignment);
        assert!(
            (region.score - full.score).abs() < 1e-3,
            "scores differ: region={} full={}",
            region.score,
            full.score
        );
        assert_eq!(region.query_start, full.query_start);
        assert_eq!(region.query_end, full.query_end);
    }

    #[test]
    fn test_align_region_fr1_only() {
        // FR1 is IMGT 1..=26, so positions[0..26].
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let mut buf = AlignBuffer::new();
        let fr1_slice = &matrix.positions[0..26];

        let params = RegionAlignParams {
            free_query_prefix: true,
            free_query_suffix: true,
            free_cons_prefix: true,
            free_cons_suffix: false,
            suffix_clip_threshold: 0.0,
        };

        let result = align_region(
            FULL_IGH.as_bytes(),
            0,
            (26 + 16).min(FULL_IGH.len()),
            fr1_slice,
            params,
            &mut buf,
        );

        assert!(result.has_alignment);
        // FR1 of FULL_IGH should align starting near IMGT pos 1 and reach pos 26.
        assert!(
            result.score > 0.0,
            "FR1 should score positive on a real IGH"
        );
        // First aligned position should be Aligned(1..=N) and last should be Aligned(26).
        let last_aligned = result
            .positions
            .iter()
            .rev()
            .find_map(|p| match p {
                AlignedPosition::Aligned(n) => Some(*n),
                _ => None,
            })
            .expect("expected at least one aligned position");
        assert_eq!(last_aligned, 26, "FR1 should reach IMGT pos 26");
    }

    #[test]
    fn test_partial_alignment_full_igh() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let mut buf = AlignBuffer::new();
        let fr1 = align_fr1(FULL_IGH, &matrix, &mut buf);
        let result = align_remaining_frs(FULL_IGH, &matrix, &fr1, &mut buf).expect("should align");

        assert_eq!(result.query_start, 0, "no prefix expected");
        assert_eq!(result.query_end, FULL_IGH.len() - 1, "no suffix expected");
        assert_eq!(result.positions.len(), FULL_IGH.len());
        assert_eq!(result.cons_start, 1);
        assert_eq!(result.cons_end, 128);
        assert!(result.score > 0.0);
        let conf = result.confidence_score / result.max_confidence_score;
        assert!(
            conf > 0.5,
            "confidence should be high on a real IGH: got {}",
            conf
        );
    }

    #[test]
    fn test_partial_alignment_with_prefix_and_suffix() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let mut buf = AlignBuffer::new();
        let prefix = "AAAAAA";
        let suffix = "AAAAAAA";
        let sequence = format!("{prefix}{FULL_IGH}{suffix}");
        let fr1 = align_fr1(&sequence, &matrix, &mut buf);
        let result =
            align_remaining_frs(&sequence, &matrix, &fr1, &mut buf).expect("should align");

        assert_eq!(result.query_start, prefix.len());
        assert_eq!(result.query_end, prefix.len() + FULL_IGH.len() - 1);
        assert_eq!(result.positions.len(), sequence.len());
    }

    #[test]
    fn test_partial_alignment_rejects_non_antibody() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let mut buf = AlignBuffer::new();
        // Random non-antibody sequence: FR1 should score below FR1_CONFIDENCE_THRESHOLD,
        // so align_remaining_frs returns None without aligning FR2-FR4.
        let junk = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let fr1 = align_fr1(junk, &matrix, &mut buf);
        let result = align_remaining_frs(junk, &matrix, &fr1, &mut buf);
        assert!(
            result.is_none(),
            "non-antibody junk should fail FR1 threshold"
        );
    }

    #[test]
    fn test_align_region_empty_window() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let mut buf = AlignBuffer::new();
        let result = align_region(
            FULL_IGH.as_bytes(),
            10,
            10,
            &matrix.positions[0..26],
            all_free(),
            &mut buf,
        );
        assert!(!result.has_alignment);
        assert!(result.positions.is_empty());
    }

    #[test]
    fn test_both_flanking() {
        let matrix = ScoringMatrix::load(Chain::IGH).unwrap();
        let prefix = "AAAAAAA";
        let suffix = "AAAAAAA";
        let sequence = format!("{prefix}{FULL_IGH}{suffix}");
        let result = test_align(&sequence, &matrix.positions);

        assert_eq!(result.query_start, prefix.len());
        assert_eq!(result.query_end, prefix.len() + FULL_IGH.len() - 1);
        assert_eq!(result.positions.len(), sequence.len());
    }
}
