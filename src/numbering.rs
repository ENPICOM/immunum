pub mod imgt;
pub mod kabat;

use crate::alignment::AlignedPosition;
use crate::types::{
    InsertionStyle, NumberingRegion, NumberingRegionType, Position, NumberingConfig,
};

/// Generate positions for a CDR region based on its length and config
///
/// Handles both deletions (len < base) and insertions (len > base).
pub fn number_with_config(len: usize, config: &NumberingConfig) -> Vec<Position> {
    if len == 0 {
        return Vec::new();
    }

    let base_len = config.base_positions.len();
    let mut result = Vec::with_capacity(len);

    if len <= base_len {
        // Deletions: select which base positions to keep
        let to_remove = base_len - len;

        match config.deletion_order {
            None => {
                // Simple case: keep first `len` positions
                result.extend(config.base_positions[..len].iter().map(|&p| Position::new(p)));
            }
            Some(order) => {
                // Build a mask of positions to skip
                let skip_set: &[u8] = &order[..to_remove];
                for &pos in config.base_positions {
                    if !skip_set.contains(&pos) {
                        result.push(Position::new(pos));
                    }
                }
            }
        }
    } else {
        // Insertions: all base positions plus extra with letters
        let extra = len - base_len;

        match config.insertion_style {
            InsertionStyle::AfterPosition(insertion_pos) => {
                for &pos in config.base_positions {
                    result.push(Position::new(pos));
                    if pos == insertion_pos {
                        // Add insertions after this position
                        for i in 0..extra {
                            result.push(Position::with_insertion(pos, (b'A' + i as u8) as char));
                        }
                    }
                }
            }
            InsertionStyle::Palindromic { left, right } => {
                let insertions_left = extra / 2;
                let insertions_right = extra.div_ceil(2);

                for &pos in config.base_positions {
                    if pos == left {
                        result.push(Position::new(pos));
                        // Left insertions: A, B, C...
                        for i in 0..insertions_left {
                            result.push(Position::with_insertion(left, (b'A' + i as u8) as char));
                        }
                        // Right insertions: ...B, A (reverse order)
                        for i in (0..insertions_right).rev() {
                            result.push(Position::with_insertion(right, (b'A' + i as u8) as char));
                        }
                    } else if pos == right {
                        // Right base position comes after insertions
                        result.push(Position::new(pos));
                    } else {
                        result.push(Position::new(pos));
                    }
                }
            }
        }
    }

    result
}

// =============================================================================
// Apply Numbering Function
// =============================================================================

/// Apply numbering directly to alignment results using numbering regions
///
/// Processes alignment in two phases:
/// 1. Extract consensus positions from alignment (skipping gaps)
/// 2. Apply region-based numbering scheme
pub fn apply_numbering(
    aligned_positions: &[AlignedPosition],
    regions: &[NumberingRegion],
) -> Vec<Position> {
    let consensus_positions = extract_consensus_positions(aligned_positions);
    number_by_regions(&consensus_positions, regions)
}

/// Extract consensus position for each residue from alignment
///
/// Returns a position for each residue (non-gap) in the alignment.
/// Insertions inherit the position of the preceding aligned residue.
fn extract_consensus_positions(aligned: &[AlignedPosition]) -> Vec<u8> {
    let mut positions = Vec::with_capacity(aligned.len());
    let mut last_pos: u8 = 0;

    for ap in aligned {
        match ap {
            AlignedPosition::QueryGap => {} // No residue, skip
            AlignedPosition::Aligned(pos) => {
                positions.push(*pos);
                last_pos = *pos;
            }
            AlignedPosition::Insertion => {
                positions.push(last_pos); // Inherit from previous
            }
        }
    }

    positions
}

/// Apply region-based numbering to consensus positions
fn number_by_regions(consensus_positions: &[u8], regions: &[NumberingRegion]) -> Vec<Position> {
    let mut result = Vec::with_capacity(consensus_positions.len());
    let mut idx = 0;

    for region in regions {
        let region_start = idx;

        // Find all positions belonging to this region
        while idx < consensus_positions.len() && region.contains(consensus_positions[idx]) {
            idx += 1;
        }

        let region_len = idx - region_start;
        if region_len == 0 {
            continue;
        }

        let region_positions = &consensus_positions[region_start..idx];

        match &region.region_type {
            NumberingRegionType::Offset {
                src_start,
                dst_start,
            } => {
                number_with_offset(region_positions, *src_start, *dst_start, &mut result);
            }
            NumberingRegionType::WithConfig(numbering_config) => {
                let numbered = number_with_config(region_len, numbering_config);
                result.extend(numbered);
            }
        }
    }

    result
}

/// Apply offset numbering to a region, handling insertions
fn number_with_offset(positions: &[u8], src_start: u8, dst_start: u8, result: &mut Vec<Position>) {
    let mut last_src_pos: Option<u8> = None;
    let mut insertion_count = 0u8;

    for &pos in positions {
        if Some(pos) == last_src_pos {
            // Same position as previous = insertion
            insertion_count += 1;
            let dst_pos = pos.wrapping_sub(src_start).wrapping_add(dst_start);
            result.push(Position::with_insertion(
                dst_pos,
                (b'A' + insertion_count - 1) as char,
            ));
        } else {
            // New position
            let dst_pos = pos.wrapping_sub(src_start).wrapping_add(dst_start);
            result.push(Position::new(dst_pos));
            last_src_pos = Some(pos);
            insertion_count = 0;
        }
    }
}
