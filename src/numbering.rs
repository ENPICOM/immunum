//! Numbering schemes for antibody and TCR sequences
//!
//! This module provides position numbering for different schemes (IMGT, Kabat, etc.).
//! Each scheme has its own submodule with scheme-specific configurations.
//!
//! The core abstraction is `RenumberConfig` which can express both:
//! - Palindromic patterns (IMGT): deletions via custom order, insertions split between two positions
//! - Sequential patterns (Kabat): deletions from end or custom order, insertions at single position

pub mod imgt;
pub mod kabat;

use crate::types::{InsertionStyle, Position, RenumberConfig};

/// Renumber positions based on configuration
///
/// This is the core renumbering function used by both IMGT and Kabat schemes.
/// It handles both deletions (when len < base) and insertions (when len > base).
pub fn renumber(positions: &[Position], config: &RenumberConfig) -> Vec<Position> {
    let len = positions.len();
    let base_len = config.base_positions.len();

    if len == 0 {
        return Vec::new();
    }

    if len <= base_len {
        renumber_deletions(len, config)
    } else {
        renumber_insertions(len, config)
    }
}

/// Handle deletions: select which base positions to use
fn renumber_deletions(len: usize, config: &RenumberConfig) -> Vec<Position> {
    let base_len = config.base_positions.len();

    let positions_to_use: Vec<u8> = match config.deletion_order {
        None => {
            // Default: use first `len` positions (delete from end)
            config.base_positions[..len].to_vec()
        }
        Some(order) => {
            // Remove positions in specified order
            let mut available: Vec<u8> = config.base_positions.to_vec();
            let to_remove = base_len - len;
            for &pos in order.iter().take(to_remove) {
                available.retain(|&p| p != pos);
            }
            available.sort_unstable();
            available
        }
    };

    positions_to_use.into_iter().map(Position::new).collect()
}

/// Handle insertions: all base positions plus insertion letters
fn renumber_insertions(len: usize, config: &RenumberConfig) -> Vec<Position> {
    let base_len = config.base_positions.len();
    let extra = len - base_len;
    let mut result = Vec::with_capacity(len);

    match config.insertion_style {
        InsertionStyle::AfterPosition(insertion_pos) => {
            // Find index of insertion position in base_positions
            let insertion_idx = config
                .base_positions
                .iter()
                .position(|&p| p == insertion_pos)
                .expect("Insertion position must be in base positions");

            // Add positions up to and including insertion_pos
            for &pos in &config.base_positions[..=insertion_idx] {
                result.push(Position::new(pos));
            }

            // Add insertions: A, B, C, etc.
            for i in 0..extra {
                result.push(Position::with_insertion(
                    insertion_pos,
                    (b'A' + i as u8) as char,
                ));
            }

            // Add remaining base positions
            for &pos in &config.base_positions[insertion_idx + 1..] {
                result.push(Position::new(pos));
            }
        }
        InsertionStyle::Palindromic { left, right } => {
            // Find indices for left and right insertion positions
            let left_idx = config
                .base_positions
                .iter()
                .position(|&p| p == left)
                .expect("Left insertion position must be in base positions");
            let right_idx = config
                .base_positions
                .iter()
                .position(|&p| p == right)
                .expect("Right insertion position must be in base positions");

            // Split insertions: first goes to right, then alternates
            let insertions_left = extra / 2;
            let insertions_right = extra.div_ceil(2);

            // Add positions up to and including left insertion point
            for &pos in &config.base_positions[..=left_idx] {
                result.push(Position::new(pos));
            }

            // Add left insertions (forward order: A, B, C...)
            for i in 0..insertions_left {
                result.push(Position::with_insertion(left, (b'A' + i as u8) as char));
            }

            // Add right insertions (reverse order: B, A for 2 insertions)
            for i in (0..insertions_right).rev() {
                result.push(Position::with_insertion(right, (b'A' + i as u8) as char));
            }

            // Add positions from right insertion point to end
            for &pos in &config.base_positions[right_idx..] {
                result.push(Position::new(pos));
            }
        }
    }

    result
}

/// Simple offset-based renumbering for framework regions
pub fn renumber_offset(positions: &[Position], src_start: u8, dst_start: u8) -> Vec<Position> {
    positions
        .iter()
        .map(|p| Position {
            number: p.number.wrapping_sub(src_start).wrapping_add(dst_start),
            insertion: p.insertion,
        })
        .collect()
}
