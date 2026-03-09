//! Kabat-specific numbering via direct position mapping from IMGT alignment
//!
//! Kabat uses sequential insertion patterns:
//! - Insertions: all at a single position (35A, 35B, 35C, ...)
//! - Deletions: from end or a custom order

use crate::types::NumberingRule;
use crate::Insertion;

// =============================================================================
// Kabat Heavy Chain Numbering Rules
// =============================================================================

/// Kabat heavy chain numbering rules
pub const KABAT_HEAVY_RULES: &[NumberingRule] = &[
    NumberingRule::fr(1, 9),
    NumberingRule::offset(11, 23, -1),
    NumberingRule::variable(
        24,
        40,
        23,
        35,
        &[35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25],
        Insertion::Sequential(35),
    ),
    NumberingRule::offset(41, 53, -5),
    NumberingRule::variable(
        54,
        65,
        49,
        57,
        &[52, 51, 50, 53, 54, 55, 56, 57, 49],
        Insertion::Sequential(52),
    ),
    NumberingRule::offset(66, 72, -8),
    NumberingRule::offset(74, 90, -9),
    NumberingRule::variable(91, 94, 82, 82, &[], Insertion::Sequential(82)),
    NumberingRule::offset(95, 105, -12),
    NumberingRule::variable(
        106,
        117,
        94,
        102,
        &[100, 99, 98, 97, 96, 95, 94],
        Insertion::Sequential(100),
    ),
    NumberingRule::offset(118, 128, -15),
];

// =============================================================================
// Kabat Light Chain Numbering Rules
// =============================================================================

/// Kabat light chain numbering rules (shared by kappa and lambda)
pub const KABAT_LIGHT_RULES: &[NumberingRule] = &[
    NumberingRule::fr(1, 25),
    NumberingRule::variable(
        26,
        38,
        26,
        32,
        &[28, 29, 30, 31, 27, 26],
        Insertion::Sequential(27),
    ),
    NumberingRule::offset(39, 55, -6),
    NumberingRule::variable(56, 65, 50, 52, &[52], Insertion::Sequential(52)),
    NumberingRule::offset(66, 72, -13),
    NumberingRule::offset(74, 79, -14),
    NumberingRule::variable(80, 82, 66, 66, &[], Insertion::Sequential(66)),
    NumberingRule::offset(83, 104, -16),
    NumberingRule::variable(
        105,
        116,
        89,
        96,
        &[95, 94, 93, 92, 91, 90, 89],
        Insertion::Sequential(95),
    ),
    NumberingRule::offset(117, 127, -20),
];
