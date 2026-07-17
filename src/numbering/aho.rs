//! AHo numbering via direct position mapping from IMGT alignment
//!
//! Derived from the ANARCI `number_aho` scheme (see `anarci_source_files/schemes.py`)
//! and validated against the consensus fixtures in `fixtures/validation/ab_*_aho.csv`.
//!
//! AHo spans positions 1-149 and, unlike Kabat/Chothia/Martin, gaps its CDRs
//! symmetrically around a centre with chain-specific deletion orders:
//! - FR1 indel on 8
//! - CDR1 insertions on 36 (chain-specific deletion order)
//! - CDR2 insertions on 63
//! - FR3 indel/insertions on 85
//! - CDR3 insertions on 123 (symmetric deletion around 123/124)
//!
//! Per the project decision, each chain has its own rule table rather than extending
//! `NumberingRule` with a per-chain deletion order.

use crate::types::NumberingRule;
use crate::Insertion;

/// CDR3 deletion order shared by all chains (ANARCI `number_aho` CDR3 ordered_deletions),
/// expressed in AHo numbering. Deletions grow symmetrically outward from 123/124.
const AHO_CDR3_DELETIONS: &[u8] = &[
    123, 124, 122, 125, 121, 126, 120, 127, 119, 128, 118, 129, 117, 130, 116, 131, 115, 132, 114,
    133, 113, 134, 112, 135, 111, 136, 110, 137, 109, 138, 108, 107,
];

// =============================================================================
// AHo Heavy Chain Numbering Rules
// =============================================================================

/// AHo heavy chain numbering rules
///
/// Boundaries derived empirically from the paired IMGT/AHo consensus fixtures.
pub const AHO_HEAVY_RULES: &[NumberingRule] = &[
    // FR1: IMGT 1-8 -> AHo 1-9, indel on 8 (position 8 deleted first, then from front).
    NumberingRule::variable(1, 8, 1, 9, &[8, 1, 2, 3, 4, 5, 6, 7, 9], Insertion::Sequential(8)),
    // FR1 continued: IMGT 9-26 -> AHo 10-... (9 -> 10 offset +1, then 11 -> 11 offset 0)
    NumberingRule::offset(9, 9, 1),
    NumberingRule::offset(11, 26, 0),
    // CDR1: IMGT 27-36 -> AHo 27-38, insertions on 36 (VH deletion order).
    NumberingRule::variable(
        27,
        36,
        27,
        38,
        &[28, 36, 35, 37, 34, 38, 27, 33, 32, 29, 30, 31],
        Insertion::Sequential(36),
    ),
    // FR2: IMGT 37-58 -> AHo 39-60 (offset +2)
    NumberingRule::offset(37, 58, 2),
    // CDR2: IMGT 59-62 -> AHo 61-65, insertions on 63.
    NumberingRule::variable(59, 62, 61, 65, &[63, 62, 64, 61, 65], Insertion::Sequential(63)),
    // FR3: IMGT 63-72 -> AHo 66-75 (offset +3), then 74-104 -> 76-106 (offset +2)
    NumberingRule::offset(63, 72, 3),
    NumberingRule::offset(74, 104, 2),
    // CDR3: IMGT 105-117 -> AHo 107-138, insertions on 123, symmetric deletion.
    NumberingRule::variable(105, 117, 107, 138, AHO_CDR3_DELETIONS, Insertion::Sequential(123)),
    // FR4: IMGT 118-128 -> AHo 139-149 (offset +21)
    NumberingRule::offset(118, 128, 21),
];
