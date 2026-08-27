//! AHo numbering via direct position mapping from IMGT alignment
//!
//! Derived from the ANARCI `number_aho` scheme (see `anarci_source_files/schemes.py`,
//! `number_aho` + `_number_regions`) and validated against the consensus fixtures in
//! `fixtures/validation/ab_*_aho.csv`.
//!
//! AHo spans positions 1-149. Unlike Kabat/Chothia/Martin it gaps its CDRs symmetrically
//! around a centre, with chain-specific CDR1 deletion orders. The rule windows follow
//! ANARCI's region boundaries exactly (wide, framework-anchored on both sides) so that the
//! residue count fed to each CDR renumbering is stable regardless of interior gap placement:
//! - FR1 indel on 8
//! - CDR1 insertions on 36 (chain-specific deletion order: VH / VK / VL)
//! - CDR2 insertions on 63
//! - FR3 indel on 85 (deletions move onto 86 then 85)
//! - CDR3 insertions on 123 (symmetric deletion around 123/124)
//!
//! Light chains additionally carry one extra C-terminal position (AHo 149) that has no IMGT
//! state (light chains end at IMGT 127 -> AHo 148); that residue is appended in
//! `Annotator::number`, matching ANARCI's `number_aho` tail rule.
//!
//! Per the project decision, each chain has its own rule table rather than extending
//! `NumberingRule` with a per-chain deletion order.

use crate::types::{NumberingRule, RegionDefinition};
use crate::Insertion;

/// AHo regions, tiling all 149 positions.
/// Unverified, unlike the other schemes, and the sources conflict. Honegger & Plueckthun give
/// CDR-H1 27-42 / CDR-H2 57-76 and CDR-L1 24-42 / CDR-L2 58-72, both with CDR3 107-138. The CDR3
/// start below (109) matches neither, so check this table against the Honegger paper before relying
/// on AHo segmentation.
pub const AHO_REGIONS: RegionDefinition = RegionDefinition {
    fr1_end: 24,
    cdr1_end: 42,
    fr2_end: 56,
    cdr2_end: 77,
    fr3_end: 108,
    cdr3_end: 138,
    fr4_end: 149,
};

/// FR1 rule for heavy and lambda: IMGT 1-10 -> AHo 1-10 with the scheme indel on 8. Heavy and
/// lambda AHo reserve position 8 as an empty gap (IMGT 8 shifts to AHo 9), so the count-based
/// deletion order drops 8 first. Kappa does NOT gap 8 (its FR1 is a true 1:1), so it uses
/// `fr(1,10)` instead.
const AHO_FR1_INDEL8: NumberingRule = NumberingRule::variable(
    1,
    10,
    1,
    10,
    &[8, 1, 2, 3, 4, 5, 6, 7, 9, 10],
    Insertion::Sequential(8),
);

/// CDR2 rule shared by all antibody chains: IMGT 56-75 -> AHo 58-77, insertions on 63,
/// deletions symmetric around 63.
const AHO_CDR2: NumberingRule = NumberingRule::variable(
    56,
    75,
    58,
    77,
    &[
        63, 62, 64, 61, 65, 60, 66, 59, 67, 58, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77,
    ],
    Insertion::Sequential(63),
);

/// FR3 rule shared by all chains: IMGT 76-91 -> AHo 78-93, deletions onto 86 then 85,
/// insertions on 85.
const AHO_FR3: NumberingRule = NumberingRule::variable(
    76,
    91,
    78,
    93,
    &[
        86, 85, 87, 84, 88, 83, 89, 82, 90, 81, 91, 80, 92, 79, 93, 78,
    ],
    Insertion::Sequential(85),
);

/// CDR3 deletion order shared by all chains (ANARCI `number_aho` CDR3 ordered_deletions,
/// AHo numbering). Deletions grow symmetrically outward from 123/124.
const AHO_CDR3_DELETIONS: &[u8] = &[
    123, 124, 122, 125, 121, 126, 120, 127, 119, 128, 118, 129, 117, 130, 116, 131, 115, 132, 114,
    133, 113, 134, 112, 135, 111, 136, 110, 137, 109, 138, 108, 107,
];

/// CDR3 rule shared by all chains: IMGT 105-117 -> AHo 107-138, insertions on 123.
const AHO_CDR3: NumberingRule = NumberingRule::variable(
    105,
    117,
    107,
    138,
    AHO_CDR3_DELETIONS,
    Insertion::Sequential(123),
);

// =============================================================================
// AHo Heavy Chain Numbering Rules
// =============================================================================

/// AHo heavy chain numbering rules
pub const AHO_HEAVY_RULES: &[NumberingRule] = &[
    AHO_FR1_INDEL8,
    // FR1 cont. (region C): IMGT 11-24 -> AHo 11-24 (1:1).
    NumberingRule::fr(11, 24),
    // CDR1 (region D): IMGT 25-40 -> AHo 25-42, insertions on 36, VH deletion order.
    NumberingRule::variable(
        25,
        40,
        25,
        42,
        &[
            28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42,
        ],
        Insertion::Sequential(36),
    ),
    // FR2 (region E): IMGT 41-55 -> AHo 43-57 (offset +2).
    NumberingRule::offset(41, 55, 2),
    AHO_CDR2,
    AHO_FR3,
    // FR3 cont. (region I): IMGT 92-104 -> AHo 94-106 (offset +2).
    NumberingRule::offset(92, 104, 2),
    AHO_CDR3,
    // FR4 (region K): IMGT 118-128 -> AHo 139-149 (offset +21).
    NumberingRule::offset(118, 128, 21),
];

// =============================================================================
// AHo Kappa Chain Numbering Rules
// =============================================================================

/// AHo kappa chain numbering rules (VK CDR1 deletion order: two-residue gap at 27/28).
///
/// Kappa FR1 is a true 1:1 map (no gap at AHo 8), so it uses `fr` rather than the indel-on-8
/// rule used by heavy/lambda; this also keeps N-terminally truncated kappa sequences correct.
pub const AHO_KAPPA_RULES: &[NumberingRule] = &[
    NumberingRule::fr(1, 10),
    NumberingRule::fr(11, 24),
    // CDR1 (region D): IMGT 25-40 -> AHo 25-42, insertions on 36, VK deletion order.
    NumberingRule::variable(
        25,
        40,
        25,
        42,
        &[
            28, 27, 36, 35, 37, 34, 38, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42,
        ],
        Insertion::Sequential(36),
    ),
    NumberingRule::offset(41, 55, 2),
    AHO_CDR2,
    AHO_FR3,
    NumberingRule::offset(92, 104, 2),
    AHO_CDR3,
    NumberingRule::offset(118, 128, 21),
];

// =============================================================================
// AHo Lambda Chain Numbering Rules
// =============================================================================

/// AHo lambda chain numbering rules (VL CDR1 deletion order).
pub const AHO_LAMBDA_RULES: &[NumberingRule] = &[
    AHO_FR1_INDEL8,
    NumberingRule::fr(11, 24),
    // CDR1 (region D): IMGT 25-40 -> AHo 25-42, insertions on 36, VL deletion order.
    NumberingRule::variable(
        25,
        40,
        25,
        42,
        &[
            28, 36, 35, 37, 34, 38, 27, 29, 33, 39, 32, 40, 26, 30, 25, 31, 41, 42,
        ],
        Insertion::Sequential(36),
    ),
    NumberingRule::offset(41, 55, 2),
    AHO_CDR2,
    AHO_FR3,
    NumberingRule::offset(92, 104, 2),
    AHO_CDR3,
    NumberingRule::offset(118, 128, 21),
];
