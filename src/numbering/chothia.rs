//! Chothia numbering via direct position mapping from IMGT alignment
//!
//! Derived from the ANARCI `number_chothia_heavy` / `number_chothia_light` schemes
//! (see `anarci_source_files/schemes.py`) and validated against the consensus
//! fixtures in `fixtures/validation/ab_*_chothia.csv`.
//!
//! Chothia uses sequential insertion patterns:
//! - CDR1 insertions on 31 (heavy) / 30 (light)
//! - CDR2 insertions on 52
//! - FR3 insertions on 82 (heavy) / 68 (light)
//! - CDR3 insertions on 100 (heavy) / 95 (light)

use crate::types::{NumberingRule, RegionDefinition};
use crate::Insertion;

/// Chothia region definition
/// Chothia CDR spans: the structural loops rather than Kabat's sequence-variability windows.
///
/// Source: table of CDR definitions, Martin group, UCL -- <http://bioinf.org.uk/abs/info.html>
///
/// | loop | Chothia definition |
/// |------|--------------------|
/// | H1   | H26-H32            |
/// | H2   | H52-H56            |
/// | H3   | H96-H101           |
/// | L1   | L26-L32            |
/// | L2   | L50-L52            |
/// | L3   | L91-L96            |
///
/// These are the *consensus* values that page adopted on 28 June 2021, after surveying the Chothia
/// papers -- which disagree with each other, as its per-paper table shows. It deliberately narrowed
/// the earlier, wider set (H3 95-102; light 24-34 / 50-56 / 89-97). AbNumber and SADIE still ship
/// that wider set, where Chothia light happens to coincide with Kabat light, so Chothia segments
/// from this crate are tighter than those tools report.
pub const CHOTHIA_HEAVY_REGIONS: RegionDefinition = RegionDefinition {
    fr1_end: 25,
    cdr1_end: 32,
    fr2_end: 51,
    cdr2_end: 56,
    fr3_end: 95,
    cdr3_end: 101,
    fr4_end: 113,
};

pub const CHOTHIA_LIGHT_REGIONS: RegionDefinition = RegionDefinition {
    fr1_end: 25,
    cdr1_end: 32,
    fr2_end: 49,
    cdr2_end: 52,
    fr3_end: 90,
    cdr3_end: 96,
    fr4_end: 107,
};

// =============================================================================
// Chothia Heavy Chain Numbering Rules
// =============================================================================

/// Chothia heavy chain numbering rules
///
/// Boundaries derived empirically from the paired IMGT/Chothia consensus fixtures.
pub const CHOTHIA_HEAVY_RULES: &[NumberingRule] = &[
    // FR1: IMGT 1-9 map 1:1 to Chothia 1-9
    NumberingRule::fr(1, 9),
    // FR1 continued: IMGT 11-26 -> Chothia 10-25 (offset -1)
    NumberingRule::offset(11, 26, -1),
    // CDR1: IMGT 27-38 -> Chothia 26-33, insertions on 31.
    NumberingRule::variable(
        27,
        38,
        26,
        33,
        &[31, 30, 29, 28, 27, 26, 32, 33],
        Insertion::Sequential(31),
    ),
    // FR2: IMGT 39-53 -> Chothia 34-48 (offset -5)
    NumberingRule::offset(39, 53, -5),
    // CDR2: IMGT 54-66 -> Chothia 49-58, insertions on 52.
    NumberingRule::variable(
        54,
        66,
        49,
        58,
        &[52, 51, 50, 53, 54, 55, 56, 57, 58, 49],
        Insertion::Sequential(52),
    ),
    // FR3: IMGT 67-88 -> Chothia 59-79 (offset -8 then -9 across the 73 gap)
    NumberingRule::offset(67, 72, -8),
    NumberingRule::offset(74, 88, -9),
    // FR3 insertions on 82: IMGT 89-95 -> Chothia 80-83, insertions on 82.
    NumberingRule::variable(89, 95, 80, 83, &[83, 82, 81, 80], Insertion::Sequential(82)),
    // FR3 continued: IMGT 96-103 -> Chothia 84-91 (offset -12)
    NumberingRule::offset(96, 103, -12),
    // CDR3: IMGT 104-117 -> Chothia 92-102, insertions on 100.
    NumberingRule::variable(
        104,
        117,
        92,
        102,
        &[100, 99, 98, 97, 96, 95, 101, 102, 94, 93, 92],
        Insertion::Sequential(100),
    ),
    // FR4: IMGT 118-128 -> Chothia 103-113 (offset -15)
    NumberingRule::offset(118, 128, -15),
];

// =============================================================================
// Chothia Light Chain Numbering Rules
// =============================================================================

/// Chothia light chain numbering rules (shared by kappa and lambda; also used by Martin light)
///
/// Boundaries derived empirically from the paired IMGT/Chothia consensus fixtures.
pub const CHOTHIA_LIGHT_RULES: &[NumberingRule] = &[
    // FR1: IMGT 1-26 map 1:1 to Chothia 1-26
    NumberingRule::fr(1, 26),
    // CDR1: IMGT 27-40 -> Chothia 27-34, insertions on 30.
    // Deletion order follows ANARCI (the Chothia fixtures); note AntPack-generated
    // Martin light gaps a few short CDR1s the other way, which we accept as a
    // tool-level fixture difference.
    NumberingRule::variable(
        27,
        40,
        27,
        34,
        &[31, 32, 33, 30, 29, 28, 27, 34],
        Insertion::Sequential(30),
    ),
    // FR2: IMGT 41-55 -> Chothia 35-49 (offset -6)
    NumberingRule::offset(41, 55, -6),
    // CDR2: IMGT 56-65 -> Chothia 50-52, insertions on 52.
    NumberingRule::variable(56, 65, 50, 52, &[52, 51, 50], Insertion::Sequential(52)),
    // FR3: IMGT 66-72 -> Chothia 53-59 (offset -13), then 74-80 -> 60-66 (offset -14)
    NumberingRule::offset(66, 72, -13),
    NumberingRule::offset(74, 80, -14),
    // FR3 insertions on 68: IMGT 81-84 -> Chothia 67-68, insertions on 68.
    NumberingRule::variable(81, 84, 67, 68, &[68, 67], Insertion::Sequential(68)),
    // FR3 continued: IMGT 85-104 -> Chothia 69-88 (offset -16)
    NumberingRule::offset(85, 104, -16),
    // CDR3: IMGT 105-117 -> Chothia 89-97, insertions on 95.
    NumberingRule::variable(
        105,
        117,
        89,
        97,
        &[95, 94, 93, 92, 91, 90, 96, 97, 89],
        Insertion::Sequential(95),
    ),
    // FR4: IMGT 118-127 -> Chothia 98-107 (offset -20)
    NumberingRule::offset(118, 127, -20),
];
