//! Martin (extended Chothia) numbering via direct position mapping from IMGT alignment
//!
//! Derived from the ANARCI `number_martin_heavy` scheme (see
//! `anarci_source_files/schemes.py`) and validated against the consensus fixtures
//! in `fixtures/validation/ab_*_martin.csv`.
//!
//! Martin heavy is identical to Chothia heavy except that FR3 insertions are placed
//! explicitly on position 72 (rather than 82 as in Chothia).
//!
//! Martin light is defined identically to Chothia light (ANARCI implements
//! `number_martin_light` by delegating to `number_chothia_light`), so we reuse
//! `CHOTHIA_LIGHT_RULES`.

use crate::types::{NumberingRule, RegionDefinition};
use crate::Insertion;

/// Martin region definition (identical to Chothia's).
pub const MARTIN_REGIONS: RegionDefinition = RegionDefinition {
    fr1_end: 25,
    cdr1_end: 32,
    fr2_end: 51,
    cdr2_end: 56,
    fr3_end: 95,
    cdr3_end: 101,
    fr4_end: 113,
};

// =============================================================================
// Martin Heavy Chain Numbering Rules
// =============================================================================

/// Martin heavy chain numbering rules
///
/// Boundaries derived empirically from the paired IMGT/Martin consensus fixtures.
pub const MARTIN_HEAVY_RULES: &[NumberingRule] = &[
    // FR1: IMGT 1-9 map 1:1 to Martin 1-9
    NumberingRule::fr(1, 9),
    // FR1 continued: IMGT 11-26 -> Martin 10-25 (offset -1)
    NumberingRule::offset(11, 26, -1),
    // CDR1: IMGT 27-38 -> Martin 26-33, insertions on 31.
    NumberingRule::variable(
        27,
        38,
        26,
        33,
        &[31, 30, 29, 28, 27, 26, 32, 33],
        Insertion::Sequential(31),
    ),
    // FR2: IMGT 39-53 -> Martin 34-48 (offset -5)
    NumberingRule::offset(39, 53, -5),
    // CDR2: IMGT 54-66 -> Martin 49-58, insertions on 52.
    NumberingRule::variable(
        54,
        66,
        49,
        58,
        &[52, 51, 50, 53, 54, 55, 56, 57, 58, 49],
        Insertion::Sequential(52),
    ),
    // FR3: IMGT 67-80 -> Martin 59-71 (offset -8 then -9 across the 73 gap)
    NumberingRule::offset(67, 72, -8),
    NumberingRule::offset(74, 80, -9),
    // FR3 insertions on 72: IMGT 81-84 -> Martin 72, insertions on 72.
    NumberingRule::variable(81, 84, 72, 72, &[72], Insertion::Sequential(72)),
    // FR3 continued: IMGT 85-103 -> Martin 73-91 (offset -12)
    NumberingRule::offset(85, 103, -12),
    // CDR3: IMGT 104-117 -> Martin 92-102, insertions on 100.
    NumberingRule::variable(
        104,
        117,
        92,
        102,
        &[100, 99, 98, 97, 96, 95, 101, 102, 94, 93, 92],
        Insertion::Sequential(100),
    ),
    // FR4: IMGT 118-128 -> Martin 103-113 (offset -15)
    NumberingRule::offset(118, 128, -15),
];

// =============================================================================
// Martin Light Chain Numbering Rules
// =============================================================================

/// Martin light chain numbering rules
///
/// ANARCI defines Martin light identically to Chothia light, and every region here
/// matches `CHOTHIA_LIGHT_RULES` except the CDR1 deletion order: the AntPack-generated
/// Martin fixtures gap short CDR1 loops front-first (keeping 32-34), whereas the
/// ANARCI-generated Chothia fixtures gap the other way. We keep a distinct table so
/// each scheme reproduces its own reference exactly.
pub const MARTIN_LIGHT_RULES: &[NumberingRule] = &[
    // FR1: IMGT 1-26 map 1:1 to Martin 1-26
    NumberingRule::fr(1, 26),
    // CDR1: IMGT 27-40 -> Martin 27-34, insertions on 30 (front-first deletion order).
    NumberingRule::variable(
        27,
        40,
        27,
        34,
        &[31, 30, 29, 28, 27, 32, 33, 34],
        Insertion::Sequential(30),
    ),
    // FR2: IMGT 41-55 -> Martin 35-49 (offset -6)
    NumberingRule::offset(41, 55, -6),
    // CDR2: IMGT 56-65 -> Martin 50-52, insertions on 52.
    NumberingRule::variable(56, 65, 50, 52, &[52, 51, 50], Insertion::Sequential(52)),
    // FR3: IMGT 66-72 -> Martin 53-59 (offset -13), then 74-80 -> 60-66 (offset -14)
    NumberingRule::offset(66, 72, -13),
    NumberingRule::offset(74, 80, -14),
    // FR3 insertions on 68: IMGT 81-84 -> Martin 67-68, insertions on 68.
    NumberingRule::variable(81, 84, 67, 68, &[68, 67], Insertion::Sequential(68)),
    // FR3 continued: IMGT 85-104 -> Martin 69-88 (offset -16)
    NumberingRule::offset(85, 104, -16),
    // CDR3: IMGT 105-117 -> Martin 89-97, insertions on 95.
    NumberingRule::variable(
        105,
        117,
        89,
        97,
        &[95, 94, 93, 92, 91, 90, 96, 97, 89],
        Insertion::Sequential(95),
    ),
    // FR4: IMGT 118-127 -> Martin 98-107 (offset -20)
    NumberingRule::offset(118, 127, -20),
];
