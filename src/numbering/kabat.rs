//! Kabat-specific numbering via direct position mapping from IMGT alignment
//!
//! Kabat uses sequential patterns:
//! - Deletions: from end or custom order
//! - Insertions: all at single position (35A, 35B, 35C, ...)

use crate::alignment::{AlignedPosition, Alignment};
use crate::numbering::apply_numbering;
use crate::types::{Chain, NumberingRegion, Position, NumberingConfig};

// =============================================================================
// Kabat Heavy Chain Numbering Configuration
// =============================================================================

// CDR1: positions 23-35, insertions after 35
const HCDR1_CONFIG: NumberingConfig = NumberingConfig::sequential(
    &[23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
    35,
    None,
);

// CDR2: positions 49-57, insertions at 52, special deletion order
const HCDR2_CONFIG: NumberingConfig = NumberingConfig::sequential(
    &[49, 50, 51, 52, 53, 54, 55, 56, 57],
    52,
    Some(&[52, 51, 50, 53, 54, 55, 56, 57, 49]),
);

// FR3 special: positions 82, 82A, 82B, 82C
const HFR3_82_CONFIG: NumberingConfig = NumberingConfig::sequential(&[82], 82, None);

// CDR3: positions 94-102, insertions after 100
const HCDR3_CONFIG: NumberingConfig = NumberingConfig::sequential(
    &[94, 95, 96, 97, 98, 99, 100, 101, 102],
    100,
    Some(&[100, 99, 98, 97, 96, 95, 94]),
);

/// Kabat Heavy chain numbering regions
pub const KABAT_HEAVY_REGIONS: &[NumberingRegion] = &[
    NumberingRegion::offset(1, 9, 1),                     // FR1 part 1
    NumberingRegion::offset(11, 23, 10),                  // FR1 part 2
    NumberingRegion::with_config(24, 40, HCDR1_CONFIG),   // CDR1
    NumberingRegion::offset(41, 53, 36),                  // FR2
    NumberingRegion::with_config(54, 65, HCDR2_CONFIG),   // CDR2
    NumberingRegion::offset(66, 72, 58),                  // FR3 part 1
    NumberingRegion::offset(74, 90, 65),                  // FR3 part 2
    NumberingRegion::with_config(91, 94, HFR3_82_CONFIG), // FR3 82 region
    NumberingRegion::offset(95, 105, 83),                 // FR3 part 3
    NumberingRegion::with_config(106, 117, HCDR3_CONFIG), // CDR3
    NumberingRegion::offset(118, 128, 103),               // FR4
];

// =============================================================================
// Kabat Light Chain Numbering Configuration
// =============================================================================

// CDR1: positions 24-34, insertions at 27, delete forward from 28
const LCDR1_CONFIG: NumberingConfig = NumberingConfig::sequential(
    &[24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34],
    27,
    Some(&[28, 29, 30, 31, 32, 33, 34]),
);

// CDR2: positions 51-54, insertions at 52
// ANARCI comment: "How to gap L2 in Chothia/Kabat is unclear so we let the alignment do it"
const LCDR2_CONFIG: NumberingConfig = NumberingConfig::sequential(
    &[51, 52, 53, 54],
    52,
    None,
);

// CDR3: positions 89-95, insertions at 95, delete backward from 95
const LCDR3_CONFIG: NumberingConfig = NumberingConfig::sequential(
    &[89, 90, 91, 92, 93, 94, 95],
    95,
    Some(&[95, 94, 93, 92, 91, 90, 89]),
);

/// Kabat Light chain numbering regions
pub const KABAT_LIGHT_REGIONS: &[NumberingRegion] = &[
    NumberingRegion::offset(1, 23, 1),                     // FR1
    NumberingRegion::with_config(24, 40, LCDR1_CONFIG),    // CDR1
    NumberingRegion::offset(41, 56, 35),                   // FR2
    NumberingRegion::with_config(57, 67, LCDR2_CONFIG),    // CDR2
    NumberingRegion::offset(68, 72, 55),                   // FR3 part 1
    NumberingRegion::offset(74, 80, 60),                   // FR3 part 2 (skip IMGT 73)
    NumberingRegion::offset(83, 104, 67),                  // FR3 part 3 (skip IMGT 81-82)
    NumberingRegion::with_config(105, 117, LCDR3_CONFIG),  // CDR3
    NumberingRegion::offset(118, 128, 98),                 // FR4
];

/// Get Kabat-specific numbering: FR regions mapped from IMGT, CDR regions from Kabat rules
pub fn get_kabat_numbering(alignment: &Alignment, chain: Chain) -> Vec<Position> {
    match chain {
        Chain::IGH => apply_numbering(&alignment.positions, KABAT_HEAVY_REGIONS),
        Chain::IGK | Chain::IGL => apply_numbering(&alignment.positions, KABAT_LIGHT_REGIONS),
        _ => panic!("Kabat numbering only supported for antibody chains"),
    }
}

/// Get Kabat numbering from aligned positions directly
pub fn number_kabat(aligned_positions: &[AlignedPosition], chain: Chain) -> Vec<Position> {
    match chain {
        Chain::IGH => apply_numbering(aligned_positions, KABAT_HEAVY_REGIONS),
        Chain::IGK | Chain::IGL => apply_numbering(aligned_positions, KABAT_LIGHT_REGIONS),
        _ => panic!("Kabat numbering only supported for antibody chains"),
    }
}

