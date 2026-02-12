//! Kabat-specific numbering via direct position mapping from IMGT alignment
//!
//! Kabat uses sequential patterns:
//! - Deletions: from end or custom order
//! - Insertions: all at single position (35A, 35B, 35C, ...)

use super::{renumber, renumber_offset};
use crate::types::{Position, RenumberConfig};

/// Region type for IMGT to Kabat conversion
#[derive(Debug, Clone, Copy)]
enum RegionType {
    /// Simple offset: imgt_pos - imgt_start + kabat_start
    Offset { imgt_start: u32, kabat_start: u32 },
    /// Renumbering with insertions/deletions
    WithInsertions(RenumberConfig),
}

/// Region definition for IMGT to Kabat conversion
#[derive(Debug, Clone, Copy)]
struct Region {
    /// IMGT position range (inclusive)
    imgt_range: (u32, u32),
    /// How to renumber this region
    region_type: RegionType,
}

impl Region {
    const fn offset(imgt_start: u32, imgt_end: u32, kabat_start: u32) -> Self {
        Self {
            imgt_range: (imgt_start, imgt_end),
            region_type: RegionType::Offset {
                imgt_start,
                kabat_start,
            },
        }
    }

    const fn with_insertions(imgt_start: u32, imgt_end: u32, config: RenumberConfig) -> Self {
        Self {
            imgt_range: (imgt_start, imgt_end),
            region_type: RegionType::WithInsertions(config),
        }
    }

    fn contains(&self, pos: u32) -> bool {
        pos >= self.imgt_range.0 && pos <= self.imgt_range.1
    }
}

// CDR1: positions 23-35, insertions after 35
const HCDR1_CONFIG: RenumberConfig = RenumberConfig::sequential(
    &[23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],
    35,
    None,
);

// CDR2: positions 49-57, insertions at 52, special deletion order
const HCDR2_CONFIG: RenumberConfig = RenumberConfig::sequential(
    &[49, 50, 51, 52, 53, 54, 55, 56, 57],
    52,
    Some(&[52, 51, 50, 53, 54, 55, 56, 57, 49]),
);

// FR3 special: positions 82, 82A, 82B, 82C
const HFR3_82_CONFIG: RenumberConfig = RenumberConfig::sequential(&[82], 82, None);

// CDR3: positions 94-102, insertions after 100
const HCDR3_CONFIG: RenumberConfig = RenumberConfig::sequential(
    &[94, 95, 96, 97, 98, 99, 100, 101, 102],
    100,
    Some(&[100, 99, 98, 97, 96, 95, 94]),
);

/// Region mapping for Heavy chains (IGH)
const HEAVY_REGIONS: &[Region] = &[
    Region::offset(1, 9, 1),
    Region::offset(11, 23, 10),
    Region::with_insertions(24, 40, HCDR1_CONFIG),
    Region::offset(41, 53, 36),
    Region::with_insertions(54, 65, HCDR2_CONFIG),
    Region::offset(66, 72, 58),
    Region::offset(74, 90, 65),
    Region::with_insertions(91, 94, HFR3_82_CONFIG),
    Region::offset(95, 105, 83),
    Region::with_insertions(106, 117, HCDR3_CONFIG),
    Region::offset(118, 128, 103),
];

/// Convert a sequence of IMGT positions to Kabat positions for Heavy chains (IGH)
pub fn imgt_to_kabat_heavy(imgt_positions: &[Position]) -> Vec<Position> {
    let mut result = Vec::with_capacity(imgt_positions.len());
    let mut idx = 0;

    for region in HEAVY_REGIONS {
        // Find positions in this region
        let start_idx = idx;
        while idx < imgt_positions.len() && region.contains(imgt_positions[idx].number) {
            idx += 1;
        }
        let region_positions = &imgt_positions[start_idx..idx];

        if region_positions.is_empty() {
            continue;
        }

        let converted = match &region.region_type {
            RegionType::Offset {
                imgt_start,
                kabat_start,
            } => renumber_offset(region_positions, *imgt_start, *kabat_start),
            RegionType::WithInsertions(config) => renumber(region_positions, config),
        };

        result.extend(converted);
    }

    result
}

// // CDR1: positions 23-35, insertions after 35
// const LCDR1_CONFIG: RenumberConfig = RenumberConfig::sequential(
//     &[24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34],
//     27,
//     None,
// );

// // CDR2: positions 49-57, insertions at 52, special deletion order
// const LCDR2_CONFIG: RenumberConfig = RenumberConfig::sequential(
//     &[50, 51, 52, 53, 54],
//     52,
//     Some(&[52, 51, 50, 53, 54, 55, 56, 57, 49]),
// );

// // FR3 special: positions 82, 82A, 82B, 82C
// const LFR3_82_CONFIG: RenumberConfig = RenumberConfig::sequential(&[82], 82, None);

// // CDR3: positions 94-102, insertions after 100
// const LCDR3_CONFIG: RenumberConfig = RenumberConfig::sequential(
//     &[94, 95, 96, 97, 98, 99, 100, 101, 102],
//     95,
//     Some(&[100, 99, 98, 97, 96, 95, 94]),
// );


// /// Region mapping for Light chains (IGK, IGL)
// const LIGHT_REGIONS: &[Region] = &[
//     Region::offset(1, 24, 1),
//     Region::with_insertions(24, 38, LCDR1_CONFIG),
//     Region::offset(41, 53, 36),
//     Region::with_insertions(54, 65, LCDR2_CONFIG),
//     Region::offset(66, 72, 58),
//     Region::offset(74, 90, 65),
//     Region::with_insertions(91, 94, LFR3_82_CONFIG),
//     Region::offset(95, 105, 83),
//     Region::with_insertions(106, 117, LCDR3_CONFIG),
//     Region::offset(118, 128, 103),
// ];

/// Convert a sequence of IMGT positions to Kabat positions for Light chains (IGK, IGL)
pub fn imgt_to_kabat_light(imgt_positions: &[Position]) -> Vec<Position> {
    imgt_positions.to_vec()
}
