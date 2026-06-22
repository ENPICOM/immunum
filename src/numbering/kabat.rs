//! Kabat-specific numbering via direct position mapping from IMGT alignment
//!
//! Kabat uses sequential insertion patterns:
//! - Insertions: all at a single position (35A, 35B, 35C, ...)
//! - Deletions: from end or a custom order

use crate::types::{NumberingRule, Region};
use crate::Insertion;

// =============================================================================
// Kabat Region Boundaries (numbering space, shared by heavy and light chains)
// =============================================================================

/// Maps each functional region to its inclusive numbering-position range under
/// Kabat. Ranges are contiguous and cover the full numbered range 1..=113.
///
/// Note: these region boundaries are independent of the rule `num_start`/
/// `num_end` boundaries — a single rule's numbering range may straddle two
/// regions (e.g. the heavy `var(24, 40, 23, 35, ...)` rule spans FR1 and CDR1).
pub const KABAT_REGIONS: &[(Region, u8, u8)] = &[
    (Region::FR1, 1, 25),
    (Region::CDR1, 26, 35),
    (Region::FR2, 36, 50),
    (Region::CDR2, 51, 57),
    (Region::FR3, 58, 92),
    (Region::CDR3, 93, 100),
    (Region::FR4, 101, 113),
];

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbering::region_for_position;
    use crate::Scheme;

    #[test]
    fn test_regions_are_contiguous() {
        let mut expected = KABAT_REGIONS[0].1;
        for (region, start, end) in KABAT_REGIONS {
            assert_eq!(
                *start, expected,
                "gap or overlap before region {region:?}: expected start {expected}, got {start}"
            );
            assert!(start <= end, "region {region:?} has inverted range");
            expected = end + 1;
        }
    }

    #[test]
    fn test_rule_numbering_ranges_have_no_gaps() {
        for (label, rules) in [("heavy", KABAT_HEAVY_RULES), ("light", KABAT_LIGHT_RULES)] {
            let mut ranges: Vec<(u8, u8)> =
                rules.iter().map(|r| (r.num_start, r.num_end)).collect();
            ranges.sort();
            let mut expected = ranges[0].0;
            for (start, end) in ranges {
                assert_eq!(
                    start, expected,
                    "gap in {label} rule numbering ranges at {expected}"
                );
                expected = end + 1;
            }
        }
    }

    #[test]
    fn test_every_rule_position_has_a_region() {
        for (label, rules) in [("heavy", KABAT_HEAVY_RULES), ("light", KABAT_LIGHT_RULES)] {
            for rule in rules {
                for pos in rule.num_start..=rule.num_end {
                    assert!(
                        region_for_position(pos, Scheme::Kabat).is_some(),
                        "{label} numbering position {pos} maps to no Kabat region"
                    );
                }
            }
        }
    }
}
