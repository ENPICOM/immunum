use crate::types::{NumberingRule, Region};
use crate::Insertion;

// =============================================================================
// IMGT Region Boundaries (numbering space, same for all chain types)
// =============================================================================

/// Maps each functional region to its inclusive numbering-position range under
/// IMGT. Ranges are contiguous and cover the full numbered range 1..=128.
pub const IMGT_REGIONS: &[(Region, u8, u8)] = &[
    (Region::FR1, 1, 26),
    (Region::CDR1, 27, 38),
    (Region::FR2, 39, 55),
    (Region::CDR2, 56, 65),
    (Region::FR3, 66, 104),
    (Region::CDR3, 105, 117),
    (Region::FR4, 118, 128),
];

// =============================================================================
// IMGT Chain Numbering Rules (same for all chain types)
// =============================================================================

pub const IMGT_RULES: &[NumberingRule] = &[
    NumberingRule::fr(1, 26),
    NumberingRule::variable(
        27,
        38,
        27,
        38,
        &[33, 32, 34, 31, 35, 30, 36, 29, 37, 28, 38],
        Insertion::Symmetric {
            left: 32,
            right: 33,
        },
    ),
    NumberingRule::fr(39, 55),
    NumberingRule::variable(
        56,
        65,
        56,
        65,
        &[61, 60, 62, 59, 63, 58, 64, 57, 65],
        Insertion::Symmetric {
            left: 60,
            right: 61,
        },
    ),
    NumberingRule::fr(66, 104),
    NumberingRule::variable(
        105,
        117,
        105,
        117,
        &[111, 112, 110, 113, 109, 114, 108, 115, 107, 116, 106, 117],
        Insertion::Symmetric {
            left: 111,
            right: 112,
        },
    ),
    NumberingRule::fr(118, 128),
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbering::number_with_rules;

    #[test]
    fn test_regions_are_contiguous() {
        let mut expected = IMGT_REGIONS[0].1;
        for (region, start, end) in IMGT_REGIONS {
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
        let mut ranges: Vec<(u8, u8)> =
            IMGT_RULES.iter().map(|r| (r.num_start, r.num_end)).collect();
        ranges.sort();
        let mut expected = ranges[0].0;
        for (start, end) in ranges {
            assert_eq!(start, expected, "gap in rule numbering ranges at {expected}");
            expected = end + 1;
        }
    }

    #[test]
    fn test_every_rule_position_has_a_region() {
        use crate::numbering::region_for_position;
        use crate::Scheme;
        for rule in IMGT_RULES {
            for pos in rule.num_start..=rule.num_end {
                assert!(
                    region_for_position(pos, Scheme::IMGT).is_some(),
                    "numbering position {pos} maps to no IMGT region"
                );
            }
        }
    }

    #[test]
    fn test_cdr1_numbering() {
        let rule = &IMGT_RULES[1]; // CDR1
        for (length, output_strings) in vec![
            (2, vec!["27", "38"]),
            (3, vec!["27", "28", "38"]),
            (4, vec!["27", "28", "37", "38"]),
            (5, vec!["27", "28", "29", "37", "38"]),
            (6, vec!["27", "28", "29", "36", "37", "38"]),
            (7, vec!["27", "28", "29", "30", "36", "37", "38"]),
            (8, vec!["27", "28", "29", "30", "35", "36", "37", "38"]),
            (
                9,
                vec!["27", "28", "29", "30", "31", "35", "36", "37", "38"],
            ),
            (
                10,
                vec!["27", "28", "29", "30", "31", "34", "35", "36", "37", "38"],
            ),
            (
                11,
                vec![
                    "27", "28", "29", "30", "31", "32", "34", "35", "36", "37", "38",
                ],
            ),
            (
                12,
                vec![
                    "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38",
                ],
            ),
            // Symmetric insertions: first goes to right (33A), pattern: 32A, 33A, 33
            (
                13,
                vec![
                    "27", "28", "29", "30", "31", "32", "33A", "33", "34", "35", "36", "37", "38",
                ],
            ),
            (
                14,
                vec![
                    "27", "28", "29", "30", "31", "32", "32A", "33A", "33", "34", "35", "36", "37",
                    "38",
                ],
            ),
            (
                15,
                vec![
                    "27", "28", "29", "30", "31", "32", "32A", "33B", "33A", "33", "34", "35",
                    "36", "37", "38",
                ],
            ),
            (
                16,
                vec![
                    "27", "28", "29", "30", "31", "32", "32A", "32B", "33B", "33A", "33", "34",
                    "35", "36", "37", "38",
                ],
            ),
        ] {
            let positions = number_with_rules(length, rule);
            assert_eq!(positions.len(), length);
            let pos_strings: Vec<String> = positions.iter().map(|p| p.to_string()).collect();
            assert_eq!(pos_strings, output_strings, "Failed for length {}", length);
        }
    }

    #[test]
    fn test_cdr2_numbering() {
        let rule = &IMGT_RULES[3]; // CDR2
        for (length, output_strings) in vec![
            (2, vec!["56", "65"]),
            (3, vec!["56", "57", "65"]),
            (4, vec!["56", "57", "64", "65"]),
            (5, vec!["56", "57", "58", "64", "65"]),
            (6, vec!["56", "57", "58", "63", "64", "65"]),
            (7, vec!["56", "57", "58", "59", "63", "64", "65"]),
            (8, vec!["56", "57", "58", "59", "62", "63", "64", "65"]),
            (
                9,
                vec!["56", "57", "58", "59", "60", "62", "63", "64", "65"],
            ),
            (
                10,
                vec!["56", "57", "58", "59", "60", "61", "62", "63", "64", "65"],
            ),
            (
                11,
                vec![
                    "56", "57", "58", "59", "60", "61A", "61", "62", "63", "64", "65",
                ],
            ),
            (
                12,
                vec![
                    "56", "57", "58", "59", "60", "60A", "61A", "61", "62", "63", "64", "65",
                ],
            ),
            (
                13,
                vec![
                    "56", "57", "58", "59", "60", "60A", "61B", "61A", "61", "62", "63", "64", "65",
                ],
            ),
            (
                14,
                vec![
                    "56", "57", "58", "59", "60", "60A", "60B", "61B", "61A", "61", "62", "63",
                    "64", "65",
                ],
            ),
            (
                15,
                vec![
                    "56", "57", "58", "59", "60", "60A", "60B", "61C", "61B", "61A", "61", "62",
                    "63", "64", "65",
                ],
            ),
            (
                16,
                vec![
                    "56", "57", "58", "59", "60", "60A", "60B", "60C", "61C", "61B", "61A", "61",
                    "62", "63", "64", "65",
                ],
            ),
        ] {
            let positions = number_with_rules(length, rule);
            assert_eq!(positions.len(), length);
            let pos_strings: Vec<String> = positions.iter().map(|p| p.to_string()).collect();
            assert_eq!(pos_strings, output_strings, "Failed for length {}", length);
        }
    }

    #[test]
    fn test_cdr3_numbering() {
        let rule = &IMGT_RULES[5]; // CDR3
        for (length, output_strings) in vec![
            (2, vec!["105", "117"]),
            (3, vec!["105", "106", "117"]),
            (4, vec!["105", "106", "116", "117"]),
            (5, vec!["105", "106", "107", "116", "117"]),
            (6, vec!["105", "106", "107", "115", "116", "117"]),
            (7, vec!["105", "106", "107", "108", "115", "116", "117"]),
            (
                8,
                vec!["105", "106", "107", "108", "114", "115", "116", "117"],
            ),
            (
                9,
                vec![
                    "105", "106", "107", "108", "109", "114", "115", "116", "117",
                ],
            ),
            (
                10,
                vec![
                    "105", "106", "107", "108", "109", "113", "114", "115", "116", "117",
                ],
            ),
            (
                11,
                vec![
                    "105", "106", "107", "108", "109", "110", "113", "114", "115", "116", "117",
                ],
            ),
            (
                12,
                vec![
                    "105", "106", "107", "108", "109", "110", "112", "113", "114", "115", "116",
                    "117",
                ],
            ),
            (
                13,
                vec![
                    "105", "106", "107", "108", "109", "110", "111", "112", "113", "114", "115",
                    "116", "117",
                ],
            ),
            (
                14,
                vec![
                    "105", "106", "107", "108", "109", "110", "111", "112A", "112", "113", "114",
                    "115", "116", "117",
                ],
            ),
            (
                15,
                vec![
                    "105", "106", "107", "108", "109", "110", "111", "111A", "112A", "112", "113",
                    "114", "115", "116", "117",
                ],
            ),
            (
                16,
                vec![
                    "105", "106", "107", "108", "109", "110", "111", "111A", "112B", "112A", "112",
                    "113", "114", "115", "116", "117",
                ],
            ),
            (
                17,
                vec![
                    "105", "106", "107", "108", "109", "110", "111", "111A", "111B", "112B",
                    "112A", "112", "113", "114", "115", "116", "117",
                ],
            ),
            (
                18,
                vec![
                    "105", "106", "107", "108", "109", "110", "111", "111A", "111B", "112C",
                    "112B", "112A", "112", "113", "114", "115", "116", "117",
                ],
            ),
            (
                19,
                vec![
                    "105", "106", "107", "108", "109", "110", "111", "111A", "111B", "111C",
                    "112C", "112B", "112A", "112", "113", "114", "115", "116", "117",
                ],
            ),
            (
                20,
                vec![
                    "105", "106", "107", "108", "109", "110", "111", "111A", "111B", "111C",
                    "112D", "112C", "112B", "112A", "112", "113", "114", "115", "116", "117",
                ],
            ),
            (
                25,
                vec![
                    "105", "106", "107", "108", "109", "110", "111", "111A", "111B", "111C",
                    "111D", "111E", "111F", "112F", "112E", "112D", "112C", "112B", "112A", "112",
                    "113", "114", "115", "116", "117",
                ],
            ),
            (
                30,
                vec![
                    "105", "106", "107", "108", "109", "110", "111", "111A", "111B", "111C",
                    "111D", "111E", "111F", "111G", "111H", "112I", "112H", "112G", "112F", "112E",
                    "112D", "112C", "112B", "112A", "112", "113", "114", "115", "116", "117",
                ],
            ),
        ] {
            let positions = number_with_rules(length, rule);
            assert_eq!(positions.len(), length);
            let pos_strings: Vec<String> = positions.iter().map(|p| p.to_string()).collect();
            assert_eq!(pos_strings, output_strings, "Failed for length {}", length);
        }
    }
}
