use crate::alignment::{AlignedPosition, Alignment};
use crate::numbering::apply_numbering;
use crate::types::{NumberingRegion, Position, NumberingConfig};

// =============================================================================
// IMGT Chain Numbering Configuration
// =============================================================================

/// IMGT numbering regions (same for all chain types)
pub const IMGT_REGIONS: &[NumberingRegion] = &[
    NumberingRegion::offset(1, 26, 1),                   // FR1
    NumberingRegion::with_config(27, 38, CDR1_CONFIG),   // CDR1
    NumberingRegion::offset(39, 55, 39),                 // FR2
    NumberingRegion::with_config(56, 65, CDR2_CONFIG),   // CDR2
    NumberingRegion::offset(66, 104, 66),                // FR3
    NumberingRegion::with_config(105, 117, CDR3_CONFIG), // CDR3
    NumberingRegion::offset(118, 128, 118),              // FR4
];

/// Get IMGT-specific numbering: FR regions from alignment, CDR regions from IMGT rules
pub fn get_imgt_numbering(alignment: &Alignment) -> Vec<Position> {
    apply_numbering(&alignment.positions, IMGT_REGIONS)
}

// =============================================================================
// CDR Numbering Configurations
// =============================================================================

// CDR1: positions 27-38 (12 base), deletions from center outward, insertions at 32/33
// Fill order: 27, 38, 28, 37, 29, 36, 30, 35, 31, 34, 32, 33
// Deletion order (reverse of fill): 33, 32, 34, 31, 35, 30, 36, 29, 37, 28, 38
const CDR1_CONFIG: NumberingConfig = NumberingConfig::palindromic(
    &[27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38],
    &[33, 32, 34, 31, 35, 30, 36, 29, 37, 28, 38],
    32,
    33,
);

// CDR2: positions 56-65 (10 base), deletions from center outward, insertions at 60/61
// Fill order: 56, 65, 57, 64, 58, 63, 59, 62, 60, 61
// Deletion order (reverse of fill): 61, 60, 62, 59, 63, 58, 64, 57, 65
const CDR2_CONFIG: NumberingConfig = NumberingConfig::palindromic(
    &[56, 57, 58, 59, 60, 61, 62, 63, 64, 65],
    &[61, 60, 62, 59, 63, 58, 64, 57, 65],
    60,
    61,
);

// CDR3: positions 105-117 (13 base), deletions from center outward, insertions at 111/112
// Fill order: 105, 117, 106, 116, 107, 115, 108, 114, 109, 113, 110, 112, 111
// Deletion order (reverse of fill): 111, 112, 110, 113, 109, 114, 108, 115, 107, 116, 106, 117
const CDR3_CONFIG: NumberingConfig = NumberingConfig::palindromic(
    &[
        105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117,
    ],
    &[111, 112, 110, 113, 109, 114, 108, 115, 107, 116, 106, 117],
    111,
    112,
);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbering::number_with_config;

    #[test]
    fn test_cdr1_numbering() {
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
            // Palindromic insertions: first goes to right (33A), pattern: 32A, 33A, 33
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
            let positions = number_with_config(length, &CDR1_CONFIG);
            assert_eq!(positions.len(), length);
            let pos_strings: Vec<String> = positions.iter().map(|p| p.to_string()).collect();
            assert_eq!(pos_strings, output_strings, "Failed for length {}", length);
        }
    }

    #[test]
    fn test_cdr2_numbering() {
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
            let positions = number_with_config(length, &CDR2_CONFIG);
            assert_eq!(positions.len(), length);
            let pos_strings: Vec<String> = positions.iter().map(|p| p.to_string()).collect();
            assert_eq!(pos_strings, output_strings, "Failed for length {}", length);
        }
    }

    #[test]
    fn test_cdr3_numbering() {
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
            let positions = number_with_config(length, &CDR3_CONFIG);
            assert_eq!(positions.len(), length);
            let pos_strings: Vec<String> = positions.iter().map(|p| p.to_string()).collect();
            assert_eq!(pos_strings, output_strings, "Failed for length {}", length);
        }
    }
}
