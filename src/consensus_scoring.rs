use crate::constants::{RegionRanges, ScoringParams, ENCODED_RESIDUES_MAP};
use crate::constants::{ACCEPTED_RESIDUES, BLOSUM62};
use crate::scoring_matrix::ScoringMatrix;
use crate::types::{Chain, Scheme};
use std::collections::{HashMap, HashSet};

/// Check if a position is in any framework region
fn is_framework_position(position: u32, region_ranges: &RegionRanges) -> bool {
    region_ranges.fr1.positions().contains(&position)
        || region_ranges.fr2.positions().contains(&position)
        || region_ranges.fr3.positions().contains(&position)
        || region_ranges.fr4.positions().contains(&position)
}

/// Check if a position is in any CDR region
fn is_cdr_position(position: u32, region_ranges: &RegionRanges) -> bool {
    region_ranges.cdr1.positions().contains(&position)
        || region_ranges.cdr2.positions().contains(&position)
        || region_ranges.cdr3.positions().contains(&position)
}

/// Calculate base penalty based on position type
fn calculate_base_penalty(
    position: u32,
    scoring: &ScoringParams,
    config: &GapPenaltyConfig,
) -> f64 {
    if config.conserved_positions_set.contains(&position) {
        scoring.gap_pen_cp
    } else if is_framework_position(position, &config.region_ranges) {
        scoring.gap_pen_fr
    } else if is_cdr_position(position, &config.region_ranges) {
        scoring.gap_pen_cdr
    } else {
        scoring.gap_pen_other
    }
}

/// Apply IMGT-specific CDR penalty adjustments
fn apply_imgt_cdr_penalty_adjustments(
    position: u32,
    scoring: &ScoringParams,
    config: &GapPenaltyConfig,
    consensus_gap_penalty: f64,
) -> f64 {
    use crate::constants::insertion_points;

    let mut penalty = consensus_gap_penalty;

    if config.region_ranges.cdr1.start <= position
        && position < config.region_ranges.cdr1.end
        && position != insertion_points::CDR1_IMGT
    {
        // cdr1
        penalty += scoring.pen_leap_insertion_point_imgt
            + scoring.cdr_increase
                * (position as isize - insertion_points::CDR1_IMGT as isize).abs() as f64;
        penalty += if position > insertion_points::CDR1_IMGT {
            0.1
        } else {
            0.0
        };
    } else if config.region_ranges.cdr2.start <= position
        && position < config.region_ranges.cdr2.end
        && position != insertion_points::CDR2_IMGT
    {
        // cdr2
        penalty += scoring.pen_leap_insertion_point_imgt
            + scoring.cdr_increase
                * (position as isize - insertion_points::CDR2_IMGT as isize).abs() as f64;
        penalty += if position > insertion_points::CDR2_IMGT {
            0.1
        } else {
            0.0
        };
    } else if config.region_ranges.cdr3.start <= position
        && position < config.region_ranges.cdr3.end
        && position != insertion_points::CDR3_IMGT
    {
        // cdr3
        penalty += scoring.pen_leap_insertion_point_imgt
            + scoring.cdr_increase
                * (position as isize - insertion_points::CDR3_IMGT as isize).abs() as f64;
        penalty += if position < insertion_points::CDR3_IMGT {
            0.1
        } else {
            0.0
        };
    }

    penalty
}

/// Apply KABAT-specific CDR penalty adjustments
fn apply_kabat_cdr_penalty_adjustments(
    position: u32,
    scoring: &ScoringParams,
    config: &GapPenaltyConfig,
    consensus_gap_penalty: f64,
) -> f64 {
    use crate::constants::insertion_points;

    let mut penalty = consensus_gap_penalty;

    let (cdr1_insertion_position, cdr2_insertion_position, cdr3_insertion_position) =
        match config.chain_type {
            Chain::IGH => (
                insertion_points::CDR1_KABAT_HEAVY,
                insertion_points::CDR2_KABAT_HEAVY,
                insertion_points::CDR3_KABAT_HEAVY,
            ),
            Chain::IGL | Chain::IGK => (
                insertion_points::CDR1_KABAT_LIGHT,
                insertion_points::CDR2_KABAT_LIGHT,
                insertion_points::CDR3_KABAT_LIGHT,
            ),
            _ => panic!("This scheme is not implemented yet (calculating gap penalty)"),
        };

    if config.region_ranges.cdr1.start <= position
        && position < config.region_ranges.cdr1.end
        && position != cdr1_insertion_position
    {
        if config.chain_type == Chain::IGH {
            // cdr1
            penalty += scoring.pen_leap_insertion_point_kabat
                + (scoring.cdr_increase
                    * (position as isize - cdr1_insertion_position as isize).abs() as f64);
        } else {
            // Exception in cdr 1 of kabat light schemes to ensure correct placement
            let relative_position: i32 = position as i32 - cdr1_insertion_position as i32;
            let penalty_addition: f64 = if relative_position < 0 {
                scoring.pen_leap_insertion_point_kabat
                    + (scoring.cdr_increase
                        * (position as isize - cdr1_insertion_position as isize).abs() as f64)
            } else if relative_position < 3 {
                // positions 29 and 30 have lower penalty
                scoring.gap_pen_cdr + (0.1 * relative_position as f64)
            } else {
                // positions 31-34
                // -1 times cdr_increase to line up penalties of 26 and 31
                scoring.pen_leap_insertion_point_kabat
                    + (scoring.cdr_increase
                        * (((position as isize - cdr1_insertion_position as isize).abs() as f64)
                            - 1.0))
            };
            penalty += penalty_addition
        }
    } else if config.region_ranges.cdr2.start <= position
        && position < config.region_ranges.cdr2.end
        && position != cdr2_insertion_position
    {
        // cdr2
        // Exception in heavy scheme:
        if config.chain_type == Chain::IGH && position > 54 {
            penalty = scoring.gap_pen_fr;
        } else {
            penalty += scoring.pen_leap_insertion_point_kabat
                + (scoring.cdr_increase
                    * (position as isize - cdr2_insertion_position as isize).abs() as f64);
        }
    } else if config.region_ranges.cdr3.start <= position
        && position < config.region_ranges.cdr3.end
        && position != cdr3_insertion_position
    {
        // cdr3
        penalty += scoring.pen_leap_insertion_point_kabat
            + (scoring.cdr_increase
                * (position as isize - cdr3_insertion_position as isize).abs() as f64);
        if position > cdr3_insertion_position {
            // higher penalty after insertion in cdr3
            penalty += 4.0 * scoring.cdr_increase;
        }
    }

    penalty
}

/// Apply start position penalty adjustment (same for all schemes)
fn apply_start_position_penalty(position: u32, consensus_gap_penalty: f64) -> f64 {
    if position < 18 {
        // Increase from 2 to 11, then add 0.1 until position 18
        // TODO does not have effect yet, because these are restricted sites
        1.0 + (if position > 10 { 10.0 } else { position as f64 })
            + (if position > 10 {
                0.1 * (position as f64 - 10.0)
            } else {
                0.0
            })
    } else {
        consensus_gap_penalty
    }
}

/// Apply insertion and gap position adjustments
fn apply_insertion_gap_adjustments(
    position: u32,
    scoring: &ScoringParams,
    config: &GapPenaltyConfig,
    query_gap_penalty: f64,
    consensus_gap_penalty: f64,
) -> (f64, f64) {
    let mut query_penalty = query_gap_penalty;
    let mut consensus_penalty = consensus_gap_penalty;

    // Adapt only query or consensus gap for insertion and gap positions
    if config.insertion_positions.contains(&position) {
        query_penalty = scoring.gap_pen_ip;
    }

    if config.gap_positions.contains(&position) {
        consensus_penalty = scoring.gap_pen_op;
    }

    (query_penalty, consensus_penalty)
}

/// Configuration for gap penalty calculation
pub struct GapPenaltyConfig {
    pub conserved_positions_set: HashSet<u32>,
    pub scheme_type: Scheme,
    pub chain_type: Chain,
    pub insertion_positions: Vec<u32>,
    pub gap_positions: Vec<u32>,
    pub region_ranges: RegionRanges,
}

/// Calculate scoring matrix from consensus data and scoring parameters
pub fn calculate_scoring_matrix(
    consensus: &HashMap<u32, Vec<u8>>,
    scoring_params: &ScoringParams,
    config: &GapPenaltyConfig,
) -> ScoringMatrix {
    let number_accepted_residues = ACCEPTED_RESIDUES.len();
    let consensus_length = consensus.len();

    let mut matrix = ScoringMatrix::zeros(consensus_length, number_accepted_residues + 2);

    // Fill in scores
    for consensus_position in 0..consensus_length {
        for (residue_index, &residue) in ACCEPTED_RESIDUES
            .iter()
            .enumerate()
            .take(number_accepted_residues)
        {
            let score =
                best_score_consensus((consensus_position + 1) as u32, residue, consensus) as f64;
            matrix.set(consensus_position, residue_index, score);
        }

        // Fill in gap penalties
        let (query_gap, consensus_gap) =
            calculate_gap_penalty((consensus_position + 1) as u32, scoring_params, config);

        // Query gap
        matrix.set(consensus_position, number_accepted_residues, query_gap);
        // Consensus gap
        matrix.set(
            consensus_position,
            number_accepted_residues + 1,
            consensus_gap,
        );
    }

    matrix
}

/// Finds highest score according to blosum62
fn best_score_consensus(position: u32, residue: u8, consensus: &HashMap<u32, Vec<u8>>) -> i32 {
    // Initialize to minimal score
    let to_check_residues: &Vec<u8> = consensus
        .get(&position)
        .expect("Position outside of consensus");

    // when any is allowed, best score is perfect match
    if to_check_residues.iter().all(|&c| c == b'-') {
        return blosum_lookup(&residue, &residue);
    }
    // iterate over residues to get max score
    blosum_lookup(
        &residue,
        to_check_residues
            .iter()
            .max_by_key(|&c| blosum_lookup(&residue, c))
            .unwrap(),
    )
}

/// Utility to easily find blosum62 score for two residues
fn blosum_lookup(residue1: &u8, residue2: &u8) -> i32 {
    if residue1 < residue2 {
        let lookup: &[u8; 2] = &[*residue1, *residue2];
        *BLOSUM62.get(lookup).unwrap()
    } else {
        let lookup: &[u8; 2] = &[*residue2, *residue1];
        *BLOSUM62.get(lookup).unwrap()
    }
}

/// Encodes sequence to indexes for the scoring matrix
pub(crate) fn encode_sequence(input: &[u8]) -> Vec<u8> {
    // Create a lookup table (once, could be static)
    input
        .iter()
        .map(|&residue| *ENCODED_RESIDUES_MAP.get(&residue).unwrap_or(&128))
        .collect()
}

/// Calculate gap penalty for a given position based on scheme configuration
pub fn calculate_gap_penalty(
    position: u32,
    scoring: &ScoringParams,
    config: &GapPenaltyConfig,
) -> (f64, f64) {
    // Step 1: Calculate base penalty based on position type
    let base_penalty = calculate_base_penalty(position, scoring, config);
    let mut query_gap_penalty = base_penalty;
    let mut consensus_gap_penalty = base_penalty;

    // Step 2: Apply scheme-specific CDR penalty adjustments
    if is_cdr_position(position, &config.region_ranges) {
        match config.scheme_type {
            Scheme::IMGT => {
                consensus_gap_penalty = apply_imgt_cdr_penalty_adjustments(
                    position,
                    scoring,
                    config,
                    consensus_gap_penalty,
                );
            }
            Scheme::KABAT => {
                consensus_gap_penalty = apply_kabat_cdr_penalty_adjustments(
                    position,
                    scoring,
                    config,
                    consensus_gap_penalty,
                );
            }
        }
        // For CDR positions, set query gap penalty to conserved position penalty
        query_gap_penalty = scoring.gap_pen_cp;
    }

    // Step 3: Apply start position penalty adjustment (same for all schemes)
    consensus_gap_penalty = apply_start_position_penalty(position, consensus_gap_penalty);

    // Step 4: Apply insertion and gap position adjustments
    let (final_query_penalty, final_consensus_penalty) = apply_insertion_gap_adjustments(
        position,
        scoring,
        config,
        query_gap_penalty,
        consensus_gap_penalty,
    );

    (final_query_penalty, final_consensus_penalty)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lookup() {
        // working
        assert_eq!(blosum_lookup(&b'B', &b'N'), 3);
        assert_eq!(blosum_lookup(&b'W', &b'Y'), 2);
        assert_eq!(blosum_lookup(&b'X', &b'E'), -1);
        // same either way of inputting arguments
        assert_eq!(blosum_lookup(&b'C', &b'S'), blosum_lookup(&b'S', &b'C'));
    }

    #[test]
    fn amino_acid_encoding() {
        assert_eq!(
            encode_sequence("ABCDEFGHIKLMNPQRSTVWXYZ".as_bytes()),
            (0..23).collect::<Vec<u8>>()
        )
    }
}
