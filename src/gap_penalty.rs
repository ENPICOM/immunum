use crate::constants::{insertion_points, ScoringParams};
use crate::schemes::SchemeConfig;
use crate::types::{Chain, RegionRange, Scheme};
use std::collections::HashSet;

/// Parameters for gap penalty calculation
pub struct GapPenaltyParams<'a> {
    pub position: u32,
    pub scheme_type: Scheme,
    pub chain_type: Chain,
    pub conserved_positions: &'a HashSet<u32>,
    pub insertion_positions: &'a [u32],
    pub gap_positions: &'a [u32],
    pub fr1: &'a RegionRange,
    pub fr2: &'a RegionRange,
    pub fr3: &'a RegionRange,
    pub fr4: &'a RegionRange,
    pub cdr1: &'a RegionRange,
    pub cdr2: &'a RegionRange,
    pub cdr3: &'a RegionRange,
    pub scoring: &'a ScoringParams,
}

/// Calculates gap penalty according to position and scheme parameters
///
/// This is a core function of the immunological numbering system that determines
/// the penalty for introducing gaps at specific positions during sequence alignment.
/// The penalty varies based on:
/// - Position type (conserved, framework, CDR, other)
/// - Numbering scheme (IMGT vs KABAT)
/// - Chain type (heavy, kappa, lambda, TCR)
/// - Proximity to insertion points
/// - Regional constraints
///
/// # Arguments
/// * `params` - All parameters needed for gap penalty calculation
///
/// # Returns
/// A tuple of (query_gap_penalty, consensus_gap_penalty)
pub fn calculate_gap_penalty(
    position: u32,
    scheme_type: Scheme,
    chain_type: Chain,
    config: &SchemeConfig,
    scoring: &ScoringParams,
) -> (f64, f64) {
    // Helper function to get all framework positions
    let framework_positions: Vec<u32> = config
        .fr1
        .positions()
        .chain(config.fr2.positions())
        .chain(config.fr3.positions())
        .chain(config.fr4.positions())
        .collect();

    // Helper function to get all CDR positions
    let cdr_positions: Vec<u32> = config
        .cdr1
        .positions()
        .chain(config.cdr2.positions())
        .chain(config.cdr3.positions())
        .collect();

    // Set initial gap penalties based on position type
    let penalty = match () {
        _ if config.conserved_positions.contains(&position) => scoring.gap_pen_cp,
        _ if framework_positions.contains(&position) => scoring.gap_pen_fr,
        _ if cdr_positions.contains(&position) => scoring.gap_pen_cdr,
        _ => scoring.gap_pen_other,
    };

    let mut query_gap_penalty = penalty;
    let mut consensus_gap_penalty = penalty;

    // Apply scheme-specific CDR penalty adjustments
    if cdr_positions.contains(&position) {
        apply_cdr_penalties(
            position,
            scheme_type,
            chain_type,
            &mut consensus_gap_penalty,
            &config.cdr1,
            &config.cdr2,
            &config.cdr3,
            scoring,
        );
        query_gap_penalty = scoring.gap_pen_cp; // TODO Maybe change to different variable
    }

    // Apply start penalty (same for all schemes)
    if position < 18 {
        apply_start_penalty(position, &mut consensus_gap_penalty);
    }

    // Apply position-specific penalties for insertion and gap positions
    if config.insertion_positions.contains(&position) {
        query_gap_penalty = scoring.gap_pen_ip;
    }

    if config.gap_positions.contains(&position) {
        consensus_gap_penalty = scoring.gap_pen_op;
    }

    (query_gap_penalty, consensus_gap_penalty)
}

/// Apply CDR-specific penalties based on the numbering scheme
fn apply_cdr_penalties(
    position: u32,
    scheme_type: Scheme,
    chain_type: Chain,
    consensus_gap_penalty: &mut f64,
    cdr1: &RegionRange,
    cdr2: &RegionRange,
    cdr3: &RegionRange,
    scoring: &ScoringParams,
) {
    match scheme_type {
        Scheme::IMGT => {
            apply_imgt_cdr_penalties(position, consensus_gap_penalty, cdr1, cdr2, cdr3, scoring)
        }
        Scheme::KABAT => apply_kabat_cdr_penalties(
            position,
            chain_type,
            consensus_gap_penalty,
            cdr1,
            cdr2,
            cdr3,
            scoring,
        ),
    }
}

/// Apply IMGT-specific CDR penalties
fn apply_imgt_cdr_penalties(
    position: u32,
    consensus_gap_penalty: &mut f64,
    cdr1: &RegionRange,
    cdr2: &RegionRange,
    cdr3: &RegionRange,
    scoring: &ScoringParams,
) {
    if cdr1.start <= position && position < cdr1.end && position != insertion_points::CDR1_IMGT {
        // CDR1 penalties
        *consensus_gap_penalty += scoring.pen_leap_insertion_point_imgt
            + scoring.cdr_increase
                * (position as isize - insertion_points::CDR1_IMGT as isize).abs() as f64;
        *consensus_gap_penalty += if position > insertion_points::CDR1_IMGT {
            0.1
        } else {
            0.0
        };
    } else if cdr2.start <= position
        && position < cdr2.end
        && position != insertion_points::CDR2_IMGT
    {
        // CDR2 penalties
        *consensus_gap_penalty += scoring.pen_leap_insertion_point_imgt
            + scoring.cdr_increase
                * (position as isize - insertion_points::CDR2_IMGT as isize).abs() as f64;
        *consensus_gap_penalty += if position > insertion_points::CDR2_IMGT {
            0.1
        } else {
            0.0
        };
    } else if cdr3.start <= position
        && position < cdr3.end
        && position != insertion_points::CDR3_IMGT
    {
        // CDR3 penalties
        *consensus_gap_penalty += scoring.pen_leap_insertion_point_imgt
            + scoring.cdr_increase
                * (position as isize - insertion_points::CDR3_IMGT as isize).abs() as f64;
        *consensus_gap_penalty += if position < insertion_points::CDR3_IMGT {
            0.1
        } else {
            0.0
        };
    }
}

/// Apply KABAT-specific CDR penalties
fn apply_kabat_cdr_penalties(
    position: u32,
    chain_type: Chain,
    consensus_gap_penalty: &mut f64,
    cdr1: &RegionRange,
    cdr2: &RegionRange,
    cdr3: &RegionRange,
    scoring: &ScoringParams,
) {
    let (cdr1_insertion_position, cdr2_insertion_position, cdr3_insertion_position) =
        match chain_type {
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

    if cdr1.start <= position && position < cdr1.end && position != cdr1_insertion_position {
        apply_kabat_cdr1_penalties(
            position,
            chain_type,
            consensus_gap_penalty,
            cdr1_insertion_position,
            scoring,
        );
    } else if cdr2.start <= position && position < cdr2.end && position != cdr2_insertion_position {
        apply_kabat_cdr2_penalties(
            position,
            chain_type,
            consensus_gap_penalty,
            cdr2_insertion_position,
            scoring,
        );
    } else if cdr3.start <= position && position < cdr3.end && position != cdr3_insertion_position {
        apply_kabat_cdr3_penalties(
            position,
            consensus_gap_penalty,
            cdr3_insertion_position,
            scoring,
        );
    }
}

/// Apply KABAT CDR1-specific penalties
fn apply_kabat_cdr1_penalties(
    position: u32,
    chain_type: Chain,
    consensus_gap_penalty: &mut f64,
    cdr1_insertion_position: u32,
    scoring: &ScoringParams,
) {
    if chain_type == Chain::IGH {
        // Heavy chain CDR1
        *consensus_gap_penalty += scoring.pen_leap_insertion_point_kabat
            + (scoring.cdr_increase
                * (position as isize - cdr1_insertion_position as isize).abs() as f64);
    } else {
        // Light chain CDR1 - special exception for correct placement
        let relative_position: i32 = position as i32 - cdr1_insertion_position as i32;
        let penalty_addition: f64 = if relative_position < 0 {
            scoring.pen_leap_insertion_point_kabat
                + (scoring.cdr_increase
                    * (position as isize - cdr1_insertion_position as isize).abs() as f64)
        } else if relative_position < 3 {
            // Positions 29 and 30 have lower penalty
            scoring.gap_pen_cdr + (0.1 * relative_position as f64)
        } else {
            // Positions 31-34 - subtract 1 times cdr_increase to align penalties of 26 and 31
            scoring.pen_leap_insertion_point_kabat
                + (scoring.cdr_increase
                    * (((position as isize - cdr1_insertion_position as isize).abs() as f64) - 1.0))
        };
        *consensus_gap_penalty += penalty_addition;
    }
}

/// Apply KABAT CDR2-specific penalties
fn apply_kabat_cdr2_penalties(
    position: u32,
    chain_type: Chain,
    consensus_gap_penalty: &mut f64,
    cdr2_insertion_position: u32,
    scoring: &ScoringParams,
) {
    // Exception in heavy scheme: positions after 54 get framework penalty
    if chain_type == Chain::IGH && position > 54 {
        *consensus_gap_penalty = scoring.gap_pen_fr;
    } else {
        *consensus_gap_penalty += scoring.pen_leap_insertion_point_kabat
            + (scoring.cdr_increase
                * (position as isize - cdr2_insertion_position as isize).abs() as f64);
    }
}

/// Apply KABAT CDR3-specific penalties
fn apply_kabat_cdr3_penalties(
    position: u32,
    consensus_gap_penalty: &mut f64,
    cdr3_insertion_position: u32,
    scoring: &ScoringParams,
) {
    *consensus_gap_penalty += scoring.pen_leap_insertion_point_kabat
        + (scoring.cdr_increase
            * (position as isize - cdr3_insertion_position as isize).abs() as f64);

    if position > cdr3_insertion_position {
        // Higher penalty after insertion in CDR3
        *consensus_gap_penalty += 4.0 * scoring.cdr_increase;
    }
}

/// Apply start penalty for early positions (same for all schemes)
fn apply_start_penalty(position: u32, consensus_gap_penalty: &mut f64) {
    // Increase from 2 to 11, then add 0.1 until position 18
    // TODO: This doesn't have effect yet because these are restricted sites
    *consensus_gap_penalty = 1.0
        + (if position > 10 { 10.0 } else { position as f64 })
        + (if position > 10 {
            0.1 * (position as f64 - 10.0)
        } else {
            0.0
        });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::get_scoring_params;

    fn create_test_regions() -> (
        RegionRange,
        RegionRange,
        RegionRange,
        RegionRange,
        RegionRange,
        RegionRange,
        RegionRange,
    ) {
        let fr1 = RegionRange { start: 1, end: 27 };
        let cdr1 = RegionRange { start: 27, end: 39 };
        let fr2 = RegionRange { start: 39, end: 56 };
        let cdr2 = RegionRange { start: 56, end: 66 };
        let fr3 = RegionRange {
            start: 66,
            end: 105,
        };
        let cdr3 = RegionRange {
            start: 105,
            end: 118,
        };
        let fr4 = RegionRange {
            start: 118,
            end: 129,
        };
        (fr1, cdr1, fr2, cdr2, fr3, cdr3, fr4)
    }

    #[test]
    fn test_basic_gap_penalty() {
        let (fr1, cdr1, fr2, cdr2, fr3, cdr3, fr4) = create_test_regions();
        let conserved_positions = vec![23, 41, 104, 118, 119, 121];
        let insertion_positions = vec![32, 60, 111];
        let gap_positions = vec![10, 73];
        let scoring_params = get_scoring_params();
        let config = SchemeConfig {
            conserved_positions,
            insertion_positions,
            gap_positions,
            fr1,
            cdr1,
            fr2,
            cdr2,
            fr3,
            cdr3,
            fr4,
        };

        // Test conserved position
        let (query_pen, consensus_pen) =
            calculate_gap_penalty(23, Scheme::IMGT, Chain::IGH, &config, &scoring_params);
        assert_eq!(query_pen, scoring_params.gap_pen_cp);
        assert_eq!(consensus_pen, scoring_params.gap_pen_cp);

        // Test framework position
        let (query_pen, consensus_pen) =
            calculate_gap_penalty(25, Scheme::IMGT, Chain::IGH, &config, &scoring_params);
        assert_eq!(query_pen, scoring_params.gap_pen_fr);
        assert_eq!(consensus_pen, scoring_params.gap_pen_fr);

        // Test gap position
        let (query_pen, consensus_pen) =
            calculate_gap_penalty(10, Scheme::IMGT, Chain::IGH, &config, &scoring_params);
        assert_eq!(query_pen, scoring_params.gap_pen_fr);
        assert_eq!(consensus_pen, scoring_params.gap_pen_op);
    }

    #[test]
    fn test_insertion_position() {
        let (fr1, cdr1, fr2, cdr2, fr3, cdr3, fr4) = create_test_regions();
        let conserved_positions: Vec<u32> = vec![];
        let insertion_positions = vec![32, 60, 111];
        let gap_positions = vec![10, 73];
        let scoring_params = get_scoring_params();
        let config = SchemeConfig {
            conserved_positions,
            insertion_positions,
            gap_positions,
            fr1,
            cdr1,
            fr2,
            cdr2,
            fr3,
            cdr3,
            fr4,
        };

        let (query_pen, _) =
            calculate_gap_penalty(32, Scheme::IMGT, Chain::IGH, &config, &scoring_params);
        assert_eq!(query_pen, scoring_params.gap_pen_ip);
    }
}
