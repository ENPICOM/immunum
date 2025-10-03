use crate::constants::{insertion_points, ScoringParams};
use crate::insertion_numbering::name_insertions;
use crate::needleman_wunsch::needleman_wunsch_consensus;
use crate::scoring_matrix::ScoringMatrix;
use crate::types::{CdrDefinitions, Chain, NumberingPosition, RegionRange, Scheme};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct NumberingScheme {
    pub scheme_type: Scheme,
    pub chain_type: Chain,
    pub cdr_definition: CdrDefinitions,
    pub insertion_positions: Vec<u32>,
    pub gap_positions: Vec<u32>,
    pub consensus_amino_acids: HashMap<u32, Vec<u8>>,
    pub scoring_matrix: ScoringMatrix,
    pub fr1: RegionRange,
    pub fr2: RegionRange,
    pub fr3: RegionRange,
    pub fr4: RegionRange,
    pub cdr1: RegionRange,
    pub cdr2: RegionRange,
    pub cdr3: RegionRange,
    pub conserved_positions_set: HashSet<u32>,
    pub restricted_sites_set: HashSet<u32>,
    // K-mer set for prefiltering (computed during initialization)
    pub kmer_set: HashSet<String>,
}

impl NumberingScheme {
    /// All framework positions
    pub fn framework_positions(&self) -> Vec<u32> {
        self.fr1
            .positions()
            .chain(self.fr2.positions())
            .chain(self.fr3.positions())
            .chain(self.fr4.positions())
            .collect()
    }
    /// All cdr positions
    pub fn cdr_positions(&self) -> Vec<u32> {
        self.cdr1
            .positions()
            .chain(self.cdr2.positions())
            .chain(self.cdr3.positions())
            .collect()
    }
    /// Calculates gap penalty according to position and scheme
    pub fn gap_penalty(&self, position: u32, scoring: &ScoringParams) -> (f64, f64) {
        // Set initial gap penalties
        let penalty = if self.conserved_positions_set.contains(&position) {
            scoring.gap_pen_cp
        } else if self.framework_positions().contains(&position) {
            scoring.gap_pen_fr
        } else if self.cdr_positions().contains(&position) {
            scoring.gap_pen_cdr
        } else {
            scoring.gap_pen_other
        };

        let mut query_gap_penalty = penalty;
        let mut consensus_gap_penalty = penalty;

        // ADAPT CDR PENALTIES, different for every scheme
        // IMGT
        if self.cdr_positions().contains(&position) {
            if self.scheme_type == Scheme::IMGT {
                if self.cdr1.start <= position
                    && position < self.cdr1.end
                    && position != insertion_points::CDR1_IMGT
                {
                    // cdr1
                    consensus_gap_penalty += scoring.pen_leap_insertion_point_imgt
                        + scoring.cdr_increase
                            * (position as isize - insertion_points::CDR1_IMGT as isize).abs()
                                as f64;
                    consensus_gap_penalty += if position > insertion_points::CDR1_IMGT {
                        0.1
                    } else {
                        0.0
                    };
                } else if self.cdr2.start <= position
                    && position < self.cdr2.end
                    && position != insertion_points::CDR2_IMGT
                {
                    // cdr2
                    consensus_gap_penalty += scoring.pen_leap_insertion_point_imgt
                        + scoring.cdr_increase
                            * (position as isize - insertion_points::CDR2_IMGT as isize).abs()
                                as f64;
                    consensus_gap_penalty += if position > insertion_points::CDR2_IMGT {
                        0.1
                    } else {
                        0.0
                    };
                } else if self.cdr3.start <= position
                    && position < self.cdr3.end
                    && position != insertion_points::CDR3_IMGT
                {
                    // cdr3
                    consensus_gap_penalty += scoring.pen_leap_insertion_point_imgt
                        + scoring.cdr_increase
                            * (position as isize - insertion_points::CDR3_IMGT as isize).abs()
                                as f64;
                    consensus_gap_penalty += if position < insertion_points::CDR3_IMGT {
                        0.1
                    } else {
                        0.0
                    };
                }
            }
            // KABAT
            else if self.scheme_type == Scheme::KABAT {
                let (cdr1_insertion_position, cdr2_insertion_position, cdr3_insertion_position) =
                    match self.chain_type {
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

                if self.cdr1.start <= position
                    && position < self.cdr1.end
                    && position != cdr1_insertion_position
                {
                    if self.chain_type == Chain::IGH {
                        // cdr1
                        consensus_gap_penalty += scoring.pen_leap_insertion_point_kabat
                            + (scoring.cdr_increase
                                * (position as isize - cdr1_insertion_position as isize).abs()
                                    as f64);
                    } else {
                        // Exception in cdr 1 of kabat light schemes to ensure correct placement
                        let relative_position: i32 =
                            position as i32 - cdr1_insertion_position as i32;
                        let penalty_addition: f64 = if relative_position < 0 {
                            scoring.pen_leap_insertion_point_kabat
                                + (scoring.cdr_increase
                                    * (position as isize - cdr1_insertion_position as isize).abs()
                                        as f64)
                        } else if relative_position < 3 {
                            // positions 29 and 30 have lower penalty
                            scoring.gap_pen_cdr + (0.1 * relative_position as f64)
                        } else {
                            // positions 31-34
                            // -1 times cdr_increase to line up penalties of 26 and 31
                            scoring.pen_leap_insertion_point_kabat
                                + (scoring.cdr_increase
                                    * (((position as isize - cdr1_insertion_position as isize).abs()
                                        as f64)
                                        - 1.0))
                        };
                        consensus_gap_penalty += penalty_addition
                    }
                } else if self.cdr2.start <= position
                    && position < self.cdr2.end
                    && position != cdr2_insertion_position
                {
                    // cdr2
                    // Exception in heavy scheme:
                    if self.chain_type == Chain::IGH && position > 54 {
                        consensus_gap_penalty = scoring.gap_pen_fr;
                    } else {
                        consensus_gap_penalty += scoring.pen_leap_insertion_point_kabat
                            + (scoring.cdr_increase
                                * (position as isize - cdr2_insertion_position as isize).abs()
                                    as f64);
                    }
                } else if self.cdr3.start <= position
                    && position < self.cdr3.end
                    && position != cdr3_insertion_position
                {
                    // cdr3
                    consensus_gap_penalty += scoring.pen_leap_insertion_point_kabat
                        + (scoring.cdr_increase
                            * (position as isize - cdr3_insertion_position as isize).abs() as f64);
                    if position > cdr3_insertion_position {
                        // higher penalty after insertion in cdr3
                        consensus_gap_penalty += 4.0 * scoring.cdr_increase;
                    }
                }
                // if self.chain_type != Chain::IGH && position == cdr1_insertion_position {
                //     consensus_gap_penalty = scoring.gap_pen_other;
                // }
            }
            query_gap_penalty = scoring.gap_pen_cp // TODO Maybe change to different variable
        }

        // Handle start penalty, same for all schemes
        if position < 18 {
            // Increase from 2 to 11, then add 0.1 until position 18
            // TODO does not have effect yet, because these are restriced sites
            consensus_gap_penalty = 1.0
                + (if position > 10 { 10.0 } else { position as f64 })
                + (if position > 10 {
                    0.1 * (position as f64 - 10.0)
                } else {
                    0.0
                });
        }

        // Adapt only query or consensus gap for insertion and gap positions
        if self.insertion_positions.contains(&position) {
            query_gap_penalty = scoring.gap_pen_ip
        }

        if self.gap_positions.contains(&position) {
            consensus_gap_penalty = scoring.gap_pen_op;
        }
        (query_gap_penalty, consensus_gap_penalty)
    }

    /// Number a sequence by alignment, returns numbers and score
    pub(crate) fn number_sequence(
        &self,
        query_sequence: &[u8],
    ) -> (Vec<NumberingPosition>, f64, usize, usize) {
        let (mut numbering, identity) = needleman_wunsch_consensus(query_sequence, self);

        // give gap positions correct names as defined by the numbering scheme
        name_insertions(&mut numbering, &self.scheme_type);

        // find start and end index
        // TODO do we need to handle numberings with no positions? (only '-')
        let start = numbering
            .iter()
            .position(|s| s != "-")
            .expect("No positions numbered");

        // Find the last non-dash element (searching from the end)
        let end = numbering
            .iter()
            .rposition(|s| s != "-")
            .expect("No positions numbered");

        (numbering, identity, start, end)
    }
}
