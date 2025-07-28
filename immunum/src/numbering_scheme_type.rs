use crate::constants::{insertion_points, scoring};
use crate::insertion_numbering::name_insertions;
use crate::needleman_wunsch::needleman_wunsch_consensus;
use crate::types::{Chain, RegionRange, Scheme};
use ndarray::Array2;
use std::collections::HashMap;

#[derive(Debug)]
pub struct NumberingScheme {
    pub name: String,
    pub description: String,
    pub scheme_type: Scheme,
    pub chain_type: Chain,
    pub conserved_positions: Vec<u32>,
    pub insertion_positions: Vec<u32>,
    pub gap_positions: Vec<u32>,
    pub consensus_amino_acids: HashMap<u32, Vec<u8>>,
    pub scoring_matrix: Array2<f64>,
    pub fr1: RegionRange,
    pub fr2: RegionRange,
    pub fr3: RegionRange,
    pub fr4: RegionRange,
    pub cdr1: RegionRange,
    pub cdr2: RegionRange,
    pub cdr3: RegionRange,
}

#[derive(Debug, Clone)]
pub struct NumberingOutput<'a> {
    pub scheme: &'a NumberingScheme,
    pub sequence: &'a [u8],
    pub numbering: Vec<String>,
    pub identity: f64,
    pub start: u32,
    pub end: u32,
}

impl NumberingScheme {
    /// Return restricted sites according to Antpack definition (non '-' positions)
    pub fn restricted_sites(&self) -> Vec<u32> {
        let mut sites = Vec::new();
        for (&key, value) in &self.consensus_amino_acids {
            if !value.contains(&b'-') {
                sites.push(key);
            }
        }
        sites
    }
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
    pub fn gap_penalty(&self, position: u32) -> (f64, f64) {
        // Set initial gap penalties
        let mut penalty = match () {
            _ if self.conserved_positions.contains(&position) => scoring::GAP_PEN_CP,
            _ if self.framework_positions().contains(&position) => scoring::GAP_PEN_FR,
            _ if self.cdr_positions().contains(&position) => scoring::GAP_PEN_CDR,
            _ => scoring::GAP_PEN_OTHER,
        };

        // ADAPT CDR PENALTIES, different for every scheme
        // IMGT
        if self.cdr_positions().contains(&position) {
            if self.scheme_type == Scheme::IMGT {
                if self.cdr1.start <= position
                    && position < self.cdr1.end
                    && position != insertion_points::CDR1_IMGT
                {
                    // cdr1
                    penalty += scoring::PEN_LEAP_FROM_INSERTION_POINT_IMGT
                        + scoring::CDR_INCREASE
                            * (position as isize - insertion_points::CDR1_IMGT as isize).abs()
                                as f64;
                    penalty += if position > insertion_points::CDR1_IMGT {
                        0.1
                    } else {
                        0.0
                    };
                } else if self.cdr2.start <= position
                    && position < self.cdr2.end
                    && position != insertion_points::CDR2_IMGT
                {
                    // cdr2
                    penalty += scoring::PEN_LEAP_FROM_INSERTION_POINT_IMGT
                        + scoring::CDR_INCREASE
                            * (position as isize - insertion_points::CDR2_IMGT as isize).abs()
                                as f64;
                    penalty += if position > insertion_points::CDR2_IMGT {
                        0.1
                    } else {
                        0.0
                    };
                } else if self.cdr3.start <= position
                    && position < self.cdr3.end
                    && position != insertion_points::CDR3_IMGT
                {
                    // cdr3
                    penalty += scoring::PEN_LEAP_FROM_INSERTION_POINT_IMGT
                        + scoring::CDR_INCREASE
                            * (position as isize - insertion_points::CDR3_IMGT as isize).abs()
                                as f64;
                    penalty += if position < insertion_points::CDR3_IMGT {
                        0.1
                    } else {
                        0.0
                    };
                }
            }
            // KABAT
            else if self.scheme_type == Scheme::KABAT {
                let cdr1_insertion_position = if self.chain_type == Chain::IGH {
                    insertion_points::CDR1_KABAT_HEAVY
                } else {
                    insertion_points::CDR1_KABAT_LIGHT
                };

                let cdr2_insertion_position = if self.chain_type == Chain::IGH {
                    insertion_points::CDR2_KABAT_HEAVY
                } else {
                    insertion_points::CDR2_KABAT_LIGHT
                };

                let cdr3_insertion_position = if self.chain_type == Chain::IGH {
                    insertion_points::CDR3_KABAT_HEAVY
                } else {
                    insertion_points::CDR3_KABAT_LIGHT
                };

                if self.cdr1.start <= position
                    && position < self.cdr1.end
                    && position != cdr1_insertion_position
                {
                    // cdr1
                    penalty += scoring::PEN_LEAP_INSERTION_POINT_KABAT
                        + (scoring::CDR_INCREASE
                            * (position as isize - cdr1_insertion_position as isize).abs() as f64);
                } else if self.cdr2.start <= position
                    && position < self.cdr2.end
                    && position != cdr2_insertion_position
                {
                    // cdr2
                    // Exception in heavy scheme:
                    if self.chain_type == Chain::IGH && position > 54 {
                        penalty = scoring::GAP_PEN_FR;
                    } else {
                        penalty += scoring::PEN_LEAP_INSERTION_POINT_KABAT
                            + (scoring::CDR_INCREASE
                                * (position as isize - cdr2_insertion_position as isize).abs()
                                    as f64);
                    }
                } else if self.cdr3.start <= position
                    && position < self.cdr3.end
                    && position != cdr3_insertion_position
                {
                    // cdr3
                    penalty += scoring::PEN_LEAP_INSERTION_POINT_KABAT
                        + (scoring::CDR_INCREASE
                            * (position as isize - cdr3_insertion_position as isize).abs() as f64);
                    if position > cdr3_insertion_position {
                        // higher penalty after insertion in cdr3
                        penalty += 5.0;
                    }
                }

                if self.chain_type != Chain::IGH && position == cdr1_insertion_position {
                    penalty = scoring::GAP_PEN_OTHER;
                }
            }
        }

        let mut query_gap_penalty = penalty;
        let mut consensus_gap_penalty = penalty;

        // Handle start penalty, same for all schemes
        if position < 18 {
            // Increase from 2 to 11, then add 0.1 until position 18
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
            query_gap_penalty = scoring::GAP_PEN_IP;
        }

        if self.gap_positions.contains(&position) {
            consensus_gap_penalty = scoring::GAP_PEN_OP;
            // TODO set gap penalty to other?
            // query_gap_penalty = GAP_PEN_OTHER;
        }

        (query_gap_penalty, consensus_gap_penalty)
    }
    /// numbers sequence, returns NumberingOutput
    pub(crate) fn number_sequence<'a>(&'a self, query_sequence: &'a [u8]) -> NumberingOutput<'a> {
        let (mut numbering, identity) = needleman_wunsch_consensus(query_sequence, self);

        name_insertions(&mut numbering, &self.scheme_type);

        // find start and end index
        // TODO handle numberings with no positions (only '-')
        let start = numbering
            .iter()
            .position(|s| s != "-")
            .expect("No positions numbered");

        // Find the last non-dash element (searching from the end)
        let end = numbering
            .iter()
            .rposition(|s| s != "-")
            .expect("No positions numbered");

        NumberingOutput {
            scheme: self,
            sequence: query_sequence,
            numbering,
            identity,
            start: start as u32,
            end: end as u32,
        }
    }
}
