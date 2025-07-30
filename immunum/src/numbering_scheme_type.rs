use crate::constants::{insertion_points, ScoringParams};
use crate::insertion_numbering::name_insertions;
use crate::needleman_wunsch::needleman_wunsch_consensus;
use crate::types::{Chain, RegionRange, Scheme};
use ndarray::{s, Array2};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct NumberingScheme {
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

impl NumberingOutput<'_> {
    /// Get slice of sequence corresponding to a framework or cdr region
    pub fn get_query_region(&self, region: &RegionRange) -> &str {
        let mut start_position = region.start;
        let mut end_position = region.end - 1;

        // get final position
        let mut final_position: usize = 0;
        for p in 0..self.numbering.len() {
            if self.numbering[p] != "-" {
                final_position = p
            }
        }

        let mut start_sequence: Option<usize> = None;
        let mut end_sequence: Option<usize> = None;

        // find start position in sequence
        while start_position <= end_position {
            let pos_option: Option<usize> = self
                .numbering
                .iter()
                .position(|num| num == &start_position.to_string());
            match pos_option {
                Some(i) => {
                    start_sequence = Some(i);
                    break;
                }
                None => {
                    start_position += 1;
                }
            }
        }
        // check if start is found, if not, region is not present
        let start_index = match start_sequence {
            Some(i) => i,
            None => return "",
        };

        // find end position
        while end_position < self.scheme.consensus_amino_acids.len() as u32 {
            let pos_option: Option<usize> = self
                .numbering
                .iter()
                .position(|num| num == &end_position.to_string());
            match pos_option {
                Some(i) => {
                    end_sequence = Some(i);
                    break;
                }
                None => {
                    end_position += 1;
                }
            }
        }
        let end_index = end_sequence.unwrap_or(final_position);

        std::str::from_utf8(&self.sequence[start_index..=end_index])
            .expect("Non-UTF8 character in sequence")
    }

    pub fn get_output_string(&self) -> String {
        let mut output_str = String::new();
        output_str.push_str(
            std::str::from_utf8(self.sequence)
                .expect("Non-UTF8 character in sequence"),
        );
        output_str.push('\t');
        output_str.push_str(&self.numbering.join(","));
        output_str.push('\t');
        output_str.push_str(&format!("{}", self.identity));
        output_str.push('\t');
        output_str.push_str(match self.scheme.chain_type {
            Chain::IGH => "H",
            Chain::IGK => "K",
            Chain::IGL => "L",
            Chain::TRA => "A",
            Chain::TRB => "B",
            Chain::TRD => "D",
            Chain::TRG => "G",
        });
        output_str.push('\t');

        output_str.push_str(self.get_query_region(&self.scheme.cdr1));
        output_str.push('\t');
        output_str.push_str(self.get_query_region(&self.scheme.cdr2));
        output_str.push('\t');
        output_str.push_str(self.get_query_region(&self.scheme.cdr3));
        output_str.push('\t');
        output_str.push_str(self.get_query_region(&self.scheme.fr1));
        output_str.push('\t');
        output_str.push_str(self.get_query_region(&self.scheme.fr2));
        output_str.push('\t');
        output_str.push_str(self.get_query_region(&self.scheme.fr3));
        output_str.push('\t');
        output_str.push_str(self.get_query_region(&self.scheme.fr4));

        output_str.push('\t');
        output_str.push_str(&format!("{}", self.start));

        output_str.push('\t');
        output_str.push_str(&format!("{}", self.end));

        output_str.push('\n');
        output_str
    }
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
    pub fn gap_penalty(&self, position: u32, scoring:&ScoringParams) -> (f64, f64) {
        // Set initial gap penalties
        let mut penalty = match () {
            _ if self.conserved_positions.contains(&position) => scoring.gap_pen_cp,
            _ if self.framework_positions().contains(&position) => scoring.gap_pen_fr,
            _ if self.cdr_positions().contains(&position) => scoring.gap_pen_cdr,
            _ => scoring.gap_pen_other,
        };

        // ADAPT CDR PENALTIES, different for every scheme
        // IMGT
        let (mut query_gap_penalty, mut consensus_gap_penalty) =
            if self.cdr_positions().contains(&position) {
            if self.scheme_type == Scheme::IMGT {
                if self.cdr1.start <= position
                    && position < self.cdr1.end
                    && position != insertion_points::CDR1_IMGT
                {
                    // cdr1
                    penalty += scoring.pen_leap_insertion_point_imgt
                        + scoring.cdr_increase
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
                    penalty += scoring.pen_leap_insertion_point_imgt
                        + scoring.cdr_increase
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
                    penalty += scoring.pen_leap_insertion_point_imgt
                        + scoring.cdr_increase
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
                    penalty += scoring.pen_leap_insertion_point_kabat
                        + (scoring.cdr_increase
                            * (position as isize - cdr1_insertion_position as isize).abs() as f64);
                } else if self.cdr2.start <= position
                    && position < self.cdr2.end
                    && position != cdr2_insertion_position
                {
                    // cdr2
                    // Exception in heavy scheme:
                    if self.chain_type == Chain::IGH && position > 54 {
                        penalty = scoring.gap_pen_fr;
                    } else {
                        penalty += scoring.pen_leap_insertion_point_kabat
                            + (scoring.cdr_increase
                                * (position as isize - cdr2_insertion_position as isize).abs()
                                    as f64);
                    }
                } else if self.cdr3.start <= position
                    && position < self.cdr3.end
                    && position != cdr3_insertion_position
                {
                    // cdr3
                    penalty += scoring.pen_leap_insertion_point_kabat
                        + (scoring.cdr_increase
                            * (position as isize - cdr3_insertion_position as isize).abs() as f64);
                    if position > cdr3_insertion_position {
                        // higher penalty after insertion in cdr3
                        penalty += 5.0 * scoring.cdr_increase;
                    }
                }

                if self.chain_type != Chain::IGH && position == cdr1_insertion_position {
                    penalty = scoring.gap_pen_other;
                }
            }
            (scoring.gap_pen_cp, penalty) //TODO maybe make seperate variable, for now just high
        }else {(penalty,penalty)};


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
            query_gap_penalty = scoring.gap_pen_ip;
        }

        if self.gap_positions.contains(&position) {
            consensus_gap_penalty = scoring.gap_pen_op;
            // TODO set gap penalty to other?
            // query_gap_penalty = GAP_PEN_OTHER;
        }

        (query_gap_penalty, consensus_gap_penalty)
    }
    /// numbers sequence, returns NumberingOutput
    pub(crate) fn number_sequence<'a>(&'a self, query_sequence: &'a [u8]) -> NumberingOutput<'a> {
        let (mut numbering, identity) = 
            needleman_wunsch_consensus(query_sequence, self);

        // give gap positions correct names as defined by the numbering scheme
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
    pub fn to_terminal_schemes(&self, terminal_length: u8) -> (NumberingScheme, NumberingScheme) {
        // Create N-terminal scheme (positions 1 to terminal_length)
        let mut n_terminal_consensus = HashMap::new();
        for i in 1..=terminal_length as u32 {
            if let Some(amino_acids) = self.consensus_amino_acids.get(&i) {
                n_terminal_consensus.insert(i, amino_acids.clone());
            }
        }

        let n_terminal_scheme = NumberingScheme {
            conserved_positions: vec![], // No conserved positions at N-terminal
            insertion_positions: vec![],
            gap_positions: vec![10], // IMGT Position 10
            consensus_amino_acids: n_terminal_consensus,
            scoring_matrix: self
                .scoring_matrix
                .slice(s![0..terminal_length as usize, ..])
                .to_owned(), // TODO wrong
            ..self.clone() // Use the remaining fields from the original
        };

        // Create C-terminal scheme
        let c_term_start: u32 = self.fr4.start;
        let mut c_terminal_consensus = HashMap::new();
        for i in 0..terminal_length as u32 {
            let original_position = c_term_start + i;
            if let Some(amino_acids) = self.consensus_amino_acids.get(&original_position) {
                c_terminal_consensus.insert(i + 1, amino_acids.clone());
            }
        }

        let fmwk4_start = self.fr4.start - 1;

        let c_terminal_scheme = NumberingScheme {
            conserved_positions: vec![1, 2, 4], // IMGT position 118, 119 and 121 (mapped to 1, 2, 4)
            insertion_positions: vec![],
            gap_positions: vec![], // No gap positions at c-terminal
            consensus_amino_acids: c_terminal_consensus,
            scoring_matrix: self
                .scoring_matrix
                .slice(s![
                    fmwk4_start as usize..(fmwk4_start + terminal_length as u32) as usize,
                    ..
                ])
                .to_owned(),
            ..self.clone() // Use the remaining fields from the original
        };

        (n_terminal_scheme, c_terminal_scheme)
    }
}
