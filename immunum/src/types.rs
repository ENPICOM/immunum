#![allow(clippy::upper_case_acronyms)]

use std::collections::HashMap;
use clap::ValueEnum;
use immunum_macros::ParseFromString;
use std::ops::Range;
use crate::constants::{CDR1_INSERTION_POSITION_IMGT, CDR1_INSERTION_POSITION_KABAT_HEAVY, CDR1_INSERTION_POSITION_KABAT_LIGHT, CDR2_INSERTION_POSITION_IMGT, CDR2_INSERTION_POSITION_KABAT_HEAVY, CDR2_INSERTION_POSITION_KABAT_LIGHT, CDR3_INSERTION_POSITION_IMGT, CDR3_INSERTION_POSITION_KABAT_HEAVY, CDR3_INSERTION_POSITION_KABAT_LIGHT, CDR_INCREASE, GAP_PEN_CDR, GAP_PEN_CP, GAP_PEN_FR, GAP_PEN_IP, GAP_PEN_OP, GAP_PEN_OTHER, PEN_LEAP_FROM_INSERTION_POINT_IMGT, PEN_LEAP_INSERTION_POINT_KABAT};
#[derive(Debug)]
pub struct RegionRange {
    pub(crate) start: u32,
    pub(crate) end: u32

}
impl RegionRange {
    pub fn positions(&self) -> Range<u32> {
        self.start..self.end
    }
}

#[derive(Debug)]
pub struct NumberingOutput<'a>{
    pub(crate) scheme: &'a NumberingScheme,
    pub(crate) sequence: String,
    pub(crate) numbering: Vec<String>,
    pub(crate) identity: f64,
    pub(crate) score: f64,
    pub(crate) start: i32,
    pub(crate) end: i32,
}

#[derive(Debug)]
pub struct NumberingScheme {
    pub(crate) name: String,
    pub(crate) description: String,
    pub(crate) scheme_type: Scheme,
    pub(crate) chain_type: Chain,
    pub(crate) conserved_positions: Vec<u32>,
    pub(crate) insertion_positions: Vec<u32>,
    pub(crate) gap_positions: Vec<u32>,
    pub(crate) consensus_amino_acids: HashMap<u32, Vec<char>>,
    pub(crate) fr1: RegionRange,
    pub(crate) fr2: RegionRange,
    pub(crate) fr3: RegionRange,
    pub(crate) fr4: RegionRange,
    pub(crate) cdr1: RegionRange,
    pub(crate) cdr2: RegionRange,
    pub(crate) cdr3: RegionRange,
    // TODO scoring matrix // framework_positions: //cdr_positions:
}

impl NumberingScheme {
    pub(crate) fn restricted_sites(&self) -> Vec<u32>{
        let mut sites = Vec::new();
        for (&key, value) in &self.consensus_amino_acids{
            if !value.contains(&'-'){
                sites.push(key);
            }
        }
        sites
    }

    pub(crate) fn framework_positions(&self) -> Vec<u32>{
        // TODO
        self.fr1.positions().chain(
            self.fr2.positions()).chain(
            self.fr3.positions()).chain(
            self.fr4.positions()).collect()
    }

    pub(crate) fn cdr_positions(&self) -> Vec<u32>{
        //TODO
        self.cdr1.positions().chain(
            self.cdr2.positions()).chain(
            self.cdr3.positions()).collect()
    }

    pub(crate) fn gap_penalty(&self, position: u32) -> (f64, f64) {
        // Set initial gap penalties
        let mut penalty = if self.conserved_positions.contains(&position) {
            GAP_PEN_CP
        } else if self.framework_positions().contains(&position) {
            GAP_PEN_FR
        } else if self.cdr_positions().contains(&position) {
            GAP_PEN_CDR
        } else {
            GAP_PEN_OTHER
        };

        // ADAPT CDR PENALTIES, different for every scheme
        // IMGT
        if self.cdr_positions().contains(&position) {
            if self.scheme_type == Scheme::IMGT {
                if self.cdr1.start <= position && position < self.cdr1.end
                    && position != CDR1_INSERTION_POSITION_IMGT
                {
                    // cdr1
                    penalty += PEN_LEAP_FROM_INSERTION_POINT_IMGT + CDR_INCREASE * (position as isize - CDR1_INSERTION_POSITION_IMGT as isize).abs() as f64;
                    penalty += if position > CDR1_INSERTION_POSITION_IMGT { 0.1 } else { 0.0 };
                } else if self.cdr2.start <= position && position < self.cdr2.end
                    && position != CDR2_INSERTION_POSITION_IMGT
                {
                    // cdr2
                    penalty += PEN_LEAP_FROM_INSERTION_POINT_IMGT + CDR_INCREASE * (position as isize - CDR2_INSERTION_POSITION_IMGT as isize).abs() as f64;
                    penalty += if position > CDR2_INSERTION_POSITION_IMGT { 0.1 } else { 0.0 };
                } else if self.cdr3.start <= position && position < self.cdr3.end
                    && position != CDR3_INSERTION_POSITION_IMGT
                {
                    // cdr3
                    penalty += PEN_LEAP_FROM_INSERTION_POINT_IMGT + CDR_INCREASE * (position as isize - CDR3_INSERTION_POSITION_IMGT as isize).abs() as f64;
                    penalty += if position < CDR3_INSERTION_POSITION_IMGT { 0.1 } else { 0.0 };
                }
            }
            // KABAT
            else if self.scheme_type == Scheme::KABAT {
                let cdr1_insertion_position = if self.chain_type == Chain::IGH {
                    CDR1_INSERTION_POSITION_KABAT_HEAVY
                } else {
                    CDR1_INSERTION_POSITION_KABAT_LIGHT
                };

                let cdr2_insertion_position = if self.chain_type == Chain::IGH {
                    CDR2_INSERTION_POSITION_KABAT_HEAVY
                } else {
                    CDR2_INSERTION_POSITION_KABAT_LIGHT
                };

                let cdr3_insertion_position = if self.chain_type == Chain::IGH {
                    CDR3_INSERTION_POSITION_KABAT_HEAVY
                } else {
                    CDR3_INSERTION_POSITION_KABAT_LIGHT
                };

                if self.cdr1.start <= position && position < self.cdr1.end
                    && position != cdr1_insertion_position
                {
                    // cdr1
                    penalty += PEN_LEAP_INSERTION_POINT_KABAT +
                        (CDR_INCREASE * (position as isize - cdr1_insertion_position as isize).abs() as f64);
                } else if self.cdr2.start <= position && position < self.cdr2.end
                    && position != cdr2_insertion_position
                {
                    // cdr2
                    // Exception in heavy scheme:
                    if self.chain_type == Chain::IGH && position > 54 {
                        penalty = GAP_PEN_FR;
                    } else {
                        penalty += PEN_LEAP_INSERTION_POINT_KABAT +
                            (CDR_INCREASE * (position as isize - cdr2_insertion_position as isize).abs() as f64);
                    }
                } else if self.cdr3.start <= position && position < self.cdr3.end
                    && position != cdr3_insertion_position
                {
                    // cdr3
                    penalty += PEN_LEAP_INSERTION_POINT_KABAT +
                        (CDR_INCREASE * (position as isize - cdr3_insertion_position as isize).abs() as f64);
                    if position > cdr3_insertion_position {
                        // higher penalty after insertion in cdr3
                        penalty += 5.0;
                    }
                }

                if self.chain_type != Chain::IGH && position == cdr1_insertion_position {
                    penalty = GAP_PEN_OTHER; // 5, 11 in antpack, here set to 11
                }
            }
        }

        let mut query_gap_penalty = penalty;
        let mut consensus_gap_penalty = penalty;

        // Handle start penalty, same for all schemes
        if position < 18 {
            // Increase from 2 to 11, then add 0.1 until position 18
            consensus_gap_penalty = 1.0 +
                (if position > 10 { 10.0 } else { position as f64 }) +
                (if position > 10 { 0.1 * (position as f64 - 10.0) } else { 0.0 });
        }

        // Adapt only query or consensus gap for insertion and gap positions
        if self.insertion_positions.contains(&position) {
            query_gap_penalty = GAP_PEN_IP;
        }

        if self.gap_positions.contains(&position) {
            consensus_gap_penalty = GAP_PEN_OP;
            // TODO set gap penalty to other?
            // query_gap_penalty = GAP_PEN_OTHER;
        }

        (query_gap_penalty, consensus_gap_penalty)
    }
}

/// Numbering schemes for immunoglobulin sequences
#[derive(Clone, Debug, PartialEq, ValueEnum, ParseFromString)]
pub enum Scheme {
    /// IMGT numbering scheme
    #[value(alias = "I")]
    IMGT,
    /// Kabat numbering scheme
    #[value(alias = "K")]
    KABAT,
}

/// Immunoglobulin and T-cell receptor chain types
#[derive(Clone, Debug, PartialEq, ValueEnum, ParseFromString)]
pub enum Chain {
    // IG Heavy chain variants
    #[value(alias = "Heavy", alias = "H")]
    IGH,

    // IG Kappa chain variants
    #[value(alias = "Kappa", alias = "K")]
    IGK,

    // IG Lambda chain variants
    #[value(alias = "Lambda", alias = "L")]
    IGL,

    // T-cell receptor Alpha chain variants
    #[value(alias = "Alpha", alias = "A")]
    TRA,

    // T-cell receptor Beta chain variants
    #[value(alias = "Beta", alias = "B")]
    TRB,

    // T-cell receptor Gamma chain variants
    #[value(alias = "Gamma", alias = "G")]
    TRG,

    // T-cell receptor Delta chain variants
    #[value(alias = "Delta", alias = "D")]
    TRD,
}
