#![allow(clippy::upper_case_acronyms)]

use std::collections::HashMap;
use clap::ValueEnum;
use immunum_macros::ParseFromString;
use std::ops::Range;

#[derive(Debug)]
pub struct RegionRange {
    pub(crate) start: i32,
    pub(crate) end: i32

}
impl RegionRange {
    pub fn positions(&self) -> Range<i32> {
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
    pub(crate) conserved_positions: Vec<i32>,
    pub(crate) insertion_positions: Vec<i32>,
    pub(crate) gap_positions: Vec<i32>,
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
