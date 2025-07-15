#![allow(clippy::upper_case_acronyms)]

use clap::ValueEnum;
use immunum_macros::ParseFromString;
use ndarray::Array2;
use std::collections::HashMap;
use std::ops::Range;
#[derive(Debug)]
pub struct RegionRange {
    pub start: u32,
    pub end: u32,
}
impl RegionRange {
    pub fn positions(&self) -> Range<u32> {
        self.start..self.end
    }
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
    // TODO scoring matrix // framework_positions: //cdr_positions:
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
