#![allow(clippy::upper_case_acronyms)]

use clap::ValueEnum;
use std::ops::Range;
use std::str::FromStr;

#[derive(Debug, Clone)]
pub struct RegionRange {
    pub start: u32,
    pub end: u32,
}
impl RegionRange {
    pub fn positions(&self) -> Range<u32> {
        self.start..self.end
    }
}

/// Struct for output of prefiltering, containing identity of terminals, start and end
pub struct PrefilterOutput {
    pub identity: f64,
    pub predicted_start: u32,
    pub predicted_end: u32,
}

/// Numbering schemes for immunoglobulin sequences
#[derive(Clone, Debug, PartialEq, ValueEnum)]
pub enum Scheme {
    /// IMGT numbering scheme
    #[value(alias = "I")]
    IMGT,
    /// Kabat numbering scheme
    #[value(alias = "K")]
    KABAT,
}

impl FromStr for Scheme {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "imgt" | "i" => Ok(Scheme::IMGT),
            "kabat" | "k" => Ok(Scheme::KABAT),
            _ => Err(format!(
                "Scheme not supported: '{}', use any of: IMGT (I), KABAT (K)",
                s
            )),
        }
    }
}


/// Immunoglobulin and T-cell receptor chain types
#[derive(Clone, Copy, Debug, PartialEq, Hash, Eq, ValueEnum)]
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

impl FromStr for Chain {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "igh" | "heavy" | "h" => Ok(Chain::IGH),
            "igk" | "kappa" | "k" => Ok(Chain::IGK),
            "igl" | "lambda" | "l" => Ok(Chain::IGL),
            "tra" | "alpha" | "a" => Ok(Chain::TRA),
            "trb" | "beta" | "b" => Ok(Chain::TRB),
            "trg" | "gamma" | "g" => Ok(Chain::TRG),
            "trd" | "delta" | "d" => Ok(Chain::TRD),
            _ => Err(format!(
                "Chain not supported: '{}', use any of: IGH (Heavy, H), IGK (Kappa, K), IGL (Lambda, L), TRA (Alpha, A), TRB (Beta, B), TRG (Gamma, G), TRD (Delta, D)",
                s
            )),
        }
    }
}

