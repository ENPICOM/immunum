#![allow(clippy::upper_case_acronyms)]

use clap::ValueEnum;
use std::ops::Range;

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
    pub _predicted_start: usize,
    pub _predicted_end: usize,
}

/// Numbering schemes for immunoglobulin sequences
#[derive(Clone, Copy, Debug, PartialEq, ValueEnum)]
#[cfg_attr(feature = "python", pyo3::pyclass)]
#[cfg_attr(feature = "wasm", wasm_bindgen::prelude::wasm_bindgen)]
pub enum Scheme {
    /// IMGT numbering scheme
    #[value(alias = "I")]
    IMGT,
    /// Kabat numbering scheme
    #[value(alias = "K")]
    KABAT,
}

/// Immunoglobulin and T-cell receptor chain types
#[derive(Clone, Copy, Debug, PartialEq, Hash, Eq, ValueEnum)]
#[cfg_attr(feature = "python", pyo3::pyclass)]
#[cfg_attr(feature = "wasm", wasm_bindgen::prelude::wasm_bindgen)]
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

impl Chain {
    pub fn to_short(self) -> &'static str {
        match self {
            Chain::IGH => "H",
            Chain::IGK => "K",
            Chain::IGL => "L",
            Chain::TRA => "A",
            Chain::TRB => "B",
            Chain::TRG => "G",
            Chain::TRD => "D",
        }
    }
}

#[derive(Debug)]
pub struct ChainNumbering {
    pub numbers: Vec<String>,
    pub identity: f64,
    pub scheme: Scheme,
    pub chain: Chain,
    pub start: usize,
    pub end: usize,
}

impl ChainNumbering {
    pub fn to_json_string(&self, name: String) -> String {
        format!("{{\"name\": \"{}\", \"numbers\": {:?}, \"identity\": {}, \"scheme\": {:?}, \"chain\": {:?}, \"start\": {}, \"end\": {}}}", name, self.numbers, self.identity, self.scheme, self.chain, self.start, self.end)
    }
}
