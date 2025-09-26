#![allow(clippy::upper_case_acronyms)]

use clap::ValueEnum;
use serde::Serialize;
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
#[derive(Debug)]
pub struct PrefilterOutput {
    pub identity: f64,
    pub _predicted_start: usize,
    pub _predicted_end: usize,
}

/// Numbering schemes for immunoglobulin sequences
#[derive(Clone, Copy, Debug, PartialEq, ValueEnum, Serialize)]
#[cfg_attr(feature = "python", pyo3::pyclass(eq, eq_int))]
#[cfg_attr(feature = "wasm", wasm_bindgen::prelude::wasm_bindgen)]
pub enum Scheme {
    /// IMGT numbering scheme
    #[value(alias = "I")]
    IMGT,
    /// Kabat numbering scheme
    #[value(alias = "K")]
    KABAT,
}

/// CDR definition schemes for determining CDR boundaries
#[derive(Clone, Copy, Debug, PartialEq, ValueEnum, Serialize)]
#[cfg_attr(feature = "python", pyo3::pyclass(eq, eq_int))]
#[cfg_attr(feature = "wasm", wasm_bindgen::prelude::wasm_bindgen)]
pub enum CdrDefinitions {
    /// IMGT CDR definitions (default for IMGT scheme)
    #[value(alias = "I")]
    IMGT,
    /// Kabat CDR definitions (default for KABAT scheme)
    #[value(alias = "K")]
    KABAT,
    /// North CDR definitions
    #[value(alias = "N")]
    NORTH,
}

impl CdrDefinitions {
    /// Get the default CDR definition for a given numbering scheme
    pub fn from_scheme(scheme: Scheme) -> Self {
        match scheme {
            Scheme::IMGT => CdrDefinitions::IMGT,
            Scheme::KABAT => CdrDefinitions::KABAT,
        }
    }
}

/// Immunoglobulin and T-cell receptor chain types
#[derive(Clone, Copy, Debug, PartialEq, Hash, Eq, ValueEnum, Serialize)]
#[cfg_attr(feature = "python", pyo3::pyclass(eq, eq_int))]
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

/// Efficient representation of a numbering position
#[derive(Debug, Clone, Serialize, PartialEq)]
#[serde(untagged)]
pub enum NumberingPosition {
    /// Gap in alignment
    Gap,
    /// Simple numeric position
    Number(u32),
    /// Position with insertion (e.g., "32A")
    Insertion { position: u32, insertion: char },
}

impl std::fmt::Display for NumberingPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NumberingPosition::Gap => write!(f, "-"),
            NumberingPosition::Number(n) => write!(f, "{}", n),
            NumberingPosition::Insertion {
                position,
                insertion,
            } => {
                write!(f, "{}{}", position, insertion)
            }
        }
    }
}

impl NumberingPosition {
    /// Create from string (for backward compatibility)
    pub fn from_string(s: &str) -> Self {
        if s == "-" {
            NumberingPosition::Gap
        } else if let Ok(num) = s.parse::<u32>() {
            NumberingPosition::Number(num)
        } else {
            // Try to parse as insertion (e.g., "32A")
            if s.len() > 1 {
                if let Ok(position) = s[..s.len() - 1].parse::<u32>() {
                    if let Some(insertion) = s.chars().last() {
                        return NumberingPosition::Insertion {
                            position,
                            insertion,
                        };
                    }
                }
            }
            // Fallback to treating as number 0 if can't parse
            NumberingPosition::Number(0)
        }
    }
}

impl PartialEq<str> for NumberingPosition {
    fn eq(&self, other: &str) -> bool {
        match self {
            NumberingPosition::Gap => other == "-",
            NumberingPosition::Number(n) => other == n.to_string(),
            NumberingPosition::Insertion {
                position,
                insertion,
            } => other == format!("{}{}", position, insertion),
        }
    }
}

impl PartialEq<&str> for NumberingPosition {
    fn eq(&self, other: &&str) -> bool {
        self == *other
    }
}

impl PartialEq<String> for NumberingPosition {
    fn eq(&self, other: &String) -> bool {
        self == other.as_str()
    }
}

impl PartialEq<NumberingPosition> for String {
    fn eq(&self, other: &NumberingPosition) -> bool {
        other == self.as_str()
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct RegionInfo {
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct Regions {
    pub fr1: RegionInfo,
    pub cdr1: RegionInfo,
    pub fr2: RegionInfo,
    pub cdr2: RegionInfo,
    pub fr3: RegionInfo,
    pub cdr3: RegionInfo,
    pub fr4: RegionInfo,
}

#[derive(Debug, Clone, Serialize)]
pub struct ChainNumbering {
    #[serde(serialize_with = "serialize_numbering_positions")]
    pub numbers: Vec<NumberingPosition>,
    pub identity: f64,
    pub scheme: Scheme,
    pub chain: Chain,
    pub cdr_definition: CdrDefinitions,
    pub start: usize,
    pub end: usize,
    pub regions: Regions,
}

/// Custom serializer to maintain backward compatibility with JSON output
fn serialize_numbering_positions<S>(
    positions: &[NumberingPosition],
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    let strings: Vec<String> = positions.iter().map(|p| p.to_string()).collect();
    strings.serialize(serializer)
}

#[derive(Debug, Clone, Serialize)]
pub struct SequenceResult {
    pub sequence_id: String,
    pub chains: Vec<ChainNumbering>,
}

impl SequenceResult {
    pub fn new(sequence_id: String, chains: Vec<ChainNumbering>) -> Self {
        Self {
            sequence_id,
            chains,
        }
    }
}
