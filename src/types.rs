//! Core types for sequence numbering

use crate::error::{Error, Result};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;
use strum_macros::{Display, EnumString};

#[cfg(feature = "python")]
use pyo3::prelude::*;

#[cfg_attr(feature = "python", pyclass(get_all))]
#[derive(Debug, EnumString, Display, PartialEq, Serialize, Deserialize, Clone, Copy)]
pub enum Chain {
    #[strum(
        serialize = "IGH",
        serialize = "H",
        serialize = "heavy",
        ascii_case_insensitive
    )]
    IGH,
    #[strum(
        serialize = "IGK",
        serialize = "K",
        serialize = "kappa",
        ascii_case_insensitive
    )]
    IGK,
    #[strum(
        serialize = "IGL",
        serialize = "L",
        serialize = "lambda",
        ascii_case_insensitive
    )]
    IGL,
    #[strum(
        serialize = "TRA",
        serialize = "A",
        serialize = "alpha",
        ascii_case_insensitive
    )]
    TRA,
    #[strum(
        serialize = "TRB",
        serialize = "B",
        serialize = "beta",
        ascii_case_insensitive
    )]
    TRB,
    #[strum(
        serialize = "TRG",
        serialize = "G",
        serialize = "gamma",
        ascii_case_insensitive
    )]
    TRG,
    #[strum(
        serialize = "TRD",
        serialize = "D",
        serialize = "delta",
        ascii_case_insensitive
    )]
    TRD,
}

impl Chain {
    /// Returns true if this is an immunoglobulin chain
    pub fn is_immunoglobulin(&self) -> bool {
        matches!(self, Chain::IGH | Chain::IGK | Chain::IGL)
    }

    /// Returns true if this is a T-cell receptor chain
    pub fn is_tcr(&self) -> bool {
        matches!(self, Chain::TRA | Chain::TRB | Chain::TRG | Chain::TRD)
    }
}

/// Numbering schemes for output
#[cfg_attr(feature = "python", pyclass(get_all))]
#[derive(Debug, EnumString, Display, PartialEq, Serialize, Deserialize, Clone, Copy)]
pub enum Scheme {
    /// IMGT numbering (canonical internal representation)
    #[strum(serialize = "IMGT", ascii_case_insensitive)]
    IMGT,
    /// Kabat numbering (derived from IMGT)
    #[strum(serialize = "Kabat", ascii_case_insensitive)]
    Kabat,
}

/// Position in a numbered sequence
/// Can be a simple number or a number with an insertion letter (e.g., "111A")
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Position {
    /// The numeric part of the position (max 128 for IMGT)
    pub number: u8,
    /// Optional insertion letter (for IMGT: A, B, C, etc.)
    pub insertion: Option<char>,
}

impl Position {
    /// Create a new position with just a number
    pub fn new(number: u8) -> Self {
        Self {
            number,
            insertion: None,
        }
    }

    /// Create a new position with a number and insertion letter
    pub fn with_insertion(number: u8, insertion: char) -> Self {
        Self {
            number,
            insertion: Some(insertion),
        }
    }
}

impl fmt::Display for Position {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(ins) = self.insertion {
            write!(f, "{}{}", self.number, ins)
        } else {
            write!(f, "{}", self.number)
        }
    }
}

impl FromStr for Position {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        let s = s.trim();
        if s.is_empty() {
            return Err(Error::InvalidPosition("empty string".to_string()));
        }

        // Find where digits end
        let digit_end = s
            .chars()
            .position(|c| !c.is_ascii_digit())
            .unwrap_or(s.len());

        if digit_end == 0 {
            return Err(Error::InvalidPosition(format!("no numeric part: {}", s)));
        }

        let number: u8 = s[..digit_end]
            .parse()
            .map_err(|_| Error::InvalidPosition(format!("invalid number: {}", s)))?;

        // Parse insertion letter if present
        let insertion = match &s[digit_end..] {
            "" => None,
            rest if rest.len() == 1 && rest.chars().next().unwrap().is_alphabetic() => {
                Some(rest.chars().next().unwrap())
            }
            _ => {
                return Err(Error::InvalidPosition(format!(
                    "invalid insertion part: {}",
                    s
                )))
            }
        };

        Ok(Self { number, insertion })
    }
}

/// Functional regions in a numbered sequence
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Region {
    FR1,
    CDR1,
    FR2,
    CDR2,
    FR3,
    CDR3,
    FR4,
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// A rule mapping a range of alignment positions to numbering positions
///
/// Defines how to handle insertions and deletions when the alignment length doesn't match the numbering range.
#[derive(Debug, Clone, Copy)]
pub struct NumberingRule {
    /// First alignment position (inclusive)
    pub align_start: u8,
    /// Last alignment position (inclusive)
    pub align_end: u8,
    /// First numbering position (inclusive)
    pub num_start: u8,
    /// Last numbering position (inclusive)
    pub num_end: u8,
    /// Order to delete positions when alignment is shorter than numbering range (for variable regions)
    pub deletion_order: &'static [u8],
    /// How to handle insertions when alignment is longer than numbering range (for variable regions)
    pub insertion: Insertion,
}

impl NumberingRule {
    /// Framework-like region with direct 1:1 mapping (alignment positions equal numbering positions)
    pub const fn fr(start: u8, end: u8) -> Self {
        Self {
            align_start: start,
            align_end: end,
            num_start: start,
            num_end: end,
            deletion_order: &[],
            insertion: Insertion::None,
        }
    }

    /// Framework region with simple offset mapping (alignment positions map to numbering positions with a fixed offset)
    pub const fn offset(align_start: u8, align_end: u8, offset: i8) -> Self {
        let num_start = (align_start as i16 + offset as i16) as u8;
        Self {
            align_start,
            align_end,
            num_start,
            num_end: num_start + (align_end - align_start),
            deletion_order: &[],
            insertion: Insertion::None,
        }
    }

    /// Variable region: CDR or other variable length region with custom deletion/insertion rules
    /// and explicit align and numbering ranges
    pub const fn variable(
        align_start: u8,
        align_end: u8,
        num_start: u8,
        num_end: u8,
        deletion_order: &'static [u8],
        insertion: Insertion,
    ) -> Self {
        Self {
            align_start,
            align_end,
            num_start,
            num_end,
            deletion_order,
            insertion,
        }
    }

    /// Check if a position falls within this rule's source range
    #[inline]
    pub const fn contains(&self, pos: u8) -> bool {
        pos >= self.align_start && pos <= self.align_end
    }
}
/// How insertions are handled when a variable region exceeds its base length
#[derive(Debug, Clone, Copy)]
pub enum Insertion {
    /// Simple offset arithmetic — no insertions possible (framework regions)
    None,
    /// All insertions after a single position: 35A, 35B, 35C (Kabat style)
    Sequential(u8),
    /// Insertions split symmetrically between two positions: 111A, 112A, 111B, 112B (IMGT style)
    Symmetric { left: u8, right: u8 },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chain_parsing() {
        assert_eq!("IGH".parse::<Chain>().unwrap(), Chain::IGH);
        assert_eq!("H".parse::<Chain>().unwrap(), Chain::IGH);
        assert_eq!("igh".parse::<Chain>().unwrap(), Chain::IGH);
        assert_eq!("TRA".parse::<Chain>().unwrap(), Chain::TRA);
        assert_eq!("A".parse::<Chain>().unwrap(), Chain::TRA);
        assert!("invalid".parse::<Chain>().is_err());
    }

    #[test]
    fn test_position_parsing() {
        let pos = "111".parse::<Position>().unwrap();
        assert_eq!(pos.number, 111);
        assert_eq!(pos.insertion, None);

        let pos = "111A".parse::<Position>().unwrap();
        assert_eq!(pos.number, 111);
        assert_eq!(pos.insertion, Some('A'));

        assert!("".parse::<Position>().is_err());
        assert!("A".parse::<Position>().is_err());
        assert!("111AB".parse::<Position>().is_err());
    }

    #[test]
    fn test_position_display() {
        assert_eq!(Position::new(111).to_string(), "111");
        assert_eq!(Position::with_insertion(111, 'A').to_string(), "111A");
    }

    #[test]
    fn test_chain_types() {
        assert!(Chain::IGH.is_immunoglobulin());
        assert!(!Chain::IGH.is_tcr());
        assert!(Chain::TRA.is_tcr());
        assert!(!Chain::TRA.is_immunoglobulin());
    }
}
