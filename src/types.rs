//! Core types for sequence numbering

use crate::error::{Error, Result};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

/// Chain types for immunoglobulins and T-cell receptors
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Chain {
    /// Immunoglobulin Heavy chain
    IGH,
    /// Immunoglobulin Kappa light chain
    IGK,
    /// Immunoglobulin Lambda light chain
    IGL,
    /// T-cell receptor Alpha chain
    TRA,
    /// T-cell receptor Beta chain
    TRB,
    /// T-cell receptor Gamma chain
    TRG,
    /// T-cell receptor Delta chain
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

impl fmt::Display for Chain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl FromStr for Chain {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        match s.to_uppercase().as_str() {
            "IGH" | "H" => Ok(Chain::IGH),
            "IGK" | "K" => Ok(Chain::IGK),
            "IGL" | "L" => Ok(Chain::IGL),
            "TRA" | "A" => Ok(Chain::TRA),
            "TRB" | "B" => Ok(Chain::TRB),
            "TRG" | "G" => Ok(Chain::TRG),
            "TRD" | "D" => Ok(Chain::TRD),
            _ => Err(Error::InvalidChain(s.to_string())),
        }
    }
}

/// Numbering schemes for output
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Scheme {
    /// IMGT numbering (canonical internal representation)
    IMGT,
    /// Kabat numbering (derived from IMGT)
    Kabat,
}

impl fmt::Display for Scheme {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl FromStr for Scheme {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        match s.to_uppercase().as_str() {
            "IMGT" => Ok(Scheme::IMGT),
            "KABAT" => Ok(Scheme::Kabat),
            _ => Err(Error::InvalidScheme(s.to_string())),
        }
    }
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

/// How insertions are handled when sequence length > base positions
#[derive(Debug, Clone, Copy)]
pub enum InsertionStyle {
    /// All insertions after a single position: 100A, 100B, 100C (Kabat style)
    AfterPosition(u8),
    /// Insertions split between two positions palindromically: 111A, 112A, 111B, 112B (IMGT style)
    Palindromic { left: u8, right: u8 },
}

/// Configuration for CDR-style renumbering with insertions and deletions
#[derive(Debug, Clone, Copy)]
pub struct RenumberConfig {
    /// Base positions to use (in order)
    pub base_positions: &'static [u8],
    /// Order to delete positions when len < base (None = delete from end)
    pub deletion_order: Option<&'static [u8]>,
    /// How to handle insertions
    pub insertion_style: InsertionStyle,
}

impl RenumberConfig {
    /// Create a config for sequential (Kabat-style) numbering
    pub const fn sequential(
        base_positions: &'static [u8],
        insertion_pos: u8,
        deletion_order: Option<&'static [u8]>,
    ) -> Self {
        Self {
            base_positions,
            deletion_order,
            insertion_style: InsertionStyle::AfterPosition(insertion_pos),
        }
    }

    /// Create a config for palindromic (IMGT-style) numbering
    pub const fn palindromic(
        base_positions: &'static [u8],
        deletion_order: &'static [u8],
        insertion_left: u8,
        insertion_right: u8,
    ) -> Self {
        Self {
            base_positions,
            deletion_order: Some(deletion_order),
            insertion_style: InsertionStyle::Palindromic {
                left: insertion_left,
                right: insertion_right,
            },
        }
    }
}

/// Functional regions in a numbered sequence
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Region {
    /// Framework region 1
    FR1,
    /// Complementarity-determining region 1
    CDR1,
    /// Framework region 2
    FR2,
    /// Complementarity-determining region 2
    CDR2,
    /// Framework region 3
    FR3,
    /// Complementarity-determining region 3
    CDR3,
    /// Framework region 4
    FR4,
}

/// How to number a region in the output scheme
#[derive(Debug, Clone, Copy)]
pub enum NumberingRegionType {
    /// Simple offset: src_pos - src_start + dst_start
    Offset { src_start: u8, dst_start: u8 },
    /// Renumbering with insertions/deletions using RenumberConfig
    WithConfig(RenumberConfig),
}

/// A region definition for numbering conversion
///
/// Maps a range of source positions to output positions using either
/// simple offset arithmetic or full renumbering with insertions/deletions.
#[derive(Debug, Clone, Copy)]
pub struct NumberingRegion {
    /// Source position range (inclusive)
    pub src_range: (u8, u8),
    /// How to number this region
    pub region_type: NumberingRegionType,
}

impl NumberingRegion {
    /// Create an offset-based region (for framework regions)
    pub const fn offset(src_start: u8, src_end: u8, dst_start: u8) -> Self {
        Self {
            src_range: (src_start, src_end),
            region_type: NumberingRegionType::Offset {
                src_start,
                dst_start,
            },
        }
    }

    /// Create a region with full renumbering config (for CDR regions)
    pub const fn with_config(src_start: u8, src_end: u8, config: RenumberConfig) -> Self {
        Self {
            src_range: (src_start, src_end),
            region_type: NumberingRegionType::WithConfig(config),
        }
    }

    /// Check if a position falls within this region
    #[inline]
    pub const fn contains(&self, pos: u8) -> bool {
        pos >= self.src_range.0 && pos <= self.src_range.1
    }
}

/// Full chain numbering configuration
///
/// Contains all regions needed to number a complete chain sequence.
/// Used by both IMGT and Kabat schemes with appropriate region configs.
#[derive(Debug, Clone, Copy)]
pub struct ChainNumberingConfig {
    /// Ordered list of regions covering the full sequence
    pub regions: &'static [NumberingRegion],
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
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
