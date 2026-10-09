//! Core types for sequence numbering

use crate::error::{Error, Result};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;
use strum::{EnumMessage, IntoEnumIterator};
use strum_macros::{Display, EnumIter, EnumMessage, EnumString};

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// A chain type. Parses case-insensitively from its locus (`IGH`), letter (`H`) or name (`heavy`).
#[cfg_attr(feature = "python", pyclass(get_all))]
#[derive(
    Debug, EnumString, EnumMessage, Display, PartialEq, Serialize, Deserialize, Clone, Copy,
)]
#[strum(parse_err_ty = Error, parse_err_fn = unknown_chain)]
pub enum Chain {
    #[strum(
        serialize = "IGH",
        to_string = "H",
        serialize = "heavy",
        ascii_case_insensitive
    )]
    IGH,
    #[strum(
        serialize = "IGK",
        to_string = "K",
        serialize = "kappa",
        ascii_case_insensitive
    )]
    IGK,
    #[strum(
        serialize = "IGL",
        to_string = "L",
        serialize = "lambda",
        ascii_case_insensitive
    )]
    IGL,
    #[strum(
        serialize = "TRA",
        to_string = "A",
        serialize = "alpha",
        ascii_case_insensitive
    )]
    TRA,
    #[strum(
        serialize = "TRB",
        to_string = "B",
        serialize = "beta",
        ascii_case_insensitive
    )]
    TRB,
    #[strum(
        serialize = "TRG",
        to_string = "G",
        serialize = "gamma",
        ascii_case_insensitive
    )]
    TRG,
    #[strum(
        serialize = "TRD",
        to_string = "D",
        serialize = "delta",
        ascii_case_insensitive
    )]
    TRD,
}

/// All chain variants
pub const ALL_CHAINS: &[Chain] = &[
    Chain::IGH,
    Chain::IGK,
    Chain::IGL,
    Chain::TRA,
    Chain::TRB,
    Chain::TRG,
    Chain::TRD,
];

/// All immunoglobulin chains
pub const IG_CHAINS: &[Chain] = &[Chain::IGH, Chain::IGK, Chain::IGL];

/// All T-cell receptor chains
pub const TCR_CHAINS: &[Chain] = &[Chain::TRA, Chain::TRB, Chain::TRG, Chain::TRD];

/// Groups of chains that [`Chain::parse_names`] accepts besides single chains
const CHAIN_GROUPS: [(&str, &[Chain]); 3] =
    [("all", ALL_CHAINS), ("ig", IG_CHAINS), ("tcr", TCR_CHAINS)];

fn unknown_chain(name: &str) -> Error {
    unknown_chain_among(name, &[])
}

// An unknown chain name, with every name that would have been accepted: each chain's names as
// `Chain` parses them, and `groups` where groups are accepted too
fn unknown_chain_among(name: &str, groups: &[(&str, &[Chain])]) -> Error {
    let mut options: Vec<String> = ALL_CHAINS.iter().map(accepted_names).collect();
    if !groups.is_empty() {
        let groups: Vec<&str> = groups.iter().map(|(group, _)| *group).collect();
        options.push(format!("or the groups {}", groups.join(", ")));
    }
    Error::InvalidChain(format!(
        "unknown chain '{name}' (options: {})",
        options.join(", ")
    ))
}

impl Chain {
    /// Parse chain names, each a chain (see [`Chain`]) or a group of chains: `ig`, `tcr` or `all`.
    /// Case-insensitive and trimmed. A chain named more than once is kept once, where first named.
    pub fn parse_names<'a>(names: impl IntoIterator<Item = &'a str>) -> Result<Vec<Chain>> {
        let mut chains = Vec::new();
        for name in names {
            let name = name.trim();
            let group = CHAIN_GROUPS
                .into_iter()
                .find(|(group, _)| name.eq_ignore_ascii_case(group));
            let named = match group {
                Some((_, group)) => group,
                None => &[name
                    .parse::<Chain>()
                    .map_err(|_| unknown_chain_among(name, &CHAIN_GROUPS))?][..],
            };
            for &chain in named {
                if !chains.contains(&chain) {
                    chains.push(chain);
                }
            }
        }
        Ok(chains)
    }
}

/// Numbering schemes for output. Parses case-insensitively from its name (`kabat`) or initial (`k`).
#[cfg_attr(feature = "python", pyclass(get_all))]
#[derive(
    Debug,
    EnumString,
    EnumMessage,
    EnumIter,
    Display,
    PartialEq,
    Serialize,
    Deserialize,
    Clone,
    Copy,
)]
#[strum(parse_err_ty = Error, parse_err_fn = unknown_scheme)]
pub enum Scheme {
    /// IMGT numbering (canonical internal representation)
    #[strum(to_string = "IMGT", serialize = "i", ascii_case_insensitive)]
    IMGT,
    /// Kabat numbering (derived from IMGT)
    #[strum(to_string = "Kabat", serialize = "k", ascii_case_insensitive)]
    Kabat,
    /// Chothia numbering (derived from IMGT)
    #[strum(to_string = "Chothia", serialize = "c", ascii_case_insensitive)]
    Chothia,
    /// Martin / extended Chothia numbering (derived from IMGT)
    #[strum(to_string = "Martin", serialize = "m", ascii_case_insensitive)]
    Martin,
    /// AHo numbering (derived from IMGT)
    #[strum(to_string = "Aho", serialize = "a", ascii_case_insensitive)]
    Aho,
}

fn unknown_scheme(name: &str) -> Error {
    let options: Vec<String> = Scheme::iter()
        .map(|scheme| accepted_names(&scheme))
        .collect();
    Error::InvalidScheme(format!(
        "unknown scheme '{name}' (options: {})",
        options.join(", ")
    ))
}

// Every name a variant parses from, as strum matches them, joined by `/` for an error message: its
// display name first, as results report it, then the rest
fn accepted_names<T: EnumMessage + fmt::Display>(variant: &T) -> String {
    let display = variant.to_string();
    let others = variant
        .get_serializations()
        .iter()
        .filter(|&&name| name != display);
    std::iter::once(display.as_str())
        .chain(others.copied())
        .collect::<Vec<_>>()
        .join("/")
}

impl Scheme {
    /// Whether this scheme numbers `chain`. IMGT numbers every chain; Kabat, Chothia, Martin and
    /// AHo rules are derived for antibody chains only. AHo is defined for TCR chains too, but
    /// immunum does not ship TCR AHo rules yet.
    pub fn supports(self, chain: Chain) -> bool {
        self == Scheme::IMGT || !TCR_CHAINS.contains(&chain)
    }

    /// An error unless this scheme [`supports`](Self::supports) `chain`.
    pub fn validate_chain(self, chain: Chain) -> Result<()> {
        if !self.supports(chain) {
            return Err(Error::InvalidScheme(format!(
                "{self} scheme only supported for antibody chains (IGH, IGK, IGL)"
            )));
        }
        Ok(())
    }
}

/// Whether `scheme` numbers `chain` (see [`Scheme::supports`]), both by the names users write:
/// see [`Scheme`] and [`Chain`]. An error for an unknown name.
pub fn scheme_supports_chain(scheme: &str, chain: &str) -> Result<bool> {
    Ok(scheme.parse::<Scheme>()?.supports(chain.parse()?))
}

/// Position in a numbered sequence
/// Can be a simple number or a number with an insertion letter (e.g., "111A")
#[cfg_attr(feature = "python", pyclass(get_all))]
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, EnumString, Display)]
pub enum Region {
    FR1,
    CDR1,
    FR2,
    CDR2,
    FR3,
    CDR3,
    FR4,
}

impl Region {
    /// The region's name as every interface writes it: `fr1`, `cdr1`, ... `fr4`.
    /// [`SEGMENT_NAMES`](crate::numbering::SEGMENT_NAMES) takes its region names from here.
    pub const fn name(self) -> &'static str {
        match self {
            Region::FR1 => "fr1",
            Region::CDR1 => "cdr1",
            Region::FR2 => "fr2",
            Region::CDR2 => "cdr2",
            Region::FR3 => "fr3",
            Region::CDR3 => "cdr3",
            Region::FR4 => "fr4",
        }
    }
}

/// Region definition for a numbering scheme.
///
/// The seven regions are contiguous starting at position 1, so a scheme's region layout is
/// fully described by the last position number of each region: FR1 = `1..=fr1_end`,
/// CDR1 = `fr1_end+1..=cdr1_end`, and so on. Positions of 0 (prefix) or beyond `fr4_end`
/// (postfix) are outside the numbered range. Each scheme defines its own in its rule module.
#[derive(Debug, Clone, Copy)]
pub struct RegionDefinition {
    pub fr1_end: u8,
    pub cdr1_end: u8,
    pub fr2_end: u8,
    pub cdr2_end: u8,
    pub fr3_end: u8,
    pub cdr3_end: u8,
    pub fr4_end: u8,
}

impl RegionDefinition {
    /// Region for a position number, or `None` if outside the numbered range.
    pub const fn region(&self, pos: u8) -> Option<Region> {
        if pos == 0 {
            None
        } else if pos <= self.fr1_end {
            Some(Region::FR1)
        } else if pos <= self.cdr1_end {
            Some(Region::CDR1)
        } else if pos <= self.fr2_end {
            Some(Region::FR2)
        } else if pos <= self.cdr2_end {
            Some(Region::CDR2)
        } else if pos <= self.fr3_end {
            Some(Region::FR3)
        } else if pos <= self.cdr3_end {
            Some(Region::CDR3)
        } else if pos <= self.fr4_end {
            Some(Region::FR4)
        } else {
            None
        }
    }

    /// The seven regions as inclusive `(start, end)` position pairs, N- to C-terminal.
    pub const fn spans(&self) -> [(Region, (u8, u8)); 7] {
        [
            (Region::FR1, (1, self.fr1_end)),
            (Region::CDR1, (self.fr1_end + 1, self.cdr1_end)),
            (Region::FR2, (self.cdr1_end + 1, self.fr2_end)),
            (Region::CDR2, (self.fr2_end + 1, self.cdr2_end)),
            (Region::FR3, (self.cdr2_end + 1, self.fr3_end)),
            (Region::CDR3, (self.fr3_end + 1, self.cdr3_end)),
            (Region::FR4, (self.cdr3_end + 1, self.fr4_end)),
        ]
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
        assert_eq!("igh".parse::<Chain>().unwrap(), Chain::IGH);
        assert_eq!("H".parse::<Chain>().unwrap(), Chain::IGH);
        assert_eq!("heavy".parse::<Chain>().unwrap(), Chain::IGH);
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
    fn parse_names_expands_each_group() {
        let ig = Chain::parse_names(["ig"]).unwrap();
        assert_eq!(ig, vec![Chain::IGH, Chain::IGK, Chain::IGL]);

        let tcr = Chain::parse_names(["tcr"]).unwrap();
        assert_eq!(tcr, vec![Chain::TRA, Chain::TRB, Chain::TRG, Chain::TRD]);

        let all = Chain::parse_names(["all"]).unwrap();
        assert_eq!(all, ALL_CHAINS);
    }

    // Groups were only accepted by the CLI; every surface now parses names through `parse_names`.
    #[test]
    fn parse_names_mixes_chains_and_groups_and_drops_repeats() {
        let chains = Chain::parse_names([" Heavy", "ig", "TCR", "b"]).unwrap();
        assert_eq!(
            chains,
            vec![
                Chain::IGH,
                Chain::IGK,
                Chain::IGL,
                Chain::TRA,
                Chain::TRB,
                Chain::TRG,
                Chain::TRD
            ]
        );
        assert!(Chain::parse_names(["ig", "IGX"]).is_err());
    }

    #[test]
    fn an_unknown_chain_lists_the_names_that_would_parse() {
        let single = "IGX".parse::<Chain>().unwrap_err().to_string();
        let listed = Chain::parse_names(["IGX"]).unwrap_err().to_string();
        for chain in ALL_CHAINS {
            for name in chain.get_serializations() {
                assert_eq!(name.parse::<Chain>().unwrap(), *chain);
                assert!(
                    single.contains(name) && listed.contains(name),
                    "{name} not listed"
                );
            }
        }
        // Groups are offered only where they're accepted
        assert!(listed.ends_with("or the groups all, ig, tcr)"), "{listed}");
        assert!(!single.contains("groups"), "{single}");
    }

    #[test]
    fn an_unknown_scheme_lists_the_names_that_would_parse() {
        let message = "Z".parse::<Scheme>().unwrap_err().to_string();
        for scheme in Scheme::iter() {
            for name in scheme.get_serializations() {
                assert_eq!(name.parse::<Scheme>().unwrap(), scheme);
                assert!(
                    message.contains(&format!("{scheme}/")),
                    "{scheme} not listed first"
                );
                assert!(message.contains(name), "{name} not listed");
            }
        }
    }

    /// A definition stores only the region ends, so the starts are arithmetic: every start is the
    /// previous end plus one, and FR1 starts at 1. Values are the IMGT table.
    #[test]
    fn spans_reconstruct_starts_from_ends() {
        let imgt = RegionDefinition {
            fr1_end: 26,
            cdr1_end: 38,
            fr2_end: 55,
            cdr2_end: 65,
            fr3_end: 104,
            cdr3_end: 117,
            fr4_end: 128,
        };

        assert_eq!(
            imgt.spans(),
            [
                (Region::FR1, (1, 26)),
                (Region::CDR1, (27, 38)),
                (Region::FR2, (39, 55)),
                (Region::CDR2, (56, 65)),
                (Region::FR3, (66, 104)),
                (Region::CDR3, (105, 117)),
                (Region::FR4, (118, 128)),
            ]
        );
    }

    #[test]
    fn imgt_numbers_every_chain() {
        for &chain in ALL_CHAINS {
            assert!(
                Scheme::IMGT.validate_chain(chain).is_ok(),
                "IMGT should number {chain}"
            );
        }
    }

    #[test]
    fn schemes_without_tcr_rules_reject_tcr_chains() {
        for scheme in [Scheme::Kabat, Scheme::Chothia, Scheme::Martin, Scheme::Aho] {
            for &chain in IG_CHAINS {
                assert!(
                    scheme.validate_chain(chain).is_ok(),
                    "{scheme} should number {chain}"
                );
            }
            for &chain in TCR_CHAINS {
                assert!(
                    scheme.validate_chain(chain).is_err(),
                    "{scheme} has no rules for {chain}"
                );
            }
        }
    }
}
