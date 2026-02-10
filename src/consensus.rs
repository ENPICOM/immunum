//! Consensus sequence parsing and representation

use crate::error::{Error, Result};
use crate::types::Chain;
use std::collections::HashSet;

/// A single position in a consensus sequence
#[derive(Debug, Clone)]
pub struct ConsensusPosition {
    /// Position number (IMGT numbering)
    pub position: u32,
    /// Allowed amino acids at this position (empty for gaps)
    pub allowed_aa: Vec<char>,
    /// Whether this is a gap position (marked with '-')
    pub is_gap: bool,
}

/// Metadata about a consensus sequence
#[derive(Debug, Clone, Default)]
pub struct ConsensusMetadata {
    /// Conserved positions (critical for alignment scoring)
    pub conserved: HashSet<u32>,
    /// Insertion point positions (where CDR insertions occur)
    pub insertion_points: HashSet<u32>,
}

/// A complete consensus sequence for a chain type
#[derive(Debug, Clone)]
pub struct Consensus {
    /// The chain type this consensus represents
    pub chain: Chain,
    /// Positions in the consensus (ordered by position number)
    pub positions: Vec<ConsensusPosition>,
    /// Metadata about conserved and insertion positions
    pub metadata: ConsensusMetadata,
}

impl Consensus {
    /// Parse a consensus file from text
    pub fn from_text(chain: Chain, text: &str) -> Result<Self> {
        let mut positions = Vec::new();
        let mut metadata = ConsensusMetadata::default();
        let mut in_metadata = false;
        let mut current_section = String::new();

        for line in text.lines() {
            let line = line.trim();

            // Skip empty lines
            if line.is_empty() {
                continue;
            }

            // End of file marker
            if line == "//" {
                break;
            }

            // Check for metadata sections
            if line.starts_with('#') {
                let content = line.trim_start_matches('#').trim();

                // Check for section headers
                if content.starts_with('[') && content.ends_with(']') {
                    in_metadata = true;
                    current_section = content
                        .trim_start_matches('[')
                        .trim_end_matches(']')
                        .to_string();
                    continue;
                }

                // Parse metadata content
                if in_metadata && !content.is_empty() {
                    Self::parse_metadata_line(&mut metadata, &current_section, content)?;
                }
                continue;
            }

            // If we hit a non-comment line, we're done with metadata
            in_metadata = false;

            // Parse position line: "position,AA1,AA2,..."
            let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
            if parts.is_empty() {
                continue;
            }

            let position: u32 = parts[0].parse().map_err(|_| {
                Error::ConsensusParseError(format!("Invalid position number: {}", parts[0]))
            })?;

            // Check if this is a gap position
            let is_gap = parts.len() == 2 && parts[1] == "-";

            let allowed_aa: Vec<char> = if is_gap {
                Vec::new()
            } else {
                parts[1..]
                    .iter()
                    .filter(|s| !s.is_empty())
                    .map(|s| {
                        if s.len() != 1 {
                            return Err(Error::ConsensusParseError(format!(
                                "Invalid amino acid: {}",
                                s
                            )));
                        }
                        Ok(s.chars().next().unwrap())
                    })
                    .collect::<Result<Vec<char>>>()?
            };

            positions.push(ConsensusPosition {
                position,
                allowed_aa,
                is_gap,
            });
        }

        if positions.is_empty() {
            return Err(Error::ConsensusParseError(
                "No positions found in consensus file".to_string(),
            ));
        }

        // Sort positions by position number
        positions.sort_by_key(|p| p.position);

        Ok(Consensus {
            chain,
            positions,
            metadata,
        })
    }

    /// Parse a metadata line based on the current section
    fn parse_metadata_line(
        metadata: &mut ConsensusMetadata,
        section: &str,
        content: &str,
    ) -> Result<()> {
        let positions: Vec<u32> = content
            .split(',')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(|s| {
                s.parse::<u32>().map_err(|_| {
                    Error::ConsensusParseError(format!("Invalid position in metadata: {}", s))
                })
            })
            .collect::<Result<Vec<u32>>>()?;

        match section {
            "conserved" => {
                metadata.conserved.extend(positions);
            }
            "insertion_points" => {
                metadata.insertion_points.extend(positions);
            }
            _ => {
                // Unknown section, ignore
            }
        }

        Ok(())
    }

    /// Get a position by its number
    pub fn get_position(&self, pos: u32) -> Option<&ConsensusPosition> {
        self.positions.iter().find(|p| p.position == pos)
    }

    /// Check if a position is conserved
    pub fn is_conserved(&self, pos: u32) -> bool {
        self.metadata.conserved.contains(&pos)
    }

    /// Check if a position is an insertion point
    pub fn is_insertion_point(&self, pos: u32) -> bool {
        self.metadata.insertion_points.contains(&pos)
    }

    /// Get the length of the consensus (number of positions)
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Check if the consensus is empty
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple_consensus() {
        let text = r#"
# Universal consensus for IGH
#
# [conserved]
# 23,41,104,118
#
# [insertion_points]
# 111,112
#
# [positions]
1,Q,E,D
2,V,Q
23,C
27,-
41,W
104,C
111,-
118,W
//
"#;

        let consensus = Consensus::from_text(Chain::IGH, text).unwrap();
        assert_eq!(consensus.chain, Chain::IGH);
        assert_eq!(consensus.positions.len(), 8);

        // Check first position
        assert_eq!(consensus.positions[0].position, 1);
        assert_eq!(consensus.positions[0].allowed_aa, vec!['Q', 'E', 'D']);
        assert!(!consensus.positions[0].is_gap);

        // Check gap position
        let gap_pos = consensus.get_position(27).unwrap();
        assert!(gap_pos.is_gap);
        assert!(gap_pos.allowed_aa.is_empty());

        // Check conserved positions
        assert!(consensus.is_conserved(23));
        assert!(consensus.is_conserved(41));
        assert!(!consensus.is_conserved(1));

        // Check insertion points
        assert!(consensus.is_insertion_point(111));
        assert!(!consensus.is_insertion_point(23));
    }

    #[test]
    fn test_parse_without_metadata() {
        let text = r#"
1,Q,E,D
2,V,Q
23,C
//
"#;

        let consensus = Consensus::from_text(Chain::IGH, text).unwrap();
        assert_eq!(consensus.positions.len(), 3);
        assert!(consensus.metadata.conserved.is_empty());
        assert!(consensus.metadata.insertion_points.is_empty());
    }
}
