pub mod aho;
pub mod chothia;
pub mod imgt;
pub mod kabat;
pub mod martin;

use std::collections::HashMap;

use crate::aho::AHO_HEAVY_RULES;
use crate::alignment::AlignedPosition;
use crate::chothia::{CHOTHIA_HEAVY_RULES, CHOTHIA_LIGHT_RULES};
use crate::imgt::IMGT_RULES;
use crate::kabat::{KABAT_HEAVY_RULES, KABAT_LIGHT_RULES};
use crate::martin::{MARTIN_HEAVY_RULES, MARTIN_LIGHT_RULES};
use crate::types::{Chain, Insertion, NumberingRule, Position, Region, Scheme};

/// Get the region for a position number under the given scheme, or None if outside numbered range
pub fn region_for_position(pos: u8, scheme: Scheme) -> Option<Region> {
    match (scheme, pos) {
        (Scheme::IMGT, 1..=26) => Some(Region::FR1),
        (Scheme::IMGT, 27..=38) => Some(Region::CDR1),
        (Scheme::IMGT, 39..=55) => Some(Region::FR2),
        (Scheme::IMGT, 56..=65) => Some(Region::CDR2),
        (Scheme::IMGT, 66..=104) => Some(Region::FR3),
        (Scheme::IMGT, 105..=117) => Some(Region::CDR3),
        (Scheme::IMGT, 118..=128) => Some(Region::FR4),
        (Scheme::Kabat, 1..=25) => Some(Region::FR1),
        (Scheme::Kabat, 26..=35) => Some(Region::CDR1),
        (Scheme::Kabat, 36..=50) => Some(Region::FR2),
        (Scheme::Kabat, 51..=57) => Some(Region::CDR2),
        (Scheme::Kabat, 58..=92) => Some(Region::FR3),
        (Scheme::Kabat, 93..=100) => Some(Region::CDR3),
        (Scheme::Kabat, 101..=113) => Some(Region::FR4),
        // Chothia heavy/light region boundaries (Chothia CDR definitions)
        (Scheme::Chothia, 1..=25) => Some(Region::FR1),
        (Scheme::Chothia, 26..=32) => Some(Region::CDR1),
        (Scheme::Chothia, 33..=51) => Some(Region::FR2),
        (Scheme::Chothia, 52..=56) => Some(Region::CDR2),
        (Scheme::Chothia, 57..=95) => Some(Region::FR3),
        (Scheme::Chothia, 96..=101) => Some(Region::CDR3),
        (Scheme::Chothia, 102..=113) => Some(Region::FR4),
        // Martin uses the same region boundaries as Chothia
        (Scheme::Martin, 1..=25) => Some(Region::FR1),
        (Scheme::Martin, 26..=32) => Some(Region::CDR1),
        (Scheme::Martin, 33..=51) => Some(Region::FR2),
        (Scheme::Martin, 52..=56) => Some(Region::CDR2),
        (Scheme::Martin, 57..=95) => Some(Region::FR3),
        (Scheme::Martin, 96..=101) => Some(Region::CDR3),
        (Scheme::Martin, 102..=113) => Some(Region::FR4),
        // AHo region boundaries (1-149 space; AHo CDR definitions)
        (Scheme::Aho, 1..=24) => Some(Region::FR1),
        (Scheme::Aho, 25..=42) => Some(Region::CDR1),
        (Scheme::Aho, 43..=56) => Some(Region::FR2),
        (Scheme::Aho, 57..=77) => Some(Region::CDR2),
        (Scheme::Aho, 78..=108) => Some(Region::FR3),
        (Scheme::Aho, 109..=138) => Some(Region::CDR3),
        (Scheme::Aho, 139..=149) => Some(Region::FR4),
        _ => None,
    }
}

/// Segment a numbered sequence into its constituent regions
///
/// Returns a HashMap with keys for all Region variants plus "Prefix" and "Postfix".
/// Prefix collects residues before the numbered region, Postfix those after.
/// All keys are always present, with empty strings for absent regions.
pub fn segment(positions: &[Position], sequence: &str, scheme: Scheme) -> HashMap<String, String> {
    let mut segments: HashMap<String, String> = [
        "prefix", "fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4", "postfix",
    ]
    .iter()
    .map(|&s| (s.to_string(), String::new()))
    .collect();

    for (position, ch) in positions.iter().zip(sequence.chars()) {
        let key = match region_for_position(position.number, scheme) {
            Some(region) => region.to_string(),
            None if position.number == 0 => "Prefix".to_string(),
            None => "Postfix".to_string(),
        };
        segments.get_mut(&key.to_lowercase()).unwrap().push(ch);
    }

    segments
}

/// Apply a numbering scheme to aligned positions
pub fn apply_numbering(
    aligned_positions: &[AlignedPosition],
    scheme: Scheme,
    chain: Chain,
) -> Vec<Position> {
    let consensus_positions = extract_consensus_positions(aligned_positions);

    let rules = match (scheme, chain) {
        (Scheme::IMGT, _) => IMGT_RULES,
        (Scheme::Kabat, Chain::IGH) => KABAT_HEAVY_RULES,
        (Scheme::Kabat, Chain::IGK) | (Scheme::Kabat, Chain::IGL) => KABAT_LIGHT_RULES,
        (Scheme::Chothia, Chain::IGH) => CHOTHIA_HEAVY_RULES,
        (Scheme::Chothia, Chain::IGK) | (Scheme::Chothia, Chain::IGL) => CHOTHIA_LIGHT_RULES,
        (Scheme::Martin, Chain::IGH) => MARTIN_HEAVY_RULES,
        (Scheme::Martin, Chain::IGK) | (Scheme::Martin, Chain::IGL) => MARTIN_LIGHT_RULES,
        (Scheme::Aho, Chain::IGH) => AHO_HEAVY_RULES,
        _ => unreachable!("invalid scheme/chain combination should be prevented by Annotator"),
    };
    number_by_rules(&consensus_positions, rules)
}

/// Extract consensus position for each residue from alignment
///
/// Returns a position for each residue (non-gap) in the alignment.
/// Insertions inherit the position of the preceding aligned residue.
pub(crate) fn extract_consensus_positions(aligned: &[AlignedPosition]) -> Vec<u8> {
    let mut positions = Vec::with_capacity(aligned.len());
    let mut last_pos: u8 = 0;

    for ap in aligned {
        match ap {
            AlignedPosition::Aligned(pos) => {
                positions.push(*pos);
                last_pos = *pos;
            }
            AlignedPosition::Insertion() => {
                positions.push(last_pos); // Inherit from previous
            }
        }
    }

    positions
}

/// Apply rule-based numbering to consensus positions
fn number_by_rules(consensus_positions: &[u8], rules: &[NumberingRule]) -> Vec<Position> {
    let mut numbered_positions = Vec::with_capacity(consensus_positions.len());
    let mut idx = 0;

    for rule in rules {
        // Skip consensus positions that fall in an inter-rule gap (a scheme "hole",
        // e.g. IMGT position 10 or 73 which some schemes omit). These are treated as
        // insertions of the previously numbered position so that a residue present in
        // the alignment is never silently dropped, which would otherwise stall the
        // scanner and truncate the whole tail of the numbering.
        while idx < consensus_positions.len()
            && consensus_positions[idx] < rule.align_start
            && !rule.contains(consensus_positions[idx])
        {
            if let Some(last) = numbered_positions.last().copied() {
                let last: Position = last;
                let next_letter = match last.insertion {
                    Some(c) => (c as u8 + 1) as char,
                    None => 'A',
                };
                numbered_positions.push(Position::with_insertion(last.number, next_letter));
            }
            idx += 1;
        }

        let rule_start = idx;

        // Find all positions belonging to this rule's source range
        while idx < consensus_positions.len() && rule.contains(consensus_positions[idx]) {
            idx += 1;
        }

        let rule_len = idx - rule_start;
        if rule_len == 0 {
            continue;
        }

        match rule.insertion {
            Insertion::None => {
                numbered_positions.extend(number_with_offset(
                    &consensus_positions[rule_start..idx],
                    rule.align_start,
                    rule.num_start,
                ));
            }
            _ => {
                numbered_positions.extend(number_with_rules(rule_len, rule));
            }
        }
    }

    numbered_positions
}

/// Apply offset numbering to a region, handling insertions
fn number_with_offset(positions: &[u8], align_start: u8, num_start: u8) -> Vec<Position> {
    let mut result = Vec::with_capacity(positions.len());
    let mut last_align_pos: Option<u8> = None;
    let mut insertion_count = 0u8;

    for &pos in positions {
        if Some(pos) == last_align_pos {
            // Same position as previous = insertion
            insertion_count += 1;
            let dst_pos = pos.wrapping_sub(align_start).wrapping_add(num_start);
            result.push(Position::with_insertion(
                dst_pos,
                (b'A' + insertion_count - 1) as char,
            ));
        } else {
            // New position
            let dst_pos = pos.wrapping_sub(align_start).wrapping_add(num_start);
            result.push(Position::new(dst_pos));
            last_align_pos = Some(pos);
            insertion_count = 0;
        }
    }

    result
}

/// Generate positions for a CDR-like region based on its length and rules
///
/// Handles both deletions (len < base range) and insertions (len > base range).
pub fn number_with_rules(len: usize, rule: &NumberingRule) -> Vec<Position> {
    if len == 0 {
        return Vec::new();
    }

    let base_len = (rule.num_end - rule.num_start + 1) as usize;
    let mut result = Vec::with_capacity(len);

    if len <= base_len {
        // Deletions: select which base positions to keep
        let to_remove = base_len - len;

        // Build a mask of positions to skip
        let skip_set: &[u8] = &rule.deletion_order[..to_remove];
        for pos in rule.num_start..=rule.num_end {
            if !skip_set.contains(&pos) {
                result.push(Position::new(pos));
            }
        }
    } else {
        // Insertions: all base positions plus extra with letters
        let extra = len - base_len;

        match rule.insertion {
            Insertion::Sequential(insertion_pos) => {
                for pos in rule.num_start..=rule.num_end {
                    result.push(Position::new(pos));
                    if pos == insertion_pos {
                        for i in 0..extra {
                            result.push(Position::with_insertion(pos, (b'A' + i as u8) as char));
                        }
                    }
                }
            }
            Insertion::Symmetric { left, right } => {
                let insertions_left = extra / 2;
                let insertions_right = extra.div_ceil(2);

                for pos in rule.num_start..=rule.num_end {
                    if pos == left {
                        result.push(Position::new(pos));
                        // Left insertions: A, B, C...
                        for i in 0..insertions_left {
                            result.push(Position::with_insertion(left, (b'A' + i as u8) as char));
                        }
                        // Right insertions: ...B, A (reverse order)
                        for i in (0..insertions_right).rev() {
                            result.push(Position::with_insertion(right, (b'A' + i as u8) as char));
                        }
                    } else {
                        result.push(Position::new(pos));
                    }
                }
            }
            Insertion::None => unreachable!("offset rules handled separately"),
        }
    }

    result
}
