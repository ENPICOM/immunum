pub mod aho;
pub mod chothia;
pub mod imgt;
pub mod kabat;
pub mod martin;

use std::collections::HashMap;

use crate::aho::{AHO_HEAVY_RULES, AHO_KAPPA_RULES, AHO_LAMBDA_RULES, AHO_REGIONS};
use crate::alignment::AlignedPosition;
use crate::chothia::{
    CHOTHIA_HEAVY_REGIONS, CHOTHIA_HEAVY_RULES, CHOTHIA_LIGHT_REGIONS, CHOTHIA_LIGHT_RULES,
};
use crate::imgt::{IMGT_REGIONS, IMGT_RULES};
use crate::kabat::{
    KABAT_HEAVY_REGIONS, KABAT_HEAVY_RULES, KABAT_LIGHT_REGIONS, KABAT_LIGHT_RULES,
};
use crate::martin::{
    MARTIN_HEAVY_REGIONS, MARTIN_HEAVY_RULES, MARTIN_LIGHT_REGIONS, MARTIN_LIGHT_RULES,
};
use crate::types::{Chain, Insertion, NumberingRule, Position, Region, RegionDefinition, Scheme};

/// Region layout for a scheme and chain.
///
/// Kabat, Chothia and Martin define their CDRs differently for heavy and light chains -- CDR1 is
/// 31-35 on a heavy chain but 24-34 on a light one, and light numbering stops at 107 rather than
/// 113 -- so the chain is required, not incidental. IMGT and AHo are deliberately chain-agnostic:
/// both align positions so that one number means one structural position across chain types.
pub fn regions_for(scheme: Scheme, chain: Chain) -> RegionDefinition {
    let heavy = matches!(chain, Chain::IGH);
    match (scheme, heavy) {
        (Scheme::IMGT, _) => IMGT_REGIONS,
        (Scheme::Aho, _) => AHO_REGIONS,
        (Scheme::Kabat, true) => KABAT_HEAVY_REGIONS,
        (Scheme::Kabat, false) => KABAT_LIGHT_REGIONS,
        (Scheme::Chothia, true) => CHOTHIA_HEAVY_REGIONS,
        (Scheme::Chothia, false) => CHOTHIA_LIGHT_REGIONS,
        (Scheme::Martin, true) => MARTIN_HEAVY_REGIONS,
        (Scheme::Martin, false) => MARTIN_LIGHT_REGIONS,
    }
}

/// Get the region for a position number under the given scheme and chain, or None if outside the
/// numbered range.
pub fn region_for_position(pos: u8, scheme: Scheme, chain: Chain) -> Option<Region> {
    regions_for(scheme, chain).region(pos)
}

/// Segment a numbered sequence into its constituent regions
///
/// Returns a HashMap with keys for all Region variants plus "Prefix" and "Postfix".
/// Prefix collects residues before the numbered region, Postfix those after.
/// All keys are always present, with empty strings for absent regions.
pub fn segment(
    positions: &[Position],
    sequence: &str,
    scheme: Scheme,
    chain: Chain,
) -> HashMap<String, String> {
    let mut segments: HashMap<String, String> = [
        "prefix", "fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4", "postfix",
    ]
    .iter()
    .map(|&s| (s.to_string(), String::new()))
    .collect();

    for (position, ch) in positions.iter().zip(sequence.chars()) {
        let key = match region_for_position(position.number, scheme, chain) {
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
        (Scheme::Aho, Chain::IGK) => AHO_KAPPA_RULES,
        (Scheme::Aho, Chain::IGL) => AHO_LAMBDA_RULES,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbering::imgt::IMGT_RULES;
    use crate::numbering::kabat::{KABAT_HEAVY_RULES, KABAT_LIGHT_RULES};

    /// Published CDR boundaries, per scheme and chain, from the Martin group's table of CDR
    /// definitions <http://bioinf.org.uk/abs/info.html>: the Kabat columns, the post-June-2021
    /// Chothia consensus, and the AbM column for Martin. That table names no "Martin" definition --
    /// "Martin (Enhanced Chothia)" is a numbering scheme there -- so Martin numbering is paired with
    /// AbM, whose H1 row for that numbering convention reads H26-H35. The Kabat and Chothia rows
    /// were cross-checked against AbNumber's `SCHEME_BORDERS` (github.com/prihoda/AbNumber). These
    /// were chain-agnostic and wrong for Kabat in both chains before: CDR1 was 26-35 (the AbM
    /// definition, correct for Martin but not Kabat) and CDR-H2 was 51-57 against Kabat's 50-65,
    /// which put nine CDR-H2 residues in FR3.
    #[test]
    fn region_boundaries_match_published_definitions() {
        /// Scheme, chain, the three CDR spans as (first, last), and the last numbered position.
        type Case = (Scheme, Chain, [(u8, u8); 3], u8);

        let cases: &[Case] = &[
            (
                Scheme::IMGT,
                Chain::IGH,
                [(27, 38), (56, 65), (105, 117)],
                128,
            ),
            (
                Scheme::IMGT,
                Chain::IGK,
                [(27, 38), (56, 65), (105, 117)],
                128,
            ),
            (
                Scheme::Kabat,
                Chain::IGH,
                [(31, 35), (50, 65), (95, 102)],
                113,
            ),
            (
                Scheme::Kabat,
                Chain::IGK,
                [(24, 34), (50, 56), (89, 97)],
                107,
            ),
            (
                Scheme::Kabat,
                Chain::IGL,
                [(24, 34), (50, 56), (89, 97)],
                107,
            ),
            (
                Scheme::Chothia,
                Chain::IGH,
                [(26, 32), (52, 56), (96, 101)],
                113,
            ),
            (
                Scheme::Chothia,
                Chain::IGK,
                [(26, 32), (50, 52), (91, 96)],
                107,
            ),
            (
                Scheme::Martin,
                Chain::IGH,
                [(26, 35), (50, 58), (95, 102)],
                113,
            ),
            (
                Scheme::Martin,
                Chain::IGL,
                [(24, 34), (50, 56), (89, 97)],
                107,
            ),
        ];

        for &(scheme, chain, cdrs, fr4_end) in cases {
            let regions = regions_for(scheme, chain);
            let expected = [Region::CDR1, Region::CDR2, Region::CDR3];
            for (&(first, last), region) in cdrs.iter().zip(expected) {
                assert_eq!(
                    region_for_position(first, scheme, chain),
                    Some(region),
                    "{scheme:?} {chain} position {first} should open {region:?}"
                );
                assert_eq!(
                    region_for_position(last, scheme, chain),
                    Some(region),
                    "{scheme:?} {chain} position {last} should close {region:?}"
                );
                // The residue either side must fall in the flanking framework, which is what
                // catches an off-by-one at a boundary.
                assert_ne!(
                    region_for_position(first - 1, scheme, chain),
                    Some(region),
                    "{scheme:?} {chain} {region:?} starts too early: {} is inside it",
                    first - 1
                );
                assert_ne!(
                    region_for_position(last + 1, scheme, chain),
                    Some(region),
                    "{scheme:?} {chain} {region:?} ends too late: {} is inside it",
                    last + 1
                );
            }
            assert_eq!(regions.fr4_end, fr4_end, "{scheme:?} {chain} FR4 end");
            assert_eq!(region_for_position(fr4_end + 1, scheme, chain), None);
            assert_eq!(region_for_position(0, scheme, chain), None);
        }
    }

    /// The regions tile `1..=fr4_end` with no gap and no overlap, for every scheme and chain, and
    /// the tiles sum to the scheme's canonical domain length.
    ///
    /// Numbering is 1-based: position 1 opens FR1 and 0 is not a numbered position at all. Zero is
    /// the sentinel for a residue ahead of the domain -- `extract_consensus_positions` seeds its
    /// carry at 0, and `segment` routes 0 to the prefix bucket -- so it must keep mapping to `None`.
    #[test]
    fn regions_tile_the_numbered_range() {
        // Canonical domain length per scheme and chain class: IMGT 128, AHo 149, and 113/107 for
        // the Kabat-derived schemes on heavy/light.
        let expected_len = |scheme: Scheme, heavy: bool| match (scheme, heavy) {
            (Scheme::IMGT, _) => 128,
            (Scheme::Aho, _) => 149,
            (_, true) => 113,
            (_, false) => 107,
        };

        for scheme in [
            Scheme::IMGT,
            Scheme::Kabat,
            Scheme::Chothia,
            Scheme::Martin,
            Scheme::Aho,
        ] {
            for chain in [Chain::IGH, Chain::IGK, Chain::IGL] {
                let regions = regions_for(scheme, chain);

                assert_eq!(
                    region_for_position(1, scheme, chain),
                    Some(Region::FR1),
                    "{scheme:?} {chain} numbering must be 1-based: position 1 opens FR1"
                );
                assert_eq!(
                    region_for_position(0, scheme, chain),
                    None,
                    "{scheme:?} {chain} position 0 is the prefix sentinel, not a numbered position"
                );

                for pos in 1..=regions.fr4_end {
                    assert!(
                        region_for_position(pos, scheme, chain).is_some(),
                        "{scheme:?} {chain} position {pos} falls in no region"
                    );
                }
                assert_eq!(
                    region_for_position(regions.fr4_end + 1, scheme, chain),
                    None,
                    "{scheme:?} {chain} position past FR4 must be postfix"
                );
                assert_eq!(
                    regions.fr4_end,
                    expected_len(scheme, matches!(chain, Chain::IGH)),
                    "{scheme:?} {chain} regions should tile the scheme's full domain length"
                );
            }
        }
    }

    /// Heavy and light must not share a table where the scheme distinguishes them, which is the
    /// bug this replaced: one table served both.
    ///
    /// Compared on CDR2 and the last numbered position rather than CDR1, because Chothia places
    /// CDR1 at 26-32 in *both* chains -- for it only CDR2, CDR3 and the FR4 end diverge. Kabat's
    /// CDR1 does differ (31-35 heavy, 24-34 light), as does Martin's AbM-derived one (26-35 heavy,
    /// 24-34 light), and the boundary test above covers both.
    #[test]
    fn chain_specific_schemes_differ_between_heavy_and_light() {
        for scheme in [Scheme::Kabat, Scheme::Chothia, Scheme::Martin] {
            let heavy = regions_for(scheme, Chain::IGH);
            let light = regions_for(scheme, Chain::IGK);
            assert_ne!(
                (heavy.cdr2_end, heavy.fr4_end),
                (light.cdr2_end, light.fr4_end),
                "{scheme:?} CDR2/FR4 should differ between heavy and light"
            );
        }
        // IMGT and AHo align positions structurally across chain types, so sharing is correct.
        for scheme in [Scheme::IMGT, Scheme::Aho] {
            let heavy = regions_for(scheme, Chain::IGH);
            let light = regions_for(scheme, Chain::IGK);
            assert_eq!(
                heavy.cdr1_end, light.cdr1_end,
                "{scheme:?} should be chain-agnostic"
            );
            assert_eq!(
                heavy.fr4_end, light.fr4_end,
                "{scheme:?} should be chain-agnostic"
            );
        }
    }

    /// Only rules with an insertion strategy reach `number_with_rules`; `number_by_rules` sends the
    /// framework rules down the offset path instead, where `deletion_order` is never consulted.
    fn variable_rules() -> Vec<(&'static str, &'static NumberingRule)> {
        [
            ("KABAT_HEAVY_RULES", KABAT_HEAVY_RULES),
            ("KABAT_LIGHT_RULES", KABAT_LIGHT_RULES),
            ("IMGT_RULES", IMGT_RULES),
        ]
        .into_iter()
        .flat_map(|(name, rules)| {
            rules
                .iter()
                .filter(|rule| !matches!(rule.insertion, Insertion::None))
                .map(move |rule| (name, rule))
        })
        .collect()
    }

    /// A region one residue long asks to delete every base position but one, so `deletion_order`
    /// needs `base_len - 1` entries. Three Kabat rules were short of that and sliced out of bounds
    /// (`&rule.deletion_order[..to_remove]`), which surfaced as a panic on truncated reads.
    #[test]
    fn deletion_order_covers_the_shortest_possible_region() {
        for (table, rule) in variable_rules() {
            let base_len = (rule.num_end - rule.num_start + 1) as usize;
            assert!(
                rule.deletion_order.len() >= base_len - 1,
                "{table} rule {}-{} has {} deletion positions but needs {} to represent a \
                 single-residue region",
                rule.num_start,
                rule.num_end,
                rule.deletion_order.len(),
                base_len - 1,
            );
        }
    }

    /// The contract of `number_with_rules`: one position per residue, for every length a region can
    /// take. Exercises the deletion path (below base length) and the insertion path (above it).
    #[test]
    fn every_region_length_yields_one_position_per_residue() {
        for (table, rule) in variable_rules() {
            let base_len = (rule.num_end - rule.num_start + 1) as usize;
            for len in 1..=base_len + 4 {
                let positions = number_with_rules(len, rule);
                assert_eq!(
                    positions.len(),
                    len,
                    "{table} rule {}-{} returned {} positions for a region of {len} residues",
                    rule.num_start,
                    rule.num_end,
                    positions.len(),
                );
            }
        }
    }

    #[test]
    fn empty_region_yields_no_positions() {
        for (_, rule) in variable_rules() {
            assert!(number_with_rules(0, rule).is_empty());
        }
    }
}
