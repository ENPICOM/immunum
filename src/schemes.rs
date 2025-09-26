use crate::consensus_scoring::calculate_scoring_matrix;
use crate::constants::{get_consensus, get_region_ranges, get_scoring_params};
use crate::numbering_scheme_type::NumberingScheme;
use crate::scoring_matrix::ScoringMatrix;
use crate::types::{CdrDefinitions, Chain, Scheme};

struct SchemeConfig {
    conserved_positions: Vec<u32>,
    insertion_positions: Vec<u32>,
    gap_positions: Vec<u32>,
}

fn get_scheme_config(scheme: &Scheme, chain: &Chain) -> SchemeConfig {
    match (scheme, chain) {
        (Scheme::IMGT, Chain::IGH) => SchemeConfig {
            conserved_positions: vec![23, 41, 104, 118, 119, 121],
            insertion_positions: vec![32, 60, 111],
            gap_positions: vec![10, 73],
        },
        (Scheme::IMGT, Chain::IGK) => SchemeConfig {
            conserved_positions: vec![23, 41, 104, 118, 119, 121],
            insertion_positions: vec![32, 60, 111],
            gap_positions: vec![10, 73, 81, 82],
        },
        (Scheme::IMGT, Chain::IGL) => SchemeConfig {
            conserved_positions: vec![23, 41, 104, 118, 119, 121],
            insertion_positions: vec![32, 60, 111],
            gap_positions: vec![10, 73, 81, 82],
        },
        (Scheme::KABAT, Chain::IGH) => SchemeConfig {
            conserved_positions: vec![22, 36, 92, 103, 104, 106],
            insertion_positions: vec![6, 35, 52, 82, 100],
            gap_positions: vec![40, 41, 42, 43, 44, 72, 73, 74],
        },
        (Scheme::KABAT, Chain::IGK) => SchemeConfig {
            conserved_positions: vec![23, 35, 88, 98, 99, 101],
            insertion_positions: vec![27, 52, 95],
            gap_positions: vec![10],
        },
        (Scheme::KABAT, Chain::IGL) => SchemeConfig {
            conserved_positions: vec![23, 35, 88, 98, 99, 101],
            insertion_positions: vec![27, 52, 95],
            gap_positions: vec![10],
        },
        // For unsupported T-cell receptor chains, default to heavy chain patterns
        (Scheme::IMGT, _) => get_scheme_config(&Scheme::IMGT, &Chain::IGH),
        (Scheme::KABAT, _) => get_scheme_config(&Scheme::KABAT, &Chain::IGH),
    }
}

pub fn get_scheme_with_cdr_definition(
    scheme: Scheme,
    chain: Chain,
    cdr_definitions: CdrDefinitions,
) -> NumberingScheme {
    let scoring_params = get_scoring_params();
    let consensus_amino_acids = get_consensus(scheme, chain);
    let config = get_scheme_config(&scheme, &chain);
    let region_ranges = get_region_ranges(cdr_definitions, scheme, chain);

    // Create a temporary scheme to access the gap_penalty method
    let temp_conserved_set: std::collections::HashSet<u32> =
        config.conserved_positions.iter().cloned().collect();
    let temp_restricted_set: std::collections::HashSet<u32> = consensus_amino_acids
        .iter()
        .filter_map(|(&key, value)| {
            if !value.contains(&b'-') {
                Some(key)
            } else {
                None
            }
        })
        .collect();

    let temp_scheme = NumberingScheme {
        scheme_type: scheme,
        chain_type: chain,
        cdr_definition: cdr_definitions,
        insertion_positions: config.insertion_positions.clone(),
        gap_positions: config.gap_positions.clone(),
        consensus_amino_acids: consensus_amino_acids.clone(),
        scoring_matrix: ScoringMatrix::zeros(1, 1), // Temporary
        fr1: region_ranges.fr1.clone(),
        cdr1: region_ranges.cdr1.clone(),
        fr2: region_ranges.fr2.clone(),
        cdr2: region_ranges.cdr2.clone(),
        fr3: region_ranges.fr3.clone(),
        cdr3: region_ranges.cdr3.clone(),
        fr4: region_ranges.fr4.clone(),
        conserved_positions_set: temp_conserved_set.clone(),
        restricted_sites_set: temp_restricted_set.clone(),
    };

    let scoring_matrix =
        calculate_scoring_matrix(&consensus_amino_acids, &scoring_params, |pos, params| {
            temp_scheme.gap_penalty(pos, params)
        });

    // Re-use the cached HashSets from temp_scheme
    let conserved_positions_set = temp_conserved_set;
    let restricted_sites_set = temp_restricted_set;

    NumberingScheme {
        scheme_type: scheme,
        chain_type: chain,
        cdr_definition: cdr_definitions,
        insertion_positions: config.insertion_positions,
        gap_positions: config.gap_positions,
        consensus_amino_acids,
        scoring_matrix,
        fr1: region_ranges.fr1,
        cdr1: region_ranges.cdr1,
        fr2: region_ranges.fr2,
        cdr2: region_ranges.cdr2,
        fr3: region_ranges.fr3,
        cdr3: region_ranges.cdr3,
        fr4: region_ranges.fr4,
        conserved_positions_set,
        restricted_sites_set,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::get_scoring_params;
    use crate::types::CdrDefinitions;

    #[test]
    fn scheme_creation() {
        let scheme = get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGH, CdrDefinitions::IMGT);
        let scoring_params = get_scoring_params();
        assert_eq!(scheme.gap_positions, vec![10, 73]);
        assert_eq!(scheme.consensus_amino_acids.len(), 128);
        assert_eq!(scheme.consensus_amino_acids[&1], vec![b'Q', b'E', b'D']);
        assert_eq!(
            scheme.gap_penalty(25, &scoring_params),
            (scoring_params.gap_pen_fr, scoring_params.gap_pen_fr)
        );
        assert_eq!(
            scheme.gap_penalty(200, &scoring_params),
            (scoring_params.gap_pen_other, scoring_params.gap_pen_other)
        );
        assert_eq!(
            scheme.gap_penalty(10, &scoring_params),
            (scoring_params.gap_pen_fr, scoring_params.gap_pen_op)
        );

        for i in 1..129 {
            println!("{:?},", scheme.gap_penalty(i as u32, &scoring_params))
        }
    }
}
