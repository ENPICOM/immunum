use crate::consensus_scoring::{calculate_scoring_matrix, GapPenaltyConfig};
use crate::constants::{get_consensus, get_region_ranges, get_scoring_params};
use crate::kmer_prefiltering::generate_scheme_kmers;
use crate::numbering_scheme_type::NumberingScheme;
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

    // Create cached HashSets for efficient lookups
    let conserved_positions_set: std::collections::HashSet<u32> =
        config.conserved_positions.iter().cloned().collect();
    let restricted_sites_set: std::collections::HashSet<u32> = consensus_amino_acids
        .iter()
        .filter_map(|(&key, value)| {
            if !value.contains(&b'-') {
                Some(key)
            } else {
                None
            }
        })
        .collect();

    // Create gap penalty configuration
    let gap_penalty_config = GapPenaltyConfig {
        conserved_positions_set: conserved_positions_set.clone(),
        scheme_type: scheme,
        chain_type: chain,
        insertion_positions: config.insertion_positions,
        gap_positions: config.gap_positions,
        region_ranges,
    };

    // Use the gap penalty config directly
    let scoring_matrix =
        calculate_scoring_matrix(&consensus_amino_acids, &scoring_params, &gap_penalty_config);

    NumberingScheme {
        scheme_type: scheme,
        chain_type: chain,
        cdr_definition: cdr_definitions,
        consensus_amino_acids: consensus_amino_acids.clone(),
        scoring_matrix,
        conserved_positions_set,
        restricted_sites_set,
        kmer_set: generate_scheme_kmers(consensus_amino_acids.clone()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use crate::constants::get_scoring_params;
    use crate::types::CdrDefinitions;

    #[test]
    fn scheme_creation() {
        let scheme = get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGH, CdrDefinitions::IMGT);
        assert_eq!(scheme.consensus_amino_acids.len(), 128);
        assert_eq!(scheme.consensus_amino_acids[&1], vec![b'Q', b'E', b'D']);
    }
}
