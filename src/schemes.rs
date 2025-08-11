use crate::consensus_scoring::calculate_scoring_matrix;
use crate::constants::{get_consensus, get_scoring_params, ScoringParams};
use crate::gap_penalty::calculate_gap_penalty;
use crate::numbering_scheme::NumberingScheme;
use crate::types::{Chain, RegionRange, Scheme};

pub struct SchemeConfig {
    pub conserved_positions: Vec<u32>,
    pub insertion_positions: Vec<u32>,
    pub gap_positions: Vec<u32>,
    pub fr1: RegionRange,
    pub cdr1: RegionRange,
    pub fr2: RegionRange,
    pub cdr2: RegionRange,
    pub fr3: RegionRange,
    pub cdr3: RegionRange,
    pub fr4: RegionRange,
}

fn get_scheme_config(scheme: &Scheme, chain: &Chain) -> SchemeConfig {
    match (scheme, chain) {
        (Scheme::IMGT, Chain::IGH) => SchemeConfig {
            conserved_positions: vec![23, 41, 104, 118, 119, 121],
            insertion_positions: vec![32, 60, 111],
            gap_positions: vec![10, 73],
            fr1: RegionRange { start: 1, end: 27 },
            cdr1: RegionRange { start: 27, end: 39 },
            fr2: RegionRange { start: 39, end: 56 },
            cdr2: RegionRange { start: 56, end: 66 },
            fr3: RegionRange {
                start: 66,
                end: 105,
            },
            cdr3: RegionRange {
                start: 105,
                end: 118,
            },
            fr4: RegionRange {
                start: 118,
                end: 129,
            },
        },
        (Scheme::IMGT, Chain::IGK) => SchemeConfig {
            conserved_positions: vec![23, 41, 104, 118, 119, 121],
            insertion_positions: vec![32, 60, 111],
            gap_positions: vec![10, 73, 81, 82],
            fr1: RegionRange { start: 1, end: 27 },
            cdr1: RegionRange { start: 27, end: 39 },
            fr2: RegionRange { start: 39, end: 56 },
            cdr2: RegionRange { start: 56, end: 66 },
            fr3: RegionRange {
                start: 66,
                end: 105,
            },
            cdr3: RegionRange {
                start: 105,
                end: 118,
            },
            fr4: RegionRange {
                start: 118,
                end: 128,
            },
        },
        (Scheme::IMGT, Chain::IGL) => SchemeConfig {
            conserved_positions: vec![23, 41, 104, 118, 119, 121],
            insertion_positions: vec![32, 60, 111],
            gap_positions: vec![10, 73, 81, 82],
            fr1: RegionRange { start: 1, end: 27 },
            cdr1: RegionRange { start: 27, end: 39 },
            fr2: RegionRange { start: 39, end: 56 },
            cdr2: RegionRange { start: 56, end: 66 },
            fr3: RegionRange {
                start: 66,
                end: 105,
            },
            cdr3: RegionRange {
                start: 105,
                end: 118,
            },
            fr4: RegionRange {
                start: 118,
                end: 129,
            },
        },
        (Scheme::KABAT, Chain::IGH) => SchemeConfig {
            conserved_positions: vec![22, 36, 92, 103, 104, 106],
            insertion_positions: vec![6, 82],
            gap_positions: vec![40, 41, 42, 43, 44, 72, 73, 74],
            fr1: RegionRange { start: 1, end: 31 },
            cdr1: RegionRange { start: 31, end: 36 },
            fr2: RegionRange { start: 36, end: 50 },
            cdr2: RegionRange { start: 50, end: 66 },
            fr3: RegionRange { start: 66, end: 95 },
            cdr3: RegionRange {
                start: 95,
                end: 103,
            },
            fr4: RegionRange {
                start: 103,
                end: 114,
            },
        },
        (Scheme::KABAT, Chain::IGK) => SchemeConfig {
            conserved_positions: vec![23, 35, 88, 98, 99, 101],
            insertion_positions: vec![27],
            gap_positions: vec![10],
            fr1: RegionRange { start: 1, end: 24 },
            cdr1: RegionRange { start: 24, end: 35 },
            fr2: RegionRange { start: 35, end: 50 },
            cdr2: RegionRange { start: 50, end: 57 },
            fr3: RegionRange { start: 57, end: 89 },
            cdr3: RegionRange { start: 89, end: 98 },
            fr4: RegionRange {
                start: 98,
                end: 108,
            },
        },
        (Scheme::KABAT, Chain::IGL) => SchemeConfig {
            conserved_positions: vec![23, 35, 88, 98, 99, 101],
            insertion_positions: vec![27],
            gap_positions: vec![10],
            fr1: RegionRange { start: 1, end: 24 },
            cdr1: RegionRange { start: 24, end: 35 },
            fr2: RegionRange { start: 35, end: 50 },
            cdr2: RegionRange { start: 50, end: 57 },
            fr3: RegionRange { start: 57, end: 89 },
            cdr3: RegionRange { start: 89, end: 98 },
            fr4: RegionRange {
                start: 98,
                end: 108,
            },
        },
        // For unsupported T-cell receptor chains, default to heavy chain patterns
        (Scheme::IMGT, _) => get_scheme_config(&Scheme::IMGT, &Chain::IGH),
        (Scheme::KABAT, _) => get_scheme_config(&Scheme::KABAT, &Chain::IGH),
    }
}

pub fn get_scheme(scheme: Scheme, chain: Chain, params: Option<ScoringParams>) -> NumberingScheme {
    let scoring_params = params.unwrap_or_else(get_scoring_params);
    let consensus_amino_acids_map = get_consensus(scheme, chain);
    let config = get_scheme_config(&scheme, &chain);

    // Convert HashMap to Vec for indexed access
    let max_position = consensus_amino_acids_map.keys().max().copied().unwrap_or(0) as usize;
    let mut consensus_amino_acids = vec![Vec::new(); max_position + 1];

    // Populate the indexed vector with consensus amino acids
    for (&pos, amino_acids) in &consensus_amino_acids_map {
        let index = pos as usize;
        if index < consensus_amino_acids.len() {
            consensus_amino_acids[index] = amino_acids.clone();
        }
    }

    // Pre-compute restricted sites (non '-' positions)
    let mut restricted_sites = std::collections::HashSet::new();
    for (&pos, amino_acids) in &consensus_amino_acids_map {
        if !amino_acids.contains(&b'-') {
            restricted_sites.insert(pos);
        }
    }

    // Convert conserved_positions to HashSet for O(1) lookup
    let conserved_positions: std::collections::HashSet<u32> =
        config.conserved_positions.iter().copied().collect();

    let scoring_matrix = calculate_scoring_matrix(
        &consensus_amino_acids_map,
        &scoring_params,
        |pos, params| calculate_gap_penalty(pos, scheme, chain, &config, params),
    );

    NumberingScheme {
        scheme_type: scheme,
        chain_type: chain,
        conserved_positions,
        consensus_amino_acids,
        restricted_sites,
        scoring_matrix,
        fr1: config.fr1,
        cdr1: config.cdr1,
        fr2: config.fr2,
        cdr2: config.cdr2,
        fr3: config.fr3,
        cdr3: config.cdr3,
        fr4: config.fr4,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::get_scoring_params;

    #[test]
    fn scheme_creation() {
        let scheme = get_scheme(Scheme::IMGT, Chain::IGH, None);
        let config = get_scheme_config(&Scheme::IMGT, &Chain::IGH);
        let scoring_params = get_scoring_params();
        assert_eq!(scheme.consensus_amino_acids.len(), 129); // Vec is 0-indexed plus max position
        assert_eq!(scheme.consensus_amino_acids[1], vec![b'Q', b'E', b'D']);
        assert_eq!(
            calculate_gap_penalty(
                25,
                scheme.scheme_type,
                scheme.chain_type,
                &config,
                &scoring_params,
            ),
            (scoring_params.gap_pen_fr, scoring_params.gap_pen_fr)
        );
        assert_eq!(
            calculate_gap_penalty(
                200,
                scheme.scheme_type,
                scheme.chain_type,
                &config,
                &scoring_params,
            ),
            (scoring_params.gap_pen_other, scoring_params.gap_pen_other)
        );
        assert_eq!(
            calculate_gap_penalty(
                10,
                scheme.scheme_type,
                scheme.chain_type,
                &config,
                &scoring_params,
            ),
            (scoring_params.gap_pen_fr, scoring_params.gap_pen_op)
        );

        for i in 1..129 {
            println!(
                "{:?},",
                calculate_gap_penalty(
                    i as u32,
                    scheme.scheme_type,
                    scheme.chain_type,
                    &config,
                    &scoring_params,
                )
            )
        }
    }
}
