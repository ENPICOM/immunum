use crate::consensus_scoring::calculate_scoring_matrix;
use crate::constants::{get_consensus, get_scoring_params};
use crate::numbering_scheme_type::NumberingScheme;
use crate::scoring_matrix::ScoringMatrix;
use crate::types::{Chain, RegionRange, Scheme};

struct SchemeConfig {
    conserved_positions: Vec<u32>,
    insertion_positions: Vec<u32>,
    gap_positions: Vec<u32>,
    fr1: RegionRange,
    cdr1: RegionRange,
    fr2: RegionRange,
    cdr2: RegionRange,
    fr3: RegionRange,
    cdr3: RegionRange,
    fr4: RegionRange,
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
            insertion_positions: vec![6, 35, 52, 82, 100],
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
            insertion_positions: vec![27, 52, 95],
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
            insertion_positions: vec![27, 52, 95],
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

pub fn get_scheme(scheme: Scheme, chain: Chain) -> NumberingScheme {
    let scoring_params = get_scoring_params();
    let consensus_amino_acids = get_consensus(scheme, chain);
    let config = get_scheme_config(&scheme, &chain);

    // Create a temporary scheme to access the gap_penalty method
    let temp_scheme = NumberingScheme {
        scheme_type: scheme,
        chain_type: chain,
        conserved_positions: config.conserved_positions.clone(),
        insertion_positions: config.insertion_positions.clone(),
        gap_positions: config.gap_positions.clone(),
        consensus_amino_acids: consensus_amino_acids.clone(),
        scoring_matrix: ScoringMatrix::zeros(1, 1), // Temporary
        fr1: config.fr1.clone(),
        cdr1: config.cdr1.clone(),
        fr2: config.fr2.clone(),
        cdr2: config.cdr2.clone(),
        fr3: config.fr3.clone(),
        cdr3: config.cdr3.clone(),
        fr4: config.fr4.clone(),
    };

    let scoring_matrix =
        calculate_scoring_matrix(&consensus_amino_acids, &scoring_params, |pos, params| {
            temp_scheme.gap_penalty(pos, params)
        });

    NumberingScheme {
        scheme_type: scheme,
        chain_type: chain,
        conserved_positions: config.conserved_positions,
        insertion_positions: config.insertion_positions,
        gap_positions: config.gap_positions,
        consensus_amino_acids,
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
        let scheme = get_scheme(Scheme::IMGT, Chain::IGH);
        let scoring_params = get_scoring_params();
        assert_eq!(scheme.gap_positions, vec![10, 73]);
        assert_eq!(scheme.restricted_sites().len(), 88);
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
