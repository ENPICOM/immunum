use crate::consensus_scoring::{read_consensus_file, read_scoring_matrix};
use crate::numbering_scheme_type::NumberingScheme;
use crate::types::{Chain, RegionRange, Scheme};
use std::path::PathBuf;

pub fn get_imgt_heavy_scheme() -> NumberingScheme {
    NumberingScheme {
        scheme_type: Scheme::IMGT,
        chain_type: Chain::IGH,
        conserved_positions: vec![23, 41, 104, 118, 119, 121],
        insertion_positions: vec![32, 60, 111],
        gap_positions: vec![10, 73],
        consensus_amino_acids: read_consensus_file(
            PathBuf::from("resources")
                .join("consensus")
                .join("IMGT_CONSENSUS_H.txt"),
        ),
        scoring_matrix: read_scoring_matrix(
            PathBuf::from("resources")
                .join("consensus")
                .join("IMGT_CONSENSUS_H.npy"),
        ),
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
    }
}

pub fn get_imgt_kappa_scheme() -> NumberingScheme {
    NumberingScheme {
        scheme_type: Scheme::IMGT,
        chain_type: Chain::IGK,
        conserved_positions: vec![23, 41, 104, 118, 119, 121],
        insertion_positions: vec![32, 60, 111],
        gap_positions: vec![10, 73, 81, 82],
        consensus_amino_acids: read_consensus_file(
            PathBuf::from("resources")
                .join("consensus")
                .join("IMGT_CONSENSUS_K.txt"),
        ),
        scoring_matrix: read_scoring_matrix(
            PathBuf::from("resources")
                .join("consensus")
                .join("IMGT_CONSENSUS_K.npy"),
        ),
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
    }
}

pub fn get_imgt_lambda_scheme() -> NumberingScheme {
    NumberingScheme {
        scheme_type: Scheme::IMGT,
        chain_type: Chain::IGL,
        conserved_positions: vec![23, 41, 104, 118, 119, 121],
        insertion_positions: vec![32, 60, 111],
        gap_positions: vec![10, 73, 81, 82],
        consensus_amino_acids: read_consensus_file(
            PathBuf::from("resources")
                .join("consensus")
                .join("IMGT_CONSENSUS_L.txt"),
        ),
        scoring_matrix: read_scoring_matrix(
            PathBuf::from("resources")
                .join("consensus")
                .join("IMGT_CONSENSUS_L.npy"),
        ),
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
    }
}

pub fn get_kabat_heavy_scheme() -> NumberingScheme {
    NumberingScheme {
        // name: "KABAT Heavy".to_string(),
        // description: "KABAT numbering scheme for heavy chains".to_string(),
        scheme_type: Scheme::KABAT,
        chain_type: Chain::IGH,
        conserved_positions: vec![22, 36, 92, 103, 104, 106],
        insertion_positions: vec![6, 82],
        gap_positions: vec![40, 41, 42, 43, 44, 72, 73, 74],
        consensus_amino_acids: read_consensus_file(
            PathBuf::from("resources")
                .join("consensus")
                .join("KABAT_CONSENSUS_H.txt"),
        ),
        scoring_matrix: read_scoring_matrix(
            PathBuf::from("resources")
                .join("consensus")
                .join("KABAT_CONSENSUS_H.npy"),
        ),
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
    }
}

pub fn get_kabat_kappa_scheme() -> NumberingScheme {
    NumberingScheme {
        scheme_type: Scheme::KABAT,
        chain_type: Chain::IGK,
        conserved_positions: vec![23, 35, 88, 98, 99, 101],
        insertion_positions: vec![27],
        gap_positions: vec![10],
        consensus_amino_acids: read_consensus_file(
            PathBuf::from("resources")
                .join("consensus")
                .join("KABAT_CONSENSUS_K.txt"),
        ),
        scoring_matrix: read_scoring_matrix(
            PathBuf::from("resources")
                .join("consensus")
                .join("KABAT_CONSENSUS_K.npy"),
        ),
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
    }
}

pub fn get_kabat_lambda_scheme() -> NumberingScheme {
    NumberingScheme {
        scheme_type: Scheme::KABAT,
        chain_type: Chain::IGL,
        conserved_positions: vec![23, 35, 88, 98, 99, 101],
        insertion_positions: vec![27],
        gap_positions: vec![10],
        consensus_amino_acids: read_consensus_file(
            PathBuf::from("resources")
                .join("consensus")
                .join("KABAT_CONSENSUS_L.txt"),
        ),
        scoring_matrix: read_scoring_matrix(
            PathBuf::from("resources")
                .join("consensus")
                .join("KABAT_CONSENSUS_L.npy"),
        ),
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
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{get_scoring_params};

    #[test]
    fn scheme_creation() {
        let scheme = get_imgt_heavy_scheme();
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
        //TODO add more tests
    }
}
