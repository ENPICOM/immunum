use crate::consensus_scoring::{read_partial_consensus_file, read_partial_scoring_matrix};
use crate::numbering_scheme_type::NumberingScheme;
use crate::schemes::{get_imgt_heavy_scheme, get_imgt_kappa_scheme, get_imgt_lambda_scheme};
use crate::types::Chain;
use std::path::PathBuf;

/// Create c-terminal NumberingScheme from a full length scheme (based on IMGT consensus)
pub fn get_terminal_schemes(
    chain: Chain,
    terminal_length: u8,
) -> (NumberingScheme, NumberingScheme) {
    let original_scheme = match chain {
        Chain::IGH => get_imgt_heavy_scheme(),
        Chain::IGK => get_imgt_kappa_scheme(),
        Chain::IGL => get_imgt_lambda_scheme(),
        other => !panic!("No scheme implemented for chain {:?} yet", other),
    };
    let file_name = original_scheme.file_name.clone();

    // Create N-terminal scheme
    let n_terminal_scheme = NumberingScheme {
        conserved_positions: vec![], // No conserved positions at N-terminal
        insertion_positions: vec![],
        gap_positions: vec![10], // IMGT Position 10
        consensus_amino_acids: read_partial_consensus_file(
            PathBuf::from("resources")
                .join("consensus")
                .join(file_name.clone() + ".txt"),
            1,
            (terminal_length) as usize,
        ),
        scoring_matrix: read_partial_scoring_matrix(
            PathBuf::from("resources")
                .join("consensus")
                .join(file_name.clone() + ".npy"),
            1,
            terminal_length as usize,
        ),
        ..original_scheme.clone() // Use the remaining fields from the original
    };

    // Calculate C-terminal positions
    let c_term_start = 118; // skip first line, start at position 118
    let c_term_end = c_term_start + terminal_length as usize - 1;

    // Create C-terminal scheme
    let c_terminal_scheme = NumberingScheme {
        conserved_positions: vec![1, 2, 4], // IMGT position 118, 119 and 121
        insertion_positions: vec![],
        gap_positions: vec![], // No gap positions at c-terminal
        consensus_amino_acids: read_partial_consensus_file(
            PathBuf::from("resources")
                .join("consensus")
                .join(file_name.clone() + ".txt"),
            c_term_start,
            c_term_end,
        ),
        scoring_matrix: read_partial_scoring_matrix(
            PathBuf::from("resources")
                .join("consensus")
                .join(file_name.clone() + ".npy"),
            c_term_start,
            c_term_end,
        ),
        ..original_scheme // Use the remaining fields from the original
    };

    // Return the two new schemes
    (n_terminal_scheme, c_terminal_scheme)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::schemes::{get_imgt_heavy_scheme, get_imgt_kappa_scheme};
    use ndarray::Ix;

    #[test]
    fn terminal_schemes_heavy() {
        let terminal_length = 10;
        let original_scheme = get_imgt_heavy_scheme();
        let (n_term_scheme, c_term_scheme) = get_terminal_schemes(Chain::IGH, terminal_length);
        for i in 0..terminal_length {
            assert_eq!(
                n_term_scheme.scoring_matrix.row(i as Ix),
                original_scheme.scoring_matrix.row(i as Ix)
            );

            assert_eq!(
                c_term_scheme.scoring_matrix.row(i as Ix),
                original_scheme.scoring_matrix.row(117 + i as Ix)
            );
        }

        for i in 1..terminal_length {
            assert_eq!(
                n_term_scheme.consensus_amino_acids[&(i as u32)],
                original_scheme.consensus_amino_acids[&(i as u32)]
            );
            assert_eq!(
                c_term_scheme.consensus_amino_acids[&(i as u32)],
                original_scheme.consensus_amino_acids[&(i as u32 + 117)]
            );
        }
    }

    #[test]
    fn terminal_schemes_kappa() {
        let terminal_length = 10;
        let original_scheme = get_imgt_kappa_scheme();
        let (n_term_scheme, c_term_scheme) = get_terminal_schemes(Chain::IGK, terminal_length);

        for i in 0..terminal_length {
            assert_eq!(
                n_term_scheme.scoring_matrix.row(i as Ix),
                original_scheme.scoring_matrix.row(i as Ix)
            );

            assert_eq!(
                c_term_scheme.scoring_matrix.row(i as Ix),
                original_scheme.scoring_matrix.row(117 + i as Ix)
            );
        }

        for i in 1..terminal_length {
            assert_eq!(
                n_term_scheme.consensus_amino_acids[&(i as u32)],
                original_scheme.consensus_amino_acids[&(i as u32)]
            );
            assert_eq!(
                c_term_scheme.consensus_amino_acids[&(i as u32)],
                original_scheme.consensus_amino_acids[&(i as u32 + 117)]
            );
        }
    }

    #[test]
    fn terminal_schemes_lambda() {
        let terminal_length = 10;
        let original_scheme = get_imgt_lambda_scheme();
        let (n_term_scheme, c_term_scheme) = get_terminal_schemes(Chain::IGL, terminal_length);

        for i in 0..terminal_length {
            assert_eq!(
                n_term_scheme.scoring_matrix.row(i as Ix),
                original_scheme.scoring_matrix.row(i as Ix)
            );

            assert_eq!(
                c_term_scheme.scoring_matrix.row(i as Ix),
                original_scheme.scoring_matrix.row(117 + i as Ix)
            );
        }

        for i in 1..terminal_length {
            assert_eq!(
                n_term_scheme.consensus_amino_acids[&(i as u32)],
                original_scheme.consensus_amino_acids[&(i as u32)]
            );
            assert_eq!(
                c_term_scheme.consensus_amino_acids[&(i as u32)],
                original_scheme.consensus_amino_acids[&(i as u32 + 117)]
            );
        }
    }
}
