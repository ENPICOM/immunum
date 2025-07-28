use crate::constants::{PRE_FILTER_TERMINAL_LENGTH, WITHIN_IDENTITY_RANGE};
use crate::numbering_scheme_type::{NumberingOutput, NumberingScheme};
use crate::types::{Chain, PrefilterOutput};
use std::collections::HashMap;

/// Create c-terminal NumberingScheme from a full length scheme (based on IMGT consensus)
pub fn get_terminal_schemes(
    schemes: &Vec<NumberingScheme>,
) -> Vec<(NumberingScheme, NumberingScheme)> {
    let mut terminal_schemes: Vec<(NumberingScheme, NumberingScheme)> = Vec::new();
    for scheme in schemes {
        terminal_schemes.push(scheme.to_terminal_schemes(PRE_FILTER_TERMINAL_LENGTH));
    }
    terminal_schemes
}

///Run pre scan to find likely chains present using c and n terminal identity
pub fn run_pre_scan(
    query_sequence: &[u8],
    all_terminal_schemes: &Vec<(NumberingScheme, NumberingScheme)>,
) -> (HashMap<Chain, PrefilterOutput>, f64) {
    let mut chain_identity_map: HashMap<Chain, PrefilterOutput> = HashMap::new();

    let mut highest_identity: f64 = 0.0;

    for terminal_schemes in all_terminal_schemes {
        let (n_terminal, c_terminal) = terminal_schemes;
        // run alignment for n and c terminal
        let n_terminal_output: NumberingOutput = n_terminal.number_sequence(query_sequence);
        let c_terminal_output: NumberingOutput = c_terminal.number_sequence(query_sequence);

        // calculate combined identity

        let combined_identity: f64 =
            (n_terminal_output.identity + c_terminal_output.identity) / 2.0;

        // Store if best match
        if combined_identity > highest_identity {
            highest_identity = combined_identity;
        }

        // Store score and predicted end and start positions
        chain_identity_map.insert(
            n_terminal.chain_type,
            PrefilterOutput {
                identity: combined_identity,
                predicted_start: n_terminal_output.start,
                predicted_end: c_terminal_output.end,
            },
        );
    }
    (chain_identity_map, highest_identity)
}

pub fn select_chains_from_pre_scan(
    pre_scan_output: &HashMap<Chain, PrefilterOutput>,
    highest_score: f64,
) -> Vec<Chain> {
    pre_scan_output
        .iter()
        .filter(|(_, output)| output.identity >= highest_score - WITHIN_IDENTITY_RANGE)
        .map(|(chain, _)| *chain)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::schemes::{get_imgt_heavy_scheme, get_imgt_kappa_scheme, get_imgt_lambda_scheme};
    use ndarray::Ix;

    #[test]
    fn heavy_chain_pre_scan() {
        let sequence_h = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS".as_bytes();
        let sequence_l = "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTRVFGTGTKVTVL".as_bytes();
        let sequence_k = "DIQMTQSPSSLSASVGDRVTITCRASQSISSWLAWYQQKPGKAPKLLIYKASSLESGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYNSYPFTFGQGTKVEIK".as_bytes();

        let terminal_schemes = get_terminal_schemes(&vec![
            get_imgt_heavy_scheme(),
            get_imgt_lambda_scheme(),
            get_imgt_kappa_scheme(),
        ]);

        let (pre_scan_output_h, _highest_score_h) = run_pre_scan(sequence_h, &terminal_schemes);
        let (pre_scan_output_l, _highest_score_l) = run_pre_scan(sequence_l, &terminal_schemes);
        let (pre_scan_output_k, _highest_score_k) = run_pre_scan(sequence_k, &terminal_schemes);

        // check if correct chain has highest identity
        assert!(
            (pre_scan_output_h[&Chain::IGH].identity > pre_scan_output_h[&Chain::IGK].identity)
                & (pre_scan_output_h[&Chain::IGH].identity
                    > pre_scan_output_h[&Chain::IGL].identity)
        );
        assert!(
            (pre_scan_output_l[&Chain::IGL].identity > pre_scan_output_l[&Chain::IGK].identity)
                & (pre_scan_output_l[&Chain::IGL].identity
                    > pre_scan_output_l[&Chain::IGH].identity)
        );
        assert!(
            (pre_scan_output_k[&Chain::IGK].identity > pre_scan_output_k[&Chain::IGH].identity)
                & (pre_scan_output_k[&Chain::IGK].identity
                    > pre_scan_output_k[&Chain::IGL].identity)
        );
    }

    #[test]
    fn terminal_schemes_heavy() {
        let terminal_length = 10;
        let original_scheme = get_imgt_heavy_scheme();
        let (n_term_scheme, c_term_scheme) =
            &get_terminal_schemes(&vec![get_imgt_heavy_scheme()])[0];
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
        let (n_term_scheme, c_term_scheme) =
            &get_terminal_schemes(&vec![get_imgt_kappa_scheme()])[0];
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
        let (n_term_scheme, c_term_scheme) =
            &get_terminal_schemes(&vec![get_imgt_lambda_scheme()])[0];
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
