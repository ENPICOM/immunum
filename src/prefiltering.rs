use crate::constants::WITHIN_IDENTITY_RANGE;
use crate::needleman_wunsch::{needleman_wunsch_consensus, MatrixPool};
use crate::numbering_scheme::NumberingScheme;
// use crate::result::AnnotationResult; // no longer needed here
use crate::types::{Chain, PrefilterOutput};
use std::collections::HashMap;

/// Prefilter using precomputed terminal schemes (cached) and thread-local matrix pool.
/// Returns the list of promising chains.
pub fn prefilter_schemes(
    sequence: &[u8],
    terminal_schemes: &[(NumberingScheme, NumberingScheme)],
) -> Vec<Chain> {
    let (terminal_number_output, highest_score) =
        terminal_number_sequence(sequence, terminal_schemes);
    select_best_chains(&terminal_number_output, highest_score)
}

/// Faster terminal scan using thread-local matrix pool and avoiding full AnnotationResult work.
pub fn terminal_number_sequence(
    query_sequence: &[u8],
    terminal_schemes: &[(NumberingScheme, NumberingScheme)],
) -> (HashMap<Chain, PrefilterOutput>, f64) {
    let mut chain_identity_map: HashMap<Chain, PrefilterOutput> = HashMap::new();
    let mut highest_identity: f64 = 0.0;

    thread_local! {
        static MATRIX_POOL: std::cell::RefCell<MatrixPool> = std::cell::RefCell::new(MatrixPool::new());
    }

    for (n_terminal, c_terminal) in terminal_schemes.iter() {
        let (n_identity, n_start) = MATRIX_POOL.with(|pool_cell| {
            let mut pool = pool_cell.borrow_mut();
            let (numbering, identity) =
                needleman_wunsch_consensus(query_sequence, n_terminal, &mut pool, -50.0);
            let start = numbering.iter().position(|s| s != "-").unwrap_or(0) as u32;
            (identity, start)
        });

        let (c_identity, c_end) = MATRIX_POOL.with(|pool_cell| {
            let mut pool = pool_cell.borrow_mut();
            let (numbering, identity) =
                needleman_wunsch_consensus(query_sequence, c_terminal, &mut pool, -50.0);
            let end = numbering
                .iter()
                .rposition(|s| s != "-")
                .unwrap_or_else(|| numbering.len().saturating_sub(1)) as u32;
            (identity, end)
        });

        let combined_identity = (n_identity + c_identity) / 2.0;
        if combined_identity > highest_identity {
            highest_identity = combined_identity;
        }

        chain_identity_map.insert(
            n_terminal.chain_type,
            PrefilterOutput {
                identity: combined_identity,
                _predicted_start: n_start,
                _predicted_end: c_end,
            },
        );
    }

    (chain_identity_map, highest_identity)
}

pub fn select_best_chains(
    terminal_numbering_output: &HashMap<Chain, PrefilterOutput>,
    highest_score: f64,
) -> Vec<Chain> {
    terminal_numbering_output
        .iter()
        .filter(|(_, output)| output.identity >= highest_score - WITHIN_IDENTITY_RANGE)
        .map(|(chain, _)| *chain)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        constants::PRE_FILTER_TERMINAL_LENGTH, schemes::get_scheme, types::{Chain, Scheme}
    };

    #[test]
    fn heavy_chain_pre_scan() {
        let sequence_h = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS".as_bytes();
        let sequence_l = "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTRVFGTGTKVTVL".as_bytes();
        let sequence_k = "DIQMTQSPSSLSASVGDRVTITCRASQSISSWLAWYQQKPGKAPKLLIYKASSLESGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYNSYPFTFGQGTKVEIK".as_bytes();

        let schemes = vec![
            get_scheme(Scheme::IMGT, Chain::IGH, None),
            get_scheme(Scheme::IMGT, Chain::IGL, None),
            get_scheme(Scheme::IMGT, Chain::IGK, None),
        ];
        let terminal_schemes: Vec<(NumberingScheme, NumberingScheme)> = schemes
            .iter()
            .map(|scheme| scheme.to_terminal_schemes(PRE_FILTER_TERMINAL_LENGTH))
            .collect();

        let (pre_scan_output_h, _highest_score_h) =
            terminal_number_sequence(sequence_h, &terminal_schemes);
        let (pre_scan_output_l, _highest_score_l) =
            terminal_number_sequence(sequence_l, &terminal_schemes);
        let (pre_scan_output_k, _highest_score_k) =
            terminal_number_sequence(sequence_k, &terminal_schemes);

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
    fn test_apply_prefiltering() {
        let heavy_chain: &[u8] = "QVQLVQSGAVIKTPGSSVKISCRASGYNFRDYSIHWVRLIPDKGFEWIGWIKPLWGAVSYARQL\
        QGRVSMTRQLSQDPDDPDWGVAYMEFSGLTPADTAEYFCVRRGSCDYCGDFPWQYWCQGTVVVVSSASTKGPSVFPLAPSSGGTAALGCLV\
        KDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPK"
            .as_bytes();

        let schemes = vec![
            get_scheme(Scheme::IMGT, Chain::IGH, None),
            get_scheme(Scheme::IMGT, Chain::IGK, None),
            get_scheme(Scheme::IMGT, Chain::IGL, None),
        ];

        let terminal_schemes: Vec<(NumberingScheme, NumberingScheme)> = schemes
            .iter()
            .map(|scheme| scheme.to_terminal_schemes(PRE_FILTER_TERMINAL_LENGTH))
            .collect();

        let filtered_chains = prefilter_schemes(heavy_chain, &terminal_schemes);

        // Should return at least one chain (preferably heavy)
        assert!(!filtered_chains.is_empty());
        assert!(filtered_chains.contains(&Chain::IGH));
    }

    #[test]
    fn test_apply_prefiltering_empty_schemes() {
        let sequence: &[u8] = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS".as_bytes();
        let terminal_schemes: Vec<(NumberingScheme, NumberingScheme)> = vec![];

        let filtered = prefilter_schemes(sequence, &terminal_schemes);
        assert!(filtered.is_empty());
    }
}
