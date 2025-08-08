use crate::constants::{PRE_FILTER_TERMINAL_LENGTH, WITHIN_IDENTITY_RANGE};
use crate::numbering_scheme::NumberingScheme;
use crate::result::AnnotationResult;
use crate::types::{Chain, PrefilterOutput};
use std::collections::HashMap;

/// Apply prefiltering to select only promising chains based on terminal identity
pub fn prefilter_schemes<'a>(
    sequence: &'a [u8],
    schemes: &'a [NumberingScheme],
) -> Vec<&'a NumberingScheme> {
    if schemes.is_empty() {
        return vec![];
    }

    // get pre-filter schemes
    let terminal_schemes = schemes
        .iter()
        .map(|scheme| scheme.to_terminal_schemes(PRE_FILTER_TERMINAL_LENGTH))
        .collect();

    // Select models to run using terminal numbering
    let (terminal_number_output, highest_score) =
        terminal_number_sequence(sequence, &terminal_schemes);
    let best_chains = select_best_chains(&terminal_number_output, highest_score);

    // Filter schemes based on prefiltering results
    schemes
        .iter()
        .filter(|scheme| best_chains.contains(&scheme.chain_type))
        .collect()
}

///Run pre scan to find likely chains present using c and n terminal identity
pub fn terminal_number_sequence(
    query_sequence: &[u8],
    terminal_schemes: &Vec<(NumberingScheme, NumberingScheme)>,
) -> (HashMap<Chain, PrefilterOutput>, f64) {
    let mut chain_identity_map: HashMap<Chain, PrefilterOutput> = HashMap::new();

    let mut highest_identity: f64 = 0.0;

    for (n_terminal, c_terminal) in terminal_schemes {
        // run alignment for n and c terminal
        let n_terminal_output: AnnotationResult = n_terminal.number_sequence(query_sequence);
        let c_terminal_output: AnnotationResult = c_terminal.number_sequence(query_sequence);

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
                _predicted_start: n_terminal_output.start,
                _predicted_end: c_terminal_output.end,
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
        schemes::get_scheme,
        types::{Chain, Scheme},
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

        let filtered_schemes = prefilter_schemes(heavy_chain, &schemes);

        // Should return at least one scheme (preferably the heavy chain one)
        assert!(!filtered_schemes.is_empty());
        // Heavy chain should be among the filtered schemes since it's the best match
        assert!(filtered_schemes.iter().any(|s| s.chain_type == Chain::IGH));
    }

    #[test]
    fn test_apply_prefiltering_empty_schemes() {
        let sequence: &[u8] = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS".as_bytes();
        let schemes: Vec<NumberingScheme> = vec![];

        let filtered_schemes = prefilter_schemes(sequence, &schemes);
        assert!(filtered_schemes.is_empty());
    }
}
