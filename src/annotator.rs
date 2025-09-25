use crate::constants::{MINIMAL_CHAIN_IDENTITY, MINIMAL_CHAIN_LENGTH};
use crate::numbering_scheme_type::NumberingScheme;
use crate::prefiltering::apply_prefiltering;
use crate::schemes::get_scheme;
use crate::sequence::SequenceRecord;
use crate::types::{Chain, ChainNumbering, Scheme};
use rayon::prelude::*;

/// Main annotator struct that consolidates all numbering functionality
pub struct Annotator {
    schemes: Vec<NumberingScheme>,
    use_prefiltering: bool,
}

impl Annotator {
    /// Create a new Annotator with the specified scheme, chains, and optional scoring parameters
    pub fn new(scheme: Scheme, chains: Vec<Chain>, use_prefiltering: bool) -> Result<Self, String> {
        // Pre-build all required schemes for performance
        let schemes: Vec<NumberingScheme> = chains
            .iter()
            .map(|&chain| get_scheme(scheme, chain))
            .collect();

        if schemes.is_empty() {
            return Err("No valid schemes could be created for the specified chains".to_string());
        }

        Ok(Annotator {
            schemes,
            use_prefiltering,
        })
    }

    /// Numbers sequences against the setup numbering schemes in the annotator class.  
    /// Returns a Vec with the same order as the input sequences, for each the name a Result with vec of ChainNumbering.
    /// This is the only public class on Annotator
    pub fn number(
        &self,
        sequences: Vec<SequenceRecord>,
        max_chains: Option<usize>,
    ) -> Vec<(String, Result<Vec<ChainNumbering>, String>)> {
        sequences
            .into_par_iter()
            .map(|record| {
                let name = String::from_utf8_lossy(&record.name).to_string();
                match self.number_sequence(&record, max_chains) {
                    Ok(results) => (name, Ok(results)),
                    Err(e) => (name, Err(e)),
                }
            })
            .collect()
    }

    fn number_sequence(
        &self,
        sequence: &SequenceRecord,
        max_chains: Option<usize>,
    ) -> Result<Vec<ChainNumbering>, String> {
        if sequence.sequence.is_empty() {
            return Err("Empty sequence provided".to_string());
        }

        // Apply prefiltering if enabled, otherwise use all schemes
        let schemes: Vec<&NumberingScheme> = if self.use_prefiltering {
            apply_prefiltering(&sequence.sequence, &self.schemes)
        } else {
            self.schemes.iter().collect()
        };

        // Use find_all_chains to find multiple chains in the sequence
        find_all_chains(&sequence.sequence, schemes.as_slice(), max_chains)
    }

}

/// Runs sequence alignment on given schemes, filters the output and return the numbering with the highest confidence
pub fn find_best_chain(
    query_sequence: &[u8],
    numbering_schemes: &[&NumberingScheme],
) -> Result<ChainNumbering, String> {
    numbering_schemes
        .iter()
        .map(|&scheme| {
            let (numbers, identity, start, end) = scheme.number_sequence(query_sequence);
            ChainNumbering {
                numbers,
                identity,
                scheme: scheme.scheme_type,
                chain: scheme.chain_type,
                start,
                end,
            }
        })
        .filter(|chain| chain.identity > MINIMAL_CHAIN_IDENTITY)
        .max_by(|numbers_a, numbers_b| {
            numbers_a
                .identity
                .partial_cmp(&numbers_b.identity)
                // This can only fail when an identity is NaN, which should never happen
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .ok_or_else(|| "No valid numbering found".to_string())
}

/// Attempts to find all antibody chains in a sequence
fn find_all_chains(
    query_sequence: &[u8],
    numbering_schemes: &[&NumberingScheme],
    max_chains: Option<usize>,
) -> Result<Vec<ChainNumbering>, String> {
    let max_chains = max_chains.unwrap_or(99); // Default is "unlimited"

    let mut chains_found: Vec<ChainNumbering> = Vec::new();

    let full_query_length = query_sequence.len();
    let mut sequence_list: Vec<(&[u8], usize, usize)> = Vec::new();
    sequence_list.push((query_sequence, 0, full_query_length - 1));

    while let Some((current_sequence, current_start, current_end)) = sequence_list.pop() {
        // Early termination if we've found enough chains
        if chains_found.len() >= max_chains {
            break;
        }

        let numbering_result = find_best_chain(current_sequence, numbering_schemes);

        // Break early for single chain optimization
        if max_chains == 1 {
            return match numbering_result {
                Ok(best_chain) => Ok(vec![best_chain]),
                Err(e) => Err(e.to_string()),
            };
        }

        match numbering_result {
            Ok(mut best_chain) => {
                // split the remaining sequence, add sequences that are long enough to the list
                let front_sequence: &[u8] = &current_sequence[0..best_chain.start];
                let end_sequence: &[u8] = &current_sequence[(best_chain.end + 1)..];

                if front_sequence.len() > MINIMAL_CHAIN_LENGTH {
                    sequence_list.push((
                        front_sequence,
                        current_start,
                        (current_start + best_chain.start - 1),
                    ))
                }
                if end_sequence.len() > MINIMAL_CHAIN_LENGTH {
                    sequence_list.push((
                        end_sequence,
                        (current_start + best_chain.end + 1),
                        current_end,
                    ))
                }

                // set sequence to full original sequence
                // best_chain.sequence = query_sequence.to_vec();
                best_chain.start += current_start;
                best_chain.end += current_start;

                // add gaps to front and end to match with length of original sequence
                let mut start_padding = vec![String::from("-"); current_start];
                let end_padding = vec![String::from("-"); full_query_length - current_end - 1];
                // TODO improve naming
                start_padding.extend(best_chain.numbers);
                start_padding.extend(end_padding);

                best_chain.numbers = start_padding;
                chains_found.push(best_chain);
            }
            Err(_) => {
                // Continue processing other segments even if current segment doesn't match
                continue;
            }
        }
    }

    Ok(chains_found)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_annotator_creation() {
        let annotator = Annotator::new(Scheme::IMGT, vec![Chain::IGH], false);
        assert!(annotator.is_ok());
    }

    #[test]
    fn test_empty_chains() {
        let annotator = Annotator::new(Scheme::IMGT, vec![], false);
        assert!(annotator.is_err());
    }

    #[test]
    fn test_light_chain() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            false,
        )
        .unwrap();

        // Heavy chain sequence that should be correctly identified with prefiltering
        let record = SequenceRecord {
                name: b"light".to_vec(),
                sequence: b"DIVMTQSPDSLAVSLGERATINCKASQSVTNDVAWYQQKPGQPPKLLIYYASNRYTGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQDYSSPYTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGE".to_vec(),
            };

        let mut results = annotator.number(vec![record], Some(1));
        assert_eq!(results.len(), 1);

        let (name, annotations) = results.pop().unwrap();
        assert!(annotations.is_ok());
        assert_eq!(name, "light");

        if let Ok(annotations) = annotations {
            assert_eq!(annotations.len(), 1);
            assert_eq!(annotations[0].chain, Chain::IGK);
        }
    }

    #[test]
    fn test_prefiltering_with_heavy_chain() {
        let annotator =
            Annotator::new(Scheme::IMGT, vec![Chain::IGH, Chain::IGK, Chain::IGL], true).unwrap();

        // Heavy chain sequence that should be correctly identified with prefiltering
        let heavy_chain_record = SequenceRecord {
                name: b"heavy".to_vec(),
                sequence: b"QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS".to_vec(),
            };

        let mut results = annotator.number(vec![heavy_chain_record], Some(1));
        assert_eq!(results.len(), 1);

        let (name, annotations) = results.pop().unwrap();
        assert!(annotations.is_ok());
        assert_eq!(name, "heavy");

        if let Ok(annotations) = annotations {
            assert_eq!(annotations.len(), 1);
            assert_eq!(annotations[0].chain, Chain::IGH);
        }
    }

    #[test]
    fn test_number_paired_sequence() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            false, // Disable prefiltering for this test
        )
        .unwrap();

        // Test sequence with multiple chains (concatenated heavy and light chain)
        let paired_record = SequenceRecord {
                name: b"heavy_light".to_vec(),
                sequence: b"QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSSDIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELK".to_vec()
        };

        let mut results = annotator.number(vec![paired_record], Some(2));
        assert_eq!(results.len(), 1);

        let (_, annotations) = results.pop().unwrap();
        assert!(annotations.is_ok());
        if let Ok(annotations) = annotations {
            // Should find two chains in the paired sequence
            assert_eq!(annotations.len(), 2);
        }
    }

    #[test]
    fn test_number_sequences() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            false,
        )
        .unwrap();

        let records = vec![
            SequenceRecord::new("seq_1".to_string(), "QVQLVQSGAEVKKPGASVKVSCKAS".to_string()),
            SequenceRecord::new("seq_2".to_string(), "DIQMTQSPSSLSASVGDRVTITC".to_string()),
            SequenceRecord::new("seq_3".to_string(), "EVQLLESGGGLVQPGGSLRLSCAAS".to_string()),
        ];

        let results = annotator.number(records, Some(1));
        assert_eq!(results.len(), 3);
    }

    #[test]
    fn test_correct_chain_identification() {
        let heavy_chain: &[u8] = "QVQLVQSGAVIKTPGSSVKISCRASGYNFRDYSIHWVRLIPDKGFEWIGWIKPLWGAVSYARQL\
        QGRVSMTRQLSQDPDDPDWGVAYMEFSGLTPADTAEYFCVRRGSCDYCGDFPWQYWCQGTVVVVSSASTKGPSVFPLAPSSGGTAALGCLV\
        KDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPK"
            .as_bytes();
        let lambda_chain: &[u8] = "SALTQPPSASGSLGQSVTISCTGTSSDVGGYNYVSWYQQHAGKAPKVIIYEVNKRPSGVPDRF\
        SGSKSGNTASLTVSGLQAEDEADYYCSSYEGSDNFVFGTGTKVTVLGQPKANPTVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWK\
        ADGSPVKAGVETTKPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS"
            .as_bytes();
        let kappa_chain: &[u8] = "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTG\
        SGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSE\
        RQNGVLNSATDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC"
            .as_bytes();

        let heavy_scheme = get_scheme(Scheme::IMGT, Chain::IGH);
        let lambda_scheme = get_scheme(Scheme::IMGT, Chain::IGL);
        let kappa_scheme = get_scheme(Scheme::KABAT, Chain::IGK);
        let schemes: Vec<&NumberingScheme> = vec![&heavy_scheme, &lambda_scheme, &kappa_scheme];

        assert_eq!(
            find_best_chain(heavy_chain, &schemes)
                .expect("Did not find match")
                .chain,
            Chain::IGH
        );

        assert_eq!(
            find_best_chain(kappa_chain, &schemes)
                .expect("Did not find match")
                .chain,
            Chain::IGK
        );

        assert_eq!(
            find_best_chain(lambda_chain, &schemes)
                .expect("Did not find match")
                .chain,
            Chain::IGL
        );
    }

    #[test]
    fn test_correct_paired_chain_identification() {
        let heavy_chain: &[u8] = "QVQLVQSGAVIKTPGSSVKISCRASGYNFRDYSIHWVRLIPDKGFEWIGWIKPLWGAVSYARQL\
        QGRVSMTRQLSQDPDDPDWGVAYMEFSGLTPADTAEYFCVRRGSCDYCGDFPWQYWCQGTVVVVSSASTKGPSVFPLAPSSGGTAALGCLV\
        KDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPK"
            .as_bytes();
        let lambda_chain: &[u8] = "SALTQPPSASGSLGQSVTISCTGTSSDVGGYNYVSWYQQHAGKAPKVIIYEVNKRPSGVPDRF\
        SGSKSGNTASLTVSGLQAEDEADYYCSSYEGSDNFVFGTGTKVTVLGQPKANPTVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWK\
        ADGSPVKAGVETTKPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS"
            .as_bytes();
        let kappa_chain: &[u8] = "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTG\
        SGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSE\
        RQNGVLNSATDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC"
            .as_bytes();
        let linker: &[u8] = "GGGGGGG".as_bytes();

        let heavy_scheme = get_scheme(Scheme::IMGT, Chain::IGH);
        let lambda_scheme = get_scheme(Scheme::IMGT, Chain::IGL);
        let kappa_scheme = get_scheme(Scheme::KABAT, Chain::IGK);
        let schemes: Vec<&NumberingScheme> = vec![&heavy_scheme, &lambda_scheme, &kappa_scheme];

        let heavy_kappa = [heavy_chain, linker, kappa_chain].concat();
        let hk_results = find_all_chains(&heavy_kappa, &schemes, Some(2));
        assert_eq!(hk_results.unwrap().len(), 2);

        let heavy_lambda = [heavy_chain, linker, lambda_chain].concat();
        let hl_results = find_all_chains(&heavy_lambda, &schemes, Some(2));
        assert_eq!(hl_results.unwrap().len(), 2);

        let kappa_lambda = [kappa_chain, linker, lambda_chain].concat();
        let kl_results = find_all_chains(&kappa_lambda, &schemes, Some(2));
        assert_eq!(kl_results.unwrap().len(), 2);

        let heavy_heavy = [heavy_chain, linker, heavy_chain].concat();
        let hh_results = find_all_chains(&heavy_heavy, &schemes, Some(2));
        assert_eq!(hh_results.unwrap().len(), 2);
    }
}
