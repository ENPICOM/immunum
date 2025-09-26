use crate::constants::{get_region_ranges, MINIMAL_CHAIN_IDENTITY, MINIMAL_CHAIN_LENGTH};
use crate::numbering_scheme_type::NumberingScheme;
use crate::prefiltering::apply_prefiltering;
use crate::schemes::get_scheme_with_cdr_definition;
use crate::sequence::SequenceRecord;
use crate::types::{
    CdrDefinition, Chain, ChainNumbering, NumberingPosition, RegionInfo, Regions, Scheme,
};

/// Main annotator struct that consolidates all numbering functionality
pub struct Annotator {
    schemes: Vec<NumberingScheme>,
    disable_prefiltering: bool,
    min_confidence: f64,
}

impl Annotator {
    /// Create a new Annotator with the specified scheme, chains, and optional scoring parameters
    pub fn new(
        scheme: Scheme,
        chains: Vec<Chain>,
        disable_prefiltering: bool,
        min_confidence: Option<f64>,
    ) -> Result<Self, String> {
        Self::new_with_cdr_definition(
            scheme,
            chains,
            CdrDefinition::from_scheme(scheme),
            disable_prefiltering,
            min_confidence,
        )
    }

    /// Create a new Annotator with custom CDR definitions
    pub fn new_with_cdr_definition(
        scheme: Scheme,
        chains: Vec<Chain>,
        cdr_definition: CdrDefinition,
        disable_prefiltering: bool,
        min_confidence: Option<f64>,
    ) -> Result<Self, String> {
        // Pre-build all required schemes for performance
        let schemes: Vec<NumberingScheme> = chains
            .iter()
            .map(|&chain| get_scheme_with_cdr_definition(scheme, chain, cdr_definition))
            .collect();

        if schemes.is_empty() {
            return Err("No valid schemes could be created for the specified chains".to_string());
        }

        Ok(Annotator {
            schemes,
            disable_prefiltering,
            min_confidence: min_confidence.unwrap_or(MINIMAL_CHAIN_IDENTITY),
        })
    }

    /// Numbers a sequence against the (prefiltered) numbering schemes in the annotator class finding at most max_chains
    pub fn number_sequence(
        &self,
        sequence: &SequenceRecord,
        max_chains: Option<usize>,
    ) -> Result<Vec<ChainNumbering>, String> {
        if sequence.sequence.is_empty() {
            return Err("Empty sequence provided".to_string());
        }

        // Apply prefiltering if enabled, otherwise use all schemes
        let schemes: Vec<&NumberingScheme> = if !self.disable_prefiltering {
            apply_prefiltering(&sequence.sequence, &self.schemes)
        } else {
            self.schemes.iter().collect()
        };

        // Use find_all_chains to find multiple chains in the sequence
        find_all_chains(
            &sequence.sequence,
            schemes.as_slice(),
            self.min_confidence,
            max_chains,
        )
    }
}

/// Runs sequence alignment on given schemes, filters the output and return the numbering with the highest confidence
pub fn find_best_chain(
    query_sequence: &[u8],
    numbering_schemes: &[&NumberingScheme],
    min_confidence: f64,
) -> Result<ChainNumbering, String> {
    numbering_schemes
        .iter()
        .map(|&scheme| {
            let (numbers, identity, start, end) = scheme.number_sequence(query_sequence);
            let regions = extract_regions(
                query_sequence,
                &numbers,
                scheme.scheme_type,
                scheme.chain_type,
                scheme.cdr_definition,
                start,
            );
            ChainNumbering {
                numbers,
                identity,
                scheme: scheme.scheme_type,
                chain: scheme.chain_type,
                cdr_definition: scheme.cdr_definition,
                start,
                end,
                regions,
            }
        })
        .filter(|chain| chain.identity > min_confidence)
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
    min_confidence: f64,
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

        let numbering_result = find_best_chain(current_sequence, numbering_schemes, min_confidence);

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
                let mut padded_numbering = vec![NumberingPosition::Gap; current_start];
                padded_numbering.extend(best_chain.numbers);
                padded_numbering.extend(vec![
                    NumberingPosition::Gap;
                    full_query_length - current_end - 1
                ]);

                best_chain.numbers = padded_numbering;

                // Re-extract regions using the full sequence and padded numbering
                best_chain.regions = extract_regions(
                    query_sequence,
                    &best_chain.numbers,
                    best_chain.scheme,
                    best_chain.chain,
                    best_chain.cdr_definition,
                    0, // start from 0 since we're using the full sequence
                );

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

/// Extract region sequences and positions from numbering and sequence
fn extract_regions(
    sequence: &[u8],
    numbers: &[NumberingPosition],
    scheme: Scheme,
    chain: Chain,
    cdr_definition: CdrDefinition,
    sequence_start: usize,
) -> Regions {
    let region_ranges = get_region_ranges(cdr_definition, scheme, chain);

    // Helper function to extract sequence for a region
    let extract_region_sequence = |start_pos: u32, end_pos: u32| -> RegionInfo {
        let mut sequence_start_idx = None;
        let mut sequence_end_idx = None;
        let mut region_sequence = String::new();

        // Find the actual sequence positions corresponding to the numbering positions
        for (seq_idx, number_pos) in numbers.iter().enumerate() {
            let position_matches = match number_pos {
                NumberingPosition::Gap => false,
                NumberingPosition::Number(n) => *n >= start_pos && *n < end_pos,
                NumberingPosition::Insertion { position, .. } => {
                    *position >= start_pos && *position < end_pos
                }
            };

            if position_matches {
                if sequence_start_idx.is_none() {
                    sequence_start_idx = Some(seq_idx);
                }
                sequence_end_idx = Some(seq_idx);

                // Add the amino acid if it's within sequence bounds
                if seq_idx < sequence.len() {
                    region_sequence.push(sequence[seq_idx] as char);
                }
            }
        }

        // Calculate absolute positions in the full sequence
        let abs_start = sequence_start_idx.map_or(0, |idx| sequence_start + idx);
        let abs_end = sequence_end_idx.map_or(abs_start, |idx| sequence_start + idx);

        RegionInfo {
            start: abs_start,
            end: abs_end,
            sequence: region_sequence,
        }
    };

    Regions {
        fr1: extract_region_sequence(region_ranges.fr1.start, region_ranges.fr1.end),
        cdr1: extract_region_sequence(region_ranges.cdr1.start, region_ranges.cdr1.end),
        fr2: extract_region_sequence(region_ranges.fr2.start, region_ranges.fr2.end),
        cdr2: extract_region_sequence(region_ranges.cdr2.start, region_ranges.cdr2.end),
        fr3: extract_region_sequence(region_ranges.fr3.start, region_ranges.fr3.end),
        cdr3: extract_region_sequence(region_ranges.cdr3.start, region_ranges.cdr3.end),
        fr4: extract_region_sequence(region_ranges.fr4.start, region_ranges.fr4.end),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_annotator_creation() {
        let annotator = Annotator::new(Scheme::IMGT, vec![Chain::IGH], false, None);
        assert!(annotator.is_ok());
    }

    #[test]
    fn test_empty_chains() {
        let annotator = Annotator::new(Scheme::IMGT, vec![], false, None);
        assert!(annotator.is_err());
    }

    #[test]
    fn test_light_chain() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            false,
            None,
        )
        .unwrap();

        // Heavy chain sequence that should be correctly identified with prefiltering
        let record = SequenceRecord {
                name: b"light".to_vec(),
                sequence: b"DIVMTQSPDSLAVSLGERATINCKASQSVTNDVAWYQQKPGQPPKLLIYYASNRYTGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQDYSSPYTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGE".to_vec(),
            };

        let annotations = annotator.number_sequence(&record, Some(1));
        assert!(annotations.is_ok());

        if let Ok(annotations) = annotations {
            assert_eq!(annotations.len(), 1);
            assert_eq!(annotations[0].chain, Chain::IGK);
        }
    }

    #[test]
    fn test_prefiltering_with_heavy_chain() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            true,
            None,
        )
        .unwrap();

        // Heavy chain sequence that should be correctly identified with prefiltering
        let heavy_chain_record = SequenceRecord {
                name: b"heavy".to_vec(),
                sequence: b"QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS".to_vec(),
            };

        let annotations = annotator.number_sequence(&heavy_chain_record, Some(1));
        assert!(annotations.is_ok());

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
            None,
        )
        .unwrap();

        // Test sequence with multiple chains (concatenated heavy and light chain)
        let paired_record = SequenceRecord {
                name: b"heavy_light".to_vec(),
                sequence: b"QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSSDIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELK".to_vec()
        };

        let annotations = annotator.number_sequence(&paired_record, Some(2));
        assert!(annotations.is_ok());
        if let Ok(annotations) = annotations {
            // Should find two chains in the paired sequence
            assert_eq!(annotations.len(), 2);
        }
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
            find_best_chain(heavy_chain, &schemes, MINIMAL_CHAIN_IDENTITY)
                .expect("Did not find match")
                .chain,
            Chain::IGH
        );

        assert_eq!(
            find_best_chain(kappa_chain, &schemes, MINIMAL_CHAIN_IDENTITY)
                .expect("Did not find match")
                .chain,
            Chain::IGK
        );

        assert_eq!(
            find_best_chain(lambda_chain, &schemes, MINIMAL_CHAIN_IDENTITY)
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
        let hk_results = find_all_chains(&heavy_kappa, &schemes, MINIMAL_CHAIN_IDENTITY, Some(2));
        assert_eq!(hk_results.unwrap().len(), 2);

        let heavy_lambda = [heavy_chain, linker, lambda_chain].concat();
        let hl_results = find_all_chains(&heavy_lambda, &schemes, MINIMAL_CHAIN_IDENTITY, Some(2));
        assert_eq!(hl_results.unwrap().len(), 2);

        let kappa_lambda = [kappa_chain, linker, lambda_chain].concat();
        let kl_results = find_all_chains(&kappa_lambda, &schemes, MINIMAL_CHAIN_IDENTITY, Some(2));
        assert_eq!(kl_results.unwrap().len(), 2);

        let heavy_heavy = [heavy_chain, linker, heavy_chain].concat();
        let hh_results = find_all_chains(&heavy_heavy, &schemes, MINIMAL_CHAIN_IDENTITY, Some(2));
        assert_eq!(hh_results.unwrap().len(), 2);
    }
}
