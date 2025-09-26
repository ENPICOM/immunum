use crate::constants::{
    get_region_ranges, MINIMAL_CHAIN_IDENTITY, MINIMAL_CHAIN_LENGTH, PRE_FILTER_TERMINAL_LENGTH,
};
use crate::numbering_scheme_type::NumberingScheme;
use crate::prefiltering::apply_prefiltering;
use crate::schemes::get_scheme_with_cdr_definition;
use crate::sequence::SequenceRecord;
use crate::types::{
    CdrDefinitions, Chain, ChainNumbering, NumberingPosition, RegionInfo, Regions, Scheme,
};

/// Main annotator struct that consolidates all numbering functionality
pub struct Annotator {
    schemes: Vec<NumberingScheme>,
    disable_prefiltering: bool,
    min_confidence: f64,
    // Cache terminal schemes for prefiltering performance
    terminal_schemes: Vec<(NumberingScheme, NumberingScheme)>,
}

impl Annotator {
    /// Create a new Annotator with the specified scheme, chains, and optional scoring parameters
    pub fn new(
        scheme: Scheme,
        chains: Vec<Chain>,
        cdr_definitions: Option<CdrDefinitions>,
        disable_prefiltering: bool,
        min_confidence: Option<f64>,
    ) -> Result<Self, String> {
        let cdr_definitions = cdr_definitions.unwrap_or(CdrDefinitions::from_scheme(scheme));
        // Pre-build all required schemes for performance
        let schemes: Vec<NumberingScheme> = chains
            .iter()
            .map(|&chain| get_scheme_with_cdr_definition(scheme, chain, cdr_definitions))
            .collect();

        if schemes.is_empty() {
            return Err("No valid schemes could be created for the specified chains".to_string());
        }

        // Pre-build terminal schemes for prefiltering performance
        let terminal_schemes: Vec<(NumberingScheme, NumberingScheme)> = schemes
            .iter()
            .map(|scheme| scheme.to_terminal_schemes(PRE_FILTER_TERMINAL_LENGTH))
            .collect();

        Ok(Annotator {
            schemes,
            disable_prefiltering,
            min_confidence: min_confidence.unwrap_or(MINIMAL_CHAIN_IDENTITY),
            terminal_schemes,
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
            apply_prefiltering(&sequence.sequence, &self.schemes, &self.terminal_schemes)
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

/// Extract region positions from numbering (optimized - positions only, no sequence extraction)
fn extract_regions(
    numbers: &[NumberingPosition],
    scheme: Scheme,
    chain: Chain,
    cdr_definition: CdrDefinitions,
    sequence_start: usize,
) -> Regions {
    let region_ranges = get_region_ranges(cdr_definition, scheme, chain);

    // Pre-allocate region boundary trackers
    let mut region_boundaries = [None; 14]; // 7 regions * 2 (start, end)

    // Single pass through numbering positions to find all region boundaries
    for (seq_idx, number_pos) in numbers.iter().enumerate() {
        let position = match number_pos {
            NumberingPosition::Gap => continue,
            NumberingPosition::Number(n) => *n,
            NumberingPosition::Insertion { position, .. } => *position,
        };

        // Check each region and update boundaries if this position falls within it
        let regions = [
            (region_ranges.fr1.start, region_ranges.fr1.end, 0),
            (region_ranges.cdr1.start, region_ranges.cdr1.end, 2),
            (region_ranges.fr2.start, region_ranges.fr2.end, 4),
            (region_ranges.cdr2.start, region_ranges.cdr2.end, 6),
            (region_ranges.fr3.start, region_ranges.fr3.end, 8),
            (region_ranges.cdr3.start, region_ranges.cdr3.end, 10),
            (region_ranges.fr4.start, region_ranges.fr4.end, 12),
        ];

        for &(start_pos, end_pos, boundary_idx) in &regions {
            if position >= start_pos && position < end_pos {
                let start_idx = boundary_idx;
                let end_idx = boundary_idx + 1;

                // Update start boundary (first occurrence)
                if region_boundaries[start_idx].is_none() {
                    region_boundaries[start_idx] = Some(sequence_start + seq_idx);
                }
                // Always update end boundary (last occurrence)
                region_boundaries[end_idx] = Some(sequence_start + seq_idx);
            }
        }
    }

    // Helper to create RegionInfo from boundary indices
    let create_region = |start_idx: usize, end_idx: usize| -> RegionInfo {
        RegionInfo {
            start: region_boundaries[start_idx].unwrap_or(0),
            end: region_boundaries[end_idx].unwrap_or(0),
        }
    };

    Regions {
        fr1: create_region(0, 1),
        cdr1: create_region(2, 3),
        fr2: create_region(4, 5),
        cdr2: create_region(6, 7),
        fr3: create_region(8, 9),
        cdr3: create_region(10, 11),
        fr4: create_region(12, 13),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::schemes::get_scheme_with_cdr_definition;

    #[test]
    fn test_annotator_creation() {
        let annotator = Annotator::new(Scheme::IMGT, vec![Chain::IGH], None, false, None);
        assert!(annotator.is_ok());
    }

    #[test]
    fn test_empty_chains() {
        let annotator = Annotator::new(Scheme::IMGT, vec![], None, false, None);
        assert!(annotator.is_err());
    }

    #[test]
    fn test_light_chain() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            None,
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
            None,
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
            None,
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

        let heavy_scheme =
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGH, CdrDefinitions::IMGT);
        let lambda_scheme =
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGL, CdrDefinitions::IMGT);
        let kappa_scheme =
            get_scheme_with_cdr_definition(Scheme::KABAT, Chain::IGK, CdrDefinitions::KABAT);
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

        let heavy_scheme =
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGH, CdrDefinitions::IMGT);
        let lambda_scheme =
            get_scheme_with_cdr_definition(Scheme::IMGT, Chain::IGL, CdrDefinitions::IMGT);
        let kappa_scheme =
            get_scheme_with_cdr_definition(Scheme::KABAT, Chain::IGK, CdrDefinitions::KABAT);
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
