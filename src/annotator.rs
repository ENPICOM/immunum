use crate::constants::{get_scoring_params, ScoringParams, MINIMAL_CHAIN_IDENTITY, MINIMAL_CHAIN_LENGTH};
use crate::numbering_scheme_type::NumberingScheme;
use crate::prefiltering::apply_prefiltering;
use crate::result::AnnotationResult;
use crate::schemes::get_scheme;
use crate::sequence::{from_path, SequenceRecord};
use crate::types::{Chain, Scheme};
use rayon::prelude::*;

/// Type alias for complex file processing results
type FileProcessingResult = Result<Vec<(String, Vec<Result<AnnotationResult, String>>)>, String>;

/// Main annotator struct that consolidates all numbering functionality
pub struct Annotator {
    // TODO! These fields are currently not used, we might remove them in the future
    _scheme: Scheme,
    _chains: Vec<Chain>,
    _scoring_params: ScoringParams,
    schemes: Vec<NumberingScheme>,
    use_prefiltering: bool,
}

impl Annotator {
    /// Create a new Annotator with the specified scheme, chains, and optional scoring parameters
    pub fn new(
        scheme: Scheme,
        chains: Vec<Chain>,
        scoring_params: Option<ScoringParams>,
        use_prefiltering: Option<bool>,
    ) -> Result<Self, String> {
        let params = scoring_params.unwrap_or_else(get_scoring_params);
        let enable_prefiltering = use_prefiltering.unwrap_or(false);

        // Pre-build all required schemes for performance
        let schemes: Vec<NumberingScheme> = chains
            .iter()
            .map(|&chain| get_scheme(scheme, chain, Some(params.clone())))
            .collect();

        if schemes.is_empty() {
            return Err("No valid schemes could be created for the specified chains".to_string());
        }

        Ok(Annotator {
            _scheme: scheme,
            _chains: chains,
            _scoring_params: params,
            schemes,
            use_prefiltering: enable_prefiltering,
        })
    }

    /// Number a paired sequence using the pre-configured schemes to find multiple chains
    pub fn number_paired_sequence(
        &self,
        sequence: &str,
        sequence_id: String,
    ) -> Vec<Result<AnnotationResult, String>> {
        if sequence.is_empty() {
            return vec![Err("Empty sequence provided".to_string())];
        }

        let sequence_bytes = sequence.as_bytes();

        // Apply prefiltering if enabled, otherwise use all schemes
        let scheme_refs: Vec<&NumberingScheme> = if self.use_prefiltering {
            apply_prefiltering(sequence_bytes, &self.schemes)
        } else {
            self.schemes.iter().collect()
        };

        // Use find_all_chains to find multiple chains in the sequence
        self.find_all_chains(sequence_bytes, scheme_refs, sequence_id)
            .into_iter()
            .map(|result| result.map_err(|e| e.to_string()))
            .collect()
    }

    /// Number a single sequence using the pre-configured schemes
    pub fn number_sequence(
        &self,
        sequence: &str,
        sequence_id: String,
    ) -> Result<AnnotationResult, String> {
        if sequence.is_empty() {
            return Err("Empty sequence provided".to_string());
        }

        let sequence_bytes = sequence.as_bytes();

        // Apply prefiltering if enabled, otherwise use all schemes
        let scheme_refs: Vec<&NumberingScheme> = if self.use_prefiltering {
            apply_prefiltering(sequence_bytes, &self.schemes)
        } else {
            self.schemes.iter().collect()
        };

        match self.find_highest_identity_chain(sequence_bytes, &scheme_refs, sequence_id) {
            Ok(output) => Ok(output),
            Err(e) => Err(e.to_string()),
        }
    }

    /// Number multiple sequences, optionally in parallel
    pub fn number_sequences(
        &self,
        sequences: &[String],
        parallel: bool,
    ) -> Vec<Result<AnnotationResult, String>> {
        if parallel {
            // Use rayon for parallel processing
            sequences
                .par_iter()
                .enumerate()
                .map(|(i, seq)| self.number_sequence(seq, format!("sequence_{}", i + 1)))
                .collect()
        } else {
            sequences
                .iter()
                .enumerate()
                .map(|(i, seq)| self.number_sequence(seq, format!("sequence_{}", i + 1)))
                .collect()
        }
    }

    /// Process sequences from a FASTA/FASTQ file and return results with sequence names
    pub fn number_file(&self, file_path: &str, parallel: bool) -> FileProcessingResult {
        // Check if input file exists
        if !std::path::Path::new(file_path).exists() {
            return Err(format!("Input file not found: {}", file_path));
        }

        // Read sequences from file
        let reader = from_path(file_path).map_err(|e| format!("Error reading file: {}", e))?;

        let records: Vec<SequenceRecord> = reader
            .collect::<Result<Vec<SequenceRecord>, _>>()
            .map_err(|e| format!("Error parsing file: {}", e))?;

        if records.is_empty() {
            return Err("No sequences found in input file".to_string());
        }

        // Process all sequences and return results with sequence names
        let results: Vec<(String, Vec<Result<AnnotationResult, String>>)> = if parallel {
            records
                .into_par_iter()
                .map(|record| {
                    let result = self.number_paired_sequence(&record.sequence, record.name.clone());
                    (record.name, result)
                })
                .collect()
        } else {
            records
                .into_iter()
                .map(|record| {
                    let result = self.number_paired_sequence(&record.sequence, record.name.clone());
                    (record.name, result)
                })
                .collect()
        };

        Ok(results)
    }

    /// Runs alignment of given schemes on sequence and selects one with highest identity
    fn find_highest_identity_chain(
        &self,
        query_sequence: &[u8],
        numbering_schemes: &Vec<&NumberingScheme>,
        sequence_id: String,
    ) -> Result<AnnotationResult, &'static str> {
        let mut highest_identity: f64 = -0.1;
        let mut best_output: Result<AnnotationResult, &'static str> =
            Err("No numbering schemes passed");

        for scheme in numbering_schemes {
            let output: AnnotationResult = scheme.number_sequence(query_sequence, sequence_id.clone());
            if output.identity > highest_identity {
                highest_identity = output.identity;
                best_output = Ok(output);
            }
        }
        best_output
    }

    /// Attempts to find all antibody chains in a sequence
    fn find_all_chains(
        &self,
        query_sequence: &[u8],
        numbering_schemes: Vec<&NumberingScheme>,
        sequence_id: String,
    ) -> Vec<Result<AnnotationResult, &'static str>> {
        let mut chains_found: Vec<Result<AnnotationResult, &str>> = Vec::new();

        let full_query_length: u32 = query_sequence.len() as u32;
        let mut sequence_list: Vec<(&[u8], u32, u32)> = Vec::new();
        sequence_list.push((query_sequence, 0u32, full_query_length - 1));

        while let Some(item) = sequence_list.pop() {
            let (current_sequence, current_start, current_end) = item;

            let numbering_result: Result<AnnotationResult, &str> =
                self.find_highest_identity_chain(current_sequence, &numbering_schemes, sequence_id.clone());

            match numbering_result {
                Ok(mut best_chain) => {
                    if best_chain.identity > MINIMAL_CHAIN_IDENTITY {
                        // split the remaining sequence, add sequences that are long enough to the list
                        let front_sequence: &[u8] = &current_sequence[0..best_chain.start as usize];
                        let end_sequence: &[u8] = &current_sequence[(best_chain.end as usize + 1)..];

                        if front_sequence.len() > MINIMAL_CHAIN_LENGTH as usize {
                            sequence_list.push((
                                front_sequence,
                                current_start,
                                (current_start + best_chain.start - 1),
                            ))
                        }
                        if end_sequence.len() > MINIMAL_CHAIN_LENGTH as usize {
                            sequence_list.push((
                                end_sequence,
                                (current_start + best_chain.end + 1),
                                current_end,
                            ))
                        }

                        // set sequence to full original sequence
                        best_chain.sequence = query_sequence.to_vec();
                        best_chain.start += current_start;
                        best_chain.end += current_start;

                        // add gaps to front and end to match with length of original sequence
                        let mut start_addition = vec![String::from("-"); current_start as usize];
                        let end_addition =
                            vec![String::from("-"); (full_query_length - current_end - 1) as usize];
                        start_addition.extend(best_chain.numbers);
                        start_addition.extend(end_addition);

                        best_chain.numbers = start_addition;
                        chains_found.push(Ok(best_chain));
                    }
                }
                Err(e) => chains_found.push(Err::<AnnotationResult, &str>(e)),
            }
        }

        chains_found
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_annotator_creation() {
        let annotator = Annotator::new(Scheme::IMGT, vec![Chain::IGH], None, None);
        assert!(annotator.is_ok());
    }

    #[test]
    fn test_empty_chains() {
        let annotator = Annotator::new(Scheme::IMGT, vec![], None, None);
        assert!(annotator.is_err());
    }

    #[test]
    fn test_empty_sequence() {
        let annotator = Annotator::new(Scheme::IMGT, vec![Chain::IGH], None, None).unwrap();

        let result = annotator.number_sequence("", "".to_string());
        assert!(result.is_err());
    }

    #[test]
    fn test_prefiltering_enabled() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            None,
            Some(true),
        );
        assert!(annotator.is_ok());
    }

    #[test]
    fn test_prefiltering_disabled() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            None,
            Some(false),
        );
        assert!(annotator.is_ok());
    }

    #[test]
    fn test_prefiltering_with_heavy_chain() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            None,
            Some(true),
        )
        .unwrap();

        // Heavy chain sequence that should be correctly identified with prefiltering
        let heavy_chain_seq = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS";

        let result = annotator.number_sequence(heavy_chain_seq, "".to_string());
        assert!(result.is_ok());
        if let Ok(annotation) = result {
            assert_eq!(annotation.chain, Chain::IGH);
        }
    }

    #[test]
    fn test_number_paired_sequence() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            None,
            Some(false), // Disable prefiltering for this test
        )
        .unwrap();

        // Test sequence with multiple chains (concatenated heavy and light chain)
        let paired_seq = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSSDIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELK";

        let results = annotator.number_paired_sequence(paired_seq, "".to_string());

        // Should find at least one chain in the paired sequence
        assert!(!results.is_empty());

        // Check that we get valid results
        let successful_results: Vec<_> = results.into_iter().filter_map(|r| r.ok()).collect();
        assert!(!successful_results.is_empty());

        // At least one result should have reasonable identity
        let has_good_identity = successful_results.iter().any(|r| r.identity > 0.5);
        assert!(has_good_identity);
    }

    #[test]
    fn test_number_paired_sequence_empty() {
        let annotator = Annotator::new(Scheme::IMGT, vec![Chain::IGH], None, None).unwrap();

        let results = annotator.number_paired_sequence("", "".to_string());
        assert_eq!(results.len(), 1);
        assert!(results[0].is_err());
        assert!(results[0].as_ref().unwrap_err().contains("Empty sequence"));
    }

    #[test]
    fn test_parallel_number_sequences() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            None,
            None,
        )
        .unwrap();

        let sequences = vec![
            "QVQLVQSGAEVKKPGASVKVSCKAS".to_string(),
            "DIQMTQSPSSLSASVGDRVTITC".to_string(),
            "EVQLLESGGGLVQPGGSLRLSCAAS".to_string(),
        ];

        // Test sequential processing
        let sequential_results = annotator.number_sequences(&sequences, false);
        assert_eq!(sequential_results.len(), 3);

        // Test parallel processing
        let parallel_results = annotator.number_sequences(&sequences, true);
        assert_eq!(parallel_results.len(), 3);

        // Results should be the same (order might differ in more complex cases, but for this simple test they should match)
        assert_eq!(sequential_results.len(), parallel_results.len());
    }

    #[test]
    fn test_parallel_number_file() {
        use std::fs::File;
        use std::io::Write;

        // Create a temporary test file
        let test_file_path = "test_parallel.fasta";
        let mut file = File::create(test_file_path).unwrap();
        writeln!(file, ">seq1").unwrap();
        writeln!(file, "QVQLVQSGAEVKKPGASVKVSCKAS").unwrap();
        writeln!(file, ">seq2").unwrap();
        writeln!(file, "DIQMTQSPSSLSASVGDRVTITC").unwrap();
        drop(file);

        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            None,
            None,
        )
        .unwrap();

        // Test sequential file processing
        let sequential_results = annotator.number_file(test_file_path, false).unwrap();
        assert_eq!(sequential_results.len(), 2);

        // Test parallel file processing
        let parallel_results = annotator.number_file(test_file_path, true).unwrap();
        assert_eq!(parallel_results.len(), 2);

        // Clean up
        std::fs::remove_file(test_file_path).ok();
    }

    #[test]
    fn multi_sequence_finall_all() {
        let seq = "PLIIITGDAWWGSNQPPKDRQQQHCRVKMDELANPVASEQSYVLTQPPSVSVSPGQTARITCSGDKLGDKYASWYQQKPGQSPVLVIYQDNKRPSEIPARFSGSNSGNTATLTISGAQAMDEADYYCQAWDSNTGVFGTGTKLTVLGGSSRSSSSGGGGSGGGGVLTQPPSVSAAPGQKVTISCSGSSSNIGNNFVSWYQQRPGTAPSLLIYETNKRPSGIPDRFSGSKSATSATLAITGLQTGDEADYYCATWAASLVFGTGTKVIVSGGSSRSSSSGGGGSGGGGYELTQPPSVSVSPGQTATITCSGDKVASKNVCWYQVKPGQSPEVVMYENYKRPSGIPDRFSGSKSGSTATLTIRGTQATDEADYYCQVWDSFSTFVFGSGTQVTVLGGSSRSSSSGGGGSGGGGYELTQPPSVSVSPGQTASITCSGDKLGNKFTSWYQRKPGQSPVLVIYQDTKRPSGIPERFSGSTSGNTATLTISGTQAMDEADYYCQAWDSSTAWVFGGGTKLEVGGSSRSSSSGGGGSGGGGSYVLTQPPSASGTPGQRVAISCSGSNSNIGSNTVHWYQQLPGAAPKLLIYSNNQRPSGVPDRFSGSNSGTSASLAISRLQSEDEADYYCAAWDDSLNGVVFGGGTKVTVLQIDCEC".as_bytes();

        let annotator = Annotator::new(
            Scheme::KABAT,
            vec![Chain::IGL, Chain::IGK, Chain::IGH],
            None,
            None,
        )
        .unwrap();

        let results = annotator.find_all_chains(seq, annotator.schemes.iter().collect(), "test_multi_sequence".to_string());
        for o in results {
            let o = o.unwrap();
            println!(
                "{} {} {} {:?}\n{}",
                o.start,
                o.end,
                o.numbers.len(),
                o.chain,
                o.sequence_string()
            );
        }
    }

    #[test]
    fn single_sequence_find_all() {
        let seq = "VLTQSPGTLSLSPGETAIISCRTSQYGSLAWYQQRPGQAPRLVIYSGSTRAAGIPDRFSGSRWGPDYNLTISNLESGDFGVYYCQQYEFFGQGTKVQVDIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLRSPVTKSFNRGEC".as_bytes();

        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGL],
            None,
            None,
        )
        .unwrap();

        let annotator_kabat = Annotator::new(
            Scheme::KABAT,
            vec![Chain::IGK],
            None,
            None,
        )
        .unwrap();

        let annotator_heavy = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH],
            None,
            None,
        )
        .unwrap();

        // Collect schemes from different annotators
        let mut schemes: Vec<&NumberingScheme> = Vec::new();
        schemes.extend(annotator.schemes.iter());
        schemes.extend(annotator_kabat.schemes.iter());
        schemes.extend(annotator_heavy.schemes.iter());

        let results = annotator.find_all_chains(seq, schemes, "test_multi_sequence".to_string());
        for o in results {
            println!("{:?}", o);
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

        let annotator_lambda = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGL],
            None,
            None,
        )
        .unwrap();

        let annotator_kappa = Annotator::new(
            Scheme::KABAT,
            vec![Chain::IGK],
            None,
            None,
        )
        .unwrap();

        let annotator_heavy = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH],
            None,
            None,
        )
        .unwrap();

        // Collect schemes from different annotators
        let mut schemes: Vec<&NumberingScheme> = Vec::new();
        schemes.extend(annotator_lambda.schemes.iter());
        schemes.extend(annotator_kappa.schemes.iter());
        schemes.extend(annotator_heavy.schemes.iter());

        assert_eq!(
            annotator_heavy.find_highest_identity_chain(heavy_chain, &schemes, "test_heavy".to_string())
                .expect("")
                .chain,
            Chain::IGH
        );

        assert_eq!(
            annotator_kappa.find_highest_identity_chain(kappa_chain, &schemes, "test_kappa".to_string())
                .expect("")
                .chain,
            Chain::IGK
        );

        assert_eq!(
            annotator_lambda.find_highest_identity_chain(lambda_chain, &schemes, "test_lambda".to_string())
                .expect("")
                .chain,
            Chain::IGL
        );
    }
}
