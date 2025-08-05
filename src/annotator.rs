use crate::annotation::{find_all_chains, find_highest_identity_chain};
use crate::constants::{get_scoring_params, ScoringParams};
use crate::fastx::{from_path, FastxRecord};
use crate::numbering_scheme_type::NumberingScheme;
use crate::prefiltering::apply_prefiltering;
use crate::result::AnnotationResult;
use crate::schemes::get_scheme;
use crate::types::{Chain, Scheme};

/// Type alias for complex file processing results
type FileProcessingResult = Result<Vec<(String, Result<AnnotationResult, String>)>, String>;

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
    pub fn number_paired_sequence(&self, sequence: &str) -> Vec<Result<AnnotationResult, String>> {
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
        find_all_chains(sequence_bytes, scheme_refs)
            .into_iter()
            .map(|result| result.map_err(|e| e.to_string()))
            .collect()
    }

    /// Number a single sequence using the pre-configured schemes
    pub fn number_sequence(&self, sequence: &str) -> Result<AnnotationResult, String> {
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

        match find_highest_identity_chain(sequence_bytes, &scheme_refs) {
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
            // TODO: Implement parallel processing using rayon
            // For now, process sequentially
            sequences
                .iter()
                .map(|seq| self.number_sequence(seq))
                .collect()
        } else {
            sequences
                .iter()
                .map(|seq| self.number_sequence(seq))
                .collect()
        }
    }

    /// Process sequences from a FASTA/FASTQ file and return results with sequence names
    pub fn number_file(
        &self,
        file_path: &str,
    ) -> FileProcessingResult {
        // Check if input file exists
        if !std::path::Path::new(file_path).exists() {
            return Err(format!("Input file not found: {}", file_path));
        }

        // Read sequences from file
        let reader = from_path(file_path).map_err(|e| format!("Error reading file: {}", e))?;

        let records: Vec<FastxRecord> = reader
            .collect::<Result<Vec<FastxRecord>, _>>()
            .map_err(|e| format!("Error parsing file: {}", e))?;

        if records.is_empty() {
            return Err("No sequences found in input file".to_string());
        }

        // Process all sequences and return results with sequence names
        let results: Vec<(String, Result<AnnotationResult, String>)> = records
            .into_iter()
            .map(|record| {
                let result = self.number_sequence(&record.sequence);
                (record._name, result)
            })
            .collect();

        Ok(results)
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

        let result = annotator.number_sequence("");
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

        let result = annotator.number_sequence(heavy_chain_seq);
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

        let results = annotator.number_paired_sequence(paired_seq);

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

        let results = annotator.number_paired_sequence("");
        assert_eq!(results.len(), 1);
        assert!(results[0].is_err());
        assert!(results[0].as_ref().unwrap_err().contains("Empty sequence"));
    }
}
