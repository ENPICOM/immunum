use crate::annotation::{apply_prefiltering, find_highest_identity_chain};
use crate::constants::{ScoringParams, get_scoring_params};
use crate::fastx::{from_path, FastxRecord};
use crate::numbering_scheme_type::NumberingScheme;
use crate::result::AnnotationResult;
use crate::schemes::get_scheme;
use crate::types::{Chain, Scheme};

/// Main annotator struct that consolidates all numbering functionality
pub struct Annotator {
    scheme: Scheme,
    chains: Vec<Chain>,
    scoring_params: ScoringParams,
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
            .map(|&chain| get_scheme(scheme.clone(), chain, Some(params.clone())))
            .collect();
        
        if schemes.is_empty() {
            return Err("No valid schemes could be created for the specified chains".to_string());
        }
        
        Ok(Annotator {
            scheme,
            chains,
            scoring_params: params,
            schemes,
            use_prefiltering: enable_prefiltering,
        })
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
            Ok(output) => {
                AnnotationResult::from_numbering_output(output, sequence.to_string())
            }
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
    
    /// Process sequences from a FASTA/FASTQ file
    pub fn number_file(
        &self,
        file_path: &str,
        output_path: Option<&str>,
    ) -> Result<String, String> {
        // Check if input file exists
        if !std::path::Path::new(file_path).exists() {
            return Err(format!("Input file not found: {}", file_path));
        }
        
        // Generate output path if not provided
        let output_file = match output_path {
            Some(path) => path.to_string(),
            None => format!("{}.numbered.txt", file_path),
        };
        
        // Read sequences from file
        let reader = from_path(file_path)
            .map_err(|e| format!("Error reading file: {}", e))?;
        
        let records: Vec<FastxRecord> = reader
            .collect::<Result<Vec<FastxRecord>, _>>()
            .map_err(|e| format!("Error parsing file: {}", e))?;
        
        if records.is_empty() {
            return Err("No sequences found in input file".to_string());
        }
        
        // Process all sequences
        let mut output_lines = Vec::new();
        
        // Add header
        output_lines.push("Name\tSequence\tNumbering\tScore\tChain\tcdr1\tcdr2\tcdr3\tfmwk1\tfmwk2\tfmwk3\tfmwk4\tStart\tEnd".to_string());
        
        for record in records {
            match self.number_sequence(&record.sequence) {
                Ok(result) => {
                    let line = format!(
                        "{}\t{}\t{}\t{:.3}\t{:?}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        record._name,
                        result.sequence,
                        result.numbers.join(","),
                        result.identity,
                        result.chain,
                        result.regions.get("cdr1").map_or("", |_| ""), // TODO: Extract actual region sequences
                        result.regions.get("cdr2").map_or("", |_| ""),
                        result.regions.get("cdr3").map_or("", |_| ""),
                        result.regions.get("fr1").map_or("", |_| ""),
                        result.regions.get("fr2").map_or("", |_| ""),
                        result.regions.get("fr3").map_or("", |_| ""),
                        result.regions.get("fr4").map_or("", |_| ""),
                        0, // start position - TODO: implement
                        result.sequence.len()
                    );
                    output_lines.push(line);
                }
                Err(e) => {
                    let line = format!(
                        "{}\t{}\tFailed numbering: {}\t0\tUnknown\t\t\t\t\t\t\t\t0\t0",
                        record._name, record.sequence, e
                    );
                    output_lines.push(line);
                }
            }
        }
        
        // Write results to file
        let output_content = output_lines.join("\n") + "\n";
        std::fs::write(&output_file, output_content)
            .map_err(|e| format!("Error writing output file: {}", e))?;
        
        Ok(output_file)
    }
    
    /// Get the scheme used by this annotator
    pub fn scheme(&self) -> &Scheme {
        &self.scheme
    }
    
    /// Get the chains used by this annotator
    pub fn chains(&self) -> &[Chain] {
        &self.chains
    }
    
    /// Get the scoring parameters used by this annotator
    pub fn scoring_params(&self) -> &ScoringParams {
        &self.scoring_params
    }
    
    /// Check if prefiltering is enabled
    pub fn uses_prefiltering(&self) -> bool {
        self.use_prefiltering
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_annotator_creation() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH],
            None,
            None,
        );
        assert!(annotator.is_ok());
        assert!(!annotator.unwrap().uses_prefiltering());
    }
    
    #[test]
    fn test_empty_chains() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![],
            None,
            None,
        );
        assert!(annotator.is_err());
    }
    
    #[test]
    fn test_empty_sequence() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH],
            None,
            None,
        ).unwrap();
        
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
        assert!(annotator.unwrap().uses_prefiltering());
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
        assert!(!annotator.unwrap().uses_prefiltering());
    }
    
    #[test]
    fn test_prefiltering_with_heavy_chain() {
        let annotator = Annotator::new(
            Scheme::IMGT,
            vec![Chain::IGH, Chain::IGK, Chain::IGL],
            None,
            Some(true),
        ).unwrap();
        
        // Heavy chain sequence that should be correctly identified with prefiltering
        let heavy_chain_seq = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARWGGRGSYAMDYWGQGTLVTVSS";
        
        let result = annotator.number_sequence(heavy_chain_seq);
        assert!(result.is_ok());
        if let Ok(annotation) = result {
            assert_eq!(annotation.chain, Chain::IGH);
        }
    }
}