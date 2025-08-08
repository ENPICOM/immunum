use crate::constants::{get_scoring_params, ScoringParams};
use crate::sequence_io::{from_path, SequenceRecord};
use crate::needleman_wunsch::{needleman_wunsch_consensus, MatrixPool};
use crate::numbering_scheme::NumberingScheme;
use crate::prefiltering::prefilter_schemes;
use crate::result::AnnotationResult;
use crate::schemes::get_scheme;
use crate::types::{Chain, Scheme};
use rayon::prelude::*;

/// Type alias for complex file processing results
type FileProcessingResult = Result<Vec<(String, Vec<Result<AnnotationResult, String>>)>, String>;

/// Main annotator struct that consolidates all numbering functionality
pub struct Annotator {
    _scheme: Scheme,
    _chains: Vec<Chain>,
    _scoring_params: ScoringParams,
    schemes: Vec<NumberingScheme>,
    use_prefiltering: bool,
    early_termination_threshold: f64,
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
            early_termination_threshold: -50.0, // Default early termination threshold
        })
    }

    /// Number a single sequence using optimized algorithm
    pub fn number_sequence(&self, sequence: &str) -> Result<AnnotationResult, String> {
        if sequence.is_empty() {
            return Err("Empty sequence provided".to_string());
        }

        let sequence_bytes = sequence.as_bytes();

        // Apply prefiltering if enabled
        let scheme_refs: Vec<&NumberingScheme> = if self.use_prefiltering {
            prefilter_schemes(sequence_bytes, &self.schemes)
        } else {
            self.schemes.iter().collect()
        };

        self.find_best_match(sequence_bytes, &scheme_refs)
    }

    /// Find highest identity chain using optimized algorithm
    fn find_best_match(
        &self,
        query_sequence: &[u8],
        numbering_schemes: &[&NumberingScheme],
    ) -> Result<AnnotationResult, String> {
        let mut highest_identity: f64 = -0.1;
        let mut best_output: Result<AnnotationResult, String> =
            Err("No numbering schemes passed".to_string());

        // Use thread-local matrix pool for better performance
        thread_local! {
            static MATRIX_POOL: std::cell::RefCell<MatrixPool> =
                std::cell::RefCell::new(MatrixPool::new());
        }

        for &scheme in numbering_schemes {
            let output = MATRIX_POOL.with(|pool_cell| {
                let mut pool = pool_cell.borrow_mut();
                let (mut numbering, identity) = needleman_wunsch_consensus(
                    query_sequence,
                    scheme,
                    &mut pool,
                    self.early_termination_threshold,
                );

                // Apply insertion naming according to the scheme (KABAT vs IMGT)
                crate::insertion_naming::name_insertions(&mut numbering, &scheme.scheme_type);

                // Create AnnotationResult
                let start = numbering.iter().position(|s| s != "-").unwrap_or(0);

                let end = numbering
                    .iter()
                    .rposition(|s| s != "-")
                    .unwrap_or(numbering.len() - 1);

                AnnotationResult {
                    sequence: query_sequence.to_vec(),
                    numbers: numbering,
                    scheme: scheme.scheme_type,
                    chain: scheme.chain_type,
                    identity,
                    start: start as u32,
                    end: end as u32,
                    cdr1: scheme.cdr1.clone(),
                    cdr2: scheme.cdr2.clone(),
                    cdr3: scheme.cdr3.clone(),
                    fr1: scheme.fr1.clone(),
                    fr2: scheme.fr2.clone(),
                    fr3: scheme.fr3.clone(),
                    fr4: scheme.fr4.clone(),
                }
            });

            if output.identity > highest_identity {
                highest_identity = output.identity;
                best_output = Ok(output);
            }
        }

        best_output
    }

    /// Number multiple sequences, optionally in parallel
    pub fn number_sequences(
        &self,
        sequences: &[String],
        parallel: bool,
    ) -> Vec<Result<AnnotationResult, String>> {
        if parallel {
            sequences
                .par_iter()
                .map(|seq| self.number_sequence(seq))
                .collect()
        } else {
            sequences
                .iter()
                .map(|seq| self.number_sequence(seq))
                .collect()
        }
    }

    /// Number a sequence and try to find all possible chains (for paired sequences)
    pub fn number_paired_sequence(&self, sequence: &str) -> Vec<Result<AnnotationResult, String>> {
        if sequence.is_empty() {
            return vec![Err("Empty sequence provided".to_string())];
        }

        let sequence_bytes = sequence.as_bytes();

        // Apply prefiltering if enabled
        let scheme_refs: Vec<&NumberingScheme> = if self.use_prefiltering {
            prefilter_schemes(sequence_bytes, &self.schemes)
        } else {
            self.schemes.iter().collect()
        };

        // Try all schemes and return all successful results
        let mut results = Vec::new();
        thread_local! {
            static MATRIX_POOL: std::cell::RefCell<MatrixPool> =
                std::cell::RefCell::new(MatrixPool::new());
        }

        for &scheme in &scheme_refs {
            let result = MATRIX_POOL.with(|pool_cell| {
                let mut pool = pool_cell.borrow_mut();
                let (mut numbering, identity) = needleman_wunsch_consensus(
                    sequence_bytes,
                    scheme,
                    &mut pool,
                    self.early_termination_threshold,
                );

                // Apply insertion naming according to the scheme (KABAT vs IMGT)
                crate::insertion_naming::name_insertions(&mut numbering, &scheme.scheme_type);

                // Create AnnotationResult
                let start = numbering.iter().position(|s| s != "-").unwrap_or(0);

                let end = numbering
                    .iter()
                    .rposition(|s| s != "-")
                    .unwrap_or(numbering.len() - 1);

                AnnotationResult {
                    sequence: sequence_bytes.to_vec(),
                    numbers: numbering,
                    scheme: scheme.scheme_type,
                    chain: scheme.chain_type,
                    identity,
                    start: start as u32,
                    end: end as u32,
                    cdr1: scheme.cdr1.clone(),
                    cdr2: scheme.cdr2.clone(),
                    cdr3: scheme.cdr3.clone(),
                    fr1: scheme.fr1.clone(),
                    fr2: scheme.fr2.clone(),
                    fr3: scheme.fr3.clone(),
                    fr4: scheme.fr4.clone(),
                }
            });

            // Only include results with reasonable identity scores
            if result.identity > 0.1 {
                results.push(Ok(result));
            }
        }

        // If no good results, return the best single result
        if results.is_empty() {
            match self.number_sequence(sequence) {
                Ok(result) => vec![Ok(result)],
                Err(e) => vec![Err(e)],
            }
        } else {
            results
        }
    }

    /// Number all sequences in a FASTA/FASTQ file
    pub fn number_file(&self, file_path: &str, parallel: bool) -> FileProcessingResult {
        let reader = from_path(file_path)
            .map_err(|e| format!("Failed to open file '{}': {}", file_path, e))?;

        let records: Vec<SequenceRecord> = reader
            .collect::<Result<Vec<SequenceRecord>, _>>()
            .map_err(|e| format!("Failed to read sequences from file: {}", e))?;

        let results: Vec<(String, Vec<Result<AnnotationResult, String>>)> = if parallel {
            records
                .par_iter()
                .map(|record| {
                    let results = vec![self.number_sequence(&record.sequence)];
                    (record._name.clone(), results)
                })
                .collect()
        } else {
            records
                .iter()
                .map(|record| {
                    let results = vec![self.number_sequence(&record.sequence)];
                    (record._name.clone(), results)
                })
                .collect()
        };

        Ok(results)
    }
}
