use crate::insertion_numbering::name_insertions;
use crate::needleman_wunsch::needleman_wunsch_consensus;
use crate::scoring_matrix::ScoringMatrix;
use crate::types::{CdrDefinitions, Chain, NumberingPosition, Scheme};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct NumberingScheme {
    pub scheme_type: Scheme,
    pub chain_type: Chain,
    pub cdr_definition: CdrDefinitions,
    pub consensus_amino_acids: HashMap<u32, Vec<u8>>,
    pub scoring_matrix: ScoringMatrix,
    pub conserved_positions_set: HashSet<u32>,
    pub restricted_sites_set: HashSet<u32>,
    // K-mer set for prefiltering (computed during initialization)
    pub kmer_set: HashSet<String>,
}

impl NumberingScheme {
    /// Number a sequence by alignment, returns numbers and score
    pub(crate) fn number_sequence(
        &self,
        query_sequence: &[u8],
    ) -> (Vec<NumberingPosition>, f64, usize, usize) {
        let (mut numbering, identity) = needleman_wunsch_consensus(query_sequence, self);

        // give gap positions correct names as defined by the numbering scheme
        name_insertions(&mut numbering, &self.scheme_type);

        // find start and end index
        // If these reach the unwrap_or, identity score should be (nearing) 0 and will fail further on
        let start = numbering.iter().position(|s| s != "-").unwrap_or(0);
        let end = numbering.iter().rposition(|s| s != "-").unwrap_or(0);

        (numbering, identity, start, end)
    }
}
