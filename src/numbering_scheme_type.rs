use crate::insertion_numbering::name_insertions;
use crate::needleman_wunsch::needleman_wunsch_consensus;
use crate::result::AnnotationResult;
use crate::scoring_matrix::ScoringMatrix;
use crate::types::{Chain, RegionRange, Scheme};
use std::collections::HashSet;

#[derive(Debug, Clone)]
pub struct NumberingScheme {
    pub scheme_type: Scheme,
    pub chain_type: Chain,
    pub conserved_positions: HashSet<u32>,
    pub insertion_positions: Vec<u32>,
    pub gap_positions: Vec<u32>,
    pub consensus_amino_acids: Vec<Vec<u8>>, // Indexed by position
    pub restricted_sites: HashSet<u32>,      // Pre-computed for performance
    pub scoring_matrix: ScoringMatrix,
    pub fr1: RegionRange,
    pub fr2: RegionRange,
    pub fr3: RegionRange,
    pub fr4: RegionRange,
    pub cdr1: RegionRange,
    pub cdr2: RegionRange,
    pub cdr3: RegionRange,
}

impl NumberingScheme {
    /// Return restricted sites according to Antpack definition (non '-' positions)
    pub fn restricted_sites(&self) -> &HashSet<u32> {
        &self.restricted_sites
    }
    /// All framework positions
    pub fn framework_positions(&self) -> Vec<u32> {
        self.fr1
            .positions()
            .chain(self.fr2.positions())
            .chain(self.fr3.positions())
            .chain(self.fr4.positions())
            .collect()
    }
    /// All cdr positions
    pub fn cdr_positions(&self) -> Vec<u32> {
        self.cdr1
            .positions()
            .chain(self.cdr2.positions())
            .chain(self.cdr3.positions())
            .collect()
    }
    /// numbers sequence, returns AnnotationResult
    pub(crate) fn number_sequence(&self, query_sequence: &[u8]) -> AnnotationResult {
        // Use the optimized implementation with temporary pool
        let mut matrix_pool = crate::needleman_wunsch::MatrixPool::new();

        let (mut numbering, identity) = needleman_wunsch_consensus(
            query_sequence,
            self,
            &mut matrix_pool,
            -50.0, // Default early termination threshold
        );

        // give gap positions correct names as defined by the numbering scheme
        name_insertions(&mut numbering, &self.scheme_type);

        // find start and end index
        // TODO handle numberings with no positions (only '-')
        let start = numbering
            .iter()
            .position(|s| s != "-")
            .expect("No positions numbered");

        // Find the last non-dash element (searching from the end)
        let end = numbering
            .iter()
            .rposition(|s| s != "-")
            .expect("No positions numbered");

        // Create result with region definitions
        AnnotationResult {
            sequence: query_sequence.to_vec(),
            numbers: numbering,
            scheme: self.scheme_type,
            chain: self.chain_type,
            identity,
            start: start as u32,
            end: end as u32,
            cdr1: self.cdr1.clone(),
            cdr2: self.cdr2.clone(),
            cdr3: self.cdr3.clone(),
            fr1: self.fr1.clone(),
            fr2: self.fr2.clone(),
            fr3: self.fr3.clone(),
            fr4: self.fr4.clone(),
        }
    }
    pub fn to_terminal_schemes(&self, terminal_length: u8) -> (NumberingScheme, NumberingScheme) {
        // Create N-terminal scheme (positions 1 to terminal_length)
        let mut n_terminal_consensus = vec![Vec::new(); terminal_length as usize + 1];
        for i in 1..=terminal_length as u32 {
            if (i as usize) < self.consensus_amino_acids.len()
                && !self.consensus_amino_acids[i as usize].is_empty()
            {
                n_terminal_consensus[i as usize] = self.consensus_amino_acids[i as usize].clone();
            }
        }

        // Calculate restricted sites for N-terminal (non '-' positions)
        let mut n_restricted_sites = HashSet::new();
        for i in 1..=terminal_length as u32 {
            let index = i as usize;
            if index < n_terminal_consensus.len()
                && !n_terminal_consensus[index].is_empty()
                && !n_terminal_consensus[index].contains(&b'-')
            {
                n_restricted_sites.insert(i);
            }
        }

        let n_terminal_scheme = NumberingScheme {
            conserved_positions: HashSet::new(), // No conserved positions at N-terminal
            insertion_positions: vec![],
            gap_positions: vec![10], // IMGT Position 10
            consensus_amino_acids: n_terminal_consensus,
            restricted_sites: n_restricted_sites,
            scoring_matrix: self
                .scoring_matrix
                .slice(0..terminal_length as usize, 0..self.scoring_matrix.ncols()),
            ..self.clone() // Use the remaining fields from the original
        };

        // Create C-terminal scheme
        let c_term_start: u32 = self.fr4.start;
        let mut c_terminal_consensus = vec![Vec::new(); terminal_length as usize + 1];
        for i in 0..terminal_length as u32 {
            let original_position = c_term_start + i;
            if (original_position as usize) < self.consensus_amino_acids.len()
                && !self.consensus_amino_acids[original_position as usize].is_empty()
            {
                let target_index = (i + 1) as usize;
                if target_index < c_terminal_consensus.len() {
                    c_terminal_consensus[target_index] =
                        self.consensus_amino_acids[original_position as usize].clone();
                }
            }
        }

        // Calculate restricted sites for C-terminal (non '-' positions)
        let mut c_restricted_sites = HashSet::new();
        for i in 1..=terminal_length as u32 {
            let index = i as usize;
            if index < c_terminal_consensus.len()
                && !c_terminal_consensus[index].is_empty()
                && !c_terminal_consensus[index].contains(&b'-')
            {
                c_restricted_sites.insert(i);
            }
        }

        let fmwk4_start = self.fr4.start - 1;

        let c_terminal_scheme = NumberingScheme {
            conserved_positions: [1, 2, 4].iter().cloned().collect(), // IMGT position 118, 119 and 121 (mapped to 1, 2, 4)
            insertion_positions: vec![],
            gap_positions: vec![], // No gap positions at c-terminal
            consensus_amino_acids: c_terminal_consensus,
            restricted_sites: c_restricted_sites,
            scoring_matrix: self.scoring_matrix.slice(
                fmwk4_start as usize..(fmwk4_start + terminal_length as u32) as usize,
                0..self.scoring_matrix.ncols(),
            ),
            ..self.clone() // Use the remaining fields from the original
        };

        (n_terminal_scheme, c_terminal_scheme)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        schemes::get_scheme,
        scoring_matrix::ScoringMatrix,
        types::{Chain, Scheme},
    };

    impl ScoringMatrix {
        /// Get a slice of a specific row
        pub fn row(&self, row: usize) -> Vec<f64> {
            let mut row_data = Vec::with_capacity(self.ncols());
            for col in 0..self.ncols() {
                row_data.push(self.get(row, col));
            }
            row_data
        }
    }

    #[test]
    fn terminal_schemes_heavy() {
        let terminal_length = 10;
        let original_scheme = get_scheme(Scheme::IMGT, Chain::IGH, None);
        let (n_term_scheme, c_term_scheme) = original_scheme.to_terminal_schemes(terminal_length);

        for i in 0..terminal_length {
            assert_eq!(
                n_term_scheme.scoring_matrix.row(i as usize),
                original_scheme.scoring_matrix.row(i as usize)
            );

            assert_eq!(
                c_term_scheme.scoring_matrix.row(i as usize),
                original_scheme.scoring_matrix.row(117 + i as usize)
            );
        }

        for i in 1..terminal_length {
            let pos = i as usize;
            if pos < n_term_scheme.consensus_amino_acids.len()
                && pos < original_scheme.consensus_amino_acids.len()
            {
                assert_eq!(
                    n_term_scheme.consensus_amino_acids[pos],
                    original_scheme.consensus_amino_acids[pos]
                );
            }

            let c_pos = (i as u32 + 117) as usize;
            if pos < c_term_scheme.consensus_amino_acids.len()
                && c_pos < original_scheme.consensus_amino_acids.len()
            {
                assert_eq!(
                    c_term_scheme.consensus_amino_acids[pos],
                    original_scheme.consensus_amino_acids[c_pos]
                );
            }
        }
    }

    #[test]
    fn terminal_schemes_kappa() {
        let terminal_length = 10;
        let original_scheme = get_scheme(Scheme::IMGT, Chain::IGK, None);
        let (n_term_scheme, c_term_scheme) = original_scheme.to_terminal_schemes(terminal_length);

        for i in 0..terminal_length {
            assert_eq!(
                n_term_scheme.scoring_matrix.row(i as usize),
                original_scheme.scoring_matrix.row(i as usize)
            );

            assert_eq!(
                c_term_scheme.scoring_matrix.row(i as usize),
                original_scheme.scoring_matrix.row(117 + i as usize)
            );
        }

        for i in 1..terminal_length {
            let pos = i as usize;
            if pos < n_term_scheme.consensus_amino_acids.len()
                && pos < original_scheme.consensus_amino_acids.len()
            {
                assert_eq!(
                    n_term_scheme.consensus_amino_acids[pos],
                    original_scheme.consensus_amino_acids[pos]
                );
            }

            let c_pos = (i as u32 + 117) as usize;
            if pos < c_term_scheme.consensus_amino_acids.len()
                && c_pos < original_scheme.consensus_amino_acids.len()
            {
                assert_eq!(
                    c_term_scheme.consensus_amino_acids[pos],
                    original_scheme.consensus_amino_acids[c_pos]
                );
            }
        }
    }

    #[test]
    fn terminal_schemes_lambda() {
        let terminal_length = 10;
        let original_scheme = get_scheme(Scheme::IMGT, Chain::IGL, None);
        let (n_term_scheme, c_term_scheme) = original_scheme.to_terminal_schemes(terminal_length);

        for i in 0..terminal_length {
            assert_eq!(
                n_term_scheme.scoring_matrix.row(i as usize),
                original_scheme.scoring_matrix.row(i as usize)
            );

            assert_eq!(
                c_term_scheme.scoring_matrix.row(i as usize),
                original_scheme.scoring_matrix.row(117 + i as usize)
            );
        }

        for i in 1..terminal_length {
            let pos = i as usize;
            if pos < n_term_scheme.consensus_amino_acids.len()
                && pos < original_scheme.consensus_amino_acids.len()
            {
                assert_eq!(
                    n_term_scheme.consensus_amino_acids[pos],
                    original_scheme.consensus_amino_acids[pos]
                );
            }

            let c_pos = (i as u32 + 117) as usize;
            if pos < c_term_scheme.consensus_amino_acids.len()
                && c_pos < original_scheme.consensus_amino_acids.len()
            {
                assert_eq!(
                    c_term_scheme.consensus_amino_acids[pos],
                    original_scheme.consensus_amino_acids[c_pos]
                );
            }
        }
    }
}
