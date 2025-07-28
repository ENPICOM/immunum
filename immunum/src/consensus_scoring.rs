use crate::constants::{ScoringParams, ENCODED_RESIDUES_MAP};
use crate::constants::{ACCEPTED_RESIDUES, BLOSUM62};
use crate::numbering_scheme_type::NumberingScheme;
use crate::schemes::{
    get_imgt_heavy_scheme, get_imgt_kappa_scheme, get_imgt_lambda_scheme, get_kabat_heavy_scheme,
    get_kabat_kappa_scheme, get_kabat_lambda_scheme,
};
use ndarray::Array2;
use ndarray_npy::{read_npy, write_npy};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

/// Utility to read chain consensus file
pub(crate) fn read_consensus_file(path: PathBuf) -> HashMap<u32, Vec<u8>> {
    let content = fs::read_to_string(path).expect("Error in reading consensus file");
    let mut consensus_aas: HashMap<u32, Vec<u8>> = HashMap::new();
    // Loop over every line of content
    let total_lines = content.lines().count();
    // Skip first and last line
    for line in content.lines().skip(1).take(total_lines - 2) {
        let split_line: Vec<&str> = line.split(',').collect();
        consensus_aas.insert(
            split_line[0].parse::<u32>().unwrap(),
            split_line[1..].join("").into_bytes(),
        );
    }
    consensus_aas
}

/// Finds highest score according to blosum62
fn best_score_consensus(position: u32, residue: u8, consensus: &HashMap<u32, Vec<u8>>) -> i32 {
    // Initialize to minimal score
    let to_check_residues: &Vec<u8> = consensus
        .get(&position)
        .expect("Position outside of consensus");

    // when any is allowed, best score is perfect match
    if to_check_residues.iter().all(|&c| c == b'-') {
        return blosum_lookup(&residue, &residue);
    }
    // iterate over residues to get max score
    blosum_lookup(
        &residue,
        to_check_residues
            .iter()
            .max_by_key(|&c| blosum_lookup(&residue, c))
            .unwrap(),
    )
}

/// Utility to easily find blosum62 score for two residues
fn blosum_lookup(residue1: &u8, residue2: &u8) -> i32 {
    if residue1 < residue2 {
        let lookup: &[u8; 2] = &[*residue1, *residue2];
        *BLOSUM62.get(lookup).unwrap()
    } else {
        let lookup: &[u8; 2] = &[*residue2, *residue1];
        *BLOSUM62.get(lookup).unwrap()
    }
}

/// Calculates and writes scoring matrix for consensus sequence
fn write_scoring_matrix(path: PathBuf, scheme: NumberingScheme, scoring_params: &ScoringParams) {
    let number_accepted_residues = ACCEPTED_RESIDUES.len();
    let consensus_length = scheme.consensus_amino_acids.len();

    let mut matrix = Array2::<f64>::zeros((consensus_length, number_accepted_residues + 2));

    // fill in scores
    for consensus_position in 0..consensus_length {
        for residue_index in 0..number_accepted_residues {
            let residue: u8 = ACCEPTED_RESIDUES[residue_index];
            matrix[[consensus_position, residue_index]] = best_score_consensus(
                (consensus_position + 1) as u32,
                residue,
                &scheme.consensus_amino_acids,
            ) as f64
        }
        // fill in gap penalties
        let (query_gap, consensus_gap) =
            scheme.gap_penalty((consensus_position + 1) as u32, scoring_params);
        // Query gap
        matrix[[consensus_position, number_accepted_residues]] = query_gap;
        // Consensus gap
        matrix[[consensus_position, number_accepted_residues + 1]] = consensus_gap;
    }
    // WRITE MATRIX TO FILE
    write_npy(path, &matrix).expect("Writing scoring matrix failed");
}

/// encodes sequence to indexes for the scoring matrix
pub(crate) fn encode_sequence(input: &[u8]) -> Vec<u8> {
    // Create a lookup table (once, could be static)
    input
        .iter()
        .map(|&residue| *ENCODED_RESIDUES_MAP.get(&residue).unwrap_or(&128))
        .collect()
}

/// Utility to read scoring matrix
pub fn read_scoring_matrix(path: PathBuf) -> Array2<f64> {
    read_npy(path).expect("Error reading scoring matrix")
}

/// Recalculates and writes all scoring matrices (IMGT and KABAT for now)
pub fn write_all_scoring_matrices(scoring_params: &ScoringParams) {
    let schemes: Vec<(NumberingScheme, PathBuf)> = vec![
        (
            get_imgt_heavy_scheme(),
            PathBuf::from("resources")
                .join("consensus")
                .join("IMGT_CONSENSUS_H.npy"),
        ),
        (
            get_imgt_kappa_scheme(),
            PathBuf::from("resources")
                .join("consensus")
                .join("IMGT_CONSENSUS_K.npy"),
        ),
        (
            get_imgt_lambda_scheme(),
            PathBuf::from("resources")
                .join("consensus")
                .join("IMGT_CONSENSUS_L.npy"),
        ),
        (
            get_kabat_heavy_scheme(),
            PathBuf::from("resources")
                .join("consensus")
                .join("KABAT_CONSENSUS_H.npy"),
        ),
        (
            get_kabat_kappa_scheme(),
            PathBuf::from("resources")
                .join("consensus")
                .join("KABAT_CONSENSUS_K.npy"),
        ),
        (
            get_kabat_lambda_scheme(),
            PathBuf::from("resources")
                .join("consensus")
                .join("KABAT_CONSENSUS_L.npy"),
        ),
    ];
    for (scheme, file_path) in schemes {
        write_scoring_matrix(file_path, scheme, scoring_params);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use crate::constants::get_scoring_params;
    //TODO tests for this part

    #[test]
    fn test_lookup() {
        // working
        assert_eq!(blosum_lookup(&b'B', &b'N'), 3);
        assert_eq!(blosum_lookup(&b'W', &b'Y'), 2);
        assert_eq!(blosum_lookup(&b'X', &b'E'), -1);
        // same either way of inputting arguments
        assert_eq!(blosum_lookup(&b'C', &b'S'), blosum_lookup(&b'S', &b'C'));
    }

    #[test]
    fn test_write_scoring_matrix() {
        write_all_scoring_matrices(&get_scoring_params());
    }

    #[test]
    fn test_read_scoring_matrix() {
        let scoring_matrix = read_scoring_matrix(
            PathBuf::from("resources")
                .join("consensus")
                .join("IMGT_CONSENSUS_H.npy"),
        );
        println! {"{scoring_matrix:?}"};
    }
    #[test]
    fn amino_acid_encoding() {
        assert_eq!(
            encode_sequence("ABCDEFGHIKLMNPQRSTVWXYZ".as_bytes()),
            (0..23).collect::<Vec<u8>>()
        )
    }
}
