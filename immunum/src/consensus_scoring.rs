#[allow(dead_code)]
use crate::constants::{ACCEPTED_RESIDUES, BLOSUM62};
use crate::schemes::{
    get_imgt_heavy_scheme, get_imgt_kappa_scheme, get_imgt_lambda_scheme, get_kabat_heavy_scheme,
    get_kabat_kappa_scheme, get_kabat_lambda_scheme,
};
use crate::types::NumberingScheme;
use ndarray::Array2;
use ndarray_npy::{read_npy, write_npy};
use std::collections::HashMap;
use std::fs;

pub(crate) fn read_consensus_file(path: &str) -> HashMap<u32, Vec<char>> {
    let content = fs::read_to_string(path);
    let content = content.unwrap_or("".to_string());
    let mut consensus_aas: HashMap<u32, Vec<char>> = HashMap::new();
    // Loop over every line of content
    let total_lines = content.lines().count();
    // Skip first and last line
    for line in content.lines().skip(1).take(total_lines - 2) {
        let split_line: Vec<&str> = line.split(',').collect();
        consensus_aas.insert(
            split_line[0].parse::<u32>().unwrap(),
            split_line[1..].iter().flat_map(|s| s.chars()).collect(),
        );
    }
    consensus_aas
}

fn best_score_consensus(position: u32, residue: char, consensus: &HashMap<u32, Vec<char>>) -> i32 {
    // Initialize to minimal score
    let to_check_residues: &Vec<char> = consensus
        .get(&position)
        .expect("Position outside of consensus");

    // when any is allowed, best score is perfect match
    if to_check_residues.iter().all(|&c| c == '-') {
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

fn blosum_lookup(residue1: &char, residue2: &char) -> i32 {
    if residue1 < residue2 {
        let lookup: String = format!("{}{}", residue1, residue2);
        *BLOSUM62.get(&lookup).unwrap()
    } else {
        let lookup: String = format!("{}{}", residue2, residue1);
        *BLOSUM62.get(&lookup).unwrap()
    }
}

fn write_scoring_matrix(path: &str, scheme: NumberingScheme) {
    let number_accepted_residues = ACCEPTED_RESIDUES.len();
    let consensus_length = scheme.consensus_amino_acids.len();

    let mut matrix = Array2::<f64>::zeros((consensus_length, number_accepted_residues + 2));

    // fill in scores
    for consensus_position in 0..consensus_length {
        for residue_index in 0..number_accepted_residues {
            let residue: char = ACCEPTED_RESIDUES[residue_index];
            matrix[[consensus_position, residue_index]] = best_score_consensus(
                (consensus_position + 1) as u32,
                residue,
                &scheme.consensus_amino_acids,
            ) as f64
        }
        // fill in gap penalties
        let (query_gap, consensus_gap) = scheme.gap_penalty((consensus_position + 1) as u32);
        // Query gap
        matrix[[consensus_position, number_accepted_residues]] = query_gap;
        // Consensus gap
        matrix[[consensus_position, number_accepted_residues + 1]] = consensus_gap;
    }
    // WRITE MATRIX TO FILE
    write_npy(path, &matrix).expect("Writing scoring matrix failed");
}

pub(crate) fn encode_sequence(input: &str) -> Vec<u8> {
    // Create a lookup table (once, could be static)
    let mut lookup = [255u8; 128]; // Assuming ASCII characters
                                   //const ACCEPTED_RESIDUES: &str = "ACDEFGHIKLMNPQRSTVWY";

    for (i, c) in ACCEPTED_RESIDUES.iter().enumerate() {
        lookup[*c as usize] = i as u8;
    }

    // Use the lookup table for conversion
    input.chars().map(|c| lookup[c as usize]).collect()
}

pub fn read_scoring_matrix(path: &str) -> Array2<f64> {
    read_npy(path).expect("Error reading scoring matrix")
}

fn write_all_scoring_matrices() {
    let schemes: Vec<(NumberingScheme, &str)> = vec![
        (
            get_imgt_heavy_scheme(),
            r"src\consensus\IMGT_HEAVY_CONSENSUS.npy",
        ),
        (
            get_imgt_kappa_scheme(),
            r"src\consensus\IMGT_KAPPA_CONSENSUS.npy",
        ),
        (
            get_imgt_lambda_scheme(),
            r"src\consensus\IMGT_LAMBDA_CONSENSUS.npy",
        ),
        (
            get_kabat_heavy_scheme(),
            r"src\consensus\KABAT_HEAVY_CONSENSUS.npy",
        ),
        (
            get_kabat_kappa_scheme(),
            r"src\consensus\KABAT_KAPPA_CONSENSUS.npy",
        ),
        (
            get_kabat_lambda_scheme(),
            r"src\consensus\KABAT_LAMBDA_CONSENSUS.npy",
        ),
    ];
    for (scheme, file_path) in schemes {
        write_scoring_matrix(file_path, scheme);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    //TODO tests for this part

    #[test]
    fn test_lookup() {
        // working
        assert_eq!(blosum_lookup(&'B', &'N'), 3);
        assert_eq!(blosum_lookup(&'W', &'Y'), 2);
        assert_eq!(blosum_lookup(&'X', &'E'), -1);
        // same either way of inputting arguments
        assert_eq!(blosum_lookup(&'C', &'S'), blosum_lookup(&'S', &'C'));
    }

    #[test]
    fn finding_best_scores() {
        let consensus =
            read_consensus_file(r"C:\Anti_Num\numbering\consensus\IMGT_CONSENSUS_H.txt");

        // correct scores
        assert_eq!(best_score_consensus(1, 'A', &consensus), -1);
        assert_eq!(best_score_consensus(1, 'B', &consensus), 4);
        assert_eq!(best_score_consensus(111, 'Y', &consensus), 7);
        assert_eq!(best_score_consensus(1, 'C', &consensus), -3);

        assert_eq!(best_score_consensus(128, 'X', &consensus), 0);
        assert_eq!(best_score_consensus(128, 'Y', &consensus), -2);
        assert_eq!(best_score_consensus(128, 'Z', &consensus), 0);
    }
    #[test]
    #[should_panic(expected = "Position outside of consensus")]
    fn look_outside_consensus() {
        let consensus =
            read_consensus_file(r"C:\Anti_Num\numbering\consensus\IMGT_CONSENSUS_H.txt");
        // test wrong positions
        best_score_consensus(200, 'A', &consensus);
    }

    #[test]
    fn test_write_scoring_matrix() {
        write_all_scoring_matrices();
    }

    #[test]
    fn test_read_scoring_matrix() {
        let scoring_matrix = read_scoring_matrix(r"src\consensus\IMGT_HEAVY_CONSENSUS.npy");
        println! {"{:?}", scoring_matrix};
    }
    #[test]
    fn amino_acid_encoding() {
        assert_eq!(
            encode_sequence("ABCDEFGHIKLMNPQRSTVWXYZ"),
            (0..23).collect::<Vec<u8>>()
        )
    }
}
