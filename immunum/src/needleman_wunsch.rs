use crate::constants::{traceback_directions, scoring};
use crate::schemes::*;
use crate::types::{NumberingOutput, NumberingScheme};

pub fn needleman_wunsch_consensus(
    query_sequence: String,
    scheme: &NumberingScheme,
) -> (Vec<String>, u32) {
    let num_positions_consensus: usize = scheme.consensus_amino_acids.len();
    let len_query_sequence: usize = query_sequence.len();

    let mut dynamic_matrix: Vec<Vec<f64>> =
        vec![vec![0.0; len_query_sequence + 1]; num_positions_consensus + 1];
    let mut traceback_matrix: Vec<Vec<u8>> =
        vec![vec![0; len_query_sequence + 1]; num_positions_consensus + 1];

    for i in 0..len_query_sequence + 1 {
        dynamic_matrix[0][i] = 0.0 - (scoring::GAP_PEN_START * i as f64); // Formula: 2 times the index
        traceback_matrix[0][i] = traceback_directions::FROM_LEFT;
    }
    for i in 0..num_positions_consensus + 1 {
        dynamic_matrix[i][0] = 0.0 - (scoring::GAP_PEN_START * i as f64); // Formula: 2 times the index
        traceback_matrix[i][0] = traceback_directions::FROM_TOP;
    }

    for consensus_position in 1..num_positions_consensus {
        for query_position in 1..len_query_sequence {
            let query_gap_penalty: f64 = 20.0; // TODO get real values
            let consensus_gap_penalty: f64 = 20.0; // TODO get real values

            let top_value: f64 =
                dynamic_matrix[consensus_position - 1][query_position] - consensus_gap_penalty;
            let left_value: f64 =
                dynamic_matrix[consensus_position][query_position - 1] - query_gap_penalty;
            let mut match_value: f64 = dynamic_matrix[consensus_position - 1][query_position - 1];

            // NOTE: sequence index is query_position - 1
            let mut best_score: f64 = 2.0; // TODO Get scoring value

            if scheme
                .conserved_positions
                .contains(&(consensus_position as u32))
            {
                best_score *= scoring::MATCH_CP_MULTIPLIER
            }

            match_value += best_score;

            let mut max_score = match_value;
            let mut transfer = traceback_directions::FROM_DIAG;
            if top_value > max_score {
                max_score = top_value;
                transfer = traceback_directions::FROM_TOP;
            }
            if left_value > max_score {
                max_score = left_value;
                transfer = traceback_directions::FROM_LEFT
            }

            dynamic_matrix[consensus_position][query_position] = max_score;
            traceback_matrix[consensus_position][query_position] = transfer;

            let seq_char = query_sequence.chars().nth(query_position - 1).unwrap();
            if transfer == traceback_directions::FROM_DIAG
                && scheme.consensus_amino_acids[&(consensus_position as u32)].contains(&seq_char)
                && scheme
                    .restricted_sites()
                    .contains(&(consensus_position as u32))
            {
                traceback_matrix[consensus_position][query_position] = traceback_directions::PERFECT_MATCH;
            }
        }
    }
    // Traceback
    let (numbering, matches) = traceback_alignment(
        len_query_sequence,
        num_positions_consensus,
        &traceback_matrix,
    );

    let identity = matches / (scheme.restricted_sites().len() as u32);
    (numbering, identity)
}

fn traceback_alignment(
    sequence_length: usize,
    consensus_length: usize,
    traceback_matrix: &Vec<Vec<u8>>,
) -> (Vec<String>, u32) {
    let mut matches: u32 = 0;
    let mut numbering = Vec::new();

    let mut i = consensus_length;
    let mut j = sequence_length;

    while i != 0 || j != 0 {
        // Figure out if it was top, left or top left
        if traceback_matrix[i][j] == traceback_directions::FROM_LEFT {
            numbering.insert(0, "-".to_string());
            j -= 1;
        } else if traceback_matrix[i][j] == traceback_directions::FROM_TOP {
            i -= 1;
        } else {
            // FROM_DIAG or PERFECT_MATCH
            numbering.insert(0, i.to_string());
            if traceback_matrix[i][j] == traceback_directions::PERFECT_MATCH {
                matches += 1; // Used to calculate identity
            }
            i -= 1;
            j -= 1;
        }
    }

    (numbering, matches)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_needleman_wunsch_consensus() {
        let heavy_chain: String =
            "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKF\
        QGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPE\
        PVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSTWPSETVTCNVAHPASSTKVDKKIVPRD"
                .to_string();
        let scheme = get_imgt_heavy_scheme();
        let output = needleman_wunsch_consensus(heavy_chain, &scheme);
        println!("{:?}", output);
    }

    #[test]
    fn test_traceback() {
        let traceback1: Vec<Vec<u8>> = vec![
            vec![1, 1, 1, 1, 1, 1, 1, 1],
            vec![2, 0, 1, 1, 1, 1, 1, 1],
            vec![2, 1, 0, 1, 1, 1, 1, 1],
            vec![2, 1, 0, 0, 1, 1, 1, 1],
            vec![2, 1, 1, 1, 3, 1, 1, 1],
            vec![2, 1, 0, 1, 1, 0, 1, 1],
            vec![2, 1, 0, 1, 1, 1, 3, 1],
            vec![2, 1, 0, 1, 1, 1, 1, 0],
        ];

        let traceback2: Vec<Vec<u8>> = vec![
            vec![1, 1, 1, 1, 1, 1, 1, 1],
            vec![2, 1, 1, 1, 1, 1, 0, 1],
            vec![2, 1, 1, 1, 1, 0, 1, 1],
            vec![2, 1, 0, 1, 1, 1, 0, 1],
            vec![2, 0, 1, 1, 1, 1, 1, 1],
            vec![2, 1, 0, 1, 1, 1, 1, 1],
            vec![2, 1, 0, 1, 3, 1, 1, 1],
            vec![2, 1, 0, 1, 1, 1, 0, 1],
        ];

        let output1 = traceback_alignment(7, 7, &traceback1);
        let output2 = traceback_alignment(7, 7, &traceback2);
        assert_eq!(output1.0, vec!["1", "2", "3", "4", "5", "6", "7"]);
        assert_eq!(output1.1, 2);
        assert_eq!(output2.0, vec!["4", "5", "-", "6", "-", "7", "-"]);
        assert_eq!(output2.1, 1);
    }
}
