use crate::types::{NumberingScheme, NumberingOutput};

pub fn needleman_wunsch_consensus<'a>(query_sequence:String, scheme: &'a NumberingScheme)
                                      -> &'a NumberingOutput {

    let num_positions_consensus: usize = scheme.consensus_amino_acids.len();
    let len_query_sequence: usize = query_sequence.len();

    (0..=num_items).fold(
        vec![vec![0; len_query_sequence + 1]; num_positions_consensus + 1],
        |mut matrix, consensus_position| {
            (0..=len_query_sequence).for_each(|query_position| {
                println!("{}{}", consensus_position, query_position);
            });
            // TODO fill fields
            NumberingOutput{};
        },
    )
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_needleman_wunsch_consensus() {
        needleman_wunsch_consensus();
    }
}
