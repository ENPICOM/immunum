use crate::annotation::find_highest_identity_chain;
use crate::schemes::get_scheme;
use crate::types::{Chain, Scheme};

/// Process a single sequence and return the numbering result
pub fn number_sequence(sequence: &str, scheme: &Scheme, chains: &[Chain]) -> String {
    // Get the appropriate schemes based on the selected scheme type and chains
    let schemes = get_schemes_for_numbering(scheme, chains);
    let scheme_refs: Vec<&_> = schemes.iter().collect();
    
    // Convert sequence to bytes for processing
    let sequence_bytes = sequence.as_bytes();
    
    // Find the best matching scheme
    match find_highest_identity_chain(sequence_bytes, &scheme_refs) {
        Ok(output) => {
            format!("Identity: {:.2}%, Numbering: {:?}", output.identity * 100.0, output.numbering)
        }
        Err(err) => {
            format!("Error: {}", err)
        }
    }
}

/// Helper function to get the appropriate schemes based on scheme type and chains
fn get_schemes_for_numbering(scheme: &Scheme, chains: &[Chain]) -> Vec<crate::numbering_scheme_type::NumberingScheme> {
    let mut schemes = Vec::new();
    
    for chain in chains {
        let numbering_scheme = get_scheme(scheme.clone(), *chain, None);
        schemes.push(numbering_scheme);
    }
    
    schemes
}
