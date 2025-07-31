use crate::annotation::find_highest_identity_chain;
use crate::schemes::{
    get_imgt_heavy_scheme, get_imgt_kappa_scheme, get_imgt_lambda_scheme,
    get_kabat_heavy_scheme, get_kabat_kappa_scheme, get_kabat_lambda_scheme,
};
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
        let numbering_scheme = match (scheme, chain) {
            (Scheme::IMGT, Chain::IGH) => get_imgt_heavy_scheme(),
            (Scheme::IMGT, Chain::IGK) => get_imgt_kappa_scheme(),
            (Scheme::IMGT, Chain::IGL) => get_imgt_lambda_scheme(),
            (Scheme::KABAT, Chain::IGH) => get_kabat_heavy_scheme(),
            (Scheme::KABAT, Chain::IGK) => get_kabat_kappa_scheme(),
            (Scheme::KABAT, Chain::IGL) => get_kabat_lambda_scheme(),
            // For other chains, default to heavy chain schemes for now
            (Scheme::IMGT, _) => get_imgt_heavy_scheme(),
            (Scheme::KABAT, _) => get_kabat_heavy_scheme(),
        };
        schemes.push(numbering_scheme);
    }
    
    schemes
}
