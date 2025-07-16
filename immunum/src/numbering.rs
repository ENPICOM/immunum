// use crate::annotation::find_highest_identity_chain;
// use crate::numbering_scheme_type::NumberingScheme;
// use crate::schemes::{
//     get_imgt_heavy_scheme, get_imgt_kappa_scheme, get_imgt_lambda_scheme, get_kabat_heavy_scheme,
//     get_kabat_kappa_scheme, get_kabat_lambda_scheme,
// };
use crate::types::{Chain, Scheme};

/// Mock function for numbering sequences
/// This is a placeholder that returns mock numbering results
pub fn number_sequence(_sequence: &str, _scheme: &Scheme, _chains: &[Chain]) -> String {
    // let schemes: Vec<NumberingScheme> = match scheme {
    //     Scheme::IMGT => vec![
    //         get_imgt_heavy_scheme(),
    //         get_imgt_kappa_scheme(),
    //         get_imgt_lambda_scheme(),
    //     ],
    //     Scheme::KABAT => vec![
    //         get_kabat_heavy_scheme(),
    //         get_kabat_kappa_scheme(),
    //         get_kabat_lambda_scheme(),
    //     ],
    // };
    //
    // let schemes: Vec<NumberingScheme> = schemes
    //     .into_iter()
    //     .filter(|scheme| chains.contains(&scheme.chain_type))
    //     .collect();
    //
    // let output = find_highest_identity_chain(sequence.as_bytes(), &schemes);
    // format!(
    //     "Numbering output: {0:?}{1:?}{2}{3}",
    //     output.numbering, output.scheme.chain_type, output.start, output.end
    // )
    "To be implemented".to_string()
}
