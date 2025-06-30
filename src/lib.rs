use pyo3::prelude::*;

use crate::types::{Chain, Scheme};

mod numbering;
mod types;

/// Perform numbering on a given sequence using the specified scheme and chains
/// # Arguments
/// * `sequence` - The input sequence to be numbered
/// * `scheme` - The numbering scheme to use (e.g., "IMGT" or "KABAT")
/// * `chains` - A vector of chain types to consider (e.g., "IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
/// # Returns
/// A string containing the numbered sequence
#[pyfunction]
fn number_sequence(sequence: &str, scheme: &str, chains: Vec<String>) -> PyResult<String> {
    let valid_scheme = match scheme.to_lowercase().as_str() {
        "imgt" | "i" => Scheme::IMGT,
        "kabat" | "k" => Scheme::KABAT,
        _ => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Invalid scheme. Valid schemes are: 'imgt', 'i', 'kabat', 'k'",
            ))
        }
    };

    let mut valid_chains = Vec::new();
    for chain in chains {
        let parsed_chain = match chain.to_lowercase().as_str() {
            "igh" | "h" => Chain::IGH,
            "igk" | "k" => Chain::IGK,
            "igl" | "l" => Chain::IGL,
            "tra" | "a" => Chain::TRA,
            "trb" | "b" => Chain::TRB,
            "trg" | "g" => Chain::TRG,
            "trd" | "d" => Chain::TRD,
            _ => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    format!("Invalid chain '{}'. Valid chains are: 'igh', 'h', 'igk', 'k', 'igl', 'l', 'tra', 'a', 'trb', 'b', 'trg', 'g', 'trd', 'd'", chain),
                ))
            }
        };
        valid_chains.push(parsed_chain);
    }

    Ok(numbering::number_sequence(
        sequence,
        &valid_scheme,
        &valid_chains,
    ))
}

/// A Python module implemented in Rust.
#[pymodule]
fn antinum(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(number_sequence, m)?)?;
    Ok(())
}
