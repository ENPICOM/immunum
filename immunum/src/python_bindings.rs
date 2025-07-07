#[cfg(feature = "python")]
use pyo3::prelude::*;

use crate::numbering;
use crate::types::Chain;
use crate::types::Scheme;

/// Perform numbering on a given sequence using the specified scheme and chains
/// # Arguments
/// * `sequence` - The input sequence to be numbered
/// * `scheme` - The numbering scheme to use (e.g., "IMGT" or "KABAT")
/// * `chains` - A vector of chain types to consider (e.g., "IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
/// # Returns
/// A string containing the numbered sequence
#[cfg(feature = "python")]
#[pyfunction]
pub fn number_sequence(sequence: &str, scheme: &str, chains: Vec<String>) -> PyResult<String> {
    let valid_scheme = Scheme::parse_from_string(scheme)
        .map_err(PyErr::new::<pyo3::exceptions::PyValueError, _>)?;

    let valid_chains = Chain::parse_vec_from_strings(chains)
        .map_err(PyErr::new::<pyo3::exceptions::PyValueError, _>)?;

    Ok(numbering::number_sequence(
        sequence,
        &valid_scheme,
        &valid_chains,
    ))
}

/// Immunum Python module
#[cfg(feature = "python")]
#[pymodule]
pub fn immunum(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(number_sequence, m)?)?;
    Ok(())
}
