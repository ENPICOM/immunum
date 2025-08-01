#[cfg(feature = "python")]
use pyo3::prelude::*;

use crate::annotation::{find_highest_identity_chain, number_sequences_and_write_output};
use crate::constants::{get_scoring_params, MINIMAL_CHAIN_IDENTITY};
use crate::consensus_scoring::write_all_scoring_matrices;
use crate::prefiltering::{get_terminal_schemes, run_pre_scan, select_chains_from_pre_scan};
use crate::schemes::{
    get_imgt_heavy_scheme, get_imgt_kappa_scheme, get_imgt_lambda_scheme, get_kabat_heavy_scheme,
    get_kabat_kappa_scheme, get_kabat_lambda_scheme,
};
use crate::types::{Chain, Scheme};

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

    // Handle empty sequence case
    if sequence.is_empty() {
        return Ok("Empty sequence provided".to_string());
    }

    // Convert sequence to bytes
    let sequence_bytes = sequence.as_bytes();

    // Update scoring matrices
    let scoring_params = get_scoring_params();
    write_all_scoring_matrices(&scoring_params);

    // Get schemes based on selected numbering method
    let schemes = match valid_scheme {
        Scheme::IMGT => vec![
            get_imgt_heavy_scheme(),
            get_imgt_kappa_scheme(),
            get_imgt_lambda_scheme(),
        ],
        Scheme::KABAT => vec![
            get_kabat_heavy_scheme(),
            get_kabat_kappa_scheme(),
            get_kabat_lambda_scheme(),
        ],
    };

    // Only include schemes for selected chains
    let schemes: Vec<_> = schemes
        .into_iter()
        .filter(|scheme| valid_chains.contains(&scheme.chain_type))
        .collect();

    if schemes.is_empty() {
        return Ok("No valid schemes found for selected chains".to_string());
    }

    // Get pre-filter schemes
    let terminal_schemes = get_terminal_schemes(&schemes);

    // Select models to run using pre-scan
    let (pre_scan_output, highest_score) = run_pre_scan(sequence_bytes, &terminal_schemes);
    let pre_filter_chains = select_chains_from_pre_scan(&pre_scan_output, highest_score);
    let filtered_schemes: Vec<_> = schemes
        .iter()
        .filter(|scheme| pre_filter_chains.contains(&scheme.chain_type))
        .collect();

    if filtered_schemes.is_empty() {
        return Ok("No chains passed pre-filtering".to_string());
    }

    // Find highest identity chain
    match find_highest_identity_chain(sequence_bytes, &filtered_schemes) {
        Ok(output) => {
            if output.identity > MINIMAL_CHAIN_IDENTITY {
                Ok(output.get_output_string())
            } else {
                Ok(format!("Low identity score: {:.3}", output.identity))
            }
        }
        Err(e) => Ok(format!("Failed numbering: {}", e)),
    }
}

/// Process sequences from a FASTA or FASTQ file using the specified scheme and chains
/// # Arguments
/// * `input_file` - Path to the FASTA or FASTQ file (supports .gz compression)
/// * `scheme` - The numbering scheme to use (e.g., "IMGT" or "KABAT")
/// * `chains` - A vector of chain types to consider (e.g., "IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
/// * `output_file` - Optional output file path (defaults to input_file + ".numbered.txt")
/// # Returns
/// Path to the output file containing numbered sequences
#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(signature = (input_file, scheme, chains, output_file = None))]
pub fn number_file(
    input_file: &str,
    scheme: &str,
    chains: Vec<String>,
    output_file: Option<String>,
) -> PyResult<String> {
    let valid_scheme = Scheme::parse_from_string(scheme)
        .map_err(PyErr::new::<pyo3::exceptions::PyValueError, _>)?;

    let valid_chains = Chain::parse_vec_from_strings(chains)
        .map_err(PyErr::new::<pyo3::exceptions::PyValueError, _>)?;

    // Generate output file path if not provided
    let output_path = match output_file {
        Some(path) => path,
        None => format!("{}.numbered.txt", input_file),
    };

    // Check if input file exists
    if !std::path::Path::new(input_file).exists() {
        return Err(PyErr::new::<pyo3::exceptions::PyFileNotFoundError, _>(
            format!("Input file not found: {}", input_file),
        ));
    }

    // Use existing file processing functionality
    let scoring_params = get_scoring_params();
    
    // Call the existing function, but handle any panics
    let result = std::panic::catch_unwind(|| {
        number_sequences_and_write_output(
            input_file,
            valid_scheme,
            &valid_chains,
            &output_path,
            true,
            &scoring_params,
        );
    });

    match result {
        Ok(_) => {
            // Check if output file was created successfully
            if std::path::Path::new(&output_path).exists() {
                Ok(output_path)
            } else {
                Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                    "Failed to create output file",
                ))
            }
        }
        Err(_) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
            "Error processing file - check file format and chains",
        )),
    }
}

/// Immunum Python module
#[cfg(feature = "python")]
#[pymodule]
pub fn immunum(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(number_sequence, m)?)?;
    m.add_function(wrap_pyfunction!(number_file, m)?)?;
    Ok(())
}
