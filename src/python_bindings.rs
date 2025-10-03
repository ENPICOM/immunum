use clap::ValueEnum;
use pyo3::prelude::*;
use rayon::prelude::*;

use crate::annotator::Annotator as RustAnnotator;
use crate::sequence::SequenceRecord;
use crate::types::{Chain, Scheme};

/// Python wrapper for Annotator - the main entry point
#[pyclass(name = "Annotator")]
pub struct Annotator {
    inner: RustAnnotator,
    thread_pool: rayon::ThreadPool,
}

//PyResult<Vec<Vec<(Vec<String>, f64, Chain)>>>
type PyChainNumbering = (Vec<String>, f64, Chain);
type PySequenceNumbering = Vec<PyChainNumbering>;

/// Helper function to parse Scheme from Python object (string or enum)
fn parse_scheme(obj: &Bound<'_, PyAny>) -> PyResult<Scheme> {
    // Try to extract as Scheme enum first
    if let Ok(scheme) = obj.extract::<Scheme>() {
        return Ok(scheme);
    }

    // Try to extract as string
    if let Ok(s) = obj.extract::<String>() {
        Scheme::from_str(&s, true).map_err(|_| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Invalid scheme '{}'. Valid values: IMGT (I), KABAT (K)",
                s
            ))
        })
    } else {
        Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
            "scheme must be either Scheme enum or string",
        ))
    }
}

/// Helper function to parse Chain from Python object (string or enum)
fn parse_chain(obj: &Bound<'_, PyAny>) -> PyResult<Chain> {
    // Try to extract as Chain enum first
    if let Ok(chain) = obj.extract::<Chain>() {
        return Ok(chain);
    }

    // Try to extract as string
    if let Ok(s) = obj.extract::<String>() {
        Chain::from_str(&s, true)
            .map_err(|_| PyErr::new::<pyo3::exceptions::PyValueError, _>(
                format!("Invalid chain '{}'. Valid values: IGH/Heavy/H, IGK/Kappa/K, IGL/Lambda/L, TRA/Alpha/A, TRB/Beta/B, TRG/Gamma/G, TRD/Delta/D", s)
            ))
    } else {
        Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
            "chain must be either Chain enum or string",
        ))
    }
}

#[pymethods]
impl Annotator {
    #[new]
    #[pyo3(signature = (scheme=None, chains=None, disable_prefiltering=false, threads=None, min_confidence=0.7, min_kmer_overlap=None))]
    pub fn new(
        scheme: Option<Bound<'_, PyAny>>,
        chains: Option<Bound<'_, PyAny>>,
        disable_prefiltering: bool,
        threads: Option<usize>,
        min_confidence: f64,
        min_kmer_overlap: Option<f64>,
    ) -> PyResult<Self> {
        // Parse scheme (default to IMGT if not provided)
        let scheme = if let Some(scheme_obj) = scheme {
            parse_scheme(&scheme_obj)?
        } else {
            Scheme::IMGT
        };

        // Parse chains (default to IGH, IGK, IGL if not provided)
        let chains = if let Some(chains_obj) = chains {
            // Try to extract as list
            if let Ok(chain_list) = chains_obj.extract::<Vec<Bound<'_, PyAny>>>() {
                chain_list
                    .iter()
                    .map(|item| parse_chain(item))
                    .collect::<PyResult<Vec<Chain>>>()?
            } else {
                return Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                    "chains must be a list of Chain enums or strings",
                ));
            }
        } else {
            vec![Chain::IGH, Chain::IGK, Chain::IGL]
        };

        let threads = threads.unwrap_or_else(num_cpus::get);

        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("Failed to create thread pool");

        // TODO enable this in the python API
        let cdr_definitions = None;

        // Create the Rust annotator
        match RustAnnotator::new(
            scheme,
            chains,
            cdr_definitions,
            disable_prefiltering,
            Some(min_confidence),
            min_kmer_overlap,
        ) {
            Ok(annotator) => Ok(Annotator {
                inner: annotator,
                thread_pool,
            }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e)),
        }
    }

    /// Number sequences
    ///
    /// Args:
    ///     sequences: List of sequences. Each can be either:
    ///         - str: Just the sequence
    ///         - tuple(str, str): (sequence_id, sequence)
    ///     max_chains: Maximum number of chains to find per sequence (default: 2)
    ///
    /// Returns:
    ///     List of results for each input sequence. Each result is a list of tuples:
    ///     [(numbers: List[str], confidence: float, chain: Chain), ...]
    #[pyo3(signature = (sequences, max_chains=2))]
    pub fn number_sequences(
        &self,
        sequences: Vec<PyObject>,
        max_chains: usize,
        py: Python,
    ) -> PyResult<Vec<PySequenceNumbering>> {
        // Parse input sequences
        let parsed_sequences: Result<Vec<SequenceRecord>, PyErr> = sequences
            .iter()
            .enumerate()
            .map(|(i, seq_obj)| {
                if let Ok(sequence_str) = seq_obj.extract::<String>(py) {
                    // Just a string - create default ID
                    Ok(SequenceRecord {
                        name: format!("sequence_{}", i).into_bytes(),
                        sequence: sequence_str.into_bytes(),
                    })
                } else if let Ok((seq_id, sequence_str)) = seq_obj.extract::<(String, String)>(py) {
                    // Tuple of (id, sequence)
                    Ok(SequenceRecord {
                        name: seq_id.into_bytes(),
                        sequence: sequence_str.into_bytes(),
                    })
                } else {
                    Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(format!(
                        "Sequence at index {} must be either str or tuple(str, str)",
                        i
                    )))
                }
            })
            .collect();

        let parsed_sequences = parsed_sequences?;

        let results: Vec<PySequenceNumbering> = self.thread_pool.install(|| {
            parsed_sequences
                .par_iter()
                .map(|sequence| {
                    match self.inner.number_sequence(sequence, Some(max_chains)) {
                        Ok(chains) => chains
                            .into_iter()
                            .map(|chain| {
                                (
                                    chain.numbers.iter().map(|n| n.to_string()).collect(),
                                    chain.confidence,
                                    chain.chain,
                                )
                            })
                            .collect(),
                        Err(_) => {
                            // Return empty vector on error
                            Vec::new()
                        }
                    }
                })
                .collect()
        });

        Ok(results)
    }
}

/// Immunum Python module configuration
#[pymodule]
#[pyo3(name = "immunum")]
pub fn immunum(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Main API classes
    m.add_class::<Annotator>()?; // Exported as "Annotator"
    m.add_class::<crate::types::Scheme>()?; // Exported as "Scheme"
    m.add_class::<crate::types::Chain>()?; // Exported as "Chain"

    Ok(())
}
