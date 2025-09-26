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

#[pymethods]
impl Annotator {
    #[new]
    #[pyo3(signature = (scheme=Scheme::IMGT, chains=None, disable_prefiltering=false, threads=None, min_confidence=0.7))]
    pub fn new(
        scheme: Scheme,
        chains: Option<Vec<Chain>>,
        disable_prefiltering: bool,
        threads: Option<usize>,
        min_confidence: f64,
    ) -> PyResult<Self> {
        // Default chains if not provided
        let chains = chains.unwrap_or_else(|| vec![Chain::IGH, Chain::IGK, Chain::IGL]);

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
                                    chain.identity,
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
