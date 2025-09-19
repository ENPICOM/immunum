use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::annotator::Annotator;
use crate::constants::{get_scoring_params, ScoringParams};
use crate::result::AnnotationResult;
use crate::types::{Chain, Scheme};

/// Python wrapper for AnnotationResult
#[pyclass(name = "AnnotationResult")]
pub struct PyAnnotationResult {
    inner: AnnotationResult,
}

#[pymethods]
impl PyAnnotationResult {
    #[getter]
    pub fn sequence_id(&self) -> String {
        self.inner.sequence_id.clone()
    }

    #[getter]
    pub fn sequence(&self) -> String {
        self.inner.sequence_string()
    }

    #[getter]
    pub fn numbers(&self) -> Vec<String> {
        self.inner.numbers.clone()
    }

    #[getter]
    pub fn scheme(&self) -> Scheme {
        self.inner.scheme
    }

    #[getter]
    pub fn chain(&self) -> Chain {
        self.inner.chain
    }

    #[getter]
    pub fn identity(&self) -> f64 {
        self.inner.identity
    }

    #[getter]
    pub fn regions(&self, py: Python) -> PyResult<PyObject> {
        let dict = PyDict::new(py);
        let region_names = ["cdr1", "cdr2", "cdr3", "fr1", "fr2", "fr3", "fr4"];
        let region_ranges = [
            &self.inner.cdr1,
            &self.inner.cdr2,
            &self.inner.cdr3,
            &self.inner.fr1,
            &self.inner.fr2,
            &self.inner.fr3,
            &self.inner.fr4,
        ];

        for (name, range) in region_names.iter().zip(region_ranges.iter()) {
            dict.set_item(*name, (range.start, range.end))?;
        }
        Ok(dict.into())
    }

    #[getter]
    pub fn start(&self) -> u32 {
        self.inner.start
    }

    #[getter]
    pub fn end(&self) -> u32 {
        self.inner.end
    }

    pub fn get_region_sequence(&self, region_name: &str) -> Option<String> {
        self.inner.get_region_sequence(region_name)
    }
}

/// Python wrapper for Annotator - the main entry point
#[pyclass(name = "Annotator")]
pub struct PyAnnotator {
    inner: Annotator,
}

#[pymethods]
impl PyAnnotator {
    #[new]
    #[pyo3(signature = (scheme, chains, scoring_params=None, use_prefiltering=None))]
    pub fn new(
        scheme: Scheme,
        chains: PyObject,
        scoring_params: Option<ScoringParams>,
        use_prefiltering: Option<bool>,
        py: Python,
    ) -> PyResult<Self> {
        // Handle both single chain and list of chains
        let rust_chains: Vec<Chain> = if let Ok(chain_list) = chains.extract::<Vec<Chain>>(py) {
            chain_list
        } else if let Ok(single_chain) = chains.extract::<Chain>(py) {
            vec![single_chain]
        } else {
            return Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "chains must be a Chain or list of Chains",
            ));
        };

        match Annotator::new(scheme, rust_chains, scoring_params, use_prefiltering) {
            Ok(annotator) => Ok(PyAnnotator { inner: annotator }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e)),
        }
    }

    #[pyo3(signature = (sequence, sequence_id="input".to_string()))]
    pub fn number_sequence(&self, sequence: &str, sequence_id: String) -> PyResult<PyAnnotationResult> {
        match self.inner.number_sequence(sequence, sequence_id) {
            Ok(result) => Ok(PyAnnotationResult { inner: result }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e)),
        }
    }

    #[pyo3(signature = (sequences, parallel=false))]
    pub fn number_sequences(
        &self,
        sequences: Vec<String>,
        parallel: bool,
    ) -> PyResult<Vec<PyAnnotationResult>> {
        let results = self.inner.number_sequences(&sequences, parallel);
        let mut py_results = Vec::new();

        for result in results {
            match result {
                Ok(annotation_result) => {
                    py_results.push(PyAnnotationResult {
                        inner: annotation_result,
                    });
                }
                Err(e) => {
                    return Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e));
                }
            }
        }

        Ok(py_results)
    }

    #[pyo3(signature = (file_path, parallel=false))]
    pub fn number_file(
        &self,
        file_path: &str,
        parallel: bool,
    ) -> PyResult<Vec<(String, PyAnnotationResult)>> {
        match self.inner.number_file(file_path, parallel) {
            Ok(results) => {
                let mut py_results = Vec::new();
                for (name, multi_result) in results {
                    for result in multi_result {
                        match result {
                            Ok(annotation_result) => {
                                py_results.push((
                                    name.clone(),
                                    PyAnnotationResult {
                                        inner: annotation_result,
                                    },
                                ));
                            }
                            Err(e) => {
                                return Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                                    format!("Error processing sequence '{}': {}", name, e),
                                ));
                            }
                        }
                    }
                }
                Ok(py_results)
            }
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e)),
        }
    }
}

/// Get default scoring parameters
#[pyfunction]
pub fn default_scoring_params() -> ScoringParams {
    get_scoring_params()
}

/// Immunum Python module configuration
#[pymodule]
pub fn immunum(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Main API classes
    m.add_class::<PyAnnotator>()?; // Exported as "Annotator"
    m.add_class::<PyAnnotationResult>()?; // Exported as "AnnotationResult"
    m.add_class::<ScoringParams>()?; // Exported as "ScoringParams"
    m.add_class::<crate::types::Scheme>()?;
    m.add_class::<crate::types::Chain>()?;

    // Utility functions
    m.add_function(wrap_pyfunction!(default_scoring_params, m)?)?;

    Ok(())
}
