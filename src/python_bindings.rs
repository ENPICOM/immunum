use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::annotator::Annotator;
use crate::constants::{get_scoring_params, ScoringParams};
use crate::result::AnnotationResult;
use crate::types::{Chain, Scheme};

/// Python wrapper for ScoringParams
#[pyclass(name = "ScoringParams")]
#[derive(Clone)]
pub struct PyScoringParams {
    inner: ScoringParams,
}

#[pymethods]
impl PyScoringParams {
    #[new]
    #[pyo3(signature = (gap_pen_cp=None, gap_pen_fr=None, gap_pen_ip=None, gap_pen_op=None, gap_pen_cdr=None, gap_pen_other=None, cdr_increase=None, pen_leap_insertion_point_imgt=None, pen_leap_insertion_point_kabat=None))]
    pub fn new(
        gap_pen_cp: Option<f64>,
        gap_pen_fr: Option<f64>,
        gap_pen_ip: Option<f64>,
        gap_pen_op: Option<f64>,
        gap_pen_cdr: Option<f64>,
        gap_pen_other: Option<f64>,
        cdr_increase: Option<f64>,
        pen_leap_insertion_point_imgt: Option<f64>,
        pen_leap_insertion_point_kabat: Option<f64>,
    ) -> Self {
        let default_params = ScoringParams::default();
        PyScoringParams {
            inner: ScoringParams {
                gap_pen_cp: gap_pen_cp.unwrap_or(default_params.gap_pen_cp),
                gap_pen_fr: gap_pen_fr.unwrap_or(default_params.gap_pen_fr),
                gap_pen_ip: gap_pen_ip.unwrap_or(default_params.gap_pen_ip),
                gap_pen_op: gap_pen_op.unwrap_or(default_params.gap_pen_op),
                gap_pen_cdr: gap_pen_cdr.unwrap_or(default_params.gap_pen_cdr),
                gap_pen_other: gap_pen_other.unwrap_or(default_params.gap_pen_other),
                cdr_increase: cdr_increase.unwrap_or(default_params.cdr_increase),
                pen_leap_insertion_point_imgt: pen_leap_insertion_point_imgt
                    .unwrap_or(default_params.pen_leap_insertion_point_imgt),
                pen_leap_insertion_point_kabat: pen_leap_insertion_point_kabat
                    .unwrap_or(default_params.pen_leap_insertion_point_kabat),
            },
        }
    }

    #[getter]
    pub fn gap_pen_cp(&self) -> f64 {
        self.inner.gap_pen_cp
    }

    #[setter]
    pub fn set_gap_pen_cp(&mut self, value: f64) {
        self.inner.gap_pen_cp = value;
    }

    #[getter]
    pub fn gap_pen_op(&self) -> f64 {
        self.inner.gap_pen_op
    }

    #[setter]
    pub fn set_gap_pen_op(&mut self, value: f64) {
        self.inner.gap_pen_op = value;
    }
}

/// Python wrapper for AnnotationResult
#[pyclass(name = "AnnotationResult")]
pub struct PyAnnotationResult {
    inner: AnnotationResult,
}

#[pymethods]
impl PyAnnotationResult {
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
        scoring_params: Option<PyScoringParams>,
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

        let scoring_params_inner = scoring_params.map(|p| p.inner);

        match Annotator::new(scheme, rust_chains, scoring_params_inner, use_prefiltering) {
            Ok(annotator) => Ok(PyAnnotator { inner: annotator }),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e)),
        }
    }

    pub fn number_sequence(&self, sequence: &str) -> PyResult<PyAnnotationResult> {
        match self.inner.number_sequence(sequence) {
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

    #[pyo3(signature = (file_path))]
    pub fn number_file(&self, file_path: &str) -> PyResult<Vec<(String, PyAnnotationResult)>> {
        match self.inner.number_file(file_path) {
            Ok(results) => {
                let mut py_results = Vec::new();
                for (name, result) in results {
                    match result {
                        Ok(annotation_result) => {
                            py_results.push((
                                name,
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
                Ok(py_results)
            }
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e)),
        }
    }
}

/// Get default scoring parameters
#[pyfunction]
pub fn default_scoring_params() -> PyScoringParams {
    PyScoringParams {
        inner: get_scoring_params(),
    }
}

/// Immunum Python module configuration
#[pymodule]
pub fn immunum(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Main API classes
    m.add_class::<PyAnnotator>()?; // Exported as "Annotator"
    m.add_class::<PyAnnotationResult>()?; // Exported as "AnnotationResult"
    m.add_class::<PyScoringParams>()?; // Exported as "ScoringParams"
    m.add_class::<crate::types::Scheme>()?;
    m.add_class::<crate::types::Chain>()?;

    // Utility functions
    m.add_function(wrap_pyfunction!(default_scoring_params, m)?)?;

    Ok(())
}
