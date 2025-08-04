#[cfg(feature = "python")]
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::annotator::Annotator;
use crate::constants::{ScoringParams, get_scoring_params};
use crate::result::AnnotationResult;
use crate::types::{Chain, Scheme};

/// Python enum for numbering schemes
#[cfg(feature = "python")]
#[pyclass(name = "Scheme")]
#[derive(Clone)]
pub enum PyScheme {
    /// IMGT numbering scheme
    IMGT,
    /// Kabat numbering scheme
    KABAT,
}

impl From<PyScheme> for Scheme {
    fn from(py_scheme: PyScheme) -> Self {
        match py_scheme {
            PyScheme::IMGT => Scheme::IMGT,
            PyScheme::KABAT => Scheme::KABAT,
        }
    }
}

/// Python enum for chain types
#[cfg(feature = "python")]
#[pyclass(name = "Chain")]
#[derive(Clone)]
pub enum PyChain {
    /// Immunoglobulin Heavy chain
    IGH,
    /// Immunoglobulin Kappa chain
    IGK,
    /// Immunoglobulin Lambda chain
    IGL,
    /// T-cell receptor Alpha chain
    TRA,
    /// T-cell receptor Beta chain
    TRB,
    /// T-cell receptor Gamma chain
    TRG,
    /// T-cell receptor Delta chain
    TRD,
}

impl From<PyChain> for Chain {
    fn from(py_chain: PyChain) -> Self {
        match py_chain {
            PyChain::IGH => Chain::IGH,
            PyChain::IGK => Chain::IGK,
            PyChain::IGL => Chain::IGL,
            PyChain::TRA => Chain::TRA,
            PyChain::TRB => Chain::TRB,
            PyChain::TRG => Chain::TRG,
            PyChain::TRD => Chain::TRD,
        }
    }
}

/// Python wrapper for ScoringParams
#[cfg(feature = "python")]
#[pyclass(name = "ScoringParams")]
#[derive(Clone)]
pub struct PyScoringParams {
    inner: ScoringParams,
}

#[cfg(feature = "python")]
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
#[cfg(feature = "python")]
#[pyclass(name = "AnnotationResult")]
pub struct PyAnnotationResult {
    inner: AnnotationResult,
}

#[cfg(feature = "python")]
#[pymethods]
impl PyAnnotationResult {
    #[getter]
    pub fn sequence(&self) -> String {
        self.inner.sequence.clone()
    }

    #[getter]
    pub fn numbers(&self) -> Vec<String> {
        self.inner.numbers.clone()
    }

    #[getter]
    pub fn scheme(&self) -> PyScheme {
        match self.inner.scheme {
            Scheme::IMGT => PyScheme::IMGT,
            Scheme::KABAT => PyScheme::KABAT,
        }
    }

    #[getter]
    pub fn chain(&self) -> PyChain {
        match self.inner.chain {
            Chain::IGH => PyChain::IGH,
            Chain::IGK => PyChain::IGK,
            Chain::IGL => PyChain::IGL,
            Chain::TRA => PyChain::TRA,
            Chain::TRB => PyChain::TRB,
            Chain::TRG => PyChain::TRG,
            Chain::TRD => PyChain::TRD,
        }
    }

    #[getter]
    pub fn identity(&self) -> f64 {
        self.inner.identity
    }

    #[getter]
    pub fn regions(&self, py: Python) -> PyResult<PyObject> {
        let dict = PyDict::new(py);
        for (key, (start, end)) in &self.inner.regions {
            dict.set_item(key, (start, end))?;
        }
        Ok(dict.into())
    }

    pub fn summary(&self) -> String {
        self.inner.summary()
    }
}

/// Python wrapper for Annotator - the main entry point
#[cfg(feature = "python")]
#[pyclass(name = "Annotator")]
pub struct PyAnnotator {
    inner: Annotator,
}

#[cfg(feature = "python")]
#[pymethods]
impl PyAnnotator {
    #[new]
    #[pyo3(signature = (scheme, chains, scoring_params=None))]
    pub fn new(
        scheme: PyScheme,
        chains: PyObject,
        scoring_params: Option<PyScoringParams>,
        py: Python,
    ) -> PyResult<Self> {
        let rust_scheme: Scheme = scheme.into();
        
        // Handle both single chain and list of chains
        let rust_chains: Vec<Chain> = if let Ok(chain_list) = chains.extract::<Vec<PyChain>>(py) {
            chain_list.into_iter().map(|c| c.into()).collect()
        } else if let Ok(single_chain) = chains.extract::<PyChain>(py) {
            vec![single_chain.into()]
        } else {
            return Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "chains must be a Chain or list of Chains"
            ));
        };

        let scoring_params_inner = scoring_params.map(|p| p.inner);

        match Annotator::new(rust_scheme, rust_chains, scoring_params_inner, None) {
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
                    py_results.push(PyAnnotationResult { inner: annotation_result });
                }
                Err(e) => {
                    return Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e));
                }
            }
        }
        
        Ok(py_results)
    }

    #[pyo3(signature = (file_path, output_path=None))]
    pub fn number_file(
        &self,
        file_path: &str,
        output_path: Option<&str>,
    ) -> PyResult<String> {
        match self.inner.number_file(file_path, output_path) {
            Ok(output_file) => Ok(output_file),
            Err(e) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(e)),
        }
    }
}

/// Get default scoring parameters
#[cfg(feature = "python")]
#[pyfunction]
pub fn default_scoring_params() -> PyScoringParams {
    PyScoringParams {
        inner: get_scoring_params(),
    }
}

/// Immunum Python module configuration
#[cfg(feature = "python")]
pub fn immunum(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Main API classes
    m.add_class::<PyAnnotator>()?; // Exported as "Annotator"
    m.add_class::<PyAnnotationResult>()?; // Exported as "AnnotationResult"
    m.add_class::<PyScoringParams>()?; // Exported as "ScoringParams"
    m.add_class::<PyScheme>()?; // Exported as "Scheme"
    m.add_class::<PyChain>()?; // Exported as "Chain"
    
    // Utility functions
    m.add_function(wrap_pyfunction!(default_scoring_params, m)?)?;
    
    Ok(())
}
