use std::str::FromStr;

use serde::{Deserialize, Serialize};

use postcard::{from_bytes, to_allocvec};
use pyo3::prelude::*;

use crate::annotator::{Annotator, NumberingResult};
use crate::types::{Chain, Scheme};

#[pymethods]
impl Annotator {
    // python methods
    #[new]
    #[pyo3(signature = (chains , scheme))]
    pub fn init(chains: Vec<String>, scheme: String) -> PyResult<Self> {
        let parsed_chains = chains
            .iter()
            .map(|chain| {
                Chain::from_str(&chain).map_err(|_| {
                    PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                        "Invalid chain: {}",
                        chain
                    ))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let parsed_scheme = Scheme::from_str(&scheme).map_err(|_| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Invalid scheme: {}", scheme))
        })?;
        let annotator =
            Annotator::new(parsed_chains.as_slice(), parsed_scheme.clone()).map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Failed to initialize annotator: {}",
                    e
                ))
            })?;

        Ok(annotator)
    }

    #[pyo3(signature = (sequence), name = "number")]
    pub fn _number(&self, sequence: &str) -> PyResult<NumberingResult> {
        let result = self
            .number(sequence)
            .map_err(|_| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Invalid sequence: {}",
                    sequence
                ))
            })
            .unwrap();
        Ok(result)
    }

    pub fn __setstate__(
        &mut self,
        state: &pyo3::Bound<'_, pyo3::types::PyBytes>,
    ) -> pyo3::PyResult<()> {
        let wrapper: AnnotatorSerializationWrapper = from_bytes(state.as_bytes()).unwrap();
        self.matrices = wrapper.inner.matrices;
        self.scheme = wrapper.inner.scheme;
        Ok(())
    }

    pub fn __getstate__<'py>(
        &self,
        py: pyo3::Python<'py>,
    ) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyBytes>> {
        let wrapper = AnnotatorSerializationWrapper::from_annotator(self);
        Ok(pyo3::types::PyBytes::new(
            py,
            &to_allocvec(&wrapper.inner).unwrap(),
        ))
    }
    pub fn __getnewargs__(&self) -> pyo3::PyResult<(Vec<String>, String)> {
        let wrapper = AnnotatorSerializationWrapper::from_annotator(self);
        Ok((
            wrapper
                .inner
                .chains
                .clone()
                .iter()
                .map(move |s| s.to_string())
                .collect(),
            wrapper.inner.scheme.to_string(),
        ))
    }
}

#[derive(Serialize, Deserialize)]
pub(crate) struct AnnotatorSerializationWrapper {
    pub(crate) inner: Annotator,
}

impl AnnotatorSerializationWrapper {
    pub(crate) fn from_annotator(annotator: &Annotator) -> Self {
        AnnotatorSerializationWrapper {
            inner: annotator.clone(),
        }
    }
}

#[pymodule]
fn _internal(_py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<Annotator>()?;
    Ok(())
}
