use postcard::{from_bytes, to_allocvec};
use std::str::FromStr;

use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::annotator::Annotator;
use crate::numbering::segment;
use crate::types::{Chain, Scheme};

#[pymethods]
impl Annotator {
    // python methods
    #[new]
    #[pyo3(signature = (chains, scheme, min_confidence=None))]
    pub fn init(
        chains: Vec<String>,
        scheme: String,
        min_confidence: Option<f32>,
    ) -> PyResult<Self> {
        let parsed_chains = chains
            .iter()
            .map(|chain| {
                Chain::from_str(chain).map_err(|_| {
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
        let annotator = Annotator::new(parsed_chains.as_slice(), parsed_scheme, min_confidence)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Failed to initialize annotator: {}",
                    e
                ))
            })?;

        Ok(annotator)
    }

    #[pyo3(signature = (sequence), name = "number")]
    pub fn _number<'py>(&self, py: Python<'py>, sequence: &str) -> PyResult<Bound<'py, PyDict>> {
        let result = self.number(sequence).map_err(|_| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Invalid sequence: {}",
                sequence
            ))
        })?;

        let numbering = PyDict::new(py);
        let aligned_seq = &sequence[result.query_start..=result.query_end];
        for (pos, ch) in result.positions.iter().zip(aligned_seq.chars()) {
            numbering.set_item(pos.to_string(), ch.to_string())?;
        }

        let dict = PyDict::new(py);
        dict.set_item("chain", result.chain.to_string())?;
        dict.set_item("scheme", result.scheme.to_string())?;
        dict.set_item("confidence", result.confidence)?;
        dict.set_item("numbering", numbering)?;
        Ok(dict)
    }

    #[pyo3(signature = (sequence), name = "segment")]
    pub fn _segment<'py>(&self, py: Python<'py>, sequence: &str) -> PyResult<Bound<'py, PyDict>> {
        let result = self.number(sequence).map_err(|_| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Invalid sequence: {}",
                sequence
            ))
        })?;

        let aligned_seq = &sequence[result.query_start..=result.query_end];
        let segments = segment(&result.positions, aligned_seq, result.scheme);

        let dict = PyDict::new(py);
        for (region, seq) in segments {
            dict.set_item(region, seq)?;
        }
        Ok(dict)
    }

    pub fn __setstate__(
        &mut self,
        state: &pyo3::Bound<'_, pyo3::types::PyBytes>,
    ) -> pyo3::PyResult<()> {
        let annotator: Annotator = from_bytes(state.as_bytes()).unwrap();
        self.matrices = annotator.matrices;
        self.scheme = annotator.scheme;
        self.chains = annotator.chains;
        self.min_confidence = annotator.min_confidence;
        Ok(())
    }

    pub fn __getstate__<'py>(
        &self,
        py: pyo3::Python<'py>,
    ) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::types::PyBytes>> {
        Ok(pyo3::types::PyBytes::new(py, &to_allocvec(&self).unwrap()))
    }

    pub fn __getnewargs__(&self) -> pyo3::PyResult<(Vec<String>, String)> {
        Ok((
            self.chains
                .clone()
                .iter()
                .map(move |s| s.to_string())
                .collect(),
            self.scheme.to_string().clone(),
        ))
    }
}

#[pymodule]
fn _internal(_py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<Annotator>()?;
    Ok(())
}
