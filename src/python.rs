use postcard::{from_bytes, to_allocvec};

use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::annotator::Annotator;
use crate::numbering::{regions_for, SEGMENT_NAMES};
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
        Annotator::from_names(chains.iter().map(String::as_str), &scheme, min_confidence)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
    }

    #[pyo3(signature = (sequence), name = "number")]
    pub fn _number<'py>(&self, py: Python<'py>, sequence: &str) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        match self.number(sequence) {
            Ok(result) => {
                let numbering = PyDict::new(py);
                for (pos, ch) in result
                    .residues(sequence)
                    .expect("`result` numbered `sequence`")
                {
                    numbering.set_item(pos.to_string(), ch.to_string())?;
                }
                dict.set_item("chain", result.chain.to_string())?;
                dict.set_item("scheme", result.scheme.to_string())?;
                dict.set_item("confidence", result.confidence)?;
                dict.set_item("numbering", numbering)?;
                dict.set_item("query_start", result.query_start)?;
                dict.set_item("query_end", result.query_end)?;
                dict.set_item("error", py.None())?;
            }
            Err(e) => {
                dict.set_item("chain", py.None())?;
                dict.set_item("scheme", py.None())?;
                dict.set_item("confidence", py.None())?;
                dict.set_item("numbering", py.None())?;
                dict.set_item("query_start", py.None())?;
                dict.set_item("query_end", py.None())?;
                dict.set_item("error", e.to_string())?;
            }
        }
        Ok(dict)
    }

    #[pyo3(signature = (sequence), name = "segment")]
    pub fn _segment<'py>(&self, py: Python<'py>, sequence: &str) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        match self.segment(sequence) {
            Ok(s) => {
                for (name, residues) in s.regions() {
                    dict.set_item(name, residues)?;
                }
                dict.set_item("error", py.None())?;
            }
            Err(e) => {
                for name in SEGMENT_NAMES {
                    dict.set_item(name, py.None())?;
                }
                dict.set_item("error", e.to_string())?;
            }
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

/// Region boundaries as `{region: (start, end)}`, both ends inclusive, keyed by the lowercase
/// region name that `segment` uses.
#[pyfunction]
fn _regions_for<'py>(py: Python<'py>, scheme: &str, chain: &str) -> PyResult<Bound<'py, PyDict>> {
    let invalid = |e: crate::Error| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string());
    let parsed_chain = chain.parse::<Chain>().map_err(invalid)?;
    let parsed_scheme = scheme.parse::<Scheme>().map_err(invalid)?;
    parsed_scheme
        .validate_chain(parsed_chain)
        .map_err(invalid)?;

    let dict = PyDict::new(py);
    for (region, span) in regions_for(parsed_scheme, parsed_chain).spans() {
        dict.set_item(region.to_string().to_lowercase(), span)?;
    }
    Ok(dict)
}

#[pymodule]
fn _internal(_py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<Annotator>()?;
    m.add_function(wrap_pyfunction!(_regions_for, m)?)?;
    Ok(())
}
