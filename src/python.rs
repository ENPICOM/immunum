use postcard::{from_bytes, to_allocvec};

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::annotator::{per_domain, Annotator, NumberingResult, SegmentResult};
use crate::numbering::{region_spans, SEGMENT_NAMES};
use crate::types::scheme_supports_chain;

// What `Annotator.number` returns for `sequence`: the numbering, or the error with every other field
// None
fn numbering_dict<'py>(
    py: Python<'py>,
    sequence: &str,
    result: crate::Result<NumberingResult>,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    match result {
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

// What `Annotator.segment` returns: the segments, or the error with every segment None
fn segment_dict(
    py: Python<'_>,
    result: crate::Result<SegmentResult>,
) -> PyResult<Bound<'_, PyDict>> {
    let dict = PyDict::new(py);
    match result {
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
        numbering_dict(py, sequence, self.number(sequence))
    }

    /// One `number` result per domain, or a single error result when `sequence` is invalid
    #[pyo3(signature = (sequence), name = "number_domains")]
    pub fn _number_domains<'py>(
        &self,
        py: Python<'py>,
        sequence: &str,
    ) -> PyResult<Bound<'py, PyList>> {
        let dicts = per_domain(self.number_domains(sequence))
            .into_iter()
            .map(|result| numbering_dict(py, sequence, result))
            .collect::<PyResult<Vec<_>>>()?;
        PyList::new(py, dicts)
    }

    #[pyo3(signature = (sequence), name = "segment")]
    pub fn _segment<'py>(&self, py: Python<'py>, sequence: &str) -> PyResult<Bound<'py, PyDict>> {
        segment_dict(py, self.segment(sequence))
    }

    /// One `segment` result per domain, or a single error result when `sequence` is invalid
    #[pyo3(signature = (sequence), name = "segment_domains")]
    pub fn _segment_domains<'py>(
        &self,
        py: Python<'py>,
        sequence: &str,
    ) -> PyResult<Bound<'py, PyList>> {
        let dicts = per_domain(self.segment_domains(sequence))
            .into_iter()
            .map(|result| segment_dict(py, result))
            .collect::<PyResult<Vec<_>>>()?;
        PyList::new(py, dicts)
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

fn invalid(e: crate::Error) -> PyErr {
    PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string())
}

/// Region boundaries as `{region: (start, end)}`, both ends inclusive, keyed by the region names
/// `segment` uses.
#[pyfunction]
fn _regions_for<'py>(py: Python<'py>, scheme: &str, chain: &str) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for (region, span) in region_spans(scheme, chain).map_err(invalid)? {
        dict.set_item(region, span)?;
    }
    Ok(dict)
}

/// Whether `scheme` numbers `chain`
#[pyfunction]
fn _scheme_supports_chain(scheme: &str, chain: &str) -> PyResult<bool> {
    scheme_supports_chain(scheme, chain).map_err(invalid)
}

#[pymodule]
fn _internal(_py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<Annotator>()?;
    m.add_function(wrap_pyfunction!(_regions_for, m)?)?;
    m.add_function(wrap_pyfunction!(_scheme_supports_chain, m)?)?;
    Ok(())
}
