use postcard::{from_bytes, to_allocvec};

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::annotator::{Annotator, NumberingResult, SegmentResult};
use crate::numbering::{regions_for, SEGMENT_NAMES};
use crate::types::{Chain, Scheme};

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

// One dict per domain, or a single error dict when the sequence couldn't be searched
fn per_domain<'py, T>(
    py: Python<'py>,
    results: crate::Result<Vec<T>>,
    dict: impl Fn(crate::Result<T>) -> PyResult<Bound<'py, PyDict>>,
) -> PyResult<Bound<'py, PyList>> {
    let dicts = match results {
        Ok(results) => results
            .into_iter()
            .map(|result| dict(Ok(result)))
            .collect::<PyResult<Vec<_>>>()?,
        Err(e) => vec![dict(Err(e))?],
    };
    PyList::new(py, dicts)
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
        per_domain(py, self.number_domains(sequence), |result| {
            numbering_dict(py, sequence, result)
        })
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
        per_domain(py, self.segment_domains(sequence), |result| {
            segment_dict(py, result)
        })
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
