use std::str::FromStr;

use serde::{Deserialize, Serialize};

use postcard::{from_bytes, to_allocvec};
use pyo3::prelude::*;

use crate::annotator::{AnnotationResult, Annotator};
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

    #[pyo3(signature = (sequence))]
    pub fn _annotate(&self, sequence: &str) -> PyResult<AnnotationResult> {
        let result = self
            .annotate(sequence)
            .map_err(|_| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                    "Invalid sequence: {}",
                    sequence
                ))
            })
            .unwrap();
        Ok(result)
    }
}
