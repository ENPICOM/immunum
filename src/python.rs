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
struct AnnotatorSerializationWrapper {
    inner: Annotator,
}

impl AnnotatorSerializationWrapper {
    fn from_annotator(annotator: &Annotator) -> Self {
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

// // the code below has been taken from: https://github.com/MarcoGorelli/polars-plugins-tutorial/issues/75
// impl serde::Serialize for Annotator {
//     fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
//     where
//         S: serde::Serializer,
//     {
//         let wrapper = AnnotatorSerializationWrapper::from_annotator(self);
//         let bytes = to_allocvec(&wrapper.inner).map_err(serde::ser::Error::custom)?;
//         serializer.serialize_bytes(&bytes)
//     }
// }
//
// impl<'de> serde::Deserialize<'de> for Annotator {
//     fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
//     where
//         D: serde::Deserializer<'de>,
//     {
//         struct FieldVisitor;
//
//         impl<'de> serde::de::Visitor<'de> for FieldVisitor {
//             type Value = Annotator;
//             fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
//                 formatter.write_str(concat!(
//                     "a byte array containing bincode-serialized ",
//                     stringify!(Encoder),
//                     " data"
//                 ))
//             }
//             fn visit_bytes<E>(self, v: &[u8]) -> Result<Self::Value, E>
//             where
//                 E: serde::de::Error,
//             {
//                 let wrapper: AnnotatorSerializationWrapper =
//                     from_bytes(v).map_err(serde::de::Error::custom)?;
//
//                 Ok(Annotator {
//                     matrices: wrapper.inner.matrices,
//                     scheme: wrapper.inner.scheme,
//                     chains: wrapper.inner.chains,
//                 })
//             }
//
//             fn visit_byte_buf<E>(self, v: Vec<u8>) -> Result<Self::Value, E>
//             where
//                 E: serde::de::Error,
//             {
//                 self.visit_bytes(&v)
//             }
//         }
//         deserializer.deserialize_bytes(FieldVisitor)
//     }
// }
//
