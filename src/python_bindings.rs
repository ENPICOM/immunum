#[cfg(feature = "python")]
use pyo3::prelude::*;

use crate::constants::ScoringParams;
use crate::numbering;
use crate::numbering_scheme_type::NumberingScheme;
use crate::schemes::get_scheme;
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
}

/// Python wrapper for NumberingScheme
#[cfg(feature = "python")]
#[pyclass(name = "NumberingScheme")]
pub struct PyNumberingScheme {
    inner: NumberingScheme,
}

#[cfg(feature = "python")]
#[pymethods]
impl PyNumberingScheme {
    #[staticmethod]
    #[pyo3(signature = (scheme, chain, params=None))]
    pub fn get_scheme(
        scheme: PyScheme,
        chain: PyChain,
        params: Option<PyScoringParams>,
    ) -> PyResult<Self> {
        let rust_scheme: Scheme = scheme.into();
        let rust_chain: Chain = chain.into();
        let scoring_params = params.map(|p| p.inner);

        let numbering_scheme = get_scheme(rust_scheme, rust_chain, scoring_params);
        Ok(PyNumberingScheme {
            inner: numbering_scheme,
        })
    }

    pub fn number_sequence(&self, sequence: &str) -> PyResult<String> {
        let output = self.inner.number_sequence(sequence.as_bytes());
        Ok(format!(
            "Identity: {:.2}%, Numbering: {:?}",
            output.identity * 100.0,
            output.numbering
        ))
    }
}

/// Batch processing function for multiple sequences
#[cfg(feature = "python")]
#[pyfunction]
pub fn number_sequences_batch(
    sequences: Vec<String>,
    scheme: PyScheme,
    chains: Vec<PyChain>,
) -> PyResult<Vec<String>> {
    let rust_scheme: Scheme = scheme.into();
    let rust_chains: Vec<Chain> = chains.into_iter().map(|c| c.into()).collect();

    let results: Vec<String> = sequences
        .iter()
        .map(|seq| numbering::number_sequence(seq, &rust_scheme, &rust_chains))
        .collect();

    Ok(results)
}

/// Immunum Python module
#[cfg(feature = "python")]
#[pymodule]
pub fn immunum(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(number_sequences_batch, m)?)?;
    m.add_class::<PyScoringParams>()?; // Exported as "ScoringParams"
    m.add_class::<PyNumberingScheme>()?; // Exported as "NumberingScheme"
    m.add_class::<PyScheme>()?; // Exported as "Scheme"
    m.add_class::<PyChain>()?; // Exported as "Chain"
    Ok(())
}
