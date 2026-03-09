//! immunum - High-performance antibody and TCR sequence numbering
//!
//! This library provides tools for aligning immunoglobulin (antibody) and
//! T-cell receptor sequences to consensus sequences and numbering them
//! according to various schemes (IMGT, Kabat).

pub mod alignment;
pub mod annotator;
pub mod error;
pub mod io;
pub mod numbering;
pub mod scoring;
pub mod types;
pub mod validation;

// Re-export for backwards compatibility
pub use numbering::imgt;
pub use numbering::kabat;

pub use alignment::{align, Alignment};
pub use annotator::{Annotator, NumberingResult};
pub use error::{Error, Result};
pub use io::{read_fasta, read_input, NumberedRecord, OutputFormat, Record};
pub use scoring::ScoringMatrix;
pub use types::{
    Chain, Insertion, NumberingRule, Position, Region, Scheme, ALL_CHAINS, IG_CHAINS, TCR_CHAINS,
};
pub use validation::{load_validation_csv, validate_entry, ValidationEntry, ValidationResult};

#[cfg(feature = "python")]
mod python;
#[cfg(feature = "python")]
use pyo3::prelude::*;

#[cfg(feature = "python")]
#[pymodule]
fn _internal(_py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_class::<annotator::Annotator>()?;
    Ok(())
}
