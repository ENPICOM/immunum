pub mod annotation;
pub mod annotator;
pub mod consensus_scoring;
pub mod constants;
pub mod fastx;
pub mod insertion_numbering;
pub mod needleman_wunsch;
pub mod numbering;
pub mod numbering_scheme_type;
pub mod prefiltering;
pub mod result;
pub mod schemes;
pub mod scoring_matrix;
pub mod sequence_stream;
pub mod types;
// Binding modules
#[cfg(feature = "python")]
pub mod python_bindings;
#[cfg(feature = "wasm")]
pub mod wasm_bindings;

// Re-export public API functions for convenience
pub use schemes::get_scheme;

// New primary API
pub use annotator::Annotator;
pub use result::AnnotationResult;

pub use constants::{ScoringParams, get_scoring_params};
pub use numbering_scheme_type::{NumberingScheme, NumberingOutput};
pub use scoring_matrix::ScoringMatrix;
pub use types::{Chain, Scheme, RegionRange};

// Make the Python module available as the main entry point for the cdylib
#[cfg(feature = "python")]
use pyo3::prelude::*;

#[cfg(feature = "python")]
#[pymodule]
fn immunum(m: &Bound<'_, PyModule>) -> PyResult<()> {
    python_bindings::immunum(m)
}
