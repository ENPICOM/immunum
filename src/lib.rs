pub mod annotation;
pub mod annotator;
pub mod consensus_scoring;
pub mod constants;
pub mod sequence;
pub mod insertion_numbering;
pub mod needleman_wunsch;
pub mod numbering_scheme_type;
pub mod prefiltering;
pub mod result;
pub mod schemes;
pub mod scoring_matrix;
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
pub use result::{AnnotationResult, OutputFormat};

pub use constants::{get_scoring_params, ScoringParams};
pub use numbering_scheme_type::NumberingScheme;
pub use scoring_matrix::ScoringMatrix;
pub use types::{Chain, RegionRange, Scheme};
