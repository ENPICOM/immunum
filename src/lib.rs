pub mod annotation;
pub mod consensus_scoring;
pub mod constants;
pub mod fastx;
pub mod insertion_numbering;
pub mod needleman_wunsch;
pub mod numbering;
pub mod numbering_scheme_type;
pub mod prefiltering;
pub mod schemes;
pub mod scoring_matrix;
pub mod sequence_stream;
pub mod types;
// Binding modules
#[cfg(feature = "python")]
pub mod python_bindings;
#[cfg(feature = "wasm")]
pub mod wasm_bindings;

// Removed proc-macro dependency - now using manual FromStr implementations

// Re-export public API functions for convenience
pub use schemes::get_scheme;

pub use constants::{ScoringParams, get_scoring_params};
pub use numbering_scheme_type::{NumberingScheme, NumberingOutput};
pub use scoring_matrix::ScoringMatrix;
pub use types::{Chain, Scheme, RegionRange};

// Re-export the Python module for external use
#[cfg(feature = "python")]
pub use python_bindings::immunum;
