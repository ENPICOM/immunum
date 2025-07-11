pub mod consensus_scoring;
pub mod constants;
pub mod fastx;
pub mod insertion_numbering;
pub mod needleman_wunsch;
pub mod numbering;
pub mod schemes;
pub mod sequence_stream;
pub mod types;

pub mod annotation;

// Binding modules
mod numbering_scheme_methods;
#[cfg(feature = "python")]
pub mod python_bindings;
#[cfg(feature = "wasm")]
pub mod wasm_bindings;

// Import the procedural macro from the separate crate
pub use immunum_macros::ParseFromString;

// Re-export the Python module for external use
#[cfg(feature = "python")]
pub use python_bindings::immunum;
