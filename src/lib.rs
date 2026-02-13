//! Immunum2 - High-performance antibody and TCR sequence numbering
//!
//! This library provides tools for aligning immunoglobulin (antibody) and
//! T-cell receptor sequences to consensus sequences and numbering them
//! according to various schemes (IMGT, Kabat).

pub mod alignment;
pub mod annotator;
pub mod error;
pub mod numbering;
pub mod scoring;
pub mod types;
pub mod validation;

// Re-export for backwards compatibility
pub use numbering::imgt;
pub use numbering::kabat;

pub use alignment::{align, Alignment};
pub use annotator::{AnnotationResult, Annotator};
pub use error::{Error, Result};
pub use scoring::ScoringMatrix;
pub use types::{Chain, InsertionStyle, Position, Region, NumberingConfig, Scheme};
pub use validation::{load_validation_csv, validate_entry, ValidationEntry, ValidationResult};
