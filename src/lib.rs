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
pub use annotator::{Annotator, NumberingResult, DEFAULT_MIN_CONFIDENCE};
pub use error::{Error, Result};
pub use io::{read_fasta, read_input, NumberedRecord, OutputFormat, Record};
pub use scoring::ScoringMatrix;
pub use types::{
    Chain, Insertion, NumberingRule, Position, Region, Scheme, ALL_CHAINS, IG_CHAINS, TCR_CHAINS,
};
pub use validation::{load_validation_csv, validate_entry, ValidationEntry, ValidationResult};

#[cfg(any(feature = "python", feature = "polars"))]
mod python;

#[cfg(feature = "polars")]
mod polars;
