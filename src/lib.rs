//! High-performance antibody and TCR sequence numbering.
//!
//! `immunum` numbers antibody and T-cell receptor (TCR) variable domain sequences
//! using Needleman-Wunsch semi-global alignment against position-specific scoring
//! matrices (PSSM) built from consensus sequences with BLOSUM62 substitution scores.
//! Chain type is detected automatically by aligning against all requested chains and
//! selecting the best match.
//!
//! # Supported chains and schemes
//!
//! | Chain | Type | Schemes |
//! |-------|------|---------|
//! | IGH | Antibody heavy | IMGT, Kabat |
//! | IGK | Antibody kappa | IMGT, Kabat |
//! | IGL | Antibody lambda | IMGT, Kabat |
//! | TRA | TCR alpha | IMGT |
//! | TRB | TCR beta | IMGT |
//! | TRG | TCR gamma | IMGT |
//! | TRD | TCR delta | IMGT |
//!
//! # Quick start
//!
//! ```rust
//! use immunum::{Annotator, Chain, Scheme};
//!
//! // Create an annotator for all antibody chains with IMGT numbering
//! let annotator = Annotator::new(&[Chain::IGH, Chain::IGK, Chain::IGL], Scheme::IMGT, None).unwrap();
//!
//! let sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";
//! let result = annotator.number(sequence).unwrap();
//!
//! println!("Chain:      {}", result.chain);       // IGH
//! println!("Confidence: {:.2}", result.confidence);
//!
//! // Iterate over (amino acid, IMGT position) pairs
//! for (aa, pos) in sequence.chars().zip(result.positions.iter()) {
//!     println!("{aa} -> {pos}");
//! }
//! ```
//!
//! # Key types
//!
//! - [`Annotator`] — main entry point; holds loaded scoring matrices and numbers sequences
//! - [`NumberingResult`] — the output of [`Annotator::number`]: positions, chain, confidence
//! - [`Position`] — an IMGT/Kabat position such as `111` or `111A`
//! - [`Chain`] — chain type (`IGH`, `IGK`, `IGL`, `TRA`, `TRB`, `TRG`, `TRD`)
//! - [`Scheme`] — numbering scheme (`IMGT` or `Kabat`)
//!
//! # Modules
//!
//! - [`annotator`] — high-level annotation API
//! - [`alignment`] — Needleman-Wunsch semi-global alignment
//! - [`numbering`] — IMGT and Kabat numbering rules
//! - [`scoring`] — position-specific scoring matrices
//! - [`io`] — FASTA parsing and TSV/JSON/JSONL output
//! - [`types`] — core domain types
//! - [`error`] — error types

pub mod alignment;
pub mod annotator;
pub mod error;
pub mod numbering;
pub mod scoring;
pub mod types;

#[cfg(not(target_arch = "wasm32"))]
pub mod io;
#[cfg(not(target_arch = "wasm32"))]
pub mod validation;

pub use numbering::imgt;
pub use numbering::kabat;

pub use alignment::{align_full, Alignment};
pub use annotator::{Annotator, NumberingResult, SegmentResult, DEFAULT_MIN_CONFIDENCE};
pub use error::{Error, Result};
pub use scoring::ScoringMatrix;
pub use types::{Chain, Insertion, NumberingRule, Position, Region, Scheme};

#[cfg(not(target_arch = "wasm32"))]
pub use io::{read_fasta, read_input, NumberedRecord, OutputFormat, Record};
#[cfg(not(target_arch = "wasm32"))]
pub use validation::{load_validation_csv, validate_entry, ValidationEntry, ValidationResult};

#[cfg(any(feature = "python", feature = "polars"))]
mod python;

#[cfg(feature = "polars")]
mod polars;

#[cfg(feature = "wasm")]
mod wasm;
