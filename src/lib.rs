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
//! | IGH | Antibody heavy | IMGT, Kabat, Chothia, Martin, AHo |
//! | IGK | Antibody kappa | IMGT, Kabat, Chothia, Martin, AHo |
//! | IGL | Antibody lambda | IMGT, Kabat, Chothia, Martin, AHo |
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
//! println!("Chain:      {}", result.chain);       // H
//! println!("Confidence: {:.2}", result.confidence);
//!
//! // Iterate over (IMGT position, amino acid) pairs
//! for (pos, aa) in result.residues(sequence).unwrap() {
//!     println!("{pos} -> {aa}");
//! }
//! ```
//!
//! # Key types
//!
//! - [`Annotator`] — main entry point; holds loaded scoring matrices and numbers sequences
//! - [`NumberingResult`] — the output of [`Annotator::number`]: positions, chain, confidence
//! - [`Position`] — a position such as `111` or `111A`
//! - [`Chain`] — chain type (`IGH`, `IGK`, `IGL`, `TRA`, `TRB`, `TRG`, `TRD`)
//! - [`Scheme`] — numbering scheme (`IMGT`, `Kabat`, `Chothia`, `Martin` or `Aho`)
//! - [`Error`] — a setup mistake, such as an unknown chain name; [`SequenceError`] — a sequence
//!   that couldn't be numbered
//!
//! # Modules
//!
//! - [`annotator`] — high-level annotation API
//! - [`alignment`] — Needleman-Wunsch semi-global alignment
//! - [`numbering`] — IMGT, Kabat, Chothia, Martin and AHo numbering rules
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

pub use numbering::aho;
pub use numbering::chothia;
pub use numbering::imgt;
pub use numbering::kabat;
pub use numbering::martin;

pub use alignment::{align, Alignment};
pub use annotator::{
    per_domain, Annotator, Domain, NumberingResult, SegmentResult, DEFAULT_MIN_CONFIDENCE,
};
pub use error::{Error, Result, SequenceError};
pub use scoring::ScoringMatrix;
pub use types::{scheme_supports_chain, Chain, Insertion, NumberingRule, Position, Region, Scheme};

#[cfg(not(target_arch = "wasm32"))]
pub use io::{read_fasta, read_input, NumberedRecord, OutputFormat, Record, SegmentedRecord};
#[cfg(not(target_arch = "wasm32"))]
pub use validation::{load_validation_csv, validate_entry, ValidationEntry, ValidationResult};

#[cfg(any(feature = "python", feature = "polars"))]
mod python;

#[cfg(feature = "polars")]
mod polars;

#[cfg(feature = "wasm")]
mod wasm;
