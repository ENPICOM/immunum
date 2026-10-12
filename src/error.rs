//! Error types for immunum
//!
//! immunum fails in two ways, and every interface reports them differently on purpose:
//!
//! - [`Error`]: immunum was set up or called wrongly, such as an unknown chain name. Nothing can be
//!   done until the caller fixes it, so every interface raises it: Python raises `immunum.Error`,
//!   JavaScript throws an `Error`, Polars fails when the expression is built and the CLI exits with
//!   status 1.
//! - [`SequenceError`]: one sequence couldn't be numbered. A batch shouldn't stop for it, so every
//!   interface returns it in that sequence's result, as `error` with its [`SequenceError::kind`] as
//!   `error_kind`.

use thiserror::Error;

use crate::types::{Chain, Scheme};

/// Result of an immunum call: [`Error`] unless given, `Result<T, SequenceError>` for one sequence
pub type Result<T, E = Error> = std::result::Result<T, E>;

/// immunum was set up or called wrongly. Every interface raises it.
#[derive(Debug, Error, PartialEq)]
#[non_exhaustive]
pub enum Error {
    /// An unknown chain name, or no chains at all
    #[error("{0}")]
    InvalidChain(String),

    /// An unknown scheme name
    #[error("{0}")]
    InvalidScheme(String),

    /// A scheme asked to number a chain it has no rules for
    #[error("{scheme} scheme only supported for antibody chains (IGH, IGK, IGL), not {chain:?}")]
    UnsupportedChain { scheme: Scheme, chain: Chain },

    /// A minimum confidence outside `[0, 1]`
    #[error("min_confidence must be in [0, 1], got {0}")]
    InvalidMinConfidence(f32),

    /// A position that doesn't read as a number with an optional insertion letter
    #[error("invalid position: {0}")]
    InvalidPosition(String),

    /// A numbering paired with a sequence other than the one it numbered
    #[error("numbered residues {start}..={end} lie outside a sequence of length {length}")]
    WrongSequence {
        start: usize,
        end: usize,
        length: usize,
    },
}

impl Error {
    /// What went wrong, as a stable code every interface reports alongside the message
    pub fn kind(&self) -> &'static str {
        match self {
            Error::InvalidChain(_) => "invalid_chain",
            Error::InvalidScheme(_) => "invalid_scheme",
            Error::UnsupportedChain { .. } => "unsupported_chain",
            Error::InvalidMinConfidence(_) => "invalid_min_confidence",
            Error::InvalidPosition(_) => "invalid_position",
            Error::WrongSequence { .. } => "wrong_sequence",
        }
    }
}

/// One sequence couldn't be numbered. Every interface returns it in that sequence's result.
#[derive(Debug, Error, PartialEq)]
#[non_exhaustive]
pub enum SequenceError {
    /// Too short, too long, or holding a character that isn't a letter
    #[error("{0}")]
    InvalidSequence(String),

    /// The best alignment's confidence is below the annotator's minimum
    #[error("alignment confidence {confidence:.4} is below min_confidence {threshold:.4}")]
    LowConfidence { confidence: f32, threshold: f32 },

    /// A search for every domain found none: the best alignment is confident but covers fewer
    /// residues than a domain must
    #[error("domain length {length} is below minimum {minimum}")]
    DomainTooShort { length: usize, minimum: usize },
}

impl SequenceError {
    /// What went wrong, as a stable code every interface reports as `error_kind`
    pub fn kind(&self) -> &'static str {
        match self {
            SequenceError::InvalidSequence(_) => "invalid_sequence",
            SequenceError::LowConfidence { .. } => "low_confidence",
            SequenceError::DomainTooShort { .. } => "domain_too_short",
        }
    }
}
