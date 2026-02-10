//! Error types for Immunum2

use thiserror::Error;

/// Result type alias for Immunum2 operations
pub type Result<T> = std::result::Result<T, Error>;

/// Errors that can occur during sequence numbering and alignment
#[derive(Debug, Error)]
pub enum Error {
    #[error("Invalid chain type: {0}")]
    InvalidChain(String),

    #[error("Invalid numbering scheme: {0}")]
    InvalidScheme(String),

    #[error("Consensus file parsing error: {0}")]
    ConsensusParseError(String),

    #[error("Alignment failed: {0}")]
    AlignmentError(String),

    #[error("Position mapping error: {0}")]
    PositionMappingError(String),

    #[error("Invalid position format: {0}")]
    InvalidPosition(String),

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Invalid sequence: {0}")]
    InvalidSequence(String),
}
