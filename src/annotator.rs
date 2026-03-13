//! High-level API for sequence annotation and chain detection
use std::cell::RefCell;

use crate::alignment::{align, AlignBuffer, Alignment};
use crate::error::{Error, Result};
use crate::numbering::apply_numbering;
use crate::scoring::ScoringMatrix;
use crate::types::{Chain, Position, Scheme, TCR_CHAINS};

#[cfg(feature = "python")]
use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

/// Result of numbering a sequence
#[cfg_attr(feature = "python", pyclass(get_all))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NumberingResult {
    /// Detected chain type
    pub chain: Chain,
    /// Numbering scheme used
    pub scheme: Scheme,
    /// Numbered positions
    pub positions: Vec<Position>,
    /// Start position in consensus
    pub start: usize,
    /// End position in consensus
    pub end: usize,
    /// Confidence score (normalized alignment score)
    pub confidence: f32,
}

/// Annotator for numbering sequences
#[cfg_attr(
    feature = "python",
    pyclass(name = "_Annotator", module = "immunum._internal", unsendable)
)]
#[derive(Serialize, Deserialize)]
pub struct Annotator {
    pub(crate) matrices: Vec<(Chain, ScoringMatrix)>,
    pub(crate) scheme: Scheme,
    pub(crate) chains: Vec<Chain>,
    /// Reusable alignment buffer to avoid per-alignment allocation
    #[serde(skip)]
    align_buf: RefCell<AlignBuffer>,
}

impl Clone for Annotator {
    fn clone(&self) -> Self {
        Self {
            matrices: self.matrices.clone(),
            scheme: self.scheme,
            chains: self.chains.clone(),
            align_buf: RefCell::new(AlignBuffer::new()),
        }
    }
}

impl Annotator {
    pub fn new(chains: &[Chain], scheme: Scheme) -> Result<Self> {
        if chains.is_empty() {
            return Err(Error::InvalidChain("chains cannot be empty".to_string()));
        }

        // Validate: Kabat only supported for antibody chains
        if scheme == Scheme::Kabat && chains.iter().any(|c| TCR_CHAINS.contains(c)) {
            return Err(Error::InvalidScheme(
                "Kabat scheme only supported for antibody chains (IGH, IGK, IGL)".to_string(),
            ));
        }

        let mut matrices = Vec::new();
        for &chain in chains {
            let matrix = ScoringMatrix::load(chain)?;
            matrices.push((chain, matrix));
        }

        Ok(Self {
            matrices,
            scheme,
            chains: chains.to_vec(),
            align_buf: RefCell::new(AlignBuffer::new()),
        })
    }

    /// Number a sequence by aligning to the configured chain types and applying the numbering scheme
    pub fn number(&self, sequence: &str) -> Result<NumberingResult> {
        if sequence.is_empty() {
            return Err(Error::InvalidSequence("empty sequence".to_string()));
        }

        let (chain, alignment) = self.get_best_alignment(sequence)?;
        let positions = apply_numbering(&alignment.positions, self.scheme, chain);
        let confidence = alignment.score / (sequence.len() as f32);

        Ok(NumberingResult {
            chain,
            scheme: self.scheme,
            positions,
            start: alignment.start_pos as usize,
            end: alignment.end_pos as usize,
            confidence,
        })
    }

    /// Align the sequence to all loaded chain types and return the best match
    /// If multiple chains were provided during initialization, this will align to all
    /// of them and return the best match. If only one chain was provided, it will
    /// align to that chain directly.
    pub fn get_best_alignment(&self, sequence: &str) -> Result<(Chain, Alignment)> {
        let mut buf = self.align_buf.borrow_mut();
        // Align to all loaded chain types and find best match by raw alignment score
        let mut best: Option<(Chain, Alignment)> = None;
        for (chain, matrix) in &self.matrices {
            if let Ok(alignment) = align(sequence, &matrix.positions, Some(&mut *buf)) {
                let is_better = match &best {
                    Some((_, prev)) => alignment.score > prev.score,
                    None => true,
                };
                if is_better {
                    best = Some((*chain, alignment));
                }
            }
        }
        best.ok_or_else(|| Error::AlignmentError("failed to align to any chain type".to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::ALL_CHAINS;

    #[test]
    fn test_create_annotator() {
        let annotator = Annotator::new(ALL_CHAINS, Scheme::IMGT).unwrap();
        assert_eq!(annotator.matrices.len(), 7);
    }

    #[test]
    fn test_create_annotator_with_chains() {
        let annotator = Annotator::new(&[Chain::IGH, Chain::IGK], Scheme::IMGT).unwrap();
        assert_eq!(annotator.matrices.len(), 2);
    }

    #[test]
    fn test_number_igh_sequence() {
        let annotator = Annotator::new(ALL_CHAINS, Scheme::IMGT).unwrap();

        // Known IGH sequence
        let sequence =
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN";

        let result = annotator.number(sequence).unwrap();

        // Should detect as IGH
        assert_eq!(result.chain, Chain::IGH);
        assert_eq!(result.scheme, Scheme::IMGT);
        assert!(result.confidence != 0.0);
        assert_eq!(result.positions.len(), sequence.len());
    }

    #[test]
    fn test_number_with_single_chain() {
        let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT).unwrap();
        let sequence =
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN";

        let result = annotator.number(sequence).unwrap();
        assert_eq!(result.chain, Chain::IGH);
    }

    #[test]
    fn test_empty_sequence() {
        let annotator = Annotator::new(ALL_CHAINS, Scheme::IMGT).unwrap();
        let result = annotator.number("");
        assert!(result.is_err());
    }
}
