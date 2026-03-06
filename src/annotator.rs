//! High-level API for sequence annotation and chain detection
use crate::alignment::{align, Alignment};
use crate::error::{Error, Result};
use crate::numbering::{imgt::get_imgt_numbering, kabat::get_kabat_numbering};
use crate::scoring::ScoringMatrix;
use crate::types::{Chain, Position, Scheme};
use serde::{Deserialize, Serialize};

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Result of annotating a sequence with chain detection
#[cfg_attr(feature = "python", pyclass(get_all))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotationResult {
    /// Detected chain type
    pub chain: Chain,
    /// Alignment result
    pub alignment: Alignment,
    /// Confidence score (normalized alignment score)
    pub confidence: f32,
}

impl AnnotationResult {
    /// Get numbering for a specific scheme
    pub fn numbering(&self, scheme: Scheme) -> Vec<Position> {
        match scheme {
            Scheme::IMGT => get_imgt_numbering(&self.alignment),
            Scheme::Kabat => get_kabat_numbering(&self.alignment, self.chain),
        }
    }
}

/// Annotator for numbering sequences
#[cfg_attr(feature = "python", pyclass)]
#[derive(Serialize, Deserialize)]
pub struct Annotator {
    matrices: Vec<(Chain, ScoringMatrix)>,
}

impl Annotator {
    pub fn new(chains: &[Chain], scheme: Scheme) -> Result<Self> {
        // Validate: Kabat only supported for antibody chains
        if scheme == Scheme::Kabat && chains.iter().any(|c| c.is_tcr()) {
            return Err(Error::InvalidScheme(
                "Kabat scheme only supported for antibody chains (IGH, IGK, IGL)".to_string(),
            ));
        }

        let mut matrices = Vec::new();
        for &chain in chains {
            let matrix = ScoringMatrix::load(chain)?;
            matrices.push((chain, matrix));
        }

        Ok(Self { matrices })
    }

    /// Annotate a sequence using the configured numbering scheme
    ///
    /// If multiple chains were provided during initialization, this will align to all
    /// of them and return the best match. If only one chain was provided, it will
    /// align to that chain directly.
    ///
    /// The numbering returned will use the scheme specified during annotator creation.
    pub fn annotate(&self, sequence: &str) -> Result<AnnotationResult> {
        if sequence.is_empty() {
            return Err(Error::InvalidSequence("empty sequence".to_string()));
        }

        // Align to all loaded chain types and find best match
        self.matrices
            .iter()
            .filter_map(|(chain, matrix)| {
                align(sequence, &matrix.positions).ok().map(|alignment| {
                    let confidence = self.calculate_confidence(&alignment, sequence.len());
                    (*chain, alignment, confidence)
                })
            })
            .max_by(|(_, _, conf1), (_, _, conf2)| {
                conf1
                    .partial_cmp(conf2)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(chain, alignment, confidence)| AnnotationResult {
                chain,
                alignment,
                confidence,
            })
            .ok_or_else(|| Error::AlignmentError("failed to align to any chain type".to_string()))
    }

    fn calculate_confidence(&self, alignment: &Alignment, seq_len: usize) -> f32 {
        // Simple confidence: normalize score by sequence length
        // Positive score per residue = good alignment
        alignment.score / (seq_len as f32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ALL_CHAINS: &[Chain] = &[
        Chain::IGH,
        Chain::IGK,
        Chain::IGL,
        Chain::TRA,
        Chain::TRB,
        Chain::TRG,
        Chain::TRD,
    ];

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
    fn test_annotate_igh_sequence() {
        let annotator = Annotator::new(ALL_CHAINS, Scheme::IMGT).unwrap();

        // Known IGH sequence
        let sequence =
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN";

        let result = annotator.annotate(sequence).unwrap();

        // Should detect as IGH
        assert_eq!(result.chain, Chain::IGH);
        assert!(result.confidence != 0.0);

        let numbering = result.numbering(Scheme::IMGT);
        assert_eq!(numbering.len(), sequence.len());
    }

    #[test]
    fn test_annotate_with_single_chain() {
        let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT).unwrap();
        let sequence =
            "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN";

        let result = annotator.annotate(sequence).unwrap();
        assert_eq!(result.chain, Chain::IGH);
    }

    #[test]
    fn test_empty_sequence() {
        let annotator = Annotator::new(ALL_CHAINS, Scheme::IMGT).unwrap();
        let result = annotator.annotate("");
        assert!(result.is_err());
    }
}
