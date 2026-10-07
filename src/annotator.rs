//! High-level API for sequence annotation and chain detection
use std::cell::RefCell;
use std::ops::Range;

use crate::alignment::{align, AlignBuffer, AlignedPosition, Alignment};
use crate::error::{Error, Result};
use crate::numbering::{apply_numbering, segment as segment_positions};
use crate::scoring::ScoringMatrix;
use crate::types::{Chain, Position, Scheme};

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
    /// Numbered positions for the aligned region only (length == query_end - query_start + 1)
    pub positions: Vec<Position>,
    /// First aligned consensus position
    pub cons_start: usize,
    /// Last aligned consensus position
    pub cons_end: usize,
    /// Confidence score (normalized alignment score)
    pub confidence: f32,
    /// 0-based index of the first antibody residue in the query (0 when no prefix)
    pub query_start: usize,
    /// 0-based index of the last antibody residue in the query (query.len()-1 when no suffix)
    pub query_end: usize,
}

/// Result of segmenting a sequence into FR/CDR regions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SegmentResult {
    pub prefix: String,
    pub fr1: String,
    pub cdr1: String,
    pub fr2: String,
    pub cdr2: String,
    pub fr3: String,
    pub cdr3: String,
    pub fr4: String,
    pub postfix: String,
}

/// Default minimum confidence threshold for accepting a numbering result.
///
/// Based on empirical analysis of validated sequences:
/// - Antibody sequences (IGH/IGK/IGL): min ~0.51, median ~0.78-0.85
/// - TCR sequences (TRA/TRB/TRG/TRD): min ~0.28, median ~0.62-0.83
///
/// A threshold of 0.5 filters non-immunoglobulin sequences while retaining
/// all validated antibody sequences. Some low-scoring TCR sequences (notably
/// TCR-A p5=0.39, TCR-B p5=0.49) may fall below this threshold due to less
/// complete consensus data. Set to 0.0 to disable filtering.
pub const DEFAULT_MIN_CONFIDENCE: f32 = 0.5;

/// Minimum allowed input sequence length.
pub const MIN_SEQUENCE_LENGTH: usize = 30;

/// Maximum allowed input sequence length.
pub const MAX_SEQUENCE_LENGTH: usize = 10000;

/// Validate that `sequence` contains only standard amino acid characters
/// (case-insensitive) and that its length is within the allowed bounds.
fn validate_sequence(sequence: &str) -> Result<()> {
    let len = sequence.len();
    if len < MIN_SEQUENCE_LENGTH {
        return Err(Error::InvalidSequence(format!(
            "sequence length {} is below minimum {}",
            len, MIN_SEQUENCE_LENGTH
        )));
    }
    if len > MAX_SEQUENCE_LENGTH {
        return Err(Error::InvalidSequence(format!(
            "sequence length {} exceeds maximum {}",
            len, MAX_SEQUENCE_LENGTH
        )));
    }
    for (i, b) in sequence.bytes().enumerate() {
        if !b.is_ascii_alphabetic() {
            return Err(Error::InvalidSequence(format!(
                "invalid character {:?} at position {i}",
                b as char
            )));
        }
    }
    Ok(())
}

thread_local! {
    static ALIGN_BUFFER: RefCell<AlignBuffer> = RefCell::new(AlignBuffer::new());
}

/// Annotator for numbering sequences
#[cfg_attr(
    feature = "python",
    pyclass(name = "_Annotator", module = "immunum._internal", unsendable)
)]
#[cfg_attr(feature = "wasm", wasm_bindgen::prelude::wasm_bindgen(skip_typescript))]
#[derive(Clone, Serialize, Deserialize)]
pub struct Annotator {
    pub(crate) matrices: Vec<(Chain, ScoringMatrix)>,
    pub(crate) scheme: Scheme,
    pub(crate) chains: Vec<Chain>,
    pub(crate) min_confidence: f32,
}

impl Annotator {
    pub fn new(chains: &[Chain], scheme: Scheme, min_confidence: Option<f32>) -> Result<Self> {
        if chains.is_empty() {
            return Err(Error::InvalidChain("chains cannot be empty".to_string()));
        }

        for &chain in chains {
            scheme.validate_chain(chain)?;
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
            min_confidence: min_confidence.unwrap_or(DEFAULT_MIN_CONFIDENCE),
        })
    }

    /// Number a sequence by aligning to the configured chain types and applying the numbering scheme
    pub fn number(&self, sequence: &str) -> Result<NumberingResult> {
        validate_sequence(sequence)?;
        self.best_domain(sequence)?.number(self.scheme)
    }

    /// Segment a sequence into FR/CDR regions
    pub fn segment(&self, sequence: &str) -> Result<SegmentResult> {
        self.number(sequence)?.segment(sequence)
    }

    /// Every variable domain in `sequence`, ordered by position.
    ///
    /// Aligns the sequence and keeps the best alignment's domain, then searches the residues before and
    /// after it the same way, each on its own, so every domain is the best alignment of the residues
    /// around it that no other domain holds. A part stops being searched once it is shorter than
    /// [`MIN_SEQUENCE_LENGTH`] or its best alignment falls below the minimum confidence. A domain
    /// shorter than [`MIN_SEQUENCE_LENGTH`] is not reported, but the residues around it are searched.
    pub fn domains(&self, sequence: &str) -> Result<Vec<Domain>> {
        validate_sequence(sequence)?;
        let mut domains = Vec::new();
        let mut parts = Vec::new();
        parts.push(0..sequence.len());
        while let Some(part) = parts.pop() {
            if part.len() < MIN_SEQUENCE_LENGTH {
                continue;
            }
            let domain = self.aligned_domain(sequence, part.clone())?;
            if domain.confidence < self.min_confidence {
                continue;
            }
            parts.push(part.start..domain.query_start);
            parts.push(domain.query_end + 1..part.end);
            if domain.query_end + 1 - domain.query_start >= MIN_SEQUENCE_LENGTH {
                domains.push(domain);
            }
        }
        domains.sort_by_key(|domain| domain.query_start);
        for i in 1..domains.len() {
            domains[i - 1].tail_limit = domains[i].query_start;
        }
        Ok(domains)
    }

    fn best_domain(&self, sequence: &str) -> Result<Domain> {
        let domain = self.aligned_domain(sequence, 0..sequence.len())?;
        if domain.confidence < self.min_confidence {
            return Err(Error::LowConfidence {
                confidence: domain.confidence,
                threshold: self.min_confidence,
            });
        }
        Ok(domain)
    }

    // The domain of the best alignment of `sequence[part]`, whatever its confidence, placed in `sequence`
    fn aligned_domain(&self, sequence: &str, part: Range<usize>) -> Result<Domain> {
        let (chain, alignment) = self.get_best_alignment(&sequence[part.clone()])?;
        let confidence = if alignment.max_confidence_score > 0.0 {
            (alignment.confidence_score / alignment.max_confidence_score).clamp(0.0, 1.0)
        } else {
            0.0
        };
        Ok(Domain::new(
            chain,
            alignment,
            confidence,
            part.start,
            sequence.len(),
        ))
    }

    /// Align the sequence to all loaded chain types and return the best match
    /// If multiple chains were provided during initialization, this will align to all
    /// of them and return the best match. If only one chain was provided, it will
    /// align to that chain directly.
    fn get_best_alignment(&self, sequence: &str) -> Result<(Chain, Alignment)> {
        ALIGN_BUFFER.with_borrow_mut(|buf| {
            // Align to all loaded chain types and find best match by raw alignment score
            let mut best: Option<(Chain, Alignment)> = None;
            for (chain, matrix) in &self.matrices {
                let alignment = align(sequence, &matrix.positions, Some(&mut *buf));
                let is_better = match &best {
                    Some((_, prev)) => alignment.score > prev.score,
                    None => true,
                };
                if is_better {
                    best = Some((*chain, alignment));
                }
            }
            best.ok_or_else(|| {
                Error::AlignmentError("failed to align to any chain type".to_string())
            })
        })
    }
}

/// A variable domain found in a sequence: the alignment that found it, before a numbering scheme is
/// applied. Number it under as many schemes as needed; none aligns again.
#[derive(Debug, Clone)]
pub struct Domain {
    chain: Chain,
    confidence: f32,
    query_start: usize,
    query_end: usize,
    cons_start: usize,
    cons_end: usize,
    // One alignment state per residue in `query_start..=query_end`
    states: Vec<AlignedPosition>,
    // First index the AHo tail residue may not reach: the next domain's start, or the sequence end.
    tail_limit: usize,
}

impl Domain {
    // `alignment` aligned the residues of the searched sequence from `offset` on
    fn new(
        chain: Chain,
        alignment: Alignment,
        confidence: f32,
        offset: usize,
        tail_limit: usize,
    ) -> Self {
        Self {
            chain,
            confidence,
            query_start: offset + alignment.query_start,
            query_end: offset + alignment.query_end,
            cons_start: alignment.cons_start as usize,
            cons_end: alignment.cons_end as usize,
            states: alignment.positions[alignment.query_start..=alignment.query_end].to_vec(),
            tail_limit,
        }
    }

    pub fn chain(&self) -> Chain {
        self.chain
    }

    /// Confidence of this domain's alignment
    pub fn confidence(&self) -> f32 {
        self.confidence
    }

    /// 0-based index of the domain's first residue in the searched sequence
    pub fn query_start(&self) -> usize {
        self.query_start
    }

    /// 0-based index of the domain's last residue
    pub fn query_end(&self) -> usize {
        self.query_end
    }

    /// First aligned consensus position
    pub fn cons_start(&self) -> usize {
        self.cons_start
    }

    /// Last aligned consensus position
    pub fn cons_end(&self) -> usize {
        self.cons_end
    }

    /// This domain numbered under `scheme`, from the alignment that found it.
    pub fn number(&self, scheme: Scheme) -> Result<NumberingResult> {
        scheme.validate_chain(self.chain)?;

        let mut positions = apply_numbering(&self.states, scheme, self.chain);
        let mut query_end = self.query_end;

        // AHo light chains carry one extra C-terminal position (149) beyond the IMGT-numbered
        // region: IMGT ends light chains at 127 -> AHo 148, so the 149 residue has no IMGT
        // state and is appended here when a residue follows, matching ANARCI's number_aho tail
        // rule. Heavy chains populate IMGT 128 -> AHo 149 directly and need no append.
        if scheme == Scheme::Aho
            && matches!(self.chain, Chain::IGK | Chain::IGL)
            && positions.last() == Some(&Position::new(148))
            && query_end + 1 < self.tail_limit
        {
            positions.push(Position::new(149));
            query_end += 1;
        }

        Ok(NumberingResult {
            chain: self.chain,
            scheme,
            positions,
            cons_start: self.cons_start,
            cons_end: self.cons_end,
            confidence: self.confidence,
            query_start: self.query_start,
            query_end,
        })
    }
}

impl NumberingResult {
    /// The FR/CDR split of the residues this result numbered, without numbering again. `sequence` is
    /// the whole sequence that was numbered or searched, not a domain's own residues.
    pub fn segment(&self, sequence: &str) -> Result<SegmentResult> {
        let aligned_seq = sequence
            .get(self.query_start..=self.query_end)
            .ok_or_else(|| {
                Error::InvalidSequence(format!(
                    "numbered residues {}..={} lie outside a sequence of length {}",
                    self.query_start,
                    self.query_end,
                    sequence.len()
                ))
            })?;

        let mut map = segment_positions(&self.positions, aligned_seq, self.scheme, self.chain);

        Ok(SegmentResult {
            prefix: map.remove("prefix").unwrap_or_default(),
            fr1: map.remove("fr1").unwrap_or_default(),
            cdr1: map.remove("cdr1").unwrap_or_default(),
            fr2: map.remove("fr2").unwrap_or_default(),
            cdr2: map.remove("cdr2").unwrap_or_default(),
            fr3: map.remove("fr3").unwrap_or_default(),
            cdr3: map.remove("cdr3").unwrap_or_default(),
            fr4: map.remove("fr4").unwrap_or_default(),
            postfix: map.remove("postfix").unwrap_or_default(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::ALL_CHAINS;

    #[test]
    fn test_create_annotator() {
        let annotator = Annotator::new(ALL_CHAINS, Scheme::IMGT, None).unwrap();
        assert_eq!(annotator.matrices.len(), 7);
    }

    #[test]
    fn test_create_annotator_with_chains() {
        let annotator = Annotator::new(&[Chain::IGH, Chain::IGK], Scheme::IMGT, None).unwrap();
        assert_eq!(annotator.matrices.len(), 2);
    }

    #[test]
    fn test_number_igh_sequence() {
        let annotator = Annotator::new(ALL_CHAINS, Scheme::IMGT, None).unwrap();

        // Known IGH sequence
        let sequence =
            "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";

        let result = annotator.number(sequence).unwrap();

        // Should detect as IGH
        assert_eq!(result.chain, Chain::IGH);
        assert_eq!(result.scheme, Scheme::IMGT);
        assert!(result.confidence > 0.0 && result.confidence <= 1.0);
        assert_eq!(
            result.positions.len(),
            result.query_end - result.query_start + 1
        );
    }

    #[test]
    fn test_number_with_single_chain() {
        let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT, None).unwrap();
        let sequence =
            "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";

        let result = annotator.number(sequence).unwrap();
        assert_eq!(result.chain, Chain::IGH);
    }

    #[test]
    fn test_empty_sequence() {
        let annotator = Annotator::new(ALL_CHAINS, Scheme::IMGT, None).unwrap();
        let result = annotator.number("");
        assert!(result.is_err());
    }

    // Full IGH from the task description (FR1 through FR4)
    const FULL_IGH: &str = "EVQLVESGGGLVQPGGSLRLSCAASGFNVSYSSIHWVRQAPGKGLEWVAYIYPSSGYTSYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARSYSTKLAMDYWGQGTLVTVSS";

    // PDB 1EFQ chain A, from fixtures/validation/ab_K_imgt.csv.
    const KAPPA: &str = "DIVMTQSPDSLAVSLGERATINCKSSQSVLYSSNSKNYLAWYQDKPGQPPKLLIYWASTRESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQYYSTPYSFGQGTKLEIK";
    // The human kappa constant region. AHo numbers the residue after a kappa domain as position 149.
    const KAPPA_CONSTANT: &str = "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC";
    const LINKER: &str = "GGGGSGGGGSGGGGS";
    const ANTIBODY_CHAINS: &[Chain] = &[Chain::IGH, Chain::IGK, Chain::IGL];

    fn spans(sequence: &str) -> Vec<(usize, usize, Chain)> {
        Annotator::new(ANTIBODY_CHAINS, Scheme::IMGT, None)
            .unwrap()
            .domains(sequence)
            .unwrap()
            .iter()
            .map(|domain| (domain.query_start, domain.query_end, domain.chain))
            .collect()
    }

    #[test]
    fn a_domain_numbers_like_number_in_every_scheme() {
        let kappa_with_constant = format!("{KAPPA}{KAPPA_CONSTANT}");
        let heavy_with_one = format!("{FULL_IGH}A");
        let kappa_with_one = format!("{KAPPA}A");
        for sequence in [
            FULL_IGH,
            KAPPA,
            kappa_with_constant.as_str(),
            heavy_with_one.as_str(),
            kappa_with_one.as_str(),
        ] {
            for scheme in [
                Scheme::IMGT,
                Scheme::Kabat,
                Scheme::Chothia,
                Scheme::Martin,
                Scheme::Aho,
            ] {
                let annotator = Annotator::new(ANTIBODY_CHAINS, scheme, None).unwrap();
                let expected = annotator.number(sequence).unwrap();
                let domains = annotator.domains(sequence).unwrap();
                assert_eq!(domains.len(), 1, "{scheme} {sequence}");
                let got = domains[0].number(scheme).unwrap();
                assert_eq!(got.positions, expected.positions, "{scheme} {sequence}");
                assert_eq!(
                    (got.query_start, got.query_end),
                    (expected.query_start, expected.query_end)
                );
                assert_eq!(
                    (got.chain, got.cons_start, got.cons_end),
                    (expected.chain, expected.cons_start, expected.cons_end)
                );
                assert_eq!(got.confidence, expected.confidence);
            }
        }
    }

    #[test]
    fn aho_numbers_the_residue_after_a_kappa_domain() {
        let sequence = format!("{KAPPA}{KAPPA_CONSTANT}");
        let annotator = Annotator::new(ANTIBODY_CHAINS, Scheme::IMGT, None).unwrap();
        let domain = &annotator.domains(&sequence).unwrap()[0];
        assert_eq!(domain.number(Scheme::IMGT).unwrap().query_end, 112);
        assert_eq!(domain.number(Scheme::Aho).unwrap().query_end, 113);
    }

    #[test]
    fn finds_both_domains_of_an_scfv() {
        assert_eq!(
            spans(&format!("{FULL_IGH}{LINKER}{KAPPA}")),
            vec![(0, 118, Chain::IGH), (134, 246, Chain::IGK)]
        );
    }

    #[test]
    fn a_kappa_domain_keeps_its_aho_tail_out_of_the_next_domain() {
        let sequence = format!("{KAPPA}{FULL_IGH}");
        let annotator = Annotator::new(ANTIBODY_CHAINS, Scheme::IMGT, None).unwrap();
        let domains = annotator.domains(&sequence).unwrap();
        assert_eq!(domains.len(), 2);
        let first = domains[0].number(Scheme::Aho).unwrap();
        assert!(first.query_end < domains[1].query_start);
    }

    // The residues after a domain are aligned on their own, so a domain that lost its N-terminal
    // residues can skip the consensus positions it lacks, as it can at the start of a sequence.
    #[test]
    fn an_n_truncated_domain_after_another_numbers_like_on_its_own() {
        let annotator = Annotator::new(ANTIBODY_CHAINS, Scheme::IMGT, None).unwrap();
        for start in [40, 60, 80] {
            let truncated = &FULL_IGH[start..];
            let domains = annotator
                .domains(&format!("{FULL_IGH}{truncated}"))
                .unwrap();
            assert_eq!(domains.len(), 2, "FULL_IGH[{start}..]");
            let got = domains[1].number(Scheme::IMGT).unwrap();
            let expected = annotator.number(truncated).unwrap();
            assert_eq!(got.positions, expected.positions, "FULL_IGH[{start}..]");
            assert_eq!(
                (got.query_start, got.query_end),
                (
                    FULL_IGH.len() + expected.query_start,
                    FULL_IGH.len() + expected.query_end
                )
            );
            assert_eq!(got.confidence, expected.confidence);
        }
    }

    #[test]
    fn finds_abutting_domains() {
        assert_eq!(
            spans(&format!("{FULL_IGH}{KAPPA}")),
            vec![(0, 118, Chain::IGH), (119, 231, Chain::IGK)]
        );
    }

    #[test]
    fn keeps_a_truncated_domain_beside_a_full_one() {
        assert_eq!(
            spans(&format!("{FULL_IGH}{}", &KAPPA[..57])),
            vec![(0, 118, Chain::IGH), (119, 175, Chain::IGK)]
        );
    }

    #[test]
    fn an_x_inside_a_domain_keeps_it_whole() {
        let sequence = format!("{}X{}", &FULL_IGH[..60], &FULL_IGH[61..]);
        assert_eq!(spans(&sequence), vec![(0, 118, Chain::IGH)]);
    }

    #[test]
    fn a_constant_region_holds_no_domain() {
        assert!(spans(KAPPA_CONSTANT).is_empty());
    }

    #[test]
    fn an_x_run_holds_no_domain() {
        assert!(spans(&"X".repeat(300)).is_empty());
    }

    #[test]
    fn segment_of_a_numbering_matches_segment() {
        let annotator = Annotator::new(ANTIBODY_CHAINS, Scheme::Kabat, None).unwrap();
        let from_numbering = annotator
            .number(FULL_IGH)
            .unwrap()
            .segment(FULL_IGH)
            .unwrap();
        let segmented = annotator.segment(FULL_IGH).unwrap();
        assert_eq!(
            (from_numbering.fr1, from_numbering.cdr3, from_numbering.fr4),
            (segmented.fr1, segmented.cdr3, segmented.fr4)
        );
    }

    #[test]
    fn segment_of_a_later_domain_needs_the_whole_sequence() {
        let scfv = format!("{FULL_IGH}{LINKER}{KAPPA}");
        let annotator = Annotator::new(ANTIBODY_CHAINS, Scheme::IMGT, None).unwrap();
        let kappa = annotator.domains(&scfv).unwrap()[1]
            .number(Scheme::IMGT)
            .unwrap();
        assert!(kappa.segment(&scfv).is_ok());
        assert!(kappa.segment(KAPPA).is_err());
    }

    #[test]
    fn an_annotator_is_shared_across_threads() {
        fn assert_sync<T: Send + Sync>() {}
        assert_sync::<Annotator>();
    }

    #[test]
    fn test_number_no_flanking_has_zero_query_start_end() {
        let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT, None).unwrap();
        let result = annotator.number(FULL_IGH).unwrap();
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, FULL_IGH.len() - 1);
        assert_eq!(result.positions.len(), FULL_IGH.len());
    }

    #[test]
    fn test_number_with_prefix() {
        let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT, None).unwrap();
        let prefix = "AAAAAA";
        let sequence = format!("{prefix}{FULL_IGH}");
        let result = annotator.number(&sequence).unwrap();
        assert_eq!(result.chain, Chain::IGH);
        assert_eq!(result.query_start, prefix.len());
        assert_eq!(result.query_end, sequence.len() - 1);
        assert_eq!(result.positions.len(), FULL_IGH.len());
    }

    #[test]
    fn test_number_with_suffix() {
        let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT, None).unwrap();
        let suffix = "AAAAAAA";
        let sequence = format!("{FULL_IGH}{suffix}");
        let result = annotator.number(&sequence).unwrap();
        assert_eq!(result.chain, Chain::IGH);
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, FULL_IGH.len() - 1);
        assert_eq!(result.positions.len(), FULL_IGH.len());
    }

    // Expected values from ANARCI: one residue past the domain is never part of it, except that AHo
    // numbers the residue after a light chain as position 149.
    #[test]
    fn test_one_trailing_residue_numbers_like_none() {
        let chains = [Chain::IGH, Chain::IGK, Chain::IGL];
        for scheme in [
            Scheme::IMGT,
            Scheme::Kabat,
            Scheme::Chothia,
            Scheme::Martin,
            Scheme::Aho,
        ] {
            let annotator = Annotator::new(&chains, scheme, None).unwrap();
            for domain in [FULL_IGH, KAPPA] {
                let bare = annotator.number(domain).unwrap();
                for residue in ["A", "G", "R", "W"] {
                    let sequence = format!("{domain}{residue}");
                    let with_tail = annotator.number(&sequence).unwrap();
                    let aho_light_tail = scheme == Scheme::Aho && bare.chain == Chain::IGK;
                    let mut expected = bare.positions.clone();
                    if aho_light_tail {
                        expected.push(Position::new(149));
                    }

                    assert_eq!(
                        with_tail.positions, expected,
                        "{scheme} {:?}+{residue}",
                        bare.chain
                    );
                    assert_eq!(
                        with_tail.query_end,
                        bare.query_start + expected.len() - 1,
                        "{scheme} {:?}+{residue}",
                        bare.chain
                    );
                    if !aho_light_tail {
                        assert_eq!(
                            annotator.segment(&sequence).unwrap().fr4,
                            annotator.segment(domain).unwrap().fr4,
                            "{scheme} {:?}+{residue}",
                            bare.chain
                        );
                    }
                }
            }
        }
    }

    /// Kabat segmentation, heavy and light. Guards the chain-specific region tables end to end:
    /// under Kabat, CDR-H2 is 50-65 (16 positions) while CDR-L2 is 50-56 (7), and light numbering
    /// stops at 107. A single shared table cannot produce both, which is what this catches.
    #[test]
    fn test_segment_kabat_heavy_and_light_differ() {
        let heavy_seq = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";
        let heavy = Annotator::new(&[Chain::IGH], Scheme::Kabat, None)
            .unwrap()
            .segment(heavy_seq)
            .unwrap();

        // Every residue lands in exactly one region, and nothing spills into prefix/postfix.
        let rebuilt = format!(
            "{}{}{}{}{}{}{}",
            heavy.fr1, heavy.cdr1, heavy.fr2, heavy.cdr2, heavy.fr3, heavy.cdr3, heavy.fr4
        );
        assert_eq!(
            rebuilt, heavy_seq,
            "Kabat heavy segments must reconstruct the input"
        );
        assert!(heavy.prefix.is_empty() && heavy.postfix.is_empty());

        // CDR-H1 is Kabat's five-residue 31-35, not the ten-residue AbM 26-35.
        assert!(
            heavy.cdr1.len() <= 7,
            "Kabat CDR-H1 should be ~5 residues (31-35, plus any 35A/35B), got {} in {:?}",
            heavy.cdr1.len(),
            heavy.cdr1
        );
        // CDR-H2 spans 50-65, so it is far longer than the seven-residue light CDR2.
        assert!(
            heavy.cdr2.len() >= 14,
            "Kabat CDR-H2 spans 50-65, expected >=14 residues, got {} in {:?}",
            heavy.cdr2.len(),
            heavy.cdr2
        );

        let light_seq = "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPPTFGQGTKVEIK";
        let light = Annotator::new(&[Chain::IGK], Scheme::Kabat, None)
            .unwrap()
            .segment(light_seq)
            .unwrap();
        let rebuilt = format!(
            "{}{}{}{}{}{}{}",
            light.fr1, light.cdr1, light.fr2, light.cdr2, light.fr3, light.cdr3, light.fr4
        );
        assert_eq!(
            rebuilt, light_seq,
            "Kabat light segments must reconstruct the input"
        );

        // CDR-L1 is 24-34: eleven positions, so clearly longer than CDR-H1.
        assert!(
            light.cdr1.len() >= 9,
            "Kabat CDR-L1 spans 24-34, expected >=9 residues, got {} in {:?}",
            light.cdr1.len(),
            light.cdr1
        );
        assert!(
            light.cdr2.len() <= 8,
            "Kabat CDR-L2 spans 50-56, expected <=8 residues, got {} in {:?}",
            light.cdr2.len(),
            light.cdr2
        );
        assert!(
            heavy.cdr2.len() > light.cdr2.len(),
            "Kabat CDR-H2 (50-65) must be longer than CDR-L2 (50-56)"
        );
    }

    /// Martin segmentation uses the AbM CDR definition, not Chothia's. On heavy chains AbM widens
    /// both loops -- H1 26-35 against Chothia's 26-32, H2 50-58 against 52-56 -- so the extra
    /// residues are exactly the ones Chothia hands to the flanking frameworks. On light chains AbM
    /// coincides with Kabat (24-34 / 50-56 / 89-97), which is what the light half checks.
    #[test]
    fn test_segment_martin_follows_abm_definition() {
        let heavy_seq = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";
        let martin = Annotator::new(&[Chain::IGH], Scheme::Martin, None)
            .unwrap()
            .segment(heavy_seq)
            .unwrap();
        let chothia = Annotator::new(&[Chain::IGH], Scheme::Chothia, None)
            .unwrap()
            .segment(heavy_seq)
            .unwrap();

        let rebuilt = format!(
            "{}{}{}{}{}{}{}{}{}",
            martin.prefix,
            martin.fr1,
            martin.cdr1,
            martin.fr2,
            martin.cdr2,
            martin.fr3,
            martin.cdr3,
            martin.fr4,
            martin.postfix
        );
        assert_eq!(
            rebuilt, heavy_seq,
            "Martin heavy segments must reconstruct the input"
        );

        // CDR-H1: AbM 26-35 = Chothia 26-32 plus 33, 34, 35, the first three Chothia FR2 residues.
        assert_eq!(
            martin.cdr1,
            format!("{}{}", chothia.cdr1, &chothia.fr2[..3]),
            "Martin CDR-H1 should extend Chothia's 26-32 to AbM's 26-35"
        );
        // CDR-H2: AbM 50-58 = Chothia 52-56 plus 50, 51 in front and 57, 58 behind.
        assert_eq!(
            martin.cdr2,
            format!(
                "{}{}{}",
                &chothia.fr2[chothia.fr2.len() - 2..],
                chothia.cdr2,
                &chothia.fr3[..2]
            ),
            "Martin CDR-H2 should span AbM's 50-58, not Chothia's 52-56"
        );
        // CDR-H3: AbM 95-102 opens one residue earlier than Chothia's 96-101 and closes one later.
        assert_eq!(
            martin.cdr3,
            format!(
                "{}{}{}",
                &chothia.fr3[chothia.fr3.len() - 1..],
                chothia.cdr3,
                &chothia.fr4[..1]
            ),
            "Martin CDR-H3 should span AbM's 95-102"
        );

        // AbM light is Kabat light, so both schemes must cut the same light chain identically.
        let light_seq = "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPPTFGQGTKVEIK";
        let martin_light = Annotator::new(&[Chain::IGK], Scheme::Martin, None)
            .unwrap()
            .segment(light_seq)
            .unwrap();
        let kabat_light = Annotator::new(&[Chain::IGK], Scheme::Kabat, None)
            .unwrap()
            .segment(light_seq)
            .unwrap();
        assert_eq!(
            (
                martin_light.cdr1.as_str(),
                martin_light.cdr2.as_str(),
                martin_light.cdr3.as_str()
            ),
            (
                kabat_light.cdr1.as_str(),
                kabat_light.cdr2.as_str(),
                kabat_light.cdr3.as_str()
            ),
            "AbM light (24-34 / 50-56 / 89-97) coincides with Kabat light"
        );
    }

    #[test]
    fn test_segment_igh_sequence() {
        let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT, None).unwrap();
        let sequence =
            "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";
        let segments = annotator.segment(sequence).unwrap();
        assert_eq!(segments.fr1, "QVQLVQSGAEVKRPGSSVTVSCKAS");
        assert_eq!(segments.cdr1, "GGSFSTYA");
        assert_eq!(segments.cdr3, "AREGTTGKPIGAFAH");
        assert_eq!(segments.fr4, "WGQGTLVTVSS");
        assert!(segments.prefix.is_empty());
        assert!(segments.postfix.is_empty());
    }

    #[test]
    fn test_number_with_both_flanking() {
        let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT, None).unwrap();
        let prefix = "AAAAAA";
        let suffix = "AAAAAAA";
        let sequence = format!("{prefix}{FULL_IGH}{suffix}");
        let result = annotator.number(&sequence).unwrap();
        assert_eq!(result.chain, Chain::IGH);
        assert_eq!(result.query_start, prefix.len());
        assert_eq!(result.query_end, prefix.len() + FULL_IGH.len() - 1);
        assert_eq!(result.positions.len(), FULL_IGH.len());
    }

    /// Truncated but productive camel VHH reads from the Observed Antibody Space (Li et al. 2017,
    /// bactrian camel, run SRR3544217). Each aligns such that a Kabat heavy CDR1 or CDR3 collapses to
    /// a single residue, which asks `number_with_rules` to delete every base position but one. Three
    /// Kabat `deletion_order` tables were one entry short of that and sliced out of bounds, so these
    /// panicked at `numbering.rs:175` under Kabat while numbering fine under IMGT.
    const TRUNCATED_VHH_READS: &[(&str, &str)] = &[
        ("CDR1", "GWFRQAPGKEREGGAYIYTSDGIARYSDSVKGRFTISVDGVKKILFLQMNELKAEDTATYYCASTGRSNDCGPAQKLLLHSARGGRDFGIWGQGTQVTVS"),
        ("CDR3", "QLVESGGGLVQPGGSLRLSCAATGFTFSNNWMHWVRQAPGKGLEWVASISRSGGNTDYADSVKGRFTISRDNAKNTLYLHLNSLKPEDTAMYYCTNWGQGTQVTVS"),
        ("CDR1", "GWFRQAPGKEREGVAFISSEGAPTYADSVQGRFTISRNVLPERLSLQMTRLKAEDTAMYYCALDPSWDGRRIVLHGTFAAWECPREERQAFGVWGLGTQVTVS"),
        ("CDR3", "QLVESGGGLVQPGGSLRLSCAASGLTFSSHAMSWVRQAPGKGLEWVSGITGGGTSYYADPVKGRFTISRDNAKNSVYLQLNSLKAEDSAMYYCAKWGQGTQVTVS"),
        ("CDR1", "TWVRQAPGKGLEWVSTINSGGDSTYYADSVKGRFTISQDSAKNILYLQMRSLKPEDTAMYYCAARSVGWCPLFEHWLGKRAYTPGGYFANWGQGTQVTVS"),
        ("CDR1", "GWFRQAPGKEREGVAVIHKNIYVASNTPGAVFYADSVKGRFTISRDSAKNTLYLQMNSLKPEDAAMYSCAADSRYASCGWLLDRFRDFAYRGQGTQVTVS"),
    ];

    #[test]
    fn numbers_truncated_reads_under_kabat() {
        let annotator = Annotator::new(&[Chain::IGH], Scheme::Kabat, None).unwrap();

        for (region, sequence) in TRUNCATED_VHH_READS {
            let result = annotator
                .number(sequence)
                .unwrap_or_else(|err| panic!("{region} read failed to number: {err}"));

            assert_eq!(result.scheme, Scheme::Kabat);
            assert_eq!(
                result.positions.len(),
                result.query_end - result.query_start + 1,
                "{region} read got {} positions for {} aligned residues",
                result.positions.len(),
                result.query_end - result.query_start + 1,
            );
        }
    }

    #[test]
    fn numbers_truncated_reads_under_imgt_too() {
        let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT, None).unwrap();
        for (region, sequence) in TRUNCATED_VHH_READS {
            annotator
                .number(sequence)
                .unwrap_or_else(|err| panic!("{region} read failed to number under IMGT: {err}"));
        }
    }
}
