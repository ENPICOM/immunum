//! High-level API for sequence annotation and chain detection
use std::cell::RefCell;
use std::ops::Range;

use crate::alignment::{align, AlignBuffer, AlignedPosition, Alignment};
use crate::error::{Error, Result, SequenceError};
use crate::numbering::{apply_numbering, segment as segment_positions, SEGMENT_NAMES};
use crate::scoring::ScoringMatrix;
use crate::types::{Chain, Position, Scheme};

#[cfg(feature = "python")]
use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

#[cfg_attr(feature = "python", pyclass(eq, eq_int))]
#[derive(Debug, Serialize, Clone, Deserialize, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward = 0,
    Reverse = 1,
}

#[cfg_attr(feature = "python", pyclass(eq))]
#[derive(Debug, Serialize, Clone, Deserialize, Copy, PartialEq, Eq)]
/// Describes frameshift including how many bases were consumed for it
pub enum FrameshiftType {
    Insertion(u8),
    Deletion(u8),
}

#[cfg_attr(feature = "python", pyclass(get_all))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Frameshift {
    pub span: (usize, usize), //in the query sequence
    pub type_frameshift: FrameshiftType,
    pub position: Position, //IMGT number
    pub repaired: bool,     //if it was repaired
}

#[cfg_attr(feature = "python", pyclass(get_all))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NucleotideNumbering {
    pub strand: Strand,
    pub frameshifts: Vec<Frameshift>,
    pub sequence_aa: String,
}

/// Result of numbering a sequence
#[cfg_attr(feature = "python", pyclass(get_all))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NumberingResult {
    /// Detected chain type
    pub chain: Chain,
    /// Numbering scheme used
    pub scheme: Scheme,
    /// Numbered positions for the aligned region only
    pub positions: Vec<Position>,
    /// First aligned consensus position
    pub cons_start: usize,
    /// Last aligned consensus position
    pub cons_end: usize,
    /// Confidence score (normalized alignment score)
    pub confidence: f32,
    /// 0-based index of the first antibody position in the query (0 when no prefix)
    pub query_start: usize,
    /// 0-based index of the last antibody position in the query (query.len()-1 when no suffix)
    pub query_end: usize,

    // NT additions
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nucleotide: Option<NucleotideNumbering>,
}

/// Same as SegmentResult without post/prefix
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SegmentResultFixless {
    pub fr1: String,
    pub cdr1: String,
    pub fr2: String,
    pub cdr2: String,
    pub fr3: String,
    pub cdr3: String,
    pub fr4: String,
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

    /// The regions' amino acids, frameshifts repaired. Nucleotide input only.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub translated: Option<SegmentResultFixless>,
}

impl SegmentResult {
    /// Each segment's residues with its name from [`SEGMENT_NAMES`], in sequence order
    pub fn regions(&self) -> [(&'static str, &str); 9] {
        let segments = [
            &self.prefix,
            &self.fr1,
            &self.cdr1,
            &self.fr2,
            &self.cdr2,
            &self.fr3,
            &self.cdr3,
            &self.fr4,
            &self.postfix,
        ];
        std::array::from_fn(|i| (SEGMENT_NAMES[i], segments[i].as_str()))
    }
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

const NT_VOCABULARY: [char; 6] = ['A', 'C', 'T', 'G', 'U', 'N'];

/// Validate that `sequence` contains only standard amino acid characters
/// (case-insensitive) and that its length is within the allowed bounds.
fn validate_aminoacid_sequence(sequence: &str) -> Result<(), SequenceError> {
    let len = sequence.len();
    if len < MIN_SEQUENCE_LENGTH {
        return Err(SequenceError::InvalidSequence(format!(
            "sequence length {} is below minimum {}",
            len, MIN_SEQUENCE_LENGTH
        )));
    }
    if len > MAX_SEQUENCE_LENGTH {
        return Err(SequenceError::InvalidSequence(format!(
            "sequence length {} exceeds maximum {}",
            len, MAX_SEQUENCE_LENGTH
        )));
    }
    for (i, b) in sequence.bytes().enumerate() {
        if !b.is_ascii_alphabetic() {
            return Err(SequenceError::InvalidSequence(format!(
                "invalid character {:?} at position {i}",
                b as char
            )));
        }
    }
    Ok(())
}

fn validate_nucleotide_sequence(sequence: &str) -> Result<(), SequenceError> {
    let len = sequence.len();

    for (i, b) in sequence.chars().enumerate() {
        if !NT_VOCABULARY.contains(&b) {
            return Err(SequenceError::InvalidSequence(format!(
                "invalid character {:?} at position {i}",
                b
            )));
        }
    }
    if len < MIN_SEQUENCE_LENGTH * 3 {
        return Err(SequenceError::InvalidSequence(format!(
            "sequence length {} is below minimum {}",
            len, MIN_SEQUENCE_LENGTH
        )));
    }
    if len > MAX_SEQUENCE_LENGTH * 3 {
        return Err(SequenceError::InvalidSequence(format!(
            "sequence length {} exceeds maximum {}",
            len, MAX_SEQUENCE_LENGTH
        )));
    }

    Ok(())
}

thread_local! {
    // Reused by every alignment on this thread. Kept between calls only up to `KEPT_ALIGN_CELLS`.
    static ALIGN_BUFFER: RefCell<AlignBuffer> = RefCell::new(AlignBuffer::new());
}

/// The most alignment matrix cells a thread keeps between calls: what a 1,000-residue query needs
/// against a 128-position consensus, about 650 KB. A longer query still aligns, but its matrices are
/// freed when its call returns.
const KEPT_ALIGN_CELLS: usize = 1_001 * 129;

/// Per-domain results as every interface returns them: one per domain from
/// [`Annotator::number_domains`] or [`Annotator::segment_domains`], or, when the sequence couldn't
/// be searched, its error as the single result.
pub fn per_domain<T>(results: Result<Vec<T>, SequenceError>) -> Vec<Result<T, SequenceError>> {
    match results {
        Ok(results) => results.into_iter().map(Ok).collect(),
        Err(e) => vec![Err(e)],
    }
}

/// Annotator for numbering sequences
#[cfg_attr(
    feature = "python",
    pyclass(name = "_Annotator", module = "immunum._internal")
)]
#[cfg_attr(feature = "wasm", wasm_bindgen::prelude::wasm_bindgen(skip_typescript))]
#[derive(Clone, Serialize, Deserialize)]
pub struct Annotator {
    pub(crate) matrices: Vec<(Chain, ScoringMatrix)>,
    pub(crate) scheme: Scheme,
    pub(crate) chains: Vec<Chain>,
    pub(crate) min_confidence: f32,
    pub(crate) nucletide_input: bool,
}

impl Annotator {
    pub fn new(chains: &[Chain], scheme: Scheme, min_confidence: Option<f32>) -> Result<Self> {
        Annotator::new_nt(chains, scheme, false, min_confidence)
    }

    pub fn new_nt(
        chains: &[Chain],
        scheme: Scheme,
        nucletide_input: bool,
        min_confidence: Option<f32>,
    ) -> Result<Self> {
        if chains.is_empty() {
            return Err(Error::InvalidChain("no chains given".to_string()));
        }
        if let Some(confidence) = min_confidence.filter(|c| !(0.0..=1.0).contains(c)) {
            return Err(Error::InvalidMinConfidence(confidence));
        }

        for &chain in chains {
            scheme.validate_chain(chain)?;
        }

        let matrices = chains
            .iter()
            .map(|&chain| (chain, ScoringMatrix::load(chain)))
            .collect();

        Ok(Self {
            matrices,
            scheme,
            chains: chains.to_vec(),
            min_confidence: min_confidence.unwrap_or(DEFAULT_MIN_CONFIDENCE),
            nucletide_input,
        })
    }

    /// An annotator from chain and scheme names, as a user writes them: see [`Chain::parse_names`]
    /// for the chains and [`Scheme`] for the scheme.
    pub fn from_names<'a>(
        chains: impl IntoIterator<Item = &'a str>,
        scheme: &str,
        min_confidence: Option<f32>,
    ) -> Result<Self> {
        Self::new(
            &Chain::parse_names(chains)?,
            scheme.parse()?,
            min_confidence,
        )
    }

    /// Number a sequence by aligning to the configured chain types and applying the numbering scheme
    pub fn number(&self, sequence: &str) -> Result<NumberingResult, SequenceError> {
        if self.nucletide_input {
            validate_nucleotide_sequence(sequence)?;
        } else {
            validate_aminoacid_sequence(sequence)?;
        }
        Ok(self.best_domain(sequence)?.numbered_as(self.scheme))
    }

    /// Segment a sequence into FR/CDR regions
    pub fn segment(&self, sequence: &str) -> Result<SegmentResult, SequenceError> {
        Ok(self
            .number(sequence)?
            .segment_within(sequence, 0..sequence.len()))
    }

    /// Every domain in `sequence` numbered under this annotator's scheme, ordered by position; see
    /// [`Annotator::domains`] for how they're found, and for the error when there is none.
    pub fn number_domains(&self, sequence: &str) -> Result<Vec<NumberingResult>, SequenceError> {
        Ok(self
            .domains(sequence)?
            .iter()
            .map(|domain| domain.numbered_as(self.scheme))
            .collect())
    }

    /// The FR/CDR split of every domain in `sequence`, ordered by position; see
    /// [`Annotator::domains`] for how they're found. Every residue lands in exactly one domain's
    /// segments: a domain's prefix holds the residues since the previous domain (or the start of the
    /// sequence), and only the last domain has the residues after it as its postfix. The domains'
    /// segments in order therefore rebuild `sequence`. See [`Annotator::domains`] for the error when
    /// there is no domain.
    pub fn segment_domains(&self, sequence: &str) -> Result<Vec<SegmentResult>, SequenceError> {
        let domains = self.number_domains(sequence)?;
        let mut start = 0;
        Ok(domains
            .iter()
            .enumerate()
            .map(|(i, domain)| {
                let end = if i + 1 == domains.len() {
                    sequence.len()
                } else {
                    domain.query_end + 1
                };
                let segments = domain.segment_within(sequence, start..end);
                start = end;
                segments
            })
            .collect())
    }

    /// Every variable domain in `sequence`, ordered by position; never empty.
    ///
    /// Aligns the sequence and keeps the best alignment's domain, then searches the residues before and
    /// after it the same way, each on its own, so every domain is the best alignment of the residues
    /// around it that no other domain holds. A part stops being searched once it is shorter than
    /// [`MIN_SEQUENCE_LENGTH`] or its best alignment falls below the minimum confidence. A domain
    /// shorter than [`MIN_SEQUENCE_LENGTH`] is not reported, but the residues around it are searched.
    ///
    /// When there is no domain, the error says why, from the best alignment of the whole sequence:
    /// [`SequenceError::LowConfidence`] when it falls below the minimum confidence, as
    /// [`Annotator::number`] reports it, or [`SequenceError::DomainTooShort`] when it is confident but
    /// too short and the residues around it hold no domain either.
    pub fn domains(&self, sequence: &str) -> Result<Vec<Domain>, SequenceError> {
        if self.nucletide_input {
            validate_nucleotide_sequence(sequence)?;
        } else {
            validate_aminoacid_sequence(sequence)?;
        }
        let best = self.best_domain(sequence)?;
        let best_length = best.query_end + 1 - best.query_start;

        // Every confident alignment found, with the part of the sequence it was the best of
        let mut found = vec![(best, 0..sequence.len())];
        let mut domains = Vec::new();
        while let Some((domain, part)) = found.pop() {
            for side in [
                part.start..domain.query_start,
                domain.query_end + 1..part.end,
            ] {
                if side.len() >= MIN_SEQUENCE_LENGTH {
                    let next = self.aligned_domain(sequence, side.clone());
                    if next.confidence >= self.min_confidence {
                        found.push((next, side));
                    }
                }
            }
            if domain.query_end + 1 - domain.query_start >= MIN_SEQUENCE_LENGTH {
                domains.push(domain);
            }
        }
        if domains.is_empty() {
            return Err(SequenceError::DomainTooShort {
                length: best_length,
                minimum: MIN_SEQUENCE_LENGTH,
            });
        }

        domains.sort_by_key(|domain| domain.query_start);
        for i in 1..domains.len() {
            domains[i - 1].tail_limit = domains[i].query_start;
        }
        Ok(domains)
    }

    fn best_domain(&self, sequence: &str) -> Result<Domain, SequenceError> {
        let domain = self.aligned_domain(sequence, 0..sequence.len());
        if domain.confidence < self.min_confidence {
            return Err(SequenceError::LowConfidence {
                confidence: domain.confidence,
                threshold: self.min_confidence,
            });
        }
        Ok(domain)
    }

    // The domain of the best alignment of `sequence[part]`, whatever its confidence, placed in `sequence`
    fn aligned_domain(&self, sequence: &str, part: Range<usize>) -> Domain {
        let (chain, alignment) = self.get_best_alignment(&sequence[part.clone()]);
        let confidence = if alignment.max_confidence_score > 0.0 {
            (alignment.confidence_score / alignment.max_confidence_score).clamp(0.0, 1.0)
        } else {
            0.0
        };
        Domain::new(chain, alignment, confidence, part.start, sequence.len())
    }

    /// Align the sequence to all loaded chain types and return the best match
    /// If multiple chains were provided during initialization, this will align to all
    /// of them and return the best match. If only one chain was provided, it will
    /// align to that chain directly.
    fn get_best_alignment(&self, sequence: &str) -> (Chain, Alignment) {
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
            buf.release_above(KEPT_ALIGN_CELLS);
            best.expect("`Annotator::new` rejects an empty chain list")
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

    /// This domain numbered under `scheme`, from the alignment that found it. An error when `scheme`
    /// doesn't number this domain's chain.
    pub fn number(&self, scheme: Scheme) -> Result<NumberingResult> {
        scheme.validate_chain(self.chain)?;
        Ok(self.numbered_as(scheme))
    }

    // This domain numbered under `scheme`, which numbers its chain
    fn numbered_as(&self, scheme: Scheme) -> NumberingResult {
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

        NumberingResult {
            chain: self.chain,
            scheme,
            positions,
            cons_start: self.cons_start,
            cons_end: self.cons_end,
            confidence: self.confidence,
            query_start: self.query_start,
            query_end,
            nucleotide: None,
        }
    }
}

impl NumberingResult {
    /// Each numbered position with its residue. `sequence` is the whole sequence that was numbered,
    /// flanking residues included. For nucleotide input the residues are the amino acids the
    /// codons translate to, frameshifts repaired, not the bases.
    pub fn residues<'a>(
        &'a self,
        sequence: &'a str,
    ) -> Result<impl Iterator<Item = (&'a Position, char)> + 'a> {
        let residues = match &self.nucleotide {
            Some(nucleotide) => nucleotide.sequence_aa.as_str(),
            None => self.numbered(sequence)?,
        };
        Ok(self.positions.iter().zip(residues.chars()))
    }

    /// The FR/CDR split of `sequence`, the whole sequence that was numbered. The residues before the
    /// numbered ones open the prefix and the residues after them close the postfix, so the regions in
    /// order rebuild `sequence`.
    pub fn segment(&self, sequence: &str) -> Result<SegmentResult> {
        self.numbered(sequence)?;
        Ok(self.segment_within(sequence, 0..sequence.len()))
    }

    // The FR/CDR split of `sequence[span]`, which holds the residues this result numbered: the
    // residues of `span` before them open the prefix and those after them close the postfix
    fn segment_within(&self, sequence: &str, span: Range<usize>) -> SegmentResult {
        let numbered = &sequence[self.query_start..=self.query_end];
        let mut map = segment_positions(&self.positions, numbered, self.scheme, self.chain);
        let [mut prefix, fr1, cdr1, fr2, cdr2, fr3, cdr3, fr4, mut postfix] =
            SEGMENT_NAMES.map(|name| map.remove(name).unwrap_or_default());
        prefix.insert_str(0, &sequence[span.start..self.query_start]);
        postfix.push_str(&sequence[self.query_end + 1..span.end]);

        SegmentResult {
            prefix,
            fr1,
            cdr1,
            fr2,
            cdr2,
            fr3,
            cdr3,
            fr4,
            postfix,
            translated: None,
        }
    }

    // The residues of `sequence` this result numbered
    fn numbered<'a>(&self, sequence: &'a str) -> Result<&'a str> {
        sequence
            .get(self.query_start..=self.query_end)
            .ok_or(Error::WrongSequence {
                start: self.query_start,
                end: self.query_end,
                length: sequence.len(),
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

    // The range was only checked by the Python wrapper; JavaScript, Polars and the CLI accepted any value.
    #[test]
    fn min_confidence_must_lie_in_zero_to_one() {
        for accepted in [None, Some(0.0), Some(1.0)] {
            assert!(Annotator::new(&[Chain::IGH], Scheme::IMGT, accepted).is_ok());
        }
        for rejected in [-0.01, 1.01, f32::NAN] {
            assert!(matches!(
                Annotator::new(&[Chain::IGH], Scheme::IMGT, Some(rejected)),
                Err(Error::InvalidMinConfidence(_))
            ));
        }
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

    // Why `sequence` holds no domain
    fn no_domain(sequence: &str) -> &'static str {
        Annotator::new(ANTIBODY_CHAINS, Scheme::IMGT, None)
            .unwrap()
            .domains(sequence)
            .unwrap_err()
            .kind()
    }

    #[test]
    fn a_constant_region_holds_no_domain() {
        assert_eq!(no_domain(KAPPA_CONSTANT), "low_confidence");
    }

    #[test]
    fn an_x_run_holds_no_domain() {
        assert_eq!(no_domain(&"X".repeat(300)), "low_confidence");
    }

    #[test]
    fn segment_domains_puts_every_residue_in_exactly_one_domain() {
        let annotator = Annotator::new(ANTIBODY_CHAINS, Scheme::IMGT, None).unwrap();
        let (leader, tag) = ("MKYLL", "HHHHHH");
        let scfv = format!("{leader}{FULL_IGH}{LINKER}{KAPPA}{tag}");
        let segments = annotator.segment_domains(&scfv).unwrap();

        let rebuilt: String = segments
            .iter()
            .flat_map(|s| s.regions())
            .map(|(_, residues)| residues)
            .collect();
        assert_eq!(rebuilt, scfv);
        let flanks: Vec<_> = segments
            .iter()
            .map(|s| (s.prefix.as_str(), s.postfix.as_str()))
            .collect();
        assert_eq!(flanks, [(leader, ""), (LINKER, tag)]);

        // Between the flanks, each domain splits as it does on its own
        for (segments, domain) in segments.iter().zip([FULL_IGH, KAPPA]) {
            let alone = annotator.segment(domain).unwrap();
            assert_eq!(segments.regions()[1..8], alone.regions()[1..8]);
        }
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

    #[test]
    fn residues_skip_the_flanking_residues() {
        let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT, None).unwrap();
        let sequence = format!("MGWSCIILFLVATATGVHSX{FULL_IGH}AAAAAAA");
        let flanked = annotator.number(&sequence).unwrap();
        let bare = annotator.number(FULL_IGH).unwrap();
        let got: Vec<_> = flanked.residues(&sequence).unwrap().collect();
        let expected: Vec<_> = bare.residues(FULL_IGH).unwrap().collect();
        assert_eq!(got, expected);
    }

    // Issue #58: residues the aligner leaves out were dropped instead of landing in prefix/postfix.
    #[test]
    fn segment_puts_flanking_residues_in_prefix_and_postfix() {
        let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT, None).unwrap();
        let sequence = "AAAAAQVQLQESGGGLVQPGGSLRLSCAASGFTFSNYKMNWVRQAPGKGLEWVSDISQSGASISYTGSVKGRFTISRDNAKNTLYLQMNSLKPEDTAVYYCARCPAPFTRDCFDVTSTTYAYRGQGTQVTVSSHHHHHHEPEA";
        let s = annotator.segment(sequence).unwrap();
        assert_eq!(s.prefix, "AAAAA");
        assert_eq!(s.postfix, "HHHHHHEPEA");
        assert_eq!(s.fr4, "RGQGTQVTVSS");
        let rebuilt = [
            &s.prefix, &s.fr1, &s.cdr1, &s.fr2, &s.cdr2, &s.fr3, &s.cdr3, &s.fr4, &s.postfix,
        ]
        .map(String::as_str)
        .concat();
        assert_eq!(rebuilt, sequence);
    }

    #[test]
    fn a_numbering_needs_the_sequence_it_numbered() {
        let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT, None).unwrap();
        let result = annotator.number(&format!("AAAAAA{FULL_IGH}")).unwrap();
        assert!(result.residues(FULL_IGH).is_err());
        assert!(result.segment(FULL_IGH).is_err());
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

    fn imgt_segments(sequence: &str) -> SegmentResult {
        Annotator::new(ANTIBODY_CHAINS, Scheme::IMGT, None)
            .unwrap()
            .segment(sequence)
            .unwrap()
    }

    // Mouse IGHV2-9*02 with G54V, a CDR3 and IGHJ4. IGHV2 has a 7-residue CDR2.
    #[test]
    fn a_variable_cdr2_residue_keeps_fr2_whole() {
        for residue in ["I", "N", "H", "Q", "W", "G"] {
            let sequence = format!("QVQLKESGPGLVAPSQSLSITCTVSGFSLTSYGVHWVRQPPGKGLEWLVV{residue}WAGGSTNYNSALMSRLSISKDNSKSQVFLKMNSLQTDDTAMYYCARDRGYYAMDYWGQGTSVTVSS");
            let segments = imgt_segments(&sequence);
            assert_eq!(segments.fr2, "VHWVRQPPGKGLEWLVV", "{residue}");
            assert_eq!(segments.cdr2, format!("{residue}WAGGST"), "{residue}");
        }
    }

    /// The only heavy chains in `fixtures/validation/ab_H_imgt.csv` with a real FR2 deletion,
    /// with ANARCI's IMGT numbering of FR2 (39-55). 5gjs_H has the same sequence as 5gjt_H.
    const FR2_DELETION_HEAVY_CHAINS: &[(&str, &str, &str)] = &[
        (
            "5gjt_H",
            "QVQLQESGPGLVKPSETLSLTCSVSGASISSYYWIWIRQPGKGLEWIGRFYTSGSPNYNPSLRSRVTMSVDTSKNQFSLKLTSVTAADTAVYYCAREEHITFGGVIVRYWGQGTLVTVSP",
            "39W 40I 41W 42I 43R 44Q 46P 47G 48K 49G 50L 51E 52W 53I 54G 55R",
        ),
        (
            "1mfe_H",
            "EVQVQQSGTVLARPGASVKMSCKASGYTFTNYWMHWIKQRPGGLEWIGAIYPGNSATFYNHKFRAKTKLTAVTSTITAYMELSSLTNEDSAVYYCTRGGHGYYGDYWGQGASLTVSS",
            "39M 40H 41W 42I 43K 44Q 45R 46P 47G 49G 50L 51E 52W 53I 54G 55A",
        ),
        (
            "1yc8_B",
            "VQLVESGGGSVQAGGSLRLSCAVSGSTYSPCTTGWVRQAPGGLEWVSSISSPGTIYYQDSVKGRFTISRDNAKNTVYLQMNSLQREDTGMYYCQIQCRSIREYWGQGTQVTVS",
            "39T 40G 41W 42V 43R 44Q 45A 46P 47G 49G 50L 51E 52W 53V 54S 55S",
        ),
    ];

    // Antibody FR2 has a conservation multiplier floor (see `conservation_multiplier_floor` in
    // build.rs), so a gap anywhere in it costs at least as much as at a conserved position. That
    // keeps FR2 whole when CDR2 is mutated, but it also makes the real FR2 deletions in these chains
    // more expensive: 5gjt_H's gap sits at 45, the one antibody FR2 position the floor raised once
    // the occupancy check was dropped. If the floor grows, these deletions are the first to leave
    // FR2 for a neighbouring CDR (at a floor of 5.0, 5gjt_H's moves to CDR1). The IGH validation
    // test would not notice: it asserts a 99% perfect rate, which leaves room for about 20
    // regressions.
    #[test]
    fn a_real_fr2_deletion_stays_in_fr2() {
        let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT, None).unwrap();
        for (header, sequence, expected_fr2) in FR2_DELETION_HEAVY_CHAINS {
            let result = annotator.number(sequence).unwrap();
            let fr2 = result
                .residues(sequence)
                .unwrap()
                .filter(|(position, _)| (39..=55).contains(&position.number))
                .map(|(position, residue)| format!("{position}{residue}"))
                .collect::<Vec<_>>()
                .join(" ");
            assert_eq!(fr2, *expected_fr2, "{header}");
        }
    }

    // Human IGKV1-39*01 and IGKJ1 with an extra S in FR1
    #[test]
    fn an_extra_fr1_residue_keeps_residue_1() {
        let segments = imgt_segments("DIQMTQSPSSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPWTFGQGTKVEIK");
        assert_eq!(segments.prefix, "");
        assert_eq!(segments.fr1, "DIQMTQSPSSSLSASVGDRVTITCRAS");
    }

    // Alpaca IGHV3-3*01 up to C104; the FR4 that follows each CDR3 lacks two J residues.
    const VHH_TO_CDR3: &str = "QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYC";
    const SHORT_FR4: &str = "WGQGTVSS";

    #[test]
    fn a_tag_after_a_short_fr4_stays_out_of_the_domain() {
        let cdr3 = "AASGYYSGTSYGPSAS";
        let segments = imgt_segments(&format!("{VHH_TO_CDR3}{cdr3}{SHORT_FR4}GGGGSHHHHHH"));
        assert_eq!(segments.cdr3, cdr3);
        assert!(segments.fr4.starts_with("WGQG"), "{}", segments.fr4);
        assert!(segments.postfix.ends_with("GGGGSHHHHHH"));
    }

    #[test]
    fn a_long_cdr3_before_a_tail_keeps_its_fr4() {
        let cdr3 = "AADPTRGLSYSGQWERPDYAS";
        let segments = imgt_segments(&format!(
            "{VHH_TO_CDR3}{cdr3}{SHORT_FR4}{LINKER}{KAPPA_CONSTANT}"
        ));
        assert_eq!(segments.cdr3, cdr3);
        assert!(segments.fr4.starts_with("WGQG"), "{}", segments.fr4);
    }

    #[test]
    fn validate_nucleotide_sequence_correctly_rejects_aa_sequence() {
        let cdr3 = "AADPTRGLSYSGQWERPDYAS";

        let Err(error) = validate_nucleotide_sequence(cdr3) else {
            panic!("expected an amino-acid sequence to be rejected");
        };

        let SequenceError::InvalidSequence(message) = error else {
            panic!("expected InvalidSequence, got {error:?}");
        };
        assert!(message.contains("invalid character"), "{message}");
    }
}

/// Nucleotide input
///
/// The reads are synthetic naive human antibodies: IMGT germline V, D and J genes (the alleles the
/// consensus is built from) joined as in V(D)J recombination, with their ends trimmed and a few
/// random non-templated bases (TdT's N nucleotides) added at each join, from codon 1 to the end of
/// FR4. Each comes with its translation. A clean read must number exactly as protein mode numbers
/// that translation, with `query_start`/`query_end` counted in bases.
///
/// The errors go where antibody genes have homopolymers, runs of 4+ identical bases, at the same
/// IMGT positions in most genes: Nanopore, Ion Torrent and 454 miscount them by one base, and today
/// such a frameshift makes immunum confidently wrong. Losing or doubling one base of a run can
/// be repaired exactly, so the read must number like its clean read, ending that many bases
/// earlier or later, at a lower but still accepted confidence. Inside the CDR3 a frameshift
/// may be a real out-of-frame join instead, so there it must not be repaired.
///
/// Every frameshift found is reported in `NumberingResult::frameshifts`: repaired ones, and the
/// ones left as they are inside the CDR3. Its `span` is taken as the bases the shifted codon
/// consumed, 0-based and inclusive like `query_start`/`query_end`: 2 for a lost base, 4 for a
/// doubled one. A homopolymer error can't be placed within its run, so the span only has to
/// touch the run (one base either side allowed, for placement at a codon boundary), and its
/// `position` has to be a codon the run touches.
#[cfg(test)]
mod nt_tests_synthetic {
    use super::*;
    use crate::types::IG_CHAINS;

    // Synthetic naive antibodies, codon 1 through FR4 (bases: 0-based, end exclusive)
    // Heavy chain: IGHV3-23*01, IGHD3-16*01, IGHJ4*01. Runs: GGGGG 21..26 (FR1 8-9), GGGGGG 42..48
    // (FR1 16-17), GGGG 128..132 (FR2 48-49), GGGGG 310..315 (CDR3), GGGG 346..350 (118-119).
    const HEAVY_DNA: &str = "GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAGCTATTAGTGGTAGTGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGTGCGAAAGGGCATGATTACGTTTGGGGGAGTTATGCTTATACCAGGACAATTGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCA";
    const HEAVY_PROTEIN: &str = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKGHDYVWGSYAYTRTIDYWGQGTLVTVSS";
    // Heavy chain: IGHV1-3*01, IGHD5-5*01, IGHJ6*01. Runs: CCCCC 118..123 (FR2 45-46), GGGGG
    // 346..351 (118-119).
    const HEAVY_IGHJ6_DNA: &str = "CAGGTCCAGCTTGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTTTCCTGCAAGGCTTCTGGATACACCTTCACTAGCTATGCTATGCATTGGGTGCGCCAGGCCCCCGGACAAAGGCTTGAGTGGATGGGATGGATCAACGCTGGCAATGGTAACACAAAATATTCACAGAAGTTCCAGGGCAGAGTCACCATTACCAGGGACACATCCGCGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAAGACACGGCTGTGTATTACTGTGCGAGAGACGAAGATACAGCTATGGTTAGCTACTACTACTACTACGGTATGGACGTCTGGGGGCAAGGGACCACGGTCACCGTCTCCTCA";
    const HEAVY_IGHJ6_PROTEIN: &str = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYAMHWVRQAPGQRLEWMGWINAGNGNTKYSQKFQGRVTITRDTSASTAYMELSSLRSEDTAVYYCARDEDTAMVSYYYYYGMDVWGQGTTVTVSS";
    // Kappa chain: IGKV1-12*01, IGKJ2*01. Runs: CCCC 127..131 (FR2 49-50), GGGG 168..172 (FR3
    // 70-71), TTTT 245..249 (FR3 98-99), TTTT 290..294 (CDR3 end and 118), GGGG 299..303 (120-121).
    const KAPPA_DNA: &str = "GACATCCAGATGACCCAGTCTCCATCTTCCGTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGTCGGGCGAGTCAGGGTATTAGCAGCTGGTTAGCCTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTATGCTGCATCCAGTTTGCAAAGTGGGGTCCCATCAAGGTTCAGCGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGCCTGCAGCCTGAAGATTTTGCAACTTACTATTGTCAACAGGCTAACAGTTTCCCTCACACTTTTGGCCAGGGGACCAAGCTGGAGATCAAA";
    const KAPPA_PROTEIN: &str = "DIQMTQSPSSVSASVGDRVTITCRASQGISSWLAWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQANSFPHTFGQGTKLEIK";
    // Lambda chain: IGLV1-40*01, IGLJ2*01. Runs: CCCC 37..41 (FR1 14-15), GGGG 87..91 (CDR1),
    // CCCCC 133..138 (FR2 49-50), GGGG 174..178 (FR3 70-71).
    const LAMBDA_DNA: &str = "CAGTCTGTGCTGACGCAGCCGCCCTCAGTGTCTGGGGCCCCAGGGCAGAGGGTCACCATCTCCTGCACTGGGAGCAGCTCCAACATCGGGGCAGGTTATGATGTACACTGGTACCAGCAGCTTCCAGGAACAGCCCCCAAACTCCTCATCTATGGTAACAGCAATCGGCCCTCAGGGGTCCCTGACCGATTCTCTGGCTCCAAGTCTGGCACCTCAGCCTCCCTGGCCATCACTGGGCTCCAGGCTGAGGATGAGGCTGATTATTACTGCCAGTCCTATGACAGCAGCCTGAGTGGTTCGGTATTCGGCGGAGGGACCAAGCTGACCGTCCTA";
    const LAMBDA_PROTEIN: &str = "QSVLTQPPSVSGAPGQRVTISCTGSSSNIGAGYDVHWYQQLPGTAPKLLIYGNSNRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYCQSYDSSLSGSVFGGGTKLTVL";
    // IGHV3-23*01, IGHD6-6*01, IGHJ4*01 joined out of frame: 353 bases.
    const OUT_OF_FRAME_HEAVY_DNA: &str = "GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAGCTATTAGTGGTAGTGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGTGCGAAAGCAGGTATAGCAGCTCCTTTGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCA";
    // The first 60 bases of IGHM CH1 and of IGKC, in frame after FR4
    const IGHM_CH1_START: &str = "GGGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTGAGAATTCCCCGTCGGAT";
    const IGKC_START: &str = "CGAACTGTGGCTGCACCATCTGTCTTCATCTTCCCGCCATCTGATGAGCAGTTGAAATCT";
    // 360 random bases
    const RANDOM_DNA: &str = "CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCTACATCACTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTCGATGTCAAACCCCGGGGGGAGCTCAGATATCCGATACAGGGATGAAGAAATAACCTCATCCCATTGGTGACGAAAGGTTGTAAGTAGCTGGCCGCCGAGATAGCTGAGCGGCGAACCACTAGAAAAGGTTCAGACCCCGGAGCCCAGCCGTCACGATTGT";

    const FIXTURES: [(&str, &str, &str); 4] = [
        (HEAVY_DNA, HEAVY_PROTEIN, "AKGHDYVWGSYAYTRTIDY"),
        (HEAVY_IGHJ6_DNA, HEAVY_IGHJ6_PROTEIN, "ARDEDTAMVSYYYYYGMDV"),
        (KAPPA_DNA, KAPPA_PROTEIN, "QQANSFPHT"),
        (LAMBDA_DNA, LAMBDA_PROTEIN, "QSYDSSLSGSV"),
    ];

    // Every confidence is reported (minimum 0), so a test sees the value instead of LowConfidence
    fn nt(scheme: Scheme) -> Annotator {
        Annotator::new_nt(IG_CHAINS, scheme, true, Some(0.0)).unwrap()
    }

    fn protein_numbering(protein: &str, scheme: Scheme) -> NumberingResult {
        Annotator::new(IG_CHAINS, scheme, Some(0.0))
            .unwrap()
            .number(protein)
            .unwrap()
    }

    // The standard genetic code, codons in TCAG order; a stop reads as X, as protein mode takes it
    fn translate(dna: &str) -> String {
        const CODE: &[u8; 64] = b"FFLLSSSSYYXXCCXWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
        let index = |b: u8| match b {
            b'T' => 0,
            b'C' => 1,
            b'A' => 2,
            _ => 3,
        };
        dna.as_bytes()
            .chunks_exact(3)
            .map(|codon| {
                // CODE lists the 64 codons in TCAG order: the first base picks a block of 16,
                // the second a block of 4 within it, the third the codon within that. Read as
                // a base-4 number with T=0, C=1, A=2, G=3, the codon is its own index into CODE:
                // TTT = 0 (F), TCA = 0*16 + 1*4 + 2 = 6 (S), GGG = 3*16 + 3*4 + 3 = 63 (G).
                let first = index(codon[0]);
                let second = index(codon[1]);
                let third = index(codon[2]);
                let amino_acid = CODE[16 * first + 4 * second + third];
                char::from(amino_acid)
            })
            .collect()
    }

    fn reverse_complement(dna: &str) -> String {
        dna.bytes()
            .rev()
            .map(|b| match b {
                b'A' => 'T',
                b'T' => 'A',
                b'C' => 'G',
                _ => 'C',
            })
            .collect()
    }

    fn deleted(dna: &str, at: usize, length: usize) -> String {
        format!("{}{}", &dna[..at], &dna[at + length..])
    }

    fn inserted(dna: &str, at: usize, bases: &str) -> String {
        format!("{}{bases}{}", &dna[..at], &dna[at..])
    }

    // `got` numbers like protein mode's `expected`, its domain spanning bases `span` of the
    // forward read
    #[track_caller]
    fn assert_numbered_like(
        got: &NumberingResult,
        expected: &NumberingResult,
        span: (usize, usize),
    ) {
        assert_eq!(got.chain, expected.chain);
        assert_eq!(got.positions, expected.positions);
        assert_eq!(
            (got.query_start, got.query_end),
            span,
            "domain span in bases"
        );
        assert_eq!(nucleotide(got).strand, Strand::Forward);
        assert_one_residue_per_position(got);
    }

    // The nucleotide part of a result, which nucleotide input always gives
    #[track_caller]
    fn nucleotide(got: &NumberingResult) -> &NucleotideNumbering {
        got.nucleotide
            .as_ref()
            .expect("nucleotide input gives a nucleotide part")
    }

    // Nucleotide input always gives its amino acids, one for each numbered position
    #[track_caller]
    fn assert_one_residue_per_position(got: &NumberingResult) {
        let residues = &nucleotide(got).sequence_aa;
        assert_eq!(residues.chars().count(), got.positions.len(), "{residues}");
    }

    // The domain's amino acids are `protein`
    #[track_caller]
    fn assert_residues(got: &NumberingResult, protein: &str) {
        assert_eq!(nucleotide(got).sequence_aa, protein);
    }

    // Read in one frame from start to end: no frameshift reported
    #[track_caller]
    fn assert_in_frame(got: &NumberingResult) {
        let frameshifts = &nucleotide(got).frameshifts;
        assert!(frameshifts.is_empty(), "{frameshifts:?}");
    }

    // The IMGT positions of the clean read's codons that bases `bases` (of the clean read) touch,
    // one base either side included
    fn codon_positions(clean: &NumberingResult, bases: Range<usize>) -> &[Position] {
        let first = bases.start.saturating_sub(1) / 3;
        let last = (bases.end / 3).min(clean.positions.len() - 1);
        &clean.positions[first..=last]
    }

    // `got` is one frameshift of `kind`, `repaired` or not, its span the 2 or 4 bases it consumed
    // and touching bases `near` of the read (one base either side allowed), numbered at one of
    // `positions`
    #[track_caller]
    fn assert_frameshift(
        got: &Frameshift,
        kind: FrameshiftType,
        repaired: bool,
        near: Range<usize>,
        positions: &[Position],
    ) {
        assert_eq!(got.type_frameshift, kind, "{got:?}");
        assert_eq!(got.repaired, repaired, "{got:?}");
        let (start, end) = got.span;
        let consumed = match kind {
            FrameshiftType::Deletion(n) => 3 - usize::from(n),
            FrameshiftType::Insertion(n) => 3 + usize::from(n),
        };
        assert_eq!(end + 1 - start, consumed, "bases consumed: {got:?}");
        assert!(
            start <= near.end && end + 1 >= near.start,
            "{got:?} isn't next to bases {near:?}"
        );
        assert!(
            positions.contains(&got.position),
            "{got:?} isn't numbered at one of {positions:?}"
        );
    }

    // `read` is `clean_dna` with sequencing errors: it must number like the clean protein, its
    // domain still spanning the whole read, and score below the clean read without falling under
    // the default minimum confidence. Returns the read's numbering and the clean one, for the
    // caller to check the frameshifts.
    #[track_caller]
    fn assert_repaired(
        read: &str,
        clean_dna: &str,
        clean_protein: &str,
    ) -> (NumberingResult, NumberingResult) {
        let clean = protein_numbering(clean_protein, Scheme::IMGT);
        assert_eq!(3 * clean_protein.len(), clean_dna.len());
        let got = nt(Scheme::IMGT).number(read).unwrap();
        assert_numbered_like(&got, &clean, (0, read.len() - 1));
        assert!(
            got.confidence < clean.confidence,
            "a repair must cost confidence: {} vs clean {}",
            got.confidence,
            clean.confidence
        );
        assert!(
            got.confidence >= DEFAULT_MIN_CONFIDENCE,
            "{}",
            got.confidence
        );
        (got, clean)
    }

    // One base of the homopolymer `dna[run]` lost, and one doubled: both reads are repaired, each
    // with one repaired frameshift at the run. Every base of a run reads the same when lost or
    // doubled, so the first one stands for all. The base the repair restores (or drops) is the
    // run's own, so the amino acids come back exactly as in the clean protein, and `segment`
    // gives the clean protein's regions in `translated` while its bases still rebuild the read.
    #[track_caller]
    fn assert_run_errors_repaired(dna: &str, protein: &str, run: Range<usize>) {
        let base = &dna[run.start..run.start + 1];
        assert_eq!(
            dna[run.clone()],
            base.repeat(run.len()),
            "the fixture's run"
        );
        assert_ne!(
            &dna[run.end..run.end + 1],
            base,
            "the run ends at {}",
            run.end
        );
        for (read, kind, run_in_read) in [
            (
                deleted(dna, run.start, 1),
                FrameshiftType::Deletion(1),
                run.start..run.end - 1,
            ),
            (
                inserted(dna, run.start, base),
                FrameshiftType::Insertion(1),
                run.start..run.end + 1,
            ),
        ] {
            let (got, clean) = assert_repaired(&read, dna, protein);
            assert_residues(&got, protein);
            let segments = nt(Scheme::IMGT).segment(&read).unwrap();
            assert_segmented(&segments, &read, protein);
            let [frameshift] = nucleotide(&got).frameshifts.as_slice() else {
                panic!(
                    "one frameshift expected: {:?}",
                    nucleotide(&got).frameshifts
                );
            };
            let positions = codon_positions(&clean, run.clone());
            assert_frameshift(frameshift, kind, true, run_in_read, positions);
        }
    }

    // `segment` on nucleotide input: the regions are bases and rebuild `dna` (the domain's
    // strand), and `translated` holds each region's amino acids exactly as protein mode cuts
    // `protein`
    #[track_caller]
    fn assert_segmented(got: &SegmentResult, dna: &str, protein: &str) {
        let bases: String = got.regions().map(|(_, region)| region).concat();
        assert_eq!(bases, dna, "the regions rebuild the read");
        let expected = Annotator::new(IG_CHAINS, Scheme::IMGT, Some(0.0))
            .unwrap()
            .segment(protein)
            .unwrap();
        let aa = got
            .translated
            .as_ref()
            .expect("nucleotide input gives translated regions");
        let regions = [
            ("fr1", &aa.fr1, &expected.fr1),
            ("cdr1", &aa.cdr1, &expected.cdr1),
            ("fr2", &aa.fr2, &expected.fr2),
            ("cdr2", &aa.cdr2, &expected.cdr2),
            ("fr3", &aa.fr3, &expected.fr3),
            ("cdr3", &aa.cdr3, &expected.cdr3),
            ("fr4", &aa.fr4, &expected.fr4),
        ];
        for (name, got, expected) in regions {
            assert_eq!(got, expected, "{name}");
        }
    }

    fn numbers(positions: &[Position], number: u8) -> bool {
        positions.contains(&Position::new(number))
    }

    // The fixtures themselves, checked in protein mode (no nucleotide input needed): each DNA
    // translates to its protein, and protein mode numbers every residue of it with the expected CDR3.
    // This is the reference every other test compares with.
    #[test]
    fn nt_fixtures_are_valid() {
        for (dna, protein, cdr3) in FIXTURES {
            assert_eq!(translate(dna), protein);
            let numbered = protein_numbering(protein, Scheme::IMGT);
            assert_eq!(
                (numbered.query_start, numbered.query_end),
                (0, protein.len() - 1)
            );
            let segments = numbered.segment(protein).unwrap();
            assert_eq!(segments.cdr3, cdr3);
        }
        // Out of frame: 353 bases, not a multiple of three
        assert_ne!(OUT_OF_FRAME_HEAVY_DNA.len() % 3, 0);
    }

    // DNA that holds no antibody isn't numbered at the default minimum confidence: random bases, and
    // a heavy chain read backwards (not complemented), which has its composition but no codons
    #[test]
    fn nt_dna_without_an_antibody_is_rejected() {
        let annotator = Annotator::new_nt(IG_CHAINS, Scheme::IMGT, true, None).unwrap();
        let backwards: String = HEAVY_DNA.chars().rev().collect();
        for read in [RANDOM_DNA, &backwards] {
            let error = annotator.number(read).unwrap_err();
            assert_eq!(error.kind(), "low_confidence", "{error}");
        }
    }

    // A clean read numbers exactly like its translation: same chain and positions, the domain
    // spanning every base, and the same confidence, for heavy, kappa and lambda. The homopolymers
    // in these reads are read as they are: no frameshift is reported.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_clean_reads_number_like_their_protein() {
        for (dna, protein, _) in FIXTURES {
            let expected = protein_numbering(protein, Scheme::IMGT);
            let got = nt(Scheme::IMGT).number(dna).unwrap();
            assert_numbered_like(&got, &expected, (0, dna.len() - 1));
            assert_in_frame(&got);
            assert_residues(&got, protein);
            assert_eq!(
                (got.cons_start, got.cons_end),
                (expected.cons_start, expected.cons_end)
            );
            assert!((got.confidence - expected.confidence).abs() < 1e-6);
        }
    }

    // Every scheme numbers a read as it numbers its translation
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_every_scheme_numbers_like_the_protein() {
        for scheme in [
            Scheme::IMGT,
            Scheme::Kabat,
            Scheme::Chothia,
            Scheme::Martin,
            Scheme::Aho,
        ] {
            for (dna, protein, _) in [FIXTURES[0], FIXTURES[2]] {
                let got = nt(scheme).number(dna).unwrap();
                assert_numbered_like(
                    &got,
                    &protein_numbering(protein, scheme),
                    (0, dna.len() - 1),
                );
                assert_in_frame(&got);
                assert_residues(&got, protein);
            }
        }
    }

    // Bases before the domain are free, which is how the alignment picks the reading frame: after
    // 1, 2, 3 or 7 extra bases the domain numbers the same, starting that many bases later. A
    // frame chosen at the start isn't a frameshift.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_reading_frame_is_found_after_leading_bases() {
        for (dna, protein, _) in [FIXTURES[0], FIXTURES[3]] {
            let expected = protein_numbering(protein, Scheme::IMGT);
            for leader in ["A", "CA", "TCA", "GATTACA"] {
                let read = format!("{leader}{dna}");
                let got = nt(Scheme::IMGT).number(&read).unwrap();
                assert_numbered_like(&got, &expected, (leader.len(), read.len() - 1));
                assert_in_frame(&got);
                assert_residues(&got, protein);
            }
        }
    }

    // A read from the other strand, as primed in the constant region, numbers like the forward read
    // and says it was numbered on the reverse strand. Not the coordinates: which strand they refer
    // to is undecided.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_reverse_strand_numbers_like_the_forward_read() {
        for (dna, protein, _) in FIXTURES {
            let expected = protein_numbering(protein, Scheme::IMGT);
            let got = nt(Scheme::IMGT).number(&reverse_complement(dna)).unwrap();
            assert_eq!(nucleotide(&got).strand, Strand::Reverse);
            assert_eq!(got.chain, expected.chain);
            assert_eq!(got.positions, expected.positions);
            assert_in_frame(&got);
            assert_residues(&got, protein);
            assert!((got.confidence - expected.confidence).abs() < 1e-6);
        }
    }

    // An amplicon read continues from FR4 into the constant region (here the first 60 bases of IGHM
    // CH1 and of IGKC, in frame). The constant region stays out of the domain, as in protein mode.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_constant_region_after_fr4_stays_out_of_the_domain() {
        for read in [
            format!("{HEAVY_DNA}{IGHM_CH1_START}"),
            format!("{KAPPA_DNA}{IGKC_START}"),
        ] {
            let translation = translate(&read);
            let expected = protein_numbering(&translation, Scheme::IMGT);
            let span = (3 * expected.query_start, 3 * expected.query_end + 2);
            let got = nt(Scheme::IMGT).number(&read).unwrap();
            assert_numbered_like(&got, &expected, span);
            assert_in_frame(&got);
            let domain = expected.query_start..=expected.query_end;
            assert_residues(&got, &translation[domain]);
        }
    }

    // Lower case and RNA (U for T) are the same sequence, and give upper-case amino acids
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_lower_case_and_rna_number_like_dna() {
        let expected = protein_numbering(HEAVY_PROTEIN, Scheme::IMGT);
        for read in [HEAVY_DNA.to_ascii_lowercase(), HEAVY_DNA.replace('T', "U")] {
            let got = nt(Scheme::IMGT).number(&read).unwrap();
            assert_numbered_like(&got, &expected, (0, read.len() - 1));
            assert_in_frame(&got);
            assert_residues(&got, HEAVY_PROTEIN);
        }
    }

    // An scFv as DNA, heavy chain, (G4S)3 linker, kappa chain, holds two domains, each numbered as
    // on its own and spanning its own bases
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_scfv_holds_both_domains() {
        let linker = "GGTGGAGGCGGTTCA".repeat(3);
        let scfv = format!("{HEAVY_DNA}{linker}{KAPPA_DNA}");
        let domains = nt(Scheme::IMGT).domains(&scfv).unwrap();
        assert_eq!(domains.len(), 2);
        let kappa_start = HEAVY_DNA.len() + linker.len();
        for (domain, protein, span) in [
            (&domains[0], HEAVY_PROTEIN, (0, HEAVY_DNA.len() - 1)),
            (&domains[1], KAPPA_PROTEIN, (kappa_start, scfv.len() - 1)),
        ] {
            let got = domain.number(Scheme::IMGT).unwrap();
            assert_numbered_like(&got, &protein_numbering(protein, Scheme::IMGT), span);
            assert_in_frame(&got);
            assert_residues(&got, protein);
        }
    }

    // IGHV3 FR1: TCT GGG GGA at 7-9 holds GGGGG over 8-9, in most IGHV3 genes (about 60% of IGHV).
    // Today one G lost there leaves a single numbered residue at confidence 1.000.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_heavy_fr1_8_9_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(HEAVY_DNA, HEAVY_PROTEIN, 21..26);
    }

    // IGHV3-23 FR1: CCT GGG GGG TCC holds a run of six Gs over 16-17 (FR1 16-17 in half of IGHV)
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_heavy_fr1_16_17_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(HEAVY_DNA, HEAVY_PROTEIN, 42..48);
    }

    // IGHV3-23 FR2: GGGG over 48-49. Today one G lost there truncates the domain at 49.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_heavy_fr2_48_49_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(HEAVY_DNA, HEAVY_PROTEIN, 128..132);
    }

    // IGHV1-3 FR2: CCCCC over 45-46 (FR2 45-49 holds a run in about 40% of IGHV)
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_heavy_fr2_45_46_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(HEAVY_IGHJ6_DNA, HEAVY_IGHJ6_PROTEIN, 118..123);
    }

    // Every IGHJ reads TGGGGC over the J-TRP 118 and G 119, so every heavy chain ends its CDR3 on a G
    // run. Today one G lost reads FR4 as WAREPWSPSPQ, one doubled as WGPGNPGHRLL, both confidently.
    // Next to the CDR3, but a sequencing error all the same: repaired. IGHJ4 here.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_heavy_j_anchor_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(HEAVY_DNA, HEAVY_PROTEIN, 346..350);
    }

    // The same run in IGHJ6, five Gs (TGGGGGC)
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_ighj6_anchor_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(HEAVY_IGHJ6_DNA, HEAVY_IGHJ6_PROTEIN, 346..351);
    }

    // IGKV FR2: CCCC over 49-50 (about 40% of IGKV)
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_kappa_fr2_49_50_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(KAPPA_DNA, KAPPA_PROTEIN, 127..131);
    }

    // IGKV FR3: GGGG over 70-71, the most common kappa hotspot (65% of IGKV)
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_kappa_fr3_70_71_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(KAPPA_DNA, KAPPA_PROTEIN, 168..172);
    }

    // IGKV FR3: TTTT over 98-99, just before the Cys 104 (63% of IGKV)
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_kappa_fr3_98_99_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(KAPPA_DNA, KAPPA_PROTEIN, 245..249);
    }

    // IGKJ2 reads ACT TTT GGC: a run of four Ts over the end of the CDR3 and the J-PHE 118. Like the
    // heavy-chain G run, it overlaps the conserved 118: a sequencing error, repaired.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_kappa_j_anchor_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(KAPPA_DNA, KAPPA_PROTEIN, 290..294);
    }

    // IGKJ2 FR4: CAG GGG ACC holds GGGG over 120-121
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_kappa_fr4_120_121_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(KAPPA_DNA, KAPPA_PROTEIN, 299..303);
    }

    // IGLV FR1: CCCC over 14-15 (40% of IGLV)
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_lambda_fr1_14_15_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(LAMBDA_DNA, LAMBDA_PROTEIN, 37..41);
    }

    // IGLV1-40 CDR1: GGGG over 31-35. A CDR has a weaker profile than a framework, but an error in
    // CDR1 or CDR2 is still most likely a sequencing error: repaired.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_lambda_cdr1_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(LAMBDA_DNA, LAMBDA_PROTEIN, 87..91);
    }

    // IGLV FR2: CCCCC over 49-50, the most common lambda hotspot (63-74% of IGLV)
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_lambda_fr2_49_50_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(LAMBDA_DNA, LAMBDA_PROTEIN, 133..138);
    }

    // IGLV FR3: GGGG over 70-71 (about half of IGLV)
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_lambda_fr3_70_71_homopolymer_error_is_repaired() {
        assert_run_errors_repaired(LAMBDA_DNA, LAMBDA_PROTEIN, 174..178);
    }

    // Two errors in one read, a G lost in FR1 8-9 and one in FR2 48-49, 40 codons apart: both
    // repaired to the clean protein, and both reported in read order. In the read the FR2 run
    // starts a base earlier.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_errors_in_two_frameworks_are_both_repaired() {
        let read = deleted(&deleted(HEAVY_DNA, 128, 1), 21, 1);
        let (got, clean) = assert_repaired(&read, HEAVY_DNA, HEAVY_PROTEIN);
        assert_residues(&got, HEAVY_PROTEIN);
        let [fr1, fr2] = nucleotide(&got).frameshifts.as_slice() else {
            panic!(
                "two frameshifts expected: {:?}",
                nucleotide(&got).frameshifts
            );
        };
        let lost = FrameshiftType::Deletion(1);
        assert_frameshift(fr1, lost, true, 21..25, codon_positions(&clean, 21..26));
        assert_frameshift(fr2, lost, true, 127..130, codon_positions(&clean, 128..132));
    }

    // Nanopore can miscount a long run by two: two Gs lost from the six over FR1 16-17. Reported as
    // one event of two bases. Next to the run, a skipped position may compete with a second
    // one-base repair, so only the domain as a whole is asserted: numbered from 1 through 128 with
    // its conserved Cys 23, Cys 104 and Trp 118, and spanning the whole read.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_two_bases_lost_in_a_long_run_keep_the_domain_whole() {
        let read: String = deleted(HEAVY_DNA, 42, 2);
        let got = nt(Scheme::IMGT).number(&read).unwrap();
        assert_eq!(got.chain, Chain::IGH);
        for number in [1, 23, 104, 118, 128] {
            assert!(numbers(&got.positions, number), "{number} isn't numbered");
        }
        assert_eq!((got.query_start, got.query_end), (0, read.len() - 1));
        assert_one_residue_per_position(&got);
        assert!(
            got.confidence >= DEFAULT_MIN_CONFIDENCE,
            "{}",
            got.confidence
        );
        let [frameshift] = nucleotide(&got).frameshifts.as_slice() else {
            panic!(
                "one two-base frameshift expected: {:?}",
                nucleotide(&got).frameshifts
            );
        };
        let clean = protein_numbering(HEAVY_PROTEIN, Scheme::IMGT);
        let positions = codon_positions(&clean, 42..48);
        let lost = FrameshiftType::Deletion(2);
        assert_frameshift(frameshift, lost, true, 42..46, positions);
    }

    // An error outside a homopolymer: the A of GAC (Asp 81, FR3) lost, or a T inserted before it.
    // There's one place the error can be, but the repair has to guess the base it rebuilds the
    // codon with: the amino acids are the clean protein's except at 81 or a neighbouring codon.
    // One repaired frameshift is reported next to base 217, at one of those codons.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_single_base_error_outside_a_homopolymer_is_repaired() {
        assert_eq!(&HEAVY_DNA[216..219], "GAC");
        for (read, kind) in [
            (deleted(HEAVY_DNA, 217, 1), FrameshiftType::Deletion(1)),
            (inserted(HEAVY_DNA, 217, "T"), FrameshiftType::Insertion(1)),
        ] {
            let (got, clean) = assert_repaired(&read, HEAVY_DNA, HEAVY_PROTEIN);
            let residues = nucleotide(&got).sequence_aa.as_str();
            let guessed = 71..=73; // codon 72 (bases 216-218) and its neighbours
            for (i, (got_aa, clean_aa)) in residues.chars().zip(HEAVY_PROTEIN.chars()).enumerate() {
                assert!(
                    got_aa == clean_aa || guessed.contains(&i),
                    "residue {i} is {got_aa}, clean {clean_aa}: {residues}"
                );
            }
            let [frameshift] = nucleotide(&got).frameshifts.as_slice() else {
                panic!(
                    "one frameshift expected: {:?}",
                    nucleotide(&got).frameshifts
                );
            };
            let positions = codon_positions(&clean, 217..218);
            assert_frameshift(frameshift, kind, true, 217..218, positions);
        }
    }

    // A whole codon lost (AGT, Ser 59 in CDR2) keeps the frame: a deletion in the protein, as B-cell
    // mutation makes them, not three frameshifts. The read numbers like its own translation,
    // at its translation's confidence.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_codon_deletion_is_not_a_frameshift() {
        assert_eq!(&HEAVY_DNA[159..162], "AGT");
        let read = deleted(HEAVY_DNA, 159, 3);
        let expected = protein_numbering(&translate(&read), Scheme::IMGT);
        let got = nt(Scheme::IMGT).number(&read).unwrap();
        assert_numbered_like(&got, &expected, (0, read.len() - 1));
        assert_in_frame(&got);
        assert_residues(&got, &translate(&read));
        assert!((got.confidence - expected.confidence).abs() < 1e-6);
    }

    // A stop codon in the CDR3 (TAG for the Ile of ...TIDY) is kept, read as an unknown residue, and
    // not avoided by inventing two frameshifts around it: the read numbers like its
    // translation, X for the stop in its amino acids. (The stop should also be reported, once
    // results can say so.)
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_stop_codon_is_kept_not_avoided() {
        let codon = 3 * HEAVY_PROTEIN.find("IDYWGQG").unwrap();
        let read = format!("{}TAG{}", &HEAVY_DNA[..codon], &HEAVY_DNA[codon + 3..]);
        let expected = protein_numbering(&translate(&read), Scheme::IMGT);
        let got = nt(Scheme::IMGT).number(&read).unwrap();
        assert_numbered_like(&got, &expected, (0, read.len() - 1));
        assert_in_frame(&got);
        assert_residues(&got, &translate(&read));
    }

    // A G lost from the GGGGG inside the CDR3 (from the D gene). In one read this can't be told from
    // an out-of-frame join, so it isn't repaired: the V gene is still numbered (1-104 as in the
    // clean read), but FR4 isn't numbered as if the read were in frame. The frameshift is still
    // reported, as a lost base next to the run, not repaired.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_frameshift_inside_cdr3_is_not_repaired() {
        assert_eq!(&HEAVY_DNA[310..315], "GGGGG");
        let frameshift = assert_not_repaired(&deleted(HEAVY_DNA, 310, 1));
        let clean = protein_numbering(HEAVY_PROTEIN, Scheme::IMGT);
        let positions = codon_positions(&clean, 310..315);
        let lost = FrameshiftType::Deletion(1);
        assert_frameshift(&frameshift, lost, false, 310..314, positions);
    }

    // A real out-of-frame rearrangement: IGHV3-23 and IGHJ4 joined over 353 bases, so FR4 lies in
    // another frame than the V gene. Most B cells carry one beside their productive chain. It must
    // not be repaired into an antibody that doesn't exist, the same as above. Whether it's reported
    // as a lost base or two extra ones is up to the alignment.
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_out_of_frame_rearrangement_is_not_repaired() {
        let frameshift = assert_not_repaired(OUT_OF_FRAME_HEAVY_DNA);
        assert!(
            matches!(
                frameshift.type_frameshift,
                FrameshiftType::Deletion(1) | FrameshiftType::Insertion(2)
            ),
            "353 bases are a base short of the frame: {frameshift:?}"
        );
    }

    // Found as a heavy chain with the V positions and residues of the clean heavy chain (same
    // IGHV3-23), and no J-TRP 118 numbered. The one frameshift reported is left unrepaired inside the CDR3 (IMGT
    // 105-117), and returned for the caller to check.
    #[track_caller]
    fn assert_not_repaired(read: &str) -> Frameshift {
        let clean = protein_numbering(HEAVY_PROTEIN, Scheme::IMGT);
        let got = nt(Scheme::IMGT).number(read).unwrap();
        assert_eq!(got.chain, Chain::IGH);
        assert_eq!(nucleotide(&got).strand, Strand::Forward);
        for position in clean.positions.iter().filter(|p| p.number <= 104) {
            assert!(
                got.positions.contains(position),
                "{position} isn't numbered"
            );
        }
        assert!(!numbers(&got.positions, 118), "FR4 numbered as if repaired");
        assert_one_residue_per_position(&got);
        let v_residues = clean.positions.iter().filter(|p| p.number <= 104).count();
        let residues = nucleotide(&got).sequence_aa.as_str();
        assert!(
            residues.starts_with(&HEAVY_PROTEIN[..v_residues]),
            "{residues}"
        );
        let [frameshift] = nucleotide(&got).frameshifts.as_slice() else {
            panic!(
                "one frameshift expected: {:?}",
                nucleotide(&got).frameshifts
            );
        };
        assert!(!frameshift.repaired, "{frameshift:?}");
        assert!(
            (105..=117).contains(&frameshift.position.number),
            "outside the CDR3: {frameshift:?}"
        );
        frameshift.clone()
    }

    // A clean read segments into its bases, which rebuild it, with each region's amino acids in
    // `translated`, as protein mode segments its translation
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_segment_gives_bases_and_amino_acids() {
        for (dna, protein, cdr3) in FIXTURES {
            let got = nt(Scheme::IMGT).segment(dna).unwrap();
            assert_segmented(&got, dna, protein);
            assert_eq!(got.translated.unwrap().cdr3, cdr3);
            assert_eq!((got.prefix.as_str(), got.postfix.as_str()), ("", ""));
        }
    }

    // Bases around the domain land in the prefix and postfix as bases; they aren't translated
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_segment_keeps_flanking_bases_untranslated() {
        let read = format!("GATTACA{HEAVY_DNA}{IGHM_CH1_START}");
        let got = nt(Scheme::IMGT).segment(&read).unwrap();
        assert_eq!(got.prefix, "GATTACA");
        assert_eq!(got.postfix, IGHM_CH1_START);
        assert_segmented(&got, &read, HEAVY_PROTEIN);
    }

    // A reverse-strand read segments in the domain's orientation: its regions read 5' to 3' as
    // on the forward strand and rebuild the reverse complement of the read
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_segment_of_a_reverse_read_is_in_the_domains_orientation() {
        for (dna, protein, _) in FIXTURES {
            let got = nt(Scheme::IMGT).segment(&reverse_complement(dna)).unwrap();
            assert_segmented(&got, dna, protein);
        }
    }

    // An scFv's domains each segment into their own bases and amino acids; all regions in order
    // rebuild the read, linker in the kappa domain's prefix
    #[test]
    #[ignore = "nucleotide input not implemented yet"]
    fn nt_segment_domains_of_an_scfv() {
        let linker = "GGTGGAGGCGGTTCA".repeat(3);
        let scfv = format!("{HEAVY_DNA}{linker}{KAPPA_DNA}");
        let domains = nt(Scheme::IMGT).segment_domains(&scfv).unwrap();
        let [heavy, kappa] = domains.as_slice() else {
            panic!("two domains expected: {domains:?}");
        };
        assert_segmented(heavy, HEAVY_DNA, HEAVY_PROTEIN);
        assert_eq!(kappa.prefix, linker);
        assert_segmented(kappa, &format!("{linker}{KAPPA_DNA}"), KAPPA_PROTEIN);
    }

    // Protein segmentation has no translated part, and doesn't serialize one
    #[test]
    fn nt_protein_segment_has_no_translated_regions() {
        let got = Annotator::new(IG_CHAINS, Scheme::IMGT, None)
            .unwrap()
            .segment(HEAVY_PROTEIN)
            .unwrap();
        assert!(got.translated.is_none());
        let json = serde_json::to_value(&got).unwrap();
        assert!(json.get("translated").is_none(), "{json}");
    }

    // Protein results serialize as before: no nucleotide part. This keeps existing JSON
    // output and its consumers unchanged.
    #[test]
    fn nt_protein_result_serializes_without_nucleotide_fields() {
        let result = protein_numbering(HEAVY_PROTEIN, Scheme::IMGT);
        let json = serde_json::to_value(&result).unwrap();
        assert!(json.get("nucleotide").is_none(), "{json}");
        let back: NumberingResult = serde_json::from_value(json).unwrap();
        assert!(back.nucleotide.is_none());
    }

    // A nucleotide result keeps its strand, frameshifts and amino acids through JSON
    #[test]
    fn nt_result_with_frameshifts_round_trips_through_json() {
        let mut result = protein_numbering(HEAVY_PROTEIN, Scheme::IMGT);
        result.nucleotide = Some(NucleotideNumbering {
            strand: Strand::Reverse,
            frameshifts: vec![Frameshift {
                span: (21, 22),
                type_frameshift: FrameshiftType::Deletion(1),
                position: Position::new(8),
                repaired: true,
            }],
            sequence_aa: HEAVY_PROTEIN.to_string(),
        });
        let json = serde_json::to_string(&result).unwrap();
        let back: NumberingResult = serde_json::from_str(&json).unwrap();
        let back = nucleotide(&back);
        assert_eq!(back.strand, Strand::Reverse);
        let [frameshift] = back.frameshifts.as_slice() else {
            panic!("{json}");
        };
        assert_eq!(frameshift.span, (21, 22));
        assert_eq!(frameshift.type_frameshift, FrameshiftType::Deletion(1));
        assert_eq!(frameshift.position, Position::new(8));
        assert!(frameshift.repaired);
        assert_eq!(back.sequence_aa, HEAVY_PROTEIN);
    }
}
