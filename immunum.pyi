"""
Python stub file for immunum - A high-performance library for numbering antibody and T-cell receptor sequences.

This module provides functionality for numbering immunoglobulin and T-cell receptor sequences
using standard immunological numbering schemes (IMGT and KABAT).
"""

from typing import Dict, List, Optional, Tuple, Union

class Scheme:
    """Numbering schemes for immunoglobulin sequences."""

    IMGT: "Scheme"
    """IMGT numbering scheme"""
    KABAT: "Scheme"
    """Kabat numbering scheme"""

class Chain:
    """Immunoglobulin and T-cell receptor chain types."""

    IGH: "Chain"
    """Heavy chain"""
    IGK: "Chain"
    """Kappa light chain"""
    IGL: "Chain"
    """Lambda light chain"""
    TRA: "Chain"
    """T-cell receptor Alpha chain"""
    TRB: "Chain"
    """T-cell receptor Beta chain"""
    TRG: "Chain"
    """T-cell receptor Gamma chain"""
    TRD: "Chain"
    """T-cell receptor Delta chain"""

class ScoringParams:
    """Scoring parameters for sequence alignment and numbering."""

    def __init__(self) -> None:
        """Create scoring parameters with default values."""
        ...

    @property
    def gap_pen_cp(self) -> float:
        """Gap penalty for conserved positions."""
        ...

    @gap_pen_cp.setter
    def gap_pen_cp(self, value: float) -> None:
        """Set gap penalty for conserved positions."""
        ...

    @property
    def gap_pen_fr(self) -> float:
        """Gap penalty for framework regions."""
        ...

    @gap_pen_fr.setter
    def gap_pen_fr(self, value: float) -> None:
        """Set gap penalty for framework regions."""
        ...

    @property
    def gap_pen_ip(self) -> float:
        """Gap penalty for insertion points."""
        ...

    @gap_pen_ip.setter
    def gap_pen_ip(self, value: float) -> None:
        """Set gap penalty for insertion points."""
        ...

    @property
    def gap_pen_op(self) -> float:
        """Gap penalty for opening positions."""
        ...

    @gap_pen_op.setter
    def gap_pen_op(self, value: float) -> None:
        """Set gap penalty for opening positions."""
        ...

    @property
    def gap_pen_cdr(self) -> float:
        """Gap penalty for CDR regions."""
        ...

    @gap_pen_cdr.setter
    def gap_pen_cdr(self, value: float) -> None:
        """Set gap penalty for CDR regions."""
        ...

    @property
    def gap_pen_other(self) -> float:
        """Gap penalty for other regions."""
        ...

    @gap_pen_other.setter
    def gap_pen_other(self, value: float) -> None:
        """Set gap penalty for other regions."""
        ...

    @property
    def cdr_increase(self) -> float:
        """CDR increase factor."""
        ...

    @cdr_increase.setter
    def cdr_increase(self, value: float) -> None:
        """Set CDR increase factor."""
        ...

    @property
    def pen_leap_insertion_point_imgt(self) -> float:
        """Penalty for leap insertion points in IMGT."""
        ...

    @pen_leap_insertion_point_imgt.setter
    def pen_leap_insertion_point_imgt(self, value: float) -> None:
        """Set penalty for leap insertion points in IMGT."""
        ...

    @property
    def pen_leap_insertion_point_kabat(self) -> float:
        """Penalty for leap insertion points in KABAT."""
        ...

    @pen_leap_insertion_point_kabat.setter
    def pen_leap_insertion_point_kabat(self, value: float) -> None:
        """Set penalty for leap insertion points in KABAT."""
        ...

class AnnotationResult:
    """Result of sequence annotation and numbering."""

    @property
    def sequence(self) -> str:
        """The original input sequence."""
        ...

    @property
    def numbers(self) -> List[str]:
        """List of position numbers corresponding to each residue."""
        ...

    @property
    def scheme(self) -> Scheme:
        """The numbering scheme used."""
        ...

    @property
    def chain(self) -> Chain:
        """The chain type identified."""
        ...

    @property
    def identity(self) -> float:
        """Identity score of the alignment."""
        ...

    @property
    def regions(self) -> Dict[str, Tuple[int, int]]:
        """Dictionary mapping region names to (start, end) positions."""
        ...

    @property
    def start(self) -> int:
        """Start position of the numbering."""
        ...

    @property
    def end(self) -> int:
        """End position of the numbering."""
        ...

    def get_region_sequence(self, region_name: str) -> Optional[str]:
        """
        Get the sequence for a specific region.

        Args:
            region_name: Name of the region (e.g., 'cdr1', 'cdr2', 'cdr3', 'fr1', 'fr2', 'fr3', 'fr4')

        Returns:
            The sequence for the specified region, or None if not found.
        """
        ...

class Annotator:
    """Main class for sequence annotation and numbering."""

    def __init__(
        self,
        scheme: Scheme,
        chains: Union[Chain, List[Chain]],
        scoring_params: Optional[ScoringParams] = None,
        disable_prefiltering: Optional[bool] = None,
    ) -> None:
        """
        Create a new annotator instance.

        Args:
            scheme: The numbering scheme to use (IMGT or KABAT)
            chains: Single chain or list of chains to annotate
            scoring_params: Optional custom scoring parameters
            disable_prefiltering: Whether to disable prefiltering optimization (default: False, prefiltering enabled)

        Raises:
            RuntimeError: If annotator creation fails
        """
        ...

    def number_sequence(self, sequence: str, *, all_chains: bool = False) -> List[AnnotationResult]:
        """
        Number a single sequence.

        Args:
            sequence: The amino acid sequence to number
            all_chains: Whether to return results for all chains

        Returns:
            List of annotation results for the sequence

        Raises:
            RuntimeError: If sequence numbering fails
        """
        ...

    def number_sequences(
        self, sequences: List[str], *, all_chains: bool = False, parallel: bool = False
    ) -> List[List[AnnotationResult]]:
        """
        Number multiple sequences.

        Args:
            sequences: List of amino acid sequences to number
            all_chains: Whether to return results for all chains
            parallel: Whether to use parallel processing

        Returns:
            List of lists of annotation results per input sequence

        Raises:
            RuntimeError: If any sequence numbering fails
        """
        ...

    def number_file(
        self, file_path: str, *, all_chains: bool = False, parallel: bool = False
    ) -> List[Tuple[str, AnnotationResult]]:
        """
        Number sequences from a FASTA or FASTQ file.

        Args:
            file_path: Path to the input file (supports .fasta, .fastq, .gz)
            parallel: Whether to use parallel processing for better performance

        Returns:
            Flattened list of tuples containing (sequence_name, annotation_result)

        Raises:
            RuntimeError: If file processing fails
        """
        ...

def default_scoring_params() -> ScoringParams:
    """
    Get default scoring parameters.

    Returns:
        Default ScoringParams instance
    """
    ...
