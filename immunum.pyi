"""
Python stub file for immunum - A high-performance library for numbering antibody and T-cell receptor sequences.

This module provides functionality for numbering immunoglobulin and T-cell receptor sequences 
using standard immunological numbering schemes (IMGT and KABAT).
"""

from typing import Dict, List, Optional, Tuple, Union

class Scheme:
    """Numbering schemes for immunoglobulin sequences."""
    IMGT: 'Scheme'
    """IMGT numbering scheme"""
    KABAT: 'Scheme'
    """Kabat numbering scheme"""

class Chain:
    """Immunoglobulin and T-cell receptor chain types."""
    IGH: 'Chain'
    """Heavy chain"""
    IGK: 'Chain'
    """Kappa light chain"""
    IGL: 'Chain'
    """Lambda light chain"""
    TRA: 'Chain'
    """T-cell receptor Alpha chain"""
    TRB: 'Chain'
    """T-cell receptor Beta chain"""
    TRG: 'Chain'
    """T-cell receptor Gamma chain"""
    TRD: 'Chain'
    """T-cell receptor Delta chain"""

class ScoringParams:
    """Scoring parameters for sequence alignment and numbering."""
    
    def __init__(
        self,
        gap_pen_cp: Optional[float] = None,
        gap_pen_fr: Optional[float] = None,
        gap_pen_ip: Optional[float] = None,
        gap_pen_op: Optional[float] = None,
        gap_pen_cdr: Optional[float] = None,
        gap_pen_other: Optional[float] = None,
        cdr_increase: Optional[float] = None,
        pen_leap_insertion_point_imgt: Optional[float] = None,
        pen_leap_insertion_point_kabat: Optional[float] = None,
    ) -> None:
        """
        Create scoring parameters with optional custom values.
        
        Args:
            gap_pen_cp: Gap penalty for conserved positions
            gap_pen_fr: Gap penalty for framework regions
            gap_pen_ip: Gap penalty for insertion points
            gap_pen_op: Gap penalty for opening positions
            gap_pen_cdr: Gap penalty for CDR regions
            gap_pen_other: Gap penalty for other regions
            cdr_increase: CDR increase factor
            pen_leap_insertion_point_imgt: Penalty for leap insertion points in IMGT
            pen_leap_insertion_point_kabat: Penalty for leap insertion points in KABAT
        """
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
    
    def summary(self) -> str:
        """Get a summary string of the annotation result."""
        ...
    
    def get_region_sequence(self, region_name: str) -> Optional[str]:
        """
        Get the sequence for a specific region.
        
        Args:
            region_name: Name of the region (e.g., 'CDR1', 'FR1', etc.)
            
        Returns:
            The sequence for the specified region, or None if not found.
        """
        ...
    
    def get_cdr_sequences(self) -> Dict[str, str]:
        """
        Get all CDR sequences.
        
        Returns:
            Dictionary mapping CDR names to their sequences.
        """
        ...
    
    def get_framework_sequences(self) -> Dict[str, str]:
        """
        Get all framework sequences.
        
        Returns:
            Dictionary mapping framework names to their sequences.
        """
        ...
    
    def is_high_confidence(self, threshold: float) -> bool:
        """
        Check if the annotation result meets a confidence threshold.
        
        Args:
            threshold: Minimum identity score required
            
        Returns:
            True if identity score >= threshold
        """
        ...

class Annotator:
    """Main class for sequence annotation and numbering."""
    
    def __init__(
        self,
        scheme: Scheme,
        chains: Union[Chain, List[Chain]],
        scoring_params: Optional[ScoringParams] = None,
        use_prefiltering: Optional[bool] = None,
    ) -> None:
        """
        Create a new annotator instance.
        
        Args:
            scheme: The numbering scheme to use (IMGT or KABAT)
            chains: Single chain or list of chains to annotate
            scoring_params: Optional custom scoring parameters
            use_prefiltering: Whether to use prefiltering for performance
            
        Raises:
            RuntimeError: If annotator creation fails
        """
        ...
    
    def number_sequence(self, sequence: str) -> AnnotationResult:
        """
        Number a single sequence.
        
        Args:
            sequence: The amino acid sequence to number
            
        Returns:
            Annotation result containing numbering and analysis
            
        Raises:
            RuntimeError: If sequence numbering fails
        """
        ...
    
    def number_sequences(
        self, 
        sequences: List[str], 
        parallel: bool = False
    ) -> List[AnnotationResult]:
        """
        Number multiple sequences.
        
        Args:
            sequences: List of amino acid sequences to number
            parallel: Whether to use parallel processing
            
        Returns:
            List of annotation results
            
        Raises:
            RuntimeError: If any sequence numbering fails
        """
        ...
    
    def number_file(self, file_path: str) -> List[Tuple[str, AnnotationResult]]:
        """
        Number sequences from a FASTA or FASTQ file.
        
        Args:
            file_path: Path to the input file (supports .fasta, .fastq, .gz)
            
        Returns:
            List of tuples containing (sequence_name, annotation_result)
            
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