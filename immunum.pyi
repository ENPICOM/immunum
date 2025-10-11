"""
Python stub file for immunum - A high-performance library for numbering antibody and T-cell receptor sequences.

This module provides functionality for numbering immunoglobulin and T-cell receptor sequences
using standard immunological numbering schemes (IMGT and KABAT).
"""

from typing import List, Optional, Tuple, Union

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

class Annotator:
    """Main class for sequence annotation and numbering."""

    def __init__(
        self,
        scheme: Scheme = Scheme.IMGT,
        chains: Optional[List[Chain]] = None,
        disable_prefiltering: bool = False,
        threads: Optional[int] = None,
        min_confidence: float = 0.7,
        min_kmer_overlap: Optional[float] = None,
    ) -> None:
        """
        Create a new annotator instance.

        Args:
            scheme: The numbering scheme to use (IMGT or KABAT). Defaults to IMGT.
            chains: List of chains to detect. Defaults to [IGH, IGK, IGL] if not provided.
            disable_prefiltering: Whether to disable k-mer prefiltering. Defaults to False.
            threads: Number of threads for parallel processing. Defaults to CPU cores.
            min_confidence: Minimum alignment confidence threshold. Defaults to 0.7.
            min_kmer_overlap: K-mer overlap threshold for prefiltering. Defaults to 0.2.

        Raises:
            RuntimeError: If annotator creation fails
        """
        ...

    def number_sequences(
        self, sequences: List[Union[str, Tuple[str, str]]], max_chains: int = 2
    ) -> List[List[Tuple[List[str], float, Chain, int, int]]]:
        """
        Number multiple sequences with automatic parallel processing.

        Args:
            sequences: List of sequences. Each can be either:
                - str: Just the sequence string
                - tuple(str, str): (sequence_id, sequence) pair
            max_chains: Maximum number of chains to find per sequence. Defaults to 2.

        Returns:
            List of results for each input sequence. Each result is a list of tuples:
            [(numbers, confidence, chain, start, end), ...]
            where:
            - numbers: List of position strings according to numbering scheme
            - confidence: Alignment confidence score (0.0-1.0)
            - chain: Detected chain type
            - start: First non gap position
            - end: Last non gap position

        Example:
            >>> annotator = Annotator(scheme=Scheme.IMGT, chains=[Chain.IGH])
            >>> sequences = ["QVQLVQSGAEVKKPGASVKVSCKAS...", "DIQMTQSPSSLSASVGDRVTITC..."]
            >>> results = annotator.number_sequences(sequences)
            >>> for i, sequence_chains in enumerate(results):
            ...     print(f"Sequence {i+1}: Found {len(sequence_chains)} chains")
            ...     for numbers, confidence, chain, start, end in sequence_chains:
            ...         print(f"  Chain {chain}: {confidence:.2f} confidence")
        """
        ...
