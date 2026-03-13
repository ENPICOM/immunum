from typing import TypedDict

__version__: str

class NumberingResult(TypedDict):
    """Result of numbering a single antibody or TCR sequence.

    Attributes:
        chain: Detected chain type (e.g. ``"H"``, ``"K"``, ``"L"``, ``"A"``, ``"B"``).
        scheme: Numbering scheme applied (e.g. ``"imgt"``, ``"kabat"``).
        positions: Position labels for each residue
            (e.g. ``["1", "2", ..., "111A", "111B", ...]``).
        residues: Amino acid at each position, one character per entry.
        confidence: Alignment confidence score — higher is better.
    """

    chain: str
    scheme: str
    positions: list[str]
    residues: list[str]
    confidence: float

class _Annotator:
    """Annotates antibody and TCR sequences with IMGT or Kabat position numbers.

    Wraps the Rust alignment engine. Create one instance per chain/scheme
    combination and reuse it across many sequences for best performance.
    """

    def __init__(
        self,
        chains: list[str],
        scheme: str,
    ) -> None:
        """Create an _Annotator for the given chains and numbering scheme.

        Args:
            chains: Chain types to consider during auto-detection.
                Antibody: ``"H"`` (heavy), ``"K"`` (kappa), ``"L"`` (lambda).
                TCR: ``"A"``, ``"B"``, ``"G"``, ``"D"``.
                Pass multiple to enable auto-detection across them.
            scheme: Numbering scheme. Supported: ``"imgt"``, ``"kabat"``.
        """
        ...

    def number(self, sequence: str) -> NumberingResult:
        """Number a single amino acid sequence.

        Args:
            sequence: Single-letter amino acid sequence string.

        Returns:
            A [NumberingResult][immunum._internal.NumberingResult] with
            positions, residues, chain, scheme, and confidence.

        Raises:
            ValueError: If the sequence cannot be aligned to any of the
                specified chains.
        """
        ...
