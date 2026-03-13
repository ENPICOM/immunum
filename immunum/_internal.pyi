__version__: str

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

    def number(self, sequence: str): ...
    def segment(self, sequence: str): ...
