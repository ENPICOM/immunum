__version__: str

class _Annotator:
    def __init__(
        self,
        chains: list[str],
        scheme: str,
        min_confidence: float | None = None,
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
