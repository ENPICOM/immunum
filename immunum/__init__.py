from immunum._internal import _Annotator  # noqa: F401
from dataclasses import dataclass


@dataclass(frozen=True)
class SegmenationResult:
    fr1: str
    cdr1: str
    fr2: str
    cdr2: str
    fr3: str
    cdr3: str
    fr4: str
    prefix: str
    postfix: str

    def as_dict(self) -> dict[str, str]:
        return {
            "fr1": self.fr1,
            "cdr1": self.cdr1,
            "fr2": self.fr2,
            "cdr2": self.cdr2,
            "fr3": self.fr3,
            "cdr3": self.cdr3,
            "fr4": self.fr4,
            "prefix": self.prefix,
            "postfix": self.postfix,
        }


@dataclass(frozen=True)
class NumberingResult:
    chain: str
    scheme: str
    confidence: float
    numbering: dict[str, str]


class Annotator:
    """Annotates antibody and T-cell receptor sequences with IMGT or Kabat position numbers.

    Args:
        chains: Chain types to consider during auto-detection. Each entry is a
            case-insensitive string. Accepted values:

            - Antibody heavy chain: ``"IGH"`` / ``"H"`` / ``"heavy"``
            - Antibody kappa chain: ``"IGK"`` / ``"K"`` / ``"kappa"``
            - Antibody lambda chain: ``"IGL"`` / ``"L"`` / ``"lambda"``
            - TCR alpha chain:       ``"TRA"`` / ``"A"`` / ``"alpha"``
            - TCR beta chain:        ``"TRB"`` / ``"B"`` / ``"beta"``
            - TCR gamma chain:       ``"TRG"`` / ``"G"`` / ``"gamma"``
            - TCR delta chain:       ``"TRD"`` / ``"D"`` / ``"delta"``

            Pass all chains you want to consider; the annotator scores each and picks the
            best-matching one. To consider every supported chain pass all seven values.

        scheme: Numbering scheme to use for output positions. Accepted values
            (case-insensitive):

            - ``"IMGT"`` / ``"i"`` — IMGT numbering (recommended; used internally)
            - ``"Kabat"`` / ``"k"`` — Kabat numbering (derived from IMGT)

            Note: Kabat is only supported for antibody chains (IGH, IGK, IGL).

        min_confidence: Minimum alignment confidence threshold in the range ``[0, 1]``.
            Sequences scoring below this value raise a ``ValueError``. Defaults to
            ``0.5``, which filters non-immunoglobulin sequences while retaining all
            validated antibody sequences. Pass ``0.0`` to disable filtering.
    """

    def __init__(
        self,
        chains: list[str],
        scheme: str,
        min_confidence: float | None = None,
    ):
        self._annotator = _Annotator(
            chains=chains, scheme=scheme, min_confidence=min_confidence
        )

    def number(self, sequence: str) -> NumberingResult:
        return NumberingResult(**self._annotator.number(sequence))

    def segment(self, sequence: str) -> SegmenationResult:
        return SegmenationResult(**self._annotator.segment(sequence))
