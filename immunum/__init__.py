from immunum._internal import _Annotator  # noqa: F401
from dataclasses import dataclass

_CHAIN_ALIASES: dict[str, str] = {
    "igh": "IGH",
    "h": "IGH",
    "heavy": "IGH",
    "igk": "IGK",
    "k": "IGK",
    "kappa": "IGK",
    "igl": "IGL",
    "l": "IGL",
    "lambda": "IGL",
    "tra": "TRA",
    "a": "TRA",
    "alpha": "TRA",
    "trb": "TRB",
    "b": "TRB",
    "beta": "TRB",
    "trg": "TRG",
    "g": "TRG",
    "gamma": "TRG",
    "trd": "TRD",
    "d": "TRD",
    "delta": "TRD",
}

_SCHEME_ALIASES: dict[str, str] = {
    "imgt": "IMGT",
    "i": "IMGT",
    "kabat": "Kabat",
    "k": "Kabat",
}


def _normalize_chains(chains: list[str]) -> list[str]:
    result = []
    for chain in chains:
        normalized = _CHAIN_ALIASES.get(chain.lower())
        if normalized is None:
            valid = sorted(set(_CHAIN_ALIASES.values()))
            raise ValueError(f"Unknown chain {chain!r}. Valid chains: {valid}")
        result.append(normalized)
    return result


def _normalize_scheme(scheme: str) -> str:
    normalized = _SCHEME_ALIASES.get(scheme.lower())
    if normalized is None:
        valid = sorted(set(_SCHEME_ALIASES.values()))
        raise ValueError(f"Unknown scheme {scheme!r}. Valid schemes: {valid}")
    return normalized


@dataclass(frozen=True)
class SegmenationResult:
    """
    Python dataclass containing numbering results. Allows for direct atribute access
    via `results.fr1`, and also for iterating through segmentation results via `as_dict()`:

    ```python
    from immunum import Annotator

    annotator = Annotator(chains=["H", "K", "L"], scheme="imgt")

    sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"

    result = annotator.segment(sequence)
    assert result.fr1 == 'QVQLVQSGAEVKRPGSSVTVSCKAS'
    assert result.cdr1 == 'GGSFSTYA'
    assert result.fr2 == 'LSWVRQAPGRGLEWMGG'
    assert result.cdr2 == 'VIPLLTIT'
    assert result.fr3 == 'NYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYC'
    assert result.cdr3 == 'AREGTTGKPIGAFAH'
    assert result.fr4 == 'WGQGTLVTVSS'

    for segment, aminoacids in result.as_dict().items():
        print(f'{segment}: {aminoacids}')

    # fr1: QVQLVQSGAEVKRPGSSVTVSCKAS
    # cdr1: GGSFSTYA
    # fr2: LSWVRQAPGRGLEWMGG
    # cdr2: VIPLLTIT
    # fr3: NYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYC
    # cdr3: AREGTTGKPIGAFAH
    # fr4: WGQGTLVTVSS
    # prefix:
    # postfix:
    ```
    """

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
        """Return dict mapping segment names to sequences

        Returns:
            dict[str, str]: dict mapping ['fr1', 'fr2', ...] to their aminoacid sequences
        """
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
    """Python dataclass containing numbering results. Allows for direct attribute access
    via `result.chain`, `result.numbering`, etc.:

    ```python
    from immunum import Annotator

    annotator = Annotator(chains=["H", "K", "L"], scheme="imgt")

    sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"

    result = annotator.number(sequence)
    assert result.chain == 'H'
    assert result.scheme == 'IMGT'
    assert isinstance(result.confidence, float)
    assert result.numbering['1'] == 'Q'

    for position, amino_acid in result.numbering.items():
        print(f'{position}: {amino_acid}')

    # 1: Q
    # 2: V
    # 3: Q
    # ...
    ```
    """

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
        """Create an Annotator.

        Args:
            chains: Chain types to consider. See class docstring for accepted values.
            scheme: Numbering scheme — ``"imgt"`` (default) or ``"kabat"``.
            min_confidence: Reject sequences with alignment confidence below this
                threshold. Defaults to ``0.5``; pass ``0.0`` to disable.

        Raises:
            ValueError: If any chain or scheme value is unrecognised, if Kabat is
                requested for TCR chains, or if ``min_confidence`` is outside ``[0, 1]``.
        """
        if min_confidence is not None and not (0 <= min_confidence <= 1.0):
            raise ValueError(
                f"min_confidence should be in [0, 1], got {min_confidence=}"
            )
        self._annotator = _Annotator(
            chains=_normalize_chains(chains),
            scheme=_normalize_scheme(scheme),
            min_confidence=min_confidence,
        )

    def number(self, sequence: str) -> NumberingResult:
        """Assign IMGT or Kabat position numbers to every residue in a sequence.

        Args:
            sequence: Amino-acid sequence string (single-letter codes).

        Returns:
            A `NumberingResult` with the detected chain, scheme, confidence score,
            and a ``{position: residue}`` numbering dict.

        Raises:
            ValueError: If the sequence is empty or scores below ``min_confidence``.
        """
        return NumberingResult(**self._annotator.number(sequence))

    def segment(self, sequence: str) -> SegmenationResult:
        """Split a sequence into FR/CDR regions.

        Args:
            sequence: Amino-acid sequence string (single-letter codes).

        Returns:
            A `SegmenationResult` with ``fr1``–``fr4``, ``cdr1``–``cdr3``,
            and any unaligned ``prefix``/``postfix`` residues.

        Raises:
            ValueError: If the sequence is empty or scores below ``min_confidence``.
        """
        return SegmenationResult(**self._annotator.segment(sequence))
