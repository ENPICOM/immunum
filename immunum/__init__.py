from immunum._internal import (  # noqa: F401
    Error,
    _Annotator,
    _regions_for,
    _scheme_supports_chain,
)
from dataclasses import dataclass, fields
from typing import Optional


@dataclass(frozen=True)
class SegmenationResult:
    """
    Python dataclass containing numbering results. Allows for direct atribute access
    via `results.fr1`, and also for iterating through segmentation results via `as_dict()`:

    ```python
    from immunum import Annotator

    annotator = Annotator(
        chains=["H", "K", "L"],
        scheme="imgt",
    )

    sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"

    result = annotator.segment(sequence)
    assert (
        result.fr1
        == "QVQLVQSGAEVKRPGSSVTVSCKAS"
    )
    assert result.cdr1 == "GGSFSTYA"
    assert result.fr2 == "LSWVRQAPGRGLEWMGG"
    assert result.cdr2 == "VIPLLTIT"
    assert (
        result.fr3
        == "NYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYC"
    )
    assert result.cdr3 == "AREGTTGKPIGAFAH"
    assert result.fr4 == "WGQGTLVTVSS"

    for (
        segment,
        aminoacids,
    ) in result.as_dict().items():
        print(f"{segment}: {aminoacids}")

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

    fr1: Optional[str]
    cdr1: Optional[str]
    fr2: Optional[str]
    cdr2: Optional[str]
    fr3: Optional[str]
    cdr3: Optional[str]
    fr4: Optional[str]
    prefix: Optional[str]
    postfix: Optional[str]
    error: Optional[str]
    """Why the sequence couldn't be segmented, or ``None`` on success."""
    error_kind: Optional[str]
    """What went wrong, as a stable code, or ``None`` on success: ``"invalid_sequence"``
    (too short, too long or not amino acids), ``"low_confidence"`` (no alignment reached
    ``min_confidence``) or, from ``segment_domains`` only, ``"domain_too_short"`` (the best
    alignment is confident but too short to be a domain)."""

    def as_dict(self) -> dict[str, Optional[str]]:
        """Return dict mapping segment names to sequences (excludes the error fields)

        Returns:
            dict[str, str | None]: dict mapping ['fr1', 'fr2', ...] to their aminoacid sequences
        """
        return {
            f.name: getattr(self, f.name)
            for f in fields(self)
            if f.name not in ("error", "error_kind")
        }


@dataclass(frozen=True)
class NumberingResult:
    """Python dataclass containing numbering results. Allows for direct attribute access
    via `result.chain`, `result.numbering`, etc.:

    ```python
    from immunum import Annotator

    annotator = Annotator(
        chains=["H", "K", "L"],
        scheme="imgt",
    )

    sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"

    result = annotator.number(sequence)
    assert result.chain == "H"
    assert result.scheme == "IMGT"
    assert isinstance(
        result.confidence, float
    )
    assert result.numbering["1"] == "Q"

    for (
        position,
        amino_acid,
    ) in result.numbering.items():
        print(f"{position}: {amino_acid}")

    # 1: Q
    # 2: V
    # 3: Q
    # ...
    ```
    """

    chain: Optional[str]
    scheme: Optional[str]
    confidence: Optional[float]
    numbering: Optional[dict[str, str]]
    query_start: Optional[int]
    query_end: Optional[int]
    error: Optional[str]
    """Why the sequence couldn't be numbered, or ``None`` on success."""
    error_kind: Optional[str]
    """What went wrong, as a stable code, or ``None`` on success: ``"invalid_sequence"``
    (too short, too long or not amino acids), ``"low_confidence"`` (no alignment reached
    ``min_confidence``) or, from ``number_domains`` only, ``"domain_too_short"`` (the best
    alignment is confident but too short to be a domain)."""


class Annotator:
    """Annotates antibody and T-cell receptor sequences with scheme-specific position numbers.

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

            A group of chains is accepted too: ``"ig"`` (IGH, IGK, IGL), ``"tcr"``
            (TRA, TRB, TRG, TRD) or ``"all"``. Pass all chains you want to consider;
            the annotator scores each and picks the best-matching one.

        scheme: Numbering scheme to use for output positions. Accepted values
            (case-insensitive):

            - ``"IMGT"`` / ``"i"`` — IMGT numbering (recommended; used internally)
            - ``"Kabat"`` / ``"k"`` — Kabat numbering (derived from IMGT)
            - ``"Chothia"`` / ``"c"`` — Chothia numbering (derived from IMGT)
            - ``"Martin"`` / ``"m"`` — Martin / extended Chothia numbering (derived from IMGT)
            - ``"Aho"`` / ``"a"`` — AHo numbering (derived from IMGT)

            Note: only IMGT supports TCR chains. Kabat, Chothia, Martin and AHo are
            restricted to antibody chains (IGH, IGK, IGL).

        min_confidence: Minimum alignment confidence threshold in the range ``[0, 1]``.
            Sequences scoring below this value get a result with ``error`` set and
            ``error_kind`` ``"low_confidence"``. Defaults to ``0.5``, which filters
            non-immunoglobulin sequences while retaining all validated antibody
            sequences. Pass ``0.0`` to disable filtering.

    Raises:
        immunum.Error: A ``ValueError`` raised when the arguments are wrong. Its
            ``kind`` attribute names what went wrong: ``"invalid_chain"`` (an unknown
            chain, or no chains), ``"invalid_scheme"``, ``"unsupported_chain"`` (a
            non-IMGT scheme for TCR chains) or ``"invalid_min_confidence"``.
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
            scheme: Numbering scheme — ``"imgt"``, ``"kabat"``, ``"chothia"``,
                ``"martin"`` or ``"aho"``. See class docstring for aliases.
            min_confidence: Reject sequences with alignment confidence below this
                threshold. Defaults to ``0.5``; pass ``0.0`` to disable.

        Raises:
            immunum.Error: If any chain or scheme value is unrecognised, if a
                non-IMGT scheme is requested for TCR chains, or if
                ``min_confidence`` is outside ``[0, 1]``. See the class docstring
                for its ``kind`` values.
        """
        self._annotator = _Annotator(
            chains=chains, scheme=scheme, min_confidence=min_confidence
        )

    def number(self, sequence: str) -> NumberingResult:
        """Assign scheme-specific position numbers to every residue in a sequence.

        Args:
            sequence: Amino-acid sequence string (single-letter codes).

        Returns:
            A `NumberingResult` with the detected chain, scheme, confidence score,
            and a ``{position: residue}`` numbering dict. On failure, ``error`` and
            ``error_kind`` are set and all other fields are ``None``.
        """
        return NumberingResult(**self._annotator.number(sequence))

    def number_domains(self, sequence: str) -> list[NumberingResult]:
        """Number every variable domain in a sequence, such as both domains of an scFv.

        ```python
        from immunum import Annotator

        heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
        kappa = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
        scfv = heavy + "GGGGSGGGGSGGGGS" + kappa

        annotator = Annotator(
            chains=["ig"], scheme="imgt"
        )
        domains = annotator.number_domains(scfv)
        assert [d.chain for d in domains] == [
            "H",
            "K",
        ]
        assert (
            domains[1].query_start
            == len(heavy) + 15
        )
        assert domains[1].numbering["1"] == "D"
        ```

        Args:
            sequence: Amino-acid sequence string (single-letter codes).

        Returns:
            One `NumberingResult` per domain, ordered by position, each what `number`
            returns for that domain; never empty. A sequence without a domain gives a
            single result with ``error`` and ``error_kind`` set: ``"low_confidence"``
            when no alignment reaches ``min_confidence``, as `number` reports it, or
            ``"domain_too_short"`` when the best alignment is confident but shorter
            than a domain must be. So does an invalid sequence (too short, too long or
            not amino acids), with ``"invalid_sequence"``.

        A domain that lacks its first IMGT positions (a light chain starting at position 2,
        say) and directly follows other residues, such as a linker, can have the residue just
        before it numbered as its first position, without a known germline or source, numbering
        can't tell a linker residue from the domain's own first residue.
        """
        return [NumberingResult(**d) for d in self._annotator.number_domains(sequence)]

    def segment(self, sequence: str) -> SegmenationResult:
        """Split a sequence into FR/CDR regions.

        Args:
            sequence: Amino-acid sequence string (single-letter codes).

        Returns:
            A `SegmenationResult` with ``fr1``–``fr4``, ``cdr1``–``cdr3``,
            and the residues before and after the domain as ``prefix``/``postfix``,
            so the regions in order rebuild the sequence. On failure, ``error`` and
            ``error_kind`` are set and all region fields are ``None``.
        """
        return SegmenationResult(**self._annotator.segment(sequence))

    def segment_domains(self, sequence: str) -> list[SegmenationResult]:
        """Split every variable domain in a sequence into FR/CDR regions.

        Every residue lands in exactly one domain's regions: a domain's ``prefix`` holds
        the residues since the previous domain (or the start of the sequence), and only
        the last domain has the residues after it as its ``postfix``. All domains'
        regions in order rebuild the sequence.

        ```python
        from immunum import Annotator

        heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
        kappa = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
        linker = "GGGGSGGGGSGGGGS"

        annotator = Annotator(
            chains=["ig"], scheme="imgt"
        )
        heavy_regions, kappa_regions = (
            annotator.segment_domains(
                heavy + linker + kappa
            )
        )
        assert (
            heavy_regions.cdr3
            == "AREGTTGKPIGAFAH"
        )
        assert kappa_regions.prefix == linker
        assert kappa_regions.cdr3 == "QQHYTTPPT"
        ```

        Args:
            sequence: Amino-acid sequence string (single-letter codes).

        Returns:
            One `SegmenationResult` per domain, ordered by position; never empty. A
            sequence without a domain or an invalid one gives a single result with
            ``error`` and ``error_kind`` set, as `number_domains` describes.
        """
        return [
            SegmenationResult(**d) for d in self._annotator.segment_domains(sequence)
        ]


def regions_for(scheme: str, chain: str) -> dict[str, tuple[int, int]]:
    """Look up the FR/CDR region boundaries a scheme uses for a chain.

    ```python
    from immunum import regions_for

    kabat_heavy = regions_for("kabat", "H")
    assert kabat_heavy["cdr1"] == (31, 35)
    assert kabat_heavy["fr4"] == (103, 113)
    ```

    Args:
        scheme: Numbering scheme — ``"imgt"``, ``"kabat"``, ``"chothia"``,
            ``"martin"`` or ``"aho"``. See `Annotator` for aliases.
        chain: Chain type — ``"IGH"`` / ``"H"`` / ``"heavy"`` and so on. See
            `Annotator` for accepted values.

    Returns:
        dict[str, tuple[int, int]]: ``{region: (start, end)}`` with both bounds
            inclusive, in N- to C-terminal order (``fr1``, ``cdr1``, … ``fr4``).
            IMGT and AHo number every chain alike; Kabat, Chothia and Martin
            place their CDRs differently on heavy and light chains.

    Raises:
        immunum.Error: If the scheme or chain is unrecognised (``kind``
            ``"invalid_scheme"`` or ``"invalid_chain"``), or if the scheme has no
            rules for the chain (``"unsupported_chain"``) — only IMGT covers TCR
            chains, so there is no Kabat, Chothia, Martin or AHo table to return for one.
    """
    return _regions_for(scheme=scheme, chain=chain)


def scheme_supports_chain(scheme: str, chain: str) -> bool:
    """Tell whether a scheme numbers a chain, before building an `Annotator` for them.

    ```python
    from immunum import scheme_supports_chain

    assert scheme_supports_chain("kabat", "H")
    assert not scheme_supports_chain("kabat", "B")  # only IMGT numbers TCR chains
    ```

    Args:
        scheme: Numbering scheme, as for `Annotator`.
        chain: A single chain, as for `Annotator`.

    Returns:
        bool: ``True`` when `Annotator` accepts the pair. IMGT numbers every chain;
            Kabat, Chothia, Martin and AHo number antibody chains (IGH, IGK, IGL) only.

    Raises:
        immunum.Error: If the scheme or chain is unrecognised (``kind``
            ``"invalid_scheme"`` or ``"invalid_chain"``).
    """
    return _scheme_supports_chain(scheme=scheme, chain=chain)
