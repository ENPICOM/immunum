from immunum._internal import (  # noqa: F401
    Error,
    __version__,
    _Annotator,
    _regions_for,
    _scheme_supports_chain,
)
from dataclasses import dataclass, fields
from typing import Optional


@dataclass(frozen=True)
class SegmentationResult:
    """
    FR/CDR regions of a sequence, plus the residues before and after the domain. Read a
    region as an attribute, such as `result.fr1`, or all of them in sequence order with
    `as_dict()`:

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
    assert (
        "".join(result.as_dict().values())
        == sequence
    )

    for (
        segment,
        aminoacids,
    ) in result.as_dict().items():
        print(f"{segment}: {aminoacids}")

    # prefix:
    # fr1: QVQLVQSGAEVKRPGSSVTVSCKAS
    # cdr1: GGSFSTYA
    # fr2: LSWVRQAPGRGLEWMGG
    # cdr2: VIPLLTIT
    # fr3: NYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYC
    # cdr3: AREGTTGKPIGAFAH
    # fr4: WGQGTLVTVSS
    # postfix:
    ```
    """

    # The regions are declared in sequence order, which `as_dict` keeps, as every interface does.
    prefix: Optional[str]
    """Residues before FR1, such as a leader or the linker after a previous domain; often
    empty. ``None`` on failure."""
    fr1: Optional[str]
    """Framework region 1. ``None`` on failure."""
    cdr1: Optional[str]
    """CDR1. ``None`` on failure."""
    fr2: Optional[str]
    """Framework region 2. ``None`` on failure."""
    cdr2: Optional[str]
    """CDR2. ``None`` on failure."""
    fr3: Optional[str]
    """Framework region 3. ``None`` on failure."""
    cdr3: Optional[str]
    """CDR3. ``None`` on failure."""
    fr4: Optional[str]
    """Framework region 4. ``None`` on failure."""
    postfix: Optional[str]
    """Residues after FR4, such as a constant region; often empty. ``None`` on failure."""
    error: Optional[str]
    """Why the sequence couldn't be segmented, or ``None`` on success."""
    error_kind: Optional[str]
    """What went wrong, as a stable code, or ``None`` on success: ``"invalid_sequence"``
    (shorter than 30 or longer than 10,000 residues, or holding a character that isn't a
    letter), ``"low_confidence"`` (no alignment reached ``min_confidence``) or, from
    ``segment_domains`` only, ``"domain_too_short"`` (the best alignment is confident but
    too short to be a domain)."""

    def as_dict(self) -> dict[str, Optional[str]]:
        """Map each region name to its residues, in sequence order.

        Returns:
            ``prefix``, ``fr1`` … ``fr4`` and ``postfix`` with their residues; the error
            fields are left out. On success, the values joined in order are the sequence.
        """
        return {
            f.name: getattr(self, f.name)
            for f in fields(self)
            if f.name not in ("error", "error_kind")
        }


@dataclass(frozen=True)
class NumberingResult:
    """The numbering of a sequence. Read a field as an attribute, such as `result.chain`
    or `result.numbering`:

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
    """Detected chain: ``"H"``, ``"K"``, ``"L"``, ``"A"``, ``"B"``, ``"G"`` or ``"D"``.
    ``None`` on failure."""
    scheme: Optional[str]
    """Scheme the positions are in: ``"IMGT"``, ``"Kabat"``, ``"Chothia"``, ``"Martin"``
    or ``"Aho"``. ``None`` on failure."""
    confidence: Optional[float]
    """Alignment confidence, between 0 and 1. ``None`` on failure."""
    numbering: Optional[dict[str, str]]
    """Position to residue, in sequence order. Positions are strings because insertions
    carry a letter, such as ``"111A"``. ``None`` on failure."""
    query_start: Optional[int]
    """0-based index of the first numbered residue in the sequence. ``None`` on failure."""
    query_end: Optional[int]
    """0-based index of the last numbered residue in the sequence, inclusive: the
    numbered residues are ``sequence[query_start : query_end + 1]``. ``None`` on
    failure."""
    error: Optional[str]
    """Why the sequence couldn't be numbered, or ``None`` on success."""
    error_kind: Optional[str]
    """What went wrong, as a stable code, or ``None`` on success: ``"invalid_sequence"``
    (shorter than 30 or longer than 10,000 residues, or holding a character that isn't a
    letter), ``"low_confidence"`` (no alignment reached ``min_confidence``) or, from
    ``number_domains`` only, ``"domain_too_short"`` (the best alignment is confident but
    too short to be a domain)."""


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
            (TRA, TRB, TRG, TRD) or ``"all"`` (every chain the scheme numbers: IGH, IGK
            and IGL with schemes other than IMGT). Pass all chains you want to consider;
            the annotator aligns against each and picks the one with the highest confidence.

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
            ``error_kind`` ``"low_confidence"``. Defaults to ``0.4``, which keeps every
            validated antibody and TCR domain and rejects sequences that aren't immune
            receptor domains. Fragments score lower; pass ``0.0`` to disable filtering.

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

        A domain that lacks its first IMGT positions (a light chain starting at position 2,
        say) and directly follows other residues, such as a linker, can have the residue just
        before it numbered as its first position: IMGT position 1 is so variable that the
        sequence alone can't tell a linker residue from the domain's own first residue.

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
            than a domain must be. So does an invalid sequence, with
            ``"invalid_sequence"``.
        """
        return [NumberingResult(**d) for d in self._annotator.number_domains(sequence)]

    def segment(self, sequence: str) -> SegmentationResult:
        """Split a sequence into FR/CDR regions.

        Args:
            sequence: Amino-acid sequence string (single-letter codes).

        Returns:
            A `SegmentationResult` with the residues before the domain as ``prefix``,
            ``fr1``–``fr4`` and ``cdr1``–``cdr3``, and the residues after it as
            ``postfix``, so the regions in order rebuild the sequence. On failure,
            ``error`` and ``error_kind`` are set and all region fields are ``None``.
        """
        return SegmentationResult(**self._annotator.segment(sequence))

    def segment_domains(self, sequence: str) -> list[SegmentationResult]:
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
            One `SegmentationResult` per domain, ordered by position; never empty. A
            sequence without a domain or an invalid one gives a single result with
            ``error`` and ``error_kind`` set, as `number_domains` describes.
        """
        return [
            SegmentationResult(**d) for d in self._annotator.segment_domains(sequence)
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
        chain: A single chain (not a group such as ``"ig"``) — ``"IGH"`` / ``"H"`` /
            ``"heavy"`` and so on. See `Annotator` for accepted values.

    Returns:
        ``{region: (start, end)}`` with both bounds inclusive, in N- to C-terminal order
        (``fr1``, ``cdr1``, … ``fr4``). IMGT and AHo number every chain alike; Kabat,
        Chothia and Martin place their CDRs differently on heavy and light chains.

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
    from immunum import (
        scheme_supports_chain,
    )

    assert scheme_supports_chain(
        "kabat", "H"
    )
    assert not scheme_supports_chain(
        "kabat", "B"
    )  # only IMGT numbers TCR chains
    ```

    Args:
        scheme: Numbering scheme, as for `Annotator`.
        chain: A single chain, as for `Annotator`.

    Returns:
        ``True`` when `Annotator` accepts the pair. IMGT numbers every chain; Kabat,
        Chothia, Martin and AHo number antibody chains (IGH, IGK, IGL) only.

    Raises:
        immunum.Error: If the scheme or chain is unrecognised (``kind``
            ``"invalid_scheme"`` or ``"invalid_chain"``).
    """
    return _scheme_supports_chain(scheme=scheme, chain=chain)
