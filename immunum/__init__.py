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
