from immunum._internal import _Annotator  # noqa: F401
from dataclasses import dataclass


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
    ):
        self._annotator = _Annotator(chains=chains, scheme=scheme)

    def number(self, sequence: str) -> NumberingResult:
        return NumberingResult(**self._annotator.number(sequence))
