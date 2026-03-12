from typing import TypedDict

__version__: str

class NumberingResult(TypedDict):
    chain: str
    scheme: str
    positions: list[str]
    residues: list[str]
    confidence: float

class _Annotator:
    def __init__(
        self,
        chains: list[str],
        scheme: str,
    ): ...
    def number(self, sequence: str): ...
