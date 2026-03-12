# Annotator

::: immunum._internal.Annotator
    options:
      show_root_heading: false

## NumberingResult

::: immunum._internal.NumberingResult
    options:
      show_root_heading: false

## Example

```python
from immunum import Annotator

annotator = Annotator(chains=["H", "K", "L"], scheme="imgt")
result = annotator.number(
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYWMSWVRQAPGKGLEWVANIKQDGSEKYYVDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARDSGYGAMDYWGQGTLVTVSS"
)

assert result["chain"] == "H"
assert result["scheme"] == "IMGT"
assert isinstance(result["positions"], list)
assert isinstance(result["confidence"], float)
```
