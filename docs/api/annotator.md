# Annotator

::: immunum._internal._Annotator
    options:
      show_root_heading: false

## Example

```python
from immunum import Annotator

annotator = Annotator(chains=["H", "K", "L"], scheme="imgt")
result = annotator.number(
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYWMSWVRQAPGKGLEWVANIKQDGSEKYYVDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARDSGYGAMDYWGQGTLVTVSS"
)

assert result.chain == "H"
assert result.scheme == "IMGT"
```
