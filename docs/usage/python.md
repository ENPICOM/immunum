# Python module

Here, interface is very simple: first you create an `Annotator` object with fixed chain types and numbering scheme you need (and optional `min_confidence` value), then call `number()` or `segment()` on each sequence.

**Numbering** assigns an IMGT or Kabat position label to every residue:

```python
from immunum import Annotator

annotator = Annotator(chains=["H", "K", "L"], scheme="imgt")

result = annotator.number(
    "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
)
print(result.chain)       # "H"
print(result.scheme)      # "IMGT"
print(result.numbering["1"])  # "Q"
```

**Segmentation** splits the sequence into FR1–FR4 and CDR1–CDR3 regions plus prefix/postfix:

```python
result = annotator.segment(
    "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
)
print(result.cdr3)  # "AREGTTGKPIGAFAH"
print(result.fr4)   # "WGQGTLVTVSS"
```

By default, sequences with an alignment confidence below `0.5` raise a `ValueError`.
Pass `min_confidence=0.0` to disable this check, or raise the threshold to filter
non-immunoglobulin sequences more aggressively.
