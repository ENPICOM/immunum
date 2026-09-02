# Python module

Here, interface is very simple: first you create an `Annotator` object with fixed chain types and numbering scheme you need (and optional `min_confidence` value), then call `number()` or `segment()` on each sequence.

**Numbering** assigns a position label to every residue, in the scheme you selected (`imgt`,
`kabat`, `chothia`, `martin` or `aho`; only `imgt` covers TCR chains):

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
from immunum import Annotator

annotator = Annotator(chains=["H", "K", "L"], scheme="imgt")
result = annotator.segment(
    "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
)
print(result.cdr3)  # "AREGTTGKPIGAFAH"
print(result.fr4)   # "WGQGTLVTVSS"
```

By default, sequences with an alignment confidence below `0.5` raise a `ValueError`.
Pass `min_confidence=0.0` to disable this check, or raise the threshold to filter
non-immunoglobulin sequences more aggressively.

**Region boundaries** come from the same tables numbering assigns residues to, and can be read
without numbering a sequence:

```python
from immunum import regions_for

print(regions_for("kabat", "H")["cdr1"])  # (31, 35)
print(regions_for("imgt", "H")["cdr3"])  # (105, 117)
```

Both ends are inclusive. IMGT and AHo number every chain alike; Kabat, Chothia and Martin place
their CDRs differently on heavy and light chains, and only IMGT covers TCR chains.
