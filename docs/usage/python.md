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

By default, sequences with an alignment confidence below `0.5` get a result with `error` set
instead of a numbering. Pass `min_confidence=0.0` to disable this check, or raise the threshold
to filter non-immunoglobulin sequences more aggressively.

**Region boundaries** come from the same tables numbering assigns residues to, and can be read
without numbering a sequence:

```python
from immunum import regions_for

print(regions_for("kabat", "H")["cdr1"])  # (31, 35)
print(regions_for("imgt", "H")["cdr3"])  # (105, 117)
```

Both ends are inclusive. IMGT and AHo number every chain alike; Kabat, Chothia and Martin place
their CDRs differently on heavy and light chains, and only IMGT covers TCR chains.

`scheme_supports_chain` tells whether a scheme numbers a chain before you build an `Annotator`:

```python
from immunum import scheme_supports_chain

assert scheme_supports_chain("kabat", "H")
assert not scheme_supports_chain("kabat", "B")
```

## Errors

Mistakes in how you set immunum up raise `immunum.Error`, a `ValueError`: an unknown chain or
scheme, a scheme that doesn't number a chain, or a `min_confidence` outside `[0, 1]`. Its `kind`
attribute names what went wrong as a stable code (`invalid_chain`, `invalid_scheme`,
`unsupported_chain` or `invalid_min_confidence`); the message says it for people.

A sequence that can't be numbered doesn't raise, so a loop over many sequences never stops for
one bad one. Its result has every field `None` except `error`, the message, and `error_kind`:
`invalid_sequence` (too short, too long or not amino acids) or `low_confidence` (no alignment
reached `min_confidence`). Both are `None` on success. `number_domains` and `segment_domains`
never return an empty list: a sequence without a domain gives that one result, with
`low_confidence`, or with `domain_too_short` when the best alignment is confident but shorter
than a domain must be.

```python
import immunum

try:
    immunum.Annotator(chains=["IGX"], scheme="imgt")
except immunum.Error as e:
    print(e.kind)  # "invalid_chain"
    print(e)  # "unknown chain 'IGX' (options: ...)"

annotator = immunum.Annotator(chains=["ig"], scheme="imgt")
result = annotator.number("AAAA")
assert result.chain is None
print(result.error_kind)  # "invalid_sequence"
print(result.error)  # "sequence length 4 is below minimum 30"
```
