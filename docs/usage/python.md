# Python

Build an `Annotator` once, with the chains to consider, the scheme and optionally
`min_confidence`, then call it on each sequence. Install with `pip install immunum`
([Getting started](../getting-started.md#python)); every parameter and field is in the
[API reference](../api/immunum.md).

## Creating an annotator

```python
from immunum import Annotator

antibodies = Annotator(chains=["H", "K", "L"], scheme="imgt")
same = Annotator(chains=["ig"], scheme="IMGT")              # a group, any case
tcrs = Annotator(chains=["tcr"], scheme="imgt", min_confidence=0.6)
```

- `chains`: the chains each sequence is scored against; the best match wins. Codes (`"H"`), names
  (`"IGH"`, `"heavy"`) and the groups `"ig"`, `"tcr"` and `"all"` are accepted, in any case; see
  [Chains](../concepts.md#chains). Each chain costs one alignment per sequence, so pass only those
  you expect.
- `scheme`: `"imgt"`, `"kabat"`, `"chothia"`, `"martin"` or `"aho"`, or their first letter. Only
  IMGT numbers TCR chains; see [Numbering schemes](../concepts.md#numbering-schemes).
- `min_confidence`: sequences whose best alignment scores below it, `0.5` by default, get a
  `low_confidence` result instead of a numbering. `0.0` turns the check off; see
  [Confidence](../concepts.md#confidence).

## Numbering

`number` gives every residue of the domain a position label:

```python
from immunum import Annotator

annotator = Annotator(chains=["ig"], scheme="imgt")
sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"

result = annotator.number(sequence)
print(result.chain)           # H
print(result.scheme)          # IMGT
print(result.numbering["1"])  # Q
print(result.numbering["111A"])  # K: insertions carry a letter, so positions are strings
print(result.query_start, result.query_end)  # 0 121: the numbered residues, 0-based inclusive
```

## Segmentation

`segment` splits the sequence into FR1–FR4 and CDR1–CDR3, with the residues before and after the
domain as `prefix` and `postfix`, so the regions always join back into the sequence:

```python
from immunum import Annotator

annotator = Annotator(chains=["ig"], scheme="imgt")
sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"

result = annotator.segment(sequence)
assert result.fr1 == "QVQLVQSGAEVKRPGSSVTVSCKAS"
assert result.cdr1 == "GGSFSTYA"
assert result.cdr2 == "VIPLLTIT"
assert result.cdr3 == "AREGTTGKPIGAFAH"
assert result.fr4 == "WGQGTLVTVSS"
assert "".join(result.as_dict().values()) == sequence
```

The regions follow the scheme you chose: Kabat puts CDR-H1 at positions 31–35, so the same
sequence gives a shorter CDR1:

```python
from immunum import Annotator

kabat = Annotator(chains=["ig"], scheme="kabat")
sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
print(kabat.segment(sequence).cdr1)  # TYALS
```

## Multiple domains

`number` and `segment` work on the best-scoring domain. `number_domains` and `segment_domains` do
the same for every variable domain in a sequence, such as both domains of an scFv, in sequence
order. Each result is what `number` or `segment` returns for that domain:

```python
from immunum import Annotator

heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
kappa = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
linker = "GGGGSGGGGSGGGGS"

annotator = Annotator(chains=["ig"], scheme="imgt")
domains = annotator.number_domains(heavy + linker + kappa)
assert [d.chain for d in domains] == ["H", "K"]

heavy_regions, kappa_regions = annotator.segment_domains(heavy + linker + kappa)
assert kappa_regions.prefix == linker
assert kappa_regions.cdr3 == "QQHYTTPPT"
```

Every residue lands in exactly one domain's regions: a domain's `prefix` holds the residues since
the previous domain, and only the last domain has a `postfix`. How domains are found, and a caveat
about IMGT position 1 after a linker, are in [Multiple domains](../concepts.md#multiple-domains).

## Region boundaries

`regions_for` gives a scheme's FR/CDR boundaries for a chain, both ends inclusive, without numbering
a sequence; `scheme_supports_chain` tells whether a scheme numbers a chain before you build an
`Annotator`:

```python
from immunum import regions_for, scheme_supports_chain

print(regions_for("kabat", "H")["cdr1"])  # (31, 35)
print(regions_for("imgt", "H")["cdr3"])   # (105, 117)

assert scheme_supports_chain("kabat", "H")
assert not scheme_supports_chain("kabat", "B")  # only IMGT numbers TCR chains
```

## Batches, threads and processes

An `Annotator` is safe to share between threads and can be pickled with its settings, so it works
with `concurrent.futures` thread and process pools and with `multiprocessing`. Numbering holds the
GIL, so a thread pool doesn't number faster than one thread; a process pool does, and for large
datasets [Polars](polars.md) runs on every core without managing a pool yourself.

```python
import pickle
from concurrent.futures import ThreadPoolExecutor

from immunum import Annotator

annotator = Annotator(chains=["ig"], scheme="imgt")
sequences = [
    "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS",
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
]

with ThreadPoolExecutor() as pool:
    chains = [r.chain for r in pool.map(annotator.number, sequences)]
assert chains == ["H", "K"]

copy = pickle.loads(pickle.dumps(annotator))
assert copy.number(sequences[0]).chain == "H"
```

## Errors

A setup mistake raises `immunum.Error`; a sequence that can't be numbered returns a result with
`error` and `error_kind` set and every other field `None`. See [Errors](errors.md#python).
