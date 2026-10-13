# Upgrading to 2.0

immunum 2.0 makes every interface report errors the same way, finds every domain in a sequence,
and fixes how flanking residues are segmented. Most code keeps working; this page lists what can
break, by interface. The [changelog](../changelog.md) has every change.

## Every interface

**Error messages lost their category prefix.** The category is now a stable code: `error_kind` on
a result that failed (`errorKind` in JavaScript), and `kind` on an error that's raised or thrown.
Match on the code, not the message:

| 1.x message | 2.0 message | 2.0 code |
|---|---|---|
| `Invalid sequence: sequence length 4 is below minimum 30` | `sequence length 4 is below minimum 30` | `invalid_sequence` |
| `Low confidence: 0.0416 < threshold 0.5000` | `alignment confidence 0.0416 is below min_confidence 0.5000` | `low_confidence` |
| `Invalid chain type: unknown chain 'IGX' …` | `unknown chain 'IGX' …` | `invalid_chain` |

The codes are listed under [Errors](errors.md#error-kinds).

**`segment` keeps every residue.** Residues the aligner leaves out before or after the domain
used to be dropped; they're now in `prefix` and `postfix`, so the regions always join back into
the sequence. A sequence with a leader or a constant region can get a longer `prefix` or
`postfix` than in 1.x, but the same FR and CDR regions.

**`min_confidence` must be in `[0, 1]`** in every interface; only Python checked it before.

## Python

```python
import immunum

annotator = immunum.Annotator(chains=["ig"], scheme="imgt")
result = annotator.number("AAAA")

# 1.x: result.error.startswith("Invalid sequence")
assert result.error_kind == "invalid_sequence"

try:
    immunum.Annotator(chains=["IGX"], scheme="imgt")
except ValueError as e:  # immunum.Error is a ValueError, so 1.x handlers still catch it
    assert isinstance(e, immunum.Error)
    assert e.kind == "invalid_chain"
```

- `SegmenationResult` is renamed `SegmentationResult`.
- `SegmentationResult` lists its regions in sequence order, `prefix` first and `postfix` last, so
  `as_dict()` and the dataclass fields put `prefix` before `fr1`. In 1.x they came after `fr4`.
- Polars is an extra: install with `pip install "immunum[polars]"` to use `immunum.polars`.

## Polars

**`number` returns a `numbering` list** of `{position, residue}` structs, with the fields Python's
`Annotator.number` returns: `chain`, `scheme`, `confidence`, `numbering`, `query_start`,
`query_end`, `error` and `error_kind`. 1.x returned `positions` and `residues` lists instead of
`numbering`. Code that read those lists can rebuild them:

```python
import polars as pl
import immunum.polars as imp

df = pl.DataFrame({"sequence": [
    "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
]})
numbering = df.select(
    imp.number("sequence", chains=["ig"], scheme="imgt").struct.field("numbering")
)
lists = numbering.select(
    pl.col("numbering").list.eval(pl.element().struct.field("position")).alias("positions"),
    pl.col("numbering").list.eval(pl.element().struct.field("residue")).alias("residues"),
)
assert lists["positions"][0][:3].to_list() == ["1", "2", "3"]
assert lists["residues"][0][:3].to_list() == ["Q", "V", "Q"]
```

Or explode `numbering` and unnest it for one row per residue.

**Numbering after a leader is fixed.** For a sequence with residues before the domain, such as a
leader, 1.x Polars `number` paired each position with the residue that many places too early
(position 1 got the leader's first residue). 2.0 numbers it as every other interface does, so
these rows change.

- `numbering_method(expr, annotator=...)` and `segmentation_method(expr, annotator=...)` are
  deprecated and warn; use `number(expr, annotator=...)` and `segment(expr, annotator=...)`, which
  return the same.
- An unknown chain or scheme name raises `immunum.Error` when you build the expression, not when
  the query runs.
- Results carry `error_kind` next to `error`.

## JavaScript

```js
const { Annotator } = require("immunum");
const annotator = new Annotator(["ig"], "imgt");
const sequence =
  "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";

const result = annotator.number(sequence);
result.queryStart; // 1.x: result.query_start
result.queryEnd;   // 1.x: result.query_end
result.errorKind;  // new: the code for a failed result

try {
  new Annotator(["IGX"], "imgt");
} catch (e) {
  // 1.x threw a string; 2.0 throws an Error with a `kind`
  console.log(e.kind);    // "invalid_chain"
  console.log(e.message); // "unknown chain 'IGX' …"
}
```

- Results are camelCase throughout: `queryStart`, `queryEnd`, `errorKind`. Every other field name
  is the same as in the other interfaces.
- A failed `segment` result sets every region to `null`; 1.x left the region fields out.

## CLI

- TSV output has an `error_kind` column after `error`, and JSON records an `error_kind` field. The
  columns before it are unchanged.
- An argument is read as a sequence only when it's all letters. A path to a file that doesn't
  exist is an error (exit status 1); 1.x numbered it as a sequence and reported an invalid
  sequence.

## Rust

Errors are split by how they're reported:

- `immunum::Error` is a setup mistake: `InvalidChain`, `InvalidScheme`, `UnsupportedChain` (was
  `InvalidScheme`), `InvalidMinConfidence`, `InvalidPosition`, `WrongSequence` (was
  `InvalidSequence`).
- `immunum::SequenceError` is one sequence that couldn't be numbered: `InvalidSequence`,
  `LowConfidence`, `DomainTooShort`. `Annotator::number`, `segment`, `domains`, `number_domains`
  and `segment_domains` return it.

```rust
use immunum::{Annotator, Chain, Scheme};

let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT, None)?; // immunum::Error
match annotator.number("AAAA") {
    Ok(result) => println!("{:?}", result.chain),
    Err(e) => println!("{}: {e}", e.kind()), // invalid_sequence: sequence length 4 …
}
```

- Both enums are `#[non_exhaustive]`: add a `_` arm when you match on them.
- `AlignmentError`, `ConsensusParseError`, `PositionMappingError` and `Io` are gone.
- `immunum::Result<T, E = Error>` takes the error type as a second parameter.
- `NumberedRecord` and `SegmentedRecord` hold `result: Result<_, SequenceError>` and are built with
  `new`, instead of `success`/`failure` and separate `result`/`error` fields.
- `Chain` and `Scheme` parse errors are `immunum::Error` instead of `strum::ParseError`.
- `Chain::parse_chain_spec(spec)` is removed; use `Chain::parse_names(spec.split(','))`.
- The `cli` feature no longer brings `chrono` and `toml`, and the crate needs Rust 1.89 or newer.
