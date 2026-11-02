# Errors

immunum fails in two ways, and every interface reports each the same way:

- **A setup mistake is raised** before any sequence is numbered, such as an unknown chain name.
- **A sequence that can't be numbered is returned**, so a batch never stops for it: its result
  has `error` (the message) and `error_kind` (a stable code) set, and every other field empty.
  On success both are empty.

Match on the code, not the message: messages are for people and can be reworded.

## Error kinds

| Kind                     | Raised or returned | When                                                         |
| ------------------------ | ------------------ | ------------------------------------------------------------ |
| `invalid_chain`          | raised             | an unknown chain name, or no chains                          |
| `invalid_scheme`         | raised             | an unknown scheme name                                       |
| `unsupported_chain`      | raised             | a scheme other than IMGT asked to number a TCR chain         |
| `invalid_min_confidence` | raised             | `min_confidence` outside `[0, 1]`                            |
| `invalid_sequence`       | returned           | a sequence shorter than 30 or longer than 10,000 residues, or holding a non-letter |
| `low_confidence`         | returned           | no alignment reaches `min_confidence` (see [Confidence](../concepts.md#confidence)) |
| `domain_too_short`       | returned           | from `*_domains` only: the best alignment is confident but shorter than a domain (30 residues), and no other domain is found |

`number_domains` and `segment_domains` never return an empty list: a sequence without a domain, or
an invalid one, gives a single result with the error set.

## Python

A setup mistake raises `immunum.Error`, a `ValueError`, with the code as `kind`. A failed result
has every field `None` except `error` and `error_kind`:

```python
import immunum

try:
    immunum.Annotator(chains=["IGX"], scheme="imgt")
except immunum.Error as e:
    print(e.kind)  # invalid_chain
    print(e)       # unknown chain 'IGX' (options: ...)

annotator = immunum.Annotator(chains=["ig"], scheme="imgt")
result = annotator.number("AAAA")
assert result.chain is None
print(result.error_kind)  # invalid_sequence
print(result.error)       # sequence length 4 is below minimum 30

failed = [r for r in map(annotator.number, ["AAAA", "A" * 40]) if r.error_kind is not None]
assert [r.error_kind for r in failed] == ["invalid_sequence", "low_confidence"]
```

## Polars

A setup mistake raises `immunum.Error` when you build the expression, before any query runs. A
sequence that can't be numbered never fails the query: its struct has every field null except
`error` and `error_kind`, so `error_kind.is_not_null()` finds the failures.

A **null sequence** isn't numbered and isn't a failure either: its struct has every field null,
`error` and `error_kind` too, and `number_domains` and `segment_domains` give it a null list. So
`error.is_null()` is true for a null input as well as for a success. Test success with a field that
only success fills, such as `chain` from `number` or `fr1` from `segment`:

```python
import polars as pl
import immunum
import immunum.polars as imp

try:
    imp.number("sequence", chains=["IGX"], scheme="imgt")
except immunum.Error as e:
    print(e.kind)  # invalid_chain

heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
df = pl.DataFrame({"sequence": [heavy, "AAAA", "A" * 40, None]}).with_columns(
    imp.number("sequence", chains=["ig"], scheme="imgt").alias("numbered")
)
print(
    df.select(
        succeeded=pl.col("numbered").struct.field("chain").is_not_null(),
        failed=pl.col("numbered").struct.field("error_kind").is_not_null(),
        error_kind=pl.col("numbered").struct.field("error_kind"),
    )
)
# shape: (4, 3)
# ┌───────────┬────────┬──────────────────┐
# │ succeeded ┆ failed ┆ error_kind       │
# │ ---       ┆ ---    ┆ ---              │
# │ bool      ┆ bool   ┆ str              │
# ╞═══════════╪════════╪══════════════════╡
# │ true      ┆ false  ┆ null             │
# │ false     ┆ true   ┆ invalid_sequence │
# │ false     ┆ true   ┆ low_confidence   │
# │ false     ┆ false  ┆ null             │
# └───────────┴────────┴──────────────────┘
```

## JavaScript

A setup mistake throws an `Error` with the code as `kind` (`ImmunumError` in the TypeScript types).
A failed result has every field `null` except `error` and `errorKind`:

```js
const { Annotator } = require("immunum");

try {
  new Annotator(["IGX"], "imgt");
} catch (err) {
  console.log(err.kind);    // "invalid_chain"
  console.log(err.message); // "unknown chain 'IGX' (options: ...)"
}

const failed = new Annotator(["ig"], "imgt").number("AAAA");
console.log(failed.errorKind); // "invalid_sequence"
console.log(failed.error);     // "sequence length 4 is below minimum 30"
```

## Rust

Setting up returns `immunum::Error`; numbering one sequence returns `immunum::SequenceError`, so a
sequence's result can only fail with `InvalidSequence` or `LowConfidence`, plus `DomainTooShort`
from the domain search. Both have `kind()`, the code every other interface reports, and both are
`#[non_exhaustive]`:

```rust
use immunum::{Annotator, Chain, Scheme};

let annotator = Annotator::new(&[Chain::IGH], Scheme::IMGT, None)?; // immunum::Error
match annotator.number("AAAA") {
    Ok(result) => println!("{}", result.chain),
    Err(e) => println!("{}: {e}", e.kind()), // invalid_sequence: sequence length 4 …
}
```

## Command line

A setup or input mistake prints `error: …` to stderr and exits with status 1. A sequence that can't
be numbered gets a record with the `error` and `error_kind` columns (or fields) set, the run
carries on, a one-line summary of how many failed goes to stderr, and the exit status stays 0. See
[Exit status](cli.md#exit-status).
