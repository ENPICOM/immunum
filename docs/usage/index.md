# Introduction

As advertised, immunum annotates antibody and T-cell receptor sequences with IMGT or Kabat position
numbers (`number` functions), or segments them into FR/CDR regions plus optional prefixes (if sequence contains more than `immunum` is confident about). Segmentation is guaranteed to cover the whole sequence, i.e. string sum of segments always equals input sequence.

Two interfaces are provided:

- **Python** (`immunum.Annotator`) — for single sequences or small batches
- **Polars** (`immunum.polars`) — for large datasets; runs as a lazy expression inside
  a Polars query plan

When multiple chains are provided, `immunum` scores each sequence against every chain type you pass and picks the best-matching one. The more chains you include, the more comparisons are made per
sequence -— so pass only the chains you expect in your data. E.g. including 6 chains in `chains=...` argument will make numbering/segmentation six times slower.

--8<-- "docs/usage/python.md"

--8<-- "docs/usage/polars.md"
