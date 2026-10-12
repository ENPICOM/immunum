# Scripts

Python helper scripts. All commands should be run with `uv run` from the project root.

## generate_consensus.py

Builds the consensus CSVs in `resources/consensus/`, which `build.rs` compiles into the scoring matrices, from
the IMGT/GENE-DB germline reference sequences. It downloads the IMGT data into `.cache/imgt/`, numbers every V-J
germline combination per species and chain on the IMGT numbering, and writes per-position amino acid frequencies,
occupancy and region labels, plus `germlines.tsv` with the alleles used. The method and IMGT citations are in
[`resources/consensus/README.md`](../resources/consensus/README.md).

```bash
uv run --no-sync python -m scripts.generate_consensus
```
