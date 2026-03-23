# Scripts

Python helper scripts for generating consensus data used by the immunum-rs alignment engine. All commands should be run with `uv run` from the project root.

## generate_consensus.py

Builds consensus CSV files from numbered antibody/TCR sequences. These consensus files drive the scoring matrices used at compile time by `build.rs`.

Accepts either pre-numbered CSVs or IMGT-aligned FASTA files (128 positions per sequence) and outputs per-position amino acid frequencies, occupancy, and region labels.

```bash
uv run python -m scripts.generate_consensus csv input.csv output.csv
uv run python -m scripts.generate_consensus fasta input.fasta output.csv
```

## generate_vj_combinations.py

Generates all V+J gene protein combinations per chain type from IMGT nucleotide FASTA files. Translates V/J genes, concatenates every V+J pair, numbers them with AntPack, and writes IMGT-aligned FASTA output (128 positions). These can then be fed into `generate_consensus.py`.

```bash
uv run python -m scripts.generate_vj_combinations v_genes.fasta j_genes.fasta -o output_dir
```

## utils.py

Shared utilities: FASTA parsing, nucleotide-to-protein translation (standard codon table), and the canonical list of 128 IMGT positions.
