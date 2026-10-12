# Consensus files

One CSV per chain (`IGH`, `IGK`, `IGL`, `TRA`, `TRB`, `TRG`, `TRD`) with the amino acid frequencies at each
position of the IMGT unique numbering. `build.rs` turns them into the position-specific scoring matrices the
aligner uses.

## Source

The germline genes come from [IMGT®, the international ImMunoGeneTics information system®](https://www.imgt.org):
the IMGT/GENE-DB amino acid reference sequences with IMGT gaps,
[`IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP`](https://www.imgt.org/download/GENE-DB/),
**release 202641-7**, downloaded 12 October 2026.

IMGT data are licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). These files are derived
from them: they hold per-position statistics of germline V-J combinations, built as described below, not the
IMGT sequences themselves. `germlines.tsv` lists every IMGT allele used.

## Regenerating

```bash
uv run --no-sync python -m scripts.generate_consensus
```

This downloads the current IMGT/GENE-DB release into `.cache/imgt/` (delete it to fetch a newer one) and rewrites
the CSVs and `germlines.tsv`. The script only uses the Python standard library.

## How they're generated

1. **Germline alleles.** Take the IMGT/GENE-DB `V-REGION` and `J-REGION` entries that are functional (`F`) and
   complete (not partial), from these species:

   | Chains        | Species                                                                                       |
   |---------------|-----------------------------------------------------------------------------------------------|
   | IGH, IGK, IGL | human, mouse, rat, rabbit, rhesus macaque, pig, alpaca (IGH only), cow                        |
   | TRA–TRD       | human, mouse                                                                                  |

   These are the species ANARCI builds its HMMs from. Each gene contributes one allele: the lowest-numbered one
   that passes the checks below. TRAV/DV genes count towards TRA only.

2. **V genes on the IMGT numbering.** IMGT gaps the V genes of each species and locus in one alignment, so their
   positions 1–104 are taken as IMGT gives them. Where an allele has an insertion, though, IMGT adds a column for
   it to the whole alignment, which shifts every later position (rhesus IGHV, IGKV and IGLV, mouse TRAV and rat
   IGHV in this release). Each species' alignment is therefore profile-aligned to the human alignment of the
   same locus (BLOSUM62 and occupancy), keeping all 104 human positions and dropping the extra columns. Alleles
   with a residue in a dropped column, or without the 1st-CYS at 23 and the 2nd-CYS at 104, are skipped. The
   V-REGION residues after 104 open the CDR3.

3. **J genes.** The J-TRP or J-PHE of the `[FW]G.G` motif is IMGT position 118 (falling back to `[FW]..G`, as in
   human TRAJ16); FR4 is 118–128 and the residues before it close the CDR3. J genes with neither motif are skipped.

4. **V-J combinations.** Every V gene is joined to every J gene of the same species and chain. The CDR3 (V
   residues after 104, then J residues before 118) is placed on 105–117 with the IMGT rule: a CDR3 shorter than 13
   leaves positions empty in the order 111, 112, 110, 113, 109, 114, …; a longer one fills 105–111 and 112–117
   from either end, and the residues in between are IMGT insertions with no column here. Germline CDR3s carry no D
   gene or N region, so they are shorter than rearranged ones.

5. **Statistics.** All combinations of a chain are pooled, without weighting by species. Per position the CSV
   holds every residue (and the gap, `-`) at 5% or more, the fraction of combinations with a residue (occupancy),
   and the IMGT region.

6. **Scheme holes.** The Kabat-derived schemes leave IMGT 10 and 73 unnumbered on IGH and IGL, and 73 on IGK
   (`src/numbering/*.rs`). These positions are left out of the antibody CSVs, and V genes with a residue there
   are skipped (human IGLV(I)-70 in this release).

## Output format

| Column          | Description                                         |
|-----------------|-----------------------------------------------------|
| `position`      | IMGT position number                                |
| `consensus_aas` | Pipe-separated residues by descending frequency     |
| `frequencies`   | Their pipe-separated frequencies                    |
| `occupancy`     | Fraction of combinations with a residue here        |
| `region`        | IMGT region (FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4)  |

## Citing IMGT

- Sanou G, Zeitoun G, Manso T, et al. IMGT® at scale: FAIR, dynamic and automated tools for immune locus
  analysis. *Nucleic Acids Res.* 2025; gkaf1024. [doi:10.1093/nar/gkaf1024](https://doi.org/10.1093/nar/gkaf1024)
- Giudicelli V, Chaume D, Lefranc M-P. IMGT/GENE-DB: a comprehensive database for human and mouse immunoglobulin
  and T cell receptor genes. *Nucleic Acids Res.* 2005; 33:D256–D261.
  [doi:10.1093/nar/gki010](https://doi.org/10.1093/nar/gki010)
- Lefranc M-P, Pommié C, Ruiz M, et al. IMGT unique numbering for immunoglobulin and T cell receptor variable
  domains and Ig superfamily V-like domains. *Dev. Comp. Immunol.* 2003; 27:55–77.
  [doi:10.1016/s0145-305x(02)00039-3](https://doi.org/10.1016/s0145-305x(02)00039-3)
