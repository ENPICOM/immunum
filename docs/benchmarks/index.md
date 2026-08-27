# Benchmarks

immunum is compared against three established antibody numbering tools:
[antpack~=0.2.7](https://github.com/jlparkI/AntPack),
[ANARCI](https://github.com/oxpig/ANARCI), and
[anarcii2](https://github.com/oxpig/ANARCII).

## How benchmarks are run

Two separate benchmark runs are combined to produce the plots on this page:

- **Accuracy** (`task benchmark-accuracy`): each table from `fixtures/validation/*.csv` is
sampled to 1,000 sequences and annotated 7 rounds with every tool. Accuracy is measured
as residue-level correctness per IMGT region (FR1–FR4, CDR1–CDR3).

- **Speed** (`task benchmark-speed`): batch sizes from 100 to 1,000,000 sequences (10× steps)
are timed 3 rounds each using IGH/IMGT. Only the annotation step is timed -— annotator
construction and result extraction are excluded.

Plots are generated from the saved CSV files with `task plots`.

### Ground-truth labelling

The validation fixtures are derived from unique PDB structures where the
IMGT numbering assigned by antpack and ANARCI agrees between the two tools. As a
consequence, **antpack and ANARCI will by definition always achieve 100% correctness
on these fixtures** — they are the source of the ground truth. This is a known
limitation; we intend to fix that in the next release.

### Known issues

The runs for ANARCI, ANARCII, and Antpack may not reflect best-case performance due
to how they are invoked in the benchmark harness (for instance, `antpack_parallel` scaling is very weird). Fixes are tracked in:

- ANARCI: [#0](https://github.com/ENPICOM/immunum-rs/issues)
- ANARCII: [#0](https://github.com/ENPICOM/immunum-rs/issues)
- Antpack: [#0](https://github.com/ENPICOM/immunum-rs/issues)

If you want to add another tool to the comparison, please open an issue on our
[issue tracker](https://github.com/ENPICOM/immunum/issues).

## Correctness by segment

Per-region residue accuracy (%) for each tool, averaged across IGH/IGK/IGL chains.
Only single-threaded variants are shown. FR and CDR columns correspond to the seven
IMGT regions (FR1–FR4 and CDR1–CDR3). Because the ground truth is derived from
antpack/ANARCI agreement, those two tools score 100% by construction — see the note
above.

<iframe src="../assets/benchmark_plot2_correctness.html" width="100%" height="700px" frameborder="0" scrolling="yes"></iframe>

### SAbDab2-paired: full PDB chains

`fixtures/validation/sabdab2paired_{H,K,L}_imgt.csv` numbers full PDB SEQRES chains from
paired antibody structures in
[SAbDab](http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/)
(`resources/opig/sabdab2_paired_and_vhh_pdb_sequences.csv`, `chain_type` "heavy"/"light",
sourced from `sabdab2_ab_split.csv`), numbered with legacy/original ANARCI.
Unlike the fixtures above, these sequences are **not** pre-trimmed to the variable
domain — most rows include the constant region (CH1/CL) appended directly after FR4, so
a tool has to locate the domain boundary itself rather than being handed just the domain.
Light chains are split into kappa/lambda based on ANARCI's own chain-type call, not known
ahead of time.

Ground truth here comes from legacy/original ANARCI alone (`anarci_original`), not from
antpack/ANARCI agreement like the fixtures above — there's no cross-tool agreement check,
so **ANARCI scores 100% by construction on this dataset** (it is the ground truth), the
same limitation noted above but for a single tool rather than two.

<iframe src="../assets/benchmark_plot2_correctness_sabdab2paired.html" width="100%" height="500px" frameborder="0" scrolling="yes"></iframe>

### SAbDab2-nano: full PDB chains, VHH

`fixtures/validation/sabdab2nano_H_imgt.csv` numbers full PDB SEQRES chains from
single-domain (nanobody) structures in SAbDab
(`resources/opig/sabdab2_paired_and_vhh_pdb_sequences.csv`, `chain_type` "VHH", sourced
from `sabdab_summary_all_single_domain`). Same "not pre-trimmed to the variable domain"
caveat as SAbDab2-paired above. A VHH has no paired light chain but is otherwise a
heavy-chain fold, so it's numbered as heavy (`--chain H`) both by ANARCI (ground truth)
and by every tool below; the row is labelled "VHH" in the plot purely for display.

<iframe src="../assets/benchmark_plot2_correctness_sabdab2nano.html" width="100%" height="250px" frameborder="0" scrolling="yes"></iframe>

## Throughput at fixed batch size

Sequences annotated per second at a fixed batch size of 10,000 IGH sequences, shown
separately for single-threaded and multi-threaded execution. Error bars span at standard
deviation across rounds. Higher is better.

<iframe src="../assets/benchmark_plot1_performance.html" width="100%" height="440px" frameborder="0" scrolling="no"></iframe>

## Scaling with batch size

Mean wall-clock time (log scale) as a function of batch size, from 100 to 1,000,000
sequences. Separate panels for single-threaded and multi-threaded execution.
Deviations (especially in the multi-threaded panel) reflect parallelisation overhead and process-spawning costs.

<iframe src="../assets/benchmark_plot3_scaling.html" width="100%" height="1000px" frameborder="0" scrolling="no"></iframe>
