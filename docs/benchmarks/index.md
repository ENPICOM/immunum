# Benchmarks

!!! note "These results predate immunum 2.0"
    The speed and scaling plots were produced with immunum 0.9.2 (March 2026) and the
    correctness plot with immunum 1.0.0 (April 2026). They haven't been re-run for 2.0.

immunum is compared against three established antibody numbering tools:
[antpack~=0.2.7](https://github.com/jlparkI/AntPack),
[ANARCI](https://github.com/oxpig/ANARCI), and
[ANARCII](https://github.com/oxpig/ANARCII).

## How benchmarks are run

Two separate benchmark runs are combined to produce the plots on this page:

- **Accuracy** (`task benchmark-accuracy`): each table from `fixtures/validation/*_imgt.csv` is
sampled to 1,000 sequences and annotated 7 rounds with every tool. Accuracy is measured
as residue-level correctness per IMGT region (FR1–FR4, CDR1–CDR3). antpack 0.2.7 and ANARCII
are run on antibody chains only, so the TCR chains compare immunum with ANARCI.

- **Speed** (`task benchmark-speed`): batch sizes from 100 to 1,000,000 sequences (100, 500,
1,000, 5,000 and so on) are timed 3 rounds each using IGH/IMGT. For most tools only the
annotation step is timed; annotator construction and result extraction are excluded. The
multi-threaded and parallel variants are the exception: multi-threaded immunum builds its
annotator with the Polars expression, ANARCII's parallel variant builds its model, and
antpack's starts a process pool with an annotator per worker, all inside the timed call.

Plots are generated from the saved CSV files with `task plots`.

### Ground-truth labelling

The antibody validation fixtures are derived from unique PDB structures where the
IMGT numbering assigned by antpack and ANARCI agrees between the two tools. As a
consequence, **antpack and ANARCI will by definition always achieve 100% correctness
on the antibody fixtures** — they are the source of the ground truth. This is a known
limitation, tracked in [#33](https://github.com/ENPICOM/immunum/issues/33).

### Known issues

The runs for ANARCI, ANARCII, and Antpack may not reflect best-case performance due
to how they are invoked in the benchmark harness (for instance, `antpack_parallel` scaling is very weird). Fixes are tracked in:

- ANARCI: [#30](https://github.com/ENPICOM/immunum/issues/30)
- ANARCII: [#31](https://github.com/ENPICOM/immunum/issues/31)
- Antpack: [#32](https://github.com/ENPICOM/immunum/issues/32)

If you want to add another tool to the comparison, please open an issue on our
[issue tracker](https://github.com/ENPICOM/immunum/issues).

## Correctness by segment

Per-region residue accuracy (%) for each tool, one row per chain: the antibody chains IGH, IGK
and IGL, and the TCR chains TRA, TRB, TRG and TRD. Only single-threaded variants are shown. FR
and CDR columns correspond to the seven IMGT regions (FR1–FR4 and CDR1–CDR3). Because the
antibody ground truth is derived from antpack/ANARCI agreement, those two tools score 100% on
IGH, IGK and IGL by construction — see the note above. On the TCR chains ANARCI isn't the ground
truth: it scores 96.9% on TRA, 99.1% on TRB, 99.9% on TRG and 93.0% on TRD. The TRD fixture
(all TRDV1) follows IMGT's gapping of TRDV1, with the CDR2 at 56, 57 and 65 and residues at 69–72;
ANARCI puts the CDR2 at 56–59 and 63–65 and leaves 69–73 empty.

<iframe src="../assets/benchmark_plot2_correctness.html" width="100%" height="700px" frameborder="0" scrolling="yes"></iframe>

## Throughput at fixed batch size

Sequences annotated per second at a fixed batch size of 10,000 IGH sequences, shown
separately for single-threaded and multi-threaded execution. Error bars span half a standard
deviation above and below the mean across rounds. Higher is better.

<iframe src="../assets/benchmark_plot1_performance.html" width="100%" height="440px" frameborder="0" scrolling="no"></iframe>

## Scaling with batch size

Mean wall-clock time (log scale) as a function of batch size, from 100 to 1,000,000
sequences. Separate panels for single-threaded and multi-threaded execution.
Deviations (especially in the multi-threaded panel) reflect parallelisation overhead and process-spawning costs.

<iframe src="../assets/benchmark_plot3_scaling.html" width="100%" height="1000px" frameborder="0" scrolling="no"></iframe>
