---
hide:
  # Previous/next would lead out of the tab, to the neighbouring tabs' pages
  - footer
---

# Benchmarks

## Current accuracy

immunum's own validation, re-run for every change to the numbering (`cargo run --release
--features dev-tools --bin benchmark`, which writes `BENCHMARKS.toml`). Each chain's reference
set in `fixtures/validation` is numbered under each scheme and compared residue by residue with
its reference numbering: "Numbered exactly" counts sequences with every position right, "Residues
correct" counts residues. Only IMGT numbers TCR chains (A, B, G, D). For how the antibody
references were made, see [Ground-truth labelling](#ground-truth-labelling).

<!-- benchmark-accuracy -->

## Comparison with other tools

immunum (the development version after 2.0.0, with its faster alignment) is compared with three
established numbering tools:
[AntPack](https://github.com/jlparkI/AntPack) 0.3.8.6.2,
[ANARCI](https://github.com/oxpig/ANARCI) 2026.2.13.2 and
[ANARCII](https://github.com/oxpig/ANARCII) 2.0.8. AntPack 0.3.8.6.2 is the release its author
recommends for benchmarking: the last one under the GPL, before AntPack's licence changed to
academic, non-commercial use only.

All results on this page were measured on 2026-10-15 on two otherwise idle cloud VMs of the same
size, 32 cores and 125 GiB of memory, running Ubuntu 26.04 and Python 3.12:

- **x86**: AMD EPYC 7543 (Zen 3)
- **ARM**: Arm Neoverse-N1, the core design of AWS Graviton2

### Speed

The time to number 10,000 antibody sequences from a sequencing run
(`fixtures/SRR_HKL_10k.fasta`: heavy, kappa and lambda chains mixed) in IMGT, with every tool
allowed the H, K and L chains. Mean of 5 rounds, each a fresh sample (with replacement) of the
file's 10,000 sequences. Lower is better. "Numbered" counts the sequences each tool returned a
numbering for; immunum turned away the few that fell below its default `min_confidence` of 0.5.

#### x86

<!-- benchmark-speed -->

On 1 core immunum numbers about 13,400 sequences a second: 1.9 times AntPack's speed, about 370
times ANARCI's and 340 times ANARCII's. On 32 cores it is 1.6 times as fast as AntPack.

<iframe src="../assets/benchmark_plot1_performance.html" width="100%" height="440px" frameborder="0" scrolling="no"></iframe>

#### ARM

<!-- benchmark-speed: arm/results_headline.csv -->

Each Neoverse-N1 core is slower than an EPYC core, so every tool takes longer on 1 core here, but
the order is the same: immunum is 1.3 times as fast as AntPack on 1 core and 2.0 times on 32
cores, and about 150 times ANARCI's and 300 times ANARCII's speed on 1 core.

### How benchmarks are run

`benches/bench_comparison.py` runs every tool through its Python API on the same sequences, once
on 1 core and once on all 32, with each tool's default settings:

| Tool | 1 core | 32 cores |
|---|---|---|
| immunum | `Annotator.number` in a loop | the `immunum.polars` expression |
| AntPack | `SingleChainAnnotator.analyze_seqs` | the same, split over a process pool of 32 workers |
| ANARCI | `anarci(..., ncpu=1)` | `run_anarci(..., ncpu=32)`, ANARCI's own process pool |
| ANARCII | `Anarcii(mode="speed")`, torch on 1 thread | the same, torch on 16 threads (x86) or 32 (ARM) |

Only the numbering is timed. Building annotators, loading ANARCII's model and starting AntPack's
process pool happen before the clock starts; `run_anarci` starts its pool inside the timed call,
as it does for every caller. ANARCII's all-core thread count is the fastest measured on each
machine (`--anarcii-threads`): numbering 1,000 sequences took 13.5 s on 16 threads on the x86 VM,
against 14.6 s on 8 and 15.6 s on 32, and 32 s on 32 threads on the ARM VM, against 82 s on 8
and 41 s on 1.

Three runs make the results on this page; `task benchmark` runs all three and makes the plots
(`task plots` makes the plots alone):

- **Speed** (`task benchmark-headline`): the 10,000 sequences above, 5 rounds.
- **Accuracy** (`task benchmark-accuracy`): each table in `fixtures/validation/*_imgt.csv`
  sampled to 1,000 sequences, 7 rounds, scored residue by residue per IMGT region (FR1–FR4,
  CDR1–CDR3). AntPack and ANARCII are run on the antibody chains only, so the TCR chains compare
  immunum with ANARCI. Every tool numbered identically on both machines.
- **Scaling** (`task benchmark-speed`): batches of 100 to 1,000,000 IGH sequences, 3 rounds each.
  A tool that takes longer than 10 minutes on a batch is left out of the larger ones.

On the ARM VM the scaling run's last batch, 1,000,000 sequences, ran out of memory: AntPack's 32
forked workers each copied the harness's memory. That batch was re-run on the same samples once
the harness froze Python's garbage collector before forking (`gc.freeze()`), which keeps that
memory shared; nothing else changed. The results are in `resources/benchmark_results` (x86) and
`resources/benchmark_results/arm` (ARM).

#### Ground-truth labelling

The antibody validation fixtures are derived from unique PDB structures where the
IMGT numbering assigned by antpack and ANARCI agrees between the two tools. As a
consequence, **antpack and ANARCI will by definition always achieve 100% correctness
on the antibody fixtures** — they are the source of the ground truth. This is a known
limitation, tracked in [#33](https://github.com/ENPICOM/immunum/issues/33).

If you want to add another tool to the comparison, please open an issue on our
[issue tracker](https://github.com/ENPICOM/immunum/issues).

### Correctness by segment

Per-region residue accuracy (%) for each tool, one row per chain: the antibody chains IGH, IGK
and IGL, and the TCR chains TRA, TRB, TRG and TRD. Only the single-core runs are shown; the
all-core runs number identically. FR and CDR columns correspond to the seven IMGT regions (FR1–FR4
and CDR1–CDR3). Because the antibody ground truth is derived from antpack/ANARCI agreement, those
two tools score 100% on IGH, IGK and IGL by construction — see the note above. A sequence a tool
doesn't number counts every one of its residues as wrong. That is most of what immunum loses on
TRA and TRB: its default `min_confidence` of 0.5 turns away about 9% of the TRA and 6% of the TRB
references, which `BENCHMARKS.toml`, run without a threshold, numbers at 99.4% and 99.8% of
residues correct. Bars are the mean over the 7 rounds; error bars span one standard deviation
above and below it.

<iframe src="../assets/benchmark_plot2_correctness.html" width="100%" height="700px" frameborder="0" scrolling="yes"></iframe>

### Scaling with batch size

Mean wall-clock time (log scale) as a function of batch size, from 100 to 1,000,000 IGH
sequences on the x86 VM, one panel for 1 core and one for 32 cores. Error bars span half a
standard deviation above and below the mean across rounds. On small batches the 32-core runs are
dominated by fixed costs, such as ANARCI starting its process pool. The ARM VM's run is in
`resources/benchmark_results/arm/results_speed.csv`.

<iframe src="../assets/benchmark_plot3_scaling.html" width="100%" height="1000px" frameborder="0" scrolling="no"></iframe>
