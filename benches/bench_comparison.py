"""
Speed and accuracy comparison: immunum against antpack, ANARCI and ANARCII.

Every tool runs through its Python API on the same sequences, once on a single core and once on
every core. Only the numbering is timed: building annotators, loading models and starting
antpack's process pool happen before the clock starts. ANARCI's parallel mode is the exception,
because `run_anarci` starts its own process pool on every call. Every tool keeps its default
settings, immunum's `min_confidence` of 0.5 included.

Single core / every core, per tool:

- immunum: `Annotator.number` in a loop / the `immunum.polars` expression
- antpack: `SingleChainAnnotator.analyze_seqs` / the same, split over a process pool
- ANARCI: `anarci(..., ncpu=1)` / `run_anarci(..., ncpu=<cores>)`
- ANARCII: `Anarcii(mode="speed")` with torch held to one thread / to `--anarcii-threads` (8)

Run with:
    task benchmark-headline      # 10,000 NGS sequences, mixed H/K/L (the README figure)
    task benchmark-accuracy      # residue accuracy per chain
    task benchmark-speed         # time against batch size, IGH
    task benchmark-comparison    # the same runners under pytest-benchmark
"""

from __future__ import annotations

import gc
import math
import multiprocessing
import os
import platform
import signal
import time
from collections.abc import Callable, Iterator
from contextlib import AbstractContextManager, contextmanager
from dataclasses import dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Annotated, Any, Optional

import polars
import pytest
import typer

import immunum.polars as imp
from immunum import Annotator

FIXTURES = Path(__file__).parent.parent / "fixtures" / "validation"
SEED = 42

META_COLS = {"header", "sequence", "species"}

AB_CHAINS = {"H", "K", "L"}  # the chains antpack and ANARCII are compared on
TCR_CHAINS = {"A", "B", "G", "D"}

# Cores the benchmark may use: the affinity mask, not the machine, so `taskset` limits it
NCPU = (
    len(os.sched_getaffinity(0))
    if hasattr(os, "sched_getaffinity")
    else os.cpu_count() or 1
)

# IMGT region boundaries (numeric part of position string)
IMGT_REGIONS: dict[str, tuple[int, int]] = {
    "FR1": (1, 26),
    "CDR1": (27, 38),
    "FR2": (39, 55),
    "CDR2": (56, 65),
    "FR3": (66, 104),
    "CDR3": (105, 117),
    "FR4": (118, 128),
}

Numbering = Optional[dict[str, str]]


def _pos_num(pos: str) -> int:
    return int("".join(c for c in pos if c.isdigit()))


def _region(pos: str) -> Optional[str]:
    n = _pos_num(pos)
    for region, (lo, hi) in IMGT_REGIONS.items():
        if lo <= n <= hi:
            return region
    return None


def _correctness(
    predictions: list[Numbering], gt: list[dict[str, str]]
) -> dict[str, float]:
    """Return per-region accuracy percentages (and 'overall')."""
    stats: dict[str, dict[str, int]] = {
        r: {"correct": 0, "total": 0} for r in IMGT_REGIONS
    }
    for pred, expected in zip(predictions, gt):
        for pos, aa in expected.items():
            region = _region(pos)
            if region is None:
                continue
            stats[region]["total"] += 1
            if pred is not None and pred.get(pos) == aa:
                stats[region]["correct"] += 1
    result = {
        r: round(100 * v["correct"] / v["total"], 2) if v["total"] else 0.0
        for r, v in stats.items()
    }
    total_c = sum(v["correct"] for v in stats.values())
    total_t = sum(v["total"] for v in stats.values())
    result["overall"] = round(100 * total_c / total_t, 2) if total_t else 0.0
    return result


def _anarci_style_numbering(numbering: list) -> dict[str, str]:
    """ANARCI and ANARCII number as `((position, insertion), residue)`, `-` for a gap."""
    return {f"{num}{ins.strip()}": aa for (num, ins), aa in numbering if aa != "-"}


# ---------------------------------------------------------------------------
# Runners
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Runner:
    name: str
    """`<tool>_singlethreaded` / `<tool>` on one core, `_multithreaded` / `_parallel` on all."""
    prepare: Callable[[list[str]], AbstractContextManager[Callable[[], Any]]]
    """Untimed set-up for one batch; yields the call that is timed."""
    extract: Callable[[Any, list[str]], list[Numbering]]
    """Turns the timed call's result into a numbering per sequence, `None` where it failed."""


def _immunum_runners(chains: set[str]) -> list[Runner]:
    chain_list = sorted(chains)

    @contextmanager
    def single(seqs: list[str]) -> Iterator[Callable[[], Any]]:
        annotator = Annotator(chains=chain_list, scheme="IMGT")
        yield lambda: [annotator.number(s) for s in seqs]

    @contextmanager
    def multi(seqs: list[str]) -> Iterator[Callable[[], Any]]:
        df = polars.DataFrame({"sequence": seqs})
        yield lambda: df.select(
            imp.number(polars.col("sequence"), chains=chain_list, scheme="IMGT").alias(
                "n"
            )
        )

    def extract_multi(result: polars.DataFrame, _seqs: list[str]) -> list[Numbering]:
        # `numbering` is a list of {position, residue} structs, null when the sequence failed
        return [
            {r["position"]: r["residue"] for r in numbering}
            if numbering is not None
            else None
            for numbering in result.get_column("n").struct.field("numbering").to_list()
        ]

    return [
        Runner(
            "immunum_singlethreaded",
            single,
            lambda results, _seqs: [r.numbering for r in results],
        ),
        Runner("immunum_multithreaded", multi, extract_multi),
    ]


_antpack_annotator: Any = None


def _antpack_init(chains: list[str]) -> None:
    from antpack import SingleChainAnnotator  # type: ignore[import]

    global _antpack_annotator
    _antpack_annotator = SingleChainAnnotator(chains, scheme="imgt")


def _antpack_number(seqs: list[str]) -> list:
    return _antpack_annotator.analyze_seqs(seqs)


@contextmanager
def _frozen_heap() -> Iterator[None]:
    """Keep the harness's objects out of a forked worker's garbage collection, so the worker shares
    their memory instead of copying it: a collection touches every object it tracks. Without this,
    AntPack's 32 workers on a million sequences ran a 125 GB machine out of memory."""
    gc.collect()
    gc.freeze()
    try:
        yield
    finally:
        gc.unfreeze()


def _antpack_runners(chains: set[str]) -> list[Runner]:
    from antpack import SingleChainAnnotator  # type: ignore[import]

    chain_list = sorted(chains)

    @contextmanager
    def single(seqs: list[str]) -> Iterator[Callable[[], Any]]:
        annotator = SingleChainAnnotator(chain_list, scheme="imgt")
        yield lambda: annotator.analyze_seqs(seqs)

    @contextmanager
    def parallel(seqs: list[str]) -> Iterator[Callable[[], Any]]:
        # Four chunks per worker, so the fast cores of a hybrid CPU don't wait on the slow ones
        size = max(1, math.ceil(len(seqs) / (NCPU * 4)))
        chunks = [seqs[i : i + size] for i in range(0, len(seqs), size)]
        # Python's default start method, fork on Linux, as AntPack's users get it: forkserver
        # workers averaged 0.19 s on the 10,000 headline sequences against fork's 0.13 s. A forked
        # worker inherits the harness's heap, which is why `main` frees each tool's results
        # before the next tool runs.
        with (
            _frozen_heap(),
            multiprocessing.Pool(
                NCPU, initializer=_antpack_init, initargs=(chain_list,)
            ) as pool,
        ):
            pool.map(
                _antpack_number, [["EVQLVESGGG"]] * NCPU * 4
            )  # workers up and ready
            yield lambda: [
                r for chunk in pool.map(_antpack_number, chunks) for r in chunk
            ]

    def extract(raw: list, seqs: list[str]) -> list[Numbering]:
        # One label per residue of the input, "-" where the residue isn't numbered
        return [
            None
            if err or chain not in chains
            else {pos: seq[i] for i, pos in enumerate(labels) if pos != "-"}
            for seq, (labels, _identity, chain, err) in zip(seqs, raw)
        ]

    return [
        Runner("antpack", single, extract),
        Runner("antpack_parallel", parallel, extract),
    ]


def _anarci_runners(chains: set[str]) -> list[Runner]:
    from anarci import anarci, run_anarci  # type: ignore[import]

    @contextmanager
    def single(seqs: list[str]) -> Iterator[Callable[[], Any]]:
        named = [(str(i), s) for i, s in enumerate(seqs)]
        # `ncpu` here is hmmscan's thread count
        yield lambda: anarci(named, scheme="imgt", allow=chains, ncpu=1)[0]

    @contextmanager
    def parallel(seqs: list[str]) -> Iterator[Callable[[], Any]]:
        named = [(str(i), s) for i, s in enumerate(seqs)]
        # Splits the sequences over a process pool it starts itself, one hmmscan per worker
        with _frozen_heap():
            yield lambda: run_anarci(named, scheme="imgt", allow=chains, ncpu=NCPU)[1]

    def extract(numbered: list, _seqs: list[str]) -> list[Numbering]:
        return [_anarci_style_numbering(n[0][0]) if n else None for n in numbered]

    return [
        Runner("anarci", single, extract),
        Runner("anarci_parallel", parallel, extract),
    ]


# How many torch threads ANARCII gets on all cores. More threads can be slower, by a lot: numbering
# 1,000 sequences took 16 s on 1 thread, 5 s on 8 and 111 s on 24 on a 24-core Core Ultra 9 285HX;
# on a 32-core AMD EPYC 7543 VM, 8 threads took 14.6 s and 16 took 13.5 s; on a 32-core Ampere
# Neoverse-N1 VM, 8 threads took 82 s and 32 took 32 s. Measure a machine's best and pass it with
# `--anarcii-threads`.
DEFAULT_ANARCII_THREADS = min(8, NCPU)


def _anarcii_runners(chains: set[str], threads_parallel: int) -> list[Runner]:
    import torch
    from anarcii import Anarcii  # type: ignore[import]

    # Built once: loading the model takes seconds. `ncpu` only sets torch's thread count, which
    # is global, so each run sets it again rather than trusting the model it was built with.
    model = Anarcii(seq_type="antibody", mode="speed", ncpu=1)

    def runner(
        threads: int,
    ) -> Callable[[list[str]], AbstractContextManager[Callable[[], Any]]]:
        @contextmanager
        def prepare(seqs: list[str]) -> Iterator[Callable[[], Any]]:
            torch.set_num_threads(threads)
            yield lambda: model.number(seqs)

        return prepare

    def extract(numbered: dict, _seqs: list[str]) -> list[Numbering]:
        return [
            None
            if n.get("error")
            or n.get("chain_type") not in chains
            or not n.get("numbering")
            else _anarci_style_numbering(n["numbering"])
            for n in numbered.values()
        ]

    return [
        Runner("anarcii2", runner(1), extract),
        Runner("anarcii2_parallel", runner(threads_parallel), extract),
    ]


def make_runners(
    chains: set[str],
    tools: Optional[list[str]] = None,
    anarcii_threads: int = DEFAULT_ANARCII_THREADS,
) -> list[Runner]:
    """Every runner for `chains` (one-letter codes) whose tool is installed and supports them.

    `tools` keeps the runners named exactly so or `<tool>_*`; `anarcii_threads` is the torch thread
    count of ANARCII's all-core run."""
    factories: list[tuple[str, Callable[[set[str]], list[Runner]]]] = [
        ("immunum", _immunum_runners),
        ("anarci", _anarci_runners),
    ]
    if chains <= AB_CHAINS:
        factories += [
            ("antpack", _antpack_runners),
            ("anarcii", lambda c: _anarcii_runners(c, anarcii_threads)),
        ]
    runners: list[Runner] = []
    for module, factory in factories:
        try:
            runners += factory(chains)
        except ImportError:
            typer.echo(f"  {module} not installed, skipping", err=True)
    if tools:
        runners = [
            r
            for r in runners
            if any(r.name == t or r.name.startswith(t + "_") for t in tools)
        ]
    return runners


def _versions() -> dict[str, str]:
    rv = {"python": platform.python_version(), "cores": str(NCPU)}
    for package in ["immunum", "antpack", "anarci", "anarcii", "torch", "polars"]:
        try:
            rv[package] = version(package)
        except PackageNotFoundError:
            pass
    return rv


# ---------------------------------------------------------------------------
# pytest-benchmark: 10,000 IGH sequences from the validation set
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def sample_df() -> polars.DataFrame:
    return polars.read_csv(FIXTURES / "ab_H_imgt.csv", infer_schema=False).sample(
        n=10_000, seed=SEED, with_replacement=True
    )


@pytest.mark.parametrize("runner", make_runners({"H"}), ids=lambda r: r.name)
def test_comparison(benchmark, sample_df: polars.DataFrame, runner: Runner) -> None:
    seqs = sample_df["sequence"].to_list()
    position_cols = [c for c in sample_df.columns if c not in META_COLS]
    gt = [
        {p: row[p] for p in position_cols if row[p]}
        for row in sample_df.iter_rows(named=True)
    ]
    with runner.prepare(seqs) as run:
        raw = benchmark.pedantic(run, rounds=3, iterations=1)
    benchmark.extra_info.update(_correctness(runner.extract(raw, seqs), gt))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

_REGIONS = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4", "overall"]
_COL_W = [24, 5, 10, 9] + [7] * len(_REGIONS)  # tool, round, time_s, numbered, regions

app = typer.Typer()


def _echo_row(vals: list[str]) -> None:
    typer.echo("  ".join(v.ljust(w) for v, w in zip(vals, _COL_W)).rstrip())


def _parse_duration(s: str) -> float:
    """Parse human-readable duration string to seconds (e.g. '1h', '30m', '90s')."""
    s = s.strip()
    for suffix, factor in (("h", 3600), ("m", 60), ("s", 1)):
        if s.endswith(suffix):
            return float(s[:-1]) * factor
    return float(s)


def _timed(fn: Callable[[], Any], timeout_s: float) -> tuple[float, Any]:
    def _handler(signum: int, frame: Any) -> None:
        raise TimeoutError

    old = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(math.ceil(timeout_s))
    try:
        t0 = time.perf_counter()
        result = fn()
        elapsed = time.perf_counter() - t0
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)
    return elapsed, result


def _read_fasta(path: Path) -> list[str]:
    seqs: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            seqs.append("")
        elif seqs:
            seqs[-1] += line.strip()
    return seqs


@app.command()
def main(
    chain: Annotated[
        str, typer.Option(help="Validation set to sample: H, K, L, A, B, G or D")
    ] = "H",
    fasta: Annotated[
        Optional[Path],
        typer.Option(
            help=(
                "Sample this FASTA of antibody sequences instead, with H, K and L allowed. It "
                "has no reference numbering, so only the time and numbered count are recorded."
            )
        ),
    ] = None,
    output: Annotated[Path, typer.Option(help="Output file (.csv or .parquet)")] = Path(
        "benchmark_results.csv"
    ),
    rounds: Annotated[int, typer.Option(help="Repeats per tool/size")] = 3,
    sizes: Annotated[
        Optional[list[int]], typer.Option(help="Sample sizes (repeatable)")
    ] = None,
    timeout: Annotated[
        str,
        typer.Option(
            help="Per-run timeout (e.g. 1h, 30m, 90s); a tool that times out skips larger sizes"
        ),
    ] = "1h",
    tools: Annotated[
        Optional[list[str]],
        typer.Option(
            help=(
                "Tools to run (repeatable). Matches exact name or '<tool>_*' prefix. "
                "E.g. --tools immunum --tools anarci. Available: immunum, antpack, anarci, "
                "anarcii2, and their _singlethreaded/_multithreaded or _parallel variants."
            )
        ),
    ] = None,
    no_correctness_check: Annotated[
        bool,
        typer.Option(
            "--no-correctness-check",
            help="Skip extract and correctness computation; only time the run step.",
        ),
    ] = False,
    anarcii_threads: Annotated[
        int,
        typer.Option(
            help="Torch threads for ANARCII's all-core run; the fastest count differs per machine"
        ),
    ] = DEFAULT_ANARCII_THREADS,
) -> None:
    if fasta is not None:
        chains, label = set(AB_CHAINS), "".join(sorted(AB_CHAINS))
        df_full = polars.DataFrame({"sequence": _read_fasta(fasta)})
        check = False
    else:
        chain = chain.upper()
        chains, label = {chain}, chain
        fixture = (
            FIXTURES / f"{'tcr' if chain in TCR_CHAINS else 'ab'}_{chain}_imgt.csv"
        )
        if not fixture.exists():
            typer.echo(f"Fixture not found: {fixture}", err=True)
            raise typer.Exit(1)
        df_full = polars.read_csv(fixture, infer_schema=False)
        check = not no_correctness_check

    typer.echo(
        ", ".join(f"{k} {v}" for k, v in _versions().items())
        + f", ANARCII all-core threads {anarcii_threads}"
    )
    timeout_s = _parse_duration(timeout)
    runners = make_runners(chains, tools, anarcii_threads)
    position_cols = [c for c in df_full.columns if c not in META_COLS]
    # Without scoring, a sample needs only the sequences, not ~150 reference columns per row
    sample = df_full if check else df_full.select("sequence")
    rows: list[dict] = []
    timed_out: set[str] = set()

    for size in sorted(sizes or [len(df_full)]):
        typer.echo(f"\n=== size={size} ===")
        _echo_row(["tool", "round", "time_s", "numbered"] + (_REGIONS if check else []))
        _echo_row(["-" * w for w in _COL_W[: 4 + (len(_REGIONS) if check else 0)]])

        for r in range(rounds):
            df = sample.sample(n=size, seed=SEED + r, with_replacement=True)
            seqs = df["sequence"].to_list()
            gt = (
                [
                    {p: row[p] for p in position_cols if row[p]}
                    for row in df.iter_rows(named=True)
                ]
                if check
                else []
            )
            for runner in runners:
                if runner.name in timed_out:
                    continue
                with runner.prepare(seqs) as run:
                    try:
                        elapsed, raw = _timed(run, timeout_s)
                    except TimeoutError:
                        _echo_row([runner.name, str(r), f"TIMEOUT (>{timeout})"])
                        timed_out.add(runner.name)
                        continue
                if no_correctness_check:
                    numbered, correctness = None, {}
                else:
                    numbering = runner.extract(raw, seqs)
                    numbered = sum(n is not None for n in numbering)
                    correctness = _correctness(numbering, gt) if check else {}
                    del numbering
                # Freed before the next tool runs: a million sequences' results take gigabytes
                del raw
                rows.append(
                    {
                        "tool": runner.name,
                        "chain": label,
                        "sample_size": size,
                        "round": r,
                        "time_s": elapsed,
                        "numbered": numbered,
                        **correctness,
                    }
                )
                _echo_row(
                    [
                        runner.name,
                        str(r),
                        f"{elapsed:.3f}s",
                        "-" if numbered is None else str(numbered),
                    ]
                    + [
                        f"{correctness[reg]:.2f}"
                        for reg in _REGIONS
                        if reg in correctness
                    ]
                )

        # After every size, so a run that dies keeps the sizes it finished
        result_df = polars.DataFrame(rows)
        if output.suffix == ".parquet":
            result_df.write_parquet(output)
        else:
            result_df.write_csv(output)
    typer.echo(f"\nSaved {len(rows)} rows to {output}")


if __name__ == "__main__":
    app()
