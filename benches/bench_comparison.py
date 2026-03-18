"""
Comparison benchmark: immunum vs antpack vs anarci.

Speed: 10000 IGH sequences timed via pytest-benchmark (3 rounds).
Correctness: residue-level accuracy per IMGT region embedded in benchmark.extra_info.

Run with:
    task benchmark-comparison
    uv run --group benchmark pytest tests/bench_comparison.py -v
"""

from __future__ import annotations

import math
import multiprocessing
import os
import signal
import time
from pathlib import Path
from typing import Annotated, Any, Optional

import polars
import pytest
import typer

import immunum.polars as imp
from immunum import Annotator

FIXTURES = Path(__file__).parent.parent / "fixtures" / "validation"
FIXTURE = FIXTURES / "ab_H_imgt.csv"
SAMPLE_SIZE = 10_000
SEED = 42

META_COLS = {"header", "sequence", "species"}

TCR_CHAINS = {"A", "B", "G", "D"}

IMMUNUM_CHAIN_MAP = {
    "H": "IGH",
    "K": "IGK",
    "L": "IGL",
    "A": "TRA",
    "B": "TRB",
    "G": "TRG",
    "D": "TRD",
}
RIOT_LOCUS_MAP = {
    "H": "igh",
    "K": "igk",
    "L": "igl",
    "A": "tra",
    "B": "trb",
    "G": "trg",
    "D": "trd",
}

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


def _pos_num(pos: str) -> int:
    return int("".join(c for c in pos if c.isdigit()))


def _region(pos: str) -> Optional[str]:
    n = _pos_num(pos)
    for region, (lo, hi) in IMGT_REGIONS.items():
        if lo <= n <= hi:
            return region
    return None


def _correctness(
    predictions: list[Optional[dict[str, str]]],
    gt: list[dict[str, str]],
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


# ---------------------------------------------------------------------------
# Runners (annotator/setup created outside the timed section)
# ---------------------------------------------------------------------------


def _run_immunum_multithreaded(
    df: polars.DataFrame, chains: list[str], scheme: str
) -> list[dict[str, str]]:
    match scheme:
        case "IMGT":
            scheme = "IMGT"
        case "KABAT":
            scheme = "Kabat"
        case _:
            raise ValueError(f"Invalid {scheme=}")
    result = df.select(
        imp.number(polars.col("sequence"), chains=chains, scheme=scheme)
        .alias("n")
        .struct.unnest()
    )

    out = [
        dict(zip(pos, res))
        for pos, res in zip(
            result.get_column("positions").to_list(),
            result.get_column("residues").to_list(),
        )
    ]

    return out


def _run_immunum_singlethreaded(
    sequences: list[tuple[str, str]], annotator: Annotator
) -> list[dict[str, str]]:
    out = []
    for _header, seq in sequences:
        n = annotator.number(seq)
        out.append(n.numbering)
    return out


def _run_antpack(
    sequences: list[tuple[str, str]], annotator, chain_letter: str
) -> list[Optional[dict[str, str]]]:
    results: list[Optional[dict[str, str]]] = []
    for _header, seq in sequences:
        numbering, _confidence, chain_type, err = annotator.analyze_seq(seq)
        if err or chain_type != chain_letter:
            results.append(None)
            continue
        d = {
            pos: seq[i]
            for i, pos in enumerate(numbering)
            if pos != "-" and i < len(seq)
        }
        results.append(d)
    return results


def _antpack_worker(
    args: tuple,
) -> list[Optional[dict[str, str]]]:
    sequences_chunk, chain_letter, scheme = args
    from antpack import SingleChainAnnotator

    annotator = SingleChainAnnotator([chain_letter], scheme=scheme)
    return _run_antpack(sequences_chunk, annotator, chain_letter)


def _run_antpack_parallel(
    sequences: list[tuple[str, str]], chain_letter: str, scheme: str
) -> list[Optional[dict[str, str]]]:
    ncpu = os.cpu_count() or 1
    chunk_size = max(1, (len(sequences) + ncpu - 1) // ncpu)
    chunks = [
        sequences[i : i + chunk_size] for i in range(0, len(sequences), chunk_size)
    ]
    with multiprocessing.Pool(ncpu) as pool:
        results = pool.map(_antpack_worker, [(c, chain_letter, scheme) for c in chunks])
    return [item for chunk_result in results for item in chunk_result]


def _run_anarci(
    sequences: list[tuple[str, str]], kwargs: dict
) -> list[Optional[dict[str, str]]]:
    from anarci import anarci  # type: ignore[missing-import]

    numbered_list, _, _ = anarci(sequences, **kwargs)
    results: list[Optional[dict[str, str]]] = []
    for numbered in numbered_list:
        if not numbered:
            results.append(None)
            continue
        domain_numbering, _, _ = numbered[0]
        d = {}
        for (num, ins), aa in domain_numbering:
            if aa != "-":
                d[str(num) + (ins.strip() if ins.strip() else "")] = aa
        results.append(d)
    return results


def _anarci_worker(
    args: tuple,
) -> list[Optional[dict[str, str]]]:
    sequences_chunk, scheme, chain_letter = args
    return _run_anarci(
        sequences_chunk, {"scheme": scheme, "allow": {chain_letter}, "ncpu": 1}
    )


def _run_anarci_parallel(
    sequences: list[tuple[str, str]], scheme: str, chain_letter: str
) -> list[Optional[dict[str, str]]]:
    ncpu = os.cpu_count() or 1
    chunk_size = max(1, (len(sequences) + ncpu - 1) // ncpu)
    chunks = [
        sequences[i : i + chunk_size] for i in range(0, len(sequences), chunk_size)
    ]
    with multiprocessing.Pool(ncpu) as pool:
        results = pool.map(_anarci_worker, [(c, scheme, chain_letter) for c in chunks])
    return [item for chunk_result in results for item in chunk_result]


def _run_anarcii2(
    sequences: list[tuple[str, str]], model, chain_letter: str
) -> list[Optional[dict[str, str]]]:
    numbered_list = model.number(sequences)
    results: list[Optional[dict[str, str]]] = []
    items = numbered_list.values() if isinstance(numbered_list, dict) else numbered_list
    for numbered in items:
        if (
            numbered is None
            or numbered.get("error")
            or numbered.get("chain_type") != chain_letter
        ):
            results.append(None)
            continue
        numbering = numbered.get("numbering")
        if not numbering:
            results.append(None)
            continue
        d = {}
        for (num, ins), aa in numbering:
            if aa != "-":
                d[str(num) + (ins.strip() if ins.strip() else "")] = aa
        results.append(d)
    return results


def _anarcii2_worker(
    args: tuple,
) -> list[Optional[dict[str, str]]]:
    sequences_chunk, chain_letter = args
    from anarcii import Anarcii  # type: ignore[missing-import]

    model = Anarcii(seq_type="antibody", mode="speed", ncpu=1)
    return _run_anarcii2(sequences_chunk, model, chain_letter)


def _run_anarcii2_parallel(
    sequences: list[tuple[str, str]], chain_letter: str
) -> list[Optional[dict[str, str]]]:
    ncpu = os.cpu_count() or 1
    chunk_size = max(1, (len(sequences) + ncpu - 1) // ncpu)
    chunks = [
        sequences[i : i + chunk_size] for i in range(0, len(sequences), chunk_size)
    ]
    with multiprocessing.Pool(ncpu) as pool:
        results = pool.map(_anarcii2_worker, [(c, chain_letter) for c in chunks])
    return [item for chunk_result in results for item in chunk_result]


def _run_riot(
    sequences: list[tuple[str, str]], riot, locus: str, riot_scheme
) -> list[Optional[dict[str, str]]]:
    results: list[Optional[dict[str, str]]] = []
    for header, seq in sequences:
        r = riot.run_on_sequence(header, seq, scheme=riot_scheme)
        if r is None or r.exc is not None or r.locus != locus:
            results.append(None)
            continue
        results.append({str(k): v for k, v in r.scheme_residue_mapping.items()})
    return results


def _riot_worker(
    args: tuple,
) -> list[Optional[dict[str, str]]]:
    sequences_chunk, locus, scheme_name = args
    from riot_na import create_riot_aa, Organism, Scheme  # type: ignore[import]

    riot = create_riot_aa(allowed_species=(Organism.HOMO_SAPIENS,))
    riot_scheme = getattr(Scheme, scheme_name.upper())
    return _run_riot(sequences_chunk, riot, locus, riot_scheme)


def _run_riot_parallel(
    sequences: list[tuple[str, str]], locus: str, scheme_name: str
) -> list[Optional[dict[str, str]]]:
    ncpu = os.cpu_count() or 1
    chunk_size = max(1, (len(sequences) + ncpu - 1) // ncpu)
    chunks = [
        sequences[i : i + chunk_size] for i in range(0, len(sequences), chunk_size)
    ]
    with multiprocessing.Pool(ncpu) as pool:
        results = pool.map(_riot_worker, [(c, locus, scheme_name) for c in chunks])
    return [item for chunk_result in results for item in chunk_result]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module", params=[FIXTURE])
def sample_df(request) -> polars.DataFrame:
    return polars.read_csv(request.param, infer_schema=False).sample(
        n=SAMPLE_SIZE, seed=SEED, with_replacement=True
    )


@pytest.fixture(scope="module")
def sample_seqs(sample_df: polars.DataFrame) -> list[tuple[str, str]]:
    return list(zip(sample_df["header"].to_list(), sample_df["sequence"].to_list()))


@pytest.fixture(scope="module")
def sample_gt(sample_df: polars.DataFrame) -> list[dict[str, str]]:
    position_cols = [c for c in sample_df.columns if c not in META_COLS]
    return [
        {pos: row[pos] for pos in position_cols if row[pos]}
        for row in sample_df.iter_rows(named=True)
    ]


@pytest.fixture(scope="module", params=[{"scheme": "imgt", "allow": {"H"}, "ncpu": 1}])
def anarci_kwargs(request) -> dict:
    pytest.importorskip("anarci")
    return request.param


@pytest.fixture(scope="module", params=[{"chains": ["IGH"], "scheme": "IMGT"}])
def immunum_annotator(request) -> Annotator:
    return Annotator(**request.param)


@pytest.fixture(scope="module")
def riot_annotator():
    pytest.importorskip("riot_na")
    from riot_na import create_riot_aa, Organism  # type: ignore[import]

    return create_riot_aa(allowed_species=(Organism.HOMO_SAPIENS,))


@pytest.fixture(
    scope="module", params=[{"seq_type": "antibody", "mode": "speed", "ncpu": 1}]
)
def anarcii2_annotator(request):
    pytest.importorskip("anarcii")
    from anarcii import Anarcii  # type: ignore[missing-import]

    return Anarcii(**request.param)


@pytest.fixture(scope="module", params=[{"chains": ["H"], "scheme": "imgt"}])
def antpack_annotator(request):
    pytest.importorskip("antpack")
    from antpack import SingleChainAnnotator

    return SingleChainAnnotator(**request.param)


# ---------------------------------------------------------------------------
# Benchmarks (speed + correctness embedded in extra_info)
# ---------------------------------------------------------------------------


def test_immunum_multithreaded(benchmark, sample_df, sample_gt):
    result = benchmark.pedantic(
        _run_immunum_multithreaded,
        args=(sample_df, ["IGH"], "IMGT"),
        rounds=3,
        iterations=1,
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_immunum_singlethreaded(benchmark, sample_seqs, sample_gt, immunum_annotator):
    result = benchmark.pedantic(
        _run_immunum_singlethreaded,
        args=(sample_seqs, immunum_annotator),
        rounds=3,
        iterations=1,
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_antpack(benchmark, sample_seqs, sample_gt, antpack_annotator):
    result = benchmark.pedantic(
        _run_antpack, args=(sample_seqs, antpack_annotator, "H"), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_antpack_parallel(benchmark, sample_seqs, sample_gt):
    pytest.importorskip("antpack")
    result = benchmark.pedantic(
        _run_antpack_parallel, args=(sample_seqs, "H", "imgt"), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_anarci(benchmark, sample_seqs, sample_gt, anarci_kwargs):
    result = benchmark.pedantic(
        _run_anarci, args=(sample_seqs, anarci_kwargs), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_anarci_parallel(benchmark, sample_seqs, sample_gt):
    pytest.importorskip("anarci")
    result = benchmark.pedantic(
        _run_anarci_parallel, args=(sample_seqs, "imgt", "H"), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_anarcii2(benchmark, sample_seqs, sample_gt, anarcii2_annotator):
    result = benchmark.pedantic(
        _run_anarcii2,
        args=(sample_seqs, anarcii2_annotator, "H"),
        rounds=3,
        iterations=1,
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_anarcii2_parallel(benchmark, sample_seqs, sample_gt):
    pytest.importorskip("anarcii")
    result = benchmark.pedantic(
        _run_anarcii2_parallel, args=(sample_seqs, "H"), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_riot(benchmark, sample_seqs, sample_gt, riot_annotator):
    from riot_na import Scheme  # type: ignore[import]

    result = benchmark.pedantic(
        _run_riot,
        args=(sample_seqs, riot_annotator, "igh", Scheme.IMGT),
        rounds=3,
        iterations=1,
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_riot_parallel(benchmark, sample_seqs, sample_gt):
    pytest.importorskip("riot_na")
    result = benchmark.pedantic(
        _run_riot_parallel, args=(sample_seqs, "igh", "IMGT"), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


# ---------------------------------------------------------------------------
# CLI entrypoint
# ---------------------------------------------------------------------------

BENCHMARK_SIZES = [100, 1_000, 10_000, 100_000]
_REGIONS = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4", "overall"]
_COL_W = [26, 5, 9] + [7] * len(_REGIONS)  # tool, round, time_s, regions

app = typer.Typer()


def _table_header() -> None:
    headers = ["tool", "round", "time_s"] + _REGIONS
    typer.echo("  ".join(h.ljust(w) for h, w in zip(headers, _COL_W)))
    typer.echo("  ".join("-" * w for w in _COL_W))


def _table_row(tool: str, r: int, elapsed: float, correctness: dict) -> None:
    vals = [tool, str(r), f"{elapsed:.3f}s"] + [
        f"{correctness[reg]:.2f}" for reg in _REGIONS
    ]
    typer.echo("  ".join(v.ljust(w) for v, w in zip(vals, _COL_W)))


def _parse_duration(s: str) -> float:
    """Parse human-readable duration string to seconds (e.g. '1h', '30m', '90s')."""
    s = s.strip()
    if s.endswith("h"):
        return float(s[:-1]) * 3600
    if s.endswith("m"):
        return float(s[:-1]) * 60
    if s.endswith("s"):
        return float(s[:-1])
    return float(s)


def _timed(fn: Any, args: tuple, timeout_s: float) -> tuple[float, Any]:
    def _handler(signum: int, frame: Any) -> None:
        raise TimeoutError

    old = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(math.ceil(timeout_s))
    try:
        t0 = time.perf_counter()
        result = fn(*args)
        elapsed = time.perf_counter() - t0
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)
    return elapsed, result


AB_CHAINS = {"H", "K", "L"}  # chains supported by antpack / anarcii2 / riot


@app.command()
def main(
    chain: Annotated[str, typer.Option(help="Chain letter: H, K, L, A, B, G, D")] = "H",
    output: Annotated[Path, typer.Option(help="Output file (.csv or .parquet)")] = Path(
        "benchmark_results.csv"
    ),
    rounds: Annotated[int, typer.Option(help="Repeats per tool/size")] = 3,
    sizes: Annotated[
        Optional[list[int]], typer.Option(help="Sample sizes (repeatable)")
    ] = None,
    timeout: Annotated[
        str, typer.Option(help="Per-function timeout (e.g. 1h, 30m, 90s)")
    ] = "1h",
    tools: Annotated[
        Optional[list[str]],
        typer.Option(
            help=(
                "Tools to run (repeatable). Matches exact name or '<tool>_*' prefix. "
                "E.g. --tools immunum --tools anarci. "
                "Available: immunum, immunum_multithreaded, immunum_singlethreaded, "
                "antpack, anarci, anarcii2, riot (each also has a _parallel variant)."
            )
        ),
    ] = None,
) -> None:
    chain = chain.upper()
    mol_type = "tcr" if chain in TCR_CHAINS else "ab"
    fixture = FIXTURES / f"{mol_type}_{chain}_imgt.csv"
    if not fixture.exists():
        typer.echo(f"Fixture not found: {fixture}", err=True)
        raise typer.Exit(1)

    immunum_chain = IMMUNUM_CHAIN_MAP[chain]
    locus = RIOT_LOCUS_MAP[chain]
    is_ab = chain in AB_CHAINS

    timeout_s = _parse_duration(timeout)
    run_sizes = sizes or BENCHMARK_SIZES
    df_full = polars.read_csv(fixture, infer_schema=False)
    rows: list[dict] = []

    for size in run_sizes:
        typer.echo(f"\n=== size={size} ===")
        _table_header()

        # Build static annotators once per size (outside the timed section).
        # Each entry: (name, fn, data_args_fn) where data_args_fn(df, seqs) -> args tuple.
        RunnerSpec = tuple[str, Any, Any]  # (name, fn, data_args_fn)
        runner_specs: list[RunnerSpec] = [
            (
                "immunum_multithreaded",
                _run_immunum_multithreaded,
                lambda df, seqs, _c=[immunum_chain]: (df, _c, "IMGT"),
            ),
            (
                "immunum_singlethreaded",
                _run_immunum_singlethreaded,
                (lambda ann: lambda df, seqs: (seqs, ann))(
                    Annotator(chains=[immunum_chain], scheme="IMGT")
                ),
            ),
        ]

        if is_ab:
            try:
                from antpack import SingleChainAnnotator  # type: ignore[import]

                ann = SingleChainAnnotator([chain], scheme="imgt")
                runner_specs += [
                    (
                        "antpack",
                        _run_antpack,
                        (lambda a, cl: lambda df, seqs: (seqs, a, cl))(ann, chain),
                    ),
                    (
                        "antpack_parallel",
                        _run_antpack_parallel,
                        lambda df, seqs, _c=chain: (seqs, _c, "imgt"),
                    ),
                ]
            except ImportError:
                typer.echo("  antpack not installed, skipping")

        try:
            import anarci as _anarci_mod  # type: ignore[import]  # noqa: F401

            anarci_kw: dict = {"scheme": "imgt", "allow": {chain}, "ncpu": 1}
            runner_specs += [
                (
                    "anarci",
                    _run_anarci,
                    (lambda kw: lambda df, seqs: (seqs, kw))(anarci_kw),
                ),
                (
                    "anarci_parallel",
                    _run_anarci_parallel,
                    lambda df, seqs, _c=chain: (seqs, "imgt", _c),
                ),
            ]
        except ImportError:
            typer.echo("  anarci not installed, skipping")

        if is_ab:
            try:
                from anarcii import Anarcii  # type: ignore[import]

                model = Anarcii(seq_type="antibody", mode="speed", ncpu=1)
                runner_specs += [
                    (
                        "anarcii2",
                        _run_anarcii2,
                        (lambda m, cl: lambda df, seqs: (seqs, m, cl))(model, chain),
                    ),
                    (
                        "anarcii2_parallel",
                        _run_anarcii2_parallel,
                        lambda df, seqs, _c=chain: (seqs, _c),
                    ),
                ]
            except ImportError:
                typer.echo("  anarcii not installed, skipping")

            try:
                from riot_na import create_riot_aa, Organism, Scheme  # type: ignore[import]

                riot = create_riot_aa(allowed_species=(Organism.HOMO_SAPIENS,))
                runner_specs += [
                    (
                        "riot",
                        _run_riot,
                        (lambda rt, lo, rs: lambda df, seqs: (seqs, rt, lo, rs))(
                            riot, locus, Scheme.IMGT
                        ),
                    ),
                    (
                        "riot_parallel",
                        _run_riot_parallel,
                        lambda df, seqs, _lo=locus: (seqs, _lo, "imgt"),
                    ),
                ]
            except ImportError:
                typer.echo("  riot-na not installed, skipping")

        if tools:
            runner_specs = [
                spec
                for spec in runner_specs
                if any(spec[0] == t or spec[0].startswith(t + "_") for t in tools)
            ]

        timed_out: set[str] = set()
        position_cols: list[str] = []

        for r in range(rounds):
            df = df_full.sample(n=size, seed=r, with_replacement=True)
            seqs = list(zip(df["header"].to_list(), df["sequence"].to_list()))
            if not position_cols:
                position_cols = [c for c in df.columns if c not in META_COLS]
            gt = [
                {pos: row[pos] for pos in position_cols if row[pos]}
                for row in df.iter_rows(named=True)
            ]

            for name, fn, data_args_fn in runner_specs:
                if name in timed_out:
                    continue
                args = data_args_fn(df, seqs)
                try:
                    elapsed, result = _timed(fn, args, timeout_s)
                except TimeoutError:
                    vals = (
                        [name, str(r)]
                        + [f"TIMEOUT (>{timeout})"]
                        + [""] * (len(_COL_W) - 3)
                    )
                    typer.echo("  ".join(v.ljust(w) for v, w in zip(vals, _COL_W)))
                    timed_out.add(name)
                    continue
                correctness = _correctness(result, gt)
                rows.append(
                    {
                        "tool": name,
                        "chain": chain,
                        "sample_size": size,
                        "round": r,
                        "time_s": elapsed,
                        **correctness,
                    }
                )
                _table_row(name, r, elapsed, correctness)

    result_df = polars.DataFrame(rows)
    if str(output).endswith(".parquet"):
        result_df.write_parquet(output)
    else:
        result_df.write_csv(output)
    typer.echo(f"\nSaved {len(rows)} rows to {output}")


if __name__ == "__main__":
    app()
