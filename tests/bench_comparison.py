"""
Comparison benchmark: immunum vs antpack vs anarci.

Speed: 10000 IGH sequences timed via pytest-benchmark (3 rounds).
Correctness: residue-level accuracy per IMGT region embedded in benchmark.extra_info.

Run with:
    task benchmark-comparison
    uv run --group benchmark pytest tests/bench_comparison.py -v
"""

from __future__ import annotations

import multiprocessing
import os
from pathlib import Path
from typing import Optional

import polars
import pytest

import immunum.polars as imp
from immunum import Annotator

FIXTURES = Path(__file__).parent.parent / "fixtures" / "validation"
FIXTURE = FIXTURES / "ab_H_imgt.csv"
SAMPLE_SIZE = 10_000
SEED = 42

META_COLS = {"header", "sequence", "species"}

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


def _run_immunum_multithreaded(df: polars.DataFrame) -> list[dict[str, str]]:
    result = df.select(
        imp.number(polars.col("sequence"), chains=["IGH"], scheme="IMGT")
        .alias("n")
        .struct.unnest()
    )

    out = [
        {"positions": pos, "residues": res}
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
    sequences: list[tuple[str, str]], annotator
) -> list[Optional[dict[str, str]]]:
    results: list[Optional[dict[str, str]]] = []
    for _header, seq in sequences:
        numbering, _confidence, chain_type, err = annotator.analyze_seq(seq)
        if err or chain_type != "H":
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
    sequences_chunk: list[tuple[str, str]],
) -> list[Optional[dict[str, str]]]:
    from antpack import SingleChainAnnotator

    annotator = SingleChainAnnotator(["H"], scheme="imgt")
    return _run_antpack(sequences_chunk, annotator)


def _run_antpack_parallel(
    sequences: list[tuple[str, str]],
) -> list[Optional[dict[str, str]]]:
    ncpu = os.cpu_count() or 1
    chunk_size = max(1, (len(sequences) + ncpu - 1) // ncpu)
    chunks = [
        sequences[i : i + chunk_size] for i in range(0, len(sequences), chunk_size)
    ]
    with multiprocessing.Pool(ncpu) as pool:
        results = pool.map(_antpack_worker, chunks)
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
    sequences_chunk: list[tuple[str, str]],
) -> list[Optional[dict[str, str]]]:
    return _run_anarci(sequences_chunk, {"scheme": "imgt", "allow": {"H"}, "ncpu": 1})


def _run_anarci_parallel(
    sequences: list[tuple[str, str]],
) -> list[Optional[dict[str, str]]]:
    ncpu = os.cpu_count() or 1
    chunk_size = max(1, (len(sequences) + ncpu - 1) // ncpu)
    chunks = [
        sequences[i : i + chunk_size] for i in range(0, len(sequences), chunk_size)
    ]
    with multiprocessing.Pool(ncpu) as pool:
        results = pool.map(_anarci_worker, chunks)
    return [item for chunk_result in results for item in chunk_result]


def _run_anarcii2(
    sequences: list[tuple[str, str]], model
) -> list[Optional[dict[str, str]]]:
    numbered_list = model.number(sequences)
    results: list[Optional[dict[str, str]]] = []
    items = numbered_list.values() if isinstance(numbered_list, dict) else numbered_list
    for numbered in items:
        if (
            numbered is None
            or numbered.get("error")
            or numbered.get("chain_type") != "H"
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
    sequences_chunk: list[tuple[str, str]],
) -> list[Optional[dict[str, str]]]:
    from anarcii import Anarcii  # type: ignore[missing-import]

    model = Anarcii(seq_type="antibody", mode="speed", ncpu=1)
    return _run_anarcii2(sequences_chunk, model)


def _run_anarcii2_parallel(
    sequences: list[tuple[str, str]],
) -> list[Optional[dict[str, str]]]:
    ncpu = os.cpu_count() or 1
    chunk_size = max(1, (len(sequences) + ncpu - 1) // ncpu)
    chunks = [
        sequences[i : i + chunk_size] for i in range(0, len(sequences), chunk_size)
    ]
    with multiprocessing.Pool(ncpu) as pool:
        results = pool.map(_anarcii2_worker, chunks)
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
        _run_immunum_multithreaded, args=(sample_df,), rounds=3, iterations=1
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
        _run_antpack, args=(sample_seqs, antpack_annotator), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_antpack_parallel(benchmark, sample_seqs, sample_gt):
    pytest.importorskip("antpack")
    result = benchmark.pedantic(
        _run_antpack_parallel, args=(sample_seqs,), rounds=3, iterations=1
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
        _run_anarci_parallel, args=(sample_seqs,), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_anarcii2(benchmark, sample_seqs, sample_gt, anarcii2_annotator):
    result = benchmark.pedantic(
        _run_anarcii2, args=(sample_seqs, anarcii2_annotator), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_anarcii2_parallel(benchmark, sample_seqs, sample_gt):
    pytest.importorskip("anarcii")
    result = benchmark.pedantic(
        _run_anarcii2_parallel, args=(sample_seqs,), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))
