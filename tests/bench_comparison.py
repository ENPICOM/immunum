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
    from antpack import SingleChainAnnotator  # type: ignore

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


def _run_anarci(sequences: list[tuple[str, str]]) -> list[Optional[dict[str, str]]]:
    from anarci import anarci  # type: ignore

    numbered_list, _, _ = anarci(sequences, scheme="imgt", allow={"H"}, ncpu=1)
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
    return _run_anarci(sequences_chunk)


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


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def sample_df() -> polars.DataFrame:
    return polars.read_csv(FIXTURE, infer_schema=False).sample(
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


# ---------------------------------------------------------------------------
# Benchmarks (speed + correctness embedded in extra_info)
# ---------------------------------------------------------------------------


def test_immunum_multithreaded(benchmark, sample_df, sample_gt):
    result = benchmark.pedantic(
        _run_immunum_multithreaded, args=(sample_df,), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_immunum_singlethreaded(benchmark, sample_seqs, sample_gt):
    annotator = Annotator(chains=["IGH"], scheme="IMGT")
    result = benchmark.pedantic(
        _run_immunum_singlethreaded,
        args=(sample_seqs, annotator),
        rounds=3,
        iterations=1,
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_antpack(benchmark, sample_seqs, sample_gt):
    pytest.importorskip("antpack")
    from antpack import SingleChainAnnotator  # type: ignore

    annotator = SingleChainAnnotator(["H"], scheme="imgt")
    result = benchmark.pedantic(
        _run_antpack, args=(sample_seqs, annotator), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_antpack_parallel(benchmark, sample_seqs, sample_gt):
    pytest.importorskip("antpack")
    result = benchmark.pedantic(
        _run_antpack_parallel, args=(sample_seqs,), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


def test_anarci(benchmark, sample_seqs, sample_gt):
    pytest.importorskip("anarci")
    result = benchmark.pedantic(
        _run_anarci, args=(sample_seqs,), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))


@pytest.mark.skip
def test_anarci_parallel(benchmark, sample_seqs, sample_gt):
    pytest.importorskip("anarci")
    result = benchmark.pedantic(
        _run_anarci_parallel, args=(sample_seqs,), rounds=3, iterations=1
    )
    benchmark.extra_info.update(_correctness(result, sample_gt))
