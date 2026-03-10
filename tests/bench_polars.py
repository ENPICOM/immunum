from __future__ import annotations

from pathlib import Path

import pytest
import polars

import immunum.polars as imp

FIXTURES = Path(__file__).parent.parent / "fixtures" / "validation"

# (fixture_stem, chains, scheme)
FIXTURE_PARAMS = [
    ("ab_H_imgt", ["IGH"], "IMGT"),
    ("ab_K_imgt", ["IGK"], "IMGT"),
    ("ab_L_imgt", ["IGL"], "IMGT"),
    ("ab_H_kabat", ["IGH"], "Kabat"),
    ("ab_K_kabat", ["IGK"], "Kabat"),
    ("ab_L_kabat", ["IGL"], "Kabat"),
    ("tcr_A_imgt", ["TRA"], "IMGT"),
    ("tcr_B_imgt", ["TRB"], "IMGT"),
    ("tcr_G_imgt", ["TRG"], "IMGT"),
    ("tcr_D_imgt", ["TRD"], "IMGT"),
]


def load_sequences(stem: str) -> polars.DataFrame:
    return polars.read_csv(FIXTURES / f"{stem}.csv", columns=["sequence"])


def run_numbering(
    df: polars.DataFrame, chains: list[str], scheme: str
) -> polars.DataFrame:
    return df.select(
        imp.number(polars.col("sequence"), chains=chains, scheme=scheme).alias(
            "numbered"
        )
    )


@pytest.mark.parametrize(
    "stem,chains,scheme",
    FIXTURE_PARAMS,
    ids=[s for s, *_ in FIXTURE_PARAMS],
)
def test_benchmark_fixture(benchmark, stem, chains, scheme):
    df = load_sequences(stem)
    benchmark.pedantic(run_numbering, args=(df, chains, scheme), rounds=5, iterations=1)


LARGE_SAMPLE_PARAMS = [
    pytest.param(["IGH"], "IMGT", id="IGH-IMGT"),
    pytest.param(["IGH", "IGK", "IGL"], "IMGT", id="all_IG-IMGT"),
    pytest.param(["TRA", "TRB", "TRG", "TRD"], "IMGT", id="all_TR-IMGT"),
    pytest.param(
        ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"], "IMGT", id="all-IMGT"
    ),
]


@pytest.fixture(scope="module")
def large_sample():
    all_sequences = polars.concat([load_sequences(stem) for stem, *_ in FIXTURE_PARAMS])
    return all_sequences.sample(n=10_000, with_replacement=True, seed=42)


@pytest.mark.parametrize("chains,scheme", LARGE_SAMPLE_PARAMS)
def test_benchmark_1m_sequences(benchmark, large_sample, chains, scheme):
    benchmark.pedantic(
        run_numbering, args=(large_sample, chains, scheme), rounds=3, iterations=1
    )
