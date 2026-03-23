from __future__ import annotations

import rtoml
from pathlib import Path

import pytest

polars = pytest.importorskip("polars")

import immunum.polars as imp  # noqa: E402

FIXTURES = Path(__file__).parent.parent / "fixtures" / "validation"
BENCHMARKS = Path(__file__).parent.parent / "BENCHMARKS.toml"
META_COLS = {"header", "sequence", "species"}

# (fixture_stem, chains, scheme, benchmark_key)
VALIDATION_FIXTURES = [
    ("ab_H_imgt", ["IGH"], "IMGT", "imgt.H"),
    ("ab_K_imgt", ["IGK"], "IMGT", "imgt.K"),
    ("ab_L_imgt", ["IGL"], "IMGT", "imgt.L"),
    ("ab_H_kabat", ["IGH"], "Kabat", "kabat.H"),
    ("ab_K_kabat", ["IGK"], "Kabat", "kabat.K"),
    ("ab_L_kabat", ["IGL"], "Kabat", "kabat.L"),
    # ("tcr_A_imgt", ["TRA"], "IMGT", "imgt.A"),
    # ("tcr_B_imgt", ["TRB"], "IMGT", "imgt.B"),
    ("tcr_G_imgt", ["TRG"], "IMGT", "imgt.G"),
    ("tcr_D_imgt", ["TRD"], "IMGT", "imgt.D"),
]


def get_benchmark_threshold(benchmark_key: str) -> float:
    """Return the known perfect_pct from BENCHMARKS.toml for the given key."""
    with open(BENCHMARKS) as f:
        data = rtoml.load(f)
    section, chain = benchmark_key.split(".")
    return data[section][chain]["perfect_pct"]


def compare_fixture(csv_path: Path, chains: list[str], scheme: str) -> tuple[int, int]:
    # no cover: start
    """Returns (mismatches, total) for a validation fixture."""
    df = polars.read_csv(csv_path, infer_schema=False)
    position_cols = [c for c in df.columns if c not in META_COLS]

    result = df.select(
        [
            "header",
            "sequence",
            *position_cols,
            imp.number(
                polars.col("sequence"), chains=chains, scheme=scheme, min_confidence=0.0
            ).alias("numbered"),
        ]
    )

    mismatches = 0
    for row in result.iter_rows(named=True):
        expected = {pos: aa for pos in position_cols if (aa := row[pos])}
        numbered = row["numbered"]
        if numbered is None or numbered["positions"] is None:
            mismatches += 1
            continue
        got = dict(zip(numbered["positions"], numbered["residues"]))
        if got != expected:
            mismatches += 1

    # no cover: stop
    return mismatches, result.height


SEQ = "SALTQPPAVSGTPGQRVTISCSGSDIGRRSVNWYQQFPGTAPKLLIYSNDQRPSVVPDRFSGSKSGTSASLAISGLQSEDEAEYYCAAWDDSLAVFGGGTQLTVGQPKA"
IGH_SEQ = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"


class TestPolarsNumber:
    def test_number_returns_expr(self):
        expr = imp.number(polars.col("sequence"), chains=["IGH"], scheme="IMGT")
        assert isinstance(expr, polars.Expr)

    def test_number_on_dataframe(self):
        df = polars.DataFrame({"sequence": [IGH_SEQ]})
        result = df.select(
            imp.number(polars.col("sequence"), chains=["IGH"], scheme="IMGT").alias(
                "numbered"
            )
        )
        assert "numbered" in result.columns
        assert result.height == 1

    def test_number_struct_fields(self):
        df = polars.DataFrame({"sequence": [IGH_SEQ]})
        result = df.select(
            imp.number(polars.col("sequence"), chains=["IGH"], scheme="IMGT").alias(
                "numbered"
            )
        ).unnest("numbered")
        assert "chain" in result.columns
        assert "positions" in result.columns
        assert "residues" in result.columns
        assert "error" in result.columns

    def test_number_error_field_null_on_success(self):
        df = polars.DataFrame({"sequence": [IGH_SEQ]})
        result = df.select(
            imp.number(polars.col("sequence"), chains=["IGH"], scheme="IMGT").alias(
                "numbered"
            )
        ).unnest("numbered")
        assert result["error"][0] is None

    def test_number_error_field_set_on_failure(self):
        df = polars.DataFrame({"sequence": ["AAAAAAAAAAAAAAAA"]})
        result = df.select(
            imp.number(polars.col("sequence"), chains=["IGH"], scheme="IMGT").alias(
                "numbered"
            )
        ).unnest("numbered")
        assert result["error"][0] is not None
        assert result["chain"][0] is None

    def test_number_multiple_sequences(self):
        df = polars.DataFrame({"sequence": [IGH_SEQ, SEQ]})
        result = df.select(
            imp.number(
                polars.col("sequence"), chains=["IGH", "IGK", "IGL"], scheme="IMGT"
            ).alias("numbered")
        )
        assert result.height == 2


@pytest.mark.slow
@pytest.mark.parametrize(
    "stem,chains,scheme,benchmark_key",
    VALIDATION_FIXTURES,
    ids=[s for s, *_ in VALIDATION_FIXTURES],
)
class TestValidationFixtures:
    def test_accuracy(self, stem, chains, scheme, benchmark_key):
        csv_path = FIXTURES / f"{stem}.csv"
        mismatches, total = compare_fixture(csv_path, chains, scheme)
        perfect = total - mismatches
        perfect_pct = 100 * perfect / total
        threshold = get_benchmark_threshold(benchmark_key)
        assert round(perfect_pct, 2) >= threshold, (
            f"{stem}: {mismatches}/{total} mismatched "
            f"({perfect_pct:.2f}% perfect, expected >= {threshold}%)"
        )


class TestPolarsSegment:
    def test_segment_returns_expr(self):
        expr = imp.segment(polars.col("sequence"), chains=["IGH"], scheme="IMGT")
        assert isinstance(expr, polars.Expr)

    def test_segment_on_dataframe(self):
        df = polars.DataFrame({"sequence": [IGH_SEQ]})
        result = df.select(
            imp.segment(polars.col("sequence"), chains=["IGH"], scheme="IMGT").alias(
                "segmented"
            )
        )
        assert "segmented" in result.columns
        assert result.height == 1

    def test_segment_struct_fields(self):
        df = polars.DataFrame({"sequence": [IGH_SEQ]})
        result = df.select(
            imp.segment(polars.col("sequence"), chains=["IGH"], scheme="IMGT").alias(
                "segmented"
            )
        ).unnest("segmented")
        expected_fields = {
            "fr1", "fr2", "fr3", "fr4",
            "cdr1", "cdr2", "cdr3",
            "prefix", "postfix", "error",
        }
        assert expected_fields.issubset(set(result.columns))

    def test_segment_error_field_null_on_success(self):
        df = polars.DataFrame({"sequence": [IGH_SEQ]})
        result = df.select(
            imp.segment(polars.col("sequence"), chains=["IGH"], scheme="IMGT").alias(
                "segmented"
            )
        ).unnest("segmented")
        assert result["error"][0] is None

    def test_segment_error_field_set_on_failure(self):
        df = polars.DataFrame({"sequence": ["AAAAAAAAAAAAAAAA"]})
        result = df.select(
            imp.segment(polars.col("sequence"), chains=["IGH"], scheme="IMGT").alias(
                "segmented"
            )
        ).unnest("segmented")
        assert result["error"][0] is not None
        assert result["fr1"][0] is None

    def test_segment_multiple_sequences(self):
        df = polars.DataFrame({"sequence": [IGH_SEQ, SEQ]})
        result = df.select(
            imp.segment(
                polars.col("sequence"), chains=["IGH", "IGK", "IGL"], scheme="IMGT"
            ).alias("segmented")
        )
        assert result.height == 2


class TestPolarsNumberingMethod:
    def test_segmentation_method_returns_expr(self):
        from immunum import Annotator

        annotator = Annotator(["IGH"], "IMGT")
        expr = imp.segmentation_method(polars.col("sequence"), annotator=annotator)
        assert isinstance(expr, polars.Expr)

    def test_segmentation_method_on_dataframe(self):
        from immunum import Annotator

        annotator = Annotator(["IGH"], "IMGT")
        df = polars.DataFrame({"sequence": [IGH_SEQ]})
        result = df.select(
            imp.segmentation_method(polars.col("sequence"), annotator=annotator).alias(
                "numbered"
            )
        )
        assert "numbered" in result.columns
        assert result.height == 1

    def test_numbering_method_returns_expr(self):
        from immunum import Annotator

        annotator = Annotator(["IGH"], "IMGT")
        expr = imp.numbering_method(polars.col("sequence"), annotator=annotator)
        assert isinstance(expr, polars.Expr)

    def test_numbering_method_on_dataframe(self):
        from immunum import Annotator

        annotator = Annotator(["IGH"], "IMGT")
        df = polars.DataFrame({"sequence": [IGH_SEQ]})
        result = df.select(
            imp.numbering_method(polars.col("sequence"), annotator=annotator).alias(
                "numbered"
            )
        )
        assert "numbered" in result.columns
        assert result.height == 1

    def test_numbering_method_struct_fields(self):
        from immunum import Annotator

        annotator = Annotator(["IGH"], "IMGT")
        df = polars.DataFrame({"sequence": [IGH_SEQ]})
        result = df.select(
            imp.numbering_method(polars.col("sequence"), annotator=annotator).alias(
                "numbered"
            )
        ).unnest("numbered")
        assert "chain" in result.columns
        assert "scheme" in result.columns
        assert "confidence" in result.columns
        assert "numbering" in result.columns
        assert "error" in result.columns
