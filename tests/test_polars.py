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
    ("ab_H_chothia", ["IGH"], "Chothia", "chothia.H"),
    ("ab_K_chothia", ["IGK"], "Chothia", "chothia.K"),
    ("ab_L_chothia", ["IGL"], "Chothia", "chothia.L"),
    ("ab_H_martin", ["IGH"], "Martin", "martin.H"),
    ("ab_K_martin", ["IGK"], "Martin", "martin.K"),
    ("ab_L_martin", ["IGL"], "Martin", "martin.L"),
    ("ab_H_aho", ["IGH"], "Aho", "aho.H"),
    ("ab_K_aho", ["IGK"], "Aho", "aho.K"),
    ("ab_L_aho", ["IGL"], "Aho", "aho.L"),
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
        if numbered is None or numbered["numbering"] is None:
            mismatches += 1
            continue
        got = {r["position"]: r["residue"] for r in numbered["numbering"]}
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

    def test_number_accepts_aliases_and_groups(self):
        df = polars.DataFrame({"sequence": [IGH_SEQ]})

        def numbering(chains, scheme):
            row = df.select(
                imp.number(polars.col("sequence"), chains=chains, scheme=scheme).alias(
                    "n"
                )
            ).unnest("n")
            return {r["position"]: r["residue"] for r in row["numbering"][0]}

        canonical = numbering(["IGH", "IGK", "IGL"], "IMGT")
        assert numbering(["ig"], "i") == canonical
        assert numbering(["heavy", "k", "lambda"], "imgt") == canonical

    def test_number_unknown_chain_raises(self):
        df = polars.DataFrame({"sequence": [IGH_SEQ]})
        with pytest.raises(polars.exceptions.ComputeError):
            df.select(imp.number(polars.col("sequence"), chains=["IGX"], scheme="IMGT"))


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
            "fr1",
            "fr2",
            "fr3",
            "fr4",
            "cdr1",
            "cdr2",
            "cdr3",
            "prefix",
            "postfix",
            "error",
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


# A signal peptide before the domain and a tag after it: residues the aligner leaves out.
FLANKED_SEQ = "MGWSCIILFLVATATGVHSX" + IGH_SEQ + "HHHHHHEPEA"
REGIONS = ("prefix", "fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4", "postfix")


def as_numbering_result(row: dict) -> dict:
    """A Polars numbering row in the shape of `Annotator.number`'s result."""
    numbering = row["numbering"]
    if numbering is not None:
        numbering = {r["position"]: r["residue"] for r in numbering}
    return {**row, "numbering": numbering}


class TestPolarsMatchesAnnotator:
    """Issues #53 and #58: every Polars expression returns what `Annotator` returns for the same
    sequence, flanking residues included, with the same fields."""

    SEQUENCES = [FLANKED_SEQ, "A" * 40]

    @pytest.fixture
    def annotator(self):
        from immunum import Annotator

        return Annotator(["IGH"], "IMGT")

    def assert_matches_annotator(self, expr, annotator):
        from dataclasses import asdict, fields
        from immunum import NumberingResult

        df = polars.DataFrame({"sequence": self.SEQUENCES})
        rows = df.select(expr.alias("n")).unnest("n").to_dicts()
        assert list(rows[0]) == [f.name for f in fields(NumberingResult)]
        assert [as_numbering_result(row) for row in rows] == [
            asdict(annotator.number(s)) for s in self.SEQUENCES
        ]

    def test_number(self, annotator):
        expr = imp.number(polars.col("sequence"), chains=["IGH"], scheme="IMGT")
        self.assert_matches_annotator(expr, annotator)

    def test_numbering_method(self, annotator):
        expr = imp.numbering_method(polars.col("sequence"), annotator=annotator)
        self.assert_matches_annotator(expr, annotator)

    def test_segment(self, annotator):
        df = polars.DataFrame({"sequence": [FLANKED_SEQ]})
        row = df.select(
            imp.segment(polars.col("sequence"), chains=["IGH"], scheme="IMGT").alias(
                "s"
            )
        ).unnest("s")
        expected = annotator.segment(FLANKED_SEQ)
        assert {r: row[r][0] for r in REGIONS} == expected.as_dict()
        assert "".join(row[r][0] for r in REGIONS) == FLANKED_SEQ

    def test_segmentation_method(self, annotator):
        df = polars.DataFrame({"sequence": [FLANKED_SEQ]})
        row = df.select(
            imp.segmentation_method(polars.col("sequence"), annotator=annotator).alias(
                "s"
            )
        ).unnest("s")
        expected = annotator.segment(FLANKED_SEQ)
        assert {r: row[r][0] for r in REGIONS} == expected.as_dict()
        assert "".join(row[r][0] for r in REGIONS) == FLANKED_SEQ


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
