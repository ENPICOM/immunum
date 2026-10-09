from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

try:
    import polars as pl
    from polars.plugins import register_plugin_function
except ImportError as e:
    raise ImportError(
        "polars is required to use immunum.polars. Install it with: pip install polars"
    ) from e

from immunum._internal import _Annotator  # noqa: F401
from immunum import Annotator

if TYPE_CHECKING:
    from immunum.typing import IntoExprColumn


LIB = Path(__file__).parent


def number(
    expr: IntoExprColumn,
    *,
    chains: list[str],
    scheme: str,
    min_confidence: float | None = None,
) -> pl.Expr:
    """Number sequences as a Polars expression.

    Each row gets the fields `Annotator.number` returns: `chain`, `scheme`, `confidence`,
    `numbering` (a list of `{position, residue}` structs), `query_start`, `query_end` and
    `error`. On failure, `error` is set and every other field is null.

    The annotator is built from `chains`, `scheme` and `min_confidence` when the query runs, so
    an unknown name or an out-of-range `min_confidence` raises a `ComputeError` then.
    `numbering_method` takes a prebuilt `Annotator` instead and returns the same fields.

    Example:

    ```python
    import polars as pl
    import immunum.polars as imp

    df = pl.DataFrame(
        {
            "sequence": [
                "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS",
                "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
            ]
        }
    ).select(
        imp.number(
            "sequence",
            chains=["h", "k"],
            scheme="imgt",
        ).alias("numbering")
    )
    assert df.dtypes == [
        pl.Struct(
            {
                "chain": pl.String,
                "scheme": pl.String,
                "confidence": pl.Float32,
                "numbering": pl.List(
                    pl.Struct(
                        {
                            "position": pl.String,
                            "residue": pl.String,
                        }
                    )
                ),
                "query_start": pl.UInt32,
                "query_end": pl.UInt32,
                "error": pl.String,
            }
        )
    ]
    print(
        df.select(
            pl.col(
                "numbering"
            ).struct.unnest()
        )
    )

    # shape: (2, 7)
    # ┌───────┬────────┬────────────┬─────────────────────────────────┬─────────────┬───────────┬───────┐
    # │ chain ┆ scheme ┆ confidence ┆ numbering                       ┆ query_start ┆ query_end ┆ error │
    # │ ---   ┆ ---    ┆ ---        ┆ ---                             ┆ ---         ┆ ---       ┆ ---   │
    # │ str   ┆ str    ┆ f32        ┆ list[struct[2]]                 ┆ u32         ┆ u32       ┆ str   │
    # ╞═══════╪════════╪════════════╪═════════════════════════════════╪═════════════╪═══════════╪═══════╡
    # │ H     ┆ IMGT   ┆ 0.784515   ┆ [{"1","Q"}, {"2","V"}, … {"128… ┆ 0           ┆ 121       ┆ null  │
    # │ K     ┆ IMGT   ┆ 0.878814   ┆ [{"1","D"}, {"2","I"}, … {"127… ┆ 0           ┆ 106       ┆ null  │
    # └───────┴────────┴────────────┴─────────────────────────────────┴─────────────┴───────────┴───────┘

    # One row per numbered residue
    print(
        df.select(
            pl.col(
                "numbering"
            ).struct.field(
                "chain",
                "numbering",
            )
        )
        .explode("numbering")
        .unnest("numbering")
        .head(3)
    )

    # shape: (3, 3)
    # ┌───────┬──────────┬─────────┐
    # │ chain ┆ position ┆ residue │
    # │ ---   ┆ ---      ┆ ---     │
    # │ str   ┆ str      ┆ str     │
    # ╞═══════╪══════════╪═════════╡
    # │ H     ┆ 1        ┆ Q       │
    # │ H     ┆ 2        ┆ V       │
    # │ H     ┆ 3        ┆ Q       │
    # └───────┴──────────┴─────────┘
    ```

    Args:
        expr (IntoExprColumn): input polars expression (e.g. `pl.col('sequence')`)
        chains (list[str]): list of chains to use for initialized `Annotator`
        scheme (str): scheme to use for initialized `Annotator`
        min_confidence (float | None, optional): confidence to use for initialized `Annotator`. Defaults to None (corresponds to 0.5)

    Returns:
        pl.Expr: numbering expression
    """
    return register_plugin_function(
        args=[expr],
        plugin_path=LIB,
        function_name="numbering_struct_expr",
        is_elementwise=True,
        kwargs={
            "chains": chains,
            "scheme": scheme,
            "min_confidence": min_confidence,
        },
    )


def segment(
    expr: IntoExprColumn,
    *,
    chains: list[str],
    scheme: str,
    min_confidence: float | None = None,
) -> pl.Expr:
    """Split sequences into FR/CDR regions as a Polars expression.

    Each row gets `prefix`, `fr1`, `cdr1`, `fr2`, `cdr2`, `fr3`, `cdr3`, `fr4`, `postfix` and
    `error`. The segments join back into the sequence. On failure, `error` is set and every
    segment is null.

    The annotator is built from `chains`, `scheme` and `min_confidence` when the query runs, so
    an unknown name or an out-of-range `min_confidence` raises a `ComputeError` then.
    `segmentation_method` takes a prebuilt `Annotator` instead and returns the same fields.

    Example:

    ```python
    import polars as pl
    import immunum.polars as imp

    df = pl.DataFrame(
        {
            "sequence": [
                "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS",
                "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
            ]
        }
    ).select(
        imp.segment(
            "sequence",
            chains=["h", "k", "l"],
            scheme="imgt",
            min_confidence=0.0,
        ).alias("segmentation")
    )
    assert df[
        "segmentation"
    ].dtype == pl.Struct(
        {
            "prefix": pl.String,
            "fr1": pl.String,
            "cdr1": pl.String,
            "fr2": pl.String,
            "cdr2": pl.String,
            "fr3": pl.String,
            "cdr3": pl.String,
            "fr4": pl.String,
            "postfix": pl.String,
            "error": pl.String,
        }
    )
    print(
        df.select(
            pl.col(
                "segmentation"
            ).struct.unnest()
        )
    )

    # shape: (2, 9)
    # ┌────────────────────────┬──────────┬───┬─────────────┬────────┬─────────┐
    # │ fr1                    ┆ cdr1     ┆ … ┆ fr4         ┆ prefix ┆ postfix │
    # │ ---                    ┆ ---      ┆   ┆ ---         ┆ ---    ┆ ---     │
    # │ str                    ┆ str      ┆   ┆ str         ┆ str    ┆ str     │
    # ╞════════════════════════╪══════════╪═══╪═════════════╪════════╪═════════╡
    # │ QVQLVQSGAEVKRPGSSVTVS… ┆ GGSFSTYA ┆ … ┆ WGQGTLVTVSS ┆        ┆         │
    # │ DIQMTQSPSSLSASVGDRVTI… ┆ RASQDVNT ┆ … ┆ FGQGTKVEIK  ┆        ┆         │
    # └────────────────────────┴──────────┴───┴─────────────┴────────┴─────────┘
    ```

    Args:
        expr (IntoExprColumn): input polars expression (e.g. `pl.col('sequence')`)
        chains (list[str]): list of chains to use for initialized `Annotator`
        scheme (str): scheme to use for initialized `Annotator`
        min_confidence (float | None, optional): confidence to use for initialized `Annotator`. Defaults to None (corresponds to 0.5)

    Returns:
        pl.Expr: segmentation expression
    """
    return register_plugin_function(
        args=[expr],
        plugin_path=LIB,
        function_name="segmentation_struct_expr",
        is_elementwise=True,
        kwargs={
            "chains": chains,
            "scheme": scheme,
            "min_confidence": min_confidence,
        },
    )


def numbering_method(expr: IntoExprColumn, *, annotator: Annotator) -> pl.Expr:
    """Number sequences with a prebuilt `Annotator`.

    Returns exactly what `number` returns, and runs as fast. Use it when your code already holds
    an `Annotator`: its chains, scheme and `min_confidence` were checked when it was built, so a
    mistake raises `ValueError` there instead of when the query runs. The annotator travels with
    the query and is rebuilt from it on every call, so it saves no set-up work over `number`.

    Example:

    ```python
    import polars as pl
    import immunum
    import immunum.polars as imp

    annotator = immunum.Annotator(
        chains=["H", "K", "L"],
        scheme="imgt",
    )

    df = pl.DataFrame(
        {
            "sequence": [
                "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS",
            ]
        }
    ).select(
        imp.numbering_method(
            pl.col("sequence"),
            annotator=annotator,
        ).alias("numbering")
    )
    assert df[
        "numbering"
    ].struct.fields == [
        "chain",
        "scheme",
        "confidence",
        "numbering",
        "query_start",
        "query_end",
        "error",
    ]
    ```

    Args:
        expr (IntoExprColumn): input polars expression (e.g. `pl.col('sequence')`)
        annotator (Annotator): pre-built `Annotator` instance

    Returns:
        pl.Expr: numbering expression
    """
    return register_plugin_function(
        args=[expr],
        plugin_path=LIB,
        function_name="numbering_class_struct_expr",
        is_elementwise=True,
        kwargs={"annotator": annotator._annotator},
    )


def segmentation_method(expr: IntoExprColumn, *, annotator: Annotator) -> pl.Expr:
    """Segment sequences with a prebuilt `Annotator`.

    Returns exactly what `segment` returns, and runs as fast. Use it when your code already holds
    an `Annotator`: its chains, scheme and `min_confidence` were checked when it was built, so a
    mistake raises `ValueError` there instead of when the query runs. The annotator travels with
    the query and is rebuilt from it on every call, so it saves no set-up work over `segment`.

    Example:

    ```python
    import polars as pl
    import immunum
    import immunum.polars as imp

    annotator = immunum.Annotator(
        chains=["H", "K", "L"],
        scheme="imgt",
    )

    df = pl.DataFrame(
        {
            "sequence": [
                "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS",
                "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
            ]
        }
    ).select(
        imp.segmentation_method(
            "sequence", annotator=annotator
        ).alias("segmentation")
    )
    assert df[
        "segmentation"
    ].dtype == pl.Struct(
        {
            "prefix": pl.String,
            "fr1": pl.String,
            "cdr1": pl.String,
            "fr2": pl.String,
            "cdr2": pl.String,
            "fr3": pl.String,
            "cdr3": pl.String,
            "fr4": pl.String,
            "postfix": pl.String,
            "error": pl.String,
        }
    )
    ```

    Args:
        expr (IntoExprColumn): input polars expression (e.g. `pl.col('sequence')`)
        annotator (Annotator): pre-built `Annotator` instance

    Returns:
        pl.Expr: segmentation expression
    """
    return register_plugin_function(
        args=[expr],
        plugin_path=LIB,
        function_name="segmentation_class_struct_expr",
        is_elementwise=True,
        kwargs={"annotator": annotator._annotator},
    )
