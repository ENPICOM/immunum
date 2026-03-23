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
from immunum import _normalize_chains, _normalize_scheme, Annotator

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
    """Segment a polars expr with immunum.

    Annotator object will be initialized with chains, scheme, min_confidence.
    Results are returned as `Struct({'chain': String, 'scheme': String, 'positions': List(String), 'residues': List(String)})`

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
            chains=["h"],
            scheme="imgt",
            min_confidence=0.0,
        ).alias("numbering")
    )
    assert df.dtypes == [
        pl.Struct(
            {
                "chain": pl.String,
                "scheme": pl.String,
                "positions": pl.List(
                    pl.String
                ),
                "residues": pl.List(
                    pl.String
                ),
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

    # shape: (2, 4)
    # ┌───────┬────────┬─────────────────────┬───────────────────┐
    # │ chain ┆ scheme ┆ positions           ┆ residues          │
    # │ ---   ┆ ---    ┆ ---                 ┆ ---               │
    # │ str   ┆ str    ┆ list[str]           ┆ list[str]         │
    # ╞═══════╪════════╪═════════════════════╪═══════════════════╡
    # │ H     ┆ IMGT   ┆ ["1", "2", … "128"] ┆ ["Q", "V", … "S"] │
    # │ H     ┆ IMGT   ┆ ["1", "2", … "127"] ┆ ["D", "I", … "K"] │
    # └───────┴────────┴─────────────────────┴───────────────────┘
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
            "chains": _normalize_chains(chains),
            "scheme": _normalize_scheme(scheme),
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

    Annotator object will be initialized with chains, scheme, min_confidence.
    Results are returned as `Struct({'fr1': String, 'cdr1': String, 'fr2': String, 'cdr2': String, 'fr3': String, 'cdr3': String, 'fr4': String, 'prefix': String, 'postfix': String})`.

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
            "chains": _normalize_chains(chains),
            "scheme": _normalize_scheme(scheme),
            "min_confidence": min_confidence,
        },
    )


def numbering_method(expr: IntoExprColumn, *, annotator: Annotator) -> pl.Expr:
    """Number sequences using a pre-built `Annotator` instance.

    Prefer this over `number` when reusing the same annotator across multiple expressions,
    as it avoids re-initializing the annotator on each call. Returns the same struct shape
    as `number`.

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
    assert set(
        df["numbering"].struct.fields
    ) == {
        "chain",
        "scheme",
        "confidence",
        "numbering",
    }
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
    """Segment sequences using a pre-built `Annotator` instance.

    Prefer this over `segment` when reusing the same annotator across multiple expressions,
    as it avoids re-initializing the annotator on each call. Returns the same struct shape
    as `segment`.

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
