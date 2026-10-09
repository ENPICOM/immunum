from __future__ import annotations

import warnings
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


def _plugin(
    expr: IntoExprColumn,
    function: str,
    chains: list[str] | None,
    scheme: str | None,
    min_confidence: float | None,
    annotator: Annotator | None,
) -> pl.Expr:
    """The plugin `function` with `annotator`, or with one built from the names now.

    Building it here raises `immunum.Error` for a bad name when the expression is built, not
    when the query runs.
    """
    if annotator is not None:
        if chains is not None or scheme is not None or min_confidence is not None:
            raise TypeError(
                "pass either `annotator` or `chains`, `scheme` and `min_confidence`, not both"
            )
    elif chains is None or scheme is None:
        raise TypeError("pass `chains` and `scheme`, or a prebuilt `annotator`")
    else:
        annotator = Annotator(chains, scheme, min_confidence)
    return register_plugin_function(
        args=[expr],
        plugin_path=LIB,
        function_name=f"{function}_class_struct_expr",
        is_elementwise=True,
        kwargs={"annotator": annotator._annotator},
    )


def number(
    expr: IntoExprColumn,
    *,
    chains: list[str] | None = None,
    scheme: str | None = None,
    min_confidence: float | None = None,
    annotator: Annotator | None = None,
) -> pl.Expr:
    """Number sequences as a Polars expression.

    Each row gets the fields `Annotator.number` returns: `chain`, `scheme`, `confidence`,
    `numbering` (a list of `{position, residue}` structs), `query_start`, `query_end`, `error`
    and `error_kind`. On failure, `error` and `error_kind` are set and every other field is null.

    Pass either `chains` and `scheme` (and optionally `min_confidence`), or a prebuilt
    `annotator`; both give the same result and run as fast. Either way the names are checked
    when the expression is built, so a bad one raises `immunum.Error` straight away.

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
                "error_kind": pl.String,
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

    # shape: (2, 8)
    # ┌───────┬────────┬────────────┬─────────────────────────────────┬─────────────┬───────────┬───────┬────────────┐
    # │ chain ┆ scheme ┆ confidence ┆ numbering                       ┆ query_start ┆ query_end ┆ error ┆ error_kind │
    # │ ---   ┆ ---    ┆ ---        ┆ ---                             ┆ ---         ┆ ---       ┆ ---   ┆ ---        │
    # │ str   ┆ str    ┆ f32        ┆ list[struct[2]]                 ┆ u32         ┆ u32       ┆ str   ┆ str        │
    # ╞═══════╪════════╪════════════╪═════════════════════════════════╪═════════════╪═══════════╪═══════╪════════════╡
    # │ H     ┆ IMGT   ┆ 0.784515   ┆ [{"1","Q"}, {"2","V"}, … {"128… ┆ 0           ┆ 121       ┆ null  ┆ null       │
    # │ K     ┆ IMGT   ┆ 0.878814   ┆ [{"1","D"}, {"2","I"}, … {"127… ┆ 0           ┆ 106       ┆ null  ┆ null       │
    # └───────┴────────┴────────────┴─────────────────────────────────┴─────────────┴───────────┴───────┴────────────┘

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
        chains (list[str] | None): chains to consider, as for `Annotator`. Required unless
            `annotator` is given.
        scheme (str | None): numbering scheme, as for `Annotator`. Required unless
            `annotator` is given.
        min_confidence (float | None, optional): minimum alignment confidence, as for
            `Annotator`. Defaults to None (corresponds to 0.5).
        annotator (Annotator | None, optional): a prebuilt `Annotator` to use instead of
            `chains`, `scheme` and `min_confidence`.

    Returns:
        pl.Expr: numbering expression

    Raises:
        immunum.Error: If `chains`, `scheme` or `min_confidence` is wrong, as for `Annotator`.
        TypeError: If both or neither of `annotator` and the names are given.
    """
    return _plugin(expr, "numbering", chains, scheme, min_confidence, annotator)


def number_domains(
    expr: IntoExprColumn,
    *,
    chains: list[str] | None = None,
    scheme: str | None = None,
    min_confidence: float | None = None,
    annotator: Annotator | None = None,
) -> pl.Expr:
    """Number every variable domain in each sequence, such as both domains of an scFv.

    Each row gets a list with one struct per domain, ordered by position, each with the fields
    `number` returns for that domain. The list is empty when no domain aligns with enough
    confidence. When the sequence itself is invalid (too short, too long or not amino acids), it
    holds a single struct with `error` and `error_kind` set and every other field null.

    A domain that lacks its first IMGT positions (a light chain starting at position 2, say) and
    directly follows other residues, such as a linker, can have the residue just before it
    numbered as its first position: IMGT position 1 is so variable that the sequence alone can't
    tell a linker residue from the domain's own first residue.

    Pass either `chains` and `scheme`, or a prebuilt `annotator`, as for `number`.

    Example:

    ```python
    import polars as pl
    import immunum.polars as imp

    heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
    kappa = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"

    df = pl.DataFrame(
        {
            "sequence": [
                heavy
                + "GGGGSGGGGSGGGGS"
                + kappa,
                heavy,
            ]
        }
    ).select(
        imp.number_domains(
            "sequence",
            chains=["ig"],
            scheme="imgt",
        ).alias("domains")
    )

    # One row per domain
    print(
        df.with_row_index("sequence")
        .explode("domains")
        .unnest("domains")
        .select(
            "sequence",
            "chain",
            "confidence",
            "query_start",
            "query_end",
        )
    )
    ```

    Args:
        expr (IntoExprColumn): input polars expression (e.g. `pl.col('sequence')`)
        chains (list[str] | None): chains to consider, as for `Annotator`. Required unless
            `annotator` is given.
        scheme (str | None): numbering scheme, as for `Annotator`. Required unless
            `annotator` is given.
        min_confidence (float | None, optional): minimum alignment confidence, as for
            `Annotator`. Defaults to None (corresponds to 0.5).
        annotator (Annotator | None, optional): a prebuilt `Annotator` to use instead of
            `chains`, `scheme` and `min_confidence`.

    Returns:
        pl.Expr: domains expression

    Raises:
        immunum.Error: If `chains`, `scheme` or `min_confidence` is wrong, as for `Annotator`.
        TypeError: If both or neither of `annotator` and the names are given.
    """
    return _plugin(expr, "number_domains", chains, scheme, min_confidence, annotator)


def segment(
    expr: IntoExprColumn,
    *,
    chains: list[str] | None = None,
    scheme: str | None = None,
    min_confidence: float | None = None,
    annotator: Annotator | None = None,
) -> pl.Expr:
    """Split sequences into FR/CDR regions as a Polars expression.

    Each row gets `prefix`, `fr1`, `cdr1`, `fr2`, `cdr2`, `fr3`, `cdr3`, `fr4`, `postfix`,
    `error` and `error_kind`. The segments join back into the sequence. On failure, `error` and
    `error_kind` are set and every segment is null.

    Pass either `chains` and `scheme`, or a prebuilt `annotator`, as for `number`.

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
            "error_kind": pl.String,
        }
    )
    print(
        df.select(
            pl.col(
                "segmentation"
            ).struct.unnest()
        )
    )
    ```

    Args:
        expr (IntoExprColumn): input polars expression (e.g. `pl.col('sequence')`)
        chains (list[str] | None): chains to consider, as for `Annotator`. Required unless
            `annotator` is given.
        scheme (str | None): numbering scheme, as for `Annotator`. Required unless
            `annotator` is given.
        min_confidence (float | None, optional): minimum alignment confidence, as for
            `Annotator`. Defaults to None (corresponds to 0.5).
        annotator (Annotator | None, optional): a prebuilt `Annotator` to use instead of
            `chains`, `scheme` and `min_confidence`.

    Returns:
        pl.Expr: segmentation expression

    Raises:
        immunum.Error: If `chains`, `scheme` or `min_confidence` is wrong, as for `Annotator`.
        TypeError: If both or neither of `annotator` and the names are given.
    """
    return _plugin(expr, "segmentation", chains, scheme, min_confidence, annotator)


def segment_domains(
    expr: IntoExprColumn,
    *,
    chains: list[str] | None = None,
    scheme: str | None = None,
    min_confidence: float | None = None,
    annotator: Annotator | None = None,
) -> pl.Expr:
    """Split every variable domain in each sequence into FR/CDR regions.

    Each row gets a list with one struct per domain, ordered by position, each with the fields
    `segment` returns. Every residue lands in exactly one domain's regions: a domain's `prefix`
    holds the residues since the previous domain (or the start of the sequence), and only the
    last domain has the residues after it as its `postfix`, so all domains' regions in order
    rebuild the sequence. The list is empty when no domain aligns with enough confidence; for an
    invalid sequence it holds a single struct with `error` and `error_kind` set and every segment
    null.

    Pass either `chains` and `scheme`, or a prebuilt `annotator`, as for `number`.

    Example:

    ```python
    import polars as pl
    import immunum.polars as imp

    heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
    kappa = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"

    df = pl.DataFrame(
        {
            "sequence": [
                heavy
                + "GGGGSGGGGSGGGGS"
                + kappa
            ]
        }
    ).select(
        imp.segment_domains(
            "sequence",
            chains=["ig"],
            scheme="imgt",
        ).alias("domains")
    )

    # One row per domain
    print(
        df.explode("domains")
        .unnest("domains")
        .select(
            "prefix",
            "cdr1",
            "cdr2",
            "cdr3",
        )
    )
    ```

    Args:
        expr (IntoExprColumn): input polars expression (e.g. `pl.col('sequence')`)
        chains (list[str] | None): chains to consider, as for `Annotator`. Required unless
            `annotator` is given.
        scheme (str | None): numbering scheme, as for `Annotator`. Required unless
            `annotator` is given.
        min_confidence (float | None, optional): minimum alignment confidence, as for
            `Annotator`. Defaults to None (corresponds to 0.5).
        annotator (Annotator | None, optional): a prebuilt `Annotator` to use instead of
            `chains`, `scheme` and `min_confidence`.

    Returns:
        pl.Expr: domains expression

    Raises:
        immunum.Error: If `chains`, `scheme` or `min_confidence` is wrong, as for `Annotator`.
        TypeError: If both or neither of `annotator` and the names are given.
    """
    return _plugin(expr, "segment_domains", chains, scheme, min_confidence, annotator)


def numbering_method(expr: IntoExprColumn, *, annotator: Annotator) -> pl.Expr:
    """Deprecated: use `number(expr, annotator=annotator)`, which returns the same."""
    warnings.warn(
        "numbering_method is deprecated; use number(expr, annotator=annotator)",
        DeprecationWarning,
        stacklevel=2,
    )
    return number(expr, annotator=annotator)


def segmentation_method(expr: IntoExprColumn, *, annotator: Annotator) -> pl.Expr:
    """Deprecated: use `segment(expr, annotator=annotator)`, which returns the same."""
    warnings.warn(
        "segmentation_method is deprecated; use segment(expr, annotator=annotator)",
        DeprecationWarning,
        stacklevel=2,
    )
    return segment(expr, annotator=annotator)
