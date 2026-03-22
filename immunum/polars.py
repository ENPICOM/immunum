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
from immunum import _normalize_chains, _normalize_scheme

if TYPE_CHECKING:
    from immunum.typing import IntoExprColumn


LIB = Path(__file__).parent


def numbering_method(expr: IntoExprColumn, *, annotator: _Annotator) -> pl.Expr:
    return register_plugin_function(
        args=[expr],
        plugin_path=LIB,
        function_name="numbering_class_struct_expr",
        is_elementwise=True,
        kwargs={"annotator": annotator},
    )


def number(
    expr: IntoExprColumn,
    *,
    chains: list[str],
    scheme: str,
    min_confidence: float | None = None,
) -> pl.Expr:
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
