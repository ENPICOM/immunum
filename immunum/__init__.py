from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING  # noqa: F401

import polars as pl
from polars.plugins import register_plugin_function

from immunum._internal import Annotator  # noqa: F401

if TYPE_CHECKING:
    from immunum.typing import IntoExprColumn


LIB = Path(__file__).parent


def numbering_method(expr: IntoExprColumn, *, annotator: Annotator) -> pl.Expr:
    return register_plugin_function(
        args=[expr],
        plugin_path=LIB,
        function_name="numbering_class_struct_expr",
        is_elementwise=True,
        kwargs={"annotator": annotator},
    )


def number(expr: IntoExprColumn, *, chains: list[str], scheme: str) -> pl.Expr:
    return register_plugin_function(
        args=[expr],
        plugin_path=LIB,
        function_name="numbering_struct_expr",
        is_elementwise=True,
        kwargs={"chains": chains, "scheme": scheme},
    )
