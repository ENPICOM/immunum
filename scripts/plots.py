# /// script
# requires-python = "==3.13"
# dependencies = [
#   "altair>=6.0.0",
#   "polars>1.30.0",
#   "pyarrow>=23.0.1",
#   "pandas<3.0.0",
#   "vl-convert-python>=1.9.0.post1",
# ]
# ///

import glob
import math

import polars as pl
import altair as alt


# ── Load data ────────────────────────────────────────────────────────────────
accuracy_files = glob.glob("resources/benchmark_results/results_*_imgt.csv")
acc_df = pl.concat([pl.read_csv(f) for f in accuracy_files])

speed_df = pl.read_csv("resources/benchmark_results/results_speed.csv")
headline_df = pl.read_csv("resources/benchmark_results/results_headline.csv")

# Derive available values from data
acc_sizes = sorted(acc_df["sample_size"].unique().to_list())
acc_size_label = ", ".join(f"n={s:,}" for s in acc_sizes)

segments = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]

# ── Plots 1 & 3: Performance at a fixed batch size ───────────────────────────
single_tools = [
    "immunum_singlethreaded",
    "anarci",
    "anarcii2",
    "antpack",
]
single_label_map = {
    "immunum_singlethreaded": "immunum",
    "anarci": "ANARCI",
    "anarcii2": "ANARCII",
    "antpack": "AntPack",
}
mt_tools = [
    "immunum_multithreaded",
    "antpack_parallel",
    "anarci_parallel",
    "anarcii2_parallel",
]
mt_label_map = {
    "immunum_multithreaded": "immunum",
    "antpack_parallel": "AntPack",
    "anarci_parallel": "ANARCI",
    "anarcii2_parallel": "ANARCII",
}

# The headline run: one batch of NGS antibody sequences, H, K and L mixed
HEADLINE_SIZE = headline_df["sample_size"].max()
# The benchmark machines' core count, which the results files don't record; scripts/mkdocs_pages.py
# has the same constant for the docs tables
CORES = 32
SINGLE_CORE = "1 core"
ALL_CORES = f"{CORES} cores"
perf_base = headline_df.with_columns(
    (pl.col("sample_size") / pl.col("time_s")).alias("throughput")
)

p1_data = perf_base.filter(pl.col("tool").is_in(single_tools)).with_columns(
    pl.col("tool").replace(single_label_map).alias("method")
)
p3_data = perf_base.filter(pl.col("tool").is_in(mt_tools)).with_columns(
    pl.col("tool").replace(mt_label_map).alias("method")
)


def make_perf_barplot(data: pl.DataFrame, title: str) -> alt.Chart:
    summary = (
        data.group_by("method")
        .agg(
            pl.col("throughput").mean().alias("mean_throughput"),
            pl.col("throughput").std().alias("std_throughput"),
        )
        .sort("mean_throughput", descending=True)
    )
    method_order = summary["method"].to_list()
    bars = (
        alt.Chart()
        .mark_bar(opacity=0.7, stroke="black", strokeWidth=1.5)
        .encode(
            x=alt.X("method:N", title="Method", sort=method_order),
            y=alt.Y("mean_throughput:Q", title="Sequences per second"),
            color=alt.Color("method:N", legend=None),
            tooltip=[
                "method",
                alt.Tooltip("mean_throughput:Q", format=",.0f", title="Mean seq/s"),
                alt.Tooltip("std_throughput:Q", format=",.0f", title="Std seq/s"),
            ],
        )
    )
    errorbars = (
        alt.Chart()
        .mark_rule(strokeWidth=2)
        .encode(
            x=alt.X("method:N", sort=method_order),
            y=alt.Y("errbar_min:Q"),
            y2=alt.Y2("errbar_max:Q"),
        )
        .transform_calculate(
            errbar_min="datum.mean_throughput - datum.std_throughput / 2",
            errbar_max="datum.mean_throughput + datum.std_throughput / 2",
        )
    )
    return alt.layer(bars, errorbars, data=summary).properties(
        title=title, width=300, height=350
    )


plot_perf = (
    make_perf_barplot(p1_data, f"{SINGLE_CORE} (n={HEADLINE_SIZE:,})")
    | make_perf_barplot(p3_data, f"{ALL_CORES} (n={HEADLINE_SIZE:,})")
).resolve_scale(y="shared")


# ── Plot 3: Scaling — time vs batch size ──────────────────────────────────────
def make_scaling_plot(
    data: pl.DataFrame, tool_label_map: dict, title: str
) -> alt.Chart:
    df = (
        data.with_columns(pl.col("tool").replace(tool_label_map).alias("method"))
        .group_by(["method", "sample_size"])
        .agg(
            pl.col("time_s").mean().alias("mean_time"),
            pl.col("time_s").std().alias("std_time"),
        )
        .sort(["method", "sample_size"])
    )
    size_order = sorted(df["sample_size"].unique().to_list())
    pdf = df.to_pandas()
    pdf["sample_size"] = pdf["sample_size"].astype(str)
    size_order_str = [str(s) for s in size_order]

    lines = (
        alt.Chart()
        .mark_line(point=True)
        .encode(
            x=alt.X(
                "sample_size:O",
                title="Batch size",
                sort=size_order_str,
                axis=alt.Axis(labelAngle=45),
            ),
            y=alt.Y(
                "mean_time:Q",
                title="Time (s)",
                scale=alt.Scale(type="log", base=10),
            ),
            color=alt.Color("method:N", title="Method"),
            tooltip=[
                "method",
                alt.Tooltip("sample_size:O", title="Batch size"),
                alt.Tooltip("mean_time:Q", format=".3f", title="Mean time (s)"),
                alt.Tooltip("std_time:Q", format=".3f", title="Std (s)"),
            ],
        )
    )
    errorbars = (
        alt.Chart()
        .mark_rule()
        .encode(
            x=alt.X("sample_size:O", sort=size_order_str),
            y=alt.Y("errbar_min:Q"),
            y2=alt.Y2("errbar_max:Q"),
            color=alt.Color("method:N"),
        )
        .transform_calculate(
            errbar_min="datum.mean_time - datum.std_time / 2",
            errbar_max="datum.mean_time + datum.std_time / 2",
        )
    )
    return alt.layer(lines, errorbars, data=pdf).properties(
        title=title, width=400, height=350
    )


plot_scaling = (
    make_scaling_plot(
        speed_df.filter(pl.col("tool").is_in(single_tools)),
        single_label_map,
        SINGLE_CORE,
    )
    & make_scaling_plot(
        speed_df.filter(pl.col("tool").is_in(mt_tools)),
        mt_label_map,
        ALL_CORES,
    )
).resolve_scale(y="shared", color="shared")

# ── Plot 2: Correctness by chain and segment ──────────────────────────────────
p2_data = (
    acc_df.filter(
        ~pl.col("tool").str.ends_with("_parallel")
        & ~pl.col("tool").str.ends_with("_multithreaded")
    )
    .with_columns(pl.col("tool").replace(single_label_map).alias("method"))
    .drop("tool")
    .unpivot(
        index=["method", "chain", "round", "sample_size"],
        on=segments,
        variable_name="segment",
        value_name="pct_correct",
    )
)

CHAIN_ORDER = ["H", "K", "L", "A", "B", "G", "D"]
# Full chain names as row labels, so the row header title is free to carry the one shared y title
CHAIN_LABEL_EXPR = (
    "{'H': 'IGH', 'K': 'IGK', 'L': 'IGL', 'A': 'TRA', 'B': 'TRB', 'G': 'TRG', 'D': 'TRD'}"
    "[datum.value]"
)
METHOD_ORDER = ["immunum", "AntPack", "ANARCI", "ANARCII"]
# The colours Vega's default scheme gives the tools in plots 1 and 3, pinned so this plot's legend
# and bars can follow METHOD_ORDER instead of the alphabetical order that colour assignment uses
METHOD_COLORS = ["#72b7b2", "#e45756", "#4c78a8", "#f58518"]

p2_summary = (
    p2_data.group_by("method", "chain", "segment")
    .agg(
        pl.col("pct_correct").mean().alias("mean_pct"),
        pl.col("pct_correct").std().fill_null(0).alias("std_pct"),
    )
    .with_columns(
        (pl.col("mean_pct") - pl.col("std_pct")).alias("err_lo"),
        # A standard deviation above a mean near 100% would pass 100%
        pl.min_horizontal(pl.col("mean_pct") + pl.col("std_pct"), 100).alias("err_hi"),
    )
)
# Zoomed in on the top of the axis, where every mean sits: on 0–100% a 97% bar is
# indistinguishable from a 100% one. The floor is the lowest error bar's end, rounded down
PCT_FLOOR = math.floor(p2_summary["err_lo"].min())
PCT_SCALE = alt.Scale(domain=[PCT_FLOOR, 100], nice=False)

_p2_base = alt.Chart()
# Outlined like plot 1's bars. Clipped: with the axis starting above 0 a bar runs down past the frame
_p2_bars = _p2_base.mark_bar(
    opacity=0.7, stroke="black", strokeWidth=1.5, clip=True
).encode(
    x=alt.X(
        "segment:N",
        title="Region",
        sort=segments,
        axis=alt.Axis(labelAngle=0),
        scale=alt.Scale(paddingInner=0.25),
    ),
    xOffset=alt.XOffset("method:N", sort=METHOD_ORDER),
    y=alt.Y(
        "mean_pct:Q",
        title=None,
        scale=PCT_SCALE,
        axis=alt.Axis(tickCount=4),
    ),
    color=alt.Color(
        "method:N",
        title="Method",
        scale=alt.Scale(domain=METHOD_ORDER, range=METHOD_COLORS),
    ),
    tooltip=[
        "method",
        "chain",
        "segment",
        alt.Tooltip("mean_pct:Q", format=".2f", title="Mean %"),
        alt.Tooltip("std_pct:Q", format=".2f", title="Std %"),
    ],
)
_p2_errorbars = _p2_base.mark_rule(color="black", strokeWidth=2, clip=True).encode(
    x=alt.X("segment:N", sort=segments),
    xOffset=alt.XOffset("method:N", sort=METHOD_ORDER),
    y=alt.Y("err_lo:Q", scale=PCT_SCALE),
    y2="err_hi:Q",
)

plot2 = (
    alt.layer(_p2_bars, _p2_errorbars, data=p2_summary)
    .properties(width=480, height=70)
    .facet(
        row=alt.Row(
            "chain:N",
            sort=CHAIN_ORDER,
            header=alt.Header(
                labelAngle=0,
                labelAlign="right",
                labelExpr=CHAIN_LABEL_EXPR,
                # Shared y title: one label left of all rows instead of one per row
                title="% Correct",
            ),
        ),
        spacing=14,
        title=f"Correctness by segment ({acc_size_label})",
    )
    # TCR rows only have immunum and ANARCI: their bars take the width of the tools missing there
    .resolve_scale(xOffset="independent")
)

# ── Save ─────────────────────────────────────────────────────────────────────
BLACK = "#000000"
WHITE = "#ffffff"


def _configure_for_export(chart):
    return (
        chart.configure(background=WHITE)
        .configure_title(color=BLACK, anchor="middle")
        .configure_legend(labelColor=BLACK, titleColor=BLACK)
        .configure_view(stroke=BLACK)
    )


# ── README figure: time for the headline batch, on 1 core and on all ─────────
def make_headline_chart(foreground: str, muted: str) -> alt.VConcatChart:
    """Horizontal bars of mean wall time per tool, fastest first, for a light or dark page."""

    def panel(data: pl.DataFrame, title: str) -> alt.LayerChart:
        summary = (
            data.group_by("method")
            .agg(pl.col("time_s").mean().alias("mean_time"))
            .sort("mean_time")
            .with_columns(
                pl.col("mean_time")
                .map_elements(
                    lambda t: f"{t:.2f}s" if t < 10 else f"{t:.0f}s",
                    return_dtype=pl.String,
                )
                .alias("label")
            )
        )
        order = summary["method"].to_list()
        base = alt.Chart(summary).encode(
            y=alt.Y(
                "method:N",
                sort=order,
                title=None,
                axis=alt.Axis(
                    labelFontSize=13,
                    labelFontWeight="bold",
                    labelPadding=8,
                    ticks=False,
                    domain=False,
                ),
            ),
            x=alt.X("mean_time:Q", title=None, axis=None),
        )
        bars = base.mark_bar(height=18, cornerRadiusEnd=3).encode(
            color=alt.condition(
                alt.datum.method == "immunum",
                alt.value("#e2711d"),
                alt.value(muted),
            )
        )
        labels = base.mark_text(
            align="left", dx=6, fontSize=12, color=foreground
        ).encode(text="label:N")
        return alt.layer(bars, labels).properties(
            width=420,
            height=110,
            title=alt.TitleParams(title, anchor="start", fontSize=13),
        )

    return (
        alt.vconcat(
            panel(p1_data, SINGLE_CORE),
            panel(p3_data, ALL_CORES),
            spacing=18,
        )
        .resolve_scale(x="independent", y="independent")
        .configure(background="transparent", font="Helvetica, Arial, sans-serif")
        .configure_axis(labelColor=foreground)
        .configure_title(color=foreground)
        .configure_view(stroke=None)
    )


for theme, foreground, muted in [
    ("light", "#1f2328", "#9aa4b2"),
    ("dark", "#e6edf3", "#6e7781"),
]:
    svg_path = f"docs/assets/benchmark_headline_{theme}.svg"
    make_headline_chart(foreground, muted).save(svg_path)
    print(f"Saved {svg_path}")


for name, p in [
    ("plot1_performance", plot_perf),
    ("plot2_correctness", _configure_for_export(plot2)),
    ("plot3_scaling", _configure_for_export(plot_scaling)),
]:
    svg_path = f"docs/assets/benchmark_{name}.svg"
    p.save(svg_path)
    print(f"Saved {svg_path}")
    html_path = f"docs/assets/benchmark_{name}.html"
    p.save(html_path)
    print(f"Saved {html_path}")
