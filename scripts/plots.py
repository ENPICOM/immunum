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

import polars as pl
import altair as alt


# ── Load data ────────────────────────────────────────────────────────────────
accuracy_files = glob.glob("resources/benchmark_results/results_ab_*_imgt.csv")
acc_df = pl.concat([pl.read_csv(f) for f in accuracy_files])

speed_df = pl.read_csv("resources/benchmark_results/results_speed.csv")

# Derive available values from data
acc_sizes = sorted(acc_df["sample_size"].unique().to_list())
acc_size_label = ", ".join(f"n={s:,}" for s in acc_sizes)

# Summarise accuracy data (mean per tool/chain/segment)
segments = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
acc_summary = (
    acc_df.group_by(["tool", "chain"])
    .agg([pl.col(s).mean() for s in segments])
    .sort(["chain", "tool"])
)

# ── Plots 1 & 3: Performance boxplots at fixed batch size ────────────────────
BATCH_SIZE = 10_000

single_tools = [
    "immunum_singlethreaded",
    "anarci",
    "anarcii2",
    "antpack",
]
single_label_map = {
    "immunum_singlethreaded": "immunum",
    "anarci": "anarci",
    "anarcii2": "anarcii2",
    "antpack": "antpack",
}
mt_tools = [
    "immunum_multithreaded",
    "antpack_parallel",
    "anarci_parallel",
    "anarcii2_parallel",
]
mt_label_map = {
    "immunum_multithreaded": "immunum",
    "antpack_parallel": "antpack",
    "anarci_parallel": "anarci",
    "anarcii2_parallel": "anarcii2",
}

METHODS = ["immunum", "anarci", "anarcii2", "antpack"]

perf_base = speed_df.filter(pl.col("sample_size") == BATCH_SIZE).with_columns(
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
    make_perf_barplot(p1_data, f"Single-threaded (n={int(BATCH_SIZE):,})")
    | make_perf_barplot(p3_data, f"Multi-threaded (n={int(BATCH_SIZE):,})")
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
        "Single-threaded",
    )
    & make_scaling_plot(
        speed_df.filter(pl.col("tool").is_in(mt_tools)),
        mt_label_map,
        "Multi-threaded",
    )
).resolve_scale(y="shared", color="shared")

# ── Plot 2: Correctness by chain and segment ──────────────────────────────────
single_thread_mask = ~(
    acc_summary["tool"].str.ends_with("_parallel")
    | acc_summary["tool"].str.ends_with("_multithreaded")
)
p2_data = (
    acc_df.filter(
        ~pl.col("tool").str.ends_with("_parallel")
        & ~pl.col("tool").str.ends_with("_multithreaded")
    )
    .group_by("tool")
    .agg([pl.col(s).mean() for s in segments])
    .with_columns(pl.col("tool").replace(single_label_map).alias("method"))
    .drop("tool")
    .unpivot(
        index=["method"], on=segments, variable_name="segment", value_name="pct_correct"
    )
)

plot2 = (
    alt.Chart(p2_data.to_pandas())
    .mark_bar(opacity=0.7, stroke="black", strokeWidth=1.5, size=20)
    .encode(
        y=alt.Y("pct_correct:Q", title="% Correct", scale=alt.Scale(domain=[0, 100])),
        x=alt.X(
            "segment:N",
            title="Segment",
            sort=segments,
            axis=alt.Axis(labelAngle=-45),
            scale=alt.Scale(paddingInner=0.3),
        ),
        color=alt.Color("method:N", title="Tool"),
        xOffset=alt.XOffset("method:N", scale=alt.Scale(paddingInner=0.01)),
        tooltip=["method", "segment", alt.Tooltip("pct_correct:Q", format=".2f")],
    )
    .properties(
        title=f"Correctness by segment ({acc_size_label}, averaged across chains)",
        width=1000,
        height=400,
    )
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


chart = (
    alt.vconcat(plot_perf, plot2)
    .resolve_scale(color="independent")
    .configure(background=WHITE)
    .configure_axis(
        labelColor=BLACK,
        titleColor=BLACK,
        gridColor="#dddddd",
        domainColor=BLACK,
        tickColor=BLACK,
    )
    .configure_title(color=BLACK)
    .configure_legend(labelColor=BLACK, titleColor=BLACK)
    .configure_view(stroke=BLACK)
)

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
