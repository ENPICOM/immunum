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

# Summarise speed data (mean ± std per tool/sample_size)
speed_summary = (
    speed_df.group_by(["tool", "sample_size"])
    .agg(
        pl.col("time_s").mean().alias("mean_time"),
        pl.col("time_s").std().alias("std_time"),
    )
    .sort(["tool", "sample_size"])
)

# Summarise accuracy data (mean per tool/chain/segment)
segments = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
acc_summary = (
    acc_df.group_by(["tool", "chain"])
    .agg([pl.col(s).mean() for s in segments])
    .sort(["chain", "tool"])
)

# ── Plot 1: Single-threaded performance ──────────────────────────────────────
single_tools = [
    "immunum_singlethreaded",
    "anarci",
    "anarcii2",
    "antpack",
    "riot",
]
label_map = {
    "immunum_singlethreaded": "immunum",
    "anarci": "anarci",
    "anarcii2": "anarcii2",
    "antpack": "antpack",
    "riot": "riot",
}

p1_data = speed_summary.filter(pl.col("tool").is_in(single_tools)).with_columns(
    pl.col("tool").replace(label_map).alias("method"),
    pl.col("sample_size").cast(pl.Utf8).alias("n_sequences"),
)

plot1 = (
    alt.Chart(p1_data.to_pandas())
    .mark_bar()
    .encode(
        x=alt.X("mean_time:Q", title="Time (seconds)"),
        y=alt.Y(
            "method:N",
            title="Method",
            sort=alt.EncodingSortField(field="mean_time", op="mean", order="ascending"),
        ),
        color=alt.Color("n_sequences:N", title="Number of sequences"),
        yOffset="n_sequences:N",
        tooltip=["method", "n_sequences", alt.Tooltip("mean_time:Q", format=".4f")],
    )
    .properties(title="Single-threaded performance", width=500, height=300)
)

# ── Plot 2: Correctness by chain and segment ──────────────────────────────────
single_thread_tools = ~(
    acc_summary["tool"].str.ends_with("_parallel")
    | acc_summary["tool"].str.ends_with("_multithreaded")
)
p2_data = (
    acc_summary.filter(single_thread_tools)
    .with_columns((pl.col("tool") + " (" + pl.col("chain") + ")").alias("method"))
    .unpivot(
        index=["method"], on=segments, variable_name="segment", value_name="pct_correct"
    )
)

plot2 = (
    alt.Chart(p2_data.to_pandas())
    .mark_bar()
    .encode(
        x=alt.X("pct_correct:Q", title="% Correct", scale=alt.Scale(domain=[0, 100])),
        y=alt.Y("segment:N", title="Segment", sort=segments),
        color=alt.Color("method:N", title="Tool (chain)"),
        yOffset="method:N",
        tooltip=["method", "segment", alt.Tooltip("pct_correct:Q", format=".2f")],
    )
    .properties(
        title=f"Correctness by segment and chain ({acc_size_label})",
        width=500,
        height=400,
    )
)

# ── Plot 3: Multi-threaded performance (line + error bars) ───────────────────
mt_tools = [
    "immunum_multithreaded",
    "antpack_parallel",
    "anarci_parallel",
    "anarcii2_parallel",
    "riot_parallel",
]
mt_label_map = {
    "immunum_multithreaded": "immunum (multithreaded)",
    "antpack_parallel": "antpack (parallel)",
    "anarci_parallel": "anarci (parallel)",
    "anarcii2_parallel": "anarcii2 (parallel)",
    "riot_parallel": "riot (parallel)",
}

p3_data = speed_summary.filter(pl.col("tool").is_in(mt_tools)).with_columns(
    pl.col("tool").replace(mt_label_map).alias("method"),
    (pl.col("mean_time") - pl.col("std_time")).alias("lower"),
    (pl.col("mean_time") + pl.col("std_time")).alias("upper"),
)

base = alt.Chart(p3_data.to_pandas()).encode(
    x=alt.X(
        "sample_size:Q",
        title="Number of sequences",
        scale=alt.Scale(type="log", base=10),
        axis=alt.Axis(format="~s"),
    ),
    color=alt.Color("method:N", title="Method"),
)

lines = base.mark_line(point=True).encode(
    y=alt.Y("mean_time:Q", title="Time (seconds)"),
    tooltip=["method", "sample_size", alt.Tooltip("mean_time:Q", format=".3f")],
)

error_bars = base.mark_errorband().encode(
    y=alt.Y("lower:Q", title="Time (seconds)"),
    y2="upper:Q",
)

plot3 = (error_bars + lines).properties(
    title="Multi-threaded / parallel performance", width=500, height=350
)

# ── Save ─────────────────────────────────────────────────────────────────────
chart = alt.vconcat(plot1, plot2, plot3).resolve_scale(color="independent")
chart.save("benchmark_plots.html")
print("Saved benchmark_plots.html")

for name, p in [
    ("plot1_singlethreaded", plot1),
    ("plot2_correctness", plot2),
    ("plot3_parallel", plot3),
]:
    path = f"benchmark_{name}.svg"
    p.save(path)
    print(f"Saved {path}")
