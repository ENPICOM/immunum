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


@alt.theme.register("enpicom", enable=True)
def altair_theme() -> dict:
    # ENPICOM Brand Colors from Guidebook
    # Primary Colors
    dark_blue = "#231E60"  # Primary for titles, body text, and key elements
    purple = "#7D5BFF"  # Primary for emphasis and hierarchy

    # Secondary Colors
    blue = "#3229C1"  # For subheadings, fine lines, or subtle outlines
    off_white = "#F0F4FC"  # For backgrounds

    # Accent Colors
    periwinkle = "#C7D6F4"
    lavander = "#D8CCFC"
    yellow = "#F9E7A4"
    mint = "#DCF7F0"

    # Font settings
    font = "Albert Sans"
    background_color = off_white
    base_size = 16
    header_font = base_size * 1.25  # Renamed lg_font to header_font
    sm_font = base_size * 0.8  # st.table size
    title_font = base_size * 1.75  # For main chart titles

    # Glassmorphism effect - adding opacity to certain elements
    purple_glass = purple + "CC"  # ~80% opacity
    blue_glass = blue + "CC"  # ~80% opacity

    config = {
        "config": {
            "view": {"fill": background_color},
            "arc": {"fill": purple},
            "area": {"fill": purple_glass},
            "circle": {"fill": purple, "stroke": dark_blue, "strokeWidth": 0.5},
            "line": {"stroke": blue},
            "path": {"stroke": blue},
            "point": {"stroke": purple},
            "rect": {"fill": purple_glass},
            "shape": {"stroke": blue},
            "symbol": {"fill": purple},
            "title": {
                "font": font,
                "fontWeight": "SemiBold",
                "color": dark_blue,
                "fontSize": title_font,  # Using the renamed variable
                "anchor": "start",
            },
            "axis": {
                "titleFont": font,
                "titleFontWeight": "SemiBold",
                "titleColor": purple,
                "titleFontSize": sm_font,
                "labelFont": font,
                "labelFontWeight": "Light",
                "labelColor": blue,
                "labelFontSize": sm_font,
                "grid": True,
                "gridColor": blue_glass,
                "gridOpacity": 0.3,
                "domain": False,
                # "domainColor": dark_blue,
                "tickColor": blue,
            },
            "header": {
                "labelFont": font,
                "labelFontWeight": "Light",
                "titleFont": font,
                "titleFontWeight": "SemiBold",
                "labelFontSize": base_size,
                "titleFontSize": header_font,  # Using the renamed header_font variable
            },
            "legend": {
                "titleFont": font,
                "titleFontWeight": "SemiBold",
                "titleColor": purple,
                "titleFontSize": sm_font,
                "labelFont": font,
                "labelFontWeight": "Light",
                "labelColor": dark_blue,
                "labelFontSize": sm_font,
                "fillOpacity": 0.8,
                "strokeOpacity": 0.8,
                "symbolOpacity": 0.9,
            },
            "range": {
                # Primary category colors using ENPICOM brand palette
                "category": [
                    purple,  # Primary accent
                    dark_blue,  # Primary dark
                    blue,  # Secondary
                    periwinkle,  # Accent
                    lavander,  # Accent
                    yellow,  # Accent
                    mint,  # Accent
                    purple + "99",  # Primary with 60% opacity
                    blue + "99",  # Secondary with 60% opacity
                ],
                "diverging": [
                    dark_blue,
                    "#2B2580",  # Slightly lighter dark blue
                    "#342AA0",  # Mid tone between dark blue and blue
                    blue,
                    "#584BE0",  # Mid tone between blue and purple
                    purple,
                    "#9A82FF",  # Slightly lighter purple
                    "#B8A9FF",  # Even lighter purple
                    lavander,
                ],
                "heatmap": [
                    mint,
                    "#E3F9F3",
                    "#EAF2FA",
                    off_white,
                    "#DED6FC",
                    lavander,
                    "#9F88FF",
                    purple,
                    dark_blue,
                ],
                "ramp": [
                    mint,
                    "#E3F9F3",
                    "#EAF2FA",
                    off_white,
                    periwinkle,
                    "#A2B7ED",
                    "#7D91D5",
                    "#584BE0",
                    blue,
                    dark_blue,
                ],
                "ordinal": [
                    dark_blue,
                    blue,
                    purple,
                    periwinkle,
                    lavander,
                    yellow,
                    mint,
                ],
            },
        }
    }
    return config


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
    "riot",
]
single_label_map = {
    "immunum_singlethreaded": "immunum",
    "anarci": "anarci",
    "anarcii2": "anarcii2",
    "antpack": "antpack",
    "riot": "riot",
}
mt_tools = [
    "immunum_multithreaded",
    "antpack_parallel",
    "anarci_parallel",
    "anarcii2_parallel",
    "riot_parallel",
]
mt_label_map = {
    "immunum_multithreaded": "immunum",
    "antpack_parallel": "antpack",
    "anarci_parallel": "anarci",
    "anarcii2_parallel": "anarcii2",
    "riot_parallel": "riot",
}

METHODS = ["immunum", "anarci", "anarcii2", "antpack", "riot"]
_enpicom_category = [
    "#7D5BFF",
    "#231E60",
    "#3229C1",
    "#C7D6F4",
    "#D8CCFC",
    "#F9E7A4",
    "#DCF7F0",
]
color_scale = alt.Scale(domain=METHODS, range=_enpicom_category[: len(METHODS)])

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
        .agg(pl.col("throughput").mean())
        .sort("throughput", descending=True)
    )
    return (
        alt.Chart(summary.to_pandas())
        .mark_bar(opacity=0.7, stroke="black", strokeWidth=1.5)
        .encode(
            x=alt.X("method:N", title="Method", sort="-y"),
            y=alt.Y(
                "throughput:Q",
                title="Sequences per second",
                # scale=alt.Scale(type="log", base=10),
            ),
            color=alt.Color("method:N", scale=color_scale, legend=None),
            tooltip=["method", alt.Tooltip("throughput:Q", format=",.0f")],
        )
        .properties(title=title, width=300, height=350)
    )


plot_perf = (
    make_perf_barplot(p1_data, "Single-threaded")
    | make_perf_barplot(p3_data, "Multi-threaded")
).resolve_scale(y="shared")

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
        color=alt.Color("method:N", title="Tool", scale=color_scale),
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
    ("plot2_correctness", plot2),
]:
    path = f"assets/benchmark_{name}.svg"
    p.save(path)
    print(f"Saved {path}")
