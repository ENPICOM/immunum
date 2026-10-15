"""MkDocs hook: fill in pages from files the build already has, so the site can't drift from them.

- `<!-- typescript-declarations -->` becomes the TypeScript declarations of the npm package: the
  web bundle's `immunum.d.ts` (built by `task build-wasm-web`, which `task docs-build` runs first)
  without wasm-bindgen's loader functions, which the Node.js package doesn't have.
- `<!-- benchmark-accuracy -->` becomes a table of the per-scheme accuracy in BENCHMARKS.toml.
- `<!-- benchmark-speed -->` becomes a table of the mean times in the headline comparison,
  resources/benchmark_results/results_headline.csv.
- `<!-- immunum-version -->` becomes the version the site documents, which mkdocs_version.py
  reads from Cargo.toml; that hook must run first.
- A page that only includes one file through pymdownx.snippets, such as the changelog, gets an
  "edit this page" link to that file instead of to its one-line wrapper.
"""

import csv
import re
import sys
from pathlib import Path

from mkdocs.config.defaults import MkDocsConfig
from mkdocs.structure.files import Files
from mkdocs.structure.pages import Page

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

ROOT = Path(__file__).resolve().parents[1]
TYPESCRIPT = ROOT / "docs" / "assets" / "wasm" / "immunum.d.ts"
BENCHMARKS = ROOT / "BENCHMARKS.toml"
HEADLINE = ROOT / "resources" / "benchmark_results" / "results_headline.csv"
# Tool, then its runner on one core and on all cores in bench_comparison.py
HEADLINE_TOOLS = [
    ("immunum", "immunum_singlethreaded", "immunum_multithreaded"),
    ("AntPack", "antpack", "antpack_parallel"),
    ("ANARCI", "anarci", "anarci_parallel"),
    ("ANARCII", "anarcii2", "anarcii2_parallel"),
]

# wasm-bindgen appends its loader after the declarations from src/wasm.rs; it starts here
LOADER_START = "export type InitInput"
SNIPPET = re.compile(r'^--8<-- "([^":]+)(?::[^"]*)?"$')

SCHEMES = {
    "imgt": "IMGT",
    "kabat": "Kabat",
    "chothia": "Chothia",
    "martin": "Martin",
    "aho": "AHo",
}
CHAINS = "HKLABGD"


def typescript_declarations() -> str:
    if not TYPESCRIPT.exists():
        raise FileNotFoundError(
            f"{TYPESCRIPT} is missing; run `task build-wasm-web` first"
        )
    declarations = TYPESCRIPT.read_text()
    start = declarations.find(LOADER_START)
    if start == -1:
        raise ValueError(
            f"{TYPESCRIPT}: no `{LOADER_START}` to cut wasm-bindgen's loader at"
        )
    body = declarations[:start]
    # wasm-bindgen's lint pragmas mean nothing to a reader
    body = re.sub(r"^/\* (tslint|eslint)-disable \*/\n", "", body, flags=re.MULTILINE)
    return f"```ts\n{body.strip()}\n```"


def benchmark_accuracy() -> str:
    results = tomllib.loads(BENCHMARKS.read_text())
    rows = [
        "| Scheme | Chain | Sequences | Numbered exactly | Residues correct |",
        "|---|---|---:|---:|---:|",
    ]
    for scheme, label in SCHEMES.items():
        for chain in CHAINS:
            r = results.get(scheme, {}).get(chain)
            if r is None:
                continue
            rows.append(
                f"| {label} | {chain} | {r['total_sequences']:,} "
                f"| {r['perfect_pct']:.2f}% ({r['perfect_sequences']:,}) "
                f"| {r['overall_accuracy']:.2f}% |"
            )
    return f"Last run: {results['last_updated']}.\n\n" + "\n".join(rows)


def benchmark_speed() -> str:
    times: dict[str, list[float]] = {}
    numbered: dict[str, list[int]] = {}
    with HEADLINE.open() as f:
        rows = list(csv.DictReader(f))
    size = int(rows[0]["sample_size"])
    for row in rows:
        times.setdefault(row["tool"], []).append(float(row["time_s"]))
        numbered.setdefault(row["tool"], []).append(int(row["numbered"]))

    def cell(runner: str) -> str:
        t = sum(times[runner]) / len(times[runner])
        return f"{t:,.2f} s | {size / t:,.0f}"

    table = [
        "| Tool | One core | Sequences/s | All cores | Sequences/s | Numbered |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for label, single, multi in HEADLINE_TOOLS:
        if single in times and multi in times:
            table.append(
                f"| {label} | {cell(single)} | {cell(multi)} "
                f"| {min(numbered[single]):,} of {size:,} |"
            )
    return "\n".join(table)


def on_page_markdown(
    markdown: str, page: Page, config: MkDocsConfig, files: Files
) -> str:
    if "<!-- typescript-declarations -->" in markdown:
        markdown = markdown.replace(
            "<!-- typescript-declarations -->", typescript_declarations()
        )
    if "<!-- immunum-version -->" in markdown:
        markdown = markdown.replace(
            "<!-- immunum-version -->", config.extra["immunum_version"]
        )
    if "<!-- benchmark-accuracy -->" in markdown:
        markdown = markdown.replace("<!-- benchmark-accuracy -->", benchmark_accuracy())
    if "<!-- benchmark-speed -->" in markdown:
        markdown = markdown.replace("<!-- benchmark-speed -->", benchmark_speed())

    lines = [line for line in markdown.splitlines() if line.strip()]
    included = {m.group(1) for line in lines if (m := SNIPPET.match(line.strip()))}
    if len(included) == 1 and all(SNIPPET.match(line.strip()) for line in lines):
        (source,) = included
        page.edit_url = f"{config.repo_url}/edit/main/{source}"
    return markdown
