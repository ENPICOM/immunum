# Polars

`immunum.polars` numbers and segments sequences as Polars expressions. They run inside a Polars
query on every core, work in lazy and streaming queries, and only compute the fields you keep. For
hundreds of thousands of sequences or more, this is the interface to use.

It needs Polars: install with `pip install "immunum[polars]"`
([Getting started](../getting-started.md#python)). Every parameter is in the
[API reference](../api/polars.md).

## Expressions

| Expression | Returns per sequence |
|---|---|
| `number` | a struct with the fields `Annotator.number` returns: `chain`, `scheme`, `confidence`, `numbering` (a list of `{position, residue}` structs), `query_start`, `query_end`, `error`, `error_kind` |
| `segment` | a struct with `prefix`, `fr1`, `cdr1`, `fr2`, `cdr2`, `fr3`, `cdr3`, `fr4`, `postfix`, `error`, `error_kind` |
| `number_domains`, `segment_domains` | a list of those structs, one per domain, such as both domains of an scFv |

Every expression takes either `chains`, `scheme` and optionally `min_confidence`, as
[`Annotator`](python.md#creating-an-annotator) does, or a prebuilt `Annotator` as `annotator=`.
Both give the same result and run as fast. Either way the annotator is built when the expression
is, so an unknown name raises `immunum.Error` right there, before any query runs.
`numbering_method()` and `segmentation_method()` are deprecated; use `number(..., annotator=...)`
and `segment(..., annotator=...)`.

## Example: regions, residues and domains

```python
import polars as pl
import immunum.polars as imp

heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
kappa = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
df = pl.DataFrame({
    "id": ["heavy", "kappa", "scfv"],
    "sequence": [heavy, kappa, heavy + "GGGGSGGGGSGGGGS" + kappa],
})

# One column per region; `segment` takes the best-scoring domain, here the scFv's kappa domain
regions = df.select(
    "id",
    imp.segment("sequence", chains=["ig"], scheme="imgt").struct.unnest(),
).select("id", "cdr1", "cdr2", "cdr3", "error_kind")
print(regions)
# shape: (3, 5)
# ┌───────┬──────────┬──────────┬─────────────────┬────────────┐
# │ id    ┆ cdr1     ┆ cdr2     ┆ cdr3            ┆ error_kind │
# │ ---   ┆ ---      ┆ ---      ┆ ---             ┆ ---        │
# │ str   ┆ str      ┆ str      ┆ str             ┆ str        │
# ╞═══════╪══════════╪══════════╪═════════════════╪════════════╡
# │ heavy ┆ GGSFSTYA ┆ VIPLLTIT ┆ AREGTTGKPIGAFAH ┆ null       │
# │ kappa ┆ QDVNTA   ┆ SAS      ┆ QQHYTTPPT       ┆ null       │
# │ scfv  ┆ QDVNTA   ┆ SAS      ┆ QQHYTTPPT       ┆ null       │
# └───────┴──────────┴──────────┴─────────────────┴────────────┘
assert regions["cdr3"].to_list() == ["AREGTTGKPIGAFAH", "QQHYTTPPT", "QQHYTTPPT"]

# One row per numbered residue: explode `numbering` and unnest it
residues = (
    df.filter(pl.col("id") == "heavy")
    .select("id", imp.number("sequence", chains=["ig"], scheme="imgt").alias("numbered"))
    .unnest("numbered")
    .explode("numbering")
    .unnest("numbering")
    .select("id", "chain", "position", "residue")
    .filter(pl.col("position").is_in(["110", "111", "111A", "112A", "112"]))
)
print(residues)
# shape: (5, 4)
# ┌───────┬───────┬──────────┬─────────┐
# │ id    ┆ chain ┆ position ┆ residue │
# │ ---   ┆ ---   ┆ ---      ┆ ---     │
# │ str   ┆ str   ┆ str      ┆ str     │
# ╞═══════╪═══════╪══════════╪═════════╡
# │ heavy ┆ H     ┆ 110      ┆ T       │
# │ heavy ┆ H     ┆ 111      ┆ G       │
# │ heavy ┆ H     ┆ 111A     ┆ K       │
# │ heavy ┆ H     ┆ 112A     ┆ P       │
# │ heavy ┆ H     ┆ 112      ┆ I       │
# └───────┴───────┴──────────┴─────────┘
assert residues["position"].to_list() == ["110", "111", "111A", "112A", "112"]

# Every domain, one row each, in a lazy query; the linker is the kappa domain's prefix
domains = (
    df.lazy()
    .select("id", imp.segment_domains("sequence", chains=["ig"], scheme="imgt").alias("domain"))
    .explode("domain")
    .unnest("domain")
    .select("id", "prefix", "cdr3")
    .collect()
)
print(domains)
# shape: (4, 3)
# ┌───────┬─────────────────┬─────────────────┐
# │ id    ┆ prefix          ┆ cdr3            │
# │ ---   ┆ ---             ┆ ---             │
# │ str   ┆ str             ┆ str             │
# ╞═══════╪═════════════════╪═════════════════╡
# │ heavy ┆                 ┆ AREGTTGKPIGAFAH │
# │ kappa ┆                 ┆ QQHYTTPPT       │
# │ scfv  ┆                 ┆ AREGTTGKPIGAFAH │
# │ scfv  ┆ GGGGSGGGGSGGGGS ┆ QQHYTTPPT       │
# └───────┴─────────────────┴─────────────────┘
assert domains["cdr3"].to_list() == ["AREGTTGKPIGAFAH", "QQHYTTPPT", "AREGTTGKPIGAFAH", "QQHYTTPPT"]
```

## Performance

- Select only the fields you need, and as late as possible: on big data, queries get
  **significantly** faster. The first example below runs in under a second on 8M rows (200,000 of
  them paired sequences), where materialising every field would take about 10 times longer.
- Pass only the chains you expect: each chain is one more alignment per sequence.
- Use `pl.scan_*` and `.collect(engine="streaming")` for data larger than memory.

## Errors

A sequence that can't be numbered never fails the query: its struct has `error` and `error_kind`
set and every other field null. A null sequence gives a struct that's null throughout, `error`
included, so test success with `chain` (from `number`) or `fr1` (from `segment`) rather than with
`error` being null. Details and an example: [Errors](errors.md#polars).

## Examples on large datasets

These examples read parquet files of paired sequences that aren't part of the repository, so they
don't run in the documentation's tests; their output is from a real run.

### Example: calculating unique cdr3_light/heavy pairs

Segment heavy and light chains in parallel, then deduplicate on their CDR3 sequences.
Note that `.unique()` is called before materializing any other columns — this avoids
loading the full FR region strings for every row.

```python
import pytest
import polars as pl
import immunum.polars as imp


def unique_cdr3_pairs(source: str) -> pl.DataFrame:
    return (
        pl.scan_parquet(source)
        .filter(
            pl.all_horizontal(
                pl.col("sequence_alignment_aa_heavy").is_not_null(),
                pl.col("sequence_alignment_aa_light").is_not_null(),
            )
        )
        .select(
            imp.segment("sequence_alignment_aa_heavy", chains=["IGH"], scheme="IMGT")
            .name.suffix_fields("_heavy")
            .struct.unnest(),
            imp.segment("sequence_alignment_aa_light", chains=["IGL", "IGK"], scheme="IMGT")
            .name.suffix_fields("_light")
            .struct.unnest(),
        )
        .unique(subset=["cdr3_heavy", "cdr3_light"])
        .select("cdr3_heavy", "cdr3_light")
        .head(20)
        .collect(engine="streaming")
    )


pytest.skip("requires parquet data files")
print(unique_cdr3_pairs("source.parquet"))

# shape: (20, 2)
# ┌──────────────────────┬──────────────┐
# │ cdr3_heavy           ┆ cdr3_light   │
# │ ---                  ┆ ---          │
# │ str                  ┆ str          │
# ╞══════════════════════╪══════════════╡
# │ ARDLSQGYFDY          ┆ QQYYSTPYT    │
# │ ASLRGITGTTDS         ┆ AAWDDSLKGVV  │
# │ ARERPQALCFDP         ┆ QQYSSSLYT    │
# │ AKDLEGKHHYFDY        ┆ QQYGSSLWT    │
# │ ARAGSWNKDYYDSSGPLDY  ┆ QQYNSYFTWT   │
# │ …                    ┆ …            │
# │ ALTSGYSSGWPFDY       ┆ QQYYSTPLT    │
# │ ARHRSIAARPAIYYYYYMDV ┆ QSYDSSLSGFYV │
# │ ARNLRGYSYGPDAFDI     ┆ QQYDNLPYT    │
# │ VRSPGWSFDF           ┆ QQYGSSPSPMYT │
# │ ARQKVGSRSPNWYFDL     ┆ QQYNNWPRT    │
# └──────────────────────┴──────────────┘
```

### Example: matching with another dataframe on framework sets

For each light+heavy framework combination the source and reference datasets share, count the
matching source/reference sequence pairs (a combination found twice in each counts 4), then
return the 20 most common combinations.

```python
import pytest
import polars as pl
import immunum.polars as imp

HEAVY = "sequence_alignment_aa_heavy"
LIGHT = "sequence_alignment_aa_light"
FRAMEWORK_COLS = ["fr1_heavy", "fr2_heavy", "fr3_heavy", "fr1_light", "fr2_light", "fr3_light"]


def top_framework_groups(source: str, target: str, n: int = 10_000) -> pl.DataFrame:
    def _segment(path: str) -> pl.LazyFrame:
        return (
            pl.scan_parquet(path)
            .filter(
                pl.all_horizontal(
                    pl.col(HEAVY).is_not_null(),
                    pl.col(LIGHT).is_not_null(),
                )
            )
            .head(n)
            .select(
                imp.segment(HEAVY, chains=["IGH"], scheme="IMGT")
                .name.suffix_fields("_heavy")
                .struct.unnest(),
                imp.segment(LIGHT, chains=["IGK", "IGL"], scheme="IMGT")
                .name.suffix_fields("_light")
                .struct.unnest(),
            )
        )

    return (
        _segment(source)
        .join(_segment(target), on=FRAMEWORK_COLS, how="inner")
        .group_by(FRAMEWORK_COLS)
        .agg(pl.len().alias("count"))
        .sort("count", descending=True)
        .head(20)
        .collect(engine="streaming")
    )


pytest.skip("requires parquet data files")
print(top_framework_groups("source.parquet", "target.parquet"))

# shape: (20, 7)
# ┌───────────────────────────┬───────────────────┬─────────────────────────────────┬────────────────────────────┬───────────────────┬─────────────────────────────────┬───────┐
# │ fr1_heavy                 ┆ fr2_heavy         ┆ fr3_heavy                       ┆ fr1_light                  ┆ fr2_light         ┆ fr3_light                       ┆ count │
# │ ---                       ┆ ---               ┆ ---                             ┆ ---                        ┆ ---               ┆ ---                             ┆ ---   │
# │ str                       ┆ str               ┆ str                             ┆ str                        ┆ str               ┆ str                             ┆ u32   │
# ╞═══════════════════════════╪═══════════════════╪═════════════════════════════════╪════════════════════════════╪═══════════════════╪═════════════════════════════════╪═══════╡
# │ EVQLLESGGGLVQPGGSLRLSCAAS ┆ MSWVRQAPGKGLEWVSA ┆ YYADSVKGRFTISRDNSKNTLYLQMNSLRA… ┆ EIVLTQSPGTLSLSPGERATLSCRAS ┆ LAWYQQKPGQAPRLLIY ┆ SRATGIPDRFSGSGSGTDFTLTISRLEPED… ┆ 3248  │
# │ QVQLVESGGGVVQPGRSLRLSCAAS ┆ MHWVRQAPGKGLEWVAV ┆ YYADSVKGRFTISRDNSKNTLYLQMNSLRA… ┆ DIQMTQSPSSLSASVGDRVTITCRAS ┆ LNWYQQKPGKAPKLLIY ┆ SLQSGVPSRFSGSGSGTDFTLTISSLQPED… ┆ 2025  │
# │ QVQLQESGPGLVKPSETLSLTCTVS ┆ WSWIRQPPGKGLEWIGY ┆ NYNPSLKSRVTISVDTSKNQFSLKLSSVTA… ┆ DIQMTQSPSSLSASVGDRVTITCRAS ┆ LNWYQQKPGKAPKLLIY ┆ SLQSGVPSRFSGSGSGTDFTLTISSLQPED… ┆ 1998  │
# │ EVQLLESGGGLVQPGGSLRLSCAAS ┆ MSWVRQAPGKGLEWVSA ┆ YYADSVKGRFTISRDNSKNTLYLQMNSLRA… ┆ EIVMTQSPATLSVSPGERATLSCRAS ┆ LAWYQQKPGQAPRLLIY ┆ TRATGIPARFSGSGSGTEFTLTISSLQSED… ┆ 1961  │
# │ QVQLQESGPGLVKPSETLSLTCTVS ┆ WSWIRQPPGKGLEWIGY ┆ NYNPSLKSRVTISVDTSKNQFSLKLSSVTA… ┆ EIVLTQSPGTLSLSPGERATLSCRAS ┆ LAWYQQKPGQAPRLLIY ┆ SRATGIPDRFSGSGSGTDFTLTISRLEPED… ┆ 1872  │
# │ …                         ┆ …                 ┆ …                               ┆ …                          ┆ …                 ┆ …                               ┆ …     │
# │ EVQLVESGGGLVQPGRSLRLSCAAS ┆ MHWVRQAPGKGLEWVSG ┆ GYADSVKGRFTISRDNAKNSLYLQMNSLRA… ┆ EIVLTQSPGTLSLSPGERATLSCRAS ┆ LAWYQQKPGQAPRLLIY ┆ SRATGIPDRFSGSGSGTDFTLTISRLEPED… ┆ 990   │
# │ EVQLLESGGGLVQPGGSLRLSCAAS ┆ MSWVRQAPGKGLEWVSA ┆ YYADSVKGRFTISRDNSKNTLYLQMNSLRA… ┆ EIVLTQSPATLSLSPGERATLSCRAS ┆ LAWYQQKPGQAPRLLIY ┆ NRATGIPARFSGSGSGTDFTLTISSLEPED… ┆ 920   │
# │ QVQLVESGGGVVQPGRSLRLSCAAS ┆ MHWVRQAPGKGLEWVAV ┆ YYADSVKGRFTISRDNSKNTLYLQMNSLRA… ┆ DIQMTQSPSTLSASVGDRVTITCRAS ┆ LAWYQQKPGKAPKLLIY ┆ SLESGVPSRFSGSGSGTEFTLTISSLQPDD… ┆ 900   │
# │ QVQLVQSGAEVKKPGSSVKVSCKAS ┆ ISWVRQAPGQGLEWMGG ┆ NYAQKFQGRVTITADESTSTAYMELSSLRS… ┆ EIVLTQSPGTLSLSPGERATLSCRAS ┆ LAWYQQKPGQAPRLLIY ┆ SRATGIPDRFSGSGSGTDFTLTISRLEPED… ┆ 899   │
# │ EVQLVESGGGLVQPGGSLRLSCAAS ┆ MSWVRQAPGKGLEWVAN ┆ YYVDSVKGRFTISRDNAKNSLYLQMNSLRA… ┆ DIQMTQSPSSLSASVGDRVTITCRAS ┆ LNWYQQKPGKAPKLLIY ┆ SLQSGVPSRFSGSGSGTDFTLTISSLQPED… ┆ 868   │
# └───────────────────────────┴───────────────────┴─────────────────────────────────┴────────────────────────────┴───────────────────┴─────────────────────────────────┴───────┘
```

## Polars integrations

Ecosystem of `polars` has many handy tools for data engineering -- for instance, you can check out [awesome-polars](https://github.com/ddotta/awesome-polars) list. Here we provide examples of two packages that we use actively: [`polars-distance`](https://github.com/ion-elgreco/polars-distance) providing convenient (and fast!) interface for string distance computations, and [`polars-bio`](https://biodatageeks.org/polars-bio/) providing direct interface for reading biological data formats (fasta/fastq/...) into polars dataframes.

### Example: annotating on-the-fly with polars-bio

[polars-bio](https://github.com/biodatageeks/polars-bio) can scan biological sequence files
(FASTA, FASTQ) directly as a `LazyFrame`. Chain that with `imp.segment` to annotate without
an intermediate file:

```python
import pytest
import polars as pl
pb = pytest.importorskip("polars_bio")
import immunum.polars as imp


def annotate_fasta(path: str) -> pl.DataFrame:
    return (
        pb.scan_fasta(path)
        .select(
            pl.col("name"),
            imp.segment("sequence", chains=["IGH", "IGK", "IGL"], scheme="IMGT")
            .struct.unnest(),
        )
        .collect(engine="streaming")
    )


pytest.skip("requires polars_bio and a FASTA file")
print(annotate_fasta("sequences.fasta"))
```

### Example: computing minimal distance between CDR3s in two datasets

Cross-join two CDR3 sets and keep the nearest neighbour by Levenshtein distance.
Pre-filter to unique CDR3s and, if the datasets are large, narrow the search to
equal-length sequences first to reduce the cross-product size.

```python
import pytest
import polars as pl
pld = pytest.importorskip("polars_distance")
import immunum.polars as imp


def nearest_cdr3s(source: str, target: str) -> pl.DataFrame:
    def _extract(path: str, alias: str) -> pl.DataFrame:
        return (
            pl.scan_parquet(path)
            .select(
                imp.segment("sequence_aa", chains=["IGH"], scheme="IMGT")
                .struct.field("cdr3")
                .alias(alias)
            )
            .filter(pl.col(alias).is_not_null())
            .unique(alias)
            .collect()
        )

    return (
        _extract(source, "cdr3_a")
        .join(_extract(target, "cdr3_b"), how="cross")
        .select(
            pl.col("cdr3_a"),
            pl.col("cdr3_b"),
            pld.col("cdr3_a").dist_str.levenshtein("cdr3_b").alias("dist"),
        )
        .group_by("cdr3_a")
        .agg(
            pl.col("dist").min().alias("min_dist"),
            pl.col("cdr3_b").sort_by("dist").first().alias("nearest_cdr3"),
        )
        .sort("min_dist")
        .head(20)
    )


pytest.skip("requires polars_distance and parquet data files")
print(nearest_cdr3s("source.parquet", "target.parquet"))

# shape: (20, 3)
# ┌───────────────────┬──────────┬───────────────────┐
# │ cdr3_a            ┆ min_dist ┆ nearest_cdr3      │
# │ ---               ┆ ---      ┆ ---               │
# │ str               ┆ u32      ┆ str               │
# ╞═══════════════════╪══════════╪═══════════════════╡
# │ ARDRDDPMADYHPLFDS ┆ 0        ┆ ARDRDDPMADYHPLFDS │
# │ ARESPPRLGHWYFDL   ┆ 0        ┆ ARESPPRLGHWYFDL   │
# │ ARESGRGVVSPYFDP   ┆ 0        ┆ ARESGRGVVSPYFDP   │
# │ ARHRGSTINIPYFDY   ┆ 0        ┆ ARHRGSTINIPYFDY   │
# │ ARDGGYSGSPWYYFDY  ┆ 0        ┆ ARDGGYSGSPWYYFDY  │
# │ …                 ┆ …        ┆ …                 │
# │ AKERGSTGSAINY     ┆ 0        ┆ AKERGSTGSAINY     │
# │ AATTRDWFDP        ┆ 0        ┆ AATTRDWFDP        │
# │ ARDPDTSNKIDY      ┆ 0        ┆ ARDPDTSNKIDY      │
# │ ARDRSSDY          ┆ 0        ┆ ARDRSSDY          │
# │ ATHWDWRFDN        ┆ 0        ┆ ATHWDWRFDN        │
# └───────────────────┴──────────┴───────────────────┘
```
