# Introduction

As advertised, immunum annotates antibody and T-cell receptor sequences with IMGT or Kabat position
numbers (`number` functions), or segments them into FR/CDR regions plus optional prefixes (if sequence contains more than `immunum` is confident about). Segmentation is guaranteed to cover the whole sequence, i.e. string sum of segments always equals input sequence.

Two interfaces are provided:

- **Python** (`immunum.Annotator`) — for single sequences or small batches
- **Polars** (`immunum.polars`) — for large datasets; runs as a lazy expression inside
  a Polars query plan

When multiple chains are provided, `immunum` scores each sequence against every chain type you pass and picks the best-matching one. The more chains you include, the more comparisons are made per
sequence -— so pass only the chains you expect in your data. E.g. including 6 chains in `chains=...` argument will make numbering/segmentation six times slower.

# Python module

Create an `Annotator` with the chain types and numbering scheme you need, then call
`number()` or `segment()` on each sequence.

**Numbering** assigns an IMGT or Kabat position label to every residue:

```python
from immunum import Annotator

annotator = Annotator(chains=["H", "K", "L"], scheme="imgt")

result = annotator.number(
    "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
)
print(result.chain)       # "H"
print(result.scheme)      # "IMGT"
print(result.numbering["1"])  # "Q"
```

**Segmentation** splits the sequence into FR1–FR4 and CDR1–CDR3 regions:

```python
result = annotator.segment(
    "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
)
print(result.cdr3)  # "AREGTTGKPIGAFAH"
print(result.fr4)   # "WGQGTLVTVSS"
```

By default, sequences with an alignment confidence below `0.5` raise a `ValueError`.
Pass `min_confidence=0.0` to disable this check, or raise the threshold to filter
non-immunoglobulin sequences more aggressively.

# Polars module
Polars module uses [`polars`](https://pola.rs) dataframe library, relying on its built-in multiprocessing and query optimization engine to allow you to run your queries blazingly fast. If you have a big of data (hundreds of thousands and more), it's strongly preferred you use this interface.

It's preferred that you use functional interface for the `polars` module, i.e. call `immunum.polars.number()` and `immunum.polars.segment()` directly instead of calling `immunum.polars.numbering_method(...)`. Those methods provide the same functionality, and will be used in future for runtime-buit annotators, but currently this syntax is just less convenient and requires building annotators in advance.

Also, note that if you don't need a column, don't materialize it afterwards -- on big dataframes, queries will get **significantly** faster. For instance, example below runs in under a second for dataset with 8M entries (200,000 of them are actual paired sequences), while full numbering would take around 10 times more.

## Example: calculating unique cdr3_light/heavy pairs

Segment heavy and light chains in parallel, then deduplicate on their CDR3 sequences.
Note that `.unique()` is called before materializing any other columns — this avoids
loading the full FR region strings for every row.

```python
import polars as pl
import immunum.polars as imp

PARQUET = "test.parquet"

unique_pairs = (
    pl.scan_parquet(PARQUET)
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
print(unique_pairs)

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

## Example: matching with another dataframe on framework sets

For each unique light+heavy framework combination in the source dataset, count how many
sequences in the reference dataset share the same frameworks, then return the 20 most
common groups.

```python
import polars as pl
import immunum.polars as imp

heavy_col = "sequence_alignment_aa_heavy"
light_col = "sequence_alignment_aa_light"
N = 10_000
SOURCE = (
    pl.scan_parquet("test.parquet")
    .filter(
        pl.all_horizontal(
            pl.col(heavy_col).is_not_null(),
            pl.col(light_col).is_not_null(),
        )
    )
    .head(N)
)
TARGET = (
    pl.scan_parquet("validate.parquet")
    .filter(
        pl.all_horizontal(
            pl.col(heavy_col).is_not_null(),
            pl.col(light_col).is_not_null(),
        )
    )
    .head(N)
)

FRAMEWORK_COLS = [
    "fr1_heavy",
    "fr2_heavy",
    "fr3_heavy",
    "fr1_light",
    "fr2_light",
    "fr3_light",
]

top_groups = (
    SOURCE.select(
        imp.segment(heavy_col, chains=["IGH"], scheme="IMGT")
        .name.suffix_fields("_heavy")
        .struct.unnest(),
        imp.segment(light_col, chains=["IGK", "IGL"], scheme="IMGT")
        .name.suffix_fields("_light")
        .struct.unnest(),
    )
    .join(
        TARGET.select(
            imp.segment(heavy_col, chains=["IGH"], scheme="IMGT")
            .name.suffix_fields("_heavy")
            .struct.unnest(),
            imp.segment(light_col, chains=["IGK", "IGL"], scheme="IMGT")
            .name.suffix_fields("_light")
            .struct.unnest(),
        ),
        on=FRAMEWORK_COLS,
        how="inner",
    )
    .group_by(FRAMEWORK_COLS)
    .agg(pl.len().alias("count"))
    .sort("count", descending=True)
    .head(20)
    .collect(engine="streaming")
)
print(top_groups)

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

# Polars integrations

Ecosystem of `polars` has many handy tools for data engineering -- for instance, you can check out [awesome-polars](https://github.com/ddotta/awesome-polars) list. Here we provide examples of two packages that we use actively: [`polars-distance`](https://github.com/ion-elgreco/polars-distance) providing convenient (and fast!) interface for string distance computations, and [`polars-bio`](https://biodatageeks.org/polars-bio/) providing direct interface for reading biological data formats (fasta/fastq/...) into polars dataframes.

## Example: annotating on-the-fly with polars-bio

[polars-bio](https://github.com/biodatageeks/polars-bio) can scan biological sequence files
(FASTA, FASTQ) directly as a `LazyFrame`. Chain that with `imp.segment` to annotate without
an intermediate file:

```python
import polars as pl
import polars_bio as pb
import immunum.polars as imp

print(
    pb.scan_fasta("sequences.fasta")
    .select(
        pl.col("name"),
        imp.segment("sequence", chains=["IGH", "IGK", "IGL"], scheme="IMGT")
        .struct.unnest(),
    )
    .collect(engine='streaming')
)
```

## Example: computing minimal distance between CDR3s in two datasets

Cross-join two CDR3 sets and keep the nearest neighbour by Levenshtein distance.
Pre-filter to unique CDR3s and, if the datasets are large, narrow the search to
equal-length sequences first to reduce the cross-product size.

```python
import polars as pl
import immunum.polars as imp
import polars_distance as pld

result = (
    SOURCE
    .select(
        imp.segment(heavy_col, chains=["IGH"], scheme="IMGT")
        .struct.field("cdr3")
        .alias("cdr3_a")
    )
    .filter(pl.col("cdr3_a").is_not_null())
    .unique("cdr3_a")
    .collect()
    .join(
        TARGET
        .select(
            imp.segment(heavy_col, chains=["IGH"], scheme="IMGT")
            .struct.field("cdr3")
            .alias("cdr3_b")
        )
        .filter(pl.col("cdr3_b").is_not_null())
        .unique("cdr3_b")
        .collect(),
        how="cross",
    )
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
    .collect(engine='streaming')
)
print(result)

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
