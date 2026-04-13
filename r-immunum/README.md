# immunum (R)

R bindings for the [`immunum`](https://github.com/ENPICOM/immunum) Rust crate:
fast numbering of antibody and T-cell receptor variable-domain sequences using
IMGT and Kabat schemes.

## Installation

```r
# install.packages("remotes")
remotes::install_github("ENPICOM/immunum", subdir = "r-immunum", build = FALSE)
```

## Build from source

You need a working R installation, the [Rust toolchain](https://rustup.rs/),
and a C compiler. From the repository root:

```bash
R CMD INSTALL r-immunum
```

Or via the project Taskfile:

```bash
task r:build
```

## Quick start

```r
library(immunum)

ann <- Annotator$new(chains = c("H", "K", "L"), scheme = "imgt")

result <- ann$number(
  "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
)
result$chain        # "H"
result$confidence   # 0.78
result$numbering    # named character vector "1"="Q", "2"="V", ...
```

See the package vignettes (`browseVignettes("immunum")`) for batch processing,
the polars integration, and benchmarks.

## License

MIT, copyright 2026 ENPICOM. See `LICENSE`.
