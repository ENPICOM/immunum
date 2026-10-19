# Rust

The `immunum` crate is the core every other interface wraps. Add it with `cargo add immunum`
([Getting started](../getting-started.md#rust)); the full API reference is on
[docs.rs/immunum](https://docs.rs/immunum).

## Features

| Feature | Default | What it adds |
|---|---|---|
| `cli` | yes | The `immunum` command-line tool (`immunum::cli`) and its `clap` dependency |
| `python`, `polars`, `wasm` | no | The Python, Polars and WebAssembly bindings; used to build those packages |
| `dev-tools` | no | The validation tools in `src/bin`, run from a checkout |

As a library dependency, leave out the CLI: `cargo add immunum --no-default-features`.

## Numbering and segmentation

```rust
use immunum::{Annotator, Chain, Scheme};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Chains to consider, the scheme, and the minimum confidence (None: 0.4)
    let annotator = Annotator::new(&[Chain::IGH, Chain::IGK, Chain::IGL], Scheme::IMGT, None)?;
    // Or from names as a user writes them, groups included
    let same = Annotator::from_names(["ig"], "imgt", None)?;

    let sequence = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";

    let result = annotator.number(sequence)?;
    println!("{} {:.2}", result.chain, result.confidence); // H 0.78
    for (position, residue) in result.residues(sequence)? {
        println!("{position} {residue}"); // 1 Q, 2 V, ... 111A K, ...
    }

    let segments = same.segment(sequence)?;
    println!("{}", segments.cdr3); // AREGTTGKPIGAFAH
    for (name, residues) in segments.regions() {
        println!("{name}: {residues}"); // prefix, fr1, cdr1, ... fr4, postfix in sequence order
    }
    Ok(())
}
```

`NumberingResult::segment(sequence)` splits a numbering you already have into regions without
aligning again.

## Multiple domains

`number_domains` and `segment_domains` return one result per variable domain, in sequence order,
such as both domains of an scFv. `domains` returns the found domains themselves; `Domain::number`
numbers one under any scheme from the alignment that found it, without aligning again:

```rust
use immunum::{Annotator, Scheme};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";
    let kappa = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK";
    let scfv = format!("{heavy}GGGGSGGGGSGGGGS{kappa}");

    let annotator = Annotator::from_names(["ig"], "imgt", None)?;
    for domain in annotator.domains(&scfv)? {
        let kabat = domain.number(Scheme::Kabat)?;
        println!("{} at {}..={}", domain.chain(), domain.query_start(), domain.query_end());
        println!("Kabat CDR3: {}", kabat.segment(&scfv)?.cdr3);
    }
    Ok(())
}
```

How domains are found, and how the residues between them are assigned, is in
[Multiple domains](../concepts.md#multiple-domains).

## Sharing an annotator

`Annotator` is `Send`, `Sync` and `Clone`, so one can be shared between threads (for example behind
an `Arc`, or by reference in `std::thread::scope` or rayon). It also implements serde's
`Serialize` and `Deserialize`.

## Errors

`Annotator::new` and `from_names` return `immunum::Error`; numbering one sequence returns
`immunum::SequenceError`. Both have `kind()`, the code every other interface reports. See
[Errors](errors.md#rust).
