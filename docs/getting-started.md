# Getting started

## Which interface should I use?

Every interface runs the same Rust core and gives the same numbering; they differ in how you feed
it sequences.

| You have | Use | Why |
|---|---|---|
| A FASTA file, or sequences in a shell pipeline | [Command line](usage/cli.md) | No code: TSV or JSON out, one record per residue or per sequence |
| A few sequences up to tens of thousands, in Python | [Python](usage/python.md) (`immunum.Annotator`) | Plain objects per sequence; thread-safe and picklable, so it fits thread and process pools |
| Hundreds of thousands of sequences or more, or data already in a DataFrame | [Polars](usage/polars.md) (`immunum.polars`) | Runs inside a Polars query on every core, streams larger-than-memory files, and only computes the columns you keep |
| A Node.js service or tool | [JavaScript](usage/javascript.md) (npm `immunum`) | WebAssembly build of the core with TypeScript types |
| A Rust program | [Rust](usage/rust.md) (crate `immunum`) | The core library itself, with no Python or JavaScript runtime |
| One sequence to look at | [Demo](demo.md) | Runs in your browser; nothing to install |

--8<-- "README.md:install"

## Next steps

- [Concepts and FAQ](concepts.md): which chains and schemes to pick, what confidence means, and
  how flanking residues are handled.
- The usage guide for your interface, from the table above.
- [Errors](usage/errors.md): what's raised and what's returned when a sequence can't be numbered.
- Coming from immunum 1.x? Read [Upgrading to 2.0](usage/upgrading.md).
