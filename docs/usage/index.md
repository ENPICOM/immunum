# Usage

Task guides for each interface. They share one model:

- **Number** a sequence: give every residue a position in the scheme you chose (`number`).
- **Segment** it: split it into FR1–FR4 and CDR1–CDR3, plus the residues before and after the
  domain as `prefix` and `postfix`, so the pieces always join back into the input (`segment`).
- Do either for **every domain** in a sequence, such as both domains of an scFv
  (`number_domains`, `segment_domains`, or `--all-domains` on the command line).

You set the chains to consider, the scheme and the minimum confidence once, when you build an
annotator; [Concepts](../concepts.md) explains each choice. A sequence that can't be numbered gives
a result with `error` set rather than stopping the batch; see [Errors](errors.md).

| Guide | For |
|---|---|
| [Command line](cli.md) | FASTA or plain-text files, shell pipelines; TSV and JSON output |
| [Python](python.md) | `immunum.Annotator`, one sequence at a time |
| [Polars](polars.md) | `immunum.polars`, large datasets inside Polars queries |
| [JavaScript](javascript.md) | The npm package for Node.js |
| [Rust](rust.md) | The `immunum` crate |
| [Errors](errors.md) | What every interface raises or returns when something goes wrong |
| [Upgrading to 2.0](upgrading.md) | What can break when coming from 1.x |

Installing: [Getting started](../getting-started.md). Every parameter and field: the
[API reference](../api/index.md).
