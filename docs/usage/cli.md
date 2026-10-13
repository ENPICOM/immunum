# Command line

The `immunum` command numbers sequences and splits them into framework and CDR regions, reading FASTA, plain text or a single sequence and writing TSV, JSON or JSON Lines.

## Installing

`uvx immunum` runs it without installing; `pip install immunum` or `cargo install immunum` installs it. See [Getting started](../getting-started.md#command-line) for details.

## Commands

```text
immunum number  [OPTIONS] [INPUT] [OUTPUT]
immunum segment [OPTIONS] [INPUT] [OUTPUT]
immunum --version
immunum --help
```

- `number` writes one TSV row per numbered residue, or one JSON record per sequence.
- `segment` writes one record per sequence, with a column or field per region.
- `immunum --version` prints the version (`immunum 2.0.0`); `immunum <command> --help` lists the options below.

`INPUT` is a sequence, a file, or `-` for stdin (the default). `OUTPUT` is a file path; without it, results go to stdout.

## Options

Both commands take the same options.

| Option | Values | Default |
| --- | --- | --- |
| `-s`, `--scheme` | `imgt` (`i`), `kabat` (`k`), `chothia` (`c`), `martin` (`m`), `aho` (`a`); case-insensitive | `imgt` |
| `-c`, `--chain` | Chains to try: `h`, `k`, `l`, `a`, `b`, `g`, `d`, or a group: `ig` (H, K, L), `tcr` (A, B, G, D), `all`. Comma-separate several (`-c h,k`). Names such as `heavy` or `igh` also work; see [Chains](../concepts.md#chains) | `ig` |
| `-f`, `--format` | `tsv`, `json`, `jsonl` | `tsv` |
| `--min-confidence` | Number in [0, 1]. A sequence whose alignment confidence is below it gets an error record with `error_kind` `low_confidence` | `0.5` |
| `--all-domains` | Flag. Number or segment every variable domain in each sequence (for example both domains of an scFv), one record per domain, with a 0-based `domain` column or field. With `segment`, every residue lands in exactly one domain's regions. A sequence without a domain gets one error record | off |

## Input

`INPUT` is read as:

| `INPUT` | Read as |
| --- | --- |
| omitted or `-` | stdin |
| an existing path | a file, FASTA or plain text (auto-detected) |
| only ASCII letters, not an existing path | one sequence, with id `seq_1` |
| anything else | an error: `error: cannot open 'EVQL-VES': no such file (a raw sequence may only contain letters)` |

So a mistyped file name such as `sequnces.fasta` is an error, not a sequence.

Files and stdin are auto-detected: if the first non-blank line starts with `>`, the input is FASTA; otherwise it is plain text.

- **FASTA**: the id is the first whitespace-separated word of the header (`>scfv1 my construct` gives `scfv1`). A sequence split over several lines is joined.
- **Plain text**: one sequence per line. Lines are trimmed and blank lines skipped; ids are `seq_1`, `seq_2`, … in order of the remaining lines.

```bash
immunum number EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS
immunum number sequences.fasta
immunum number sequences.txt
cat sequences.fasta | immunum number
immunum number - < sequences.fasta
```

If `INPUT` is omitted and stdin is a terminal, `immunum` exits with `error: no input provided`.

## Output

Results go to stdout, or to `OUTPUT` if given:

```bash
immunum number sequences.fasta results.tsv
immunum number -f json sequences.fasta results.json
```

`json` writes one JSON array of records; `jsonl` writes one JSON object per line. A sequence that can't be numbered gets a record with `error` and `error_kind` set and every other field empty (`null` in JSON), and the run carries on.

### number

TSV has a header row, then one row per numbered residue:

| # | Column | Content |
| --- | --- | --- |
| 1 | `sequence_id` | id of the input sequence |
| 2 | `chain` | chain code, such as `H`, `K`, `A` |
| 3 | `scheme` | `IMGT`, `Kabat`, `Chothia`, `Martin` or `Aho` |
| 4 | `confidence` | alignment confidence, 4 decimals |
| 5 | `position` | position number, with insertion letter (`111A`) |
| 6 | `residue` | amino acid at that position |
| 7 | `error` | empty, or the error message |
| 8 | `error_kind` | empty, or the error kind |

A failed sequence gives one row with only `sequence_id`, `error` and `error_kind` filled.

```text
$ immunum number EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS | head -4
sequence_id	chain	scheme	confidence	position	residue	error	error_kind
seq_1	H	IMGT	0.9344	1	E
seq_1	H	IMGT	0.9344	2	V
seq_1	H	IMGT	0.9344	3	Q

$ immunum number AAAAAAAAAA
sequence_id	chain	scheme	confidence	position	residue	error	error_kind
seq_1						sequence length 10 is below minimum 30	invalid_sequence
warning: 1 of 1 sequences had errors; see the error and error_kind columns
```

JSON and JSONL have one record per sequence:

| Field | Content |
| --- | --- |
| `sequence_id` | id of the input sequence |
| `chain` | chain code |
| `scheme` | scheme name |
| `confidence` | alignment confidence, full precision |
| `numbering` | object mapping each position (`"111A"`) to its residue |
| `query_start` | index in the input of the first numbered residue, 0-based |
| `query_end` | index of the last numbered residue, 0-based inclusive |
| `error` | `null`, or the error message |
| `error_kind` | `null`, or the error kind |

```text
$ immunum number -f jsonl EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS
{"sequence_id":"seq_1","chain":"H","scheme":"IMGT","confidence":0.93440181016922,"numbering":{"1":"E","2":"V","3":"Q",…},"query_start":0,"query_end":34,"error":null,"error_kind":null}
```

### segment

TSV and JSON alike have one record per sequence, with these columns or fields in order: `sequence_id`, `prefix`, `fr1`, `cdr1`, `fr2`, `cdr2`, `fr3`, `cdr3`, `fr4`, `postfix`, `error`, `error_kind`. `prefix` and `postfix` hold residues before and after the numbered domain; see [Prefix and postfix](../concepts.md#prefix-and-postfix).

```text
$ immunum segment fixtures/ig.fasta | head -3
sequence_id	prefix	fr1	cdr1	fr2	cdr2	fr3	cdr3	fr4	postfix	error	error_kind
4qo1_A|Heavy|A		VQLQESGGGLVQAGGSLRLSCAAS	ERTFSTYA	MGWFRQAPGREREFLAQ	INWSGTTT	YYAESVKDRTTISRDNAKNTVYLEMNNLNADDTGIYFC	AAHPQRGWGSTLGWTY	WGQGTQVTVSS
4olu_L|Light|L		VLTQSPGTLSLSPGETAIISCRTS	QYGS	LAWYQQRPGQAPRLVIY	SGS	TRAAGIPDRFSGSRWGPDYTLTISNLESGDFGVYYC	QQYEF	FGQGTKVQVD

$ immunum segment -f jsonl fixtures/ig.fasta | head -1
{"sequence_id":"4qo1_A|Heavy|A","prefix":"","fr1":"VQLQESGGGLVQAGGSLRLSCAAS","cdr1":"ERTFSTYA","fr2":"MGWFRQAPGREREFLAQ","cdr2":"INWSGTTT","fr3":"YYAESVKDRTTISRDNAKNTVYLEMNNLNADDTGIYFC","cdr3":"AAHPQRGWGSTLGWTY","fr4":"WGQGTQVTVSS","postfix":"","error":null,"error_kind":null}
```

### With --all-domains

`--all-domains` inserts a `domain` column (0-based) as column 2 of the TSV output of both commands, and a `domain` field after `sequence_id` in JSON. **Every later TSV column moves one to the right**: in `number` output, `chain` is `$3` and `position` is `$6` instead of `$5` in awk. An error record leaves `domain` empty (`null` in JSON).

```text
$ immunum number --all-domains scfvs.fasta | head -2
sequence_id	domain	chain	scheme	confidence	position	residue	error	error_kind
scfv1	0	H	IMGT	0.7012	2	V

$ immunum segment --all-domains scfvs.fasta | cut -f1,2,9
sequence_id	domain	cdr3
scfv1	0	AAHPQRGWGSTLGWTY
scfv1	1	QQYEF
```

See [Multiple domains](../concepts.md#multiple-domains).

## Exit status

| Status | When | stderr |
| --- | --- | --- |
| 0 | Success, **even if some sequences failed**: they get error records | `warning: N of M sequences had errors; see the error and error_kind columns` when any failed |
| 0 | The reader closed the pipe early, as in `immunum number seqs.fasta \| head` | nothing |
| 1 | Setup or input error: unknown scheme, chain or format, `--min-confidence` outside [0, 1], missing input file, output file can't be created, no input | `error: …`, for example `error: unknown chain 'q' (options: …)` |
| 2 | Usage error: unknown option, invalid option value (such as `--min-confidence x`) | clap's usage message |

To detect failed sequences, check the `error_kind` column or field, not the exit status.

## Examples

```bash
# Kabat scheme, JSON output
immunum number -s kabat -f json EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMS

# All chains (antibody + TCR), JSONL output
immunum number -c all -f jsonl sequences.fasta

# TCR sequences only, save to file
immunum number -c tcr tcr_sequences.fasta output.tsv

# Both domains of each scFv, one JSONL record per domain
immunum number --all-domains -f jsonl scfvs.fasta

# FR/CDR regions of each sequence, one TSV row per sequence
immunum segment sequences.fasta

# CDR3 of every domain of each scFv
immunum segment --all-domains scfvs.fasta | cut -f1,2,9

# Extract sequences from a TSV column and pipe in (see fixtures/ig.tsv)
tail -n +2 fixtures/ig.tsv | cut -f2 | immunum number
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="sequence") c=i} NR>1{print $c}' fixtures/ig.tsv | immunum number

# Filter TSV output to CDR3 positions (105-117 in IMGT); `+0` reads insertions like 111A as 111
immunum number sequences.fasta | awk -F'\t' 'NR==1 || ($5+0 >= 105 && $5+0 <= 117)'

# The same with --all-domains: position is now column 6
immunum number --all-domains scfvs.fasta | awk -F'\t' 'NR==1 || ($6+0 >= 105 && $6+0 <= 117)'

# Filter to heavy chain results only
immunum number -c all sequences.fasta | awk -F'\t' 'NR==1 || $2=="H"'

# CDR3 of each sequence with jq
immunum segment -f json sequences.fasta | jq -r '.[] | [.sequence_id, .cdr3] | @tsv'
```

## Errors

A sequence that can't be numbered gets an `error` message and an `error_kind`; [Error kinds](errors.md#error-kinds) lists the kinds and what causes each.
