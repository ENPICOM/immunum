# Concepts and FAQ

How immunum picks a chain, numbers it under each scheme, scores the result and handles input.
The examples use Python; every interface behaves the same way.

## Chains

immunum numbers seven chains. Results report a chain by its one-letter code. Every name in the
table is accepted, in any case; spaces around a name are ignored, and a chain named twice is
kept once.

| Chain             | Code | Also accepted   |
|-------------------|------|-----------------|
| Antibody heavy    | `H`  | `IGH`, `heavy`  |
| Antibody kappa    | `K`  | `IGK`, `kappa`  |
| Antibody lambda   | `L`  | `IGL`, `lambda` |
| T-cell receptor α | `A`  | `TRA`, `alpha`  |
| T-cell receptor β | `B`  | `TRB`, `beta`   |
| T-cell receptor γ | `G`  | `TRG`, `gamma`  |
| T-cell receptor δ | `D`  | `TRD`, `delta`  |

Wherever a list of chains is expected, a group can stand in for several chains:

| Group | Chains             |
|-------|--------------------|
| `ig`  | H, K, L            |
| `tcr` | A, B, G, D         |
| `all` | H, K, L, A, B, G, D |

Each list entry is one name: `["ig", "B"]` asks for H, K, L and B, while `"igh,igk"` is rejected
as an unknown chain. (The command line splits its `--chain` value on commas itself; see
[Command line](usage/cli.md).) Groups are only accepted where several chains are, not where one
chain is expected, such as `regions_for`.

### Chain detection

The annotator aligns the sequence to the consensus profile of every chain you asked for and keeps
the chain whose alignment scores highest. The [confidence](#confidence) is then computed for that
alignment. Each requested chain is one more alignment per sequence, so time grows roughly in
proportion to the number of chains: `["all"]` takes about five times as long as `["H"]`. Pass only the
chains you expect in your data.

## Numbering schemes

| Scheme  | Accepted names (any case) | Reported as | Chains          |
|---------|---------------------------|-------------|-----------------|
| IMGT    | `imgt`, `i`               | `IMGT`      | all seven       |
| Kabat   | `kabat`, `k`              | `Kabat`     | H, K, L         |
| Chothia | `chothia`, `c`            | `Chothia`   | H, K, L         |
| Martin  | `martin`, `m`             | `Martin`    | H, K, L         |
| AHo     | `aho`, `a`                | `Aho`       | H, K, L         |

Only IMGT numbers T-cell receptor chains. Asking another scheme for one raises an
`unsupported_chain` error when the annotator is built; `scheme_supports_chain(scheme, chain)`
tells you beforehand.

Every sequence is aligned once against IMGT-numbered consensus profiles (IMGT positions 1–128).
The other schemes are derived from that alignment by per-scheme rules that map IMGT positions to
the scheme's positions. The detected chain, the confidence and the residues numbered are
therefore the same under every scheme; only the position labels and region boundaries change.

### Region boundaries

CDR boundaries per scheme, inclusive, from `immunum.regions_for(scheme, chain)`. The framework
regions fill the gaps: FR1 runs from 1 to just before CDR1, FR4 from just after CDR3 to the last
position.

| Scheme  | Chain    | CDR1    | CDR2  | CDR3    | Positions |
|---------|----------|---------|-------|---------|-----------|
| IMGT    | any      | 27–38   | 56–65 | 105–117 | 1–128     |
| AHo     | H, K, L  | 25–42   | 57–77 | 109–138 | 1–149     |
| Kabat   | H        | 31–35   | 50–65 | 95–102  | 1–113     |
| Kabat   | K, L     | 24–34   | 50–56 | 89–97   | 1–107     |
| Chothia | H        | 26–32   | 52–56 | 96–101  | 1–113     |
| Chothia | K, L     | 26–32   | 50–52 | 91–96   | 1–107     |
| Martin  | H        | 26–35   | 50–58 | 95–102  | 1–113     |
| Martin  | K, L     | 24–34   | 50–56 | 89–97   | 1–107     |

IMGT and AHo use the same positions for every chain, so one number means the same place on a
heavy and a light chain. Kabat, Chothia and Martin define heavy and light CDRs differently.
Martin's CDRs follow the AbM definition.

The same sequence therefore gets different CDRs. For this heavy chain:

| Scheme  | CDR-H1         | CDR-H2               | CDR-H3            |
|---------|----------------|----------------------|-------------------|
| IMGT    | `GGSFSTYA`     | `VIPLLTIT`           | `AREGTTGKPIGAFAH` |
| Kabat   | `TYALS`        | `GVIPLLTITNYAPRFQG`  | `EGTTGKPIGAFAH`   |
| Chothia | `GGSFSTY`      | `IPLLTI`             | `GTTGKPIGAFA`     |
| Martin  | `GGSFSTYALS`   | `GVIPLLTITN`         | `EGTTGKPIGAFAH`   |
| AHo     | `ASGGSFSTYALS` | `GVIPLLTITNYAPRFQGR` | `EGTTGKPIGAFAH`   |

```python
from immunum import Annotator, regions_for

heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"

assert regions_for("kabat", "H")["cdr1"] == (31, 35)
assert regions_for("chothia", "H")["cdr1"] == (26, 32)
assert regions_for("kabat", "K")["cdr1"] == (24, 34)
assert regions_for("imgt", "H") == regions_for("imgt", "L")

kabat = Annotator(chains=["H"], scheme="kabat").segment(heavy)
chothia = Annotator(chains=["H"], scheme="chothia").segment(heavy)
assert kabat.cdr1 == "TYALS"
assert chothia.cdr1 == "GGSFSTY"
```

### Insertions and deletions

A loop longer than a scheme's positions gets insertion letters; positions are strings for that
reason, such as `"111A"`. Results list positions in sequence order.

- **IMGT** splits insertions between the two centre positions of a loop: 32/33 in CDR1, 60/61 in
  CDR2 and 111/112 in CDR3. The heavy chain above has a 15-residue CDR3, numbered with `111A`
  and `112A`; four residues more give `111A`, `111B`, `111C`, `112C`, `112B`, `112A`. A loop
  shorter than its positions leaves out the centre ones first: the 6-residue CDR1 `QDVNTA` of
  the kappa chain below is numbered 27, 28, 29, 36, 37, 38.
- **Kabat** and **Chothia** put every insertion after one position. For the heavy chain above
  both give `52A`, `82A`–`82C` and `100A`–`100E`.
- **Martin** heavy chains are numbered like Chothia, except that framework 3 insertions go on 72
  instead of 82: `72A`–`72C` for the heavy chain above.
- **AHo** has positions for longer loops (1–149); none of the sequences on this page get an
  insertion letter under AHo.

```python
from immunum import Annotator

heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
kappa = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"


def insertions(scheme, sequence):
    numbering = Annotator(chains=["ig"], scheme=scheme).number(sequence).numbering
    return [position for position in numbering if position[-1].isalpha()]


assert insertions("imgt", heavy) == ["111A", "112A"]
assert insertions("kabat", heavy) == [
    "52A", "82A", "82B", "82C", "100A", "100B", "100C", "100D", "100E"
]
assert insertions("martin", heavy) == [
    "52A", "72A", "72B", "72C", "100A", "100B", "100C", "100D", "100E"
]

kappa_numbering = Annotator(chains=["ig"], scheme="imgt").number(kappa).numbering
cdr1 = [p for p in kappa_numbering if 27 <= int(p) <= 38]
assert cdr1 == ["27", "28", "29", "36", "37", "38"]
```

## Confidence

Every result carries a `confidence` between 0 and 1: how well the sequence matches the consensus
profile of the detected chain. It is computed from the alignment:

1. Each consensus position scores every amino acid by BLOSUM62, averaged over the residues seen
   at that position in the chain's reference sequences and weighted by their frequency. Letters
   other than the 20 standard amino acids score −4.
2. Only positions occupied in more than 90% of the reference sequences count. In the heavy-chain
   profile that is every framework position but only 21 of the 35 CDR positions: the centre of
   each loop, where shorter loops have no residue, doesn't count.
3. Over the counted positions between the first and last consensus position the alignment
   reaches, immunum sums the score of the residue aligned there, or the position's gap penalty
   when the alignment skips it.
4. That sum is divided by the best possible sum over the same positions and clamped to 0–1.

Inserted residues, and residues before or after the domain, don't enter the score. Because only
the reached positions count, a fragment can score as high as a full domain: the first 30
residues of the heavy chain above score 0.77, the whole domain 0.78. Changing CDR3 residues
barely moves the score; changing a conserved framework residue does.

### Minimum confidence

An annotator has a `min_confidence`, 0.5 unless you pass another value in `[0, 1]`. When the best
alignment of a sequence scores below it, the sequence isn't numbered: its result has `error_kind`
`low_confidence` and the score in `error`. In the [domain search](#multiple-domains) it also
decides when to stop looking for further domains.

```python
from immunum import Annotator

heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
gfp = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"

annotator = Annotator(chains=["ig"], scheme="imgt")
print(round(annotator.number(heavy).confidence, 3))  # 0.785

result = annotator.number(gfp)
print(result.error_kind)  # low_confidence
print(result.error)  # alignment confidence 0.1525 is below min_confidence 0.5000

unfiltered = Annotator(chains=["ig"], scheme="imgt", min_confidence=0.0)
print(round(unfiltered.number(gfp).confidence, 3))  # 0.153
```

### Choosing a value

Scores of the bundled validation sequences (`fixtures/validation`), numbered with
`chains=["all"]`:

| Chain | Sequences | Lowest | Median | Below 0.5 |
|-------|-----------|--------|--------|-----------|
| H     | 2452      | 0.550  | 0.803  | 0         |
| K     | 1485      | 0.583  | 0.850  | 0         |
| L     | 371       | 0.553  | 0.821  | 0         |
| A     | 865       | 0.285  | 0.645  | 77        |
| B     | 934       | 0.269  | 0.765  | 57        |
| G     | 25        | 0.567  | 0.849  | 0         |
| D     | 23        | 0.766  | 0.825  | 0         |

Green fluorescent protein scores 0.153 and hen egg-white lysozyme 0.082 against the antibody
chains.

- Keep 0.5 for antibodies: it kept every antibody in the validation sets.
- Lower it for TCR α and β chains, where the default rejects some real sequences: at 0.25 every
  TRA and TRB validation sequence passes.
- Raise it to reject sequences that align poorly, such as heavily mutated or engineered domains,
  at the cost of losing some real ones.
- `0.0` disables the filter. Every valid sequence then gets numbered, even one that isn't an
  antibody: GFP gets four residues numbered as a kappa chain.

## Prefix and postfix

Segmentation accounts for every residue: joining `prefix`, `fr1` … `fr4` and `postfix` in order
gives back the input sequence.

- `prefix` holds the residues before FR1, such as a signal peptide, or the linker after a
  previous domain.
- `postfix` holds the residues after FR4, such as the start of a constant region.

Both are often empty. `number` leaves these residues out of `numbering`; `query_start` and
`query_end` (0-based, inclusive) give the span it numbered.

```python
from immunum import Annotator

leader = "MEFGLSWVFLVALLRGVQC"
heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
constant = "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV"
sequence = leader + heavy + constant

annotator = Annotator(chains=["ig"], scheme="imgt")

regions = annotator.segment(sequence)
assert regions.prefix == leader
assert regions.fr1 == "QVQLVQSGAEVKRPGSSVTVSCKAS"
assert regions.postfix == constant
assert "".join(regions.as_dict().values()) == sequence

result = annotator.number(sequence)
assert (result.query_start, result.query_end) == (19, 140)
assert sequence[result.query_start : result.query_end + 1] == heavy
```

With several domains (`segment_domains`), every residue still lands in exactly one domain's
regions: a domain's `prefix` holds the residues since the previous domain (or the start of the
sequence), and only the last domain has a `postfix`. The domains' regions in order rebuild the
sequence.

## Multiple domains

`number` and `segment` return one domain: the best-scoring alignment of the whole sequence. For
an scFv with a heavy and a kappa domain, `number` returns the kappa domain, because it scores
higher. `number_domains` and `segment_domains` return every domain, ordered by position. On the
command line, `--all-domains` does the same.

The domain search:

1. Align the whole sequence and keep the best alignment as a domain. If it scores below
   `min_confidence`, the result is `low_confidence`, as for `number`.
2. Search the residues before that domain and the residues after it the same way, each on its
   own, and keep each part's best alignment if it reaches `min_confidence`.
3. Repeat on the residues around every domain found. A part shorter than 30 residues isn't
   searched.
4. A domain shorter than 30 residues isn't reported, but the residues around it are still
   searched.

The lists are never empty. When no domain is found, the list holds a single result with `error`
set: `low_confidence` when the best alignment scores too low, or `domain_too_short` when it is
confident but shorter than 30 residues. An invalid sequence gives a single `invalid_sequence`
result. `number` has no length minimum for the domain it returns, so a fragment can be numbered
by `number` and still give `domain_too_short` from `number_domains`.

```python
from immunum import Annotator

heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
kappa = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
linker = "GGGGSGGGGSGGGGS"
scfv = heavy + linker + kappa

annotator = Annotator(chains=["ig"], scheme="imgt")

assert annotator.number(scfv).chain == "K"
domains = annotator.number_domains(scfv)
assert [d.chain for d in domains] == ["H", "K"]
assert domains[1].query_start == len(heavy + linker)

heavy_regions, kappa_regions = annotator.segment_domains(scfv)
assert heavy_regions.postfix == ""
assert kappa_regions.prefix == linker

# The last 28 residues of the heavy chain, then 40 glycines
fragment = heavy[-28:] + "G" * 40
assert annotator.number(fragment).error is None
(result,) = annotator.number_domains(fragment)
assert result.error_kind == "domain_too_short"
print(result.error)  # domain length 28 is below minimum 30
```

## Input rules

- **Length**: 30 to 10,000 characters. Shorter or longer sequences get an `invalid_sequence`
  error.
- **Characters**: every character must be an ASCII letter, `A`–`Z` or `a`–`z`. Digits, spaces,
  newlines, `*`, `-`, `.` and non-ASCII characters give an `invalid_sequence` error that names the
  character and its 0-based position.
- **Non-standard letters** such as `X`, `B` or `Z` are accepted and numbered, but score as a
  mismatch (−4), so they lower the [confidence](#confidence) at counted positions: one `X` at
  IMGT position 6 of the heavy chain above lowers it from 0.785 to 0.767.
- **Case**: lower and upper case score the same, and residues come back as you gave them.
- **Whitespace**: the Python, JavaScript and Rust APIs take the sequence as given, so a trailing
  newline or a space inside the sequence is an invalid character; strip it first. The command
  line trims each input line and joins the lines of a FASTA record; whitespace inside a line is
  still an error. The [demo](demo.md) removes all whitespace before numbering.

```python
from immunum import Annotator

heavy = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
annotator = Annotator(chains=["H"], scheme="imgt")

assert annotator.number(heavy.lower()).numbering["1"] == "q"

result = annotator.number(heavy + "\n")
print(result.error_kind)  # invalid_sequence
print(result.error)  # invalid character '\n' at position 122

assert annotator.number("".join((heavy + "\n").split())).error is None

with_x = annotator.number(heavy[:5] + "X" + heavy[6:])
assert with_x.numbering["6"] == "X"

print(annotator.number(heavy[:29]).error)  # sequence length 29 is below minimum 30
```

See [Errors](usage/errors.md#error-kinds) for every error kind.

## FAQ

### Why is my sequence `low_confidence`?

Its best alignment scored below `min_confidence` (0.5 by default). The `error` message gives the
score. Common causes:

- It isn't an antibody or TCR variable domain.
- You didn't ask for its chain. The TCR β chain 1AO7_E from the validation set scores 0.30
  against `["ig"]` and 0.76 against `["tcr"]`.
- It is a TCR α or β chain that scores low; see [Choosing a value](#choosing-a-value).
- It is a nucleotide sequence; see [below](#how-do-i-number-nucleotide-sequences).
- It holds many non-standard letters such as `X` at conserved positions.

To see the score of a rejected sequence, number it with `min_confidence=0.0`.

### Which scheme should I use?

Use the scheme your downstream tools or data expect. Otherwise:

- **IMGT** is the only one that numbers TCR chains, uses the same positions for every chain, and
  is what immunum aligns to internally.
- **AHo** also uses the same positions for every chain, with room for longer loops (1–149).
- **Kabat**, **Chothia** and **Martin** have separate heavy and light layouts and differ mostly in
  where the CDRs start and end; see [Region boundaries](#region-boundaries).

The numbered residues and the confidence don't depend on the scheme, so you can number a sequence
under several schemes and compare.

### Why does my numbering differ from ANARCI for some sequences?

immunum doesn't run ANARCI. It aligns to its own consensus profiles and derives every scheme from
that IMGT alignment, so where a numbering is ambiguous it can choose differently. The
[benchmarks](benchmarks/index.md) measure per-region agreement with the validation data, which
for antibodies is where ANARCI and antpack agree.

### Is the Annotator thread-safe? Can I pickle it?

Yes. In Python an `Annotator` can be shared between threads and pickled with
its settings, so it works with thread and process pools; see
[Batches, threads and processes](usage/python.md#batches-threads-and-processes). In Rust it is
`Send`, `Sync` and `Clone`; see [Sharing an annotator](usage/rust.md#sharing-an-annotator).

### How do I number nucleotide sequences?

Translate them to amino acids first, in the right frame; immunum numbers protein sequences only.
A DNA sequence is all letters, so it isn't rejected as invalid, but it is scored as amino acids
and comes back `low_confidence`. Remove stop codons (`*`) from the translation: `*` is an invalid
character.

```python
from immunum import Annotator

# The first 147 nucleotides of a human IGHV1 gene
dna = "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGTCCTCGGTGAAGGTCTCCTGCAAGGCTTCTGGAGGCACCTTCAGCAGCTATGCTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGA"

result = Annotator(chains=["ig"], scheme="imgt").number(dna)
assert result.error_kind == "low_confidence"
```
