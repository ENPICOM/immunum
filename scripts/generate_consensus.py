"""Generate resources/consensus/*.csv from IMGT/GENE-DB germline reference sequences.

Usage:
    uv run python -m scripts.generate_consensus [--cache-dir .cache/imgt] [--output-dir resources/consensus]

Downloads the IMGT/GENE-DB amino acid reference sequences with IMGT gaps (CC BY 4.0) into the
cache directory, builds every V x J germline combination per species and chain, and writes
per-position amino acid frequencies on the IMGT unique numbering (positions 1-128). The method is
described in resources/consensus/README.md; the alleles used are listed in
resources/consensus/germlines.tsv.
"""

import argparse
import re
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path
from typing import NamedTuple

IMGT_DOWNLOAD_URL = "https://www.imgt.org/download/GENE-DB/"
IMGT_FASTA = "IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP"
IMGT_RELEASE = "RELEASE"

# Species are named as in the IMGT/GENE-DB FASTA header. The antibody set is the one ANARCI builds
# its HMMs from; TCR germlines are only taken from human and mouse.
IG_SPECIES = (
    "Homo sapiens",
    "Mus musculus",
    "Rattus norvegicus",
    "Oryctolagus cuniculus",
    "Macaca mulatta",
    "Sus scrofa",
    "Vicugna pacos",
    "Bos taurus",
)
TR_SPECIES = ("Homo sapiens", "Mus musculus")
# Each species' V alignment is mapped onto the human one of the same V group (see map_columns).
REFERENCE_SPECIES = "Homo sapiens"

CHAINS = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
GAP = "."

# IMGT V-domain landmarks: 1st-CYS 23 and 2nd-CYS 104 (1-based)
FIRST_CYS, SECOND_CYS = 23, 104
# J-TRP/J-PHE: the F or W of the J motif [FW]G.G sits at IMGT 118; FR4 runs 118-128.
J_ANCHOR_POSITION, LAST_POSITION = 118, 128
J_MOTIF = re.compile(r"[FW]G.G")
# Fallback for J genes with a non-canonical motif, e.g. human TRAJ16 FARG, TRGJP1 FAEG
J_MOTIF_FALLBACK = re.compile(r"[FW]..G")
# CDR3-IMGT is 105-117. A CDR3 shorter than 13 residues leaves these positions empty, in order.
CDR3_FIRST, CDR3_LAST = 105, 117
CDR3_DELETION_ORDER = (111, 112, 110, 113, 109, 114, 108, 115, 107, 116, 106, 117)

# Positions the Kabat-derived schemes leave unnumbered on antibody chains, kept out of the
# antibody consensus so an alignment never places a residue there (see src/numbering/*.rs).
SCHEME_HOLES = {"IGH": (10, 73), "IGK": (73,), "IGL": (10, 73)}

# Profile alignment weights for map_columns
OCCUPANCY_MISMATCH_WEIGHT = 4.0
DELETION_WEIGHT = 8.0

MIN_FREQUENCY = 0.05

REGION_BOUNDS = (
    (26, "FR1"),
    (38, "CDR1"),
    (55, "FR2"),
    (65, "CDR2"),
    (104, "FR3"),
    (117, "CDR3"),
)

BLOSUM62 = (
    (4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2),
    (0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2),
    (-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3),
    (-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2),
    (-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3),
    (0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3),
    (-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2),
    (-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1),
    (-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2),
    (-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1),
    (-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1),
    (-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2),
    (-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3),
    (-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1),
    (-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2),
    (1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2),
    (0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2),
    (0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1),
    (-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2),
    (-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7),
)
AA_INDEX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}


class Allele(NamedTuple):
    name: str  # e.g. IGHV1-18*01
    species: str
    sequence: str  # V: IMGT-gapped; J: ungapped

    @property
    def gene(self) -> str:
        return self.name.split("*")[0]

    @property
    def allele_number(self) -> int:
        return int(self.name.split("*")[1])


class VGene(NamedTuple):
    allele: str
    framework: str  # IMGT positions 1-104, GAP where empty
    cdr3: str  # V-encoded residues after 2nd-CYS


class JGene(NamedTuple):
    allele: str
    cdr3: str  # J-encoded residues before J-TRP/J-PHE
    framework: str  # residues from IMGT 118 onward


def download_imgt(cache_dir: Path) -> tuple[Path, str]:
    """Fetch the IMGT/GENE-DB FASTA and its release tag unless already cached."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    for name in (IMGT_RELEASE, IMGT_FASTA):
        path = cache_dir / name
        if not path.exists():
            url = IMGT_DOWNLOAD_URL + urllib.parse.quote(name)
            print(f"Downloading {url}")
            urllib.request.urlretrieve(url, path)
    release = (cache_dir / IMGT_RELEASE).read_text().strip()
    return cache_dir / IMGT_FASTA, release


def read_alleles(fasta: Path, region: str, species: set[str]) -> list[Allele]:
    """Functional, complete alleles of one IMGT label (V-REGION or J-REGION)."""
    alleles = []
    header, chunks = None, []
    with open(fasta) as f:
        for line in [*f, ">"]:
            line = line.strip()
            if not line.startswith(">"):
                chunks.append(line)
                continue
            if header is not None:
                # 15 '|'-separated fields, see IMGT/GENE-DB README.txt
                fields = header.split("|")
                allele_species = fields[2].split("_")[0]  # drop strain/breed suffix
                sequence = "".join(chunks)
                if (
                    fields[4] == region
                    and fields[3] == "F"
                    and not fields[13].strip()  # not partial
                    and allele_species in species
                    and set(sequence) <= set(AMINO_ACIDS + GAP)
                ):
                    alleles.append(Allele(fields[1], allele_species, sequence))
            header, chunks = line[1:], []
    return alleles


def species_for(chain: str) -> tuple[str, ...]:
    return IG_SPECIES if chain.startswith("IG") else TR_SPECIES


def column_profile(
    sequences: list[str], length: int
) -> tuple[list[dict[int, float]], list[float]]:
    """Per column: amino acid frequencies among residues, and fraction of sequences occupied."""
    counts = [Counter() for _ in range(length)]
    for seq in sequences:
        for i, c in enumerate(seq[:length]):
            if c != GAP:
                counts[i][AA_INDEX[c]] += 1
    freqs, occupancy = [], []
    for column in counts:
        total = sum(column.values())
        freqs.append({aa: n / total for aa, n in column.items()})
        occupancy.append(total / len(sequences))
    return freqs, occupancy


def map_columns(
    sequences: list[str], reference: list[str]
) -> tuple[dict[int, int], list[int]]:
    """Map the columns of one species' IMGT-gapped V alignment onto IMGT positions.

    IMGT/GENE-DB gaps every species' V genes in one alignment per V group, but where an allele
    carries an insertion IMGT adds a column for it to the whole group (rhesus IGHV, mouse TRAV,
    rat IGHV, ...), shifting every later position. Profile-align the group to the human alignment
    of the same V group: every reference position 1-104 is matched, extra columns are deleted,
    and the columns after the one matched to 104 become 105, 106, ...

    Returns (0-based column -> IMGT position, deleted 0-based columns).
    """
    n_cols = max(map(len, sequences))
    freqs, occupancy = column_profile(sequences, n_cols)
    ref_freqs, ref_occupancy = column_profile(reference, SECOND_CYS)

    def match(i: int, j: int) -> float:
        expected = sum(
            fa * fb * BLOSUM62[a][b]
            for a, fa in freqs[i].items()
            for b, fb in ref_freqs[j].items()
        )
        o, ro = occupancy[i], ref_occupancy[j]
        return expected * min(o, ro) - OCCUPANCY_MISMATCH_WEIGHT * abs(o - ro)

    def delete(i: int) -> float:
        # Small constant: of two equal alignments, prefer deleting the later column
        return -DELETION_WEIGHT * occupancy[i] - 0.01

    # score[i][j]: best alignment of the first i columns to the first j reference positions
    neg = float("-inf")
    score = [[neg] * (SECOND_CYS + 1) for _ in range(n_cols + 1)]
    matched = [[False] * (SECOND_CYS + 1) for _ in range(n_cols + 1)]
    score[0][0] = 0.0
    for i in range(1, n_cols + 1):
        score[i][0] = score[i - 1][0] + delete(i - 1)
        for j in range(1, SECOND_CYS + 1):
            m = score[i - 1][j - 1] + match(i - 1, j - 1)
            d = score[i - 1][j] + delete(i - 1)
            score[i][j], matched[i][j] = (m, True) if m >= d else (d, False)

    end = max(range(n_cols + 1), key=lambda i: score[i][SECOND_CYS])
    mapping = {col: SECOND_CYS + 1 + k for k, col in enumerate(range(end, n_cols))}
    deleted = []
    i, j = end, SECOND_CYS
    while i > 0:
        if j > 0 and matched[i][j]:
            mapping[i - 1] = j
            j -= 1
        else:
            deleted.append(i - 1)
        i -= 1
    return mapping, deleted


def build_v_genes(
    alleles: list[Allele],
) -> tuple[dict[tuple[str, str], list[VGene]], list[str]]:
    """IMGT-numbered V genes per (species, chain), one allele per gene."""
    groups: dict[tuple[str, str], list[Allele]] = defaultdict(list)
    for allele in alleles:
        groups[(allele.species, allele.name[:4])].append(allele)

    by_chain: dict[tuple[str, str], list[VGene]] = {}
    log = []
    for (species, group), members in sorted(groups.items()):
        reference = [a.sequence for a in groups[(REFERENCE_SPECIES, group)]]
        mapping, deleted = map_columns([a.sequence for a in members], reference)
        if deleted:
            columns = ", ".join(str(c + 1) for c in sorted(deleted))
            log.append(f"{species} {group}: removed IMGT/GENE-DB columns {columns}")

        chosen: dict[str, VGene] = {}
        rejected = 0
        for allele in sorted(members, key=lambda a: (a.gene, a.allele_number)):
            if allele.gene in chosen:
                continue
            if any(
                allele.sequence[c] != GAP for c in deleted if c < len(allele.sequence)
            ):
                rejected += 1  # residue in an insertion column: no IMGT position for it
                continue
            numbered = dict.fromkeys(range(1, SECOND_CYS + 1), GAP)
            cdr3 = []
            for col, residue in enumerate(allele.sequence):
                if col in deleted or residue == GAP:
                    continue
                if mapping[col] <= SECOND_CYS:
                    numbered[mapping[col]] = residue
                else:
                    cdr3.append(residue)
            if numbered[FIRST_CYS] != "C" or numbered[SECOND_CYS] != "C":
                rejected += 1
                continue
            chosen[allele.gene] = VGene(
                allele.name, "".join(numbered.values()), "".join(cdr3)
            )
        if rejected:
            log.append(
                f"{species} {group}: skipped {rejected} alleles without IMGT positions"
            )
        # TRAV/DV genes are named and gapped as TRAV and count towards TRA only
        by_chain[(species, group[:3])] = list(chosen.values())
    return by_chain, log


def build_j_genes(
    alleles: list[Allele],
) -> tuple[dict[tuple[str, str], list[JGene]], list[str]]:
    """J genes split at J-TRP/J-PHE per (species, chain), one allele per gene."""
    by_chain: dict[tuple[str, str], dict[str, JGene]] = defaultdict(dict)
    log = []
    for allele in sorted(alleles, key=lambda a: (a.gene, a.allele_number)):
        genes = by_chain[(allele.species, allele.name[:3])]
        if allele.gene in genes:
            continue
        fr4_length = LAST_POSITION - J_ANCHOR_POSITION + 1
        anchors = [m.start() for m in J_MOTIF.finditer(allele.sequence)] or [
            m.start()
            for m in J_MOTIF_FALLBACK.finditer(allele.sequence)
            if fr4_length - 1 <= len(allele.sequence) - m.start() <= fr4_length + 1
        ]
        if not anchors:
            log.append(f"{allele.species} {allele.name}: no J-TRP/J-PHE motif, skipped")
            continue
        anchor = anchors[-1]
        seq = allele.sequence
        genes[allele.gene] = JGene(
            allele.name, seq[:anchor], seq[anchor : anchor + fr4_length]
        )
    return {key: list(genes.values()) for key, genes in by_chain.items()}, log


def number_combination(v: VGene, j: JGene) -> str:
    """IMGT-gapped 128-position sequence of a V-J germline.

    A CDR3 over 13 residues fills 105-111 and 112-117 from either end; the residues in between
    are IMGT insertions (111.1, 112.1, ...), which have no column here.
    """
    cdr3 = v.cdr3 + j.cdr3
    n_positions = CDR3_LAST - CDR3_FIRST + 1
    if len(cdr3) > n_positions:
        left = CDR3_DELETION_ORDER[0] - CDR3_FIRST + 1
        cdr3 = cdr3[:left] + cdr3[left - n_positions :]
    empty = set(CDR3_DELETION_ORDER[: n_positions - len(cdr3)])
    cdr3_positions = [p for p in range(CDR3_FIRST, CDR3_LAST + 1) if p not in empty]
    numbered = dict(zip(cdr3_positions, cdr3))
    cdr3_gapped = "".join(
        numbered.get(p, GAP) for p in range(CDR3_FIRST, CDR3_LAST + 1)
    )
    fr4 = j.framework.ljust(LAST_POSITION - J_ANCHOR_POSITION + 1, GAP)
    return v.framework + cdr3_gapped + fr4


def region_of(position: int) -> str:
    for upper, name in REGION_BOUNDS:
        if position <= upper:
            return name
    return "FR4"


def write_consensus(chain: str, sequences: list[str], output: Path) -> None:
    holes = SCHEME_HOLES.get(chain, ())
    lines = ["position,consensus_aas,frequencies,occupancy,region"]
    for position in range(1, LAST_POSITION + 1):
        column = Counter(seq[position - 1] for seq in sequences)
        occupancy = 1 - column[GAP] / len(sequences)
        if position in holes:
            assert occupancy == 0, f"{chain} has residues at scheme hole {position}"
            continue
        frequent = [
            ("-" if aa == GAP else aa, n / len(sequences))
            for aa, n in sorted(column.items(), key=lambda item: (-item[1], item[0]))
            if n / len(sequences) >= MIN_FREQUENCY
        ]
        lines.append(
            f"{position},{'|'.join(aa for aa, _ in frequent)},"
            f"{'|'.join(f'{f:.3f}' for _, f in frequent)},{occupancy:.3f},{region_of(position)}"
        )
    output.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--cache-dir", type=Path, default=Path(".cache/imgt"))
    parser.add_argument("--output-dir", type=Path, default=Path("resources/consensus"))
    args = parser.parse_args()

    fasta, release = download_imgt(args.cache_dir)
    print(f"IMGT/GENE-DB release {release}")
    alleles = {
        region: [
            a
            for a in read_alleles(fasta, region, set(IG_SPECIES) | set(TR_SPECIES))
            if a.species in species_for(a.name[:3])
        ]
        for region in ("V-REGION", "J-REGION")
    }
    v_genes, v_log = build_v_genes(alleles["V-REGION"])
    j_genes, j_log = build_j_genes(alleles["J-REGION"])
    for line in v_log + j_log:
        print(f"  {line}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    germlines = [f"# IMGT/GENE-DB release {release}", "chain\tspecies\tsegment\tallele"]
    for chain in CHAINS:
        holes = SCHEME_HOLES.get(chain, ())
        sequences = []
        for sp in species_for(chain):
            vs, js = v_genes.get((sp, chain), []), j_genes.get((sp, chain), [])
            for v in vs:
                if any(v.framework[p - 1] != GAP for p in holes):
                    print(
                        f"  {sp} {v.allele}: residue at a {chain} scheme hole, skipped"
                    )
            vs = [v for v in vs if all(v.framework[p - 1] == GAP for p in holes)]
            if not vs or not js:
                continue
            germlines += [f"{chain}\t{sp}\tV\t{v.allele}" for v in vs]
            germlines += [f"{chain}\t{sp}\tJ\t{j.allele}" for j in js]
            sequences += [number_combination(v, j) for v, j in product(vs, js)]
            print(f"{chain} {sp}: {len(vs)} V x {len(js)} J")
        write_consensus(chain, sequences, args.output_dir / f"{chain}.csv")
        print(
            f"{chain}: {len(sequences)} germline combinations -> {args.output_dir / chain}.csv"
        )
    (args.output_dir / "germlines.tsv").write_text("\n".join(germlines) + "\n")


if __name__ == "__main__":
    main()
