"""Generate all V+J protein combinations per chain type from nucleotide FASTA files.

Usage:
    uv run python scripts/generate_vj_combinations.py v_genes.fasta j_genes.fasta -o output_dir

Input:
  V and J gene nucleotide FASTA files with IMGT-style headers: gene_name|functionality[|...|frame=N]
  Only functional (F) genes are included by default; use --include-pseudo for P/ORF genes.

Output:
  Per-chain FASTA files ({chain}_VJ_combinations.fasta) with IMGT-numbered sequences
  (128 positions, '_' for gaps), one file per chain type found in both V and J inputs.
"""

import argparse
import re
from itertools import product
from pathlib import Path
from typing import NamedTuple

from antpack import SingleChainAnnotator

from .utils import parse_fasta, translate, IMGT_POSITIONS

CHAIN_PATTERN = re.compile(r"^(IG[HKL]|TR[ABDG])")


class Gene(NamedTuple):
    name: str
    protein: str


class ChainInfo(NamedTuple):
    letter: str
    allowed: set[str]


CHAIN_MAP: dict[str, ChainInfo] = {
    "IGH": ChainInfo("H", {"H", "K", "L"}),
    "IGK": ChainInfo("K", {"H", "K", "L"}),
    "IGL": ChainInfo("L", {"H", "K", "L"}),
    "TRA": ChainInfo("A", {"A", "B", "D", "G"}),
    "TRB": ChainInfo("B", {"A", "B", "D", "G"}),
    "TRD": ChainInfo("D", {"A", "B", "D", "G"}),
    "TRG": ChainInfo("G", {"A", "B", "D", "G"}),
}


def get_chain(gene_name: str) -> str | None:
    """Extract chain type (e.g. IGH, TRA) from an IMGT gene name."""
    m = CHAIN_PATTERN.match(gene_name)
    return m.group(1) if m else None


def translate_v(seq: str) -> str:
    """Translate V gene from start, trimming incomplete codon at end."""
    remainder = len(seq) % 3
    if remainder:
        seq = seq[: len(seq) - remainder]
    return translate(seq)


def translate_j(seq: str, frame: int | None = None) -> str:
    """Translate J gene based on the frame or from the end, trimming incomplete codons"""
    if frame is not None:
        seq = seq[frame:]
    else:
        remainder = len(seq) % 3
        if remainder:
            seq = seq[remainder:]
    return translate(seq)


def parse_genes(
    fasta_path: Path, translate_fn, allowed_functionality: set[str] | None
) -> dict[str, list[Gene]]:
    """Parse and translate V or J genes from FASTA, grouped by chain type."""
    by_chain: dict[str, list[Gene]] = {}
    for header, seq in parse_fasta(fasta_path):
        parts = header.split("|")
        gene_name = parts[0]
        functionality = parts[1] if len(parts) > 1 else ""
        if allowed_functionality and functionality not in allowed_functionality:
            continue
        chain = get_chain(gene_name)
        if chain is None:
            print(f"  Skipping gene with unknown chain: {gene_name}")
            continue

        if translate_fn == translate_j:
            frame = int(parts[3].split("=")[1]) if len(parts) > 3 else None
            protein = translate_j(seq.upper(), frame)
        else:
            protein = translate_v(seq.upper())

        by_chain.setdefault(chain, []).append(Gene(gene_name, protein))
    return by_chain


def number_sequences(
    sequences: list[Gene],
    chain_prefix: str,
    min_identity: float = 0.7,
) -> list[Gene]:
    """Number sequences with AntPack and return gapped 128-position strings.

    Sequences with wrong chain type or percent identity below min_identity are filtered out.
    """
    chain_letter = chain_prefix[2]
    annotator = SingleChainAnnotator(chains=[chain_letter], scheme="imgt")

    gapped: list[Gene] = []
    failed = wrong_chain = low_identity = 0

    for gene in sequences:
        numbering, pct_identity, chain_type, error = annotator.analyze_seq(gene.protein)

        if error:
            failed += 1
            continue
        if chain_type != chain_letter:
            wrong_chain += 1
            continue
        if pct_identity < min_identity:
            low_identity += 1
            continue

        pos_dict: dict[str, str] = {}
        for i, pos in enumerate(numbering):
            if pos != "-":
                pos_dict[pos] = gene.protein[i]

        gapped_seq = "".join(pos_dict.get(pos, "_") for pos in IMGT_POSITIONS)
        gapped.append(Gene(gene.name, gapped_seq))

    if failed:
        print(f"  Failed to number {failed} sequences")
    if wrong_chain:
        print(f"  Filtered {wrong_chain} sequences with wrong chain type")
    if low_identity:
        print(f"  Filtered {low_identity} sequences with identity < {min_identity}")
    return gapped


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("v_fasta", type=Path, help="V region FASTA file")
    parser.add_argument("j_fasta", type=Path, help="J region FASTA file")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Output directory (default: current directory)",
    )
    parser.add_argument(
        "--include-pseudo",
        action="store_true",
        help="Include pseudogenes (P) and ORFs (default: only functional F genes)",
    )
    parser.add_argument(
        "--min_identity",
        type=float,
        default=0.7,
        help="Minimum percent identity to keep a sequence (default: 0.7)",
    )
    args = parser.parse_args()

    allowed_func = None if args.include_pseudo else {"F"}

    v_by_chain = parse_genes(args.v_fasta, translate_v, allowed_func)
    j_by_chain = parse_genes(args.j_fasta, translate_j, allowed_func)

    common_chains = sorted(set(v_by_chain) & set(j_by_chain))
    if not common_chains:
        print("No matching chains found between V and J files.")
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for chain in common_chains:
        v_genes, j_genes = v_by_chain[chain], j_by_chain[chain]
        combos = [
            Gene(f"{v.name}_{j.name}", v.protein + j.protein)
            for v, j in product(v_genes, j_genes)
        ]
        print(
            f"{chain}: {len(v_genes)} V x {len(j_genes)} J = {len(combos)} combinations"
        )

        print("  Numbering with ANARCI...")
        gapped = number_sequences(combos, chain, min_identity=args.min_identity)
        print(f"  {len(gapped)} sequences numbered successfully")

        out_path = args.output_dir / f"{chain}_VJ_combinations.fasta"
        with open(out_path, "w") as f:
            for gene in gapped:
                f.write(f">{gene.name}\n{gene.protein}\n")
        print(f"  -> {out_path}")


if __name__ == "__main__":
    main()
