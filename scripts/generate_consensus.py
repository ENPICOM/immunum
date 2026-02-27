"""Generate consensus CSV files from numbered sequences.

Usage:
    uv run python scripts/generate_consensus.py csv input.csv output.csv
    uv run python scripts/generate_consensus.py fasta input1.fasta [input2.fasta ...] output.csv

Input formats:
  csv:   Numbered CSV with columns: 1, 2, ..., 128
         Each position column contains a single amino acid or is empty (gap).
  fasta: IMGT-aligned FASTA where each sequence is exactly 128 characters,
         one per IMGT position. Non-amino-acid characters are treated as gaps.

Output columns: position, consensus_aas, frequencies, occupancy, region
List values in consensus_aas and frequencies are pipe-separated (e.g. "A|G|S").
"""

import argparse
from pathlib import Path

import pandas as pd

from .utils import IMGT_POSITIONS

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# IMGT region boundaries
REGION_BOUNDS = [
    (26, "FR1"),
    (38, "CDR1"),
    (55, "FR2"),
    (65, "CDR2"),
    (104, "FR3"),
    (117, "CDR3"),
]


def position_to_region(pos: int) -> str:
    """Map IMGT position number to region name."""
    for upper, name in REGION_BOUNDS:
        if pos <= upper:
            return name
    return "FR4"


def build_consensus(
    df: pd.DataFrame, output_path: Path, min_freq: float = 0.05
) -> None:
    """Build consensus CSV from a DataFrame with IMGT position columns."""
    pos_cols = [c for c in df.columns if c.isdigit()]
    assert len(pos_cols) == 128, f"Expected 128 position columns, found {len(pos_cols)}"

    total = len(df)
    rows = []
    for pos in pos_cols:
        col = df[pos].replace("", "-")
        freqs = col.value_counts() / total
        freqs = freqs[freqs >= min_freq].sort_values(ascending=False)

        if freqs.empty:
            continue

        occupancy = (df[pos] != "").sum() / total
        region = position_to_region(int("".join(c for c in pos if c.isdigit())))

        rows.append(
            {
                "position": pos,
                "consensus_aas": "|".join(freqs.index),
                "frequencies": "|".join(f"{f:.3f}" for f in freqs),
                "occupancy": f"{occupancy:.3f}",
                "region": region,
            }
        )

    pd.DataFrame(rows).to_csv(output_path, index=False)
    print(f"  Saved {output_path}")


def load_aligned_fasta(filepath: Path) -> pd.DataFrame:
    """Load IMGT-aligned FASTA (128 chars/seq) into a DataFrame with position columns."""
    headers, sequences = [], []
    with open(filepath) as f:
        header = None
        seq_lines: list[str] = []
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    sequences.append("".join(seq_lines))
                header = line[1:]
                headers.append(header)
                seq_lines = []
            elif header is not None:
                seq_lines.append(line)
        if header is not None:
            sequences.append("".join(seq_lines))

    rows = []
    for hdr, seq in zip(headers, sequences):
        if len(seq) != 128:
            print(f"  Warning: skipping {hdr} (length {len(seq)}, expected 128)")
            continue
        row = {
            "header": hdr,
            "sequence": "".join(c for c in seq if c.upper() in AMINO_ACIDS),
        }
        for i, pos in enumerate(IMGT_POSITIONS):
            c = seq[i].upper()
            row[pos] = c if c in AMINO_ACIDS else ""
        rows.append(row)

    return pd.DataFrame(rows).fillna("")


def generate_from_csv(input_path: Path, output_path: Path) -> None:
    """Generate consensus from a numbered CSV file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(input_path, dtype=str).fillna("")
    print(f"  {len(df)} sequences from {input_path}")
    build_consensus(df, output_path)


def generate_from_fasta(input_paths: list[Path], output_path: Path) -> None:
    """Generate consensus from one or more IMGT-aligned FASTA files."""
    dfs = []
    for path in input_paths:
        print(f"  Loading {path}")
        df = load_aligned_fasta(path)
        print(f"  {len(df)} sequences")
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True)
    print(f"  Total: {len(df)} sequences")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    build_consensus(df, output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate consensus CSV files.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    csv_parser = subparsers.add_parser("csv", help="From a numbered CSV file")
    csv_parser.add_argument("input", type=Path, help="Input numbered CSV file")
    csv_parser.add_argument("output", type=Path, help="Output consensus CSV file")

    fasta_parser = subparsers.add_parser(
        "fasta", help="From IMGT-aligned FASTA file(s)"
    )
    fasta_parser.add_argument("input", type=Path, nargs="+", help="Input FASTA file(s)")
    fasta_parser.add_argument("output", type=Path, help="Output consensus CSV file")

    args = parser.parse_args()

    if args.command == "csv":
        generate_from_csv(args.input, args.output)
    elif args.command == "fasta":
        generate_from_fasta(args.input, args.output)


if __name__ == "__main__":
    main()
