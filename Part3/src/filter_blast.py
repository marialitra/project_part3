#!/usr/bin/env python3
"""Filter BLAST outfmt 6 results by E-value and per-query top-N bit scores."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

Hit = Tuple[float, float, List[str]]  # (bitscore, evalue, row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter BLAST tabular (outfmt 6) results: keep rows with E <= max-evalue, "
            "sort per query by descending bit score, and retain top N per query."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to BLAST outfmt 6 TSV file (no header).",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path to write filtered TSV results.",
    )
    parser.add_argument(
        "-n",
        "--top",
        type=int,
        required=True,
        help="Keep only the top N hits per query (sorted by bit score).",
    )
    parser.add_argument(
        "--max-evalue",
        type=float,
        default=1e-3,
        help="Maximum E-value to keep (default: 1e-3).",
    )
    return parser.parse_args()


def filter_hits(input_path: Path, max_evalue: float) -> Dict[str, List[Hit]]:
    hits_by_query: Dict[str, List[Hit]] = defaultdict(list)
    with input_path.open(newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 12:
                # Skip malformed rows.
                continue
            try:
                evalue = float(row[10])
                bitscore = float(row[11])
            except ValueError:
                # Skip rows with non-numeric evalue/bitscore.
                continue
            if evalue > max_evalue:
                continue
            query_id = row[0]
            hits_by_query[query_id].append((bitscore, evalue, row))
    return hits_by_query


def write_top_hits(hits_by_query: Dict[str, List[Hit]], output_path: Path, top_n: int) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        for query_id, hits in hits_by_query.items():
            # Sort by descending bitscore, then ascending evalue, then subject id for determinism.
            hits.sort(key=lambda h: (-h[0], h[1], h[2][1] if len(h[2]) > 1 else ""))
            for _, _, row in hits[:top_n]:
                writer.writerow(row)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    if args.top <= 0:
        raise ValueError("--top must be a positive integer")

    hits_by_query = filter_hits(input_path, args.max_evalue)
    write_top_hits(hits_by_query, output_path, args.top)


if __name__ == "__main__":
    main()
