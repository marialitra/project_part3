#!/usr/bin/env python3
"""Filter BLAST outfmt 6 results by E-value and per-query top-N bit scores."""

import libraries
from libraries import Dict, List, Tuple, defaultdict, Path

# A BLAST hit is represented as:
# (bitscore, evalue, original row fields)
Hit = Tuple[float, float, List[str]]

def parse_args() -> libraries.argparse.Namespace:
    """
        Parse and validate command-line arguments.
    """

    parser = libraries.argparse.ArgumentParser(
        description=(
            "Filter BLAST tabular (outfmt 6) results: keep rows with E <= max-evalue, "
            "sort per query by descending bit score, and retain top N per query."
        )
    )
    parser.add_argument("-i", "--input", required=True, help="Path to BLAST outfmt 6 TSV file (no header).")
    parser.add_argument("-o", "--output", required=True, help="Path to write filtered TSV results.")
    parser.add_argument("-n", "--top", type=int, required=True, help="Keep only the top N hits per query (sorted by bit score).")
    parser.add_argument("--max-evalue", type=float, default=1e-3, help="Maximum E-value to keep (default: 1e-3).")
    return parser.parse_args()


def filter_hits(input_path: Path, max_evalue: float) -> Dict[str, List[Hit]]:
    """
        Read a BLAST outfmt 6 file and filter hits by E-value.
        For each query sequence, all hits with an E-value less than or equal
        to `max_evalue` are collected.
    """
        
    hits_by_query: Dict[str, List[Hit]] = defaultdict(list)
    with input_path.open(newline="") as f:
        reader = libraries.csv.reader(f, delimiter="\t")

        for row in reader:
            if not row or row[0].startswith("#"):
                continue

            # BLAST outfmt 6 should have at least 12 columns
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
    
    """
        Write the top N unique protein hits per query to an output file.
        Duplicate protein IDs are removed.
    """
    
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="") as f:
        writer = libraries.csv.writer(f, delimiter="\t")

        for query_id, hits in hits_by_query.items():
            # Sort by descending bitscore, then ascending evalue, then subject id for determinism.
            hits.sort(key=lambda h: (-h[0], h[1], h[2][1] if len(h[2]) > 1 else ""))
            
            # Track seen proteins to avoid duplicates
            seen_proteins = set()
            unique_count = 0

            for _, _, row in hits:
                if unique_count >= top_n:
                    break
                
                # Extract protein ID from subject field (sp|ID|...)
                subject_field = row[1] if len(row) > 1 else ""
                
                if subject_field.startswith("sp|"):
                    protein_id = subject_field.split("|")[1]
                else:
                    protein_id = subject_field
                
                # Skip if we've already seen this protein
                if protein_id in seen_proteins:
                    continue
                
                seen_proteins.add(protein_id)
                writer.writerow(row)
                unique_count += 1


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    # -----------------------------
    # Parse command-line arguments
    # -----------------------------
    args = parse_args()
    input_path = Path(args.input)
    
    output_path = Path(args.output)
    if args.top <= 0:
        raise ValueError("--top must be a positive integer")

    # --------------------------------
    # Filter BLAST hits by E-value
    # --------------------------------
    hits_by_query = filter_hits(input_path, args.max_evalue)

    
    # -----------------------------------------------------
    # Writee top N unique hits per query to output
    # -----------------------------------------------------
    write_top_hits(hits_by_query, output_path, args.top)


if __name__ == "__main__":
    main()
