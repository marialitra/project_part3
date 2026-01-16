#!/usr/bin/env python3
"""Filter BLAST outfmt 6 results by E-value and per-query top-N bit scores."""

import libraries
from libraries import Path

# -----------------------------
# Main
# -----------------------------
def main() -> None:
    # -----------------------------
    # Parse command-line arguments
    # -----------------------------
    args = libraries.parse_args_blast()
    input_path = Path(args.input)
    
    output_path = Path(args.output)
    if args.top <= 0:
        raise ValueError("--top must be a positive integer")

    # --------------------------------
    # Filter BLAST hits by E-value
    # --------------------------------
    hits_by_query = libraries.filter_hits(input_path, args.max_evalue)

    
    # -----------------------------------------------------
    # Write top N unique hits per query to output
    # -----------------------------------------------------
    libraries.write_top_hits(hits_by_query, output_path, args.top)

if __name__ == "__main__":
    main()