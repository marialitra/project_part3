from libraries import argparse, List, Tuple, Dict, defaultdict, re, Optional


def parse_args_embed() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate protein embeddings using ESM2")
    parser.add_argument('-i', '--input', type=str, required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output vectors.dat file')
    parser.add_argument('-model', '--model', type=str, default='esm2_t6_8M_UR50D')
    parser.add_argument('--batch-size', type=int, default=64)
    parser.add_argument('--max-len', type=int, default=1022)

    return parser.parse_args()

def parse_args_blast() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=("Filter BLAST tabular (outfmt 6) results: keep rows with E <= max-evalue,"
            "sort per query by descending bit score, and retain top N per query."))
    parser.add_argument("-i", "--input", required=True, help="Path to BLAST outfmt 6 TSV file (no header).")
    parser.add_argument("-o", "--output", required=True, help="Path to write filtered TSV results.")
    parser.add_argument("-n", "--top", type=int, required=True, help="Keep only the top N hits per query (sorted by bit score).")
    parser.add_argument("--max-evalue", type=float, default=1e-3, help="Maximum E-value to keep (default: 1e-3).")

    return parser.parse_args()

def parse_args_search() -> argparse.Namespace:
    parser = argparse.ArgumentParser("Protein ANN search wrapper")
    
    parser.add_argument("-d", type=str, required=True, help="Path to base protein vectors .dat (float32 N x 320)")
    parser.add_argument("-q", type=str, required=True, help="Path to FASTA with query sequences")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output neighbors file (text)")
    parser.add_argument("-method", type=str, required=True, choices=["lsh", "hypercube", "ivfflat", "ivfpq", "neural", "all"],
						help="ANN algorithm to use (or 'all' to run all methods)")
	
	# Global Parameters
    parser.add_argument("-N", type=int, default=10, help="Number of nearest neighbors")
    parser.add_argument("-R", type=float, default=0.5, help="Range search radius (for range search mode)")
    parser.add_argument("-range", type=bool, default=False, help="Flag to enable Range Search")
    parser.add_argument("-seed", type=int, default=42, help="Random seed")
	
	# LSH params
    parser.add_argument("-k", type=int, default=2, help="Hash functions per table (LSH)")
    parser.add_argument("-L", type=int, default=5, help="Number of hash tables (LSH)")
    parser.add_argument("--lsh-w", type=float, default=20.0, help="Window width (LSH)")
	
	# Hypercube params
    parser.add_argument("-kproj", type=int, default=12, help="Number of projections (Hypercube)")
    parser.add_argument("--hyper-w", type=float, default=20.0, help="Window width (Hypercube)")
    parser.add_argument("--hyper-M", type=int, default=10000, help="Max candidates to check (Hypercube)")
    parser.add_argument("-probes", type=int, default=100, help="Vertices to examine (Hypercube)")
	
	# IVFFlat params
    parser.add_argument("--flat-kclusters", type=int, default=200, help="Number of clusters (IVF-Flat)")
    parser.add_argument("--flat-nprobe", type=int, default=100, help="Number of probes (IVF-Flat)")
	
	# IVFPQ params
    parser.add_argument("--pq-kclusters", type=int, default=500, help="Number of clusters (IVFPQ)")
    parser.add_argument("--pq-nprobe", type=int, default=100, help="Number of probes (IVFPQ)")
    parser.add_argument("--pq-M", type=int, default=16, help="Number of subvectors (IVFPQ, must divide dimension)")
    parser.add_argument("-nbits", type=int, default=8, help="Bits per subspace (IVFPQ)")

	# NLSH specific - search phase
    parser.add_argument("--nlsh-index", type=str, default=None, help="Directory for NLSH index (default: alongside output)")
    parser.add_argument("--nlsh-T", type=int, default=1000, help="Number of bins to probe (NLSH)")

	# NLSH specific - build phase
    parser.add_argument("--nlsh-m", type=int, default=1800, help="Number of parts for KaHIP (NLSH build)")
    parser.add_argument("--nlsh-imbalance", type=float, default=0.1, help="KaHIP imbalance (NLSH build)")
    parser.add_argument("--nlsh-kahip-mode", type=int, default=0, help="KaHIP mode (NLSH build)")
    parser.add_argument("--nlsh-layers", type=int, default=5, help="MLP layers (NLSH build)")
    parser.add_argument("--nlsh-nodes", type=int, default=128, help="MLP hidden units (NLSH build)")
    parser.add_argument("--nlsh-epochs", type=int, default=8, help="Training epochs (NLSH build)")
    parser.add_argument("--nlsh-batch-size", type=int, default=512, help="Batch size (NLSH build)")
    parser.add_argument("--nlsh-lr", type=float, default=1e-3, help="Learning rate (NLSH build)")

    return parser.parse_args()


def parse_neighbor_results(results_txt: str, topN: int) -> Dict[str, List[Tuple[str, float]]]:
	results = defaultdict(list)
	current_query = None
	
	# Match both formats: "Nearest neighbor-N: ID, 0.123" and "Nearest neighbor-N: ID, Distance: 0.123"
	neighbor_re = re.compile(r"Nearest neighbor-\d+\s*:\s*([^,\s]+)\s*,\s*(?:Distance:\s*)?([\d.]+)")

	with open(results_txt, "r") as f:
		for line in f:
			line = line.strip()
			
			if line.startswith("Query:"):
				current_query = line.split("Query:")[1].strip()
			elif line.startswith("Nearest neighbor") and current_query:
				match = neighbor_re.search(line)
				
				if match and len(results[current_query]) < topN:
					neighbor_id = match.group(1).strip()
					distance = float(match.group(2))
					results[current_query].append((neighbor_id, distance))

	return results

def parse_blast_tsv(blast_tsv_path, topN):
    gt = defaultdict(list)

    with open(blast_tsv_path, "r") as f:
        for line in f:
            if not line.strip():
                continue

            cols = line.strip().split("\t")
            query_id = cols[0]

            # Extract protein ID from sp|ID|...
            target_field = cols[1]
            if target_field.startswith("sp|"):
                target_id = target_field.split("|")[1]
            else:
                continue

            if len(gt[query_id]) < topN:
                gt[query_id].append(target_id)

    # Convert lists to sets
    return {q: set(v) for q, v in gt.items()}


def _extract_neighbor_id(line: str) -> Optional[str]:
	"""
	    Return the neighbor token from a result line, ignoring distance.
	"""
	
	# Match both formats
	match = re.search(r"Nearest neighbor-\d+\s*:\s*([^,\s]+)", line)
	if match:
		return match.group(1).strip()
	
	return None

def parse_ann_txt(ann_txt_path, topN):
	results = defaultdict(list)
	current_query = None

	with open(ann_txt_path, "r") as f:
		for raw_line in f:
			line = raw_line.strip()

			if line.startswith("Query:"):
				current_query = line.split("Query:")[1].strip()
			elif line.startswith("Nearest neighbor") and current_query:
				neighbor_id = _extract_neighbor_id(line)
				
				if neighbor_id and len(results[current_query]) < topN:
					results[current_query].append(neighbor_id)

	return results