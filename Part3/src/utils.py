import libraries
from libraries import Path, Tuple, List, Dict, defaultdict, csv, esm, torch, SeqIO, Any, np, os, re

# A BLAST hit is represented as:
# (bitscore, evalue, original row fields)
Hit = Tuple[float, float, List[str]]

def filter_hits(input_path: Path, max_evalue: float) -> Dict[str, List[Hit]]:
    """
        Read a BLAST outfmt 6 file and filter hits by E-value.
        For each query sequence, all hits with an E-value less than or equal
        to `max_evalue` are collected.
    """
        
    hits_by_query: Dict[str, List[Hit]] = defaultdict(list)
    with input_path.open(newline="") as f:
        reader = csv.reader(f, delimiter="\t")

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
        writer = csv.writer(f, delimiter="\t")

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






def load_model(DEVICE: Any) -> Tuple[torch.nn.Module, esm.data.Alphabet, callable]:
    """
        Load a pretrained ESM model and move it to the specified device.
    """
        
    model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    model = model.to(DEVICE)
    model.eval()
    batch_converter = alphabet.get_batch_converter()

    return model, alphabet, batch_converter

def load_sequences(FASTA: str, MAX_LEN: int) -> Tuple[list, int]:
    """
        Load protein sequences from a FASTA file and truncate them to a maximum length.
        Sequences are sorted by descending length to improve batching efficiency.
    """

    data = []
    for record in SeqIO.parse(FASTA, "fasta"):
        # Extract UniProt accession from >sp|A7IC07|RBFA_XANP2
        parts = record.id.split('|')
        protein_id = parts[1] if len(parts) == 3 else record.id

        seq = str(record.seq)[:MAX_LEN]
        data.append((protein_id, seq))

    # Sort by length for efficient batching
    data.sort(key=lambda x: len(x[1]), reverse=True)
    N = len(data)

    return data, N

def save_output(VECTORS_FILE: str, N: int, EMBED_DIM: int, IDS_FILE: str, BATCH_SIZE: int, data: list,
                batch_converter: callable, DEVICE: Any, model: torch.nn.Module, alphabet: esm.data.Alphabet,
                vectors: np.memmap) -> None:
    """
        Generate embeddings for protein sequences and save them to disk.
        Embeddings are written to a memory-mapped NumPy array, and the corresponding
        protein IDs are written to a separate text file.
    """

    with open(IDS_FILE, 'w', encoding='ascii', errors='replace') as id_file_handle:
        with torch.no_grad():
            idx = 0
            for i in range(0, N, BATCH_SIZE):
                batch = data[i:i + BATCH_SIZE]

                _, _, tokens = batch_converter(batch)
                tokens = tokens.to(DEVICE)

                results = model(tokens, repr_layers=[6])
                reps = results["representations"][6]
                for j, tok in enumerate(tokens):
                    valid = (tok != alphabet.padding_idx) & \
                            (tok != alphabet.cls_idx) & \
                            (tok != alphabet.eos_idx)

                    emb = reps[j][valid].mean(dim=0).cpu().numpy()

                    vectors[idx] = emb
                    id_file_handle.write(batch[j][0] + '\n')
                    idx += 1

    vectors.flush()

    print(f"Saved embeddings to: {VECTORS_FILE}")
    print(f"Saved ID mapping to: {IDS_FILE}")
    print(f"Shape: ({N}, {EMBED_DIM})")



# FOR SEARCH CODE:


def remap_output_ids(output_txt: str, base_ids_txt: str, query_ids_txt: str):
	"""
		Remap numeric indices in output to actual protein IDs (distance-safe).
	"""

	if not os.path.exists(output_txt):
		return
	
	# Load ID mappings
	base_ids = []
	if os.path.exists(base_ids_txt):
		with open(base_ids_txt, "r") as f:
			base_ids = [line.strip() for line in f if line.strip()]
	
	query_ids = []
	if os.path.exists(query_ids_txt):
		with open(query_ids_txt, "r") as f:
			query_ids = [line.strip() for line in f if line.strip()]
	
	if not base_ids and not query_ids:
		return
	
	with open(output_txt, "r") as f:
		lines = f.readlines()
	
	remapped_lines = []
	neighbor_colon_re = re.compile(r"^(Nearest neighbor-\d+)\s*:\s*(\d+)(.*)$")

	for raw_line in lines:
		line = raw_line.rstrip("\n")

		if line.startswith("Query:"):
			parts = line.split()
			if len(parts) >= 2:
				try:
					idx = int(parts[1])
					if 0 <= idx < len(query_ids):
						line = f"Query: {query_ids[idx]}"
				except (ValueError, IndexError):
					pass
		elif line.startswith("Nearest neighbor"):
			match = neighbor_colon_re.match(line)
			if match:
				prefix, idx_str, rest = match.groups()
				try:
					idx = int(idx_str)
					if 0 <= idx < len(base_ids):
						line = f"{prefix}: {base_ids[idx]}{rest}"
				except ValueError:
					pass
		remapped_lines.append(line + "\n")
	
	with open(output_txt, "w") as f:
		f.writelines(remapped_lines)

	print(f"[protein_search] Remapped indices to protein IDs in {output_txt}")


def print_recall(name, mean_recall):
    if name == "lsh":
        print(f"  {'Euclidean LSH':15s}: {mean_recall:.4f}")
    elif name == "hypercube":
        print(f"  {'Hypercube':15s}: {mean_recall:.4f}")
    elif name == "ivfflat":
        print(f"  {'IVF-Flat':15s}: {mean_recall:.4f}")
    elif name == "ivfpq":
        print(f"  {'IVF-PQ':15s}: {mean_recall:.4f}")
    elif name == "nlsh":
        print(f"  {'Neural LSH':15s}: {mean_recall:.4f}")
    else:
        print(f"  {name.capitalize():15s}: {mean_recall:.4f}")

def calculate_recall(blast_tsv, ann_txt, topN):
    blast_gt = libraries.parse_blast_tsv(blast_tsv, topN)
    ann_res = libraries.parse_ann_txt(ann_txt, topN)

    recalls = []
    per_query = {}

    for query, gt_neighbors in blast_gt.items():
        ann_neighbors = set(ann_res.get(query, []))

        if not gt_neighbors:
            continue

        intersection = gt_neighbors & ann_neighbors
        recall = len(intersection) / len(gt_neighbors)

        recalls.append(recall)
        per_query[query] = recall

    mean_recall = sum(recalls) / len(recalls) if recalls else 0.0
    return mean_recall, per_query

def compute_recall(args, method, current_methods, blast_results, query_ids, all_qps, blast_identity, blast_time_per_query, blast_qps):
	# Compute Recall against BLAST results (topN.tsv) vs the results.txt
	# If method.lower == "all" then compare all the results_algo.txt vs the BLAST_results.tsv
	if method == "all":
		base_output = os.path.splitext(args.output)[0]

		print("\nRecall results (average):")
		method_recall = {}
		method_results = {}

		for algo in current_methods:
			ann_txt = f"{base_output}_{algo}.txt"
			mean_recall, _ = calculate_recall(blast_tsv=blast_results, ann_txt=ann_txt, topN=args.N)
			method_recall[algo] = mean_recall

			# Parse neighbor results for formatted report
			method_results[algo] = libraries.parse_neighbor_results(ann_txt, args.N)
				
			# Print Recall
			print_recall(algo, mean_recall)

		# Parse BLAST top-N
		blast_gt = libraries.parse_blast_tsv(blast_results, args.N)

		# Generate a single consolidated report for all methods
		# report_path = base_output + "_REPORT.txt"
		libraries.generate_all_methods_report(
			output_report=args.output,
			query_ids=query_ids,
			topN=args.N,
			method_qps=all_qps,
			method_results=method_results,
			methods=current_methods,
			blast_results_topn=blast_gt,
			blast_results_identity=blast_identity,
			blast_time_per_query = blast_time_per_query,
			blast_qps = blast_qps
		)
	else:
		mean_recall, _ = calculate_recall(blast_tsv=blast_results, ann_txt=args.output, topN=args.N)

		# Print Recall
		print("\nRecall result (average):")
		print_recall(method, mean_recall)

		# Parse neighbor results
		method_results = libraries.parse_neighbor_results(args.output, args.N)

		# Parse BLAST top-N
		blast_gt = libraries.parse_blast_tsv(blast_results, args.N)

		# Determine method display name
		method_display_map = {
			"lsh": "Euclidean LSH",
			"hypercube": "Hypercube",
			"ivfflat": "IVF-Flat",
			"ivfpq": "IVF-PQ",
			"nlsh": "Neural LSH",
		}
		method_display = method_display_map.get(method, method.capitalize())

		# Generate per-query report
		libraries.generate_per_query_report(
			output_report=args.output,
			query_ids=query_ids,
			topN=args.N,
			method_name=method_display,
			method_key=method,
			qps=all_qps,
			method_results=method_results,
			blast_results_topn=blast_gt,
			blast_results_identity=blast_identity,
			blast_time_per_query = blast_time_per_query,
			blast_qps = blast_qps
		)


def print_QPS_dictionary(all_qps, method, current_methods) -> list[str]:
	# Print All QPS status if it is a dictionary
	if method == "all":
		print("\nQPS results:")

		for name, qps in all_qps.items():
			if qps == None:
				current_methods.remove(name)
				continue

			if name == "lsh":
				name = "Euclidean LSH"
			elif name == "hypercube":
				name = "Hypercube"
			elif name == "ivfflat":
				name = "IVF-Flat"
			elif name == "ivfpq":
				name = "IVF-PQ"
			elif name == "nlsh":
				name = "Neural LSH"
			else:
				name = name.capitalize() # Everything else
			
			print(f"  {name:15s}: {qps}")

		return current_methods
	else:
		print("Wrong method called, method = all")
		return []

def print_QPS(all_qps, method):
		# Print All QPS status
		if method == "lsh":
			method = "Euclidean LSH"
		elif method == "hypercube":
			method = "Hypercube"
		elif method == "ivfflat":
			method = "IVF-Flat"
		elif method == "ivfpq":
			method = "IVF-PQ"
		elif method == "nlsh":
			method = "Neural LSH"
		else:
			method = method.capitalize() # Everything else

		print("\nQPS result:")
		print(f"  {method:10s}: {all_qps}")


def delete_files(args, method, all_methods):
	# Deleting unecessary files
	if method == "all":
		base_output = os.path.splitext(args.output)[0]
		
		for algo in all_methods:
			if os.path.exists(f"{base_output}_{algo}.txt"):
				os.remove(f"{base_output}_{algo}.txt")

			if os.path.exists(f"{base_output}_{algo}.queries_ids.txt"):
				os.remove(f"{base_output}_{algo}.queries_ids.txt")

			if os.path.exists(f"{base_output}_{algo}.queries.dat"):
				os.remove(f"{base_output}_{algo}.queries.dat")
	else:
		base_output = os.path.splitext(args.output)[0]

		if os.path.exists(f"{base_output}.queries.dat"):
			os.remove(f"{base_output}.queries.dat")
		
		if os.path.exists(f"{base_output}.queries_ids.txt"):
			os.remove(f"{base_output}.queries_ids.txt")