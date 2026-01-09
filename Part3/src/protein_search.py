import libraries
from libraries import os, subprocess, Optional, Dict, List, Tuple, defaultdict

# Local imports
libraries.sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "Algorithms", "src"))
from runSearchExe import build_executable, run_algorithm


def blast_executable(neighbors: int, output: str):
    """
        Runs 'make blast' to run the BLAST method and its filtering.
    """

    try:
        start = libraries.time.perf_counter()
        build_process = subprocess.run(["make", "blast", f"N={neighbors}", f"out={output}"], capture_output=True, text=True, check=True)
        end = libraries.time.perf_counter()

        total_time = end - start
        
        print("\n\nBuild complete: BLAST results are ready.")
        return total_time
    except subprocess.CalledProcessError as e:
        print("\n\n--- ERROR: Build failed. ---")
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        return False
    except FileNotFoundError:
         print("\n\n--- ERROR: 'make' command not found. Is it installed? ---")
         return False


def parse_blast_tsv(blast_tsv_path, topN):
    """
    Returns:
        dict: {query_id: set(topN target_ids)}
    """
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
	"""Return the neighbor token from a result line, ignoring distance."""
	# Match both formats
	match = libraries.re.search(r"Nearest neighbor-\d+\s*:\s*([^,\s]+)", line)
	if match:
		return match.group(1).strip()
	return None


def parse_ann_txt(ann_txt_path, topN):
	"""
	Returns:
		dict: {query_id: list(topN neighbor_ids)}
	"""
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

def compute_recall(blast_tsv, ann_txt, topN):
    blast_gt = parse_blast_tsv(blast_tsv, topN)
    ann_res = parse_ann_txt(ann_txt, topN)

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


def parse_blast_results_with_identity(blast_tsv_path: str) -> Dict[str, List[Tuple[str, float]]]:
	"""
	Parse BLAST results.tsv and return {query_id: [(target_id, identity%), ...]}
	Identity is extracted from column 2 (0-indexed).
	"""
	blast_data = defaultdict(list)

	with open(blast_tsv_path, "r") as f:
		for line in f:
			if not line.strip():
				continue
			cols = line.strip().split("\t")
			if len(cols) < 3:
				continue

			query_id = cols[0]
			target_field = cols[1]
			identity_str = cols[2]

			# Extract protein ID from sp|ID|... format
			if target_field.startswith("sp|"):
				target_id = target_field.split("|")[1]
			else:
				target_id = target_field

			try:
				identity = float(identity_str)
			except ValueError:
				continue

			blast_data[query_id].append((target_id, identity))

	return blast_data


def parse_neighbor_results(results_txt: str, topN: int) -> Dict[str, List[Tuple[str, float]]]:
	"""Parse results file and return {query_id: [(neighbor_id, distance), ...]}"""
	results = defaultdict(list)
	current_query = None
	# Match both formats: "Nearest neighbor-N: ID, 0.123" and "Nearest neighbor-N: ID, Distance: 0.123"
	neighbor_re = libraries.re.compile(r"Nearest neighbor-\d+\s*:\s*([^,\s]+)\s*,\s*(?:Distance:\s*)?([\d.]+)")

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



def generate_per_query_report(
	output_report: str,
	query_ids: list,
	topN: int,
	method_name: str,
	method_key: str,
	qps: Optional[float],
	method_results: Dict[str, List[Tuple[str, float]]],
	blast_results_topn: Dict[str, set],
	blast_results_identity: Dict[str, List[Tuple[str, float]]],
	blast_time_per_query: float,
	blast_qps: float
):
	"""
	Generate comprehensive per-query report with individual recall.

	Args:
		output_report: Path to write report
		query_fasta: Path to query FASTA file
		topN: Number of neighbors in Top-N
		method_name: Display name (e.g., "Euclidean LSH")
		method_key: Internal key (e.g., "lsh")
		qps: Queries per second
		method_results: {query_id: [(neighbor_id, distance), ...]}
		blast_results_topn: {query_id: set(neighbor_ids)}
		blast_results_identity: {query_id: [(neighbor_id, identity%), ...]}
	"""

	with open(output_report, "w") as f:
		# Iterate through each query
		for query_id in query_ids:
			f.write("\n" + "=" * 110 + "\n")
			f.write(f"Query: {query_id}\n")
			f.write(f"N = {topN}\n")
			f.write("=" * 110 + "\n\n")

			# Get neighbors for this query from results
			neighbors = method_results.get(query_id, [])
			blast_top_n = blast_results_topn.get(query_id, set())
			blast_identities = {tid: ident for tid, ident in blast_results_identity.get(query_id, [])}

			# Compute recall for this query
			if blast_top_n:
				neighbors_in_blast = sum(1 for nid, _ in neighbors[:topN] if nid in blast_top_n)
				recall = neighbors_in_blast / len(blast_top_n)
			else:
				recall = 0.0

			# [1] Method comparison (this method vs BLAST)
			f.write("[1] Brief Method Comparison\n")
			f.write("-" * 110 + "\n")
			f.write(f"{'Method':<20} | {'Time/query (s)':<15} | {'QPS':<10} | {'Recall@N vs BLAST Top-N':<25}\n")
			f.write("-" * 110 + "\n")

			if qps is not None and qps > 0:
				time_per_query = 1.0 / qps
				f.write(f"{method_name:<20} | {time_per_query:<15.3f} | {qps:<10.3f} | {recall:<25.2f}\n")
			else:
				f.write(f"{method_name:<20} | {'N/A':<15s} | {'N/A':<10s} | {recall:<25.2f}\n")

			f.write(f"{'BLAST (Ref)':<20} | {blast_time_per_query:<15.3f} | {blast_qps:<10.3f} | {1.00:<25.2f}\n")
			f.write("-" * 110 + "\n\n")

			# [2] Top-N neighbors
			f.write("[2] Top-N neighbors\n")
			f.write("-" * 110 + "\n")
			f.write("\n" + "Method: " + method_name + "\n")
			f.write(f"{'Rank':<6} | {'Neighbor ID':<15} | {'Distance':<12} | {'BLAST Identity':<17} | {'In BLAST Top-N?':<17} | {'Bio Comment':<30}\n")
			f.write("-" * 110 + "\n")

			# Table rows
			for rank, (neighbor_id, distance) in enumerate(neighbors[:topN], 1):
				in_blast = neighbor_id in blast_top_n
				blast_in_str = "Yes" if in_blast else "No"
				blast_id_val = blast_identities.get(neighbor_id)
				blast_id_str = f"{f'{blast_id_val:.2f}%' : <17}" if blast_id_val is not None else "undetected".ljust(17)

				# Bio comment logic
				if in_blast and blast_id_val is not None and blast_id_val > 30:
					bio_comment = "Homolog"
				elif in_blast and blast_id_val is not None and 20 < blast_id_val <= 30:
					bio_comment = "Remote homolog"
				elif not in_blast:
					bio_comment = "Possible false positive"
				else:
					bio_comment = ""

				f.write(
					f"{rank:<6} | {neighbor_id:<15} | {distance:<12.2f} | {blast_id_str} | {blast_in_str:<17} | {bio_comment:<30}\n"
				)

	print(f"[protein_search] Per-query report written to {output_report}")


def generate_all_methods_report(
	output_report: str,
	query_ids: list,
	topN: int,
	method_qps: Dict[str, Optional[float]],
	method_results: Dict[str, Dict[str, List[Tuple[str, float]]]],
	blast_results_topn: Dict[str, set],
	blast_results_identity: Dict[str, List[Tuple[str, float]]],
	blast_time_per_query: float,
	blast_qps:float
):
	"""
	Generate a single consolidated report for method="all" with per-query summaries and neighbors across methods.
	Per-query recall is computed per method. The neighbor tables reuse the current per-method format.
	"""

	method_display = {
		"lsh": "Euclidean LSH",
		"hypercube": "Hypercube",
		"ivfflat": "IVF-Flat",
		"ivfpq": "IVF-PQ",
		"nlsh": "Neural LSH",
	}
	method_order = ["lsh", "hypercube", "ivfflat", "ivfpq", "nlsh"]

	with open(output_report, "w") as f:
		f.write("=" * 110 + "\n")
		f.write("ANN Evaluation Report (method = all)\n")
		f.write(f"Top-N = {topN}\n")
		f.write("=" * 110 + "\n\n")

		for query_id in query_ids:
			f.write("\n" + "=" * 110 + "\n")
			f.write(f"Query: {query_id}\n")
			f.write(f"N = {topN}\n")
			f.write("=" * 110 + "\n\n")

			# [1] Per-query method comparison
			f.write("[1] Brief Method Comparison\n")
			f.write("-" * 110 + "\n")
			f.write(f"{'Method':<20} | {'Time/query (s)':<15} | {'QPS':<10} | {'Recall@N vs BLAST Top-N':<25}\n")
			f.write("-" * 110 + "\n")

			blast_top_n = blast_results_topn.get(query_id, set())

			for m in method_order:
				name = method_display[m]
				qps_val = method_qps.get(m)
				neighbors = method_results.get(m, {}).get(query_id, [])
				if blast_top_n:
					hits = sum(1 for nid, _ in neighbors[:topN] if nid in blast_top_n)
					recall_q = hits / len(blast_top_n)
				else:
					recall_q = 0.0

				if qps_val is not None and qps_val > 0:
					time_per_query = 1.0 / qps_val
					f.write(f"{name:<20} | {time_per_query:<15.3f} | {qps_val:<10.3f} | {recall_q:<25.2f}\n")
				else:
					f.write(f"{name:<20} | {'N/A':<14s} | {'N/A':<10s} | {recall_q:<25.2f}\n")

			# BLAST reference row
			f.write(f"{'BLAST (Ref)':<20} | {blast_time_per_query:<15.3f} | {blast_qps:<10.3f} | {1.00:<25.2f}\n")
			f.write("-" * 110 + "\n\n")

			# [2] Top-N neighbors per method
			f.write("[2] Top-N neighbors per method\n")
			for m in method_order:
				name = method_display[m]
				neighbors = method_results.get(m, {}).get(query_id, [])
				blast_identities = {tid: ident for tid, ident in blast_results_identity.get(query_id, [])}

				f.write("\n" + "Method: " + name + "\n")
				f.write(f"{'Rank':<6} | {'Neighbor ID':<15} | {'Distance':<12} | {'BLAST Identity':<17} | {'In BLAST Top-N?':<17} | {'Bio Comment':<30}\n")
				f.write("-" * 110 + "\n")

				for rank, (neighbor_id, distance) in enumerate(neighbors[:topN], 1):
					in_blast = neighbor_id in blast_top_n
					blast_in_str = "Yes" if in_blast else "No"
					blast_id_val = blast_identities.get(neighbor_id)
					blast_id_str = f"{f'{blast_id_val:.0f}%' : <17}" if blast_id_val is not None else "undetected".ljust(17)

					if in_blast and blast_id_val is not None and blast_id_val > 30:
						bio_comment = "Homolog"
					elif in_blast and blast_id_val is not None and 20 < blast_id_val <= 30:
						bio_comment = "Remote homolog"
					elif not in_blast:
						bio_comment = "Possible false positive"
					else:
						bio_comment = ""

					f.write(
						f"{rank:<6} | {neighbor_id:<15} | {distance:<12.2f} | {blast_id_str} | {blast_in_str:<17} | {bio_comment:<30}\n"
					)


	print(f"[protein_search] Consolidated all-methods report written to {output_report}")




def remap_output_ids(output_txt: str, base_ids_txt: str, query_ids_txt: str):
	"""Remap numeric indices in output to actual protein IDs (distance-safe)."""
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
	neighbor_colon_re = libraries.re.compile(r"^(Nearest neighbor-\d+)\s*:\s*(\d+)(.*)$")

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

def run_nlsh(
	base_dat: str,
	query_dat: str,
	output_txt: str,
	N: int,
	R: float,
	seed: int,
	nlsh_index: Optional[str],
	nlsh_T: int,
	nlsh_m: int,
	nlsh_imbalance: float,
	nlsh_kahip_mode: int,
	nlsh_layers: int,
	nlsh_nodes: int,
	nlsh_epochs: int,
	nlsh_batch_size: int,
	nlsh_lr: float,
):
	"""Build (if missing) and run NLSH on protein .dat data."""
	alg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "Algorithms"))
	build_py = os.path.join(alg_root, "src", "nlsh_build.py")
	search_py = os.path.join(alg_root, "src", "nlsh_search.py")

	index_dir = nlsh_index or os.path.join(os.path.dirname(output_txt), "nlsh_index_protein")
	os.makedirs(index_dir, exist_ok=True)

	model_path = os.path.join(index_dir, "model.pth")
	inverted_path = os.path.join(index_dir, "inverted_file.npy")

	need_build = not (os.path.exists(model_path) and os.path.exists(inverted_path))
	if need_build:
		build_cmd = [
			libraries.sys.executable,
			build_py,
			"-d", os.path.abspath(base_dat),
			"-i", os.path.abspath(index_dir),
			"--type", "protein",
			"--knn", str(N),
			"-m", str(nlsh_m),
			"--imbalance", str(nlsh_imbalance),
			"--kahip_mode", str(nlsh_kahip_mode),
			"--layers", str(nlsh_layers),
			"--nodes", str(nlsh_nodes),
			"--epochs", str(nlsh_epochs),
			"--batch_size", str(nlsh_batch_size),
			"--lr", str(nlsh_lr),
			"--seed", str(seed),
		]

		print(f"[protein_search] Building NLSH index at {index_dir} ...")
		subprocess.run(build_cmd, check=True, cwd=alg_root)
	else:
		print(f"[protein_search] Reusing existing NLSH index at {index_dir} (model + inverted_file found)")

	search_cmd = [
		libraries.sys.executable,
		search_py,
		"-d", os.path.abspath(base_dat),
		"-q", os.path.abspath(query_dat),
		"-i", os.path.abspath(index_dir),
		"-o", os.path.abspath(output_txt),
		"-type", "protein",
		"-N", str(N),
		"-R", str(R),
		"-T", str(nlsh_T),
		"-range", "false",
	]

	print("[protein_search] Running NLSH search on protein data ...")
	subprocess.run(search_cmd, check=True, cwd=alg_root)
	
	# Remap output IDs
	base_ids_txt = base_dat.replace(".dat", "_ids.txt")
	if not os.path.exists(base_ids_txt):
		raise RuntimeError("Failed to find the base ids text")

	query_ids_txt = query_dat.replace(".dat", "_ids.txt")
	if not os.path.exists(query_ids_txt):
		raise RuntimeError("Failed to find the queries ids text")
	
	remap_output_ids(output_txt, base_ids_txt, query_ids_txt)

	qps = None
	output_lines = []

	with open(output_txt, "r") as f:
		for line in f:
			if line.startswith("QPS:"):
				qps = float(line.split(":", 1)[1].strip())
			else:
				output_lines.append(line)

	with open(output_txt, "w") as f:
		f.writelines(output_lines)

	return qps

def run_protein_search(
	base_dat: str,
	query_fasta: str,
	output_txt: str,
	method: str = "ivfflat",
	N: int = 10,
	R: float = 0.5,
	seed: int = 42,
	# IVFFlat / IVFPQ params
	kclusters: int = 200,
	nprobe: int = 80,
	# LSH params
	k: int = 4,
	L: int = 70,
	w: float = 4.0,
	# Hypercube params
	kproj: int = 14,
	M: int = 2000,
	probes: int = 100,
	# IVFPQ specific
	nbits: int = 8,
	ivfpq_M: int = 16,
	# NLSH specific
	nlsh_index: Optional[str] = None,
	nlsh_T: int = 300,
	nlsh_m: int = 2000,
	nlsh_imbalance: float = 0.1,
	nlsh_kahip_mode: int = 0,
	nlsh_layers: int = 3,
	nlsh_nodes: int = 256,
	nlsh_epochs: int = 5,
	nlsh_batch_size: int = 512,
	nlsh_lr: float = 1e-3,
):
	"""
	1) Build query embeddings via Part3/src/protein_embed.py -> query.dat
	2) Build C executable if needed
	3) Run selected ANN algorithm on protein data (cosine), writing results to output_txt
	
	Args:
		method: ANN algorithm to use ('lsh', 'hypercube', 'ivfflat', 'ivfpq')
	"""

	# Handle 'all' method by recursively calling for each algorithm
	if method == "all":
		all_methods = ["lsh", "hypercube", "ivfflat", "ivfpq", "nlsh"]
		all_qps = {}
		base_output = os.path.splitext(output_txt)[0]
		
		for algo in all_methods:
			algo_output = f"{base_output}_{algo}.txt"
			print(f"\n{'='*60}")
			print(f"[protein_search] Running method: {algo.upper()}")
			print(f"{'='*60}")
			
			try:
				if algo == "nlsh":

					# 1. Create query .dat using protein_embed
					embed_py = os.path.join(os.path.dirname(__file__), "protein_embed.py")
					query_dat = os.path.splitext(algo_output)[0] + ".queries.dat"

					os.makedirs(os.path.dirname(algo_output), exist_ok=True)

					embed_cmd = [
						libraries.sys.executable,
						embed_py,
						"-i", query_fasta,
						"-o", query_dat,
					]

					print("[protein_search] Generating query embeddings (.dat)...")
					subprocess.run(embed_cmd, check=True)

					# 2. Build executable if missing (needed for C-based methods and for NLSH graph build)
					alg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "Algorithms"))
					alg_part1 = os.path.join(alg_root, "AlgorithmsPart1")
					exe_path = os.path.join(alg_part1, "search")
					if not os.path.exists(exe_path):
						print("[protein_search] Building C executable ...")
						prev_cwd = os.getcwd()
						try:
							os.chdir(alg_root)
							if not build_executable():
								raise RuntimeError("Failed to build C executable")
						except Exception:
							os.chdir(prev_cwd)
							raise
						os.chdir(prev_cwd)

					qps = run_nlsh(
						base_dat=base_dat,
						query_dat=query_dat,
						output_txt=algo_output,
						N=N,
						R=R,
						seed=seed,
						nlsh_index=nlsh_index,
						nlsh_T=nlsh_T,
						nlsh_m=nlsh_m,
						nlsh_imbalance=nlsh_imbalance,
						nlsh_kahip_mode=nlsh_kahip_mode,
						nlsh_layers=nlsh_layers,
						nlsh_nodes=nlsh_nodes,
						nlsh_epochs=nlsh_epochs,
						nlsh_batch_size=nlsh_batch_size,
						nlsh_lr=nlsh_lr,
					)
					all_qps[algo] = qps

					print(f"[protein_search] Done. Results at: {algo_output}")

				else:
					qps = run_protein_search(
						base_dat=base_dat,
						query_fasta=query_fasta,
						output_txt=algo_output,
						method=algo,
						N=N,
						R=R,
						seed=seed,
						kclusters=kclusters,
						nprobe=nprobe,
						k=k,
						L=L,
						w=w,
						kproj=kproj,
						M=M,
						probes=probes,
						nbits=nbits,
						ivfpq_M=ivfpq_M,
						nlsh_index=nlsh_index,
						nlsh_T=nlsh_T,
						nlsh_m=nlsh_m,
						nlsh_imbalance=nlsh_imbalance,
						nlsh_kahip_mode=nlsh_kahip_mode,
						nlsh_layers=nlsh_layers,
						nlsh_nodes=nlsh_nodes,
						nlsh_epochs=nlsh_epochs,
						nlsh_batch_size=nlsh_batch_size,
						nlsh_lr=nlsh_lr,
					)
					all_qps[algo] = qps
			except Exception as e:
				print(f"[protein_search] ERROR running {algo.upper()}: {e}")
				all_qps[algo] = None
				continue

		print(f"\n{'='*60}")
		print(f"[protein_search] All methods complete!")
		print(f"{'='*60}")

		return all_qps
	
	# Handle 'nlsh' method separately
	if method == "nlsh":
		# 1. Create query .dat using protein_embed
		embed_py = os.path.join(os.path.dirname(__file__), "protein_embed.py")
		query_dat = os.path.splitext(output_txt)[0] + ".queries.dat"

		os.makedirs(os.path.dirname(output_txt), exist_ok=True)

		embed_cmd = [
			libraries.sys.executable,
			embed_py,
			"-i", query_fasta,
			"-o", query_dat,
		]

		print("[protein_search] Generating query embeddings (.dat)...")
		subprocess.run(embed_cmd, check=True)

		# 2. Build executable if missing (needed for NLSH graph build)
		alg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "Algorithms"))
		alg_part1 = os.path.join(alg_root, "AlgorithmsPart1")
		exe_path = os.path.join(alg_part1, "search")
		if not os.path.exists(exe_path):
			print("[protein_search] Building C executable ...")
			prev_cwd = os.getcwd()
			try:
				os.chdir(alg_root)
				if not build_executable():
					raise RuntimeError("Failed to build C executable")
			except Exception:
				os.chdir(prev_cwd)
				raise
			os.chdir(prev_cwd)

		qps = run_nlsh(
			base_dat=base_dat,
			query_dat=query_dat,
			output_txt=output_txt,
			N=N,
			R=R,
			seed=seed,
			nlsh_index=nlsh_index,
			nlsh_T=nlsh_T,
			nlsh_m=nlsh_m,
			nlsh_imbalance=nlsh_imbalance,
			nlsh_kahip_mode=nlsh_kahip_mode,
			nlsh_layers=nlsh_layers,
			nlsh_nodes=nlsh_nodes,
			nlsh_epochs=nlsh_epochs,
			nlsh_batch_size=nlsh_batch_size,
			nlsh_lr=nlsh_lr,
		)
		
		print(f"[protein_search] Done. Results at: {output_txt}")
		return qps
	
	# 1. Create query .dat using protein_embed
	embed_py = os.path.join(os.path.dirname(__file__), "protein_embed.py")
	query_dat = os.path.splitext(output_txt)[0] + ".queries.dat"

	os.makedirs(os.path.dirname(output_txt), exist_ok=True)

	embed_cmd = [
		libraries.sys.executable,
		embed_py,
		"-i", query_fasta,
		"-o", query_dat,
	]

	print("[protein_search] Generating query embeddings (.dat)...")
	subprocess.run(embed_cmd, check=True)

	# 2. Build executable if missing (needed for C-based methods and for NLSH graph build)
	alg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "Algorithms"))
	alg_part1 = os.path.join(alg_root, "AlgorithmsPart1")
	exe_path = os.path.join(alg_part1, "search")
	if not os.path.exists(exe_path):
		print("[protein_search] Building C executable ...")
		prev_cwd = os.getcwd()
		try:
			os.chdir(alg_root)
			if not build_executable():
				raise RuntimeError("Failed to build C executable")
		except Exception:
			os.chdir(prev_cwd)
			raise
		os.chdir(prev_cwd)

	# 3. Run selected ANN algorithm (cosine) on protein data
	# We invoke the C binary exactly like Part1 examples, but with -type protein
	cmd = [
		exe_path,
		"-d", os.path.abspath(base_dat),
		"-q", os.path.abspath(query_dat),
		"-o", os.path.abspath(output_txt),
		"-N", str(N),
		"-R", str(R),
		"-type", "protein",
		"-range", "false",
	]

	if method == "lsh":
		cmd.extend(["-k", str(k), "-L", str(L), "-w", str(w), "-lsh"])
		algo = "lsh"
	elif method == "hypercube":
		cmd.extend(["-kproj", str(kproj), "-w", str(w), "-M", str(M), "-probes", str(probes), "-hypercube"])
		algo = "hypercube"
	elif method == "ivfflat":
		cmd.extend(["-kclusters", str(kclusters), "-nprobe", str(nprobe), "-ivfflat"])
		algo = "ivfflat"
	elif method == "ivfpq":
		cmd.extend(["-kclusters", str(kclusters), "-nprobe", str(nprobe), "-M", str(ivfpq_M), "-nbits", str(nbits), "-ivfpq"])
		algo = "ivfpq"
	else:
		raise ValueError(f"Unknown method: {method}. Choose from: lsh, hypercube, ivfflat, ivfpq, nlsh")

	cmd.extend(["-seed", str(seed)])

	print(f"[protein_search] Running {method.upper()} on protein data ...")
	error_flag = run_algorithm(cmd)
	
	if not error_flag:	
		base_ids_txt = base_dat.replace(".dat", "_ids.txt")
		
		if not os.path.exists(base_ids_txt):
			raise ValueError(f"Unknown file: {base_ids_txt}.")
		
		query_ids_txt = query_dat.replace(".dat", "_ids.txt")
		remap_output_ids(output_txt, base_ids_txt, query_ids_txt)
		
		qps = None
		output_lines = []

		with open(output_txt, "r") as f:
			for line in f:
				if line.startswith("QPS:"):
					qps = float(line.split(":", 1)[1].strip())
				else:
					output_lines.append(line)

		with open(output_txt, "w") as f:
			f.writelines(output_lines)
			
		print(f"[protein_search] Done. Results at: {output_txt}")

		return qps
	
	if error_flag:
		return None


def main():
	parser = libraries.argparse.ArgumentParser("Protein ANN search wrapper")
	parser.add_argument("-d", required=True, help="Path to base protein vectors .dat (float32 N x 320)")
	parser.add_argument("-q", required=True, help="Path to FASTA with query sequences")
	parser.add_argument("-o", "--output", required=True, help="Output neighbors file (text)")
	parser.add_argument("-method", type=str, default="ivfflat", choices=["lsh", "hypercube", "ivfflat", "ivfpq", "neural", "all"],
						help="ANN algorithm to use (or 'all' to run all methods)")
	parser.add_argument("-N", type=int, default=10, help="Number of nearest neighbors")
	parser.add_argument("-R", type=float, default=0.5, help="Range search radius (for range search mode)")
	parser.add_argument("--seed", type=int, default=42, help="Random seed")
	
	# IVFFlat / IVFPQ params
	parser.add_argument("--kclusters", type=int, default=200, help="Number of clusters (ivfflat/ivfpq)")
	parser.add_argument("--nprobe", type=int, default=80, help="Number of probes (ivfflat/ivfpq)")
	
	# LSH params
	parser.add_argument("-k", type=int, default=4, help="Hash functions per table (lsh)")
	parser.add_argument("-L", type=int, default=70, help="Number of hash tables (lsh)")
	parser.add_argument("-w", type=float, default=4.0, help="Window width (lsh/hypercube)")
	
	# Hypercube params
	parser.add_argument("--kproj", type=int, default=14, help="Number of projections (hypercube)")
	parser.add_argument("-M", type=int, default=2000, help="Max candidates to check (hypercube)")
	parser.add_argument("--probes", type=int, default=100, help="Vertices to examine (hypercube)")
	
	# IVFPQ specific
	parser.add_argument("--nbits", type=int, default=8, help="Bits per subspace (ivfpq)")
	parser.add_argument("--ivfpq-M", type=int, default=16, help="Number of subvectors (ivfpq, must divide dimension)")

	# NLSH specific
	parser.add_argument("--nlsh-index", type=str, default=None, help="Directory for NLSH index (default: alongside output)")
	parser.add_argument("--nlsh-T", type=int, default=1500, help="Number of bins to probe (nlsh)")
	parser.add_argument("--nlsh-m", type=int, default=2000, help="Number of parts for KaHIP (nlsh build)")
	parser.add_argument("--nlsh-imbalance", type=float, default=0.1, help="KaHIP imbalance (nlsh build)")
	parser.add_argument("--nlsh-kahip-mode", type=int, default=0, help="KaHIP mode (nlsh build)")
	parser.add_argument("--nlsh-layers", type=int, default=10, help="MLP layers (nlsh build)")
	parser.add_argument("--nlsh-nodes", type=int, default=256, help="MLP hidden units (nlsh build)")
	parser.add_argument("--nlsh-epochs", type=int, default=8, help="Training epochs (nlsh build)")
	parser.add_argument("--nlsh-batch-size", type=int, default=512, help="Batch size (nlsh build)")
	parser.add_argument("--nlsh-lr", type=float, default=1e-3, help="Learning rate (nlsh build)")

	args = parser.parse_args()

	method_lower = args.method.lower()

	# The user will depict Neural LSH as neural
	# But we, throughout the code, we call it nlsh
	# For simplicity across the many programs
	# So let's change it!
	if method_lower == "neural":
		method_lower = "nlsh"

	print(method_lower)

	all_qps = run_protein_search(
		base_dat=args.d,
		query_fasta=args.q,
		output_txt=args.output,
		method=method_lower,
		N=args.N,
		R=args.R,
		seed=args.seed,
		kclusters=args.kclusters,
		nprobe=args.nprobe,
		k=args.k,
		L=args.L,
		w=args.w,
		kproj=args.kproj,
		M=args.M,
		probes=args.probes,
		nbits=args.nbits,
		ivfpq_M=args.ivfpq_M,
		nlsh_index=args.nlsh_index,
		nlsh_T=args.nlsh_T,
		nlsh_m=args.nlsh_m,
		nlsh_imbalance=args.nlsh_imbalance,
		nlsh_kahip_mode=args.nlsh_kahip_mode,
		nlsh_layers=args.nlsh_layers,
		nlsh_nodes=args.nlsh_nodes,
		nlsh_epochs=args.nlsh_epochs,
		nlsh_batch_size=args.nlsh_batch_size,
		nlsh_lr=args.nlsh_lr,
	)

	# Print All QPS status
	if isinstance(all_qps, dict):
		if method_lower == "all":
			print("\nQPS results:")

			for method, qps in all_qps.items():
				if method == "lsh":
					name = method.upper()  # LSH
				elif method == "ivfpq":
					name = method[:3].upper() + "-" + method[3:].upper()  # IVF-PQ
				elif method == "ivfflat":
					name = method[:3].upper() + "-" + method[3:].capitalize()  #IVF-Flat
				elif method == "nlsh":
					name = "Neural LSH"
				else:
					name = method.capitalize() # everything else
				
				print(f"  {name:10s}: {qps}")
	else:
		# Set the method's name in order to be printed correctly
		if method_lower == "lsh":
			name = method_lower.upper()  # LSH
		elif method_lower == "ivfpq":
			name = method_lower[:3].upper() + "-" + method_lower[3:].upper()  # IVF-PQ
		elif method_lower == "ivfflat":
			name = method_lower[:3].upper() + "-" + method_lower[3:].capitalize()  #IVF-Flat
		elif method_lower == "nlsh":
			name = "Neural LSH"
		else:
			name = method_lower.capitalize() # everything else

		print("\nQPS result:")
		print(f"  {name:10s}: {all_qps}")


	# Running BLAST command
	blast_results = f"output/blast/topN/blast_results_top{args.N}.tsv"
	blast_time = blast_executable(args.N, blast_results)

	# Read all query IDs from FASTA and
	# Find the total number of queries
	query_ids = []
	with open(args.q, "r") as f:
		for line in f:
			if line.startswith(">"):
				query_id = line.strip()[1:].split()[0]
				query_ids.append(query_id)

	total_proteins = len(query_ids)
	blast_time_per_query = blast_time / total_proteins
	blast_qps = total_proteins / blast_time

	# Parse BLAST results with identity
	blast_results_path = "output/blast/search/blast_results.tsv"
	blast_identity = parse_blast_results_with_identity(blast_results_path) if os.path.exists(blast_results_path) else {}

	# Compute Recall against BLAST results (topN.tsv) vs the results.txt
	# If method.lower == "all" then compare all the results_algo.txt vs the BLAST_results.tsv
	if method_lower == "all":
		all_methods = ["lsh", "hypercube", "ivfflat", "ivfpq", "nlsh"]
		all_methods_display = ["Euclidean LSH", "Hypercube", "IVF-Flat", "IVF-PQ", "Neural LSH"]
		base_output = os.path.splitext(args.output)[0]

		print("\nRecall results (average, for logging only):")
		method_recall = {}
		method_results = {}

		for algo, algo_display in zip(all_methods, all_methods_display):
			ann_txt = f"{base_output}_{algo}.txt"
			mean_recall, _ = compute_recall(blast_tsv=blast_results, ann_txt=ann_txt, topN=args.N)
			method_recall[algo] = mean_recall

			# Parse neighbor results for formatted report
			method_results[algo] = parse_neighbor_results(ann_txt, args.N)

			if algo == "lsh" or algo == "nlsh" or algo == "ivfpq":
				print(f"  {algo.upper():10s}: {mean_recall:.4f}")  # LSH / NLSH / IVFPQ
			elif algo == "ivfflat":
				algo_display_print = algo[:4].upper() + algo[4:]  # IVFFlat
				print(f"  {algo_display_print:10s}: {mean_recall:.4f}")
			else:
				print(f"  {algo.capitalize():10s}: {mean_recall:.4f}")

		# Parse BLAST top-N
		blast_gt = parse_blast_tsv(blast_results, args.N)

		# Generate a single consolidated report for all methods
		# report_path = base_output + "_REPORT.txt"
		generate_all_methods_report(
			output_report=args.output,
			query_ids=query_ids,
			topN=args.N,
			method_qps=all_qps,
			method_results=method_results,
			blast_results_topn=blast_gt,
			blast_results_identity=blast_identity,
			blast_time_per_query = blast_time_per_query,
			blast_qps = blast_qps
		)
	else:
		mean_recall, _ = compute_recall(blast_tsv=blast_results, ann_txt=args.output, topN=args.N)
		print(f"\nRecall@{args.N} (average): {mean_recall:.4f}")

		# Parse neighbor results
		method_results = parse_neighbor_results(args.output, args.N)

		# Parse BLAST top-N
		blast_gt = parse_blast_tsv(blast_results, args.N)

		# Determine method display name
		method_display_map = {
			"lsh": "Euclidean LSH",
			"hypercube": "Hypercube",
			"ivfflat": "IVF-Flat",
			"ivfpq": "IVF-PQ",
			"nlsh": "Neural LSH",
		}
		method_display = method_display_map.get(method_lower, method_lower.capitalize())

		# Generate per-query report
		# report_path = os.path.splitext(args.output)[0] + "_REPORT.txt"
		generate_per_query_report(
			output_report=args.output,
			query_ids=query_ids,
			topN=args.N,
			method_name=method_display,
			method_key=method_lower,
			qps=all_qps,
			method_results=method_results,
			blast_results_topn=blast_gt,
			blast_results_identity=blast_identity,
			blast_time_per_query = blast_time_per_query,
			blast_qps = blast_qps
		)


	# Deleting unecessary files
	if method_lower == "all":
		all_methods = ["lsh", "hypercube", "ivfflat", "ivfpq", "nlsh"]
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
	


if __name__ == "__main__":
	main()