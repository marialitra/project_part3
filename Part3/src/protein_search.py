import os
import argparse
import subprocess
import sys
from typing import Optional

# Local imports
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "Algorithms", "src"))
from runSearchExe import build_executable, run_algorithm


def blast_executable(neighbors: int, output: str):
    """
        Runs 'make blast' to run the BLAST method and its filtering.
    """

    try:
        build_process = subprocess.run(["make", "blast", f"N={neighbors}", f"out={output}"], capture_output=True, text=True, check=True)
        print("\n\nBuild complete: BLAST results are ready.")
        return True
    except subprocess.CalledProcessError as e:
        print("\n\n--- ERROR: Build failed. ---")
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        return False
    except FileNotFoundError:
         print("\n\n--- ERROR: 'make' command not found. Is it installed? ---")
         return False

# def compute_recall():


def remap_output_ids(output_txt: str, base_ids_txt: str, query_ids_txt: str):
	"""Remap numeric indices in output to actual protein IDs."""
	if not os.path.exists(output_txt):
		return
	
	# Load ID mappings
	base_ids = []
	if os.path.exists(base_ids_txt):
		with open(base_ids_txt, 'r') as f:
			base_ids = [line.strip() for line in f if line.strip()]
	
	query_ids = []
	if os.path.exists(query_ids_txt):
		with open(query_ids_txt, 'r') as f:
			query_ids = [line.strip() for line in f if line.strip()]
	
	if not base_ids and not query_ids:
		return
	
	# Read and remap output
	with open(output_txt, 'r') as f:
		lines = f.readlines()
	
	remapped_lines = []
	for line in lines:
		if line.startswith("Query:"):
			# Extract query index and remap to ID
			parts = line.split()
			if len(parts) >= 2:
				try:
					idx = int(parts[1])
					if 0 <= idx < len(query_ids):
						line = f"Query: {query_ids[idx]}\n"
				except (ValueError, IndexError):
					pass
		elif line.startswith("Nearest neighbor"):
			# Extract neighbor index and remap to ID
			parts = line.split()
			if len(parts) >= 2:
				try:
					idx = int(parts[-1])
					if 0 <= idx < len(base_ids):
						prefix = " ".join(parts[:-1])
						line = f"{prefix} {base_ids[idx]}\n"
				except (ValueError, IndexError):
					pass
		remapped_lines.append(line)
	
	# Write back
	with open(output_txt, 'w') as f:
		f.writelines(remapped_lines)
	
	print(f"[protein_search] Remapped indices to protein IDs in {output_txt}")

def run_nlsh_pipeline(
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
			sys.executable,
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
		sys.executable,
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
	kclusters: int = 1024,
	nprobe: int = 10,
	# LSH params
	k: int = 4,
	L: int = 5,
	w: float = 4.0,
	# Hypercube params
	kproj: int = 12,
	M: int = 5000,
	probes: int = 10,
	# IVFPQ specific
	nbits: int = 8,
	# NLSH specific
	nlsh_index: Optional[str] = None,
	nlsh_T: int = 50,
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
		method: ANN algorithm to use ('lsh', 'hypercube', 'ivfflat', 'ivfpq', 'nlsh')
	"""

	# Handle 'all' method by recursively calling for each algorithm
	if method.lower() == "all":
		all_methods = ["lsh", "hypercube", "ivfflat", "ivfpq", "nlsh"]
		all_qps = {}
		base_output = os.path.splitext(output_txt)[0]
		
		for algo in all_methods:
			algo_output = f"{base_output}_{algo}.txt"
			print(f"\n{'='*60}")
			print(f"[protein_search] Running method: {algo.upper()}")
			print(f"{'='*60}")
			
			try:
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
	
	# 1. Create query .dat using protein_embed
	embed_py = os.path.join(os.path.dirname(__file__), "protein_embed.py")
	query_dat = os.path.splitext(output_txt)[0] + ".queries.dat"

	os.makedirs(os.path.dirname(output_txt), exist_ok=True)

	embed_cmd = [
		sys.executable,
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
	method_lower = method.lower()
	if method_lower == "nlsh":
		qps = run_nlsh_pipeline(
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

	if method_lower == "lsh":
		cmd.extend(["-k", str(k), "-L", str(L), "-w", str(w), "-lsh"])
		algo = "lsh"
	elif method_lower == "hypercube":
		cmd.extend(["-kproj", str(kproj), "-w", str(w), "-M", str(M), "-probes", str(probes), "-hypercube"])
		algo = "hypercube"
	elif method_lower == "ivfflat":
		cmd.extend(["-kclusters", str(kclusters), "-nprobe", str(nprobe), "-ivfflat"])
		algo = "ivfflat"
	elif method_lower == "ivfpq":
		cmd.extend(["-kclusters", str(kclusters), "-nprobe", str(nprobe), "-M", str(M), "-nbits", str(nbits), "-ivfpq"])
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
	parser = argparse.ArgumentParser("Protein ANN search wrapper")
	parser.add_argument("-d", required=True, help="Path to base protein vectors .dat (float32 N x 320)")
	parser.add_argument("-q", required=True, help="Path to FASTA with query sequences")
	parser.add_argument("-o", "--output", required=True, help="Output neighbors file (text)")
	parser.add_argument("-method", type=str, default="ivfflat", choices=["lsh", "hypercube", "ivfflat", "ivfpq", "nlsh", "all"],
						help="ANN algorithm to use (or 'all' to run all methods)")
	parser.add_argument("-N", type=int, default=5, help="Number of nearest neighbors")
	parser.add_argument("-R", type=float, default=0.5, help="Range search radius (for range search mode)")
	parser.add_argument("--seed", type=int, default=42, help="Random seed")
	
	# IVFFlat / IVFPQ params
	parser.add_argument("--kclusters", type=int, default=50, help="Number of clusters (ivfflat/ivfpq)")
	parser.add_argument("--nprobe", type=int, default=5, help="Number of probes (ivfflat/ivfpq)")
	
	# LSH params
	parser.add_argument("-k", type=int, default=4, help="Hash functions per table (lsh)")
	parser.add_argument("-L", type=int, default=5, help="Number of hash tables (lsh)")
	parser.add_argument("-w", type=float, default=4.0, help="Window width (lsh/hypercube)")
	
	# Hypercube params
	parser.add_argument("--kproj", type=int, default=12, help="Number of projections (hypercube)")
	parser.add_argument("-M", type=int, default=500, help="Max candidates to check (hypercube) or subvectors (ivfpq)")
	parser.add_argument("--probes", type=int, default=2, help="Vertices to examine (hypercube)")
	
	# IVFPQ specific
	parser.add_argument("--nbits", type=int, default=8, help="Bits per subspace (ivfpq)")

	# NLSH specific
	parser.add_argument("--nlsh-index", type=str, default=None, help="Directory for NLSH index (default: alongside output)")
	parser.add_argument("--nlsh-T", type=int, default=50, help="Number of bins to probe (nlsh)")
	parser.add_argument("--nlsh-m", type=int, default=2000, help="Number of parts for KaHIP (nlsh build)")
	parser.add_argument("--nlsh-imbalance", type=float, default=0.1, help="KaHIP imbalance (nlsh build)")
	parser.add_argument("--nlsh-kahip-mode", type=int, default=0, help="KaHIP mode (nlsh build)")
	parser.add_argument("--nlsh-layers", type=int, default=3, help="MLP layers (nlsh build)")
	parser.add_argument("--nlsh-nodes", type=int, default=256, help="MLP hidden units (nlsh build)")
	parser.add_argument("--nlsh-epochs", type=int, default=5, help="Training epochs (nlsh build)")
	parser.add_argument("--nlsh-batch-size", type=int, default=512, help="Batch size (nlsh build)")
	parser.add_argument("--nlsh-lr", type=float, default=1e-3, help="Learning rate (nlsh build)")

	args = parser.parse_args()
	all_qps = run_protein_search(
		base_dat=args.d,
		query_fasta=args.q,
		output_txt=args.output,
		method=args.method,
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

	# Deleting unecessary files
	method = args.method
	
	if method.lower() == "all":
		all_methods = ["lsh", "hypercube", "ivfflat", "ivfpq", "nlsh"]
		base_output = os.path.splitext(args.output)[0]
		
		for algo in all_methods:
			os.remove(f"{base_output}_{algo}.queries_ids.txt")
			os.remove(f"{base_output}_{algo}.queries.dat")
	else:
		os.remove(os.path.splitext(args.output)[0] + ".queries.dat")
		os.remove(os.path.splitext(args.output)[0] + ".queries_ids.txt")
	

	# Print All QPS status
	if isinstance(all_qps, dict):
		if method.lower() == "all":
			print("\nQPS results:")

			for method, qps in all_qps.items():
				print(f"  {method.upper():10s}: {qps}")
	else:
		print("\nQPS result:")
		
		method = method.capitalize()
		print(f"{method}: {all_qps}")


	# Running BLAST command
	output = f"output/blast/topN/blast_results_top{args.N}.tsv"
	blast_executable(args.N, output)

	# Compute Recall against BLAST results (topN.tsv) vs the results.txt
	# If method.lower == "all" then compare all the results_algo.txt vs the BLAST_results.tsv
	
	# compute_recall()


if __name__ == "__main__":
	main()