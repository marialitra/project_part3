import libraries
from libraries import os, subprocess

# Local imports from the second prt of the project
libraries.sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "Algorithms", "src"))
from runSearchExe import build_executable, run_algorithm


def run_nlsh(
	base_dat: str,
	query_dat: str,
	output_txt: str,
	N: int,
	R: float,
	range: bool,
	seed: int,
	nlsh_index: str,
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
	"""
		Build (if missing) and run NLSH on protein .dat data.
	"""

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
		"-range", str(range).lower(),
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
	
	libraries.remap_output_ids(output_txt, base_ids_txt, query_ids_txt)

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
		method: str,
		N: int,
		R: float,
		range: bool,
		seed: int,
		k: int,
		L: int,
		lsh_w: float,
		kproj: int,
		hyper_w: float,
		hyper_M: int,
		probes: int,
		flat_kclusters: int,
		flat_nprobe: int,
		pq_kclusters: int,
		pq_nprobe: int,
		pq_M: int,
		nbits: int,
		nlsh_index: str,
		nlsh_T: int,
		nlsh_m: int,
		nlsh_imbalance: float,
		nlsh_kahip_mode: int,
		nlsh_layers: int,
		nlsh_nodes: int,
		nlsh_epochs: int,
		nlsh_batch_size: int,
		nlsh_lr: float
):
	"""
		1) Build query embeddings via Part3/src/protein_embed.py -> query.dat
		2) Build C executable if needed
		3) Run selected ANN algorithm on protein data, writing results to output_txt
	"""

	# Handle 'all' method by recursively calling for each algorithm
	if method == "all":
		all_methods = ["lsh", "hypercube", "ivfflat", "ivfpq", "nlsh"]
		all_qps = {}
		base_output = os.path.splitext(output_txt)[0]
		
		for algo in all_methods:
			algo_output = f"{base_output}_{algo}.txt"
			print(f"\n{'='*75}")
			print(f"[protein_search] Running method: {algo.upper()}")
			print(f"{'='*75}")
			
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
						range=range,
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
						range=range,
						seed=seed,
						k=k,
						L=L,
						lsh_w=lsh_w,
						kproj=kproj,
						hyper_w=hyper_w,
						hyper_M=hyper_M,
						probes=probes,
						flat_kclusters=flat_kclusters,
						flat_nprobe=flat_nprobe,
						pq_kclusters=pq_kclusters,
						pq_nprobe=pq_nprobe,
						pq_M=pq_M,
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
						nlsh_lr=nlsh_lr
					)
					all_qps[algo] = qps
						
			except Exception as e:
				print(f"[protein_search] ERROR running {algo.upper()}: {e}")
				all_qps[algo] = None
				continue

		print(f"\n{'='*75}")
		flag = 0
		for name, qps in all_qps.items():
			if qps == None:
				print(f"[protein_search] Not Completed. Method: {name}, did not run.")
				flag = 1
			else:
				print(f"[protein_search] Done. Method: {name}, run successfully!")


		if flag == 0:
			print(f"[protein_search] All methods complete!")

		print(f"{'='*75}")

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
			output_txt=os.path.abspath(output_txt),
			N=N,
			R=R,
			range=range,
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
		
		print(f"\n{'='*75}")
		if qps == None:
			print(f"[protein_search] Not Completed. Method: {method}, did not run.")
		else:
			print(f"[protein_search] Done. Results at: {output_txt}")		
		print(f"{'='*75}")
	
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
		"-range", str(range).lower(),
		"-seed", str(seed),
	]

	if method == "lsh":
		cmd.extend(["-k", str(k), "-L", str(L), "-w", str(lsh_w), "-lsh"])
		algo = "lsh"
	elif method == "hypercube":
		cmd.extend(["-kproj", str(kproj), "-w", str(hyper_w), "-M", str(hyper_M), "-probes", str(probes), "-hypercube"])
		algo = "hypercube"
	elif method == "ivfflat":
		cmd.extend(["-kclusters", str(flat_kclusters), "-nprobe", str(flat_nprobe), "-ivfflat"])
		algo = "ivfflat"
	elif method == "ivfpq":
		cmd.extend(["-kclusters", str(pq_kclusters), "-nprobe", str(pq_nprobe), "-M", str(pq_M), "-nbits", str(nbits), "-ivfpq"])
		algo = "ivfpq"
	else:
		raise ValueError(f"Unknown method: {method}. Choose from: lsh, hypercube, ivfflat, ivfpq, nlsh")

	print(f"[protein_search] Running {method.upper()} on protein data ...")
	error_flag = run_algorithm(cmd)
	
	if not error_flag:	
		base_ids_txt = base_dat.replace(".dat", "_ids.txt")
		
		if not os.path.exists(base_ids_txt):
			raise ValueError(f"Unknown file: {base_ids_txt}.")
		
		query_ids_txt = query_dat.replace(".dat", "_ids.txt")
		libraries.remap_output_ids(output_txt, base_ids_txt, query_ids_txt)
		
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
			
		print(f"\n{'='*75}")
		print(f"[protein_search] Done. Results at: {output_txt}")
		print(f"{'='*75}")

		return qps
	
	if error_flag:
		print(f"\n{'='*75}")
		print(f"[protein_search] Not Completed. Method: {method}, did not run.")
		print(f"{'='*75}")

		return None




def run_method(args, method, current_methods):
	# -------------------------------------
	# Run the method specified by the user
	# -------------------------------------
	all_qps = run_protein_search(
		base_dat=args.d,
		query_fasta=args.q,
		output_txt=args.output,
		method=method,
		N=args.N,
		R=args.R,
		range=args.range,
		seed=args.seed,
		k=args.k,
		L=args.L,
		lsh_w=args.lsh_w,
		kproj=args.kproj,
		hyper_w=args.hyper_w,
		hyper_M=args.hyper_M,
		probes=args.probes,
		flat_kclusters=args.flat_kclusters,
		flat_nprobe=args.flat_nprobe,
		pq_kclusters=args.pq_kclusters,
		pq_nprobe=args.pq_nprobe,
		pq_M=args.pq_M,
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

	# If results.txt exist, keep the old verion
	# As the algorithm did not run, delete it!
	if all_qps == None:
		base_output = os.path.splitext(args.output)[0]

		if os.path.exists(f"{base_output}.queries.dat"):
			os.remove(f"{base_output}.queries.dat")
		
		if os.path.exists(f"{base_output}.queries_ids.txt"):
			os.remove(f"{base_output}.queries_ids.txt")

		if os.path.exists(f"{args.output}"):
			os.remove(f"{args.output}")
	
	
	# Print QPS status in terminal
	answer = []
	if isinstance(all_qps, dict):
		answer = libraries.print_QPS_dictionary(all_qps, method, current_methods)
		if answer == []:
			print("QPS: []")
	else:
		libraries.print_QPS(all_qps, method)

	return all_qps, answer
		


        



