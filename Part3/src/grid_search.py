#!/usr/bin/env python3
"""
Grid search for ANN hyperparameter tuning - Version 2
Runs all method/parameter combinations and saves results to CSV.
Supports appending to existing results.
"""

import os
import sys
import csv
import subprocess
import argparse
from itertools import product
from typing import Dict, List, Tuple, Optional, Set


# Local imports
sys.path.append(os.path.dirname(__file__))
from libraries import run_protein_search, calculate_recall



def ensure_blast_results(N_vals: List[int]):
	"""Ensure BLAST top-N results exist for all N values."""
	print(f"\n{'='*80}")
	print("ENSURING BLAST RESULTS EXIST")
	print(f"{'='*80}")
	
	for N in N_vals:
		blast_tsv = f"output/blast/topN/blast_results_top{N}.tsv"
		if os.path.exists(blast_tsv):
			print(f"✓ BLAST top-{N} results exist: {blast_tsv}")
		else:
			print(f"✗ BLAST top-{N} results missing. Generating...")
			try:
				subprocess.run(
					["make", "blast", f"N={N}"],
					check=True,
					cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
				)
				print(f"✓ Generated BLAST top-{N} results")
			except subprocess.CalledProcessError as e:
				print(f"✗ ERROR generating BLAST top-{N}: {e}")
				raise


def load_existing_results(csv_path: str) -> Tuple[List[Dict], Set[Tuple]]:
	"""Load existing CSV results and return (results, completed_param_sets)."""
	if not os.path.exists(csv_path):
		return [], set()
	
	results = []
	completed = set()
	
	with open(csv_path, "r", newline="") as f:
		reader = csv.DictReader(f)
		for row in results:
			results.append(row)
			# Create tuple of all params except recall/qps/time_per_query
			param_tuple = tuple(
				(k, v) for k, v in sorted(row.items())
				if k not in ["recall", "qps", "time_per_query"]
			)
			completed.add(param_tuple)
	
	return results, completed


def get_lsh_grid():
	"""Generate LSH parameter grid: k, L, w"""
	k_vals = [2, 4, 6]
	L_vals = [5, 10, 15]
	w_vals = [20, 40, 80, 120, 300, 500, 1000]
	
	for k, L, w in product(k_vals, L_vals, w_vals):
		yield {
			"k": k,
			"L": L,
			"lsh_w": w,
		}


def get_hypercube_grid():
	"""Generate Hypercube parameter grid: kproj, w, M, probes"""
	kproj_vals = [12, 14, 16, 20]
	w_vals = [20, 40, 80, 120, 300, 500, 1000]
	M_vals = [2000, 5000, 10000]
	probes_vals = [100, 250, 500]
	
	for kproj, w, M, probes in product(kproj_vals, w_vals, M_vals, probes_vals):
		yield {
			"kproj": kproj,
			"hyper_w": w,
			"hyper_M": M,
			"probes": probes,
		}


def get_ivfflat_grid():
	"""Generate IVF-Flat parameter grid: kclusters, nprobe"""
	kclusters_vals = [200, 500, 1000, 5000]
	nprobe_vals = [50, 100, 250, 500, 1500]
	
	for kclusters, nprobe in product(kclusters_vals, nprobe_vals):
		# Only yield if nprobe < kclusters
		if nprobe < kclusters:
			yield {
				"flat_kclusters": kclusters,
				"flat_nprobe": nprobe,
			}


def get_ivfpq_grid():
	"""Generate IVF-PQ parameter grid: kclusters, nprobe, M (subvectors)"""
	kclusters_vals = [200, 500, 1000, 5000]
	nprobe_vals = [50, 100, 250, 500, 1500]
	M_vals = [8, 16]
	
	for kclusters, nprobe, M in product(kclusters_vals, nprobe_vals, M_vals):
		# Only yield if nprobe < kclusters
		if nprobe < kclusters:
			yield {
				"pq_kclusters": kclusters,
				"pq_nprobe": nprobe,
				"pq_M": M,
			}


def get_nlsh_grid():
	"""Generate NLSH parameter grid: T, m, layers, nodes"""
	T_vals = [1000, 1500]
	m_vals = [1800, 2000, 2200]
	layers_vals = [5, 10, 15]
	nodes_vals = [128, 256, 512]
	
	for T, m, layers, nodes in product(T_vals, m_vals, layers_vals, nodes_vals):
		yield {
			"nlsh_T": T,
			"nlsh_m": m,
			"nlsh_layers": layers,
			"nlsh_nodes": nodes,
		}


def run_grid_search(
	base_dat: str,
	query_fasta: str,
	output_dir: str,
	method: str = "all",
	seed: int = 42,
	nbits: int = 8,
	nlsh_epochs: int = 8,
	nlsh_batch_size: int = 512,
	nlsh_lr: float = 1e-3,
	append: bool = True,
):
	"""
	Run grid search for specified method(s).
	
	Args:
		base_dat: Path to base protein vectors .dat
		query_fasta: Path to query FASTA
		output_dir: Directory to save results and CSVs
		method: "all" or specific method name
		seed: Random seed
		nbits: Bits per subspace (IVF-PQ)
		nlsh_epochs: Training epochs (NLSH)
		nlsh_batch_size: Batch size (NLSH)
		nlsh_lr: Learning rate (NLSH)
		append: If True, append to existing CSVs. If False, overwrite.
	"""
	
	os.makedirs(output_dir, exist_ok=True)
	
	# BLAST reference results
	blast_results_path = "output/blast/topN"
	
	N_vals = [1, 10, 50]
	seed_val = seed
	
	# Ensure BLAST results exist for all N values
	ensure_blast_results(N_vals)
	
	# Dictionary to store results per method
	all_results = {
		"lsh": [],
		"hypercube": [],
		"ivfflat": [],
		"ivfpq": [],
		"nlsh": [],
	}
	
	# Load existing results if appending
	if append:
		for m in all_results.keys():
			csv_path = os.path.join(output_dir, f"grid_search_{m}.csv")
			if os.path.exists(csv_path):
				with open(csv_path, "r", newline="") as f:
					reader = csv.DictReader(f)
					all_results[m] = list(reader)
				print(f"Loaded {len(all_results[m])} existing results for {m}")
	
	methods_to_run = [method] if method != "all" else ["lsh", "hypercube", "ivfflat", "ivfpq", "nlsh"]
	
	run_count = 0

	# LSH grid search
	if "lsh" in methods_to_run:
		print(f"\n{'='*80}")
		print("LSH GRID SEARCH")
		print(f"{'='*80}")
		
		for N in N_vals:
			for param_set in get_lsh_grid():
				run_count += 1
				
				# Check if already exists
				already_exists = False
				if append:
					for existing in all_results["lsh"]:
						if (int(existing["N"]) == N and 
							int(existing["seed"]) == seed_val and
							int(existing["k"]) == param_set["k"] and
							int(existing["L"]) == param_set["L"] and
							float(existing["lsh_w"]) == param_set["lsh_w"]):
							already_exists = True
							break
				
				if already_exists:
					print(f"\n[{run_count}] LSH: N={N}, k={param_set['k']}, L={param_set['L']}, lsh_w={param_set['lsh_w']} - SKIPPING (already exists)")
					continue
				
				param_str = f"N{N}_k{param_set['k']}_L{param_set['L']}_lsh_w{int(param_set['lsh_w'])}"
				output_txt = os.path.join(output_dir, f"lsh_{param_str}.txt")
				
				print(f"\n[{run_count}] LSH: N={N}, k={param_set['k']}, L={param_set['L']}, lsh_w={param_set['lsh_w']}")
				
				try:
					qps = run_protein_search(
						base_dat=base_dat,
						query_fasta=query_fasta,
						output_txt=output_txt,
						method="lsh",
						N=N,
						seed=seed_val,
						k=param_set['k'],
						L=param_set['L'],
						lsh_w=param_set['lsh_w'],
						R=0.5,
						range=False,
						kproj=0,
						hyper_w=0.0,
						hyper_M=0,
						probes=0,
						flat_kclusters=0,
						flat_nprobe=0,
						pq_kclusters=0,
						pq_nprobe=0,
						pq_M=0,
						nbits=0,
						nlsh_index="",
						nlsh_T=0,
						nlsh_m=0,
						nlsh_imbalance=0.0,
						nlsh_kahip_mode=0,
						nlsh_layers=0,
						nlsh_nodes=0,
						nlsh_epochs=0,
						nlsh_batch_size=0,
						nlsh_lr=0.0
					)
					
					# Compute recall
					blast_tsv = os.path.join(blast_results_path, f"blast_results_top{N}.tsv")
					if os.path.exists(blast_tsv) and os.path.exists(output_txt):
						mean_recall, _ = calculate_recall(blast_tsv, output_txt, N)
						
						result = {
							"N": N,
							"seed": seed_val,
							"k": param_set['k'],
							"L": param_set['L'],
							"lsh_w": param_set['lsh_w'],
							"recall": mean_recall,
							"qps": qps if qps else 0.0,
							"time_per_query": 1.0/qps if qps and qps > 0 else float('inf'),
						}
						all_results["lsh"].append(result)
						
						#print(f"  Recall@{N}: {mean_recall:.4f}, QPS: {qps:.2f if qps else 0}")
						# Chnged the baove line with the 2 above this sentence due to:
						# ERROR: Invalid Format Specifier
						qps_val = qps if qps else 0.0
						print(f"  Recall@{N}: {mean_recall:.4f}, QPS: {qps_val:.2f}")
					
					# Clean up
					if os.path.exists(output_txt):
						os.remove(output_txt)
					query_dat = output_txt.replace(".txt", ".queries.dat")
					query_ids = output_txt.replace(".txt", ".queries_ids.txt")
					if os.path.exists(query_dat):
						os.remove(query_dat)
					if os.path.exists(query_ids):
						os.remove(query_ids)
					
				except Exception as e:
					print(f"  ERROR: {e}")
					continue
	
	# Hypercube grid search
	if "hypercube" in methods_to_run:
		print(f"\n{'='*80}")
		print("HYPERCUBE GRID SEARCH")
		print(f"{'='*80}")
		
		for N in N_vals:
			for param_set in get_hypercube_grid():
				run_count += 1
				
				# Check if already exists
				already_exists = False
				if append:
					for existing in all_results["hypercube"]:
						if (int(existing["N"]) == N and 
							int(existing["seed"]) == seed_val and
							int(existing["kproj"]) == param_set["kproj"] and
							float(existing["hyper_w"]) == param_set["hyper_w"] and
							int(existing["hyper_M"]) == param_set["hyper_M"] and
							int(existing["probes"]) == param_set["probes"]):
							already_exists = True
							break
				
				if already_exists:
					print(f"\n[{run_count}] Hypercube: N={N}, kproj={param_set['kproj']}, hyper_w={param_set['hyper_w']}, hyper_M={param_set['hyper_M']}, probes={param_set['probes']} - SKIPPING (already exists)")
					continue
				
				param_str = f"N{N}_kproj{param_set['kproj']}_hyper_w{int(param_set['hyper_w'])}_hyper_M{param_set['hyper_M']}_probes{param_set['probes']}"
				output_txt = os.path.join(output_dir, f"hypercube_{param_str}.txt")
				
				print(f"\n[{run_count}] Hypercube: N={N}, kproj={param_set['kproj']}, hyper_w={param_set['hyper_w']}, hyper_M={param_set['hyper_M']}, probes={param_set['probes']}")
				
				try:
					qps = run_protein_search(
						base_dat=base_dat,
						query_fasta=query_fasta,
						output_txt=output_txt,
						method="hypercube",
						N=N,
						seed=seed_val,
						kproj=param_set['kproj'],
						hyper_w=param_set['hyper_w'],
						hyper_M=param_set['hyper_M'],
						probes=param_set['probes'],
						R=0.5,
						range=False,
						k=0,
						L=0,
						lsh_w=0.0,
						flat_kclusters=0,
						flat_nprobe=0,
						pq_kclusters=0,
						pq_nprobe=0,
						pq_M=0,
						nbits=0,
						nlsh_index="",
						nlsh_T=0,
						nlsh_m=0,
						nlsh_imbalance=0.0,
						nlsh_kahip_mode=0,
						nlsh_layers=0,
						nlsh_nodes=0,
						nlsh_epochs=0,
						nlsh_batch_size=0,
						nlsh_lr=0.0
					)
					
					# Compute recall
					blast_tsv = os.path.join(blast_results_path, f"blast_results_top{N}.tsv")
					if os.path.exists(blast_tsv) and os.path.exists(output_txt):
						mean_recall, _ = calculate_recall(blast_tsv, output_txt, N)
						
						result = {
							"N": N,
							"seed": seed_val,
							"kproj": param_set['kproj'],
							"hyper_w": param_set['hyper_w'],
							"hyper_M": param_set['hyper_M'],
							"probes": param_set['probes'],
							"recall": mean_recall,
							"qps": qps if qps else 0.0,
							"time_per_query": 1.0/qps if qps and qps > 0 else float('inf'),
						}
						all_results["hypercube"].append(result)
						
						#print(f"  Recall@{N}: {mean_recall:.4f}, QPS: {qps:.2f if qps else 0}")
						# Chnged the baove line with the 2 above this sentence due to:
						# ERROR: Invalid Format Specifier
						qps_val = qps if qps else 0.0
						print(f"  Recall@{N}: {mean_recall:.4f}, QPS: {qps_val:.2f}")
					
					# Clean up
					if os.path.exists(output_txt):
						os.remove(output_txt)
					query_dat = output_txt.replace(".txt", ".queries.dat")
					query_ids = output_txt.replace(".txt", ".queries_ids.txt")
					if os.path.exists(query_dat):
						os.remove(query_dat)
					if os.path.exists(query_ids):
						os.remove(query_ids)
					
				except Exception as e:
					print(f"  ERROR: {e}")
					continue
	
	# IVF-Flat grid search
	if "ivfflat" in methods_to_run:
		print(f"\n{'='*80}")
		print("IVF-FLAT GRID SEARCH")
		print(f"{'='*80}")
		
		for N in N_vals:
			for param_set in get_ivfflat_grid():
				run_count += 1
				
				# Check if already exists
				already_exists = False
				if append:
					for existing in all_results["ivfflat"]:
						if (int(existing["N"]) == N and 
							int(existing["seed"]) == seed_val and
							int(existing["flat_kclusters"]) == param_set["flat_kclusters"] and
							int(existing["flat_nprobe"]) == param_set["flat_nprobe"]):
							already_exists = True
							break
				
				if already_exists:
					print(f"\n[{run_count}] IVF-Flat: N={N}, flat_kclusters={param_set['flat_kclusters']}, flat_nprobe={param_set['flat_nprobe']} - SKIPPING (already exists)")
					continue
				
				param_str = f"N{N}_flat_kclusters{param_set['flat_kclusters']}_flat_nprobe{param_set['flat_nprobe']}"
				output_txt = os.path.join(output_dir, f"ivfflat_{param_str}.txt")
				
				print(f"\n[{run_count}] IVF-Flat: N={N}, flat_kclusters={param_set['flat_kclusters']}, flat_nprobe={param_set['flat_nprobe']}")
				
				try:
					qps = run_protein_search(
						base_dat=base_dat,
						query_fasta=query_fasta,
						output_txt=output_txt,
						method="ivfflat",
						N=N,
						seed=seed_val,
						flat_kclusters=param_set['flat_kclusters'],
						flat_nprobe=param_set['flat_nprobe'],
						R=0.5,
						range=False,
						k=0,
						L=0,
						lsh_w=0.0,
						kproj=0,
						hyper_w=0.0,
						hyper_M=0,
						probes=0,
						pq_kclusters=0,
						pq_nprobe=0,
						pq_M=0,
						nbits=0,
						nlsh_index="",
						nlsh_T=0,
						nlsh_m=0,
						nlsh_imbalance=0.0,
						nlsh_kahip_mode=0,
						nlsh_layers=0,
						nlsh_nodes=0,
						nlsh_epochs=0,
						nlsh_batch_size=0,
						nlsh_lr=0.0
					)
					
					# Compute recall
					blast_tsv = os.path.join(blast_results_path, f"blast_results_top{N}.tsv")
					if os.path.exists(blast_tsv) and os.path.exists(output_txt):
						mean_recall, _ = calculate_recall(blast_tsv, output_txt, N)
						
						result = {
							"N": N,
							"seed": seed_val,
							"flat_kclusters": param_set['flat_kclusters'],
							"flat_nprobe": param_set['flat_nprobe'],
							"recall": mean_recall,
							"qps": qps if qps else 0.0,
							"time_per_query": 1.0/qps if qps and qps > 0 else float('inf'),
						}
						all_results["ivfflat"].append(result)
						
						#print(f"  Recall@{N}: {mean_recall:.4f}, QPS: {qps:.2f if qps else 0}")
						# Chnged the baove line with the 2 above this sentence due to:
						# ERROR: Invalid Format Specifier
						qps_val = qps if qps else 0.0
						print(f"  Recall@{N}: {mean_recall:.4f}, QPS: {qps_val:.2f}")
					
					# Clean up
					if os.path.exists(output_txt):
						os.remove(output_txt)
					query_dat = output_txt.replace(".txt", ".queries.dat")
					query_ids = output_txt.replace(".txt", ".queries_ids.txt")
					if os.path.exists(query_dat):
						os.remove(query_dat)
					if os.path.exists(query_ids):
						os.remove(query_ids)
					
				except Exception as e:
					print(f"  ERROR: {e}")
					continue
	
	# IVF-PQ grid search
	if "ivfpq" in methods_to_run:
		print(f"\n{'='*80}")
		print("IVF-PQ GRID SEARCH")
		print(f"{'='*80}")
		
		for N in N_vals:
			for param_set in get_ivfpq_grid():
				run_count += 1
				
				# Check if already exists
				already_exists = False
				if append:
					for existing in all_results["ivfpq"]:
						if (int(existing["N"]) == N and 
							int(existing["seed"]) == seed_val and
							int(existing["pq_kclusters"]) == param_set["pq_kclusters"] and
							int(existing["pq_nprobe"]) == param_set["pq_nprobe"] and
							int(existing["pq_M"]) == param_set["pq_M"]):
							already_exists = True
							break
				
				if already_exists:
					print(f"\n[{run_count}] IVF-PQ: N={N}, pq_kclusters={param_set['pq_kclusters']}, pq_nprobe={param_set['pq_nprobe']}, pq_M={param_set['pq_M']} - SKIPPING (already exists)")
					continue
				
				param_str = f"N{N}_pq_kclusters{param_set['pq_kclusters']}_pq_nprobe{param_set['pq_nprobe']}_pq_M{param_set['pq_M']}"
				output_txt = os.path.join(output_dir, f"ivfpq_{param_str}.txt")
				
				print(f"\n[{run_count}] IVF-PQ: N={N}, pq_kclusters={param_set['pq_kclusters']}, pq_nprobe={param_set['pq_nprobe']}, pq_M={param_set['pq_M']}")
				
				try:
					qps = run_protein_search(
						base_dat=base_dat,
						query_fasta=query_fasta,
						output_txt=output_txt,
						method="ivfpq",
						N=N,
						seed=seed_val,
						pq_kclusters=param_set['pq_kclusters'],
						pq_nprobe=param_set['pq_nprobe'],
						nbits=nbits,
						pq_M=param_set['pq_M'],
						R=0.5,
						range=False,
						k=0,
						L=0,
						lsh_w=0.0,
						kproj=0,
						hyper_w=0.0,
						hyper_M=0,
						probes=0,
						flat_kclusters=0,
						flat_nprobe=0,
						nlsh_index="",
						nlsh_T=0,
						nlsh_m=0,
						nlsh_imbalance=0.0,
						nlsh_kahip_mode=0,
						nlsh_layers=0,
						nlsh_nodes=0,
						nlsh_epochs=0,
						nlsh_batch_size=0,
						nlsh_lr=0.0
					)
					
					# Compute recall
					blast_tsv = os.path.join(blast_results_path, f"blast_results_top{N}.tsv")
					if os.path.exists(blast_tsv) and os.path.exists(output_txt):
						mean_recall, _ = calculate_recall(blast_tsv, output_txt, N)
						
						result = {
							"N": N,
							"seed": seed_val,
							"pq_kclusters": param_set['pq_kclusters'],
							"pq_nprobe": param_set['pq_nprobe'],
							"pq_M": param_set['pq_M'],
							"nbits": nbits,
							"recall": mean_recall,
							"qps": qps if qps else 0.0,
							"time_per_query": 1.0/qps if qps and qps > 0 else float('inf'),
						}
						all_results["ivfpq"].append(result)
						
						#print(f"  Recall@{N}: {mean_recall:.4f}, QPS: {qps:.2f if qps else 0}")
						# Chnged the baove line with the 2 above this sentence due to:
						# ERROR: Invalid Format Specifier
						qps_val = qps if qps else 0.0
						print(f"  Recall@{N}: {mean_recall:.4f}, QPS: {qps_val:.2f}")
					
					# Clean up
					if os.path.exists(output_txt):
						os.remove(output_txt)
					query_dat = output_txt.replace(".txt", ".queries.dat")
					query_ids = output_txt.replace(".txt", ".queries_ids.txt")
					if os.path.exists(query_dat):
						os.remove(query_dat)
					if os.path.exists(query_ids):
						os.remove(query_ids)
					
				except Exception as e:
					print(f"  ERROR: {e}")
					continue
	
	# NLSH grid search
	if "nlsh" in methods_to_run:
		print(f"\n{'='*80}")
		print("NLSH GRID SEARCH")
		print(f"{'='*80}")
		
		for N in N_vals:
			for param_set in get_nlsh_grid():
				run_count += 1
				
				# Check if already exists
				already_exists = False
				if append:
					for existing in all_results["nlsh"]:
						if (int(existing["N"]) == N and 
							int(existing["seed"]) == seed_val and
							int(existing["nlsh_T"]) == param_set["nlsh_T"] and
							int(existing["nlsh_m"]) == param_set["nlsh_m"] and
							int(existing["nlsh_layers"]) == param_set["nlsh_layers"] and
							int(existing["nlsh_nodes"]) == param_set["nlsh_nodes"]):
							already_exists = True
							break
				
				if already_exists:
					print(f"\n[{run_count}] NLSH: N={N}, T={param_set['nlsh_T']}, m={param_set['nlsh_m']}, layers={param_set['nlsh_layers']}, nodes={param_set['nlsh_nodes']} - SKIPPING (already exists)")
					continue
				
				param_str = f"N{N}_T{param_set['nlsh_T']}_m{param_set['nlsh_m']}_layers{param_set['nlsh_layers']}_nodes{param_set['nlsh_nodes']}"
				output_txt = os.path.join(output_dir, f"nlsh_{param_str}.txt")
				
				print(f"\n[{run_count}] NLSH: N={N}, T={param_set['nlsh_T']}, m={param_set['nlsh_m']}, layers={param_set['nlsh_layers']}, nodes={param_set['nlsh_nodes']}")
				
				try:
					qps = run_protein_search(
						base_dat=base_dat,
						query_fasta=query_fasta,
						output_txt=output_txt,
						method="nlsh",
						N=N,
						seed=seed_val,
						nlsh_T=param_set['nlsh_T'],
						nlsh_m=param_set['nlsh_m'],
						nlsh_layers=param_set['nlsh_layers'],
						nlsh_nodes=param_set['nlsh_nodes'],
						nlsh_epochs=nlsh_epochs,
						nlsh_batch_size=nlsh_batch_size,
						nlsh_lr=nlsh_lr,
						R=0.5,
						range=False,
						k=0,
						L=0,
						lsh_w=0.0,
						kproj=0,
						hyper_w=0.0,
						hyper_M=0,
						probes=0,
						flat_kclusters=0,
						flat_nprobe=0,
						pq_kclusters=0,
						pq_nprobe=0,
						pq_M=0,
						nbits=0,
						nlsh_index="",
						nlsh_imbalance=0.1,
						nlsh_kahip_mode=0,
					)

					# Compute recall
					blast_tsv = os.path.join(blast_results_path, f"blast_results_top{N}.tsv")
					if os.path.exists(blast_tsv) and os.path.exists(output_txt):
						mean_recall, _ = calculate_recall(blast_tsv, output_txt, N)
						
						result = {
							"N": N,
							"seed": seed_val,
							"nlsh_T": param_set['nlsh_T'],
							"nlsh_m": param_set['nlsh_m'],
							"nlsh_layers": param_set['nlsh_layers'],
							"nlsh_nodes": param_set['nlsh_nodes'],
							"nlsh_epochs": nlsh_epochs,
							"nlsh_batch_size": nlsh_batch_size,
							"nlsh_lr": nlsh_lr,
							"recall": mean_recall,
							"qps": qps if qps else 0.0,
							"time_per_query": 1.0/qps if qps and qps > 0 else float('inf'),
						}
						all_results["nlsh"].append(result)
						
						#print(f"  Recall@{N}: {mean_recall:.4f}, QPS: {qps:.2f if qps else 0}")
						# Chnged the baove line with the 2 above this sentence due to:
						# ERROR: Invalid Format Specifier
						qps_val = qps if qps else 0.0
						print(f"  Recall@{N}: {mean_recall:.4f}, QPS: {qps_val:.2f}")

					# Clean up
					if os.path.exists(output_txt):
						os.remove(output_txt)
					query_dat = output_txt.replace(".txt", ".queries.dat")
					query_ids = output_txt.replace(".txt", ".queries_ids.txt")
					if os.path.exists(query_dat):
						os.remove(query_dat)
					if os.path.exists(query_ids):
						os.remove(query_ids)
					
				except Exception as e:
					print(f"  ERROR: {e}")
					continue
	
	# Save results to CSV files
	print(f"\n{'='*80}")
	print("SAVING RESULTS TO CSV")
	print(f"{'='*80}")
	
	for method, results in all_results.items():
		if results:
			csv_path = os.path.join(output_dir, f"grid_search_{method}.csv")
			
			# Get fieldnames from first result
			fieldnames = list(results[0].keys())
			
			with open(csv_path, "w", newline="") as f:
				writer = csv.DictWriter(f, fieldnames=fieldnames)
				writer.writeheader()
				writer.writerows(results)
			
			print(f"Saved {len(results)} results to {csv_path}")


def main():#
	parser = argparse.ArgumentParser("Grid search for ANN hyperparameter tuning - V2 with append support")
	parser.add_argument("-d", required=True, help="Path to base protein vectors .dat")
	parser.add_argument("-q", required=True, help="Path to query FASTA")
	parser.add_argument("--output-dir", default="output/grid_search", help="Output directory for results")
	parser.add_argument("--method", default="all", choices=["lsh", "hypercube", "ivfflat", "ivfpq", "neural", "all"],
						help="Method to grid search")
	parser.add_argument("--seed", type=int, default=42, help="Random seed")
	parser.add_argument("--nbits", type=int, default=8, help="Bits per subspace (IVF-PQ)")
	parser.add_argument("--nlsh-epochs", type=int, default=8, help="NLSH training epochs")
	parser.add_argument("--nlsh-batch-size", type=int, default=512, help="NLSH batch size")
	parser.add_argument("--nlsh-lr", type=float, default=1e-3, help="NLSH learning rate")
	parser.add_argument("--no-append", action="store_true", help="Overwrite existing CSVs instead of appending")
	
	args = parser.parse_args()
	method = args.method
	method = method.lower()
	if method == "neural":
		method = "nlsh"
	
	run_grid_search(
		base_dat=args.d,
		query_fasta=args.q,
		output_dir=args.output_dir,
		method=method,
		seed=args.seed,
		nbits=args.nbits,
		nlsh_epochs=args.nlsh_epochs,
		nlsh_batch_size=args.nlsh_batch_size,
		nlsh_lr=args.nlsh_lr,
		append=not args.no_append,
	)


if __name__ == "__main__":
	main()
