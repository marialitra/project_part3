from libraries import Optional, List, Tuple, Dict

def generate_per_query_report(
	output_report: str,
	query_ids: list,
	topN: int,
	method_name: str,
	method_key: str,
	qps: Optional[float],
	method_results: Dict[str, List[Tuple[str, float]]],
	blast_results_topn: Dict[str, set],
	# blast_results_identity: Dict[str, List[Tuple[str, float]]],
	blast_results_identity: Dict[str, Dict[str, float]],
	blast_time_per_query: float,
	blast_qps: float
):
	"""
		Generate comprehensive per-query report with individual recall.
	"""

	with open(output_report, "w") as f:
		f.write("=" * 110 + "\n")
		f.write(f"ANN Evaluation Report (method = {method_name})\n")
		f.write("=" * 110 + "\n\n")

		# Iterate through each query
		for query_id in query_ids:
			f.write("\n" + "=" * 110 + "\n")
			f.write(f"Query: {query_id}\n")
			f.write(f"N = {topN}\n")
			f.write("=" * 110 + "\n\n")

			# Get neighbors for this query from results
			neighbors = method_results.get(query_id, [])
			blast_top_n = blast_results_topn.get(query_id, set())
			blast_identities = blast_results_identity.get(query_id, {})

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
				f.write(f"{method_name:<20} | {time_per_query:<15.3f} | {qps:<10.2f} | {recall:<25.2f}\n")
			else:
				f.write(f"{method_name:<20} | {'N/A':<15s} | {'N/A':<10s} | {recall:<25.2f}\n")

			f.write(f"{'BLAST (Ref)':<20} | {blast_time_per_query:<15.3f} | {blast_qps:<10.2f} | {1.00:<25.2f}\n")
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
				blast_id_str = f"{f'{blast_id_val:.0f}%' : <17}" if blast_id_val is not None else "undetected".ljust(17)

				# Bio comment logic
				if in_blast and blast_id_val is not None and blast_id_val > 30:
					bio_comment = "Homolog"
				elif in_blast and blast_id_val is not None and 20 < blast_id_val <= 30:
					bio_comment = "Remote homolog"
				elif in_blast and blast_id_val is not None and blast_id_val <= 20:
					bio_comment = "Weak similarity"
				elif not in_blast:
					bio_comment = "Possible false positive"
				else:
					bio_comment = ""

				f.write(f"{rank:<6} | {neighbor_id:<15} | {distance:<12.2f} | {blast_id_str} | {blast_in_str:<17} | {bio_comment:<30}\n")
				
	print(f"\n{'='*75}")
	print(f"[protein_search] Per-query report written to {output_report}")
	print(f"{'='*75}")


def generate_all_methods_report(
	output_report: str,
	query_ids: list,
	topN: int,
	method_qps: Dict[str, Optional[float]],
	method_results: Dict[str, Dict[str, List[Tuple[str, float]]]],
	methods: List[str],
	blast_results_topn: Dict[str, set],
	# blast_results_identity: Dict[str, List[Tuple[str, float]]],
	blast_results_identity: Dict[str, Dict[str, float]],
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

	with open(output_report, "w") as f:
		f.write("=" * 110 + "\n")
		f.write("ANN Evaluation Report (method = all)\n")
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

			for m in methods:
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
					f.write(f"{name:<20} | {time_per_query:<15.3f} | {qps_val:<10.2f} | {recall_q:<25.2f}\n")
				else:
					f.write(f"{name:<20} | {'N/A':<14s} | {'N/A':<10s} | {recall_q:<25.2f}\n")

			# BLAST reference row
			f.write(f"{'BLAST (Ref)':<20} | {blast_time_per_query:<15.3f} | {blast_qps:<10.2f} | {1.00:<25.2f}\n")
			f.write("-" * 110 + "\n\n")

			# [2] Top-N neighbors per method
			f.write("[2] Top-N neighbors per method\n")
			for m in methods:
				name = method_display[m]
				neighbors = method_results.get(m, {}).get(query_id, [])
				# blast_identities = {tid: ident for tid, ident in blast_results_identity.get(query_id, [])}
				blast_identities = blast_results_identity.get(query_id, {})

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
					elif in_blast and blast_id_val is not None and blast_id_val <= 20:
						bio_comment = "Weak similarity"
					elif not in_blast:
						bio_comment = "Possible false positive"
					else:
						bio_comment = ""

					f.write(f"{rank:<6} | {neighbor_id:<15} | {distance:<12.2f} | {blast_id_str} | {blast_in_str:<17} | {bio_comment:<30}\n")
					
	print(f"\n{'='*75}")
	print(f"[protein_search] Consolidated all-methods report written to {output_report}")
	print(f"{'='*75}")

