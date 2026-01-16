import libraries
from libraries import Optional, Dict, List, Tuple, defaultdict


# -----------------------------
# Main
# -----------------------------
def main():
	# -----------------------------
	# Parse Arguments
	# -----------------------------
	args = libraries.parse_args_search()

	method = args.method.lower()

	# We know for sure that:
	all_methods = ["lsh", "hypercube", "ivfflat", "ivfpq", "nlsh"]
	current_methods = ["lsh", "hypercube", "ivfflat", "ivfpq", "nlsh"]

	# The user will depict Neural LSH as neural
	# But we, throughout the code, we call it nlsh
	# For simplicity across the many programs
	# So let's change it!
	if method == "neural":
		method = "nlsh"
		
			
	# -----------------------------
	# Run Given Method
	# -----------------------------
	all_qps, answer = libraries.run_method(args, method, current_methods)

	# If method did not run we end the execution here!

	# But first if it is not all, we have already deleted uneccessary files
	if method != "all" and all_qps == None:
		print("\nRecall result (average):")
		libraries.print_recall(method, 0)
		return
	# If it is all, so we have NOT already deleted uneccessary file, we must do so, NOW
	elif method == "all" and all_qps == None:
		# Delete uneccessary files
		libraries.delete_files(args, method, all_methods)

		print("\nRecall result (average):")
		libraries.print_recall(method, 0)
		return


	# -----------------------------
	# Running BLAST command
	# -----------------------------
	query_ids, blast_results, blast_identity, blast_time_per_query, blast_qps = libraries.running_blast(args)

	# -----------------------------
	# Compute the Recall@N, using as ground truth the rsults BLAST returned
	# -----------------------------
	libraries.compute_recall(args, method, answer, blast_results, query_ids, all_qps, blast_identity, blast_time_per_query, blast_qps)

	# -----------------------------
	# Delete uneccessary files
	# -----------------------------
	libraries.delete_files(args, method, all_methods)

if __name__ == "__main__":
	main()