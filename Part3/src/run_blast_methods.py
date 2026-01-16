from libraries import os, time, subprocess, Dict, defaultdict

def blast_executable(neighbors: int, output: str):
    """
        Runs 'make blast' to run the BLAST method and its filtering.
    """

    try:
        start = time.perf_counter()
        build_process = subprocess.run(["make", "blast", f"N={neighbors}", f"out={output}"], capture_output=True, text=True, check=True)
        end = time.perf_counter()

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


def parse_blast_results_with_identity(blast_tsv_path: str) -> Dict[str, Dict[str, float]]:
    blast_data = defaultdict(dict)

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

            if target_field.startswith("sp|"):
                target_id = target_field.split("|")[1]
            else:
                continue

            try:
                identity = float(identity_str)
            except ValueError:
                continue

            prev = blast_data[query_id].get(target_id)

            if prev is None or identity > prev:
                blast_data[query_id][target_id] = identity

    return blast_data





def running_blast(args):
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
     
	return query_ids, blast_results, blast_identity, blast_time_per_query, blast_qps


