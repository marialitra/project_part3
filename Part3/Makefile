emb:
	python3 src/protein_embed.py \
	-i Data/subset_1k.fasta \
	-o output/embeddings/vectors.dat \
	-model esm2_t6_8M_UR50D

emb_50k:
	python3 src/protein_embed.py \
	-i Data/swissprot_50k.fasta \
	-o output/embeddings/vectors_50k_opt.dat \
	--batch-size 64 \
	-model esm2_t6_8M_UR50D

search:
	python3 src/protein_search.py \
	-d output/embeddings/vectors_50k_opt.dat -q Data/targets.fasta \
	-o output/search/results.txt \
	-method lsh
# Possible Inputs:
# -method <all|lsh|hypercube|ivfflat|ivfpq|neural>

# BLAST: build SwissProt DB and run blastp on targets
blast_db:
	@mkdir -p output/blast/database
	@mkdir -p output/blast/search
	@command -v makeblastdb >/dev/null || (echo "makeblastdb not found. Install NCBI BLAST+ (apt: ncbi-blast+, conda: bioconda/blast) or add it to PATH." && exit 1)
	makeblastdb -in Data/swissprot_50k.fasta -dbtype prot -out output/blast/database/swissprot_db

blast_search: blast_db
	@command -v blastp >/dev/null || (echo "blastp not found. Install NCBI BLAST+ (apt: ncbi-blast+, conda: bioconda/blast) or add it to PATH." && exit 1)
	blastp -db output/blast/database/swissprot_db -query Data/targets.fasta -outfmt 6 -out output/blast/search/blast_results.tsv

.PHONY: blast_db blast_search

# Filter BLAST results (outfmt 6) to top N per query by bit score with E <= max-evalue
N = 10
out = output/blast/topN/blast_results_top${N}.tsv
EVALUE = 1e-3 # Hardcoded it is the filter we choose to reduce noise

blast: blast_search
	python3 src/filter_blast.py -i output/blast/search/blast_results.tsv -o ${out} -n ${N} --max-evalue ${EVALUE}

.PHONY: blast_filter

# Grid search for hyperparameter tuning
# Usage: make grid_search [method=<all|lsh|hypercube|ivfflat|ivfpq|neural>]
method = all
grid_search:
	python3 src/grid_search.py \
	-d output/embeddings/vectors_50k_opt.dat \
	-q Data/targets.fasta \
	--output-dir output/grid_search \
	--method ${method}

.PHONY: grid_search
