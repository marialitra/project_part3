emb:
	python3 src/protein_embed.py \
	-i Data/subset_1k.fasta \
	-o output/vectors.dat \
	-model esm2_t6_8M_UR50D

emb_50k:
	python3 src/protein_embed_danai.py \
	-i Data/swissprot_50k.fasta \
	-o output/vectors_50k_opt.dat \
	--batch-size 64 \
	-model esm2_t6_8M_UR50D

search:
	python src/protein_search.py \
	-d output/vectors.dat -q Data/targets.fasta \
	-o output/results.txt \
	-method ivfflat

# BLAST: build SwissProt DB and run blastp on targets
blast_db:
	@mkdir -p output/blast
	@command -v makeblastdb >/dev/null || (echo "makeblastdb not found. Install NCBI BLAST+ (apt: ncbi-blast+, conda: bioconda/blast) or add it to PATH." && exit 1)
	makeblastdb -in Data/swissprot_50k.fasta -dbtype prot -out output/blast/swissprot_db

blast_search: blast_db
	@command -v blastp >/dev/null || (echo "blastp not found. Install NCBI BLAST+ (apt: ncbi-blast+, conda: bioconda/blast) or add it to PATH." && exit 1)
	blastp -db output/blast/swissprot_db -query Data/targets.fasta -outfmt 6 -out output/blast/blast_results.tsv

.PHONY: blast_db blast_search

# Filter BLAST results (outfmt 6) to top N per query by bit score with E <= max-evalue
N = 5
EVALUE = 1e-3

blast_filter: blast_search
	python3 src/filter_blast.py -i output/blast/blast_results.tsv -o output/blast/blast_results_top${N}.tsv -n ${N} --max-evalue ${EVALUE}

.PHONY: blast_filter

# Possible Inputs:
# -method <all|lsh|hypercube|ivfflat|ivfpq>