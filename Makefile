emb:
	python3 src/protein_embed_danai.py \
	-i Data/swissprot_50k.fasta \
	-o output/vectors.dat \
	-model esm2_t6_8M_UR50D