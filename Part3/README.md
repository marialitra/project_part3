# project_part3

1.)

ALSO CHANGE THE FOLLOWING IN PART 1 POF THE PROJECT WHICH WE WILL UPLOAD IN GITHUB!!!!!
IMPORTANT CHANGE!!!!!!!!!!!!!!!

changes in ivfpq c source code, in ivfpq_init about
    // index->k = k_clusters; line 335

    and change it into:
    
    if(ivf_temp->k <= 0)
    {
        fprintf(stderr, "Error in kclusters\n");
        exit(1);
    }
    
    index->k = ivf_temp->k;
    k_clusters = index->k;

    lines 351 - 358

2.)

not called ids.txt but verctors_ids.txt, problem?

3.) in query.c of the first exercise we re entered the lines about the count of the QPS
    and we write it at the end of the file so that python part3 handles it as it wants to

4.) in nlsh_core.py from part2 we inserted the things about computing the QPS (there were there before)

5.) problem with M in pq and hypercube share the same in defaults so we shall see

6.) BLAST takes the same N as search will take

7.) it was running for nlsh but it is not an  option so i fixed it not to run for it

now lets do grid search to compute hyperparameters. N=[1,10,50], seed =42, kculsters = [200, 500, 1000, 5000], nprobe =[50, 100, 250, 500, 1500] only if nprobe < kclusters, k =[2,4,6], L = [5, 10, 15], w = [20, 40, 80, 120, 300, 500, 1000], kproj =[12, 14, 16, 20], M= [2000, 5000, 10000], probes=[100, 250, 500], nbits=8, ivfpq-M =[8,16], nlsh-T = [1000, 1500], nlsh-m=[1800, 2000, 2200], nlsh-layers =[5,10,15], nlsh-nodes = [128, 256, 512], epochs = 8, nlsh-batch-size = 512

8.) we thought that the extra command about the nueral lsh will be neural as it says in the pdf of the assignment (in the lab's pdf it does not say so)


9.)
BLAST may return fewer than N hits after filtering.

so instead of:
recall = neighbors_in_blast / args.N
we do:
recall = neighbors_in_blast / len(blast_top_n)


10.) IMPORTANT
Why the same protein appears multiple times in BLAST

Let’s take your example:

sp|P11961|ODP2_GEOSE


This is one unique UniProt protein.

Yet BLAST reports it three times, with different stats:

33% identity over 439 aa
38% identity over 78 aa
33% identity over 108 aa

Why?

Because BLAST does local alignment, not global alignment.

What BLAST actually reports (important concept)

BLAST reports HSPs (High-scoring Segment Pairs).

That means:

One protein sequence can align to multiple regions of the query

Each alignment is reported as a separate row

Each row can have:

Different start/end positions

Different alignment length

Different % identity

Different bit scores

Different E-values

This is 100% expected behavior.


“BLAST reports multiple local alignments per protein; we collapse hits by UniProt ID and evaluate recall on unique target proteins, which is standard practice.”
Also, in the blast results in the report of the output
we chose to take the results being shown from the best match of the protein in blast's list,
the one with the grater blast identity.
|-> “When multiple BLAST HSPs map to the same target protein, we retain the maximum percent 
    identity to represent the strongest local similarity.”


11.) we hardcoded the dimension of the data! as we will only take these as the assignment states

12.) use of cosine and eucleidian tell in report (use of eucleidian only in ivfpq)
“Cosine similarity was used for all embedding-based ANN methods except IVF-PQ. Product quantization relies on Euclidean geometry and residual encoding, which is incompatible with cosine similarity. Therefore, IVF-PQ was evaluated using Euclidean distance, consistent with standard practice. Neural LSH inherits the cosine similarity metric from its IVF-Flat search backend.”
13.) need of normalization in order to use the cosine distance metric
“All protein embeddings were L2-normalized prior to indexing so that cosine similarity could be approximated using Euclidean distance across all ANN methods.”


14.) normalization is happening inside of cosine, if you do normalization before cosine then you must only do the dot product, no need for cosine
    cos(a,b) = a dot b / norm(a) * norm(b)\

15.) if not method - all then if the algorithm did not run delete the results.txt if there was any
if it is method = all and one algorithm did not run then just skip it from the results .txt!!!
and if none of them runs when the method is all then we delete the requirements.txt