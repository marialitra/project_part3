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